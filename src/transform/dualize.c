#include "dualize.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

static const double MARGIN = 0.5; /* how much smaller the dual must look */

static int length(const arc_csc *a, int j)
{
    return a->p[j + 1] - a->p[j];
}

/* A slack: plain, one entry. */
static int is_slack(const arc_lp *lp, int j)
{
    return j >= lp->n_free && j < arc_lp_first_bounded(lp)
        && length(&lp->a, j) == 1;
}

/* Rows with a slack are dense in no normal matrix; everything else is. */
int arc_dual_is_cheaper(const arc_lp *lp)
{
    const arc_csc *a = &lp->a;
    int *row_len = calloc((size_t)a->m + 1, sizeof(int));
    if (row_len == NULL)
    {
        return 0;
    }
    double primal = 0.0;
    double dual = 0.0;
    int kept = 0;
    for (int j = 0; j < a->n; j++)
    {
        const double len = length(a, j);
        primal += len * len;
        if (!is_slack(lp, j))
        {
            kept++;
            for (int k = a->p[j]; k < a->p[j + 1]; k++)
            {
                row_len[a->i[k]]++;
            }
        }
    }
    for (int i = 0; i < a->m; i++)
    {
        dual += (double)row_len[i] * row_len[i];
    }
    free(row_len);
    return kept > 0 && kept < a->m && dual + kept < MARGIN * primal;
}

/* The dual's sizes, known before any array exists. */
static void measure(arc_dual *d)
{
    const arc_lp *lp = d->primal;
    d->cols = lp->a.m;
    for (int j = 0; j < lp->a.n; j++)
    {
        if (!is_slack(lp, j))
        {
            const int columns = (j >= lp->n_free)
                              + (j >= arc_lp_first_bounded(lp));
            d->n_kept++;
            d->cols += columns;
            d->nnz += length(&lp->a, j) + columns;
        }
    }
}

static void layout(arc_arena *ar, void *owner)
{
    arc_dual *d = owner;
    const size_t m = (size_t)d->primal->a.m;
    const size_t cols = (size_t)d->cols;
    const size_t nnz = (size_t)d->nnz;
    d->row_of = arc_arena_take(ar, (size_t)d->primal->a.n, sizeof(int));
    d->kept = arc_arena_take(ar, (size_t)d->n_kept, sizeof(int));
    double **by_column[] = {&d->c, &d->lower, &d->upper};
    arc_arena_take_vectors(ar, by_column, 3, cols);
    d->b = arc_arena_take(ar, (size_t)d->n_kept, sizeof(double));
    d->x = arc_arena_take(ar, nnz, sizeof(double));
    d->left = arc_arena_take(ar, m, sizeof(double));
    d->p = arc_arena_take(ar, cols + 1, sizeof(int));
    d->i = arc_arena_take(ar, nnz, sizeof(int));
    d->pos = arc_arena_take(ar, nnz, sizeof(int));
    d->next = arc_arena_take(ar, m, sizeof(int));
}

static void number_kept(arc_dual *d)
{
    int t = 0;
    for (int j = 0; j < d->primal->a.n; j++)
    {
        d->row_of[j] = -1;
        if (!is_slack(d->primal, j))
        {
            d->row_of[j] = t;
            d->kept[t++] = j;
        }
    }
}

/* The bound a slack puts on y_i: a y_i <= c. */
static void bound_by_slack(arc_dual *d, int j)
{
    const arc_lp *lp = d->primal;
    const int i = lp->a.i[lp->a.p[j]];
    const double value = lp->a.x[lp->a.p[j]];
    const double limit = lp->c[j] / value;
    if (value > 0.0)
    {
        d->upper[i] = fmin(d->upper[i], limit);
    }
    else
    {
        d->lower[i] = fmax(d->lower[i], limit);
    }
}

/* Columns y_i of the dual: the kept entries of row i of A. */
static void fill_y(arc_dual *d)
{
    const arc_lp *lp = d->primal;
    const int m = lp->a.m;
    arc_csc_index_rows(&lp->a, NULL, d->row_of, d->p, d->i, d->pos, d->next);
    for (int k = 0; k < d->p[m]; k++)
    {
        d->i[k] = d->row_of[d->i[k]];
        d->x[k] = lp->a.x[d->pos[k]];
    }
    for (int i = 0; i < m; i++)
    {
        d->c[i] = -lp->b[i];
        d->lower[i] = -INFINITY;
        d->upper[i] = INFINITY;
    }
}

static void push_unit(arc_dual *d, int *col, int row, double value,
                      double cost)
{
    const int slot = d->p[*col];
    d->i[slot] = row;
    d->x[slot] = value;
    d->c[*col] = cost;
    d->lower[*col] = 0.0;
    d->upper[*col] = INFINITY;
    d->p[++*col] = slot + 1;
}

/* Columns s (one per kept column with a sign constraint), then z. */
static void fill_slacks(arc_dual *d)
{
    const arc_lp *lp = d->primal;
    const int first = arc_lp_first_bounded(lp);
    int col = lp->a.m;
    for (int t = 0; t < d->n_kept; t++)
    {
        d->b[t] = lp->c[d->kept[t]];
        if (d->kept[t] >= lp->n_free)
        {
            push_unit(d, &col, t, 1.0, 0.0);
        }
    }
    for (int t = 0; t < d->n_kept; t++)
    {
        if (d->kept[t] >= first)
        {
            push_unit(d, &col, t, -1.0, lp->upper[d->kept[t] - first]);
        }
    }
}

int arc_dual_create(const arc_lp *lp, arc_dual *d)
{
    memset(d, 0, sizeof(*d));
    d->primal = lp;
    measure(d);
    if (arc_arena_build(&d->arena, layout, d) != 0)
    {
        return -1;
    }
    const arc_csc view = {d->n_kept, d->cols, d->p, d->i, d->x};
    d->a = view;
    number_kept(d);
    fill_y(d);
    fill_slacks(d);
    for (int j = 0; j < lp->a.n; j++)
    {
        if (d->row_of[j] < 0)
        {
            bound_by_slack(d, j);
        }
    }
    return 0;
}

void arc_dual_free(arc_dual *d)
{
    arc_arena_release(&d->arena);
}

/* What the kept columns leave of row i goes to the slack of that row with
 * the matching sign and the tightest bound on y_i. */
static void recover_slacks(const arc_dual *d, double *x)
{
    const arc_lp *lp = d->primal;
    double *left = d->left;
    for (int j = 0; j < lp->a.n; j++)
    {
        if (d->row_of[j] >= 0)
        {
            continue;
        }
        const int k = lp->a.p[j];
        const int i = lp->a.i[k];
        const double limit = lp->c[j] / lp->a.x[k];
        const int tightest = lp->a.x[k] > 0.0 ? limit <= d->upper[i]
                                              : limit >= d->lower[i];
        const double value = left[i] / lp->a.x[k];
        x[j] = tightest && value > 0.0 ? value : 0.0;
        if (x[j] > 0.0)
        {
            left[i] = 0.0;
        }
    }
}

void arc_dual_recover(const arc_dual *d, const double *multipliers,
                      double *x)
{
    const arc_lp *lp = d->primal;
    const int first = arc_lp_first_bounded(lp);
    memcpy(d->left, lp->b, sizeof(double) * (size_t)lp->a.m);
    for (int t = 0; t < d->n_kept; t++)
    {
        const int j = d->kept[t];
        double value = -multipliers[t];
        if (j >= lp->n_free)
        {
            value = fmax(value, 0.0);
        }
        if (j >= first)
        {
            value = fmin(value, lp->upper[j - first]);
        }
        x[j] = value;
        for (int k = lp->a.p[j]; k < lp->a.p[j + 1]; k++)
        {
            d->left[lp->a.i[k]] -= lp->a.x[k] * value;
        }
    }
    recover_slacks(d, x);
}
