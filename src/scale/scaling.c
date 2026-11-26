#include "scaling.h"

#include <math.h>
#include <string.h>

enum
{
    PASSES = 8
};

static void largest_entries(const arc_csc *a, const double *values,
                            double *row_max, double *col_max)
{
    memset(row_max, 0, sizeof(double) * (size_t)a->m);
    for (int j = 0; j < a->n; j++)
    {
        col_max[j] = 0.0;
        for (int k = a->p[j]; k < a->p[j + 1]; k++)
        {
            const double v = fabs(values[k]);
            col_max[j] = v > col_max[j] ? v : col_max[j];
            row_max[a->i[k]] = v > row_max[a->i[k]] ? v : row_max[a->i[k]];
        }
    }
}

static double inverse_root(double largest)
{
    return largest > 0.0 ? 1.0 / sqrt(largest) : 1.0;
}

/* One Ruiz pass: divide every row and column by the square root of its
 * largest absolute entry. */
static void equilibrate(const arc_csc *a, double *values, arc_scaling *s)
{
    largest_entries(a, values, s->row_max, s->col_max);
    for (int i = 0; i < a->m; i++)
    {
        s->row_max[i] = inverse_root(s->row_max[i]);
        s->row[i] *= s->row_max[i];
    }
    for (int j = 0; j < a->n; j++)
    {
        const double f = inverse_root(s->col_max[j]);
        s->col[j] *= f;
        for (int k = a->p[j]; k < a->p[j + 1]; k++)
        {
            values[k] *= f * s->row_max[a->i[k]];
        }
    }
}

static void layout(arc_arena *ar, void *owner)
{
    arc_scaling *s = owner;
    double **by_row[] = {&s->row, &s->row_unit, &s->row_max};
    double **by_col[] = {&s->col, &s->col_unit, &s->col_max};
    arc_arena_take_vectors(ar, by_row, 3, (size_t)s->m);
    arc_arena_take_vectors(ar, by_col, 3, (size_t)s->n);
}

int arc_scale_problem(arc_lp *lp, arc_scaling *s)
{
    const arc_csc *a = &lp->a;
    const size_t m = (size_t)a->m;
    const size_t n = (size_t)a->n;
    s->m = a->m;
    s->n = a->n;
    if (arc_arena_build(&s->arena, layout, s) != 0)
    {
        return -1;
    }
    for (size_t i = 0; i < m; i++)
    {
        s->row[i] = 1.0;
    }
    for (size_t j = 0; j < n; j++)
    {
        s->col[j] = 1.0;
    }
    for (int pass = 0; pass < PASSES; pass++)
    {
        equilibrate(a, lp->values, s);
    }
    for (size_t i = 0; i < m; i++)
    {
        lp->b[i] *= s->row[i];
        s->row_unit[i] = 1.0 / s->row[i];
    }
    for (size_t j = 0; j < n; j++)
    {
        lp->c[j] *= s->col[j];
        s->col_unit[j] = 1.0 / s->col[j];
    }
    const int first = arc_lp_first_bounded(lp);
    for (int k = 0; k < lp->n_bounded; k++)
    {
        lp->upper[k] *= s->col_unit[first + k];
    }
    return 0;
}

void arc_unscale(const arc_scaling *s, const arc_lp *lp, double *x,
                 double *y)
{
    for (int j = 0; j < lp->a.n; j++)
    {
        x[j] *= s->col[j];
    }
    for (int i = 0; y != NULL && i < lp->a.m; i++)
    {
        y[i] *= s->row[i];
    }
}

void arc_scaling_free(arc_scaling *s)
{
    arc_arena_release(&s->arena);
}
