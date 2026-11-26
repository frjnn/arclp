/* Compresses the surviving rows and columns into the solver's LP form. */
#include "presolve.h"

#include <math.h>

/* 0 free, 1 plain, 2 bounded: the column order of arc_lp. */
static int kind(const arc_presolve *ps, int j)
{
    return ps->is_free[j] ? 0 : isinf(ps->upper[j]) ? 1 : 2;
}

static void measure(arc_presolve *ps)
{
    const ps_matrix *a = &ps->a;
    arc_lp *lp = &ps->reduced;
    lp->a.m = lp->a.n = 0;
    lp->n_free = lp->n_bounded = 0;
    ps->reduced_nnz = 0;
    for (int i = 0; i < a->m; i++)
    {
        if (a->row_alive[i])
        {
            lp->a.m++;
            ps->reduced_nnz += a->row_len[i];
        }
    }
    for (int j = 0; j < a->n; j++)
    {
        if (a->col_alive[j])
        {
            lp->a.n++;
            lp->n_free += kind(ps, j) == 0;
            lp->n_bounded += kind(ps, j) == 2;
        }
    }
}

static void layout(arc_arena *ar, void *owner)
{
    arc_presolve *ps = owner;
    arc_lp *lp = &ps->reduced;
    const size_t nnz = (size_t)ps->reduced_nnz;
    const size_t n = (size_t)lp->a.n;
    lp->values = arc_arena_take(ar, nnz, sizeof(double));
    lp->b = arc_arena_take(ar, (size_t)lp->a.m, sizeof(double));
    lp->c = arc_arena_take(ar, n, sizeof(double));
    lp->upper = arc_arena_take(ar, (size_t)lp->n_bounded, sizeof(double));
    ps->reduced_p = arc_arena_take(ar, n + 1, sizeof(int));
    ps->reduced_i = arc_arena_take(ar, nnz, sizeof(int));
    ps->reduced_col = arc_arena_take(ar, n, sizeof(int));
    const arc_csc view = {lp->a.m, lp->a.n, ps->reduced_p, ps->reduced_i,
                          lp->values};
    lp->a = view;
}

static void place(arc_presolve *ps, int j, int k)
{
    const int first = arc_lp_first_bounded(&ps->reduced);
    ps->reduced_col[k] = j;
    ps->reduced.c[k] = ps->c[j];
    ps->reduced_p[k + 1] = ps->a.col_len[j];
    if (k >= first)
    {
        ps->reduced.upper[k - first] = ps->upper[j];
    }
}

/* New indices of surviving columns (others -1) and their data. */
static void number_columns(arc_presolve *ps, int *new_col)
{
    int next[3] = {0, ps->reduced.n_free, arc_lp_first_bounded(&ps->reduced)};
    for (int j = 0; j < ps->a.n; j++)
    {
        new_col[j] = -1;
        if (ps->a.col_alive[j])
        {
            new_col[j] = next[kind(ps, j)]++;
            place(ps, j, new_col[j]);
        }
    }
    ps->reduced_p[0] = 0;
    for (int k = 0; k < ps->reduced.a.n; k++)
    {
        ps->reduced_p[k + 1] += ps->reduced_p[k];
    }
}

/* Visiting rows in order leaves every column sorted by row. */
static void fill_entries(arc_presolve *ps, const int *new_col, int *next)
{
    const ps_matrix *a = &ps->a;
    for (int k = 0; k < ps->reduced.a.n; k++)
    {
        next[k] = ps->reduced_p[k];
    }
    int row = 0;
    for (int i = 0; i < a->m; i++)
    {
        if (!a->row_alive[i])
        {
            continue;
        }
        PS_FOR_ROW(a, k, i)
        {
            const int slot = next[new_col[a->e[k].col]]++;
            ps->reduced_i[slot] = row;
            ps->reduced.values[slot] = a->e[k].value;
        }
        ps->reduced.b[row++] = ps->b[i];
    }
}

int arc_extract_reduced(arc_presolve *ps)
{
    measure(ps);
    if (arc_arena_build(&ps->arena, layout, ps) != 0)
    {
        return -1;
    }
    int *new_col = ps->scatter; /* n entries, freed after this */
    int *next = ps->row_cols;
    number_columns(ps, new_col);
    fill_entries(ps, new_col, next);
    return 0;
}
