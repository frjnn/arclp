/* Compresses the surviving rows and columns into a CSC problem. */
#include "presolve.h"

#include <stdlib.h>

typedef struct
{
    int m, n, nnz;
} dimensions;

static dimensions count_live(const ps_matrix *a)
{
    dimensions d = {0, 0, 0};
    for (int i = 0; i < a->m; i++)
    {
        if (a->row_alive[i])
        {
            d.m++;
            d.nnz += a->row_len[i];
        }
    }
    for (int j = 0; j < a->n; j++)
    {
        d.n += a->col_alive[j];
    }
    return d;
}

static int allocate(arc_presolve *ps, dimensions d)
{
    const size_t doubles = (size_t)d.nnz + (size_t)d.m + (size_t)d.n;
    const size_t ints = (size_t)d.n + 1 + (size_t)d.nnz + (size_t)d.n;
    double *block = malloc(sizeof(double) * doubles + sizeof(int) * ints);
    if (block == NULL)
    {
        return -1;
    }
    ps->reduced_x = block;
    ps->reduced_b = ps->reduced_x + d.nnz;
    ps->reduced_c = ps->reduced_b + d.m;
    ps->reduced_p = (int *)(void *)(ps->reduced_c + d.n);
    ps->reduced_i = ps->reduced_p + d.n + 1;
    ps->reduced_col = ps->reduced_i + d.nnz;
    const arc_csc view = {d.m, d.n, ps->reduced_p, ps->reduced_i,
                          ps->reduced_x};
    ps->reduced = view;
    return 0;
}

/* New indices of surviving columns (others -1) and the reduced costs. */
static void number_columns(arc_presolve *ps, int *new_col)
{
    int k = 0;
    for (int j = 0; j < ps->a.n; j++)
    {
        new_col[j] = ps->a.col_alive[j] ? k : -1;
        if (ps->a.col_alive[j])
        {
            ps->reduced_col[k] = j;
            ps->reduced_c[k] = ps->c[j];
            ps->reduced_p[k + 1] = ps->a.col_len[j];
            k++;
        }
    }
    ps->reduced_p[0] = 0;
    for (k = 0; k < ps->reduced.n; k++)
    {
        ps->reduced_p[k + 1] += ps->reduced_p[k];
    }
}

/* Visiting rows in order leaves every column sorted by row. */
static void fill_entries(arc_presolve *ps, const int *new_col, int *next)
{
    const ps_matrix *a = &ps->a;
    for (int k = 0; k < ps->reduced.n; k++)
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
        for (int k = a->row_head[i]; k != -1; k = a->e[k].next_in_row)
        {
            if (ps_live(a, k))
            {
                const int slot = next[new_col[a->e[k].col]]++;
                ps->reduced_i[slot] = row;
                ps->reduced_x[slot] = a->e[k].value;
            }
        }
        ps->reduced_b[row++] = ps->b[i];
    }
}

int arc_extract_reduced(arc_presolve *ps)
{
    if (allocate(ps, count_live(&ps->a)) != 0)
    {
        return -1;
    }
    int *new_col = ps->scatter;      /* n entries, freed after this */
    int *next = ps->row_cols;
    number_columns(ps, new_col);
    fill_entries(ps, new_col, next);
    return 0;
}
