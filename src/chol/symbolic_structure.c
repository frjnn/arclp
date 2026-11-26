#include "symbolic_steps.h"

#include <stdlib.h>

static int compare_ints(const void *a, const void *b)
{
    const int x = *(const int *)a;
    const int y = *(const int *)b;
    return (x > y) - (x < y);
}

static void prefix_offsets(arc_symbolic *s, const int *srows)
{
    s->row_ptr[0] = 0;
    s->val_ptr[0] = 0;
    s->max_rows = 0;
    for (int k = 0; k < s->nsuper; k++)
    {
        s->row_ptr[k + 1] = s->row_ptr[k] + srows[k];
        s->val_ptr[k + 1] = s->val_ptr[k] + (long)srows[k] * arc_sn_cols(s, k);
        s->max_rows = srows[k] > s->max_rows ? srows[k] : s->max_rows;
        for (int j = s->first[k]; j < s->first[k + 1]; j++)
        {
            s->owner[j] = k;
        }
    }
}

/* Lists the supernodal children of every supernode. */
static void link_children(const arc_symbolic *s, const int *parent,
                          int *head, int *next)
{
    for (int k = 0; k < s->nsuper; k++)
    {
        head[k] = -1;
    }
    for (int k = s->nsuper - 1; k >= 0; k--)
    {
        const int up = parent[s->first[k + 1] - 1];
        if (up != -1)
        {
            next[k] = head[s->owner[up]];
            head[s->owner[up]] = k;
        }
    }
}

typedef struct
{
    int *out;
    int count;
    int capacity;
    int last_col;
    int stamp;
    int *mark;
} row_set;

static int add_row(row_set *set, int i)
{
    if (i <= set->last_col || set->mark[i] == set->stamp)
    {
        return 0;
    }
    if (set->count == set->capacity)
    {
        return -1;
    }
    set->mark[i] = set->stamp;
    set->out[set->count++] = i;
    return 0;
}

/* Rows of supernode k: its columns, the entries of P A A' P' below them,
 * and the off-diagonal rows of its children. */
static int collect_rows(const arc_graph *g, arc_symbolic *s, int k,
                        const int *head, const int *next, int *mark)
{
    const int f = s->first[k];
    row_set set = {s->rows + s->row_ptr[k], 0, arc_sn_rows(s, k),
                   s->first[k + 1] - 1, k, mark};
    for (int j = f; j <= set.last_col; j++)
    {
        set.out[set.count++] = j;
        mark[j] = k;
    }
    int status = 0;
    for (int j = f; j <= set.last_col; j++)
    {
        const int v = s->perm[j];
        for (int t = g->p[v]; t < g->p[v + 1]; t++)
        {
            status |= add_row(&set, s->pinv[g->i[t]]);
        }
    }
    for (int c = head[k]; c != -1; c = next[c])
    {
        const int *crows = s->rows + s->row_ptr[c];
        for (int r = arc_sn_cols(s, c); r < arc_sn_rows(s, c); r++)
        {
            status |= add_row(&set, crows[r]);
        }
    }
    const int ncols = arc_sn_cols(s, k);
    qsort(set.out + ncols, (size_t)(set.count - ncols), sizeof(int),
          compare_ints);
    return (status != 0 || set.count != set.capacity) ? -1 : 0;
}

/* Largest block L_d(p1:end, :) * L_d(p1:p2, :)' over all updates. */
static void measure_updates(arc_symbolic *s)
{
    s->max_update = 0;
    for (int d = 0; d < s->nsuper; d++)
    {
        const int *rows = s->rows + s->row_ptr[d];
        const int nrows = arc_sn_rows(s, d);
        for (int p1 = arc_sn_cols(s, d), p2; p1 < nrows; p1 = p2)
        {
            const int end = s->first[s->owner[rows[p1]] + 1];
            for (p2 = p1; p2 < nrows && rows[p2] < end; p2++)
            {
            }
            const long size = (long)(p2 - p1) * (nrows - p1);
            s->max_update = size > s->max_update ? size : s->max_update;
        }
    }
}

int arc_build_structure(const arc_graph *g, const int *parent,
                        const int *srows, arc_symbolic *s, int *work)
{
    int *mark = work;
    int *head = work + s->m;
    int *next = work + 2 * s->m;
    prefix_offsets(s, srows);
    link_children(s, parent, head, next);
    for (int i = 0; i < s->m; i++)
    {
        mark[i] = -1;
    }
    for (int k = 0; k < s->nsuper; k++)
    {
        if (collect_rows(g, s, k, head, next, mark) != 0)
        {
            return -1;
        }
    }
    measure_updates(s);
    return 0;
}
