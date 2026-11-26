#include "ps_matrix.h"

#include "util/grow.h"

#include <string.h>

static int link(ps_matrix *a, int row, int col, double value)
{
    if (a->count == a->capacity)
    {
        const int capacity = arc_grown_capacity(a->capacity, a->count + 1, 64);
        if (arc_resize(&a->e, (size_t)capacity, sizeof(*a->e)) != 0)
        {
            return -1;
        }
        a->capacity = capacity;
    }
    const int k = a->count++;
    const ps_entry entry = {row, col, a->row_head[row], a->col_head[col],
                            value};
    a->e[k] = entry;
    a->row_head[row] = k;
    a->col_head[col] = k;
    a->row_len[row]++;
    a->col_len[col]++;
    return 0;
}

int ps_matrix_init(ps_matrix *a, const arc_csc *src)
{
    const size_t m = (size_t)src->m;
    const size_t n = (size_t)src->n;
    memset(a, 0, sizeof(*a));
    a->m = src->m;
    a->n = src->n;
    a->row_head = malloc(sizeof(int) * (2 * m + 2 * n + 1));
    a->row_alive = malloc(m + n + 1);
    a->capacity = src->p[src->n] + src->p[src->n] / 4 + 64;
    a->e = malloc(sizeof(ps_entry) * (size_t)a->capacity);
    if (a->row_head == NULL || a->row_alive == NULL || a->e == NULL)
    {
        ps_matrix_free(a);
        return -1;
    }
    a->col_head = a->row_head + m;
    a->row_len = a->col_head + n;
    a->col_len = a->row_len + m;
    a->col_alive = a->row_alive + m;
    memset(a->row_head, 0xff, sizeof(int) * (m + n));
    memset(a->row_len, 0, sizeof(int) * (m + n));
    memset(a->row_alive, 1, m + n);
    /* Columns are linked in reverse so every list ends up in index order. */
    for (int j = src->n - 1; j >= 0; j--)
    {
        for (int t = src->p[j + 1] - 1; t >= src->p[j]; t--)
        {
            if (src->x[t] != 0.0 && link(a, src->i[t], j, src->x[t]) != 0)
            {
                ps_matrix_free(a);
                return -1;
            }
        }
    }
    return 0;
}

void ps_matrix_free(ps_matrix *a)
{
    free(a->e);
    free(a->row_head);
    free(a->row_alive);
    memset(a, 0, sizeof(*a));
}

int ps_matrix_insert(ps_matrix *a, int row, int col, double value)
{
    return link(a, row, col, value);
}

void ps_matrix_set(ps_matrix *a, int k, double value)
{
    ps_entry *e = &a->e[k];
    const int delta = (value != 0.0) - (e->value != 0.0);
    a->row_len[e->row] += delta;
    a->col_len[e->col] += delta;
    e->value = value;
}

void ps_matrix_drop_row(ps_matrix *a, int i)
{
    PS_FOR_ROW(a, k, i)
    {
        a->col_len[a->e[k].col]--;
    }
    a->row_alive[i] = 0;
}

void ps_matrix_drop_col(ps_matrix *a, int j)
{
    PS_FOR_COL(a, k, j)
    {
        a->row_len[a->e[k].row]--;
    }
    a->col_alive[j] = 0;
}
