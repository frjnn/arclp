/*
 * Sparse matrix for presolve: every entry is linked into its row and its
 * column, so whole rows and columns can be deleted in O(length) and fill-in
 * can be inserted.  Deleted rows/columns are flagged rather than unlinked;
 * an entry is live while its value is nonzero and its row and column live.
 */
#ifndef ARC_PS_MATRIX_H
#define ARC_PS_MATRIX_H

#include "sparse/csc.h"

typedef struct
{
    int row, col;
    int next_in_row, next_in_col;
    double value;
} ps_entry;

typedef struct
{
    int m, n;
    ps_entry *e;
    int count, capacity;
    int *row_head, *col_head;
    int *row_len, *col_len; /* live entries per row / column */
    unsigned char *row_alive, *col_alive;
} ps_matrix;

int ps_matrix_init(ps_matrix *a, const arc_csc *src);
void ps_matrix_free(ps_matrix *a);

/* Inserts a new entry (fill-in); returns 0 on success. */
int ps_matrix_insert(ps_matrix *a, int row, int col, double value);

/* Changes a value, keeping the live counts consistent. */
void ps_matrix_set(ps_matrix *a, int k, double value);

void ps_matrix_drop_row(ps_matrix *a, int i);
void ps_matrix_drop_col(ps_matrix *a, int j);

static inline int ps_live(const ps_matrix *a, int k)
{
    const ps_entry *e = &a->e[k];
    return e->value != 0.0 && a->row_alive[e->row] && a->col_alive[e->col];
}

/* First live entry at or after k along a row (column) list, -1 at the end. */
static inline int ps_live_in_row(const ps_matrix *a, int k)
{
    while (k != -1 && !ps_live(a, k))
    {
        k = a->e[k].next_in_row;
    }
    return k;
}

static inline int ps_live_in_col(const ps_matrix *a, int k)
{
    while (k != -1 && !ps_live(a, k))
    {
        k = a->e[k].next_in_col;
    }
    return k;
}

/* Visits the live entries of a row / column.  The body may delete the rows
 * and columns it meets. */
#define PS_FOR_ROW(a, k, i)                                                 \
    for (int k = ps_live_in_row(a, (a)->row_head[i]); k != -1;              \
         k = ps_live_in_row(a, (a)->e[k].next_in_row))
#define PS_FOR_COL(a, k, j)                                                 \
    for (int k = ps_live_in_col(a, (a)->col_head[j]); k != -1;              \
         k = ps_live_in_col(a, (a)->e[k].next_in_col))

#endif
