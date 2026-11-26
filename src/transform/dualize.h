/*
 * The dual of an internal LP as another LP.  A problem with many rows and
 * few columns is cheaper to solve as its dual, whose normal matrix has one
 * row per column.  With A x = b, 0 <= x <= u (free: no bound):
 *
 *     min  -b'y + u'z   s.t.   A_j'y + s_j - z_j = c_j   (kept columns)
 *
 * Plain singleton columns (slacks of inequalities) are not kept: they only
 * bound y, so their rows disappear.
 */
#ifndef ARC_DUALIZE_H
#define ARC_DUALIZE_H

#include "sparse/lp.h"
#include "util/arena.h"

typedef struct
{
    const arc_lp *primal;
    arc_csc a;            /* n_kept rows; columns y | s | z           */
    double *b, *c, *lower, *upper;
    int *kept;            /* primal column of each dual row           */
    int n_kept;
    int *row_of;          /* dual row of each primal column, -1 slack */
    int *p, *i;
    double *x;            /* writable storage behind a                */
    double *left;         /* primal rows: what the kept columns leave */
    int *pos, *next;      /* scratch                                  */
    int cols, nnz;        /* sizes, known before the arrays           */
    arc_arena arena;
} arc_dual;

/* 1 when the dual has the smaller normal matrix, judged by the sums of
 * squared column lengths that bound the two sparsity patterns. */
int arc_dual_is_cheaper(const arc_lp *lp);

int arc_dual_create(const arc_lp *lp, arc_dual *d);
void arc_dual_free(arc_dual *d);

/* Primal solution from the multipliers of the dual's rows. */
void arc_dual_recover(const arc_dual *d, const double *multipliers,
                      double *x);

#endif
