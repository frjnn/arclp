/*
 * Symbolic supernodal Cholesky factorization of P (A A') P'.
 *
 * Supernode s owns the columns first[s] .. first[s+1]-1 of L and stores them
 * as one dense column-major block of nrows(s) x ncols(s) entries, whose row
 * labels are rows[row_ptr[s] .. row_ptr[s+1]): the supernode's own columns
 * first, then the remaining rows in increasing order.
 */
#ifndef ARC_SYMBOLIC_H
#define ARC_SYMBOLIC_H

#include "sparse/csc.h"
#include "util/arena.h"
#include "supernodes.h"

typedef struct
{
    int m;
    int *perm;          /* perm[k] = row of A eliminated k-th            */
    int *pinv;
    int nsuper;
    int *first;         /* nsuper + 1                                    */
    int *owner;         /* supernode owning each column                  */
    int *row_ptr;       /* nsuper + 1                                    */
    int *rows;
    long *val_ptr;      /* nsuper + 1 offsets of the dense blocks        */
    int max_rows;       /* largest nrows(s)                              */
    long max_update;    /* largest descendant update block               */
    /* A with rows relabelled by pinv, each column sorted by label, and its
     * row-wise index: (column, position in a_perm) of each entry.         */
    arc_csc a_perm;
    int *ap_p, *ap_i;   /* writable storage behind a_perm */
    double *ap_x;
    int *by_row_ptr;    /* m + 1 */
    int *by_row_col;
    int *by_row_pos;
    arc_arena arena;
} arc_symbolic;

int arc_symbolic_analyse(const arc_csc *a, const arc_relax_policy *policy,
                         arc_symbolic *sym);

void arc_symbolic_free(arc_symbolic *sym);

static inline int arc_sn_cols(const arc_symbolic *s, int k)
{
    return s->first[k + 1] - s->first[k];
}

static inline int arc_sn_rows(const arc_symbolic *s, int k)
{
    return s->row_ptr[k + 1] - s->row_ptr[k];
}

#endif
