/*
 * Linear program in the solver's internal form:
 *     min c'x  s.t.  A x = b,
 * with the columns ordered as  free | nonnegative | bounded, the last
 * n_bounded ones satisfying 0 <= x <= upper.
 */
#ifndef ARC_LP_H
#define ARC_LP_H

#include "csc.h"

typedef struct
{
    arc_csc a;
    double *values; /* writable storage behind a.x    */
    double *b, *c;
    double *upper;  /* n_bounded entries              */
    int n_free, n_bounded;
} arc_lp;

static inline int arc_lp_first_bounded(const arc_lp *lp)
{
    return lp->a.n - lp->n_bounded;
}

#endif
