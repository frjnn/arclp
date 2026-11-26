/*
 * Refinement of normal-equation solves by conjugate gradients on
 * A diag(d) A', preconditioned with the Cholesky factorization.  One or two
 * steps recover the accuracy that replaced pivots and ill-conditioned
 * scalings cost the factorization late in a solve.  It is switched on only
 * when needed: checking the residual alone costs two products with A.
 */
#ifndef ARC_REFINE_H
#define ARC_REFINE_H

#include "numeric.h"
#include "util/arena.h"

typedef struct
{
    const arc_csc *a;
    const double *d;    /* scaling of the last factorization */
    int enabled;        /* 0: only the plain solve           */
    double *rhs;        /* m x ARC_MAX_RHS copy of the right-hand sides */
    double *r, *p, *q, *z; /* m each */
    double *best;          /* m, iterate with the smallest residual */
    double *t;          /* n         */
    arc_arena arena;
} arc_refine;

int arc_refine_create(arc_refine *f, const arc_csc *a);
void arc_refine_free(arc_refine *f);

/* Overwrites r (m x nrhs) with the refined solution of A D A' x = r. */
void arc_refine_solve(arc_refine *f, arc_numeric *num, double *r, int nrhs);

#endif
