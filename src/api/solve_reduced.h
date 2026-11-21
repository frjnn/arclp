/* Interior-point solve of the presolved problem (use-case layer). */
#ifndef ARC_SOLVE_REDUCED_H
#define ARC_SOLVE_REDUCED_H

#include <arclp/arclp.h>

#include "sparse/csc.h"

/* Solves min c'x, A x = b, x >= 0 for the reduced problem; x has a->n
 * entries.  Fills the iteration, factorization and timing fields of info. */
arclp_status arc_solve_reduced(const arc_csc *a, const double *b,
                               const double *c, const arclp_options *options,
                               double *x, arclp_info *info);

#endif
