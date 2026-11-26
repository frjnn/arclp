/* Interior-point solve of the presolved problem (use-case layer). */
#ifndef ARC_SOLVE_REDUCED_H
#define ARC_SOLVE_REDUCED_H

#include <arclp/arclp.h>

#include "sparse/lp.h"

/* Solves the reduced problem, which is scaled in place; x has a.n entries.
 * Fills the iteration, factorization and timing fields of info. */
arclp_status arc_solve_reduced(arc_lp *lp, const arclp_options *options,
                               double *x, arclp_info *info);

#endif
