/* How far an iterate is from optimal, and proofs that no optimum exists. */
#ifndef ARC_CONVERGENCE_H
#define ARC_CONVERGENCE_H

#include "ipm.h"

typedef struct
{
    double primal; /* c'x        */
    double dual;   /* b'y - u'z  */
} arc_ipm_objectives;

arc_ipm_objectives arc_ipm_objective(const arc_ipm *ipm);

/* Measures errors in the original units of a scaled problem (NULL: not
 * scaled): 1 / row scaling, 1 / column scaling and the column scaling. */
void arc_ipm_set_units(arc_ipm *ipm, const double *row_unit,
                       const double *col_unit, const double *col);

/* Relative residuals and duality gap of the current iterate. */
arc_ipm_errors arc_ipm_measure(const arc_ipm *ipm);

static inline double arc_ipm_worst(arc_ipm_errors e)
{
    const double residual = e.primal > e.dual ? e.primal : e.dual;
    return residual > e.gap ? residual : e.gap;
}

/* Complementarity x's + w'z relative to the objective: unlike the duality
 * gap, it cannot pass through zero far from the optimum. */
double arc_ipm_complementarity(const arc_ipm *ipm);

/* ARC_IPM_PRIMAL_INFEASIBLE or ARC_IPM_UNBOUNDED when the iterate proves
 * it, ARC_IPM_STALLED when it proves nothing. */
arc_ipm_outcome arc_ipm_certified(const arc_ipm *ipm);

#endif
