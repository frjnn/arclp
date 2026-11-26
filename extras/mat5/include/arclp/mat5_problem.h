/* A loaded .mat problem as the solver's input; links nothing. */
#ifndef ARCLP_MAT5_PROBLEM_H
#define ARCLP_MAT5_PROBLEM_H

#include <arclp/arclp.h>
#include <arclp/mat5.h>

static inline arclp_problem mat_lp_as_problem(const mat_lp *lp)
{
    const arclp_problem problem = {lp->m, lp->n, lp->col_ptr, lp->row_idx,
                                   lp->values, lp->b, lp->c, NULL, NULL};
    return problem;
}

#endif
