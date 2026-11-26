/*
 * Newton systems shared by the step rules.  With theta = s/x on the pairs
 * of a column (plus z/w when it is bounded) and D = 1/theta, the normal
 * matrix is M = A D A'.  Each step is split so that the caller combines
 * the two matrix products as it likes:
 *
 *   first order     M dy = A g - b          arc_newton_first_rhs
 *                   (dx, ds) from A'dy      arc_newton_first_recover
 *   second order    M dy = -A w             arc_newton_second_rhs
 *   x ds + s dx = target, A dx = 0          arc_newton_second_recover
 *
 * Directions are subtracted from the iterate (first order) or added
 * (second order), as in arcLP.m.  Arrays of pairs have n_total entries.
 */
#ifndef ARC_NEWTON_H
#define ARC_NEWTON_H

#include "ipm.h"

/* d = D for the current iterate. */
void arc_newton_scaling(arc_ipm *ipm);

void arc_newton_first_rhs(const arc_ipm *ipm, double *g);

/* ds holds A'dy on entry. */
void arc_newton_first_recover(const arc_ipm *ipm, double *dx, double *ds);

void arc_newton_second_rhs(const arc_ipm *ipm, const double *target,
                           double *w);

/* ds holds A'dy on entry; w comes from arc_newton_second_rhs. */
void arc_newton_second_recover(const arc_ipm *ipm, const double *target,
                               const double *w, double *dx, double *ds);

#endif
