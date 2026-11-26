/* Stages of the Mehrotra rule; internal to mehrotra_*.c. */
#ifndef ARC_MEHROTRA_STEPS_H
#define ARC_MEHROTRA_STEPS_H

#include "step_rule.h"

typedef struct
{
    arc_step_rule base;     /* must stay first */
    double *dx, *dy, *ds;   /* affine-scaling direction          */
    double *ddx, *ddy, *dds;/* centering-corrector direction     */
    double *work, *target;
    double step_x, step_s;
    int iteration;
} mehrotra_rule;

/* Affine direction; same Newton system as the arc's first derivative. */
void mehrotra_affine(mehrotra_rule *r, arc_ipm *ipm);

/* Centering parameter (mu_aff / mu)^3 and the corrector direction. */
void mehrotra_corrector(mehrotra_rule *r, arc_ipm *ipm);

/* Step keeping v + a (dd - d) > 0: bisection to 1e-6, then the
 * (1 -/+ shrink) perturbation of mehrotra.m. */
double mehrotra_step_length(int len, const double *v, const double *d,
                            const double *dd, double shrink);

#endif
