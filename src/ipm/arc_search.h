/*
 * Arc-search rule of Yang (2018), Algorithm 5.1.  The iterate moves along
 * the ellipse
 *     v(alpha) = v - v' sin(alpha) + v''(sigma) (1 - cos(alpha)),
 * where v' solves the Newton system for the residuals and
 * v''(sigma) = sigma p + q is the second-order term, linear in the
 * centering parameter sigma.
 */
#ifndef ARC_SEARCH_H
#define ARC_SEARCH_H

#include "step_rule.h"

typedef struct
{
    double *x, *y, *s;
} arc_direction;

typedef struct
{
    arc_step_rule base; /* must stay first */
    arc_direction dot;  /* first derivative of the arc            */
    arc_direction p, q; /* second derivative is sigma p + q       */
    double *work;       /* n scratch                              */
    double *work_p;     /* n scratch                              */
    double *target;     /* complementarity targets, all pairs     */
    double nu;          /* nu_k = prod (1 - sin alpha_j)          */
    double alpha;       /* last accepted angle                    */
} arc_search_rule;

typedef struct
{
    double alpha;
    double sigma;
} arc_angle;

/* Computes dot, p and q from the factorized normal equations. */
void arc_search_directions(arc_search_rule *rule, arc_ipm *ipm);

/* Algorithm 6.1 plus the duality-measure safeguard and rescaling. */
arc_angle arc_search_angle(const arc_search_rule *rule, const arc_ipm *ipm);

/* Largest angle in (0, pi/2] keeping
 *     gap - v1 sin(a) + v2 (1 - cos(a)) >= 0,
 * for gap > 0 (the five cases of arcLP.m). */
double arc_component_angle(double gap, double v1, double v2);

#endif
