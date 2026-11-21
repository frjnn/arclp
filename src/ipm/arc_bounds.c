/*
 * Step-angle bounds for one centering parameter sigma (Algorithm 6.1,
 * Step 2).  Each running minimum keeps its trigonometry, so a component
 * that provably cannot lower it is skipped without evaluating arcsines;
 * the minimum only changes O(log n) times per pass.
 */
#include "arc_bounds.h"

#include <math.h>

static const double HALF_PI = 1.57079632679489661923;

static inline void lower(arc_bound *b, double angle)
{
    if (angle < b->angle)
    {
        b->angle = angle;
        b->trig = arc_trig_of(angle);
    }
}

/* dp is the sigma-derivative of the component's second-order term. */
static inline void consider(arc_bounds *b, double gap, double v1, double v2,
                            double dp)
{
    arc_bound *relevant = dp < 0.0   ? &b->minus
                          : dp > 0.0 ? &b->plus
                                     : &b->all;
    if (!arc_angle_may_be_below(gap, v1, v2, relevant->trig))
    {
        return;
    }
    const double angle = arc_component_angle(gap, v1, v2);
    lower(&b->all, angle);
    if (dp < 0.0)
    {
        lower(&b->minus, angle);
    }
    else if (dp > 0.0)
    {
        lower(&b->plus, angle);
    }
}

arc_bounds arc_bounds_at(const arc_search_rule *r, const arc_ipm *ipm,
                         double sigma, double phi, double psi)
{
    const arc_bound quarter = {HALF_PI, {1.0, 0.0, 1.0}};
    arc_bounds b = {quarter, quarter, quarter};
    for (int i = 0; i < ipm->n; i++)
    {
        consider(&b, ipm->x[i] - phi, r->dot.x[i],
                 sigma * r->p.x[i] + r->q.x[i], r->p.x[i]);
        consider(&b, ipm->s[i] - psi, r->dot.s[i],
                 sigma * r->p.s[i] + r->q.s[i], r->p.s[i]);
    }
    return b;
}
