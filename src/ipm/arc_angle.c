/*
 * Step selection of arcLP.m.  For each centering parameter sigma the
 * largest angle keeping x >= phi and s >= psi follows in closed form per
 * component; sigma is bisected (Algorithm 6.1) to balance the components
 * whose second derivative decreases with sigma against those that grow.
 * Each side then takes its own bound at that sigma as its angle, without
 * the halving of Step 2 of Algorithm 5.1 that keeps mu falling: on an
 * infeasible problem mu must rise while the dual moves off towards the
 * certificate.  A second derivative too large for any angle falls back to
 * a step along the first derivative.
 */
#include "arc_bounds.h"

#include "util/vec.h"

#include <math.h>

static const double HALF_PI = 1.57079632679489661923;
static const double RHO = 0.1;
static const double SIGMA_MIN = 1e-6;
static const double SIGMA_MAX = 0.3;
static const double SIGMA_WIDTH = 1e-4; /* sqrt(1e-8) in arcLP.m */
static const double ANGLE_SHRINK = 0.999;
static const double ANGLE_CAP_FRACTION = 0.99;

/* An angle a little inside the boundary, and short of a quarter turn. */
static double usable(double angle)
{
    return arc_fmin(ANGLE_SHRINK * angle, ANGLE_CAP_FRACTION * HALF_PI);
}

/* The bound a side keeps: a share of its smallest entry, at most nu. */
static double floor_of(const arc_search_rule *r, const arc_ipm *ipm,
                       const double *v)
{
    return arc_fmin(RHO * arc_min(ipm->n_pairs, v + ipm->n_free), r->nu);
}

/* The largest angles along the first derivative alone. */
static arc_angle first_order(const arc_search_rule *r, const arc_ipm *ipm,
                             double phi, double psi)
{
    double ax = HALF_PI;
    double as = HALF_PI;
    for (int i = ipm->n_free; i < ipm->n_total; i++)
    {
        ax = arc_fmin(ax, arc_component_angle(ipm->x[i] - phi, r->dot.x[i],
                                              0.0));
        as = arc_fmin(as, arc_component_angle(ipm->s[i] - psi, r->dot.s[i],
                                              0.0));
    }
    const arc_angle result = {usable(ax), usable(as), 0.0, 0.0};
    return result;
}

arc_angle arc_search_angle(const arc_search_rule *r, const arc_ipm *ipm)
{
    const double phi = floor_of(r, ipm, ipm->x);
    const double psi = floor_of(r, ipm, ipm->s);
    double lower = SIGMA_MIN;
    double upper = SIGMA_MAX;
    double sigma = lower;
    const arc_bound quarter = {HALF_PI, {1.0, 0.0, 1.0}};
    arc_bounds b = {quarter, quarter, quarter, quarter, quarter};
    while (upper - lower > SIGMA_WIDTH)
    {
        sigma = lower + 0.5 * (upper - lower);
        b = arc_bounds_at(r, ipm, sigma, phi, psi);
        if (b.minus.angle > b.plus.angle)
        {
            lower = sigma;
        }
        else
        {
            upper = sigma;
        }
    }
    const arc_angle result = {usable(b.x.angle), usable(b.s.angle), sigma,
                              1.0};
    if (arc_fmax(result.alpha_x, result.alpha_s) < ARC_STALL_ANGLE)
    {
        return first_order(r, ipm, phi, psi);
    }
    return result;
}
