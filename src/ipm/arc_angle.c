/*
 * Step selection of arcLP.m.  For each centering parameter sigma the
 * largest angle keeping x >= phi and s >= psi follows in closed form per
 * component; sigma is bisected (Algorithm 6.1) to balance the components
 * whose second derivative decreases with sigma against those that grow.
 */
#include "arc_bounds.h"

#include "util/vec.h"

#include <math.h>

static const double HALF_PI = 1.57079632679489661923;
static const double RHO = 0.1;
static const double SIGMA_MIN = 1e-6;
static const double SIGMA_MAX = 0.3;
static const double SIGMA_WIDTH = 1e-4;   /* sqrt(1e-8) in arcLP.m */
static const double ANGLE_SHRINK = 0.999;
static const double ANGLE_CAP_FRACTION = 0.99;
static const double SMALLEST_ANGLE = 1e-300;

/* Halves alpha until mu(sigma, alpha) < mu (Step 2 of Algorithm 5.1). */
static double reduce_for_mu(const arc_search_rule *r, const arc_ipm *ipm,
                            double alpha, double sigma)
{
    const int n = ipm->n;
    const double a1 = arc_dot(n, r->dot.x, r->dot.s) / n;
    const double a2 = (arc_dot(n, r->dot.x, r->q.s)
                       + arc_dot(n, r->dot.s, r->q.x)) / n;
    const double a3 = (arc_dot(n, r->dot.x, r->p.s)
                       + arc_dot(n, r->dot.s, r->p.x)) / n;
    const double mu = ipm->mu;
    while (alpha > SMALLEST_ANGLE) {
        const double sa = sin(alpha);
        const double ca = 1.0 - cos(alpha);
        const double mu_alpha =
            mu * (1.0 - sa) + ca * ((mu - a3 * sa) * sigma - a1 * ca - a2 * sa);
        if (mu > mu_alpha) {
            break;
        }
        alpha *= 0.5;
    }
    return alpha;
}

arc_angle arc_search_angle(const arc_search_rule *r, const arc_ipm *ipm)
{
    const double phi = arc_fmin(RHO * arc_min(ipm->n, ipm->x), r->nu);
    const double psi = arc_fmin(r->nu, RHO * arc_min(ipm->n, ipm->s));
    double lower = SIGMA_MIN;
    double upper = SIGMA_MAX;
    double sigma = lower;
    const arc_bound quarter = {HALF_PI, {1.0, 0.0, 1.0}};
    arc_bounds b = {quarter, quarter, quarter};
    while (upper - lower > SIGMA_WIDTH) {
        sigma = lower + 0.5 * (upper - lower);
        b = arc_bounds_at(r, ipm, sigma, phi, psi);
        if (b.minus.angle > b.plus.angle) {
            lower = sigma;
        } else {
            upper = sigma;
        }
    }
    double alpha = arc_fmin(arc_fmin(b.minus.angle, b.plus.angle),
                            b.all.angle);
    alpha = reduce_for_mu(r, ipm, alpha, sigma);
    const arc_angle result = {arc_fmin(ANGLE_SHRINK * alpha,
                                         ANGLE_CAP_FRACTION * HALF_PI), sigma};
    return result;
}
