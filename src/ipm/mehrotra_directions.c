#include "mehrotra_steps.h"

#include "newton.h"
#include "util/vec.h"

#include <math.h>

static const double MAX_CENTERING = 0.5;

/* Largest step in (0, 1] keeping v - a d >= 0. */
static double ratio_test(int n, const double *v, const double *d)
{
    double step = 1.0;
    for (int i = 0; i < n; i++)
    {
        if (d[i] > 0.0)
        {
            step = arc_fmin(step, v[i] / d[i]);
        }
    }
    return step;
}

void mehrotra_affine(mehrotra_rule *r, arc_ipm *ipm)
{
    arc_newton_first_rhs(ipm, r->work);
    arc_csc_mul(ipm->a, r->work, r->dy);
    arc_subtract(ipm->m, ipm->b, r->dy);
    arc_normal_solve(ipm->solver, r->dy, 1);
    arc_csc_mul_t(ipm->a, r->dy, r->ds);
    arc_newton_first_recover(ipm, r->dx, r->ds);
}

/* (mu_aff / mu)^3 mu, with mu_aff after the largest affine steps. */
static double centering_target(const mehrotra_rule *r, const arc_ipm *ipm)
{
    const int f = ipm->n_free;
    const int pairs = ipm->n_pairs;
    const double ax = ratio_test(pairs, ipm->x + f, r->dx + f);
    const double as = ratio_test(pairs, ipm->s + f, r->ds + f);
    double gap = 0.0;
    for (int k = f; k < ipm->n_total; k++)
    {
        gap += (ipm->x[k] - ax * r->dx[k]) * (ipm->s[k] - as * r->ds[k]);
    }
    const double ratio = gap / pairs / ipm->mu;
    return arc_fmin(ratio * ratio * ratio, MAX_CENTERING) * ipm->mu;
}

void mehrotra_corrector(mehrotra_rule *r, arc_ipm *ipm)
{
    const double target = centering_target(r, ipm);
    for (int k = 0; k < ipm->n_total; k++)
    {
        r->target[k] = target - 2.0 * r->dx[k] * r->ds[k];
    }
    double *w = r->work;
    arc_newton_second_rhs(ipm, r->target, w);
    arc_csc_mul(ipm->a, w, r->ddy);
    arc_negate(ipm->m, r->ddy);
    arc_normal_solve(ipm->solver, r->ddy, 1);
    arc_csc_mul_t(ipm->a, r->ddy, r->dds);
    arc_newton_second_recover(ipm, r->target, w, r->ddx, r->dds);
}
