#include "mehrotra_steps.h"

#include "util/vec.h"

#include <math.h>

static const double MAX_CENTERING = 0.5;

/* Largest step in (0, 1] keeping v - a d >= 0. */
static double ratio_test(int n, const double *v, const double *d)
{
    double step = 1.0;
    for (int i = 0; i < n; i++) {
        if (d[i] > 0.0) {
            step = arc_fmin(step, v[i] / d[i]);
        }
    }
    return step;
}

void mehrotra_affine(mehrotra_rule *r, arc_ipm *ipm)
{
    double *w = r->work;
    for (int j = 0; j < ipm->n; j++) {
        w[j] = ipm->d[j] * ipm->rc[j];
    }
    arc_csc_mul(ipm->a, w, r->dy);
    arc_subtract(ipm->m, ipm->b, r->dy);
    arc_normal_solve(ipm->solver, r->dy, 1);
    arc_csc_mul_t(ipm->a, r->dy, r->ds);
    for (int j = 0; j < ipm->n; j++) {
        r->ds[j] = ipm->rc[j] - r->ds[j];
        r->dx[j] = ipm->x[j] - ipm->d[j] * r->ds[j];
    }
}

void mehrotra_corrector(mehrotra_rule *r, arc_ipm *ipm)
{
    const int n = ipm->n;
    const double ax = ratio_test(n, ipm->x, r->dx);
    const double as = ratio_test(n, ipm->s, r->ds);
    double gap = 0.0;
    for (int j = 0; j < n; j++) {
        gap += (ipm->x[j] - ax * r->dx[j]) * (ipm->s[j] - as * r->ds[j]);
    }
    const double ratio = gap / n / ipm->mu;
    const double target = arc_fmin(ratio * ratio * ratio, MAX_CENTERING)
                        * ipm->mu;
    double *w = r->work;
    for (int j = 0; j < n; j++) {
        w[j] = (2.0 * r->dx[j] * r->ds[j] - target) / ipm->s[j];
    }
    arc_csc_mul(ipm->a, w, r->ddy);
    arc_normal_solve(ipm->solver, r->ddy, 1);
    arc_csc_mul_t(ipm->a, r->ddy, r->dds);
    for (int j = 0; j < n; j++) {
        r->dds[j] = -r->dds[j];
        r->ddx[j] = (target - ipm->x[j] * r->dds[j]
                     - 2.0 * r->dx[j] * r->ds[j]) / ipm->s[j];
    }
}
