#include "newton.h"

#include "util/vec.h"

static const double FREE_FLOOR = 1.0;   /* |x| a free column is given at least */
static const double SMALLEST_MU = 1e-12;

/* The scaling a pair x+ = x- = t of a split free variable would have,
 * t^2 / (mu / 2 + mu / 2): a proximal term that vanishes with mu. */
static double free_scaling(double x, double mu)
{
    const double t = arc_fmax(FREE_FLOOR, fabs(x));
    return 2.0 * t * t / arc_fmax(mu, SMALLEST_MU);
}

static double bounded_scaling(const arc_ipm *ipm, int j)
{
    const int k = arc_ipm_slack(ipm, j);
    return 1.0 / (ipm->s[j] / ipm->x[j] + ipm->s[k] / ipm->x[k]);
}

void arc_newton_scaling(arc_ipm *ipm)
{
    double *d = ipm->d;
    ARC_BY_KIND(ipm, j,
                d[j] = free_scaling(ipm->x[j], ipm->mu),
                d[j] = ipm->x[j] / ipm->s[j],
                d[j] = bounded_scaling(ipm, j));
}

/* Slope of the slack dual minus the column's own: s - z + z rw / w. */
static inline double bound_pull(const arc_ipm *ipm, int j)
{
    const int k = arc_ipm_slack(ipm, j);
    return ipm->s[j] - ipm->s[k]
         + ipm->s[k] * ipm->rw[k - ipm->n] / ipm->x[k];
}

void arc_newton_first_rhs(const arc_ipm *ipm, double *g)
{
    const double *x = ipm->x;
    const double *d = ipm->d;
    const double *rc = ipm->rc;
    ARC_BY_KIND(ipm, j,
                g[j] = x[j] + d[j] * rc[j],
                g[j] = d[j] * rc[j],
                g[j] = x[j] + d[j] * (rc[j] - bound_pull(ipm, j)));
}

/* Directions of a bounded column and its slack once dx_j is known: x + w
 * changes by slack_sum and each dual is base - (s / x) dx. */
typedef struct
{
    double slack_sum, base_s, base_z;
} pair_targets;

static inline void close_pair(const arc_ipm *ipm, int j, double dxj,
                              pair_targets t, double *dx, double *ds)
{
    const int k = arc_ipm_slack(ipm, j);
    dx[j] = dxj;
    dx[k] = t.slack_sum - dxj;
    ds[j] = t.base_s - ipm->s[j] / ipm->x[j] * dx[j];
    ds[k] = t.base_z - ipm->s[k] / ipm->x[k] * dx[k];
}

static inline void first_free(const arc_ipm *ipm, int j, double *dx,
                              double *ds)
{
    dx[j] = ipm->d[j] * (ds[j] - ipm->rc[j]);
    ds[j] = 0.0;
}

static inline void first_plain(const arc_ipm *ipm, int j, double *dx,
                               double *ds)
{
    ds[j] = ipm->rc[j] - ds[j];
    dx[j] = ipm->x[j] - ipm->d[j] * ds[j];
}

static inline void first_bounded(const arc_ipm *ipm, int j, double *dx,
                                 double *ds)
{
    const int k = arc_ipm_slack(ipm, j);
    const pair_targets t = {ipm->rw[k - ipm->n], ipm->s[j], ipm->s[k]};
    close_pair(ipm, j, ipm->d[j] * (ds[j] - ipm->rc[j] + bound_pull(ipm, j)),
               t, dx, ds);
}

void arc_newton_first_recover(const arc_ipm *ipm, double *dx, double *ds)
{
    ARC_BY_KIND(ipm, j, first_free(ipm, j, dx, ds),
                first_plain(ipm, j, dx, ds), first_bounded(ipm, j, dx, ds));
}

static inline double bounded_weight(const arc_ipm *ipm, const double *target,
                                    int j)
{
    const int k = arc_ipm_slack(ipm, j);
    return ipm->d[j] * (target[j] / ipm->x[j] - target[k] / ipm->x[k]);
}

void arc_newton_second_rhs(const arc_ipm *ipm, const double *target,
                           double *w)
{
    ARC_BY_KIND(ipm, j, w[j] = 0.0, w[j] = target[j] / ipm->s[j],
                w[j] = bounded_weight(ipm, target, j));
}

static inline void second_free(const arc_ipm *ipm, int j, double *dx,
                               double *ds)
{
    dx[j] = ipm->d[j] * ds[j];
    ds[j] = 0.0;
}

static inline void second_plain(const arc_ipm *ipm, const double *w, int j,
                                double *dx, double *ds)
{
    dx[j] = w[j] + ipm->d[j] * ds[j];
    ds[j] = -ds[j];
}

static inline void second_bounded(const arc_ipm *ipm, const double *target,
                                  const double *w, int j, double *dx,
                                  double *ds)
{
    const int k = arc_ipm_slack(ipm, j);
    const pair_targets t = {0.0, target[j] / ipm->x[j],
                            target[k] / ipm->x[k]};
    close_pair(ipm, j, ipm->d[j] * ds[j] + w[j], t, dx, ds);
}

void arc_newton_second_recover(const arc_ipm *ipm, const double *target,
                               const double *w, double *dx, double *ds)
{
    ARC_BY_KIND(ipm, j, second_free(ipm, j, dx, ds),
                second_plain(ipm, w, j, dx, ds),
                second_bounded(ipm, target, w, j, dx, ds));
}
