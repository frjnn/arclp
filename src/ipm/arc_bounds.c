/*
 * Step-angle bounds for one centering parameter sigma (Algorithm 6.1,
 * Step 2).  Each running minimum keeps its trigonometry, so a component
 * that provably cannot lower it is skipped without evaluating arcsines;
 * the minimum only changes O(log n) times per pass.
 */
#include "arc_bounds.h"

#include "util/parallel.h"

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

static ARC_ALWAYS_INLINE void scan(arc_bounds *b, const arc_search_rule *r,
                                   const arc_ipm *ipm, double sigma,
                                   double phi, double psi, int lo, int hi)
{
    for (int i = lo; i < hi; i++)
    {
        consider(b, ipm->x[i] - phi, r->dot.x[i],
                 sigma * r->p.x[i] + r->q.x[i], r->p.x[i]);
        consider(b, ipm->s[i] - psi, r->dot.s[i],
                 sigma * r->p.s[i] + r->q.s[i], r->p.s[i]);
    }
}

static const arc_bound QUARTER = {HALF_PI, {1.0, 0.0, 1.0}};

#if defined(ARCLP_USE_OPENMP)
static void merge(arc_bounds *into, const arc_bounds *part)
{
    lower(&into->all, part->all.angle);
    lower(&into->minus, part->minus.angle);
    lower(&into->plus, part->plus.angle);
}

/* One slice of components per thread, each with its own minima, merged at
 * the end: the minima do not depend on the split, so the result is the
 * same for any number of threads. */
static arc_bounds parallel_bounds(const arc_search_rule *r,
                                  const arc_ipm *ipm, double sigma,
                                  double phi, double psi)
{
    arc_bounds b = {QUARTER, QUARTER, QUARTER};
#pragma omp parallel
    {
        const long threads = omp_get_num_threads();
        const long t = omp_get_thread_num();
        arc_bounds part = {QUARTER, QUARTER, QUARTER};
        scan(&part, r, ipm, sigma, phi, psi, (int)(ipm->n * t / threads),
             (int)(ipm->n * (t + 1) / threads));
#pragma omp critical(arc_bounds_merge)
        merge(&b, &part);
    }
    return b;
}
#endif

arc_bounds arc_bounds_at(const arc_search_rule *r, const arc_ipm *ipm,
                         double sigma, double phi, double psi)
{
#if defined(ARCLP_USE_OPENMP)
    if (16L * ipm->n >= ARC_PARALLEL_WORK)   /* ~16 units per component */
    {
        return parallel_bounds(r, ipm, sigma, phi, psi);
    }
#endif
    arc_bounds b = {QUARTER, QUARTER, QUARTER};
    scan(&b, r, ipm, sigma, phi, psi, 0, ipm->n);
    return b;
}
