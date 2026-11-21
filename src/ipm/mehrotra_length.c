/*
 * Step length of mehrotra.m: bisection on [0, 1] to width 1e-6 for the
 * largest a keeping v + a (dd - d) > 0, then a (1 -/+ shrink) perturbation.
 *
 * min_i (v_i + a g_i) is concave and positive at a = 0, so the feasible
 * steps form an interval [0, a*) whose end a* is a ratio test.  The
 * bisection's path is therefore decided by comparing each midpoint with a*,
 * without touching the vectors; only a midpoint within rounding distance of
 * a* is checked componentwise, exactly as the reference does.  Two passes
 * over the data replace the reference's twenty-two.
 */
#include "mehrotra_steps.h"

#include "util/vec.h"

static const double BISECTION_WIDTH = 1e-6;
static const double SMALLEST_STEP = 1e-8;
static const double BOUNDARY_FUZZ = 1e-12;

/* v_i + a (dd_i - d_i), evaluated as in the reference. */
static inline double moved(double v, double d, double dd, double a)
{
    return v - d * a + dd * a;
}

static double lowest_at(int len, const double *v, const double *d,
                        const double *dd, double a)
{
    double lowest = INFINITY;
    for (int i = 0; i < len; i++)
    {
        lowest = arc_fmin(lowest, moved(v[i], d[i], dd[i], a));
    }
    return lowest;
}

/* a* = min over decreasing components of v_i / (d_i - dd_i). */
static double boundary(int len, const double *v, const double *d,
                       const double *dd)
{
    double limit = INFINITY;
    for (int i = 0; i < len; i++)
    {
        const double fall = d[i] - dd[i];
        if (fall > 0.0)
        {
            limit = arc_fmin(limit, v[i] / fall);
        }
    }
    return limit;
}

static int feasible(int len, const double *v, const double *d,
                    const double *dd, double a, double limit)
{
    if (limit == INFINITY || fabs(a - limit) > BOUNDARY_FUZZ * limit)
    {
        return a < limit;
    }
    return lowest_at(len, v, d, dd, a) > 0.0;
}

typedef struct
{
    double step, shorter, longer;
} three_minima;

/* The three minima the final perturbation compares, in one pass. */
static three_minima lowest_at_three(int len, const double *v, const double *d,
                                    const double *dd, three_minima a)
{
    three_minima low = {INFINITY, INFINITY, INFINITY};
    for (int i = 0; i < len; i++)
    {
        low.step = arc_fmin(low.step, moved(v[i], d[i], dd[i], a.step));
        low.shorter =
            arc_fmin(low.shorter, moved(v[i], d[i], dd[i], a.shorter));
        low.longer = arc_fmin(low.longer, moved(v[i], d[i], dd[i], a.longer));
    }
    return low;
}

double mehrotra_step_length(int len, const double *v, const double *d,
                            const double *dd, double shrink)
{
    const double limit = boundary(len, v, d, dd);
    double lo = 0.0;
    double hi = 1.0;
    while (hi - lo > BISECTION_WIDTH)
    {
        const double mid = 0.5 * (lo + hi);
        if (feasible(len, v, d, dd, mid, limit))
        {
            lo = mid;
        }
        else
        {
            hi = mid;
        }
    }
    const three_minima steps = {lo,
                                arc_fmax((1.0 - shrink) * lo, SMALLEST_STEP),
                                arc_fmin((1.0 + shrink) * lo, 1.0)};
    const three_minima low = lowest_at_three(len, v, d, dd, steps);
    double step = steps.step;
    if (low.shorter > low.step)
    {
        step = steps.shorter;
    }
    if (low.longer > arc_fmin(low.step, low.shorter))
    {
        step = steps.longer;
    }
    return step;
}
