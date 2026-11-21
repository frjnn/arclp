/*
 * Cheap conservative test run for every component in every bisection pass,
 * hence inline.  g(a) = gap - v1 sin(a) + v2 (1 - cos(a)) is a sinusoid,
 * so its minimum over [0, t] (t <= pi/2) lies at t or at the interior
 * critical point tan(a*) = v1 / v2, a minimum inside (0, t) only for
 * v1, v2 > 0 and v1 cos(t) < v2 sin(t), where g(a*) = gap + v2 - r.
 * The margin keeps the test conservative under rounding.
 */
#ifndef ARC_FILTER_H
#define ARC_FILTER_H

#include <math.h>

/* Precomputed trigonometry of a threshold angle t. */
typedef struct
{
    double sin, cos, versine; /* versine = 1 - cos(t) */
} arc_trig;

static inline arc_trig arc_trig_of(double t)
{
    const arc_trig trig = {sin(t), cos(t), 1.0 - cos(t)};
    return trig;
}

/* Returns 0 only when arc_component_angle(gap, v1, v2) >= t. */
static inline int arc_angle_may_be_below(double gap, double v1, double v2,
                                         arc_trig t)
{
    const double margin = 1e-9 * (gap + fabs(v1) + fabs(v2));
    if (gap - v1 * t.sin + v2 * t.versine < margin)
    {
        return 1;
    }
    if (v1 > 0.0 && v2 > 0.0 && v1 * t.cos < v2 * t.sin)
    {
        return gap + v2 < sqrt(v1 * v1 + v2 * v2) + margin;
    }
    return 0;
}

#endif
