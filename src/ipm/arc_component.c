/*
 * Closed-form step angle for one component (Yang 2018, Section 6; cases
 * 1-5 of arcLP.m).  The boundary condition
 *     gap - v1 sin(a) + v2 (1 - cos(a)) = 0
 * is  v1 sin(a) + v2 cos(a) = gap + v2,  i.e.  r sin(a + beta) = gap + v2
 * with r = sqrt(v1^2 + v2^2) and sin(beta) = v2 / r.
 */
#include "arc_search.h"

#include <math.h>

static const double PI = 3.14159265358979323846;

/* First boundary crossing on (0, pi]; pi/2 when none is closer. */
static double crossing(double gap, double v1, double v2)
{
    const double target = gap + v2;
    if (v1 == 0.0) {
        return target < 0.0 ? acos(target / v2) : 0.5 * PI;
    }
    if (v2 == 0.0) {
        return v1 > gap ? asin(gap / v1) : 0.5 * PI;
    }
    const double r = sqrt(v1 * v1 + v2 * v2);   /* as in arcLP.m */
    if (v1 > 0.0 && target < r) {
        return asin(target / r) - asin(v2 / r);
    }
    if (v1 < 0.0 && v2 < 0.0 && target < 0.0) {
        return PI - asin(-target / r) - asin(-v2 / r);
    }
    return 0.5 * PI;
}

double arc_component_angle(double gap, double v1, double v2)
{
    return fmin(crossing(gap, v1, v2), 0.5 * PI);
}
