/* The closed-form step angle keeps the component feasible up to the angle
 * and, when it is below pi/2, hits the boundary exactly there. */
#include "check.h"

#include "ipm/arc_filter.h"
#include "ipm/arc_search.h"

static const double HALF_PI = 1.57079632679489661923;

static double along_arc(double gap, double v1, double v2, double a)
{
    return gap - v1 * sin(a) + v2 * (1.0 - cos(a));
}

static void check_case(double gap, double v1, double v2)
{
    const double a = arc_component_angle(gap, v1, v2);
    CHECK(a > 0.0 && a <= HALF_PI + 1e-15);
    const double scale = fabs(gap) + fabs(v1) + fabs(v2);
    for (int k = 0; k <= 100; k++) {
        CHECK(along_arc(gap, v1, v2, a * k / 100.0) >= -1e-12 * scale);
    }
    if (a < HALF_PI) {
        CHECK_NEAR(along_arc(gap, v1, v2, a), 0.0, 1e-10 * scale);
    }
}

/* The filter may only skip components whose angle is not below t. */
static void check_filter(double gap, double v1, double v2, double t)
{
    if (arc_component_angle(gap, v1, v2) < t) {
        CHECK(arc_angle_may_be_below(gap, v1, v2, arc_trig_of(t)));
    }
}

int main(void)
{
    check_case(1.0, 0.0, -3.0);     /* case 1: no first-order motion   */
    check_case(1.0, 4.0, 0.0);      /* case 2: no second-order motion  */
    check_case(1.0, 2.0, 0.5);      /* case 3                          */
    check_case(1.0, 2.0, -0.5);     /* case 4                          */
    check_case(1.0, -0.5, -3.0);    /* case 5                          */
    check_case(1.0, -1.0, 2.0);     /* moving away: full quarter turn  */
    CHECK_NEAR(arc_component_angle(1.0, -1.0, 2.0), HALF_PI, 0.0);
    unsigned seed = 3;
    for (int t = 0; t < 2000; t++) {
        check_case(check_random(&seed) + 1e-3, 4.0 * check_random(&seed) - 2.0,
                   4.0 * check_random(&seed) - 2.0);
    }
    for (int t = 0; t < 200000; t++) {
        const double gap = check_random(&seed) + 1e-6;
        const double v1 = 6.0 * check_random(&seed) - 3.0;
        const double v2 = 6.0 * check_random(&seed) - 3.0;
        check_filter(gap, v1, v2, HALF_PI * check_random(&seed));
        check_filter(gap, t % 7 == 0 ? 0.0 : v1, t % 5 == 0 ? 0.0 : v2,
                     HALF_PI * check_random(&seed));
    }
    return check_report("arc component angle");
}
