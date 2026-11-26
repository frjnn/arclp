/* The closed-form step angle keeps the component feasible up to the angle
 * and, when it is below pi/2, hits the boundary exactly there; a collapsed
 * arc falls back to the first derivative. */
#include "check.h"

#include "ipm/arc_filter.h"
#include "ipm/arc_search.h"

#include <string.h>

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
    for (int k = 0; k <= 100; k++)
    {
        CHECK(along_arc(gap, v1, v2, a * k / 100.0) >= -1e-12 * scale);
    }
    if (a < HALF_PI)
    {
        CHECK_NEAR(along_arc(gap, v1, v2, a), 0.0, 1e-10 * scale);
    }
}

/* The filter may only skip components whose angle is not below t. */
static void check_filter(double gap, double v1, double v2, double t)
{
    if (arc_component_angle(gap, v1, v2) < t)
    {
        CHECK(arc_angle_may_be_below(gap, v1, v2, arc_trig_of(t)));
    }
}

/* Two pairs at x = s = 1 with first derivative 1/2 and second derivative
 * q; p is zero.  A side whose q is huge has no usable angle, and only when
 * both have none does the step drop the second derivative. */
static void collapsed_arc_uses_first_derivative(void)
{
    double x[] = {1.0, 1.0}, s[] = {1.0, 1.0};
    double half[] = {0.5, 0.5}, zero[] = {0.0, 0.0};
    double huge[] = {-1e30, -1e30};
    arc_ipm ipm;
    memset(&ipm, 0, sizeof(ipm));
    ipm.x = x;
    ipm.s = s;
    ipm.n_total = 2;
    ipm.n_pairs = 2;
    arc_search_rule r;
    memset(&r, 0, sizeof(r));
    r.nu = 1.0;
    r.dot.x = half;
    r.dot.s = half;
    r.p.x = r.p.s = r.q.x = r.q.s = zero;
    arc_angle a = arc_search_angle(&r, &ipm);
    CHECK(a.curvature == 1.0 && a.alpha_x > 1.0 && a.alpha_s > 1.0);
    r.q.x = huge;
    a = arc_search_angle(&r, &ipm);
    CHECK(a.curvature == 1.0 && a.alpha_x < ARC_STALL_ANGLE
          && a.alpha_s > 1.0);
    r.q.s = huge;
    a = arc_search_angle(&r, &ipm);
    CHECK(a.curvature == 0.0 && a.alpha_x > 1.0 && a.alpha_s > 1.0);
}

int main(void)
{
    collapsed_arc_uses_first_derivative();
    check_case(1.0, 0.0, -3.0);     /* case 1: no first-order motion   */
    check_case(1.0, 4.0, 0.0);      /* case 2: no second-order motion  */
    check_case(1.0, 2.0, 0.5);      /* case 3                          */
    check_case(1.0, 2.0, -0.5);     /* case 4                          */
    check_case(1.0, -0.5, -3.0);    /* case 5                          */
    check_case(1.0, -1.0, 2.0);     /* moving away: full quarter turn  */
    CHECK_NEAR(arc_component_angle(1.0, -1.0, 2.0), HALF_PI, 0.0);
    unsigned seed = 3;
    for (int t = 0; t < 2000; t++)
    {
        check_case(check_random(&seed) + 1e-3, 4.0 * check_random(&seed) - 2.0,
                   4.0 * check_random(&seed) - 2.0);
    }
    for (int t = 0; t < 200000; t++)
    {
        const double gap = check_random(&seed) + 1e-6;
        const double v1 = 6.0 * check_random(&seed) - 3.0;
        const double v2 = 6.0 * check_random(&seed) - 3.0;
        check_filter(gap, v1, v2, HALF_PI * check_random(&seed));
        check_filter(gap, t % 7 == 0 ? 0.0 : v1, t % 5 == 0 ? 0.0 : v2,
                     HALF_PI * check_random(&seed));
    }
    return check_report("arc component angle");
}
