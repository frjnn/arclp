/* End-to-end solves through the public API. */
#include "check.h"

#include <arclp/arclp.h>
#include <arclp/mat5_problem.h>

#include "small_lp.h"

#include <string.h>

static const arclp_method METHODS[] = {ARCLP_METHOD_ARC_SEARCH,
                                       ARCLP_METHOD_MEHROTRA};

static arclp_options options(arclp_method method, int presolve)
{
    arclp_options o;
    arclp_default_options(&o);
    o.method = method;
    o.presolve = presolve;
    return o;
}

static void small_optimum(void)
{
    const arclp_problem lp = {SMALL_M, SMALL_N, SMALL_P, SMALL_I, SMALL_V,
                              SMALL_B, SMALL_C, NULL, NULL};
    for (int presolve = 0; presolve <= 1; presolve++)
    {
        for (int k = 0; k < 2; k++)
        {
            const arclp_options o = options(METHODS[k], presolve);
            double x[4];
            arclp_info info;
            CHECK(arclp_solve(&lp, &o, x, &info) == ARCLP_OPTIMAL);
            CHECK_NEAR(info.objective, -5.0, 1e-7);
            CHECK_NEAR(x[0], 3.0, 1e-6);
            CHECK_NEAR(x[1], 1.0, 1e-6);
        }
    }
}

/* The two-variable example of infeasibleexample.m (the paper's Figure 2):
 * min x1  s.t.  x1 + x2 = 5,  optimum x = (0, 5). */
static void figure_two_example(void)
{
    const int p[] = {0, 1, 2};
    const int i[] = {0, 0};
    const double v[] = {1, 1}, b[] = {5}, c[] = {1, 0};
    const arclp_problem lp = {1, 2, p, i, v, b, c, NULL, NULL};
    for (int k = 0; k < 2; k++)
    {
        const arclp_options o = options(METHODS[k], 0);
        double x[2];
        CHECK(arclp_solve(&lp, &o, x, NULL) == ARCLP_OPTIMAL);
        CHECK_NEAR(x[0], 0.0, 1e-7);
        CHECK_NEAR(x[1], 5.0, 1e-7);
    }
}

/* The status of lp for every method, with and without presolve. */
static void expect_status(const arclp_problem *lp, arclp_status expected)
{
    for (int presolve = 0; presolve <= 1; presolve++)
    {
        for (int k = 0; k < 2; k++)
        {
            const arclp_options o = options(METHODS[k], presolve);
            double x[4];
            CHECK(arclp_solve(lp, &o, x, NULL) == expected);
        }
    }
}

static void infeasible_and_unbounded(void)
{
    /* x1 + x2 = 1 and x1 + x2 = 3 */
    const int p1[] = {0, 2, 4};
    const int i1[] = {0, 1, 0, 1};
    const double v1[] = {1, 1, 1, 1}, b1[] = {1, 3}, c1[] = {1, 1};
    const arclp_problem infeasible = {2, 2, p1, i1, v1, b1, c1, NULL, NULL};
    /* min -x1  s.t.  x1 - x2 = 1 */
    const int p2[] = {0, 1, 2};
    const int i2[] = {0, 0};
    const double v2[] = {1, -1}, b2[] = {1}, c2[] = {-1, 0};
    const arclp_problem unbounded = {1, 2, p2, i2, v2, b2, c2, NULL, NULL};
    /* min -x1  s.t.  x2 + x3 = -1: x1 lowers the cost without limit, but
     * no x >= 0 is feasible. */
    const int p3[] = {0, 0, 1, 2};
    const int i3[] = {0, 0};
    const double v3[] = {1, 1}, b3[] = {-1}, c3[] = {-1, 0, 0};
    const arclp_problem ray_but_infeasible = {1, 3, p3, i3, v3, b3, c3, NULL,
                                              NULL};
    expect_status(&infeasible, ARCLP_INFEASIBLE);
    expect_status(&unbounded, ARCLP_UNBOUNDED);
    expect_status(&ray_but_infeasible, ARCLP_INFEASIBLE);
}

/* A 1 x 1 problem with the given data is rejected. */
static void rejects(int row, double value, double cost)
{
    const int p[] = {0, 1};
    const double b[] = {1};
    const arclp_problem lp = {1, 1, p, &row, &value, b, &cost, NULL, NULL};
    double x[1];
    CHECK(arclp_solve(&lp, NULL, x, NULL) == ARCLP_INVALID_INPUT);
}

static void invalid_input(void)
{
    rejects(5, 1.0, 1.0);           /* row out of range */
    rejects(0, INFINITY, 1.0);
    rejects(0, 1.0, NAN);
    double x[1];
    CHECK(arclp_solve(NULL, NULL, x, NULL) == ARCLP_INVALID_INPUT);
}

typedef struct
{
    const char *name;
    double objective;
} netlib_case;

static void netlib(const char *dir)
{
    static const netlib_case cases[] = {
        {"afiro", -4.647531428571e+02}, {"adlittle", 2.254949631624e+05},
        {"blend", -3.081214984583e+01}, {"share2b", -4.157322407414e+02},
        {"degen2", -1.435178000000e+03}, {"scagr7", -2.331389824331e+06}};
    for (size_t t = 0; t < sizeof(cases) / sizeof(*cases); t++)
    {
        mat_lp m;
        if (mat_lp_load_named(dir, cases[t].name, &m) != 0)
        {
            continue;
        }
        const arclp_problem lp = mat_lp_as_problem(&m);
        double *x = malloc(sizeof(double) * (size_t)m.n);
        for (int k = 0; k < 2; k++)
        {
            const arclp_options o = options(METHODS[k], 1);
            arclp_info info;
            CHECK(arclp_solve(&lp, &o, x, &info) == ARCLP_OPTIMAL);
            const double ref = cases[t].objective;
            CHECK_NEAR(info.objective, ref, 1e-7 * fabs(ref));
        }
        free(x);
        mat_lp_free(&m);
    }
}

int main(void)
{
    small_optimum();
    figure_two_example();
    infeasible_and_unbounded();
    invalid_input();
    const char *dir = getenv("NETLIB_DIR");
    netlib(dir != NULL ? dir : "data/arcLP");
    return check_report("solver end to end");
}
