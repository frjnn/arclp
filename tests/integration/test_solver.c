/* End-to-end solves through the public API. */
#include "check.h"

#include <arclp/arclp.h>
#include <arclp/mat5.h>

#include <string.h>

static const arclp_method METHODS[] = {ARCLP_METHOD_ARC_SEARCH,
                                       ARCLP_METHOD_MEHROTRA};

/* min -x1 - 2 x2  s.t.  x1 + x2 + s1 = 4,  x1 + 3 x2 + s2 = 6  ->  -5 */
static void small_optimum(void)
{
    const int p[] = {0, 2, 4, 5, 6};
    const int i[] = {0, 1, 0, 1, 0, 1};
    const double v[] = {1, 1, 1, 3, 1, 1};
    const double b[] = {4, 6};
    const double c[] = {-1, -2, 0, 0};
    const arclp_problem lp = {2, 4, p, i, v, b, c};
    for (int presolve = 0; presolve <= 1; presolve++) {
        for (int k = 0; k < 2; k++) {
            arclp_options o;
            arclp_default_options(&o);
            o.method = METHODS[k];
            o.presolve = presolve;
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
    const arclp_problem lp = {1, 2, p, i, v, b, c};
    for (int k = 0; k < 2; k++) {
        arclp_options o;
        arclp_default_options(&o);
        o.method = METHODS[k];
        o.presolve = 0;
        double x[2];
        CHECK(arclp_solve(&lp, &o, x, NULL) == ARCLP_OPTIMAL);
        CHECK_NEAR(x[0], 0.0, 1e-7);
        CHECK_NEAR(x[1], 5.0, 1e-7);
    }
}

static arclp_status solve_raw(const arclp_problem *lp, arclp_method method)
{
    arclp_options o;
    arclp_default_options(&o);
    o.method = method;
    o.presolve = 0;
    double x[4];
    return arclp_solve(lp, &o, x, NULL);
}

static void infeasible_and_unbounded(void)
{
    /* x1 + x2 = 1 and x1 + x2 = 3 */
    const int p1[] = {0, 2, 4};
    const int i1[] = {0, 1, 0, 1};
    const double v1[] = {1, 1, 1, 1}, b1[] = {1, 3}, c1[] = {1, 1};
    const arclp_problem infeasible = {2, 2, p1, i1, v1, b1, c1};
    /* min -x1  s.t.  x1 - x2 = 1 */
    const int p2[] = {0, 1, 2};
    const int i2[] = {0, 0};
    const double v2[] = {1, -1}, b2[] = {1}, c2[] = {-1, 0};
    const arclp_problem unbounded = {1, 2, p2, i2, v2, b2, c2};
    for (int k = 0; k < 2; k++) {
        CHECK(solve_raw(&infeasible, METHODS[k]) != ARCLP_OPTIMAL);
        CHECK(solve_raw(&unbounded, METHODS[k]) != ARCLP_OPTIMAL);
    }
}

static void non_finite_input(void)
{
    const int p[] = {0, 1};
    const int i[] = {0};
    const double v[] = {1}, b[] = {1}, c[] = {NAN};
    const arclp_problem lp = {1, 1, p, i, v, b, c};
    double x[1];
    CHECK(arclp_solve(&lp, NULL, x, NULL) == ARCLP_INVALID_INPUT);
    const double vi[] = {INFINITY}, ci[] = {1};
    const arclp_problem inf = {1, 1, p, i, vi, b, ci};
    CHECK(arclp_solve(&inf, NULL, x, NULL) == ARCLP_INVALID_INPUT);
}

static void invalid_input(void)
{
    const int p[] = {0, 1};
    const int i[] = {5};               /* row out of range */
    const double v[] = {1}, b[] = {1}, c[] = {1};
    const arclp_problem lp = {1, 1, p, i, v, b, c};
    double x[1];
    CHECK(arclp_solve(&lp, NULL, x, NULL) == ARCLP_INVALID_INPUT);
    CHECK(arclp_solve(NULL, NULL, x, NULL) == ARCLP_INVALID_INPUT);
}

typedef struct {
    const char *name;
    double objective;
} netlib_case;

static void netlib(const char *dir)
{
    static const netlib_case cases[] = {
        {"afiro", -4.647531428571e+02}, {"adlittle", 2.254949631624e+05},
        {"blend", -3.081214984583e+01}, {"share2b", -4.157322407414e+02},
        {"degen2", -1.435178000000e+03}, {"scagr7", -2.331389824331e+06}};
    for (size_t t = 0; t < sizeof(cases) / sizeof(*cases); t++) {
        char path[4096], err[256];
        snprintf(path, sizeof(path), "%s/lp_%s.mat", dir, cases[t].name);
        mat_lp m;
        if (mat_lp_load(path, &m, err, sizeof(err)) != 0) {
            printf("  (skipping %s: %s)\n", cases[t].name, err);
            continue;
        }
        const arclp_problem lp = {m.m, m.n, m.col_ptr, m.row_idx, m.values,
                                  m.b, m.c};
        double *x = malloc(sizeof(double) * (size_t)m.n);
        for (int k = 0; k < 2; k++) {
            arclp_options o;
            arclp_default_options(&o);
            o.method = METHODS[k];
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
    non_finite_input();
    const char *dir = getenv("NETLIB_DIR");
    netlib(dir != NULL ? dir : "data/arcLP");
    return check_report("solver end to end");
}
