/* Each presolve rule on a small problem, plus postsolve consistency. */
#include "check.h"

#include "presolve/presolve.h"

typedef struct
{
    int m, n;
    int p[8];
    int i[16];
    double x[16];
    double b[4];
    double c[6];
} small_lp;

static arc_presolve_status presolve(const small_lp *lp, const double *lower,
                                    const double *upper, int rules,
                                    arc_presolve *ps)
{
    const arc_csc a = {lp->m, lp->n, lp->p, lp->i, lp->x};
    return arc_presolve_run(ps, &a, lp->b, lp->c, lower, upper, rules);
}

/* Rule 5 fixes x0 = 2 from row 0; rule 9 then substitutes x1 = 3 + x2 out
 * of row 1; row 2, an inequality on one variable, becomes 0 <= x2 <= 1, and
 * the empty column x2 goes to the bound its remaining cost prefers. */
static void singleton_substitution_and_bound(void)
{
    /* row0: 2 x0 = 4;  row1: x0 + x1 - x2 = 5;  row2: x2 + x3 = 1 */
    const small_lp lp = {3, 4, {0, 2, 3, 5, 6}, {0, 1, 1, 1, 2, 2},
                         {2, 1, 1, -1, 1, 1}, {4, 5, 1}, {1, 1, 1, 3}};
    arc_presolve ps;
    CHECK(presolve(&lp, NULL, NULL, 1, &ps)
          == ARC_PRESOLVE_REDUCED);
    CHECK(ps.reduced.a.m == 0 && ps.reduced.a.n == 0);
    const double xr[1] = {0.0};
    double x[4];
    arc_postsolve(&ps, xr, x);
    CHECK_NEAR(x[0], 2.0, 1e-15);
    CHECK_NEAR(x[0] + x[1] - x[2], 5.0, 1e-14);
    CHECK_NEAR(x[2] + x[3], 1.0, 1e-14);
    CHECK_NEAR(x[2], 1.0, 1e-14);
    arc_presolve_free(&ps);
}

/* Bounds -3 <= x0 <= -1, x1 <= 4, x2 free are mapped to [0, 2], [0, inf)
 * and free, and the original values come back in postsolve. */
static void bounds_are_normalised(void)
{
    const small_lp lp = {1, 3, {0, 1, 2, 3}, {0, 0, 0}, {1, 1, 1}, {5},
                         {1, 2, 3}};
    const double lower[] = {-3, -INFINITY, -INFINITY};
    const double upper[] = {-1, 4, INFINITY};
    arc_presolve ps;
    CHECK(presolve(&lp, lower, upper, 0, &ps)
          == ARC_PRESOLVE_REDUCED);
    CHECK(ps.reduced.n_free == 1 && ps.reduced.n_bounded == 1);
    CHECK_NEAR(ps.reduced.upper[0], 2.0, 0.0);
    double xr[3], x[3];
    xr[0] = 1.0;   /* x2, free          */
    xr[1] = 2.0;   /* x1 = 4 - xr[1]    */
    xr[2] = 0.5;   /* x0 = -3 + xr[2]   */
    arc_postsolve(&ps, xr, x);
    CHECK_NEAR(x[0], -2.5, 1e-15);
    CHECK_NEAR(x[1], 2.0, 1e-15);
    CHECK_NEAR(x[2], 1.0, 1e-15);
    arc_presolve_free(&ps);
}

static void forcing_row_fixes_zero(void)
{
    /* row0: x0 + x1 = 0 forces x0 = x1 = 0;  row1: x1 + x2 = 3 */
    const small_lp lp = {2, 3, {0, 1, 3, 4}, {0, 0, 1, 1}, {1, 1, 1, 1},
                         {0, 3}, {1, 1, 1}};
    arc_presolve ps;
    CHECK(presolve(&lp, NULL, NULL, 1, &ps) == ARC_PRESOLVE_REDUCED);
    double x[3];
    const double none[1] = {0.0};
    arc_postsolve(&ps, none, x);
    CHECK(x[0] == 0.0 && x[1] == 0.0);
    CHECK_NEAR(x[2], 3.0, 1e-14);
    arc_presolve_free(&ps);
}

static void detects_infeasible_and_unbounded(void)
{
    /* x0 + x1 = -1 with x >= 0 is infeasible (rule 7). */
    const small_lp infeasible = {1, 2, {0, 1, 2}, {0, 0}, {1, 1}, {-1},
                                 {1, 1}};
    arc_presolve ps;
    CHECK(presolve(&infeasible, NULL, NULL, 1, &ps)
          == ARC_PRESOLVE_INFEASIBLE);
    arc_presolve_free(&ps);
    /* Column 1 is empty with negative cost: unbounded (rule 3). */
    const small_lp unbounded = {1, 2, {0, 1, 1}, {0}, {1}, {1}, {1, -1}};
    CHECK(presolve(&unbounded, NULL, NULL, 1, &ps)
          == ARC_PRESOLVE_UNBOUNDED);
    arc_presolve_free(&ps);
}

/* Row 0 is a tenth of row 1, up to rounding: what eliminating the free x0
 * leaves of it is not a row that fixes x1 (at -1, hence "infeasible"). */
static void cancelled_entries_vanish(void)
{
    const small_lp lp = {2, 2, {0, 2, 4}, {0, 1, 0, 1}, {0.3, 3, 0.1, 1},
                         {-0.1, -1}, {0, 1}};
    const double lower[] = {-INFINITY, 0};
    arc_presolve ps;
    CHECK(presolve(&lp, lower, NULL, 1, &ps) == ARC_PRESOLVE_REDUCED);
    const double none[1] = {0.0};
    double x[2];
    arc_postsolve(&ps, none, x);
    CHECK_NEAR(x[0], -1.0 / 3.0, 1e-15);
    CHECK(x[1] == 0.0);
    arc_presolve_free(&ps);
}

static void disabled_rules_copy_problem(void)
{
    const small_lp lp = {2, 3, {0, 1, 3, 4}, {0, 0, 1, 1}, {1, 2, 3, 4},
                         {1, 2}, {1, 2, 3}};
    arc_presolve ps;
    CHECK(presolve(&lp, NULL, NULL, 0, &ps) == ARC_PRESOLVE_REDUCED);
    CHECK(ps.reduced.a.m == 2 && ps.reduced.a.n == 3);
    CHECK(ps.reduced.a.p[3] == 4 && ps.reduced.a.x[2] == 3.0);
    arc_presolve_free(&ps);
}

int main(void)
{
    singleton_substitution_and_bound();
    bounds_are_normalised();
    forcing_row_fixes_zero();
    detects_infeasible_and_unbounded();
    cancelled_entries_vanish();
    disabled_rules_copy_problem();
    return check_report("presolve");
}
