/* Each presolve rule on a small problem, plus postsolve consistency. */
#include "check.h"

#include "presolve/presolve.h"

typedef struct {
    int m, n;
    int p[8];
    int i[16];
    double x[16];
    double b[4];
    double c[6];
} small_lp;

static arc_csc view(const small_lp *lp)
{
    const arc_csc a = {lp->m, lp->n, lp->p, lp->i, lp->x};
    return a;
}

/* Rule 5 fixes x0 = 2 from row 0; rule 9 then substitutes x1 = 3 + x2
 * out of row 1 (its only positive coefficient in a b > 0 row), leaving
 * the single row x2 + x3 = 1. */
static void singleton_and_substitution(void)
{
    /* row0: 2 x0 = 4;  row1: x0 + x1 - x2 = 5;  row2: x2 + x3 = 1 */
    const small_lp lp = {3, 4, {0, 2, 3, 5, 6}, {0, 1, 1, 1, 2, 2},
                         {2, 1, 1, -1, 1, 1}, {4, 5, 1}, {1, 1, 1, 1}};
    const arc_csc a = view(&lp);
    arc_presolve ps;
    CHECK(arc_presolve_run(&ps, &a, lp.b, lp.c, 1) == ARC_PRESOLVE_REDUCED);
    CHECK(ps.reduced.m == 1 && ps.reduced.n == 2);
    CHECK_NEAR(ps.reduced_c[0], 2.0, 1e-15);   /* c2 absorbed c1 */
    const double xr[2] = {0.5, 0.5};
    double x[4];
    arc_postsolve(&ps, xr, x);
    CHECK_NEAR(x[0], 2.0, 1e-15);
    CHECK_NEAR(x[1], 3.5, 1e-15);
    CHECK_NEAR(x[0] + x[1] - x[2], 5.0, 1e-14);
    arc_presolve_free(&ps);
}

static void forcing_row_fixes_zero(void)
{
    /* row0: x0 + x1 = 0 forces x0 = x1 = 0;  row1: x1 + x2 = 3 */
    const small_lp lp = {2, 3, {0, 1, 3, 4}, {0, 0, 1, 1}, {1, 1, 1, 1},
                         {0, 3}, {1, 1, 1}};
    const arc_csc a = view(&lp);
    arc_presolve ps;
    CHECK(arc_presolve_run(&ps, &a, lp.b, lp.c, 1) == ARC_PRESOLVE_REDUCED);
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
    arc_csc a = view(&infeasible);
    arc_presolve ps;
    CHECK(arc_presolve_run(&ps, &a, infeasible.b, infeasible.c, 1)
          == ARC_PRESOLVE_INFEASIBLE);
    arc_presolve_free(&ps);
    /* Column 1 is empty with negative cost: unbounded (rule 3). */
    const small_lp unbounded = {1, 2, {0, 1, 1}, {0}, {1}, {1}, {1, -1}};
    a = view(&unbounded);
    CHECK(arc_presolve_run(&ps, &a, unbounded.b, unbounded.c, 1)
          == ARC_PRESOLVE_UNBOUNDED);
    arc_presolve_free(&ps);
}

static void disabled_rules_copy_problem(void)
{
    const small_lp lp = {2, 3, {0, 1, 3, 4}, {0, 0, 1, 1}, {1, 2, 3, 4},
                         {1, 2}, {1, 2, 3}};
    const arc_csc a = view(&lp);
    arc_presolve ps;
    CHECK(arc_presolve_run(&ps, &a, lp.b, lp.c, 0) == ARC_PRESOLVE_REDUCED);
    CHECK(ps.reduced.m == 2 && ps.reduced.n == 3);
    CHECK(ps.reduced.p[3] == 4 && ps.reduced.x[2] == 3.0);
    arc_presolve_free(&ps);
}

int main(void)
{
    singleton_and_substitution();
    forcing_row_fixes_zero();
    detects_infeasible_and_unbounded();
    disabled_rules_copy_problem();
    return check_report("presolve");
}
