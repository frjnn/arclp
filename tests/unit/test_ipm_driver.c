/* The driver against a fake step rule: a rule that produces non-finite
 * iterates must end in a numerical error, never in "converged". */
#include "check.h"

#include "chol/chol_solver.h"
#include "ipm/step_rule.h"

typedef struct
{
    arc_step_rule base;
    int steps;
} poison_rule;

static int poison_step(arc_step_rule *self, arc_ipm *ipm)
{
    ((poison_rule *)self)->steps++;
    ipm->x[0] = NAN;
    return 0;
}

static int never_stalls(const arc_step_rule *self, const arc_ipm *ipm)
{
    (void)self;
    (void)ipm;
    return 0;
}

static void no_destroy(arc_step_rule *self)
{
    (void)self;
}

static const arc_step_rule_ops POISON = {poison_step, never_stalls,
                                         no_destroy};

int main(void)
{
    /* x1 + x2 + s1 = 4,  x1 + 3 x2 + s2 = 6 */
    const int p[] = {0, 2, 4, 5, 6};
    const int i[] = {0, 1, 0, 1, 0, 1};
    const double v[] = {1, 1, 1, 3, 1, 1}, b[] = {4, 6}, c[] = {-1, -2, 0, 0};
    const arc_csc a = {2, 4, p, i, v};
    arc_normal_solver *solver = arc_chol_solver_create(&a);
    arc_ipm ipm;
    CHECK(solver != NULL && arc_ipm_init(&ipm, &a, b, c, solver) == 0);
    CHECK(arc_ipm_start_point(&ipm) == 0);
    poison_rule rule = {{&POISON}, 0};
    const arc_ipm_settings settings = {1e-8, 100, 0};
    const arc_ipm_result r = arc_ipm_run(&ipm, &rule.base, &settings);
    CHECK(r.outcome == ARC_IPM_NUMERICAL_ERROR);
    CHECK(rule.steps == 1);
    arc_ipm_free(&ipm);
    arc_normal_destroy(solver);
    return check_report("ipm driver");
}
