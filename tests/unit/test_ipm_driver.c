/* The start point, and the driver against fake step rules: non-finite
 * iterates must end in a numerical error, never in "converged", and a
 * falling mu must not be taken for a stall. */
#include "check.h"

#include "chol/chol_solver.h"
#include "ipm/convergence.h"
#include "ipm/step_rule.h"
#include "small_lp.h"
#include "util/vec.h"

#include <string.h>

typedef struct
{
    arc_step_rule base;
    int steps;
} poison_rule;

static void poison_step(arc_step_rule *self, arc_ipm *ipm)
{
    ((poison_rule *)self)->steps++;
    ipm->x[0] = NAN;
}

static int never_stalls(const arc_step_rule *self, const arc_ipm *ipm)
{
    (void)self;
    (void)ipm;
    return 0;
}

static const arc_step_rule_ops POISON = {poison_step, never_stalls, NULL,
                                         NULL, 1.0, sizeof(poison_rule)};

static void scale(int n, double factor, double *v)
{
    for (int i = 0; i < n; i++)
    {
        v[i] *= factor;
    }
}

/* Shrinks x and y alike: mu falls by a fifth per step, and so do both
 * objectives. */
static void shrink_step(arc_step_rule *self, arc_ipm *ipm)
{
    (void)self;
    scale(ipm->n_total, 0.8, ipm->x);
    scale(ipm->m, 0.8, ipm->y);
}

static const arc_step_rule_ops SHRINK = {shrink_step, never_stalls, NULL,
                                         NULL, 1.0, sizeof(arc_step_rule)};

/* The small problem with given b and c, at its start point.  The iteration
 * and the solver keep pointers into lp. */
typedef struct
{
    arc_lp lp;
    arc_ipm ipm;
    arc_normal_solver *solver;
} fixture;

static void start(fixture *f, double *b, double *c)
{
    const arc_lp lp = {{SMALL_M, SMALL_N, SMALL_P, SMALL_I, SMALL_V},
                       NULL, b, c, NULL, 0, 0};
    f->lp = lp;
    f->solver = arc_chol_solver_create(&f->lp.a);
    CHECK(f->solver != NULL
          && arc_ipm_init(&f->ipm, &f->lp, f->solver) == 0);
    arc_ipm_set_units(&f->ipm, NULL, NULL, NULL);
    CHECK(arc_ipm_start_point(&f->ipm, 1.0) == 0);
}

static void stop(fixture *f)
{
    arc_ipm_free(&f->ipm);
    arc_normal_destroy(f->solver);
}

static void poisoned_iterates_are_errors(double *b, double *c)
{
    fixture f;
    start(&f, b, c);
    poison_rule rule = {{&POISON, 0, 0, {0}}, 0};
    const arc_ipm_settings settings = {1e-8, 100, 0};
    const arc_ipm_result r = arc_ipm_run(&f.ipm, &rule.base, &settings);
    CHECK(r.outcome == ARC_IPM_NUMERICAL_ERROR);
    CHECK(rule.steps == 1);
    stop(&f);
}

/* A falling mu is progress even when the objectives fall as fast, which
 * keeps mu relative to them constant. */
static void falling_mu_is_progress(double *b, double *c)
{
    fixture f;
    start(&f, b, c);
    arc_step_rule rule = {&SHRINK, 0, 0, {0}};
    const arc_ipm_settings settings = {1e-8, 40, 0};
    const arc_ipm_result r = arc_ipm_run(&f.ipm, &rule, &settings);
    CHECK(r.outcome == ARC_IPM_ITERATION_LIMIT);
    stop(&f);
}

/* Neither a zero b nor a constant cost may start on the boundary. */
static void start_is_interior(double *b, double *c)
{
    fixture f;
    start(&f, b, c);
    CHECK(arc_min(SMALL_N, f.ipm.x) > 1e-6);
    CHECK(arc_min(SMALL_N, f.ipm.s) > 1e-6);
    stop(&f);
}

int main(void)
{
    double b[SMALL_M], c[SMALL_N];
    double zero[SMALL_N] = {0};
    double constant[] = {2, 4, 1, 1};   /* both rows added up */
    memcpy(b, SMALL_B, sizeof(b));
    memcpy(c, SMALL_C, sizeof(c));
    poisoned_iterates_are_errors(b, c);
    double big_b[SMALL_M], big_c[SMALL_N];   /* objectives far above 1 */
    memcpy(big_b, SMALL_B, sizeof(big_b));
    memcpy(big_c, SMALL_C, sizeof(big_c));
    scale(SMALL_M, 1e6, big_b);
    scale(SMALL_N, 1e6, big_c);
    falling_mu_is_progress(big_b, big_c);
    start_is_interior(zero, c);
    start_is_interior(b, constant);
    start_is_interior(b, zero);
    start_is_interior(zero, zero);
    return check_report("ipm driver");
}
