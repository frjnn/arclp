/* Iteration driver: stopping rules and safeguards shared by all rules. */
#include "step_rule.h"

#include "util/vec.h"

#include <stdio.h>

static const double DIVERGENCE_LIMIT = 1e10;
static const double GROWTH_FACTOR = 10.0;
static const double GROWTH_FLOOR = 1e-6;
static const double INFEASIBLE_RESIDUAL = 0.1;

typedef struct {
    double rb, rc;
} residual_norms;

static residual_norms norms(const arc_ipm *ipm)
{
    const residual_norms r = {arc_norm2(ipm->m, ipm->rb),
                              arc_norm2(ipm->n, ipm->rc)};
    return r;
}

/* LIPSOL-style relative measure used by arcLP.m. */
static double optimality_error(const arc_ipm *ipm, residual_norms r)
{
    const double nb = ipm->norm_b;
    const double nc = ipm->norm_c;
    const double gap = fabs(arc_dot(ipm->n, ipm->c, ipm->x)
                            - arc_dot(ipm->m, ipm->b, ipm->y));
    return r.rb / arc_fmax(1.0, nb) + r.rc / arc_fmax(1.0, nc)
         + gap / arc_fmax(nb, arc_fmax(1.0, nc));
}

static int grew(double before, double after)
{
    return GROWTH_FACTOR * before < after && after > GROWTH_FLOOR;
}

/* arcLP.m restores x and s only; restoring y too keeps the final iterate
 * consistent, so the residuals classify() sees are meaningful. */
static void restore_previous(arc_ipm *ipm)
{
    arc_copy(ipm->n, ipm->x_old, ipm->x);
    arc_copy(ipm->n, ipm->s_old, ipm->s);
    arc_copy(ipm->m, ipm->y_old, ipm->y);
    ipm->mu = arc_dot(ipm->n, ipm->x, ipm->s) / ipm->n;
}

static void remember(arc_ipm *ipm)
{
    arc_copy(ipm->n, ipm->x, ipm->x_old);
    arc_copy(ipm->n, ipm->s, ipm->s_old);
    arc_copy(ipm->m, ipm->y, ipm->y_old);
}

static int prepare_scaling(arc_ipm *ipm)
{
    for (int j = 0; j < ipm->n; j++) {
        ipm->d[j] = ipm->x[j] / ipm->s[j];
    }
    return arc_normal_factor(ipm->solver, ipm->d);
}

/* Blow-up of x signals an unbounded primal, of y an infeasible one. */
static int diverged(const arc_ipm *ipm, arc_ipm_result *result, double rb)
{
    if (arc_max(ipm->n, ipm->x) > DIVERGENCE_LIMIT) {
        result->outcome = ARC_IPM_UNBOUNDED;
    } else if (arc_max(ipm->m, ipm->y) > DIVERGENCE_LIMIT) {
        result->outcome = ARC_IPM_PRIMAL_INFEASIBLE;
    } else {
        return 0;
    }
    result->infeasibility = rb;
    return 1;
}

static void trace(const arc_ipm *ipm, int k, double error)
{
    printf("%4d  mu %.3e  |rb| %.3e  |rc| %.3e  err %.3e\n", k, ipm->mu,
           arc_norm2(ipm->m, ipm->rb), arc_norm2(ipm->n, ipm->rc), error);
}

/* Non-finite iterates end as a numerical error (arcLP.m would report them
 * optimal: every comparison with NaN is false).  arcLP.m only tests the
 * primal residual.  In addition, an unconverged dual residual together
 * with a primal objective far below the dual one -- weak duality violated,
 * the objective running off to -infinity -- is reported as unbounded. */
static arc_ipm_outcome classify(arc_ipm *ipm, arc_ipm_result *result,
                                const arc_ipm_settings *settings)
{
    arc_ipm_residuals(ipm);
    const residual_norms r = norms(ipm);
    result->infeasibility = r.rb;
    if (!(isfinite(r.rb) && isfinite(r.rc) && isfinite(ipm->mu))) {
        return ARC_IPM_NUMERICAL_ERROR;
    }
    if (result->iterations == settings->max_iterations) {
        return ARC_IPM_ITERATION_LIMIT;
    }
    if (r.rb > INFEASIBLE_RESIDUAL) {
        return ARC_IPM_PRIMAL_INFEASIBLE;
    }
    const double dual = r.rc / arc_fmax(1.0, ipm->norm_c);
    const double by = arc_dot(ipm->m, ipm->b, ipm->y);
    const double below = by - arc_dot(ipm->n, ipm->c, ipm->x);
    const int unbounded = dual > INFEASIBLE_RESIDUAL
                       && below > INFEASIBLE_RESIDUAL * arc_fmax(1.0, fabs(by));
    return unbounded ? ARC_IPM_UNBOUNDED : ARC_IPM_CONVERGED;
}

arc_ipm_result arc_ipm_run(arc_ipm *ipm, arc_step_rule *rule,
                           const arc_ipm_settings *settings)
{
    arc_ipm_result result = {ARC_IPM_FAILURE, 0, 0.0};
    residual_norms now = norms(ipm);
    remember(ipm);
    while (optimality_error(ipm, now) > settings->tolerance
           && result.iterations < settings->max_iterations) {
        if (settings->verbose) {
            trace(ipm, result.iterations, optimality_error(ipm, now));
        }
        if (prepare_scaling(ipm) != 0 || rule->ops->step(rule, ipm) != 0) {
            return result;
        }
        const residual_norms before = now;
        arc_ipm_residuals(ipm);
        now = norms(ipm);
        if (diverged(ipm, &result, now.rb)) {
            return result;
        }
        result.iterations++;
        if (rule->ops->stalled(rule, ipm)) {
            break;
        }
        if (grew(before.rb, now.rb) || grew(before.rc, now.rc)) {
            restore_previous(ipm);
            break;
        }
        remember(ipm);
    }
    result.outcome = classify(ipm, &result, settings);
    return result;
}
