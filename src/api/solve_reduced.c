#include "solve_reduced.h"

#include "chol/chol_solver.h"
#include "ipm/step_rule.h"
#include "util/clock.h"
#include "util/vec.h"

typedef arc_step_rule *(*rule_factory)(int m, int n);

static const rule_factory RULES[] = {
    [ARCLP_METHOD_ARC_SEARCH] = arc_search_rule_create,
    [ARCLP_METHOD_MEHROTRA] = arc_mehrotra_rule_create,
};

static const double EMPTY_RESIDUAL = 1e-8;

/* Problems the interior-point method cannot represent: no rows or no
 * columns left.  Returns -1 when the problem is not degenerate. */
static int solve_trivial(const arc_csc *a, const double *b, const double *c,
                         double *x)
{
    if (a->n == 0) {
        return arc_norm2(a->m, b) <= EMPTY_RESIDUAL ? ARCLP_OPTIMAL
                                                    : ARCLP_INFEASIBLE;
    }
    if (a->m == 0) {
        for (int j = 0; j < a->n; j++) {
            x[j] = 0.0;
        }
        return arc_min(a->n, c) < 0.0 ? ARCLP_UNBOUNDED : ARCLP_OPTIMAL;
    }
    return -1;
}

static arclp_status status_of(arc_ipm_outcome outcome)
{
    switch (outcome) {
    case ARC_IPM_CONVERGED: return ARCLP_OPTIMAL;
    case ARC_IPM_PRIMAL_INFEASIBLE: return ARCLP_INFEASIBLE;
    case ARC_IPM_UNBOUNDED: return ARCLP_UNBOUNDED;
    case ARC_IPM_ITERATION_LIMIT: return ARCLP_MAX_ITERATIONS;
    case ARC_IPM_NUMERICAL_ERROR: return ARCLP_NUMERICAL_ERROR;
    default: return ARCLP_OUT_OF_MEMORY;
    }
}

static arclp_status iterate(arc_ipm *ipm, arc_step_rule *rule,
                            const arclp_options *options, arclp_info *info)
{
    const arc_ipm_settings settings = {options->tolerance,
                                       options->max_iterations,
                                       options->verbose};
    if (arc_ipm_start_point(ipm) != 0) {
        return ARCLP_OUT_OF_MEMORY;
    }
    const arc_ipm_result result = arc_ipm_run(ipm, rule, &settings);
    info->iterations = result.iterations;
    info->infeasibility = result.infeasibility;
    arc_normal_stats stats;
    ipm->solver->ops->stats(ipm->solver, &stats);
    info->factor_nonzeros = stats.factor_nonzeros;
    info->tiny_pivots = stats.tiny_pivots;
    return status_of(result.outcome);
}

arclp_status arc_solve_reduced(const arc_csc *a, const double *b,
                               const double *c, const arclp_options *options,
                               double *x, arclp_info *info)
{
    const int trivial = solve_trivial(a, b, c, x);
    if (trivial >= 0) {
        return (arclp_status)trivial;
    }
    const double start = arc_seconds();
    arc_normal_solver *solver = arc_chol_solver_create(a);
    arc_step_rule *rule = RULES[options->method](a->m, a->n);
    arc_ipm ipm;
    const int ready = solver != NULL && rule != NULL
                   && arc_ipm_init(&ipm, a, b, c, solver) == 0;
    const double analysed = arc_seconds();
    info->seconds_analyse = analysed - start;
    arclp_status status = ARCLP_OUT_OF_MEMORY;
    if (ready) {
        status = iterate(&ipm, rule, options, info);
        arc_copy(a->n, ipm.x, x);
        arc_ipm_free(&ipm);
    }
    info->seconds_iterate = arc_seconds() - analysed;
    arc_step_rule_destroy(rule);
    arc_normal_destroy(solver);
    return status;
}
