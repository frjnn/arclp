#include "solve_reduced.h"

#include "chol/chol_solver.h"
#include "ipm/convergence.h"
#include "ipm/step_rule.h"
#include "presolve/presolve.h"
#include "scale/scaling.h"
#include "transform/dualize.h"
#include "transform/stretch.h"
#include "util/clock.h"
#include "util/vec.h"

#include <stdlib.h>

static const arc_step_rule_ops *const RULES[] = {
    [ARCLP_METHOD_ARC_SEARCH] = &ARC_SEARCH_OPS,
    [ARCLP_METHOD_MEHROTRA] = &ARC_MEHROTRA_OPS,
};

static const double EMPTY_RESIDUAL = 1e-8;

/* Cheapest value of a column on its own: 0, upper, or unbounded. */
static int solve_separable(const arc_lp *lp, double *x)
{
    for (int j = 0; j < lp->a.n; j++)
    {
        const int k = j - arc_lp_first_bounded(lp);
        const int improves = lp->c[j] < 0.0 || (j < lp->n_free
                                                && lp->c[j] != 0.0);
        if (improves && k < 0)
        {
            return ARCLP_UNBOUNDED;
        }
        x[j] = improves ? lp->upper[k] : 0.0;
    }
    return ARCLP_OPTIMAL;
}

/* Problems the interior-point method cannot represent: no rows or no
 * columns left.  Returns -1 when the problem is not degenerate. */
static int solve_trivial(const arc_lp *lp, double *x)
{
    if (lp->a.n == 0)
    {
        return arc_norm2(lp->a.m, lp->b) <= EMPTY_RESIDUAL ? ARCLP_OPTIMAL
                                                           : ARCLP_INFEASIBLE;
    }
    return lp->a.m == 0 ? solve_separable(lp, x) : -1;
}

static arclp_status status_of(arc_ipm_outcome outcome)
{
    switch (outcome)
    {
    case ARC_IPM_CONVERGED: return ARCLP_OPTIMAL;
    case ARC_IPM_PRIMAL_INFEASIBLE: return ARCLP_INFEASIBLE;
    case ARC_IPM_UNBOUNDED: return ARCLP_UNBOUNDED;
    case ARC_IPM_STALLED: return ARCLP_STALLED;
    case ARC_IPM_ITERATION_LIMIT: return ARCLP_MAX_ITERATIONS;
    case ARC_IPM_NUMERICAL_ERROR:
    case ARC_IPM_FACTOR_FAILED: return ARCLP_NUMERICAL_ERROR;
    }
    return ARCLP_NUMERICAL_ERROR;
}

static arclp_status iterate(arc_ipm *ipm, arc_step_rule *rule,
                            const arclp_options *options, arclp_info *info)
{
    const arc_ipm_settings settings = {options->tolerance,
                                       options->max_iterations,
                                       options->verbose};
    if (arc_ipm_start_point(ipm, rule->ops->start_spread) != 0)
    {
        return ARCLP_OUT_OF_MEMORY;
    }
    const arc_ipm_result result = arc_ipm_run(ipm, rule, &settings);
    info->iterations = result.iterations;
    info->infeasibility = result.errors.primal * (1.0 + ipm->norm_b);
    info->primal_error = result.errors.primal;
    info->dual_error = result.errors.dual;
    info->gap = result.errors.gap;
    arc_normal_stats stats;
    ipm->solver->ops->stats(ipm->solver, &stats);
    info->factor_nonzeros = stats.factor_nonzeros;
    info->tiny_pivots = stats.tiny_pivots;
    return status_of(result.outcome);
}

/* x has a.n entries; y (a.m) may be NULL. */
typedef struct
{
    double *x, *y;
} solution;

static arclp_status solve_scaled(const arc_lp *lp, const arc_scaling *scale,
                                 const arclp_options *options,
                                 solution out, arclp_info *info)
{
    const double start = arc_seconds();
    arc_normal_solver *solver = arc_chol_solver_create(&lp->a);
    arc_step_rule *rule =
        arc_step_rule_create(RULES[options->method], lp->a.m,
                             lp->a.n + lp->n_bounded);
    arc_ipm ipm;
    const int ready = solver != NULL && rule != NULL
                   && arc_ipm_init(&ipm, lp, solver) == 0;
    const double analysed = arc_seconds();
    info->seconds_analyse = analysed - start;
    arclp_status status = ARCLP_OUT_OF_MEMORY;
    if (ready)
    {
        arc_ipm_set_units(&ipm, scale->row_unit, scale->col_unit,
                          scale->col);
        status = iterate(&ipm, rule, options, info);
        arc_copy(lp->a.n, ipm.x, out.x);
        if (out.y != NULL)
        {
            arc_copy(lp->a.m, ipm.y, out.y);
        }
        arc_ipm_free(&ipm);
    }
    info->seconds_iterate = arc_seconds() - analysed;
    arc_step_rule_destroy(rule);
    arc_normal_destroy(solver);
    return status;
}

/* Scales lp in place and solves it. */
static arclp_status solve_direct(arc_lp *lp, const arclp_options *options,
                                 solution out, arclp_info *info)
{
    arc_scaling scale;
    if (arc_scale_problem(lp, &scale) != 0)
    {
        return ARCLP_OUT_OF_MEMORY;
    }
    const arclp_status status = solve_scaled(lp, &scale, options, out, info);
    arc_unscale(&scale, lp, out.x, out.y);
    arc_scaling_free(&scale);
    return status;
}

/* Stretches the dense columns, if any, and solves. */
static arclp_status solve_stretched(arc_lp *lp, const arclp_options *options,
                                    solution out, arclp_info *info)
{
    arc_stretched s;
    const int stretched = arc_stretch(lp, &s);
    if (stretched < 0)
    {
        return ARCLP_OUT_OF_MEMORY;
    }
    if (stretched == 0)
    {
        return solve_direct(lp, options, out, info);
    }
    const size_t size = (size_t)s.lp.a.n + (size_t)s.lp.a.m;
    double *block = malloc(sizeof(double) * size);
    arclp_status status = ARCLP_OUT_OF_MEMORY;
    if (block != NULL)
    {
        const solution inner = {block, block + s.lp.a.n};
        status = solve_direct(&s.lp, options, inner, info);
        for (int j = 0; j < lp->a.n; j++)
        {
            out.x[j] = inner.x[s.first[j]];
        }
        if (out.y != NULL)
        {
            arc_copy(lp->a.m, inner.y, out.y);
        }
    }
    free(block);
    arc_stretched_free(&s);
    return status;
}

static arclp_status dual_status(arclp_status s)
{
    return s == ARCLP_INFEASIBLE  ? ARCLP_UNBOUNDED
         : s == ARCLP_UNBOUNDED ? ARCLP_INFEASIBLE
                                : s;
}

/* Solves the dual; its row multipliers give x. */
static arclp_status solve_dual(const arc_lp *lp, const arclp_options *options,
                               double *x, arclp_info *info)
{
    arc_dual d;
    if (arc_dual_create(lp, &d) != 0)
    {
        return ARCLP_OUT_OF_MEMORY;
    }
    arc_presolve ps;
    arclp_status status = ARCLP_OUT_OF_MEMORY;
    const arc_presolve_status bounds =
        arc_presolve_run(&ps, &d.a, d.b, d.c, d.lower, d.upper, 0);
    double *block = malloc(sizeof(double) * ((size_t)d.a.m + (size_t)d.a.n + 1));
    if (bounds == ARC_PRESOLVE_REDUCED && block != NULL)
    {
        const solution inner = {block, block + d.a.n};
        status = dual_status(solve_stretched(&ps.reduced, options, inner, info));
        arc_dual_recover(&d, inner.y, x);
    }
    else if (bounds != ARC_PRESOLVE_OUT_OF_MEMORY && block != NULL)
    {
        status = bounds == ARC_PRESOLVE_INFEASIBLE ? ARCLP_UNBOUNDED
                                                   : ARCLP_INFEASIBLE;
    }
    const double primal = info->primal_error;
    info->primal_error = info->dual_error;
    info->dual_error = primal;
    free(block);
    arc_presolve_free(&ps);
    arc_dual_free(&d);
    return status;
}

arclp_status arc_solve_reduced(arc_lp *lp, const arclp_options *options,
                               double *x, arclp_info *info)
{
    const int trivial = solve_trivial(lp, x);
    if (trivial >= 0)
    {
        return (arclp_status)trivial;
    }
    if (arc_dual_is_cheaper(lp))
    {
        return solve_dual(lp, options, x, info);
    }
    const solution out = {x, NULL};
    return solve_stretched(lp, options, out, info);
}
