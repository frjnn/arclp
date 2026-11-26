/* Iteration driver: stopping rules and safeguards shared by all rules. */
#include "convergence.h"
#include "newton.h"
#include "step_rule.h"

#include "util/vec.h"

#include <stdio.h>

static const double PROGRESS = 0.9;

enum
{
    PATIENCE = 15 /* iterations without progress before giving up */
};

static int prepare_scaling(arc_ipm *ipm)
{
    arc_newton_scaling(ipm);
    return arc_normal_factor(ipm->solver, ipm->d);
}

static void trace(const arc_ipm *ipm, int k, arc_ipm_errors e)
{
    printf("%4d  mu %.3e  primal %.3e  dual %.3e  gap %.3e\n", k, ipm->mu,
           e.primal, e.dual, e.gap);
}

/* Progress means that the primal residual, the dual residual or mu falls
 * clearly; amounts below the tolerance no longer count.  They are tracked
 * one by one, because the largest of them can stand still while the others
 * converge.  mu counts as it is: relative to the objective it rises while
 * an objective passes near zero, even as mu falls. */
enum
{
    MEASURES = 3
};

typedef struct
{
    double best;                /* worst error of the remembered iterate */
    double reference[MEASURES]; /* each measure at its last improvement  */
    int waited;
} progress;

/* Keeps the best iterate; returns 1 after PATIENCE iterations in which no
 * measure improved. */
static int no_progress(progress *p, arc_ipm *ipm, arc_ipm_errors now,
                       double tolerance)
{
    const double measures[MEASURES] = {now.primal, now.dual, ipm->mu};
    const double floors[MEASURES] = {tolerance, tolerance,
                                     arc_ipm_mu_floor(ipm, tolerance)};
    const double error = arc_ipm_worst(now);
    if (error < p->best)
    {
        p->best = error;
        arc_ipm_remember(ipm);
    }
    p->waited++;
    for (int k = 0; k < MEASURES; k++)
    {
        if (arc_fmax(measures[k], floors[k]) < PROGRESS * p->reference[k])
        {
            p->reference[k] = measures[k];
            p->waited = 0;
        }
    }
    return p->waited >= PATIENCE;
}

/* Refinement starts once a step has clearly raised the primal residual:
 * exact steps never do. */
static const double RESIDUAL_GROWTH = 1.1;
static const double RESIDUAL_FLOOR = 0.01; /* of the tolerance */

static void watch_residual(arc_ipm *ipm, double primal, double tolerance,
                           double *lowest)
{
    if (primal > RESIDUAL_GROWTH * *lowest
        && primal > RESIDUAL_FLOOR * tolerance)
    {
        arc_normal_refine(ipm->solver, 1);
    }
    *lowest = arc_fmin(*lowest, primal);
}

static arc_ipm_outcome iterate(arc_ipm *ipm, arc_step_rule *rule,
                               const arc_ipm_settings *settings,
                               int *iterations)
{
    arc_ipm_errors now = arc_ipm_measure(ipm);
    progress p = {arc_ipm_worst(now), {INFINITY, INFINITY, INFINITY}, 0};
    double lowest = now.primal;
    arc_normal_refine(ipm->solver, 0);
    arc_ipm_remember(ipm);
    while (arc_ipm_worst(now) > settings->tolerance)
    {
        if (settings->verbose)
        {
            trace(ipm, *iterations, now);
        }
        if (*iterations == settings->max_iterations)
        {
            return ARC_IPM_ITERATION_LIMIT;
        }
        if (prepare_scaling(ipm) != 0)
        {
            return ARC_IPM_FACTOR_FAILED;
        }
        rule->ops->step(rule, ipm);
        (*iterations)++;
        arc_ipm_residuals(ipm);
        now = arc_ipm_measure(ipm);
        watch_residual(ipm, now.primal, settings->tolerance, &lowest);
        if (!(isfinite(arc_ipm_worst(now)) && isfinite(ipm->mu)))
        {
            return ARC_IPM_NUMERICAL_ERROR;
        }
        const arc_ipm_outcome proof = arc_ipm_certified(ipm);
        if (proof != ARC_IPM_STALLED)
        {
            return proof;
        }
        if (no_progress(&p, ipm, now, settings->tolerance)
            || rule->ops->stalled(rule, ipm))
        {
            return ARC_IPM_STALLED;
        }
    }
    return ARC_IPM_CONVERGED;
}

/* When the iterations end without converging, the best iterate is the one
 * returned; it may still meet the tolerance. */
arc_ipm_result arc_ipm_run(arc_ipm *ipm, arc_step_rule *rule,
                           const arc_ipm_settings *settings)
{
    arc_ipm_result result = {ARC_IPM_FACTOR_FAILED, 0, {0.0, 0.0, 0.0}};
    result.outcome = iterate(ipm, rule, settings, &result.iterations);
    if (result.outcome != ARC_IPM_CONVERGED)
    {
        arc_ipm_restore(ipm);
    }
    result.errors = arc_ipm_measure(ipm);
    return result;
}
