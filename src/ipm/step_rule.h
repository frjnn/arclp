/*
 * Strategy that turns a factorized iterate into the next one.  The driver
 * factorizes A diag(x./s) A' before calling step(); implementations reuse
 * that factorization for all their right-hand sides.
 */
#ifndef ARC_STEP_RULE_H
#define ARC_STEP_RULE_H

#include "ipm.h"

typedef struct arc_step_rule arc_step_rule;

typedef struct {
    /* Replaces (x, y, s) by the next iterate; returns 0 on success. */
    int (*step)(arc_step_rule *self, arc_ipm *ipm);
    /* Rule-specific stagnation test, evaluated after the residual update. */
    int (*stalled)(const arc_step_rule *self, const arc_ipm *ipm);
    void (*destroy)(arc_step_rule *self);
} arc_step_rule_ops;

struct arc_step_rule {
    const arc_step_rule_ops *ops;
};

/* Yang's arc-search rule (Algorithm 5.1 with Algorithm 6.1). */
arc_step_rule *arc_search_rule_create(int m, int n);

/* Mehrotra's predictor-corrector rule as implemented in mehrotra.m. */
arc_step_rule *arc_mehrotra_rule_create(int m, int n);

static inline void arc_step_rule_destroy(arc_step_rule *rule)
{
    if (rule != 0) {
        rule->ops->destroy(rule);
    }
}

/* Runs iterations from the current iterate until a stopping rule fires. */
arc_ipm_result arc_ipm_run(arc_ipm *ipm, arc_step_rule *rule,
                           const arc_ipm_settings *settings);

#endif
