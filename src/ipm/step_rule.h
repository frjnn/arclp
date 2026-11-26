/*
 * Strategy that turns a factorized iterate into the next one.  The driver
 * factorizes A diag(x./s) A' before calling step(); implementations reuse
 * that factorization for all their right-hand sides.
 */
#ifndef ARC_STEP_RULE_H
#define ARC_STEP_RULE_H

#include "ipm.h"

typedef struct arc_step_rule arc_step_rule;

typedef struct
{
    /* Replaces (x, y, s) by the next iterate; returns 0 on success. */
    int (*step)(arc_step_rule *self, arc_ipm *ipm);
    /* Rule-specific stagnation test, evaluated after the residual update. */
    int (*stalled)(const arc_step_rule *self, const arc_ipm *ipm);
    /* Takes the rule's vectors from the arena; owner is the rule. */
    arc_arena_layout layout;
    /* Sets the fields that do not start at zero; may be NULL. */
    void (*init)(arc_step_rule *self);
    double start_spread; /* multiplies the start point's shifts */
    size_t size;         /* bytes of the concrete rule */
} arc_step_rule_ops;

/* Concrete rules start with this struct. */
struct arc_step_rule
{
    const arc_step_rule_ops *ops;
    int m, n; /* n counts the pair entries: n_total */
    arc_arena arena;
};

/* Yang's arc-search rule (Algorithm 5.1 with Algorithm 6.1). */
extern const arc_step_rule_ops ARC_SEARCH_OPS;

/* Mehrotra's predictor-corrector rule as implemented in mehrotra.m. */
extern const arc_step_rule_ops ARC_MEHROTRA_OPS;

/* A zero-initialized rule of the given kind, or NULL. */
arc_step_rule *arc_step_rule_create(const arc_step_rule_ops *ops, int m,
                                    int n);

void arc_step_rule_destroy(arc_step_rule *rule);

/* Runs iterations from the current iterate until a stopping rule fires. */
arc_ipm_result arc_ipm_run(arc_ipm *ipm, arc_step_rule *rule,
                           const arc_ipm_settings *settings);

#endif
