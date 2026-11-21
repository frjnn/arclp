/*
 * Infeasible primal-dual interior-point iteration for
 *     min c'x  s.t.  A x = b, x >= 0     (dual: A'y + s = c, s >= 0).
 * The driver owns the iterate and the stopping logic; how the next iterate
 * is produced is delegated to an arc_step_rule.
 */
#ifndef ARC_IPM_H
#define ARC_IPM_H

#include "sparse/csc.h"
#include "util/arena.h"
#include "normal_solver.h"

typedef struct {
    const arc_csc *a;
    const double *b, *c;
    int m, n;
    arc_normal_solver *solver;
    double *x, *y, *s;        /* current iterate (y is lambda)            */
    double *rb, *rc;          /* rb = A x - b,  rc = A'y + s - c          */
    double *d;                /* scaling x ./ s of the normal equations   */
    double mu;                /* x's / n                                  */
    double norm_b, norm_c;    /* ||b||, ||c||, fixed for the whole solve  */
    double *x_old, *y_old, *s_old; /* previous iterate, restored when the
                                    * residuals grow                     */
    double *work_n[2];        /* scratch for the start point              */
    double *work_m[2];
    arc_arena arena;
} arc_ipm;

typedef struct {
    double tolerance;
    int max_iterations;
    int verbose;
} arc_ipm_settings;

typedef enum {
    ARC_IPM_CONVERGED,
    ARC_IPM_PRIMAL_INFEASIBLE,
    ARC_IPM_UNBOUNDED,
    ARC_IPM_ITERATION_LIMIT,
    ARC_IPM_NUMERICAL_ERROR,
    ARC_IPM_FAILURE
} arc_ipm_outcome;

typedef struct {
    arc_ipm_outcome outcome;
    int iterations;
    double infeasibility;     /* ||A x - b|| */
} arc_ipm_result;

int arc_ipm_init(arc_ipm *ipm, const arc_csc *a, const double *b,
                 const double *c, arc_normal_solver *solver);
void arc_ipm_free(arc_ipm *ipm);

/* rb, rc and mu from the current (x, y, s). */
void arc_ipm_residuals(arc_ipm *ipm);

/* Better of Mehrotra's and Lustig's starting points (as in arcLP.m). */
int arc_ipm_start_point(arc_ipm *ipm);

#endif
