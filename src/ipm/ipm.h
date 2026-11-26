/*
 * Infeasible primal-dual interior-point iteration for
 *     min c'x  s.t.  A x = b,  x + w = u,  x, w >= 0   (bounded columns)
 *     max b'y - u'z  s.t.  A'y + s - z = c,  s, z >= 0,
 * with x >= 0 only on plain columns and x free, s = 0 on free ones.
 *
 * Columns are ordered free | plain | bounded, and x, s hold w, z behind
 * them, so the complementary pairs are the contiguous range
 * [n_free, n_total):
 *     x:  free | plain | bounded | w
 *     s:   0   |   s   |    s    | z
 * The driver owns the iterate and the stopping logic; how the next iterate
 * is produced is delegated to an arc_step_rule.
 */
#ifndef ARC_IPM_H
#define ARC_IPM_H

#include "sparse/lp.h"
#include "util/arena.h"
#include "normal_solver.h"

typedef struct
{
    const arc_csc *a;
    const double *b, *c, *upper;
    int m, n;
    int n_free, n_bounded;
    int n_total;                   /* n + n_bounded                            */
    int n_pairs;                   /* complementary pairs: n_total - n_free    */
    const double *row_unit;        /* 1 / row and column scaling, NULL when    */
    const double *col_unit;        /* the problem is not scaled                */
    const double *bound_unit;      /* scaling of the bounded columns           */
    arc_normal_solver *solver;
    double *x, *y, *s;             /* current iterate (y is lambda)            */
    double *rb, *rc;               /* rb = A x - b,  rc = A'y + s - z - c      */
    double *rw;                    /* x + w - u                                */
    double *d;                     /* scaling of the normal equations          */
    double mu;                     /* x's over the pairs / n_pairs             */
    double norm_b, norm_c;         /* of (b, u) and c in original units        */
    double *x_best, *y_best, *s_best; /* iterate with the smallest error     */
    double *work_n[2];             /* scratch for the start point              */
    double *work_m[2];
    arc_arena arena;
} arc_ipm;

typedef struct
{
    double tolerance;
    int max_iterations;
    int verbose;
} arc_ipm_settings;

typedef enum
{
    ARC_IPM_CONVERGED,
    ARC_IPM_PRIMAL_INFEASIBLE,
    ARC_IPM_UNBOUNDED,
    ARC_IPM_STALLED, /* stopped improving short of the tolerance */
    ARC_IPM_ITERATION_LIMIT,
    ARC_IPM_NUMERICAL_ERROR,
    ARC_IPM_FAILURE
} arc_ipm_outcome;

typedef struct
{
    double primal; /* ||A x - b|| / (1 + ||b||)           */
    double dual;   /* ||A'y + s - c|| / (1 + ||c||)       */
    double gap;    /* |c'x - b'y| / (1 + the smaller one) */
} arc_ipm_errors;

typedef struct
{
    arc_ipm_outcome outcome;
    int iterations;
    arc_ipm_errors errors; /* of the returned iterate */
} arc_ipm_result;

int arc_ipm_init(arc_ipm *ipm, const arc_lp *lp, arc_normal_solver *solver);
void arc_ipm_free(arc_ipm *ipm);

/* Runs one statement per column j, chosen by the column's kind. */
#define ARC_BY_KIND(ipm, j, free_column, plain_column, bounded_column)      \
    do {                                                                    \
        for (int j = 0; j < (ipm)->n_free; j++) {                           \
            free_column;                                                    \
        }                                                                   \
        for (int j = (ipm)->n_free; j < arc_ipm_bounded(ipm); j++) {        \
            plain_column;                                                   \
        }                                                                   \
        for (int j = arc_ipm_bounded(ipm); j < (ipm)->n; j++) {             \
            bounded_column;                                                 \
        }                                                                   \
    } while (0)

/* Bounded columns are [arc_ipm_bounded(ipm), n). */
static inline int arc_ipm_bounded(const arc_ipm *ipm)
{
    return ipm->n - ipm->n_bounded;
}

/* Index of the pair (w, z) of bounded column j. */
static inline int arc_ipm_slack(const arc_ipm *ipm, int j)
{
    return j + ipm->n_bounded;
}

/* rb, rc, rw and mu from the current (x, y, s). */
void arc_ipm_residuals(arc_ipm *ipm);

/* Keeps the iterate, and puts the kept one back with its residuals. */
void arc_ipm_remember(arc_ipm *ipm);
void arc_ipm_restore(arc_ipm *ipm);

/* Better of Mehrotra's and Lustig's starting points (as in arcLP.m). */
int arc_ipm_start_point(arc_ipm *ipm, double spread);

#endif
