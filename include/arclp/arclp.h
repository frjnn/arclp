/*
 * arclp — arc-search infeasible interior-point method for linear programs
 *
 *     minimize c'x   subject to   A x = b,  lower <= x <= upper.
 *
 * C implementation of Y. Yang's arcLP (Numerical Algorithms 79, 2018,
 * Algorithm 5.1 with the step selection of Algorithm 6.1).  Mehrotra's
 * predictor-corrector method is provided as an alternative step rule so the
 * two can be compared on identical presolve, start point and stopping rules.
 */
#ifndef ARCLP_H
#define ARCLP_H

#include <arclp/export.h>

#ifdef __cplusplus
extern "C"
{
#endif

    /* Read-only view of a problem; A is in compressed sparse column form. */
    typedef struct
    {
        int m;                /* rows of A                                  */
        int n;                /* columns of A                               */
        const int *col_ptr;   /* n + 1 column pointers                      */
        const int *row_idx;   /* row index of each nonzero                  */
        const double *values; /* value of each nonzero                      */
        const double *b;      /* right-hand side, length m                  */
        const double *c;      /* objective, length n                        */
        const double *lower;  /* bounds, length n, -INFINITY allowed; NULL   */
        const double *upper;  /* means 0 for lower, +INFINITY for upper      */
    } arclp_problem;

    typedef enum
    {
        ARCLP_OPTIMAL = 0,
        ARCLP_INFEASIBLE = 1,
        ARCLP_UNBOUNDED = 2, /* feasible, with a cost that falls without
                              * limit                                  */
        ARCLP_MAX_ITERATIONS = 3,
        ARCLP_INVALID_INPUT = 4, /* malformed matrix or non-finite data   */
        ARCLP_OUT_OF_MEMORY = 5,
        ARCLP_NUMERICAL_ERROR = 6, /* the iterates stopped being finite    */
        ARCLP_STALLED = 7 /* stopped improving short of the tolerance: x is
                           * the best iterate, see the errors in info     */
    } arclp_status;

    typedef enum
    {
        ARCLP_METHOD_ARC_SEARCH = 0, /* Yang's arcLP (default)        */
        ARCLP_METHOD_MEHROTRA = 1    /* Mehrotra predictor-corrector  */
    } arclp_method;

    typedef struct
    {
        arclp_method method;
        double tolerance; /* relative residual + gap stopping threshold  */
        int max_iterations;
        int presolve; /* nonzero: apply the reduction rules of arcLP */
        int verbose;  /* nonzero: print one line per iteration       */
    } arclp_options;

    typedef struct
    {
        arclp_status status;
        int iterations;       /* also those of the zero-cost solve that
                               * checks a ray for feasibility               */
        double objective;     /* c'x of the original problem                 */
        double infeasibility; /* ||A x - b|| of the reduced problem          */
        int reduced_rows;     /* size of the problem after presolve          */
        int reduced_cols;
        long factor_nonzeros; /* stored entries of the Cholesky factor       */
        int tiny_pivots;      /* pivots replaced to cope with rank loss      */
        double seconds_presolve;
        double seconds_analyse; /* ordering and symbolic factorization        */
        double seconds_iterate; /* interior-point iterations                  */
        double seconds_total;
        double primal_error; /* relative errors of the returned iterate on  */
        double dual_error;   /* the reduced problem; all are at most        */
        double gap;          /* options.tolerance when status is optimal    */
    } arclp_info;

    ARCLP_API void arclp_default_options(arclp_options *options);

    /* Solves the problem; x (length n) receives the primal solution.
     * options may be NULL for the defaults; info may be NULL. */
    ARCLP_API arclp_status arclp_solve(const arclp_problem *problem,
                                       const arclp_options *options,
                                       double *x, arclp_info *info);

    ARCLP_API const char *arclp_status_name(arclp_status status);

#ifdef __cplusplus
}
#endif

#endif
