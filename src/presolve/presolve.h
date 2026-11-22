/*
 * Presolve of arcLP.m: the reduction rules 1, 3, 5, 7 and 9 applied until
 * none fires, the reduced problem handed to the solver, and the postsolve
 * that maps its solution back to the original variables.
 */
#ifndef ARC_PRESOLVE_H
#define ARC_PRESOLVE_H

#include "sparse/csc.h"
#include "ps_matrix.h"

typedef enum
{
    ARC_PRESOLVE_REDUCED,
    ARC_PRESOLVE_INFEASIBLE,
    ARC_PRESOLVE_UNBOUNDED,
    ARC_PRESOLVE_OUT_OF_MEMORY
} arc_presolve_status;

/* x_p = rhs - sum coef_k x_{col_k}: a variable eliminated by rule 9. */
typedef struct
{
    int *pivot;
    double *rhs;
    int *start; /* count + 1 offsets into term_col / term_coef */
    int *term_col;
    double *term_coef;
    int count, capacity, terms, term_capacity;
} arc_substitutions;

typedef struct
{
    int n;           /* columns of the original problem       */
    ps_matrix a;     /* working matrix, freed once reduced is built */
    double *b, *c;   /* working right-hand side and costs     */
    double *fixed_x; /* values of the eliminated variables    */
    arc_substitutions subs;
    int *scatter;  /* n, -1 outside use                     */
    int *row_cols; /* n, scratch for one row                */
    double *row_vals;
    double abs_ratio; /* max |a_ij| / min |a_ij| of the input  */
    /* Reduced problem */
    arc_csc reduced; /* view of the storage below            */
    int *reduced_p, *reduced_i;
    double *reduced_x;
    double *reduced_b, *reduced_c;
    int *reduced_col; /* original column of each reduced one   */
} arc_presolve;

arc_presolve_status arc_presolve_run(arc_presolve *ps, const arc_csc *a,
                                     const double *b, const double *c,
                                     int apply_rules);

/* x (original length) from the reduced solution xr. */
void arc_postsolve(const arc_presolve *ps, const double *xr, double *x);

void arc_presolve_free(arc_presolve *ps);

/* Individual rules; each returns 1 when it changed the problem, 0 when not,
 * or a negative arc_presolve_status when it proved the outcome. */
int arc_rule_empty_rows(arc_presolve *ps);
int arc_rule_empty_cols(arc_presolve *ps);
int arc_rule_singleton_rows(arc_presolve *ps);
int arc_rule_forcing_rows(arc_presolve *ps);
int arc_rule_implied_free(arc_presolve *ps);

/* Builds the reduced CSC problem from the surviving rows and columns. */
int arc_extract_reduced(arc_presolve *ps);

#endif
