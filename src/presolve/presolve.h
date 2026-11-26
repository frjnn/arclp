/*
 * Presolve: the reduction rules 1, 3, 5, 7 and 9 of arcLP.m, extended to
 * bounded and free variables, and the bounds implied by single-variable
 * inequalities, applied until none fires; the reduced problem handed to the
 * solver; and the postsolve that maps its solution back.
 */
#ifndef ARC_PRESOLVE_H
#define ARC_PRESOLVE_H

#include "sparse/lp.h"
#include "ps_matrix.h"
#include "util/arena.h"

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
    double *upper;   /* 0 <= x <= upper unless is_free        */
    unsigned char *is_free;
    arc_substitutions subs;
    int *scatter;  /* n, -1 outside use                     */
    int *row_cols; /* n, scratch for one row                */
    double *row_vals;
    double abs_ratio; /* max |a_ij| / min |a_ij| of the input  */
    arc_lp reduced;  /* the problem handed to the solver      */
    int *reduced_p, *reduced_i; /* writable storage behind it  */
    int reduced_nnz;
    int *reduced_col; /* original column of each reduced one   */
    arc_arena arena;  /* storage of the reduced problem        */
} arc_presolve;

typedef struct
{
    double lower, upper;
} arc_ps_interval;

arc_presolve_status arc_presolve_run(arc_presolve *ps, const arc_csc *a,
                                     const double *b, const double *c,
                                     const double *lower, const double *upper,
                                     int apply_rules);

/* x (original length) from the reduced solution xr. */
void arc_postsolve(const arc_presolve *ps, const double *xr, double *x);

void arc_presolve_free(arc_presolve *ps);

static inline int arc_ps_proved(arc_presolve_status status)
{
    return -(int)status;
}

/* Individual rules, applied to one live row or column.  Each returns 1 when
 * it changed the problem, 0 when not, or arc_ps_proved(status) when it
 * settled the outcome. */
int arc_rule_empty_row(arc_presolve *ps, int i);
int arc_rule_empty_col(arc_presolve *ps, int j);
int arc_rule_singleton_row(arc_presolve *ps, int i);
int arc_rule_forcing_row(arc_presolve *ps, int i);
int arc_rule_implied_free(arc_presolve *ps, int i);
int arc_rule_slack_bound(arc_presolve *ps, int i);

/* Removes column j at the given value. */
void arc_ps_fix(arc_presolve *ps, int j, double value);

/* Gives x_j the bounds [lower, upper] (in its current frame): restores the
 * [0, upper] or free form, or fixes the variable.  Returns 0 or a negative
 * arc_presolve_status. */
int arc_ps_set_bounds(arc_presolve *ps, int j, double lower, double upper);

/* Bounds of all columns at once; NULL means 0 and +infinity. */
int arc_ps_import_bounds(arc_presolve *ps, const double *lower,
                         const double *upper);

arc_ps_interval arc_ps_bounds_of(const arc_presolve *ps, int j);

/* Builds the reduced CSC problem from the surviving rows and columns. */
int arc_extract_reduced(arc_presolve *ps);

#endif
