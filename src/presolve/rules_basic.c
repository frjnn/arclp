/* Rules 1, 3, 5 and 7 of arcLP.m; rows and columns are visited from the
 * last to the first, as in the MATLAB loops. */
#include "presolve.h"

#include <math.h>

static const double EPSILON = 1e-8;
static const double COST_TOLERANCE = 100.0 * 1e-8;

static int proved(arc_presolve_status status)
{
    return -(int)status;
}

/* Rule 1: an empty row must have b_i = 0. */
int arc_rule_empty_rows(arc_presolve *ps)
{
    int changed = 0;
    for (int i = ps->a.m - 1; i >= 0; i--) {
        if (!ps->a.row_alive[i] || ps->a.row_len[i] != 0) {
            continue;
        }
        if (fabs(ps->b[i]) > EPSILON) {
            return proved(ARC_PRESOLVE_INFEASIBLE);
        }
        ps_matrix_drop_row(&ps->a, i);
        changed = 1;
    }
    return changed;
}

/* Rule 3: an empty column is fixed at zero unless its cost is negative. */
int arc_rule_empty_cols(arc_presolve *ps)
{
    int changed = 0;
    for (int j = ps->a.n - 1; j >= 0; j--) {
        if (!ps->a.col_alive[j] || ps->a.col_len[j] != 0) {
            continue;
        }
        if (ps->c[j] < -COST_TOLERANCE) {
            return proved(ARC_PRESOLVE_UNBOUNDED);
        }
        ps_matrix_drop_col(&ps->a, j);
        changed = 1;
    }
    return changed;
}

static int first_live(const ps_matrix *a, int i)
{
    int k = a->row_head[i];
    while (k != -1 && !ps_live(a, k)) {
        k = a->e[k].next_in_row;
    }
    return k;
}

/* Rule 5: a singleton row a_ij x_j = b_i fixes x_j. */
int arc_rule_singleton_rows(arc_presolve *ps)
{
    ps_matrix *a = &ps->a;
    int changed = 0;
    for (int i = a->m - 1; i >= 0; i--) {
        if (!a->row_alive[i] || a->row_len[i] != 1) {
            continue;
        }
        const ps_entry *pivot = &a->e[first_live(a, i)];
        const int j = pivot->col;
        const double xj = ps->b[i] / pivot->value;
        if (xj < 0.0) {
            return proved(ARC_PRESOLVE_INFEASIBLE);
        }
        ps->fixed_x[j] = xj;
        for (int k = a->col_head[j]; k != -1; k = a->e[k].next_in_col) {
            if (ps_live(a, k)) {
                ps->b[a->e[k].row] -= a->e[k].value * xj;
            }
        }
        ps_matrix_drop_col(a, j);
        ps_matrix_drop_row(a, i);
        changed = 1;
    }
    return changed;
}

/* Drops row i and every column whose coefficient has the given sign;
 * those variables are forced to zero. */
static void force_to_zero(ps_matrix *a, int i, int sign)
{
    for (int k = a->row_head[i]; k != -1; k = a->e[k].next_in_row) {
        if (ps_live(a, k) && (a->e[k].value > 0.0) == (sign > 0)) {
            ps_matrix_drop_col(a, a->e[k].col);
        }
    }
    ps_matrix_drop_row(a, i);
}

/* Rule 7: rows whose coefficients share one sign. */
int arc_rule_forcing_rows(arc_presolve *ps)
{
    ps_matrix *a = &ps->a;
    int changed = 0;
    for (int i = a->m - 1; i >= 0; i--) {
        if (!a->row_alive[i]) {
            continue;
        }
        int positive = 0;
        int negative = 0;
        for (int k = a->row_head[i]; k != -1; k = a->e[k].next_in_row) {
            if (ps_live(a, k)) {
                positive |= a->e[k].value > 0.0;
                negative |= a->e[k].value < 0.0;
            }
        }
        const double bi = ps->b[i];
        if ((bi < 0.0 && !negative) || (bi > 0.0 && !positive)) {
            return proved(ARC_PRESOLVE_INFEASIBLE);
        }
        if (bi == 0.0 && !(positive && negative)) {
            force_to_zero(a, i, positive ? 1 : -1);
            changed = 1;
        }
    }
    return changed;
}
