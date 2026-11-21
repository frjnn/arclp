/*
 * Rule 9 of arcLP.m.  When row i has a single coefficient a_ip whose sign
 * agrees with b_i (or, for b_i = 0, a sign shared by no other entry) then
 *     x_p = b_i / a_ip - sum_k (a_ik / a_ip) x_k
 * is nonnegative whenever the other x_k are, so x_p is substituted out of
 * every other row and the objective, and recovered in postsolve.
 */
#include "presolve.h"
#include "substitutions.h"

static const double RATIO_LIMIT = 1e6;

static int choose_pivot(const ps_matrix *a, int i, double bi)
{
    int positives = 0, negatives = 0, last_positive = -1, last_negative = -1;
    for (int k = a->row_head[i]; k != -1; k = a->e[k].next_in_row) {
        if (!ps_live(a, k)) {
            continue;
        }
        if (a->e[k].value > 0.0) {
            positives++;
            last_positive = k;
        } else {
            negatives++;
            last_negative = k;
        }
    }
    if (bi > 0.0 || (bi == 0.0 && positives <= negatives)) {
        return positives == 1 ? last_positive : -1;
    }
    return negatives == 1 ? last_negative : -1;
}

/* row_j := row_j - (a_jp / a_ip) row_i over the gathered entries of row i. */
static int eliminate(arc_presolve *ps, int j, double a_jp, int count)
{
    ps_matrix *a = &ps->a;
    for (int k = a->row_head[j]; k != -1; k = a->e[k].next_in_row) {
        if (a->col_alive[a->e[k].col]) {
            ps->scatter[a->e[k].col] = k;
        }
    }
    int status = 0;
    for (int t = 0; t < count && status == 0; t++) {
        const int col = ps->row_cols[t];
        const double delta = a_jp * ps->row_vals[t];
        const int k = ps->scatter[col];
        if (k >= 0) {
            ps_matrix_set(a, k, a->e[k].value - delta);
        } else if (delta != 0.0) {
            status = ps_matrix_insert(a, j, col, -delta);
        }
    }
    for (int k = a->row_head[j]; k != -1; k = a->e[k].next_in_row) {
        ps->scatter[a->e[k].col] = -1;
    }
    return status;
}

/* Gathers (a_ik / a_ip) for the other live entries of row i. */
static int gather_row(arc_presolve *ps, int i, int pivot_entry)
{
    const ps_matrix *a = &ps->a;
    const double a_ip = a->e[pivot_entry].value;
    int count = 0;
    for (int k = a->row_head[i]; k != -1; k = a->e[k].next_in_row) {
        if (k != pivot_entry && ps_live(a, k)) {
            ps->row_cols[count] = a->e[k].col;
            ps->row_vals[count++] = a->e[k].value / a_ip;
        }
    }
    return count;
}

static int substitute(arc_presolve *ps, int i, int pivot_entry)
{
    ps_matrix *a = &ps->a;
    const int p = a->e[pivot_entry].col;
    const double rhs = ps->b[i] / a->e[pivot_entry].value;
    const int count = gather_row(ps, i, pivot_entry);
    if (arc_substitutions_add(&ps->subs, p, rhs, ps->row_cols, ps->row_vals,
                              count) != 0) {
        return -1;
    }
    for (int k = a->col_head[p]; k != -1; k = a->e[k].next_in_col) {
        const int j = a->e[k].row;
        if (j == i || !ps_live(a, k)) {
            continue;
        }
        const double a_jp = a->e[k].value;
        ps->b[j] -= a_jp * rhs;
        if (eliminate(ps, j, a_jp, count) != 0) {
            return -1;
        }
    }
    for (int t = 0; t < count; t++) {
        ps->c[ps->row_cols[t]] -= ps->c[p] * ps->row_vals[t];
    }
    ps_matrix_drop_row(a, i);
    ps_matrix_drop_col(a, p);
    return 0;
}

int arc_rule_implied_free(arc_presolve *ps)
{
    if (!(ps->abs_ratio < RATIO_LIMIT)) {
        return 0;
    }
    int changed = 0;
    for (int i = ps->a.m - 1; i >= 0; i--) {
        if (!ps->a.row_alive[i]) {
            continue;
        }
        const int pivot_entry = choose_pivot(&ps->a, i, ps->b[i]);
        if (pivot_entry < 0) {
            continue;
        }
        if (substitute(ps, i, pivot_entry) != 0) {
            return -(int)ARC_PRESOLVE_OUT_OF_MEMORY;
        }
        changed = 1;
    }
    return changed;
}
