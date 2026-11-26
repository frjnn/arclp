/*
 * Rule 9 of arcLP.m.  When row i has a single coefficient a_ip whose sign
 * agrees with b_i (or, for b_i = 0, a sign shared by no other entry) then
 *     x_p = b_i / a_ip - sum_k (a_ik / a_ip) x_k
 * is nonnegative whenever the other x_k are, so x_p is substituted out of
 * every other row and the objective, and recovered in postsolve.
 */
#include "presolve.h"
#include "substitutions.h"

#include <math.h>

static const double RATIO_LIMIT = 1e6;
/* Share of its row's largest entry a free pivot needs: postsolve multiplies
 * the errors of the other variables by up to its inverse. */
static const double MIN_PIVOT = 0.1;
/* An entry this much smaller than what it was computed from is rounding. */
static const double CANCELLED = 1e-12;
static const long MAX_FILL = 1000;

/* A free variable can always be substituted out; a nonnegative one when
 * the row itself keeps it nonnegative.  Only columns with few entries are
 * used, so that the rule does not fill the matrix. */
static int cheap(const ps_matrix *a, int i, int j)
{
    return (long)(a->row_len[i] - 1) * (a->col_len[j] - 1) <= MAX_FILL;
}

static double largest_entry(const ps_matrix *a, int i)
{
    double largest = 0.0;
    PS_FOR_ROW(a, k, i)
    {
        largest = fmax(largest, fabs(a->e[k].value));
    }
    return largest;
}

/* The shortest free column with a large enough entry in row i, or -1. */
static int free_pivot(const arc_presolve *ps, int i, int *has_free)
{
    const ps_matrix *a = &ps->a;
    const double smallest = MIN_PIVOT * largest_entry(a, i);
    int best = -1;
    *has_free = 0;
    PS_FOR_ROW(a, k, i)
    {
        const int j = a->e[k].col;
        if (!ps->is_free[j])
        {
            continue;
        }
        *has_free = 1;
        if (fabs(a->e[k].value) >= smallest
            && (best < 0 || a->col_len[j] < a->col_len[a->e[best].col]))
        {
            best = k;
        }
    }
    return best;
}

static int choose_pivot(const arc_presolve *ps, int i)
{
    const ps_matrix *a = &ps->a;
    const double bi = ps->b[i];
    int has_free;
    const int free_entry = free_pivot(ps, i, &has_free);
    if (has_free)
    {
        return free_entry;
    }
    int positives = 0, negatives = 0, last_positive = -1, last_negative = -1;
    PS_FOR_ROW(a, k, i)
    {
        if (a->e[k].value > 0.0)
        {
            positives++;
            last_positive = k;
        }
        else
        {
            negatives++;
            last_negative = k;
        }
    }
    if (bi > 0.0 || (bi == 0.0 && positives <= negatives))
    {
        return positives == 1 ? last_positive : -1;
    }
    return negatives == 1 ? last_negative : -1;
}

/* old - delta, exactly zero when only rounding would be left: such an entry
 * must not become a coefficient of its own. */
static double subtract(double old, double delta)
{
    const double left = old - delta;
    return fabs(left) <= CANCELLED * fabs(delta) ? 0.0 : left;
}

/* row_j := row_j - (a_jp / a_ip) row_i over the gathered entries of row i. */
static int eliminate(arc_presolve *ps, int j, double a_jp, int count)
{
    ps_matrix *a = &ps->a;
    for (int k = a->row_head[j]; k != -1; k = a->e[k].next_in_row)
    {
        if (a->col_alive[a->e[k].col])
        {
            ps->scatter[a->e[k].col] = k;
        }
    }
    int status = 0;
    for (int t = 0; t < count && status == 0; t++)
    {
        const int col = ps->row_cols[t];
        const double delta = a_jp * ps->row_vals[t];
        const int k = ps->scatter[col];
        if (k >= 0)
        {
            ps_matrix_set(a, k, subtract(a->e[k].value, delta));
        }
        else if (delta != 0.0)
        {
            status = ps_matrix_insert(a, j, col, -delta);
        }
    }
    for (int k = a->row_head[j]; k != -1; k = a->e[k].next_in_row)
    {
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
    PS_FOR_ROW(a, k, i)
    {
        if (k != pivot_entry)
        {
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
                              count) != 0)
    {
        return -1;
    }
    PS_FOR_COL(a, k, p)
    {
        const int j = a->e[k].row;
        if (j == i)
        {
            continue;
        }
        const double a_jp = a->e[k].value;
        ps->b[j] -= a_jp * rhs;
        if (eliminate(ps, j, a_jp, count) != 0)
        {
            return -1;
        }
    }
    for (int t = 0; t < count; t++)
    {
        ps->c[ps->row_cols[t]] -= ps->c[p] * ps->row_vals[t];
    }
    ps_matrix_drop_row(a, i);
    ps_matrix_drop_col(a, p);
    return 0;
}

static int usable(const arc_presolve *ps, int i, int k)
{
    const int j = ps->a.e[k].col;
    return cheap(&ps->a, i, j) && (ps->is_free[j] || isinf(ps->upper[j]));
}

int arc_rule_implied_free(arc_presolve *ps, int i)
{
    if (!(ps->abs_ratio < RATIO_LIMIT))
    {
        return 0;
    }
    const int pivot_entry = choose_pivot(ps, i);
    if (pivot_entry < 0 || !usable(ps, i, pivot_entry))
    {
        return 0;
    }
    return substitute(ps, i, pivot_entry) == 0
               ? 1
               : arc_ps_proved(ARC_PRESOLVE_OUT_OF_MEMORY);
}
