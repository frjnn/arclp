/* Rules 1, 3, 5 and 7 of arcLP.m; rows and columns are visited from the
 * last to the first, as in the MATLAB loops. */
#include "presolve.h"

#include <math.h>

static const double EPSILON = 1e-8;
static const double COST_TOLERANCE = 100.0 * 1e-8;

/* Rule 1: an empty row must have b_i = 0. */
int arc_rule_empty_row(arc_presolve *ps, int i)
{
    if (ps->a.row_len[i] != 0)
    {
        return 0;
    }
    if (fabs(ps->b[i]) > EPSILON)
    {
        return arc_ps_proved(ARC_PRESOLVE_INFEASIBLE);
    }
    ps_matrix_drop_row(&ps->a, i);
    return 1;
}

/* Rule 3: an empty column sits at the bound its cost prefers. */
int arc_rule_empty_col(arc_presolve *ps, int j)
{
    if (ps->a.col_len[j] != 0)
    {
        return 0;
    }
    const int falls = ps->c[j] < -COST_TOLERANCE;
    const int rises = ps->c[j] > COST_TOLERANCE;
    if ((falls && isinf(ps->upper[j])) || (rises && ps->is_free[j]))
    {
        return arc_ps_proved(ARC_PRESOLVE_UNBOUNDED);
    }
    arc_ps_fix(ps, j, falls ? ps->upper[j] : 0.0);
    return 1;
}

/* Rule 5: a singleton row a_ij x_j = b_i fixes x_j. */
int arc_rule_singleton_row(arc_presolve *ps, int i)
{
    ps_matrix *a = &ps->a;
    if (a->row_len[i] != 1)
    {
        return 0;
    }
    const ps_entry *pivot = &a->e[ps_live_in_row(a, a->row_head[i])];
    const int j = pivot->col;
    const double xj = ps->b[i] / pivot->value;
    if (!ps->is_free[j]
        && (xj < -EPSILON || xj > ps->upper[j] + EPSILON * (1.0 + xj)))
    {
        return arc_ps_proved(ARC_PRESOLVE_INFEASIBLE);
    }
    arc_ps_fix(ps, j, ps->is_free[j] ? xj
                                     : fmin(fmax(xj, 0.0), ps->upper[j]));
    ps_matrix_drop_row(a, i);
    return 1;
}

/* Least and greatest value of row i over the bounds; -inf / +inf when a
 * term is unbounded. */
static arc_ps_interval activity(const arc_presolve *ps, int i)
{
    const ps_matrix *a = &ps->a;
    arc_ps_interval range = {0.0, 0.0};
    PS_FOR_ROW(a, k, i)
    {
        const double v = a->e[k].value;
        const arc_ps_interval x = arc_ps_bounds_of(ps, a->e[k].col);
        const double low = v > 0.0 ? x.lower : x.upper;
        const double high = v > 0.0 ? x.upper : x.lower;
        range.lower += low != 0.0 ? v * low : 0.0;
        range.upper += high != 0.0 ? v * high : 0.0;
    }
    return range;
}

/* Row i can only be met with every term at its extreme: fixes them all. */
static void force(arc_presolve *ps, int i, int at_lowest)
{
    ps_matrix *a = &ps->a;
    PS_FOR_ROW(a, k, i)
    {
        const int j = a->e[k].col;
        const arc_ps_interval x = arc_ps_bounds_of(ps, j);
        const int low = (a->e[k].value > 0.0) == at_lowest;
        arc_ps_fix(ps, j, low ? x.lower : x.upper);
    }
    ps_matrix_drop_row(a, i);
}

/* Rule 7, extended to bounds: rows whose right-hand side is out of reach
 * or exactly at the end of the reachable interval. */
int arc_rule_forcing_row(arc_presolve *ps, int i)
{
    if (ps->a.row_len[i] == 0)
    {
        return 0;
    }
    const arc_ps_interval range = activity(ps, i);
    const double margin = EPSILON * (1.0 + fabs(ps->b[i]));
    if (ps->b[i] < range.lower - margin || ps->b[i] > range.upper + margin)
    {
        return arc_ps_proved(ARC_PRESOLVE_INFEASIBLE);
    }
    const int at_lowest = ps->b[i] == range.lower;
    if (!at_lowest && ps->b[i] != range.upper)
    {
        return 0;
    }
    force(ps, i, at_lowest);
    return 1;
}
