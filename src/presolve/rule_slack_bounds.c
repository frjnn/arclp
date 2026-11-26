/*
 * Inequalities on a single variable.  A row  a x + g s = b  whose s appears
 * nowhere else only bounds x: s is replaced by b/g - (a/g) x and the row
 * disappears.
 */
#include "presolve.h"
#include "substitutions.h"

#include <math.h>

static const double MAX_RATIO = 1e6;

typedef struct
{
    int x, s; /* live entries of the row */
} row_pair;

/* Two live entries, one of them a column singleton whose variable is not
 * free; 0 when the row does not have this shape. */
static int slack_pair(const arc_presolve *ps, int i, row_pair *pair)
{
    const ps_matrix *a = &ps->a;
    if (a->row_len[i] != 2)
    {
        return 0;
    }
    int found[2];
    int count = 0;
    PS_FOR_ROW(a, k, i)
    {
        found[count++] = k;
    }
    for (int t = 0; t < 2; t++)
    {
        const int s = a->e[found[t]].col;
        if (a->col_len[s] == 1 && !ps->is_free[s])
        {
            pair->s = found[t];
            pair->x = found[1 - t];
            const double ratio = a->e[pair->x].value / a->e[pair->s].value;
            return fabs(ratio) < MAX_RATIO && fabs(ratio) > 1.0 / MAX_RATIO;
        }
    }
    return 0;
}

/* The interval of x for which s = (b - a x) / g stays in [0, upper_s]. */
static arc_ps_interval implied_range(const arc_presolve *ps, int i,
                                     row_pair p)
{
    const double a = ps->a.e[p.x].value;
    const double g = ps->a.e[p.s].value;
    const double at_zero = ps->b[i] / a;                        /* s = 0 */
    const double upper = ps->upper[ps->a.e[p.s].col];
    const double at_upper = isinf(upper)
                                ? (g / a > 0.0 ? -INFINITY : INFINITY)
                                : (ps->b[i] - g * upper) / a;
    const arc_ps_interval range = {fmin(at_zero, at_upper),
                                   fmax(at_zero, at_upper)};
    return range;
}

static int replace_row(arc_presolve *ps, int i, row_pair p)
{
    ps_matrix *a = &ps->a;
    const int x = a->e[p.x].col;
    const int s = a->e[p.s].col;
    const double coef = a->e[p.x].value / a->e[p.s].value;
    const arc_ps_interval implied = implied_range(ps, i, p);
    const arc_ps_interval own = arc_ps_bounds_of(ps, x);
    if (arc_substitutions_add(&ps->subs, s, ps->b[i] / a->e[p.s].value, &x,
                              &coef, 1) != 0)
    {
        return arc_ps_proved(ARC_PRESOLVE_OUT_OF_MEMORY);
    }
    ps->c[x] -= ps->c[s] * coef;
    ps_matrix_drop_row(a, i);
    ps_matrix_drop_col(a, s);
    const int status = arc_ps_set_bounds(ps, x, fmax(own.lower, implied.lower),
                                         fmin(own.upper, implied.upper));
    return status != 0 ? status : 1;
}

int arc_rule_slack_bound(arc_presolve *ps, int i)
{
    row_pair pair;
    return slack_pair(ps, i, &pair) ? replace_row(ps, i, pair) : 0;
}
