/*
 * Variable bounds.  During presolve every column is either free or lives in
 * [0, upper]; a change of variables that restores this form is recorded in
 * the substitution log, which postsolve replays backwards.
 */
#include "presolve.h"
#include "substitutions.h"

#include <math.h>
#include <stddef.h>

static const double EMPTY_RANGE = 1e-9;

void arc_ps_fix(arc_presolve *ps, int j, double value)
{
    ps_matrix *a = &ps->a;
    ps->fixed_x[j] = value;
    PS_FOR_COL(a, k, j)
    {
        ps->b[a->e[k].row] -= a->e[k].value * value;
    }
    ps_matrix_drop_col(a, j);
}

/* x_j := offset + sign x_j */
static int replace_variable(arc_presolve *ps, int j, double offset,
                            double sign)
{
    const double coef = -sign;
    if (arc_substitutions_add(&ps->subs, j, offset, &j, &coef, 1) != 0)
    {
        return -1;
    }
    ps_matrix *a = &ps->a;
    for (int k = a->col_head[j]; k != -1; k = a->e[k].next_in_col)
    {
        ps->b[a->e[k].row] -= a->e[k].value * offset;
        a->e[k].value *= sign;
    }
    ps->c[j] *= sign;
    return 0;
}

int arc_ps_set_bounds(arc_presolve *ps, int j, double lower, double upper)
{
    if (lower > upper + EMPTY_RANGE * (1.0 + fabs(upper)))
    {
        return arc_ps_proved(ARC_PRESOLVE_INFEASIBLE);
    }
    if (lower >= upper)
    {
        arc_ps_fix(ps, j, upper);
        return 0;
    }
    int status = 0;
    ps->is_free[j] = isinf(lower) && isinf(upper);
    ps->upper[j] = INFINITY;
    if (isfinite(lower))
    {
        ps->upper[j] = upper - lower;
        status = lower != 0.0 ? replace_variable(ps, j, lower, 1.0) : 0;
    }
    else if (isfinite(upper))
    {
        status = replace_variable(ps, j, upper, -1.0);
    }
    return status != 0 ? arc_ps_proved(ARC_PRESOLVE_OUT_OF_MEMORY) : 0;
}

int arc_ps_import_bounds(arc_presolve *ps, const double *lower,
                         const double *upper)
{
    for (int j = 0; j < ps->a.n; j++)
    {
        const int status =
            arc_ps_set_bounds(ps, j, lower != NULL ? lower[j] : 0.0,
                              upper != NULL ? upper[j] : INFINITY);
        if (status != 0)
        {
            return status;
        }
    }
    return 0;
}

arc_ps_interval arc_ps_bounds_of(const arc_presolve *ps, int j)
{
    const arc_ps_interval range = {ps->is_free[j] ? -INFINITY : 0.0,
                                   ps->upper[j]};
    return range;
}
