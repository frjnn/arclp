#include "presolve.h"

#include "substitutions.h"

#include <stdlib.h>
#include <string.h>

static int open_workspace(arc_presolve *ps, const arc_csc *a, const double *b,
                          const double *c)
{
    const size_t m = (size_t)a->m;
    const size_t n = (size_t)a->n;
    memset(ps, 0, sizeof(*ps));
    ps->n = a->n;
    if (ps_matrix_init(&ps->a, a) != 0)
    {
        return -1;
    }
    ps->b = malloc(sizeof(double) * (m + 4 * n + 1));
    ps->scatter = malloc(sizeof(int) * (2 * n + 1));
    ps->is_free = malloc(n + 1);
    if (ps->b == NULL || ps->scatter == NULL || ps->is_free == NULL)
    {
        return -1;
    }
    ps->c = ps->b + m;
    ps->fixed_x = ps->c + n;
    ps->row_vals = ps->fixed_x + n;
    ps->upper = ps->row_vals + n;
    ps->row_cols = ps->scatter + n;
    memcpy(ps->b, b, sizeof(double) * m);
    memcpy(ps->c, c, sizeof(double) * n);
    memset(ps->fixed_x, 0, sizeof(double) * n);
    memset(ps->scatter, 0xff, sizeof(int) * n);
    double lowest = 0.0, highest = 0.0;
    ps->abs_ratio = arc_csc_abs_range(a, &lowest, &highest) ? highest / lowest
                                                            : 0.0;
    return 0;
}

typedef struct
{
    int (*apply)(arc_presolve *ps, int index);
    int on_columns;
    int forward; /* first to last instead of last to first */
} presolve_rule;

/* One pass of a rule over the live rows or columns. */
static int sweep(arc_presolve *ps, const presolve_rule *rule)
{
    const unsigned char *alive = rule->on_columns ? ps->a.col_alive
                                                  : ps->a.row_alive;
    const int count = rule->on_columns ? ps->a.n : ps->a.m;
    int changed = 0;
    for (int t = 0; t < count; t++)
    {
        const int index = rule->forward ? t : count - 1 - t;
        const int outcome = alive[index] ? rule->apply(ps, index) : 0;
        if (outcome < 0)
        {
            return outcome;
        }
        changed |= outcome;
    }
    return changed;
}

/* Applies the rules in arcLP.m's order until a full sweep changes nothing. */
static arc_presolve_status reduce(arc_presolve *ps)
{
    static const presolve_rule rules[] = {
        {arc_rule_empty_row, 0, 0},    {arc_rule_empty_col, 1, 0},
        {arc_rule_singleton_row, 0, 0}, {arc_rule_slack_bound, 0, 1},
        {arc_rule_forcing_row, 0, 0},  {arc_rule_implied_free, 0, 0}};
    int changed = 1;
    while (changed)
    {
        changed = 0;
        for (size_t r = 0; r < sizeof(rules) / sizeof(*rules); r++)
        {
            const int outcome = sweep(ps, &rules[r]);
            if (outcome < 0)
            {
                return (arc_presolve_status)(-outcome);
            }
            changed |= outcome;
        }
    }
    return ARC_PRESOLVE_REDUCED;
}

arc_presolve_status arc_presolve_run(arc_presolve *ps, const arc_csc *a,
                                     const double *b, const double *c,
                                     const double *lower, const double *upper,
                                     int apply_rules)
{
    if (open_workspace(ps, a, b, c) != 0)
    {
        return ARC_PRESOLVE_OUT_OF_MEMORY;
    }
    const int bounds = arc_ps_import_bounds(ps, lower, upper);
    if (bounds != 0)
    {
        return (arc_presolve_status)(-bounds);
    }
    const arc_presolve_status status = apply_rules ? reduce(ps)
                                                   : ARC_PRESOLVE_REDUCED;
    if (status != ARC_PRESOLVE_REDUCED)
    {
        return status;
    }
    if (arc_extract_reduced(ps) != 0)
    {
        return ARC_PRESOLVE_OUT_OF_MEMORY;
    }
    /* Postsolve needs neither, and the solver runs while ps is alive. */
    ps_matrix_free(&ps->a);
    free(ps->scatter);
    ps->scatter = NULL;
    ps->row_cols = NULL;
    return ARC_PRESOLVE_REDUCED;
}

void arc_postsolve(const arc_presolve *ps, const double *xr, double *x)
{
    memcpy(x, ps->fixed_x, sizeof(double) * (size_t)ps->n);
    for (int k = 0; k < ps->reduced.a.n; k++)
    {
        x[ps->reduced_col[k]] = xr[k];
    }
    arc_substitutions_apply(&ps->subs, x);
}

void arc_presolve_free(arc_presolve *ps)
{
    ps_matrix_free(&ps->a);
    arc_substitutions_free(&ps->subs);
    free(ps->b);
    free(ps->scatter);
    free(ps->is_free);
    arc_arena_release(&ps->arena);
    memset(ps, 0, sizeof(*ps));
}
