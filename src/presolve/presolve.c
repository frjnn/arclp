#include "presolve.h"

#include "substitutions.h"

#include <stdlib.h>
#include <string.h>

typedef int (*presolve_rule)(arc_presolve *ps);

static int open_workspace(arc_presolve *ps, const arc_csc *a, const double *b,
                          const double *c)
{
    const size_t m = (size_t)a->m;
    const size_t n = (size_t)a->n;
    memset(ps, 0, sizeof(*ps));
    if (ps_matrix_init(&ps->a, a) != 0) {
        return -1;
    }
    ps->b = malloc(sizeof(double) * (m + 3 * n + 1));
    ps->scatter = malloc(sizeof(int) * (2 * n + 1));
    if (ps->b == NULL || ps->scatter == NULL) {
        return -1;
    }
    ps->c = ps->b + m;
    ps->fixed_x = ps->c + n;
    ps->row_vals = ps->fixed_x + n;
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

/* Applies the rules in arcLP.m's order until a full sweep changes nothing. */
static arc_presolve_status reduce(arc_presolve *ps)
{
    static const presolve_rule rules[] = {
        arc_rule_empty_rows, arc_rule_empty_cols, arc_rule_singleton_rows,
        arc_rule_forcing_rows, arc_rule_implied_free};
    int changed = 1;
    while (changed) {
        changed = 0;
        for (size_t r = 0; r < sizeof(rules) / sizeof(*rules); r++) {
            const int outcome = rules[r](ps);
            if (outcome < 0) {
                return (arc_presolve_status)(-outcome);
            }
            changed |= outcome;
        }
    }
    return ARC_PRESOLVE_REDUCED;
}

arc_presolve_status arc_presolve_run(arc_presolve *ps, const arc_csc *a,
                                     const double *b, const double *c,
                                     int apply_rules)
{
    if (open_workspace(ps, a, b, c) != 0) {
        return ARC_PRESOLVE_OUT_OF_MEMORY;
    }
    const arc_presolve_status status = apply_rules ? reduce(ps)
                                                   : ARC_PRESOLVE_REDUCED;
    if (status != ARC_PRESOLVE_REDUCED) {
        return status;
    }
    return arc_extract_reduced(ps) == 0 ? ARC_PRESOLVE_REDUCED
                                        : ARC_PRESOLVE_OUT_OF_MEMORY;
}

void arc_postsolve(const arc_presolve *ps, const double *xr, double *x)
{
    memcpy(x, ps->fixed_x, sizeof(double) * (size_t)ps->a.n);
    for (int k = 0; k < ps->reduced.n; k++) {
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
    free(ps->reduced_x);
    memset(ps, 0, sizeof(*ps));
}
