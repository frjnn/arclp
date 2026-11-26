/* Public facade: validation, presolve, interior-point solve, postsolve. */
#include <arclp/arclp.h>

#include "presolve/presolve.h"
#include "util/clock.h"
#include "util/vec.h"
#include "solve_reduced.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

void arclp_default_options(arclp_options *options)
{
    options->method = ARCLP_METHOD_ARC_SEARCH;
    options->tolerance = 1e-8;
    options->max_iterations = 10000;
    options->presolve = 1;
    options->verbose = 0;
}

const char *arclp_status_name(arclp_status status)
{
    static const char *const names[] = {"optimal", "infeasible", "unbounded",
                                        "iteration limit", "invalid input",
                                        "out of memory", "numerical error",
                                        "stalled"};
    const unsigned k = (unsigned)status;
    return k < sizeof(names) / sizeof(*names) ? names[k] : "unknown";
}

static int all_finite(const double *v, int len)
{
    for (int k = 0; k < len; k++)
    {
        if (!isfinite(v[k]))
        {
            return 0;
        }
    }
    return 1;
}

/* lower < +inf, upper > -inf and no NaN; crossed bounds are for presolve. */
static int valid_bounds(const arclp_problem *p)
{
    for (int j = 0; j < p->n; j++)
    {
        const double lower = p->lower != NULL ? p->lower[j] : 0.0;
        const double upper = p->upper != NULL ? p->upper[j] : INFINITY;
        if (!(lower < INFINITY && upper > -INFINITY))
        {
            return 0;
        }
    }
    return 1;
}

static int valid_matrix(const arclp_problem *p)
{
    if (p->col_ptr[0] != 0)
    {
        return 0;
    }
    for (int j = 0; j < p->n; j++)
    {
        if (p->col_ptr[j + 1] < p->col_ptr[j])
        {
            return 0;
        }
    }
    for (int k = 0; k < p->col_ptr[p->n]; k++)
    {
        if (p->row_idx[k] < 0 || p->row_idx[k] >= p->m)
        {
            return 0;
        }
    }
    return all_finite(p->values, p->col_ptr[p->n]) && all_finite(p->b, p->m)
        && all_finite(p->c, p->n) && valid_bounds(p);
}

static int valid(const arclp_problem *p, const arclp_options *o, double *x)
{
    return p != NULL && x != NULL && p->m >= 0 && p->n >= 0
        && p->col_ptr != NULL && p->b != NULL && p->c != NULL
        && (p->col_ptr[p->n] == 0 || (p->row_idx && p->values))
        && o->tolerance > 0.0 && o->max_iterations > 0
        && (o->method == ARCLP_METHOD_ARC_SEARCH
            || o->method == ARCLP_METHOD_MEHROTRA)
        && valid_matrix(p);
}

static arclp_status presolve_status(arc_presolve_status s)
{
    static const arclp_status map[] = {ARCLP_OPTIMAL, ARCLP_INFEASIBLE,
                                       ARCLP_UNBOUNDED, ARCLP_OUT_OF_MEMORY};
    return map[s];
}

static arclp_status solve(const arclp_problem *p, const arclp_options *o,
                          double *x, arclp_info *info)
{
    const double start = arc_seconds();
    const arc_csc a = {p->m, p->n, p->col_ptr, p->row_idx, p->values};
    arc_presolve ps;
    arclp_status status =
        presolve_status(arc_presolve_run(&ps, &a, p->b, p->c, p->lower,
                                         p->upper, o->presolve));
    info->seconds_presolve = arc_seconds() - start;
    if (status == ARCLP_OPTIMAL)
    {
        info->reduced_rows = ps.reduced.a.m;
        info->reduced_cols = ps.reduced.a.n;
        double *xr = malloc(sizeof(double) * (size_t)(ps.reduced.a.n + 1));
        status = xr == NULL ? ARCLP_OUT_OF_MEMORY
                            : arc_solve_reduced(&ps.reduced, o, xr, info);
        if (xr != NULL)
        {
            arc_postsolve(&ps, xr, x);
        }
        free(xr);
    }
    arc_presolve_free(&ps);
    const int usable = status == ARCLP_OPTIMAL || status == ARCLP_STALLED
                    || status == ARCLP_MAX_ITERATIONS;
    info->objective = usable ? arc_dot(p->n, p->c, x) : NAN;
    return status;
}

arclp_status arclp_solve(const arclp_problem *problem,
                         const arclp_options *options, double *x,
                         arclp_info *info)
{
    const double start = arc_seconds();
    arclp_options defaults;
    arclp_default_options(&defaults);
    const arclp_options *o = options != NULL ? options : &defaults;
    arclp_info local;
    arclp_info *out = info != NULL ? info : &local;
    memset(out, 0, sizeof(*out));
    out->objective = NAN;
    out->status = valid(problem, o, x) ? solve(problem, o, x, out)
                                       : ARCLP_INVALID_INPUT;
    out->seconds_total = arc_seconds() - start;
    return out->status;
}
