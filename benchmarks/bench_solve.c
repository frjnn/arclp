#include "bench.h"

#include <arclp/mat5_problem.h>

#include <math.h>
#include <stdlib.h>

static const double SOLVED_ERROR = 1e-6;

int bench_solve(const mat_lp *lp, arclp_method method, int repeat,
                double reference, bench_result *out)
{
    const arclp_problem problem = mat_lp_as_problem(lp);
    arclp_options options;
    arclp_default_options(&options);
    options.method = method;
    double *x = malloc(sizeof(double) * (size_t)(lp->n + 1));
    if (x == NULL)
    {
        return -1;
    }
    out->seconds = INFINITY;
    for (int r = 0; r < repeat; r++)
    {
        arclp_solve(&problem, &options, x, &out->info);
        out->seconds = fmin(out->seconds, out->info.seconds_total);
    }
    out->objective_error = fabs(out->info.objective - reference)
                         / fmax(1.0, fabs(reference));
    free(x);
    return 0;
}

void bench_accumulate(bench_totals *t, const bench_result *r,
                      int paper_iterations)
{
    t->problems++;
    t->solved += r->info.status == ARCLP_OPTIMAL
              && r->objective_error < SOLVED_ERROR;
    t->iterations += r->info.iterations;
    t->paper_iterations += paper_iterations;
    t->seconds += r->seconds;
}
