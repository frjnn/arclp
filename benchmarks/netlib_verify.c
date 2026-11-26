/*
 * Checks every returned solution against the ORIGINAL problem (before
 * presolve), with presolve enabled and disabled:
 *     primal residual  ||A x - b|| / (1 + ||b||),  min x_j,
 *     objective error  |c'x - ref| / max(1, |ref|)  (HiGHS reference).
 *     netlib_verify DATA_DIR
 * Exits with failure when a run is not within the tolerances below.
 */
#include <arclp/arclp.h>
#include <arclp/mat5.h>

#include "paper_table.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static const double MAX_RESIDUAL = 1e-6;
static const double MAX_NEGATIVE = 1e-8;
static const double MAX_OBJECTIVE_ERROR = 1e-5;

typedef struct
{
    arclp_status status;
    double residual, lowest, objective_error;
} check;

static check verify(const mat_lp *m, const double *x, arclp_status status,
                    double reference)
{
    check c = {status, 0.0, INFINITY, 0.0};
    double *ax = calloc((size_t)m->m + 1, sizeof(double));
    double nb = 0.0, objective = 0.0;
    for (int j = 0; j < m->n; j++)
    {
        for (int k = m->col_ptr[j]; k < m->col_ptr[j + 1]; k++)
        {
            ax[m->row_idx[k]] += m->values[k] * x[j];
        }
        c.lowest = fmin(c.lowest, x[j]);
        objective += m->c[j] * x[j];
    }
    for (int i = 0; i < m->m; i++)
    {
        c.residual += (ax[i] - m->b[i]) * (ax[i] - m->b[i]);
        nb += m->b[i] * m->b[i];
    }
    c.residual = sqrt(c.residual) / (1.0 + sqrt(nb));
    c.objective_error = fabs(objective - reference) / fmax(1.0, fabs(reference));
    free(ax);
    return c;
}

static int acceptable(check c)
{
    return c.status == ARCLP_OPTIMAL && c.residual <= MAX_RESIDUAL
        && c.lowest >= -MAX_NEGATIVE && c.objective_error <= MAX_OBJECTIVE_ERROR;
}

static check run(const mat_lp *m, arclp_method method, int presolve,
                 double reference)
{
    const arclp_problem lp = {m->m, m->n, m->col_ptr, m->row_idx, m->values,
                              m->b, m->c};
    arclp_options o;
    arclp_default_options(&o);
    o.method = method;
    o.presolve = presolve;
    double *x = calloc((size_t)m->n + 1, sizeof(double));
    const arclp_status s = arclp_solve(&lp, &o, x, NULL);
    const check c = verify(m, x, s, reference);
    free(x);
    return c;
}

static int report(const char *name, const char *method, int presolve, check c)
{
    const int ok = acceptable(c);
    printf("%-9s %-8s presolve=%d  %-10s  residual %8.1e  min x %9.1e  "
           "obj err %8.1e  %s\n", name, method, presolve,
           arclp_status_name(c.status), c.residual, c.lowest,
           c.objective_error, ok ? "ok" : "FAIL");
    return ok;
}

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        fprintf(stderr, "usage: %s DATA_DIR\n", argv[0]);
        return EXIT_FAILURE;
    }
    int failures = 0;
    for (int k = 0; k < PAPER_TABLE_SIZE; k++)
    {
        const paper_row *row = &PAPER_TABLE[k];
        char path[4096], err[256];
        snprintf(path, sizeof(path), "%s/lp_%s.mat", argv[1], row->name);
        mat_lp m;
        if (mat_lp_load(path, &m, err, sizeof(err)) != 0)
        {
            fprintf(stderr, "skipping %s: %s\n", row->name, err);
            continue;
        }
        for (int presolve = 1; presolve >= 0; presolve--)
        {
            failures += !report(row->name, "arc", presolve,
                                run(&m, ARCLP_METHOD_ARC_SEARCH, presolve,
                                    row->reference_objective));
            failures += !report(row->name, "mehrotra", presolve,
                                run(&m, ARCLP_METHOD_MEHROTRA, presolve,
                                    row->reference_objective));
        }
        mat_lp_free(&m);
    }
    printf("%d run(s) outside tolerance\n", failures);
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
