/* Solves running in parallel threads agree with sequential ones (to
 * rounding: concurrent solves may run BLAS with another thread count), and
 * the library leaves the caller's BLAS threading mode as it found it.  The
 * absence of data races is checked by running this under ThreadSanitizer. */
#include "check.h"

#include <arclp/arclp.h>
#include <arclp/mat5_problem.h>

#include <pthread.h>

#if defined(__APPLE__) && defined(ARCLP_TEST_ACCELERATE)
#define ACCELERATE_NEW_LAPACK
#include <Accelerate/Accelerate.h>
#endif

enum { THREADS = 4, ROUNDS = 3 };

typedef struct
{
    const arclp_problem *problem;
    arclp_method method;
    double objective;
    int matches;
} job;

static void *run(void *arg)
{
    job *j = arg;
    arclp_options o;
    arclp_default_options(&o);
    o.method = j->method;
    double *x = malloc(sizeof(double) * (size_t)j->problem->n);
    j->matches = 1;
    for (int r = 0; r < ROUNDS; r++)
    {
        arclp_info info;
        const arclp_status s = arclp_solve(j->problem, &o, x, &info);
        j->matches &= s == ARCLP_OPTIMAL
                   && fabs(info.objective - j->objective)
                          <= 1e-9 * fabs(j->objective);
    }
    free(x);
    return NULL;
}

static void check_problem(const char *dir, const char *name)
{
    mat_lp m;
    if (mat_lp_load_named(dir, name, &m) != 0)
    {
        return;
    }
    const arclp_problem lp = mat_lp_as_problem(&m);
    job jobs[THREADS];
    pthread_t threads[THREADS];
    for (int t = 0; t < THREADS; t++)
    {
        arclp_options o;
        arclp_default_options(&o);
        o.method = t % 2 ? ARCLP_METHOD_MEHROTRA : ARCLP_METHOD_ARC_SEARCH;
        double *x = malloc(sizeof(double) * (size_t)m.n);
        arclp_info info;
        CHECK(arclp_solve(&lp, &o, x, &info) == ARCLP_OPTIMAL);
        free(x);
        const job j = {&lp, o.method, info.objective, 0};
        jobs[t] = j;
    }
    for (int t = 0; t < THREADS; t++)
    {
        pthread_create(&threads[t], NULL, run, &jobs[t]);
    }
    for (int t = 0; t < THREADS; t++)
    {
        pthread_join(threads[t], NULL);
        CHECK(jobs[t].matches);
    }
    mat_lp_free(&m);
}

int main(void)
{
    const char *dir = getenv("NETLIB_DIR");
    dir = dir != NULL ? dir : "data/arcLP";
#if defined(__APPLE__) && defined(ARCLP_TEST_ACCELERATE)
    BLASSetThreading(BLAS_THREADING_SINGLE_THREADED);
#endif
    check_problem(dir, "afiro");
    check_problem(dir, "scsd8");
    check_problem(dir, "ship04s");
    check_problem(dir, "maros_r7");
#if defined(__APPLE__) && defined(ARCLP_TEST_ACCELERATE)
    CHECK(BLASGetThreading() == BLAS_THREADING_SINGLE_THREADED);
#endif
    return check_report("concurrent solves");
}
