/* Normal-equations solver: random sparse A, residual of (A D A') y = r. */
#include "check.h"

#include "chol/chol_solver.h"
#include "chol/graph.h"
#include "util/vec.h"

typedef struct
{
    int *p, *i;
    double *x;
    arc_csc view;
} test_matrix;

static void test_matrix_free(test_matrix *t)
{
    free(t->p);
    free(t->i);
    free(t->x);
}

/* m x n matrix: an identity block (full row rank) plus `extra` random
 * entries per column; when duplicate_row >= 0 that row copies row 0. */
static test_matrix random_matrix(int m, int n, int extra, int duplicate_row,
                                 unsigned seed)
{
    test_matrix a = {malloc(sizeof(int) * (size_t)(n + 1)),
                     malloc(sizeof(int) * (size_t)n * (size_t)(extra + 2)),
                     malloc(sizeof(double) * (size_t)n * (size_t)(extra + 2)),
                     {0, 0, NULL, NULL, NULL}};
    int nz = 0;
    for (int j = 0; j < n; j++)
    {
        a.p[j] = nz;
        int used_zero = 0;
        for (int t = 0; t <= extra; t++)
        {
            int row = t == 0 ? j % m : (int)(check_random(&seed) * m);
            if (row == duplicate_row)
            {
                row = 0;
            }
            int seen = 0;
            for (int q = a.p[j]; q < nz; q++)
            {
                seen |= a.i[q] == row;
            }
            if (seen)
            {
                continue;
            }
            const double v = t == 0 ? 2.0 + check_random(&seed)
                                    : check_random(&seed) - 0.5;
            a.i[nz] = row;
            a.x[nz++] = v;
            used_zero |= row == 0;
        }
        if (duplicate_row >= 0 && used_zero)
        {
            for (int q = a.p[j]; q < nz; q++)
            {
                if (a.i[q] == 0)
                {
                    a.i[nz] = duplicate_row;
                    a.x[nz++] = a.x[q];
                }
            }
        }
    }
    a.p[n] = nz;
    const arc_csc view = {m, n, a.p, a.i, a.x};
    a.view = view;
    return a;
}

static double relative_residual(const arc_csc *a, const double *d,
                                const double *y, const double *r)
{
    double *t = malloc(sizeof(double) * (size_t)a->n);
    double *q = malloc(sizeof(double) * (size_t)a->m);
    arc_csc_mul_t(a, y, t);
    for (int j = 0; j < a->n; j++)
    {
        t[j] *= d[j];
    }
    arc_csc_mul(a, t, q);
    arc_subtract(a->m, r, q);
    const double res = arc_norm2(a->m, q) / arc_norm2(a->m, r);
    free(t);
    free(q);
    return res;
}

static double solve_case(int m, int n, int extra, int duplicate_row, int nrhs,
                         int *tiny)
{
    unsigned seed = (unsigned)(m * 31 + n);
    test_matrix t = random_matrix(m, n, extra, duplicate_row, seed);
    const arc_csc a = t.view;
    double *d = malloc(sizeof(double) * (size_t)n);
    double *r = malloc(sizeof(double) * (size_t)m * (size_t)nrhs);
    double *y = malloc(sizeof(double) * (size_t)m * (size_t)nrhs);
    for (int j = 0; j < n; j++)
    {
        d[j] = pow(10.0, 4.0 * check_random(&seed) - 2.0);
    }
    for (int k = 0; k < m * nrhs; k++)
    {
        y[k] = r[k] = check_random(&seed) - 0.5;
    }
    if (duplicate_row >= 0) /* keep r in the range of A D A' */
    {
        for (int c = 0; c < nrhs; c++)
        {
            y[duplicate_row + c * m] = r[duplicate_row + c * m] = r[c * m];
        }
    }
    arc_normal_solver *s = arc_chol_solver_create(&a);
    double worst = INFINITY;
    if (s != NULL)
    {
        arc_normal_factor(s, d);
        arc_normal_solve(s, y, nrhs);
        arc_normal_stats stats;
        s->ops->stats(s, &stats);
        *tiny = stats.tiny_pivots;
        worst = 0.0;
        for (int c = 0; c < nrhs; c++)
        {
            worst = fmax(worst, relative_residual(&a, d, y + c * m, r + c * m));
        }
        arc_normal_destroy(s);
    }
    test_matrix_free(&t);
    free(d);
    free(r);
    free(y);
    return worst;
}

static void check_amd_is_permutation(void)
{
    unsigned seed = 7;
    test_matrix t = random_matrix(300, 900, 3, -1, seed);
    arc_graph g;
    CHECK(arc_graph_of_aat(&t.view, &g) == 0);
    int *perm = malloc(sizeof(int) * 300);
    int *seen = calloc(300, sizeof(int));
    CHECK(arc_amd_order(&g, perm) == 0);
    for (int k = 0; k < 300; k++)
    {
        CHECK(perm[k] >= 0 && perm[k] < 300 && !seen[perm[k]]);
        seen[perm[k]] = 1;
    }
    free(perm);
    free(seen);
    arc_graph_free(&g);
    test_matrix_free(&t);
}

int main(void)
{
    int tiny = 0;
    check_amd_is_permutation();
    CHECK(solve_case(1, 3, 0, -1, 1, &tiny) < 1e-12);
    CHECK(solve_case(50, 120, 2, -1, 1, &tiny) < 1e-10);
    CHECK(solve_case(400, 1000, 3, -1, 2, &tiny) < 1e-9);
    CHECK(solve_case(2000, 6000, 4, -1, 3, &tiny) < 1e-8);
    CHECK(solve_case(1500, 1500, 12, -1, 1, &tiny) < 1e-8);
    CHECK(tiny == 0);
    CHECK(solve_case(300, 800, 3, 17, 1, &tiny) < 1e-6);
    return check_report("cholesky");
}
