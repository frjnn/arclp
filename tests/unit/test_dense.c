/* Dense kernels against naive loops, on both sides of the size thresholds
 * where the BLAS backend switches routines or threading. */
#include "check.h"

#include "dense/dense.h"

#include <string.h>

static double *random_block(int rows, int cols, unsigned *seed)
{
    double *a = malloc(sizeof(double) * (size_t)rows * (size_t)cols);
    for (int k = 0; k < rows * cols; k++)
    {
        a[k] = check_random(seed) - 0.5;
    }
    return a;
}

/* Well-conditioned lower triangle. */
static double *random_lower(int n, unsigned *seed)
{
    double *l = random_block(n, n, seed);
    for (int j = 0; j < n; j++)
    {
        l[j + j * n] = 2.0 + check_random(seed);
        for (int i = 0; i < j; i++)
        {
            l[i + j * n] = 0.0;
        }
        for (int i = j + 1; i < n; i++)
        {
            l[i + j * n] /= n;
        }
    }
    return l;
}

static double max_diff(int len, const double *x, const double *y)
{
    double worst = 0.0;
    for (int k = 0; k < len; k++)
    {
        worst = fmax(worst, fabs(x[k] - y[k]));
    }
    return worst;
}

/* y := alpha op(A) x + beta y for `cols` columns, one entry at a time; the
 * stored A has leading dimension lda and op(A) is rows x inner. */
static void naive_product(int trans, int rows, int inner, int cols,
                          double alpha, const double *a, int lda,
                          const double *x, double beta, double *y)
{
    for (int c = 0; c < cols; c++)
    {
        for (int i = 0; i < rows; i++)
        {
            double sum = 0.0;
            for (int p = 0; p < inner; p++)
            {
                sum += (trans ? a[p + i * lda] : a[i + p * lda])
                     * x[p + c * inner];
            }
            y[i + c * rows] = alpha * sum + beta * y[i + c * rows];
        }
    }
}

static void check_product(int trans, int rows, int inner, int cols,
                          double alpha, double beta, unsigned *seed)
{
    const int lda = trans ? inner : rows;
    double *a = random_block(rows, inner, seed);
    double *x = random_block(inner, cols, seed);
    double *y = random_block(rows, cols, seed);
    double *expect = malloc(sizeof(double) * (size_t)rows * cols);
    memcpy(expect, y, sizeof(double) * (size_t)rows * cols);
    naive_product(trans, rows, inner, cols, alpha, a, lda, x, beta, expect);
    const arc_factor left = {a, lda, trans ? ARC_TRANSPOSED : ARC_PLAIN};
    const arc_factor right = {x, inner, ARC_PLAIN};
    arc_gemm(rows, cols, inner, alpha, left, right, beta, y, rows);
    CHECK(max_diff(rows * cols, y, expect) < 1e-12 * inner);
    free(a); free(x); free(y); free(expect);
}

/* L^{-1} and L^{-T} solves recover a known x. */
static void check_triangular(int n, int cols, unsigned *seed)
{
    double *l = random_lower(n, seed);
    double *x = random_block(n, cols, seed);
    double *b = calloc((size_t)n * (size_t)cols, sizeof(double));
    for (int trans = 0; trans <= 1; trans++)
    {
        naive_product(trans, n, n, cols, 1.0, l, n, x, 0.0, b);
        arc_trsm_left(trans ? ARC_TRANSPOSED : ARC_PLAIN, n, cols, l, n, b,
                      n);
        CHECK(max_diff(n * cols, b, x) < 1e-10);
    }
    free(l); free(x); free(b);
}

int main(void)
{
    unsigned seed = 11;
    arc_dense_begin();
    const int sizes[][2] = {{3, 5}, {40, 30}, {400, 300}, {700, 500}};
    const int columns[] = {1, 2, 3, 5};
    for (size_t s = 0; s < sizeof(sizes) / sizeof(*sizes); s++)
    {
        for (size_t c = 0; c < sizeof(columns) / sizeof(*columns); c++)
        {
            const int m = sizes[s][0];
            const int k = sizes[s][1];
            check_product(0, m, k, columns[c], 2.0, 0.5, &seed);
            check_product(1, k, m, columns[c], -1.0, 1.0, &seed);
            check_triangular(sizes[s][0], columns[c], &seed);
        }
    }
    arc_dense_end();
    return check_report("dense kernels");
}
