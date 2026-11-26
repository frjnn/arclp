/* Dense kernels against naive loops, on both sides of the size thresholds
 * where the BLAS backend switches routines or threading. */
#include "check.h"

#include "dense/dense.h"

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

/* y = alpha A x + beta y and y = alpha A' x + beta y, x with `cols` columns. */
static void check_products(int m, int k, int cols, unsigned *seed)
{
    double *a = random_block(m, k, seed);
    double *x = random_block(k > m ? k : m, cols, seed);
    double *y = random_block(m > k ? m : k, cols, seed);
    double *expect = malloc(sizeof(double) * (size_t)(m > k ? m : k) * cols);
    for (int c = 0; c < cols; c++)
    {
        for (int i = 0; i < m; i++)
        {
            double sum = 0.0;
            for (int p = 0; p < k; p++)
            {
                sum += a[i + p * m] * x[p + c * k];
            }
            expect[i + c * m] = 2.0 * sum + 0.5 * y[i + c * m];
        }
    }
    arc_gemm_nn(m, cols, k, 2.0, a, m, x, k, 0.5, y, m);
    CHECK(max_diff(m * cols, y, expect) < 1e-12 * k);
    for (int c = 0; c < cols; c++)
    {
        for (int j = 0; j < k; j++)
        {
            double sum = 0.0;
            for (int p = 0; p < m; p++)
            {
                sum += a[p + j * m] * x[p + c * m];
            }
            expect[j + c * k] = -sum + y[j + c * k];
        }
    }
    arc_gemm_tn(k, cols, m, -1.0, a, m, x, m, 1.0, y, k);
    CHECK(max_diff(k * cols, y, expect) < 1e-12 * m);
    free(a); free(x); free(y); free(expect);
}

/* L^{-1} and L^{-T} solves recover a known x. */
static void check_triangular(int n, int cols, unsigned *seed)
{
    double *l = random_lower(n, seed);
    double *x = random_block(n, cols, seed);
    double *b = malloc(sizeof(double) * (size_t)n * cols);
    for (int trans = 0; trans <= 1; trans++)
    {
        for (int c = 0; c < cols; c++)
        {
            for (int i = 0; i < n; i++)
            {
                double sum = 0.0;
                for (int p = 0; p < n; p++)
                {
                    sum += (trans ? l[p + i * n] : l[i + p * n]) * x[p + c * n];
                }
                b[i + c * n] = sum;
            }
        }
        if (trans)
        {
            arc_trsm_left_lt(n, cols, l, n, b, n);
        }
        else
        {
            arc_trsm_left_ln(n, cols, l, n, b, n);
        }
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
            check_products(sizes[s][0], sizes[s][1], columns[c], &seed);
            check_triangular(sizes[s][0], columns[c], &seed);
        }
    }
    arc_dense_end();
    return check_report("dense kernels");
}
