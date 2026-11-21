/* Portable implementation of the dense kernels (no BLAS dependency). */
#include "dense.h"

void arc_dense_begin(void) {}
void arc_dense_end(void) {}

static void scale(int m, int n, double beta, double *c, int ldc)
{
    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < m; i++)
        {
            c[i + j * ldc] = beta == 0.0 ? 0.0 : beta * c[i + j * ldc];
        }
    }
}

void arc_gemm_nt(int m, int n, int k, double alpha, const double *a, int lda,
                 const double *b, int ldb, double beta, double *c, int ldc)
{
    scale(m, n, beta, c, ldc);
    for (int j = 0; j < n; j++)
    {
        for (int p = 0; p < k; p++)
        {
            const double t = alpha * b[j + p * ldb];
            for (int i = 0; i < m; i++)
            {
                c[i + j * ldc] += t * a[i + p * lda];
            }
        }
    }
}

void arc_gemm_nn(int m, int n, int k, double alpha, const double *a, int lda,
                 const double *b, int ldb, double beta, double *c, int ldc)
{
    scale(m, n, beta, c, ldc);
    for (int j = 0; j < n; j++)
    {
        for (int p = 0; p < k; p++)
        {
            const double t = alpha * b[p + j * ldb];
            for (int i = 0; i < m; i++)
            {
                c[i + j * ldc] += t * a[i + p * lda];
            }
        }
    }
}

void arc_gemm_tn(int m, int n, int k, double alpha, const double *a, int lda,
                 const double *b, int ldb, double beta, double *c, int ldc)
{
    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < m; i++)
        {
            double sum = 0.0;
            for (int p = 0; p < k; p++)
            {
                sum += a[p + i * lda] * b[p + j * ldb];
            }
            const double old = beta == 0.0 ? 0.0 : beta * c[i + j * ldc];
            c[i + j * ldc] = old + alpha * sum;
        }
    }
}

void arc_syrk_ln(int n, int k, double alpha, const double *a, int lda,
                 double beta, double *c, int ldc)
{
    for (int j = 0; j < n; j++)
    {
        for (int i = j; i < n; i++)
        {
            c[i + j * ldc] = beta == 0.0 ? 0.0 : beta * c[i + j * ldc];
        }
        for (int p = 0; p < k; p++)
        {
            const double t = alpha * a[j + p * lda];
            for (int i = j; i < n; i++)
            {
                c[i + j * ldc] += t * a[i + p * lda];
            }
        }
    }
}

void arc_trsm_right_lt(int m, int n, const double *l, int ldl, double *b,
                       int ldb)
{
    for (int j = 0; j < n; j++)
    {
        for (int p = 0; p < j; p++)
        {
            const double t = l[j + p * ldl];
            for (int i = 0; i < m; i++)
            {
                b[i + j * ldb] -= t * b[i + p * ldb];
            }
        }
        const double inv = 1.0 / l[j + j * ldl];
        for (int i = 0; i < m; i++)
        {
            b[i + j * ldb] *= inv;
        }
    }
}

void arc_trsm_left_ln(int n, int nrhs, const double *l, int ldl, double *x,
                      int ldx)
{
    for (int r = 0; r < nrhs; r++)
    {
        double *v = x + r * ldx;
        for (int j = 0; j < n; j++)
        {
            v[j] /= l[j + j * ldl];
            for (int i = j + 1; i < n; i++)
            {
                v[i] -= l[i + j * ldl] * v[j];
            }
        }
    }
}

void arc_trsm_left_lt(int n, int nrhs, const double *l, int ldl, double *x,
                      int ldx)
{
    for (int r = 0; r < nrhs; r++)
    {
        double *v = x + r * ldx;
        for (int j = n - 1; j >= 0; j--)
        {
            double sum = v[j];
            for (int i = j + 1; i < n; i++)
            {
                sum -= l[i + j * ldl] * v[i];
            }
            v[j] = sum / l[j + j * ldl];
        }
    }
}
