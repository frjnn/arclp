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

/* Plain A: C gathers multiples of A's columns. */
static void gemm_by_columns(int m, int n, int k, double alpha, arc_factor a,
                            arc_factor b, double beta, double *c, int ldc)
{
    const int row_step = b.op == ARC_PLAIN ? 1 : b.ld;
    const int col_step = b.op == ARC_PLAIN ? b.ld : 1;
    scale(m, n, beta, c, ldc);
    for (int j = 0; j < n; j++)
    {
        for (int p = 0; p < k; p++)
        {
            const double t = alpha * b.x[p * row_step + j * col_step];
            for (int i = 0; i < m; i++)
            {
                c[i + j * ldc] += t * a.x[i + p * a.ld];
            }
        }
    }
}

/* Transposed A, plain B: every entry of C is a dot product. */
static void gemm_by_dots(int m, int n, int k, double alpha, arc_factor a,
                         arc_factor b, double beta, double *c, int ldc)
{
    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < m; i++)
        {
            double sum = 0.0;
            for (int p = 0; p < k; p++)
            {
                sum += a.x[p + i * a.ld] * b.x[p + j * b.ld];
            }
            const double old = beta == 0.0 ? 0.0 : beta * c[i + j * ldc];
            c[i + j * ldc] = old + alpha * sum;
        }
    }
}

void arc_gemm(int m, int n, int k, double alpha, arc_factor a, arc_factor b,
              double beta, double *c, int ldc)
{
    if (a.op == ARC_PLAIN)
    {
        gemm_by_columns(m, n, k, alpha, a, b, beta, c, ldc);
    }
    else
    {
        gemm_by_dots(m, n, k, alpha, a, b, beta, c, ldc);
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

static void solve_plain(int n, const double *l, int ldl, double *v)
{
    for (int j = 0; j < n; j++)
    {
        v[j] /= l[j + j * ldl];
        for (int i = j + 1; i < n; i++)
        {
            v[i] -= l[i + j * ldl] * v[j];
        }
    }
}

static void solve_transposed(int n, const double *l, int ldl, double *v)
{
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

void arc_trsm_left(arc_op op, int n, int nrhs, const double *l, int ldl,
                   double *x, int ldx)
{
    for (int r = 0; r < nrhs; r++)
    {
        if (op == ARC_PLAIN)
        {
            solve_plain(n, l, ldl, x + r * ldx);
        }
        else
        {
            solve_transposed(n, l, ldl, x + r * ldx);
        }
    }
}
