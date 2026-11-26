/*
 * CBLAS implementation of the dense kernels.  Products and solves with a
 * few right-hand sides against a large block run one level-2 call per
 * column, which streams the block better than level 3; everything else
 * is level 3.
 */
#include "blas_threading.h"
#include "dense.h"

#if defined(ARCLP_USE_ACCELERATE)
#define ACCELERATE_NEW_LAPACK
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#endif

static const int FEW_COLUMNS = 4;
static const double COLUMN_WORK = 1e5;

static int by_columns(int columns, double block)
{
    return columns == 1 || (columns <= FEW_COLUMNS && block >= COLUMN_WORK);
}

static enum CBLAS_TRANSPOSE transpose(arc_op op)
{
    return op == ARC_PLAIN ? CblasNoTrans : CblasTrans;
}

void arc_gemm(int m, int n, int k, double alpha, arc_factor a, arc_factor b,
              double beta, double *c, int ldc)
{
    const enum CBLAS_TRANSPOSE op = transpose(a.op);
    if (b.op == ARC_PLAIN && by_columns(n, (double)m * k))
    {
        const int rows = a.op == ARC_PLAIN ? m : k;
        const int cols = a.op == ARC_PLAIN ? k : m;
        arc_blas_threads_for((double)m * k);
        for (int col = 0; col < n; col++)
        {
            cblas_dgemv(CblasColMajor, op, rows, cols, alpha, a.x, a.ld,
                        b.x + (long)col * b.ld, 1, beta,
                        c + (long)col * ldc, 1);
        }
        return;
    }
    arc_blas_threads_for((double)m * n * k);
    cblas_dgemm(CblasColMajor, op, transpose(b.op), m, n, k, alpha, a.x, a.ld,
                b.x, b.ld, beta, c, ldc);
}

void arc_syrk_ln(int n, int k, double alpha, const double *a, int lda,
                 double beta, double *c, int ldc)
{
    arc_blas_threads_for((double)n * n * k);
    cblas_dsyrk(CblasColMajor, CblasLower, CblasNoTrans, n, k, alpha, a, lda,
                beta, c, ldc);
}

void arc_trsm_right_lt(int m, int n, const double *l, int ldl, double *b,
                       int ldb)
{
    arc_blas_threads_for((double)m * n * n);
    cblas_dtrsm(CblasColMajor, CblasRight, CblasLower, CblasTrans,
                CblasNonUnit, m, n, 1.0, l, ldl, b, ldb);
}

void arc_trsm_left(arc_op op, int n, int nrhs, const double *l, int ldl,
                   double *x, int ldx)
{
    const enum CBLAS_TRANSPOSE trans = transpose(op);
    if (by_columns(nrhs, (double)n * n))
    {
        arc_blas_threads_for((double)n * n);
        for (int r = 0; r < nrhs; r++)
        {
            cblas_dtrsv(CblasColMajor, CblasLower, trans, CblasNonUnit, n, l,
                        ldl, x + (long)r * ldx, 1);
        }
        return;
    }
    arc_blas_threads_for((double)n * n * nrhs);
    cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, trans, CblasNonUnit, n,
                nrhs, 1.0, l, ldl, x, ldx);
}
