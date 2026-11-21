/*
 * CBLAS implementation of the dense kernels.  Products and solves with a
 * few right-hand sides against a large block run one level-2 call per
 * column, which streams the block better than level 3 (EXPERIMENTS.md,
 * E23); everything else is level 3.
 */
#include "blas_threading.h"
#include "dense.h"

#ifdef __APPLE__
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

void arc_gemm_nt(int m, int n, int k, double alpha, const double *a, int lda,
                 const double *b, int ldb, double beta, double *c, int ldc)
{
    arc_blas_threads_for((double)m * n * k);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, m, n, k, alpha, a,
                lda, b, ldb, beta, c, ldc);
}

void arc_gemm_nn(int m, int n, int k, double alpha, const double *a, int lda,
                 const double *b, int ldb, double beta, double *c, int ldc)
{
    if (by_columns(n, (double)m * k)) {
        arc_blas_threads_for((double)m * k);
        for (int col = 0; col < n; col++) {
            cblas_dgemv(CblasColMajor, CblasNoTrans, m, k, alpha, a, lda,
                        b + (long)col * ldb, 1, beta, c + (long)col * ldc, 1);
        }
        return;
    }
    arc_blas_threads_for((double)m * n * k);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, a,
                lda, b, ldb, beta, c, ldc);
}

void arc_gemm_tn(int m, int n, int k, double alpha, const double *a, int lda,
                 const double *b, int ldb, double beta, double *c, int ldc)
{
    if (by_columns(n, (double)m * k)) {
        arc_blas_threads_for((double)m * k);
        for (int col = 0; col < n; col++) {
            cblas_dgemv(CblasColMajor, CblasTrans, k, m, alpha, a, lda,
                        b + (long)col * ldb, 1, beta, c + (long)col * ldc, 1);
        }
        return;
    }
    arc_blas_threads_for((double)m * n * k);
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, m, n, k, alpha, a,
                lda, b, ldb, beta, c, ldc);
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

static void trsm_left(enum CBLAS_TRANSPOSE trans, int n, int nrhs,
                      const double *l, int ldl, double *x, int ldx)
{
    if (by_columns(nrhs, (double)n * n)) {
        arc_blas_threads_for((double)n * n);
        for (int r = 0; r < nrhs; r++) {
            cblas_dtrsv(CblasColMajor, CblasLower, trans, CblasNonUnit, n, l,
                        ldl, x + (long)r * ldx, 1);
        }
        return;
    }
    arc_blas_threads_for((double)n * n * nrhs);
    cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, trans, CblasNonUnit, n,
                nrhs, 1.0, l, ldl, x, ldx);
}

void arc_trsm_left_ln(int n, int nrhs, const double *l, int ldl, double *x,
                      int ldx)
{
    trsm_left(CblasNoTrans, n, nrhs, l, ldl, x, ldx);
}

void arc_trsm_left_lt(int n, int nrhs, const double *l, int ldl, double *x,
                      int ldx)
{
    trsm_left(CblasTrans, n, nrhs, l, ldl, x, ldx);
}
