/*
 * Dense column-major kernels used by the supernodal factorization.  Two
 * interchangeable implementations exist: dense_blas.c (any CBLAS, e.g.
 * Accelerate or OpenBLAS) and dense_ref.c (portable C).  The build selects
 * one; callers never know which.
 */
#ifndef ARC_DENSE_H
#define ARC_DENSE_H

/* Bracket a sequence of kernel calls.  A backend may tune library state
 * (e.g. BLAS threading) in between; arc_dense_end restores it. */
void arc_dense_begin(void);
void arc_dense_end(void);

typedef enum
{
    ARC_PLAIN,
    ARC_TRANSPOSED
} arc_op;

typedef struct
{
    const double *x;
    int ld;
    arc_op op;
} arc_factor;

/* C := alpha op(A) op(B) + beta C;  C is m x n, the inner dimension is k. */
void arc_gemm(int m, int n, int k, double alpha, arc_factor a, arc_factor b,
              double beta, double *c, int ldc);

/* Lower triangle of C := alpha A A' + beta C;  A is n x k. */
void arc_syrk_ln(int n, int k, double alpha, const double *a, int lda,
                 double beta, double *c, int ldc);

/* B := B L^{-T};  B is m x n, L is n x n lower triangular. */
void arc_trsm_right_lt(int m, int n, const double *l, int ldl, double *b,
                       int ldb);

/* X := op(L)^{-1} X;  L is n x n lower triangular, X is n x nrhs. */
void arc_trsm_left(arc_op op, int n, int nrhs, const double *l, int ldl,
                   double *x, int ldx);

/* In-place Cholesky of the m x n trapezoid A (m >= n): the top n x n block
 * becomes L11 and the rest A21 L11^{-T}.  A pivot of column j not
 * exceeding floor[j] is replaced by `huge`, which removes that direction
 * from the solution (the classic remedy for rank-deficient A in
 * interior-point methods).  Returns the number of replaced pivots. */
int arc_dense_cholesky(int m, int n, double *a, int lda, const double *floor,
                       double huge);

#endif
