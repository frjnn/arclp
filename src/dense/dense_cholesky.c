/* Recursive blocked Cholesky of a supernode's dense trapezoid. */
#include "dense.h"

#include <math.h>

enum { LEAF = 32 };

/* Unblocked factorization of an nb x nb diagonal leaf. */
static int factor_diagonal(int nb, double *a, int lda, const double *floor,
                           double huge)
{
    int replaced = 0;
    for (int j = 0; j < nb; j++) {
        double *col = a + (long)j * lda;
        double pivot = col[j];
        if (!(pivot > floor[j])) {
            pivot = huge;
            replaced++;
        }
        const double ljj = sqrt(pivot);
        const double inv = 1.0 / ljj;
        col[j] = ljj;
        for (int i = j + 1; i < nb; i++) {
            col[i] *= inv;
        }
        for (int k = j + 1; k < nb; k++) {
            double *target = a + (long)k * lda;
            const double t = col[k];
            for (int i = k; i < nb; i++) {
                target[i] -= col[i] * t;
            }
        }
    }
    return replaced;
}

/* Recursive column split: every flop outside the leaves is BLAS-3. */
int arc_dense_cholesky(int m, int n, double *a, int lda, const double *floor,
                       double huge)
{
    if (n <= LEAF) {
        const int replaced = factor_diagonal(n, a, lda, floor, huge);
        if (m > n) {
            arc_trsm_right_lt(m - n, n, a, lda, a + n, lda);
        }
        return replaced;
    }
    const int n1 = n / 2 > 256 ? n / 2 / 256 * 256 : n / 2;
    const int n2 = n - n1;
    int replaced = arc_dense_cholesky(m, n1, a, lda, floor, huge);
    const double *left = a + n1;                 /* rows n1..m, cols 0..n1 */
    double *right = a + n1 + (long)n1 * lda;     /* rows n1..m, cols n1..n */
    arc_syrk_ln(n2, n1, -1.0, left, lda, 1.0, right, lda);
    if (m - n1 > n2) {
        arc_gemm_nt(m - n1 - n2, n2, n1, -1.0, left + n2, lda, left, lda, 1.0,
                    right + n2, lda);
    }
    replaced += arc_dense_cholesky(m - n1, n2, right, lda, floor + n1, huge);
    return replaced;
}
