/* Supernodal forward and backward substitution with several right sides. */
#include "numeric.h"

#include "dense/dense.h"

static void forward(arc_numeric *num, int nrhs)
{
    const arc_symbolic *sym = num->sym;
    const int m = sym->m;
    for (int s = 0; s < sym->nsuper; s++)
    {
        const int nrows = arc_sn_rows(sym, s);
        const int ncols = arc_sn_cols(sym, s);
        const int below = nrows - ncols;
        const double *l = num->lx + sym->val_ptr[s];
        const int *rows = sym->rows + sym->row_ptr[s] + ncols;
        double *xs = num->x + sym->first[s];
        arc_trsm_left_ln(ncols, nrhs, l, nrows, xs, m);
        if (below == 0)
        {
            continue;
        }
        double *w = num->gather;
        arc_gemm_nn(below, nrhs, ncols, 1.0, l + ncols, nrows, xs, m, 0.0, w,
                    below);
        for (int r = 0; r < nrhs; r++)
        {
            for (int i = 0; i < below; i++)
            {
                num->x[rows[i] + (long)r * m] -= w[i + (long)r * below];
            }
        }
    }
}

static void backward(arc_numeric *num, int nrhs)
{
    const arc_symbolic *sym = num->sym;
    const int m = sym->m;
    for (int s = sym->nsuper - 1; s >= 0; s--)
    {
        const int nrows = arc_sn_rows(sym, s);
        const int ncols = arc_sn_cols(sym, s);
        const int below = nrows - ncols;
        const double *l = num->lx + sym->val_ptr[s];
        const int *rows = sym->rows + sym->row_ptr[s] + ncols;
        double *xs = num->x + sym->first[s];
        if (below > 0)
        {
            double *w = num->gather;
            for (int r = 0; r < nrhs; r++)
            {
                for (int i = 0; i < below; i++)
                {
                    w[i + (long)r * below] = num->x[rows[i] + (long)r * m];
                }
            }
            arc_gemm_tn(ncols, nrhs, below, -1.0, l + ncols, nrows, w, below,
                        1.0, xs, m);
        }
        arc_trsm_left_lt(ncols, nrhs, l, nrows, xs, m);
    }
}

void arc_numeric_solve(arc_numeric *num, double *b, int nrhs)
{
    const int m = num->sym->m;
    const int *perm = num->sym->perm;
    for (int r = 0; r < nrhs; r++)
    {
        for (int k = 0; k < m; k++)
        {
            num->x[k + (long)r * m] = b[perm[k] + (long)r * m];
        }
    }
    arc_dense_begin();
    forward(num, nrhs);
    backward(num, nrhs);
    arc_dense_end();
    for (int r = 0; r < nrhs; r++)
    {
        for (int k = 0; k < m; k++)
        {
            b[perm[k] + (long)r * m] = num->x[k + (long)r * m];
        }
    }
}
