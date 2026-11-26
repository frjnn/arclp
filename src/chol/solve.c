/* Supernodal forward and backward substitution with several right sides. */
#include "numeric.h"

#include "dense/dense.h"

/* Supernode s of the factor and its part of the right-hand sides. */
typedef struct
{
    int nrows, ncols, below; /* below: rows under the diagonal block */
    const double *l;         /* nrows x ncols                        */
    const int *rows;         /* labels of the rows below             */
    double *x;               /* its entries of the right-hand sides  */
} supernode;

static supernode supernode_at(const arc_numeric *num, int s)
{
    const arc_symbolic *sym = num->sym;
    const int nrows = arc_sn_rows(sym, s);
    const int ncols = arc_sn_cols(sym, s);
    const supernode sn = {nrows, ncols, nrows - ncols,
                          num->lx + sym->val_ptr[s],
                          sym->rows + sym->row_ptr[s] + ncols,
                          num->x + sym->first[s]};
    return sn;
}

static void forward(arc_numeric *num, int nrhs)
{
    const int m = num->sym->m;
    for (int s = 0; s < num->sym->nsuper; s++)
    {
        const supernode sn = supernode_at(num, s);
        arc_trsm_left(ARC_PLAIN, sn.ncols, nrhs, sn.l, sn.nrows, sn.x, m);
        if (sn.below == 0)
        {
            continue;
        }
        double *w = num->gather;
        const arc_factor l = {sn.l + sn.ncols, sn.nrows, ARC_PLAIN};
        const arc_factor x = {sn.x, m, ARC_PLAIN};
        arc_gemm(sn.below, nrhs, sn.ncols, 1.0, l, x, 0.0, w, sn.below);
        for (int r = 0; r < nrhs; r++)
        {
            for (int i = 0; i < sn.below; i++)
            {
                num->x[sn.rows[i] + (long)r * m] -= w[i + (long)r * sn.below];
            }
        }
    }
}

static void backward(arc_numeric *num, int nrhs)
{
    const int m = num->sym->m;
    for (int s = num->sym->nsuper - 1; s >= 0; s--)
    {
        const supernode sn = supernode_at(num, s);
        if (sn.below > 0)
        {
            double *w = num->gather;
            for (int r = 0; r < nrhs; r++)
            {
                for (int i = 0; i < sn.below; i++)
                {
                    w[i + (long)r * sn.below] =
                        num->x[sn.rows[i] + (long)r * m];
                }
            }
            const arc_factor l = {sn.l + sn.ncols, sn.nrows, ARC_TRANSPOSED};
            const arc_factor below = {w, sn.below, ARC_PLAIN};
            arc_gemm(sn.ncols, nrhs, sn.below, -1.0, l, below, 1.0, sn.x, m);
        }
        arc_trsm_left(ARC_TRANSPOSED, sn.ncols, nrhs, sn.l, sn.nrows, sn.x, m);
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
