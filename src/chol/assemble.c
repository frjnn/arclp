#include "numeric.h"

/* Column k of P A D A' P' below the diagonal is the sum over the entries
 * a_kj of row k of  a_kj d_j  times the part of column j at labels >= k. */
void arc_assemble_supernode(const arc_symbolic *sym, int s, const double *d,
                            const int *map, double *block)
{
    const arc_csc *ap = &sym->a_perm;
    const int nrows = arc_sn_rows(sym, s);
    const int first = sym->first[s];
    for (int c = 0; c < arc_sn_cols(sym, s); c++) {
        const int k = first + c;
        double *col = block + (long)c * nrows;
        for (int t = sym->by_row_ptr[k]; t < sym->by_row_ptr[k + 1]; t++) {
            const int j = sym->by_row_col[t];
            const int q = sym->by_row_pos[t];
            const double w = ap->x[q] * d[j];
            for (int r = q; r < ap->p[j + 1]; r++) {
                col[map[ap->i[r]]] += w * ap->x[r];
            }
        }
    }
}
