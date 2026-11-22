#include "numeric.h"

#include "util/parallel.h"

/* Column c of supernode s: the sum over the entries a_kj of row k of
 * a_kj d_j times the part of column j at labels >= k. */
static inline void assemble_column(const arc_symbolic *sym, int k,
                                   const double *d, const int *map,
                                   double *col)
{
    const arc_csc *ap = &sym->a_perm;
    for (int t = sym->by_row_ptr[k]; t < sym->by_row_ptr[k + 1]; t++) {
        const int j = sym->by_row_col[t];
        const int q = sym->by_row_pos[t];
        const double w = ap->x[q] * d[j];
        for (int r = q; r < ap->p[j + 1]; r++) {
            col[map[ap->i[r]]] += w * ap->x[r];
        }
    }
}

/* Columns of P A D A' P' owned by supernode s, each written by itself. */
void arc_assemble_supernode(const arc_symbolic *sym, int s, const double *d,
                            const int *map, double *block)
{
    const int nrows = arc_sn_rows(sym, s);
    const int first = sym->first[s];
    const int ncols = arc_sn_cols(sym, s);
    const long terms = sym->by_row_ptr[first + ncols] - sym->by_row_ptr[first];
    ARC_FOR(c, ncols, 8 * terms,
            assemble_column(sym, first + c, d, map,
                            block + (long)c * nrows));
}
