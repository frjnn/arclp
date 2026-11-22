/* Internal stages of arc_symbolic_analyse. */
#ifndef ARC_SYMBOLIC_STEPS_H
#define ARC_SYMBOLIC_STEPS_H

#include "graph.h"
#include "symbolic.h"

/* Postordered fill-reducing permutation (AMD, or METIS nested dissection
 * when available and it gives a smaller stored factor), its elimination
 * tree and column counts.  work needs 5m + 1 integers. */
int arc_order(const arc_graph *g, const arc_relax_policy *policy,
              int *perm, int *pinv, int *parent, int *count, int *work);

/* Fills owner, row_ptr, rows, val_ptr, max_rows and max_update from the
 * partition in sym->first.  srows[s] is the row count of supernode s;
 * work needs 3m integers. */
int arc_build_structure(const arc_graph *g, const int *parent,
                        const int *srows, arc_symbolic *sym, int *work);

/* Fills a_perm and the row-wise index from A and sym->pinv.
 * work needs max(m, n) integers. */
void arc_build_assembly(const arc_csc *a, arc_symbolic *sym, int *work);

#endif
