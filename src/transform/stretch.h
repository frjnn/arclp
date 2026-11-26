/*
 * Dense columns.  One column with L entries makes an L x L block of A A'
 * dense.  Stretching replaces it by copies that hold about PIECE entries
 * each and are tied together by rows x^t - x^(t+1) = 0, which leaves small
 * cliques and a chain.
 */
#ifndef ARC_STRETCH_H
#define ARC_STRETCH_H

#include "sparse/lp.h"
#include "util/arena.h"

typedef struct
{
    arc_lp lp;   /* the stretched problem                       */
    int *p, *i;  /* writable storage behind lp.a                */
    int *first;  /* its column for each original one            */
    int rows, nnz; /* sizes, known before the arrays           */
    int original;  /* columns of the problem being stretched   */
    arc_arena arena;
} arc_stretched;

/* Returns 1 and fills s when some column was stretched, 0 when none needs
 * it, -1 when memory runs out. */
int arc_stretch(const arc_lp *lp, arc_stretched *s);
void arc_stretched_free(arc_stretched *s);

#endif
