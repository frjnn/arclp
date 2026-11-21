/* Elimination tree utilities for a symmetric pattern under a permutation. */
#ifndef ARC_ETREE_H
#define ARC_ETREE_H

#include "graph.h"

/* parent[k] of label k in the elimination tree of P G P'.  perm[k] is the
 * vertex labelled k and pinv its inverse; work holds n integers. */
void arc_etree(const arc_graph *g, const int *perm, const int *pinv,
               int *parent, int *work);

/* Appends the postorder of the tree rooted at j to post[k..] and returns
 * the new k.  Children are linked by head/next (head is consumed); stack
 * needs room for the tree height. */
int arc_tree_dfs(int j, int k, int *head, const int *next, int *post,
                 int *stack);

/* post[k] is the k-th label in a postorder of the forest; work: 3n ints. */
void arc_etree_postorder(int n, const int *parent, int *post, int *work);

/* count[k] = nonzeros in column k of the Cholesky factor (with diagonal). */
void arc_column_counts(const arc_graph *g, const int *perm, const int *pinv,
                       const int *parent, int *count, int *mark);

#endif
