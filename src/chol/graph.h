/* Symmetric sparsity pattern stored as adjacency lists, no self loops. */
#ifndef ARC_GRAPH_H
#define ARC_GRAPH_H

#include "sparse/csc.h"

typedef struct {
    int n;
    int *p;   /* n + 1 list pointers   */
    int *i;   /* neighbours            */
} arc_graph;

/* Builds the pattern of A A' (an m-vertex graph).  Returns 0 on success. */
int arc_graph_of_aat(const arc_csc *a, arc_graph *g);

void arc_graph_free(arc_graph *g);

/* Fill-reducing ordering: perm[k] is the vertex eliminated k-th. */
int arc_amd_order(const arc_graph *g, int *perm);

#if defined(ARCLP_USE_METIS)
/* Nested dissection by METIS; same convention as arc_amd_order. */
int arc_metis_order(const arc_graph *g, int *perm);
#endif

#endif
