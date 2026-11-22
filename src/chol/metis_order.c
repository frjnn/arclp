/* Nested-dissection ordering through METIS (built when ARCLP_USE_METIS). */
#include "graph.h"

#include <metis.h>
#include <stdlib.h>

int arc_metis_order(const arc_graph *g, int *perm)
{
    idx_t n = g->n;
    const int edges = g->p[g->n];
    if (n == 0 || edges == 0) {
        return -1;   /* nothing to dissect; the other orderings suffice */
    }
    idx_t *xadj = malloc(sizeof(idx_t) * ((size_t)n + 1));
    idx_t *adjncy = malloc(sizeof(idx_t) * (size_t)edges);
    idx_t *order = malloc(sizeof(idx_t) * (size_t)n * 2);
    int status = -1;
    if (xadj != NULL && adjncy != NULL && order != NULL) {
        for (int k = 0; k <= g->n; k++) {
            xadj[k] = g->p[k];
        }
        for (int k = 0; k < edges; k++) {
            adjncy[k] = g->i[k];
        }
        idx_t options[METIS_NOPTIONS];
        METIS_SetDefaultOptions(options);
        if (METIS_NodeND(&n, xadj, adjncy, NULL, options, order,
                         order + n) == METIS_OK) {
            for (int k = 0; k < g->n; k++) {
                perm[k] = (int)order[k];
            }
            status = 0;
        }
    }
    free(xadj);
    free(adjncy);
    free(order);
    return status;
}
