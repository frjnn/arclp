#include "etree.h"
#include "symbolic_steps.h"

static void invert(int m, const int *perm, int *pinv)
{
    for (int k = 0; k < m; k++) {
        pinv[perm[k]] = k;
    }
}

int arc_order(const arc_graph *g, int *perm, int *pinv, int *parent,
              int *count, int *work)
{
    const int m = g->n;
    int *amd = count;   /* count doubles as scratch until the last step */
    if (arc_amd_order(g, amd) != 0) {
        return -1;
    }
    invert(m, amd, pinv);
    arc_etree(g, amd, pinv, parent, work);
    int *post = work + 3 * m;
    arc_etree_postorder(m, parent, post, work);
    for (int k = 0; k < m; k++) {
        perm[k] = amd[post[k]];
    }
    invert(m, perm, pinv);
    arc_etree(g, perm, pinv, parent, work);
    arc_column_counts(g, perm, pinv, parent, count, work);
    return 0;
}
