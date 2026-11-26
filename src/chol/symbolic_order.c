/*
 * Fill-reducing ordering.  AMD always; when its factor is large, the other
 * available methods too, keeping the one whose stored factor (supernodes
 * after amalgamation, known before any numeric work) is smallest.
 * Factorization time follows the stored size more closely than the flop
 * count, because scattering the updates is bound by memory traffic; for
 * small factors the extra analysis would cost more than it can save.
 */
#include "etree.h"
#include "supernodes.h"
#include "symbolic_steps.h"

#include <stdlib.h>
#include <string.h>

typedef int (*ordering_method)(const arc_graph *g, int *perm);

static const long TRY_OTHERS_ABOVE = 1000000;   /* stored entries */

static const ordering_method ORDERINGS[] = {
    arc_amd_order,
#if defined(ARCLP_USE_METIS)
    arc_metis_order,
#endif
};

static void invert(int m, const int *perm, int *pinv)
{
    for (int k = 0; k < m; k++)
    {
        pinv[perm[k]] = k;
    }
}

typedef struct
{
    int *perm, *pinv, *parent, *count;
} ordering;

/* Postordered permutation, its elimination tree and column counts. */
static int postordered(const arc_graph *g, ordering_method method,
                       ordering o, int *work)
{
    const int m = g->n;
    int *raw = o.count;   /* count doubles as scratch until the last step */
    if (method(g, raw) != 0)
    {
        return -1;
    }
    invert(m, raw, o.pinv);
    arc_etree(g, raw, o.pinv, o.parent, work);
    int *post = work + 3 * m;
    arc_etree_postorder(m, o.parent, post, work);
    for (int k = 0; k < m; k++)
    {
        o.perm[k] = raw[post[k]];
    }
    invert(m, o.perm, o.pinv);
    arc_etree(g, o.perm, o.pinv, o.parent, work);
    arc_column_counts(g, o.perm, o.pinv, o.parent, o.count, work);
    return 0;
}

/* Entries the supernodal factor will store for this ordering; first and
 * rows need m + 1 integers, work 5m + 1. */
static long stored_size(int m, ordering o, const arc_relax_policy *policy,
                        int *first, int *rows, int *work)
{
    const int ns = arc_partition_supernodes(m, o.parent, o.count, policy,
                                            first, rows, work);
    long entries = 0;
    for (int k = 0; k < ns; k++)
    {
        entries += (long)rows[k] * (first[k + 1] - first[k]);
    }
    return entries;
}

static void keep(int m, ordering into, ordering from)
{
    const size_t bytes = sizeof(int) * (size_t)m;
    memcpy(into.perm, from.perm, bytes);
    memcpy(into.pinv, from.pinv, bytes);
    memcpy(into.parent, from.parent, bytes);
    memcpy(into.count, from.count, bytes);
}

int arc_order(const arc_graph *g, const arc_relax_policy *policy,
              int *perm, int *pinv, int *parent, int *count, int *work)
{
    const int m = g->n;
    const ordering best = {perm, pinv, parent, count};
    const size_t methods = sizeof(ORDERINGS) / sizeof(*ORDERINGS);
    if (postordered(g, ORDERINGS[0], best, work) != 0)
    {
        return -1;
    }
    int *block = methods > 1 ? malloc(sizeof(int) * 6 * (size_t)(m + 1))
                             : NULL;
    if (block == NULL)
    {
        return 0;   /* keep AMD */
    }
    int *first = block + 4 * (m + 1);
    int *rows = block + 5 * (m + 1);
    long smallest = stored_size(m, best, policy, first, rows, work);
    const ordering other = {block, block + m + 1, block + 2 * (m + 1),
                            block + 3 * (m + 1)};
    for (size_t k = 1; k < methods && smallest >= TRY_OTHERS_ABOVE; k++)
    {
        if (postordered(g, ORDERINGS[k], other, work) != 0)
        {
            continue;
        }
        const long size = stored_size(m, other, policy, first, rows, work);
        if (size < smallest)
        {
            keep(m, best, other);
            smallest = size;
        }
    }
    free(block);
    return 0;
}
