#include "supernodes.h"

arc_relax_policy arc_default_relax_policy(void)
{
    const arc_relax_policy policy = {4, 16, 0.8, 48, 0.1, 0.05};
    return policy;
}

typedef struct {
    long cols;
    long rows;
    long nonzeros;   /* structural nonzeros of the fundamental members */
} group;

static long trapezoid(long cols, long rows)
{
    return cols * rows - cols * (cols - 1) / 2;
}

static int accept(const arc_relax_policy *policy, const group *g)
{
    if (g->cols <= policy->small_cols) {
        return 1;
    }
    const long entries = trapezoid(g->cols, g->rows);
    const double zeros = (double)(entries - g->nonzeros) / (double)entries;
    return (g->cols <= policy->mid_cols && zeros < policy->mid_zeros)
        || (g->cols <= policy->large_cols && zeros < policy->large_zeros)
        || zeros < policy->any_zeros;
}

/* Fundamental supernodes: chains with nested structure and single child. */
static int fundamental(int n, const int *parent, const int *count,
                       int *nchild, int *ffirst)
{
    for (int j = 0; j < n; j++) {
        nchild[j] = 0;
    }
    for (int j = 0; j < n; j++) {
        if (parent[j] != -1) {
            nchild[parent[j]]++;
        }
    }
    int nf = 0;
    for (int j = 0; j < n; j++) {
        const int extends = j > 0 && parent[j - 1] == j
                         && count[j - 1] == count[j] + 1 && nchild[j] == 1;
        if (!extends) {
            ffirst[nf++] = j;
        }
    }
    ffirst[nf] = n;
    return nf;
}

static group fundamental_group(const int *ffirst, const int *count, int f)
{
    const long cols = ffirst[f + 1] - ffirst[f];
    const long rows = count[ffirst[f]];
    const group g = {cols, rows, trapezoid(cols, rows)};
    return g;
}

/* Walks the fundamental supernodes from the root end, merging each one
 * into its parent (always the next one) while the policy accepts; marks
 * group starts and records the row count of every group. */
static void relax(int nf, const int *ffirst, const int *owner,
                  const int *parent, const int *count,
                  const arc_relax_policy *policy, int *starts, int *grows)
{
    group current = fundamental_group(ffirst, count, nf - 1);
    for (int f = nf - 2; f >= 0; f--) {
        const int last = ffirst[f + 1] - 1;
        const group mine = fundamental_group(ffirst, count, f);
        const group merged = {mine.cols + current.cols,
                              mine.cols + current.rows,
                              mine.nonzeros + current.nonzeros};
        if (parent[last] != -1 && owner[parent[last]] == f + 1
            && accept(policy, &merged)) {
            starts[f + 1] = 0;
            current = merged;
        } else {
            grows[f + 1] = (int)current.rows;
            current = mine;
        }
    }
    grows[0] = (int)current.rows;
}

int arc_partition_supernodes(int n, const int *parent, const int *count,
                             const arc_relax_policy *policy,
                             int *first, int *rows, int *work)
{
    if (n == 0) {
        first[0] = 0;
        return 0;
    }
    int *ffirst = work;               /* n + 1 */
    int *owner = work + n + 1;        /* fundamental supernode of a column */
    int *starts = work + 2 * n + 1;   /* 1 if f begins a relaxed supernode */
    int *grows = work + 3 * n + 1;    /* rows of the relaxed group at f   */
    const int nf = fundamental(n, parent, count, work + 4 * n + 1, ffirst);
    for (int f = 0; f < nf; f++) {
        for (int j = ffirst[f]; j < ffirst[f + 1]; j++) {
            owner[j] = f;
        }
        starts[f] = 1;
    }
    relax(nf, ffirst, owner, parent, count, policy, starts, grows);
    int ns = 0;
    for (int f = 0; f < nf; f++) {
        if (starts[f]) {
            first[ns] = ffirst[f];
            rows[ns++] = grows[f];
        }
    }
    first[ns] = n;
    return ns;
}
