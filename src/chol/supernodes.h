/* Partition of a postordered elimination tree into (relaxed) supernodes. */
#ifndef ARC_SUPERNODES_H
#define ARC_SUPERNODES_H

/* Amalgamation limits in the style of CHOLMOD: a merged supernode with
 * `cols` columns is accepted while its fraction of explicit zeros stays
 * below the threshold of the first size class it fits in. */
typedef struct
{
    int small_cols;       /* always merge up to this many columns */
    int mid_cols;
    double mid_zeros;
    int large_cols;
    double large_zeros;
    double any_zeros;
} arc_relax_policy;

arc_relax_policy arc_default_relax_policy(void);

/* Writes the first column of each supernode to first[] (plus first[ns] = n)
 * and returns ns.  rows[s] receives the row count of supernode s.
 * work needs 5n + 1 integers. */
int arc_partition_supernodes(int n, const int *parent, const int *count,
                             const arc_relax_policy *policy,
                             int *first, int *rows, int *work);

#endif
