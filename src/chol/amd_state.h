/*
 * Approximate minimum degree ordering (Amestoy, Davis and Duff, SIMAX 1996;
 * presentation of Davis, "Direct Methods for Sparse Linear Systems", ch. 7).
 *
 * The quotient graph lives in one index array `ci`.  Each live object j
 * (variable or element) owns ci[cp[j] .. cp[j] + len[j]); for a variable the
 * first elen[j] entries are elements, the rest variables.  Absorbed objects
 * store flip(parent) in cp, building the assembly tree used to postorder.
 */
#ifndef ARC_AMD_STATE_H
#define ARC_AMD_STATE_H

typedef struct
{
    int n;
    int *cp, *ci;            /* quotient graph (cp has n + 1 entries)       */
    int capacity, used;      /* size of ci and its occupied prefix          */
    int *len, *nv, *next, *last, *head, *elen, *degree, *w, *hhead;
    int mindeg, eliminated, mark, lemax, dense;
} amd_state;

/* Current pivot element being built. */
typedef struct
{
    int k;          /* pivot variable, becomes the new element  */
    int elenk;      /* number of elements adjacent to k         */
    int nvk;        /* supervariable size of k                  */
    int begin, end; /* the element's variables: ci[begin..end)  */
    int degree;     /* |Lk|, external degree of the element     */
} amd_pivot;

static inline int amd_flip(int i) { return -i - 2; }

int amd_init(amd_state *s, int n, const int *gp, const int *gi);
void amd_free(amd_state *s);
int amd_clear_marks(amd_state *s, int mark);

/* Puts variable i first in the list of degree d (its `last` is untouched). */
void amd_file_by_degree(amd_state *s, int i, int d);
void amd_unlink_degree(amd_state *s, int i);
amd_pivot amd_select_pivot(amd_state *s);
void amd_compact(amd_state *s);
void amd_build_element(amd_state *s, amd_pivot *e);
void amd_set_differences(amd_state *s, const amd_pivot *e);
void amd_update_degrees(amd_state *s, amd_pivot *e);
void amd_merge_indistinguishable(amd_state *s, const amd_pivot *e);
void amd_finalize_element(amd_state *s, amd_pivot *e);
void amd_postorder(amd_state *s, int *perm);

#endif
