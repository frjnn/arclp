#include "amd_state.h"
#include "etree.h"
#include "graph.h"

#include <stdlib.h>

/* Orders variables by a postorder of the assembly tree. */
void amd_postorder(amd_state *s, int *perm)
{
    const int n = s->n;
    for (int i = 0; i < n; i++)
    {
        s->cp[i] = amd_flip(s->cp[i]);
    }
    for (int j = 0; j <= n; j++)
    {
        s->head[j] = -1;
    }
    for (int j = n; j >= 0; j--)
    {
        if (s->nv[j] <= 0)
        {
            s->next[j] = s->head[s->cp[j]];
            s->head[s->cp[j]] = j;
        }
    }
    for (int e = n; e >= 0; e--)
    {
        if (s->nv[e] > 0 && s->cp[e] != -1)
        {
            s->next[e] = s->head[s->cp[e]];
            s->head[s->cp[e]] = e;
        }
    }
    for (int k = 0, i = 0; i <= n; i++)
    {
        if (s->cp[i] == -1)
        {
            k = arc_tree_dfs(i, k, s->head, s->next, perm, s->w);
        }
    }
}

static void eliminate_all(amd_state *s)
{
    while (s->eliminated < s->n)
    {
        amd_pivot e = amd_select_pivot(s);
        if (e.elenk > 0 && s->used + s->mindeg >= s->capacity)
        {
            amd_compact(s);
        }
        amd_build_element(s, &e);
        amd_set_differences(s, &e);
        amd_update_degrees(s, &e);
        amd_merge_indistinguishable(s, &e);
        amd_finalize_element(s, &e);
    }
}

int arc_amd_order(const arc_graph *g, int *perm)
{
    const int n = g->n;
    if (n == 0)
    {
        return 0;
    }
    amd_state s;
    int *order = malloc(sizeof(int) * (size_t)(n + 1));
    if (order == NULL || amd_init(&s, n, g->p, g->i) != 0)
    {
        free(order);
        return -1;
    }
    eliminate_all(&s);
    amd_postorder(&s, order);
    for (int k = 0; k < n; k++)
    {
        perm[k] = order[k];
    }
    amd_free(&s);
    free(order);
    return 0;
}
