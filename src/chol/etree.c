#include "etree.h"

/* Liu's algorithm with path compression through `ancestor`. */
void arc_etree(const arc_graph *g, const int *perm, const int *pinv,
               int *parent, int *ancestor)
{
    for (int k = 0; k < g->n; k++)
    {
        const int v = perm[k];
        parent[k] = -1;
        ancestor[k] = -1;
        for (int t = g->p[v]; t < g->p[v + 1]; t++)
        {
            int i = pinv[g->i[t]];
            while (i != -1 && i < k)
            {
                const int next = ancestor[i];
                ancestor[i] = k;
                if (next == -1)
                {
                    parent[i] = k;
                }
                i = next;
            }
        }
    }
}

int arc_tree_dfs(int j, int k, int *head, const int *next, int *post,
                 int *stack)
{
    int top = 0;
    stack[0] = j;
    while (top >= 0)
    {
        const int p = stack[top];
        const int child = head[p];
        if (child == -1)
        {
            top--;
            post[k++] = p;
        }
        else
        {
            head[p] = next[child];
            stack[++top] = child;
        }
    }
    return k;
}

void arc_etree_postorder(int n, const int *parent, int *post, int *work)
{
    int *head = work;
    int *next = work + n;
    for (int j = 0; j < n; j++)
    {
        head[j] = -1;
    }
    for (int j = n - 1; j >= 0; j--)
    {
        if (parent[j] != -1)
        {
            next[j] = head[parent[j]];
            head[parent[j]] = j;
        }
    }
    for (int root = 0, k = 0; root < n; root++)
    {
        if (parent[root] == -1)
        {
            k = arc_tree_dfs(root, k, head, next, post, work + 2 * n);
        }
    }
}

/* Row k of L is the union of tree paths from each neighbour i < k to k. */
void arc_column_counts(const arc_graph *g, const int *perm, const int *pinv,
                       const int *parent, int *count, int *mark)
{
    for (int k = 0; k < g->n; k++)
    {
        count[k] = 1;
        mark[k] = -1;
    }
    for (int k = 0; k < g->n; k++)
    {
        const int v = perm[k];
        mark[k] = k;
        for (int t = g->p[v]; t < g->p[v + 1]; t++)
        {
            for (int i = pinv[g->i[t]]; i < k && mark[i] != k; i = parent[i])
            {
                count[i]++;
                mark[i] = k;
            }
        }
    }
}
