#include "graph.h"

#include <limits.h>
#include <stdlib.h>

/* Visits the neighbours of row r in A A'; writes them when out != NULL. */
static int neighbours(const arc_csc *a, const int *rp, const int *cols,
                      int r, int *mark, int *out)
{
    int count = 0;
    mark[r] = r;
    for (int t = rp[r]; t < rp[r + 1]; t++)
    {
        const int j = cols[t];
        for (int k = a->p[j]; k < a->p[j + 1]; k++)
        {
            const int q = a->i[k];
            if (mark[q] != r)
            {
                mark[q] = r;
                if (out != NULL)
                {
                    out[count] = q;
                }
                count++;
            }
        }
    }
    return count;
}

static void reset(int *mark, int m)
{
    for (int r = 0; r < m; r++)
    {
        mark[r] = -1;
    }
}

/* Row pointers of the pattern; 0 when its size exceeds an int. */
static int count_neighbours(const arc_csc *a, const int *rp, const int *cols,
                            int *mark, int *p)
{
    long total = 0;
    p[0] = 0;
    for (int r = 0; r < a->m; r++)
    {
        total += neighbours(a, rp, cols, r, mark, NULL);
        if (total > INT_MAX - 1)
        {
            return 0;
        }
        p[r + 1] = (int)total;
    }
    return 1;
}

int arc_graph_of_aat(const arc_csc *a, arc_graph *g)
{
    const int m = a->m;
    int *rp = malloc(sizeof(int) * (size_t)(2 * m + 1));
    int *cols = malloc(sizeof(int) * (size_t)(a->p[a->n] + 1));
    g->n = m;
    g->p = malloc(sizeof(int) * (size_t)(m + 1));
    g->i = NULL;
    if (rp == NULL || cols == NULL || g->p == NULL)
    {
        free(rp);
        free(cols);
        arc_graph_free(g);
        return -1;
    }
    int *mark = rp + m + 1;
    arc_csc_index_rows(a, NULL, NULL, rp, cols, NULL, mark);
    reset(mark, m);
    const int fits = count_neighbours(a, rp, cols, mark, g->p);
    g->i = fits ? malloc(sizeof(int) * (size_t)(g->p[m] + 1)) : NULL;
    if (g->i != NULL)
    {
        reset(mark, m);
        for (int r = 0; r < m; r++)
        {
            neighbours(a, rp, cols, r, mark, g->i + g->p[r]);
        }
    }
    free(rp);
    free(cols);
    return g->i == NULL ? -1 : 0;
}

void arc_graph_free(arc_graph *g)
{
    free(g->p);
    free(g->i);
    g->p = NULL;
    g->i = NULL;
}
