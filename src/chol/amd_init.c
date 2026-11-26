#include "amd_state.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

enum { AMD_ARRAYS = 10 };

static int dense_threshold(int n)
{
    int dense = (int)fmax(16.0, 10.0 * sqrt((double)n));
    return dense < n - 2 ? dense : n - 2;
}

/* Empty variables are eliminated at once; dense ones are absorbed into the
 * placeholder element n and ordered last. */
static void classify_variables(amd_state *s)
{
    const int n = s->n;
    for (int i = 0; i < n; i++)
    {
        const int d = s->degree[i];
        if (d == 0)
        {
            s->elen[i] = -2;
            s->eliminated++;
            s->cp[i] = -1;
            s->w[i] = 0;
        }
        else if (d > s->dense)
        {
            s->nv[i] = 0;
            s->elen[i] = -1;
            s->eliminated++;
            s->cp[i] = amd_flip(n);
            s->nv[n]++;
        }
        else
        {
            amd_file_by_degree(s, i, d);
        }
    }
}

static int allocate(amd_state *s, int n, int nnz)
{
    s->capacity = nnz + nnz / 5 + 2 * n + 1;
    s->ci = malloc(sizeof(int) * (size_t)s->capacity);
    s->cp = malloc(sizeof(int) * (size_t)(n + 1) * AMD_ARRAYS);
    if (s->ci == NULL || s->cp == NULL)
    {
        amd_free(s);
        return -1;
    }
    int **arrays[AMD_ARRAYS - 1] = {&s->len, &s->nv, &s->next, &s->last,
                                    &s->head, &s->elen, &s->degree, &s->w,
                                    &s->hhead};
    for (int a = 0; a < AMD_ARRAYS - 1; a++)
    {
        *arrays[a] = s->cp + (size_t)(n + 1) * (size_t)(a + 1);
    }
    return 0;
}

int amd_init(amd_state *s, int n, const int *gp, const int *gi)
{
    const int nnz = gp[n];
    memset(s, 0, sizeof(*s));
    s->n = n;
    if (allocate(s, n, nnz) != 0)
    {
        return -1;
    }
    memcpy(s->ci, gi, sizeof(int) * (size_t)nnz);
    for (int k = 0; k <= n; k++)
    {
        s->cp[k] = k < n ? gp[k] : nnz;
        s->len[k] = k < n ? gp[k + 1] - gp[k] : 0;
        s->head[k] = s->last[k] = s->next[k] = s->hhead[k] = -1;
        s->nv[k] = 1;
        s->w[k] = 1;
        s->elen[k] = 0;
        s->degree[k] = s->len[k];
    }
    s->used = nnz;
    s->dense = dense_threshold(n);
    s->mark = amd_clear_marks(s, 0);
    s->elen[n] = -2;
    s->cp[n] = -1;
    s->w[n] = 0;
    classify_variables(s);
    return 0;
}

void amd_free(amd_state *s)
{
    free(s->ci);
    free(s->cp);
    s->ci = s->cp = NULL;
}

/* Ensures every live w[] entry is below the returned mark. */
int amd_clear_marks(amd_state *s, int mark)
{
    if (mark < 2 || (long)mark + s->lemax >= INT_MAX)
    {
        for (int k = 0; k < s->n; k++)
        {
            if (s->w[k] != 0)
            {
                s->w[k] = 1;
            }
        }
        mark = 2;
    }
    return mark;
}
