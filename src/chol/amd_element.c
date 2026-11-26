#include "amd_state.h"

void amd_unlink_degree(amd_state *s, int i)
{
    if (s->next[i] != -1)
    {
        s->last[s->next[i]] = s->last[i];
    }
    if (s->last[i] != -1)
    {
        s->next[s->last[i]] = s->next[i];
    }
    else
    {
        s->head[s->degree[i]] = s->next[i];
    }
}

amd_pivot amd_select_pivot(amd_state *s)
{
    amd_pivot e = {0};
    int k = -1;
    while (s->mindeg < s->n && (k = s->head[s->mindeg]) == -1)
    {
        s->mindeg++;
    }
    if (s->next[k] != -1)
    {
        s->last[s->next[k]] = -1;
    }
    s->head[s->mindeg] = s->next[k];
    e.k = k;
    e.elenk = s->elen[k];
    e.nvk = s->nv[k];
    s->eliminated += e.nvk;
    return e;
}

/* Garbage collection: slides every live object to the front of ci. */
void amd_compact(amd_state *s)
{
    for (int j = 0; j < s->n; j++)
    {
        const int p = s->cp[j];
        if (p >= 0)
        {
            s->cp[j] = s->ci[p];
            s->ci[p] = amd_flip(j);
        }
    }
    int q = 0;
    for (int p = 0; p < s->used;)
    {
        const int j = amd_flip(s->ci[p++]);
        if (j >= 0)
        {
            s->ci[q] = s->cp[j];
            s->cp[j] = q++;
            for (int t = 0; t < s->len[j] - 1; t++)
            {
                s->ci[q++] = s->ci[p++];
            }
        }
    }
    s->used = q;
}

/* Gathers the variables of one object into the new element Lk. */
static void absorb_object(amd_state *s, amd_pivot *e, int obj, int pj, int ln)
{
    for (int t = 0; t < ln; t++)
    {
        const int i = s->ci[pj++];
        const int nvi = s->nv[i];
        if (nvi <= 0)
        {
            continue;
        }
        e->degree += nvi;
        s->nv[i] = -nvi;
        s->ci[e->end++] = i;
        amd_unlink_degree(s, i);
    }
    if (obj != e->k)
    {
        s->cp[obj] = amd_flip(e->k);
        s->w[obj] = 0;
    }
}

/* Lk = (adjacent variables of k) ∪ (variables of elements adjacent to k). */
void amd_build_element(amd_state *s, amd_pivot *e)
{
    const int k = e->k;
    int p = s->cp[k];
    s->nv[k] = -e->nvk;
    e->begin = e->elenk == 0 ? p : s->used;
    e->end = e->begin;
    e->degree = 0;
    for (int t = 0; t < e->elenk; t++)
    {
        const int obj = s->ci[p++];
        absorb_object(s, e, obj, s->cp[obj], s->len[obj]);
    }
    absorb_object(s, e, k, p, s->len[k] - e->elenk);
    if (e->elenk != 0)
    {
        s->used = e->end;
    }
    s->degree[k] = e->degree;
    s->cp[k] = e->begin;
    s->len[k] = e->end - e->begin;
    s->elen[k] = -2;
}
