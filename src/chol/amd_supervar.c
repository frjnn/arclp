#include "amd_state.h"

static int same_adjacency(const amd_state *s, int i, int j)
{
    if (s->len[j] != s->len[i] || s->elen[j] != s->elen[i])
    {
        return 0;
    }
    for (int p = s->cp[j] + 1; p < s->cp[j] + s->len[j]; p++)
    {
        if (s->w[s->ci[p]] != s->mark)
        {
            return 0;
        }
    }
    return 1;
}

/* Merges every j after i in the bucket whose adjacency equals that of i. */
static void merge_bucket_member(amd_state *s, int i)
{
    for (int p = s->cp[i] + 1; p < s->cp[i] + s->len[i]; p++)
    {
        s->w[s->ci[p]] = s->mark;
    }
    int jlast = i;
    for (int j = s->next[i]; j != -1;)
    {
        if (same_adjacency(s, i, j))
        {
            s->cp[j] = amd_flip(i);
            s->nv[i] += s->nv[j];
            s->nv[j] = 0;
            s->elen[j] = -1;
            j = s->next[j];
            s->next[jlast] = j;
        }
        else
        {
            jlast = j;
            j = s->next[j];
        }
    }
}

/* Supervariable detection through the hash buckets filled by scan 2. */
void amd_merge_indistinguishable(amd_state *s, const amd_pivot *e)
{
    for (int pk = e->begin; pk < e->end; pk++)
    {
        int i = s->ci[pk];
        if (s->nv[i] >= 0)
        {
            continue;
        }
        const int h = s->last[i];
        i = s->hhead[h];
        s->hhead[h] = -1;
        for (; i != -1 && s->next[i] != -1; i = s->next[i], s->mark++)
        {
            merge_bucket_member(s, i);
        }
    }
}

/* Restores the surviving variables of Lk and refiles them by degree. */
void amd_finalize_element(amd_state *s, amd_pivot *e)
{
    int p = e->begin;
    for (int pk = e->begin; pk < e->end; pk++)
    {
        const int i = s->ci[pk];
        const int nvi = -s->nv[i];
        if (nvi <= 0)
        {
            continue;
        }
        s->nv[i] = nvi;
        int d = s->degree[i] + e->degree - nvi;
        const int bound = s->n - s->eliminated - nvi;
        d = d < bound ? d : bound;
        amd_file_by_degree(s, i, d);
        s->last[i] = -1;
        s->mindeg = s->mindeg < d ? s->mindeg : d;
        s->degree[i] = d;
        s->ci[p++] = i;
    }
    s->nv[e->k] = e->nvk;
    s->len[e->k] = p - e->begin;
    if (s->len[e->k] == 0)
    {
        s->cp[e->k] = -1;
        s->w[e->k] = 0;
    }
    if (e->elenk != 0)
    {
        s->used = p;
    }
}
