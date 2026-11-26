#include "amd_state.h"

/* Scan 1: w[e] - mark becomes |Le \ Lk| for every element e adjacent to Lk. */
void amd_set_differences(amd_state *s, const amd_pivot *e)
{
    s->mark = amd_clear_marks(s, s->mark);
    for (int pk = e->begin; pk < e->end; pk++)
    {
        const int i = s->ci[pk];
        const int eln = s->elen[i];
        if (eln <= 0)
        {
            continue;
        }
        const int nvi = -s->nv[i];
        const int wnvi = s->mark - nvi;
        for (int p = s->cp[i]; p < s->cp[i] + eln; p++)
        {
            const int el = s->ci[p];
            if (s->w[el] >= s->mark)
            {
                s->w[el] -= nvi;
            }
            else if (s->w[el] != 0)
            {
                s->w[el] = s->degree[el] + wnvi;
            }
        }
    }
}

/* Prunes the element list of i, absorbing elements now inside Lk.
 * Returns the summed external degree; accumulates a hash of the list. */
static int prune_elements(amd_state *s, int k, int i, int *pn, unsigned *hash)
{
    const int p1 = s->cp[i];
    int d = 0;
    for (int p = p1; p < p1 + s->elen[i]; p++)
    {
        const int el = s->ci[p];
        if (s->w[el] == 0)
        {
            continue;
        }
        const int dext = s->w[el] - s->mark;
        if (dext > 0)
        {
            d += dext;
            s->ci[(*pn)++] = el;
            *hash += (unsigned)el;
        }
        else
        {
            s->cp[el] = amd_flip(k);
            s->w[el] = 0;
        }
    }
    return d;
}

/* Keeps only live variables outside Lk in the variable list of i. */
static int prune_variables(amd_state *s, int i, int old_elen, int *pn,
                           unsigned *hash)
{
    const int p1 = s->cp[i];
    int d = 0;
    for (int p = p1 + old_elen; p < p1 + s->len[i]; p++)
    {
        const int j = s->ci[p];
        const int nvj = s->nv[j];
        if (nvj <= 0)
        {
            continue;
        }
        d += nvj;
        s->ci[(*pn)++] = j;
        *hash += (unsigned)j;
    }
    return d;
}

static void mass_eliminate(amd_state *s, amd_pivot *e, int i)
{
    const int nvi = -s->nv[i];
    s->cp[i] = amd_flip(e->k);
    e->degree -= nvi;
    e->nvk += nvi;
    s->eliminated += nvi;
    s->nv[i] = 0;
    s->elen[i] = -1;
}

/* Re-links i with k as its first element and files it in a hash bucket. */
static void relink(amd_state *s, int k, int i, int p3, int pn, unsigned hash)
{
    const int p1 = s->cp[i];
    s->ci[pn] = s->ci[p3];
    s->ci[p3] = s->ci[p1];
    s->ci[p1] = k;
    s->len[i] = pn - p1 + 1;
    const int h = (int)(hash % (unsigned)s->n);
    s->next[i] = s->hhead[h];
    s->hhead[h] = i;
    s->last[i] = h;
}

/* Scan 2: approximate external degree of every variable in Lk. */
void amd_update_degrees(amd_state *s, amd_pivot *e)
{
    for (int pk = e->begin; pk < e->end; pk++)
    {
        const int i = s->ci[pk];
        unsigned hash = 0;
        const int old_elen = s->elen[i];
        int pn = s->cp[i];
        int d = prune_elements(s, e->k, i, &pn, &hash);
        const int p3 = pn;
        s->elen[i] = pn - s->cp[i] + 1;
        d += prune_variables(s, i, old_elen, &pn, &hash);
        if (d == 0)
        {
            mass_eliminate(s, e, i);
        }
        else
        {
            s->degree[i] = s->degree[i] < d ? s->degree[i] : d;
            relink(s, e->k, i, p3, pn, hash);
        }
    }
    s->degree[e->k] = e->degree;
    s->lemax = s->lemax > e->degree ? s->lemax : e->degree;
    s->mark = amd_clear_marks(s, s->mark + s->lemax);
}
