#include "symbolic_steps.h"

/* Row-wise index by new label, holding original entry positions for now. */
static void index_rows(const arc_csc *a, arc_symbolic *s, int *next)
{
    int *ptr = s->by_row_ptr;
    for (int k = 0; k <= s->m; k++)
    {
        ptr[k] = 0;
    }
    for (int t = 0; t < a->p[a->n]; t++)
    {
        ptr[s->pinv[a->i[t]] + 1]++;
    }
    for (int k = 0; k < s->m; k++)
    {
        ptr[k + 1] += ptr[k];
        next[k] = ptr[k];
    }
    for (int j = 0; j < a->n; j++)
    {
        for (int t = a->p[j]; t < a->p[j + 1]; t++)
        {
            const int slot = next[s->pinv[a->i[t]]]++;
            s->by_row_col[slot] = j;
            s->by_row_pos[slot] = t;
        }
    }
}

/* Visiting rows in label order yields columns sorted by label. */
void arc_build_assembly(const arc_csc *a, arc_symbolic *s, int *work)
{
    int *fill = work;
    index_rows(a, s, work);
    for (int j = 0; j <= a->n; j++)
    {
        s->ap_p[j] = a->p[j];
    }
    for (int j = 0; j < a->n; j++)
    {
        fill[j] = a->p[j];
    }
    for (int k = 0; k < s->m; k++)
    {
        for (int slot = s->by_row_ptr[k]; slot < s->by_row_ptr[k + 1]; slot++)
        {
            const int pos = fill[s->by_row_col[slot]]++;
            s->ap_i[pos] = k;
            s->ap_x[pos] = a->x[s->by_row_pos[slot]];
            s->by_row_pos[slot] = pos;
        }
    }
}
