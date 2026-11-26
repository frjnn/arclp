#include "symbolic_steps.h"

/* Visiting rows in label order yields columns sorted by label. */
void arc_build_assembly(const arc_csc *a, arc_symbolic *s, int *work)
{
    int *fill = work;
    arc_csc_index_rows(a, s->pinv, NULL, s->by_row_ptr, s->by_row_col,
                       s->by_row_pos, work);
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
