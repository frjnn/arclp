#include "stretch.h"

#include <string.h>

enum
{
    DENSE = 2000,
    PIECE = 100
};

static int pieces_of(const arc_csc *a, int j)
{
    const int len = a->p[j + 1] - a->p[j];
    return len > DENSE ? (len + PIECE - 1) / PIECE : 1;
}

/* The stretched problem's sizes; returns the number of links. */
static int measure(const arc_lp *lp, arc_stretched *s)
{
    arc_lp *out = &s->lp;
    int links = 0;
    out->n_free = lp->n_free;
    out->n_bounded = lp->n_bounded;
    for (int j = 0; j < lp->a.n; j++)
    {
        const int extra = pieces_of(&lp->a, j) - 1;
        links += extra;
        out->n_free += j < lp->n_free ? extra : 0;
        out->n_bounded += j >= arc_lp_first_bounded(lp) ? extra : 0;
    }
    s->original = lp->a.n;
    s->rows = lp->a.m + links;
    s->nnz = lp->a.p[lp->a.n] + 2 * links;
    out->a.m = s->rows;
    out->a.n = lp->a.n + links;
    return links;
}

static void layout(arc_arena *ar, void *owner)
{
    arc_stretched *s = owner;
    arc_lp *out = &s->lp;
    out->values = arc_arena_take(ar, (size_t)s->nnz, sizeof(double));
    out->b = arc_arena_take(ar, (size_t)s->rows, sizeof(double));
    out->c = arc_arena_take(ar, (size_t)out->a.n, sizeof(double));
    out->upper = arc_arena_take(ar, (size_t)out->n_bounded, sizeof(double));
    s->p = arc_arena_take(ar, (size_t)out->a.n + 1, sizeof(int));
    s->i = arc_arena_take(ar, (size_t)s->nnz, sizeof(int));
    s->first = arc_arena_take(ar, (size_t)s->original, sizeof(int));
    const arc_csc view = {out->a.m, out->a.n, s->p, s->i, out->values};
    out->a = view;
}

typedef struct
{
    int col, slot, link;
} cursor;

/* Piece t of column j: its share of the entries, then the two links. */
static void write_piece(arc_stretched *s, const arc_lp *lp, int j, int t,
                        int pieces, cursor *at)
{
    int *i = s->i;
    const long len = lp->a.p[j + 1] - lp->a.p[j];
    const int from = lp->a.p[j] + (int)(len * t / pieces);
    const int to = lp->a.p[j] + (int)(len * (t + 1) / pieces);
    for (int k = from; k < to; k++)
    {
        i[at->slot] = lp->a.i[k];
        s->lp.values[at->slot++] = lp->a.x[k];
    }
    if (t > 0)
    {
        i[at->slot] = at->link - 1;
        s->lp.values[at->slot++] = -1.0;
    }
    if (t < pieces - 1)
    {
        i[at->slot] = at->link++;
        s->lp.values[at->slot++] = 1.0;
    }
    s->lp.c[at->col] = lp->c[j] / pieces;
    const int bounded = at->col - arc_lp_first_bounded(&s->lp);
    if (bounded >= 0)
    {
        s->lp.upper[bounded] = lp->upper[j - arc_lp_first_bounded(lp)];
    }
    s->p[++at->col] = at->slot;
}

int arc_stretch(const arc_lp *lp, arc_stretched *s)
{
    memset(s, 0, sizeof(*s));
    if (measure(lp, s) == 0)
    {
        return 0;
    }
    if (arc_arena_build(&s->arena, layout, s) != 0)
    {
        return -1;
    }
    memcpy(s->lp.b, lp->b, sizeof(double) * (size_t)lp->a.m);
    cursor at = {0, 0, lp->a.m};
    for (int j = 0; j < lp->a.n; j++)
    {
        const int pieces = pieces_of(&lp->a, j);
        s->first[j] = at.col;
        for (int t = 0; t < pieces; t++)
        {
            write_piece(s, lp, j, t, pieces, &at);
        }
    }
    return 1;
}

void arc_stretched_free(arc_stretched *s)
{
    arc_arena_release(&s->arena);
}
