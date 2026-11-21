#include "symbolic.h"

#include "symbolic_steps.h"

#include <stdlib.h>
#include <string.h>

/* Temporary arrays of the analysis; freed when it ends. */
typedef struct {
    int *perm, *pinv, *parent, *count, *first, *srows, *work;
    int *block;
} scratch;

static int scratch_open(scratch *t, int m, int n)
{
    const size_t size = (size_t)6 * (size_t)(m + 1) + (size_t)(5 * m + 1 + n);
    t->block = malloc(sizeof(int) * size);
    if (t->block == NULL) {
        return -1;
    }
    int *cursor = t->block;
    int **parts[6] = {&t->perm, &t->pinv, &t->parent, &t->count, &t->srows,
                      &t->first};
    for (int k = 0; k < 6; k++) {
        *parts[k] = cursor;
        cursor += m + 1;
    }
    t->work = cursor;
    return 0;
}

static void layout(arc_arena *ar, arc_symbolic *s, const arc_csc *a,
                   long total_rows)
{
    const size_t m = (size_t)a->m;
    const size_t nnz = (size_t)a->p[a->n];
    const size_t ns = (size_t)s->nsuper;
    s->perm = arc_arena_take(ar, m, sizeof(int));
    s->pinv = arc_arena_take(ar, m, sizeof(int));
    s->first = arc_arena_take(ar, ns + 1, sizeof(int));
    s->owner = arc_arena_take(ar, m, sizeof(int));
    s->row_ptr = arc_arena_take(ar, ns + 1, sizeof(int));
    s->rows = arc_arena_take(ar, (size_t)total_rows, sizeof(int));
    s->val_ptr = arc_arena_take(ar, ns + 1, sizeof(long));
    s->ap_p = arc_arena_take(ar, (size_t)a->n + 1, sizeof(int));
    s->ap_i = arc_arena_take(ar, nnz, sizeof(int));
    s->ap_x = arc_arena_take(ar, nnz, sizeof(double));
    const arc_csc view = {a->m, a->n, s->ap_p, s->ap_i, s->ap_x};
    s->a_perm = view;
    s->by_row_ptr = arc_arena_take(ar, m + 1, sizeof(int));
    s->by_row_col = arc_arena_take(ar, nnz, sizeof(int));
    s->by_row_pos = arc_arena_take(ar, nnz, sizeof(int));
}

static int allocate(arc_symbolic *s, const arc_csc *a, const scratch *t)
{
    long total_rows = 0;
    for (int k = 0; k < s->nsuper; k++) {
        total_rows += t->srows[k];
    }
    arc_arena_measure(&s->arena);
    layout(&s->arena, s, a, total_rows);
    if (arc_arena_commit(&s->arena) != 0) {
        return -1;
    }
    layout(&s->arena, s, a, total_rows);
    memcpy(s->perm, t->perm, sizeof(int) * (size_t)s->m);
    memcpy(s->pinv, t->pinv, sizeof(int) * (size_t)s->m);
    memcpy(s->first, t->first, sizeof(int) * (size_t)(s->nsuper + 1));
    return 0;
}

static int analyse(const arc_csc *a, const arc_relax_policy *policy,
                   arc_symbolic *s, const arc_graph *g, scratch *t)
{
    if (arc_order(g, t->perm, t->pinv, t->parent, t->count, t->work) != 0) {
        return -1;
    }
    s->nsuper = arc_partition_supernodes(a->m, t->parent, t->count, policy,
                                         t->first, t->srows, t->work);
    if (allocate(s, a, t) != 0
        || arc_build_structure(g, t->parent, t->srows, s, t->work) != 0) {
        return -1;
    }
    arc_build_assembly(a, s, t->work);
    return 0;
}

int arc_symbolic_analyse(const arc_csc *a, const arc_relax_policy *policy,
                         arc_symbolic *s)
{
    memset(s, 0, sizeof(*s));
    s->m = a->m;
    arc_graph g = {0};
    scratch t = {0};
    int status = -1;
    if (arc_graph_of_aat(a, &g) == 0 && scratch_open(&t, a->m, a->n) == 0) {
        status = analyse(a, policy, s, &g, &t);
    }
    free(t.block);
    arc_graph_free(&g);
    if (status != 0) {
        arc_symbolic_free(s);
    }
    return status;
}

void arc_symbolic_free(arc_symbolic *s)
{
    arc_arena_release(&s->arena);
    memset(s, 0, sizeof(*s));
}
