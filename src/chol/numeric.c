#include "numeric.h"

#include "dense/dense.h"
#include "util/parallel.h"

#include <string.h>

static const double HUGE_PIVOT = 1e128;
/* Chosen on the Netlib benchmark. */
static const double DEFAULT_PIVOT_TOLERANCE = 1e-14;

static void layout(arc_arena *ar, arc_numeric *num)
{
    const arc_symbolic *s = num->sym;
    const size_t ns = (size_t)s->nsuper;
    num->lx = arc_arena_take(ar, (size_t)s->val_ptr[s->nsuper], sizeof(double));
    num->map = arc_arena_take(ar, (size_t)s->m, sizeof(int));
    num->link_head = arc_arena_take(ar, ns, sizeof(int));
    num->link_next = arc_arena_take(ar, ns, sizeof(int));
    num->next_row = arc_arena_take(ar, ns, sizeof(int));
    num->update = arc_arena_take(ar, (size_t)s->max_update, sizeof(double));
    num->x = arc_arena_take(ar, (size_t)s->m * ARC_MAX_RHS, sizeof(double));
    num->gather = arc_arena_take(ar, (size_t)s->max_rows * ARC_MAX_RHS,
                                 sizeof(double));
    num->pivot_floor = arc_arena_take(ar, (size_t)s->m, sizeof(double));
    num->relative = arc_arena_take(ar, (size_t)s->max_rows, sizeof(int));
}

int arc_numeric_create(const arc_symbolic *sym, arc_numeric *num)
{
    memset(num, 0, sizeof(*num));
    num->sym = sym;
    num->pivot_tolerance = DEFAULT_PIVOT_TOLERANCE;
    arc_arena_measure(&num->arena);
    layout(&num->arena, num);
    if (arc_arena_commit(&num->arena) != 0) {
        return -1;
    }
    layout(&num->arena, num);
    return 0;
}

void arc_numeric_free(arc_numeric *num)
{
    arc_arena_release(&num->arena);
    memset(num, 0, sizeof(*num));
}

static void push(arc_numeric *num, int d, int row_pos)
{
    const arc_symbolic *s = num->sym;
    const int target = s->owner[s->rows[s->row_ptr[d] + row_pos]];
    num->next_row[d] = row_pos;
    num->link_next[d] = num->link_head[target];
    num->link_head[target] = d;
}

/* Rows p1..p2 of descendant d fall in the columns of s; rows p2..end lie
 * below them.  C := L_d(p1:end, :) L_d(p1:p2, :)' is n2 x n1 (lower). */
typedef struct {
    int p1, p2, n1, n2;
} update_rows;

static update_rows rows_in(const arc_numeric *num, int d, int s)
{
    const arc_symbolic *sym = num->sym;
    const int *rows = sym->rows + sym->row_ptr[d];
    const int ld = arc_sn_rows(sym, d);
    update_rows u = {num->next_row[d], num->next_row[d], 0, 0};
    while (u.p2 < ld && rows[u.p2] < sym->first[s + 1]) {
        u.p2++;
    }
    u.n1 = u.p2 - u.p1;
    u.n2 = ld - u.p1;
    return u;
}

static void form_update(const arc_numeric *num, int d, update_rows u)
{
    const arc_symbolic *sym = num->sym;
    const int ld = arc_sn_rows(sym, d);
    const int k = arc_sn_cols(sym, d);
    const double *l = num->lx + sym->val_ptr[d];
    arc_syrk_ln(u.n1, k, 1.0, l + u.p1, ld, 0.0, num->update, u.n2);
    if (u.n2 > u.n1) {
        arc_gemm_nt(u.n2 - u.n1, u.n1, k, 1.0, l + u.p2, ld, l + u.p1, ld,
                    0.0, num->update + u.n1, u.n2);
    }
}

/* Column jj of the update lands in its own column of the target block. */
static inline void scatter_column(double *block, int nrows, const int *rel,
                                  const double *update, update_rows u, int jj)
{
    double *target = block + (long)rel[jj] * nrows;
    const double *src = update + (long)jj * u.n2;
    for (int ii = jj; ii < u.n2; ii++) {
        target[rel[ii]] -= src[ii];
    }
}

/* block -= C, through the positions of d's rows inside supernode s. */
static void scatter_update(arc_numeric *num, int d, int s, update_rows u,
                           double *block)
{
    const arc_symbolic *sym = num->sym;
    const int *rows = sym->rows + sym->row_ptr[d] + u.p1;
    const int nrows = arc_sn_rows(sym, s);
    int *rel = num->relative;
    for (int ii = 0; ii < u.n2; ii++) {
        rel[ii] = num->map[rows[ii]];
    }
    ARC_FOR(jj, u.n1, (long)u.n1 * u.n2,
            scatter_column(block, nrows, rel, num->update, u, jj));
}

static void apply_update(arc_numeric *num, int d, int s, double *block)
{
    const update_rows u = rows_in(num, d, s);
    form_update(num, d, u);
    scatter_update(num, d, s, u, block);
    if (u.p2 < arc_sn_rows(num->sym, d)) {
        push(num, d, u.p2);
    }
}

/* floor[k] = tolerance * (P A D A' P')_kk: a pivot below it signals a row
 * that is numerically dependent on the rows eliminated before it. */
static void pivot_floors(arc_numeric *num, const double *d)
{
    const arc_csc *ap = &num->sym->a_perm;
    double *floor = num->pivot_floor;
    for (int k = 0; k < ap->m; k++) {
        floor[k] = 0.0;
    }
    for (int j = 0; j < ap->n; j++) {
        for (int t = ap->p[j]; t < ap->p[j + 1]; t++) {
            floor[ap->i[t]] += ap->x[t] * ap->x[t] * d[j];
        }
    }
    for (int k = 0; k < ap->m; k++) {
        floor[k] *= num->pivot_tolerance;
    }
}

static void factor_supernode(arc_numeric *num, int s, const double *d)
{
    const arc_symbolic *sym = num->sym;
    const int *rows = sym->rows + sym->row_ptr[s];
    const int nrows = arc_sn_rows(sym, s);
    const int ncols = arc_sn_cols(sym, s);
    double *block = num->lx + sym->val_ptr[s];
    memset(block, 0, sizeof(double) * (size_t)nrows * (size_t)ncols);
    for (int r = 0; r < nrows; r++) {
        num->map[rows[r]] = r;
    }
    arc_assemble_supernode(sym, s, d, num->map, block);
    for (int dd = num->link_head[s], next; dd != -1; dd = next) {
        next = num->link_next[dd];
        apply_update(num, dd, s, block);
    }
    num->tiny_pivots += arc_dense_cholesky(nrows, ncols, block, nrows,
                                           num->pivot_floor + sym->first[s],
                                           HUGE_PIVOT);
    if (nrows > ncols) {
        push(num, s, ncols);
    }
}

int arc_numeric_factor(arc_numeric *num, const double *d)
{
    const arc_symbolic *sym = num->sym;
    pivot_floors(num, d);
    num->tiny_pivots = 0;
    arc_dense_begin();
    for (int s = 0; s < sym->nsuper; s++) {
        num->link_head[s] = -1;
    }
    for (int s = 0; s < sym->nsuper; s++) {
        factor_supernode(num, s, d);
    }
    arc_dense_end();
    return 0;
}
