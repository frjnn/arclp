#include "ipm.h"

#include "util/vec.h"

#include <string.h>

static void layout(arc_arena *ar, void *owner)
{
    arc_ipm *ipm = owner;
    const size_t m = (size_t)ipm->m;
    const size_t n = (size_t)ipm->n;
    const size_t total = (size_t)ipm->n_total;
    double **vectors_total[] = {&ipm->x, &ipm->s, &ipm->x_best, &ipm->s_best,
                                &ipm->work_n[0], &ipm->work_n[1]};
    double **vectors_m[] = {&ipm->y, &ipm->y_best, &ipm->rb, &ipm->work_m[0],
                            &ipm->work_m[1]};
    arc_arena_take_vectors(ar, vectors_total, 6, total);
    arc_arena_take_vectors(ar, vectors_m, 5, m);
    ipm->rc = arc_arena_take(ar, n, sizeof(double));
    ipm->d = arc_arena_take(ar, n, sizeof(double));
    ipm->rw = arc_arena_take(ar, (size_t)ipm->n_bounded, sizeof(double));
}

int arc_ipm_init(arc_ipm *ipm, const arc_lp *lp, arc_normal_solver *solver)
{
    memset(ipm, 0, sizeof(*ipm));
    ipm->a = &lp->a;
    ipm->b = lp->b;
    ipm->c = lp->c;
    ipm->upper = lp->upper;
    ipm->m = lp->a.m;
    ipm->n = lp->a.n;
    ipm->n_free = lp->n_free;
    ipm->n_bounded = lp->n_bounded;
    ipm->n_total = ipm->n + ipm->n_bounded;
    ipm->n_pairs = ipm->n_total - ipm->n_free;
    ipm->solver = solver;
    return arc_arena_build(&ipm->arena, layout, ipm);
}

void arc_ipm_free(arc_ipm *ipm)
{
    arc_arena_release(&ipm->arena);
}

static void bound_residuals(arc_ipm *ipm)
{
    const int first = arc_ipm_bounded(ipm);
    for (int k = 0; k < ipm->n_bounded; k++)
    {
        const int j = first + k;
        ipm->rw[k] = (ipm->x[j] - ipm->upper[k]) + ipm->x[ipm->n + k];
        ipm->rc[j] -= ipm->s[ipm->n + k];
    }
}

/* One pass over the columns of A yields A x, A'y and x's together. */
void arc_ipm_residuals(arc_ipm *ipm)
{
    const arc_csc *a = ipm->a;
    double *rb = ipm->rb;
    double gap = 0.0;
    arc_zero(ipm->m, rb);
    for (int j = 0; j < ipm->n; j++)
    {
        const double xj = ipm->x[j];
        double aty = 0.0;
        for (int k = a->p[j]; k < a->p[j + 1]; k++)
        {
            rb[a->i[k]] += a->x[k] * xj;
            aty += a->x[k] * ipm->y[a->i[k]];
        }
        ipm->rc[j] = aty + (ipm->s[j] - ipm->c[j]);
        gap += xj * ipm->s[j];
    }
    arc_subtract(ipm->m, ipm->b, rb);
    bound_residuals(ipm);
    gap += arc_dot(ipm->n_bounded, ipm->x + ipm->n, ipm->s + ipm->n);
    ipm->mu = ipm->n_pairs > 0 ? gap / ipm->n_pairs : 0.0;
}

void arc_ipm_remember(arc_ipm *ipm)
{
    arc_copy(ipm->n_total, ipm->x, ipm->x_best);
    arc_copy(ipm->n_total, ipm->s, ipm->s_best);
    arc_copy(ipm->m, ipm->y, ipm->y_best);
}

void arc_ipm_restore(arc_ipm *ipm)
{
    arc_copy(ipm->n_total, ipm->x_best, ipm->x);
    arc_copy(ipm->n_total, ipm->s_best, ipm->s);
    arc_copy(ipm->m, ipm->y_best, ipm->y);
    arc_ipm_residuals(ipm);
}
