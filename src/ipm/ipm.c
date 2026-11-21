#include "ipm.h"

#include "util/vec.h"

#include <string.h>

static void layout(arc_arena *ar, arc_ipm *ipm)
{
    const size_t m = (size_t)ipm->m;
    const size_t n = (size_t)ipm->n;
    double **vectors_n[] = {&ipm->x, &ipm->s, &ipm->rc, &ipm->d,
                            &ipm->x_old, &ipm->s_old, &ipm->work_n[0],
                            &ipm->work_n[1]};
    double **vectors_m[] = {&ipm->y, &ipm->y_old, &ipm->rb, &ipm->work_m[0],
                            &ipm->work_m[1]};
    for (size_t k = 0; k < sizeof(vectors_n) / sizeof(*vectors_n); k++) {
        *vectors_n[k] = arc_arena_take(ar, n, sizeof(double));
    }
    for (size_t k = 0; k < sizeof(vectors_m) / sizeof(*vectors_m); k++) {
        *vectors_m[k] = arc_arena_take(ar, m, sizeof(double));
    }
}

int arc_ipm_init(arc_ipm *ipm, const arc_csc *a, const double *b,
                 const double *c, arc_normal_solver *solver)
{
    memset(ipm, 0, sizeof(*ipm));
    ipm->a = a;
    ipm->b = b;
    ipm->c = c;
    ipm->m = a->m;
    ipm->n = a->n;
    ipm->solver = solver;
    ipm->norm_b = arc_norm2(ipm->m, b);
    ipm->norm_c = arc_norm2(ipm->n, c);
    arc_arena_measure(&ipm->arena);
    layout(&ipm->arena, ipm);
    if (arc_arena_commit(&ipm->arena) != 0) {
        return -1;
    }
    layout(&ipm->arena, ipm);
    return 0;
}

void arc_ipm_free(arc_ipm *ipm)
{
    arc_arena_release(&ipm->arena);
}

/* One pass over the columns of A yields A x, A'y and x's together. */
void arc_ipm_residuals(arc_ipm *ipm)
{
    const arc_csc *a = ipm->a;
    double *rb = ipm->rb;
    double gap = 0.0;
    for (int i = 0; i < ipm->m; i++) {
        rb[i] = 0.0;
    }
    for (int j = 0; j < ipm->n; j++) {
        const double xj = ipm->x[j];
        double aty = 0.0;
        for (int k = a->p[j]; k < a->p[j + 1]; k++) {
            rb[a->i[k]] += a->x[k] * xj;
            aty += a->x[k] * ipm->y[a->i[k]];
        }
        ipm->rc[j] = aty + (ipm->s[j] - ipm->c[j]);
        gap += xj * ipm->s[j];
    }
    arc_subtract(ipm->m, ipm->b, rb);
    ipm->mu = gap / ipm->n;
}
