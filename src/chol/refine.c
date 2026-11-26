#include "refine.h"

#include "util/vec.h"

#include <string.h>

static const double TARGET = 1e-10; /* relative residual worth no more work */
static const double PROGRESS = 0.5; /* a step that helps */
enum
{
    MAX_STEPS = 20,
    PATIENCE = 3 /* steps without progress before giving up */
};

static void layout(arc_arena *ar, void *owner)
{
    arc_refine *f = owner;
    const size_t m = (size_t)f->a->m;
    f->rhs = arc_arena_take(ar, m * ARC_MAX_RHS, sizeof(double));
    f->r = arc_arena_take(ar, m, sizeof(double));
    f->p = arc_arena_take(ar, m, sizeof(double));
    f->q = arc_arena_take(ar, m, sizeof(double));
    f->z = arc_arena_take(ar, m, sizeof(double));
    f->best = arc_arena_take(ar, m, sizeof(double));
    f->t = arc_arena_take(ar, (size_t)f->a->n, sizeof(double));
}

int arc_refine_create(arc_refine *f, const arc_csc *a)
{
    memset(f, 0, sizeof(*f));
    f->a = a;
    f->enabled = 1;
    return arc_arena_build(&f->arena, layout, f);
}

void arc_refine_free(arc_refine *f)
{
    arc_arena_release(&f->arena);
}

/* out := A D A' v */
static void apply(const arc_refine *f, const double *v, double *out)
{
    arc_csc_mul_t(f->a, v, f->t);
    for (int j = 0; j < f->a->n; j++)
    {
        f->t[j] *= f->d[j];
    }
    arc_csc_mul(f->a, f->t, out);
}

/* p := z + beta p */
static void next_direction(int n, double beta, const double *z, double *p)
{
    for (int i = 0; i < n; i++)
    {
        p[i] = z[i] + beta * p[i];
    }
}

static void precondition(const arc_refine *f, arc_numeric *num)
{
    arc_copy(f->a->m, f->r, f->z);
    arc_numeric_solve(num, f->z, 1);
}

/* r := rhs - A D A' x */
static double residual(const arc_refine *f, const double *rhs,
                       const double *x)
{
    apply(f, x, f->q);
    for (int i = 0; i < f->a->m; i++)
    {
        f->r[i] = rhs[i] - f->q[i];
    }
    return arc_norm2(f->a->m, f->r);
}

/* Conjugate gradients from x.  Their residuals need not fall at every step,
 * so the best iterate is kept: the result is never worse than the first
 * solve. */
static void improve(arc_refine *f, arc_numeric *num, const double *rhs,
                    double *x)
{
    const int m = f->a->m;
    const double limit = TARGET * arc_norm2(m, rhs);
    double norm = residual(f, rhs, x);
    double lowest = norm;
    double reference = norm; /* residual at the last progress */
    int waited = 0;
    double rho = 0.0;
    arc_copy(m, x, f->best);
    for (int step = 0;
         step < MAX_STEPS && lowest > limit && waited < PATIENCE; step++)
    {
        precondition(f, num);
        const double rz = arc_dot(m, f->r, f->z);
        next_direction(m, step > 0 ? rz / rho : 0.0, f->z, f->p);
        rho = rz;
        apply(f, f->p, f->q);
        const double curvature = arc_dot(m, f->p, f->q);
        if (!(curvature > 0.0) || !isfinite(rz))
        {
            break;
        }
        const double alpha = rz / curvature;
        arc_axpy(m, alpha, f->p, x);
        arc_axpy(m, -alpha, f->q, f->r);
        norm = arc_norm2(m, f->r);
        if (norm < lowest)
        {
            lowest = norm;
            arc_copy(m, x, f->best);
        }
        waited++;
        if (norm < PROGRESS * reference)
        {
            reference = norm;
            waited = 0;
        }
    }
    arc_copy(m, f->best, x);
}

void arc_refine_solve(arc_refine *f, arc_numeric *num, double *r, int nrhs)
{
    const int m = f->a->m;
    arc_copy(m * nrhs, r, f->rhs);
    arc_numeric_solve(num, r, nrhs);
    if (!f->enabled)
    {
        return;
    }
    for (int k = 0; k < nrhs; k++)
    {
        improve(f, num, f->rhs + (long)k * m, r + (long)k * m);
    }
}
