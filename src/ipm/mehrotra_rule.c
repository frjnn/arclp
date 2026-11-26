/*
 * Mehrotra's predictor-corrector rule, ported from mehrotra.m (Yang's
 * reference implementation) so that both methods share the presolve, the
 * starting point and the stopping rules of the comparison in the paper.
 */
#include "mehrotra_steps.h"

#include <math.h>
#include <stdlib.h>

static const double STALL_STEP = 1e-8;

static void layout(arc_arena *ar, mehrotra_rule *r)
{
    const size_t m = (size_t)r->m;
    const size_t n = (size_t)r->n;
    r->dy = arc_arena_take(ar, m, sizeof(double));
    r->ddy = arc_arena_take(ar, m, sizeof(double));
    double **vectors[] = {&r->dx, &r->ds, &r->ddx, &r->dds, &r->work};
    for (size_t k = 0; k < sizeof(vectors) / sizeof(*vectors); k++)
    {
        *vectors[k] = arc_arena_take(ar, n, sizeof(double));
    }
}

/* v := v + step (dd - d) */
static void move(int len, double *v, const double *d, const double *dd,
                 double step)
{
    for (int i = 0; i < len; i++)
    {
        v[i] += -d[i] * step + dd[i] * step;
    }
}

static int mehrotra_step(arc_step_rule *self, arc_ipm *ipm)
{
    mehrotra_rule *r = (mehrotra_rule *)self;
    mehrotra_affine(r, ipm);
    mehrotra_corrector(r, ipm);
    const double shrink = exp(-(double)(r->iteration + 2));
    r->step_x = mehrotra_step_length(ipm->n, ipm->x, r->dx, r->ddx, shrink);
    r->step_s = mehrotra_step_length(ipm->n, ipm->s, r->ds, r->dds, shrink);
    move(ipm->n, ipm->x, r->dx, r->ddx, r->step_x);
    move(ipm->n, ipm->s, r->ds, r->dds, r->step_s);
    move(ipm->m, ipm->y, r->dy, r->ddy, r->step_s);
    r->iteration++;
    return 0;
}

static int mehrotra_stalled(const arc_step_rule *self, const arc_ipm *ipm)
{
    const mehrotra_rule *r = (const mehrotra_rule *)self;
    (void)ipm;
    return r->step_x < STALL_STEP && r->step_s < STALL_STEP;
}

static void mehrotra_destroy(arc_step_rule *self)
{
    mehrotra_rule *r = (mehrotra_rule *)self;
    arc_arena_release(&r->arena);
    free(r);
}

static const arc_step_rule_ops MEHROTRA_OPS = {mehrotra_step, mehrotra_stalled,
                                               mehrotra_destroy};

arc_step_rule *arc_mehrotra_rule_create(int m, int n)
{
    mehrotra_rule *r = calloc(1, sizeof(*r));
    if (r == NULL)
    {
        return NULL;
    }
    r->base.ops = &MEHROTRA_OPS;
    r->m = m;
    r->n = n;
    arc_arena_measure(&r->arena);
    layout(&r->arena, r);
    if (arc_arena_commit(&r->arena) != 0)
    {
        free(r);
        return NULL;
    }
    layout(&r->arena, r);
    return &r->base;
}
