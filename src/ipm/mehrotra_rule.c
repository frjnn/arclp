/*
 * Mehrotra's predictor-corrector rule, ported from mehrotra.m (Yang's
 * reference implementation) so that both methods share the presolve, the
 * starting point and the stopping rules of the comparison in the paper.
 */
#include "mehrotra_steps.h"

#include <math.h>
#include <stddef.h>

static const double STALL_STEP = 1e-8;

static void layout(arc_arena *ar, void *owner)
{
    mehrotra_rule *r = owner;
    double **multipliers[] = {&r->dy, &r->ddy};
    double **vectors[] = {&r->dx, &r->ds, &r->ddx, &r->dds, &r->work,
                          &r->target};
    arc_arena_take_vectors(ar, multipliers, 2, (size_t)r->base.m);
    arc_arena_take_vectors(ar, vectors, sizeof(vectors) / sizeof(*vectors),
                           (size_t)r->base.n);
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

static void mehrotra_step(arc_step_rule *self, arc_ipm *ipm)
{
    mehrotra_rule *r = (mehrotra_rule *)self;
    mehrotra_affine(r, ipm);
    mehrotra_corrector(r, ipm);
    const double shrink = exp(-(double)(r->iteration + 2));
    const int f = ipm->n_free;
    r->step_x = mehrotra_step_length(ipm->n_pairs, ipm->x + f, r->dx + f,
                                     r->ddx + f, shrink);
    r->step_s = mehrotra_step_length(ipm->n_pairs, ipm->s + f, r->ds + f,
                                     r->dds + f, shrink);
    move(ipm->n_total, ipm->x, r->dx, r->ddx, r->step_x);
    move(ipm->n_total, ipm->s, r->ds, r->dds, r->step_s);
    move(ipm->m, ipm->y, r->dy, r->ddy, r->step_s);
    r->iteration++;
}

static int mehrotra_stalled(const arc_step_rule *self, const arc_ipm *ipm)
{
    const mehrotra_rule *r = (const mehrotra_rule *)self;
    (void)ipm;
    return r->step_x < STALL_STEP && r->step_s < STALL_STEP;
}

const arc_step_rule_ops ARC_MEHROTRA_OPS = {
    mehrotra_step, mehrotra_stalled, layout, NULL, 1.0,
    sizeof(mehrotra_rule)};
