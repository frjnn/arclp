#include "arc_search.h"

#include <math.h>
#include <stdlib.h>

static const double STALL_ANGLE = 1e-8;
static const double STALL_MU = 1e-8;

static void layout(arc_arena *ar, arc_search_rule *r)
{
    const size_t m = (size_t)r->m;
    const size_t n = (size_t)r->n;
    double *ys = arc_arena_take(ar, 2 * m, sizeof(double));
    r->dot.y = ys;
    r->p.y = ys != NULL ? ys + m : NULL; /* adjacent: solved together */
    r->q.y = arc_arena_take(ar, m, sizeof(double));
    double **vectors[] = {&r->dot.x, &r->dot.s, &r->p.x, &r->p.s, &r->q.x,
                          &r->q.s, &r->work};
    for (size_t k = 0; k < sizeof(vectors) / sizeof(*vectors); k++)
    {
        *vectors[k] = arc_arena_take(ar, n, sizeof(double));
    }
}

/* v := v - dot sin(alpha) + (sigma p + q)(1 - cos(alpha)) */
static void move(int len, double *v, const double *dot, const double *p,
                 const double *q, arc_angle a)
{
    const double sa = sin(a.alpha);
    const double ca = 1.0 - cos(a.alpha);
    for (int i = 0; i < len; i++)
    {
        v[i] += -dot[i] * sa + (a.sigma * p[i] + q[i]) * ca;
    }
}

static int arc_step(arc_step_rule *self, arc_ipm *ipm)
{
    arc_search_rule *r = (arc_search_rule *)self;
    arc_search_directions(r, ipm);
    const arc_angle a = arc_search_angle(r, ipm);
    move(ipm->m, ipm->y, r->dot.y, r->p.y, r->q.y, a);
    move(ipm->n, ipm->s, r->dot.s, r->p.s, r->q.s, a);
    move(ipm->n, ipm->x, r->dot.x, r->p.x, r->q.x, a);
    r->nu *= 1.0 - sin(a.alpha);
    r->alpha = a.alpha;
    r->sigma = a.sigma;
    return 0;
}

static int arc_stalled(const arc_step_rule *self, const arc_ipm *ipm)
{
    const arc_search_rule *r = (const arc_search_rule *)self;
    return r->alpha < STALL_ANGLE || ipm->mu < STALL_MU;
}

static void arc_destroy(arc_step_rule *self)
{
    arc_search_rule *r = (arc_search_rule *)self;
    arc_arena_release(&r->arena);
    free(r);
}

static const arc_step_rule_ops ARC_OPS = {arc_step, arc_stalled, arc_destroy};

arc_step_rule *arc_search_rule_create(int m, int n)
{
    arc_search_rule *r = calloc(1, sizeof(*r));
    if (r == NULL)
    {
        return NULL;
    }
    r->base.ops = &ARC_OPS;
    r->m = m;
    r->n = n;
    r->nu = 1.0;
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
