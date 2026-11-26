#include "arc_search.h"

#include "util/vec.h"

#include <math.h>

static void layout(arc_arena *ar, void *owner)
{
    arc_search_rule *r = owner;
    const size_t m = (size_t)r->base.m;
    double *ys = arc_arena_take(ar, 2 * m, sizeof(double));
    r->dot.y = ys;
    r->p.y = ys != NULL ? ys + m : NULL; /* adjacent: solved together */
    r->q.y = arc_arena_take(ar, m, sizeof(double));
    double **vectors[] = {&r->dot.x, &r->dot.s, &r->p.x, &r->p.s, &r->q.x,
                          &r->q.s, &r->work, &r->work_p, &r->target};
    arc_arena_take_vectors(ar, vectors, sizeof(vectors) / sizeof(*vectors),
                           (size_t)r->base.n);
}

static void init(arc_step_rule *self)
{
    ((arc_search_rule *)self)->nu = 1.0;
}

/* v := v - dot sin(alpha) + curvature (sigma p + q)(1 - cos(alpha)) */
static void move(int len, double *v, const double *dot, const double *p,
                 const double *q, double alpha, const arc_angle *a)
{
    const double sa = sin(alpha);
    const double ca = a->curvature * (1.0 - cos(alpha));
    for (int i = 0; i < len; i++)
    {
        v[i] += -dot[i] * sa + (a->sigma * p[i] + q[i]) * ca;
    }
}

static void arc_step(arc_step_rule *self, arc_ipm *ipm)
{
    arc_search_rule *r = (arc_search_rule *)self;
    arc_search_directions(r, ipm);
    const arc_angle a = arc_search_angle(r, ipm);
    move(ipm->m, ipm->y, r->dot.y, r->p.y, r->q.y, a.alpha_s, &a);
    move(ipm->n_total, ipm->s, r->dot.s, r->p.s, r->q.s, a.alpha_s, &a);
    move(ipm->n_total, ipm->x, r->dot.x, r->p.x, r->q.x, a.alpha_x, &a);
    r->nu *= 1.0 - sin(arc_fmin(a.alpha_x, a.alpha_s));
    r->alpha = arc_fmax(a.alpha_x, a.alpha_s);
}

static int arc_stalled(const arc_step_rule *self, const arc_ipm *ipm)
{
    const arc_search_rule *r = (const arc_search_rule *)self;
    (void)ipm;
    return r->alpha < ARC_STALL_ANGLE;
}

/* A start further from the boundary pays: 4 took the fewest iterations on
 * Netlib of the values tried, 1 to 8. */
enum
{
    START_SPREAD = 4
};

const arc_step_rule_ops ARC_SEARCH_OPS = {
    arc_step, arc_stalled, layout, init, START_SPREAD,
    sizeof(arc_search_rule)};
