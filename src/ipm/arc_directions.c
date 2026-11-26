/*
 * The three solves of an arc-search iteration share one factorization of
 * M = A D A' (see newton.h):
 *     dot   first derivative, from the residuals
 *     p     centering part of the second derivative, target mu e
 *     q     second-order part, target -2 (dot.x o dot.s)
 * The first two are independent, so they are solved together.
 */
#include "arc_search.h"

#include "newton.h"
#include "util/vec.h"

static void fill_targets(arc_search_rule *r, const arc_ipm *ipm)
{
    for (int k = 0; k < ipm->n_total; k++)
    {
        r->target[k] = ipm->mu;
    }
}

void arc_search_directions(arc_search_rule *r, arc_ipm *ipm)
{
    double *g = r->work;
    double *w = r->work_p;
    arc_newton_first_rhs(ipm, g);
    fill_targets(r, ipm);
    arc_newton_second_rhs(ipm, r->target, w);
    arc_csc_mul2(ipm->a, g, w, r->dot.y, r->p.y);
    arc_subtract(ipm->m, ipm->b, r->dot.y);
    arc_negate(ipm->m, r->p.y);
    arc_normal_solve(ipm->solver, r->dot.y, 2); /* p.y follows dot.y */
    arc_csc_mul_t2(ipm->a, r->dot.y, r->p.y, r->dot.s, r->p.s);
    arc_newton_first_recover(ipm, r->dot.x, r->dot.s);
    arc_newton_second_recover(ipm, r->target, w, r->p.x, r->p.s);
    for (int k = 0; k < ipm->n_total; k++)
    {
        r->target[k] = -2.0 * r->dot.x[k] * r->dot.s[k];
    }
    arc_newton_second_rhs(ipm, r->target, w);
    arc_csc_mul(ipm->a, w, r->q.y);
    arc_negate(ipm->m, r->q.y);
    arc_normal_solve(ipm->solver, r->q.y, 1);
    arc_csc_mul_t(ipm->a, r->q.y, r->q.s);
    arc_newton_second_recover(ipm, r->target, w, r->q.x, r->q.s);
}
