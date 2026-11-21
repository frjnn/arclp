/*
 * With D = X S^{-1} and M = A D A' the three solves of an arc-search
 * iteration share one factorization of M:
 *     M dot.y = A D rc - b       (first derivative)
 *     M p.y   = -A S^{-1} mu e   (centering part of the second derivative)
 *     M q.y   = 2 A S^{-1} (dot.x o dot.s)
 * The first two are independent, so they are solved together.
 */
#include "arc_search.h"

#include "util/vec.h"

static void negate(int n, double *v)
{
    for (int j = 0; j < n; j++)
    {
        v[j] = -v[j];
    }
}

static void first_derivative(arc_search_rule *r, const arc_ipm *ipm)
{
    for (int j = 0; j < ipm->n; j++)
    {
        r->dot.s[j] = ipm->rc[j] - r->dot.s[j];
        r->dot.x[j] = ipm->x[j] - ipm->d[j] * r->dot.s[j];
    }
}

static void centering_term(arc_search_rule *r, const arc_ipm *ipm)
{
    for (int j = 0; j < ipm->n; j++)
    {
        r->p.s[j] = -r->p.s[j];
        r->p.x[j] = ipm->mu / ipm->s[j] - ipm->d[j] * r->p.s[j];
    }
}

void arc_search_directions(arc_search_rule *r, arc_ipm *ipm)
{
    const int n = ipm->n;
    double *w1 = r->work;
    double *w2 = r->q.x; /* free until q is computed */
    for (int j = 0; j < n; j++)
    {
        w1[j] = ipm->d[j] * ipm->rc[j];
        w2[j] = ipm->mu / ipm->s[j];
    }
    arc_csc_mul2(ipm->a, w1, w2, r->dot.y, r->p.y);
    arc_subtract(ipm->m, ipm->b, r->dot.y);
    negate(ipm->m, r->p.y);
    arc_normal_solve(ipm->solver, r->dot.y, 2); /* p.y follows dot.y */
    arc_csc_mul_t2(ipm->a, r->dot.y, r->p.y, r->dot.s, r->p.s);
    first_derivative(r, ipm);
    centering_term(r, ipm);
    for (int j = 0; j < n; j++)
    {
        w1[j] = 2.0 * r->dot.x[j] * r->dot.s[j] / ipm->s[j];
    }
    arc_csc_mul(ipm->a, w1, r->q.y);
    arc_normal_solve(ipm->solver, r->q.y, 1);
    arc_csc_mul_t(ipm->a, r->q.y, r->q.s);
    for (int j = 0; j < n; j++)
    {
        r->q.s[j] = -r->q.s[j];
        r->q.x[j] = -ipm->d[j] * r->q.s[j] - w1[j];
    }
}
