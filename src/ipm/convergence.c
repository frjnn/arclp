#include "convergence.h"

#include "util/vec.h"

static const double CERTIFICATE_TOLERANCE = 1e-8;

/* ||u o v||, in the original units when a scaling was applied. */
static double norm_scaled(int len, const double *v, const double *unit)
{
    if (unit == NULL)
    {
        return arc_norm2(len, v);
    }
    double sum = 0.0;
    for (int i = 0; i < len; i++)
    {
        const double real = v[i] * unit[i];
        sum += real * real;
    }
    return sqrt(sum);
}

static double bound_norm(const arc_ipm *ipm, const double *v)
{
    return norm_scaled(ipm->n_bounded, v, ipm->bound_unit);
}

arc_ipm_objectives arc_ipm_objective(const arc_ipm *ipm)
{
    const arc_ipm_objectives o = {
        arc_dot(ipm->n, ipm->c, ipm->x),
        arc_dot(ipm->m, ipm->b, ipm->y)
            - arc_dot(ipm->n_bounded, ipm->upper, ipm->s + ipm->n)};
    return o;
}

void arc_ipm_set_units(arc_ipm *ipm, const double *row_unit,
                       const double *col_unit, const double *col)
{
    ipm->row_unit = row_unit;
    ipm->col_unit = col_unit;
    ipm->bound_unit = col != NULL ? col + arc_ipm_bounded(ipm) : NULL;
    ipm->norm_b = hypot(norm_scaled(ipm->m, ipm->b, ipm->row_unit),
                        bound_norm(ipm, ipm->upper));
    ipm->norm_c = norm_scaled(ipm->n, ipm->c, ipm->col_unit);
}

/* The gap is relative to the objective itself: dividing it by the norms of
 * b and c, as arcLP.m does, hides errors of any size when b is large. */
arc_ipm_errors arc_ipm_measure(const arc_ipm *ipm)
{
    const arc_ipm_objectives o = arc_ipm_objective(ipm);
    const double rb = norm_scaled(ipm->m, ipm->rb, ipm->row_unit);
    const arc_ipm_errors e = {
        hypot(rb, bound_norm(ipm, ipm->rw)) / (1.0 + ipm->norm_b),
        norm_scaled(ipm->n, ipm->rc, ipm->col_unit) / (1.0 + ipm->norm_c),
        fabs(o.primal - o.dual)
            / (1.0 + arc_fmin(fabs(o.primal), fabs(o.dual)))};
    return e;
}

double arc_ipm_complementarity(const arc_ipm *ipm)
{
    const arc_ipm_objectives o = arc_ipm_objective(ipm);
    return ipm->mu * ipm->n_pairs
         / (1.0 + arc_fmin(fabs(o.primal), fabs(o.dual)));
}

static double norm_of_sum(int len, const double *u, const double *v)
{
    double sum = 0.0;
    for (int i = 0; i < len; i++)
    {
        sum += (u[i] + v[i]) * (u[i] + v[i]);
    }
    return sqrt(sum);
}

/* Farkas: b'y - u'z > 0 with A'y + s - z = 0 proves the primal infeasible;
 * c'x < 0 with A x = 0 and x zero on the bounded columns proves it
 * unbounded.  A'y + s - z = rc + c and A x = rb + b need no product. */
arc_ipm_outcome arc_ipm_certified(const arc_ipm *ipm)
{
    const arc_ipm_objectives o = arc_ipm_objective(ipm);
    if (o.dual > 0.0 && norm_of_sum(ipm->n, ipm->rc, ipm->c)
                            <= CERTIFICATE_TOLERANCE * o.dual)
    {
        return ARC_IPM_PRIMAL_INFEASIBLE;
    }
    const double limit = -CERTIFICATE_TOLERANCE * o.primal;
    if (o.primal < 0.0 && norm_of_sum(ipm->m, ipm->rb, ipm->b) <= limit
        && arc_norm2(ipm->n_bounded,
                     ipm->x + arc_ipm_bounded(ipm)) <= limit)
    {
        return ARC_IPM_UNBOUNDED;
    }
    return ARC_IPM_STALLED;
}
