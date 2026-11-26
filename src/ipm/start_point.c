/*
 * Starting point: candidates of Mehrotra (SIAM J. Optim. 1992) and of
 * Lustig, Marsten and Shanno (ORSA J. Comput. 1992); the one with the
 * smaller max(||rb||, ||rc||, mu) wins.
 */
#include "ipm.h"

#include "util/vec.h"

static double merit(arc_ipm *ipm)
{
    arc_ipm_residuals(ipm);
    const double rb = hypot(arc_norm2(ipm->m, ipm->rb),
                            arc_norm2(ipm->n_bounded, ipm->rw));
    const double rc = arc_norm2(ipm->n, ipm->rc);
    return arc_fmax(arc_fmax(rb, rc), ipm->mu);
}

/* v := (AA')^{-1} v,  out := A'v */
static void solve_and_lift(arc_ipm *ipm, double *v, double *out)
{
    arc_normal_solve(ipm->solver, v, 1);
    arc_csc_mul_t(ipm->a, v, out);
}

/* Least-squares solutions: tx = A'(AA')^{-1} b, ty = (AA')^{-1} A c,
 * ts = c - A'ty; for bounded columns also w = u - tx and the split of ts
 * into s - z. */
static int least_squares(arc_ipm *ipm, double *tx, double *ts, double *ty,
                         double *u)
{
    for (int j = 0; j < ipm->n; j++)
    {
        ipm->d[j] = 1.0;
    }
    if (arc_normal_factor(ipm->solver, ipm->d) != 0)
    {
        return -1;
    }
    arc_csc_mul(ipm->a, ipm->c, ty);
    solve_and_lift(ipm, ty, ts);
    for (int j = 0; j < ipm->n; j++)
    {
        ts[j] = ipm->c[j] - ts[j];
    }
    arc_copy(ipm->m, ipm->b, u);
    solve_and_lift(ipm, u, tx);
    arc_zero(ipm->n_free, ts);
    for (int j = arc_ipm_bounded(ipm); j < ipm->n; j++)
    {
        const int k = arc_ipm_slack(ipm, j);
        tx[k] = ipm->upper[k - ipm->n] - tx[j];
        ts[k] = arc_fmax(-ts[j], 0.0);
        ts[j] = arc_fmax(ts[j], 0.0);
    }
    return 0;
}

/* Free variables keep their least-squares values. */
static void place_free(arc_ipm *ipm, const double *tx)
{
    for (int k = 0; k < ipm->n_free; k++)
    {
        ipm->x[k] = tx[k];
        ipm->s[k] = 0.0;
    }
}

/* A cost that is constant on the feasible set has rounding errors for its
 * least-squares s, and shifts just as small, which would start s on the
 * boundary; a zero b does the same to x.  This share of the data is below
 * every shift on Netlib. */
static const double LEAST_SHIFT = 1e-4;

static double at_least(double shift, double size)
{
    return arc_fmax(shift, LEAST_SHIFT * (1.0 + size));
}

/* Half the products x_j s_j, shared out over the other side. */
static double balance(double product, double other_sum)
{
    return other_sum > 0.0 ? 0.5 * product / other_sum : 0.0;
}

/* Mehrotra's shift: both candidates become positive on the pairs, then are
 * moved so that x and s are balanced against the products x_j s_j; spread
 * multiplies the shifts. */
static void mehrotra_candidate(arc_ipm *ipm, double spread,
                               const double *tx,
                               const double *ts, const double *ty)
{
    const int f = ipm->n_free;
    const double dx = arc_fmax(-1.5 * arc_min(ipm->n_pairs, tx + f), 0.0);
    const double ds = arc_fmax(-1.5 * arc_min(ipm->n_pairs, ts + f), 0.0);
    double sum_x = 0.0;
    double sum_s = 0.0;
    double product = 0.0;
    for (int k = f; k < ipm->n_total; k++)
    {
        sum_x += tx[k] + dx;
        sum_s += ts[k] + ds;
        product += (tx[k] + dx) * (ts[k] + ds);
    }
    const double shift_x = at_least(spread * (dx + balance(product, sum_s)),
                                    arc_max_abs(ipm->m, ipm->b));
    const double shift_s = at_least(spread * (ds + balance(product, sum_x)),
                                    arc_max_abs(ipm->n, ipm->c));
    place_free(ipm, tx);
    for (int k = f; k < ipm->n_total; k++)
    {
        ipm->x[k] = tx[k] + shift_x;
        ipm->s[k] = ts[k] + shift_s;
    }
    arc_copy(ipm->m, ty, ipm->y);
}

static void lustig_candidate(arc_ipm *ipm, const double *tx)
{
    const int f = ipm->n_free;
    const double floor_b = arc_fmax(arc_norm1(ipm->m, ipm->b) / 100.0, 100.0);
    const double floor_x =
        arc_fmax(-arc_min(ipm->n_pairs, tx + f), floor_b);
    const double floor_s = 1.0 + arc_norm1(ipm->n, ipm->c);
    place_free(ipm, tx);
    for (int k = f; k < ipm->n_total; k++)
    {
        const double c = k < ipm->n ? ipm->c[k] : 0.0;
        ipm->x[k] = arc_fmax(floor_x, tx[k]);
        if (c >= 0.0)
        {
            ipm->s[k] = c + floor_s;
        }
        else
        {
            ipm->s[k] = c < -floor_s ? -c : floor_s;
        }
    }
    arc_zero(ipm->m, ipm->y);
}

int arc_ipm_start_point(arc_ipm *ipm, double spread)
{
    double *tx = ipm->work_n[0];
    double *ts = ipm->work_n[1];
    double *ty = ipm->work_m[0];
    if (least_squares(ipm, tx, ts, ty, ipm->work_m[1]) != 0)
    {
        return -1;
    }
    mehrotra_candidate(ipm, spread, tx, ts, ty);
    const double mehrotra = merit(ipm);
    arc_ipm_remember(ipm);
    lustig_candidate(ipm, tx);
    const double lustig = merit(ipm);
    if (mehrotra < lustig)
    {
        arc_ipm_restore(ipm);
    }
    return 0;
}
