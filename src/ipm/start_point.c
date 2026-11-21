/*
 * Starting point of arcLP.m: candidates of Mehrotra (SIAM J. Optim. 1992)
 * and Lustig, Marsten and Shanno (ORSA J. Comput. 1992); the one with the
 * smaller max(||rb||, ||rc||, mu) wins.
 */
#include "ipm.h"

#include "util/vec.h"

static double merit(arc_ipm *ipm)
{
    arc_ipm_residuals(ipm);
    const double rb = arc_norm2(ipm->m, ipm->rb);
    const double rc = arc_norm2(ipm->n, ipm->rc);
    return arc_fmax(arc_fmax(rb, rc), ipm->mu);
}

/* Least-squares solutions: tx = A'(AA')^{-1} b, ty = (AA')^{-1} A c,
 * ts = c - A'ty. */
static int least_squares(arc_ipm *ipm, double *tx, double *ts, double *ty,
                         double *u)
{
    for (int j = 0; j < ipm->n; j++) {
        ipm->d[j] = 1.0;
    }
    if (arc_normal_factor(ipm->solver, ipm->d) != 0) {
        return -1;
    }
    arc_csc_mul(ipm->a, ipm->c, ty);
    arc_normal_solve(ipm->solver, ty, 1);
    arc_csc_mul_t(ipm->a, ty, ts);
    for (int j = 0; j < ipm->n; j++) {
        ts[j] = ipm->c[j] - ts[j];
    }
    arc_copy(ipm->m, ipm->b, u);
    arc_normal_solve(ipm->solver, u, 1);
    arc_csc_mul_t(ipm->a, u, tx);
    return 0;
}

/* Mehrotra's shift.  As in arcLP.m the normalising sums use only the last
 * component (the MATLAB loop assigns instead of accumulating). */
static void mehrotra_candidate(arc_ipm *ipm, const double *tx,
                               const double *ts, const double *ty)
{
    const int n = ipm->n;
    const double dx = arc_fmax(-1.1 * arc_min(n, tx), 0.0);
    const double ds = arc_fmax(-1.1 * arc_min(n, ts), 0.0);
    const double sum_x = tx[n - 1] + dx;
    const double sum_s = ts[n - 1] + ds;
    double product = 0.0;
    for (int j = 0; j < n; j++) {
        product += (tx[j] + dx) * (ts[j] + ds);
    }
    const double shift_x = dx + 0.5 * product / sum_x;
    const double shift_s = ds + 0.5 * product / sum_s;
    for (int j = 0; j < n; j++) {
        ipm->x[j] = tx[j] + shift_x;
        ipm->s[j] = ts[j] + shift_s;
    }
    arc_copy(ipm->m, ty, ipm->y);
}

static void lustig_candidate(arc_ipm *ipm, const double *tx)
{
    const int n = ipm->n;
    const double floor_b = arc_fmax(arc_norm1(ipm->m, ipm->b) / 100.0, 100.0);
    const double floor_x = arc_fmax(-arc_min(n, tx), floor_b);
    const double floor_s = 1.0 + arc_norm1(n, ipm->c);
    for (int j = 0; j < n; j++) {
        const double c = ipm->c[j];
        ipm->x[j] = arc_fmax(floor_x, tx[j]);
        if (c >= 0.0) {
            ipm->s[j] = c + floor_s;
        } else {
            ipm->s[j] = c < -floor_s ? -c : floor_s;
        }
    }
    for (int i = 0; i < ipm->m; i++) {
        ipm->y[i] = 0.0;
    }
}

int arc_ipm_start_point(arc_ipm *ipm)
{
    double *tx = ipm->work_n[0];
    double *ts = ipm->work_n[1];
    double *ty = ipm->work_m[0];
    if (least_squares(ipm, tx, ts, ty, ipm->work_m[1]) != 0) {
        return -1;
    }
    mehrotra_candidate(ipm, tx, ts, ty);
    const double mehrotra = merit(ipm);
    arc_copy(ipm->n, ipm->x, ipm->x_old);
    arc_copy(ipm->n, ipm->s, ipm->s_old);
    lustig_candidate(ipm, tx);
    const double lustig = merit(ipm);
    if (mehrotra < lustig) {
        arc_copy(ipm->n, ipm->x_old, ipm->x);
        arc_copy(ipm->n, ipm->s_old, ipm->s);
        arc_copy(ipm->m, ty, ipm->y);
        arc_ipm_residuals(ipm);
    }
    return 0;
}
