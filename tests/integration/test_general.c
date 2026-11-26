/*
 * Random problems with a known optimum: x, y, s, z are chosen
 * complementary first and b, c follow from them, so the optimal objective
 * is exact whatever mix of bounds, free variables, many rows (solved as the
 * dual) or dense columns (stretched) the problem has.
 */
#include "check.h"

#include <arclp/arclp.h>

typedef struct
{
    int rows, cols;        /* structural columns; `slacks` rows get one more */
    int slacks;
    double density;
    int dense;             /* leading columns with every entry               */
    int flat;              /* zero reduced costs: every feasible x is optimal */
} shape;

typedef struct
{
    int m, n;
    int *p, *i;
    double *v, *b, *c, *lower, *upper, *x;
    double optimum;
} instance;

static unsigned seed = 12345;

static double uniform(void)
{
    return check_random(&seed);
}

static void *take(size_t count, size_t size)
{
    void *p = calloc(count + 1, size);
    if (p == NULL)
    {
        exit(EXIT_FAILURE);
    }
    return p;
}

/* A variable at a bound or inside, with the multiplier that allows it;
 * returns the reduced cost s - z. */
static double draw_variable(instance *lp, int j)
{
    const double kind = uniform();
    const double base = 4.0 * uniform() - 2.0;
    const int inside = uniform() < 0.5;
    const double gap = 1.0 + uniform();
    lp->lower[j] = -INFINITY;
    lp->upper[j] = INFINITY;
    lp->x[j] = base;
    if (kind < 0.2)
    {
        return 0.0;                                  /* free */
    }
    if (kind < 0.6)
    {
        lp->lower[j] = inside ? base - gap : base;   /* lower only */
        return inside ? 0.0 : gap;
    }
    if (kind < 0.7)
    {
        lp->upper[j] = inside ? base + gap : base;   /* upper only */
        return inside ? 0.0 : -gap;
    }
    lp->lower[j] = inside ? base - gap : base;       /* boxed */
    lp->upper[j] = base + gap;
    return inside ? 0.0 : gap;
}

static instance generate(shape s)
{
    instance lp = {0};
    lp.m = s.rows;
    lp.n = s.cols + s.slacks;
    lp.p = take((size_t)lp.n + 1, sizeof(int));
    lp.i = take((size_t)lp.m * (size_t)lp.n, sizeof(int));
    lp.v = take((size_t)lp.m * (size_t)lp.n, sizeof(double));
    double **vectors[] = {&lp.b, &lp.c, &lp.lower, &lp.upper, &lp.x};
    for (int k = 0; k < 5; k++)
    {
        *vectors[k] = take((size_t)lp.m + (size_t)lp.n, sizeof(double));
    }
    double *y = take((size_t)lp.m, sizeof(double));
    for (int i = 0; i < lp.m; i++)
    {
        y[i] = 2.0 * uniform() - 1.0;
    }
    int nnz = 0;
    for (int j = 0; j < lp.n; j++)
    {
        const int slack = j >= s.cols;
        double reduced = slack ? 0.0 : draw_variable(&lp, j);
        if (slack)
        {
            const int tight = uniform() < 0.5;
            lp.upper[j] = INFINITY;
            lp.x[j] = tight ? 0.0 : 1.0 + uniform();
            y[j - s.cols] = tight ? -1.0 - uniform() : 0.0;
            reduced = tight ? -y[j - s.cols] : 0.0;
        }
        lp.p[j] = nnz;
        for (int i = 0; i < lp.m; i++)
        {
            const int take_entry = slack ? i == j - s.cols
                : j < s.dense || uniform() < s.density || i == j % lp.m;
            if (take_entry)
            {
                lp.i[nnz] = i;
                lp.v[nnz++] = slack ? 1.0 : 2.0 * uniform() - 1.0;
            }
        }
        lp.p[j + 1] = nnz;
        lp.c[j] = s.flat ? 0.0 : reduced;
    }
    for (int j = 0; j < lp.n; j++)
    {
        for (int k = lp.p[j]; k < lp.p[j + 1]; k++)
        {
            lp.b[lp.i[k]] += lp.v[k] * lp.x[j];
            lp.c[j] += lp.v[k] * y[lp.i[k]];
        }
        lp.optimum += lp.c[j] * lp.x[j];
    }
    free(y);
    return lp;
}

static void release(instance *lp)
{
    free(lp->p);
    free(lp->i);
    free(lp->v);
    free(lp->b);
    free(lp->c);
    free(lp->lower);
    free(lp->upper);
    free(lp->x);
}

/* Largest violation of A x = b and of the bounds. */
static double violation(const instance *lp, const double *x)
{
    double *ax = take((size_t)lp->m, sizeof(double));
    double worst = 0.0;
    for (int j = 0; j < lp->n; j++)
    {
        worst = fmax(worst, fmax(lp->lower[j] - x[j], x[j] - lp->upper[j]));
        for (int k = lp->p[j]; k < lp->p[j + 1]; k++)
        {
            ax[lp->i[k]] += lp->v[k] * x[j];
        }
    }
    for (int i = 0; i < lp->m; i++)
    {
        worst = fmax(worst, fabs(ax[i] - lp->b[i]) / (1.0 + fabs(lp->b[i])));
    }
    free(ax);
    return worst;
}

static void solve(shape s)
{
    instance lp = generate(s);
    const arclp_problem problem = {lp.m, lp.n, lp.p, lp.i, lp.v, lp.b,
                                   lp.c, lp.lower, lp.upper};
    double *x = take((size_t)lp.n, sizeof(double));
    for (int method = 0; method < 2; method++)
    {
        for (int presolve = 0; presolve < 2; presolve++)
        {
            arclp_options o;
            arclp_default_options(&o);
            o.method = method == 0 ? ARCLP_METHOD_ARC_SEARCH
                                   : ARCLP_METHOD_MEHROTRA;
            o.presolve = presolve;
            arclp_info info;
            CHECK(arclp_solve(&problem, &o, x, &info) == ARCLP_OPTIMAL);
            CHECK_NEAR(info.objective, lp.optimum,
                       1e-6 * (1.0 + fabs(lp.optimum)));
            CHECK(violation(&lp, x) < 1e-5);
        }
    }
    free(x);
    release(&lp);
}

int main(void)
{
    for (int repeat = 0; repeat < 100; repeat++)
    {
        const shape wide = {20, 50, 10, 0.2, 0, 0};
        const shape tall = {120, 20, 120, 0.2, 0, 0};   /* solved as the dual */
        const shape square = {40, 40, 20, 0.3, 0, 0};
        solve(wide);
        solve(tall);
        solve(square);
    }
    const shape dense = {2100, 2400, 0, 0.0005, 4, 0};  /* stretched columns */
    solve(dense);
    solve(dense);
    for (int repeat = 0; repeat < 10; repeat++)
    {
        const shape flat = {20, 60, 0, 0.2, 0, 1};
        solve(flat);
    }
    return check_report("general form");
}
