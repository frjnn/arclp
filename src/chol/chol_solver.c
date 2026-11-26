#include "chol_solver.h"

#include "numeric.h"
#include "refine.h"
#include "symbolic.h"

#include <stdlib.h>

typedef struct
{
    arc_normal_solver base;   /* must stay first */
    arc_symbolic sym;
    arc_numeric num;
    arc_refine refine;
} chol_solver;

static int chol_factor(arc_normal_solver *self, const double *d)
{
    chol_solver *solver = (chol_solver *)self;
    solver->refine.d = d;
    return arc_numeric_factor(&solver->num, d);
}

static void chol_solve(arc_normal_solver *self, double *r, int nrhs)
{
    chol_solver *solver = (chol_solver *)self;
    for (int done = 0; done < nrhs; done += ARC_MAX_RHS)
    {
        const int chunk = nrhs - done < ARC_MAX_RHS ? nrhs - done : ARC_MAX_RHS;
        arc_refine_solve(&solver->refine, &solver->num,
                         r + (long)done * solver->sym.m, chunk);
    }
}

static void chol_refine(arc_normal_solver *self, int enabled)
{
    ((chol_solver *)self)->refine.enabled = enabled;
}

static void chol_stats(const arc_normal_solver *self, arc_normal_stats *out)
{
    const chol_solver *solver = (const chol_solver *)self;
    out->factor_nonzeros = solver->sym.val_ptr[solver->sym.nsuper];
    out->tiny_pivots = solver->num.tiny_pivots;
}

static void chol_destroy(arc_normal_solver *self)
{
    chol_solver *solver = (chol_solver *)self;
    arc_refine_free(&solver->refine);
    arc_numeric_free(&solver->num);
    arc_symbolic_free(&solver->sym);
    free(solver);
}

static const arc_normal_solver_ops CHOL_OPS = {
    chol_factor, chol_solve, chol_refine, chol_stats, chol_destroy};

arc_normal_solver *arc_chol_solver_create(const arc_csc *a)
{
    chol_solver *solver = calloc(1, sizeof(*solver));
    if (solver == NULL)
    {
        return NULL;
    }
    solver->base.ops = &CHOL_OPS;
    const arc_relax_policy policy = arc_default_relax_policy();
    if (arc_symbolic_analyse(a, &policy, &solver->sym) != 0)
    {
        free(solver);
        return NULL;
    }
    if (arc_numeric_create(&solver->sym, &solver->num) != 0
        || arc_refine_create(&solver->refine, a) != 0)
        {
        chol_destroy(&solver->base);
        return NULL;
    }
    return &solver->base;
}
