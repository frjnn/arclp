#include "chol_solver.h"

#include "numeric.h"
#include "symbolic.h"

#include <stdlib.h>

typedef struct {
    arc_normal_solver base;   /* must stay first */
    arc_symbolic sym;
    arc_numeric num;
} chol_solver;

static int chol_factor(arc_normal_solver *self, const double *d)
{
    return arc_numeric_factor(&((chol_solver *)self)->num, d);
}

static void chol_solve(arc_normal_solver *self, double *r, int nrhs)
{
    chol_solver *solver = (chol_solver *)self;
    for (int done = 0; done < nrhs; done += ARC_MAX_RHS) {
        const int chunk = nrhs - done < ARC_MAX_RHS ? nrhs - done : ARC_MAX_RHS;
        arc_numeric_solve(&solver->num, r + (long)done * solver->sym.m, chunk);
    }
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
    arc_numeric_free(&solver->num);
    arc_symbolic_free(&solver->sym);
    free(solver);
}

static const arc_normal_solver_ops CHOL_OPS = {chol_factor, chol_solve,
                                               chol_stats, chol_destroy};

arc_normal_solver *arc_chol_solver_create(const arc_csc *a)
{
    chol_solver *solver = calloc(1, sizeof(*solver));
    if (solver == NULL) {
        return NULL;
    }
    solver->base.ops = &CHOL_OPS;
    const arc_relax_policy policy = arc_default_relax_policy();
    if (arc_symbolic_analyse(a, &policy, &solver->sym) != 0) {
        free(solver);
        return NULL;
    }
    if (arc_numeric_create(&solver->sym, &solver->num) != 0) {
        arc_symbolic_free(&solver->sym);
        free(solver);
        return NULL;
    }
    return &solver->base;
}
