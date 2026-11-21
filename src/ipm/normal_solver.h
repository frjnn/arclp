/*
 * Abstract solver for the normal equations (A D A') y = r that every
 * interior-point iteration needs.  The iteration logic depends only on this
 * interface; concrete factorizations plug in behind it.
 */
#ifndef ARC_NORMAL_SOLVER_H
#define ARC_NORMAL_SOLVER_H

typedef struct arc_normal_solver arc_normal_solver;

typedef struct
{
    long factor_nonzeros;
    int tiny_pivots; /* pivots replaced during the last factorization */
} arc_normal_stats;

typedef struct
{
    /* Factorizes A diag(d) A'; returns 0 on success. */
    int (*factor)(arc_normal_solver *self, const double *d);
    /* Overwrites the m-by-nrhs column-major block r with (A D A')^{-1} r. */
    void (*solve)(arc_normal_solver *self, double *r, int nrhs);
    void (*stats)(const arc_normal_solver *self, arc_normal_stats *out);
    void (*destroy)(arc_normal_solver *self);
} arc_normal_solver_ops;

struct arc_normal_solver
{
    const arc_normal_solver_ops *ops;
};

static inline int arc_normal_factor(arc_normal_solver *s, const double *d)
{
    return s->ops->factor(s, d);
}

static inline void arc_normal_solve(arc_normal_solver *s, double *r, int nrhs)
{
    s->ops->solve(s, r, nrhs);
}

static inline void arc_normal_destroy(arc_normal_solver *s)
{
    if (s != 0)
    {
        s->ops->destroy(s);
    }
}

#endif
