/*
 * Equilibration of min c'x, A x = b: rows and columns
 * are divided by factors that bring the largest entry of each close to 1
 * (Ruiz iterations), and the solution is mapped back by arc_unscale.
 *     A^ = R A C,   b^ = R b,   c^ = C c,   u^ = u / C,   x = C x^
 */
#ifndef ARC_SCALING_H
#define ARC_SCALING_H

#include "sparse/lp.h"
#include "util/arena.h"

typedef struct
{
    double *row, *col; /* the diagonals of R and C                        */
    double *row_unit;  /* 1 / row: turns scaled residuals into real ones  */
    double *col_unit;
    double *row_max, *col_max; /* scratch                                */
    int m, n;
    arc_arena arena;
} arc_scaling;

/* Scales the problem in place; returns -1 when memory runs out. */
int arc_scale_problem(arc_lp *lp, arc_scaling *s);

/* Maps a solution of the scaled problem back; y may be NULL. */
void arc_unscale(const arc_scaling *s, const arc_lp *lp, double *x,
                 double *y);

void arc_scaling_free(arc_scaling *s);

#endif
