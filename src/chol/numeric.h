/*
 * Numeric supernodal Cholesky factorization L L' = P (A D A') P'.  All
 * storage is allocated once in arc_numeric_create; factorizations and
 * solves then run without touching the heap.
 */
#ifndef ARC_NUMERIC_H
#define ARC_NUMERIC_H

#include "symbolic.h"

enum { ARC_MAX_RHS = 4 };

typedef struct {
    const arc_symbolic *sym;
    double *lx;          /* dense blocks of every supernode              */
    int *map;            /* row label -> position in current supernode   */
    int *relative;       /* target positions of one descendant update     */
    int *link_head;      /* descendants waiting to update a supernode    */
    int *link_next;
    int *next_row;       /* first unapplied row of each descendant       */
    double *update;      /* descendant update block                      */
    double *x;           /* permuted right-hand sides, m x ARC_MAX_RHS   */
    double *gather;      /* max_rows x ARC_MAX_RHS                       */
    double *pivot_floor;     /* per column: tolerance * diagonal entry   */
    double pivot_tolerance;
    int tiny_pivots;
    arc_arena arena;
} arc_numeric;

int arc_numeric_create(const arc_symbolic *sym, arc_numeric *num);

/* Factorizes A diag(d) A'; returns 0 on success. */
int arc_numeric_factor(arc_numeric *num, const double *d);

/* Overwrites the m x nrhs block b with the solution, nrhs <= ARC_MAX_RHS. */
void arc_numeric_solve(arc_numeric *num, double *b, int nrhs);

void arc_numeric_free(arc_numeric *num);

/* Adds the columns of P A D A' P' owned by supernode s into its block. */
void arc_assemble_supernode(const arc_symbolic *sym, int s, const double *d,
                            const int *map, double *block);

#endif
