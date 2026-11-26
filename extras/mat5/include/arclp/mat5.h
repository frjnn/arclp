/*
 * Minimal MATLAB level-5 MAT-file reader for SuiteSparse/Netlib LP files:
 * loads Problem.A (sparse, CSC), Problem.b and Problem.aux.c.  Optional
 * component (target arclp::mat5); the solver itself does not need it.
 */
#ifndef ARCLP_MAT5_H
#define ARCLP_MAT5_H

#include <stddef.h>

typedef struct
{
    int m, n;
    int *col_ptr, *row_idx;
    double *values, *b, *c;
} mat_lp;

/* Returns 0 on success; otherwise writes a message to err. */
int mat_lp_load(const char *path, mat_lp *lp, char *err, size_t err_len);

/* Loads DIR/lp_NAME.mat; on failure reports that NAME is skipped. */
int mat_lp_load_named(const char *dir, const char *name, mat_lp *lp);

void mat_lp_free(mat_lp *lp);

#endif
