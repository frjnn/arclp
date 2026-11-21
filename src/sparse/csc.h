/* Read-only view of a compressed sparse column matrix and the products
 * the solver needs.  Whoever builds a matrix keeps its own writable
 * storage; everything downstream sees this const view. */
#ifndef ARC_CSC_H
#define ARC_CSC_H

typedef struct
{
    int m;
    int n;
    const int *p;    /* n + 1 column pointers */
    const int *i;    /* row indices           */
    const double *x; /* values                */
} arc_csc;

/* y := A x */
void arc_csc_mul(const arc_csc *a, const double *x, double *y);

/* x := A' y */
void arc_csc_mul_t(const arc_csc *a, const double *y, double *x);

/* y1 := A x1 and y2 := A x2 in one pass over A. */
void arc_csc_mul2(const arc_csc *a, const double *x1, const double *x2,
                  double *y1, double *y2);

/* x1 := A' y1 and x2 := A' y2 in one pass over A. */
void arc_csc_mul_t2(const arc_csc *a, const double *y1, const double *y2,
                    double *x1, double *x2);

/* Smallest and largest absolute nonzero; returns 0 for an all-zero matrix. */
int arc_csc_abs_range(const arc_csc *a, double *lowest, double *highest);

#endif
