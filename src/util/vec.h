/* Dense vector kernels used by the interior-point iterations. */
#ifndef ARC_VEC_H
#define ARC_VEC_H

#include <math.h>

static inline double arc_dot(int n, const double *x, const double *y)
{
    double sum = 0.0;
    for (int i = 0; i < n; i++)
    {
        sum += x[i] * y[i];
    }
    return sum;
}

static inline double arc_norm2(int n, const double *x)
{
    return sqrt(arc_dot(n, x, x));
}

static inline double arc_norm1(int n, const double *x)
{
    double sum = 0.0;
    for (int i = 0; i < n; i++)
    {
        sum += fabs(x[i]);
    }
    return sum;
}

static inline double arc_min(int n, const double *x)
{
    double lowest = INFINITY;
    for (int i = 0; i < n; i++)
    {
        lowest = x[i] < lowest ? x[i] : lowest;
    }
    return lowest;
}

static inline double arc_max(int n, const double *x)
{
    double highest = -INFINITY;
    for (int i = 0; i < n; i++)
    {
        highest = x[i] > highest ? x[i] : highest;
    }
    return highest;
}

static inline void arc_copy(int n, const double *x, double *y)
{
    for (int i = 0; i < n; i++)
    {
        y[i] = x[i];
    }
}

/* y := y - x */
static inline void arc_subtract(int n, const double *x, double *y)
{
    for (int i = 0; i < n; i++)
    {
        y[i] -= x[i];
    }
}

static inline double arc_fmin(double a, double b) { return a < b ? a : b; }
static inline double arc_fmax(double a, double b) { return a > b ? a : b; }

#endif
