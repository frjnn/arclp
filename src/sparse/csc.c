#include "csc.h"

#include <math.h>

void arc_csc_mul(const arc_csc *a, const double *x, double *y)
{
    for (int r = 0; r < a->m; r++)
    {
        y[r] = 0.0;
    }
    for (int j = 0; j < a->n; j++)
    {
        const double xj = x[j];
        for (int k = a->p[j]; k < a->p[j + 1]; k++)
        {
            y[a->i[k]] += a->x[k] * xj;
        }
    }
}

void arc_csc_mul_t(const arc_csc *a, const double *y, double *x)
{
    for (int j = 0; j < a->n; j++)
    {
        double sum = 0.0;
        for (int k = a->p[j]; k < a->p[j + 1]; k++)
        {
            sum += a->x[k] * y[a->i[k]];
        }
        x[j] = sum;
    }
}

void arc_csc_mul2(const arc_csc *a, const double *x1, const double *x2,
                  double *y1, double *y2)
{
    for (int r = 0; r < a->m; r++)
    {
        y1[r] = 0.0;
        y2[r] = 0.0;
    }
    for (int j = 0; j < a->n; j++)
    {
        const double u = x1[j];
        const double v = x2[j];
        for (int k = a->p[j]; k < a->p[j + 1]; k++)
        {
            y1[a->i[k]] += a->x[k] * u;
            y2[a->i[k]] += a->x[k] * v;
        }
    }
}

void arc_csc_mul_t2(const arc_csc *a, const double *y1, const double *y2,
                    double *x1, double *x2)
{
    for (int j = 0; j < a->n; j++)
    {
        double u = 0.0;
        double v = 0.0;
        for (int k = a->p[j]; k < a->p[j + 1]; k++)
        {
            u += a->x[k] * y1[a->i[k]];
            v += a->x[k] * y2[a->i[k]];
        }
        x1[j] = u;
        x2[j] = v;
    }
}

int arc_csc_abs_range(const arc_csc *a, double *lowest, double *highest)
{
    int found = 0;
    *lowest = INFINITY;
    *highest = 0.0;
    for (int k = 0; k < a->p[a->n]; k++)
    {
        const double v = fabs(a->x[k]);
        if (v != 0.0)
        {
            found = 1;
            *lowest = v < *lowest ? v : *lowest;
            *highest = v > *highest ? v : *highest;
        }
    }
    return found;
}
