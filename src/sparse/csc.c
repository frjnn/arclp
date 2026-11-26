#include "csc.h"

#include "util/parallel.h"
#include "util/vec.h"

#include <math.h>
#include <stddef.h>

void arc_csc_index_rows(const arc_csc *a, const int *label, const int *skip,
                        int *ptr, int *cols, int *pos, int *next)
{
    for (int r = 0; r <= a->m; r++)
    {
        ptr[r] = 0;
    }
    for (int j = 0; j < a->n; j++)
    {
        for (int k = a->p[j]; k < a->p[j + 1]; k++)
        {
            const int r = label != NULL ? label[a->i[k]] : a->i[k];
            ptr[r + 1] += skip == NULL || skip[j] >= 0;
        }
    }
    for (int r = 0; r < a->m; r++)
    {
        ptr[r + 1] += ptr[r];
        next[r] = ptr[r];
    }
    for (int j = 0; j < a->n; j++)
    {
        if (skip != NULL && skip[j] < 0)
        {
            continue;
        }
        for (int k = a->p[j]; k < a->p[j + 1]; k++)
        {
            const int slot = next[label != NULL ? label[a->i[k]] : a->i[k]]++;
            cols[slot] = j;
            if (pos != NULL)
            {
                pos[slot] = k;
            }
        }
    }
}

void arc_csc_mul(const arc_csc *a, const double *x, double *y)
{
    arc_zero(a->m, y);
    for (int j = 0; j < a->n; j++)
    {
        const double xj = x[j];
        for (int k = a->p[j]; k < a->p[j + 1]; k++)
        {
            y[a->i[k]] += a->x[k] * xj;
        }
    }
}

static inline double column_dot(const arc_csc *a, int j, const double *y)
{
    double sum = 0.0;
    for (int k = a->p[j]; k < a->p[j + 1]; k++)
    {
        sum += a->x[k] * y[a->i[k]];
    }
    return sum;
}

void arc_csc_mul_t(const arc_csc *a, const double *y, double *x)
{
    ARC_FOR(j, a->n, (long)a->p[a->n], x[j] = column_dot(a, j, y));
}

void arc_csc_mul2(const arc_csc *a, const double *x1, const double *x2,
                  double *y1, double *y2)
{
    arc_zero(a->m, y1);
    arc_zero(a->m, y2);
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

static inline void column_dot2(const arc_csc *a, int j, const double *y1,
                               const double *y2, double *x1, double *x2)
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

void arc_csc_mul_t2(const arc_csc *a, const double *y1, const double *y2,
                    double *x1, double *x2)
{
    ARC_FOR(j, a->n, (long)a->p[a->n], column_dot2(a, j, y1, y2, x1, x2));
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
