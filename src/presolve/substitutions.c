#include "substitutions.h"

#include <stdlib.h>
#include <string.h>

static int reserve_records(arc_substitutions *s)
{
    if (s->count < s->capacity)
    {
        return 0;
    }
    const int capacity = s->capacity > 0 ? 2 * s->capacity : 64;
    int *pivot = realloc(s->pivot, sizeof(int) * (size_t)capacity);
    if (pivot != NULL)
    {
        s->pivot = pivot;
    }
    double *rhs = realloc(s->rhs, sizeof(double) * (size_t)capacity);
    if (rhs != NULL)
    {
        s->rhs = rhs;
    }
    int *start = realloc(s->start, sizeof(int) * (size_t)(capacity + 1));
    if (start != NULL)
    {
        s->start = start;
    }
    if (pivot == NULL || rhs == NULL || start == NULL)
    {
        return -1;
    }
    if (s->capacity == 0)
    {
        s->start[0] = 0;
    }
    s->capacity = capacity;
    return 0;
}

static int reserve_terms(arc_substitutions *s, int extra)
{
    if (s->terms + extra <= s->term_capacity)
    {
        return 0;
    }
    int capacity = s->term_capacity > 0 ? 2 * s->term_capacity : 256;
    while (capacity < s->terms + extra)
    {
        capacity *= 2;
    }
    int *cols = realloc(s->term_col, sizeof(int) * (size_t)capacity);
    if (cols != NULL)
    {
        s->term_col = cols;
    }
    double *coefs = realloc(s->term_coef, sizeof(double) * (size_t)capacity);
    if (coefs != NULL)
    {
        s->term_coef = coefs;
    }
    if (cols == NULL || coefs == NULL)
    {
        return -1;
    }
    s->term_capacity = capacity;
    return 0;
}

int arc_substitutions_add(arc_substitutions *s, int pivot, double rhs,
                          const int *cols, const double *coefs, int count)
{
    if (reserve_records(s) != 0 || reserve_terms(s, count) != 0)
    {
        return -1;
    }
    /* A constant substitution has no terms and may have no term storage. */
    if (count > 0)
    {
        memcpy(s->term_col + s->terms, cols, sizeof(int) * (size_t)count);
        memcpy(s->term_coef + s->terms, coefs, sizeof(double) * (size_t)count);
        s->terms += count;
    }
    s->pivot[s->count] = pivot;
    s->rhs[s->count] = rhs;
    s->start[++s->count] = s->terms;
    return 0;
}

void arc_substitutions_apply(const arc_substitutions *s, double *x)
{
    for (int r = s->count - 1; r >= 0; r--)
    {
        double value = s->rhs[r];
        for (int t = s->start[r]; t < s->start[r + 1]; t++)
        {
            value -= s->term_coef[t] * x[s->term_col[t]];
        }
        x[s->pivot[r]] = value;
    }
}

void arc_substitutions_free(arc_substitutions *s)
{
    free(s->pivot);
    free(s->rhs);
    free(s->start);
    free(s->term_col);
    free(s->term_coef);
    memset(s, 0, sizeof(*s));
}
