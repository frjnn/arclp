#include "substitutions.h"

#include "util/grow.h"

#include <string.h>

static int reserve_records(arc_substitutions *s)
{
    if (s->count < s->capacity)
    {
        return 0;
    }
    const int capacity = arc_grown_capacity(s->capacity, s->count + 1, 64);
    const size_t n = (size_t)capacity;
    if (arc_resize(&s->pivot, n, sizeof(int)) != 0
        || arc_resize(&s->rhs, n, sizeof(double)) != 0
        || arc_resize(&s->start, n + 1, sizeof(int)) != 0)
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
    const int needed = s->terms + extra;
    if (needed <= s->term_capacity)
    {
        return 0;
    }
    const int capacity = arc_grown_capacity(s->term_capacity, needed, 256);
    if (arc_resize(&s->term_col, (size_t)capacity, sizeof(int)) != 0
        || arc_resize(&s->term_coef, (size_t)capacity,
                      sizeof(double)) != 0)
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
