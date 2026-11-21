/* Minimal assertion helpers shared by the test programs. */
#ifndef ARC_CHECK_H
#define ARC_CHECK_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int check_failures = 0;

#define CHECK(cond)                                                \
    do                                                             \
    {                                                              \
        if (!(cond))                                               \
        {                                                          \
            fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__, \
                    __LINE__, #cond);                              \
            check_failures++;                                      \
        }                                                          \
    } while (0)

#define CHECK_NEAR(a, b, tol)                                      \
    do                                                             \
    {                                                              \
        const double check_a_ = (a), check_b_ = (b);               \
        if (!(fabs(check_a_ - check_b_) <= (tol)))                 \
        {                                                          \
            fprintf(stderr, "%s:%d: %s = %.12g, expected %.12g\n", \
                    __FILE__, __LINE__, #a, check_a_, check_b_);   \
            check_failures++;                                      \
        }                                                          \
    } while (0)

static inline int check_report(const char *suite)
{
    printf("%-24s %s\n", suite, check_failures == 0 ? "ok" : "FAILED");
    return check_failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}

/* Deterministic generator so failures are reproducible. */
static inline double check_random(unsigned *state)
{
    *state = *state * 1664525u + 1013904223u;
    return (double)(*state >> 8) / 16777216.0;
}

#endif
