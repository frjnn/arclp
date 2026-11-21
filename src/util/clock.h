/* Monotonic wall clock in seconds. */
#ifndef ARC_CLOCK_H
#define ARC_CLOCK_H

#include <time.h>

static inline double arc_seconds(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + 1e-9 * (double)t.tv_nsec;
}

#endif
