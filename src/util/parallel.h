/*
 * Optional OpenMP loop parallelism.
 *
 *     ARC_FOR(i, count, work, statement);
 *
 * runs `statement` for i = 0 .. count-1, split across threads when the
 * build has OpenMP (ARCLP_USE_OPENMP) and the loop does at least
 * ARC_PARALLEL_WORK units of work.  Smaller loops run the plain serial
 * loop without entering the OpenMP runtime at all: even a region that
 * decides not to fork costs a runtime call and an outlined loop, which
 * measurably slows the many small loops of a factorization.  Only used for
 * iterations that write disjoint data, so results do not depend on the
 * number of threads.
 */
#ifndef ARC_PARALLEL_H
#define ARC_PARALLEL_H

/* Work is counted in units of roughly one memory-bound multiply-add; call
 * sites weight their iterations accordingly. */
enum { ARC_PARALLEL_WORK = 262144 };

/* For a loop body shared by the serial path and an OpenMP region: with two
 * call sites the compiler stops inlining it, and the serial path, which
 * every small problem takes, measurably slows down. */
#if defined(__GNUC__)
#define ARC_ALWAYS_INLINE inline __attribute__((always_inline))
#else
#define ARC_ALWAYS_INLINE inline
#endif

#if defined(ARCLP_USE_OPENMP)
#include <omp.h>
#define ARC_PRAGMA(text) _Pragma(#text)
#define ARC_FOR(i, count, work, statement)                                  \
    do {                                                                    \
        const int arc_for_count_ = (count);                                 \
        if ((work) >= ARC_PARALLEL_WORK) {                                  \
            ARC_PRAGMA(omp parallel for schedule(static))                   \
            for (int i = 0; i < arc_for_count_; i++) {                      \
                statement;                                                  \
            }                                                               \
        } else {                                                            \
            for (int i = 0; i < arc_for_count_; i++) {                      \
                statement;                                                  \
            }                                                               \
        }                                                                   \
    } while (0)
#else
#define ARC_FOR(i, count, work, statement)                                  \
    do {                                                                    \
        const int arc_for_count_ = (count);                                 \
        (void)(work);                                                       \
        for (int i = 0; i < arc_for_count_; i++) {                          \
            statement;                                                      \
        }                                                                   \
    } while (0)
#endif

#endif
