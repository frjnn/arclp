/*
 * Accelerate's worker threads only pay off for large blocks; for small ones
 * the synchronisation costs more than it saves.  The mode is process-wide,
 * so it is only changed by the first of possibly several concurrent
 * begin/end brackets (the owner), which restores the caller's setting at
 * its end; solves running alongside it leave the mode alone, which affects
 * speed only.  Other CBLAS libraries are not touched.
 */
#include "blas_threading.h"
#include "dense.h"

#if defined(ARCLP_USE_ACCELERATE)
#define ACCELERATE_NEW_LAPACK
#include <Accelerate/Accelerate.h>
#include <pthread.h>

static const double PARALLEL_WORK = 1e7;

static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
static int sessions;                     /* open brackets, all threads   */
static enum BLAS_THREADING saved_mode;   /* caller's mode, owner only    */
static _Thread_local int owner;          /* this thread may set the mode */
static _Thread_local int multi = -1;     /* mode it set; -1: none        */

void arc_dense_begin(void)
{
    pthread_mutex_lock(&lock);
    owner = sessions++ == 0;
    if (owner) {
        if (__builtin_available(macOS 15.0, *)) {
            saved_mode = BLASGetThreading();
        }
    }
    pthread_mutex_unlock(&lock);
    multi = -1;
}

void arc_dense_end(void)
{
    pthread_mutex_lock(&lock);
    if (owner && multi != -1) {
        if (__builtin_available(macOS 15.0, *)) {
            BLASSetThreading(saved_mode);
        }
    }
    sessions--;
    pthread_mutex_unlock(&lock);
    owner = 0;
    multi = -1;
}

void arc_blas_threads_for(double work)
{
    const int want = work >= PARALLEL_WORK;
    if (!owner || want == multi) {
        return;
    }
    if (__builtin_available(macOS 15.0, *)) {
        BLASSetThreading(want ? BLAS_THREADING_MULTI_THREADED
                              : BLAS_THREADING_SINGLE_THREADED);
    }
    multi = want;
}

#else

void arc_dense_begin(void) {}
void arc_dense_end(void) {}
void arc_blas_threads_for(double work) { (void)work; }

#endif
