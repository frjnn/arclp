/* Threading policy of the CBLAS backend (see docs/EXPERIMENTS.md, E16). */
#ifndef ARC_BLAS_THREADING_H
#define ARC_BLAS_THREADING_H

/* Selects the library's threading for a call doing `work` multiply-adds:
 * multi-threaded only for large blocks.  Effective between
 * arc_dense_begin and arc_dense_end; a no-op where unsupported. */
void arc_blas_threads_for(double work);

#endif
