/* Supernodal Cholesky implementation of the normal-equations interface. */
#ifndef ARC_CHOL_SOLVER_H
#define ARC_CHOL_SOLVER_H

#include "ipm/normal_solver.h"
#include "sparse/csc.h"

/* Orders and analyses A A' once and allocates all numeric storage.  The
 * matrix must outlive the solver.  Returns NULL when memory runs out. */
arc_normal_solver *arc_chol_solver_create(const arc_csc *a);

#endif
