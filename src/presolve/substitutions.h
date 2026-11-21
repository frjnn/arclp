/* Growable log of rule-9 substitutions, replayed backwards in postsolve. */
#ifndef ARC_SUBSTITUTIONS_H
#define ARC_SUBSTITUTIONS_H

#include "presolve.h"

int arc_substitutions_add(arc_substitutions *s, int pivot, double rhs,
                          const int *cols, const double *coefs, int count);

/* x[pivot] = rhs - sum coef * x[col], newest substitution first. */
void arc_substitutions_apply(const arc_substitutions *s, double *x);

void arc_substitutions_free(arc_substitutions *s);

#endif
