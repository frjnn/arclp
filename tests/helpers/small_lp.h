/* min -x1 - 2 x2  s.t.  x1 + x2 + s1 = 4,  x1 + 3 x2 + s2 = 6,  x, s >= 0;
 * the optimum is x = (3, 1) with objective -5. */
#ifndef ARC_SMALL_LP_H
#define ARC_SMALL_LP_H

static const int SMALL_P[] = {0, 2, 4, 5, 6};
static const int SMALL_I[] = {0, 1, 0, 1, 0, 1};
static const double SMALL_V[] = {1, 1, 1, 3, 1, 1};
static const double SMALL_B[] = {4, 6};
static const double SMALL_C[] = {-1, -2, 0, 0};

enum { SMALL_M = 2, SMALL_N = 4 };

#endif
