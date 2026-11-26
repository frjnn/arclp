/* Smallest step angles over the components of x and s for one sigma. */
#ifndef ARC_BOUNDS_H
#define ARC_BOUNDS_H

#include "arc_filter.h"
#include "arc_search.h"

typedef struct
{
    double angle;
    arc_trig trig;
} arc_bound;

typedef struct
{
    arc_bound minus; /* components whose second derivative falls with sigma */
    arc_bound plus;  /* components whose second derivative grows with sigma */
    arc_bound all;   /* every component                                    */
    arc_bound x;     /* components of x and w                              */
    arc_bound s;     /* components of s and z                              */
} arc_bounds;

arc_bounds arc_bounds_at(const arc_search_rule *r, const arc_ipm *ipm,
                         double sigma, double phi, double psi);

#endif
