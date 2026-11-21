/*
 * Table 1 of Y. Yang, "Two computationally efficient polynomial-iteration
 * infeasible interior-point algorithms for linear programming", Numerical
 * Algorithms 79 (2018) 957-992: the 51 Netlib problems in standard form.
 * Iteration counts are those reported for arclp2.m (Algorithm 5.1, the
 * method of arcLP.m) and mehrotra.m.  The reference objective was computed
 * independently with HiGHS to 13 significant digits.
 */
#ifndef ARC_PAPER_TABLE_H
#define ARC_PAPER_TABLE_H

typedef struct {
    const char *name;
    int paper_arc_iterations;
    int paper_mehrotra_iterations;
    double reference_objective;
} paper_row;

extern const paper_row PAPER_TABLE[];
extern const int PAPER_TABLE_SIZE;

#endif
