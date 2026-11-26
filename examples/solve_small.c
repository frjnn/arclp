/*
 * Minimal use of the public API:
 *     minimize  -x1 - 2 x2
 *     s.t.       x1 +   x2 + s1      = 4
 *                x1 + 3 x2      + s2 = 6,   x, s >= 0
 * Optimum: x = (3, 1), objective -5.
 */
#include <arclp/arclp.h>

#include <stdio.h>

int main(void)
{
    /* A in compressed sparse column form, columns x1, x2, s1, s2. */
    const int col_ptr[] = {0, 2, 4, 5, 6};
    const int row_idx[] = {0, 1, 0, 1, 0, 1};
    const double values[] = {1, 1, 1, 3, 1, 1};
    const double b[] = {4, 6};
    const double c[] = {-1, -2, 0, 0};
    const arclp_problem problem = {2, 4, col_ptr, row_idx, values, b, c, NULL, NULL};

    arclp_options options;
    arclp_default_options(&options);
    double x[4];
    arclp_info info;
    const arclp_status status = arclp_solve(&problem, &options, x, &info);

    printf("status     %s\n", arclp_status_name(status));
    printf("iterations %d\n", info.iterations);
    printf("objective  %.10f\n", info.objective);
    printf("x          %.6f %.6f\n", x[0], x[1]);
    return status == ARCLP_OPTIMAL ? 0 : 1;
}
