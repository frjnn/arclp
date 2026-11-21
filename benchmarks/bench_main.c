/*
 * Reproduces the Netlib comparison of Yang (2018), Table 1:
 *     arclp_bench DATA_DIR [--problem NAME] [--method arc|mehrotra|both]
 *                          [--repeat N] [--csv FILE]
 * DATA_DIR holds the lp_<name>.mat files of github.com/yaguangyang/arcLP.
 */
#include "bench.h"

#include <stdlib.h>
#include <string.h>

static int parse(int argc, char **argv, bench_config *c)
{
    const bench_config defaults = {NULL, NULL, NULL, 1, 1, 1};
    *c = defaults;
    for (int k = 1; k < argc; k++) {
        const char *arg = argv[k];
        const char *value = k + 1 < argc ? argv[k + 1] : NULL;
        if (strcmp(arg, "--problem") == 0 && value) {
            c->only = argv[++k];
        } else if (strcmp(arg, "--csv") == 0 && value) {
            c->csv_path = argv[++k];
        } else if (strcmp(arg, "--repeat") == 0 && value) {
            c->repeat = atoi(argv[++k]);
        } else if (strcmp(arg, "--method") == 0 && value) {
            c->run_arc = strcmp(value, "mehrotra") != 0;
            c->run_mehrotra = strcmp(value, "arc") != 0;
            k++;
        } else if (c->data_dir == NULL && arg[0] != '-') {
            c->data_dir = arg;
        } else {
            return -1;
        }
    }
    return c->data_dir != NULL && c->repeat > 0 ? 0 : -1;
}

typedef struct {
    bench_totals arc, mehrotra;
    FILE *csv;
} bench_state;

static void run_problem(const bench_config *c, const paper_row *row,
                        bench_state *st)
{
    char path[4096], err[256];
    snprintf(path, sizeof(path), "%s/lp_%s.mat", c->data_dir, row->name);
    mat_lp lp;
    if (mat_lp_load(path, &lp, err, sizeof(err)) != 0) {
        fprintf(stderr, "skipping %s: %s\n", row->name, err);
        return;
    }
    bench_result arc = {0}, mehrotra = {0};
    if (c->run_arc) {
        bench_solve(&lp, ARCLP_METHOD_ARC_SEARCH, c->repeat,
                    row->reference_objective, &arc);
        bench_accumulate(&st->arc, &arc, row->paper_arc_iterations);
    }
    if (c->run_mehrotra) {
        bench_solve(&lp, ARCLP_METHOD_MEHROTRA, c->repeat,
                    row->reference_objective, &mehrotra);
        bench_accumulate(&st->mehrotra, &mehrotra,
                         row->paper_mehrotra_iterations);
    }
    bench_print_row(stdout, row, &lp, c, &arc, &mehrotra);
    if (st->csv != NULL && c->run_arc) {
        bench_csv_row(st->csv, row, "arc", row->paper_arc_iterations, &arc);
    }
    if (st->csv != NULL && c->run_mehrotra) {
        bench_csv_row(st->csv, row, "mehrotra",
                      row->paper_mehrotra_iterations, &mehrotra);
    }
    mat_lp_free(&lp);
}

int main(int argc, char **argv)
{
    bench_config config;
    if (parse(argc, argv, &config) != 0) {
        fprintf(stderr, "usage: %s DATA_DIR [--problem NAME] "
                        "[--method arc|mehrotra|both] [--repeat N] "
                        "[--csv FILE]\n", argv[0]);
        return EXIT_FAILURE;
    }
    bench_state st = {{0}, {0}, NULL};
    if (config.csv_path != NULL && (st.csv = fopen(config.csv_path, "w"))) {
        fprintf(st.csv, "problem,method,status,iterations,paper_iterations,"
                        "objective,reference,objective_error,infeasibility,"
                        "rows,cols,factor_nnz,tiny_pivots,t_presolve,"
                        "t_analyse,t_iterate,t_total\n");
    }
    bench_print_header(stdout, &config);
    for (int k = 0; k < PAPER_TABLE_SIZE; k++) {
        if (config.only == NULL || strcmp(config.only, PAPER_TABLE[k].name) == 0) {
            run_problem(&config, &PAPER_TABLE[k], &st);
        }
    }
    if (config.run_arc) {
        bench_print_totals(stdout, "arcLP", &st.arc);
    }
    if (config.run_mehrotra) {
        bench_print_totals(stdout, "Mehrotra", &st.mehrotra);
    }
    if (st.csv != NULL) {
        fclose(st.csv);
    }
    return EXIT_SUCCESS;
}
