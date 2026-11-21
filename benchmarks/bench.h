/* Benchmark harness for the paper's Netlib comparison. */
#ifndef ARC_BENCH_H
#define ARC_BENCH_H

#include <arclp/arclp.h>

#include <arclp/mat5.h>
#include "paper_table.h"

#include <stdio.h>

typedef struct {
    const char *data_dir;
    const char *only;        /* single problem name, or NULL           */
    const char *csv_path;    /* optional machine-readable output       */
    int repeat;              /* timed runs per problem; best is kept   */
    int run_arc, run_mehrotra;
} bench_config;

typedef struct {
    arclp_info info;
    double seconds;          /* best total time over the repeats       */
    double objective_error;  /* |obj - ref| / max(1, |ref|)            */
} bench_result;

typedef struct {
    int problems, solved;
    long iterations, paper_iterations;
    double seconds;
} bench_totals;

int bench_solve(const mat_lp *lp, arclp_method method, int repeat,
                double reference, bench_result *out);

void bench_print_header(FILE *out, const bench_config *config);
void bench_print_row(FILE *out, const paper_row *row, const mat_lp *lp,
                     const bench_config *config, const bench_result *arc,
                     const bench_result *mehrotra);
void bench_print_totals(FILE *out, const char *label, const bench_totals *t);
void bench_csv_row(FILE *csv, const paper_row *row, const char *method,
                   int paper_iterations, const bench_result *r);
void bench_accumulate(bench_totals *t, const bench_result *r,
                      int paper_iterations);

#endif
