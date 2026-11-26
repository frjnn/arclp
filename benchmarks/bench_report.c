#include "bench.h"

static void print_method_header(FILE *out, const char *name)
{
    fprintf(out, " | %-8s %4s %5s %8s %8s %9s", name, "it", "paper", "obj err",
            "infeas", "seconds");
}

void bench_print_header(FILE *out, const bench_config *config)
{
    fprintf(out, "%-9s %6s %7s %6s %7s", "problem", "m", "n", "m'", "n'");
    if (config->run_arc)
    {
        print_method_header(out, "arcLP");
    }
    if (config->run_mehrotra)
    {
        print_method_header(out, "Mehrotra");
    }
    fputc('\n', out);
}

static void print_method(FILE *out, const bench_result *r, int paper)
{
    const char *flag = r->info.status == ARCLP_OPTIMAL ? "" : "!";
    fprintf(out, " | %-8s %4d %5d %8.1e %8.1e %9.4f", flag,
            r->info.iterations, paper, r->objective_error,
            r->info.infeasibility, r->seconds);
}

void bench_print_row(FILE *out, const paper_row *row, const mat_lp *lp,
                     const bench_config *config, const bench_result *arc,
                     const bench_result *mehrotra)
{
    const bench_result *any = config->run_arc ? arc : mehrotra;
    fprintf(out, "%-9s %6d %7d %6d %7d", row->name, lp->m, lp->n,
            any->info.reduced_rows, any->info.reduced_cols);
    if (config->run_arc)
    {
        print_method(out, arc, row->paper_arc_iterations);
    }
    if (config->run_mehrotra)
    {
        print_method(out, mehrotra, row->paper_mehrotra_iterations);
    }
    fputc('\n', out);
    fflush(out);
}

void bench_print_totals(FILE *out, const char *label, const bench_totals *t)
{
    fprintf(out, "%-9s solved %d/%d   iterations %ld (paper %ld)   "
                 "time %.3f s\n",
            label, t->solved, t->problems, t->iterations, t->paper_iterations,
            t->seconds);
}

void bench_csv_row(FILE *csv, const paper_row *row, const char *method,
                   int paper_iterations, const bench_result *r)
{
    const arclp_info *i = &r->info;
    fprintf(csv, "%s,%s,%s,%d,%d,%.12e,%.12e,%.3e,%.3e,%d,%d,%ld,%d,"
                 "%.6f,%.6f,%.6f,%.6f\n",
            row->name, method, arclp_status_name(i->status), i->iterations,
            paper_iterations, i->objective, row->reference_objective,
            r->objective_error, i->infeasibility, i->reduced_rows,
            i->reduced_cols, i->factor_nonzeros, i->tiny_pivots,
            i->seconds_presolve, i->seconds_analyse, i->seconds_iterate,
            r->seconds);
}
