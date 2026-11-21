"""Renders the CSV written by `arclp_bench --csv` as Markdown tables."""
import csv
import sys

SOLVED = 1e-6


def load(path):
    rows = {}
    for r in csv.DictReader(open(path)):
        rows.setdefault(r["problem"], {})[r["method"]] = r
    return rows


def solved(r):
    return r["status"] == "optimal" and r["objective_error"] != "nan" \
        and float(r["objective_error"]) < SOLVED


def cell(r):
    mark = "" if solved(r) else " †"
    return f"{r['iterations']} ({r['paper_iterations']}){mark}"


def main(path):
    rows = load(path)
    print("| problem | rows×cols (presolved) | arcLP it (paper) | Mehrotra it (paper) "
          "| arcLP rel. obj. error | arcLP ‖Ax−b‖ | arcLP time [s] | Mehrotra time [s] |")
    print("|---|---|---|---|---|---|---|---|")
    for name in rows:
        a, m = rows[name]["arc"], rows[name]["mehrotra"]
        print(f"| {name} | {a['rows']}×{a['cols']} | {cell(a)} | {cell(m)} "
              f"| {float(a['objective_error']):.1e} | {float(a['infeasibility']):.1e} "
              f"| {float(a['t_total']):.4f} | {float(m['t_total']):.4f} |")
    print()
    for method, label in (("arc", "arcLP"), ("mehrotra", "Mehrotra")):
        rs = [rows[n][method] for n in rows]
        it = sum(int(r["iterations"]) for r in rs)
        paper = sum(int(r["paper_iterations"]) for r in rs)
        t = sum(float(r["t_total"]) for r in rs)
        print(f"- **{label}**: {sum(map(solved, rs))}/{len(rs)} solved to relative "
              f"objective error < {SOLVED:g}; {it} iterations in total (paper: {paper}); "
              f"{t:.2f} s in total.")


if __name__ == "__main__":
    main(sys.argv[1])
