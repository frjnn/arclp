# arclp

A C11 implementation of Yaguang Yang's **arcLP**, the arc-search infeasible
interior-point method of

> Y. Yang, *Two computationally efficient polynomial-iteration infeasible interior-point algorithms for linear programming*, Numerical Algorithms 79 (2018) 957–992

It solves

```
minimize c'x   subject to   A x = b,  x >= 0
```

with the reference's presolve, starting point, stopping rules and safeguards. Mehrotra's predictor–corrector is included as a second, swappable step rule so the paper's comparison can be reproduced on identical footing.

## Build and test

Requirements: CMake >= 3.21, a C11 compiler, `zlib` (only for the `.mat` reader used by tests and benchmarks), and optionally:

- a CBLAS (Accelerate on macOS, e.g. OpenBLAS on Linux) for the dense kernels;
- [METIS](https://github.com/KarypisLab/METIS), whose nested-dissection ordering replaces AMD when it gives a smaller factor (`brew install metis`, `apt install libmetis-dev`);
- OpenMP for parallel loops (built into GCC; on macOS `brew install libomp`).

Each optional piece is detected automatically and only changes speed: results are identical with or without OpenMP and for any thread count.

```sh
scripts/fetch_netlib.sh            # optional: Netlib data for tests/benchmark
cmake --preset release             # configure  (build/release)
cmake --build --preset release     # build
ctest --preset release             # unit + integration tests
```

| preset | what it builds |
|---|---|
| `release` | optimized, `-march=native`, system BLAS |
| `debug` | debug info, warnings as errors |
| `sanitize` | debug with AddressSanitizer + UBSan |
| `portable` | no BLAS (portable C kernels), no `-march=native` |
| `shared` | shared library exporting only the public API |

### Build options

| option | default | what it controls |
|---|---|---|
| `ARCLP_BLAS` | `auto` | Backend of the dense matrix kernels: `accelerate` (Apple), `cblas` (any CBLAS, e.g. OpenBLAS), or `none` (portable C, no dependency, slower). `auto` picks Accelerate on macOS, else a CBLAS if found, else `none`. |
| `ARCLP_METIS` | `auto` | `auto` uses METIS if found, `on` requires it, `off` never uses it. When used, analysis also tries nested dissection for factors above 10⁶ stored entries and keeps whichever ordering gives the smaller factor. |
| `ARCLP_OPENMP` | `auto` | `auto` uses OpenMP if found, `on` requires it, `off` builds single-threaded. Parallel loops run only above a work threshold; thread count via `OMP_NUM_THREADS`. |
| `ARCLP_NATIVE` | `ON` | Adds `-march=native` (optimizes for the build machine). Turn off for binaries that must run on other CPUs. |
| `ARCLP_WERROR` | `OFF` | Treats compiler warnings as errors (on in the `debug` preset). |
| `ARCLP_SANITIZE` | empty | Runtime checkers, comma-separated: e.g. `address,undefined` (memory errors, undefined behavior) or `thread` (data races). |
| `ARCLP_BUILD_TESTS` | `ON`* | Builds the test suite and the `.mat` reader it uses. |
| `ARCLP_BUILD_BENCHMARKS` | `ON`* | Builds `arclp_bench` and `netlib_verify`. |
| `ARCLP_BUILD_EXAMPLES` | `ON`* | Builds `examples/solve_small`. |
| `ARCLP_NETLIB_DIR` | `data/arcLP` | Where tests and the `benchmark` target look for the `lp_*.mat` files. |
| `BUILD_SHARED_LIBS` | `OFF` | Standard CMake switch: shared (`.dylib`/`.so`) instead of static library. |

\* `ON` when `arclp` is the top-level project, `OFF` when another project
includes it with `add_subdirectory`.

`cmake -LH <build-dir>` lists every option with its current value.

Install and use from another project:

```sh
cmake --install build/release --prefix /opt/arclp
```

```cmake
find_package(arclp 1.0 REQUIRED)            # or pkg-config: arclp
target_link_libraries(app PRIVATE arclp::arclp)
```

## Usage

```c
#include <arclp/arclp.h>

/* A in compressed sparse column form */
arclp_problem lp = {m, n, col_ptr, row_idx, values, b, c, NULL, NULL};
arclp_options opt;
arclp_default_options(&opt);          /* arc-search, tol 1e-8, presolve on */
arclp_info info;
if (arclp_solve(&lp, &opt, x, &info) == ARCLP_OPTIMAL) {
    printf("%d iterations, objective %.10g\n", info.iterations, info.objective);
}
```

`arclp_info` also reports the reduced problem size, factor size, replaced pivots and per-phase timings. A complete program is in [examples/solve_small.c](examples/solve_small.c).

## Benchmark

```sh
scripts/fetch_netlib.sh                          # .mat files into data/arcLP
cmake --build --preset release --target benchmark  # writes build/release/netlib_results.csv
python3 scripts/benchmark_report.py build/release/netlib_results.csv
```

Single problems or one method: `build/release/benchmarks/arclp_bench data/arcLP --problem qap15 --method arc`.

On the 51 Netlib problems of Table 1 from Yang 2018 (Apple M2 Max):

| | solved (rel. obj. error < 1e-6) | iterations (paper) | time |
|---|---|---|---|
| arcLP | 46/51, other 5 within 6e-6 | 1077 (1064) | 5.9 s |
| Mehrotra | 47/51 | 1073 (1099) | 5.6 s |

(with METIS and OpenMP; without them 7.6 s and 6.8 s)

Full table: [BENCHMARKS.md](BENCHMARKS.md).

## Differences from `arcLP.m`

All deliberate and measured:

1. **Rank deficiency.** The MATLAB option `d = 1` (Markowitz elimination to make A full rank, plus dropping columns with `x < 1e-6` during the iterations) is replaced by the standard interior-point remedy inside the Cholesky factorization: a pivot below `1e-14` times its own diagonal entry is replaced by `1e128`. This solves `qap15` in the paper's 27 iterations without the option.
2. **Status classification.** `arcLP.m` only tests the primal residual when the iterations stop early. A run whose dual residual never converged *and* whose primal objective lies far below the dual one (weak duality violated) is now reported as unbounded instead of optimal; non-finite iterates end as `ARCLP_NUMERICAL_ERROR` instead of "optimal", and non-finite input is rejected as `ARCLP_INVALID_INPUT`.
3. **Rule 7 of the presolve** tests "all nonzeros share a sign" instead of `min(A(i,:)) == 0`; the two differ only for completely dense rows, where the MATLAB test misses the reduction.
4. The μ-safeguard halving loop stops at `1e-300` instead of looping forever.
5. When the residual-growth safeguard restores the previous iterate, `y` is restored together with `x` and `s` (`arcLP.m` keeps the new `y`), so the final iterate is consistent; the returned `x` is unaffected.