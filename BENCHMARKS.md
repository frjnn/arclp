# Benchmarks

Three test sets, all measured on an Apple M2 Max (Accelerate BLAS, METIS ordering, OpenMP on) at tolerance 1e-8:

- [Netlib, Table 1 of Yang (2018)](#netlib-table-1-of-yang-2018): 51 small
  feasible problems, arcLP against Mehrotra's method.
- [Netlib infeasible problems](#netlib-infeasible-problems): 29 problems with
  no feasible point, which must end *infeasible*.
- [Mittelmann LPfeas benchmark](#mittelmann-lpfeas-benchmark): 49 large
  problems, against published HiGHS and MOSEK runs.

## Netlib, Table 1 of Yang (2018)

Reproduces the comparison of Yang (2018), Table 1: the 51 Netlib problems in standard form, solved by arcLP (Algorithm 5.1, the paper's `arclp2.m`) and by Mehrotra's predictor–corrector (`mehrotra.m`), with the same presolve, starting point and stopping rules. `Osa_60` is excluded as in the paper.

- **it (paper)**: iterations of this implementation, with the paper's count in parentheses. **†** marks runs that do not end *optimal* with a relative objective error below 1e-6 versus an independent HiGHS solution. All of them are `scfxm1–3` ending *stalled*, with objective errors from 4e-9 to 3e-7; the paper reports the same accuracy class for these problems, e.g. infeasibility 4.1e-4 for `scfxm3`.
- **time**: wall clock of the whole solve (presolve, ordering, symbolic analysis, iterations), best of 3.
- Reproduce with `cmake --build --preset release --target benchmark`, then `python3 scripts/benchmark_report.py build/release/netlib_results.csv`.

The iteration counts no longer follow `arclp2.m` row by row, because `x` and `(y, s)` move by separate angles (README, difference 4): 23 of 51 are within $\pm 2$ of the paper, and the total is 987 against the paper's 1064. The largest gaps are `agg` (54 vs 20) and `osa_14` (23 vs 42).

| problem | rows×cols (presolved) | arcLP it (paper) | Mehrotra it (paper) | arcLP rel. obj. error | arcLP ‖Ax−b‖ | arcLP time [s] | Mehrotra time [s] |
|---|---|---|---|---|---|---|---|
| adlittle | 53×135 | 15 (17) | 13 (15) | 8.4e-13 | 1.5e-09 | 0.0007 | 0.0003 |
| afiro | 8×32 | 9 (9) | 10 (9) | 1.1e-09 | 1.4e-11 | 0.0001 | 0.0000 |
| agg | 357×445 | 54 (20) | 32 (22) | 3.3e-11 | 6.2e-08 | 0.0063 | 0.0019 |
| agg2 | 484×725 | 20 (21) | 19 (20) | 4.0e-10 | 9.7e-07 | 0.0080 | 0.0058 |
| agg3 | 484×725 | 19 (20) | 18 (18) | 2.6e-09 | 9.7e-07 | 0.0077 | 0.0056 |
| bandm | 173×328 | 19 (20) | 17 (22) | 1.4e-11 | 1.6e-09 | 0.0026 | 0.0016 |
| beaconfd | 57×147 | 12 (11) | 11 (11) | 8.7e-12 | 5.2e-10 | 0.0010 | 0.0007 |
| blend | 47×87 | 11 (14) | 10 (14) | 1.1e-10 | 1.1e-10 | 0.0004 | 0.0002 |
| bnl1 | 393×1278 | 36 (34) | 40 (35) | 4.0e-10 | 9.4e-09 | 0.0164 | 0.0104 |
| bnl2 | 946×3005 | 33 (35) | 33 (38) | 2.8e-09 | 4.0e-08 | 0.0500 | 0.0347 |
| brandy | 102×207 | 17 (23) | 18 (19) | 3.0e-11 | 8.7e-07 | 0.0019 | 0.0013 |
| degen2 | 440×753 | 15 (19) | 13 (17) | 6.5e-10 | 3.9e-10 | 0.0079 | 0.0054 |
| degen3 | 1489×2590 | 20 (26) | 17 (22) | 4.4e-10 | 1.2e-09 | 0.0672 | 0.0513 |
| fffff800 | 476×980 | 24 (28) | 29 (31) | 2.5e-09 | 3.8e-08 | 0.0110 | 0.0090 |
| israel | 163×305 | 21 (25) | 20 (29) | 7.5e-11 | 3.1e-07 | 0.0041 | 0.0029 |
| lotfi | 108×321 | 21 (16) | 19 (18) | 3.0e-10 | 2.5e-08 | 0.0017 | 0.0008 |
| maros_r7 | 2152×7440 | 19 (20) | 18 (21) | 7.0e-13 | 2.8e-08 | 0.1398 | 0.1145 |
| osa_07 | 1081×25030 | 17 (30) | 22 (35) | 9.6e-10 | 1.5e-09 | 0.0675 | 0.0603 |
| osa_14 | 2300×54760 | 23 (42) | 23 (37) | 3.2e-09 | 5.3e-09 | 0.1790 | 0.1428 |
| osa_30 | 4313×104337 | 25 (39) | 26 (36) | 1.3e-10 | 6.5e-09 | 0.3563 | 0.2976 |
| qap12 | 3048×8712 | 21 (22) | 20 (24) | 6.4e-11 | 2.6e-11 | 0.6631 | 0.5863 |
| qap15 | 6105×22050 | 23 (27) | 23 (44) | 5.7e-09 | 1.7e-10 | 3.2491 | 3.1814 |
| qap8 | 848×1568 | 10 (12) | 9 (13) | 2.1e-10 | 1.5e-13 | 0.0265 | 0.0220 |
| sc105 | 44×102 | 11 (11) | 10 (11) | 2.8e-09 | 6.5e-11 | 0.0003 | 0.0001 |
| sc205 | 89×201 | 12 (12) | 12 (12) | 1.6e-12 | 4.5e-10 | 0.0006 | 0.0003 |
| sc50a | 19×47 | 11 (10) | 8 (9) | 1.2e-11 | 5.5e-11 | 0.0001 | 0.0000 |
| sc50b | 14×42 | 8 (10) | 8 (8) | 2.0e-09 | 1.0e-10 | 0.0001 | 0.0000 |
| scagr25 | 222×422 | 17 (19) | 16 (18) | 2.2e-09 | 3.3e-07 | 0.0020 | 0.0009 |
| scagr7 | 60×116 | 14 (17) | 15 (17) | 1.6e-09 | 1.0e-08 | 0.0005 | 0.0002 |
| scfxm1 | 217×480 | 34 (21) † | 33 (22) † | 3.7e-09 | 1.5e-04 | 0.0068 | 0.0041 |
| scfxm2 | 437×963 | 38 (24) † | 21 (26) | 5.9e-09 | 1.4e-04 | 0.0158 | 0.0049 |
| scfxm3 | 657×1446 | 36 (25) † | 20 (23) | 2.6e-07 | 8.6e-04 | 0.0250 | 0.0072 |
| scrs8 | 115×893 | 20 (27) | 19 (30) | 7.6e-10 | 1.3e-10 | 0.0052 | 0.0030 |
| scsd1 | 77×760 | 10 (11) | 10 (13) | 1.3e-11 | 3.6e-13 | 0.0016 | 0.0008 |
| scsd6 | 147×1350 | 12 (16) | 13 (16) | 4.7e-10 | 1.8e-11 | 0.0033 | 0.0016 |
| scsd8 | 397×2750 | 12 (15) | 11 (14) | 9.9e-12 | 1.5e-11 | 0.0061 | 0.0026 |
| sctap1 | 284×644 | 18 (20) | 16 (27) | 2.3e-11 | 1.2e-09 | 0.0036 | 0.0018 |
| sctap2 | 1033×2443 | 18 (20) | 16 (21) | 1.4e-10 | 6.0e-11 | 0.0138 | 0.0076 |
| sctap3 | 1408×3268 | 18 (21) | 15 (22) | 6.0e-11 | 9.9e-11 | 0.0188 | 0.0100 |
| share1b | 100×236 | 19 (26) | 20 (25) | 2.8e-09 | 3.1e-06 | 0.0013 | 0.0007 |
| share2b | 84×150 | 14 (15) | 13 (15) | 1.9e-12 | 3.6e-09 | 0.0007 | 0.0004 |
| ship04l | 292×1905 | 15 (18) | 13 (18) | 7.4e-11 | 2.7e-09 | 0.0065 | 0.0027 |
| ship04s | 216×1281 | 15 (19) | 14 (20) | 2.1e-09 | 3.1e-09 | 0.0044 | 0.0019 |
| ship08l | 470×3121 | 17 (20) | 16 (22) | 3.4e-09 | 5.8e-10 | 0.0124 | 0.0057 |
| ship08s | 274×1600 | 15 (19) | 15 (20) | 6.1e-09 | 1.0e-09 | 0.0058 | 0.0027 |
| ship12l | 610×4171 | 20 (21) | 19 (21) | 6.2e-09 | 1.5e-09 | 0.0179 | 0.0074 |
| ship12s | 340×1943 | 18 (20) | 18 (19) | 9.6e-10 | 1.5e-09 | 0.0071 | 0.0030 |
| stocfor1 | 28×76 | 12 (14) | 11 (14) | 4.9e-09 | 2.6e-09 | 0.0003 | 0.0001 |
| stocfor2 | 760×1648 | 18 (22) | 19 (22) | 2.2e-09 | 1.9e-07 | 0.0104 | 0.0067 |
| stocfor3 | 5968×12834 | 31 (37) | 31 (38) | 3.6e-11 | 6.8e-08 | 0.1338 | 0.0833 |
| truss | 1000×8806 | 20 (24) | 21 (26) | 8.6e-10 | 1.9e-09 | 0.0445 | 0.0251 |

- **arcLP**: 48/51 solved to relative objective error < 1e-06; 987 iterations in total (paper: 1064); 5.22 s in total.
- **Mehrotra**: 50/51 solved to relative objective error < 1e-06; 913 iterations in total (paper: 1099); 4.72 s in total.

## Netlib infeasible problems

The 29 problems of the Netlib `infeas` set. A run counts only when it ends *infeasible*, which needs a certificate (README, difference 2); 0 iterations means presolve found the contradiction. arcLP proves 24 infeasible and Mehrotra 25; the rest stall.

| problem | arcLP | Mehrotra |
|---|---|---|
| bgdbg1 | infeasible (0) | infeasible (0) |
| bgetam | infeasible (7) | infeasible (7) |
| bgindy | infeasible (5) | infeasible (5) |
| bgprtr | infeasible (7) | infeasible (8) |
| box1 | infeasible (5) | infeasible (4) |
| ceria3d | infeasible (8) | infeasible (8) |
| chemcom | infeasible (7) | infeasible (7) |
| cplex1 | infeasible (6) | infeasible (7) |
| cplex2 | stalled (42) | stalled (38) |
| ex72a | infeasible (5) | infeasible (4) |
| ex73a | infeasible (5) | infeasible (4) |
| forest6 | infeasible (7) | infeasible (10) |
| galenet | infeasible (0) | infeasible (0) |
| gosh | infeasible (14) | infeasible (14) |
| gran | infeasible (0) | infeasible (0) |
| greenbea | infeasible (0) | infeasible (0) |
| itest2 | infeasible (0) | infeasible (0) |
| itest6 | infeasible (0) | infeasible (0) |
| klein1 | stalled (52) | stalled (32) |
| klein2 | infeasible (45) | infeasible (34) |
| klein3 | stalled (57) | infeasible (38) |
| mondou2 | infeasible (6) | infeasible (7) |
| pang | stalled (79) | stalled (40) |
| pilot4i | infeasible (0) | infeasible (0) |
| qual | stalled (33) | stalled (37) |
| reactor | infeasible (0) | infeasible (0) |
| refinery | infeasible (14) | infeasible (19) |
| vol1 | infeasible (27) | infeasible (34) |
| woodinfe | infeasible (0) | infeasible (0) |

## Mittelmann LPfeas benchmark

The 49 instances of Hans Mittelmann's [LPfeas benchmark](https://plato.asu.edu/ftp/lpfeas.html), solved by arcLP one at a time with a limit of 900 s and 64 GB.

- **status, it, time**: arcLP's result, iterations and wall clock (one run). *time limit* means the 900 s ran out; *killed (memory)* means macOS ended the process.
- **rel. obj. error**: versus the closest of the HiGHS, MOSEK and COPT objectives in the benchmark's logs, which differ among themselves by up to 1e-4 on some instances.
- **HiGHS, MOSEK**: times from the benchmark's own logs, on a different machine (i7-11700K, 64 GB), shown only where that solver reached an optimum. They give a sense of scale, not a direct comparison.

arcLP ends *optimal* on 32 instances, all to a relative objective error of at most 1.4e-6. Of the rest, 9 stall or hit a numerical error, 6 run out of time and 2 out of memory. On the 31 instances that arcLP, HiGHS and MOSEK all solve, the shifted (by 10 s) geometric mean of the times is 30.0 s for arcLP, 66.0 s for HiGHS and 11.3 s for MOSEK; arcLP is faster than HiGHS on 25 of them and faster than MOSEK on 7.

| instance | rows×cols | status | it | rel. obj. error | time [s] | peak memory [GB] | HiGHS [s] | MOSEK [s] |
|---|---|---|---|---|---|---|---|---|
| a2864 | 22117×200787 | optimal | 18 | 1e-09 | 135.6 | 5.5 | 196.4 | 0.4 |
| bdry2 | 376500×250998 | stalled | 41 | 1e-04 | 21.1 | 0.9 | — | 16.4 |
| cont1 | 160792×40398 | stalled | 32 | 6e-07 | 2.9 | 0.2 | 6.8 | 0.3 |
| cont11 | 160792×80396 | optimal | 20 | 2e-08 | 1.6 | 0.2 | 10.4 | 0.5 |
| datt256_lp | 11077×262144 | optimal | 11 | 5e-09 | 1.3 | 0.3 | 80.3 | 1.7 |
| degme | 185501×659415 | optimal | 51 | 4e-09 | 55.8 | 2.0 | 71.8 | 364.8 |
| dlr1 | 1735470×9142907 | stalled | 137 | 2e+00 | 376.7 | 6.3 | — | 119.6 |
| dlr2 | — | time limit | — | — | > 900 | 50.9 | — | 3314.2 |
| Dual2_5000 | 30000600×33050602 | out of memory | 0 | — | 8.1 | 15.6 | — | — |
| ex10 | 69608×17680 | optimal | 24 | 0e+00 | 59.4 | 2.0 | 44.1 | 1.3 |
| fhnw-binschedule1 | 772872×1141653 | optimal | 50 | 2e-11 | 413.7 | 9.1 | 382.6 | 46.4 |
| fome13 | 48568×97840 | optimal | 44 | 5e-09 | 5.2 | 0.2 | 9.7 | 0.6 |
| graph40-40_lp | 360900×102600 | optimal | 13 | 0e+00 | 20.7 | 2.2 | 11.7 | 0.1 |
| irish-electricity | 104259×61728 | optimal | 54 | 4e-09 | 5.3 | 0.3 | — | 3.0 |
| L1_sixm1000obs | — | killed (memory) | — | — | — | 58.4 | 1225.3 | — |
| L1_sixm250obs | 986069×428032 | time limit | — | — | > 900 | 60.7 | 109.6 | — |
| L2CTA3D | — | time limit | — | — | > 900 | 13.6 | — | 2185.7 |
| Linf_520c | 93326×69004 | numerical error | 38 | 3e+00 | 78.3 | 1.2 | 512.3 | 6.0 |
| neos | 479119×36786 | stalled | 107 | 9e-08 | 35.2 | 0.8 | 86.5 | 9.1 |
| neos-3025225_lp | 91572×69846 | optimal | 63 | 5e-09 | 63.3 | 2.3 | 62.4 | 5.8 |
| neos-5052403-cygnet | 38268×32868 | optimal | 25 | 4e-10 | 19.3 | 1.2 | 9.7 | 4.0 |
| neos-5251015_lp | 486531×136971 | optimal | 14 | 0e+00 | 4.4 | 0.7 | 687.2 | 5.2 |
| neos3 | 512209×6624 | optimal | 14 | 4e-10 | 3.8 | 0.7 | 36.1 | 0.2 |
| ns1687037 | 50622×43749 | stalled | 166 | 4e-05 | 41.6 | 0.3 | 120.5 | — |
| ns1688926 | 32768×16587 | stalled | 7 | 6e+06 | 3.2 | 0.5 | — | 0.5 |
| nug08-3rd | 19728×20448 | optimal | 13 | 1e-09 | 16.2 | 1.3 | 52.5 | 0.6 |
| pds-100 | 156243×505360 | optimal | 81 | 5e-11 | 42.6 | 0.8 | 51.2 | 19.9 |
| physiciansched3-3 | 266227×79555 | stalled | 77 | 1e-08 | 22.2 | 0.6 | — | 16.6 |
| Primal2_1000 | 1299380×2559380 | optimal | 87 | 2e-08 | 649.9 | 6.0 | 134.7 | 23.4 |
| qap15 | 6330×22275 | optimal | 23 | 5e-09 | 3.3 | 0.2 | 10.4 | 0.3 |
| rail02 | 95791×270869 | optimal | 78 | 3e-08 | 30.4 | 0.6 | 57.4 | 32.6 |
| rail4284 | 4284×1092610 | optimal | 76 | 8e-09 | 22.7 | 1.2 | 95.0 | 9.4 |
| rmine15_lp | 358395×42438 | time limit | — | — | > 900 | 8.1 | 88.0 | 26.4 |
| s100 | 14733×364417 | optimal | 100 | 8e-10 | 13.3 | 0.4 | 23.0 | 3.7 |
| s250r10 | 10962×273142 | optimal | 112 | 5e-09 | 9.3 | 0.3 | 17.1 | 2.9 |
| s82 | 87878×1690631 | optimal | 169 | 5e-09 | 114.4 | 1.8 | 376.8 | 33.5 |
| savsched1 | 295989×328575 | optimal | 23 | 4e-10 | 10.3 | 0.9 | 74.8 | 8.7 |
| scpm1_lp | 5000×500000 | optimal | 27 | 4e-10 | 10.5 | 1.3 | 66.5 | 14.7 |
| set-cover-model | 10000×1102008 | optimal | 55 | 2e-09 | 95.5 | 4.3 | 384.9 | 131.3 |
| shs1023 | 133944×444625 | optimal | 81 | 1e-06 | 7.7 | 0.4 | 47.5 | 12.5 |
| square41 | 40160×62234 | optimal | 44 | 3e-09 | 13.2 | 0.8 | 15.2 | 1.8 |
| stat96v2 | 29089×957432 | stalled | 113 | 6e-06 | 26.1 | 0.8 | 498.8 | 19.3 |
| stormG2_1000 | 528185×1259121 | time limit | — | — | > 900 | 11.7 | 175.5 | 21.3 |
| stp3d | 159488×204880 | optimal | 53 | 1e-09 | 12.5 | 0.4 | 23.9 | 4.5 |
| supportcase10 | 165684×14770 | optimal | 29 | 6e-11 | 5.7 | 0.3 | 22.0 | 3.0 |
| thk_48 | — | time limit | — | — | > 900 | 42.0 | — | 5716.4 |
| thk_63 | 5694387×7701112 | optimal | 47 | 2e-10 | 548.9 | 12.7 | 1177.0 | 132.8 |
| tpl-tub-ws1617_lp | 1154615×747601 | optimal | 158 | 4e-08 | 80.0 | 1.6 | 225.4 | 21.5 |
| woodlands09 | 194599×382147 | optimal | 19 | 1e-11 | 11.1 | 1.0 | 40.1 | 3.1 |

### Start-point spread

The start point spreads the iterates by a factor (`START_SPREAD` in `src/ipm/arc_rule.c`). Spread 1 solves two more Netlib problems (50/51 instead of 48/51, gaining `scfxm2` and `scfxm3`) but takes 1031 iterations instead of 987. On the infeasible set both spreads prove 24. On the 41 Mittelmann instances that end within the limits, spread 1 solves 30 instead of 32, losing `fhnw-binschedule1` (time limit) and `square41` (stalled). On the 30 both solve it takes 1563 iterations instead of 1582, but 2220 s instead of 2060 s. The default stays at 4. Instances whose status or iteration count (by more than 10) differ:

| instance | spread 4: status, it, time [s] | spread 1: status, it, time [s] |
|---|---|---|
| dlr1 | stalled, 137, 376.7 | stalled, 21, 72.6 |
| fhnw-binschedule1 | optimal, 50, 413.7 | time limit, —, > 900 |
| ns1687037 | stalled, 166, 41.6 | stalled, 204, 48.8 |
| Primal2_1000 | optimal, 87, 649.9 | optimal, 120, 875.5 |
| rail02 | optimal, 78, 30.4 | optimal, 110, 42.3 |
| rail4284 | optimal, 76, 22.7 | optimal, 57, 18.3 |
| s250r10 | optimal, 112, 9.3 | optimal, 97, 8.3 |
| scpm1_lp | optimal, 27, 10.5 | optimal, 40, 13.7 |
| shs1023 | optimal, 81, 7.7 | optimal, 69, 7.0 |
| square41 | optimal, 44, 13.2 | stalled, 104, 25.0 |
| stat96v2 | stalled, 113, 26.1 | stalled, 75, 15.7 |
| tpl-tub-ws1617_lp | optimal, 158, 80.0 | optimal, 147, 74.3 |
