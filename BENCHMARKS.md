# Netlib benchmark

Reproduces the comparison of Yang (2018), Table 1: the 51 Netlib problems in
standard form, solved by arcLP (Algorithm 5.1, the paper's `arclp2.m`) and by
Mehrotra's predictor–corrector (`mehrotra.m`), with the same presolve,
starting point and stopping rules (tolerance 1e-8). `Osa_60` is excluded as in
the paper.

- **it (paper)**: iterations of this implementation, with the paper's count
  in parentheses. **†** marks runs whose relative objective error versus an
  independent HiGHS solution is above 1e-6 (all still end with status
  *optimal*; the paper reports the same accuracy class for these problems,
  e.g. infeasibility 4.1e-4 for `scfxm3`).

- **time**: wall clock of the whole solve (presolve, ordering, symbolic
  analysis, iterations), best of 3, Apple M2 Max, Accelerate BLAS, METIS
  ordering available, OpenMP on.
- Reproduce with `cmake --build --preset release --target benchmark`, then
  `python3 scripts/benchmark_report.py build/release/netlib_results.csv`.

Agreement with the paper's `arclp2.m` iteration counts: 34 of 51 exact, 46
within $\pm 1$, 48 within $\pm 2$ (largest gap: `osa_30`, 44 vs 39). Exact agreement on
every row is not expected.

| problem | rows×cols (presolved) | arcLP it (paper) | Mehrotra it (paper) | arcLP rel. obj. error | arcLP ‖Ax−b‖ | arcLP time [s] | Mehrotra time [s] |
|---|---|---|---|---|---|---|---|
| adlittle | 54×136 | 16 (17) | 15 (15) | 4.9e-13 | 5.4e-11 | 0.0010 | 0.0006 |
| afiro | 8×32 | 9 (9) | 9 (9) | 1.6e-09 | 1.2e-12 | 0.0001 | 0.0001 |
| agg | 391×479 | 20 (20) | 22 (22) | 4.3e-12 | 5.3e-06 | 0.0066 | 0.0048 |
| agg2 | 514×755 | 21 (21) | 20 (20) | 5.1e-11 | 2.6e-08 | 0.0085 | 0.0064 |
| agg3 | 514×755 | 20 (20) | 18 (18) | 4.0e-10 | 4.9e-08 | 0.0082 | 0.0060 |
| bandm | 192×347 | 20 (20) | 22 (22) | 2.8e-09 | 3.3e-11 | 0.0025 | 0.0020 |
| beaconfd | 57×147 | 11 (11) | 11 (11) | 4.4e-12 | 1.5e-11 | 0.0008 | 0.0006 |
| blend | 49×89 | 14 (14) | 14 (14) | 2.6e-11 | 2.0e-12 | 0.0004 | 0.0003 |
| bnl1 | 429×1314 | 34 (34) | 35 (35) | 1.5e-10 | 1.4e-09 | 0.0128 | 0.0090 |
| bnl2 | 1007×3066 | 35 (35) | 36 (38) | 9.3e-10 | 1.2e-09 | 0.0493 | 0.0388 |
| brandy | 113×218 | 23 (23) | 19 (19) | 5.3e-10 | 2.6e-07 | 0.0022 | 0.0015 |
| degen2 | 440×753 | 18 (19) | 18 (17) | 7.2e-12 | 8.8e-12 | 0.0083 | 0.0069 |
| degen3 | 1490×2591 | 24 (26) | 24 (22) | 1.7e-10 | 1.1e-10 | 0.0746 | 0.0686 |
| fffff800 | 487×991 | 31 (28) | 31 (31) | 1.3e-09 | 6.3e-09 | 0.0126 | 0.0097 |
| israel | 174×316 | 26 (25) | 21 (29) † | 5.7e-11 | 3.3e-08 | 0.0047 | 0.0032 |
| lotfi | 113×326 | 16 (16) † | 18 (18) | 5.9e-06 | 1.5e-08 | 0.0011 | 0.0007 |
| maros_r7 | 2152×7440 | 20 (20) | 21 (21) | 2.0e-13 | 1.7e-09 | 0.1420 | 0.1273 |
| osa_07 | 1081×25030 | 31 (30) | 31 (35) | 3.1e-10 | 1.0e-09 | 0.1050 | 0.0768 |
| osa_14 | 2300×54760 | 41 (42) | 31 (37) † | 7.3e-10 | 7.3e-09 | 0.2720 | 0.1629 |
| osa_30 | 4313×104337 | 44 (39) | 35 (36) | 3.1e-10 | 2.5e-08 | 0.5475 | 0.3618 |
| qap12 | 3048×8712 | 22 (22) | 25 (24) | 8.8e-07 | 2.0e-09 | 0.6958 | 0.7158 |
| qap15 | 6105×22050 | 25 (27) † | 27 (44) | 3.4e-06 | 1.6e-08 | 3.6005 | 3.7954 |
| qap8 | 848×1568 | 12 (12) | 13 (13) † | 2.6e-11 | 9.5e-12 | 0.0308 | 0.0300 |
| sc105 | 44×102 | 11 (11) | 11 (11) | 4.5e-09 | 3.3e-12 | 0.0002 | 0.0001 |
| sc205 | 89×201 | 12 (12) | 12 (12) | 1.9e-08 | 6.5e-11 | 0.0005 | 0.0003 |
| sc50a | 19×47 | 10 (10) | 9 (9) | 2.6e-08 | 5.5e-13 | 0.0001 | 0.0000 |
| sc50b | 14×42 | 10 (10) | 8 (8) | 1.2e-10 | 3.9e-13 | 0.0001 | 0.0000 |
| scagr25 | 343×543 | 19 (19) | 18 (18) | 1.1e-12 | 1.9e-08 | 0.0027 | 0.0017 |
| scagr7 | 91×147 | 17 (17) | 17 (17) | 1.5e-12 | 1.4e-09 | 0.0006 | 0.0004 |
| scfxm1 | 238×500 | 21 (21) † | 27 (22) | 1.4e-05 | 8.0e-05 | 0.0034 | 0.0029 |
| scfxm2 | 479×1003 | 25 (24) † | 26 (26) | 2.6e-06 | 6.7e-05 | 0.0081 | 0.0058 |
| scfxm3 | 720×1506 | 25 (25) † | 29 (23) † | 1.4e-06 | 4.7e-05 | 0.0124 | 0.0099 |
| scrs8 | 115×893 | 29 (27) | 30 (30) | 2.3e-08 | 2.7e-10 | 0.0056 | 0.0039 |
| scsd1 | 77×760 | 11 (11) | 13 (13) | 6.0e-08 | 9.5e-15 | 0.0014 | 0.0008 |
| scsd6 | 147×1350 | 16 (16) | 16 (16) | 1.1e-09 | 1.9e-13 | 0.0034 | 0.0017 |
| scsd8 | 397×2750 | 14 (15) | 14 (14) | 2.3e-10 | 3.2e-13 | 0.0056 | 0.0028 |
| sctap1 | 284×644 | 21 (20) | 24 (27) | 4.4e-11 | 1.4e-11 | 0.0035 | 0.0026 |
| sctap2 | 1033×2443 | 21 (20) | 21 (21) | 8.4e-11 | 1.5e-12 | 0.0136 | 0.0095 |
| sctap3 | 1408×3268 | 21 (21) | 22 (22) | 8.4e-10 | 1.3e-12 | 0.0187 | 0.0137 |
| share1b | 102×238 | 27 (26) | 25 (25) | 1.4e-09 | 1.3e-07 | 0.0014 | 0.0008 |
| share2b | 87×153 | 15 (15) | 15 (15) | 2.7e-09 | 1.1e-10 | 0.0006 | 0.0004 |
| ship04l | 292×1905 | 18 (18) | 18 (18) | 2.1e-12 | 1.3e-10 | 0.0060 | 0.0032 |
| ship04s | 216×1281 | 19 (19) | 19 (20) | 1.7e-12 | 2.4e-09 | 0.0043 | 0.0022 |
| ship08l | 470×3121 | 21 (20) | 22 (22) | 2.1e-12 | 8.8e-12 | 0.0118 | 0.0069 |
| ship08s | 274×1600 | 19 (19) | 20 (20) | 7.6e-12 | 3.5e-11 | 0.0057 | 0.0031 |
| ship12l | 610×4171 | 21 (21) | 22 (21) | 2.1e-13 | 1.7e-11 | 0.0144 | 0.0078 |
| ship12s | 340×1943 | 20 (20) | 19 (19) | 3.5e-11 | 9.8e-10 | 0.0061 | 0.0028 |
| stocfor1 | 34×82 | 14 (14) | 15 (14) | 4.8e-12 | 5.4e-11 | 0.0003 | 0.0002 |
| stocfor2 | 766×1654 | 22 (22) | 22 (22) | 9.3e-12 | 7.4e-09 | 0.0107 | 0.0075 |
| stocfor3 | 5974×12840 | 41 (37) | 38 (38) | 1.2e-11 | 6.6e-08 | 0.1450 | 0.0979 |
| truss | 1000×8806 | 24 (24) | 25 (26) | 4.1e-13 | 5.2e-11 | 0.0416 | 0.0257 |

- **arcLP**: 46/51 solved to relative objective error < 1e-06; 1077 iterations in total (paper: 1064); 5.92 s in total.
- **Mehrotra**: 47/51 solved to relative objective error < 1e-06; 1073 iterations in total (paper: 1099); 5.64 s in total.
