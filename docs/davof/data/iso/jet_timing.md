# The jet in cross flow: the timing rows (coarse mesh, 32 ranks, dt 5e-7; prof/jet-run.sbatch)

## The published scheme

| mesh | row | rep | np | steps | clock s | exec s | wall s | min alpha | max overshoot | rc |
|---|---|---|---|---|---|---|---|---|---|---|
| coarse | RING | 1 | 32 | 574 | 583 | 578.25 | 587.67 | -7.414e+124 | 7.414e+124 | 143 |
| coarse | RING | full | 32 | 574 | 571 | 565.52 | 575.45 | -7.414e+124 | 7.414e+124 | 143 |
| coarse | RDFccs | 1 | 32 | 2000 | 1889 | 1875.47 | 1891.51 | -8.899e-13 | 2.220e-16 | 0 |
| coarse | RDF | 1 | 32 | 2000 | 1926 | 1911.62 | 1928.70 | -9.427e-13 | 2.220e-16 | 0 |
| coarse | RDFccs | 3 | 32 | 2000 | 1787 | 1773.51 | 1789.22 | -8.899e-13 | 2.220e-16 | 0 |
| coarse | RDF | 3 | 32 | 2000 | 1804 | 1790.02 | 1805.98 | -9.427e-13 | 2.220e-16 | 0 |
| coarse | RDFccs | 2 | 32 | 2000 | 1863 | 1851.05 | 1865.80 | -8.899e-13 | 2.220e-16 | 0 |
| coarse | RDF | 2 | 32 | 2000 | 1866 | 1853.14 | 1868.50 | -9.427e-13 | 2.220e-16 | 0 |
| coarse | RDFccs | full | 32 | 16000 | 16731 | 16626.92 | 16733.41 | -3.845e-06 | 6.217e-15 | 0 |
| coarse | RDF | full | 32 | 16000 | 17000 | 16887.98 | 17002.37 | -2.827e-06 | 9.992e-15 | 0 |
| coarse | RINGresetNone | 1 | 32 | 597 | 566 | 561.59 | 570.36 | -5.808e+129 | 7.031e+129 | 143 |
| coarse | RINGreset | 1 | 32 | 800 | 812 | 804.31 | 817.04 | -4.551e+205 | 7.552e+205 | 143 |

## The earlier (cubic) scheme, a record only: every row diverged

| mesh | row | rep | np | steps | clock s | exec s | wall s | min alpha | max overshoot | rc |
|---|---|---|---|---|---|---|---|---|---|---|
| coarse | RING | full | 32 | 623 | 635 | 629.11 | 640.07 | -2.116e+134 | 6.567e+133 | 143 |
| coarse | RING | 1 | 32 | 623 | 646 | 640.19 | 650.46 | -2.116e+134 | 6.567e+133 | 143 |
| coarse | RDF | 1 | 32 | 2000 | 2043 | 2026.47 | 2045.54 | -4.116e-06 | 4.441e-16 | 0 |
| coarse | RDFccs | 1 | 32 | 2000 | 2076 | 2059.44 | 2078.20 | -1.856e-06 | 4.441e-16 | 0 |
| coarse | RING | 1 | 32 | 455 | 465 | 460.71 | 469.21 | -1.498e+160 | 1.379e+160 | 143 |
| coarse | RING | full | 32 | 455 | 473 | 469.45 | 477.64 | -1.498e+160 | 1.379e+160 | 143 |
| coarse | RDFccs | full | 32 | 3594 | 4092 | 4062.71 | 4097.28 | -9.589e+178 | 1.297e+179 | 136 |
| coarse | RDF | 2 | 32 | 2000 | 2047 | 2032.29 | 2049.41 | -4.116e-06 | 4.441e-16 | 0 |
| coarse | RDFccs | 2 | 32 | 2000 | 2047 | 2031.72 | 2049.80 | -1.856e-06 | 4.441e-16 | 0 |
| coarse | RDF | 3 | 32 | 2000 | 2074 | 2057.81 | 2076.81 | -4.116e-06 | 4.441e-16 | 0 |
| coarse | RDFccs | 3 | 32 | 2000 | 2073 | 2057.11 | 2075.37 | -1.856e-06 | 4.441e-16 | 0 |
| coarse | RDF | full | 32 | 4063 | 4920 | 4883.28 | 4924.95 | -6.655e+169 | 5.516e+169 | 143 |

rep `full`: the 8 ms run (16,000 steps); rep 1-3: the 1 ms windows (2,000 steps); rc 143/137: the run killed after the floating-point exception or the MPI abort; RINGreset = scale + normalAndAlphaReset, RINGresetNone = none + normalAndAlphaReset (the bf70393 build, without the surfaces function object).
