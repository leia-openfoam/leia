
### sw2_shear2D_N64 (OF Co 0.2; diagnostics on; 6 cases concurrently on 6 cores)

| variant | status | steps | first step with violation > 1e-3 V_k | run max violation before limiting (V_k) | run max violation after limiting (V_k) | S1 faces/step mean (max) | S2 faces/step mean (max) [alpha_up<0.5 / >=0.5 mean] | limiter faces/step mean (max) | limiter cells/step mean (max) | limiter volume total (cell volumes) | limiter passes max | E_shape | Emass (csv) | wall clock (s) | ranks |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| base | SIGFPE (last DIAG step 750) | 750 | 57 | 0.158 | 0.155 | 0 (0) | 0 (0) [0 / 0] | 21.1 (71) | 31.8 (100) | 172 | 3 |  |  | 74 | 1 |
| S1 | completed | 3328 | 57 | 0.245 | 0.192 | 5.25 (16) | 0 (0) [0 / 0] | 193 (428) | 244 (533) | 5.92e+03 | 3 | 0.000702 | 5.85e-10 | 491 | 1 |
| S2low | SIGFPE (last DIAG step 2727) | 2727 | 145 | 0.107 | 0.00127 | 0 (0) | 4.8e+03 (6778) [4.8e+03 / 0] | 15.6 (133) | 28.9 (216) | 216 | 3 |  |  | 745 | 1 |
| S2all | SIGFPE (last DIAG step 1292) | 1292 | 294 | 0.00155 | 0.000144 | 0 (0) | 3.51e+03 (5844) [3.19e+03 / 311] | 23.5 (313) | 32.6 (352) | 3.07 | 3 |  |  | 421 | 1 |
| S1+S2low | completed | 3328 | 145 | 0.146 | 0.00254 | 16.3 (122) | 5.01e+03 (6726) [5.01e+03 / 0] | 65 (743) | 95.2 (873) | 305 | 3 | 0.000292 | 1.88e-08 | 976 | 1 |
| S1+S2all | completed | 3328 | 294 | 0.00499 | 0.00327 | 53.6 (118) | 5.01e+03 (6334) [4.72e+03 / 297] | 160 (693) | 210 (769) | 90.8 | 3 | 0.000249 | 2.88e-08 | 1.22e+03 | 1 |

pieces stored (Co 0.2, N 64, ranks 1): status completed, E_shape 1.82e-04, wall clock 284 s

plicRDF stored (Co 0.2, N 64, ranks 1): status completed, E_shape 3.12e-04, wall clock 13 s

### sw2_shear3D_N32 (OF Co 0.2; diagnostics on; 6 cases concurrently on 6 cores)

| variant | status | steps | first step with violation > 1e-3 V_k | run max violation before limiting (V_k) | run max violation after limiting (V_k) | S1 faces/step mean (max) | S2 faces/step mean (max) [alpha_up<0.5 / >=0.5 mean] | limiter faces/step mean (max) | limiter cells/step mean (max) | limiter volume total (cell volumes) | limiter passes max | E_shape | Emass (csv) | wall clock (s) | ranks |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| base | SIGFPE (last DIAG step 236) | 236 | 5 | 0.254 | 0.183 | 0 (0) | 0 (0) [0 / 0] | 167 (2616) | 213 (2299) | 259 | 3 |  |  | 88 | 1 |
| S1 | abort step 253 | 253 | 5 | 0.525 | 0.508 | 2.62 (63) | 0 (0) [0 / 0] | 612 (10249) | 585 (8621) | 1.5e+03 | 3 |  |  | 103 | 1 |
| S2low | SIGFPE (last DIAG step 275) | 275 | 9 | 0.0563 | 4.46e-05 | 0 (0) | 2.01e+04 (34577) [2.01e+04 / 0] | 53.1 (620) | 78.7 (802) | 1.47 | 3 |  |  | 210 | 1 |
| S2all | SIGFPE (last DIAG step 276) | 276 | 36 | 0.00156 | 4.46e-05 | 0 (0) | 2.07e+04 (34749) [2.01e+04 / 583] | 51.6 (892) | 73.8 (1141) | 0.28 | 3 |  |  | 225 | 1 |
| S1+S2low | completed | 720 | 9 | 0.0563 | 0.000743 | 833 (1881) | 7.51e+04 (130695) [7.51e+04 / 0] | 1.08e+04 (32074) | 1.07e+04 (28000) | 61.4 | 3 | 0.0125 | 0.00617 | 3.04e+03 | 1 |
| S1+S2all | completed | 720 | 36 | 0.00542 | 0.00171 | 884 (2003) | 7.81e+04 (131699) [7.75e+04 / 579] | 1.15e+04 (35456) | 1.12e+04 (30313) | 59.9 | 3 | 0.0155 | 0.00924 | 3.09e+03 | 1 |

pieces stored (Co 0.2, N 32, ranks 1): status completed, E_shape 6.16e-03, wall clock 548 s

plicRDF stored (Co 0.2, N 32, ranks 1): status completed, E_shape 4.71e-03, wall clock 19 s

### sw2b_deformation3D_N32 (OF Co 0.2; diagnostics on; 6 cases concurrently on 6 cores)

| variant | status | steps | first step with violation > 1e-3 V_k | run max violation before limiting (V_k) | run max violation after limiting (V_k) | S1 faces/step mean (max) | S2 faces/step mean (max) [alpha_up<0.5 / >=0.5 mean] | limiter faces/step mean (max) | limiter cells/step mean (max) | limiter volume total (cell volumes) | limiter passes max | E_shape | Emass (csv) | wall clock (s) | ranks |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| base | abort step 236 | 236 | 1 | 0.502 | 0.553 | 0 (0) | 0 (0) [0 / 0] | 284 (7303) | 277 (5988) | 743 | 3 |  |  | 39 | 1 |
| S1 | abort step 236 | 236 | 1 | 0.502 | 0.553 | 0 (0) | 0 (0) [0 / 0] | 284 (7303) | 277 (5988) | 743 | 3 |  |  | 39 | 1 |
| S2low | completed | 1200 | 2 | 0.0876 | 0.00104 | 0 (0) | 4.73e+04 (75645) [4.73e+04 / 0] | 4.47e+03 (8709) | 6.09e+03 (11435) | 104 | 3 | 0.00989 | 1.73e-17 | 3.87e+03 | 1 |
| S2all | completed | 1200 | 74 | 0.00389 | 0.000596 | 0 (0) | 4.75e+04 (74450) [4.71e+04 / 416] | 4.47e+03 (8622) | 6.09e+03 (11445) | 74.1 | 3 | 0.0098 | 3.64e-17 | 3.88e+03 | 1 |
| S1+S2low | completed | 1200 | 2 | 0.0876 | 0.00104 | 0 (0) | 4.59e+04 (73199) [4.59e+04 / 0] | 4.47e+03 (8709) | 6.09e+03 (11435) | 104 | 3 | 0.00989 | 1.73e-17 | 3.84e+03 | 1 |
| S1+S2all | completed | 1200 | 74 | 0.00389 | 0.000596 | 0 (0) | 4.61e+04 (71971) [4.57e+04 / 416] | 4.47e+03 (8622) | 6.09e+03 (11445) | 74.1 | 3 | 0.0098 | 3.64e-17 | 3.87e+03 | 1 |

pieces stored (Co 0.2, N 32, ranks 1): status completed, E_shape 9.65e-03, wall clock 1926 s

plicRDF stored (Co 0.2, N 32, ranks 1): status completed, E_shape 7.93e-03, wall clock 37 s
