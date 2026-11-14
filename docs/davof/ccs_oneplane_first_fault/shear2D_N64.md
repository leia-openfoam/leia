## shear2D_N64 OF Co 0.2 (final build, pieces comparison every step)

steps with a DIAG line: 750; last step 750, t = 1.80288; abort/trap lines: ['trapFpe: Floating point exception trapping enabled (FOAM_SIGFPE).', '#1  Foam::sigFpe::sigHandler(int) in ~/OpenFOAM/OpenFOAM-v2606/platforms/linux64GccDPInt32Opt/lib/libOpenFOAM.so']

### Timeline

| event | step | t |
|---|---|---|
| first step with a pre-limiter violation outside [-1e-10, 1+1e-10] | 2 | 0.00480769 |
| first step with a no-plane cut face anywhere | 5 | 0.0120192 |
| first step with a violating cell that has a no-plane cut face | 164 | 0.394231 |
| first step with a violation > 1e-3 | 57 | 0.137019 |
| first step with any limiter volume moved | 2 | 0.00480769 |
| first step with limiter volume > 1e-3 cell volumes | 57 | 0.137019 |
| last step printed (crash or abort or end) | 750 | 1.80288 |

first violation step 2: 1 cells, min alpha -1.54e-05, max alpha - 1 0, summed violation 1.54e-05 cell volumes, no-plane faces uncut 152 cut 0

### The first 50 steps with violations (sums over those steps)

steps 2..56 (50 steps): violating cells 112; with a no-plane cut face 0; with a no-plane face (cut or uncut) 0; with an unprocessed face 0; with all faces with planes 112; all faces with planes and a face |one - pieces| > 0.1 V_k 0 (compare steps among them: 50); steps with max violation < 1e-3: 50; steps with a cut no-plane violating cell: 0

| step | cells out | min alpha | max alpha-1 | summed viol | cells > 1e-3 | cut no-plane | no-plane | unproc. | all planes | all planes, diff > 0.1 | limiter vol | compare Linf | compare n>0.1 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 1 | -1.54e-05 | 0 | 1.54e-05 | 0 | 0 | 0 | 0 | 1 | 0 | 1.54e-05 | 0.00101 | 0 |
| 3 | 2 | -2e-06 | 2.23e-05 | 2.43e-05 | 0 | 0 | 0 | 0 | 2 | 0 | 2.43e-05 | 0.000857 | 0 |
| 4 | 2 | -0.000198 | 1.53e-05 | 0.000213 | 0 | 0 | 0 | 0 | 2 | 0 | 0.000213 | 0.00105 | 0 |
| 5 | 4 | -0.000131 | 0.000179 | 0.000439 | 0 | 0 | 0 | 0 | 4 | 0 | 0.000484 | 0.000758 | 0 |
| 6 | 1 | 0 | 6.02e-05 | 6.02e-05 | 0 | 0 | 0 | 0 | 1 | 0 | 6.02e-05 | 0.000334 | 0 |
| 7 | 1 | -3.21e-05 | 0 | 3.21e-05 | 0 | 0 | 0 | 0 | 1 | 0 | 3.21e-05 | 0.000446 | 0 |
| 8 | 5 | -0.000452 | 0 | 0.000774 | 0 | 0 | 0 | 0 | 5 | 0 | 0.000774 | 0.0014 | 0 |
| 9 | 3 | 0 | 0.000138 | 0.000295 | 0 | 0 | 0 | 0 | 3 | 0 | 0.000295 | 0.0022 | 0 |
| 11 | 1 | -7.05e-05 | 0 | 7.05e-05 | 0 | 0 | 0 | 0 | 1 | 0 | 7.05e-05 | 0.000268 | 0 |
| 12 | 2 | -0.000275 | 0 | 0.000357 | 0 | 0 | 0 | 0 | 2 | 0 | 0.000357 | 0.000791 | 0 |
| 13 | 4 | -0.000141 | 0 | 0.000283 | 0 | 0 | 0 | 0 | 4 | 0 | 0.000307 | 0.000977 | 0 |
| 14 | 1 | 0 | 1.1e-06 | 1.1e-06 | 0 | 0 | 0 | 0 | 1 | 0 | 1.1e-06 | 0.0007 | 0 |
| 15 | 2 | -5.71e-05 | 1.01e-05 | 6.72e-05 | 0 | 0 | 0 | 0 | 2 | 0 | 6.72e-05 | 0.0016 | 0 |
| 16 | 2 | -1.26e-05 | 0 | 1.96e-05 | 0 | 0 | 0 | 0 | 2 | 0 | 1.96e-05 | 0.000982 | 0 |
| 17 | 4 | -2.72e-05 | 0.000279 | 0.00046 | 0 | 0 | 0 | 0 | 4 | 0 | 0.00046 | 0.000384 | 0 |
| 18 | 2 | -4.19e-06 | 0.000184 | 0.000189 | 0 | 0 | 0 | 0 | 2 | 0 | 0.000189 | 0.000363 | 0 |
| 19 | 1 | -0.000194 | 0 | 0.000194 | 0 | 0 | 0 | 0 | 1 | 0 | 0.000194 | 0.00141 | 0 |
| 20 | 3 | -0.000174 | 0.000128 | 0.00044 | 0 | 0 | 0 | 0 | 3 | 0 | 0.00044 | 0.000847 | 0 |
| 21 | 4 | -0.000108 | 6.04e-06 | 0.000182 | 0 | 0 | 0 | 0 | 4 | 0 | 0.000182 | 0.000705 | 0 |
| 22 | 3 | -0.000331 | 0.000106 | 0.000467 | 0 | 0 | 0 | 0 | 3 | 0 | 0.000467 | 0.0011 | 0 |
| 23 | 1 | 0 | 3.07e-06 | 3.07e-06 | 0 | 0 | 0 | 0 | 1 | 0 | 3.07e-06 | 0.000322 | 0 |
| 24 | 1 | -4.02e-06 | 0 | 4.02e-06 | 0 | 0 | 0 | 0 | 1 | 0 | 4.02e-06 | 0.0018 | 0 |
| 25 | 2 | -0.000206 | 0 | 0.000248 | 0 | 0 | 0 | 0 | 2 | 0 | 0.000248 | 0.000601 | 0 |
| 27 | 1 | -0.000145 | 0 | 0.000145 | 0 | 0 | 0 | 0 | 1 | 0 | 0.000145 | 0.000756 | 0 |
| 28 | 3 | -6.48e-05 | 4.48e-05 | 0.000154 | 0 | 0 | 0 | 0 | 3 | 0 | 0.000163 | 0.000433 | 0 |
| 30 | 3 | -0.000192 | 0.000159 | 0.000353 | 0 | 0 | 0 | 0 | 3 | 0 | 0.000353 | 0.000813 | 0 |
| 31 | 1 | 0 | 2.01e-05 | 2.01e-05 | 0 | 0 | 0 | 0 | 1 | 0 | 2.01e-05 | 0.000333 | 0 |
| 32 | 1 | 0 | 1.98e-05 | 1.98e-05 | 0 | 0 | 0 | 0 | 1 | 0 | 1.98e-05 | 0.00166 | 0 |
| 33 | 2 | -0.000164 | 6.71e-05 | 0.000231 | 0 | 0 | 0 | 0 | 2 | 0 | 0.000231 | 0.000358 | 0 |
| 34 | 2 | -6.89e-05 | 4.41e-05 | 0.000113 | 0 | 0 | 0 | 0 | 2 | 0 | 0.000113 | 0.000521 | 0 |
| 35 | 4 | -8.99e-05 | 0 | 0.000154 | 0 | 0 | 0 | 0 | 4 | 0 | 0.000176 | 0.000926 | 0 |
| 36 | 2 | 0 | 0.000245 | 0.000269 | 0 | 0 | 0 | 0 | 2 | 0 | 0.000269 | 0.00251 | 0 |
| 37 | 1 | -8.6e-05 | 0 | 8.6e-05 | 0 | 0 | 0 | 0 | 1 | 0 | 8.6e-05 | 0.000997 | 0 |
| 38 | 2 | -0.000251 | 0.000152 | 0.000404 | 0 | 0 | 0 | 0 | 2 | 0 | 0.000404 | 0.000251 | 0 |
| 39 | 1 | -2.36e-07 | 0 | 2.36e-07 | 0 | 0 | 0 | 0 | 1 | 0 | 2.36e-07 | 0.000235 | 0 |
| 41 | 1 | -6.09e-05 | 0 | 6.09e-05 | 0 | 0 | 0 | 0 | 1 | 0 | 6.09e-05 | 0.000937 | 0 |
| 42 | 2 | -4.94e-05 | 0 | 9.86e-05 | 0 | 0 | 0 | 0 | 2 | 0 | 0.000106 | 0.00024 | 0 |
| 43 | 1 | -1.5e-05 | 0 | 1.5e-05 | 0 | 0 | 0 | 0 | 1 | 0 | 1.5e-05 | 0.0015 | 0 |
| 44 | 1 | 0 | 1.96e-05 | 1.96e-05 | 0 | 0 | 0 | 0 | 1 | 0 | 1.96e-05 | 0.000615 | 0 |
| 45 | 4 | -4.79e-05 | 0.000401 | 0.000743 | 0 | 0 | 0 | 0 | 4 | 0 | 0.000743 | 0.00136 | 0 |
| 46 | 3 | -0.000462 | 0 | 0.000909 | 0 | 0 | 0 | 0 | 3 | 0 | 0.000909 | 0.000776 | 0 |
| 47 | 1 | 0 | 2.57e-06 | 2.57e-06 | 0 | 0 | 0 | 0 | 1 | 0 | 2.57e-06 | 0.00207 | 0 |
| 48 | 4 | -0.00019 | 3.38e-06 | 0.000311 | 0 | 0 | 0 | 0 | 4 | 0 | 0.000311 | 0.000295 | 0 |
| 49 | 5 | -4.4e-05 | 3.08e-05 | 0.000143 | 0 | 0 | 0 | 0 | 5 | 0 | 0.000146 | 0.000235 | 0 |
| 51 | 2 | -4.69e-05 | 0 | 7.49e-05 | 0 | 0 | 0 | 0 | 2 | 0 | 7.49e-05 | 0.00102 | 0 |
| 52 | 2 | -3.96e-05 | 2.53e-05 | 6.49e-05 | 0 | 0 | 0 | 0 | 2 | 0 | 6.49e-05 | 0.000274 | 0 |
| 53 | 1 | 0 | 3.33e-05 | 3.33e-05 | 0 | 0 | 0 | 0 | 1 | 0 | 3.33e-05 | 0.000252 | 0 |
| 54 | 1 | -7.31e-05 | 0 | 7.31e-05 | 0 | 0 | 0 | 0 | 1 | 0 | 7.31e-05 | 0.000238 | 0 |
| 55 | 2 | -2.78e-05 | 0.000389 | 0.000417 | 0 | 0 | 0 | 0 | 2 | 0 | 0.000417 | 0.00178 | 0 |
| 56 | 7 | -0.000138 | 4.92e-06 | 0.000228 | 0 | 0 | 0 | 0 | 7 | 0 | 0.000228 | 0.000759 | 0 |

### The violating cells at the first violation step (step 2), first 1 by cell id

| cell | alpha^n | alpha new | has plane | outflow/V_k | inflow/V_k | faces: plane / noPlane uncut / noPlane CUT | donor = centroid cell | max prismWet | max flux out/V_k | max alpha^n of donors | face with max |one-pieces|/V_k |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 3286 | 0.000656 | -1.54e-05 | 1 | 0.000763 | 9.1e-05 | 4 / 0 / 0 (unproc. 0) | 3 of 4 | 0.00644 | 0.000763 | 0.155 | 1.54e-05 |

### The violating cells at the first step with a violation > 1e-3 (step 57), first 1 by cell id

| cell | alpha^n | alpha new | has plane | outflow/V_k | inflow/V_k | faces: plane / noPlane uncut / noPlane CUT | donor = centroid cell | max prismWet | max flux out/V_k | max alpha^n of donors | face with max |one-pieces|/V_k |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 3681 | -1.98e-16 | -0.00165 | 0 | 0.00165 | 0 | 4 / 0 / 0 (unproc. 0) | 0 of 4 | 0.217 | 0.0011 | 0.851 | 0.0011 |

### The violating cells at the last step (step 750), first 20 by cell id

| cell | alpha^n | alpha new | has plane | outflow/V_k | inflow/V_k | faces: plane / noPlane uncut / noPlane CUT | donor = centroid cell | max prismWet | max flux out/V_k | max alpha^n of donors | face with max |one-pieces|/V_k |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 35 | -2.47e-05 | -0.00258 | 0 | 0.00265 | 9.46e-05 | 3 / 1 / 0 (unproc. 0) | 1 of 4 | 0.997 | 0.00265 | 7.59e-05 | 0.00265 |
| 99 | -2.65e-18 | -0.0128 | 0 | 0.0132 | 0.000378 | 3 / 1 / 0 (unproc. 0) | 1 of 4 | 0.999 | 0.0131 | 0.000604 | 0.0131 |
| 163 | -0.000222 | -0.0253 | 0 | 0.026 | 0.000847 | 3 / 1 / 0 (unproc. 0) | 1 of 4 | 0.999 | 0.0256 | 0.00102 | 0.0256 |
| 227 | -2.31e-05 | -0.0348 | 0 | 0.0362 | 0.0015 | 3 / 1 / 0 (unproc. 0) | 1 of 4 | 0.999 | 0.0354 | 0.00146 | 0.0354 |
| 291 | -0.000426 | -0.0449 | 0 | 0.0468 | 0.00232 | 3 / 1 / 0 (unproc. 0) | 1 of 4 | 1 | 0.0453 | 0.00146 | 0.0453 |
| 351 | -1.67e-24 | -6.25e-10 | 0 | 6.25e-10 | 1e-25 | 0 / 0 / 0 (unproc. 4) | 0 of 4 | 0 | 6.25e-10 | 0 | n/a |
| 352 | -1.04e-08 | -3.84e-08 | 0 | 2.86e-08 | 6.28e-10 | 0 / 1 / 0 (unproc. 3) | 1 of 4 | 0 | 2.86e-08 | -1.37e-31 | 2.31e-08 |
| 355 | -0.00027 | -0.0614 | 0 | 0.0611 | 0 | 2 / 2 / 0 (unproc. 0) | 2 of 4 | 1 | 0.0588 | 0.00146 | 0.0588 |
| 374 | -1.45e-10 | -1.71e-10 | 0 | 2.88e-11 | 2.75e-12 | 0 / 0 / 0 (unproc. 4) | 0 of 4 | 0 | 2.54e-11 | 0 | n/a |
| 375 | -3.13e-10 | -3.68e-10 | 0 | 6.06e-11 | 5.04e-12 | 0 / 0 / 0 (unproc. 4) | 0 of 4 | 0 | 5.45e-11 | 0 | n/a |
| 376 | -6.92e-10 | -8.06e-10 | 0 | 1.23e-10 | 9.3e-12 | 0 / 0 / 0 (unproc. 4) | 0 of 4 | 0 | 1.18e-10 | 0 | n/a |
| 377 | -7.43e-10 | -8.42e-10 | 0 | 1.07e-10 | 8.12e-12 | 0 / 0 / 0 (unproc. 4) | 0 of 4 | 0 | 1.07e-10 | 0 | n/a |
| 414 | -2.79e-22 | -1.04e-07 | 0 | 1.04e-07 | 1.98e-23 | 0 / 0 / 0 (unproc. 4) | 0 of 4 | 0 | 1.04e-07 | 0 | n/a |
| 415 | -1.5e-06 | -5.52e-06 | 0 | 4.13e-06 | 1.04e-07 | 0 / 1 / 0 (unproc. 3) | 1 of 4 | 0 | 4.13e-06 | -1.5e-06 | 1.02e-09 |
| 416 | -5.93e-05 | -0.000105 | 0 | 4.96e-05 | 4.16e-06 | 0 / 1 / 0 (unproc. 3) | 1 of 4 | 0 | 4.96e-05 | -2.3e-32 | 1.58e-35 |
| 417 | -0.000714 | -0.000678 | 0 | 1.41e-05 | 4.96e-05 | 0 / 2 / 0 (unproc. 2) | 2 of 4 | 0 | 1.41e-05 | -4e-32 | 1e-06 |
| 418 | -0.000204 | -0.000193 | 0 | 2.96e-06 | 1.41e-05 | 0 / 2 / 0 (unproc. 2) | 2 of 4 | 0 | 2.96e-06 | -7.86e-19 | 4.74e-07 |
| 419 | -4.34e-05 | -4.04e-05 | 0 | 1.67e-19 | 2.96e-06 | 0 / 2 / 0 (unproc. 2) | 2 of 4 | 0 | 1.67e-19 | -1.11e-16 | 1.39e-07 |
| 436 | -3.74e-10 | -4.6e-10 | 0 | 9.69e-11 | 1.14e-11 | 0 / 0 / 0 (unproc. 4) | 0 of 4 | 0 | 6.96e-11 | 0 | n/a |
| 437 | -1.48e-09 | -1.77e-09 | 0 | 3.26e-10 | 3.98e-11 | 0 / 0 / 0 (unproc. 4) | 0 of 4 | 0 | 2.76e-10 | 0 | n/a |

### Limiter per step (last 5 steps)

- step 745: passes 3, faces changed 62, cells touched 80, volume moved 1.27 cell volumes
- step 746: passes 3, faces changed 63, cells touched 86, volume moved 1.21 cell volumes
- step 747: passes 3, faces changed 70, cells touched 98, volume moved 1.21 cell volumes
- step 748: passes 3, faces changed 71, cells touched 100, volume moved 1.2 cell volumes
- step 749: passes 3, faces changed 70, cells touched 96, volume moved 1.19 cell volumes

### boundFlux division operands (denominator below 1e-12)

- DIAG boundFlux division by dVftot = 0 (davofSweptAdvectionTemplates.C, fluidToPassThroughFace = mag(fluidToPassOn)*mag(phi[fi]*dt)/dVftot): cell 35 iter 1 alpha -7.2054731509015428e-05 alphaOvershoot -2.1083860500335302e-06 fluidToPassOn -8.0428544999448022e-12 Vi 3.814697265625e-06 dt 0.0024038461540000001 faces to pass through 1: [face 8227 phi 0 phi*dt 0 dVfmax 5.6060001186124528e-13 dVf -5.6060001186124528e-13 numerator 0]
- (1 such lines)

### onePlane vs pieces flux difference |flux_onePlane - flux_pieces|/V_k over the faces of interface cells

| step | faces | L1 (mean) | Linf | count > 0.01 | count > 0.1 | of those: donor plane / no plane cut / no plane uncut |
|---|---|---|---|---|---|---|
| 1 | 228 | 1.13e-06 | 0.000195 | 0 | 0 | 0 / 0 / 0 |
| 2 | 228 | 1.83e-05 | 0.00101 | 0 | 0 | 0 / 0 / 0 |
| 47 | 240 | 1.58e-05 | 0.00207 | 0 | 0 | 0 / 0 / 0 |
| 57 | 238 | 1.42e-05 | 0.0011 | 0 | 0 | 0 / 0 / 0 |
| 740 | 916 | 0.000634 | 0.1 | 15 | 1 | 1 / 0 / 0 |
| 750 | 902 | 0.000543 | 0.0956 | 12 | 0 | 0 / 0 / 0 |

