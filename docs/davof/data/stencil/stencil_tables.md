
### Test 1a: the plane gate (72 cases per variant: N 40/80, Co 0.2/0.5, both arms)

| variant | completed | gate < 1e-11 | max L1(alphaf) | max L1(n) | min alpha | max excess | max |dev| % | cells dropped | fit fallbacks | = production | step / production |
|---|---|---|---|---|---|---|---|---|---|---|---|
| production | 8/8 | 8/8 | 2.842e-13 | 1.597e-13 | -4.776e-14 | 1.898e-13 | 9.752e+01 | 0 | 0 | 8/8 | 1.00 |
| twoRings | 8/8 | 8/8 | 2.833e-13 | 8.108e-14 | -4.433e-14 | 1.898e-13 | 9.752e+01 | 0 | 0 | 0/8 | 1.07 |
| kernelGaussian | 8/8 | 8/8 | 2.854e-13 | 2.803e-13 | -4.160e-14 | 1.898e-13 | 9.752e+01 | 0 | 0 | 0/8 | 0.98 |
| kernelWendland | 8/8 | 8/8 | 2.894e-13 | 6.053e-13 | -4.181e-14 | 1.898e-13 | 9.752e+01 | 0 | 0 | 0/8 | 0.99 |
| angle30 | 8/8 | 8/8 | 2.842e-13 | 1.597e-13 | -4.776e-14 | 1.898e-13 | 9.752e+01 | 0 | 0 | 8/8 | 1.02 |
| angle60 | 8/8 | 8/8 | 2.842e-13 | 1.597e-13 | -4.776e-14 | 1.898e-13 | 9.752e+01 | 0 | 0 | 8/8 | 0.97 |
| planeFitMean | 8/8 | 8/8 | 2.896e-13 | 7.837e-14 | -4.220e-14 | 1.898e-13 | 9.752e+01 | 0 | 726 | 0/8 | 1.02 |
| planeFitOwn | 8/8 | 7/8 | 1.167e-08 | 3.731e-08 | -8.122e-13 | 2.631e-13 | 9.752e+01 | 0 | 122598 | 0/8 | 1.00 |
| ringQuadric | 8/8 | 8/8 | 2.838e-13 | 1.802e-13 | -4.550e-14 | 1.898e-13 | 9.752e+01 | 0 | 3285 | 0/8 | 1.02 |

G-0, the production rows against the port's record (port_enrightKeys_plane_leia, recommended): 2 of 2 shared cells identical to the digit at every step.

### Test 1b: the seeded perturbation (1e-8 on the partial faces): g_early / g_late / L1(alphaf) after step 1 / at the end

| variant | N 40 Co 0.2 normal | N 40 Co 0.2 oblique | N 40 Co 0.5 normal | N 40 Co 0.5 oblique | N 80 Co 0.2 normal | N 80 Co 0.2 oblique | N 80 Co 0.5 normal | N 80 Co 0.5 oblique |
|---|---|---|---|---|---|---|---|---|
| production | 0.68 / 0.98 / 1.9e-10 / 6.8e-12 | 0.71 / 0.95 / 2.0e-10 / 2.3e-12 | 0.75 / 0.97 / 1.5e-10 / 1.8e-11 | 0.75 / 0.91 / 1.6e-10 / 9.9e-12 | 0.68 / 0.99 / 1.8e-10 / 4.1e-12 | 0.70 / 0.97 / 1.9e-10 / 2.9e-13 | 0.73 / 0.97 / 1.5e-10 / 8.8e-12 | 0.74 / 0.91 / 1.6e-10 / 1.3e-12 |
| twoRings | 0.67 / 0.98 / 9.2e-11 / 2.4e-12 | 0.69 / 0.93 / 9.7e-11 / 5.6e-13 | 0.70 / 0.96 / 7.5e-11 / 5.7e-12 | 0.73 / 0.89 / 7.8e-11 / 3.2e-12 | 0.68 / 0.99 / 8.5e-11 / 1.5e-12 | 0.68 / 0.99 / 8.8e-11 / 2.4e-13 | 0.70 / 0.97 / 7.1e-11 / 2.9e-12 | 0.72 / 0.91 / 7.2e-11 / 3.7e-13 |
| kernelGaussian | 0.72 / 0.99 / 3.0e-10 / 1.6e-11 | 0.73 / 0.96 / 3.2e-10 / 7.7e-12 | 0.75 / 0.97 / 2.5e-10 / 3.4e-11 | 0.76 / 0.94 / 2.6e-10 / 2.3e-11 | 0.72 / 0.99 / 2.8e-10 / 8.0e-12 | 0.72 / 0.97 / 3.1e-10 / 1.2e-12 | 0.75 / 0.98 / 2.3e-10 / 1.6e-11 | 0.76 / 0.94 / 2.5e-10 / 6.0e-12 |
| kernelWendland | 0.76 / 0.99 / 5.2e-10 / 2.8e-11 | 0.74 / 0.98 / 5.6e-10 / 2.0e-11 | 0.76 / 0.97 / 4.3e-10 / 5.9e-11 | 0.77 / 0.95 / 4.6e-10 / 4.7e-11 | 0.75 / 0.99 / 5.1e-10 / 1.4e-11 | 0.73 / 0.98 / 5.7e-10 / 7.8e-12 | 0.75 / 0.98 / 4.2e-10 / 2.9e-11 | 0.76 / 0.96 / 4.7e-10 / 2.3e-11 |
| angle30 | 0.68 / 0.98 / 1.9e-10 / 6.8e-12 | 0.71 / 0.95 / 2.0e-10 / 2.3e-12 | 0.75 / 0.97 / 1.5e-10 / 1.8e-11 | 0.75 / 0.91 / 1.6e-10 / 9.9e-12 | 0.68 / 0.99 / 1.8e-10 / 4.1e-12 | 0.70 / 0.97 / 1.9e-10 / 2.9e-13 | 0.73 / 0.97 / 1.5e-10 / 8.8e-12 | 0.74 / 0.91 / 1.6e-10 / 1.3e-12 |
| angle60 | 0.68 / 0.98 / 1.9e-10 / 6.8e-12 | 0.71 / 0.95 / 2.0e-10 / 2.3e-12 | 0.75 / 0.97 / 1.5e-10 / 1.8e-11 | 0.75 / 0.91 / 1.6e-10 / 9.9e-12 | 0.68 / 0.99 / 1.8e-10 / 4.1e-12 | 0.70 / 0.97 / 1.9e-10 / 2.9e-13 | 0.73 / 0.97 / 1.5e-10 / 8.8e-12 | 0.74 / 0.91 / 1.6e-10 / 1.3e-12 |
| planeFitMean | 0.84 / 1.00 / 1.1e-13 / 7.2e-14 | 0.72 / 1.00 / 9.5e-13 / 1.2e-13 | 1.01 / 1.02 / 4.8e-14 / 6.2e-14 | 0.68 / 1.00 / 6.6e-13 / 8.3e-14 | 1.01 / 1.01 / 1.9e-13 / 3.7e-13 | 1.00 / 1.00 / 2.1e-13 / 2.4e-13 | 1.01 / 1.01 / 1.9e-13 / 2.7e-13 | 1.01 / 1.00 / 2.1e-13 / 2.3e-13 |
| planeFitOwn | 0.67 / 0.97 / 2.6e-11 / 2.3e-13 | 0.79 / 0.98 / 3.5e-11 / 2.9e-12 | 0.68 / 0.92 / 1.6e-11 / 3.9e-13 | 0.86 / 1.00 / 3.1e-11 / 8.5e-12 | 0.71 / 1.01 / 1.2e-11 / 3.9e-13 | 0.80 / 1.00 / 1.6e-11 / 9.0e-09 | 0.72 / 1.00 / 6.7e-12 / 3.0e-13 | 0.82 / 1.09 / 1.4e-11 / 4.8e-11 |
| ringQuadric | 0.71 / 0.99 / 1.3e-10 / 4.8e-12 | 0.72 / 0.99 / 1.5e-10 / 1.0e-11 | 0.75 / 0.97 / 1.0e-10 / 1.2e-11 | 0.77 / 0.97 / 1.2e-10 / 1.7e-11 | 0.68 / 0.99 / 1.2e-10 / 2.8e-12 | 0.71 / 0.99 / 1.3e-10 / 1.4e-12 | 0.73 / 0.98 / 9.3e-11 / 5.8e-12 | 0.75 / 0.95 / 1.1e-10 / 2.8e-12 |

L1(alphaf) at step 1 (median over the cases): 1.509e-10; steps per cell: (40, '0.2', 'normal'): 50, (40, '0.2', 'oblique'): 50, (40, '0.5', 'normal'): 20, (40, '0.5', 'oblique'): 20, (80, '0.2', 'normal'): 100, (80, '0.2', 'oblique'): 100, (80, '0.5', 'normal'): 40, (80, '0.5', 'oblique'): 40.

### Test 2: the translated sphere (Co 0.2, one period; N 32/64/128 on 1/8/64 ranks; three repetitions)

| variant | completed | L1(n) N32 | L1(n) N64 | L1(n) N128 | p(n) | L1(alphaf) N128 | p(alphaf) | L1(x) N128 | p(x) | E_shape N128 | p(E) | min alpha | max excess | max |dev| % | wall N32 s (spread) | wall N64 s | wall N128 s | step/prod N32 | N64 | N128 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| production | 9/9 | 4.107e-02 | 2.239e-02 | 1.091e-02 | 0.96 | 5.560e-03 | 0.83 | 1.025e-04 | 1.86 | 2.443e-05 | 1.84 | -7.329e-13 | 7.905e-14 | 3.997e-13 | 34 (0.21) | 80 (0.45) | 814 (0.01) | 1.00 | 1.00 | 1.00 |
| twoRings | 9/9 | 5.635e-02 | 2.625e-02 | 1.092e-02 | 1.18 | 4.609e-03 | 0.92 | 8.771e-05 | 1.88 | 1.968e-05 | 1.84 | -2.871e-13 | 6.395e-14 | 4.663e-13 | 29 (0.07) | 99 (0.24) | 817 (0.01) | 0.81 | 1.29 | 1.09 |
| kernelGaussian | 9/9 | 3.788e-02 | 2.110e-02 | 9.878e-03 | 0.97 | 6.500e-03 | 0.83 | 1.153e-04 | 1.89 | 2.814e-05 | 1.85 | -5.947e-13 | 7.883e-14 | 4.774e-13 | 31 (0.13) | 97 (0.10) | 815 (0.00) | 0.86 | 1.21 | 1.01 |
| kernelWendland | 9/9 | 6.294e-02 | 3.370e-02 | 1.637e-02 | 0.97 | 8.302e-03 | 0.78 | 1.379e-04 | 1.88 | 3.386e-05 | 1.84 | -7.795e-13 | 9.170e-14 | 5.440e-13 | 31 (0.19) | 105 (0.21) | 818 (0.34) | 0.86 | 1.17 | 1.01 |
| angle30 | 9/9 | 4.194e-02 | 2.242e-02 | 1.092e-02 | 0.97 | 5.606e-03 | 0.85 | 1.029e-04 | 1.88 | 2.454e-05 | 1.86 | -7.062e-13 | 7.017e-14 | 4.330e-13 | 28 (0.07) | 88 (0.17) | 814 (0.34) | 0.81 | 1.13 | 1.01 |
| angle60 | 9/9 | 4.112e-02 | 2.242e-02 | 1.091e-02 | 0.96 | 5.573e-03 | 0.83 | 1.026e-04 | 1.86 | 2.446e-05 | 1.84 | -7.289e-13 | 8.571e-14 | 4.330e-13 | 28 (0.29) | 93 (0.08) | 676 (0.22) | 0.81 | 1.17 | 1.01 |
| planeFitMean | 9/9 | 4.468e-02 | 2.464e-02 | 1.170e-02 | 0.97 | 6.101e-03 | 0.85 | 1.112e-04 | 1.86 | 2.704e-05 | 1.87 | -5.855e-13 | 8.171e-14 | 4.885e-13 | 29 (0.10) | 92 (0.16) | 815 (0.17) | 0.76 | 1.17 | 1.02 |
| planeFitOwn | 9/9 | 4.065e-02 | 2.231e-02 | 1.096e-02 | 0.95 | 5.665e-03 | 0.82 | 1.040e-04 | 1.85 | 2.487e-05 | 1.83 | -6.930e-13 | 1.115e-13 | 4.441e-13 | 29 (0.14) | 102 (0.24) | 547 (0.48) | 0.76 | 1.29 | 1.02 |
| ringQuadric | 9/9 | 2.440e-02 | 1.150e-02 | 5.269e-03 | 1.11 | 5.619e-03 | 0.85 | 1.049e-04 | 1.89 | 2.500e-05 | 1.87 | -4.412e-13 | 7.683e-14 | 4.663e-13 | 31 (0.23) | 90 (0.06) | 573 (0.47) | 0.88 | 1.17 | 1.02 |

Repetitions and identity: production N 32: reps identical, = production True; production N 64: reps identical, = production True; production N 128: reps identical, = production True; twoRings N 32: reps identical, = production False; twoRings N 64: reps identical, = production False; twoRings N 128: reps identical, = production False; kernelGaussian N 32: reps identical, = production False; kernelGaussian N 64: reps identical, = production False; kernelGaussian N 128: reps identical, = production False; kernelWendland N 32: reps identical, = production False; kernelWendland N 64: reps identical, = production False; kernelWendland N 128: reps identical, = production False; angle30 N 32: reps identical, = production False; angle30 N 64: reps identical, = production False; angle30 N 128: reps identical, = production False; angle60 N 32: reps identical, = production False; angle60 N 64: reps identical, = production False; angle60 N 128: reps identical, = production False; planeFitMean N 32: reps identical, = production False; planeFitMean N 64: reps identical, = production False; planeFitMean N 128: reps identical, = production False; planeFitOwn N 32: reps identical, = production False; planeFitOwn N 64: reps identical, = production False; planeFitOwn N 128: reps identical, = production False; ringQuadric N 32: reps identical, = production False; ringQuadric N 64: reps identical, = production False; ringQuadric N 128: reps identical, = production False.

### Test 3a: the rotation (one turn, Co 0.2, N 32/64 serial, three repetitions)

| variant | N | completed | min alpha | max excess | max |dev| % | E_shape | L1(x) | L1(n) | L1(alphaf) | interface cells end | max | cells dropped (max/step) | fit fallbacks | wall s (spread) | step/prod |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| production | 32 | 3/3 | -3.499e-13 | 1.354e-14 | 4.996e-13 | 7.612e-04 | 2.811e-03 | 4.944e-02 | 4.164e-02 | 414 | 447 | 0 | 0 | 109 (0.07) | 1.00 |
| production | 64 | 3/3 | -7.723e-13 | 2.147e-13 | 1.998e-12 | 2.114e-04 | 7.584e-04 | 2.276e-02 | 2.323e-02 | 1739 | 1760 | 0 | 0 | 1094 (0.03) | 1.00 |
| twoRings | 32 | 3/3 | -1.238e-13 | 1.354e-14 | 3.775e-13 | 6.953e-04 | 2.484e-03 | 6.198e-02 | 3.892e-02 | 413 | 452 | 0 | 0 | 112 (0.06) | 1.08 |
| twoRings | 64 | 3/3 | -1.612e-12 | 5.256e-13 | 1.699e-12 | 2.004e-04 | 7.139e-04 | 2.786e-02 | 2.223e-02 | 1742 | 1760 | 0 | 0 | 1100 (0.12) | 1.02 |
| kernelGaussian | 32 | 3/3 | -1.203e-14 | 1.354e-14 | 4.441e-13 | 7.992e-04 | 3.082e-03 | 4.425e-02 | 4.410e-02 | 422 | 449 | 0 | 0 | 111 (0.08) | 1.00 |
| kernelGaussian | 64 | 3/3 | -8.665e-13 | 5.329e-14 | 7.327e-13 | 2.176e-04 | 7.881e-04 | 1.950e-02 | 2.400e-02 | 1744 | 1765 | 0 | 0 | 1096 (0.03) | 1.02 |
| kernelWendland | 32 | 3/3 | -8.417e-13 | 1.354e-14 | 9.992e-13 | 8.374e-04 | 3.234e-03 | 6.632e-02 | 4.608e-02 | 423 | 449 | 0 | 0 | 114 (0.15) | 1.04 |
| kernelWendland | 64 | 3/3 | -8.993e-13 | 1.983e-13 | 9.104e-13 | 2.209e-04 | 8.131e-04 | 2.835e-02 | 2.466e-02 | 1744 | 1767 | 0 | 0 | 1093 (0.09) | 1.02 |
| angle30 | 32 | 3/3 | -4.756e-13 | 1.421e-14 | 4.663e-13 | 7.679e-04 | 2.841e-03 | 5.031e-02 | 4.211e-02 | 416 | 447 | 259 | 0 | 106 (0.05) | 1.00 |
| angle30 | 64 | 3/3 | -6.812e-13 | 2.052e-13 | 2.076e-12 | 2.119e-04 | 7.600e-04 | 2.281e-02 | 2.328e-02 | 1739 | 1760 | 305 | 0 | 1140 (0.10) | 1.08 |
| angle60 | 32 | 3/3 | -1.082e-14 | 3.264e-14 | 3.997e-13 | 7.616e-04 | 2.813e-03 | 4.947e-02 | 4.167e-02 | 414 | 447 | 22 | 0 | 106 (0.05) | 0.96 |
| angle60 | 64 | 3/3 | -6.474e-13 | 2.185e-13 | 1.943e-12 | 2.115e-04 | 7.588e-04 | 2.277e-02 | 2.324e-02 | 1739 | 1760 | 112 | 0 | 1131 (0.17) | 1.08 |
| planeFitMean | 32 | 3/3 | -6.720e-14 | 1.510e-14 | 4.330e-13 | 7.674e-04 | 2.822e-03 | 4.813e-02 | 4.208e-02 | 414 | 450 | 0 | 0 | 109 (0.28) | 1.04 |
| planeFitMean | 64 | 3/3 | -9.433e-13 | 5.329e-14 | 1.976e-12 | 2.141e-04 | 7.683e-04 | 2.359e-02 | 2.355e-02 | 1737 | 1760 | 0 | 0 | 1144 (0.16) | 1.08 |
| planeFitOwn | 32 | 3/3 | -2.788e-13 | 1.443e-14 | 3.886e-13 | 7.611e-04 | 2.811e-03 | 4.912e-02 | 4.163e-02 | 414 | 447 | 0 | 514212 | 114 (0.18) | 1.08 |
| planeFitOwn | 64 | 3/3 | -7.713e-13 | 2.185e-13 | 1.787e-12 | 2.115e-04 | 7.587e-04 | 2.275e-02 | 2.323e-02 | 1739 | 1760 | 0 | 4128111 | 1122 (0.02) | 1.04 |
| ringQuadric | 32 | 3/3 | -3.099e-13 | 1.354e-14 | 3.886e-13 | 7.831e-04 | 2.932e-03 | 3.466e-02 | 4.360e-02 | 422 | 448 | 0 | 0 | 129 (0.19) | 1.23 |
| ringQuadric | 64 | 3/3 | -9.349e-13 | 5.329e-14 | 8.660e-13 | 2.146e-04 | 7.758e-04 | 1.197e-02 | 2.376e-02 | 1742 | 1760 | 0 | 9 | 1097 (0.06) | 1.02 |

Spearman rank correlation at N 32 over 9 variants: L1(n) vs E_shape -0.15, L1(n) vs L1(x) -0.15.

Spearman rank correlation at N 64 over 9 variants: L1(n) vs E_shape -0.07, L1(n) vs L1(x) -0.07.

### Test 3b: the 3D deformation (N 32, OpenFOAM Co 0.5, T = 3, three repetitions)

| variant | N | completed | min alpha | max excess | max |dev| % | E_shape | L1(x) | L1(n) | L1(alphaf) | interface cells end | max | cells dropped (max/step) | fit fallbacks | wall s (spread) | step/prod |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| production | 32 | 3/3 | -1.342e-07 | 1.003e-07 | 2.598e-12 | 9.315e-03 | 8.469e-02 | 7.216e-01 | 2.624e-01 | 688 | 1688 | 0 | 0 | 587 (0.05) | 1.00 |
| twoRings | 32 | 3/3 | -1.414e-03 | 5.996e-07 | 2.043e-12 | 1.276e-02 | 1.089e-01 | 7.218e-01 | 3.102e-01 | 741 | 1523 | 0 | 0 | 611 (0.06) | 1.01 |
| kernelGaussian | 32 | 3/3 | -1.992e-07 | 8.094e-09 | 4.774e-12 | 7.008e-03 | 6.127e-02 | 7.034e-01 | 2.122e-01 | 713 | 1458 | 0 | 0 | 575 (0.13) | 1.07 |
| kernelWendland | 32 | 3/3 | -2.173e-06 | 1.815e-09 | 4.530e-12 | 6.197e-03 | 9.433e-02 | 8.365e-01 | 1.716e-01 | 856 | 1561 | 0 | 0 | 522 (0.03) | 1.01 |
| angle30 | 32 | 3/3 | -7.397e-08 | 4.261e-13 | 2.998e-12 | 7.893e-03 | 1.210e-01 | 7.903e-01 | 2.080e-01 | 780 | 1852 | 17431 | 0 | 1288 (0.05) | 2.50 |
| angle60 | 32 | 3/3 | -9.292e-08 | 2.243e-14 | 3.397e-12 | 6.168e-03 | 2.189e-01 | 1.022e+00 | 1.573e-01 | 981 | 2329 | 15225 | 0 | 1179 (0.12) | 2.22 |
| planeFitMean | 32 | 3/3 | -3.872e-04 | 7.675e-09 | 2.820e-12 | 1.017e-02 | 3.083e-01 | 1.439e+00 | 1.071e-01 | 3698 | 3783 | 0 | 3840 | 1566 (0.15) | 2.47 |
| planeFitOwn | 32 | 3/3 | -1.322e-03 | 2.638e-07 | 3.286e-12 | 8.536e-03 | 8.876e-02 | 6.683e-01 | 2.526e-01 | 654 | 1501 | 0 | 1076682 | 595 (0.03) | 0.99 |
| ringQuadric | 32 | 3/3 | -7.937e-03 | 2.145e-08 | 2.887e-12 | 8.547e-03 | 7.347e-02 | 6.822e-01 | 2.406e-01 | 682 | 1778 | 0 | 42201 | 1250 (0.05) | 2.18 |

Spearman rank correlation at N 32 over 9 variants: L1(n) vs E_shape -0.13, L1(n) vs L1(x) 0.83.
