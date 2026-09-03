# Route B, preliminary 2 — does m stay bounded as the grid refines?

Left-endpoint discrete kernel with smallest lag delta (the paper's Definition "Exact convolution"),
so the lifts are optimised for L2(delta,T).  h = -0.4, H = 0.1, T = 1.0, kappa = 0.5.
Var[V_T] continuous = 2.254595.

## Phase 1 — the threshold to beat

This is the error the scheme already commits with the TRUE kernel, purely from cutting
time into n pieces.  A lift whose added error sits below it is not the binding constraint.

| n | delta | Var_disc(true K) | Var ratio (target 1) | threshold = relative error | delta^{H/2} for comparison |
|---|---|---|---|---|---|
| 64 | 1.562e-02 | 1.387137 | 0.615249 | **3.8475e-01** | 0.8123 |
| 128 | 7.812e-03 | 1.498128 | 0.664478 | **3.3552e-01** | 0.7846 |
| 256 | 3.906e-03 | 1.595401 | 0.707622 | **2.9238e-01** | 0.7579 |
| 512 | 1.953e-03 | 1.680407 | 0.745326 | **2.5467e-01** | 0.7320 |
| 1024 | 9.766e-04 | 1.754572 | 0.778221 | **2.2178e-01** | 0.7071 |
| 2048 | 4.883e-04 | 1.819218 | 0.806894 | **1.9311e-01** | 0.6830 |
| 4096 | 2.441e-04 | 1.875537 | 0.831873 | **1.6813e-01** | 0.6598 |
| 8192 | 1.221e-04 | 1.924585 | 0.853628 | **1.4637e-01** | 0.6373 |
| 16384 | 6.104e-05 | 1.967295 | 0.872571 | **1.2743e-01** | 0.6156 |

The threshold decays in n with log-log slope **-0.1994**, against the predicted -2H = -0.2.

## Phase 2 — m*(n), the smallest usable number of factors

| n | threshold | m* (Abi Jaber–El Euch) | m* (achievable floor) |
|---|---|---|---|
| 64 | 3.8475e-01 | 2 | **1** |
| 128 | 3.3552e-01 | 2 | **1** |
| 256 | 2.9238e-01 | 3 | **1** |
| 512 | 2.5467e-01 | 3 | **1** |
| 1024 | 2.2178e-01 | 4 | **1** |
| 2048 | 1.9311e-01 | 5 | **1** |
| 4096 | 1.6813e-01 | 5 | **2** |
| 8192 | 1.4637e-01 | 6 | **2** |
| 16384 | 1.2743e-01 | 8 | **2** |

### The full error table (achievable floor), added variance error

| n | m=1 | m=2 | m=3 | m=4 | m=5 | m=6 | m=8 | m=10 | m=12 |
|---|---|---|---|---|---|---|---|---|---|
| 64 | 9.19e-02 | 2.06e-02 | 5.12e-03 | 1.14e-03 | 2.19e-04 | 4.45e-05 | 1.24e-05 | 4.76e-06 | 9.06e-07 |
| 128 | 1.04e-01 | 2.39e-02 | 6.98e-03 | 1.93e-03 | 4.80e-04 | 1.05e-04 | 8.37e-06 | 3.39e-06 | 1.50e-06 |
| 256 | 1.19e-01 | 2.71e-02 | 8.69e-03 | 2.80e-03 | 8.43e-04 | 2.25e-04 | 4.72e-05 | 6.32e-06 | 2.48e-06 |
| 512 | 1.35e-01 | 3.04e-02 | 1.02e-02 | 3.69e-03 | 1.28e-03 | 4.12e-04 | 1.01e-04 | 2.42e-05 | 6.94e-06 |
| 1024 | 1.51e-01 | 3.41e-02 | 1.16e-02 | 4.53e-03 | 1.74e-03 | 6.41e-04 | 1.84e-04 | 3.53e-05 | 1.03e-05 |
| 2048 | 1.68e-01 | 3.84e-02 | 1.30e-02 | 5.30e-03 | 2.22e-03 | 8.99e-04 | 1.26e-04 | 7.79e-05 | 1.15e-05 |
| 4096 | 1.83e-01 | 4.34e-02 | 1.44e-02 | 6.00e-03 | 2.67e-03 | 1.17e-03 | 3.15e-04 | 7.66e-05 | 2.58e-05 |
| 8192 | 1.98e-01 | 4.90e-02 | 1.59e-02 | 6.66e-03 | 3.08e-03 | 1.44e-03 | 2.90e-04 | 7.56e-05 | 4.46e-05 |
| 16384 | 2.11e-01 | 5.51e-02 | 1.76e-02 | 7.30e-03 | 3.47e-03 | 1.70e-03 | 3.89e-04 | 1.17e-04 | 3.09e-05 |

### And the kernel L2(delta,T) error behind it

| n | m=1 | m=2 | m=3 | m=4 | m=5 | m=6 | m=8 | m=10 | m=12 |
|---|---|---|---|---|---|---|---|---|---|
| 64 | 2.14e-01 | 3.64e-02 | 6.02e-03 | 9.95e-04 | 1.65e-04 | 2.78e-05 | 9.71e-06 | 3.69e-06 | 6.20e-07 |
| 128 | 2.77e-01 | 5.78e-02 | 1.17e-02 | 2.36e-03 | 4.77e-04 | 9.67e-05 | 8.24e-06 | 3.25e-06 | 1.14e-06 |
| 256 | 3.39e-01 | 8.35e-02 | 1.99e-02 | 4.70e-03 | 1.11e-03 | 2.69e-04 | 5.32e-05 | 7.14e-06 | 2.54e-06 |
| 512 | 3.98e-01 | 1.13e-01 | 3.07e-02 | 8.26e-03 | 2.23e-03 | 6.02e-04 | 1.30e-04 | 2.89e-05 | 8.78e-06 |
| 1024 | 4.53e-01 | 1.44e-01 | 4.39e-02 | 1.32e-02 | 3.96e-03 | 1.19e-03 | 2.84e-04 | 5.35e-05 | 1.52e-05 |
| 2048 | 5.01e-01 | 1.77e-01 | 5.93e-02 | 1.96e-02 | 6.45e-03 | 2.13e-03 | 2.32e-04 | 1.17e-04 | 2.05e-05 |
| 4096 | 5.45e-01 | 2.11e-01 | 7.65e-02 | 2.74e-02 | 9.76e-03 | 3.48e-03 | 7.39e-04 | 1.55e-04 | 5.02e-05 |
| 8192 | 5.83e-01 | 2.44e-01 | 9.53e-02 | 3.66e-02 | 1.40e-02 | 5.33e-03 | 7.80e-04 | 1.79e-04 | 9.48e-05 |
| 16384 | 6.16e-01 | 2.77e-01 | 1.15e-01 | 4.70e-02 | 1.90e-02 | 7.72e-03 | 1.27e-03 | 3.24e-04 | 8.14e-05 |

## Phase 3 — confirmation with the whole covariance matrix

The terminal variance is one number and could flatter the lift; this repeats the crossing
with the relative Frobenius norm of the full n x n covariance matrix.

| n | threshold (Frobenius) | m* (floor) |
|---|---|---|
| 64 | 1.9630e-01 | **1** |
| 128 | 1.3121e-01 | **1** |
| 256 | 8.6778e-02 | **2** |
| 512 | 5.7108e-02 | **2** |

## Phase 4 — growth of m* and what it costs

Fit: `m*(n) ~ a*log(n) + b` with a = **0.2164**, b = -0.1667.
Extrapolated: m* = 2.8 at n = 10^6, m* = 9.8 at n = 10^20.

A logarithmic fit is the claim being tested: if m* grows like log n, Route B is usable
at every n; if it grows like a power of n, it is not.

| target accuracy eps | n needed (rate delta^{H/2}) | m needed | cost eps^-q with Route A' inside |
|---|---|---|---|
| 0.1 | 1e+20 | 9.8 | 221.0 |
| 0.01 | 1e+40 | 19.8 | 421.0 |
| 0.001 | 1e+60 | 29.7 | 621.0 |

Route A' alone (inconsistent, so not a competitor): eps^-41.0.

## Compute time

| phase | wall time | seconds |
|---|---|---|
| 1. the threshold: the scheme's own error with the TRUE kernel | 0s | 0.0 |
| 2. m*(n): smallest m whose added error is below the threshold | 57s | 56.83 |
| 3. confirmation with the full covariance matrix (Frobenius) | 36s | 35.98 |
| **total** | **1m32s** | **92.81** |

Machine: Darwin arm64, python 3.9.6.
