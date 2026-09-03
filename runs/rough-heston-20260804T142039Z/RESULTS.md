# Rough Heston: an independent reference, and what it costs

Model parameters {'V0': 0.09, 'theta': 0.09, 'lam': 1.5, 'nu': 0.3, 'rho': -0.7}, T = 1.0, S0 = K = 100.
Fourier price by the Lewis representation on the fractional-Riccati characteristic function
(El Euch--Rosenbaum), solved by the Diethelm--Ford--Freed predictor--corrector.

## Phase 1 --- exact controls

### (a) nu = lambda = 0 collapses to Black--Scholes

| V0 | Fourier | exact | absolute error |
|---|---|---|---|
| 0.04 | 7.9655674554 | 7.9655674554 | 2.46e-11 |
| 0.09 | 11.9235384741 | 11.9235384740 | 2.44e-11 |
| 0.25 | 19.7412651366 | 19.7412651366 | 2.39e-11 |

### (b) phi(0) = 1 and phi(-i) = 1 (the martingale property)

| H | error on phi(0) | error on phi(-i) |
|---|---|---|
| 0.05 | 0.00e+00 | 0.00e+00 |
| 0.1 | 0.00e+00 | 0.00e+00 |
| 0.3 | 0.00e+00 | 0.00e+00 |
| 0.45 | 0.00e+00 | 0.00e+00 |

### (c) nu = 0, lambda > 0: deterministic fractional variance

This is the control that matters: the variance is non-trivial and fractional, only its
randomness is switched off, so the whole fractional machinery is exercised at every H.

| H | exact | err, 200 steps | 400 | 800 | observed order | predicted min(2, 1+alpha) |
|---|---|---|---|---|---|---|
| 0.05 | 14.20081114 | 3.73e-04 | 1.24e-04 | 4.05e-05 | **1.6** | 1.55 |
| 0.1 | 14.17909844 | 2.68e-04 | 8.64e-05 | 2.75e-05 | **1.64** | 1.6 |
| 0.2 | 14.13301716 | 1.40e-04 | 4.23e-05 | 1.26e-05 | **1.73** | 1.7 |
| 0.3 | 14.08236181 | 7.31e-05 | 2.07e-05 | 5.82e-06 | **1.83** | 1.8 |
| 0.45 | 13.99477668 | 2.79e-05 | 7.17e-06 | 1.82e-06 | **1.97** | 1.95 |

## Phase 2 --- the reference's own accuracy (H = 0.10)

Taken as truth: the Fourier price at 3200 steps, **11.45567447**.

| Riccati steps | price | error | seconds |
|---|---|---|---|
| 100 | 11.45544129 | 2.33e-04 | 0.004 |
| 200 | 11.45560323 | 7.12e-05 | 0.005 |
| 400 | 11.45565253 | 2.19e-05 | 0.011 |
| 800 | 11.45566797 | 6.50e-06 | 0.028 |
| 1600 | 11.45567288 | 1.59e-06 | 0.103 |

| u_max | n_q | price | error |
|---|---|---|---|
| 15 | 100 | 11.45489099 | 7.83e-04 |
| 20 | 100 | 11.45573133 | 5.69e-05 |
| 30 | 150 | 11.45566801 | 6.46e-06 |
| 60 | 200 | 11.45566797 | 6.50e-06 |
| 60 | 400 | 11.45566797 | 6.50e-06 |
| 120 | 400 | 11.45566797 | 6.50e-06 |

## Phase 3 --- Fourier against Monte-Carlo

Monte-Carlo: 400,000 antithetic paths, 200 steps, with the nu = 0 control variate.

| H | Fourier | s | Monte-Carlo | 95% band | s | MC - Fourier | in bands | negative-variance hits | Fourier faster by |
|---|---|---|---|---|---|---|---|---|---|
| 0.05 | 11.448809 | 0.318 | 11.4920 | ±0.0219 | 5.05 | +0.0432 | 2.0x | 2,468,610 | **15.9x** |
| 0.1 | 11.455674 | 0.303 | 11.4843 | ±0.0198 | 5.0 | +0.0286 | 1.4x | 1,290,440 | **16.5x** |
| 0.2 | 11.470413 | 0.304 | 11.4865 | ±0.0164 | 4.99 | +0.0161 | 1.0x | 209,106 | **16.4x** |
| 0.3 | 11.487021 | 0.307 | 11.4966 | ±0.0140 | 5.01 | +0.0095 | 0.7x | 15,529 | **16.3x** |
| 0.45 | 11.516951 | 0.308 | 11.5191 | ±0.0117 | 5.05 | +0.0022 | 0.2x | 124 | **16.4x** |

## Phase 4 --- the Monte-Carlo's bias floor

H = 0.1, 200,000 paths throughout; Fourier reference **11.455674**.  The Euler--Volterra
scheme truncates the variance at zero, and that bias is a function of the step size, not
of the number of paths.

| steps | price | 95% band | bias vs Fourier | in bands | negative-variance hits | seconds |
|---|---|---|---|---|---|---|
| 50 | 11.4875 | ±0.0256 | +0.0318 | 1.2x | 180,886 | 0.2 |
| 100 | 11.4729 | ±0.0268 | +0.0172 | 0.6x | 336,700 | 0.67 |
| 200 | 11.4917 | ±0.0281 | +0.0360 | 1.3x | 645,660 | 2.53 |
| 400 | 11.4647 | ±0.0291 | +0.0091 | 0.3x | 1,268,848 | 11.33 |

## Phase 5 --- cost to a target accuracy (H = 0.10)

Monte-Carlo bias floor over the step counts tested: **0.0091**.  A target below that floor cannot be
reached by adding paths at all.

| target | Fourier: steps | seconds | MC: paths needed | seconds | reachable? | Fourier faster by |
|---|---|---|---|---|---|---|
| 0.01 | 100 | 0.004 | 1.57e+06 | 19.7 | yes | 4915.1x |
| 0.001 | 100 | 0.004 | 1.57e+08 | 1966.0 | **no** | 491512.5x |
| 0.0001 | 200 | 0.005 | 1.57e+10 | 196605.0 | **no** | 39320999.3x |

## Compute time

| phase | wall time | seconds |
|---|---|---|
| 1. exact controls | 4s | 3.67 |
| 2. the reference's own accuracy | 1s | 0.72 |
| 3. Fourier against Monte-Carlo, with wall time | 27s | 26.65 |
| 4. the Monte-Carlo's bias floor | 15s | 15.04 |
| **total** | **46s** | **46.09** |

Machine: Darwin arm64, python 3.9.6.
