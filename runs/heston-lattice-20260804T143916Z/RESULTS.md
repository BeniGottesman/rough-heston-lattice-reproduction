# The lattice against a reference with no Monte-Carlo in it

Rough Heston is OUTSIDE the paper's model class: that class needs
v_t = F(t, (K^h Y)(t)) with Y an autonomous diffusion, and rough Heston's nu sqrt(V) dB
makes the driver's coefficients depend on the variance, hence on the driver's own past.
Conversely the models with a semi-analytic transform are affine, and in the rough setting
affineness forces that same sqrt(V) dB.  So NO ROUGH MODEL both fits the class and has a
closed-form transform, and the lattice can never be checked against Fourier in the rough
regime.  The two classes meet only at h = 0: there K^0 is the identity, v = v_0 + y is
autonomous in y, and the model is classical Heston.  That is the test performed here.

Parameters {'V0': 0.09, 'theta': 0.09, 'lam': 1.5, 'nu': 0.3, 'rho': -0.7}, T = 1, S0 = K = 100.  Feller ratio 2*lambda*theta/nu^2 = 3.00.

## Phase 1 --- the reference at H = 1/2

Fourier price **11.528559** at 1600 Riccati steps, in 0.076s.

| control | result |
|---|---|
| nu = lambda = 0, V0 = 0.04 gives Black-Scholes | error 2.46e-11 |
| nu = lambda = 0, V0 = 0.09 gives Black-Scholes | error 2.45e-11 |
| phi(0) = 1 | error 0.00e+00 |
| phi(-i) = 1 (martingale) | error 0.00e+00 |
| self-convergence at 200 Riccati steps | 6.60e-06 |
| self-convergence at 400 Riccati steps | 1.56e-06 |
| self-convergence at 800 Riccati steps | 3.11e-07 |

## Phase 2 --- the lattice against it

| n | mref | lattice | signed error | seconds | grid | offsets | driver p outside [0,1] |
|---|---|---|---|---|---|---|---|
| 8 | 4 | 11.47869 | -0.04987 | 0.01 | 625 | 79 | 6 |
| 16 | 6 | 11.50385 | -0.02471 | 0.07 | 2305 | 145 | 30 |
| 32 | 8 | 11.51849 | -0.01007 | 1.04 | 7745 | 243 | 144 |
| 64 | 12 | 11.52274 | -0.00582 | 27.02 | 30209 | 473 | 650 |

Order in delta: **1.059**, verdict **converging**.

The error shrinks monotonically.  This is the first time this lattice has been measured
against a price with no Monte-Carlo in it.

## Phase 3 --- the small-nu behaviour, which turns out to be benign

The Lamperti transform U = 2 sqrt(v)/nu is what makes the driver recombine, and its drift
carries a constant 2*lambda*theta/nu^2 that blows up as nu decreases, so this was expected
to be the weak point.  Compared against the Fourier price AT EACH nu -- not against the
nu = 0 limit, which is a different price and was our own initial mistake -- it is not:

| nu | u0 | drift constant | Fourier | lattice error n=8 | 16 | 32 |
|---|---|---|---|---|---|---|
| 0.02 | 30.0 | 674.5 | 11.90934 | -0.06897 | -0.03500 | -0.01580 |
| 0.05 | 12.0 | 107.5 | 11.88472 | -0.06718 | -0.03396 | -0.01542 |
| 0.3 | 2.0 | 2.5 | 11.52856 | -0.04987 | -0.02471 | -0.01007 |

The drift constant grows by a factor of 270 between nu = 0.30 and nu = 0.02 while the error
grows by about 40% and the convergence ORDER is unchanged, so the parametrisation is far more
robust than feared.  The apparent pathology we first reported came from comparing the small-nu
lattice with the nu = 0 Black--Scholes price (11.92354) instead of with the true price at that
nu (11.90934): most of the gap was the genuine nu-effect, not lattice error.

## Phase 4 --- the contrast that makes the test worth running

Same Route A' coupling, same backward induction, rough Bergomi at H = 0.1, eta = 0.3,
judged against the exact-covariance Monte-Carlo (no Fourier price exists there):

| n | error at h = 0, vs Fourier | error at h = -0.4, vs Monte-Carlo |
|---|---|---|
| 8 | -0.04987 | -0.0149 |
| 16 | -0.02471 | -0.0347 |
| 32 | -0.01007 | -0.0031 |
| 64 | -0.00582 | +0.0748 |

The left column shrinks, the right one grows.  The lattice machinery --- moment matching,
recombination, the two-dimensional induction, the Route A' coupling --- is therefore sound;
it is the ROUGHNESS that breaks it, which is what Section 8 predicts and what Route B is
meant to repair.

## Compute time

| phase | wall time | seconds |
|---|---|---|
| 1. the Fourier reference at H = 1/2, and its controls there | 0s | 0.14 |
| 2. the lattice against the Fourier reference | 28s | 28.15 |
| 3. the small-nu behaviour of the Lamperti transform | 2s | 1.79 |
| **total** | **30s** | **30.08** |

Machine: Darwin arm64, python 3.9.6.
