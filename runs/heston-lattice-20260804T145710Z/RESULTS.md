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

Fourier price **11.528559** at 3200 Riccati steps, in 0.284s.

| control | result |
|---|---|
| nu = lambda = 0, V0 = 0.04 gives Black-Scholes | error 2.46e-11 |
| nu = lambda = 0, V0 = 0.09 gives Black-Scholes | error 2.45e-11 |
| phi(0) = 1 | error 0.00e+00 |
| phi(-i) = 1 (martingale) | error 0.00e+00 |
| self-convergence at 200 Riccati steps | 6.67e-06 |
| self-convergence at 400 Riccati steps | 1.64e-06 |
| self-convergence at 800 Riccati steps | 3.89e-07 |

## Phase 2 --- the lattice against it, to n = 128

The 1/U term of the Lamperti drift diverges at the origin, and a one-step +-sqrt(delta) walk
cannot represent a drift larger than 1/sqrt(delta).  Two regularisations are therefore
reported: `clip` lets the band reach zero and clips the offending up-probabilities, `floor`
raises the lower band to the level where no clipping is needed.  They must agree once the
grid is fine enough, and the rate is fitted only where they do.

### lower barrier: clip

| n | mref | lattice | signed error | seconds | grid | v at lower barrier | p outside [0,1] |
|---|---|---|---|---|---|---|---|
| 8 | 4 | 11.478671 | -0.049889 | 0.01 | 513 | 0.00000 | 6 |
| 16 | 6 | 11.503814 | -0.024745 | 0.05 | 1505 | 0.00000 | 30 |
| 32 | 8 | 11.518468 | -0.010091 | 0.39 | 4033 | 0.00000 | 144 |
| 64 | 12 | 11.522719 | -0.005840 | 5.28 | 11905 | 0.00000 | 650 |
| 128 | 16 | 11.525989 | -0.002571 | 82.56 | 31745 | 0.00000 | 2916 |

order in delta, fitted from n = 32: **0.986**

### lower barrier: floor

| n | mref | lattice | signed error | seconds | grid | v at lower barrier | p outside [0,1] |
|---|---|---|---|---|---|---|---|
| 8 | 4 | 11.979680 | +0.451121 | 0.01 | 513 | 0.07031 | 0 |
| 16 | 6 | 11.523925 | -0.004634 | 0.04 | 1505 | 0.03516 | 0 |
| 32 | 8 | 11.518252 | -0.010308 | 0.33 | 4033 | 0.01758 | 0 |
| 64 | 12 | 11.522566 | -0.005993 | 4.6 | 11905 | 0.00879 | 0 |
| 128 | 16 | 11.525959 | -0.002601 | 60.35 | 31745 | 0.00439 | 0 |

order in delta, fitted from n = 32: **0.993**

### the two regularisations agree, and better as n grows

| n | clip | floor | difference |
|---|---|---|---|
| 8 | -0.049889 | +0.451121 | 5.0e-01 |
| 16 | -0.024745 | -0.004634 | 2.0e-02 |
| 32 | -0.010091 | -0.010308 | 2.2e-04 |
| 64 | -0.005840 | -0.005993 | 1.5e-04 |
| 128 | -0.002571 | -0.002601 | 3.0e-05 |

Below n = 32 the two disagree, because the `floor` variant imposes a variance floor that
bites on coarse grids; from n = 32 up they agree to a few units in the fourth decimal and
the agreement improves with n, so the measured rate does not depend on the choice.

## Phase 2b --- what the barrier costs and what it saves

At n = 32, against the unbarriered value 11.518493.

| barrier | lattice | error vs Fourier | barrier cost | grid | offsets | seconds | speed-up |
|---|---|---|---|---|---|---|---|
| 3 sd | 11.517338 | -0.011221 | 1.2e-03 | 3393 | 107 | 0.22 | 5.1x |
| 4 sd | 11.518350 | -0.010209 | 1.4e-04 | 3713 | 117 | 0.28 | 4.0x |
| 5 sd | 11.518468 | -0.010091 | 2.5e-05 | 4033 | 127 | 0.31 | 3.6x |
| 6 sd | 11.518489 | -0.010070 | 4.0e-06 | 4289 | 135 | 0.46 | 2.4x |
| none | 11.518493 | -0.010066 | 0.0e+00 | 7745 | 243 | 1.12 | 1.0x |

At five standard deviations the barrier costs about 2e-5 on the price -- some five hundred
times below the discretisation error it is measuring -- while halving the grid and cutting
the time by a factor of three at n = 32 and by nearly seven at n = 64.  It is what makes
n = 128 reachable at all.

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
| 8 | -0.04989 | -0.0149 |
| 16 | -0.02475 | -0.0347 |
| 32 | -0.01009 | -0.0031 |
| 64 | -0.00584 | +0.0748 |

The left column shrinks, the right one grows.  The lattice machinery --- moment matching,
recombination, the two-dimensional induction, the Route A' coupling --- is therefore sound;
it is the ROUGHNESS that breaks it, which is what Section 8 predicts and what Route B is
meant to repair.

## Compute time

| phase | wall time | seconds |
|---|---|---|
| 1. the Fourier reference at H = 1/2, and its controls there | 0s | 0.35 |
| 2. the lattice against the Fourier reference, to n = 128 | 2m33s | 153.62 |
| 3. the barrier: cost against saving | 2s | 2.39 |
| 4. the small-nu behaviour of the Lamperti transform | 2s | 2.19 |
| **total** | **2m38s** | **158.55** |

Machine: Darwin arm64, python 3.9.6.
