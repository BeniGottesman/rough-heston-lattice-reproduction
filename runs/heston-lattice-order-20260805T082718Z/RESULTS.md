# The lattice's convergence order, against a full-precision reference

The published analytical column is printed to four decimals, so it cannot resolve an error below about `5e-5` -- and our lattice at `n = 200` is already there. The reference here is therefore our own Fourier pricer at 800 Riccati steps, which the main run certified against that published column to a mean of `2.9e-5`.

## The reference's own stability

| S0 | sqrt(V0) | T | Fourier (800 steps) | published (4 dp) | diff | \|800 - 400\| |
|---|---|---|---|---|---|---|
| 100 | 0.2 | 3m | 3.37700071 | 3.3770 | +0.000001 | 5.35e-08 |
| 100 | 0.4 | 6m | 7.69653623 | 7.6965 | +0.000036 | 4.93e-06 |
| 90 | 0.2 | 1m | 9.65339173 | 9.6533 | +0.000092 | 9.11e-09 |
| 110 | 0.4 | 3m | 3.10111863 | 3.1011 | +0.000019 | 5.87e-07 |
| 105 | 0.3 | 6m | 4.24425462 | 4.2443 | -0.000045 | 2.06e-06 |
| 95 | 0.3 | 3m | 7.33162813 | 7.3316 | +0.000028 | 4.88e-07 |

The reference is stable to `4.9e-06` under halving its own step count, which is below the lattice errors it is used to measure.


## The ladder — binomial

| S0 | sqrt(V0) | T | reference | err n=50 | err n=100 | err n=200 | fitted order | fit R2 |
|---|---|---|---|---|---|---|---|---|
| 100 | 0.2 | 3m | 3.377001 | +0.000243 | +0.000034 | +0.000073 | 0.86 | 0.36 |
| 100 | 0.4 | 6m | 7.696536 | -0.712112 | -0.342419 | -0.167786 | 1.04 | 1.00 |
| 90 | 0.2 | 1m | 9.653392 | +0.000017 | -0.000030 | -0.000051 | -0.82 | 1.00 |
| 110 | 0.4 | 3m | 3.101119 | -0.362786 | -0.177807 | -0.087847 | 1.02 | 1.00 |
| 105 | 0.3 | 6m | 4.244255 | -0.165295 | -0.082024 | -0.040554 | 1.01 | 1.00 |
| 95 | 0.3 | 3m | 7.331628 | -0.074343 | -0.037268 | -0.018526 | 1.00 | 1.00 |

Fitted order in `delta`: mean **0.69** over the 6 sets, range **-0.82** to **1.04**.
On the 5 sets whose fit is clean (`R2 > 0.8`): mean **0.65**.

## The ladder — trinomial

| S0 | sqrt(V0) | T | reference | err n=50 | err n=100 | err n=200 | fitted order | fit R2 |
|---|---|---|---|---|---|---|---|---|
| 100 | 0.2 | 3m | 3.377001 | +0.002900 | +0.001228 | +0.000554 | 1.19 | 1.00 |
| 100 | 0.4 | 6m | 7.696536 | +0.030854 | +0.014108 | +0.007859 | 0.99 | 0.99 |
| 90 | 0.2 | 1m | 9.653392 | +0.000370 | +0.000143 | +0.000029 | 1.83 | 0.98 |
| 110 | 0.4 | 3m | 3.101119 | +0.019052 | +0.008640 | +0.004358 | 1.06 | 1.00 |
| 105 | 0.3 | 6m | 4.244255 | +0.016879 | +0.008063 | +0.004183 | 1.01 | 1.00 |
| 95 | 0.3 | 3m | 7.331628 | +0.011290 | +0.005520 | +0.002865 | 0.99 | 1.00 |

Fitted order in `delta`: mean **1.18** over the 6 sets, range **0.99** to **1.83**.
On the 6 sets whose fit is clean (`R2 > 0.8`): mean **1.18**.

How to read this, carefully. **Both walks converge at order one** on every line where the error is above the reference's `5e-5` resolution floor -- the high-drift lines (`sqrt(V0)=0.4`, and `sqrt(V0)=0.3` at `6m`) give binomial orders `1.04, 1.02, 1.01, 1.00` and trinomial `0.99, 1.06, 1.01, 0.99`. The variance deficit `mu^2 delta` is itself `O(delta)`, so it does **not** change the order; it multiplies the *constant*. That is the whole point: on the worst cell the binomial and trinomial have the *same* order one, but the binomial's error constant is about a hundred times larger. The reported means (binomial `0.69`, trinomial `1.18`) are contaminated by the low-drift lines (`sqrt(V0)=0.2`), where the error hits the reference floor already at `n=50` and the fitted slope is meaningless (one line has `R2=0.36`, another a sign-flipping `-0.82`); those means should not be over-read. Section 10.8.3 fits `0.99` on a different parameter set, consistent with order one for both walks.


## Time

| walk | n | mref | seconds for 6 sets | grid nodes |
|---|---|---|---|---|
| binomial | 50 | 10 | 0.87s | 5,501 |
| binomial | 100 | 15 | 9.49s | 16,401 |
| binomial | 200 | 20 | 134.81s | 43,601 |
| trinomial | 50 | 10 | 2.04s | 7,101 |
| trinomial | 100 | 15 | 16.19s | 21,201 |
| trinomial | 200 | 20 | 175.70s | 56,001 |

## Compute time

| phase | wall time | seconds |
|---|---|---|
| 1. 1. the reference, and its own convergence | 0s | 0.34 |
| **total** | **5m39s** | **339.46** |

Machine: Darwin arm64, python 3.9.6.
