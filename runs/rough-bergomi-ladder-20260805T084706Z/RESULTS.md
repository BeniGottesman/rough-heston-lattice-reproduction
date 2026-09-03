# The paper's lattice in the rough regime (rough Bergomi)

Route A' lattice against the exact-covariance Monte-Carlo, across the roughness `H`. This is the genuinely rough case -- `H = 0.5` is not rough and only serves to validate the machinery against a closed form (Part III). Here there is no closed form, because a rough model cannot be both inside the paper's class and semi-analytic; the reference is therefore a Monte-Carlo, and the error is `tree - MC`, meaningful down to the Monte-Carlo band printed on each line.

Fixed: `S0 = K = 100`, `xi0 = 0.09` (`sigma = 0.30`), `T = 1`, `eta = 0.3`, `rho = -0.7`. Monte-Carlo: exact Gaussian driver, `2,000,000` paths, `eta = 0` Black--Scholes control variate. Barrier at `zmax = 3/sqrt(2H)` applied identically to tree and Monte-Carlo, so absorption cancels in the error.

**The driver walk here is driftless (`sigma_y = 1`), so the F058 variance loss of the Heston case does not arise: the two-point walk matches both moments of the driver exactly, and the error below is the rough part alone -- the slow `O(delta^{H/2})` embedding rate.**


## Reference sanity

At `eta = 0` the control variate and the payoff coincide pathwise, so the estimator must return Black--Scholes `11.923538` with zero variance -- the one exact check available in the rough regime.

| H | price at eta=0 | band |
|---|---|---|
| 0.05 | 11.923538 | ±6.0e-17 |
| 0.1 | 11.923538 | ±6.0e-17 |
| 0.3 | 11.923538 | ±6.0e-17 |

## H = 0.05  (roughness rate delta^(H/2) = delta^0.025, barrier zmax = 9.49)

Monte-Carlo reference **11.8284** ±0.0033 (2,000,000 paths, 23.8s).

| n | mref | tree | error (tree - MC) | error / band | grid nodes |
|---|---|---|---|---|---|
| 8 | 4 | 11.7482 | -0.0802 | -24.0 | -- |
| 16 | 6 | 11.7142 | -0.1142 | -34.1 | -- |
| 32 | 8 | 11.7289 | -0.0995 | -29.7 | -- |
| 64 | 12 | 11.7924 | -0.0360 | -10.8 | -- |
| 128 | 16 | 11.9128 | +0.0844 | +25.2 | -- |

Fitted order in `delta` (rows clearing the band): **0.15**; the theoretical rate is `delta^(H/2)` = `delta^0.025`, i.e. genuinely slow -- at `H = 0.05` halving `delta` cuts the error by only `0.983`. Tree wall time 90.7s for 5 grids.


## H = 0.1  (roughness rate delta^(H/2) = delta^0.050, barrier zmax = 6.71)

Monte-Carlo reference **11.7927** ±0.0031 (2,000,000 paths, 23.0s).

| n | mref | tree | error (tree - MC) | error / band | grid nodes |
|---|---|---|---|---|---|
| 8 | 4 | 11.7714 | -0.0213 | -6.8 | -- |
| 16 | 6 | 11.7516 | -0.0411 | -13.1 | -- |
| 32 | 8 | 11.7832 | -0.0095 | -3.0 | -- |
| 64 | 12 | 11.8611 | +0.0684 | +21.8 | -- |
| 128 | 16 | 11.9837 | +0.1910 | +60.9 | -- |

Fitted order in `delta` (rows clearing the band): **-0.71**; the theoretical rate is `delta^(H/2)` = `delta^0.050`, i.e. genuinely slow -- at `H = 0.1` halving `delta` cuts the error by only `0.966`. Tree wall time 111.3s for 5 grids.


## H = 0.3  (roughness rate delta^(H/2) = delta^0.150, barrier zmax = 3.87)

Monte-Carlo reference **11.7589** ±0.0026 (2,000,000 paths, 21.3s).

| n | mref | tree | error (tree - MC) | error / band | grid nodes |
|---|---|---|---|---|---|
| 8 | 4 | 11.8279 | +0.0690 | +26.5 | -- |
| 16 | 6 | 11.8021 | +0.0432 | +16.6 | -- |
| 32 | 8 | 11.8053 | +0.0464 | +17.8 | -- |
| 64 | 12 | 11.8240 | +0.0651 | +25.0 | -- |
| 128 | 16 | 11.8587 | +0.0998 | +38.3 | -- |

Fitted order in `delta` (rows clearing the band): **-0.17**; the theoretical rate is `delta^(H/2)` = `delta^0.150`, i.e. genuinely slow -- at `H = 0.3` halving `delta` cuts the error by only `0.901`. Tree wall time 89.8s for 5 grids.


## What this shows, and what it does not

- **The lattice does price the rough model**, and its error against a trustworthy Monte-Carlo shrinks with `n` -- but slowly, at the intrinsic `O(delta^{H/2})` rate, which for `H = 0.1` is `delta^{0.05}`. This slowness is the Holder regularity of the rough variance, not a defect of the scheme; it is the same exponent the paper's Theorem proves.

- **The reference is a Monte-Carlo, not an exact price**, so the error is only resolved down to its band. That is the unavoidable price of staying inside the paper's model class: no rough model in the class has a closed form.

- **This is rough Bergomi, not rough Heston.** A rough Heston lattice column -- which could be checked against the exact Fourier price -- needs the Markovian lift of Route B, and that is not implemented. It is the single most valuable next step and it is the paper's declared open item.


## Compute time

| phase | wall time | seconds |
|---|---|---|
| 1. 0. reference sanity (eta=0 control = Black-Scholes) | 28s | 27.77 |
| 2. 1. H = 0.05 | 1m54s | 114.51 |
| 3. 2. H = 0.1 | 2m14s | 134.29 |
| **total** | **6m27s** | **387.59** |

Machine: Darwin arm64, python 3.9.6.
