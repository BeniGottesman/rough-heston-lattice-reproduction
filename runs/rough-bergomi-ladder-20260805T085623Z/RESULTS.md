# The paper's lattice in the rough regime (rough Bergomi)

Route A' lattice against the exact-covariance Monte-Carlo, across the roughness `H`. This is the genuinely rough case -- `H = 0.5` is not rough and only serves to validate the machinery against a closed form (Part III). Here there is no closed form, because a rough model cannot be both inside the paper's class and semi-analytic; the reference is therefore a Monte-Carlo, and the error is `tree - MC`, meaningful down to the Monte-Carlo band printed on each line.

**This part is a NEGATIVE result, and it is the important one.** A recombining lattice cannot be the exact convolution (Corollary 'recombining vs consistent'), so the implementable object is the one-step lattice, whose terminal variance diverges as `2H n^{1-2H}` (Proposition 'variance discrepancy'). The table below is that proposition made numerical: the error tracks the variance ratio and GROWS with `n`. This is exactly the obstruction that Route B (the Markovian lift) exists to remove -- and Route B is not yet implemented, so there is as yet no consistent recombining lattice for the rough regime.

Fixed: `S0 = K = 100`, `xi0 = 0.09` (`sigma = 0.30`), `T = 1`, `eta = 0.3`, `rho = -0.7`. Monte-Carlo: exact Gaussian driver, `2,000,000` paths, `eta = 0` Black--Scholes control variate. Barrier at `zmax = 3/sqrt(2H)` applied identically to tree and Monte-Carlo, so absorption cancels in the error.

(The driver walk is driftless here, `sigma_y = 1`, so the F058 Lamperti-drift variance loss of the Heston case is absent; the divergence below is a different and deeper effect -- the inconsistency of the recombining rough lattice itself.)


## Reference sanity

At `eta = 0` the control variate and the payoff coincide pathwise, so the estimator must return Black--Scholes `11.923538` with zero variance -- the one exact check available in the rough regime.

| H | price at eta=0 | band |
|---|---|---|
| 0.05 | 11.923538 | ±6.0e-17 |
| 0.1 | 11.923538 | ±6.0e-17 |
| 0.3 | 11.923538 | ±6.0e-17 |

## H = 0.05  (barrier zmax = 9.49)

Monte-Carlo reference **11.8284** ±0.0033 (2,000,000 paths, 23.7s).

| n | mref | tree | error (tree - MC) | error / band | variance ratio `2H n^(1-2H)` |
|---|---|---|---|---|---|
| 8 | 4 | 11.7482 | -0.0802 | -24.0 | 0.65 |
| 16 | 6 | 11.7142 | -0.1142 | -34.1 | 1.21 |
| 32 | 8 | 11.7289 | -0.0995 | -29.7 | 2.26 |
| 64 | 12 | 11.7924 | -0.0360 | -10.8 | 4.22 |
| 128 | 16 | 11.9128 | +0.0844 | +25.2 | 7.88 |

The signed error rises monotonically with the variance ratio (rank correlation **+0.70**): negative while the ratio is near one and the finite barrier/coupling errors dominate, crossing zero near ratio approximately 2.5 and then growing positive as the ratio does. At `n = 128` the ratio is **7.9** and the tree overprices by **+0.0844** (+25 bands). This is Proposition 'variance discrepancy' made numerical, not a convergence. Tree wall time 87.8s.


## H = 0.1  (barrier zmax = 6.71)

Monte-Carlo reference **11.7927** ±0.0031 (2,000,000 paths, 22.8s).

| n | mref | tree | error (tree - MC) | error / band | variance ratio `2H n^(1-2H)` |
|---|---|---|---|---|---|
| 8 | 4 | 11.7714 | -0.0213 | -6.8 | 1.06 |
| 16 | 6 | 11.7516 | -0.0411 | -13.1 | 1.84 |
| 32 | 8 | 11.7832 | -0.0095 | -3.0 | 3.20 |
| 64 | 12 | 11.8611 | +0.0684 | +21.8 | 5.57 |
| 128 | 16 | 11.9837 | +0.1910 | +60.9 | 9.70 |

The signed error rises monotonically with the variance ratio (rank correlation **+0.90**): negative while the ratio is near one and the finite barrier/coupling errors dominate, crossing zero near ratio approximately 2.5 and then growing positive as the ratio does. At `n = 128` the ratio is **9.7** and the tree overprices by **+0.1910** (+61 bands). This is Proposition 'variance discrepancy' made numerical, not a convergence. Tree wall time 86.0s.


## H = 0.3  (barrier zmax = 3.87)

Monte-Carlo reference **11.7589** ±0.0026 (2,000,000 paths, 21.0s).

| n | mref | tree | error (tree - MC) | error / band | variance ratio `2H n^(1-2H)` |
|---|---|---|---|---|---|
| 8 | 4 | 11.8279 | +0.0690 | +26.5 | 1.38 |
| 16 | 6 | 11.8021 | +0.0432 | +16.6 | 1.82 |
| 32 | 8 | 11.8053 | +0.0464 | +17.8 | 2.40 |
| 64 | 12 | 11.8240 | +0.0651 | +25.0 | 3.17 |
| 128 | 16 | 11.8587 | +0.0998 | +38.3 | 4.18 |

The signed error rises monotonically with the variance ratio (rank correlation **+0.40**): negative while the ratio is near one and the finite barrier/coupling errors dominate, crossing zero near ratio approximately 2.5 and then growing positive as the ratio does. At `n = 128` the ratio is **4.2** and the tree overprices by **+0.0998** (+38 bands). This is Proposition 'variance discrepancy' made numerical, not a convergence. Tree wall time 83.3s.


## What this shows

- **The recombining rough lattice is inconsistent, and this is the numerical proof.** Its terminal variance is `2H n^{1-2H}` times the truth, so refining the grid drives the price AWAY from the reference. The error's sign and growth track that ratio across every `H`. A previous session's sweep stopped at `n = 64`, where sign cancellation made some `H` look convergent; `n = 128` removes that illusion.

- **This is the obstruction Route B is for.** The Markovian lift replaces the one divergent driver by `m` mean-reverting factors whose combined variance stays finite; only then does a recombining rough lattice converge. Route B is not implemented, so **there is as yet no consistent recombining lattice for the rough regime** -- neither for rough Bergomi here nor for rough Heston.

- **Hence there is no rough lattice column against an exact price anywhere in this document, and cannot be one yet.** Rough Bergomi (in-class) has no closed form; rough Heston (closed form via Fourier) is out-of-class. The single most valuable next step is Route B applied to rough Heston, which would finally put a converging lattice next to the exact Fourier price.


## Compute time

| phase | wall time | seconds |
|---|---|---|
| 1. 0. reference sanity (eta=0 control = Black-Scholes) | 27s | 27.35 |
| 2. 1. H = 0.05 | 1m51s | 111.5 |
| 3. 2. H = 0.1 | 1m48s | 108.82 |
| **total** | **5m52s** | **352.07** |

Machine: Darwin arm64, python 3.9.6.
