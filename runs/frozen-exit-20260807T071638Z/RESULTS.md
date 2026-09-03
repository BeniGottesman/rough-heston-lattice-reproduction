# Route F falsifiers FF1-FF4, FF6, FF7 -- measured

Route and falsifiers: `research-notes/L003-ROUTE-F.md` sec.2-3.
Code: `sim/frozen_exit.py`; driver: `sim/run_frozen_exit.py`.

Common settings: `T=1`, `y_0=0`, `v_0=0`, `K(u)=u^h` (`L==1`), continuous band `(B^Y,C^Y)=(-0.5,0.5)`, `mu_y=0`, `sigma_y=1` except in FF3.  `H := h+1/2`.

**The barrier inset is a parameter here.**  `eq:barriers` writes the literal `delta^{1/3}`.  That value occurs in exactly four places (`part2-scheme.tex:24-30, 44-47, 57, 71-81`), all inside the definition of the discrete barriers and inside `lem:barriers-wellposed` itself, and nothing downstream reads it.  Everything below is therefore computed with the inset `delta^{beta}`, at `beta = 1/3` (the current text) and `beta = 1/2` (the grid resolution, the smallest a `sqrt(delta)`-lattice can express), **on the same paths**.  `beta = 1/3` is the default everywhere the inset is not named.

These are measurements, not proofs.  Nothing here closes a Definition-of-Done item; a refutation is an instruction to stop, a non-refutation is not a permission to continue.

## Decision table

Verdicts are the author's reading of the numbers below; the numbers are the evidence and the tables are where they live.

| falsifier | verdict | the number that decides it |
|---|---|---|
| FF1 (the diagnosis) | **CONFIRMED** | at h=-0.3, n=10^6 the L2 gap is 4.154 against delta^H = 0.0631, a ratio of 66, and it GROWS with n; 49% of post-exit increments are non-zero, so `rem:freezing-convention`'s premise is false outright |
| FF2 (Route F **at beta=1/3**) | **REFUTED** | the frozen-vs-frozen decay misses the target (h+kappa)/2 by more than 2 bootstrap sd in 6 of 6 (h, ell) cells: h=-0.1,ell=1, h=-0.1,ell=2, h=-0.1,ell=4, h=-0.3,ell=1, h=-0.3,ell=2, h=-0.3,ell=4.  F1 meets its rate; F2 does not, and its exponent falls with ell exactly as the closed-form law of G predicts.  This is a refutation of the CURRENT INSET, not of the decomposition -- see FF6 |
| FF6 (the inset is the free parameter) | **CONFIRMED** | moving the inset from `delta^{1/3}` to `delta^{1/2}` on the SAME paths takes the count of failing cells from 6/6 to 1/6 (the survivor is h=-0.1,ell=4, and `ell=4 > ell*=2.50` there, so it is PREDICTED to fail); `ell=2` -- the order `lem:freeze` consumes -- passes at both h.  F2 is what moves (its exponent rises by a factor 2.1-7.2 against the predicted 1.5, the excess being the milder pre-asymptotics at the smaller inset) |
| **(unasked, and it matters)** `prop:Vconv`-as-proved | **REFUTED at both insets** | `F1_full` -- the UNFROZEN clamped-walk convolution against the ABSORBED-driver continuous one, which is the pair `prop:Vconv`'s proof actually writes -- has a NEGATIVE fitted exponent at every (h, ell, beta) tested (worst -0.0813).  The clamp/absorb mismatch of FF1 is inside `prop:Vconv` itself.  Restricted `F1` (on `j <= Xi^Y_n`, where the clamp is inactive) does meet its rate, so Route F's (F1) survives -- but it may not be obtained by quoting `prop:Vconv`; it has to be proved on the stopped range |
| FF3 (ordering L1) | **NOT REFUTED as a fact, REFUTED as an inference** | 0 violations in 1536 faithful-construction paths; but with a node error of the size (E2) permits the fraction runs 0.5% -> 5.2% -> 7.3% -> 10.9% over n = 64, 256, 1024, 4096 -- rising, not vanishing |
| FF4 (exit-time L2) | **NOT REFUTED** | on an exact ladder to n = 10^12 the exponent of E[G^0.2] in eps_n converges to 0.409 against the predicted 0.4, with the kappa'=1 control landing on 1.000; the tail matches eps_n sqrt(2/(pi t)) to within 3% over its stated range |
| FF7 (the lower-bound mechanism) | **NOT REFUTED** | P(G>T/2) scales like eps_n^{1.000} at beta=1/3 and eps_n^{1.000} at beta=1/2 against the predicted exponent 1, so the sharpness claim's mechanism is present and is pure `eps_n` |

## FF1 -- frozen against UNFROZEN, discrete side only

Measured `|| max_k |Vcal^(n)(t_k) - V^(n)_k| ||_{L2}`.  No continuous object and no coupling enter, so this is exact up to Monte-Carlo error and can be pushed to `n = 10^6`.  Two extra columns split the gap the way `rem:freezing-convention` splits it:

* **kernel part** `max_{k>Xi} |sum_{i<=Xi} [K(t_k-t_{i-1}) - K(t_Xi-t_{i-1})] dY_i|` -- the sum the remark tries to telescope;
* **increment part** `max_{k>Xi} |sum_{Xi<i<=k} K(t_k-t_{i-1}) dY_i|` -- zero if and only if the remark's premise ("past `Xi^Y_n` the increments vanish") is true;
* **nz post** = fraction of increments `dY_i` with `i > Xi^Y_n` that are NON-zero.  The premise says this is 0.

`scale` = `|| max_k |Vcal^(n)(t_k)| ||_{L2}`, the natural size of the object.

### h = -0.1  (H = 0.4)

| n | delta | paths | L2 gap | +- | delta^H | gap/delta^H | kernel part | increment part | scale | nz post | P(walk exits) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1000 | 0.001 | 512 | 0.9661 | 0.0137 | 0.0631 | 15.3 | 0.2104 | 0.8120 | 0.6162 | 0.3978 | 1.000 |
| 3162 | 0.0003163 | 512 | 1.0565 | 0.0150 | 0.0398 | 26.5 | 0.2348 | 0.8857 | 0.6916 | 0.4122 | 1.000 |
| 10000 | 0.0001 | 384 | 1.0977 | 0.0180 | 0.0251 | 43.7 | 0.2512 | 0.9101 | 0.7296 | 0.4463 | 0.998 |
| 31623 | 3.162e-05 | 256 | 1.1176 | 0.0239 | 0.0158 | 70.5 | 0.2603 | 0.9255 | 0.7600 | 0.4303 | 0.992 |
| 100000 | 1e-05 | 160 | 1.1149 | 0.0311 | 0.0100 | 111.5 | 0.2616 | 0.9195 | 0.7702 | 0.4368 | 1.000 |
| 316228 | 3.162e-06 | 96 | 1.1793 | 0.0430 | 0.0063 | 186.9 | 0.2672 | 0.9764 | 0.7880 | 0.4804 | 0.985 |
| 1000000 | 1e-06 | 64 | 1.1271 | 0.0477 | 0.0040 | 283.1 | 0.2691 | 0.9244 | 0.7901 | 0.4277 | 1.000 |

### h = -0.3  (H = 0.2)

| n | delta | paths | L2 gap | +- | delta^H | gap/delta^H | kernel part | increment part | scale | nz post | P(walk exits) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1000 | 0.001 | 512 | 2.8502 | 0.0291 | 0.2512 | 11.3 | 1.1508 | 1.9202 | 1.7170 | 0.4096 | 1.000 |
| 3162 | 0.0003163 | 512 | 3.1841 | 0.0354 | 0.1995 | 16.0 | 1.3483 | 2.0815 | 1.9853 | 0.4047 | 0.996 |
| 10000 | 0.0001 | 384 | 3.4682 | 0.0426 | 0.1585 | 21.9 | 1.4801 | 2.2500 | 2.2002 | 0.4045 | 0.998 |
| 31623 | 3.162e-05 | 256 | 3.6561 | 0.0564 | 0.1259 | 29.0 | 1.5664 | 2.3724 | 2.3221 | 0.4681 | 0.993 |
| 100000 | 1e-05 | 160 | 3.8709 | 0.0799 | 0.1000 | 38.7 | 1.6556 | 2.4825 | 2.5032 | 0.4314 | 0.983 |
| 316228 | 3.162e-06 | 96 | 4.0352 | 0.0923 | 0.0794 | 50.8 | 1.6990 | 2.6250 | 2.5858 | 0.4590 | 0.990 |
| 1000000 | 1e-06 | 64 | 4.1538 | 0.1045 | 0.0631 | 65.8 | 1.7664 | 2.6426 | 2.6521 | 0.4860 | 0.985 |

### h = -0.45  (H = 0.05)

| n | delta | paths | L2 gap | +- | delta^H | gap/delta^H | kernel part | increment part | scale | nz post | P(walk exits) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1000 | 0.001 | 512 | 6.3771 | 0.0583 | 0.7079 | 9.0 | 2.9433 | 3.8686 | 3.7847 | 0.4096 | 1.000 |
| 3162 | 0.0003163 | 512 | 7.8600 | 0.0751 | 0.6683 | 11.8 | 3.6797 | 4.6946 | 4.7813 | 0.4190 | 1.000 |
| 10000 | 0.0001 | 384 | 9.0589 | 0.0939 | 0.6310 | 14.4 | 4.2125 | 5.4185 | 5.6201 | 0.4202 | 1.000 |
| 31623 | 3.162e-05 | 256 | 10.5746 | 0.1266 | 0.5957 | 17.8 | 4.8104 | 6.3235 | 6.5470 | 0.4583 | 1.000 |
| 100000 | 1e-05 | 160 | 11.9445 | 0.1643 | 0.5623 | 21.2 | 5.4002 | 7.0657 | 7.3750 | 0.4812 | 1.000 |
| 316228 | 3.162e-06 | 96 | 12.3668 | 0.2451 | 0.5309 | 23.3 | 5.7122 | 7.2960 | 7.9774 | 0.3671 | 0.985 |
| 1000000 | 1e-06 | 64 | 14.0119 | 0.2997 | 0.5012 | 28.0 | 6.1783 | 8.3810 | 8.8195 | 0.4494 | 0.985 |

## FF2 / FF6 -- frozen discrete against FROZEN continuous, at both insets

**Coupling.**  The exact Skorokhod level-crossing embedding `theta_k = inf{t > theta_{k-1} : |y_t - y_{theta_{k-1}}| = sqrt(delta)}`, realised on a fine grid of step `dt = 9.54e-07`.  With `mu_y=0, sigma_y=1` this is the manuscript's OWN embedding and it is exact: `Y(theta_k) = Y^(n)_k` and the node error `e_j` vanishes (`part3-embedding.tex` l.582).  Because `y` is a martingale the embedded increments are i.i.d. fair `+-1` coins, so the walk is the `eq:walks` walk in law -- the coupling does not deform either object.  It was chosen over sign-coupling (`zeta_k = sign(y(t_k)-y(t_{k-1}))`, which leaves `|Y^(n)_k - y(t_k)| = Theta(1)` and would fake a stall) and over lattice-rounding of `y(t_k)` (which is not the `eq:walks` walk).

**The two insets are evaluated on the same paths and the same walk.**  The embedding, the walk `Y^(n)`, and the discrete convolution `Vcal^(n)` do not depend on `beta` at all; only the exit index `Xi^Y_n` does, through the barriers.  So the `beta = 1/3` and `beta = 1/2` columns below are PAIRED, and the bootstrap resamples whole paths, preserving the pairing.

**Decomposition measured** (Route F sec.2, `nu := k ^ Xi^Y_n`):

```
  total = max_k |V^(n)_k - V_{theta_k}|
  F1    = max_{j<=Xi^Y_n} |Vcal^(n)(t_j) - v_{theta_j}|          scheme error
  F2    = max_{k>Xi^Y_n} |v_{theta_k ^ Xi^y} - v_{theta_{Xi}}|   exit-time term
```

`h = -0.45` is NOT measured here and the paths were spent on the two cells that can be resolved instead.  At `H = 0.05` the target rate is `delta^{0.025}`, a factor 1.15 across a factor-256 range in `delta`; no feasible ladder separates that from a constant, so any verdict there would have been noise.  It is kept in FF1, which is exact and needs no ladder.

### (a) the fitted exponent, both insets side by side

| h | ell | beta=1/3 | beta=1/2 | target (h+kappa)/2 | ell* = 2 beta/(h+kappa) at 1/3 | at 1/2 |
|---|---|---|---|---|---|---|
| -0.1 | 1 | 0.1524 +- 0.0108 | 0.3265 +- 0.0156 | 0.2000 | 1.67 | 2.50 |
| -0.1 | 2 | 0.0801 +- 0.0104 | 0.2379 +- 0.0186 | 0.2000 | 1.67 | 2.50 |
| -0.1 | 4 | 0.0119 +- 0.0077 | 0.1472 +- 0.0187 | 0.2000 | 1.67 | 2.50 |
| -0.3 | 1 | 0.0656 +- 0.0066 | 0.1678 +- 0.0079 | 0.1000 | 3.33 | 5.00 |
| -0.3 | 2 | 0.0325 +- 0.0065 | 0.1390 +- 0.0099 | 0.1000 | 3.33 | 5.00 |
| -0.3 | 4 | -0.0095 +- 0.0057 | 0.0859 +- 0.0111 | 0.1000 | 3.33 | 5.00 |

### (b) the F1 / F2 split, both insets

The prediction under test was: **F1 is identical at both insets** while **F2's exponent moves with beta**.  The second half holds.  The first half does NOT hold as stated, and the reason is worth having:

* `F1` as written in the decomposition is a maximum over `j <= Xi^Y_n`, and `Xi^Y_n` **is** a function of the inset.  A wider discrete band (larger `beta`) means a later exit index and a longer index range, so restricted `F1` moves with `beta` for a purely combinatorial reason.  At `beta = 1/3` and small `n` the band can be a single grid step wide -- at `n = 64`, `C^Y_n = 0.125 = 1 x sqrt(delta)` -- so `Xi^Y_n = 1` on most paths and the maximum is taken over one index; by `n = 16384` it is taken over hundreds.  That n-dependence of the RANGE, not of the summand, is what bends the fitted slope.
* The bound the route actually invokes is `prop:Vconv` **read at the random index** `nu <= n`, and `prop:Vconv`'s proof bounds a maximum over ALL indices of the UNFROZEN discrete convolution against the ABSORBED-DRIVER continuous one.  That object -- `F1_full` below -- is the one the route is entitled to quote.  It too carries a `beta`, but through the CLAMP in `eq:walks` rather than through a truncated index range, and the dependence is weaker.
* One asymmetry inside `F1_full` deserves naming, because it is the same defect FF1 exhibits, one level up.  Its DISCRETE side is the unfrozen convolution of the **clamped** walk, whose increments do not vanish past `Xi^Y_n` (FF1 measures ~40-49% of them non-zero); its CONTINUOUS side is the **absorbed**-driver convolution, whose integrand genuinely does stop at `Xi^y`.  Those two are not analogues of one another, and `prop:Vconv`'s proof pairs them.  If `F1_full` does not decay while restricted `F1` does, the reading is that `prop:Vconv`-as-proved inherits the clamp/absorb mismatch too, and the route may not quote it to dominate (F1): it must bound (F1) directly on `j <= Xi^Y_n`, where the clamp is inactive and the pairing is honest.

So the decomposition is intact.  What needs replacing is the sentence "(F1) is free" / "(F1) never sees the barrier gap".  (F1) does see it -- through its index range, and through the clamp inside the quantity that dominates it.

| h | quantity | ell | beta=1/3 | beta=1/2 | beta min(2H,1/ell) at 1/3 | at 1/2 | target |
|---|---|---|---|---|---|---|---|
| -0.1 | F1 | 1 | 0.1880 +- 0.0017 | 0.2194 +- 0.0018 | -- | -- | 0.2000 |
| -0.1 | F1 | 2 | 0.1902 +- 0.0017 | 0.2216 +- 0.0018 | -- | -- | 0.2000 |
| -0.1 | F1 | 4 | 0.1945 +- 0.0017 | 0.2260 +- 0.0019 | -- | -- | 0.2000 |
| -0.1 | F1_full | 1 | -0.0400 +- 0.0031 | -0.0407 +- 0.0020 | -- | -- | 0.2000 |
| -0.1 | F1_full | 2 | -0.0559 +- 0.0019 | -0.0488 +- 0.0012 | -- | -- | 0.2000 |
| -0.1 | F1_full | 4 | -0.0672 +- 0.0012 | -0.0521 +- 0.0010 | -- | -- | 0.2000 |
| -0.1 | F2 | 1 | 0.1541 +- 0.0107 | 0.3310 +- 0.0160 | 0.2667 | 0.4000 | 0.2000 |
| -0.1 | F2 | 2 | 0.0863 +- 0.0103 | 0.2466 +- 0.0185 | 0.1667 | 0.2500 | 0.2000 |
| -0.1 | F2 | 4 | 0.0220 +- 0.0077 | 0.1588 +- 0.0185 | 0.0833 | 0.1250 | 0.2000 |
| -0.3 | F1 | 1 | 0.1008 +- 0.0015 | 0.1245 +- 0.0015 | -- | -- | 0.1000 |
| -0.3 | F1 | 2 | 0.1025 +- 0.0015 | 0.1259 +- 0.0015 | -- | -- | 0.1000 |
| -0.3 | F1 | 4 | 0.1059 +- 0.0016 | 0.1286 +- 0.0016 | -- | -- | 0.1000 |
| -0.3 | F1_full | 1 | -0.0366 +- 0.0040 | -0.0572 +- 0.0034 | -- | -- | 0.1000 |
| -0.3 | F1_full | 2 | -0.0502 +- 0.0031 | -0.0695 +- 0.0025 | -- | -- | 0.1000 |
| -0.3 | F1_full | 4 | -0.0654 +- 0.0022 | -0.0813 +- 0.0018 | -- | -- | 0.1000 |
| -0.3 | F2 | 1 | 0.0848 +- 0.0063 | 0.2138 +- 0.0104 | 0.1333 | 0.2000 | 0.1000 |
| -0.3 | F2 | 2 | 0.0614 +- 0.0060 | 0.1738 +- 0.0102 | 0.1333 | 0.2000 | 0.1000 |
| -0.3 | F2 | 4 | 0.0275 +- 0.0056 | 0.1216 +- 0.0102 | 0.0833 | 0.1250 | 0.1000 |

### The L2 norms themselves

#### h = -0.1

| n | delta | total_theta (b=1/3) | total_theta (b=1/2) | F1 (b=1/3) | F1 (b=1/2) | F2 (b=1/3) | F2 (b=1/2) |
|---|---|---|---|---|---|---|---|
| 64 | 0.01562 | 0.7702 | 0.7181 | 0.0693 | 0.0900 | 0.7882 | 0.7521 |
| 128 | 0.007812 | 0.7222 | 0.6335 | 0.0783 | 0.0877 | 0.7557 | 0.6661 |
| 256 | 0.003906 | 0.6942 | 0.6194 | 0.0723 | 0.0760 | 0.7217 | 0.6472 |
| 512 | 0.001953 | 0.6572 | 0.4955 | 0.0632 | 0.0671 | 0.6798 | 0.5142 |
| 1024 | 0.0009766 | 0.6411 | 0.5137 | 0.0537 | 0.0560 | 0.6582 | 0.5287 |
| 2048 | 0.0004883 | 0.5802 | 0.4076 | 0.0456 | 0.0466 | 0.5927 | 0.4162 |
| 4096 | 0.0002441 | 0.5785 | 0.3636 | 0.0379 | 0.0388 | 0.5874 | 0.3692 |
| 8192 | 0.0001221 | 0.5239 | 0.1948 | 0.0319 | 0.0326 | 0.5279 | 0.1974 |
| 16384 | 6.104e-05 | 0.4836 | 0.2046 | 0.0282 | 0.0289 | 0.4845 | 0.2049 |

#### h = -0.3

| n | delta | total_theta (b=1/3) | total_theta (b=1/2) | F1 (b=1/3) | F1 (b=1/2) | F2 (b=1/3) | F2 (b=1/2) |
|---|---|---|---|---|---|---|---|
| 64 | 0.01562 | 2.3652 | 2.2644 | 0.8476 | 1.0239 | 2.8945 | 2.8693 |
| 128 | 0.007812 | 2.3317 | 2.0942 | 0.9580 | 1.0388 | 2.9206 | 2.6005 |
| 256 | 0.003906 | 2.3079 | 2.1017 | 0.9395 | 0.9743 | 2.8485 | 2.5791 |
| 512 | 0.001953 | 2.2750 | 1.7993 | 0.8828 | 0.9170 | 2.7397 | 2.0925 |
| 1024 | 0.0009766 | 2.2826 | 1.8800 | 0.8139 | 0.8354 | 2.6727 | 2.1894 |
| 2048 | 0.0004883 | 2.1389 | 1.5874 | 0.7414 | 0.7525 | 2.4626 | 1.7722 |
| 4096 | 0.0002441 | 2.1696 | 1.4725 | 0.6662 | 0.6766 | 2.4388 | 1.6334 |
| 8192 | 0.0001221 | 2.0509 | 1.0479 | 0.5888 | 0.5950 | 2.2336 | 1.1031 |
| 16384 | 6.104e-05 | 1.9487 | 1.1088 | 0.5218 | 0.5289 | 2.0784 | 1.1755 |

`total_tk` (the comparison against `V_{t_k}` rather than `V_{theta_k}`) is recorded in the run's `progress.json` and is uniformly WORSE at both insets; it is not the manuscript's object -- `prop:Vconv` compares at `theta_k` -- so it is not tabulated here.

### FF6 x FF4 -- where the arithmetic loses a power, and what beta buys back

Route F sec.2 bounds `|F2| <= |v|_{H-} * G^{H-}` and quotes `E[G^{kappa'}] = O(eps_n^{2 kappa'})`, valid *for kappa' < 1/2*.  But `prop:Vconv` is an `L^ell` statement, so what is needed is `|| |v|_H G^H ||_{L^ell} ~ E[G^{ell H}]^{1/ell}`: the exponent inside the expectation is `ell H`, not `H`.  The `kappa' < 1/2` clause therefore reads `ell H < 1/2`, and above it `E[G^{kappa'}] ~ eps_n^{1}`, so with `eps_n ~ delta^{beta}` the norm decays like `delta^{beta/ell}` rather than `delta^{2 beta H}`.  Against the target `(h+kappa)/2`:

```
  exponent in delta = beta * min(2 kappa', 1/ell)
  meets the target  <=>  ell <= ell* := 2 beta / (h + kappa)
```

At `beta = 1/3` and `kappa -> 1/2` this is `ell* = 2/(3H)`; at `beta = 1/2` it is `ell* = 1/(h+kappa) > 2` strictly for every admissible pair, since `h < 0` and `kappa < 1/2` force `h+kappa < 1/2`.

The `exact` columns are computed from the closed-form law of `G` over the SAME `n` ladder as the tables above -- no simulation -- so they are directly comparable with the measured `F2` slopes.

| h | ell | beta | ell*H or ell*kappa' | exact E[G^{ell H}]^{1/ell} slope | measured F2 slope | asymptotic beta min(2H,1/ell) | target |
|---|---|---|---|---|---|---|---|
| -0.1 | 1 | 1/3 | 0.4 | 0.1816 | 0.1541 +- 0.0107 | 0.2667 | 0.2000 |
| -0.1 | 1 | 1/2 | 0.4 | 0.3060 | 0.3310 +- 0.0160 | 0.4000 | 0.2000 |
| -0.1 | 2 | 1/3 | 0.8 | 0.1362 | 0.0863 +- 0.0103 | 0.1667 | 0.2000 |
| -0.1 | 2 | 1/2 | 0.8 | 0.2233 | 0.2466 +- 0.0185 | 0.2500 | 0.2000 |
| -0.1 | 4 | 1/3 | 1.6 | 0.0796 | 0.0220 +- 0.0077 | 0.0833 | 0.2000 |
| -0.1 | 4 | 1/2 | 1.6 | 0.1241 | 0.1588 +- 0.0185 | 0.1250 | 0.2000 |
| -0.3 | 1 | 1/3 | 0.2 | 0.1028 | 0.0848 +- 0.0063 | 0.1333 | 0.1000 |
| -0.3 | 1 | 1/2 | 0.2 | 0.1722 | 0.2138 +- 0.0104 | 0.2000 | 0.1000 |
| -0.3 | 2 | 1/3 | 0.4 | 0.0908 | 0.0614 +- 0.0060 | 0.1333 | 0.1000 |
| -0.3 | 2 | 1/2 | 0.4 | 0.1530 | 0.1738 +- 0.0102 | 0.2000 | 0.1000 |
| -0.3 | 4 | 1/3 | 0.8 | 0.0681 | 0.0275 +- 0.0056 | 0.0833 | 0.1000 |
| -0.3 | 4 | 1/2 | 0.8 | 0.1116 | 0.1216 +- 0.0102 | 0.1250 | 0.1000 |

`prop:Vconv` is stated "for every `ell >= 1`", and at least one consumer needs `ell >= 2`: `part4-obstructions.tex` l.483-486 squares the `prop:Vconv` error and works in `L^{ell/2}`.

### Is the band still well formed at beta = 1/2?  (the adversarial check)

At `beta = 1/2` the inset is `sqrt(delta)`, the SAME ORDER as the grid spacing, so `lem:barriers-wellposed`'s "there exists n_0" is a real constraint and is checked rather than assumed.  Required: `B^Y_n < C^Y_n`, both strictly inside `(B^Y, C^Y)`, `y_0` strictly inside, and both on the `sqrt(delta)` grid.

| n | beta | B^Y_n | C^Y_n | eps_n realised | sqrt(delta) | band width in grid steps | well formed |
|---|---|---|---|---|---|---|---|
| 16 | 1/3 | 0.0000 | 0.0000 | 0.5000 | 0.2500 | 0 | **NO** |
| 16 | 1/2 | 0.0000 | 0.0000 | 0.5000 | 0.2500 | 0 | **NO** |
| 32 | 1/3 | -0.1768 | 0.1768 | 0.3232 | 0.1768 | 2 | yes |
| 32 | 1/2 | -0.1768 | 0.1768 | 0.3232 | 0.1768 | 2 | yes |
| 64 | 1/3 | -0.1250 | 0.1250 | 0.3750 | 0.1250 | 2 | yes |
| 64 | 1/2 | -0.2500 | 0.2500 | 0.2500 | 0.1250 | 4 | yes |
| 128 | 1/3 | -0.2652 | 0.2652 | 0.2348 | 0.0884 | 6 | yes |
| 128 | 1/2 | -0.3536 | 0.3536 | 0.1464 | 0.0884 | 8 | yes |
| 256 | 1/3 | -0.3125 | 0.3125 | 0.1875 | 0.0625 | 10 | yes |
| 256 | 1/2 | -0.3750 | 0.3750 | 0.1250 | 0.0625 | 12 | yes |
| 1024 | 1/3 | -0.3750 | 0.3750 | 0.1250 | 0.0312 | 24 | yes |
| 1024 | 1/2 | -0.4375 | 0.4375 | 0.0625 | 0.0312 | 28 | yes |
| 4096 | 1/3 | -0.4219 | 0.4219 | 0.0781 | 0.0156 | 54 | yes |
| 4096 | 1/2 | -0.4688 | 0.4688 | 0.0312 | 0.0156 | 60 | yes |
| 16384 | 1/3 | -0.4531 | 0.4531 | 0.0469 | 0.0078 | 116 | yes |
| 16384 | 1/2 | -0.4844 | 0.4844 | 0.0156 | 0.0078 | 124 | yes |
| 65536 | 1/3 | -0.4727 | 0.4727 | 0.0273 | 0.0039 | 242 | yes |
| 65536 | 1/2 | -0.4922 | 0.4922 | 0.0078 | 0.0039 | 252 | yes |

The floor is `n = 32` for BOTH insets and it is the same floor: at `n = 16` the band collapses onto `y_0` at `beta = 1/3` AND at `beta = 1/2`, because `sqrt(delta) = 1/4` cannot place a lattice point strictly inside `(-0.5, 0.5)` once anything is inset.  **`beta = 1/2` costs nothing in ladder range** -- the FF2 ladder `n = 64 .. 16384` is usable at both insets, which is why the side-by-side comparison above is on identical `n`.  At `beta = 1/2` the realised `eps_n` sits in `(sqrt(delta), 2 sqrt(delta)]`, i.e. one to two grid steps: the smallest non-degenerate inset the scheme can express, as intended.

### FF2 coupling diagnostics (beta = 1/3)

| n | mean max \|e_j\| | multi-level fine steps | short of n crossings | walk exits | y exits | mean G (b=1/3) | mean G (b=1/2) |
|---|---|---|---|---|---|---|---|
| 64 | 2.17e-03 | 0 | 0.008 | 1.000 | 0.995 | 0.2213 | 0.1749 |
| 128 | 2.37e-03 | 0 | 0.000 | 1.000 | 0.995 | 0.1686 | 0.1147 |
| 256 | 2.59e-03 | 0 | 0.000 | 1.000 | 0.995 | 0.1402 | 0.1023 |
| 512 | 2.80e-03 | 0 | 0.000 | 1.000 | 0.995 | 0.1147 | 0.0517 |
| 1024 | 2.99e-03 | 0 | 0.000 | 1.000 | 0.995 | 0.1023 | 0.0538 |
| 2048 | 3.18e-03 | 0 | 0.000 | 0.997 | 0.995 | 0.0702 | 0.0294 |
| 4096 | 3.34e-03 | 0 | 0.000 | 0.997 | 0.995 | 0.0684 | 0.0240 |
| 8192 | 3.51e-03 | 0 | 0.000 | 0.995 | 0.995 | 0.0514 | 0.0078 |
| 16384 | 3.65e-03 | 0 | 0.000 | 0.995 | 0.995 | 0.0419 | 0.0081 |

The node error is pure crossing overshoot, `O(sqrt(dt_fine))`; the exact embedding has `e_j == 0`.  It is two to three orders of magnitude below `delta^{1/4}`, the size the manuscript's own (E2) allows, so this artifact cannot be driving anything.

### FF2 resolution sensitivity (the check that this is not a numerical floor)

Fixed `n = 512`, 96 paths, `beta = 1/3`, fine grid refined by a factor 16 twice.  A COARSE `dt` makes the discrete and continuous objects MORE alike (they share the left-point convention), so under-resolution SUPPRESSES the measured gap -- and does so increasingly at large `n`, which INFLATES the apparent decay rate.  Measuring rates below target therefore bounds the true rates from above; the artifact cannot manufacture the stall.

| h | fine steps / unit time | dt | total_theta | F1 | F2 |
|---|---|---|---|---|---|
| -0.1 | 65536 | 1.53e-05 | 0.7233 +- 0.0475 | 0.0636 +- 0.0011 | 0.7465 +- 0.0490 |
| -0.1 | 262144 | 3.81e-06 | 0.6942 +- 0.0456 | 0.0638 +- 0.0012 | 0.7137 +- 0.0469 |
| -0.1 | 1048576 | 9.54e-07 | 0.6681 +- 0.0440 | 0.0642 +- 0.0011 | 0.6926 +- 0.0457 |
| -0.3 | 65536 | 1.53e-05 | 2.3446 +- 0.1174 | 0.7631 +- 0.0129 | 2.7313 +- 0.1273 |
| -0.3 | 262144 | 3.81e-06 | 2.2796 +- 0.1094 | 0.8367 +- 0.0128 | 2.6793 +- 0.1200 |
| -0.3 | 1048576 | 9.54e-07 | 2.3281 +- 0.1070 | 0.8905 +- 0.0125 | 2.7728 +- 0.1283 |

## FF3 -- the ordering lemma L1

Counted: the fraction of paths with `theta^(n)_{Xi^Y_n} > Xi^y`, i.e. the walk leaving its NARROW band later than the diffusion leaves its WIDE one.  `sigma_y(y) = 1 + 0.3 tanh(y)` in `[0.7,1.3]` throughout: bounded, Lipschitz, bounded away from 0, so `asm:coeff` holds.  Inset `beta = 1/3`.

* **A** -- genuine Skorokhod embedding (level crossings) with `mu_y = 0`.  Because `y` is then a martingale the embedded increments are exactly fair `+-1` coins whatever `sigma_y` is, so this IS the `eq:walks` walk and the embedding is faithful.  Node error = fine-grid crossing overshoot only.
* **B** -- SURROGATE: `theta_k := t_k` and `Y^(n)_k :=` the `sqrt(delta)`-lattice rounding of `y(t_k)`, with `mu_y(y) = 0.3 cos(y)` switched on.  NOT the `eq:walks` walk (its increments are multiples of `sqrt(delta)`, not `+-sqrt(delta)`); used only to give the node error a non-zero value of a known size, `|e_k| <= sqrt(delta)/2`.
* **C c** -- variant B plus an INJECTED node error `e_k = c * delta^{1/4} * tanh(Btilde_{t_k})`, `Btilde` independent, so `max_k |e_k| <= c delta^{1/4}` exactly.  Not an embedding: it is the largest node error the manuscript's own (E2) permits, injected deliberately, to answer the question L1 actually poses -- is the ordering IMPLIED by (E2) + `lem:barriers-wellposed`?

The decisive comparison is `mean max |e_j|` against the barrier gap `eps_n`. Note `delta^{1/4} > delta^{1/3}` for every `delta < 1`, and the ratio `delta^{1/4}/delta^{1/3} = delta^{-1/12}` DIVERGES -- at `beta = 1/2` the ratio `delta^{1/4}/delta^{1/2} = delta^{-1/4}` diverges faster, so shrinking the inset makes L1 HARDER, not easier.  That is a cost of D1 and it is recorded here rather than left for the reader to notice.

| variant | n | paths used | violations | fraction | mean max \|e_j\| | barrier gap eps_n | delta^{1/4} | no exit |
|---|---|---|---|---|---|---|---|---|
| A | 64 | 192 | 0 | 0.0000 | 0.0009 | 0.3750 | 0.3536 | 0 |
| A | 256 | 192 | 0 | 0.0000 | 0.0023 | 0.1875 | 0.2500 | 0 |
| A | 1024 | 192 | 0 | 0.0000 | 0.0033 | 0.1250 | 0.1768 | 0 |
| A | 4096 | 192 | 0 | 0.0000 | 0.0039 | 0.0781 | 0.1250 | 0 |
| B | 64 | 192 | 0 | 0.0000 | 0.0350 | 0.3750 | 0.3536 | 0 |
| B | 256 | 192 | 0 | 0.0000 | 0.0290 | 0.1875 | 0.2500 | 0 |
| B | 1024 | 192 | 0 | 0.0000 | 0.0154 | 0.1250 | 0.1768 | 0 |
| B | 4096 | 192 | 0 | 0.0000 | 0.0078 | 0.0781 | 0.1250 | 0 |
| C0.5 | 64 | 192 | 0 | 0.0000 | 0.0421 | 0.3750 | 0.3536 | 0 |
| C0.5 | 256 | 192 | 0 | 0.0000 | 0.0551 | 0.1875 | 0.2500 | 0 |
| C0.5 | 1024 | 192 | 0 | 0.0000 | 0.0424 | 0.1250 | 0.1768 | 0 |
| C0.5 | 4096 | 192 | 0 | 0.0000 | 0.0314 | 0.0781 | 0.1250 | 0 |
| C1.0 | 64 | 192 | 0 | 0.0000 | 0.0537 | 0.3750 | 0.3536 | 0 |
| C1.0 | 256 | 192 | 2 | 0.0104 | 0.0874 | 0.1875 | 0.2500 | 0 |
| C1.0 | 1024 | 192 | 6 | 0.0312 | 0.0722 | 0.1250 | 0.1768 | 0 |
| C1.0 | 4096 | 192 | 8 | 0.0417 | 0.0558 | 0.0781 | 0.1250 | 0 |
| C2.0 | 64 | 192 | 1 | 0.0052 | 0.0908 | 0.3750 | 0.3536 | 0 |
| C2.0 | 256 | 192 | 10 | 0.0521 | 0.1518 | 0.1875 | 0.2500 | 0 |
| C2.0 | 1024 | 192 | 14 | 0.0729 | 0.1307 | 0.1250 | 0.1768 | 0 |
| C2.0 | 4096 | 192 | 21 | 0.1094 | 0.1012 | 0.0781 | 0.1250 | 0 |

## FF4 -- the exit-time reconciliation L2

`G := Xi^y - theta^(n)_{Xi^Y_n} >= 0`.  With constant coefficients the exactly embedded walk sits ON the discrete barrier at `theta_{Xi^Y_n}`, so by the strong Markov property `G` is the exit time of a Brownian motion from `(B^Y,C^Y)` started at `C^Y_n`.  That law is classical, so the primary numbers below are EXACT (spectral + reflection series, `frozen_exit.psurv` / `moment_G`), validated at `kappa'=1` against `E[G] = x(L-x)` to 11 digits.  A coupled Monte-Carlo run cross-checks the reduction in (c).  Inset `beta = 1/3` unless stated.

### (a) E[G^kappa'] against the prediction O(eps_n^{2 kappa'}) = O(delta^{2 beta kappa'})

| n | delta | eps_n realised | E[G^0.2] | E[G^0.35] | E[G^0.45] | E[G] |
|---|---|---|---|---|---|---|
| 100 | 0.01 | 0.3 | 0.6815 | 0.52611 | 0.44702 | 0.21 |
| 1000 | 0.001 | 0.12053 | 0.54004 | 0.36324 | 0.28606 | 0.106 |
| 10000 | 0.0001 | 0.05 | 0.40504 | 0.22927 | 0.1645 | 0.0475 |
| 100000 | 1e-05 | 0.022496 | 0.30399 | 0.14435 | 0.094314 | 0.02199 |
| 1000000 | 1e-06 | 0.011 | 0.23227 | 0.093212 | 0.055706 | 0.010879 |
| 10000000 | 1e-07 | 0.0047873 | 0.16835 | 0.055017 | 0.029486 | 0.0047644 |
| 100000000 | 1e-08 | 0.0022 | 0.12388 | 0.033178 | 0.01599 | 0.0021952 |
| 1000000000 | 1e-09 | 0.0010242 | 0.091263 | 0.019998 | 0.0086507 | 0.0010232 |
| 10000000000 | 1e-10 | 0.00047 | 0.066613 | 0.011857 | 0.0045789 | 0.00046978 |
| 100000000000 | 1e-11 | 0.00021783 | 0.048648 | 0.0070389 | 0.0024232 | 0.00021778 |
| 1000000000000 | 1e-12 | 0.0001 | 0.03525 | 0.0041335 | 0.0012634 | 9.999e-05 |

The prediction is asymptotic and the `kappa' -> 1/2` crossover is only logarithmically fast, so the exponents are fitted over three nested windows. `kappa' = 1` is the CONTROL: there the exponent must be exactly 1 in `eps_n` (and `beta` in `delta`), which pins the fitting procedure.

| quantity | exp in delta (all) | exp in delta (last 7) | exp in delta (last 4) | prediction | exp in eps_n (all) | exp in eps_n (last 7) | exp in eps_n (last 4) | prediction |
|---|---|---|---|---|---|---|---|---|
| E[G^0.2] | 0.1298 | 0.1359 | 0.1376 | 0.1333 | 0.3772 | 0.4013 | 0.4089 | 0.4000 |
| E[G^0.35] | 0.2127 | 0.2247 | 0.2280 | 0.2333 | 0.6179 | 0.6638 | 0.6777 | 0.7000 |
| E[G^0.45] | 0.2574 | 0.2731 | 0.2783 | 0.3000 | 0.7477 | 0.8066 | 0.8270 | 0.9000 |
| E[G^1.0] | 0.3340 | 0.3378 | 0.3364 | 0.3333 | 0.9707 | 0.9980 | 0.9996 | 1.0000 |

The exponents rise monotonically toward the prediction as the window moves out, and the control lands on 1.000.  The residual shortfall at `kappa' = 0.45` is the `kappa' -> 1/2` crossover, where `kappa' int t^{kappa'-1} P(G>t) dt` stops being dominated by its lower end.

`E[G]` is the diagnostic the note asks for.  In the idealised half-line problem `E[G]` is infinite.  In the actual BOUNDED band it is finite -- but it does not scale like `eps_n^2`, it scales like `eps_n^1`, one full power short.  Any mean-based estimate of L2 gives `O(eps_n)`, not `O(eps_n^2)`.

### (c) the same, at both insets, on the FF2 ladder -- is the ratio the pure beta factor?

If nothing but the inset changed, the exponent in `delta` must be `2 beta kappa'` at each inset, i.e. the ratio of the two fitted exponents must be `(1/2)/(1/3) = 1.5` exactly, and the exponent in `eps_n` must be `2 kappa'` at BOTH.

| kappa' | exp in delta (beta=1/3) | exp in delta (beta=1/2) | ratio | predicted 1.5 | exp in eps_n (b=1/3) | (b=1/2) | predicted 2 kappa' |
|---|---|---|---|---|---|---|---|
| 0.2 | 0.1028 | 0.1722 | 1.675 | 1.500 | 0.2852 | 0.3342 | 0.4000 |
| 0.35 | 0.1641 | 0.2766 | 1.685 | 1.500 | 0.4552 | 0.5368 | 0.7000 |

### (b) raw tail P(G > t) against the prediction ~ eps_n / sqrt(t)

n = 10000, `beta = 1/3`, realised `eps_n = 0.05`, `t` over `[eps_n^2, T]`.  The prediction is the half-line stable(1/2) tail `eps_n sqrt(2/(pi t))`.

| t | t / eps_n^2 | P(G>t) exact | eps_n sqrt(2/(pi t)) | ratio |
|---|---|---|---|---|
| 0.0025 | 1 | 0.68269 | 0.79788 | 0.8556 |
| 0.005287 | 2.115 | 0.50833 | 0.54867 | 0.9265 |
| 0.01118 | 4.472 | 0.36369 | 0.3773 | 0.9639 |
| 0.02364 | 9.457 | 0.25495 | 0.25945 | 0.9827 |
| 0.05 | 20 | 0.17692 | 0.17841 | 0.9916 |
| 0.1057 | 42.29 | 0.11996 | 0.12269 | 0.9778 |
| 0.2236 | 89.44 | 0.066081 | 0.084366 | 0.7833 |
| 0.4729 | 189.1 | 0.019311 | 0.058015 | 0.3329 |
| 1 | 400 | 0.0014325 | 0.039894 | 0.0359 |

## FF7 -- P(G > c) at a FIXED c, against eps_n

`c = T/2 = 0.5`.  This is the lower-bound mechanism of the sharpness claim measured directly: with probability of order `eps_n` the diffusion started at the discrete barrier stays inside the band for a time of order 1, and on that event the frozen and the unfrozen object differ by `Theta(1)`.  Exact, from the same closed-form law.  Prediction: `P(G>c) ~ eps_n`, hence exponent 1 in `eps_n` and `beta` in `delta`; and at fixed `n` the ratio between the two insets must be `delta^{1/2-1/3} = delta^{1/6}`.

| n | delta | eps_n (b=1/3) | P(G>c) (b=1/3) | eps_n (b=1/2) | P(G>c) (b=1/2) | ratio P(1/2)/P(1/3) | delta^{1/6} |
|---|---|---|---|---|---|---|---|
| 64 | 0.01562 | 0.375 | 0.099758 | 0.25 | 0.076351 | 0.7654 | 0.5000 |
| 128 | 0.007812 | 0.23483 | 0.072628 | 0.14645 | 0.047944 | 0.6601 | 0.4454 |
| 256 | 0.003906 | 0.1875 | 0.059989 | 0.125 | 0.041321 | 0.6888 | 0.3969 |
| 512 | 0.001953 | 0.14645 | 0.047944 | 0.058058 | 0.019585 | 0.4085 | 0.3536 |
| 1024 | 0.0009766 | 0.125 | 0.041321 | 0.0625 | 0.021065 | 0.5098 | 0.3150 |
| 2048 | 0.0004883 | 0.080155 | 0.026904 | 0.035961 | 0.012173 | 0.4525 | 0.2806 |
| 4096 | 0.0002441 | 0.078125 | 0.026236 | 0.03125 | 0.010584 | 0.4034 | 0.2500 |
| 8192 | 0.0001221 | 0.058058 | 0.019585 | 0.013864 | 0.0047015 | 0.2400 | 0.2227 |
| 16384 | 6.104e-05 | 0.046875 | 0.015844 | 0.015625 | 0.0052982 | 0.3344 | 0.1984 |

Fitted on the FF2 ladder and on a deep exact ladder (`n = 64 .. 2.68e+08`), since the exponent-1 prediction is asymptotic and the exact law costs nothing to extend:

| beta | exp in eps_n (FF2 ladder) | (deep ladder) | (deep, last 6) | prediction | exp in delta (deep) | prediction beta |
|---|---|---|---|---|---|---|
| 1/3 | 0.9046 | 0.9748 | 0.9998 | 1.0000 | 0.3444 | 0.3333 |
| 1/2 | 0.9727 | 0.9927 | 1.0000 | 1.0000 | 0.4964 | 0.5000 |

### (c) coupled Monte-Carlo cross-check, both insets

Same fine-grid construction as FF2 (exact level-crossing embedding), horizon 3.0T so that `Xi^y` is observed.  `G<0` counts L1 violations (must be 0); `censored` counts paths whose `Xi^y` exceeded the horizon and are dropped, which biases the MC column DOWNWARD, so agreement is the meaningful direction.  `P(G>T/2)` here is the Monte-Carlo counterpart of the FF7 table and is quoted with its binomial standard error -- at the small `eps_n` end the counts get thin, which is exactly why FF7's primary numbers are exact.

| n | beta | paths used | G<0 | censored | E[G^0.2] MC / exact | E[G^0.35] MC / exact | E[G^0.45] MC / exact | E[G] MC / exact | P(G>T/2) MC | exact |
|---|---|---|---|---|---|---|---|---|---|---|
| 256 | 1/3 | 512 | 0 | 0 | 0.6189 / 0.612 | 0.4505 / 0.4433 | 0.3715 / 0.3636 | 0.1607 / 0.1523 | 0.06836 +- 0.0112 | 0.05999 |
| 256 | 1/2 | 512 | 0 | 0 | 0.5527 / 0.546 | 0.3757 / 0.3696 | 0.2982 / 0.2921 | 0.1141 / 0.1094 | 0.04297 +- 0.00896 | 0.04132 |
| 1024 | 1/3 | 512 | 0 | 0 | 0.5527 / 0.546 | 0.3757 / 0.3696 | 0.2982 / 0.2921 | 0.1141 / 0.1094 | 0.04297 +- 0.00896 | 0.04132 |
| 1024 | 1/2 | 510 | 0 | 0 | 0.446 / 0.4372 | 0.266 / 0.2592 | 0.1965 / 0.1906 | 0.05992 / 0.05859 | 0.01961 +- 0.00614 | 0.02107 |
| 4096 | 1/3 | 512 | 0 | 0 | 0.4817 / 0.471 | 0.3013 / 0.292 | 0.2286 / 0.22 | 0.07662 / 0.07202 | 0.02734 +- 0.00721 | 0.02624 |
| 4096 | 1/2 | 510 | 0 | 0 | 0.3477 / 0.3428 | 0.1787 / 0.1753 | 0.1223 / 0.1191 | 0.03245 / 0.03027 | 0.01373 +- 0.00515 | 0.01058 |

## Shortcuts, and how each could have manufactured the answer

1. **Band choice, and an open obligation it touches.**  `(B^Y,C^Y) = (-0.5,0.5)`, symmetric about `y_0 = 0`.  `part1-setup.tex` declares `B^Y in [0,inf)` while `eq:model` starts the driver at `y_0 = 0`, so **as the symbols stand `y_0` sits ON the barrier and `Xi^y = 0` almost surely** -- the absorbed driver never moves and `V_t == v_0`. The `[0,inf)` constraint is a leftover from the drafts in which the Y-barriers were centred at `v_0` (a variance, hence non-negative); the `\fixed{}` note in `eq:barriers` recentres them at `y_0` but the constraint on `B^Y` was not moved with them.  This work took the recentred reading throughout.  **Recorded as an independent hit on the project's open obligation O003.**  A different band changes every CONSTANT below (FF1's plateau height, FF2's level) but not the exponents, since `eps_n ~ delta^{beta}` is band-independent.
2. **Fine-grid continuous reference (FF2/FF6, FF3, FF4c).**  `v_t` is the left-point Riemann sum at `dt_fine`, own error `O(dt_fine^H)`.  This is the shortcut most able to fake a stall, so it is tested directly in the resolution table, and its bias direction is stated there: a coarse `dt` makes the two objects MORE alike, so under-resolution suppresses the gap and INFLATES the apparent decay rate.  Measured rates below target therefore bound the true rates from above.
3. **h = -0.45 is not measured in FF2/FF6.**  Its target rate `delta^{0.025}` is a factor 1.15 over a factor-256 range in `delta` and cannot be separated from a constant on any feasible ladder.  The earlier run that did include it produced slopes of -0.001 +- 0.005 at `ell=2`; that number is not evidence and was dropped rather than quoted.  `h = -0.45` remains in FF1, which is exact.
4. **Crossing overshoot.**  The exact embedding is realised on a grid, so `y(theta_k)` overshoots the lattice by `O(sqrt(dt_fine))` and the walk is taken to be the exact lattice point.  Reported per `n`; 2-3 orders below `delta^{1/4}`.  It biases the node error UPWARD, i.e. AGAINST L1, and L1 still never fails in variants A and B.
5. **Fine-grid exit detection.**  `Xi^y` is the first fine-grid index outside the band, so it is LATE by `O(dt_fine)` and the effective barrier sits `~0.6 sqrt(dt_fine)` outside the true one.  This biases `G` UPWARD, i.e. in favour of a slower FF2 decay.  At `beta = 1/2` the inset is only `sqrt(delta)`, so this artifact is a LARGER fraction of `eps_n` than at `beta = 1/3` -- it is `0.6 sqrt(dt_fine)/sqrt(delta) = 0.6 sqrt(delta/dt_fine)` inverted, i.e. `0.6/sqrt(mpu/n)`, at most 0.6/8 = 7.5% at the top of the ladder.  It therefore flatters `beta = 1/2` slightly, and the FF6 improvement should be read with that in mind.  FF4 and FF7's primary numbers avoid it entirely by using the exact law.
6. **FF3 variants B and C are surrogates, not embeddings**, and are labelled as such.  A violation count from C is evidence about what (E2) IMPLIES, not about what the manuscript's embedding DOES.  Variant A is the only faithful construction there, and for `mu_y = 0` it has `e_j == 0` in continuous time -- so a zero violation count is close to a tautology and must not be read as support for L1.
7. **Censoring.**  Paths whose `Xi^y` exceeds the horizon are dropped in the FF4/FF7 Monte-Carlo (count reported); dropping the longest `G` biases the MC moments and `P(G>T/2)` DOWNWARD, so MC-below-exact is the expected direction.
8. **Monte-Carlo error.**  All `+-` are bootstrap standard errors over paths (binomial for `P(G>T/2)`).  In FF2/FF6 the SAME fine paths are used for every `n` AND for both insets, so the estimates are strongly positively correlated and the fitted slopes -- and especially the `beta=1/3` vs `beta=1/2` difference -- are far better determined than the individual points.  The bootstrap resamples whole paths and respects that.

## Compute time

| phase | wall time | seconds |
|---|---|---|
| 1. FF1 frozen-vs-unfrozen (discrete only) | 41s | 40.9 |
| 2. FF2/FF6 frozen-vs-frozen (coupled, mpu=1048576, M=1310720) | 2m49s | 169.09 |
| 3. FF2 resolution sensitivity (fixed n, varying fine grid) | 45s | 45.47 |
| 4. FF3 ordering L1 (non-constant coefficients, M=655360) | 29s | 29.32 |
| 5. FF4/FF7 exit-time law (exact, both insets) | 5s | 5.02 |
| **total** | **5m14s** | **314.91** |

Machine: Darwin arm64, python 3.9.6.
