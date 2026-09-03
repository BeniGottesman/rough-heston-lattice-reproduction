# Route F falsifiers FF1-FF4 -- measured

Route and falsifiers: `research-notes/L003-ROUTE-F.md` sec.2-3.
Code: `sim/frozen_exit.py`; driver: `sim/run_frozen_exit.py`.

Common settings: `T=1`, `y_0=0`, `v_0=0`, `K(u)=u^h` (`L==1`), continuous band `(B^Y,C^Y)=(-0.5,0.5)`, `mu_y=0`, `sigma_y=1` except in FF3.  `H := h+1/2`.  Discrete barriers from `eq:barriers` verbatim.

These are measurements, not proofs.  Nothing here closes a Definition-of-Done item; a refutation is an instruction to stop, a non-refutation is not a permission to continue.

## Decision table

Verdicts are the author's reading of the numbers below; the numbers are the evidence and the tables are where they live.

| falsifier | verdict | the number that decides it |
|---|---|---|
| FF1 (the diagnosis) | **CONFIRMED** | at h=-0.3, n=10^6 the L2 gap is 4.154 against delta^H = 0.0631, a ratio of 66, and it GROWS with n; 49% of post-exit increments are non-zero, so the remark's premise is false outright |
| FF2 (Route F itself) | **REFUTED** | the frozen-vs-frozen decay misses the target (h+kappa)/2 by more than 2 bootstrap sd in 9 of 9 (h, ell) cells: h=-0.1,ell=1, h=-0.1,ell=2, h=-0.1,ell=4, h=-0.3,ell=1, h=-0.3,ell=2, h=-0.3,ell=4, h=-0.45,ell=1, h=-0.45,ell=2, h=-0.45,ell=4.  F1 meets its rate; F2 does not, and its exponent falls with ell exactly as the closed-form law of G predicts |
| FF3 (ordering L1) | **NOT REFUTED as a fact, REFUTED as an inference** | 0 violations in 1535 faithful-construction paths; but with a node error of the size (E2) permits the fraction runs 0.5% -> 5.7% -> 8.9% -> 17.7% over n = 64, 256, 1024, 4096 -- rising, not vanishing |
| FF4 (exit-time L2) | **NOT REFUTED** | on an exact ladder to n = 10^12 the exponent of E[G^0.2] in eps_n converges to 0.409 against the predicted 0.4, with the kappa'=1 control landing on 1.000; the tail matches eps_n sqrt(2/(pi t)) to within 3% over its stated range |

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

## FF2 -- frozen discrete against FROZEN continuous (coupled)

**Coupling.**  The exact Skorokhod level-crossing embedding `theta_k = inf{t > theta_{k-1} : |y_t - y_{theta_{k-1}}| = sqrt(delta)}`, realised on a fine grid of step `dt = 9.54e-07`.  With `mu_y=0, sigma_y=1` this is the manuscript's OWN embedding and it is exact: `Y(theta_k) = Y^(n)_k` and the node error `e_j` vanishes (`part3-embedding.tex` l.582).  Because `y` is a martingale, the embedded increments are i.i.d. fair `+-1` coins, so the walk is the `eq:walks` walk in law -- the coupling does not deform either object.  It was chosen over sign-coupling (`zeta_k = sign(y(t_k)-y(t_{k-1}))`, which leaves `|Y^(n)_k - y(t_k)| = Theta(1)` and would fake a stall) and over lattice-rounding of `y(t_k)` (which is not the `eq:walks` walk).

**Decomposition measured** (Route F sec.2, `nu := k ^ Xi^Y_n`):

```
  total = max_k |V^(n)_k - V_{theta_k}|
  F1    = max_{j<=Xi^Y_n} |Vcal^(n)(t_j) - v_{theta_j}|          scheme error
  F2    = max_{k>Xi^Y_n} |v_{theta_k ^ Xi^y} - v_{theta_{Xi}}|   exit-time term
```

`prop:Vconv` is claimed for EVERY `ell >= 1`, so every norm is reported in `L^1`, `L^2` and `L^4`.

### h = -0.1  (H = 0.4; target exponent (h+kappa)/2 -> 0.2 as kappa -> 1/2)

| n | delta | total_theta L2 | total_tk L2 | F1 L2 | F2 L2 |
|---|---|---|---|---|---|
| 64 | 0.01562 | 0.7891 +- 0.0122 | 0.7852 +- 0.0123 | 0.0691 +- 0.0007 | 0.8115 +- 0.0158 |
| 128 | 0.007812 | 0.7297 +- 0.0229 | 0.7321 +- 0.0226 | 0.0784 +- 0.0009 | 0.7636 +- 0.0250 |
| 256 | 0.003906 | 0.6919 +- 0.0247 | 0.7073 +- 0.0239 | 0.0720 +- 0.0009 | 0.7186 +- 0.0259 |
| 512 | 0.001953 | 0.6434 +- 0.0287 | 0.6805 +- 0.0256 | 0.0631 +- 0.0007 | 0.6643 +- 0.0297 |
| 1024 | 0.0009766 | 0.6284 +- 0.0299 | 0.6677 +- 0.0266 | 0.0537 +- 0.0006 | 0.6439 +- 0.0304 |
| 2048 | 0.0004883 | 0.5970 +- 0.0325 | 0.6676 +- 0.0273 | 0.0452 +- 0.0004 | 0.6096 +- 0.0329 |
| 4096 | 0.0002441 | 0.5930 +- 0.0334 | 0.6903 +- 0.0267 | 0.0377 +- 0.0003 | 0.6022 +- 0.0337 |
| 8192 | 0.0001221 | 0.5327 +- 0.0356 | 0.7006 +- 0.0244 | 0.0318 +- 0.0003 | 0.5369 +- 0.0357 |
| 16384 | 6.104e-05 | 0.5174 +- 0.0363 | 0.7422 +- 0.0221 | 0.0280 +- 0.0002 | 0.5185 +- 0.0364 |

### h = -0.3  (H = 0.2; target exponent (h+kappa)/2 -> 0.1 as kappa -> 1/2)

| n | delta | total_theta L2 | total_tk L2 | F1 L2 | F2 L2 |
|---|---|---|---|---|---|
| 64 | 0.01562 | 2.4027 +- 0.0306 | 2.3307 +- 0.0322 | 0.8452 +- 0.0090 | 2.9473 +- 0.0558 |
| 128 | 0.007812 | 2.3452 +- 0.0519 | 2.1799 +- 0.0534 | 0.9589 +- 0.0090 | 2.9272 +- 0.0683 |
| 256 | 0.003906 | 2.2992 +- 0.0564 | 2.1669 +- 0.0566 | 0.9348 +- 0.0091 | 2.8386 +- 0.0659 |
| 512 | 0.001953 | 2.2393 +- 0.0659 | 2.1612 +- 0.0613 | 0.8819 +- 0.0074 | 2.6758 +- 0.0743 |
| 1024 | 0.0009766 | 2.2439 +- 0.0714 | 2.2105 +- 0.0645 | 0.8138 +- 0.0063 | 2.6169 +- 0.0771 |
| 2048 | 0.0004883 | 2.1754 +- 0.0792 | 2.2733 +- 0.0641 | 0.7389 +- 0.0054 | 2.5011 +- 0.0830 |
| 4096 | 0.0002441 | 2.1966 +- 0.0824 | 2.4147 +- 0.0637 | 0.6646 +- 0.0047 | 2.4683 +- 0.0856 |
| 8192 | 0.0001221 | 2.0678 +- 0.0896 | 2.5430 +- 0.0585 | 0.5878 +- 0.0042 | 2.2551 +- 0.0910 |
| 16384 | 6.104e-05 | 2.0315 +- 0.0921 | 2.7345 +- 0.0524 | 0.5207 +- 0.0033 | 2.1661 +- 0.0933 |

### h = -0.45  (H = 0.05; target exponent (h+kappa)/2 -> 0.025 as kappa -> 1/2)

| n | delta | total_theta L2 | total_tk L2 | F1 L2 | F2 L2 |
|---|---|---|---|---|---|
| 64 | 0.01562 | 7.4787 +- 0.0699 | 6.9147 +- 0.0786 | 4.2824 +- 0.0459 | 10.7551 +- 0.1695 |
| 128 | 0.007812 | 7.5040 +- 0.1045 | 6.4616 +- 0.1089 | 4.9552 +- 0.0435 | 10.6768 +- 0.1889 |
| 256 | 0.003906 | 7.4605 +- 0.1133 | 6.4111 +- 0.1187 | 5.0481 +- 0.0410 | 10.5433 +- 0.1774 |
| 512 | 0.001953 | 7.4824 +- 0.1336 | 6.4989 +- 0.1313 | 4.9889 +- 0.0359 | 10.1457 +- 0.1914 |
| 1024 | 0.0009766 | 7.6169 +- 0.1485 | 6.7447 +- 0.1372 | 4.8579 +- 0.0311 | 10.0907 +- 0.1924 |
| 2048 | 0.0004883 | 7.5477 +- 0.1721 | 7.0614 +- 0.1360 | 4.6513 +- 0.0305 | 9.8024 +- 0.2065 |
| 4096 | 0.0002441 | 7.6586 +- 0.1841 | 7.6542 +- 0.1380 | 4.3972 +- 0.0271 | 9.7275 +- 0.2127 |
| 8192 | 0.0001221 | 7.4916 +- 0.2049 | 8.2602 +- 0.1293 | 4.1160 +- 0.0258 | 9.1339 +- 0.2218 |
| 16384 | 6.104e-05 | 7.4991 +- 0.2158 | 8.9898 +- 0.1149 | 3.8011 +- 0.0220 | 8.8881 +- 0.2309 |

### Fitted decay exponents (slope of log ||.||_{L^ell} against log delta)

`refined prediction` is what the route's own heuristic gives once the norm exponent is carried through the heavy tail: `||F2||_{L^ell} ~ E[G^{ell H}]^{1/ell}`, and `E[G^{kappa'}] ~ eps_n^{min(2 kappa', 1)}`, so the exponent in `delta` is `(1/3) min(2H, 1/ell)`.  It equals the route's `2H/3` only while `ell H < 1/2`.

| h | quantity | ell | fitted exponent | bootstrap sd | 95% CI | route target (h+1/2)/2 | refined prediction |
|---|---|---|---|---|---|---|---|
| -0.1 | total_theta | 1 | 0.1468 | 0.0133 | [0.1221, 0.1737] | 0.2000 | 0.2667 |
| -0.1 | total_theta | 2 | 0.0736 | 0.0123 | [0.0500, 0.0990] | 0.2000 | 0.1667 |
| -0.1 | total_theta | 4 | 0.0072 | 0.0085 | [-0.0079, 0.0242] | 0.2000 | 0.0833 |
| -0.1 | total_tk | 1 | 0.0164 | 0.0050 | [0.0060, 0.0259] | 0.2000 | 0.2667 |
| -0.1 | total_tk | 2 | 0.0108 | 0.0058 | [-0.0003, 0.0217] | 0.2000 | 0.1667 |
| -0.1 | total_tk | 4 | -0.0052 | 0.0064 | [-0.0174, 0.0078] | 0.2000 | 0.0833 |
| -0.1 | F1 | 1 | 0.1889 | 0.0020 | [0.1846, 0.1928] | 0.2000 | -- |
| -0.1 | F1 | 2 | 0.1911 | 0.0020 | [0.1866, 0.1948] | 0.2000 | -- |
| -0.1 | F1 | 4 | 0.1952 | 0.0021 | [0.1906, 0.1992] | 0.2000 | -- |
| -0.1 | F2 | 1 | 0.1487 | 0.0132 | [0.1241, 0.1755] | 0.2000 | 0.2667 |
| -0.1 | F2 | 2 | 0.0801 | 0.0122 | [0.0568, 0.1059] | 0.2000 | 0.1667 |
| -0.1 | F2 | 4 | 0.0173 | 0.0085 | [0.0023, 0.0345] | 0.2000 | 0.0833 |
| -0.3 | total_theta | 1 | 0.0625 | 0.0082 | [0.0470, 0.0786] | 0.1000 | 0.1333 |
| -0.3 | total_theta | 2 | 0.0288 | 0.0080 | [0.0135, 0.0453] | 0.1000 | 0.1333 |
| -0.3 | total_theta | 4 | -0.0133 | 0.0066 | [-0.0250, 0.0017] | 0.1000 | 0.0833 |
| -0.3 | total_tk | 1 | -0.0323 | 0.0038 | [-0.0401, -0.0247] | 0.1000 | 0.1333 |
| -0.3 | total_tk | 2 | -0.0324 | 0.0040 | [-0.0400, -0.0240] | 0.1000 | 0.1333 |
| -0.3 | total_tk | 4 | -0.0354 | 0.0046 | [-0.0439, -0.0258] | 0.1000 | 0.0833 |
| -0.3 | F1 | 1 | 0.1011 | 0.0017 | [0.0977, 0.1045] | 0.1000 | -- |
| -0.3 | F1 | 2 | 0.1026 | 0.0017 | [0.0990, 0.1061] | 0.1000 | -- |
| -0.3 | F1 | 4 | 0.1057 | 0.0019 | [0.1017, 0.1094] | 0.1000 | -- |
| -0.3 | F2 | 1 | 0.0803 | 0.0075 | [0.0656, 0.0965] | 0.1000 | 0.1333 |
| -0.3 | F2 | 2 | 0.0573 | 0.0073 | [0.0435, 0.0731] | 0.1000 | 0.1333 |
| -0.3 | F2 | 4 | 0.0235 | 0.0064 | [0.0122, 0.0377] | 0.1000 | 0.0833 |
| -0.45 | total_theta | 1 | 0.0135 | 0.0052 | [0.0037, 0.0234] | 0.0250 | 0.0333 |
| -0.45 | total_theta | 2 | -0.0011 | 0.0053 | [-0.0111, 0.0094] | 0.0250 | 0.0333 |
| -0.45 | total_theta | 4 | -0.0253 | 0.0049 | [-0.0344, -0.0151] | 0.0250 | 0.0333 |
| -0.45 | total_tk | 1 | -0.0532 | 0.0031 | [-0.0591, -0.0473] | 0.0250 | 0.0333 |
| -0.45 | total_tk | 2 | -0.0530 | 0.0032 | [-0.0592, -0.0469] | 0.0250 | 0.0333 |
| -0.45 | total_tk | 4 | -0.0531 | 0.0034 | [-0.0595, -0.0466] | 0.0250 | 0.0333 |
| -0.45 | F1 | 1 | 0.0318 | 0.0016 | [0.0287, 0.0348] | 0.0250 | -- |
| -0.45 | F1 | 2 | 0.0333 | 0.0016 | [0.0300, 0.0364] | 0.0250 | -- |
| -0.45 | F1 | 4 | 0.0361 | 0.0018 | [0.0325, 0.0397] | 0.0250 | -- |
| -0.45 | F2 | 1 | 0.0386 | 0.0050 | [0.0288, 0.0495] | 0.0250 | 0.0333 |
| -0.45 | F2 | 2 | 0.0345 | 0.0045 | [0.0263, 0.0442] | 0.0250 | 0.0333 |
| -0.45 | F2 | 4 | 0.0211 | 0.0043 | [0.0131, 0.0299] | 0.0250 | 0.0333 |

### FF2 x FF4 -- where the route's arithmetic loses a power (exact, not simulated)

Route F sec.2 bounds `|F2| <= |v|_{H-} * G^{H-}` and then quotes `E[G^{kappa'}] = O(eps_n^{2 kappa'})`, valid *for kappa' < 1/2*, with `kappa' = H`.  But `prop:Vconv` is an `L^ell` statement, so what is needed is `|| |v|_H G^H ||_{L^ell} ~ E[G^{ell H}]^{1/ell}` -- the exponent inside the expectation is `ell H`, not `H`.  The `kappa' < 1/2` restriction therefore reads `ell H < 1/2`, and above it `E[G^{kappa'}] ~ eps_n^{1}` (not `eps_n^{2 kappa'}`), so the norm decays like `eps_n^{1/ell} = delta^{1/(3 ell)}` instead of `delta^{2H/3}`.  Comparing with the target `(h+kappa)/2 -> H/2`:

```
  route's F2 cost meets the target  <=>  (1/3) min(2H, 1/ell) >= H/2
                                    <=>  ell <= 2 / (3H).
```

The column `exact E[G^{ell H}]^{1/ell} slope` below is computed from the closed-form law of `G` (no simulation) over the SAME `n` ladder as the FF2 tables, so it is directly comparable with the measured `F2` slope next to it.

| h | H | ell | ell*H | ell limit 2/(3H) | exact E[G^{ell H}]^{1/ell} slope | measured F2 slope | route 2H/3 | refined min(2H,1/ell)/3 | target (h+1/2)/2 |
|---|---|---|---|---|---|---|---|---|---|
| -0.1 | 0.4 | 1 | 0.4 | 1.67 | 0.1816 | 0.1487 +- 0.0132 | 0.2667 | 0.2667 | 0.2000 |
| -0.1 | 0.4 | 2 | 0.8 | 1.67 | 0.1362 | 0.0801 +- 0.0122 | 0.2667 | 0.1667 | 0.2000 |
| -0.1 | 0.4 | 4 | 1.6 | 1.67 | 0.0796 | 0.0173 +- 0.0085 | 0.2667 | 0.0833 | 0.2000 |
| -0.3 | 0.2 | 1 | 0.2 | 3.33 | 0.1028 | 0.0803 +- 0.0075 | 0.1333 | 0.1333 | 0.1000 |
| -0.3 | 0.2 | 2 | 0.4 | 3.33 | 0.0908 | 0.0573 +- 0.0073 | 0.1333 | 0.1333 | 0.1000 |
| -0.3 | 0.2 | 4 | 0.8 | 3.33 | 0.0681 | 0.0235 +- 0.0064 | 0.1333 | 0.0833 | 0.1000 |
| -0.45 | 0.05 | 1 | 0.05 | 13.3 | 0.0362 | 0.0386 +- 0.0050 | 0.0333 | 0.0333 | 0.0250 |
| -0.45 | 0.05 | 2 | 0.1 | 13.3 | 0.0286 | 0.0345 +- 0.0045 | 0.0333 | 0.0333 | 0.0250 |
| -0.45 | 0.05 | 4 | 0.2 | 13.3 | 0.0257 | 0.0211 +- 0.0043 | 0.0333 | 0.0333 | 0.0250 |

`prop:Vconv` is stated "for every `ell >= 1`", and at least one consumer needs `ell >= 2`: `part4-obstructions.tex` l.483-486 squares the `prop:Vconv` error and works in `L^{ell/2}`.  So the failure at `ell = 2, h = -0.1` is inside the range the manuscript actually uses, and the failure for large `ell` is inside the range the proposition claims, at EVERY `h`.

### FF2 coupling diagnostics

| n | mean max node error \|e_j\| | multi-level fine steps | short of n crossings | walk exits | y exits | mean G |
|---|---|---|---|---|---|---|
| 64 | 2.16e-03 | 0 | 0.008 | 1.000 | 1.000 | 0.2132 |
| 128 | 2.36e-03 | 0 | 0.000 | 1.000 | 1.000 | 0.1610 |
| 256 | 2.60e-03 | 0 | 0.000 | 1.000 | 1.000 | 0.1317 |
| 512 | 2.79e-03 | 0 | 0.000 | 1.000 | 1.000 | 0.1032 |
| 1024 | 3.01e-03 | 0 | 0.000 | 1.000 | 1.000 | 0.0913 |
| 2048 | 3.17e-03 | 0 | 0.000 | 1.000 | 1.000 | 0.0680 |
| 4096 | 3.36e-03 | 0 | 0.000 | 1.000 | 1.000 | 0.0660 |
| 8192 | 3.49e-03 | 0 | 0.000 | 1.000 | 1.000 | 0.0504 |
| 16384 | 3.66e-03 | 0 | 0.000 | 1.000 | 1.000 | 0.0453 |

The node error is pure crossing overshoot, `O(sqrt(dt_fine))`; the exact embedding has `e_j == 0`.  Compare it with `delta^{1/4}`, the size the manuscript's own (E2) allows: it is two to three orders of magnitude smaller, so this artifact cannot be driving anything.

### FF2 resolution sensitivity (the check that this is not a numerical floor)

Fixed `n = 512`, 96 paths, fine grid refined by a factor 16 twice.  If the reported stall were the fine-grid quadrature error of the continuous convolution (`O(dt_fine^H)`), these columns would fall.

| h | fine steps / unit time | dt | total_theta | F1 | F2 |
|---|---|---|---|---|---|
| -0.1 | 65536 | 1.53e-05 | 0.7024 +- 0.0463 | 0.0638 +- 0.0011 | 0.7264 +- 0.0481 |
| -0.1 | 262144 | 3.81e-06 | 0.6707 +- 0.0509 | 0.0641 +- 0.0012 | 0.6902 +- 0.0528 |
| -0.1 | 1048576 | 9.54e-07 | 0.6565 +- 0.0461 | 0.0637 +- 0.0011 | 0.6811 +- 0.0477 |
| -0.3 | 65536 | 1.53e-05 | 2.3022 +- 0.1098 | 0.7542 +- 0.0114 | 2.6759 +- 0.1283 |
| -0.3 | 262144 | 3.81e-06 | 2.2621 +- 0.1224 | 0.8416 +- 0.0124 | 2.6642 +- 0.1368 |
| -0.3 | 1048576 | 9.54e-07 | 2.2708 +- 0.1111 | 0.8879 +- 0.0116 | 2.7809 +- 0.1217 |
| -0.45 | 65536 | 1.53e-05 | 6.6215 +- 0.2310 | 3.6373 +- 0.0500 | 8.5495 +- 0.3230 |
| -0.45 | 262144 | 3.81e-06 | 7.0466 +- 0.2427 | 4.4085 +- 0.0566 | 9.3331 +- 0.3242 |
| -0.45 | 1048576 | 9.54e-07 | 7.5347 +- 0.2331 | 4.9877 +- 0.0521 | 10.5914 +- 0.2954 |

## FF3 -- the ordering lemma L1

Counted: the fraction of paths with `theta^(n)_{Xi^Y_n} > Xi^y`, i.e. the walk leaving its NARROW band later than the diffusion leaves its WIDE one.  `sigma_y(y) = 1 + 0.3 tanh(y)` in `[0.7,1.3]` throughout: bounded, Lipschitz, bounded away from 0, so `asm:coeff` holds.

* **A** -- genuine Skorokhod embedding (level crossings) with `mu_y = 0`.  Because `y` is then a martingale, the embedded increments are exactly fair `+-1` coins whatever `sigma_y` is, so this IS the `eq:walks` walk and the embedding is faithful.  Node error = fine-grid crossing overshoot only.
* **B** -- SURROGATE: `theta_k := t_k` and `Y^(n)_k :=` the `sqrt(delta)`-lattice rounding of `y(t_k)`, with `mu_y(y) = 0.3 cos(y)` switched on.  This is NOT the `eq:walks` walk (its increments are multiples of `sqrt(delta)`, not `+-sqrt(delta)`); it is used only to give the node error a non-zero value of a known size, `|e_k| <= sqrt(delta)/2`.
* **C c** -- variant B plus an INJECTED node error `e_k = c * delta^{1/4} * tanh(Btilde_{t_k})`, `Btilde` an independent Brownian motion, so `max_k |e_k| <= c delta^{1/4}` exactly.  This is not an embedding; it is the largest node error the manuscript's own hypothesis (E2) permits, injected deliberately.  It answers the question L1 actually poses: is the ordering IMPLIED by (E2) + `lem:barriers-wellposed`?

The decisive comparison is `mean max |e_j|` against the barrier gap `eps_n`. Note `delta^{1/4} > delta^{1/3}` for every `delta < 1`, and the ratio `delta^{1/4}/delta^{1/3} = delta^{-1/12}` DIVERGES.

| variant | n | paths used | violations | fraction | mean max \|e_j\| | barrier gap eps_n | delta^{1/4} | no exit |
|---|---|---|---|---|---|---|---|---|
| A | 64 | 192 | 0 | 0.0000 | 0.0008 | 0.3750 | 0.3536 | 0 |
| A | 256 | 192 | 0 | 0.0000 | 0.0024 | 0.1875 | 0.2500 | 0 |
| A | 1024 | 192 | 0 | 0.0000 | 0.0033 | 0.1250 | 0.1768 | 0 |
| A | 4096 | 191 | 0 | 0.0000 | 0.0040 | 0.0781 | 0.1250 | 1 |
| B | 64 | 192 | 0 | 0.0000 | 0.0357 | 0.3750 | 0.3536 | 0 |
| B | 256 | 192 | 0 | 0.0000 | 0.0289 | 0.1875 | 0.2500 | 0 |
| B | 1024 | 192 | 0 | 0.0000 | 0.0154 | 0.1250 | 0.1768 | 0 |
| B | 4096 | 192 | 0 | 0.0000 | 0.0078 | 0.0781 | 0.1250 | 0 |
| C0.5 | 64 | 192 | 0 | 0.0000 | 0.0438 | 0.3750 | 0.3536 | 0 |
| C0.5 | 256 | 192 | 0 | 0.0000 | 0.0565 | 0.1875 | 0.2500 | 0 |
| C0.5 | 1024 | 192 | 0 | 0.0000 | 0.0431 | 0.1250 | 0.1768 | 0 |
| C0.5 | 4096 | 192 | 0 | 0.0000 | 0.0315 | 0.0781 | 0.1250 | 0 |
| C1.0 | 64 | 192 | 0 | 0.0000 | 0.0573 | 0.3750 | 0.3536 | 0 |
| C1.0 | 256 | 192 | 2 | 0.0104 | 0.0872 | 0.1875 | 0.2500 | 0 |
| C1.0 | 1024 | 192 | 6 | 0.0312 | 0.0740 | 0.1250 | 0.1768 | 0 |
| C1.0 | 4096 | 192 | 6 | 0.0312 | 0.0556 | 0.0781 | 0.1250 | 0 |
| C2.0 | 64 | 192 | 1 | 0.0052 | 0.0930 | 0.3750 | 0.3536 | 0 |
| C2.0 | 256 | 192 | 11 | 0.0573 | 0.1524 | 0.1875 | 0.2500 | 0 |
| C2.0 | 1024 | 192 | 17 | 0.0885 | 0.1311 | 0.1250 | 0.1768 | 0 |
| C2.0 | 4096 | 192 | 34 | 0.1771 | 0.1055 | 0.0781 | 0.1250 | 0 |

## FF4 -- the exit-time reconciliation L2

`G := Xi^y - theta^(n)_{Xi^Y_n} >= 0`.  With constant coefficients the exactly embedded walk sits ON the discrete barrier at `theta_{Xi^Y_n}`, so by the strong Markov property `G` is the exit time of a Brownian motion from `(B^Y,C^Y)` started at `C^Y_n`.  That law is classical, so the primary numbers below are EXACT (spectral + reflection series, `frozen_exit.psurv` / `moment_G`), validated at `kappa'=1` against `E[G] = x(L-x)` to 14 digits.  A coupled Monte-Carlo run cross-checks the reduction in (c).

### (a) E[G^kappa'] against the prediction O(eps_n^{2 kappa'}) = O(delta^{2 kappa'/3})

`eps_n` is the REALISED barrier gap `C^Y - C^Y_n`, which lies in `(delta^{1/3}, delta^{1/3}+sqrt(delta)]` and carries a lattice sawtooth; the fit against `eps_n` (target `2 kappa'`) is the cleaner one, the fit against `delta` (target `2 kappa'/3`) is the one the note quotes.

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

The prediction is asymptotic and the `kappa' -> 1/2` crossover is only logarithmically fast, so the exponents are fitted over three nested windows. `kappa' = 1` is the CONTROL: there the exponent must be exactly 1 in `eps_n` (and 1/3 in `delta`), which pins the fitting procedure.

| quantity | exp in delta (all) | exp in delta (last 7) | exp in delta (last 4) | prediction | exp in eps_n (all) | exp in eps_n (last 7) | exp in eps_n (last 4) | prediction |
|---|---|---|---|---|---|---|---|---|
| E[G^0.2] | 0.1298 | 0.1359 | 0.1376 | 0.1333 | 0.3772 | 0.4013 | 0.4089 | 0.4000 |
| E[G^0.35] | 0.2127 | 0.2247 | 0.2280 | 0.2333 | 0.6179 | 0.6638 | 0.6777 | 0.7000 |
| E[G^0.45] | 0.2574 | 0.2731 | 0.2783 | 0.3000 | 0.7477 | 0.8066 | 0.8270 | 0.9000 |
| E[G^1.0] | 0.3340 | 0.3378 | 0.3364 | 0.3333 | 0.9707 | 0.9980 | 0.9996 | 1.0000 |

The exponents rise monotonically toward the prediction as the window moves out, and the control lands on 1.000.  The residual shortfall at `kappa' = 0.45` is the `kappa' -> 1/2` crossover, where the integral `kappa' int t^{kappa'-1} P(G>t) dt` stops being dominated by its lower end.

`E[G]` is the diagnostic the note asks for.  In the idealised half-line problem `E[G]` is infinite.  In the actual BOUNDED band it is finite -- but it does not scale like `eps_n^2`, it scales like `eps_n^1`, one full power short.  The heavy tail is therefore present and biting exactly where the note says: any mean-based estimate of L2 gives `O(eps_n)`, not `O(eps_n^2)`.

### (b) raw tail P(G > t) against the prediction ~ eps_n / sqrt(t)

n = 10000, realised `eps_n = 0.05`, `t` over `[eps_n^2, T]`.  The prediction plotted is the half-line stable(1/2) tail `eps_n sqrt(2/(pi t))`.

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

### (c) coupled Monte-Carlo cross-check of the strong-Markov reduction

Same fine-grid construction as FF2 (exact level-crossing embedding), horizon 3.0T so that `Xi^y` is observed.  `G<0` counts L1 violations (must be 0); `censored` counts paths whose `Xi^y` exceeded the horizon and are dropped (this biases the MC column DOWNWARD, so agreement is the meaningful direction).

| n | paths used | G<0 | censored | walk never exits | E[G^0.2] MC / exact | E[G^0.35] MC / exact | E[G^0.45] MC / exact | E[G] MC / exact |
|---|---|---|---|---|---|---|---|---|
| 256 | 512 | 0 | 0 | 0 | 0.6113 / 0.612 | 0.4405 / 0.4433 | 0.3605 / 0.3636 | 0.1481 / 0.1523 |
| 1024 | 512 | 0 | 0 | 0 | 0.5507 / 0.546 | 0.373 / 0.3696 | 0.2952 / 0.2921 | 0.1101 / 0.1094 |
| 4096 | 512 | 0 | 0 | 0 | 0.4778 / 0.471 | 0.2991 / 0.292 | 0.2277 / 0.22 | 0.07861 / 0.07202 |

## Shortcuts, and how each could have manufactured the answer

1. **Band choice.**  `(B^Y,C^Y) = (-0.5,0.5)`, symmetric about `y_0 = 0`.  `part1-setup.tex` declares `B^Y in [0,inf)` with `y_0 = 0`, which puts `y_0` on the barrier; the constraint is a leftover from the drafts where the Y-barriers were centred at `v_0`.  A different band changes every CONSTANT below (FF1's value, FF2's plateau height) but not the exponents, because the barrier gap `eps_n ~ delta^{1/3}` is band-independent.
2. **Fine-grid continuous reference (FF2, FF3, FF4c).**  `v_t` is the left-point Riemann sum at step `dt_fine`, whose own error is `O(dt_fine^H)`. This is the shortcut most likely to fake a stall, so it is tested directly in the resolution table: refining by 16 twice moves nothing outside the error bars.  It is also the reason FF2 at `h=-0.45` is reported as INCONCLUSIVE -- there `H = 0.05`, the target rate is `delta^{0.025}`, and NO feasible `n` range separates it from a constant.
3. **Crossing overshoot.**  The exact embedding is realised on a grid, so `y(theta_k)` overshoots the lattice by `O(sqrt(dt_fine))` and the walk is taken to be the exact lattice point.  Reported per `n`; it is 2-3 orders of magnitude below `delta^{1/4}`.  It biases the node error UPWARD, i.e. against L1, and L1 still never fails in variants A and B.
4. **Fine-grid exit detection.**  `Xi^y` is the first fine-grid index outside the band, so it is LATE by `O(dt_fine)` and the effective barrier sits `~0.6 sqrt(dt_fine)` outside the true one.  This biases `G` UPWARD, i.e. in favour of a slower FF2 decay.  It is `<2%` of `eps_n` at every `n` used, and FF4's primary numbers avoid it entirely by using the exact law.
5. **FF3 variants B and C are surrogates, not embeddings**, and are labelled as such in the table.  A violation count from C is evidence about what (E2) IMPLIES, not about what the manuscript's embedding DOES.  Variant A is the only faithful construction in FF3, and for `mu_y = 0` it has `e_j == 0` in continuous time -- so a zero violation count there is close to a tautology and must not be read as support for L1.
6. **Censoring.**  Paths whose `Xi^y` exceeds the horizon are dropped in the FF4 Monte-Carlo (count reported).  Dropping the longest `G` biases the MC moments DOWNWARD; the exact column is uncensored, so MC-below-exact is the expected direction and is not evidence of disagreement.
7. **Monte-Carlo error.**  All `+-` are bootstrap standard errors over paths. In FF2 the SAME fine paths are used for every `n`, so the per-`n` estimates are strongly positively correlated and the fitted slope is far better determined than the individual points; the slope bootstrap resamples whole paths and therefore respects that.

## Compute time

| phase | wall time | seconds |
|---|---|---|
| 1. FF1 frozen-vs-unfrozen (discrete only) | 43s | 42.75 |
| 2. FF2 frozen-vs-frozen (coupled, mpu=1048576, M=1310720) | 1m30s | 90.88 |
| 3. FF2 resolution sensitivity (fixed n, varying fine grid) | 36s | 35.75 |
| 4. FF3 ordering L1 (non-constant coefficients, M=655360) | 29s | 29.21 |
| 5. FF4 exit-time reconciliation (exact law) | 3s | 2.76 |
| **total** | **3m44s** | **224.92** |

Machine: Darwin arm64, python 3.9.6.
