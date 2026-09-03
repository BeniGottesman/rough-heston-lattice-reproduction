# What a lattice price costs, and what the same accuracy costs by Monte-Carlo

`ALGO_VERSION` 1.3.2 · `code_sha` 457dcce1119a

A lattice price is one deterministic number carrying a **bias**. A
Monte-Carlo price is unbiased but carries a **standard error** that
shrinks like `1/sqrt(N)`. So the only meaningful comparison matches
the accuracies first:

    N*  =  (sigma_eff / |lattice bias|)^2,

and then puts the two wall times side by side.

Parameters: rough Bergomi, `T=1`, `S0=K=100`, `xi0=0.09`, `eta=0.3`, `rho=-0.7`,
`mref = max(4, ceil(4 sqrt(n/8)))`, barrier `zmax = 3/sqrt(2H)`. The
reference is UNCAPPED, so the lattice's barrier truncation counts as
part of its error.

## 0. The answer at n = 16

| H | pricer | price | bias | lattice time | N* paths | MC time (paths only) | MC / lattice |
|---|---|---|---|---|---|---|---|
| 0.05 | one-step | 11.7142 | -0.1085 | 28ms | 587 | 4ms | 0.14x |
| 0.05 | lift m=1 | 11.7049 | -0.1178 | 17ms | 498 | 3ms | 0.19x |
| 0.05 | lift m=2 | 11.7079 | -0.1148 | 847ms | 524 | 4ms | 0.00x |
| 0.1 | one-step | 11.7516 | -0.0370 | 28ms | 4.3k | 26ms | 0.95x |
| 0.1 | lift m=1 | 11.7253 | -0.0633 | 26ms | 1.5k | 9ms | 0.35x |
| 0.1 | lift m=2 | 11.7298 | -0.0587 | 1.58s | 1.7k | 10ms | 0.01x |
| 0.3 | one-step | 11.8021 | +0.0443 | 28ms | 2.0k | 11ms | 0.40x |
| 0.3 | lift m=1 | 11.7816 | +0.0237 | 37ms | 7.0k | 40ms | 1.07x |
| 0.3 | lift m=2 | 11.7802 | +0.0223 | 3.55s | 7.9k | 45ms | 0.01x |

`MC / lattice` above 1 means the Monte-Carlo is the slower way to
reach that accuracy; below 1 means it is the faster way.

## 1. Monte-Carlo cost, measured

| H | setup (cold) | marginal rate | sigma_eff (empirical) | sigma_eff (code's formula) | ratio |
|---|---|---|---|---|---|
| 0.05 | 3.58s | 147,781 paths/s | 2.6283 | 2.4063 | 1.0923 |
| 0.1 | 3.88s | 163,469 paths/s | 2.4329 | 2.2602 | 1.0764 |
| 0.3 | 3.81s | 175,768 paths/s | 1.9857 | 1.8832 | 1.0545 |

`sigma_eff` is what the Monte-Carlo's error actually scales with, and
it is measured EXACTLY rather than estimated. With `N` paths laid out
as `N/2` antithetic pairs,

    Var(mean) = Var(sample) (1 + rho_pair) / N,

so `sigma_eff = sd(sample) * sqrt(1 + rho_pair)`, and both factors come
from one large sample with negligible uncertainty. The code's own
`stderr` is the `rho_pair = 0` case, so the last column is exactly the
factor by which that formula is wrong.

### The antithetic device and the control variate work against each other

| H | rho_pair (raw put) | rho_pair (control-adjusted) | sd (raw) | sd (control) | sigma_eff inflation |
|---|---|---|---|---|---|
| 0.05 | -0.5912 | +0.1931 | 15.255 | 2.406 | 1.0923 |
| 0.1 | -0.5829 | +0.1586 | 15.368 | 2.260 | 1.0764 |
| 0.3 | -0.5784 | +0.1119 | 15.452 | 1.883 | 1.0545 |

This is the finding of the phase, and it was not expected. On the RAW
put the antithetic device works well — a strongly negative pair
correlation, which is what antithetic sampling is for, since a put
payoff is monotone in the driver. But the `eta = 0` control variate
subtracts a second, similar put, and a DIFFERENCE of two puts is not
monotone in the driver. On that adjusted sample the pair correlation
turns POSITIVE, so the antithetic pairing is mildly counterproductive
and the quoted standard error is optimistic by the inflation factor.

The control variate is by far the bigger effect — it cuts the sample
spread several-fold — so the net estimator is still much better than a
plain one. But two consequences follow, and both matter for the rest
of the project:

1. Every `+-s.e.` quoted from this Monte-Carlo, and every 'within N MC
   bands' claim resting on one, is too generous by the inflation
   factor above. The `N*` columns here use the corrected `sigma_eff`.
2. Dropping the antithetic pairing (while keeping the control variate)
   would make the estimator slightly BETTER, not worse. That is a
   one-line change and it has not been made or tested here.

### Cross-validation of the rebuilt sample

`sigma_eff_direct` rebuilds the estimator's per-path sample, because
the pricer returns only aggregates. That duplication is a risk, so it
is checked against the real pricer at the same seed and path count.

What the check can be is limited, and worth being precise about: the
pricer draws in CHUNKS from one rng stream while the rebuild draws the
whole block at once, so the two consume the stream differently and are
**different draws from the same distribution**. The price therefore
agrees only to within a couple of standard errors — that column is
reported in s.e. units so it judges itself. What must agree tightly is
the standard error, because that estimates a property of the
distribution rather than of the particular draw.

| H | price (pricer) | price (rebuilt) | gap in s.e. | s.e. (pricer) | s.e. (rebuilt) | s.e. ratio |
|---|---|---|---|---|---|---|
| 0.05 | 11.833932 | 11.828716 | -0.97 | 0.005390 | 0.005381 | 0.9983 |
| 0.1 | 11.799883 | 11.793670 | -1.23 | 0.005065 | 0.005054 | 0.9979 |
| 0.3 | 11.767504 | 11.761072 | -1.53 | 0.004217 | 0.004211 | 0.9986 |

### The independent replication estimator, for comparison

The same quantity from the scatter of 100 independent runs of
20,000 paths. This estimator needs no theory at all, which
makes it a genuine check — but a standard deviation from
100 replications carries about
7% relative uncertainty of
its own, and that uncertainty enters `N*` SQUARED, which is why it is
not the estimator used:

| H | sigma_eff (exact) | sigma_eff (replications) | ratio |
|---|---|---|---|
| 0.05 | 2.6283 | 2.4899 | 0.947 |
| 0.1 | 2.4329 | 2.2551 | 0.927 |
| 0.3 | 1.9857 | 1.9390 | 0.977 |

Agreement to within that uncertainty is the check passing; a large
systematic gap in one direction across every `H` would mean the exact
formula is missing a correlation the replications can see. This
estimator was run at low power first — 24 replications, ±15% — and
came out about 1.3x high on every H, which looked systematic; raising
it to 100 replications (±7%) brought it back to about
1.08. The lesson is kept rather than quietly deleted: a std from two
dozen samples is not evidence of a systematic effect.

The **setup** is a one-off `nfine x nfine` covariance factorisation
(`nfine = 512`). It is charged once per process, not per path,
which is why a small Monte-Carlo can measure slower than a large one
on a cold cache. Both readings are given below.

## 2. The reference prices (the 'truth')

| H | price | s.e. | 95% CI half-width | paths | wall time |
|---|---|---|---|---|---|
| 0.05 | 11.8227 | 0.0017 | 0.0033 | 2,000,000 | 12.73s |
| 0.1 | 11.7886 | 0.0016 | 0.0031 | 2,000,000 | 11.91s |
| 0.3 | 11.7579 | 0.0013 | 0.0026 | 2,000,000 | 11.50s |

Any bias smaller than about twice these standard errors is not
resolved, and the corresponding `N*` is a band rather than a number.

## 3. Every pricer, every n

| H | pricer | n | price | bias | bias / ref s.e. | resolved | state space | best-of-3 | median |
|---|---|---|---|---|---|---|---|---|---|
| 0.05 | one-step | 16 | 11.7142 | -0.1085 | 63.6 | yes | — | 28ms | 28ms |
| 0.05 | one-step | 32 | 11.7289 | -0.0938 | 55.0 | yes | — | 222ms | 228ms |
| 0.05 | one-step | 64 | 11.7924 | -0.0303 | 17.7 | yes | — | 2.96s | 3.57s |
| 0.05 | lift m=1 | 16 | 11.7049 | -0.1178 | 69.0 | yes | 51 | 17ms | 18ms |
| 0.05 | lift m=1 | 32 | 11.6863 | -0.1364 | 80.0 | yes | 77 | 100ms | 102ms |
| 0.05 | lift m=1 | 64 | 11.6766 | -0.1461 | 85.6 | yes | 113 | 760ms | 786ms |
| 0.05 | lift m=2 | 16 | 11.7079 | -0.1148 | 67.3 | yes | 2635 | 847ms | 853ms |
| 0.05 | lift m=2 | 32 | 11.6915 | -0.1312 | 76.9 | yes | 6477 | 10.20s | 10.70s |
| 0.1 | one-step | 16 | 11.7516 | -0.0370 | 23.1 | yes | — | 28ms | 28ms |
| 0.1 | one-step | 32 | 11.7832 | -0.0054 | 3.4 | yes | — | 222ms | 223ms |
| 0.1 | one-step | 64 | 11.8611 | +0.0725 | 45.3 | yes | — | 3.05s | 3.10s |
| 0.1 | lift m=1 | 16 | 11.7253 | -0.0633 | 39.6 | yes | 67 | 26ms | 26ms |
| 0.1 | lift m=1 | 32 | 11.7092 | -0.0794 | 49.6 | yes | 99 | 141ms | 143ms |
| 0.1 | lift m=1 | 64 | 11.7015 | -0.0871 | 54.4 | yes | 145 | 1.12s | 1.12s |
| 0.1 | lift m=2 | 16 | 11.7298 | -0.0587 | 36.7 | yes | 4025 | 1.58s | 1.58s |
| 0.1 | lift m=2 | 32 | 11.7165 | -0.0721 | 45.0 | yes | 9633 | 14.85s | 15.44s |
| 0.3 | one-step | 16 | 11.8021 | +0.0443 | 33.3 | yes | — | 28ms | 30ms |
| 0.3 | one-step | 32 | 11.8053 | +0.0474 | 35.7 | yes | — | 220ms | 225ms |
| 0.3 | one-step | 64 | 11.8240 | +0.0661 | 49.7 | yes | — | 3.03s | 3.11s |
| 0.3 | lift m=1 | 16 | 11.7816 | +0.0237 | 17.9 | yes | 91 | 37ms | 37ms |
| 0.3 | lift m=1 | 32 | 11.7653 | +0.0075 | 5.6 | yes | 131 | 206ms | 207ms |
| 0.3 | lift m=1 | 64 | 11.7559 | -0.0019 | 1.5 | NO | 187 | 1.67s | 1.75s |
| 0.3 | lift m=2 | 16 | 11.7802 | +0.0223 | 16.8 | yes | 7865 | 3.55s | 3.63s |
| 0.3 | lift m=2 | 32 | 11.7678 | +0.0099 | 7.5 | yes | 13975 | 23.46s | 23.70s |

## 4. The equivalence, with the reference's uncertainty propagated

`N*` band = what `N*` becomes if the true bias is at either end of
`|bias| ± 2 · ref s.e.`. A smaller true bias needs MORE paths, so the
band is wide when the bias is barely resolved.

| H | pricer | n | N* | N* band | MC time (paths) | MC time (+setup) | MC/lattice (paths) | MC/lattice (+setup) | break-even error |
|---|---|---|---|---|---|---|---|---|---|
| 0.05 | one-step | 16 | 587 | 552–625 | 4ms | 3.59s | 0.14x | 129.21x | 0.0410 |
| 0.05 | one-step | 32 | 786 | 732–846 | 5ms | 3.59s | 0.02x | 16.19x | 0.0145 |
| 0.05 | one-step | 64 | 7.5k | 6.1k–9.6k | 51ms | 3.64s | 0.02x | 1.23x | 0.0040 |
| 0.05 | lift m=1 | 16 | 498 | 470–528 | 3ms | 3.59s | 0.19x | 205.34x | 0.0517 |
| 0.05 | lift m=1 | 32 | 371 | 353–390 | 3ms | 3.59s | 0.03x | 35.72x | 0.0216 |
| 0.05 | lift m=1 | 64 | 324 | 309–339 | 2ms | 3.59s | 0.00x | 4.72x | 0.0078 |
| 0.05 | lift m=2 | 16 | 524 | 495–557 | 4ms | 3.59s | 0.00x | 4.24x | 0.0074 |
| 0.05 | lift m=2 | 32 | 401 | 381–423 | 3ms | 3.59s | 0.00x | 0.35x | 0.0021 |
| 0.1 | one-step | 16 | 4.3k | 3.7k–5.2k | 26ms | 3.91s | 0.95x | 141.20x | 0.0362 |
| 0.1 | one-step | 32 | 204.6k | 80.4k–1.2M | 1.25s | 5.13s | 5.64x | 23.14x | 0.0128 |
| 0.1 | one-step | 64 | 1.1k | 1.0k–1.2k | 7ms | 3.89s | 0.00x | 1.28x | 0.0034 |
| 0.1 | lift m=1 | 16 | 1.5k | 1.3k–1.6k | 9ms | 3.89s | 0.35x | 149.88x | 0.0373 |
| 0.1 | lift m=1 | 32 | 940 | 868–1.0k | 6ms | 3.89s | 0.04x | 27.51x | 0.0160 |
| 0.1 | lift m=1 | 64 | 780 | 726–841 | 5ms | 3.89s | 0.00x | 3.48x | 0.0057 |
| 0.1 | lift m=2 | 16 | 1.7k | 1.5k–1.9k | 10ms | 3.89s | 0.01x | 2.47x | 0.0048 |
| 0.1 | lift m=2 | 32 | 1.1k | 1.0k–1.2k | 7ms | 3.89s | 0.00x | 0.26x | 0.0016 |
| 0.3 | one-step | 16 | 2.0k | 1.8k–2.3k | 11ms | 3.82s | 0.40x | 134.87x | 0.0281 |
| 0.3 | one-step | 32 | 1.8k | 1.6k–2.0k | 10ms | 3.82s | 0.05x | 17.34x | 0.0101 |
| 0.3 | one-step | 64 | 902 | 834–979 | 5ms | 3.82s | 0.00x | 1.26x | 0.0027 |
| 0.3 | lift m=1 | 16 | 7.0k | 5.7k–8.9k | 40ms | 3.85s | 1.07x | 103.56x | 0.0246 |
| 0.3 | lift m=1 | 32 | 70.7k | 38.5k–170.7k | 402ms | 4.21s | 1.96x | 20.48x | 0.0104 |
| 0.3 | lift m=1 | 64 | 1.0M | 186.3k–— | 5.96s | 9.77s | 3.57x | 5.86x | 0.0037 |
| 0.3 | lift m=2 | 16 | 7.9k | 6.3k–10.2k | 45ms | 3.86s | 0.01x | 1.09x | 0.0025 |
| 0.3 | lift m=2 | 32 | 40.0k | 24.9k–74.7k | 228ms | 4.04s | 0.01x | 0.17x | 0.0010 |

**Break-even error** is the accuracy at which the Monte-Carlo, given
the lattice's wall time in paths, lands. If it is SMALLER than the
lattice's bias, the Monte-Carlo wins at that budget.

## 5. What Monte-Carlo cannot match cheaply: the American price

The lattice gets the early-exercise price from the same backward
pass. Monte-Carlo needs Longstaff--Schwartz, which is more expensive
and itself biased — this project has an LSM estimator for rough
Heston, not for rough Bergomi, so no equal-accuracy comparison is
offered here and none is implied.

| H | pricer | n | European | American | premium | European time | American time |
|---|---|---|---|---|---|---|---|
| 0.05 | lift m=1 | 16 | 11.7049 | 11.7054 | +0.0004 | 17ms | 19ms |
| 0.05 | lift m=1 | 32 | 11.6863 | 11.6866 | +0.0003 | 100ms | 109ms |
| 0.1 | lift m=1 | 16 | 11.7253 | 11.7257 | +0.0005 | 26ms | 25ms |
| 0.1 | lift m=1 | 32 | 11.7092 | 11.7096 | +0.0003 | 141ms | 140ms |
| 0.3 | lift m=1 | 16 | 11.7816 | 11.7822 | +0.0006 | 37ms | 35ms |
| 0.3 | lift m=1 | 32 | 11.7653 | 11.7658 | +0.0004 | 206ms | 200ms |

## 6. How to read this — and how not to

- **This compares implementations, not methods.** The Monte-Carlo's
  inner loop is a BLAS matrix product and runs multi-threaded; the
  lattices are single-threaded Python loops over a state space, with
  no C++ port for Route B. A constant factor of tens is
  implementation, not mathematics.
- **The scalings, so a reader can rescale:** Monte-Carlo is
  `O(N x nfine)` per price; the one-step lattice is `O(n x nx)`; the
  lift at `m` factors is `O(n x nx x prod_i N_i)`, and that state
  space is what makes `m = 2` expensive.
- **A small bias is not a good thing here.** Part V established that
  the one-step lattice's error at small `n` is a CANCELLATION between
  a shrinking discretisation term and a growing variance term, with
  the crossing near `n* = (2H)^{-1/(1-2H)}` — about 13 at H = 0.05,
  7 at H = 0.1, 4 at H = 0.3. So a flattering bias at `n = 16` is
  luck, not accuracy: it is payoff-dependent and carries no error
  bound. Refining `n` makes the one-step scheme WORSE, which is
  visible in the n = 32 and n = 64 rows of section 3.
- **The lattice's case is not speed.** It is that Route B is a
  convergent scheme with a proof and an American price in the same
  pass. Nothing in this table argues it is the fast way to get a
  European number.

## Timing

| phase | wall time |
|---|---|
| setup | 15.00s |
| sigma | 53.06s |
| reference | 36.14s |
| lattice | 3m33s |
| equiv | 590µs |
| american | 534ms |
| **total** | **5m18s** |

Machine: single process. Reference 2,000,000 paths at nfine=512 per H; sigma from 100x20,000; lattice timings best-of-3.