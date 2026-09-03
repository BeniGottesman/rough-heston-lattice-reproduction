# `sim/cpp` — the C++ Monte-Carlo engine

`rheston_mc.cpp` prices rough Heston by Euler--Volterra Monte-Carlo. It exists
because the Volterra convolution costs `O(N^2)` per path: at `N = 800` steps and
200 000 paths that is `6.4e10` multiplications, which NumPy does not do in a
reasonable time and 9 threads of C++ do in a couple of seconds.

## Build

```bash
clang++ -O3 -march=native -std=c++17 -pthread -o sim/cpp/build/rheston_mc sim/cpp/rheston_mc.cpp
```

No dependencies — not even BLAS: the 8x8 least-squares solve of the
Longstaff--Schwartz step is a hand-written Cholesky, so the binary is
self-contained and the build cannot drift with a library version.

## Interface

Reads a CSV on stdin, writes a CSV on stdout, one line per parameter set, and
echoes a progress line per set on stderr. Columns in:

```
id,H,V0,theta,kappa,eta,rho,r,T,S0,K,steps,paths_eu,paths_reg,paths_val,ex_stride,seed
```

`paths_eu` is the number of paths for the European pass (half that many
antithetic pairs). `paths_reg` and `paths_val` are the Longstaff--Schwartz
regression and valuation samples; set both to `0` to skip the American pass.
`ex_stride` is the exercise spacing in steps, so the number of exercise dates is
`steps / ex_stride`. Every run is reproducible from `seed`: each antithetic pair
is seeded from `(seed, pair index)`, so results do not depend on how the pairs
were scheduled across threads.

Driven by `sim/run_rheston_tables.py` and
`sim/run_rheston_american_anchor.py`; do not call it by hand for anything that
is meant to end up in a document.

## What it is careful about

- **The kernel is never evaluated at 0**, where it is infinite for `H < 1/2`. The
  smallest lag is `dt`, matching the paper's exact-convolution convention.
- **Negative variances are counted, not hidden.** The `neg_hits` column is the
  number of truncations at zero over the whole run; at the Beliaeva--Nawalkha
  parameters (Feller ratio 24) it is a handful out of `10^8`.
- **The American price has no in-sample look-ahead.** The exercise policy is
  fitted on `paths_reg` and applied to the disjoint `paths_val`, so the reported
  figure is a genuine lower bound with an honest standard error. The in-sample
  figure is reported separately, and the two bracket the estimate.
- **The American call is a free control.** With `r >= 0` and no dividend its true
  value is the European call, so the gap measures what the exercise policy costs.
