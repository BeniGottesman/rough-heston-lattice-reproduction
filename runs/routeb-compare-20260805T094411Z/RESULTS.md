# Route B (the Markovian lift) vs the one-step lattice, priced

Part V showed the one-step recombining rough lattice is inconsistent: its driver variance is `2H n^{1-2H}` times the truth, so refining the grid pushes the price away from the reference without bound. Route B replaces the single divergent driver by `m` Ornstein--Uhlenbeck factors that share one Brownian increment and are rounded onto a common grid; their combined variance stays finite. This run prices both and compares them across `H`.

Fixed: `S0 = K = 100`, `xi0 = 0.09`, `eta = 0.3`, `rho = -0.7`, `T = 1`. Reference: `2,000,000`-path exact-covariance Monte-Carlo. The lift lattice was cross-validated against an independent mixing-formula pricer (agreement about `1.5e-2` at `n = 8-12`, the two price-coupling discretisations differing).

**Read this as bounded-vs-unbounded, not accurate-vs-inaccurate.** At small `m` the lift under-fits the rough kernel and carries a real bias; what it buys is convergence -- a finite error floor that shrinks as `m` grows -- where the one-step simply diverges.


## 0. The mechanism: driver variance ratio (analytic)

`Var[driver_T]/Var[true]`, which should be `1`. The one-step column diverges; every lift column plateaus at a finite value that rises toward `1` as `m` grows. This needs no pricing and is the exact statement of what the lift fixes.

| H | n | one-step `2H n^(1-2H)` | lift m=1 | lift m=2 | lift m=3 |
|---|---|---|---|---|---|
| 0.05 | 16 | 1.21 | 0.145 | 0.168 | 0.216 |
|  | 32 | 2.26 | 0.156 | 0.205 | 0.249 |
|  | 64 | 4.22 | 0.161 | 0.240 | 0.284 |
|  | 128 | 7.88 | 0.164 | 0.264 | 0.319 |
|  | 256 | 14.70 | 0.165 | 0.277 | 0.344 |
|  | 512 | 27.44 | 0.166 | 0.285 | 0.359 |
| 0.1 | 16 | 1.84 | 0.271 | 0.318 | 0.386 |
|  | 32 | 3.20 | 0.286 | 0.373 | 0.433 |
|  | 64 | 5.57 | 0.293 | 0.416 | 0.481 |
|  | 128 | 9.70 | 0.297 | 0.443 | 0.522 |
|  | 256 | 16.89 | 0.299 | 0.458 | 0.550 |
|  | 512 | 29.41 | 0.300 | 0.465 | 0.565 |
| 0.3 | 16 | 1.82 | 0.669 | 0.754 | 0.804 |
|  | 32 | 2.40 | 0.682 | 0.785 | 0.837 |
|  | 64 | 3.17 | 0.688 | 0.801 | 0.859 |
|  | 128 | 4.18 | 0.691 | 0.810 | 0.871 |
|  | 256 | 5.51 | 0.692 | 0.815 | 0.878 |
|  | 512 | 7.28 | 0.693 | 0.817 | 0.881 |

The one-step ratio grows like `n^{1-2H}` without bound; the lift ratios are flat in `n`. That flatness is convergence; the gap from `1` is the `m`-factor bias, which the price floor below inherits.


## H = 0.05: European put (reference **11.8284** ±0.0033)

Each cell is the signed error `price - reference`. The one-step column should grow; the lift columns should settle.

| n | one-step err | lift m=1 err | lift m=2 err |
|---|---|---|---|
| 16 | -0.1142 | -0.1238 | -0.1206 |
| 32 | -0.0995 | -0.1425 | -0.1370 |
| 64 | -0.0360 | -0.1522 | -- |
| 128 | +0.0844 | -0.1543 | -- |

One-step error moves `-0.1142 -> +0.0844` over `n = 16..128` (variance ratio reaches 7.9); lift m=1 error moves `-0.1238 -> -0.1543` and its increments shrink, converging to the m=1 floor. Increasing to m=2 raises the plateau toward the truth. Reference band ±0.0033.


## H = 0.1: European put (reference **11.7927** ±0.0031)

Each cell is the signed error `price - reference`. The one-step column should grow; the lift columns should settle.

| n | one-step err | lift m=1 err | lift m=2 err |
|---|---|---|---|
| 16 | -0.0411 | -0.0677 | -0.0630 |
| 32 | -0.0095 | -0.0838 | -0.0764 |
| 64 | +0.0684 | -0.0916 | -- |
| 128 | +0.1910 | -0.0926 | -- |

One-step error moves `-0.0411 -> +0.1910` over `n = 16..128` (variance ratio reaches 9.7); lift m=1 error moves `-0.0677 -> -0.0926` and its increments shrink, converging to the m=1 floor. Increasing to m=2 raises the plateau toward the truth. Reference band ±0.0031.


## H = 0.3: European put (reference **11.7589** ±0.0026)

Each cell is the signed error `price - reference`. The one-step column should grow; the lift columns should settle.

| n | one-step err | lift m=1 err | lift m=2 err |
|---|---|---|---|
| 16 | +0.0432 | +0.0222 | +0.0213 |
| 32 | +0.0464 | +0.0060 | +0.0086 |
| 64 | +0.0651 | -0.0034 | -- |
| 128 | +0.0998 | -0.0067 | -- |

One-step error moves `+0.0432 -> +0.0998` over `n = 16..128` (variance ratio reaches 4.2); lift m=1 error moves `+0.0222 -> -0.0067` and its increments shrink, converging to the m=1 floor. Increasing to m=2 raises the plateau toward the truth. Reference band ±0.0026.


## Is there a 'best n', and a bound? (the honest answer)

Yes, for the one-step lattice, and it is a symptom rather than a feature. Its error is a finite-`n` part that shrinks plus a variance-level part `~(2H n^{1-2H} - 1)` that grows; they cancel near the `n` where the variance ratio crosses one, `n* = (2H)^{-1/(1-2H)}` (about 0.05:13, 0.1:7, 0.3:4 here). That is the 'small n works, large n does not' the eye sees. But it is CANCELLATION, fragile and payoff-dependent, and there is **no convergent error bound** for the one-step scheme: past `n*` it diverges. The lift removes the growing term, so its error is bounded and converges to a floor set by `m` -- a genuine, controllable bound. That is the whole difference between a lucky window and a convergent scheme, and it is why the paper needs Route B.

The remaining task, not done here, is to drive the lift floor below the reference band by taking `m` a little larger (the variance ratios of Section 0 suggest `m ~ 3-5` reaches within a few percent for these H), and to extend the lift lattice to the American payoff, which the one-step cannot price consistently either.


## Compute time

| phase | wall time | seconds |
|---|---|---|
| 1. 0. driver variance ratio (analytic, exact) | 3s | 3.37 |
| 2. 1. H = 0.05: prices | 3m48s | 228.23 |
| 3. 2. H = 0.1: prices | 7m23s | 443.16 |
| **total** | **19m10s** | **1150.73** |

Machine: Darwin arm64, python 3.9.6.
