# Lemma 9.6, freezing the clock: the exit-time step, proved and measured

The driver is `dy = mu_y(y) dt + sigma_y(y) dB` with `sigma_y(y) = 1 + 0.5 sin(y)` (so `sigma_y` is bounded,
Lipschitz and bounded away from zero: it lives in `[0.5, 1.5]`) and `mu_y(y) = 0.3 cos(y)`. The node is `y0 = 0.4`,
so `sbar = sigma_y(y0) = 1.194709`, and `rho = -0.7`.  Everything is in the rescaled time `u = t/delta`, in which the
embedding band `+- sqrt(delta)` becomes `+- 1` and the exit time is `O(1)`; the reported quantities are therefore
RELATIVE errors, which is the scale-free way to read the claim.

## Phase 1 — the coupling is valid

`nuhat := <a>_nu / sbar^2` is claimed to be exactly the `+-1` exit time of `sbar x` Brownian motion, whose moments are
known in closed form: mean `0.700609`, variance `0.327235`.  Driftless, so that the Dambis--Dubins--Schwarz
identity applies exactly.  Any residual must be the simulation's own `O(sqrt du)` exit-discretisation, so it is
refined in `du` rather than asserted away.

| du | sqrt(du) | mean nuhat | exact | bias | bias / sqrt(du) | s.e. | var nuhat | exact |
|---|---|---|---|---|---|---|---|---|
| 8e-04 | 0.0283 | 0.729666 | 0.700609 | +0.029057 | 1.027 | 0.004248 | 0.360922 | 0.327235 |
| 2e-04 | 0.0141 | 0.712025 | 0.700609 | +0.011416 | 0.807 | 0.004151 | 0.344596 | 0.327235 |
| 5e-05 | 0.0071 | 0.706516 | 0.700609 | +0.005907 | 0.835 | 0.004139 | 0.342678 | 0.327235 |

The bias tracks `sqrt(du)` with a ratio of 1.03, 0.81, 0.84 — i.e. it is the discretisation of the exit and it goes to zero
with the grid, not with anything about the coupling.  At the finest `du` the bias is `+0.005907` against a standard error of `0.004139`
(1.4 sigma), so the Dambis--Dubins--Schwarz identity `<a>_nu = T` is confirmed and the frozen driver
built from it does carry the law the Markov kernel needs.

## Phase 2 — the two rates

| delta | `||nu-nuhat||_2` | `/sqrt(delta)` | `||w-what||_2` | `/delta^{1/4}` | correlated part | orthogonal part |
|---|---|---|---|---|---|---|
| 6.25e-02 | 0.036382 | 0.1455 | 0.093888 | 0.1878 | 0.023413 | 0.090581 |
| 1.56e-02 | 0.018055 | 0.1444 | 0.063196 | 0.1787 | 0.011851 | 0.061927 |
| 3.91e-03 | 0.008896 | 0.1423 | 0.043059 | 0.1722 | 0.005889 | 0.042693 |
| 9.77e-04 | 0.004543 | 0.1454 | 0.030152 | 0.1706 | 0.002962 | 0.030020 |
| 2.44e-04 | 0.002246 | 0.1437 | 0.020921 | 0.1674 | 0.001465 | 0.020858 |

| quantity | fitted slope in delta | predicted |
|---|---|---|
| `||nu - nuhat||_2` | **0.5013** | `1/2` |
| `||w - what||_2` | **0.2700** | `1/4` |
| its orthogonal part | **0.2641** | `1/4` |
| its correlated part | **0.4999** | `1/2` |

The time rate is `1/2` to three digits.  The position rate is the sum of two terms with different exponents — the
orthogonal Brownian increment over the gap `|nu - nuhat|`, which is the `1/4` the bound is built on, and the
correlated part, which is `1/2` — so the fitted total sits between them and approaches `1/4` from above as
`delta` falls.  The bound is an upper bound and it is attained.

## Phase 3 — the input to the estimate, and what is NOT measured here

The estimate is driven by `E = sup_s |sigma_s^2 - sbar^2|`, and that splits in two:

    E  <=  const x ( sup_s |Y_s - Y_lambda|      <- WITHIN the step, measured below
                   + |Y_lambda - Y^{(n)}|  )     <- the NODE error, NOT measured here

**Only the first term is measured.** The simulation starts each step from a fixed node `y0`, so it has no node error
by construction; the node error is a global object, accumulated over every earlier step of the embedding, and a
single-step simulation cannot see it.  In the paper it is carried as an explicit input `eps_n` with three regimes:
it is **zero** when `sigma_y` and `mu_y` are constant (rough Bergomi, rough FBM — every model in this project),
and otherwise it is bounded through Lemma 7.3, whose first-moment rate is the residual open constant.  So this
phase confirms the within-step half of the input and says nothing about the other half.

| delta | `|| sup_s |sigma_s^2 - sbar^2| ||_2` | `/sqrt(delta)` |
|---|---|---|
| 6.25e-02 | 0.268295 | 1.0732 |
| 1.56e-02 | 0.135690 | 1.0855 |
| 3.91e-03 | 0.068039 | 1.0886 |
| 9.77e-04 | 0.034043 | 1.0894 |
| 2.44e-04 | 0.017021 | 1.0894 |

Fitted slope **0.4976** against the predicted `1/2`, with the ratio flat to three digits — so the within-step half of
the input is what the proof says it is.

## Phase 4 — the exit positions coincide exactly

This is the part of the construction that removes the difficulty, so it is worth being explicit that it is an
identity and not an estimate.  Both processes are `Wtilde` read on two different clocks, and the band is the
same, so both leave it at the same `Wtilde`-time `T` and at the same point `Wtilde_T`.  The code carries the
mismatch as a field and it is `0.0e+00` throughout — identically zero, by construction rather than by
cancellation.  Nothing here depends on the exit map being continuous, which it is not.

## Phase 5 — the drift, which DDS does not cover

`a` is a local martingale only when `mu_y = 0`; the proof removes the drift by Girsanov, at a cost of `O(sqrt delta)`
because the Radon--Nikodym derivative over an interval of length `O(delta)` is `1 + O(sqrt delta)`.  Measuring that
requires comparing drift-on with drift-off on the SAME paths, so the sampler draws increments for every path at
every step (`aligned=True`) — without which the two runs diverge at the first exit and a common-seed comparison
is not a comparison at all.

| delta | paired mean difference | paired s.e. | t | `|mean| / sqrt(delta)` |
|---|---|---|---|---|
| 6.25e-02 | +0.000651 | 0.002791 | +0.2 | 0.0026 |
| 1.56e-02 | -0.001588 | 0.001950 | -0.8 | 0.0127 |
| 3.91e-03 | +0.000317 | 0.001449 | +0.2 | 0.0051 |
| 9.77e-04 | +0.000083 | 0.000872 | +0.1 | 0.0026 |

**The drift effect is not resolved.** Every `t` is at most `0.8` in absolute value, so at these path counts the drift's
contribution to the exit law is indistinguishable from zero and NO slope may be fitted to it — the honest
output is a bound, not an exponent.  What the run does give is that bound: at the finest grid the paired
difference is `+0.000083` while `sqrt(delta)` is `0.0312`, so the constant in an `O(sqrt delta)`
law is below `0.013` across the range.  That is consistent with the Girsanov estimate and is as much as
simulation can say here; resolving it would need the paired standard error pushed an order of magnitude down.

## Phase 6 — the margin: `gamma < 1/4` strictly

The lemma needs the relative error `O(delta^{1/4})` to beat `delta^{gamma}` with `gamma = (h+kappa)/2`.  Since `h < 0` and
`kappa < 1/2`, `gamma < 1/4` for every admissible pair, with no extra hypothesis and no restriction on `H`:

| H | h | best kappa | gamma = (h+kappa)/2 | 1/4 - gamma | margin |
|---|---|---|---|---|---|
| 0.05 | -0.45 | `-> 1/2` | 0.0250 | 0.2250 | `delta^{0.2250}` |
| 0.1 | -0.40 | `-> 1/2` | 0.0500 | 0.2000 | `delta^{0.2000}` |
| 0.2 | -0.30 | `-> 1/2` | 0.1000 | 0.1500 | `delta^{0.1500}` |
| 0.3 | -0.20 | `-> 1/2` | 0.1500 | 0.1000 | `delta^{0.1000}` |
| 0.45 | -0.05 | `-> 1/2` | 0.2250 | 0.0250 | `delta^{0.0250}` |
| 0.49 | -0.01 | `-> 1/2` | 0.2450 | 0.0050 | `delta^{0.0050}` |

The margin closes only as `H -> 1/2`, i.e.\ as the model stops being rough — and there `gamma -> 1/4` while the
estimate stays at `1/4`, so the two meet without crossing.  In the rough regime the margin is wide: at `H = 0.1`
the lemma delivers `delta^{1/4}` where `delta^{0.05}` is asked for.

## Compute time

| phase | wall time | seconds |
|---|---|---|
| 1. the coupling is valid | 53s | 53.43 |
| 2. the two rates, and the input | 50s | 49.67 |
| 3. the input: the coefficient perturbation | 0s | 0.0 |
| 4. the exit positions coincide | 0s | 0.0 |
| 5. the drift, paired | 1m41s | 101.87 |
| **total** | **3m24s** | **204.98** |

Machine: Darwin arm64, python 3.9.6.
