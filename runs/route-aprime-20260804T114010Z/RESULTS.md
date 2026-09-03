# Route A' — Level 4 redone

Scheme: Route A' (embed the driver, couple the price);  H=0.1, eta=0.3, rho=-0.7.

## Phase 1 (V8) — eta=0 against the exact CRR reference

reference = **11.923389**

| n | mref | value | abs err | negative probs | mass error |
|---|---|---|---|---|---|
| 8 | 4 | 11.983327 | 0.059938 | 0 | 3.3e-07 |
| 16 | 6 | 11.952926 | 0.029536 | 0 | 3.6e-07 |
| 32 | 8 | 11.940125 | 0.016735 | 0 | 3.6e-07 |
| 64 | 12 | 11.931107 | 0.007718 | 0 | 3.6e-07 |

observed error slope in delta: **0.969**

## Phase 2 (V3') — admissibility sweep at n=32

| barrier | 4-point value | 4-point violations | A' value | A' negative probs |
|---|---|---|---|---|
| 0.84 sd (4-point limit) | 12.697 | 0 | 11.7181 | 0 |
| 2.38 sd (9-point limit) | 12.9651 | 23424 | 11.7765 | 0 |
| 3.00 sd | 12.8444 | 25530 | 11.7834 | 0 |
| none | 10.7596 | 43930 | 11.7863 | 0 |

## Phase 3 (V13) — American put under Route A'

| barrier | n=8 | n=16 | n=32 | n=64 |
|---|---|---|---|---|
| 0.84 sd (4-point limit) | 11.76809 | 11.7303 | 11.7181 | 11.71697 |
| 2.38 sd (9-point limit) | 11.77202 | 11.75123 | 11.77647 | 11.83513 |
| 3.00 sd | 11.77204 | 11.75189 | 11.78339 | 11.86118 |
| none | 11.77204 | 11.75199 | 11.78629 | 11.88542 |

spread across barrier widths: n=8: 0.00395, n=16: 0.02169, n=32: 0.06819, n=64: 0.16845

## Phase 4 (V9) — what Route A' does NOT fix

control at eta=0: MC 11.849 +- 0.0658 vs Black-Scholes 11.9235

Monte-Carlo of the TRUE rough model: **11.7309 +- 0.0482** (95%)

| n | A' tree European | gap vs MC |
|---|---|---|
| 8 | 11.7714 | +0.0405 |
| 16 | 11.7517 | +0.0208 |
| 32 | 11.7861 | +0.0553 |
| 64 | 11.8854 | +0.1545 |

The gap GROWS with n: Route A' converges, but not to the
right limit. That is the covariance discrepancy of §8, which
Route A' leaves untouched (Remark on the scope of (A'1)).
