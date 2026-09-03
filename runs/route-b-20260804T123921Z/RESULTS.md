# Route B quantified — how many factors does the lift need?

Kernel Riemann-Liouville, K(u) = u^h / Gamma(1+h), h = -0.4, H = 0.1, T = 1.0, kappa = 0.5.
Lattice rate exponent (h+kappa)/2 = 0.04999999999999999; ||K||_(L2(0,T)) = 1.501531; Gamma(1+h) = 1.489192.

## Phase 1 — AJEE lift in L2(0,T), the norm (B1) is stated in

| m | partition | eta_max | shape | L2 error | quad check | relative | share on (0,d) n=64 | n=256 | node range |
|---|---|---|---|---|---|---|---|---|---|
| 1 | power | 6.444 | 1.0119 | 9.3585e-01 | 9.3585e-01 | 6.2326e-01 | 0.832 | 0.715 | 1.841 – 1.841 |
| 2 | power | 427.7 | 6.4819 | 6.9490e-01 | 6.9490e-01 | 4.6279e-01 | 0.751 | 0.719 | 1.367 – 146.2 |
| 3 | geometric | 8445.0 | 46.3274 | 5.5047e-01 | 5.5047e-01 | 3.6660e-01 | 0.731 | 0.716 | 1.124 – 3062.0 |
| 5 | geometric | 680200.0 | 21.852 | 3.8421e-01 | 3.8421e-01 | 2.5588e-01 | 0.722 | 0.651 | 0.8524 – 270600.0 |
| 8 | geometric | 67830000.0 | 11.7743 | 2.5938e-01 | 1.9145e-01 | 1.7275e-01 | 0.659 | 0.622 | 0.6177 – 29920000.0 |
| 12 | geometric | 4193000000.0 | 7.4974 | 1.7471e-01 | 1.2993e-01 | 1.1635e-01 | 0.675 | 0.622 | 0.2847 – 2036000000.0 |
| 20 | power | 41630000.0 | 7.2002 | 2.6118e-01 | 1.8278e-01 | 1.7394e-01 | 0.549 | 0.524 | 0.005102 – 34970000.0 |
| 30 | power | 23250000000.0 | 15.6515 | 1.4638e-01 | 1.0773e-01 | 9.7489e-02 | 0.551 | 0.509 | 5.05e-14 – 18210000000.0 |

decay in m (log-log slope): **-0.542**

The 'share' columns are the fraction of the squared error living on (0, delta), which a
lattice of time step delta never visits.

## Phase 2 — the achievable floor in L2(0,T)

| m | L2 error | relative | active factors | gain vs AJEE |
|---|---|---|---|---|
| 1 | 8.2838e-01 | 5.5169e-01 | 1 | x1.13 |
| 2 | 5.6120e-01 | 3.7375e-01 | 2 | x1.24 |
| 3 | 4.0372e-01 | 2.6887e-01 | 3 | x1.36 |
| 5 | 2.3057e-01 | 1.5356e-01 | 5 | x1.67 |
| 8 | 1.1457e-01 | 7.6301e-02 | 8 | x2.26 |
| 12 | 5.2366e-02 | 3.4875e-02 | 12 | x3.34 |
| 20 | 2.3416e-02 | 1.5595e-02 | 18 | x11.15 |
| 30 | 8.4386e-02 | 5.6200e-02 | 23 | x1.73 |

decay in m (log-log slope): **-0.97**

## Phase 3 — the norm the discrete scheme sees, L2(delta,T)

### n = 64, delta = 0.015625, ||K||_(L2(delta,T)) = 1.128373

| m | AJEE error | quad check | AJEE relative | floor error | floor relative |
|---|---|---|---|---|---|
| 1 | 3.6827e-01 | 3.6827e-01 | 3.2637e-01 | 2.1371e-01 | 1.8940e-01 |
| 2 | 1.6442e-01 | 1.6442e-01 | 1.4572e-01 | 3.6372e-02 | 3.2234e-02 |
| 3 | 9.3477e-02 | 9.3477e-02 | 8.2842e-02 | 6.0179e-03 | 5.3332e-03 |
| 5 | 4.1011e-02 | 4.1011e-02 | 3.6345e-02 | 1.6515e-04 | 1.4636e-04 |
| 8 | 1.8647e-02 | 1.8647e-02 | 1.6526e-02 | 9.7093e-06 | 8.6047e-06 |
| 12 | 9.2908e-03 | 9.2908e-03 | 8.2338e-03 | 6.1961e-07 | 5.4912e-07 |
| 20 | 3.8163e-03 | 3.8163e-03 | 3.3821e-03 | 1.7314e-07 | 1.5344e-07 |
| 30 | 1.8718e-03 | 1.8718e-03 | 1.6588e-03 | 0.0000e+00 | 0.0000e+00 |

decay slopes in m: AJEE **-1.584**, floor **-5.163**

### n = 256, delta = 0.003906, ||K||_(L2(delta,T)) = 1.229169

| m | AJEE error | quad check | AJEE relative | floor error | floor relative |
|---|---|---|---|---|---|
| 1 | 4.9701e-01 | 4.9701e-01 | 4.0435e-01 | 3.3947e-01 | 2.7618e-01 |
| 2 | 2.3676e-01 | 2.3676e-01 | 1.9262e-01 | 8.3485e-02 | 6.7920e-02 |
| 3 | 1.3807e-01 | 1.3807e-01 | 1.1233e-01 | 1.9902e-02 | 1.6192e-02 |
| 5 | 6.1251e-02 | 6.1251e-02 | 4.9832e-02 | 1.1142e-03 | 9.0647e-04 |
| 8 | 2.7732e-02 | 2.7732e-02 | 2.2562e-02 | 5.3153e-05 | 4.3243e-05 |
| 12 | 1.3692e-02 | 1.3692e-02 | 1.1139e-02 | 2.5406e-06 | 2.0670e-06 |
| 20 | 5.5477e-03 | 5.5477e-03 | 4.5134e-03 | 3.4466e-07 | 2.8040e-07 |
| 30 | 2.6912e-03 | 2.6912e-03 | 2.1895e-03 | 4.7122e-08 | 3.8336e-08 |

decay slopes in m: AJEE **-1.571**, floor **-5.031**

## Phase 4 — the crossing

lattice error = delta^{(h+kappa)/2} = delta^0.04999999999999999.

### n = 64  (delta = 0.015625, lattice error = **0.81225**)

| criterion | smallest m (AJEE) | smallest m (floor) |
|---|---|---|
| L2(0,T), absolute, K = u^h/Gamma(1+h) | 2 | 2 |
| L2(0,T), absolute, K = u^h | 5 | 3 |
| L2(0,T), relative to ||K|| | 1 | 1 |
| L2(delta,T), absolute, K = u^h/Gamma(1+h) | 1 | 1 |
| L2(delta,T), absolute, K = u^h | 1 | 1 |
| L2(delta,T), relative to ||K|| | 1 | 1 |

### n = 256  (delta = 0.003906, lattice error = **0.75786**)

| criterion | smallest m (AJEE) | smallest m (floor) |
|---|---|---|
| L2(0,T), absolute, K = u^h/Gamma(1+h) | 2 | 2 |
| L2(0,T), absolute, K = u^h | 5 | 3 |
| L2(0,T), relative to ||K|| | 1 | 1 |
| L2(delta,T), absolute, K = u^h/Gamma(1+h) | 1 | 1 |
| L2(delta,T), absolute, K = u^h | 1 | 1 |
| L2(delta,T), relative to ||K|| | 1 | 1 |

### For its own sake: m needed for a given relative accuracy

| target relative error | AJEE, L2(0,T) | floor, L2(0,T) | floor, L2(delta,T) at n=256 |
|---|---|---|---|
| 0.1 | 30 | 8 | 2 |
| 0.01 | None | None | 5 |
| 0.001 | None | None | 5 |

A `None` means no m in the tested range reaches that target.

## Phase 5 — the covariance surface, what Route B is for

Target for every ratio is 1.0 and for every error 0.

### Full grid on (0,T]^2, lift optimised for L2(0,T)

| m | max relative error | at (u,v) | relative Frobenius | Var[V_T] ratio |
|---|---|---|---|---|
| 1 | 6.4903e-01 | (0.0833, 0.0833) | 2.6881e-01 | 0.695640 |
| 2 | 4.3646e-01 | (0.0833, 1.0) | 1.4207e-01 | 0.860312 |
| 3 | 3.0468e-01 | (0.0833, 1.0) | 8.3137e-02 | 0.927709 |
| 5 | 1.3429e-01 | (0.0833, 1.0) | 3.4627e-02 | 0.976420 |
| 8 | 4.9461e-02 | (0.0833, 1.0) | 1.1893e-02 | 0.994178 |
| 12 | 2.4359e-02 | (0.0833, 1.0) | 4.3648e-03 | 0.998784 |
| 20 | 3.9465e-03 | (0.0833, 1.0) | 6.5500e-04 | 0.999757 |
| 30 | 1.1525e-02 | (0.0833, 0.3333) | 3.6284e-03 | 0.996842 |

### Restricted to u,v >= delta = 0.003906 (n=256), lift optimised for L2(delta,T)

| m | max relative error | at (u,v) | relative Frobenius | Var[V_T] ratio |
|---|---|---|---|---|
| 1 | 7.2719e-01 | (0.0833, 0.0833) | 3.0198e-01 | 0.627938 |
| 2 | 4.6919e-01 | (0.0833, 0.0833) | 2.3095e-01 | 0.711934 |
| 3 | 4.2226e-01 | (0.0833, 0.0833) | 2.0506e-01 | 0.743447 |
| 5 | 3.6746e-01 | (0.0833, 0.0833) | 1.7835e-01 | 0.776430 |
| 8 | 3.3555e-01 | (0.0833, 0.0833) | 1.6270e-01 | 0.795864 |
| 12 | 3.1532e-01 | (0.0833, 0.0833) | 1.5282e-01 | 0.808169 |
| 20 | 3.0580e-01 | (0.0833, 0.0833) | 1.4818e-01 | 0.813962 |
| 30 | 2.9814e-01 | (0.0833, 0.0833) | 1.4445e-01 | 0.818620 |

For comparison, the one-step lattice of Proposition 8.3 (target 1.0):

| n | 8 | 16 | 32 | 64 | 256 | 4096 |
|---|---|---|---|---|---|---|
| Var ratio | 1.06 | 1.84 | 3.2 | 5.57 | 16.89 | 155.21 |

## Phase 6 — cost in the target accuracy

n ~ eps^{-2/H} = eps^{-20} at H = 0.1.  Route A' alone: **eps^-41.0**.

| m | states | cost eps^-q | with Route A' inside |
|---|---|---|---|
| 1 | O(n^2) | 40.0 | 41.0 |
| 2 | O(n^3) | 60.0 | 61.0 |
| 3 | O(n^4) | 80.0 | 81.0 |
| 5 | O(n^6) | 120.0 | 121.0 |
| 8 | O(n^9) | 180.0 | 181.0 |
| 12 | O(n^13) | 260.0 | 261.0 |
| 20 | O(n^21) | 420.0 | 421.0 |
| 30 | O(n^31) | 620.0 | 621.0 |

## Compute time

| phase | wall time | seconds |
|---|---|---|
| 1. AJEE lift in L2(0,T) — the norm (B1) is stated in | 0s | 0.3 |
| 2. achievable floor in L2(0,T) — free nodes, NNLS weights | 14s | 13.79 |
| 3. the norm the discrete scheme sees: L2(delta,T) | 28s | 27.75 |
| 4. the crossing: kernel error vs lattice error | 0s | 0.0 |
| 5. the covariance surface — what Route B is for | 0s | 0.13 |
| **total** | **42s** | **41.99** |

Machine: Darwin arm64, python 3.9.6.
