# Route B quantified — how many factors does the lift need?

Riemann--Liouville kernel K(u) = u^h/Gamma(1+h) with h = -0.4, H = 0.1, T = 1.0, kappa = 0.5.
Lattice rate exponent (h+kappa)/2 = 0.05; ||K||_L2(0,T) = 1.501531; Gamma(1+h) = 1.489192.

## Phase 1 — AJEE lift in L2(0,T), the norm (B1) is stated in

| m | partition | eta_max | shape | L2 error | quad check | relative | error share on (0,d), n=64 | n=256 | node range |
|---|---|---|---|---|---|---|---|---|---|
| 1 | power | 6.444 | 1.0119 | 9.3585e-01 | 9.3585e-01 | 6.2326e-01 | 0.832 | 0.715 | 1.841 – 1.841 |
| 2 | power | 427.7 | 6.4819 | 6.9490e-01 | 6.9490e-01 | 4.6279e-01 | 0.751 | 0.719 | 1.367 – 146.2 |
| 3 | geometric | 8445.0 | 46.3274 | 5.5047e-01 | 5.5047e-01 | 3.6660e-01 | 0.731 | 0.716 | 1.124 – 3062.0 |
| 5 | geometric | 680200.0 | 21.852 | 3.8421e-01 | 3.8421e-01 | 2.5588e-01 | 0.722 | 0.651 | 0.8524 – 270600.0 |
| 8 | geometric | 67830000.0 | 11.7743 | 2.5938e-01 | 2.5938e-01 | 1.7275e-01 | 0.659 | 0.622 | 0.6177 – 29920000.0 |
| 12 | geometric | 4193000000.0 | 7.4974 | 1.7471e-01 | 1.7471e-01 | 1.1635e-01 | 0.675 | 0.622 | 0.2847 – 2036000000.0 |
| 20 | geometric | 3113000000000.0 | 4.5997 | 9.8415e-02 | 9.8415e-02 | 6.5543e-02 | 0.637 | 0.575 | 0.2276 – 1717000000000.0 |
| 30 | geometric | 615500000000000.0 | 3.3012 | 6.0447e-02 | 6.0447e-02 | 4.0257e-02 | 0.621 | 0.555 | 0.1599 – 376000000000000.0 |

Decay in m (log-log slope): **-0.814**.

The two share columns give the fraction of the SQUARED error living on (0, delta), which a
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
| 20 | 1.8139e-02 | 1.2080e-02 | 19 | x5.43 |
| 30 | 8.4096e-03 | 5.6007e-03 | 28 | x7.19 |

Decay in m (log-log slope): **-1.391**.

## Phase 3 — the norm the discrete kernel sees, L2(delta,T)

### n = 64, delta = 0.015625, ||K||_L2(delta,T) = 1.128373

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

Decay slopes in m: AJEE **-1.584**, floor **-5.163**.

### n = 256, delta = 0.003906, ||K||_L2(delta,T) = 1.229169

| m | AJEE error | quad check | AJEE relative | floor error | floor relative |
|---|---|---|---|---|---|
| 1 | 4.9701e-01 | 4.9701e-01 | 4.0435e-01 | 3.3947e-01 | 2.7618e-01 |
| 2 | 2.3676e-01 | 2.3676e-01 | 1.9262e-01 | 8.3485e-02 | 6.7920e-02 |
| 3 | 1.3807e-01 | 1.3807e-01 | 1.1233e-01 | 1.9902e-02 | 1.6192e-02 |
| 5 | 6.1251e-02 | 6.1251e-02 | 4.9832e-02 | 1.1142e-03 | 9.0647e-04 |
| 8 | 2.7732e-02 | 2.7732e-02 | 2.2562e-02 | 5.3153e-05 | 4.3243e-05 |
| 12 | 1.3692e-02 | 1.3692e-02 | 1.1139e-02 | 3.2536e-06 | 2.6470e-06 |
| 20 | 5.5477e-03 | 5.5477e-03 | 4.5134e-03 | 3.4466e-07 | 2.8040e-07 |
| 30 | 2.6912e-03 | 2.6912e-03 | 2.1895e-03 | 4.7122e-08 | 3.8336e-08 |

Decay slopes in m: AJEE **-1.571**, floor **-5.013**.

## Phase 4 — the crossing

lattice error = delta^{(h+kappa)/2} = delta^0.05.

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

### For its own sake: m needed for a given relative kernel accuracy

| target relative error | AJEE, L2(0,T) | floor, L2(0,T) | floor, L2(delta,T) at n=256 |
|---|---|---|---|
| 0.1 | 20 | 8 | 2 |
| 0.01 | None | 30 | 5 |
| 0.001 | None | None | 5 |

`None` means no m in the tested range reaches that target.

## Phase 5 — the continuous covariance surface

Lift L2(0,T)-optimised.  Target: every error 0, every ratio 1.

| m | max relative error | at (u,v) | relative Frobenius | Var[V_T] ratio |
|---|---|---|---|---|
| 1 | 6.4903e-01 | (0.0833, 0.0833) | 2.6881e-01 | 0.695640 |
| 2 | 4.3646e-01 | (0.0833, 1.0) | 1.4207e-01 | 0.860312 |
| 3 | 3.0468e-01 | (0.0833, 1.0) | 8.3137e-02 | 0.927709 |
| 5 | 1.3429e-01 | (0.0833, 1.0) | 3.4627e-02 | 0.976420 |
| 8 | 4.9461e-02 | (0.0833, 1.0) | 1.1893e-02 | 0.994178 |
| 12 | 2.4359e-02 | (0.0833, 1.0) | 4.3648e-03 | 0.998784 |
| 20 | 2.5187e-03 | (0.0833, 1.0) | 4.9456e-04 | 0.999854 |
| 30 | 1.6413e-03 | (0.0833, 0.8333) | 2.6398e-04 | 0.999969 |

For comparison, the one-step lattice of Proposition 8.3 (target 1.0):

| n | 8 | 16 | 32 | 64 | 256 | 4096 |
|---|---|---|---|---|---|---|
| Var[V_T] ratio | 1.06 | 1.84 | 3.2 | 5.57 | 16.89 | 155.21 |

## Phase 6 — the DISCRETE covariance, which decides the matter

The lattice covariance is C_kl = delta sum_{j<k^l} G((k-j)d) G((l-j)d) for the discrete
kernel G.  Two conventions for G are reported: the cell average of K, which is its L2
projection on the grid, and the naive left endpoint.  The lift's own error is the
'vs true discrete' block, and it must be judged against the discretisation error the
scheme already carries with the TRUE kernel — the line under each table.

### n = 64, delta = 0.015625

#### discrete kernel: exact convolution (cell average of K)

| lift optimised for | m | vs true discrete: max rel | Frobenius | Var[V_T] ratio | total vs continuous: Frobenius | Var ratio |
|---|---|---|---|---|---|---|
| L2(0,T) | 1 | 8.0933e-01 | 1.8334e-01 | 0.863043 | 2.0989e-01 | 0.695574 |
| L2(0,T) | 2 | 4.8780e-01 | 1.1599e-01 | 0.926708 | 1.5026e-01 | 0.746885 |
| L2(0,T) | 3 | 3.1109e-01 | 7.6731e-02 | 0.996458 | 1.0972e-01 | 0.803100 |
| L2(0,T) | 5 | 1.5607e-01 | 3.5815e-02 | 0.998017 | 8.7643e-02 | 0.804357 |
| L2(0,T) | 8 | 8.4347e-02 | 1.3705e-02 | 0.996388 | 8.2311e-02 | 0.803044 |
| L2(0,T) | 12 | 3.2634e-02 | 5.1520e-03 | 1.000868 | 8.0025e-02 | 0.806655 |
| L2(0,T) | 20 | 4.2766e-03 | 5.8937e-04 | 0.999976 | 8.0149e-02 | 0.805936 |
| L2(0,T) | 30 | 3.5563e-03 | 3.2592e-04 | 1.000037 | 8.0143e-02 | 0.805985 |
| L2(delta,T) | 1 | 8.8645e-01 | 1.5064e-01 | 0.709594 | 2.0876e-01 | 0.571900 |
| L2(delta,T) | 2 | 6.5460e-01 | 1.0163e-01 | 0.803218 | 1.6447e-01 | 0.647357 |
| L2(delta,T) | 3 | 5.1292e-01 | 7.7220e-02 | 0.846157 | 1.4332e-01 | 0.681964 |
| L2(delta,T) | 5 | 3.6816e-01 | 5.3858e-02 | 0.889537 | 1.2346e-01 | 0.716927 |
| L2(delta,T) | 8 | 3.0760e-01 | 4.4536e-02 | 0.907708 | 1.1564e-01 | 0.731572 |
| L2(delta,T) | 12 | 2.6769e-01 | 3.8513e-02 | 0.919681 | 1.1063e-01 | 0.741221 |
| L2(delta,T) | 20 | 2.5743e-01 | 3.6979e-02 | 0.922760 | 1.0936e-01 | 0.743703 |
| L2(delta,T) | 30 | 2.3643e-01 | 3.3854e-02 | 0.929062 | 1.0679e-01 | 0.748782 |

Baseline — the scheme's OWN discretisation error with the TRUE kernel at this n: Frobenius **8.0154e-02**, Var[V_T] ratio **0.805955**.  A lift whose error sits below this line is not the binding constraint.

#### discrete kernel: naive left endpoint K(l*delta)

| lift optimised for | m | vs true discrete: max rel | Frobenius | Var[V_T] ratio | total vs continuous: Frobenius | Var ratio |
|---|---|---|---|---|---|---|
| L2(0,T) | 1 | 6.7922e-01 | 2.6685e-01 | 1.092848 | 1.9658e-01 | 0.672374 |
| L2(0,T) | 2 | 6.3466e-01 | 1.2122e-01 | 0.950560 | 1.9801e-01 | 0.584831 |
| L2(0,T) | 3 | 3.7555e-01 | 7.6535e-02 | 0.957110 | 2.0275e-01 | 0.588861 |
| L2(0,T) | 5 | 2.1241e-01 | 4.1001e-02 | 1.020099 | 1.8597e-01 | 0.627615 |
| L2(0,T) | 8 | 1.1965e-01 | 1.3901e-02 | 0.994488 | 1.9709e-01 | 0.611858 |
| L2(0,T) | 12 | 3.6018e-02 | 5.1974e-03 | 1.000931 | 1.9568e-01 | 0.615822 |
| L2(0,T) | 20 | 1.0992e-02 | 8.3285e-04 | 1.001021 | 1.9581e-01 | 0.615877 |
| L2(0,T) | 30 | 3.0819e-03 | 3.4499e-04 | 0.999817 | 1.9637e-01 | 0.615137 |
| L2(delta,T) | 1 | 6.9188e-01 | 9.2266e-02 | 0.908067 | 2.1800e-01 | 0.558687 |
| L2(delta,T) | 2 | 1.9120e-01 | 1.5129e-02 | 0.979429 | 2.0382e-01 | 0.602593 |
| L2(delta,T) | 3 | 3.9522e-02 | 3.1210e-03 | 0.994884 | 1.9840e-01 | 0.612102 |
| L2(delta,T) | 5 | 1.3832e-03 | 1.2190e-04 | 0.999781 | 1.9639e-01 | 0.615114 |
| L2(delta,T) | 8 | 9.5401e-05 | 6.7077e-06 | 0.999988 | 1.9630e-01 | 0.615241 |
| L2(delta,T) | 12 | 6.5274e-06 | 4.7382e-07 | 0.999999 | 1.9630e-01 | 0.615249 |
| L2(delta,T) | 20 | 3.1710e-06 | 2.4667e-07 | 1.000000 | 1.9630e-01 | 0.615249 |
| L2(delta,T) | 30 | 5.1220e-07 | 4.2542e-08 | 1.000000 | 1.9630e-01 | 0.615249 |

Baseline — the scheme's OWN discretisation error with the TRUE kernel at this n: Frobenius **1.9630e-01**, Var[V_T] ratio **0.615249**.  A lift whose error sits below this line is not the binding constraint.

### n = 256, delta = 0.003906

#### discrete kernel: exact convolution (cell average of K)

| lift optimised for | m | vs true discrete: max rel | Frobenius | Var[V_T] ratio | total vs continuous: Frobenius | Var ratio |
|---|---|---|---|---|---|---|
| L2(0,T) | 1 | 9.3549e-01 | 1.8137e-01 | 0.815573 | 1.8804e-01 | 0.695636 |
| L2(0,T) | 2 | 4.8642e-01 | 1.1610e-01 | 0.980464 | 1.2036e-01 | 0.836278 |
| L2(0,T) | 3 | 3.1046e-01 | 7.6609e-02 | 0.987962 | 8.2986e-02 | 0.842673 |
| L2(0,T) | 5 | 1.9216e-01 | 3.5947e-02 | 0.987230 | 4.8876e-02 | 0.842049 |
| L2(0,T) | 8 | 8.5631e-02 | 1.3731e-02 | 1.000979 | 3.3995e-02 | 0.853776 |
| L2(0,T) | 12 | 3.7477e-02 | 5.1557e-03 | 0.999297 | 3.1860e-02 | 0.852341 |
| L2(0,T) | 20 | 6.2158e-03 | 6.0168e-04 | 1.000172 | 3.1298e-02 | 0.853087 |
| L2(0,T) | 30 | 3.7221e-03 | 3.2478e-04 | 1.000010 | 3.1329e-02 | 0.852950 |
| L2(delta,T) | 1 | 9.5128e-01 | 1.2616e-01 | 0.736200 | 1.4041e-01 | 0.627935 |
| L2(delta,T) | 2 | 7.5515e-01 | 5.7661e-02 | 0.834520 | 7.8836e-02 | 0.711796 |
| L2(delta,T) | 3 | 6.0071e-01 | 4.1302e-02 | 0.870892 | 6.4732e-02 | 0.742819 |
| L2(delta,T) | 5 | 4.3162e-01 | 2.8281e-02 | 0.907266 | 5.3509e-02 | 0.773844 |
| L2(delta,T) | 8 | 3.4063e-01 | 2.1862e-02 | 0.926811 | 4.8105e-02 | 0.790515 |
| L2(delta,T) | 12 | 2.8989e-01 | 1.8413e-02 | 0.937713 | 4.5258e-02 | 0.799814 |
| L2(delta,T) | 20 | 2.6485e-01 | 1.6739e-02 | 0.943094 | 4.3895e-02 | 0.804404 |
| L2(delta,T) | 30 | 2.4693e-01 | 1.5554e-02 | 0.946944 | 4.2937e-02 | 0.807688 |

Baseline — the scheme's OWN discretisation error with the TRUE kernel at this n: Frobenius **3.1329e-02**, Var[V_T] ratio **0.852941**.  A lift whose error sits below this line is not the binding constraint.

#### discrete kernel: naive left endpoint K(l*delta)

| lift optimised for | m | vs true discrete: max rel | Frobenius | Var[V_T] ratio | total vs continuous: Frobenius | Var ratio |
|---|---|---|---|---|---|---|
| L2(0,T) | 1 | 8.2232e-01 | 2.1852e-01 | 0.974777 | 1.8256e-01 | 0.689774 |
| L2(0,T) | 2 | 5.3534e-01 | 1.1547e-01 | 0.940691 | 1.2526e-01 | 0.665654 |
| L2(0,T) | 3 | 4.3136e-01 | 8.1757e-02 | 1.025918 | 9.5953e-02 | 0.725962 |
| L2(0,T) | 5 | 2.6010e-01 | 3.6442e-02 | 0.992892 | 8.8634e-02 | 0.702592 |
| L2(0,T) | 8 | 1.0062e-01 | 1.4028e-02 | 1.003842 | 8.5385e-02 | 0.710341 |
| L2(0,T) | 12 | 5.9130e-02 | 5.1241e-03 | 0.997427 | 8.7345e-02 | 0.705801 |
| L2(0,T) | 20 | 8.9679e-03 | 6.1808e-04 | 1.000420 | 8.6564e-02 | 0.707919 |
| L2(0,T) | 30 | 3.5471e-03 | 3.3088e-04 | 0.999822 | 8.6831e-02 | 0.707496 |
| L2(delta,T) | 1 | 8.6560e-01 | 1.2920e-01 | 0.881233 | 1.3977e-01 | 0.623580 |
| L2(delta,T) | 2 | 3.7818e-01 | 2.9406e-02 | 0.972889 | 9.4255e-02 | 0.688438 |
| L2(delta,T) | 3 | 1.1741e-01 | 6.0631e-03 | 0.991307 | 8.9262e-02 | 0.701471 |
| L2(delta,T) | 5 | 8.5539e-03 | 3.7239e-04 | 0.999157 | 8.7039e-02 | 0.707026 |
| L2(delta,T) | 8 | 4.7843e-04 | 1.8943e-05 | 0.999953 | 8.6793e-02 | 0.707589 |
| L2(delta,T) | 12 | 3.3016e-05 | 1.0247e-06 | 0.999997 | 8.6779e-02 | 0.707620 |
| L2(delta,T) | 20 | 5.7276e-06 | 1.8147e-07 | 0.999999 | 8.6779e-02 | 0.707622 |
| L2(delta,T) | 30 | 1.4089e-06 | 4.9455e-08 | 1.000000 | 8.6778e-02 | 0.707622 |

Baseline — the scheme's OWN discretisation error with the TRUE kernel at this n: Frobenius **8.6778e-02**, Var[V_T] ratio **0.707622**.  A lift whose error sits below this line is not the binding constraint.

## Phase 7 — cost in the target accuracy

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
| 1. AJEE lift in L2(0,T) — the norm (B1) is stated in | 0s | 0.37 |
| 2. achievable floor in L2(0,T) — free nodes, NNLS weights | 17s | 16.8 |
| 3. the norm the discrete kernel sees: L2(delta,T) | 25s | 25.47 |
| 4. the crossing: kernel error vs lattice error | 0s | 0.0 |
| 5. the continuous covariance surface | 0s | 0.08 |
| 6. the DISCRETE covariance — adjudicating the two norms | 2s | 1.91 |
| **total** | **45s** | **44.64** |

Machine: Darwin arm64, python 3.9.6.
