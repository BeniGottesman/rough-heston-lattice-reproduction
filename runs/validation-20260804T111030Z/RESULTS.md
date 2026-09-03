# Validation run — results

Example: fractional driver + constant-vol American put

## Phase 1 (V1, V2) — the rough component

| H | slope of Var[Vcheck]/Var[V] vs n | predicted 1-2H | Var[V^(n)]/Var[V] at n=4096 | Var[Vcheck]/Var[V] at n=4096 |
|---|---|---|---|---|
| 0.10 | 0.8 | 0.8 | 0.8319 | 155.21 |
| 0.20 | 0.6 | 0.6 | 0.972 | 58.81 |
| 0.30 | 0.4 | 0.4 | 0.9954 | 16.71 |
| 0.40 | 0.2 | 0.2 | 0.9993 | 4.22 |

## Phase 2 (V3, V8) — lattice vs CRR reference

CRR reference (n=20000): **11.923389**

| n | lattice | abs err | prob. violations |
|---|---|---|---|
| 8 | 11.559493 | 0.363897 | 0 |
| 16 | 11.739737 | 0.183652 | 0 |
| 32 | 11.831232 | 0.092157 | 0 |
| 64 | 11.877289 | 0.046101 | 0 |
| 128 | 11.90039 | 0.022999 | 0 |
| 256 | 11.911959 | 0.011431 | 0 |

observed error slope in delta: 0.999

## Phase 3a (V13) — rough Bergomi, NO barrier

| n | value | prob. violations |
|---|---|---|
| 8 | 11.661999 | 546 |
| 16 | 11.626183 | 5266 |
| 32 | 10.7596 | 43930 |
| 64 | 6.530714 | 351894 |
| 128 | 0.713493 | 2783030 |

observed r = -2.432 — meaningless: the value collapses, the scheme diverges.

## Phase 3b (V13) — rough Bergomi, WITH the admissible barrier

zeta clamped to |zeta| <= 1.877 (= 0.84 standard deviations of zeta_T; the 9-point kernel of Thm 9.1 would allow 5.317)

| n | value | prob. violations |
|---|---|---|
| 8 | 12.001208 | 0 |
| 16 | 12.416449 | 0 |
| 32 | 12.696974 | 0 |
| 64 | 12.891303 | 0 |
| 128 | 13.026175 | 0 |

observed r = 0.54 (theory H/2 = 0.05)

