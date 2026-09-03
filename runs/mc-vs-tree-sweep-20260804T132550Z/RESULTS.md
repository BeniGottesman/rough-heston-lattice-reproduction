# Rough Bergomi: the tree against a trustworthy Monte-Carlo, by parameter

S0 = K = 100, xi0 = 0.0900 (sigma = 0.30), T = 1.  Baseline H = 0.1, eta = 0.3, rho = -0.7.
Grids n = [8, 16, 32, 64], refinement mref = max(4, ceil(4 sqrt(n/8))).
Barrier at 3 sd of the TRUE driver, zmax = 3/sqrt(2H); the same cap is applied to the Monte-Carlo, so the scheme error and the
absorption error are separated.  Reference: exact-covariance driver + eta=0 control variate, nfine = 512, 400,000 paths.

**Reading rule: the error must SHRINK as n grows.**  A positive slope means the scheme
converges to the wrong limit.

## Phase 1 — is the reference trustworthy?

| test | price | 95% band | driver Var[V^H_T] ratio (target 1) |
|---|---|---|---|
| nfine = 256 (X time-stepping) | 11.7833 | ±0.0099 | 0.9964 |
| nfine = 512 (X time-stepping) | 11.7855 | ±0.0099 | 0.9992 |
| nfine = 1024 (X time-stepping) | 11.7794 | ±0.0100 | 1.0086 |
| old left-endpoint driver, same paths | 11.7244 | ±0.0081 | 0.7426 |

The eta = 0 control must return Black--Scholes = **11.923538** exactly, with zero variance, because the control variate
and the payoff coincide pathwise there:

| H | price at eta=0 | 95% band | absolute error vs BS |
|---|---|---|---|
| 0.05 | 11.923538 | ±1.8e-16 | 0.0e+00 |
| 0.1 | 11.923538 | ±1.8e-16 | 0.0e+00 |
| 0.2 | 11.923538 | ±1.8e-16 | 0.0e+00 |
| 0.3 | 11.923538 | ±1.8e-16 | 0.0e+00 |
| 0.45 | 11.923538 | ±1.8e-16 | 0.0e+00 |

## Sweep: H  (roughness)

### H = 0.05  (H=0.05, eta=0.3, rho=-0.7, zmax=9.49)

MC true = **11.8204** ±0.0075 · MC capped = **11.8201** ±0.0075 · absorption cost = **-0.0002**

| n | mref | tree | scheme error (tree − MC capped) | total error (tree − MC true) |
|---|---|---|---|---|
| 8 | 4 | 11.7482 | -0.0719 | -0.0722 |
| 16 | 6 | 11.7142 | -0.1059 | -0.1062 |
| 32 | 8 | 11.7289 | -0.0912 | -0.0914 |
| 64 | 12 | 11.7924 | -0.0277 | -0.0279 |

slope of |scheme error| in n: **-0.435** → **converging**

### H = 0.1  (H=0.1, eta=0.3, rho=-0.7, zmax=6.71)

MC true = **11.7863** ±0.0070 · MC capped = **11.7863** ±0.0070 · absorption cost = **-0.0000**

| n | mref | tree | scheme error (tree − MC capped) | total error (tree − MC true) |
|---|---|---|---|---|
| 8 | 4 | 11.7714 | -0.0149 | -0.0149 |
| 16 | 6 | 11.7516 | -0.0347 | -0.0347 |
| 32 | 8 | 11.7832 | -0.0031 | -0.0031 |
| 64 | 12 | 11.8611 | +0.0748 | +0.0748 |

slope of |scheme error| in n: **0.348** → **diverging**

### H = 0.2  (H=0.2, eta=0.3, rho=-0.7, zmax=4.74)

MC true = **11.7597** ±0.0063 · MC capped = **11.7597** ±0.0063 · absorption cost = **+0.0000**

| n | mref | tree | scheme error (tree − MC capped) | total error (tree − MC true) |
|---|---|---|---|---|
| 8 | 4 | 11.8071 | +0.0474 | +0.0474 |
| 16 | 6 | 11.7900 | +0.0303 | +0.0303 |
| 32 | 8 | 11.8131 | +0.0534 | +0.0534 |
| 64 | 12 | 11.8655 | +0.1057 | +0.1058 |

slope of |scheme error| in n: **0.429** → **diverging**

### H = 0.3  (H=0.3, eta=0.3, rho=-0.7, zmax=3.87)

MC true = **11.7562** ±0.0058 · MC capped = **11.7562** ±0.0058 · absorption cost = **-0.0000**

| n | mref | tree | scheme error (tree − MC capped) | total error (tree − MC true) |
|---|---|---|---|---|
| 8 | 4 | 11.8279 | +0.0717 | +0.0717 |
| 16 | 6 | 11.8021 | +0.0459 | +0.0459 |
| 32 | 8 | 11.8053 | +0.0491 | +0.0491 |
| 64 | 12 | 11.8240 | +0.0678 | +0.0677 |

slope of |scheme error| in n: **-0.015** → **flat**

### H = 0.45  (H=0.45, eta=0.3, rho=-0.7, zmax=3.16)

MC true = **11.7655** ±0.0053 · MC capped = **11.7655** ±0.0053 · absorption cost = **-0.0000**

| n | mref | tree | scheme error (tree − MC capped) | total error (tree − MC true) |
|---|---|---|---|---|
| 8 | 4 | 11.8451 | +0.0796 | +0.0796 |
| 16 | 6 | 11.8077 | +0.0423 | +0.0423 |
| 32 | 8 | 11.7922 | +0.0267 | +0.0267 |
| 64 | 12 | 11.7840 | +0.0185 | +0.0185 |

slope of |scheme error| in n: **-0.698** → **converging**


## Sweep: eta  (vol of vol)

### eta = 0.0  (H=0.1, eta=0.0, rho=-0.7, zmax=6.71)

MC true = **11.9235** ±0.0000 · MC capped = **11.9235** ±0.0000 · absorption cost = **+0.0000**

| n | mref | tree | scheme error (tree − MC capped) | total error (tree − MC true) |
|---|---|---|---|---|
| 8 | 4 | 11.9824 | +0.0589 | +0.0589 |
| 16 | 6 | 11.9525 | +0.0290 | +0.0290 |
| 32 | 8 | 11.9399 | +0.0163 | +0.0163 |
| 64 | 12 | 11.9310 | +0.0074 | +0.0074 |

slope of |scheme error| in n: **-0.978** → **converging**

### eta = 0.1  (H=0.1, eta=0.1, rho=-0.7, zmax=6.71)

MC true = **11.8877** ±0.0024 · MC capped = **11.8877** ±0.0024 · absorption cost = **+0.0000**

| n | mref | tree | scheme error (tree − MC capped) | total error (tree − MC true) |
|---|---|---|---|---|
| 8 | 4 | 11.9406 | +0.0529 | +0.0529 |
| 16 | 6 | 11.9035 | +0.0158 | +0.0158 |
| 32 | 8 | 11.8868 | -0.0009 | -0.0009 |
| 64 | 12 | 11.8778 | -0.0099 | -0.0099 |

slope of |scheme error| in n: **-1.133** → **converging**

### eta = 0.2  (H=0.1, eta=0.2, rho=-0.7, zmax=6.71)

MC true = **11.8417** ±0.0047 · MC capped = **11.8418** ±0.0047 · absorption cost = **+0.0000**

| n | mref | tree | scheme error (tree − MC capped) | total error (tree − MC true) |
|---|---|---|---|---|
| 8 | 4 | 11.8696 | +0.0278 | +0.0278 |
| 16 | 6 | 11.8363 | -0.0055 | -0.0055 |
| 32 | 8 | 11.8345 | -0.0073 | -0.0072 |
| 64 | 12 | 11.8552 | +0.0134 | +0.0134 |

slope of |scheme error| in n: **-0.275** → **converging**

### eta = 0.3  (H=0.1, eta=0.3, rho=-0.7, zmax=6.71)

MC true = **11.7863** ±0.0070 · MC capped = **11.7863** ±0.0070 · absorption cost = **-0.0000**

| n | mref | tree | scheme error (tree − MC capped) | total error (tree − MC true) |
|---|---|---|---|---|
| 8 | 4 | 11.7714 | -0.0149 | -0.0149 |
| 16 | 6 | 11.7516 | -0.0347 | -0.0347 |
| 32 | 8 | 11.7832 | -0.0031 | -0.0031 |
| 64 | 12 | 11.8611 | +0.0748 | +0.0748 |

slope of |scheme error| in n: **0.348** → **diverging**

### eta = 0.5  (H=0.1, eta=0.5, rho=-0.7, zmax=6.71)

MC true = **11.6464** ±0.0114 · MC capped = **11.6463** ±0.0114 · absorption cost = **-0.0001**

| n | mref | tree | scheme error (tree − MC capped) | total error (tree − MC true) |
|---|---|---|---|---|
| 8 | 4 | 11.5007 | -0.1456 | -0.1457 |
| 16 | 6 | 11.5309 | -0.1154 | -0.1155 |
| 32 | 8 | 11.6804 | +0.0341 | +0.0340 |
| 64 | 12 | 11.9515 | +0.3053 | +0.3051 |

slope of |scheme error| in n: **0.145** → **flat**


## Sweep: rho  (correlation)

### rho = -0.9  (H=0.1, eta=0.3, rho=-0.9, zmax=6.71)

MC true = **11.7525** ±0.0068 · MC capped = **11.7526** ±0.0068 · absorption cost = **+0.0000**

| n | mref | tree | scheme error (tree − MC capped) | total error (tree − MC true) |
|---|---|---|---|---|
| 8 | 4 | 11.7809 | +0.0284 | +0.0284 |
| 16 | 6 | 11.7189 | -0.0337 | -0.0337 |
| 32 | 8 | 11.7179 | -0.0347 | -0.0346 |
| 64 | 12 | 11.7675 | +0.0149 | +0.0149 |

slope of |scheme error| in n: **-0.275** → **converging**

### rho = -0.7  (H=0.1, eta=0.3, rho=-0.7, zmax=6.71)

MC true = **11.7863** ±0.0070 · MC capped = **11.7863** ±0.0070 · absorption cost = **-0.0000**

| n | mref | tree | scheme error (tree − MC capped) | total error (tree − MC true) |
|---|---|---|---|---|
| 8 | 4 | 11.7714 | -0.0149 | -0.0149 |
| 16 | 6 | 11.7516 | -0.0347 | -0.0347 |
| 32 | 8 | 11.7832 | -0.0031 | -0.0031 |
| 64 | 12 | 11.8611 | +0.0748 | +0.0748 |

slope of |scheme error| in n: **0.348** → **diverging**

### rho = -0.4  (H=0.1, eta=0.3, rho=-0.4, zmax=6.71)

MC true = **11.8346** ±0.0073 · MC capped = **11.8344** ±0.0072 · absorption cost = **-0.0002**

| n | mref | tree | scheme error (tree − MC capped) | total error (tree − MC true) |
|---|---|---|---|---|
| 8 | 4 | 11.7990 | -0.0354 | -0.0356 |
| 16 | 6 | 11.8084 | -0.0260 | -0.0261 |
| 32 | 8 | 11.8713 | +0.0369 | +0.0368 |
| 64 | 12 | 11.9743 | +0.1399 | +0.1397 |

slope of |scheme error| in n: **0.646** → **diverging**

### rho = 0.0  (H=0.1, eta=0.3, rho=0.0, zmax=6.71)

MC true = **11.8917** ±0.0074 · MC capped = **11.8913** ±0.0074 · absorption cost = **-0.0004**

| n | mref | tree | scheme error (tree − MC capped) | total error (tree − MC true) |
|---|---|---|---|---|
| 8 | 4 | 11.8449 | -0.0464 | -0.0468 |
| 16 | 6 | 11.8824 | -0.0089 | -0.0093 |
| 32 | 8 | 11.9706 | +0.0793 | +0.0789 |
| 64 | 12 | 12.0968 | +0.2055 | +0.2051 |

slope of |scheme error| in n: **0.959** → **diverging**

## Synthesis

| sweep | value | MC true | ±band | absorption | error n=8 | error n=64 | n=64 in bands | slope | verdict |
|---|---|---|---|---|---|---|---|---|---|
| H | 0.05 | 11.8204 | ±0.0075 | -0.0002 | -0.0719 | -0.0277 | 3.7x | -0.435 | converging |
| H | 0.1 | 11.7863 | ±0.0070 | -0.0000 | -0.0149 | +0.0748 | 10.7x | 0.348 | diverging |
| H | 0.2 | 11.7597 | ±0.0063 | +0.0000 | +0.0474 | +0.1057 | 16.8x | 0.429 | diverging |
| H | 0.3 | 11.7562 | ±0.0058 | -0.0000 | +0.0717 | +0.0678 | 11.7x | -0.015 | flat |
| H | 0.45 | 11.7655 | ±0.0053 | -0.0000 | +0.0796 | +0.0185 | 3.5x | -0.698 | converging |
| eta | 0.0 | 11.9235 | ±0.0000 | +0.0000 | +0.0589 | +0.0074 | 7437906481.4x | -0.978 | converging |
| eta | 0.1 | 11.8877 | ±0.0024 | +0.0000 | +0.0529 | -0.0099 | 4.1x | -1.133 | converging |
| eta | 0.2 | 11.8417 | ±0.0047 | +0.0000 | +0.0278 | +0.0134 | 2.8x | -0.275 | converging |
| eta | 0.3 | 11.7863 | ±0.0070 | -0.0000 | -0.0149 | +0.0748 | 10.7x | 0.348 | diverging |
| eta | 0.5 | 11.6464 | ±0.0114 | -0.0001 | -0.1456 | +0.3053 | 26.7x | 0.145 | flat |
| rho | -0.9 | 11.7525 | ±0.0068 | +0.0000 | +0.0284 | +0.0149 | 2.2x | -0.275 | converging |
| rho | -0.7 | 11.7863 | ±0.0070 | -0.0000 | -0.0149 | +0.0748 | 10.7x | 0.348 | diverging |
| rho | -0.4 | 11.8346 | ±0.0073 | -0.0002 | -0.0354 | +0.1399 | 19.3x | 0.646 | diverging |
| rho | 0.0 | 11.8917 | ±0.0074 | -0.0004 | -0.0464 | +0.2055 | 27.8x | 0.959 | diverging |

## Compute time

| phase | wall time | seconds |
|---|---|---|
| 1. reference quality: is the Monte-Carlo trustworthy? | 39s | 38.52 |
| 2. sweep H  (roughness) | 39s | 39.21 |
| 3. sweep eta  (vol of vol) | 33s | 33.11 |
| 4. sweep rho  (correlation) | 27s | 27.18 |
| **total** | **2m18s** | **138.02** |

Machine: Darwin arm64, python 3.9.6.
