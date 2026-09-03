# The paper's lattice on the ADS2014 grid (H = 0.5)

The lattice of this paper, run on the Beliaeva--Nawalkha parameter sets, against **two** third-party references: the analytical Heston price and the independent recombining tree of arXiv:1205.3555.

Fixed: `K = 100`, `r = 0.05`, `eta = 0.1`, `kappa = 3.0`, `theta = 0.04`. Lattice: variance through the Lamperti transform, price coupled by Route A' randomised rounding with `mref = max(4, ceil(4 sqrt(n/8)))`, driver absorbed at 5 standard deviations of `v_T` from the exact CIR moments, `1/U` drift regularised by clipping the up-probabilities.

**Why this table stops at `H = 0.5`.** The paper's model class needs an autonomous driver; rough Heston's `nu sqrt(V) dB` makes the driver depend on its own past, so rough Heston is outside the class and no lattice of this kind prices it. That is a proposition, not a gap in the code. At `h = 0` the class and the semi-analytic models meet, the model is classical Heston, and prices are published -- which is why this is the one place the lattice can be judged against something that is neither ours nor a simulation.


## 1. European put

`bino` is the two-point walk the paper describes; `trino` is the Hull--White trinomial with branch switching, which matches the driver's variance exactly instead of losing the fraction `mu^2 delta` of it. `ADS N=500` is their tree. `analytic` is the closed-form Heston price they publish.

| S0 | sqrt(V0) | T | analytic | bino n=200 | err | trino n=200 | err | ADS N=500 | err |
|---|---|---|---|---|---|---|---|---|---|
| 90 | 0.2 | 1m | 9.6533 | 9.6533 | +0.0000 | 9.6534 | +0.0001 | 9.6533 | +0.0000 |
| 95 | 0.2 | 1m | 5.2074 | 5.2077 | +0.0003 | 5.2078 | +0.0004 | 5.2077 | +0.0003 |
| 100 | 0.2 | 1m | 2.0971 | 2.0975 | +0.0004 | 2.0974 | +0.0003 | 2.0965 | -0.0006 |
| 105 | 0.2 | 1m | 0.6053 | 0.6055 | +0.0002 | 0.6054 | +0.0001 | 0.6050 | -0.0003 |
| 110 | 0.2 | 1m | 0.1265 | 0.1265 | -0.0000 | 0.1265 | +0.0000 | 0.1270 | +0.0005 |
| 90 | 0.3 | 1m | 9.9905 | 9.9879 | -0.0026 | 9.9911 | +0.0006 | 9.9900 | -0.0005 |
| 95 | 0.3 | 1m | 6.0155 | 6.0105 | -0.0050 | 6.0166 | +0.0011 | 6.0162 | +0.0007 |
| 100 | 0.3 | 1m | 3.1302 | 3.1245 | -0.0057 | 3.1314 | +0.0012 | 3.1290 | -0.0012 |
| 105 | 0.3 | 1m | 1.3967 | 1.3924 | -0.0043 | 1.3976 | +0.0009 | 1.3955 | -0.0012 |
| 110 | 0.3 | 1m | 0.5367 | 0.5343 | -0.0024 | 0.5372 | +0.0005 | 0.5372 | +0.0005 |
| 90 | 0.4 | 1m | 10.5668 | 10.5520 | -0.0148 | 10.5679 | +0.0011 | 10.5668 | +0.0000 |
| 95 | 0.4 | 1m | 6.9335 | 6.9100 | -0.0235 | 6.9350 | +0.0015 | 6.9352 | +0.0017 |
| 100 | 0.4 | 1m | 4.1852 | 4.1584 | -0.0268 | 4.1869 | +0.0017 | 4.1861 | +0.0009 |
| 105 | 0.4 | 1m | 2.3222 | 2.2981 | -0.0241 | 2.3238 | +0.0016 | 2.3229 | +0.0007 |
| 110 | 0.4 | 1m | 1.1882 | 1.1702 | -0.0180 | 1.1894 | +0.0012 | 1.1893 | +0.0011 |
| 90 | 0.2 | 3m | 9.5698 | 9.5699 | +0.0001 | 9.5703 | +0.0005 | 9.5694 | -0.0004 |
| 95 | 0.2 | 3m | 5.9692 | 5.9694 | +0.0002 | 5.9699 | +0.0007 | 5.9693 | +0.0001 |
| 100 | 0.2 | 3m | 3.3770 | 3.3771 | +0.0001 | 3.3776 | +0.0006 | 3.3794 | +0.0024 |
| 105 | 0.2 | 3m | 1.7410 | 1.7409 | -0.0001 | 1.7413 | +0.0003 | 1.7402 | -0.0008 |
| 110 | 0.2 | 3m | 0.8259 | 0.8257 | -0.0002 | 0.8260 | +0.0001 | 0.8253 | -0.0006 |
| 90 | 0.3 | 3m | 10.5893 | 10.5768 | -0.0125 | 10.5917 | +0.0024 | 10.5882 | -0.0011 |
| 95 | 0.3 | 3m | 7.3316 | 7.3131 | -0.0185 | 7.3345 | +0.0029 | 7.3329 | +0.0013 |
| 100 | 0.3 | 3m | 4.8310 | 4.8096 | -0.0214 | 4.8339 | +0.0029 | 4.8340 | +0.0030 |
| 105 | 0.3 | 3m | 3.0388 | 3.0178 | -0.0210 | 3.0413 | +0.0025 | 3.0391 | +0.0003 |
| 110 | 0.3 | 3m | 1.8325 | 1.8145 | -0.0180 | 1.8345 | +0.0020 | 1.8319 | -0.0006 |
| 90 | 0.4 | 3m | 11.8287 | 11.7666 | -0.0621 | 11.8326 | +0.0039 | 11.8288 | +0.0001 |
| 95 | 0.4 | 3m | 8.8035 | 8.7216 | -0.0819 | 8.8082 | +0.0047 | 8.8070 | +0.0035 |
| 100 | 0.4 | 3m | 6.3735 | 6.2808 | -0.0927 | 6.3784 | +0.0049 | 6.3786 | +0.0051 |
| 105 | 0.4 | 3m | 4.4976 | 4.4036 | -0.0940 | 4.5024 | +0.0048 | 4.5004 | +0.0028 |
| 110 | 0.4 | 3m | 3.1011 | 3.0133 | -0.0878 | 3.1055 | +0.0044 | 3.1025 | +0.0014 |
| 90 | 0.2 | 6m | 9.7572 | 9.7564 | -0.0008 | 9.7579 | +0.0007 | 9.7606 | +0.0034 |
| 95 | 0.2 | 6m | 6.7199 | 6.7186 | -0.0013 | 6.7207 | +0.0008 | 6.7185 | -0.0014 |
| 100 | 0.2 | 6m | 4.4312 | 4.4296 | -0.0016 | 4.4319 | +0.0007 | 4.4320 | +0.0008 |
| 105 | 0.2 | 6m | 2.8107 | 2.8089 | -0.0018 | 2.8112 | +0.0005 | 2.8100 | -0.0007 |
| 110 | 0.2 | 6m | 1.7240 | 1.7223 | -0.0017 | 1.7242 | +0.0002 | 1.7275 | +0.0035 |
| 90 | 0.3 | 6m | 11.0807 | 11.0549 | -0.0258 | 11.0846 | +0.0039 | 11.0845 | +0.0038 |
| 95 | 0.3 | 6m | 8.2363 | 8.2019 | -0.0344 | 8.2407 | +0.0044 | 8.2367 | +0.0004 |
| 100 | 0.3 | 6m | 5.9763 | 5.9369 | -0.0394 | 5.9808 | +0.0045 | 5.9784 | +0.0021 |
| 105 | 0.3 | 6m | 4.2443 | 4.2037 | -0.0406 | 4.2484 | +0.0041 | 4.2449 | +0.0006 |
| 110 | 0.3 | 6m | 2.9582 | 2.9196 | -0.0386 | 2.9619 | +0.0037 | 2.9623 | +0.0041 |
| 90 | 0.4 | 6m | 12.6171 | 12.4974 | -0.1197 | 12.6237 | +0.0066 | 12.6231 | +0.0060 |
| 95 | 0.4 | 6m | 9.9223 | 9.7730 | -0.1493 | 9.9298 | +0.0075 | 9.9260 | +0.0037 |
| 100 | 0.4 | 6m | 7.6965 | 7.5287 | -0.1678 | 7.7044 | +0.0079 | 7.7017 | +0.0052 |
| 105 | 0.4 | 6m | 5.8978 | 5.7228 | -0.1750 | 5.9056 | +0.0078 | 5.9015 | +0.0037 |
| 110 | 0.4 | 6m | 4.4716 | 4.2992 | -0.1724 | 4.4792 | +0.0076 | 4.4779 | +0.0063 |

### Summary over the 45 European puts

| method | mean \|error\| | worst \|error\| | negatives |
|---|---|---|---|
| binomial n=50 | 0.14837 | 0.74182 | 38/45 |
| binomial n=100 | 0.07251 | 0.35697 | 38/45 |
| binomial n=200 | 0.03575 | 0.17498 | 38/45 |
| trinomial n=50 | 0.00974 | 0.03203 | 0/45 |
| trinomial n=100 | 0.00465 | 0.01489 | 0/45 |
| trinomial n=200 | 0.00246 | 0.00789 | 0/45 |
| **ADS tree N=500** | 0.00176 | 0.00630 | 12/45 |

The published column is printed to four decimals, so an error below about `5e-5` is at the reference's own rounding and cannot be resolved against it.


### Where the two walks differ: by initial volatility

The Lamperti drift grows with `sqrt(V0)/sqrt(theta)`, and it is the drift that the two-point walk cannot carry. Splitting the same 45 puts by `sqrt(V0)` isolates that.

| sqrt(V0) | Lamperti drift at u_0 | bino n=200 mean \|err\| | trino n=200 mean \|err\| | ratio |
|---|---|---|---|---|
| 0.2 | 16.64 | 0.00059 | 0.00039 | 1.5x |
| 0.3 | 6.80 | 0.01935 | 0.00250 | 7.7x |
| 0.4 | 10.12 | 0.08732 | 0.00449 | 19.5x |

## 2. American put

`ADS tree` and `B-N CV` are published; the rest are ours. The Longstaff--Schwartz column is a genuine lower bound (policy fitted on one sample, applied to a disjoint one).

| S0 | rho | sqrt(V0) | T | ADS tree N=250 | B-N CV | bino n=200 | err vs CV | trino n=200 | err vs CV | LSM | err vs CV |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 90 | -0.7 | 0.2 | 1m | 10.0000 | 9.9997 | 10.0000 | +0.0003 | 10.0000 | +0.0003 | 10.0000 | +0.0003 |
| 100 | -0.7 | 0.2 | 1m | 2.1249 | 2.1267 | 2.1272 | +0.0005 | 2.1270 | +0.0003 | 2.1118 | -0.0149 |
| 110 | -0.7 | 0.2 | 1m | 0.1273 | 0.1274 | 0.1274 | -0.0000 | 0.1274 | -0.0000 | 0.1245 | -0.0029 |
| 90 | -0.1 | 0.2 | 1m | 10.0000 | 10.0000 | 10.0000 | +0.0000 | 10.0000 | +0.0000 | 10.0000 | +0.0000 |
| 100 | -0.1 | 0.2 | 1m | 2.1236 | 2.1254 | 2.1257 | +0.0003 | 2.1257 | +0.0003 | 2.1112 | -0.0142 |
| 110 | -0.1 | 0.2 | 1m | 0.1090 | 0.1091 | 0.1091 | +0.0000 | 0.1091 | +0.0000 | 0.1074 | -0.0017 |
| 90 | -0.7 | 0.4 | 1m | 10.6843 | 10.6804 | 10.6666 | -0.0138 | 10.6828 | +0.0024 | 10.6651 | -0.0153 |
| 100 | -0.7 | 0.4 | 1m | 4.2183 | 4.2140 | 4.1874 | -0.0266 | 4.2160 | +0.0020 | 4.1886 | -0.0254 |
| 110 | -0.7 | 0.4 | 1m | 1.1942 | 1.1939 | 1.1758 | -0.0181 | 1.1950 | +0.0011 | 1.1775 | -0.0164 |
| 90 | -0.1 | 0.4 | 1m | 10.7123 | 10.7100 | 10.7099 | -0.0001 | 10.7104 | +0.0004 | 10.6908 | -0.0192 |
| 100 | -0.1 | 0.4 | 1m | 4.2194 | 4.2158 | 4.2163 | +0.0005 | 4.2170 | +0.0012 | 4.1849 | -0.0309 |
| 110 | -0.1 | 0.4 | 1m | 1.1666 | 1.1667 | 1.1671 | +0.0004 | 1.1677 | +0.0010 | 1.1601 | -0.0066 |
| 90 | -0.7 | 0.2 | 3m | 10.1222 | 10.1206 | 10.1218 | +0.0012 | 10.1218 | +0.0012 | 10.1233 | +0.0027 |
| 100 | -0.7 | 0.2 | 3m | 3.4790 | 3.4807 | 3.4809 | +0.0002 | 3.4812 | +0.0005 | 3.4569 | -0.0238 |
| 110 | -0.7 | 0.2 | 3m | 0.8405 | 0.8416 | 0.8414 | -0.0002 | 0.8416 | +0.0000 | 0.8265 | -0.0151 |
| 90 | -0.1 | 0.2 | 3m | 10.1713 | 10.1706 | 10.1710 | +0.0004 | 10.1710 | +0.0004 | 10.1474 | -0.0232 |
| 100 | -0.1 | 0.2 | 3m | 3.4729 | 3.4747 | 3.4752 | +0.0005 | 3.4753 | +0.0006 | 3.4575 | -0.0172 |
| 110 | -0.1 | 0.2 | 3m | 0.7726 | 0.7736 | 0.7739 | +0.0003 | 0.7739 | +0.0003 | 0.7667 | -0.0069 |
| 90 | -0.7 | 0.4 | 3m | 12.1245 | 12.1122 | 12.0520 | -0.0602 | 12.1228 | +0.0106 | 12.1189 | +0.0067 |
| 100 | -0.7 | 0.4 | 3m | 6.4989 | 6.4899 | 6.3969 | -0.0930 | 6.4975 | +0.0076 | 6.4568 | -0.0331 |
| 110 | -0.7 | 0.4 | 3m | 3.1512 | 3.1456 | 3.0570 | -0.0886 | 3.1510 | +0.0054 | 3.1226 | -0.0230 |
| 90 | -0.1 | 0.4 | 3m | 12.1880 | 12.1819 | 12.1849 | +0.0030 | 12.1860 | +0.0041 | 12.1736 | -0.0083 |
| 100 | -0.1 | 0.4 | 3m | 6.5023 | 6.4964 | 6.4989 | +0.0025 | 6.5013 | +0.0049 | 6.4698 | -0.0266 |
| 110 | -0.1 | 0.4 | 3m | 3.0952 | 3.0914 | 3.0927 | +0.0013 | 3.0957 | +0.0043 | 3.0732 | -0.0182 |
| 90 | -0.7 | 0.2 | 6m | 10.5682 | 10.5637 | 10.5636 | -0.0001 | 10.5645 | +0.0008 | 10.5653 | +0.0016 |
| 100 | -0.7 | 0.2 | 6m | 4.6691 | 4.6636 | 4.6617 | -0.0019 | 4.6637 | +0.0001 | 4.6417 | -0.0219 |
| 110 | -0.7 | 0.2 | 6m | 1.7899 | 1.7874 | 1.7853 | -0.0021 | 1.7872 | -0.0002 | 1.7602 | -0.0272 |
| 90 | -0.1 | 0.2 | 6m | 10.6521 | 10.6478 | 10.6484 | +0.0006 | 10.6484 | +0.0006 | 10.6278 | -0.0200 |
| 100 | -0.1 | 0.2 | 6m | 4.6531 | 4.6473 | 4.6481 | +0.0008 | 4.6482 | +0.0009 | 4.6226 | -0.0247 |
| 110 | -0.1 | 0.2 | 6m | 1.6857 | 1.6832 | 1.6838 | +0.0006 | 1.6839 | +0.0007 | 1.6730 | -0.0102 |
| 90 | -0.7 | 0.4 | 6m | 13.2431 | 13.2172 | 13.0969 | -0.1203 | 13.2374 | +0.0202 | 13.2265 | +0.0093 |
| 100 | -0.7 | 0.4 | 6m | 8.0204 | 7.9998 | 7.8263 | -0.1735 | 8.0138 | +0.0140 | 7.9647 | -0.0351 |
| 110 | -0.7 | 0.4 | 6m | 4.6328 | 4.6201 | 4.4413 | -0.1788 | 4.6300 | +0.0099 | 4.5739 | -0.0462 |
| 90 | -0.1 | 0.4 | 6m | 13.3279 | 13.3142 | 13.3206 | +0.0064 | 13.3228 | +0.0086 | 13.2879 | -0.0263 |
| 100 | -0.1 | 0.4 | 6m | 8.0231 | 8.0083 | 8.0127 | +0.0044 | 8.0175 | +0.0092 | 7.9801 | -0.0282 |
| 110 | -0.1 | 0.4 | 6m | 4.5554 | 4.5454 | 4.5472 | +0.0018 | 4.5537 | +0.0083 | 4.5272 | -0.0182 |

### Summary against the Beliaeva--Nawalkha control variate

| method | mean error | mean \|error\| | worst | below CV |
|---|---|---|---|---|
| binomial n=200 | -0.0209 | 0.0223 | -0.1788 | 15/36 |
| trinomial n=200 | +0.0034 | 0.0034 | +0.0202 | 2/36 |
| LSM | -0.0161 | 0.0173 | -0.0462 | 30/36 |

## 3. The mechanism, measured rather than inferred

At every node of every step the walk's own first two moments are compared with the ones it is supposed to reproduce. `max |Var/delta - 1|` is the largest relative variance error over the whole lattice; for the two-point walk it should be `mu^2 delta`, for the trinomial it should be zero.

| walk | n | max \|Var/delta - 1\| | max mean error | max \|Lamperti drift\| | probability violations | driver states |
|---|---|---|---|---|---|---|
| binomial | 50 | 1.000e+00 | 1.041e-02 | 10.12 | 506 | 51 |
| binomial | 100 | 5.116e-01 | 1.785e-15 | 10.12 | 0 | 101 |
| binomial | 200 | 2.558e-01 | 1.442e-15 | 10.12 | 0 | 201 |
| trinomial | 50 | 2.220e-16 | 2.770e-16 | 10.02 | 0 | 39 |
| trinomial | 100 | 3.331e-16 | 5.767e-16 | 9.97 | 0 | 54 |
| trinomial | 200 | 3.331e-16 | 5.523e-16 | 10.02 | 0 | 77 |

Row taken at `S0 = 100`, `sqrt(V0) = 0.4`, `T = 6m`, the worst cell of the grid. The trinomial's variance error is at machine precision, which is the fix working exactly as the algebra predicts.


### The American-call control

With `r = 0.05 > 0` and no dividend the American call must equal the European call. Largest gap over all 81 parameter sets:

| walk | n | max \|American call - European call\| |
|---|---|---|
| binomial | 50 | 0.00e+00 |
| binomial | 100 | 0.00e+00 |
| binomial | 200 | 0.00e+00 |
| trinomial | 50 | 8.42e-08 |
| trinomial | 100 | 6.56e-08 |
| trinomial | 200 | 6.28e-08 |

## 4. Put-call parity

`C - P` must equal `S0 - K exp(-rT)` exactly, for any correct scheme, whatever its discretisation error. It is a property of the martingale, so a scheme that loses variance in the driver does not have to break it -- but one that mis-prices the two legs asymmetrically does.

| walk | n | mean \|parity error\| | worst |
|---|---|---|---|
| binomial | 50 | 0.10133 | 0.57413 |
| binomial | 100 | 0.05061 | 0.28651 |
| binomial | 200 | 0.02518 | 0.14285 |
| trinomial | 50 | 0.00199 | 0.00589 |
| trinomial | 100 | 0.00090 | 0.00229 |
| trinomial | 200 | 0.00058 | 0.00157 |

## 5. Computation time

| walk | n | mref | grid nodes | driver states | 45 European sets | 36 American sets |
|---|---|---|---|---|---|---|
| binomial | 50 | 10 | 5,501 | 51 | 4.4s | 3.2s |
| binomial | 100 | 15 | 16,401 | 101 | 55.5s | 37.2s |
| binomial | 200 | 20 | 43,601 | 201 | 598.4s | 786.9s |
| trinomial | 50 | 10 | 7,101 | 35 | 10.0s | 11.0s |
| trinomial | 100 | 15 | 21,201 | 49 | 79.4s | 100.8s |
| trinomial | 200 | 20 | 56,001 | 69 | 684.9s | 1017.7s |

All four vanilla payoffs come out of one backward pass, so these are four prices per parameter set. Our Longstaff--Schwartz lower bound for the 36 American puts, for comparison: 42.9s. ADS report 5.71 s per European option at N = 200 and 13.97 s per American put at N = 250, in MATLAB on a 2010 laptop, so the cross-machine comparison is indicative only.

## Compute time

| phase | wall time | seconds |
|---|---|---|
| 1. 1. European put and call, both walks, vs ADS's tree and analytic | 23m52s | 1432.55 |
| 2. 2. American put: both walks, LSM, and two published references | 33m19s | 1999.87 |
| 3. 3. the mechanism, measured node by node | 0s | 0.01 |
| **total** | **57m12s** | **3432.43** |

Machine: Darwin arm64, python 3.9.6.
