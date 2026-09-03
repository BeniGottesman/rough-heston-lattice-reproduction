"""verify_mc_band — O018: the reference Monte-Carlo band must be measured on the
antithetic PAIRS, not on the paths.

Four checks.  T1 and T2 are algebra and cost nothing.  T3 is the regression guard that
FAILS ON THE OLD FORMULA -- that is the whole point of this file.  T4 is the invariant
that proves the repair did not move the price.

    python3 verify/verify_mc_band.py
"""
from __future__ import annotations

import math
import sys

import numpy as np

sys.path.insert(0, "sim")
import mc_reference as MR                                     # noqa: E402

FAIL: list[str] = []


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"  -- {detail}" if detail else ""))
    if not ok:
        FAIL.append(name)


# ------------------------------------------------------------------ T1, algebra
# Var(A) = Var(s)(1+rho)/2 for A = (s_1 + s_2)/2, so SE on the mean of n pairs is
# SE_iid * sqrt(1+rho).  Verified on a synthetic sample with a KNOWN rho.
print("T1  the identity SE_pair = SE_iid * sqrt(1 + rho), on synthetic data")
rng = np.random.default_rng(11)
worst = 0.0
for rho_t in (-0.9, -0.5, 0.0, 0.25, 0.75):
    n = 400_000
    g = rng.standard_normal(n)
    e = rng.standard_normal(n)
    a = g
    b = rho_t * g + math.sqrt(1.0 - rho_t ** 2) * e           # corr(a, b) = rho_t
    s = np.concatenate([a, b])
    pair = 0.5 * (a + b)
    se_iid = s.std(ddof=0) / math.sqrt(s.size)
    se_pair = pair.std(ddof=0) / math.sqrt(pair.size)
    pred = se_iid * math.sqrt(1.0 + float(np.corrcoef(a, b)[0, 1]))
    worst = max(worst, abs(se_pair - pred) / pred)
check("T1 identity holds", worst < 5e-3, f"worst relative deviation {worst:.2e}")

# ------------------------------------------------------- T2, the sign of the effect
print("\nT2  the direction: rho > 0 widens the band, rho < 0 narrows it")
check("T2 rho>0 widens", math.sqrt(1 + 0.1931) > 1.0, "sqrt(1.1931) = 1.0923")
check("T2 rho<0 narrows", math.sqrt(1 - 0.5912) < 1.0, "sqrt(0.4088) = 0.6394")

# --------------------------------------- T3, the regression guard on the real pricer
# With the control variate ON, the within-pair correlation is POSITIVE, so the correct
# band is STRICTLY WIDER than the i.i.d. one.  The old code returned the i.i.d. figure,
# so `stderr == stderr_iid_SUPERSEDED` and this check FAILS on it.
print("\nT3  the pricer returns the PAIR band, not the i.i.d. band  (fails on the old code)")
r = MR.european_put_mc(H=0.1, eta=0.3, rho=-0.7, nfine=48, paths=20_000,
                       chunk=10_000, seed=7, control=True)
for k in ("stderr_iid_SUPERSEDED", "rho_pair", "band_inflation", "variance_units"):
    check(f"T3 the diagnostic {k} is reported", k in r)
if "stderr_iid_SUPERSEDED" in r:
    check("T3 rho_pair is positive with the control variate on",
          r["rho_pair"] > 0.0, f"rho_pair = {r['rho_pair']:+.4f}")
    check("T3 the reported band is strictly wider than the i.i.d. one",
          r["stderr"] > r["stderr_iid_SUPERSEDED"] * 1.01,
          f"inflation = {r['band_inflation']:.4f}")
    check("T3 the inflation equals sqrt(1 + rho_pair)",
          abs(r["band_inflation"] - math.sqrt(1.0 + r["rho_pair"])) < 1e-9,
          f"{r['band_inflation']:.9f} vs {math.sqrt(1.0+r['rho_pair']):.9f}")
    check("T3 the number of variance units is the number of pairs",
          r["variance_units"] == r["paths"] // 2,
          f"{r['variance_units']} pairs for {r['paths']} paths")

# ------------------------------------------- T4, the price did not move: eta = 0
# At eta = 0 the rough and flat payoffs coincide pathwise, so the control-variate
# estimator must return the Black--Scholes value with a band of exactly zero.  This
# invariant sees the MEAN, and it is what proves the repair touched only the band.
print("\nT4  the price is untouched: at eta = 0 the estimator is BS with a zero band")
r0 = MR.european_put_mc(H=0.1, eta=0.0, rho=-0.7, nfine=32, paths=8_000,
                        chunk=8_000, seed=3, control=True)
check("T4 price equals the exact Black-Scholes put",
      abs(r0["price"] - MR.bs_put()) < 1e-12,
      f"|{r0['price']:.12f} - {MR.bs_put():.12f}| = {abs(r0['price']-MR.bs_put()):.2e}")
# NOT bitwise zero, and the code is right: at eta = 0 the two legs compute the same
# number by DIFFERENT summation orders -- (XI0*df).sum() against XI0*T, and
# (sqrt(XI0)*dB).sum() against sqrt(XI0)*dB.sum() -- so the pathwise difference is a
# few ulp rather than 0.  The first version of this check asserted == 0.0 and failed at
# 1.03e-16; the assertion was wrong, not the pricer, and it is recorded here rather
# than quietly loosened.  The MEAN, which is what this test exists to protect, matches
# the exact Black-Scholes value to 0.00e+00.
check("T4 the band is zero to floating point", r0["stderr"] < 1e-12,
      f"stderr = {r0['stderr']:.3e}  (summation order, not a real band)")

print(f"\n{'FAILED: ' + ', '.join(FAIL) if FAIL else 'all checks passed'}")
sys.exit(1 if FAIL else 0)
