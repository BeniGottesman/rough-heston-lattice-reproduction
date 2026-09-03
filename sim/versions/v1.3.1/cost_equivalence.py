"""
cost_equivalence — what a lattice price costs, and what the same accuracy costs
by Monte-Carlo.

The question this answers.  A lattice price at n steps is a single deterministic
number carrying a BIAS.  A Monte-Carlo price is unbiased (with the exact-
covariance driver) but carries a STANDARD ERROR that shrinks like 1/sqrt(N).  So
"how do they compare" only means something once the accuracies are matched:

    N*(target) = (sigma_eff / target)^2,      target = |lattice bias|,

where sigma_eff is the per-path standard deviation of whatever estimator the
Monte-Carlo actually uses.  Then the two wall times can be put side by side.

Three things have to be measured rather than assumed, or the answer is wrong.

  sigma_eff       The reference Monte-Carlo uses ANTITHETIC pairs and an
                  eta = 0 Black--Scholes CONTROL VARIATE.  Its own reported
                  `stderr` is computed as sqrt(Var(sample)/paths), which treats
                  the antithetic pairs as independent — they are negatively
                  correlated, so that formula is conservative.  sigma_eff is
                  therefore measured EMPIRICALLY, from the scatter of R
                  independent replications, which captures both devices exactly
                  and needs no formula.  Both are reported, so the size of the
                  discrepancy is visible instead of assumed away.
  setup vs path   The Monte-Carlo pays a large ONE-OFF cost: `driver_factor`
                  factorises an nfine x nfine covariance matrix (nfine = 512
                  here).  At small N that setup dominates the run completely —
                  20k paths measured slower than 200k on a cold cache.  Setup
                  and marginal path cost are timed separately, and the
                  comparison is reported both ways.
  the bias        |lattice - truth| needs a "truth", and the truth here is
                  itself a Monte-Carlo with a standard error.  If the bias is
                  not several reference standard errors, N* is not resolved and
                  a band is reported instead of a number.

What this comparison is NOT
---------------------------
It compares THESE IMPLEMENTATIONS, not the two methods.  The Monte-Carlo is
vectorised numpy whose inner loop is a BLAS matrix product (multi-threaded); the
lattices are Python loops over a state space, single-threaded, with no C++ port
for Route B.  A constant factor of tens is implementation, not mathematics.  The
algorithmic scalings are reported alongside so a reader can rescale:

    Monte-Carlo   O(N x nfine)                     per price
    one-step      O(n x nx)                        per price
    lift, m       O(n x nx x prod_i N_i)           per price, state space grows

And the lattice's real advantage is not priced here at all: it computes an
AMERICAN price in the same backward pass, where Monte-Carlo needs
Longstaff--Schwartz — more expensive, and itself biased.  The American lattice
time is measured; no Monte-Carlo equivalent is offered, because a fair one would
require an LSM estimator for rough Bergomi, which this project has only for
rough Heston.
"""
from __future__ import annotations

import math
import time
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mc_reference as mcr                                      # noqa: E402
import route_aprime as ra                                       # noqa: E402
from route_b_lattice import LiftedLattice                        # noqa: E402


# the project's documented convention, kept for comparability with Parts V-VI
def mref_for(n: int) -> int:
    return max(4, int(math.ceil(4.0 * math.sqrt(n / 8.0))))


def zmax_for(H: float) -> float:
    return 3.0 / math.sqrt(2.0 * H)


# ------------------------------------------------------------------ timing
def timed(fn, repeats: int = 3) -> tuple[float, float, object]:
    """Run fn repeats times; return (best, median, last result).

    Best-of is the honest number for a deterministic computation: it is the one
    least polluted by other load on the machine.  The median is reported too, so
    a big gap between them is visible rather than hidden.
    """
    ts, out = [], None
    for _ in range(max(1, repeats)):
        t0 = time.perf_counter()
        out = fn()
        ts.append(time.perf_counter() - t0)
    return min(ts), float(np.median(ts)), out


# ------------------------------------------------------------- the pricers
def price_onestep(n: int, H: float, eta: float, rho: float) -> float:
    return float(ra.route_aprime_european_put(
        n, H, eta, rho, zmax=zmax_for(H), mref=mref_for(n)))


def price_lift(n: int, H: float, eta: float, rho: float, m: int,
               american: bool = False) -> dict:
    lat = LiftedLattice(n, H, eta, rho, m, mref=mref_for(n))
    r = lat.price(american=american)
    return {"price": float(r["value"]), "state_space": int(r["state_space"]),
            "nx": int(r["nx"])}


# -------------------------------------------------------- Monte-Carlo cost
def mc_setup_cost(H: float, nfine: int) -> float:
    """Wall time of the one-off covariance factorisation, measured cold.

    `mc_reference.driver_factor` caches, so this has to be timed on a cache the
    caller has not yet warmed; the runner calls this before anything else.
    """
    t0 = time.perf_counter()
    mcr.driver_factor(H, nfine)
    return time.perf_counter() - t0


def mc_path_rate(H: float, eta: float, rho: float, nfine: int, paths: int,
                 seed: int = 101) -> float:
    """Marginal paths per second, with the setup already warm."""
    mcr.driver_factor(H, nfine)                     # ensure warm
    t0 = time.perf_counter()
    mcr.european_put_mc(H, eta, rho, nfine=nfine, paths=paths, seed=seed)
    return paths / (time.perf_counter() - t0)


def mc_sigma_eff(H: float, eta: float, rho: float, nfine: int, paths: int,
                 reps: int, seed0: int = 1000) -> dict:
    """sigma_eff from the scatter of `reps` independent runs of `paths` each.

    se(N) = sigma_eff / sqrt(N) by definition of sigma_eff, so
    sigma_eff = std(replication means) * sqrt(paths).  This is the estimator's
    TRUE spread: antithetic pairing and the control variate are both inside it,
    with no formula standing in for them.  The code's own `stderr` is returned
    alongside for comparison.
    """
    means, reported = [], []
    for i in range(reps):
        r = mcr.european_put_mc(H, eta, rho, nfine=nfine, paths=paths,
                                seed=seed0 + 17 * i)
        means.append(r["price"])
        reported.append(r["stderr"])
    means = np.asarray(means)
    se_rep = float(means.std(ddof=1))               # se of ONE paths-run
    return {"sigma_eff": se_rep * math.sqrt(paths),
            "se_at_paths": se_rep,
            "sigma_eff_reported": float(np.mean(reported)) * math.sqrt(paths),
            "pooled_mean": float(means.mean()),
            "pooled_se": se_rep / math.sqrt(reps),
            "reps": reps, "paths_each": paths}


def sigma_eff_direct(H: float, eta: float, rho: float, nfine: int,
                     paths: int = 40_000, seed: int = 4242) -> dict:
    """sigma_eff EXACTLY, from one sample, via the antithetic pair correlation.

    The replication estimator above is honest but weak: a standard deviation
    estimated from R replications carries 1/sqrt(2(R-1)) relative uncertainty —
    15% at R = 24 — which propagates squared into N*.  There is a far better
    route, because the only thing the reference's own `stderr` formula gets wrong
    is the antithetic pairing.  With N paths laid out as N/2 antithetic pairs,

        Var(mean) = Var(sample) (1 + rho_pair) / N,

    so   sigma_eff = sd(sample) * sqrt(1 + rho_pair),

    and both factors are measured from a single large sample with negligible
    uncertainty.  rho_pair = 0 recovers the code's formula; rho_pair < 0 is the
    antithetic device working; rho_pair > 0 means it is COUNTERPRODUCTIVE and the
    quoted standard error is optimistic.

    The estimator's per-path sample is rebuilt here rather than extracted from
    `mc_reference.european_put_mc`, which returns only aggregates.  That
    duplication is a risk, so the runner cross-validates it.  Note what the
    check can and cannot be: the pricer draws in CHUNKS from one rng stream,
    while this rebuild draws the whole block at once, so the two consume the
    stream differently and are DIFFERENT draws from the same distribution.  The
    price therefore agrees only to within a couple of standard errors; what must
    agree tightly is `stderr_formula`, since that estimates a property of the
    distribution rather than of the particular draw.
    """
    L = mcr.driver_factor(H, nfine)
    df = mcr.T / nfine
    ts = np.arange(1, nfine + 1) * df
    half = paths // 2
    rng = np.random.default_rng(seed)
    Z = rng.standard_normal((half, nfine))
    Z = np.concatenate([Z, -Z], axis=0)
    Zp = rng.standard_normal((half, nfine))
    Zp = np.concatenate([Zp, -Zp], axis=0)

    with np.errstate(all="ignore"):
        VH = Z @ L.T
    v = mcr.XI0 * np.exp(eta * VH - 0.5 * eta ** 2 * ts[None, :] ** (2.0 * H))
    v = np.concatenate([np.full((2 * half, 1), mcr.XI0), v[:, :-1]], axis=1)
    dW = Z * math.sqrt(df)
    dB = rho * dW + math.sqrt(1.0 - rho ** 2) * Zp * math.sqrt(df)
    sB = dB.sum(axis=1)
    logST = (math.log(mcr.S0) - 0.5 * (v * df).sum(axis=1)
             + (np.sqrt(v) * dB).sum(axis=1))
    pay = np.maximum(mcr.KSTRIKE - np.exp(logST), 0.0)
    logflat = (math.log(mcr.S0) - 0.5 * mcr.XI0 * mcr.T
               + math.sqrt(mcr.XI0) * sB)
    payflat = np.maximum(mcr.KSTRIKE - np.exp(logflat), 0.0)

    out = {}
    for name, s in (("raw", pay), ("control", pay - payflat)):
        a, b = s[:half], s[half:]
        rho_pair = float(np.corrcoef(a, b)[0, 1])
        sd = float(s.std(ddof=1))
        out[name] = {"sd": sd, "rho_pair": rho_pair,
                     "inflation": math.sqrt(max(0.0, 1.0 + rho_pair)),
                     "sigma_eff": sd * math.sqrt(max(0.0, 1.0 + rho_pair))}
    return {**out, "paths": 2 * half, "nfine": nfine,
            "mean": float((pay - payflat).mean()) + mcr.bs_put(),
            "stderr_formula": float((pay - payflat).std(ddof=1)
                                    / math.sqrt(2 * half))}


def sigma_scaling_check(H: float, eta: float, rho: float, nfine: int,
                        paths_small: int, paths_large: int, reps: int,
                        seed0: int = 5000) -> dict:
    """Does se(N) really scale like 1/sqrt(N)?

    Every N* in this module is an EXTRAPOLATION from a measured sigma_eff, and
    that extrapolation is only valid if sigma_eff is the same at both ends.  It
    need not be: an antithetic estimator's variance per path depends on how the
    pairs are laid out inside a chunk, and the control variate's effectiveness
    can drift with the sample.  So it is measured at two path counts and the
    ratio reported.  A ratio far from 1 invalidates the N* column, and the run
    says so rather than quietly extrapolating.
    """
    a = mc_sigma_eff(H, eta, rho, nfine, paths_small, reps, seed0)
    b = mc_sigma_eff(H, eta, rho, nfine, paths_large, max(4, reps // 2),
                     seed0 + 7919)
    return {"paths_small": paths_small, "paths_large": paths_large,
            "sigma_small": a["sigma_eff"], "sigma_large": b["sigma_eff"],
            "ratio": (b["sigma_eff"] / a["sigma_eff"]
                      if a["sigma_eff"] > 0 else float("nan")),
            "reps_small": a["reps"], "reps_large": b["reps"],
            # relative uncertainty of a std estimated from k samples ~ 1/sqrt(2(k-1))
            "rel_unc_small": 1.0 / math.sqrt(2.0 * max(1, a["reps"] - 1)),
            "rel_unc_large": 1.0 / math.sqrt(2.0 * max(1, b["reps"] - 1))}


def paths_for_target(sigma_eff: float, target: float) -> float:
    """N* such that the Monte-Carlo standard error equals `target`."""
    if target <= 0 or not math.isfinite(target):
        return float("inf")
    return (sigma_eff / target) ** 2


def equivalence(lattice_price: float, lattice_time: float,
                ref_price: float, ref_se: float,
                sigma_eff: float, path_rate: float,
                setup_time: float) -> dict:
    """The comparison, with the reference's own uncertainty propagated.

    The bias is |lattice - reference| and the reference is itself noisy, so the
    bias is only known to +-2*ref_se.  N* is reported for the bias and for both
    ends of that band; when the bias is not resolved (|bias| < 2*ref_se) the
    point value is meaningless and `resolved` says so.
    """
    bias = lattice_price - ref_price
    ab = abs(bias)
    lo, hi = max(ab - 2.0 * ref_se, 0.0), ab + 2.0 * ref_se
    n_star = paths_for_target(sigma_eff, ab)
    out = {
        "bias": bias, "abs_bias": ab,
        "bias_in_ref_se": ab / ref_se if ref_se > 0 else float("inf"),
        "resolved": ab > 2.0 * ref_se,
        "n_star": n_star,
        "n_star_lo": paths_for_target(sigma_eff, hi),   # smaller bias -> more N
        "n_star_hi": paths_for_target(sigma_eff, lo),
        "mc_time_paths_only": n_star / path_rate if path_rate > 0 else float("inf"),
        "mc_time_with_setup": (n_star / path_rate + setup_time
                               if path_rate > 0 else float("inf")),
        "lattice_time": lattice_time,
    }
    out["ratio_paths_only"] = (out["mc_time_paths_only"] / lattice_time
                               if lattice_time > 0 else float("nan"))
    out["ratio_with_setup"] = (out["mc_time_with_setup"] / lattice_time
                               if lattice_time > 0 else float("nan"))
    return out


def breakeven_accuracy(sigma_eff: float, path_rate: float,
                       lattice_time: float) -> float:
    """The accuracy at which Monte-Carlo (paths only) costs the lattice's time.

    Below this error the Monte-Carlo is the cheaper way to get the number; above
    it the lattice is.  Reported as an absolute price error.
    """
    n = lattice_time * path_rate
    return sigma_eff / math.sqrt(n) if n > 0 else float("inf")
