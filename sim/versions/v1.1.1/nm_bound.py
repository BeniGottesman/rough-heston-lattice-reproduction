"""
nm_bound — the link between the time step and the number of lift factors.

Why the two are linked at all
-----------------------------
The one-step recombining lattice carries the roughness driver in a single
counter: up-steps minus down-steps.  For a standard process that is enough — how
much it moves depends only on how many steps there were.  A rough process is
different: its increment depends on its whole past with DECAYING weights (the
kernel K(u) = u^{H-1/2}), and a net counter knows only the balance, not who moved
when.  Refine the grid and the counter accumulates more and more independent
±1's, so its spread grows like n^{1-2H} while the true spread is fixed.  That is
Proposition 8.3, and it is why refining makes the one-step scheme WORSE.

Route B replaces the single counter by m Ornstein--Uhlenbeck factors, each one a
counter that forgets its own past at its own rate s_i.  A stiff factor (large
s_i) remembers only the very recent past; a soft one reaches far back.  With the
right weights, m factors reproduce the weighted memory the rough kernel needs.

The link is then forced by two facts:

  * the kernel's detail lives at SHORT lags — K is nearly singular at u = 0 — and
    the shortest lag a grid of n steps ever evaluates is delta = T/n.  Following
    the memory down to that lag needs a factor of stiffness s ~ 1/delta = n/T;
  * the factors have to cover every time scale from delta up to T, and each
    factor covers a bounded RATIO of scales.  Covering a range of n at a fixed
    ratio r per factor costs m ~ log(n)/log(r) factors.

So m grows like log n, not like n.  Doubling the number of time steps costs a
fraction of one extra factor.  This is the mechanism behind `F025`.

Where the exponentials come from
--------------------------------
The reason a finite sum of exponentials can do this at all is the completely
monotone (Laplace) representation used by Abi Jaber--El Euch,

    K(t) = t^h / Gamma(1+h) = int_0^infty e^{-s t} mu_h(ds),
    mu_h(ds) = c_h s^{-h-1} ds,      h = H - 1/2 in (-1/2, 0),

so an m-factor lift is exactly an m-point QUADRATURE of mu_h, and `m` is a node
count.  mu_h is integrable at s = 0 and infinite at s = infinity, so what a
finite m must sacrifice is the TOP of the s-range — the stiff factors — which is
precisely the short-lag detail the fine grid needs.  Hence the coupling.

Two constructions, two sides of the bound
-----------------------------------------
  `geometric_lift`  an EXPLICIT geometric partition of the s-range from ~1/T to
                    ~1/delta.  No optimisation: m in, error out.  This is the
                    "m is sufficient" side, and the one whose error can be
                    bounded by hand, since the ratio r between consecutive nodes
                    is the single accuracy knob and m = 1 + log(n)/log(r).
  `route_b.best_lift`  nodes free, weights non-negative least squares — the
                    achievable floor.  No m-factor lift with positive weights
                    does better (up to the local search), so the m it needs is
                    the empirical "m is necessary" side.

`route_b.ajee_optimised` sits between them: the AJEE construction with its
two-parameter partition family optimised numerically.

The error criterion
-------------------
Everywhere here the error is the RELATIVE L^2 kernel error on the lags the
scheme actually sees,

    rel(m) = || K - K^m ||_{L^2(delta,T)} / || K ||_{L^2(delta,T)},

with delta = T/n.  Relative, because both K and K^m are linear in the weights:
the ratio is unchanged by the normalisation convention (u^h/Gamma(1+h) versus
the sqrt(2H) u^h used by the rough Bergomi code), so a tolerance means the same
thing everywhere in the project.  On (delta, T) rather than (0, T), because no
lattice with step delta ever evaluates K closer to the origin than delta, and on
(0, T) the singularity that no finite sum of exponentials can follow dominates
the norm and the criterion says nothing about the scheme.

What is measured and what is proved
-----------------------------------
Measured (see `runs/nm-bound-*`): m*(n) at fixed tolerance grows like log n, for
every H tested and for both constructions, and the two constructions differ in
how the tolerance enters.  NOT proved: the inequality itself.  `m_required`
below is a fit with a certification fallback, not a theorem.  The open side is
the lower bound — that no m-factor lift with fewer factors can do it — and then
propagating a kernel-error bound to a PRICE error through Proposition (B1'),
whose constant is still unevaluated.  Do not cite `m_required` as a bound.
"""
from __future__ import annotations

import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import route_b as rb                                            # noqa: E402


M_MAX = 14          # search ceiling; every table here stays well inside it


# ----------------------------------------------------------------- the criterion
def rel_l2(w: np.ndarray, s: np.ndarray, H: float, T: float, t0: float) -> float:
    """|| K - K^m ||_{L^2(t0,T)} / || K ||_{L^2(t0,T)}, exact, scale-free."""
    h = H - 0.5
    num = rb.l2_error(w, s, h, T, t0)
    den = math.sqrt(rb.norm_K_sq(h, T, t0))
    return num / den


# ------------------------------------------------- the explicit construction
def geometric_lift(m: int, H: float, n: int, T: float = 1.0,
                   c_lo: float = 1.0, c_hi: float = 1.0
                   ) -> tuple[np.ndarray, np.ndarray, float]:
    """AJEE cell masses/barycentres on a GEOMETRIC partition of the s-range.

    The partition's upper edges are

        eta_i = eta_max r^{i-m},   i = 1..m,
        eta_max = c_hi / delta  (the stiffness the finest lag needs),
        eta_1   = c_lo / T      (the softest scale the horizon needs),

    so the ratio between consecutive nodes is fixed at

        r = (eta_max / eta_1)^{1/(m-1)} = (c_hi n / c_lo)^{1/(m-1)}.

    Inverting gives the shape of the whole question:

        m = 1 + log(c_hi n / c_lo) / log(r),

    i.e. m is logarithmic in n at fixed per-node ratio r.  Whether a fixed r
    buys a fixed accuracy — the step that makes the log law an inequality — is
    what phase `ratio` of `run_nm_bound.py` measures.

    Weights are positive by construction (each is a genuine OU variance
    loading).  Returns (w, s, r), with w normalised as K = u^h/Gamma(1+h).
    """
    if m < 1:
        raise ValueError("m >= 1")
    h = H - 0.5
    delta = T / n
    eta_max = c_hi / delta
    if m == 1:
        eta = np.array([eta_max])
        r = float("inf")
    else:
        eta_1 = c_lo / T
        r = (eta_max / eta_1) ** (1.0 / (m - 1))
        eta = eta_max * r ** (np.arange(1, m + 1, dtype=float) - m)
    w, s = rb.ajee_from_partition(eta, h)
    return w, s, float(r)


def geometric_lift_by_ratio(r: float, H: float, n: int, T: float = 1.0,
                            c_hi: float = 1.0
                            ) -> tuple[np.ndarray, np.ndarray, int]:
    """The same partition specified by the RATIO instead of the count.

    m = 1 + ceil(log(c_hi n) / log r), so the range [1/T, c_hi/delta] is covered
    with steps of at most r.  This is the form the bound is stated in: fix the
    accuracy per node, read off how many nodes the range costs.
    """
    if r <= 1.0:
        raise ValueError("r > 1")
    delta = T / n
    span = (c_hi / delta) / (1.0 / T)
    m = 1 + int(math.ceil(math.log(span) / math.log(r)))
    w, s, _ = geometric_lift(m, H, n, T, c_lo=1.0, c_hi=c_hi)
    return w, s, m


# --------------------------------------------------------------- error curves
def error_curve(n: int, H: float, T: float = 1.0, m_max: int = M_MAX,
                construction: str = "best") -> list[float]:
    """rel(m) for m = 1..m_max, on the lags the scheme sees (t0 = delta).

    `construction` is one of

      "best"       route_b.best_lift  — nodes free, NNLS weights: the achievable
                   floor, warm-started from m-1 so the curve is monotone.  This
                   is the branch the Route B lattice actually uses (`lift_for`).
      "ajee"       route_b.ajee_optimised — the AJEE two-parameter partition
                   family, optimised, warm-started likewise.
      "geometric"  the explicit `geometric_lift` — no optimisation at all.
    """
    h = H - 0.5
    delta = T / n
    out: list[float] = []
    if construction == "geometric":
        for m in range(1, m_max + 1):
            w, s, _ = geometric_lift(m, H, n, T)
            out.append(rel_l2(w, s, H, T, delta))
        return out
    if construction == "best":
        prev = None
        for m in range(1, m_max + 1):
            s0 = np.geomspace(1.0, 50.0, m) if m > 1 else np.array([2.0])
            b = rb.best_lift(m, h, T, s0, t0=delta,
                             extra_inits=[prev] if prev is not None else None)
            prev = b["s"]
            out.append(rel_l2(b["w"], b["s"], H, T, delta))
        return out
    if construction == "ajee":
        seeds: list = []
        for m in range(1, m_max + 1):
            b = rb.ajee_optimised(m, h, T, t0=delta, seeds=seeds)
            seeds = [(b["family"], b["eta_max"], b["shape"])]
            out.append(rel_l2(b["w"], b["s"], H, T, delta))
        return out
    raise ValueError(f"unknown construction {construction!r}")


def m_star_from_curve(curve: list[float], tol: float) -> int | None:
    """Smallest m whose relative error is at or below `tol`; None if none is.

    The curve is not assumed monotone — the optimisers are local searches and can
    wobble — so this is the first crossing, and a later m that fails does not
    undo it.  A wobble large enough to matter is reported by the runner.
    """
    for i, e in enumerate(curve):
        if e <= tol:
            return i + 1
    return None


def m_star(n: int, H: float, tol: float, T: float = 1.0,
           construction: str = "best", m_max: int = M_MAX) -> int | None:
    """Measured m*(n, H, tol): the smallest number of factors meeting `tol`."""
    return m_star_from_curve(error_curve(n, H, T, m_max, construction), tol)


# ------------------------------------------------------------------- the fits
def fit_log_law(ns, ms) -> dict:
    """Least squares of m* against ln n: m* ~ a ln n + b, with residuals.

    Reported with the residuals and the R^2 so the log law is TESTED rather than
    assumed; m* is integer-valued, so a fit through a handful of points can look
    good for the wrong reason and the residual spread is the honest read.
    """
    x = np.log(np.asarray(ns, float))
    y = np.asarray(ms, float)
    if len(x) < 3:
        return {"slope": float("nan"), "intercept": float("nan"),
                "r2": float("nan"), "max_resid": float("nan"), "npts": len(x)}
    a, b = np.polyfit(x, y, 1)
    pred = a * x + b
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    return {"slope": float(a), "intercept": float(b),
            "r2": (1.0 - ss_res / ss_tot) if ss_tot > 0 else float("nan"),
            "max_resid": float(np.max(np.abs(y - pred))), "npts": len(x)}


def fit_power_law(xs, ys) -> dict:
    """Least squares of ln y against ln x: y ~ C x^p.  Used to ask how the
    tolerance enters the slope — a power of 1/eps or a log of it."""
    lx = np.log(np.asarray(xs, float))
    ly = np.log(np.asarray(ys, float))
    if len(lx) < 2:
        return {"exponent": float("nan"), "const": float("nan"), "npts": len(lx)}
    p, c = np.polyfit(lx, ly, 1)
    pred = p * lx + c
    ss_res = float(np.sum((ly - pred) ** 2))
    ss_tot = float(np.sum((ly - ly.mean()) ** 2))
    return {"exponent": float(p), "const": float(math.exp(c)),
            "r2": (1.0 - ss_res / ss_tot) if ss_tot > 0 else float("nan"),
            "npts": len(lx)}


# ------------------------------------------------- the candidate law, for use
#
# Calibrated by `runs/nm-bound-20260805T*` — see ALGORITHMS.md and
# docs/N_M_BOUND.md.  The shape is the one the mechanism above predicts,
#
#     m*(n, H, eps)  ~=  A(H) * ln n * (1/eps)^p(H)  +  B(H),
#
# fitted per H on the "best" (achievable-floor) construction, which is the
# branch the Route B lattice uses.  These are FITTED CONSTANTS, not a theorem:
# `m_required` therefore verifies its answer by computation before returning it.
_LAW: dict[float, dict] = {}          # filled from the run; see set_law()
_LAW_SOURCE = "uncalibrated"


def set_law(law: dict, source: str) -> None:
    """Install fitted constants (used by the runner, and by the frozen table)."""
    global _LAW, _LAW_SOURCE
    _LAW, _LAW_SOURCE = dict(law), source


def _law_for(H: float) -> dict | None:
    if not _LAW:
        return None
    key = min(_LAW, key=lambda k: abs(k - H))
    return _LAW[key] if abs(key - H) <= 0.051 else None


def m_guess(n: int, H: float, tol: float) -> int:
    """The fitted law's guess at m*, with no verification.  Cheap, may be wrong."""
    law = _law_for(H)
    if law is None:                       # uncalibrated: the mechanism's shape
        return max(1, int(math.ceil(0.5 * math.log(n) * (0.03 / tol) ** 0.5)))
    m = law["A"] * math.log(n) * (1.0 / tol) ** law["p"] + law["B"]
    return max(1, int(math.ceil(m)))


def m_required(n: int, H: float, tol: float = 0.03, T: float = 1.0,
               construction: str = "best", verify: bool = True,
               m_max: int = M_MAX) -> dict:
    """How many lift factors a run at `n` steps needs — the usable entry point.

    Starts from the fitted law, then, if `verify`, walks up or down until the
    returned m is the smallest one that actually meets `tol` under the chosen
    construction.  So the answer is never wrong, only sometimes conservative,
    and the fit is a starting point rather than a claim.

    Returns the m, its measured relative error, and whether the fit was right,
    so a runner can record that it chose m from the bound and not by hand.
    """
    guess = min(max(1, m_guess(n, H, tol)), m_max)
    if not verify:
        return {"m": guess, "rel_error": None, "verified": False,
                "guess": guess, "source": _LAW_SOURCE}
    curve = error_curve(n, H, T, m_max, construction)
    m = m_star_from_curve(curve, tol)
    if m is None:
        return {"m": m_max, "rel_error": curve[-1], "verified": True,
                "met": False, "guess": guess, "source": _LAW_SOURCE,
                "construction": construction}
    return {"m": m, "rel_error": curve[m - 1], "verified": True, "met": True,
            "guess": guess, "guess_ok": guess == m, "source": _LAW_SOURCE,
            "construction": construction}


# --------------------------------------------------- the other criterion: variance
def variance_ratio_criterion(n: int, H: float, m: int, T: float = 1.0,
                             construction: str = "best") -> dict:
    """The DRIVER VARIANCE criterion at n steps, with the context to read it.

    The kernel L^2 error is a proxy; what actually drives the price error is the
    driver's variance, and Proposition 8.3 is stated in exactly that quantity.
    `route_b.discrete_covariance_report` computes it exactly.  Returned:

      lift_vs_true    |Var[V^m_T]/Var[V_T] - 1|, both discrete: the error the
                      lift itself adds.  This is the criterion.
      true_vs_cont    the same for the true kernel discretised versus the
                      continuous covariance: the error the SCHEME already has,
                      without any lift.  The criterion is meaningless without it
                      — there is no point driving the lift below it.
      rel_l2_delta    the L^2(delta,T) relative kernel error of the same lift.
      rel_l2_zero     the same on L^2(0,T), i.e. INCLUDING the singular region.

    The last two are what makes the comparison diagnostic.  A lattice with step
    delta never evaluates K at a lag below delta, which is why L^2(delta,T) is
    the natural criterion — but its first cell average runs over (0, delta] and
    therefore does straddle the singularity.  If m* under the variance criterion
    tracks `rel_l2_zero` rather than `rel_l2_delta`, that first cell is the
    reason, and an L^2(delta,T) tolerance alone cannot control the variance.
    """
    h = H - 0.5
    delta = T / n
    if construction == "geometric":
        w, s, _ = geometric_lift(m, H, n, T)
    elif construction == "best":
        prev, w, s = None, None, None
        for mm in range(1, m + 1):
            s0 = np.geomspace(1.0, 50.0, mm) if mm > 1 else np.array([2.0])
            b = rb.best_lift(mm, h, T, s0, t0=delta,
                             extra_inits=[prev] if prev is not None else None)
            prev, w, s = b["s"], b["w"], b["s"]
    else:
        raise ValueError(f"unknown construction {construction!r}")
    rep = rb.discrete_covariance_report(w, s, h, T, n, mode="cellavg")
    return {"lift_vs_true":
            abs(rep["lift_vs_true_discrete"]["var_T_ratio"] - 1.0),
            "true_vs_cont":
            abs(rep["true_discrete_vs_continuous"]["var_T_ratio"] - 1.0),
            "rel_l2_delta": rel_l2(w, s, H, T, delta),
            "rel_l2_zero": rel_l2(w, s, H, T, 0.0)}


def onestep_variance_ratio(n: int, H: float) -> float:
    """For contrast: the one-step scheme's ratio, 2H n^{1-2H} — it diverges."""
    return rb.onestep_variance_ratio(n, H)
