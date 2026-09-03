#!/usr/bin/env python3
"""
run_route_b_scaling — preliminary 2: does the number of factors m stay bounded
as the time grid is refined?

The first Route B run (runs/route-b-*) answered "how many factors" at n = 64 and
n = 256 and got small numbers.  But the threshold the lift has to beat is the
scheme's OWN discretisation error, and that shrinks as n grows, so m must grow
too.  If it grows fast, Route B is only cheap at grids too coarse to matter.

Preliminary 1 settled the convention: the paper's discrete kernel is the LEFT
ENDPOINT with smallest lag delta (Definition "Exact convolution", lags
t_k - t_{i-1} = (k-i+1) delta), so the norm that governs the lift is
L^2(delta,T) and the lifts are optimised for it.

Criterion.  Everything is closed form, so n can go far beyond what a covariance
matrix allows:

    baseline(n)   = |Var_disc_true(n) - Var_cont| / Var_cont
                    the error the scheme already commits with the TRUE kernel
    lift(n,m)     = |Var_disc_lift(n,m) - Var_disc_true(n)| / Var_disc_true(n)
                    the error the lift adds on top
    m*(n)         = smallest m with lift(n,m) < baseline(n)

with Var_disc_G(n) = delta * sum_{j=1}^{n} G(j delta)^2 and
Var_cont = ||K||^2_{L^2(0,T)}.

A confirmatory pass repeats the crossing with the full covariance matrix and the
relative Frobenius norm, at the n where that matrix is affordable, to check that
the terminal variance is not giving a flattering answer.

Four phases; every number lands in RESULTS.md and is registered in FINDINGS.md.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from progress import Progress            # noqa: E402
import route_b as rb                     # noqa: E402

H = 0.10
KAPPA = 0.50
h = H - 0.5
T = 1.0
NS = (64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384)
MS = (1, 2, 3, 4, 5, 6, 8, 10, 12)
NS_FROB = (64, 128, 256, 512)             # where the n x n matrix is affordable
VAR_CONT = rb.norm_K_sq(h, T)


def var_disc(g: np.ndarray, d: float) -> float:
    """delta * sum_j G(j delta)^2 — the terminal variance the lattice carries."""
    return float(d * (g ** 2).sum())


def var_disc_true(n: int) -> float:
    d = T / n
    return var_disc(rb.K_true(np.arange(1, n + 1) * d, h), d)


def phase1(pr: Progress) -> dict:
    pr.phase("the threshold: the scheme's own error with the TRUE kernel",
             total=len(NS))
    rows = []
    for n in NS:
        vt = var_disc_true(n)
        rows.append({"n": n, "delta": T / n,
                     "var_disc_true": vt, "var_continuous": VAR_CONT,
                     "var_ratio": vt / VAR_CONT,
                     "baseline_rel_error": abs(vt - VAR_CONT) / VAR_CONT,
                     "lattice_rate_delta_pow": rb.lattice_error(n, h, KAPPA, T)})
        pr.tick(note=f"n={n}")
        pr.partial.update({f"baseline_n{n}": round(rows[-1]["baseline_rel_error"], 5)})
    ns = np.array([r["n"] for r in rows], float)
    bl = np.array([r["baseline_rel_error"] for r in rows], float)
    slope = float(np.polyfit(np.log(ns), np.log(bl), 1)[0])
    return {"rows": rows, "baseline_decay_slope_in_n": round(slope, 4),
            "predicted_slope_minus_2H": -2.0 * H}


def phase2(pr: Progress, r1: dict) -> dict:
    pr.phase("m*(n): smallest m whose added error is below the threshold",
             total=len(NS) * len(MS))
    base = {r["n"]: r["baseline_rel_error"] for r in r1["rows"]}
    out = []
    for n in NS:
        d = T / n
        vt = var_disc_true(n)
        lags = np.arange(1, n + 1) * d
        seed, prev = [], None
        errs_a, errs_f, l2a, l2f = {}, {}, {}, {}
        for m in MS:
            b = rb.ajee_optimised(m, h, T, t0=d, seeds=seed)
            seed = [(b["family"], b["eta_max"], b["shape"])]
            f = rb.best_lift(m, h, T, b["s"], t0=d,
                             extra_inits=[prev] if prev is not None else None)
            prev = f["s"]
            for lift, errs, l2s in ((b, errs_a, l2a), (f, errs_f, l2f)):
                v = var_disc(rb.K_lift(lags, lift["w"], lift["s"]), d)
                errs[m] = abs(v - vt) / vt
                l2s[m] = lift["l2"]
            pr.tick(note=f"n={n} m={m}")
            pr.partial.update({f"n{n}_m{m}_floor": float(f"{errs_f[m]:.3g}")})

        def cross(errs):
            for m in MS:
                if errs[m] < base[n]:
                    return m
            return None

        out.append({"n": n, "baseline": base[n],
                    "m_star_ajee": cross(errs_a), "m_star_floor": cross(errs_f),
                    "errors_ajee": {m: errs_a[m] for m in MS},
                    "errors_floor": {m: errs_f[m] for m in MS},
                    "l2_ajee": {m: l2a[m] for m in MS},
                    "l2_floor": {m: l2f[m] for m in MS}})
        pr.result(f"m_star_n{n}", {"ajee": out[-1]["m_star_ajee"],
                                   "floor": out[-1]["m_star_floor"]})
    return {"rows": out}


def phase3(pr: Progress) -> dict:
    pr.phase("confirmation with the full covariance matrix (Frobenius)",
             total=len(NS_FROB) * len(MS))
    out = []
    for n in NS_FROB:
        d = T / n
        seed, prev = [], None
        rows, base = [], None
        for m in MS:
            b = rb.ajee_optimised(m, h, T, t0=d, seeds=seed)
            seed = [(b["family"], b["eta_max"], b["shape"])]
            f = rb.best_lift(m, h, T, b["s"], t0=d,
                             extra_inits=[prev] if prev is not None else None)
            prev = f["s"]
            rep = rb.discrete_covariance_report(f["w"], f["s"], h, T, n,
                                                mode="left")
            base = rep["true_discrete_vs_continuous"]["rel_frobenius"]
            rows.append({"m": m,
                         "frob_floor": rep["lift_vs_true_discrete"]["rel_frobenius"]})
            pr.tick(note=f"n={n} m={m}")
        mstar = next((r["m"] for r in rows if r["frob_floor"] < base), None)
        out.append({"n": n, "baseline_frobenius": base, "m_star_floor": mstar,
                    "rows": rows})
        pr.result(f"frob_m_star_n{n}", mstar)
    return {"rows": out}


def phase4(pr: Progress, r1: dict, r2: dict) -> dict:
    pr.phase("growth of m*(n) and what it costs", total=3)
    pts = [(r["n"], r["m_star_floor"]) for r in r2["rows"]
           if r["m_star_floor"] is not None]
    pr.tick()
    fit = None
    if len(pts) >= 3:
        ns = np.array([p[0] for p in pts], float)
        ms = np.array([p[1] for p in pts], float)
        a, b = np.polyfit(np.log(ns), ms, 1)          # m* ~ a log n + b
        fit = {"form": "m*(n) ~ a*log(n) + b", "a": round(float(a), 4),
               "b": round(float(b), 4),
               "m_star_at_n_1e6": round(float(a * math.log(1e6) + b), 1),
               "m_star_at_n_1e20": round(float(a * math.log(1e20) + b), 1)}
    pr.tick()
    # n needed for a target accuracy eps at the rate delta^{H/2}
    tgt = {}
    for eps in (1e-1, 1e-2, 1e-3):
        n_needed = eps ** (-2.0 / H)
        m_needed = (None if fit is None
                    else max(1.0, fit["a"] * math.log(n_needed) + fit["b"]))
        tgt[eps] = {"n_needed": f"{n_needed:.3g}",
                    "m_needed": None if m_needed is None else round(m_needed, 1),
                    "cost_exponent_q": None if m_needed is None else
                    round(rb.cost_exponent(int(math.ceil(m_needed)), H, H / 2), 1)}
    pr.tick()
    return {"points": pts, "fit": fit, "targets": tgt,
            "route_aprime_alone_cost_exponent": round((2.0 + H / 2) * 2.0 / H, 1)}


def build_results(pr, r1, r2, r3, r4) -> str:
    L = ["# Route B, preliminary 2 — does m stay bounded as the grid refines?", "",
         f"Left-endpoint discrete kernel with smallest lag delta (the paper's "
         f"Definition \"Exact convolution\"),", f"so the lifts are optimised for "
         f"L2(delta,T).  h = {h}, H = {H}, T = {T}, kappa = {KAPPA}.",
         f"Var[V_T] continuous = {VAR_CONT:.6f}.", "",
         "## Phase 1 — the threshold to beat", "",
         "This is the error the scheme already commits with the TRUE kernel, "
         "purely from cutting", "time into n pieces.  A lift whose added error "
         "sits below it is not the binding constraint.", "",
         "| n | delta | Var_disc(true K) | Var ratio (target 1) | threshold "
         "= relative error | delta^{H/2} for comparison |",
         "|---|---|---|---|---|---|"]
    for r in r1["rows"]:
        L.append(f"| {r['n']} | {r['delta']:.3e} | {r['var_disc_true']:.6f} "
                 f"| {r['var_ratio']:.6f} | **{r['baseline_rel_error']:.4e}** "
                 f"| {r['lattice_rate_delta_pow']:.4f} |")
    L += ["", f"The threshold decays in n with log-log slope "
              f"**{r1['baseline_decay_slope_in_n']}**, against the predicted "
              f"-2H = {r1['predicted_slope_minus_2H']}.", "",
          "## Phase 2 — m*(n), the smallest usable number of factors", "",
          "| n | threshold | m* (Abi Jaber–El Euch) | m* (achievable floor) |",
          "|---|---|---|---|"]
    for r in r2["rows"]:
        L.append(f"| {r['n']} | {r['baseline']:.4e} | {r['m_star_ajee']} "
                 f"| **{r['m_star_floor']}** |")
    L += ["", "### The full error table (achievable floor), added variance error",
          "", "| n | " + " | ".join(f"m={m}" for m in MS) + " |",
          "|---|" + "---|" * len(MS)]
    for r in r2["rows"]:
        L.append(f"| {r['n']} | "
                 + " | ".join(f"{r['errors_floor'][m]:.2e}" for m in MS) + " |")
    L += ["", "### And the kernel L2(delta,T) error behind it", "",
          "| n | " + " | ".join(f"m={m}" for m in MS) + " |",
          "|---|" + "---|" * len(MS)]
    for r in r2["rows"]:
        L.append(f"| {r['n']} | "
                 + " | ".join(f"{r['l2_floor'][m]:.2e}" for m in MS) + " |")
    L += ["", "## Phase 3 — confirmation with the whole covariance matrix", "",
          "The terminal variance is one number and could flatter the lift; this "
          "repeats the crossing", "with the relative Frobenius norm of the full "
          "n x n covariance matrix.", "",
          "| n | threshold (Frobenius) | m* (floor) |", "|---|---|---|"]
    for r in r3["rows"]:
        L.append(f"| {r['n']} | {r['baseline_frobenius']:.4e} "
                 f"| **{r['m_star_floor']}** |")
    L += ["", "## Phase 4 — growth of m* and what it costs", ""]
    if r4["fit"]:
        f = r4["fit"]
        L += [f"Fit: `{f['form']}` with a = **{f['a']}**, b = {f['b']}.",
              f"Extrapolated: m* = {f['m_star_at_n_1e6']} at n = 10^6, "
              f"m* = {f['m_star_at_n_1e20']} at n = 10^20.", "",
              "A logarithmic fit is the claim being tested: if m* grows like "
              "log n, Route B is usable", "at every n; if it grows like a power "
              "of n, it is not.", ""]
    else:
        L += ["Not enough crossings to fit a growth law.", ""]
    L += ["| target accuracy eps | n needed (rate delta^{H/2}) | m needed "
          "| cost eps^-q with Route A' inside |", "|---|---|---|---|"]
    for eps, t in r4["targets"].items():
        L.append(f"| {eps:g} | {t['n_needed']} | {t['m_needed']} "
                 f"| {t['cost_exponent_q']} |")
    L += ["", f"Route A' alone (inconsistent, so not a competitor): "
              f"eps^-{r4['route_aprime_alone_cost_exponent']}.", "",
          pr.timing_table_md()]
    return "\n".join(L)


def main() -> None:
    meta = {"item": "preliminary 2 — is m* bounded as n grows?",
            "convention": "left-endpoint discrete kernel, smallest lag delta "
                          "(paper Dfn 'Exact convolution')",
            "norm_for_the_lift": "L2(delta,T)",
            "H": H, "h": h, "kappa": KAPPA, "T": T,
            "n_values": list(NS), "m_values": list(MS),
            "n_values_frobenius": list(NS_FROB),
            "var_continuous": round(VAR_CONT, 6)}
    with Progress("route-b-scaling", total_phases=4, meta=meta) as pr:
        r1 = phase1(pr); pr.result("phase1_threshold", r1)
        r2 = phase2(pr, r1); pr.result("phase2_m_star", r2)
        r3 = phase3(pr); pr.result("phase3_frobenius_confirmation", r3)
        r4 = phase4(pr, r1, r2); pr.result("phase4_growth_and_cost", r4)
        pr.write_results_md(build_results(pr, r1, r2, r3, r4))


if __name__ == "__main__":
    main()
