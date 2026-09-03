#!/usr/bin/env python3
"""
run_route_b — item (N4)/(B1): quantify the Markovian lift BEFORE proving it.

The question that decides whether Route B is usable: how many factors m does the
lift need, and what does that cost?

Seven phases.

  1  AJEE lift for m = 1,2,3,5,8,12,20,30 in the norm ||.||_{L^2(0,T)} that (B1)
     of the paper is stated in: exact error, an independent quadrature check that
     resolves every boundary layer, the empirical decay rate in m, and how much
     of the error lives on (0, delta) — a neighbourhood of the origin that a
     lattice of time step delta never visits.
  2  The achievable floor in the same norm: nodes free, weights non-negative
     least-squares.  No verdict on m may rest on a poor choice of nodes.
  3  The norm ||.||_{L^2(delta,T)}, with the partition re-optimised for it.
  4  The crossing: smallest m for which the kernel error falls below the lattice
     error delta^{(h+kappa)/2}, under three normalisations and both norms.  The
     constant in (B1) is unspecified and an L^2 kernel norm is not dimensionless,
     so a single number here would be dishonest.
  5  The CONTINUOUS covariance surface: relative error of Cov(V_u, V_v).
  6  The DISCRETE covariance, which is the object the scheme's consistency
     actually depends on, measured against the quadrature error the scheme
     already carries.  This is where the two competing norms are adjudicated.
  7  Cost: eps^{-2(m+1)/H}, against the eps^{-41} of Route A' alone.

Every number lands in runs/<name>/RESULTS.md and is registered in FINDINGS.md.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from scipy.special import gamma as Gamma

sys.path.insert(0, str(Path(__file__).resolve().parent))
from progress import Progress            # noqa: E402
import route_b as rb                     # noqa: E402

H = 0.10
KAPPA = 0.50                             # so (h+kappa)/2 = H/2, the rate of Thm 7.1
h = H - 0.5
T = 1.0
MS = (1, 2, 3, 5, 8, 12, 20, 30)
NS = (64, 256)
NORM_K = math.sqrt(rb.norm_K_sq(h, T))
G1H = float(Gamma(1.0 + h))


def _slope(rows: list, key: str) -> float:
    ms = np.array([r["m"] for r in rows], float)
    er = np.array([r[key] for r in rows], float)
    ok = er > 0
    if ok.sum() < 2:
        return float("nan")
    return float(np.polyfit(np.log(ms[ok]), np.log(er[ok]), 1)[0])


def phase1(pr: Progress) -> dict:
    pr.phase("AJEE lift in L2(0,T) — the norm (B1) is stated in", total=len(MS))
    rows, store, seed = [], {}, []
    for m in MS:
        b = rb.ajee_optimised(m, h, T, seeds=seed)
        seed = [(b["family"], b["eta_max"], b["shape"])]
        rows.append({
            "m": m, "family": b["family"],
            "eta_max": float(f"{b['eta_max']:.4g}"),
            "shape": round(float(b["shape"]), 4),
            "l2_abs": b["l2"],
            "l2_quad_check": rb.l2_error_quadrature(b["w"], b["s"], h, T),
            "l2_rel": b["l2"] / NORM_K,
            "share_below_delta_n64": rb.singularity_share(b["w"], b["s"], h, T, T / 64),
            "share_below_delta_n256": rb.singularity_share(b["w"], b["s"], h, T, T / 256),
            "s_min": float(f"{b['s'].min():.4g}"),
            "s_max": float(f"{b['s'].max():.4g}")})
        store[m] = b
        pr.tick(note=f"m={m}")
        pr.partial.update({f"ajee_l2_m{m}": float(f"{b['l2']:.4g}")})
    return {"norm_K": NORM_K, "rows": rows,
            "decay_slope_in_m": round(_slope(rows, "l2_abs"), 3),
            "_store": store}


def phase2(pr: Progress, store: dict) -> dict:
    pr.phase("achievable floor in L2(0,T) — free nodes, NNLS weights",
             total=len(MS))
    rows, store2, prev = [], {}, None
    for m in MS:
        b = rb.best_lift(m, h, T, store[m]["s"],
                         extra_inits=[prev] if prev is not None else None)
        prev = b["s"]
        rows.append({"m": m, "l2_abs": b["l2"], "l2_rel": b["l2"] / NORM_K,
                     "active_factors": b["active_factors"],
                     "gain_vs_ajee": round(store[m]["l2"] / b["l2"], 2)
                     if b["l2"] > 0 else None})
        store2[m] = b
        pr.tick(note=f"m={m}")
        pr.partial.update({f"floor_l2_m{m}": float(f"{b['l2']:.4g}")})
    return {"rows": rows, "decay_slope_in_m": round(_slope(rows, "l2_abs"), 3),
            "_store": store2}


def phase3(pr: Progress) -> dict:
    pr.phase("the norm the discrete kernel sees: L2(delta,T)",
             total=len(NS) * len(MS))
    out, store = {}, {}
    for n in NS:
        d = T / n
        nk = math.sqrt(rb.norm_K_sq(h, T, d))
        rows, st, seed, prev = [], {}, [], None
        for m in MS:
            b = rb.ajee_optimised(m, h, T, t0=d, seeds=seed)
            seed = [(b["family"], b["eta_max"], b["shape"])]
            f = rb.best_lift(m, h, T, b["s"], t0=d,
                             extra_inits=[prev] if prev is not None else None)
            prev = f["s"]
            rows.append({"m": m,
                         "l2_abs_ajee": b["l2"],
                         "l2_quad_check": rb.l2_error_quadrature(
                             b["w"], b["s"], h, T, d),
                         "l2_rel_ajee": b["l2"] / nk,
                         "l2_abs_floor": f["l2"],
                         "l2_rel_floor": f["l2"] / nk})
            st[m] = f
            pr.tick(note=f"n={n} m={m}")
            pr.partial.update({f"delta{n}_l2_m{m}": float(f"{b['l2']:.4g}")})
        out[n] = {"delta": d, "norm_K_on_interval": nk, "rows": rows,
                  "decay_slope_ajee": round(_slope(rows, "l2_abs_ajee"), 3),
                  "decay_slope_floor": round(_slope(rows, "l2_abs_floor"), 3)}
        store[n] = st
    return {"by_n": out, "_store": store}


def _first_m(rows: list, key, thresh: float):
    for r in rows:
        v = key(r)
        if v is not None and v < thresh:
            return r["m"]
    return None


def phase4(pr: Progress, r1: dict, r2: dict, r3: dict) -> dict:
    pr.phase("the crossing: kernel error vs lattice error", total=len(NS))
    out = []
    for n in NS:
        lat = rb.lattice_error(n, h, KAPPA, T)
        d3 = r3["by_n"][n]
        rows = {
            "L2(0,T), absolute, K = u^h/Gamma(1+h)": {
                "ajee": _first_m(r1["rows"], lambda r: r["l2_abs"], lat),
                "floor": _first_m(r2["rows"], lambda r: r["l2_abs"], lat)},
            "L2(0,T), absolute, K = u^h": {
                "ajee": _first_m(r1["rows"], lambda r: r["l2_abs"] * G1H, lat),
                "floor": _first_m(r2["rows"], lambda r: r["l2_abs"] * G1H, lat)},
            "L2(0,T), relative to ||K||": {
                "ajee": _first_m(r1["rows"], lambda r: r["l2_rel"], lat),
                "floor": _first_m(r2["rows"], lambda r: r["l2_rel"], lat)},
            "L2(delta,T), absolute, K = u^h/Gamma(1+h)": {
                "ajee": _first_m(d3["rows"], lambda r: r["l2_abs_ajee"], lat),
                "floor": _first_m(d3["rows"], lambda r: r["l2_abs_floor"], lat)},
            "L2(delta,T), absolute, K = u^h": {
                "ajee": _first_m(d3["rows"], lambda r: r["l2_abs_ajee"] * G1H, lat),
                "floor": _first_m(d3["rows"], lambda r: r["l2_abs_floor"] * G1H, lat)},
            "L2(delta,T), relative to ||K||": {
                "ajee": _first_m(d3["rows"], lambda r: r["l2_rel_ajee"], lat),
                "floor": _first_m(d3["rows"], lambda r: r["l2_rel_floor"], lat)},
        }
        out.append({"n": n, "delta": T / n, "lattice_error": lat,
                    "crossings": rows})
        pr.tick(note=f"n={n}")
        pr.partial.update({f"lattice_err_n{n}": round(lat, 5)})
    tiers = {}
    for tgt in (1e-1, 1e-2, 1e-3):
        tiers[tgt] = {
            "m_ajee_L2_0T": _first_m(r1["rows"], lambda r: r["l2_rel"], tgt),
            "m_floor_L2_0T": _first_m(r2["rows"], lambda r: r["l2_rel"], tgt),
            "m_floor_L2_deltaT_n256": _first_m(
                r3["by_n"][256]["rows"], lambda r: r["l2_rel_floor"], tgt)}
    return {"rows": out, "relative_accuracy_tiers": tiers,
            "gamma_1_plus_h": round(G1H, 6),
            "lattice_rate_note": "lattice error = delta^{(h+kappa)/2} = "
                                 f"delta^{0.5*(h+KAPPA):.2f}"}


def phase5(pr: Progress, store0: dict) -> dict:
    pr.phase("the continuous covariance surface", total=len(MS))
    rows = []
    for m in MS:
        b = store0[m]
        rep = rb.covariance_report(b["w"], b["s"], h, T, grid=12)
        rows.append({"m": m, "max_rel_error": rep["max_rel_error"],
                     "max_rel_error_at": rep["max_rel_error_at"],
                     "rel_frobenius": rep["rel_frobenius_error"],
                     "var_T_ratio": rep["var_T_ratio"]})
        pr.tick(note=f"m={m}")
        pr.partial.update({f"cov_m{m}": float(f"{rep['max_rel_error']:.4g}")})
    return {"rows": rows,
            "lift_optimised_for": "L2(0,T)",
            "onestep_lattice_var_ratio": {
                n: round(rb.onestep_variance_ratio(n, H), 2)
                for n in (8, 16, 32, 64, 256, 4096)},
            "note": "target for every ratio is 1.0; the one-step lattice ratio "
                    "diverges like 2H n^{1-2H} (Prop 8.3)"}


MODES = (("cellavg", "exact convolution (cell average of K)"),
         ("left", "naive left endpoint K(l*delta)"))


def phase6(pr: Progress, store0: dict, store_d: dict) -> dict:
    pr.phase("the DISCRETE covariance — adjudicating the two norms",
             total=len(MODES) * 2 * len(NS) * len(MS))
    out = {}
    for n in NS:
        per_mode = {}
        for mode, mlabel in MODES:
            base, rows = None, []
            for label, st in (("L2(0,T)", store0), ("L2(delta,T)", store_d[n])):
                for m in MS:
                    b = st[m]
                    rep = rb.discrete_covariance_report(b["w"], b["s"], h, T,
                                                        n, mode=mode)
                    base = rep["true_discrete_vs_continuous"]
                    lv, lc = (rep["lift_vs_true_discrete"],
                              rep["lift_discrete_vs_continuous"])
                    rows.append({"lift": label, "m": m,
                                 "max_rel": lv["max_rel_error"],
                                 "frob": lv["rel_frobenius"],
                                 "varT": lv["var_T_ratio"],
                                 "tot_frob": lc["rel_frobenius"],
                                 "tot_varT": lc["var_T_ratio"]})
                    pr.tick(note=f"n={n} {mode} {label} m={m}")
                    pr.partial.update(
                        {f"disc_{mode}_n{n}_m{m}": float(f"{lv['rel_frobenius']:.4g}")})
            per_mode[mode] = {"label": mlabel, "rows": rows,
                              "scheme_own_discretisation_error": base}
        out[n] = {"delta": T / n, "by_mode": per_mode}
    return {"by_n": out, "modes": [m for m, _ in MODES]}


def phase7(pr: Progress) -> dict:
    pr.phase("cost in the target accuracy", total=len(MS) + 1)
    gamma = H / 2.0
    aprime = (2.0 + gamma) * 2.0 / H
    pr.result("route_aprime_alone_cost_exponent", round(aprime, 1))
    pr.tick()
    rows = []
    for m in MS:
        rows.append({"m": m, "states": f"O(n^{m+1})",
                     "cost_exponent_q": round(rb.cost_exponent(m, H), 1),
                     "with_aprime": round(rb.cost_exponent(m, H, gamma), 1)})
        pr.tick(note=f"m={m}")
    return {"n_in_eps": "n ~ eps^{-2/H} = eps^{-20} at H = 0.1",
            "route_aprime_alone": round(aprime, 1), "rows": rows}


def build_results(pr, r1, r2, r3, r4, r5, r6, r7) -> str:
    L = ["# Route B quantified — how many factors does the lift need?", "",
         f"Riemann--Liouville kernel K(u) = u^h/Gamma(1+h) with h = {h}, "
         f"H = {H}, T = {T}, kappa = {KAPPA}.",
         f"Lattice rate exponent (h+kappa)/2 = {0.5*(h+KAPPA):.2f}; "
         f"||K||_L2(0,T) = {NORM_K:.6f}; Gamma(1+h) = {G1H:.6f}.", "",
         "## Phase 1 — AJEE lift in L2(0,T), the norm (B1) is stated in", "",
         "| m | partition | eta_max | shape | L2 error | quad check | relative "
         "| error share on (0,d), n=64 | n=256 | node range |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for r in r1["rows"]:
        L.append(f"| {r['m']} | {r['family']} | {r['eta_max']} | {r['shape']} "
                 f"| {r['l2_abs']:.4e} | {r['l2_quad_check']:.4e} "
                 f"| {r['l2_rel']:.4e} | {r['share_below_delta_n64']:.3f} "
                 f"| {r['share_below_delta_n256']:.3f} "
                 f"| {r['s_min']} – {r['s_max']} |")
    L += ["", f"Decay in m (log-log slope): **{r1['decay_slope_in_m']}**.", "",
          "The two share columns give the fraction of the SQUARED error living "
          "on (0, delta), which a", "lattice of time step delta never visits.", "",
          "## Phase 2 — the achievable floor in L2(0,T)", "",
          "| m | L2 error | relative | active factors | gain vs AJEE |",
          "|---|---|---|---|---|"]
    for r in r2["rows"]:
        L.append(f"| {r['m']} | {r['l2_abs']:.4e} | {r['l2_rel']:.4e} "
                 f"| {r['active_factors']} | x{r['gain_vs_ajee']} |")
    L += ["", f"Decay in m (log-log slope): **{r2['decay_slope_in_m']}**.", "",
          "## Phase 3 — the norm the discrete kernel sees, L2(delta,T)", ""]
    for n in NS:
        d3 = r3["by_n"][n]
        L += [f"### n = {n}, delta = {d3['delta']:.6f}, "
              f"||K||_L2(delta,T) = {d3['norm_K_on_interval']:.6f}", "",
              "| m | AJEE error | quad check | AJEE relative | floor error "
              "| floor relative |", "|---|---|---|---|---|---|"]
        for r in d3["rows"]:
            L.append(f"| {r['m']} | {r['l2_abs_ajee']:.4e} "
                     f"| {r['l2_quad_check']:.4e} | {r['l2_rel_ajee']:.4e} "
                     f"| {r['l2_abs_floor']:.4e} | {r['l2_rel_floor']:.4e} |")
        L += ["", f"Decay slopes in m: AJEE **{d3['decay_slope_ajee']}**, "
                  f"floor **{d3['decay_slope_floor']}**.", ""]
    L += ["## Phase 4 — the crossing", "", f"{r4['lattice_rate_note']}.", ""]
    for row in r4["rows"]:
        L += [f"### n = {row['n']}  (delta = {row['delta']:.6f}, lattice error "
              f"= **{row['lattice_error']:.5f}**)", "",
              "| criterion | smallest m (AJEE) | smallest m (floor) |",
              "|---|---|---|"]
        for name, c in row["crossings"].items():
            L.append(f"| {name} | {c['ajee']} | {c['floor']} |")
        L.append("")
    L += ["### For its own sake: m needed for a given relative kernel accuracy",
          "", "| target relative error | AJEE, L2(0,T) | floor, L2(0,T) "
          "| floor, L2(delta,T) at n=256 |", "|---|---|---|---|"]
    for tgt, c in r4["relative_accuracy_tiers"].items():
        L.append(f"| {tgt:g} | {c['m_ajee_L2_0T']} | {c['m_floor_L2_0T']} "
                 f"| {c['m_floor_L2_deltaT_n256']} |")
    L += ["", "`None` means no m in the tested range reaches that target.", "",
          "## Phase 5 — the continuous covariance surface", "",
          f"Lift {r5['lift_optimised_for']}-optimised.  Target: every error 0, "
          "every ratio 1.", "",
          "| m | max relative error | at (u,v) | relative Frobenius | "
          "Var[V_T] ratio |", "|---|---|---|---|---|"]
    for r in r5["rows"]:
        L.append(f"| {r['m']} | {r['max_rel_error']:.4e} "
                 f"| {r['max_rel_error_at']} | {r['rel_frobenius']:.4e} "
                 f"| {r['var_T_ratio']:.6f} |")
    L += ["", "For comparison, the one-step lattice of Proposition 8.3 "
          "(target 1.0):", "",
          "| n | " + " | ".join(str(k) for k in
                                r5["onestep_lattice_var_ratio"]) + " |",
          "|---|" + "---|" * len(r5["onestep_lattice_var_ratio"]),
          "| Var[V_T] ratio | " + " | ".join(
              str(v) for v in r5["onestep_lattice_var_ratio"].values()) + " |",
          "", "## Phase 6 — the DISCRETE covariance, which decides the matter",
          "", "The lattice covariance is C_kl = delta sum_{j<k^l} G((k-j)d) "
          "G((l-j)d) for the discrete", "kernel G.  Two conventions for G are "
          "reported: the cell average of K, which is its L2", "projection on the "
          "grid, and the naive left endpoint.  The lift's own error is the",
          "'vs true discrete' block, and it must be judged against the "
          "discretisation error the", "scheme already carries with the TRUE "
          "kernel — the line under each table.", ""]
    for n in NS:
        d6 = r6["by_n"][n]
        L.append(f"### n = {n}, delta = {d6['delta']:.6f}")
        L.append("")
        for mode, _ in MODES:
            blk = d6["by_mode"][mode]
            b = blk["scheme_own_discretisation_error"]
            L += [f"#### discrete kernel: {blk['label']}", "",
                  "| lift optimised for | m | vs true discrete: max rel "
                  "| Frobenius | Var[V_T] ratio | total vs continuous: Frobenius "
                  "| Var ratio |", "|---|---|---|---|---|---|---|"]
            for r in blk["rows"]:
                L.append(f"| {r['lift']} | {r['m']} | {r['max_rel']:.4e} "
                         f"| {r['frob']:.4e} | {r['varT']:.6f} "
                         f"| {r['tot_frob']:.4e} | {r['tot_varT']:.6f} |")
            L += ["", f"Baseline — the scheme's OWN discretisation error with "
                  f"the TRUE kernel at this n: Frobenius "
                  f"**{b['rel_frobenius']:.4e}**, Var[V_T] ratio "
                  f"**{b['var_T_ratio']:.6f}**.  A lift whose error sits below "
                  "this line is not the binding constraint.", ""]
    L += ["## Phase 7 — cost in the target accuracy", "",
          f"{r7['n_in_eps']}.  Route A' alone: "
          f"**eps^-{r7['route_aprime_alone']}**.", "",
          "| m | states | cost eps^-q | with Route A' inside |",
          "|---|---|---|---|"]
    for r in r7["rows"]:
        L.append(f"| {r['m']} | {r['states']} | {r['cost_exponent_q']} "
                 f"| {r['with_aprime']} |")
    L += ["", pr.timing_table_md()]
    return "\n".join(L)


def main() -> None:
    meta = {"item": "(N4)/(B1) — quantify the Markovian lift",
            "kernel": "Riemann-Liouville, K(u) = u^h / Gamma(1+h)",
            "H": H, "h": h, "kappa": KAPPA, "T": T,
            "m_values": list(MS), "n_values": list(NS),
            "lattice_rate_exponent": 0.5 * (h + KAPPA),
            "norm_K_L2_0T": round(NORM_K, 6)}
    with Progress("route-b", total_phases=7, meta=meta) as pr:
        r1 = phase1(pr); s1 = r1.pop("_store"); pr.result("phase1_ajee_L2_0T", r1)
        r2 = phase2(pr, s1); s2 = r2.pop("_store"); pr.result("phase2_floor_L2_0T", r2)
        r3 = phase3(pr); s3 = r3.pop("_store"); pr.result("phase3_L2_deltaT", r3)
        r4 = phase4(pr, r1, r2, r3); pr.result("phase4_crossing", r4)
        r5 = phase5(pr, s2); pr.result("phase5_continuous_covariance", r5)
        r6 = phase6(pr, s2, s3); pr.result("phase6_discrete_covariance", r6)
        r7 = phase7(pr); pr.result("phase7_cost", r7)
        pr.write_results_md(build_results(pr, r1, r2, r3, r4, r5, r6, r7))


if __name__ == "__main__":
    main()
