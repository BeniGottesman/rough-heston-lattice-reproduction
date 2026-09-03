#!/usr/bin/env python3
"""
run_mc_vs_tree_sweep — the tree against a trustworthy Monte-Carlo, swept over the
model parameters, so that the error can be read as a function of them.

Model: rough Bergomi (the paper's Section "Rough Bergomi", Example 4.14), which
is what the code implements.  Rough Heston is a different model -- a fractional
CIR variance -- and would be a new implementation, not a parameter change.

Why the reference had to be rebuilt first.  The previous Monte-Carlo discretised
the fractional driver with the left-endpoint kernel, and finding F022 measures
what that costs: it loses 25% of Var[V^H_T] at H = 0.1 with 512 steps, and its
price sat 0.062 away from the exact-driver value -- more than its own error bar.
`mc_reference` samples the driver with its EXACT covariance (Cholesky of the true
covariance matrix) and adds an eta = 0 control variate, which is exact because
the control's mean is Black--Scholes in closed form.  The band drops from +-0.048
to +-0.007 and the driver's variance is right to machine precision.

The barrier.  Without a variance barrier the lattice is not merely expensive but
infeasible: at eta = 0.5, n = 64 the effective volatility ceiling is 33, and the
state space explodes.  So a barrier is imposed at three standard deviations of
the TRUE driver, zmax = 3/sqrt(2H), which is comparable across parameters -- and
the SAME cap is applied to the Monte-Carlo driver.  That splits the error in two:

    scheme error       tree - MC(capped)     like for like
    absorption error   MC(capped) - MC(true) what the barrier itself costs
    total error        tree - MC(true)       what one would actually pay

Reading rule: the error must SHRINK as n grows.  A positive slope of |error| in n
means the scheme converges to the wrong limit, which is the refutation of
F008 -- here re-measured against a reference seven times tighter, and resolved by
parameter.

Five phases; every number lands in RESULTS.md and is registered in FINDINGS.md.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from progress import Progress            # noqa: E402
import mc_reference as mc                # noqa: E402
import route_aprime as ra                # noqa: E402

NS = (8, 16, 32, 64)
CAP = 3.0                                # barrier, in sd of the true driver
NFINE = 512
PATHS = 400_000

BASE = {"H": 0.10, "eta": 0.30, "rho": -0.70}
SWEEPS = (
    ("H  (roughness)", "H", (0.05, 0.10, 0.20, 0.30, 0.45)),
    ("eta  (vol of vol)", "eta", (0.00, 0.10, 0.20, 0.30, 0.50)),
    ("rho  (correlation)", "rho", (-0.90, -0.70, -0.40, 0.00)),
)


def mref_for(n: int) -> int:
    return max(4, int(math.ceil(4.0 * math.sqrt(n / 8.0))))


def zmax_for(H: float) -> float:
    """cap sd of the true driver, in the units of the tree's discrete zc."""
    return CAP / math.sqrt(2.0 * H)


_MC: dict = {}


def mc_price(H: float, eta: float, rho: float, capped: bool) -> dict:
    key = (round(H, 6), round(eta, 6), round(rho, 6), capped)
    if key not in _MC:
        _MC[key] = mc.european_put_mc(H, eta, rho, nfine=NFINE, paths=PATHS,
                                      cap=CAP if capped else None)
    return _MC[key]


def phase1(pr: Progress) -> dict:
    """Establish that the reference deserves to be believed."""
    pr.phase("reference quality: is the Monte-Carlo trustworthy?",
             total=4 + len(SWEEPS[0][2]))
    rows = []
    for nf in (256, 512, 1024):
        r = mc.european_put_mc(BASE["H"], BASE["eta"], BASE["rho"],
                               nfine=nf, paths=200_000)
        rows.append({"test": f"nfine = {nf} (X time-stepping)",
                     "price": r["price"], "ci95": r["ci95"],
                     "driver_var_ratio": r["driver_var_ratio"]})
        pr.tick(note=f"nfine={nf}")
    old = mc.european_put_mc(BASE["H"], BASE["eta"], BASE["rho"], nfine=NFINE,
                             paths=200_000, exact_driver=False)
    rows.append({"test": "old left-endpoint driver, same paths",
                 "price": old["price"], "ci95": old["ci95"],
                 "driver_var_ratio": old["driver_var_ratio"]})
    pr.tick(note="old driver")
    zero = []
    for H in SWEEPS[0][2]:
        r = mc.european_put_mc(H, 0.0, BASE["rho"], nfine=NFINE, paths=20_000)
        zero.append({"H": H, "price": r["price"], "ci95": r["ci95"],
                     "bs_exact": mc.bs_put(),
                     "abs_err": abs(r["price"] - mc.bs_put())})
        pr.tick(note=f"eta=0 H={H}")
    return {"rows": rows, "eta0_control": zero, "bs_exact": mc.bs_put()}


def one_set(pr: Progress, H: float, eta: float, rho: float) -> dict:
    mt = mc_price(H, eta, rho, capped=False)
    mcap = mc_price(H, eta, rho, capped=True)
    z = zmax_for(H)
    rows = []
    for n in NS:
        tv = ra.route_aprime_european_put(n, H, eta, rho, zmax=z,
                                          mref=mref_for(n))
        rows.append({"n": n, "mref": mref_for(n), "tree": tv,
                     "gap_scheme": tv - mcap["price"],
                     "gap_total": tv - mt["price"]})
        pr.tick(note=f"H={H} eta={eta} rho={rho} n={n}")
        pr.partial.update({f"tree_H{H}_e{eta}_r{rho}_n{n}": round(tv, 4)})
    ns = np.array([r["n"] for r in rows], float)
    g = np.abs([r["gap_scheme"] for r in rows])
    slope = (float(np.polyfit(np.log(ns), np.log(np.maximum(g, 1e-12)), 1)[0])
             if np.all(g > 0) else float("nan"))
    return {"H": H, "eta": eta, "rho": rho, "zmax": z,
            "mc_true": mt["price"], "mc_true_ci": mt["ci95"],
            "mc_capped": mcap["price"], "mc_capped_ci": mcap["ci95"],
            "absorption": mcap["price"] - mt["price"],
            "rows": rows,
            "gap_slope_in_n": round(slope, 3),
            "verdict": ("converging" if slope < -0.15 else
                        "diverging" if slope > 0.15 else "flat")}


def sweep_phase(pr: Progress, label: str, key: str, values) -> dict:
    pr.phase(f"sweep {label}", total=len(values) * len(NS))
    out = []
    for v in values:
        p = dict(BASE)
        p[key] = v
        out.append(one_set(pr, p["H"], p["eta"], p["rho"]))
        pr.result(f"sweep_{key}_{v}",
                  {"mc_true": round(out[-1]["mc_true"], 4),
                   "slope": out[-1]["gap_slope_in_n"],
                   "verdict": out[-1]["verdict"]})
    return {"label": label, "key": key, "values": list(values), "sets": out}


def phase5(pr: Progress, sweeps: list) -> dict:
    pr.phase("synthesis", total=len(sweeps))
    tbl = []
    for s in sweeps:
        for st in s["sets"]:
            g64 = next(r for r in st["rows"] if r["n"] == 64)
            g8 = next(r for r in st["rows"] if r["n"] == 8)
            tbl.append({"sweep": s["key"], "value": st[s["key"]],
                        "H": st["H"], "eta": st["eta"], "rho": st["rho"],
                        "mc_true": st["mc_true"], "mc_ci": st["mc_true_ci"],
                        "absorption": st["absorption"],
                        "gap_n8": g8["gap_scheme"], "gap_n64": g64["gap_scheme"],
                        "gap_n64_in_ci": (abs(g64["gap_scheme"])
                                          / max(st["mc_capped_ci"], 1e-12)),
                        "slope": st["gap_slope_in_n"],
                        "verdict": st["verdict"]})
        pr.tick(note=s["key"])
    return {"table": tbl}


def build_results(pr, r1, sweeps, r5) -> str:
    L = ["# Rough Bergomi: the tree against a trustworthy Monte-Carlo, by "
         "parameter", "",
         f"S0 = K = {mc.S0:.0f}, xi0 = {mc.XI0:.4f} (sigma = "
         f"{math.sqrt(mc.XI0):.2f}), T = {mc.T:.0f}.  Baseline "
         f"H = {BASE['H']}, eta = {BASE['eta']}, rho = {BASE['rho']}.",
         f"Grids n = {list(NS)}, refinement mref = max(4, ceil(4 sqrt(n/8))).",
         f"Barrier at {CAP:.0f} sd of the TRUE driver, zmax = "
         f"{CAP:.0f}/sqrt(2H); the same cap is applied to the Monte-Carlo, so "
         "the scheme error and the",
         "absorption error are separated.  Reference: exact-covariance driver "
         f"+ eta=0 control variate, nfine = {NFINE}, {PATHS:,} paths.", "",
         "**Reading rule: the error must SHRINK as n grows.**  A positive slope "
         "means the scheme", "converges to the wrong limit.", "",
         "## Phase 1 — is the reference trustworthy?", "",
         "| test | price | 95% band | driver Var[V^H_T] ratio (target 1) |",
         "|---|---|---|---|"]
    for r in r1["rows"]:
        L.append(f"| {r['test']} | {r['price']:.4f} | ±{r['ci95']:.4f} "
                 f"| {r['driver_var_ratio']:.4f} |")
    L += ["", f"The eta = 0 control must return Black--Scholes = "
              f"**{r1['bs_exact']:.6f}** exactly, with zero variance, because "
              "the control variate", "and the payoff coincide pathwise there:",
          "", "| H | price at eta=0 | 95% band | absolute error vs BS |",
          "|---|---|---|---|"]
    for r in r1["eta0_control"]:
        L.append(f"| {r['H']} | {r['price']:.6f} | ±{r['ci95']:.1e} "
                 f"| {r['abs_err']:.1e} |")
    for s in sweeps:
        L += ["", f"## Sweep: {s['label']}", ""]
        for st in s["sets"]:
            L += [f"### {s['key']} = {st[s['key']]}  "
                  f"(H={st['H']}, eta={st['eta']}, rho={st['rho']}, "
                  f"zmax={st['zmax']:.2f})", "",
                  f"MC true = **{st['mc_true']:.4f}** ±{st['mc_true_ci']:.4f} · "
                  f"MC capped = **{st['mc_capped']:.4f}** "
                  f"±{st['mc_capped_ci']:.4f} · absorption cost = "
                  f"**{st['absorption']:+.4f}**", "",
                  "| n | mref | tree | scheme error (tree − MC capped) "
                  "| total error (tree − MC true) |", "|---|---|---|---|---|"]
            for r in st["rows"]:
                L.append(f"| {r['n']} | {r['mref']} | {r['tree']:.4f} "
                         f"| {r['gap_scheme']:+.4f} | {r['gap_total']:+.4f} |")
            L += ["", f"slope of |scheme error| in n: **{st['gap_slope_in_n']}** "
                      f"→ **{st['verdict']}**", ""]
    L += ["## Synthesis", "",
          "| sweep | value | MC true | ±band | absorption | error n=8 "
          "| error n=64 | n=64 in bands | slope | verdict |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for r in r5["table"]:
        L.append(f"| {r['sweep']} | {r['value']} | {r['mc_true']:.4f} "
                 f"| ±{r['mc_ci']:.4f} | {r['absorption']:+.4f} "
                 f"| {r['gap_n8']:+.4f} | {r['gap_n64']:+.4f} "
                 f"| {r['gap_n64_in_ci']:.1f}x | {r['slope']} "
                 f"| {r['verdict']} |")
    L += ["", pr.timing_table_md()]
    return "\n".join(L)


def main() -> None:
    meta = {"model": "rough Bergomi", "S0": mc.S0, "K": mc.KSTRIKE,
            "xi0": mc.XI0, "T": mc.T, "baseline": BASE, "n_values": list(NS),
            "barrier_sd_of_true_driver": CAP, "mc_nfine": NFINE,
            "mc_paths": PATHS,
            "reference": "exact-covariance driver + eta=0 control variate"}
    with Progress("mc-vs-tree-sweep", total_phases=2 + len(SWEEPS),
                  meta=meta) as pr:
        r1 = phase1(pr); pr.result("phase1_reference_quality", r1)
        sweeps = [sweep_phase(pr, lab, key, vals) for lab, key, vals in SWEEPS]
        for s in sweeps:
            pr.result(f"sweep_{s['key']}_sets",
                      [{k: (round(v, 5) if isinstance(v, float) else v)
                        for k, v in st.items() if k != "rows"}
                       for st in s["sets"]])
        r5 = phase5(pr, sweeps); pr.result("phase5_synthesis", r5)
        pr.write_results_md(build_results(pr, r1, sweeps, r5))


if __name__ == "__main__":
    main()
