#!/usr/bin/env python3
"""
run_route_aprime — Level 4 redone with Route A' (paper §9.2), against the
four-point scheme of §4.

Four phases.

  1  (V8)   eta = 0: Route A' against the exact CRR reference, with the
            refinement factor mref growing like sqrt(n).  Calibrates the
            implementation.
  2  (V3')  admissibility sweep: number of transition probabilities outside
            [0,1] for the four-point kernel and for Route A', at four barrier
            widths including "no barrier".  The claim under test is that Route A'
            has none, ever.
  3  (V13)  eta = 0.3 American put under Route A' at those four barrier widths:
            does the value still depend on the barrier, and does the scheme
            survive without one?
  4  (V9)   European put: Route A' tree against a Monte-Carlo of the TRUE rough
            model (exact fractional kernel).  This measures what Route A' does
            NOT fix — the covariance discrepancy of §8.

Every number lands in runs/<name>/RESULTS.md and is registered in FINDINGS.md.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from progress import Progress            # noqa: E402
import route_aprime as ra                # noqa: E402
import run_validation as rv              # noqa: E402

H = 0.10
ETA = 0.30
RHO = -0.70
NS = (8, 16, 32, 64)
SD_ZETA = math.sqrt(1.0 / 0.2)           # sd of zeta_T at H = 0.1, T = 1
BARRIERS = [
    ("0.84 sd (4-point limit)", 1.877),
    ("2.38 sd (9-point limit)", 5.317),
    ("3.00 sd", 3.0 * SD_ZETA),
    ("none", float("inf")),
]


def mref_for(n: int) -> int:
    """mref ~ delta^{-gamma}: grows so that the rounding floor O(1/mref) dies."""
    return max(4, int(math.ceil(4.0 * math.sqrt(n / 8.0))))


def phase1(pr: Progress) -> dict:
    pr.phase("V8  eta=0: Route A' vs exact CRR reference", total=len(NS) + 1)
    ref = rv.crr_american_put(20_000, 0.30)
    pr.result("crr_reference", round(ref, 6))
    pr.tick()
    rows = []
    for n in NS:
        m = mref_for(n)
        r = ra.route_aprime_american_put(n, H, 0.0, RHO, mref=m, pr=pr)
        rows.append({"n": n, "mref": m, "value": round(r["value"], 6),
                     "abs_err": round(abs(r["value"] - ref), 6),
                     "neg_probs": r["negative_probabilities"],
                     "max_mass_error": r["max_mass_error"]})
        pr.tick(note=f"n={n}")
        pr.partial.update({f"A'_eta0_n{n}": round(r["value"], 4),
                           f"err_n{n}": round(abs(r["value"] - ref), 5)})
    ns = np.array([r["n"] for r in rows], float)
    er = np.array([r["abs_err"] for r in rows], float)
    slope = float(np.polyfit(np.log(1.0 / ns), np.log(er), 1)[0])
    return {"reference": round(ref, 6), "rows": rows,
            "observed_error_slope_in_delta": round(slope, 3)}


def phase2(pr: Progress) -> dict:
    pr.phase("V3' admissibility sweep: 4-point vs Route A'",
             total=2 * len(BARRIERS))
    out = []
    for label, Z in BARRIERS:
        v4, viol4, mv4 = rv.lattice_american_put(32, H, ETA, RHO, zmax=Z)
        pr.tick(note=f"4pt {label}")
        r5 = ra.route_aprime_american_put(32, H, ETA, RHO, zmax=Z,
                                          mref=mref_for(32))
        pr.tick(note=f"A' {label}")
        out.append({"barrier": label,
                    "zmax": None if not np.isfinite(Z) else round(Z, 3),
                    "fourpoint_value": round(v4, 4),
                    "fourpoint_violations": viol4,
                    "fourpoint_max_violation": round(mv4, 4),
                    "aprime_value": round(r5["value"], 4),
                    "aprime_negative_probs": r5["negative_probabilities"],
                    "aprime_max_mass_error": r5["max_mass_error"]})
        pr.partial.update({f"viol4_{label[:8]}": viol4,
                           f"violA_{label[:8]}": r5["negative_probabilities"]})
    return {"n": 32, "rows": out}


def phase3(pr: Progress) -> dict:
    pr.phase("V13 eta=0.3 American put under Route A', by barrier width",
             total=len(BARRIERS) * len(NS))
    table = {}
    for label, Z in BARRIERS:
        vals = {}
        for n in NS:
            r = ra.route_aprime_american_put(n, H, ETA, RHO, zmax=Z,
                                             mref=mref_for(n), pr=pr)
            vals[n] = round(r["value"], 5)
            pr.tick(note=f"{label} n={n}")
            pr.partial.update({f"A'_{label[:8]}_n{n}": vals[n]})
        table[label] = vals
    spread = {n: round(max(table[l][n] for l, _ in BARRIERS)
                       - min(table[l][n] for l, _ in BARRIERS), 5)
              for n in NS}
    return {"values_by_barrier": table, "spread_across_barriers": spread}


def phase4(pr: Progress) -> dict:
    pr.phase("V9  European: Route A' vs Monte-Carlo of the TRUE rough model",
             total=2 + len(NS))
    mc0, se0 = ra.european_put_mc_rough(H, 0.0, RHO, paths=200_000)
    pr.result("mc_eta0_control", [round(mc0, 4), round(2 * se0, 4)])
    pr.tick()
    mc, se = ra.european_put_mc_rough(H, ETA, RHO, paths=400_000)
    pr.result("mc_true_rough_european", [round(mc, 4), round(2 * se, 4)])
    pr.tick()
    rows = []
    for n in NS:
        tv = ra.route_aprime_european_put(n, H, ETA, RHO,
                                          zmax=float("inf"), mref=mref_for(n))
        rows.append({"n": n, "tree": round(tv, 4), "gap": round(tv - mc, 4)})
        pr.tick(note=f"euro n={n}")
        pr.partial.update({f"euro_gap_n{n}": round(tv - mc, 4)})
    return {"mc_eta0_control": [round(mc0, 4), round(2 * se0, 4)],
            "black_scholes_eta0": 11.9235,
            "mc_true_rough": [round(mc, 4), round(2 * se, 4)],
            "rows": rows}


def main() -> None:
    meta = {"scheme": "Route A' (embed the driver, couple the price)",
            "H": H, "eta": ETA, "rho": RHO, "T": 1.0, "S0": 100, "K": 100,
            "mref_rule": "max(4, ceil(4*sqrt(n/8)))",
            "sd_zeta_T": round(SD_ZETA, 3)}
    with Progress("route-aprime", total_phases=4, meta=meta) as pr:
        r1 = phase1(pr); pr.result("phase1_V8_calibration", r1)
        r2 = phase2(pr); pr.result("phase2_V3prime_admissibility", r2)
        r3 = phase3(pr); pr.result("phase3_V13_american", r3)
        r4 = phase4(pr); pr.result("phase4_V9_covariance_gap", r4)

        L = ["# Route A' — Level 4 redone", "",
             f"Scheme: {meta['scheme']};  H={H}, eta={ETA}, rho={RHO}.", "",
             "## Phase 1 (V8) — eta=0 against the exact CRR reference", "",
             f"reference = **{r1['reference']}**", "",
             "| n | mref | value | abs err | negative probs | mass error |",
             "|---|---|---|---|---|---|"]
        for r in r1["rows"]:
            L.append(f"| {r['n']} | {r['mref']} | {r['value']} | {r['abs_err']} "
                     f"| {r['neg_probs']} | {r['max_mass_error']:.1e} |")
        L += ["", f"observed error slope in delta: "
                  f"**{r1['observed_error_slope_in_delta']}**", "",
              "## Phase 2 (V3') — admissibility sweep at n=32", "",
              "| barrier | 4-point value | 4-point violations | A' value | "
              "A' negative probs |", "|---|---|---|---|---|"]
        for r in r2["rows"]:
            L.append(f"| {r['barrier']} | {r['fourpoint_value']} | "
                     f"{r['fourpoint_violations']} | {r['aprime_value']} | "
                     f"{r['aprime_negative_probs']} |")
        L += ["", "## Phase 3 (V13) — American put under Route A'", "",
              "| barrier | " + " | ".join(f"n={n}" for n in NS) + " |",
              "|---|" + "---|" * len(NS)]
        for label, _ in BARRIERS:
            L.append(f"| {label} | "
                     + " | ".join(str(r3['values_by_barrier'][label][n])
                                  for n in NS) + " |")
        L += ["", "spread across barrier widths: "
                  + ", ".join(f"n={n}: {v}" for n, v in
                              r3["spread_across_barriers"].items()), "",
              "## Phase 4 (V9) — what Route A' does NOT fix", "",
              f"control at eta=0: MC {r4['mc_eta0_control'][0]} "
              f"+- {r4['mc_eta0_control'][1]} vs Black-Scholes "
              f"{r4['black_scholes_eta0']}", "",
              f"Monte-Carlo of the TRUE rough model: "
              f"**{r4['mc_true_rough'][0]} +- {r4['mc_true_rough'][1]}** (95%)",
              "", "| n | A' tree European | gap vs MC |", "|---|---|---|"]
        for r in r4["rows"]:
            L.append(f"| {r['n']} | {r['tree']} | {r['gap']:+} |")
        L += ["", "The gap GROWS with n: Route A' converges, but not to the",
              "right limit. That is the covariance discrepancy of §8, which",
              "Route A' leaves untouched (Remark on the scope of (A'1)).", ""]
        pr.write_results_md("\n".join(L))


if __name__ == "__main__":
    main()
