#!/usr/bin/env python3
"""
run_rough_heston — item (N6): the independent reference, and what it costs.

Every price comparison in this project so far has set our lattice against our own
Monte-Carlo, so a shared misreading of the model would cancel and stay invisible.
Rough Heston breaks that: its characteristic function solves a fractional Riccati
equation and the price follows by Fourier inversion, with no Monte-Carlo anywhere
in it.  This run establishes that our Fourier pricer is correct, then measures
Fourier against Monte-Carlo on the same model -- price, accuracy and wall time.

Five phases.

  1  Exact controls.  (a) nu = lambda = 0 collapses the model to Black--Scholes,
     which the Lewis inversion must reproduce in closed form.  (b) phi(0) = 1 and
     phi(-i) = 1, the latter being the martingale property of S.  (c) nu = 0 with
     lambda > 0 leaves a deterministic but non-trivial fractional variance, whose
     price is Black--Scholes with the integrated variance -- exact FOR EVERY H,
     so it tests the fractional machinery without degrading the roughness.
  2  The reference's own accuracy: convergence of the Riccati solver in the
     number of steps and of the Fourier rule in (u_max, n_q), on a genuinely
     rough case where no closed form exists.
  3  Fourier against Monte-Carlo across H, with wall time for each.
  4  The Monte-Carlo's bias floor: its error against the Fourier value as the
     number of steps grows at fixed paths, which no number of paths removes.
  5  Cost to a target accuracy for each method.

Every number lands in RESULTS.md and is registered in FINDINGS.md.
"""
from __future__ import annotations

import math
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from progress import Progress            # noqa: E402
import rough_heston as rh                # noqa: E402

# comparable to the project's rough Bergomi numerics: sigma = 30% at t = 0
P = {"V0": 0.09, "theta": 0.09, "lam": 1.5, "nu": 0.30, "rho": -0.70}
HS = (0.05, 0.10, 0.20, 0.30, 0.45)
T = 1.0
REF_STEPS = 3200                         # the Fourier value taken as truth
MC_STEPS = 200
MC_PATHS = 400_000


def fourier(H: float, steps: int, nu_max: float = 60.0, nq: int = 200,
            **over) -> float:
    p = dict(P)
    p.update(over)
    return rh.put_fourier(H, p["V0"], p["theta"], p["lam"], p["nu"], p["rho"],
                          T=T, steps=steps, nu_max=nu_max, nq=nq)


def timed(fn, *a, **k):
    t = time.time()
    out = fn(*a, **k)
    return out, time.time() - t


def phase1(pr: Progress) -> dict:
    pr.phase("exact controls", total=3 + len(HS) + 4)
    # (a) nu = lambda = 0 -> Black-Scholes
    a = []
    for V0 in (0.04, 0.09, 0.25):
        ex = rh.bs_put(rh.S0, rh.KSTRIKE, V0 * T)
        got = fourier(0.10, 200, V0=V0, theta=V0, lam=0.0, nu=0.0)
        a.append({"V0": V0, "fourier": got, "exact": ex, "abs_err": abs(got - ex)})
        pr.tick(note=f"BS V0={V0}")
    # (b) phi(0) = phi(-i) = 1
    b = []
    for H in (0.05, 0.10, 0.30, 0.45):
        c = rh.char_fn(np.array([0.0 + 0j, -1j]), H, P["V0"], P["theta"],
                       P["lam"], P["nu"], P["rho"], T, 400)
        b.append({"H": H, "phi0_err": abs(c[0] - 1.0),
                  "phi_minus_i_err": abs(c[1] - 1.0)})
        pr.tick(note=f"martingale H={H}")
    # (c) nu = 0, lambda > 0: deterministic fractional variance
    c_rows = []
    for H in HS:
        ex, tot = rh.put_deterministic(H, P["V0"], 0.16, P["lam"], T,
                                       steps=6000)
        errs = []
        for st in (200, 400, 800):
            got = fourier(H, st, theta=0.16, nu=0.0)
            errs.append(abs(got - ex))
        order = math.log2(errs[0] / errs[2]) / 2.0
        c_rows.append({"H": H, "exact": ex, "integrated_variance": tot,
                       "err_200": errs[0], "err_400": errs[1],
                       "err_800": errs[2], "observed_order": round(order, 2),
                       "predicted_order_min2_1plusalpha": round(min(2.0, 1.5 + H), 2)})
        pr.tick(note=f"deterministic H={H}")
    return {"bs_collapse": a, "martingale": b, "deterministic_variance": c_rows}


def phase2(pr: Progress) -> dict:
    pr.phase("the reference's own accuracy", total=6 + 6)
    ref = fourier(0.10, REF_STEPS)
    steps_rows = []
    for st in (100, 200, 400, 800, 1600):
        v, sec = timed(fourier, 0.10, st)
        steps_rows.append({"steps": st, "price": v, "err_vs_ref": abs(v - ref),
                           "seconds": round(sec, 3)})
        pr.tick(note=f"steps={st}")
    pr.tick()
    quad_rows = []
    for nu_max, nq in ((15, 100), (20, 100), (30, 150), (60, 200),
                       (60, 400), (120, 400)):
        v = fourier(0.10, 800, nu_max=nu_max, nq=nq)
        quad_rows.append({"u_max": nu_max, "n_q": nq, "price": v,
                          "err_vs_ref": abs(v - ref)})
        pr.tick(note=f"u_max={nu_max}")
    return {"reference_price": ref, "reference_steps": REF_STEPS,
            "steps_convergence": steps_rows, "quadrature": quad_rows}


def phase3(pr: Progress) -> dict:
    pr.phase("Fourier against Monte-Carlo, with wall time", total=2 * len(HS))
    rows = []
    for H in HS:
        fo, tf = timed(fourier, H, REF_STEPS)
        pr.tick(note=f"fourier H={H}")
        mc, tm = timed(rh.put_mc, H, P["V0"], P["theta"], P["lam"], P["nu"],
                       P["rho"], T=T, steps=MC_STEPS, paths=MC_PATHS,
                       chunk=25_000)
        pr.tick(note=f"mc H={H}")
        gap = mc["price"] - fo
        rows.append({"H": H, "fourier": fo, "fourier_seconds": round(tf, 3),
                     "mc": mc["price"], "mc_ci95": mc["ci95"],
                     "mc_seconds": round(tm, 2),
                     "gap_mc_minus_fourier": gap,
                     "gap_in_ci": abs(gap) / max(mc["ci95"], 1e-12),
                     "negative_variance_hits": mc["negative_variance_hits"],
                     "speedup_fourier": round(tm / max(tf, 1e-9), 1)})
        pr.result(f"H{H}", {"fourier": round(fo, 4), "mc": round(mc["price"], 4),
                            "gap": round(gap, 4)})
    return {"mc_steps": MC_STEPS, "mc_paths": MC_PATHS, "rows": rows}


def phase4(pr: Progress) -> dict:
    pr.phase("the Monte-Carlo's bias floor", total=4)
    H = 0.10
    ref = fourier(H, REF_STEPS)
    rows = []
    for st in (50, 100, 200, 400):
        mc, sec = timed(rh.put_mc, H, P["V0"], P["theta"], P["lam"], P["nu"],
                        P["rho"], T=T, steps=st, paths=200_000, chunk=25_000)
        rows.append({"steps": st, "price": mc["price"], "ci95": mc["ci95"],
                     "bias_vs_fourier": mc["price"] - ref,
                     "bias_in_ci": abs(mc["price"] - ref) / max(mc["ci95"], 1e-12),
                     "negative_variance_hits": mc["negative_variance_hits"],
                     "seconds": round(sec, 2)})
        pr.tick(note=f"steps={st}")
    return {"H": H, "fourier_reference": ref, "paths": 200_000, "rows": rows}


def phase5(pr: Progress, r2: dict, r3: dict, r4: dict) -> dict:
    pr.phase("cost to a target accuracy", total=3)
    H = 0.10
    fo_rows = r2["steps_convergence"]
    mc_row = next(r for r in r3["rows"] if r["H"] == H)
    out = []
    for tgt in (1e-2, 1e-3, 1e-4):
        # Fourier: smallest tested step count whose error is below the target
        hit = next((r for r in fo_rows if r["err_vs_ref"] < tgt), None)
        f_sec = hit["seconds"] if hit else None
        f_steps = hit["steps"] if hit else None
        # Monte-Carlo: paths needed for a 95% half-width of tgt, at fixed steps,
        # scaling as 1/sqrt(paths) -- and only if the bias floor allows it at all
        need = mc_row["mc_ci95"] ** 2 / tgt ** 2 * r3["mc_paths"]
        sec = mc_row["mc_seconds"] * need / r3["mc_paths"]
        floor = min(abs(x["bias_vs_fourier"]) for x in r4["rows"])
        out.append({"target": tgt,
                    "fourier_steps": f_steps, "fourier_seconds": f_sec,
                    "mc_paths_needed": f"{need:.3g}",
                    "mc_seconds_needed": round(sec, 1),
                    "mc_reachable": bool(floor < tgt),
                    "speedup": (None if not f_sec else round(sec / f_sec, 1))})
        pr.tick(note=f"target={tgt}")
    return {"H": H, "mc_bias_floor": min(abs(x["bias_vs_fourier"])
                                         for x in r4["rows"]), "rows": out}


def build_results(pr, r1, r2, r3, r4, r5) -> str:
    L = ["# Rough Heston: an independent reference, and what it costs", "",
         f"Model parameters {P}, T = {T}, S0 = K = {rh.S0:.0f}.",
         "Fourier price by the Lewis representation on the fractional-Riccati "
         "characteristic function", "(El Euch--Rosenbaum), solved by the "
         "Diethelm--Ford--Freed predictor--corrector.", "",
         "## Phase 1 --- exact controls", "",
         "### (a) nu = lambda = 0 collapses to Black--Scholes", "",
         "| V0 | Fourier | exact | absolute error |", "|---|---|---|---|"]
    for r in r1["bs_collapse"]:
        L.append(f"| {r['V0']} | {r['fourier']:.10f} | {r['exact']:.10f} "
                 f"| {r['abs_err']:.2e} |")
    L += ["", "### (b) phi(0) = 1 and phi(-i) = 1 (the martingale property)", "",
          "| H | error on phi(0) | error on phi(-i) |", "|---|---|---|"]
    for r in r1["martingale"]:
        L.append(f"| {r['H']} | {r['phi0_err']:.2e} | {r['phi_minus_i_err']:.2e} |")
    L += ["", "### (c) nu = 0, lambda > 0: deterministic fractional variance",
          "", "This is the control that matters: the variance is non-trivial and "
          "fractional, only its", "randomness is switched off, so the whole "
          "fractional machinery is exercised at every H.", "",
          "| H | exact | err, 200 steps | 400 | 800 | observed order "
          "| predicted min(2, 1+alpha) |", "|---|---|---|---|---|---|---|"]
    for r in r1["deterministic_variance"]:
        L.append(f"| {r['H']} | {r['exact']:.8f} | {r['err_200']:.2e} "
                 f"| {r['err_400']:.2e} | {r['err_800']:.2e} "
                 f"| **{r['observed_order']}** "
                 f"| {r['predicted_order_min2_1plusalpha']} |")
    L += ["", "## Phase 2 --- the reference's own accuracy (H = 0.10)", "",
          f"Taken as truth: the Fourier price at {r2['reference_steps']} steps, "
          f"**{r2['reference_price']:.8f}**.", "",
          "| Riccati steps | price | error | seconds |", "|---|---|---|---|"]
    for r in r2["steps_convergence"]:
        L.append(f"| {r['steps']} | {r['price']:.8f} | {r['err_vs_ref']:.2e} "
                 f"| {r['seconds']} |")
    L += ["", "| u_max | n_q | price | error |", "|---|---|---|---|"]
    for r in r2["quadrature"]:
        L.append(f"| {r['u_max']} | {r['n_q']} | {r['price']:.8f} "
                 f"| {r['err_vs_ref']:.2e} |")
    L += ["", "## Phase 3 --- Fourier against Monte-Carlo", "",
          f"Monte-Carlo: {r3['mc_paths']:,} antithetic paths, {r3['mc_steps']} "
          "steps, with the nu = 0 control variate.", "",
          "| H | Fourier | s | Monte-Carlo | 95% band | s | MC - Fourier "
          "| in bands | negative-variance hits | Fourier faster by |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for r in r3["rows"]:
        L.append(f"| {r['H']} | {r['fourier']:.6f} | {r['fourier_seconds']} "
                 f"| {r['mc']:.4f} | ±{r['mc_ci95']:.4f} | {r['mc_seconds']} "
                 f"| {r['gap_mc_minus_fourier']:+.4f} | {r['gap_in_ci']:.1f}x "
                 f"| {r['negative_variance_hits']:,} "
                 f"| **{r['speedup_fourier']}x** |")
    L += ["", "## Phase 4 --- the Monte-Carlo's bias floor", "",
          f"H = {r4['H']}, {r4['paths']:,} paths throughout; Fourier reference "
          f"**{r4['fourier_reference']:.6f}**.  The Euler--Volterra",
          "scheme truncates the variance at zero, and that bias is a function of "
          "the step size, not", "of the number of paths.", "",
          "| steps | price | 95% band | bias vs Fourier | in bands "
          "| negative-variance hits | seconds |",
          "|---|---|---|---|---|---|---|"]
    for r in r4["rows"]:
        L.append(f"| {r['steps']} | {r['price']:.4f} | ±{r['ci95']:.4f} "
                 f"| {r['bias_vs_fourier']:+.4f} | {r['bias_in_ci']:.1f}x "
                 f"| {r['negative_variance_hits']:,} | {r['seconds']} |")
    L += ["", "## Phase 5 --- cost to a target accuracy (H = 0.10)", "",
          f"Monte-Carlo bias floor over the step counts tested: "
          f"**{r5['mc_bias_floor']:.4f}**.  A target below that floor cannot be",
          "reached by adding paths at all.", "",
          "| target | Fourier: steps | seconds | MC: paths needed | seconds "
          "| reachable? | Fourier faster by |",
          "|---|---|---|---|---|---|---|"]
    for r in r5["rows"]:
        L.append(f"| {r['target']:g} | {r['fourier_steps']} "
                 f"| {r['fourier_seconds']} | {r['mc_paths_needed']} "
                 f"| {r['mc_seconds_needed']} "
                 f"| {'yes' if r['mc_reachable'] else '**no**'} "
                 f"| {r['speedup'] if r['speedup'] else '---'}x |")
    L += ["", pr.timing_table_md()]
    return "\n".join(L)


def main() -> None:
    meta = {"item": "(N6) an independent reference",
            "model": "rough Heston (El Euch-Rosenbaum)",
            "params": P, "T": T, "S0": rh.S0, "K": rh.KSTRIKE,
            "H_values": list(HS), "reference_riccati_steps": REF_STEPS,
            "mc_steps": MC_STEPS, "mc_paths": MC_PATHS,
            "riccati_scheme": "Diethelm-Ford-Freed predictor-corrector, "
                              "order min(2, 1+alpha)"}
    with Progress("rough-heston", total_phases=5, meta=meta) as pr:
        r1 = phase1(pr); pr.result("phase1_controls", r1)
        r2 = phase2(pr); pr.result("phase2_reference_accuracy", r2)
        r3 = phase3(pr); pr.result("phase3_fourier_vs_mc", r3)
        r4 = phase4(pr); pr.result("phase4_mc_bias_floor", r4)
        r5 = phase5(pr, r2, r3, r4); pr.result("phase5_cost", r5)
        pr.write_results_md(build_results(pr, r1, r2, r3, r4, r5))


if __name__ == "__main__":
    main()
