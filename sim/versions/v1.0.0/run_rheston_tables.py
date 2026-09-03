#!/usr/bin/env python3
"""
run_rheston_tables — rough Heston priced on the Akyildirim--Dolinsky--Soner grid.

This run is deliberately KEPT ASIDE from the paper.  The paper's model class needs
an autonomous driver, and rough Heston is outside it (Section 10.8), so nothing
here is a lattice result.  What it is: a reference table set for the rough Heston
model itself, laid out like the numerical section of

    E. Akyildirim, Y. Dolinsky, H. M. Soner,
    "Approximating stochastic volatility by recombinant trees",
    Ann. Appl. Probab. 24(5) 2176--2205, 2014   (arXiv:1205.3555),

whose parameter sets are those of Beliaeva and Nawalkha:

    European call and put   K = 100,  S0 in {90, 95, 100, 105, 110},
                            T in {1, 3, 6} months,  sqrt(V0) in {0.2, 0.3, 0.4},
                            r = 0.05, eta = 0.1, kappa = 3.0, theta = 0.04,
                            rho = -0.7
    American put            K = 100,  S0 in {90, 100, 110},
                            T in {1, 3, 6} months,  sqrt(V0) in {0.2, 0.4},
                            rho in {-0.1, -0.7},  r, eta, kappa, theta as above

and the roughness is added on top: H in {0.05, 0.10, 0.30} plus H = 0.50, which
IS classical Heston and is therefore the anchor -- there ADS publish the
analytical Heston prices and we must reproduce them.

Two independent engines.

  Fourier   sim/rough_heston.py: Lewis inversion of the fractional-Riccati
            characteristic function (El Euch--Rosenbaum).  No Monte-Carlo.  Only
            European.  Adapted here to r != 0 through the forward.
  MC        sim/cpp/rheston_mc.cpp: Euler--Volterra, antithetic, multithreaded
            C++, with Longstaff--Schwartz for the American problem on a
            regression sample and a disjoint valuation sample, so the American
            figure is a genuine lower bound with an honest standard error.

Five phases.

  1  The anchor: H = 0.5 against the analytical Heston prices PUBLISHED by ADS,
     for both engines.  An external reference, not one of ours.
  2  European put and call across H, Fourier against Monte-Carlo: price, standard
     error, signed error, wall time for each engine.
  3  American put across H by LSM, with the European put as a lower bound and the
     early-exercise premium isolated.
  4  The American CALL control: with r >= 0 and no dividend its true value is the
     European call, which Fourier gives independently, so the gap measures what
     the LSM exercise policy costs.
  5  The Euler--Volterra bias: refine the step count at fixed paths and watch the
     price move.  This is the error no number of paths removes.

Every number lands in RESULTS.md and is registered in FINDINGS.md.
"""
from __future__ import annotations

import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from progress import Progress            # noqa: E402
import rough_heston as rh                # noqa: E402

BIN = HERE / "cpp" / "build" / "rheston_mc"

# ------------------------------------------------------- the Beliaeva--Nawalkha set
K = 100.0
R = 0.05
ETA = 0.10                               # vol of vol
KAPPA = 3.0                              # mean reversion
THETA = 0.04                             # long-run variance
RHO = -0.70
MONTHS = {"1m": 1.0 / 12.0, "3m": 0.25, "6m": 0.5}
T_ADS = {"1m": 0.0833, "3m": 0.25, "6m": 0.5}      # as ADS print them
S0S = (90.0, 95.0, 100.0, 105.0, 110.0)
VOLS = (0.2, 0.3, 0.4)                   # sqrt(V0)
HS = (0.05, 0.10, 0.30, 0.50)

# American grid of ADS
S0S_AM = (90.0, 100.0, 110.0)
VOLS_AM = (0.2, 0.4)
RHOS_AM = (-0.1, -0.7)

STEPS = 200                              # ADS's smallest tree size
PATHS_EU = 400_000
PATHS_REG = 60_000
PATHS_VAL = 140_000
EX_STRIDE = 2                            # 100 exercise dates at STEPS = 200
SEED = 20260805

RIC_STEPS = 400                          # fractional-Riccati steps
NQ = 200
UMAX = 60.0


# --------------------------------------------------------- Fourier with a rate
def fourier_put(H, V0, T, s0, k=K, r=R, rho=RHO, eta=ETA, kappa=KAPPA,
                theta=THETA, steps=RIC_STEPS):
    """Rough Heston put with r != 0, through the forward.

    Under dS = r S dt + S sqrt(V) dW we have S_T = s0 e^{rT} M_T with M a
    martingale of the zero-rate model, so the put is e^{-rT} times the zero-rate
    put struck at k with spot s0 e^{rT}.
    """
    fwd = s0 * math.exp(r * T)
    p = rh.put_fourier(H, V0, theta, kappa, eta, rho, T=T, s0=fwd, k=k,
                       steps=steps, nu_max=UMAX, nq=NQ)
    return math.exp(-r * T) * p


def fourier_call(H, V0, T, s0, **kw):
    r = kw.get("r", R)
    k = kw.get("k", K)
    return fourier_put(H, V0, T, s0, **kw) + s0 - k * math.exp(-r * T)


# -------------------------------------------------------------- the MC harness
def mc_batch(cfgs: list[dict]) -> dict[str, dict]:
    """Run the C++ engine on a list of configs; return results keyed by id."""
    cols = ["id", "H", "V0", "theta", "kappa", "eta", "rho", "r", "T", "S0", "K",
            "steps", "paths_eu", "paths_reg", "paths_val", "ex_stride", "seed"]
    lines = [",".join(cols)]
    for c in cfgs:
        lines.append(",".join(str(c[k]) for k in cols))
    proc = subprocess.run([str(BIN)], input="\n".join(lines) + "\n",
                          capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"rheston_mc failed: {proc.stderr[-2000:]}")
    out = {}
    rows = [r for r in proc.stdout.strip().split("\n") if r]
    head = rows[0].split(",")
    for row in rows[1:]:
        f = row.split(",")
        d = dict(zip(head, f))
        rec = {k: (v if k == "id" else float(v)) for k, v in d.items()}
        out[rec["id"]] = rec
    return out


def cfg(cid, H, V0, T, s0, rho=RHO, steps=STEPS, american=True,
        paths_eu=PATHS_EU, seed=SEED):
    return {"id": cid, "H": H, "V0": V0, "theta": THETA, "kappa": KAPPA,
            "eta": ETA, "rho": rho, "r": R, "T": T, "S0": s0, "K": K,
            "steps": steps, "paths_eu": paths_eu,
            "paths_reg": PATHS_REG if american else 0,
            "paths_val": PATHS_VAL if american else 0,
            "ex_stride": EX_STRIDE, "seed": seed}


# ----------------------------------------- ADS's own published analytical column
def ads_published() -> tuple[dict, dict]:
    """Parse the Analytical-solution column of ADS Tables 1 and 2.

    Returns {(S0, vol, Tlabel): price} for the put and for the call.  These are
    Heston closed-form values published by a third party, so they are the one
    reference in this document that none of our code produced.
    """
    import re
    src = HERE.parent / "docs" / "ads2014_tables.txt"
    if not src.exists():
        return {}, {}
    t = src.read_text().replace("\n", "")
    # Tables 1 (put) and 2 (call) sit between the caption of Table 1 and Table 3;
    # slicing keeps the later American tables, which share the row format, out.
    i1 = t.find("Table1Convergence")
    i3 = t.find("Table3", i1)
    if i1 < 0 or i3 < 0:
        return {}, {}
    seg = t[i1:i3]
    pat = re.compile(
        r"(90|95|100|105|110)(0\.[234])(0:0833|0:25|0:5)"
        r"(\d+:\d{4})(\d+:\d{4})(\d+:\d{4})(\d+:\d{4})")
    tlab = {"0:0833": "1m", "0:25": "3m", "0:5": "6m"}
    ms = list(pat.finditer(seg))
    if len(ms) != 90:                    # 45 puts then 45 calls; refuse a partial parse
        raise RuntimeError(f"ADS table parse found {len(ms)} rows, expected 90")
    put, call = {}, {}
    for i, m in enumerate(ms):
        key = (float(m.group(1)), float(m.group(2)), tlab[m.group(3)])
        (put if i < 45 else call)[key] = float(m.group(7).replace(":", "."))
    return put, call


# ---------------------------------------------------------------------- helpers
def fmt(x, nd=4):
    return f"{x:.{nd}f}"


def main() -> None:
    if not BIN.exists():
        raise SystemExit(f"missing binary {BIN}; build it first (see README)")

    meta = {"K": K, "r": R, "eta": ETA, "kappa": KAPPA, "theta": THETA,
            "rho": RHO, "steps": STEPS, "paths_eu": PATHS_EU,
            "paths_reg": PATHS_REG, "paths_val": PATHS_VAL,
            "ex_stride": EX_STRIDE, "riccati_steps": RIC_STEPS,
            "H_values": list(HS), "seed": SEED,
            "reference": "arXiv:1205.3555 (ADS 2014), Beliaeva-Nawalkha parameters"}

    with Progress("rheston-tables", total_phases=5, meta=meta) as P:
        md: list[str] = []
        md.append("# Rough Heston on the Akyildirim--Dolinsky--Soner grid\n")
        md.append("Parameter sets of Beliaeva and Nawalkha as used in "
                  "arXiv:1205.3555 (Ann. Appl. Probab. 24(5) 2176--2205, 2014), "
                  "with the roughness added on top.\n")
        md.append(f"Fixed: `K = {K:.0f}`, `r = {R}`, vol of vol `eta = {ETA}`, "
                  f"mean reversion `kappa = {KAPPA}`, long-run variance "
                  f"`theta = {THETA}`, `rho = {RHO}` (European tables).\n")
        md.append(f"Monte-Carlo: Euler--Volterra, `{STEPS}` steps, "
                  f"`{PATHS_EU:,}` antithetic paths (European); "
                  f"`{PATHS_REG:,}` regression + `{PATHS_VAL:,}` valuation paths "
                  f"and `{STEPS // EX_STRIDE}` exercise dates (American). "
                  f"Fourier: fractional Riccati, `{RIC_STEPS}` steps, "
                  f"Gauss--Legendre `nq = {NQ}` on `(0, {UMAX:.0f})`.\n")

        # ---------------------------------------------------------- phase 1
        pub_put, pub_call = ads_published()
        keys = [(s0, v, tl) for tl in MONTHS for v in VOLS for s0 in S0S]
        P.phase("1. anchor at H=0.5 against the analytical Heston prices ADS publish",
                total=len(keys) + 1)

        cfgs = []
        for (s0, v, tl) in keys:
            cfgs.append(cfg(f"A_{tl}_{v}_{int(s0)}", 0.50, v * v, MONTHS[tl], s0,
                            american=False))
        t0 = time.time()
        mc = mc_batch(cfgs)
        t_mc_anchor = time.time() - t0

        rows, errs_f, errs_m = [], [], []
        t_f = 0.0
        for i, (s0, v, tl) in enumerate(keys):
            cid = f"A_{tl}_{v}_{int(s0)}"
            tt = time.time()
            fp = fourier_put(0.50, v * v, MONTHS[tl], s0)
            fc = fourier_call(0.50, v * v, MONTHS[tl], s0)
            t_f += time.time() - tt
            m = mc[cid]
            ref = pub_put.get((s0, v, tl))
            refc = pub_call.get((s0, v, tl))
            if ref is not None:
                errs_f.append(fp - ref)
                errs_m.append(m["eu_put"] - ref)
            rows.append((s0, v, tl, ref, fp, m["eu_put"], m["eu_put_se"],
                         refc, fc, m["eu_call"], m["eu_call_se"]))
            P.tick(i + 1)

        md.append("\n## 1. The anchor: H = 0.5 is classical Heston\n")
        md.append("At `H = 0.5` the kernel is constant and the model IS the Heston "
                  "model, so ADS's own **published analytical column** is an "
                  "external reference that none of our code produced. Both our "
                  "engines are measured against it.\n")
        md.append("| S0 | sqrt(V0) | T | ADS analytic put | Fourier put | err | "
                  "MC put | s.e. | err | ADS analytic call | Fourier call | "
                  "MC call | s.e. |")
        md.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for (s0, v, tl, ref, fp, mp, mse, refc, fc, mc_, mcse) in rows:
            rs = fmt(ref) if ref is not None else "--"
            rcs = fmt(refc) if refc is not None else "--"
            ef = fmt(fp - ref, 4) if ref is not None else "--"
            em = fmt(mp - ref, 4) if ref is not None else "--"
            md.append(f"| {s0:.0f} | {v} | {tl} | {rs} | {fmt(fp)} | {ef} | "
                      f"{fmt(mp)} | {mse:.1e} | {em} | {rcs} | {fmt(fc)} | "
                      f"{fmt(mc_)} | {mcse:.1e} |")
        if errs_f:
            af = np.abs(errs_f); am = np.abs(errs_m)
            md.append(f"\nAgainst the published analytical column, over "
                      f"{len(errs_f)} European puts: Fourier mean |error| "
                      f"**{af.mean():.5f}**, max **{af.max():.5f}**; "
                      f"Monte-Carlo mean |error| **{am.mean():.5f}**, max "
                      f"**{am.max():.5f}**.\n")
            P.result("anchor_fourier_mean_abs_err", float(af.mean()))
            P.result("anchor_fourier_max_abs_err", float(af.max()))
            P.result("anchor_mc_mean_abs_err", float(am.mean()))
            P.result("anchor_mc_max_abs_err", float(am.max()))
        md.append(f"Wall time: Fourier {t_f:.2f}s for {len(keys)} puts+calls, "
                  f"Monte-Carlo {t_mc_anchor:.2f}s for {len(keys)} parameter sets.\n")
        P.tick(len(keys) + 1)

        # ---------------------------------------------------------- phase 2
        P.phase("2. European put and call across H: Fourier vs Monte-Carlo",
                total=len(HS) * len(keys))
        eu_tables: dict[float, list] = {}
        eu_time: dict[float, tuple[float, float]] = {}
        done = 0
        for H in HS:
            cfgs = [cfg(f"E{H}_{tl}_{v}_{int(s0)}", H, v * v, MONTHS[tl], s0,
                        american=False) for (s0, v, tl) in keys]
            t0 = time.time()
            mc = mc_batch(cfgs)
            tm = time.time() - t0
            rows = []
            tf = 0.0
            for (s0, v, tl) in keys:
                cid = f"E{H}_{tl}_{v}_{int(s0)}"
                tt = time.time()
                fp = fourier_put(H, v * v, MONTHS[tl], s0)
                fc = fourier_call(H, v * v, MONTHS[tl], s0)
                tf += time.time() - tt
                m = mc[cid]
                rows.append((s0, v, tl, fp, m["eu_put"], m["eu_put_se"],
                             fc, m["eu_call"], m["eu_call_se"], m["neg_hits"]))
                done += 1
                P.tick(done)
            eu_tables[H] = rows
            eu_time[H] = (tf, tm)

        md.append("\n## 2. European options across the roughness\n")
        md.append("`Fourier` is the reference (no Monte-Carlo in it). `err` is "
                  "`MC - Fourier`; compare it with `s.e.` -- an error inside two "
                  "standard errors is not evidence of a discretisation bias, an "
                  "error outside it is.\n")
        for H in HS:
            tag = " (= classical Heston)" if H == 0.50 else ""
            md.append(f"\n### 2.{HS.index(H) + 1} H = {H}{tag}\n")
            md.append("| S0 | sqrt(V0) | T | Fourier put | MC put | s.e. | err | "
                      "err/s.e. | Fourier call | MC call | s.e. | err | err/s.e. |")
            md.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
            for (s0, v, tl, fp, mp, mse, fc, mc_, mcse, neg) in eu_tables[H]:
                ep, ec = mp - fp, mc_ - fc
                md.append(f"| {s0:.0f} | {v} | {tl} | {fmt(fp)} | {fmt(mp)} | "
                          f"{mse:.1e} | {ep:+.4f} | {ep / mse:+.1f} | {fmt(fc)} | "
                          f"{fmt(mc_)} | {mcse:.1e} | {ec:+.4f} | {ec / mcse:+.1f} |")
            tf, tm = eu_time[H]
            e = np.array([r[4] - r[3] for r in eu_tables[H]])
            s = np.array([r[5] for r in eu_tables[H]])
            md.append(f"\nPut error: mean `{e.mean():+.4f}`, mean |error| "
                      f"`{np.abs(e).mean():.4f}`, max |error| "
                      f"`{np.abs(e).max():.4f}`, mean |error|/s.e. "
                      f"`{np.abs(e / s).mean():.1f}`. "
                      f"Time: Fourier {tf:.2f}s, MC {tm:.2f}s "
                      f"({tm / len(keys):.3f}s per parameter set).\n")
            P.result(f"eu_put_mean_err_H{H}", float(e.mean()))
            P.result(f"eu_put_max_abs_err_H{H}", float(np.abs(e).max()))
            P.result(f"eu_time_fourier_H{H}", tf)
            P.result(f"eu_time_mc_H{H}", tm)

        # ---------------------------------------------------------- phase 3+4
        am_keys = [(s0, v, tl, rr) for rr in RHOS_AM for tl in MONTHS
                   for v in VOLS_AM for s0 in S0S_AM]
        P.phase("3-4. American put by LSM, and the American-call control",
                total=len(HS) * len(am_keys))
        am_tables: dict[float, list] = {}
        am_time: dict[float, tuple[float, float]] = {}
        done = 0
        for H in HS:
            cfgs = [cfg(f"M{H}_{tl}_{v}_{int(s0)}_{rr}", H, v * v, MONTHS[tl], s0,
                        rho=rr) for (s0, v, tl, rr) in am_keys]
            t0 = time.time()
            mc = mc_batch(cfgs)
            tm = time.time() - t0
            rows = []
            tf = 0.0
            for (s0, v, tl, rr) in am_keys:
                cid = f"M{H}_{tl}_{v}_{int(s0)}_{rr}"
                tt = time.time()
                fp = fourier_put(H, v * v, MONTHS[tl], s0, rho=rr)
                fc = fourier_call(H, v * v, MONTHS[tl], s0, rho=rr)
                tf += time.time() - tt
                m = mc[cid]
                rows.append((s0, v, tl, rr, fp, m["eu_put"], m["am_put"],
                             m["am_put_se"], m["am_put_insample"],
                             fc, m["am_call"], m["am_call_se"]))
                done += 1
                P.tick(done)
            am_tables[H] = rows
            am_time[H] = (tf, tm)

        md.append("\n## 3. American put\n")
        md.append("No closed form and no unbiased simulation of the American price "
                  "exists for this model, so there is no reference column and none "
                  "is invented. What is reported is a **lower bound**: the "
                  "Longstaff--Schwartz policy is fitted on the regression sample "
                  "and applied to a disjoint valuation sample, so the estimate "
                  "carries no in-sample look-ahead. `LSM in-sample` is the "
                  "look-ahead figure, printed only because the pair brackets where "
                  "the true value sits. `Euro` is the European put from the same "
                  "paths, a hard lower bound, and `premium` is the early-exercise "
                  "premium the policy captures.\n")
        for H in HS:
            tag = " (= classical Heston)" if H == 0.50 else ""
            md.append(f"\n### 3.{HS.index(H) + 1} H = {H}{tag}\n")
            md.append("| S0 | sqrt(V0) | T | rho | Euro put (Fourier) | "
                      "American LSM | s.e. | LSM in-sample | premium | premium % |")
            md.append("|---|---|---|---|---|---|---|---|---|---|")
            for (s0, v, tl, rr, fp, mep, ap, ase, ains, fc, ac, acse) in am_tables[H]:
                prem = ap - fp
                md.append(f"| {s0:.0f} | {v} | {tl} | {rr} | {fmt(fp)} | "
                          f"{fmt(ap)} | {ase:.1e} | {fmt(ains)} | {prem:+.4f} | "
                          f"{100 * prem / fp:+.2f}% |")
            tf, tm = am_time[H]
            pr = np.array([r[6] - r[4] for r in am_tables[H]])
            md.append(f"\nEarly-exercise premium: mean `{pr.mean():+.4f}`, max "
                      f"`{pr.max():+.4f}`, min `{pr.min():+.4f}`. "
                      f"Time: MC+LSM {tm:.2f}s for {len(am_keys)} parameter sets "
                      f"({tm / len(am_keys):.3f}s each).\n")
            P.result(f"am_premium_mean_H{H}", float(pr.mean()))
            P.result(f"am_premium_max_H{H}", float(pr.max()))
            P.result(f"am_time_mc_H{H}", tm)

        md.append("\n## 4. The American-call control\n")
        md.append("With `r >= 0` and no dividend, early exercise of a call is never "
                  "optimal, so the true American call **is** the European call, "
                  "which Fourier gives independently. The gap below is therefore "
                  "not a modelling error: it is exactly what the LSM exercise "
                  "policy and the finite exercise grid cost, measured instead of "
                  "assumed. It is the honest size of the downward bias to attach "
                  "to the American put column above.\n")
        md.append("| H | mean gap (LSM call - Fourier call) | mean gap % | "
                  "worst gap | worst gap % |")
        md.append("|---|---|---|---|---|")
        for H in HS:
            g = np.array([r[10] - r[9] for r in am_tables[H]])
            gp = np.array([100 * (r[10] - r[9]) / r[9] for r in am_tables[H]])
            w = int(np.argmin(g))
            md.append(f"| {H} | {g.mean():+.4f} | {gp.mean():+.2f}% | "
                      f"{g[w]:+.4f} | {gp[w]:+.2f}% |")
            P.result(f"lsm_call_gap_mean_H{H}", float(g.mean()))
            P.result(f"lsm_call_gap_worst_H{H}", float(g[w]))

        # ---------------------------------------------------------- phase 5
        sub = [(100.0, 0.2, "3m"), (100.0, 0.4, "6m"), (90.0, 0.2, "1m"),
               (110.0, 0.4, "3m")]
        step_grid = (100, 200, 400, 800)
        P.phase("5. the Euler--Volterra bias: refine the steps at fixed paths",
                total=len(HS) * len(step_grid))
        md.append("\n## 5. The discretisation bias of the Monte-Carlo\n")
        md.append("The standard error shrinks with the number of paths; the "
                  "Euler--Volterra bias does not. Here the paths are held fixed "
                  "and the step count refined, against the Fourier value. A column "
                  "that keeps moving in one direction as the steps double is bias, "
                  "not noise.\n")
        md.append("| H | S0 | sqrt(V0) | T | Fourier put | " +
                  " | ".join(f"MC N={n}" for n in step_grid) + " | s.e. at N=800 |")
        md.append("|---|---|---|---|---|" + "---|" * (len(step_grid) + 1))
        done = 0
        bias_rows = []
        for H in HS:
            refs = {}
            for (s0, v, tl) in sub:
                refs[(s0, v, tl)] = fourier_put(H, v * v, MONTHS[tl], s0)
            per_steps = {}
            for n in step_grid:
                cfgs = [cfg(f"B{H}_{n}_{tl}_{v}_{int(s0)}", H, v * v, MONTHS[tl],
                            s0, steps=n, american=False, paths_eu=200_000)
                        for (s0, v, tl) in sub]
                per_steps[n] = mc_batch(cfgs)
                done += 1
                P.tick(done)
            for (s0, v, tl) in sub:
                ref = refs[(s0, v, tl)]
                vals = [per_steps[n][f"B{H}_{n}_{tl}_{v}_{int(s0)}"]["eu_put"]
                        for n in step_grid]
                se8 = per_steps[step_grid[-1]][f"B{H}_{step_grid[-1]}_{tl}_{v}_{int(s0)}"]["eu_put_se"]
                md.append(f"| {H} | {s0:.0f} | {v} | {tl} | {fmt(ref)} | " +
                          " | ".join(f"{x - ref:+.4f}" for x in vals) +
                          f" | {se8:.1e} |")
                bias_rows.append((H, s0, v, tl, ref, vals, se8))
        md.append("\nEntries are `MC - Fourier`, so a column of numbers shrinking "
                  "towards zero is the discretisation error being removed.\n")
        for H in HS:
            b = [r for r in bias_rows if r[0] == H]
            first = np.mean([r[5][0] - r[4] for r in b])
            last = np.mean([r[5][-1] - r[4] for r in b])
            md.append(f"- `H = {H}`: mean signed error `{first:+.4f}` at "
                      f"`N = {step_grid[0]}`, `{last:+.4f}` at `N = {step_grid[-1]}`.")
            P.result(f"bias_first_H{H}", float(first))
            P.result(f"bias_last_H{H}", float(last))

        # ---------------------------------------------------------- timings
        md.append("\n## 6. Computation time\n")
        md.append("Machine: Apple silicon, 10 cores, the C++ engine using "
                  f"{max(1, 10 - 1)} threads; the Fourier reference is "
                  "single-threaded Python/NumPy. ADS report their own MATLAB tree "
                  "at 5.71 s per European option at `N = 200` and 13.97 s per "
                  "American put at `N = 250` on a 2010 laptop (their Table 10); "
                  "the comparison across machines and languages is indicative "
                  "only.\n")
        md.append("| H | Fourier, 45 European puts+calls | MC, 45 parameter sets "
                  "(European) | MC+LSM, 36 parameter sets (American) | "
                  "per set, European | per set, American |")
        md.append("|---|---|---|---|---|---|")
        for H in HS:
            tf, tm = eu_time[H]
            _, ta = am_time[H]
            md.append(f"| {H} | {tf:.2f}s | {tm:.2f}s | {ta:.2f}s | "
                      f"{tm / len(keys):.3f}s | {ta / len(am_keys):.3f}s |")
        md.append("")
        md.append(P.timing_table_md())
        P.write_results_md("\n".join(md))
        P.done()


if __name__ == "__main__":
    main()
