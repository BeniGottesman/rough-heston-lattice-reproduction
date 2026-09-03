#!/usr/bin/env python3
"""
run_rheston_american_anchor — how good is our American estimator, really?

The American table of the rough Heston document has no reference column, because
for the rough model none exists.  That leaves the reader with no way to judge the
Longstaff--Schwartz numbers.  This run supplies the missing calibration by
running the SAME estimator at H = 0.5, where the model is classical Heston and
where two independent American put prices are PUBLISHED:

  Tree, N = 250     Akyildirim--Dolinsky--Soner, arXiv:1205.3555, Table 3
  Control variate   Beliaeva--Nawalkha, as quoted in the same table
                    (CV = Tree American + (closed-form Euro - Tree Euro))

Both are third-party numbers on exactly the grid we use: K = 100, r = 0.05,
eta = 0.1, kappa = 3.0, theta = 0.04, S0 in {90, 100, 110},
rho in {-0.1, -0.7}, sqrt(V0) in {0.2, 0.4}, T in {1, 3, 6} months.

Whatever gap our LSM shows here is the gap to carry over to the rough tables: it
is the cost of the exercise policy and of the finite exercise grid, and it does
not go away when H drops.  The exercise grid is set to 200 dates so that the
comparison against a 250-step tree is close to like for like.
"""
from __future__ import annotations

import re
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from progress import Progress            # noqa: E402
import run_rheston_tables as RT          # noqa: E402

STEPS = 200
EX_STRIDE = 1                            # 200 exercise dates
PATHS_EU = 400_000
PATHS_REG = 50_000
PATHS_VAL = 100_000
SEED = 20260805


def ads_american() -> dict:
    """Parse ADS Table 3: {(S0, rho, sqrt(V0), Tlabel): (tree_N250, cv_N200)}."""
    src = HERE.parent / "docs" / "ads2014_tables.txt"
    t = src.read_text().replace("\n", "")
    i = t.find("Table3", t.find("Table2Convergence"))
    j = t.find("Table4", i)
    seg = t[i:j]
    pat = re.compile(
        r"(90|100|110)\x000:([17])(0\.[24])(0:0833|0:25|0:5)"
        r"(\d+:\d{4})(\d+:\d{4})")
    tlab = {"0:0833": "1m", "0:25": "3m", "0:5": "6m"}
    out = {}
    for m in pat.finditer(seg):
        key = (float(m.group(1)), -float("0." + m.group(2)),
               float(m.group(3)), tlab[m.group(4)])
        out[key] = (float(m.group(5).replace(":", ".")),
                    float(m.group(6).replace(":", ".")))
    if len(out) != 36:
        raise RuntimeError(f"ADS Table 3 parse found {len(out)} rows, expected 36")
    return out


def main() -> None:
    pub = ads_american()
    keys = sorted(pub, key=lambda k: (k[3], k[2], k[1], k[0]))

    meta = {"purpose": "calibrate the LSM American estimator at H=0.5 against "
                       "published Heston American put prices",
            "steps": STEPS, "ex_dates": STEPS // EX_STRIDE,
            "paths_eu": PATHS_EU, "paths_reg": PATHS_REG,
            "paths_val": PATHS_VAL,
            "reference": "ADS2014 Table 3 (tree N=250) and Beliaeva-Nawalkha CV"}

    with Progress("rheston-american-anchor", total_phases=2, meta=meta) as P:
        P.phase("1. American put at H=0.5 against two published references",
                total=len(keys) + 1)
        cfgs = []
        for (s0, rr, v, tl) in keys:
            cid = f"AN_{tl}_{v}_{int(s0)}_{rr}"
            c = RT.cfg(cid, 0.50, v * v, RT.MONTHS[tl], s0, rho=rr,
                       steps=STEPS, paths_eu=PATHS_EU, seed=SEED)
            c["paths_reg"] = PATHS_REG
            c["paths_val"] = PATHS_VAL
            c["ex_stride"] = EX_STRIDE
            cfgs.append(c)
        t0 = time.time()
        mc = RT.mc_batch(cfgs)
        t_mc = time.time() - t0
        P.tick(1)

        md = ["# The American estimator, calibrated at H = 0.5\n",
              "At `H = 0.5` the rough Heston model IS the Heston model, and two "
              "independent American put prices are published for exactly this "
              "grid: the recombining tree of arXiv:1205.3555 at `N = 250`, and "
              "the Beliaeva--Nawalkha control-variate value quoted in the same "
              "table. Our Longstaff--Schwartz estimator is run unchanged and "
              "measured against both.\n",
              f"Fixed: `K = 100`, `r = 0.05`, `eta = 0.1`, `kappa = 3.0`, "
              f"`theta = 0.04`. Monte-Carlo: `{STEPS}` Euler steps, "
              f"`{STEPS // EX_STRIDE}` exercise dates, `{PATHS_REG:,}` regression "
              f"+ `{PATHS_VAL:,}` valuation paths.\n",
              "| S0 | rho | sqrt(V0) | T | ADS tree N=250 | B-N control variate | "
              "our LSM | s.e. | LSM - tree | LSM - tree % | Euro put (Fourier) | "
              "premium |",
              "|---|---|---|---|---|---|---|---|---|---|---|---|"]

        rows, gaps, gpct = [], [], []
        tf = 0.0
        for i, (s0, rr, v, tl) in enumerate(keys):
            cid = f"AN_{tl}_{v}_{int(s0)}_{rr}"
            m = mc[cid]
            tree, cv = pub[(s0, rr, v, tl)]
            tt = time.time()
            fp = RT.fourier_put(0.50, v * v, RT.MONTHS[tl], s0, rho=rr)
            tf += time.time() - tt
            g = m["am_put"] - tree
            gaps.append(g)
            gpct.append(100.0 * g / tree)
            rows.append((s0, rr, v, tl, tree, cv, m["am_put"], m["am_put_se"],
                         g, fp))
            md.append(f"| {s0:.0f} | {rr} | {v} | {tl} | {tree:.4f} | {cv:.4f} | "
                      f"{m['am_put']:.4f} | {m['am_put_se']:.1e} | {g:+.4f} | "
                      f"{100 * g / tree:+.2f}% | {fp:.4f} | "
                      f"{m['am_put'] - fp:+.4f} |")
            P.tick(i + 2)

        ga = np.array(gaps)
        gp = np.array(gpct)
        # the deep-in-the-money rows are pinned at the intrinsic value 10.0 by
        # both methods, so they flatter any estimator; report with and without
        interior = np.array([abs(r[6] - max(100.0 - r[0], 0.0)) > 1e-6 for r in rows])
        md.append("")
        md.append(f"Over all {len(rows)} parameter sets: mean gap to the tree "
                  f"**{ga.mean():+.4f}** ({gp.mean():+.2f}%), mean absolute gap "
                  f"**{np.abs(ga).mean():.4f}** ({np.abs(gp).mean():.2f}%), worst "
                  f"**{ga[np.argmax(np.abs(ga))]:+.4f}** "
                  f"({gp[np.argmax(np.abs(gp))]:+.2f}%).")
        if interior.any():
            md.append(f"Excluding the {int((~interior).sum())} rows where both "
                      f"methods sit exactly on the intrinsic value: mean gap "
                      f"**{ga[interior].mean():+.4f}** "
                      f"({gp[interior].mean():+.2f}%), mean absolute gap "
                      f"**{np.abs(ga[interior]).mean():.4f}** "
                      f"({np.abs(gp[interior]).mean():.2f}%).")
        md.append("")
        md.append(f"Wall time: Monte-Carlo + LSM {t_mc:.2f}s for {len(keys)} "
                  f"parameter sets ({t_mc / len(keys):.3f}s each); Fourier "
                  f"{tf:.2f}s for {len(keys)} European puts. ADS report 13.97s "
                  f"per American put at `N = 250` in MATLAB on a 2010 laptop.")

        P.result("am_anchor_mean_gap", float(ga.mean()))
        P.result("am_anchor_mean_abs_gap", float(np.abs(ga).mean()))
        P.result("am_anchor_mean_abs_gap_pct", float(np.abs(gp).mean()))
        P.result("am_anchor_worst_gap", float(ga[np.argmax(np.abs(ga))]))
        P.result("am_anchor_time_mc", t_mc)
        if interior.any():
            P.result("am_anchor_mean_abs_gap_interior",
                     float(np.abs(ga[interior]).mean()))

        # ------------------------------------------------------------ phase 2
        P.phase("2. does the gap come from the exercise grid?", total=3)
        md.append("\n## Where the gap comes from\n")
        md.append("A Longstaff--Schwartz price is low for two separate reasons: "
                  "the exercise decision is only as good as the regression, and "
                  "exercise is only allowed on a finite grid of dates. Refining "
                  "the number of exercise dates separates them: what moves is the "
                  "grid, what stays is the policy.\n")
        sub = [(100.0, -0.7, 0.2, "3m"), (100.0, -0.7, 0.4, "6m"),
               (110.0, -0.1, 0.4, "6m")]
        grids = (50, 100, 200)
        md.append("| S0 | rho | sqrt(V0) | T | ADS tree N=250 | " +
                  " | ".join(f"{g} dates" for g in grids) + " |")
        md.append("|---|---|---|---|---|" + "---|" * len(grids))
        per = {}
        for gi, nd in enumerate(grids):
            cfgs = []
            for (s0, rr, v, tl) in sub:
                cid = f"G{nd}_{tl}_{v}_{int(s0)}_{rr}"
                c = RT.cfg(cid, 0.50, v * v, RT.MONTHS[tl], s0, rho=rr,
                           steps=STEPS, paths_eu=200_000, seed=SEED + 1)
                c["paths_reg"] = PATHS_REG
                c["paths_val"] = PATHS_VAL
                c["ex_stride"] = STEPS // nd
                cfgs.append(c)
            per[nd] = RT.mc_batch(cfgs)
            P.tick(gi + 1)
        for (s0, rr, v, tl) in sub:
            tree, _ = pub[(s0, rr, v, tl)]
            vals = [per[nd][f"G{nd}_{tl}_{v}_{int(s0)}_{rr}"]["am_put"]
                    for nd in grids]
            md.append(f"| {s0:.0f} | {rr} | {v} | {tl} | {tree:.4f} | " +
                      " | ".join(f"{x:.4f} ({x - tree:+.4f})" for x in vals) + " |")
        md.append("\nEntries are the price and, in brackets, the gap to the "
                  "published tree value.\n")
        md.append(P.timing_table_md())
        P.write_results_md("\n".join(md))
        P.done()


if __name__ == "__main__":
    main()
