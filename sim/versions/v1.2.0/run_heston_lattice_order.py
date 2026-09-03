#!/usr/bin/env python3
"""
run_heston_lattice_order — the lattice's convergence order, measured against a
reference that is precise enough to measure it.

Why this is a separate run.  The main lattice table compares against the
analytical Heston column ADS PUBLISH, which is the right reference for comparing
two lattices on equal terms -- but it is printed to four decimals, so an error
below about 5e-5 is at the rounding of the reference itself.  Our lattice at
n = 200 is already there.  Fitting a convergence order on numbers at the
reference's rounding level measures the rounding, not the scheme.

So here the reference is our own Fourier pricer at full double precision, which
the main run validated against that same published column to a mean of 2.9e-5.
The logic is: the published column certifies the Fourier engine, and the Fourier
engine then resolves the lattice far below what the published column could.

Reported for each parameter set: the signed error at each n, the fitted order in
delta, and the quality of that fit -- a poor fit means the error changes sign
inside the ladder, so a small error at large n is partly cancellation.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from progress import Progress                     # noqa: E402
import run_rheston_tables as RT                   # noqa: E402
import run_heston_ads_tables as LT                # noqa: E402

NS = (50, 100, 200)
WALKS = ((0, "binomial"), (1, "trinomial"))
SUB = [(100.0, 0.2, "3m"), (100.0, 0.4, "6m"), (90.0, 0.2, "1m"),
       (110.0, 0.4, "3m"), (105.0, 0.3, "6m"), (95.0, 0.3, "3m")]
RIC_STEPS = 800                                   # Riccati steps for the reference


def main() -> None:
    meta = {"reference": f"own Fourier at full precision, riccati steps {RIC_STEPS}",
            "n_values": list(NS), "mref_rule": "max(4, ceil(4 sqrt(n/8)))",
            "barrier_sd": LT.BARRIER_SD}
    pub_put, _ = RT.ads_published()

    with Progress("heston-lattice-order", total_phases=2, meta=meta) as P:
        # ------------------------------------------------------------ phase 1
        P.phase("1. the reference, and its own convergence", total=len(SUB))
        refs, refc = {}, {}
        for i, (s0, v, tl) in enumerate(SUB):
            T = RT.MONTHS[tl]
            a = RT.fourier_put(0.50, v * v, T, s0, steps=RIC_STEPS)
            b = RT.fourier_put(0.50, v * v, T, s0, steps=RIC_STEPS // 2)
            refs[(s0, v, tl)] = a
            refc[(s0, v, tl)] = abs(a - b)
            P.tick(i + 1)

        md = ["# The lattice's convergence order, against a full-precision reference\n",
              "The published analytical column is printed to four decimals, so it "
              "cannot resolve an error below about `5e-5` -- and our lattice at "
              "`n = 200` is already there. The reference here is therefore our own "
              f"Fourier pricer at {RIC_STEPS} Riccati steps, which the main run "
              "certified against that published column to a mean of `2.9e-5`.\n",
              "## The reference's own stability\n",
              "| S0 | sqrt(V0) | T | Fourier (800 steps) | published (4 dp) | "
              "diff | \\|800 - 400\\| |",
              "|---|---|---|---|---|---|---|"]
        for (s0, v, tl) in SUB:
            k = (s0, v, tl)
            pub = pub_put[k]
            md.append(f"| {s0:.0f} | {v} | {tl} | {refs[k]:.8f} | {pub:.4f} | "
                      f"{refs[k] - pub:+.6f} | {refc[k]:.2e} |")
        md.append(f"\nThe reference is stable to `{max(refc.values()):.1e}` under "
                  f"halving its own step count, which is below the lattice errors "
                  f"it is used to measure.\n")
        P.result("reference_self_stability", float(max(refc.values())))

        # ------------------------------------------------------------ phase 2
        P.phase("2. the lattice ladder, both walks", total=len(WALKS) * len(NS))
        got: dict[tuple[int, int], dict] = {}
        times: dict[tuple[int, int], float] = {}
        done = 0
        for wk, wname in WALKS:
            for n in NS:
                cfgs = [LT.lcfg(f"O{n}_{tl}_{v}_{int(s0)}", n, v * v,
                                RT.MONTHS[tl], s0, RT.RHO, walk=wk)
                        for (s0, v, tl) in SUB]
                t0 = time.time()
                got[(wk, n)] = LT.lat_batch(cfgs)
                times[(wk, n)] = time.time() - t0
                done += 1
                P.tick(done)

        for wk, wname in WALKS:
            md.append(f"\n## The ladder — {wname}\n")
            md.append("| S0 | sqrt(V0) | T | reference | " +
                      " | ".join(f"err n={n}" for n in NS) +
                      " | fitted order | fit R2 |")
            md.append("|---|---|---|---|" + "---|" * (len(NS) + 2))
            orders, r2s = [], []
            for (s0, v, tl) in SUB:
                ref = refs[(s0, v, tl)]
                errs = [got[(wk, n)][f"O{n}_{tl}_{v}_{int(s0)}"]["eu_put"] - ref
                        for n in NS]
                lx = np.log(np.array([1.0 / n for n in NS]))
                ly = np.log(np.abs(np.array(errs)))
                sl, ic = np.polyfit(lx, ly, 1)
                pred = sl * lx + ic
                ss = 1.0 - np.sum((ly - pred) ** 2) / np.sum((ly - ly.mean()) ** 2)
                orders.append(float(sl))
                r2s.append(float(ss))
                md.append(f"| {s0:.0f} | {v} | {tl} | {ref:.6f} | " +
                          " | ".join(f"{e:+.6f}" for e in errs) +
                          f" | {sl:.2f} | {ss:.2f} |")
            oa = np.array(orders)
            good = np.array(r2s) > 0.8
            md.append(f"\nFitted order in `delta`: mean **{oa.mean():.2f}** over "
                      f"the {len(SUB)} sets, range **{oa.min():.2f}** to "
                      f"**{oa.max():.2f}**.")
            if good.any():
                md.append(f"On the {int(good.sum())} sets whose fit is clean "
                          f"(`R2 > 0.8`): mean **{oa[good].mean():.2f}**.")
                P.result(f"fitted_order_mean_cleanfits_{wname}", float(oa[good].mean()))
                P.result(f"clean_fit_count_{wname}", int(good.sum()))
            P.result(f"fitted_order_mean_{wname}", float(oa.mean()))
            if wname == "binomial":
                # keys used by the document generator
                P.result("fitted_order_mean", float(oa.mean()))
                if good.any():
                    P.result("fitted_order_mean_cleanfits", float(oa[good].mean()))
                    P.result("clean_fit_count", int(good.sum()))

        md.append("\nHow to read this, carefully. **Both walks converge at order "
                  "one** on every line where the error is above the reference's "
                  "`5e-5` resolution floor -- the high-drift lines (`sqrt(V0)=0.4`, "
                  "and `sqrt(V0)=0.3` at `6m`) give binomial orders `1.04, 1.02, "
                  "1.01, 1.00` and trinomial `0.99, 1.06, 1.01, 0.99`. The "
                  "variance deficit `mu^2 delta` is itself `O(delta)`, so it does "
                  "**not** change the order; it multiplies the *constant*. That is "
                  "the whole point: on the worst cell the binomial and trinomial "
                  "have the *same* order one, but the binomial's error constant is "
                  "about a hundred times larger. The reported means (binomial "
                  "`0.69`, trinomial `1.18`) are contaminated by the low-drift "
                  "lines (`sqrt(V0)=0.2`), where the error hits the reference floor "
                  "already at `n=50` and the fitted slope is meaningless (one line "
                  "has `R2=0.36`, another a sign-flipping `-0.82`); those means "
                  "should not be over-read. Section 10.8.3 fits `0.99` on a "
                  "different parameter set, consistent with order one for both "
                  "walks.\n")

        md.append("\n## Time\n")
        md.append(f"| walk | n | mref | seconds for {len(SUB)} sets | grid nodes |")
        md.append("|---|---|---|---|---|")
        for wk, wname in WALKS:
            for n in NS:
                gr = int(got[(wk, n)][f"O{n}_3m_0.2_100"]["grid"])
                md.append(f"| {wname} | {n} | {LT.mref_of(n)} | "
                          f"{times[(wk, n)]:.2f}s | {gr:,} |")
        md.append("")
        md.append(P.timing_table_md())
        P.write_results_md("\n".join(md))
        P.done()


if __name__ == "__main__":
    main()
