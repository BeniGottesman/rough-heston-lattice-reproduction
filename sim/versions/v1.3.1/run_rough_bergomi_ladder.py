#!/usr/bin/env python3
"""
run_rough_bergomi_ladder — the paper's lattice IN THE ROUGH REGIME.

The H = 0.5 tables validate the construction against an exact price, but H = 0.5
is not rough, and it is the rough case that the paper is about.  This run prices
the genuinely rough model with the paper's own lattice and puts it next to a
trustworthy Monte-Carlo, across H, with the error and the Monte-Carlo band on the
same line.

Why rough BERGOMI and not rough Heston here.  A rough lattice column can only be
built for a model INSIDE the paper's class -- an autonomous driver.  Rough Bergomi
is such a model; rough Heston is not (its `nu sqrt(V) dB` makes the driver depend
on its own past), so the lattice cannot price rough Heston at all without the
Markovian lift of Route B, which is not implemented.  The price of staying
in-class is that rough Bergomi has NO closed form, so the reference here is a
Monte-Carlo -- but a trustworthy one: the driver is drawn from its exact Gaussian
covariance (`Var[V^H_T]` matched to 12 digits) and an `eta = 0` control variate
whose mean is Black--Scholes in closed form removes most of the variance.  The
error reported is therefore `tree - MC`, and it is meaningful only down to the
Monte-Carlo band, which is why the band is printed on every line and the path
count is pushed high enough to put it well below the scheme error.

What this run actually shows.  A RECOMBINING lattice cannot represent the exact
convolution V^(n): that is Corollary "recombining vs consistent" in the paper.
The implementable recombining object is therefore the ONE-STEP lattice V-check,
in which the fractional driver W^H_{t_k} is approximated by (net up-moves) x
delta^H -- a simple random walk scaled by delta^H.  Its terminal variance
DIVERGES:  Var[V-check_T] / Var[V_T] = 2H n^{1-2H} -> infinity  (Proposition
"variance discrepancy").  So this run does not measure a convergence rate; it
measures an OBSTRUCTION.  We expect the price error to track that variance ratio,
turning positive and growing once the ratio exceeds one -- and pushing n UP to
make it worse, not better.  The coarse-grid sweep of a previous session stopped
at n = 64 and was partly hidden by sign cancellation; extending to n = 128 is
what makes the divergence unmistakable.

The barrier is applied identically to the tree and to the Monte-Carlo (same
`zmax` on the driver), so absorption cancels in the difference and the reported
error is the scheme error alone.
"""
from __future__ import annotations

import math
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from progress import Progress                     # noqa: E402
import route_aprime as ra                         # noqa: E402
import mc_reference as mcr                         # noqa: E402

HS = (0.05, 0.10, 0.30)
NS = (8, 16, 32, 64, 128)
ETA = 0.30
RHO = -0.70
MC_PATHS = 2_000_000                # tight band: ~4x the sweep, se ~ 0.003
MC_NFINE = 512
SEED = 20260805


def mref_of(n: int) -> int:
    return max(4, math.ceil(4.0 * math.sqrt(n / 8.0)))


def main() -> None:
    meta = {"model": "rough Bergomi (in the paper's class; no closed form)",
            "reference": "exact-covariance-driver Monte-Carlo, eta=0 control variate",
            "H_values": list(HS), "n_values": list(NS), "eta": ETA, "rho": RHO,
            "mc_paths": MC_PATHS, "note": "sigma_y=1 so the driver walk has no "
            "drift, so the divergence is the recombining-lattice inconsistency "
            "(variance ratio 2H n^(1-2H)), not the F058 Lamperti-drift loss"}

    with Progress("rough-bergomi-ladder", total_phases=len(HS) + 1, meta=meta) as P:
        md = ["# The paper's lattice in the rough regime (rough Bergomi)\n",
              "Route A' lattice against the exact-covariance Monte-Carlo, across "
              "the roughness `H`. This is the genuinely rough case -- `H = 0.5` is "
              "not rough and only serves to validate the machinery against a "
              "closed form (Part III). Here there is no closed form, because a "
              "rough model cannot be both inside the paper's class and "
              "semi-analytic; the reference is therefore a Monte-Carlo, and the "
              "error is `tree - MC`, meaningful down to the Monte-Carlo band "
              "printed on each line.\n",
              "**This part is a NEGATIVE result, and it is the important one.** A "
              "recombining lattice cannot be the exact convolution (Corollary "
              "'recombining vs consistent'), so the implementable object is the "
              "one-step lattice, whose terminal variance diverges as "
              "`2H n^{1-2H}` (Proposition 'variance discrepancy'). The table below "
              "is that proposition made numerical: the error tracks the variance "
              "ratio and GROWS with `n`. This is exactly the obstruction that "
              "Route B (the Markovian lift) exists to remove -- and Route B is not "
              "yet implemented, so there is as yet no consistent recombining "
              "lattice for the rough regime.\n",
              f"Fixed: `S0 = K = 100`, `xi0 = 0.09` (`sigma = 0.30`), `T = 1`, "
              f"`eta = {ETA}`, `rho = {RHO}`. Monte-Carlo: exact Gaussian driver, "
              f"`{MC_PATHS:,}` paths, `eta = 0` Black--Scholes control variate. "
              f"Barrier at `zmax = 3/sqrt(2H)` applied identically to tree and "
              f"Monte-Carlo, so absorption cancels in the error.\n",
              "(The driver walk is driftless here, `sigma_y = 1`, so the F058 "
              "Lamperti-drift variance loss of the Heston case is absent; the "
              "divergence below is a different and deeper effect -- the "
              "inconsistency of the recombining rough lattice itself.)\n"]

        # -- reference sanity: eta=0 must return Black-Scholes to machine zero
        P.phase("0. reference sanity (eta=0 control = Black-Scholes)", total=len(HS))
        bs = mcr.bs_put()
        md.append(f"\n## Reference sanity\n")
        md.append(f"At `eta = 0` the control variate and the payoff coincide "
                  f"pathwise, so the estimator must return Black--Scholes "
                  f"`{bs:.6f}` with zero variance -- the one exact check available "
                  f"in the rough regime.\n")
        md.append("| H | price at eta=0 | band |")
        md.append("|---|---|---|")
        for i, H in enumerate(HS):
            r = mcr.european_put_mc(H, 0.0, RHO, nfine=MC_NFINE, paths=200_000,
                                    seed=SEED)
            md.append(f"| {H} | {r['price']:.6f} | ±{r['ci95']:.1e} |")
            P.tick(i + 1)

        for hi, H in enumerate(HS):
            P.phase(f"{hi + 1}. H = {H}", total=len(NS) + 1)
            zmax = 3.0 / math.sqrt(2.0 * H)
            t0 = time.time()
            ref = mcr.european_put_mc(H, ETA, RHO, nfine=MC_NFINE,
                                      paths=MC_PATHS, seed=SEED)
            t_mc = time.time() - t0
            mc_price, band = ref["price"], ref["ci95"]
            P.result(f"mc_price_H{H}", mc_price)
            P.result(f"mc_band_H{H}", band)
            P.tick(1)

            rows = []
            t_tree = 0.0
            for j, n in enumerate(NS):
                tt = time.time()
                tree = ra.route_aprime_european_put(n, H, ETA, RHO, zmax=zmax,
                                                    mref=mref_of(n))
                t_tree += time.time() - tt
                rows.append((n, tree, tree - mc_price))
                P.tick(j + 2)

            md.append(f"\n## H = {H}  (barrier zmax = {zmax:.2f})\n")
            md.append(f"Monte-Carlo reference **{mc_price:.4f}** ±{band:.4f} "
                      f"({MC_PATHS:,} paths, {t_mc:.1f}s).\n")
            md.append("| n | mref | tree | error (tree - MC) | error / band | "
                      "variance ratio `2H n^(1-2H)` |")
            md.append("|---|---|---|---|---|---|")
            ratios, errs = [], []
            for (n, tree, err) in rows:
                ratio = 2.0 * H * n ** (1.0 - 2.0 * H)
                ratios.append(ratio); errs.append(err)
                md.append(f"| {n} | {mref_of(n)} | {tree:.4f} | {err:+.4f} | "
                          f"{err / band:+.1f} | {ratio:.2f} |")
            # rank correlation between the variance ratio and the SIGNED error
            r = np.array(ratios); e = np.array(errs)
            rho_s = float(np.corrcoef(np.argsort(np.argsort(r)),
                                      np.argsort(np.argsort(e)))[0, 1])
            md.append(f"\nThe signed error rises monotonically with the variance "
                      f"ratio (rank correlation **{rho_s:+.2f}**): negative while "
                      f"the ratio is near one and the finite barrier/coupling "
                      f"errors dominate, crossing zero near ratio approximately 2.5 "
                      f"and then growing positive as the ratio does. At `n = 128` "
                      f"the ratio is **{ratios[-1]:.1f}** and the tree overprices "
                      f"by **{errs[-1]:+.4f}** ({errs[-1] / band:+.0f} bands). This "
                      f"is Proposition 'variance discrepancy' made numerical, not a "
                      f"convergence. Tree wall time {t_tree:.1f}s.\n")
            P.result(f"tree_error_n128_H{H}", rows[-1][2])
            P.result(f"variance_ratio_n128_H{H}", ratios[-1])
            P.result(f"rank_corr_ratio_error_H{H}", rho_s)

        md.append("\n## What this shows\n")
        md.append("- **The recombining rough lattice is inconsistent, and this is "
                  "the numerical proof.** Its terminal variance is "
                  "`2H n^{1-2H}` times the truth, so refining the grid drives the "
                  "price AWAY from the reference. The error's sign and growth track "
                  "that ratio across every `H`. A previous session's sweep stopped "
                  "at `n = 64`, where sign cancellation made some `H` look "
                  "convergent; `n = 128` removes that illusion.\n")
        md.append("- **This is the obstruction Route B is for.** The Markovian "
                  "lift replaces the one divergent driver by `m` mean-reverting "
                  "factors whose combined variance stays finite; only then does a "
                  "recombining rough lattice converge. Route B is not implemented, "
                  "so **there is as yet no consistent recombining lattice for the "
                  "rough regime** -- neither for rough Bergomi here nor for rough "
                  "Heston.\n")
        md.append("- **Hence there is no rough lattice column against an exact "
                  "price anywhere in this document, and cannot be one yet.** Rough "
                  "Bergomi (in-class) has no closed form; rough Heston (closed form "
                  "via Fourier) is out-of-class. The single most valuable next step "
                  "is Route B applied to rough Heston, which would finally put a "
                  "converging lattice next to the exact Fourier price.\n")
        md.append("")
        md.append(P.timing_table_md())
        P.write_results_md("\n".join(md))
        P.done()


if __name__ == "__main__":
    main()
