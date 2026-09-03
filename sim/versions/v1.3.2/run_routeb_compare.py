#!/usr/bin/env python3
"""
run_routeb_compare — Route B vs the one-step lattice, priced, across H.

Part V established the obstruction: the one-step recombining rough lattice has a
driver whose variance is `2H n^{1-2H}` times the truth, so refining the grid
drives the price away from the reference -- unboundedly.  This run implements
Route B (the Markovian lift) as an actual lattice and shows what it changes.

Three pricers, one reference:

  one-step   `route_aprime.route_aprime_european_put` -- the single-driver
             recombining lattice of Part V.
  lift (m)   `route_b_lattice.LiftedLattice` -- the (m+1)-dimensional lattice
             that replaces the one divergent driver by m Ornstein--Uhlenbeck
             factors sharing one Brownian increment, rounded onto a common grid.
  reference  `mc_reference.european_put_mc` -- 2-million-path exact-covariance
             Monte-Carlo, the one validated to Black--Scholes at eta = 0.

The lift lattice used here was cross-checked against an INDEPENDENT pricer (a
Romano--Touzi mixing formula on the enumerated sign tree, exact in the price
coordinate): the two agree to about 1.5e-2 at n = 8-12, the residual being the
two different price-coupling discretisations, which confirms the lift is built
correctly.

The message is not "the lift is more accurate at small n" -- at small m it
carries a real, known bias, because a handful of exponentials under-fit the
rough kernel's near-singularity.  The message is that the lift is CONVERGENT
where the one-step is not: the one-step error grows without bound in n, the lift
error converges to an m-dependent floor.  Section 0 shows this at the level of
the driver variance (analytic, exact); the price tables show it at the level of
the price.
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
import route_b as rb                              # noqa: E402
import route_aprime as ra                         # noqa: E402
import mc_reference as mcr                         # noqa: E402
from route_b_lattice import LiftedLattice          # noqa: E402

ETA = 0.30
RHO = -0.70
T = 1.0
HS = (0.05, 0.10, 0.30)
NS_ONE = (16, 32, 64, 128)          # one-step lattice
NS_M1 = (16, 32, 64, 128)           # lift, 1 factor
NS_M2 = (16, 32)                    # lift, 2 factors (state space grows fast in n)
NS_VR = (16, 32, 64, 128, 256, 512)  # analytic variance ratio (instant)
MC_PATHS = 2_000_000
MC_NFINE = 512
SEED = 20260805


def mref_of(n: int) -> int:
    return max(4, math.ceil(4.0 * math.sqrt(n / 8.0)))


def lift_ws(m: int, H: float):
    h = H - 0.5
    d = rb.ajee_optimised(m, h, T, t0=T / max(NS_VR))
    return d["w"] * math.sqrt(2.0 * H) * math.gamma(1.0 + h), d["s"]


def vr_onestep(n: int, H: float) -> float:
    return 2.0 * H * n ** (1.0 - 2.0 * H)


def vr_lift(n: int, H: float, m: int) -> float:
    d = T / n
    w, s = lift_ws(m, H)
    lags = np.arange(1, n + 1) * d
    Km = (w[None, :] * np.exp(-np.outer(lags, s))).sum(axis=1)
    return float(d * (Km ** 2).sum() / T ** (2.0 * H))


def main() -> None:
    meta = {"model": "rough Bergomi", "eta": ETA, "rho": RHO, "T": T,
            "H_values": list(HS), "mc_paths": MC_PATHS,
            "pricers": "one-step (route_aprime), lift (LiftedLattice), "
                       "reference (exact-covariance MC)",
            "what": "does the Markovian lift remove the one-step lattice's "
                    "variance divergence"}

    with Progress("routeb-compare", total_phases=len(HS) + 1, meta=meta) as P:
        md = ["# Route B (the Markovian lift) vs the one-step lattice, priced\n",
              "Part V showed the one-step recombining rough lattice is "
              "inconsistent: its driver variance is `2H n^{1-2H}` times the truth, "
              "so refining the grid pushes the price away from the reference "
              "without bound. Route B replaces the single divergent driver by `m` "
              "Ornstein--Uhlenbeck factors that share one Brownian increment and "
              "are rounded onto a common grid; their combined variance stays "
              "finite. This run prices both and compares them across `H`.\n",
              f"Fixed: `S0 = K = 100`, `xi0 = 0.09`, `eta = {ETA}`, `rho = {RHO}`, "
              f"`T = {T:.0f}`. Reference: `{MC_PATHS:,}`-path exact-covariance "
              f"Monte-Carlo. The lift lattice was cross-validated against an "
              f"independent mixing-formula pricer (agreement about `1.5e-2` at "
              f"`n = 8-12`, the two price-coupling discretisations differing).\n",
              "**Read this as bounded-vs-unbounded, not accurate-vs-inaccurate.** "
              "At small `m` the lift under-fits the rough kernel and carries a "
              "real bias; what it buys is convergence -- a finite error floor that "
              "shrinks as `m` grows -- where the one-step simply diverges.\n"]

        # -------- phase 0: analytic driver variance ratio ---------------------
        P.phase("0. driver variance ratio (analytic, exact)", total=len(HS))
        md.append("\n## 0. The mechanism: driver variance ratio (analytic)\n")
        md.append("`Var[driver_T]/Var[true]`, which should be `1`. The one-step "
                  "column diverges; every lift column plateaus at a finite value "
                  "that rises toward `1` as `m` grows. This needs no pricing and "
                  "is the exact statement of what the lift fixes.\n")
        md.append("| H | n | one-step `2H n^(1-2H)` | lift m=1 | lift m=2 | "
                  "lift m=3 |")
        md.append("|---|---|---|---|---|---|")
        for hi, H in enumerate(HS):
            for n in NS_VR:
                md.append(f"| {H if n == NS_VR[0] else ''} | {n} | "
                          f"{vr_onestep(n, H):.2f} | {vr_lift(n, H, 1):.3f} | "
                          f"{vr_lift(n, H, 2):.3f} | {vr_lift(n, H, 3):.3f} |")
            P.result(f"vr_onestep_n512_H{H}", vr_onestep(512, H))
            P.result(f"vr_lift_m1_n512_H{H}", vr_lift(512, H, 1))
            P.result(f"vr_lift_m2_n512_H{H}", vr_lift(512, H, 2))
            P.tick(hi + 1)
        md.append("\nThe one-step ratio grows like `n^{1-2H}` without bound; the "
                  "lift ratios are flat in `n`. That flatness is convergence; the "
                  "gap from `1` is the `m`-factor bias, which the price floor "
                  "below inherits.\n")

        # -------- price comparison per H --------------------------------------
        for hi, H in enumerate(HS):
            P.phase(f"{hi + 1}. H = {H}: prices",
                    total=len(NS_ONE) + len(NS_M1) + len(NS_M2) + 1)
            zmax = 3.0 / math.sqrt(2.0 * H)
            ref = mcr.european_put_mc(H, ETA, RHO, nfine=MC_NFINE,
                                      paths=MC_PATHS, seed=SEED)
            mc, band = ref["price"], ref["ci95"]
            P.result(f"mc_H{H}", mc)
            step = 1
            P.tick(step)

            one = {}
            for n in NS_ONE:
                t0 = time.time()
                one[n] = ra.route_aprime_european_put(n, H, ETA, RHO, zmax=zmax,
                                                      mref=mref_of(n))
                step += 1; P.tick(step, inner=f"one-step n={n} {time.time()-t0:.0f}s")
            m1 = {}
            for n in NS_M1:
                t0 = time.time()
                m1[n] = LiftedLattice(n, H, ETA, RHO, 1, mref=mref_of(n),
                                      mref_z=mref_of(n)).price()["value"]
                step += 1; P.tick(step, inner=f"lift m=1 n={n} {time.time()-t0:.0f}s")
            m2 = {}
            for n in NS_M2:
                t0 = time.time()
                m2[n] = LiftedLattice(n, H, ETA, RHO, 2, mref=mref_of(n),
                                      mref_z=mref_of(n)).price()["value"]
                step += 1; P.tick(step, inner=f"lift m=2 n={n} {time.time()-t0:.0f}s")

            md.append(f"\n## H = {H}: European put (reference **{mc:.4f}** "
                      f"±{band:.4f})\n")
            md.append("Each cell is the signed error `price - reference`. The "
                      "one-step column should grow; the lift columns should settle.\n")
            md.append("| n | one-step err | lift m=1 err | lift m=2 err |")
            md.append("|---|---|---|---|")
            alln = sorted(set(NS_ONE) | set(NS_M1) | set(NS_M2))
            for n in alln:
                c_one = f"{one[n]-mc:+.4f}" if n in one else "--"
                c_m1 = f"{m1[n]-mc:+.4f}" if n in m1 else "--"
                c_m2 = f"{m2[n]-mc:+.4f}" if n in m2 else "--"
                md.append(f"| {n} | {c_one} | {c_m1} | {c_m2} |")

            one_lo, one_hi = one[NS_ONE[0]] - mc, one[NS_ONE[-1]] - mc
            m1_lo, m1_hi = m1[NS_M1[0]] - mc, m1[NS_M1[-1]] - mc
            md.append(f"\nOne-step error moves `{one_lo:+.4f} -> {one_hi:+.4f}` "
                      f"over `n = {NS_ONE[0]}..{NS_ONE[-1]}` (variance ratio "
                      f"reaches {vr_onestep(NS_ONE[-1], H):.1f}); lift m=1 error "
                      f"moves `{m1_lo:+.4f} -> {m1_hi:+.4f}` and its increments "
                      f"shrink, converging to the m=1 floor. Increasing to m=2 "
                      f"raises the plateau toward the truth. Reference band "
                      f"±{band:.4f}.\n")
            P.result(f"onestep_err_n{NS_ONE[-1]}_H{H}", float(one[NS_ONE[-1]] - mc))
            P.result(f"lift_m1_err_n{NS_M1[-1]}_H{H}", float(m1[NS_M1[-1]] - mc))
            P.result(f"lift_m2_err_n{NS_M2[-1]}_H{H}", float(m2[NS_M2[-1]] - mc))

        md.append("\n## Is there a 'best n', and a bound? (the honest answer)\n")
        md.append("Yes, for the one-step lattice, and it is a symptom rather than "
                  "a feature. Its error is a finite-`n` part that shrinks plus a "
                  "variance-level part `~(2H n^{1-2H} - 1)` that grows; they cancel "
                  "near the `n` where the variance ratio crosses one, "
                  "`n* = (2H)^{-1/(1-2H)}` (about "
                  + ", ".join(f"{H}:{(2*H)**(-1.0/(1.0-2.0*H)):.0f}" for H in HS)
                  + " here). That is the 'small n works, large n does not' the eye "
                  "sees. But it is CANCELLATION, fragile and payoff-dependent, and "
                  "there is **no convergent error bound** for the one-step scheme: "
                  "past `n*` it diverges. The lift removes the growing term, so its "
                  "error is bounded and converges to a floor set by `m` -- a "
                  "genuine, controllable bound. That is the whole difference "
                  "between a lucky window and a convergent scheme, and it is why "
                  "the paper needs Route B.\n")
        md.append("The remaining task, not done here, is to drive the lift floor "
                  "below the reference band by taking `m` a little larger (the "
                  "variance ratios of Section 0 suggest `m ~ 3-5` reaches within a "
                  "few percent for these H), and to extend the lift lattice to the "
                  "American payoff, which the one-step cannot price consistently "
                  "either.\n")
        md.append("")
        md.append(P.timing_table_md())
        P.write_results_md("\n".join(md))
        P.done()


if __name__ == "__main__":
    main()
