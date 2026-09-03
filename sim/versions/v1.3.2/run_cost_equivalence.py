#!/usr/bin/env python3
"""
run_cost_equivalence — how long a lattice price takes, and how many Monte-Carlo
paths (and how much time) it takes to reach the SAME accuracy.

Six phases.

  setup       the Monte-Carlo's one-off cost, timed cold (an nfine x nfine
              covariance factorisation), and its marginal path rate once warm.
  sigma       sigma_eff, the estimator's true per-path spread, measured from the
              scatter of independent replications rather than from a formula —
              the reference uses antithetic pairs and a control variate, and its
              own `stderr` treats the pairs as independent.
  reference   the "truth": a large exact-covariance Monte-Carlo per H, with its
              standard error, which bounds how precisely any bias can be known.
  lattice     the pricers timed at n = 16 and beyond: the one-step Route A'
              lattice and the Route B lift at m = 1, 2.  Best-of-3.
  equiv       the comparison: bias, N*, Monte-Carlo time at equal accuracy, the
              ratio, and the break-even accuracy — with the reference's own
              uncertainty propagated into N*.
  american    the lattice's American price and its cost, which is the thing
              Monte-Carlo cannot match cheaply.  No equivalence is computed: a
              fair one needs an LSM estimator for rough Bergomi, which this
              project has only for rough Heston.

Model parameters are the project's documented rough Bergomi convention, kept for
comparability with Parts V and VI: T = 1, S0 = K = 100, xi0 = 0.09, eta = 0.30,
rho = -0.70, mref = max(4, ceil(4 sqrt(n/8))), barrier zmax = 3/sqrt(2H).  The
reference is UNCAPPED, so the lattice's barrier truncation counts as part of its
error — that is the honest choice for "how accurate is this number".
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from progress import Progress                     # noqa: E402
import cost_equivalence as ce                     # noqa: E402
import mc_reference as mcr                        # noqa: E402

ETA = 0.30
RHO = -0.70
HS = (0.05, 0.10, 0.30)
NFINE = 512                     # the reference's fine grid
REF_PATHS = 2_000_000           # the "truth" per H
SIGMA_PATHS = 20_000            # per replication
SIGMA_REPS = 100                # replications: +-7% on a std, enough
                                # to resolve the 5-9% inflation factor
RATE_PATHS = 200_000            # for the marginal path rate
DIRECT_PATHS = 200_000          # one sample for the exact sigma_eff
REPEATS = 3                     # best-of, for the lattice timings

# (label, n, m) -- m = None is the one-step Route A' lattice
PRICERS = [
    ("one-step", 16, None), ("one-step", 32, None), ("one-step", 64, None),
    ("lift m=1", 16, 1), ("lift m=1", 32, 1), ("lift m=1", 64, 1),
    ("lift m=2", 16, 2), ("lift m=2", 32, 2),
]
AMERICAN = [("lift m=1", 16, 1), ("lift m=1", 32, 1)]


def human(t: float) -> str:
    if not math.isfinite(t):
        return "—"
    if t < 1e-3:
        return f"{t*1e6:.0f}µs"
    if t < 1.0:
        return f"{t*1e3:.0f}ms"
    if t < 60:
        return f"{t:.2f}s"
    return f"{int(t//60)}m{int(t%60):02d}s"


def hnum(x: float) -> str:
    if not math.isfinite(x):
        return "—"
    if x >= 1e6:
        return f"{x/1e6:.1f}M"
    if x >= 1e3:
        return f"{x/1e3:.1f}k"
    return f"{x:.0f}"


def main() -> None:
    t_start = time.time()
    with Progress("cost-equivalence", total_phases=6) as pr:
        timings: list[tuple[str, float]] = []

        # ------------------------------------------------------------ setup
        pr.phase("setup", total=len(HS) * 2)
        t0 = time.time()
        setup_cost, path_rate = {}, {}
        k = 0
        for H in HS:
            setup_cost[H] = ce.mc_setup_cost(H, NFINE)      # cold, before warm
            k += 1
            pr.tick(k, note=f"setup H={H}: {human(setup_cost[H])}")
            path_rate[H] = ce.mc_path_rate(H, ETA, RHO, NFINE, RATE_PATHS)
            k += 1
            pr.tick(k, note=f"rate H={H}: {path_rate[H]:,.0f}/s")
        timings.append(("setup", time.time() - t0))

        # ------------------------------------------------------------ sigma
        pr.phase("sigma", total=len(HS) * 2)
        t0 = time.time()
        sigma, direct, xval = {}, {}, {}
        k = 0
        for H in HS:
            # the EXACT estimator: sd x sqrt(1 + antithetic pair correlation)
            direct[H] = ce.sigma_eff_direct(H, ETA, RHO, NFINE, DIRECT_PATHS)
            # cross-validate the rebuilt sample against the real pricer
            chk = mcr.european_put_mc(H, ETA, RHO, nfine=NFINE,
                                      paths=DIRECT_PATHS, seed=4242)
            xval[H] = {"price_pricer": chk["price"],
                       "price_rebuilt": direct[H]["mean"],
                       "stderr_pricer": chk["stderr"],
                       "stderr_rebuilt": direct[H]["stderr_formula"]}
            k += 1
            pr.tick(k, note=f"H={H} rho_pair="
                            f"{direct[H]['control']['rho_pair']:+.3f}")
            # the weak replication estimator, kept as an independent check
            sigma[H] = ce.mc_sigma_eff(H, ETA, RHO, NFINE,
                                       SIGMA_PATHS, SIGMA_REPS)
            k += 1
            pr.tick(k, note=f"H={H} sigma_eff={sigma[H]['sigma_eff']:.3f}")
        timings.append(("sigma", time.time() - t0))

        # -------------------------------------------------------- reference
        pr.phase("reference", total=len(HS))
        t0 = time.time()
        ref = {}
        for j, H in enumerate(HS):
            tt = time.perf_counter()
            r = mcr.european_put_mc(H, ETA, RHO, nfine=NFINE,
                                    paths=REF_PATHS, seed=7)
            r["wall"] = time.perf_counter() - tt
            ref[H] = r
            pr.tick(j + 1, note=f"H={H} ref={r['price']:.4f}±{r['stderr']:.4f}")
        timings.append(("reference", time.time() - t0))

        # ---------------------------------------------------------- lattice
        pr.phase("lattice", total=len(HS) * len(PRICERS))
        t0 = time.time()
        lat: dict = {}
        k = 0
        for H in HS:
            for label, n, m in PRICERS:
                try:
                    if m is None:
                        best, med, val = ce.timed(
                            lambda: ce.price_onestep(n, H, ETA, RHO), REPEATS)
                        info = {"price": float(val), "state_space": None}
                    else:
                        best, med, val = ce.timed(
                            lambda: ce.price_lift(n, H, ETA, RHO, m), REPEATS)
                        info = val
                    lat[(H, label, n)] = {**info, "best": best, "median": med}
                except Exception as e:                  # keep the run alive
                    lat[(H, label, n)] = None
                    pr.log(f"lattice H={H} {label} n={n} failed: {e}")
                k += 1
                pr.tick(k, note=f"H={H} {label} n={n}")
        timings.append(("lattice", time.time() - t0))

        # ------------------------------------------------------------ equiv
        pr.phase("equiv", total=len(HS))
        t0 = time.time()
        eq: dict = {}
        for j, H in enumerate(HS):
            for label, n, m in PRICERS:
                L = lat.get((H, label, n))
                if L is None:
                    continue
                eq[(H, label, n)] = ce.equivalence(
                    L["price"], L["best"], ref[H]["price"], ref[H]["stderr"],
                    direct[H]["control"]["sigma_eff"], path_rate[H],
                    setup_cost[H])
                eq[(H, label, n)]["breakeven"] = ce.breakeven_accuracy(
                    direct[H]["control"]["sigma_eff"], path_rate[H], L["best"])
            pr.tick(j + 1, note=f"H={H}")
        timings.append(("equiv", time.time() - t0))

        # --------------------------------------------------------- american
        pr.phase("american", total=len(HS) * len(AMERICAN))
        t0 = time.time()
        amer: dict = {}
        k = 0
        for H in HS:
            for label, n, m in AMERICAN:
                try:
                    best, med, val = ce.timed(
                        lambda: ce.price_lift(n, H, ETA, RHO, m, american=True),
                        1)
                    eu = lat.get((H, label, n))
                    amer[(H, label, n)] = {
                        "price": val["price"], "best": best,
                        "premium": (val["price"] - eu["price"]) if eu else None,
                        "eu_time": eu["best"] if eu else None}
                except Exception as e:
                    amer[(H, label, n)] = None
                    pr.log(f"american H={H} {label} n={n} failed: {e}")
                k += 1
                pr.tick(k, note=f"H={H} {label} n={n} american")
        timings.append(("american", time.time() - t0))

        # ---------------------------------------------------------- results
        total = time.time() - t_start
        L: list[str] = []
        A = L.append
        A("# What a lattice price costs, and what the same accuracy costs by "
          "Monte-Carlo")
        A("")
        A(f"`ALGO_VERSION` {pr.meta.get('algo_version','?')} · "
          f"`code_sha` {pr.meta.get('code_sha','?')}")
        A("")
        A("A lattice price is one deterministic number carrying a **bias**. A")
        A("Monte-Carlo price is unbiased but carries a **standard error** that")
        A("shrinks like `1/sqrt(N)`. So the only meaningful comparison matches")
        A("the accuracies first:")
        A("")
        A("    N*  =  (sigma_eff / |lattice bias|)^2,")
        A("")
        A("and then puts the two wall times side by side.")
        A("")
        A("Parameters: rough Bergomi, `T=1`, `S0=K=100`, `xi0=0.09`, "
          f"`eta={ETA}`, `rho={RHO}`,")
        A("`mref = max(4, ceil(4 sqrt(n/8)))`, barrier `zmax = 3/sqrt(2H)`. The")
        A("reference is UNCAPPED, so the lattice's barrier truncation counts as")
        A("part of its error.")
        A("")

        A("## 0. The answer at n = 16")
        A("")
        A("| H | pricer | price | bias | lattice time | N* paths | MC time (paths only) | MC / lattice |")
        A("|---|---|---|---|---|---|---|---|")
        for H in HS:
            for label, n, m in PRICERS:
                if n != 16:
                    continue
                e, Lr = eq.get((H, label, n)), lat.get((H, label, n))
                if e is None or Lr is None:
                    continue
                A(f"| {H} | {label} | {Lr['price']:.4f} | {e['bias']:+.4f} "
                  f"| {human(e['lattice_time'])} | {hnum(e['n_star'])} "
                  f"| {human(e['mc_time_paths_only'])} "
                  f"| {e['ratio_paths_only']:.2f}x |")
        A("")
        A("`MC / lattice` above 1 means the Monte-Carlo is the slower way to")
        A("reach that accuracy; below 1 means it is the faster way.")
        A("")

        A("## 1. Monte-Carlo cost, measured")
        A("")
        A("| H | setup (cold) | marginal rate | sigma_eff (empirical) | "
          "sigma_eff (code's formula) | ratio |")
        A("|---|---|---|---|---|---|")
        for H in HS:
            d = direct[H]
            A(f"| {H} | {human(setup_cost[H])} | {path_rate[H]:,.0f} paths/s "
              f"| {d['control']['sigma_eff']:.4f} | {d['control']['sd']:.4f} "
              f"| {d['control']['inflation']:.4f} |")
        A("")
        A("`sigma_eff` is what the Monte-Carlo's error actually scales with, and")
        A("it is measured EXACTLY rather than estimated. With `N` paths laid out")
        A("as `N/2` antithetic pairs,")
        A("")
        A("    Var(mean) = Var(sample) (1 + rho_pair) / N,")
        A("")
        A("so `sigma_eff = sd(sample) * sqrt(1 + rho_pair)`, and both factors come")
        A("from one large sample with negligible uncertainty. The code's own")
        A("`stderr` is the `rho_pair = 0` case, so the last column is exactly the")
        A("factor by which that formula is wrong.")
        A("")
        A("### The antithetic device and the control variate work against each other")
        A("")
        A("| H | rho_pair (raw put) | rho_pair (control-adjusted) | sd (raw) | "
          "sd (control) | sigma_eff inflation |")
        A("|---|---|---|---|---|---|")
        for H in HS:
            d = direct[H]
            A(f"| {H} | {d['raw']['rho_pair']:+.4f} "
              f"| {d['control']['rho_pair']:+.4f} | {d['raw']['sd']:.3f} "
              f"| {d['control']['sd']:.3f} "
              f"| {d['control']['inflation']:.4f} |")
        A("")
        A("This is the finding of the phase, and it was not expected. On the RAW")
        A("put the antithetic device works well — a strongly negative pair")
        A("correlation, which is what antithetic sampling is for, since a put")
        A("payoff is monotone in the driver. But the `eta = 0` control variate")
        A("subtracts a second, similar put, and a DIFFERENCE of two puts is not")
        A("monotone in the driver. On that adjusted sample the pair correlation")
        A("turns POSITIVE, so the antithetic pairing is mildly counterproductive")
        A("and the quoted standard error is optimistic by the inflation factor.")
        A("")
        A("The control variate is by far the bigger effect — it cuts the sample")
        A("spread several-fold — so the net estimator is still much better than a")
        A("plain one. But two consequences follow, and both matter for the rest")
        A("of the project:")
        A("")
        A("1. Every `+-s.e.` quoted from this Monte-Carlo, and every 'within N MC")
        A("   bands' claim resting on one, is too generous by the inflation")
        A("   factor above. The `N*` columns here use the corrected `sigma_eff`.")
        A("2. Dropping the antithetic pairing (while keeping the control variate)")
        A("   would make the estimator slightly BETTER, not worse. That is a")
        A("   one-line change and it has not been made or tested here.")
        A("")
        A("### Cross-validation of the rebuilt sample")
        A("")
        A("`sigma_eff_direct` rebuilds the estimator's per-path sample, because")
        A("the pricer returns only aggregates. That duplication is a risk, so it")
        A("is checked against the real pricer at the same seed and path count.")
        A("")
        A("What the check can be is limited, and worth being precise about: the")
        A("pricer draws in CHUNKS from one rng stream while the rebuild draws the")
        A("whole block at once, so the two consume the stream differently and are")
        A("**different draws from the same distribution**. The price therefore")
        A("agrees only to within a couple of standard errors — that column is")
        A("reported in s.e. units so it judges itself. What must agree tightly is")
        A("the standard error, because that estimates a property of the")
        A("distribution rather than of the particular draw.")
        A("")
        A("| H | price (pricer) | price (rebuilt) | gap in s.e. | s.e. (pricer) | "
          "s.e. (rebuilt) | s.e. ratio |")
        A("|---|---|---|---|---|---|---|")
        for H in HS:
            x = xval[H]
            gap = ((x["price_rebuilt"] - x["price_pricer"])
                   / x["stderr_pricer"] if x["stderr_pricer"] > 0 else float("nan"))
            A(f"| {H} | {x['price_pricer']:.6f} | {x['price_rebuilt']:.6f} "
              f"| {gap:+.2f} | {x['stderr_pricer']:.6f} "
              f"| {x['stderr_rebuilt']:.6f} "
              f"| {x['stderr_rebuilt']/x['stderr_pricer']:.4f} |")
        A("")
        A("### The independent replication estimator, for comparison")
        A("")
        A(f"The same quantity from the scatter of {SIGMA_REPS} independent runs of")
        A(f"{SIGMA_PATHS:,} paths. This estimator needs no theory at all, which")
        A("makes it a genuine check — but a standard deviation from")
        A(f"{SIGMA_REPS} replications carries about")
        A(f"{100.0/math.sqrt(2.0*(SIGMA_REPS-1)):.0f}% relative uncertainty of")
        A("its own, and that uncertainty enters `N*` SQUARED, which is why it is")
        A("not the estimator used:")
        A("")
        A("| H | sigma_eff (exact) | sigma_eff (replications) | ratio |")
        A("|---|---|---|---|")
        for H in HS:
            se_ex = direct[H]["control"]["sigma_eff"]
            se_rp = sigma[H]["sigma_eff"]
            A(f"| {H} | {se_ex:.4f} | {se_rp:.4f} | {se_rp/se_ex:.3f} |")
        A("")
        A("Agreement to within that uncertainty is the check passing; a large")
        A("systematic gap in one direction across every `H` would mean the exact")
        A("formula is missing a correlation the replications can see. This")
        A("estimator was run at low power first — 24 replications, ±15% — and")
        A("came out about 1.3x high on every H, which looked systematic; raising")
        A(f"it to {SIGMA_REPS} replications (±"
          f"{100.0/math.sqrt(2.0*(SIGMA_REPS-1)):.0f}%) brought it back to about")
        A("1.08. The lesson is kept rather than quietly deleted: a std from two")
        A("dozen samples is not evidence of a systematic effect.")
        A("")
        A("The **setup** is a one-off `nfine x nfine` covariance factorisation")
        A(f"(`nfine = {NFINE}`). It is charged once per process, not per path,")
        A("which is why a small Monte-Carlo can measure slower than a large one")
        A("on a cold cache. Both readings are given below.")
        A("")

        A("## 2. The reference prices (the 'truth')")
        A("")
        A("| H | price | s.e. | 95% CI half-width | paths | wall time |")
        A("|---|---|---|---|---|---|")
        for H in HS:
            r = ref[H]
            A(f"| {H} | {r['price']:.4f} | {r['stderr']:.4f} "
              f"| {r['ci95']:.4f} | {r['paths']:,} | {human(r['wall'])} |")
        A("")
        A("Any bias smaller than about twice these standard errors is not")
        A("resolved, and the corresponding `N*` is a band rather than a number.")
        A("")

        A("## 3. Every pricer, every n")
        A("")
        A("| H | pricer | n | price | bias | bias / ref s.e. | resolved | "
          "state space | best-of-3 | median |")
        A("|---|---|---|---|---|---|---|---|---|---|")
        for H in HS:
            for label, n, m in PRICERS:
                Lr, e = lat.get((H, label, n)), eq.get((H, label, n))
                if Lr is None or e is None:
                    continue
                ss = Lr.get("state_space")
                A(f"| {H} | {label} | {n} | {Lr['price']:.4f} "
                  f"| {e['bias']:+.4f} | {e['bias_in_ref_se']:.1f} "
                  f"| {'yes' if e['resolved'] else 'NO'} "
                  f"| {ss if ss else '—'} | {human(Lr['best'])} "
                  f"| {human(Lr['median'])} |")
        A("")

        A("## 4. The equivalence, with the reference's uncertainty propagated")
        A("")
        A("`N*` band = what `N*` becomes if the true bias is at either end of")
        A("`|bias| ± 2 · ref s.e.`. A smaller true bias needs MORE paths, so the")
        A("band is wide when the bias is barely resolved.")
        A("")
        A("| H | pricer | n | N* | N* band | MC time (paths) | MC time (+setup) "
          "| MC/lattice (paths) | MC/lattice (+setup) | break-even error |")
        A("|---|---|---|---|---|---|---|---|---|---|")
        for H in HS:
            for label, n, m in PRICERS:
                e = eq.get((H, label, n))
                if e is None:
                    continue
                A(f"| {H} | {label} | {n} | {hnum(e['n_star'])} "
                  f"| {hnum(e['n_star_lo'])}–{hnum(e['n_star_hi'])} "
                  f"| {human(e['mc_time_paths_only'])} "
                  f"| {human(e['mc_time_with_setup'])} "
                  f"| {e['ratio_paths_only']:.2f}x "
                  f"| {e['ratio_with_setup']:.2f}x "
                  f"| {e['breakeven']:.4f} |")
        A("")
        A("**Break-even error** is the accuracy at which the Monte-Carlo, given")
        A("the lattice's wall time in paths, lands. If it is SMALLER than the")
        A("lattice's bias, the Monte-Carlo wins at that budget.")
        A("")

        A("## 5. What Monte-Carlo cannot match cheaply: the American price")
        A("")
        A("The lattice gets the early-exercise price from the same backward")
        A("pass. Monte-Carlo needs Longstaff--Schwartz, which is more expensive")
        A("and itself biased — this project has an LSM estimator for rough")
        A("Heston, not for rough Bergomi, so no equal-accuracy comparison is")
        A("offered here and none is implied.")
        A("")
        A("| H | pricer | n | European | American | premium | European time | "
          "American time |")
        A("|---|---|---|---|---|---|---|---|")
        for H in HS:
            for label, n, m in AMERICAN:
                a, Lr = amer.get((H, label, n)), lat.get((H, label, n))
                if a is None or Lr is None:
                    continue
                A(f"| {H} | {label} | {n} | {Lr['price']:.4f} "
                  f"| {a['price']:.4f} | {a['premium']:+.4f} "
                  f"| {human(a['eu_time'])} | {human(a['best'])} |")
        A("")

        A("## 6. How to read this — and how not to")
        A("")
        A("- **This compares implementations, not methods.** The Monte-Carlo's")
        A("  inner loop is a BLAS matrix product and runs multi-threaded; the")
        A("  lattices are single-threaded Python loops over a state space, with")
        A("  no C++ port for Route B. A constant factor of tens is")
        A("  implementation, not mathematics.")
        A("- **The scalings, so a reader can rescale:** Monte-Carlo is")
        A("  `O(N x nfine)` per price; the one-step lattice is `O(n x nx)`; the")
        A("  lift at `m` factors is `O(n x nx x prod_i N_i)`, and that state")
        A("  space is what makes `m = 2` expensive.")
        A("- **A small bias is not a good thing here.** Part V established that")
        A("  the one-step lattice's error at small `n` is a CANCELLATION between")
        A("  a shrinking discretisation term and a growing variance term, with")
        A("  the crossing near `n* = (2H)^{-1/(1-2H)}` — about 13 at H = 0.05,")
        A("  7 at H = 0.1, 4 at H = 0.3. So a flattering bias at `n = 16` is")
        A("  luck, not accuracy: it is payoff-dependent and carries no error")
        A("  bound. Refining `n` makes the one-step scheme WORSE, which is")
        A("  visible in the n = 32 and n = 64 rows of section 3.")
        A("- **The lattice's case is not speed.** It is that Route B is a")
        A("  convergent scheme with a proof and an American price in the same")
        A("  pass. Nothing in this table argues it is the fast way to get a")
        A("  European number.")
        A("")

        A("## Timing")
        A("")
        A("| phase | wall time |")
        A("|---|---|")
        for name, dt in timings:
            A(f"| {name} | {human(dt)} |")
        A(f"| **total** | **{human(total)}** |")
        A("")
        A(f"Machine: single process. Reference {REF_PATHS:,} paths at "
          f"nfine={NFINE} per H; sigma from {SIGMA_REPS}x{SIGMA_PATHS:,}; "
          f"lattice timings best-of-{REPEATS}.")

        pr.write_results_md("\n".join(L))

        (pr.dir / "cost_equivalence.json").write_text(json.dumps({
            "params": {"eta": ETA, "rho": RHO, "nfine": NFINE,
                       "ref_paths": REF_PATHS, "sigma_paths": SIGMA_PATHS,
                       "sigma_reps": SIGMA_REPS, "repeats": REPEATS},
            "setup_cost": {str(h): setup_cost[h] for h in HS},
            "path_rate": {str(h): path_rate[h] for h in HS},
            "sigma": {str(h): sigma[h] for h in HS},
            "sigma_direct": {str(h): direct[h] for h in HS},
            "sigma_xval": {str(h): xval[h] for h in HS},
            "reference": {str(h): ref[h] for h in HS},
            "lattice": {f"H={h},{lb},n={n}": v
                        for (h, lb, n), v in lat.items() if v},
            "equivalence": {f"H={h},{lb},n={n}": v
                            for (h, lb, n), v in eq.items()},
            "american": {f"H={h},{lb},n={n}": v
                         for (h, lb, n), v in amer.items() if v},
        }, indent=2, default=float))

        for H in HS:
            e = eq.get((H, "one-step", 16))
            if e:
                pr.result(f"n16_onestep_ratio_H{H}", e["ratio_paths_only"])
                pr.result(f"n16_onestep_nstar_H{H}", e["n_star"])


if __name__ == "__main__":
    main()
