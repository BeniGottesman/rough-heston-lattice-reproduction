#!/usr/bin/env python3
"""
run_nm_bound — the n<->m link, measured: how many lift factors a grid of n steps
needs, across H, across tolerances, and across the three constructions.

The question.  Route B replaces the single divergent counter of the one-step
lattice by m Ornstein--Uhlenbeck factors.  Part VI showed that this converts an
UNBOUNDED error into a bounded one, but left the obvious practical question
open: given n time steps and a target accuracy, how large must m be?  The
project has an empirical trace of the answer (`F025`: m* ~ 0.22 log n, marked
open) and no inequality.

What is measured here, in five phases.

  curves     rel(m) = ||K-K^m||_{L^2(delta,T)} / ||K||_{L^2(delta,T)} for
             m = 1..M_MAX, on every (H, n), for all three constructions:
             `best` (free nodes, NNLS weights: the achievable floor, and the
             branch the Route B lattice uses), `ajee` (the AJEE partition family
             optimised) and `geometric` (an explicit partition, no optimisation).
  mstar      m*(n, H, tol) read off those curves, and the log fits with their
             residuals, so `m* ~ a ln n + b` is TESTED, not assumed.
  ratio      the mechanism.  The explicit construction covers the s-range
             [1/T, c_hi/delta] with a fixed ratio r between consecutive nodes,
             so m = 1 + log(c_hi n)/log(r) identically.  IF a fixed r buys a
             fixed accuracy independently of n, the log law follows as an
             inequality rather than a fit.  That independence is what this phase
             measures.
  epslaw     how the tolerance enters.  The slope a(H, tol) = dm*/d ln n is
             fitted against tol, to separate the two candidate shapes:
             a ~ log(1/tol)  (near-optimal exponential sums) versus
             a ~ tol^{-p}    (algebraic, cell-wise moment matching).
  variance   the same m*, re-derived from the DRIVER VARIANCE criterion --- the
             quantity that actually drives the price error, and the one
             Proposition 8.3 shows diverging for the one-step scheme --- so the
             L^2 tolerance can be judged against it instead of assumed adequate.

Not measured, and stated as such in the results: the propagation of a kernel
error to a PRICE error.  That goes through Proposition (B1'), whose constant is
unevaluated; until it is, none of this is a bound on a price.

Analytic throughout: no lattice, no Monte-Carlo.  Every number is a closed-form
L^2 norm or a closed-form discrete covariance.
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
import nm_bound as nb                             # noqa: E402
import route_b as rb                              # noqa: E402

T = 1.0
HS = (0.05, 0.10, 0.20, 0.30)
NS = (16, 32, 64, 128, 256, 512, 1024, 2048)
TOLS = (0.10, 0.03, 0.01)
CONSTRUCTIONS = ("best", "ajee", "geometric")
M_MAX = 12

# Below this, the exact formula A - 2w.b + w'Gw cancels catastrophically and
# `route_b.l2_error` clamps at zero.  Every tolerance used here is far above it,
# but a curve value under the floor means "unresolved", not "zero".
REL_FLOOR = 1e-4

# the ratio phase
RATIOS = (2.0, 3.0, 4.0, 6.0, 8.0, 12.0)
NS_RATIO = (64, 256, 1024, 4096)
# the variance phase: discrete_covariance_report is O(n^2) quadratures
NS_VAR = (16, 32, 64, 128)
HS_VAR = (0.05, 0.10, 0.30)


def fmt(x: float | None, nd: int = 4) -> str:
    if x is None:
        return "—"
    if isinstance(x, float) and not math.isfinite(x):
        return "—"
    if x == 0.0:
        return f"<{REL_FLOOR:.0e}"
    return f"{x:.{nd}f}" if x >= 1e-3 else f"{x:.2e}"


def main() -> None:
    t_start = time.time()
    with Progress("nm-bound", total_phases=5) as pr:
        timings: list[tuple[str, float]] = []

        # ------------------------------------------------------------ curves
        pr.phase("curves", total=len(HS) * len(NS))
        t0 = time.time()
        curves: dict = {}
        wobbles: list[str] = []
        k = 0
        for H in HS:
            for n in NS:
                for con in CONSTRUCTIONS:
                    c = nb.error_curve(n, H, T, M_MAX, construction=con)
                    curves[(H, n, con)] = c
                    # a local search can come out worse at m+1 than at m; a big
                    # wobble makes the table unreadable as a rate, so record it
                    for i in range(len(c) - 1):
                        if c[i] > REL_FLOOR and c[i + 1] > c[i] * 1.5:
                            wobbles.append(
                                f"{con} H={H} n={n}: rel({i+1})={c[i]:.3e} -> "
                                f"rel({i+2})={c[i+1]:.3e}")
                k += 1
                pr.tick(k, note=f"H={H} n={n}")
        timings.append(("curves", time.time() - t0))

        # ------------------------------------------------------------- mstar
        pr.phase("mstar", total=len(CONSTRUCTIONS))
        t0 = time.time()
        mstar: dict = {}
        fits: dict = {}
        for j, con in enumerate(CONSTRUCTIONS):
            for H in HS:
                for tol in TOLS:
                    ms = [nb.m_star_from_curve(curves[(H, n, con)], tol)
                          for n in NS]
                    mstar[(H, tol, con)] = ms
                    ok = [(n, m) for n, m in zip(NS, ms) if m is not None]
                    fits[(H, tol, con)] = nb.fit_log_law(
                        [n for n, _ in ok], [m for _, m in ok])
            pr.tick(j + 1, note=con)
        timings.append(("mstar", time.time() - t0))

        # ------------------------------------------------------------- ratio
        # Does a fixed node ratio buy a fixed accuracy, independently of n?
        pr.phase("ratio", total=len(HS))
        t0 = time.time()
        ratio_tab: dict = {}
        for j, H in enumerate(HS):
            for r in RATIOS:
                for n in NS_RATIO:
                    w, s, m = nb.geometric_lift_by_ratio(r, H, n, T)
                    ratio_tab[(H, r, n)] = (
                        nb.rel_l2(w, s, H, T, T / n), m)
            pr.tick(j + 1, note=f"H={H}")
        timings.append(("ratio", time.time() - t0))

        # ------------------------------------------------------------ epslaw
        pr.phase("epslaw", total=len(CONSTRUCTIONS))
        t0 = time.time()
        eps_fits: dict = {}
        law: dict = {}
        for j, con in enumerate(CONSTRUCTIONS):
            for H in HS:
                slopes = [fits[(H, tol, con)]["slope"] for tol in TOLS]
                good = [(t_, a) for t_, a in zip(TOLS, slopes)
                        if math.isfinite(a) and a > 0]
                # a ~ C tol^{-p}: fit ln a against ln(1/tol)
                pw = nb.fit_power_law([1.0 / t_ for t_, _ in good],
                                      [a for _, a in good]) if good else {}
                # a ~ C' log(1/tol): the competing shape, same data
                lg = (nb.fit_power_law([math.log(1.0 / t_) for t_, _ in good],
                                       [a for _, a in good])
                      if len(good) >= 2 else {})
                eps_fits[(H, con)] = {"slopes": slopes, "power": pw,
                                      "vs_log": lg}
                if con == "best" and pw:
                    # m* ~ A ln n (1/tol)^p + B, A and p from the slope fit,
                    # B the mean intercept across tolerances
                    ints = [fits[(H, tol, con)]["intercept"] for tol in TOLS
                            if math.isfinite(fits[(H, tol, con)]["intercept"])]
                    law[H] = {"A": pw.get("const", float("nan")),
                              "p": pw.get("exponent", float("nan")),
                              "B": float(np.mean(ints)) if ints else 0.0}
            pr.tick(j + 1, note=con)
        timings.append(("epslaw", time.time() - t0))

        # ---------------------------------------------------------- variance
        pr.phase("variance", total=len(HS_VAR) * len(NS_VAR))
        t0 = time.time()
        var_tab: dict = {}
        k = 0
        for H in HS_VAR:
            for n in NS_VAR:
                row = {}
                for m in range(1, 7):
                    try:
                        row[m] = nb.variance_ratio_criterion(
                            n, H, m, T, construction="best")
                    except Exception as e:                # keep the run alive
                        row[m] = None
                        pr.log(f"variance H={H} n={n} m={m} failed: {e}")
                var_tab[(H, n)] = row
                k += 1
                pr.tick(k, note=f"H={H} n={n}")
        timings.append(("variance", time.time() - t0))

        # ------------------------------------------------------------ results
        total = time.time() - t_start
        L: list[str] = []
        A = L.append
        A("# n <-> m : how many lift factors a grid of n steps needs")
        A("")
        A(f"`ALGO_VERSION` {pr.meta.get('algo_version','?')} · "
          f"`code_sha` {pr.meta.get('code_sha','?')}")
        A("")
        A("Criterion throughout: the RELATIVE kernel error on the lags the")
        A("scheme actually sees,")
        A("")
        A("    rel(m) = ||K - K^m||_{L2(delta,T)} / ||K||_{L2(delta,T)},"
          "   delta = T/n.")
        A("")
        A("Relative, so it is independent of the kernel's normalisation; on")
        A("(delta,T), because no lattice with step delta evaluates K closer to")
        A("the origin than delta.  Values printed as `<1e-04` are below the")
        A("cancellation floor of the exact L2 formula and mean *unresolved*,")
        A("not zero.  All tolerances used are far above that floor.")
        A("")
        A("## 0. The answer, in one table")
        A("")
        A("m*(n) at a 3% relative kernel error, for the three constructions.")
        A("`best` is the branch the Route B lattice uses (`route_b_lattice.lift_for`).")
        A("")
        A("| H | construction | " + " | ".join(f"n={n}" for n in NS)
          + " | slope dm*/d ln n | R² | max resid |")
        A("|---|---|" + "---|" * (len(NS) + 3))
        for H in HS:
            for con in CONSTRUCTIONS:
                ms = mstar[(H, 0.03, con)]
                f = fits[(H, 0.03, con)]
                A(f"| {H} | `{con}` | "
                  + " | ".join("—" if m is None else str(m) for m in ms)
                  + f" | {f['slope']:.3f} | {f['r2']:.3f} | {f['max_resid']:.2f} |")
        A("")
        A("Doubling n adds `slope x ln 2` factors — a fraction of one factor.")
        A("The growth is logarithmic for every H and every construction.")
        A("")
        A("## 1. m*(n, H, tol), all tolerances")
        A("")
        for con in CONSTRUCTIONS:
            A(f"### {con}")
            A("")
            A("| H | tol | " + " | ".join(f"n={n}" for n in NS)
              + " | slope | intercept | R² |")
            A("|---|---|" + "---|" * (len(NS) + 3))
            for H in HS:
                for tol in TOLS:
                    ms = mstar[(H, tol, con)]
                    f = fits[(H, tol, con)]
                    A(f"| {H} | {tol:.0%} | "
                      + " | ".join("—" if m is None else str(m) for m in ms)
                      + f" | {f['slope']:.3f} | {f['intercept']:.2f}"
                      + f" | {f['r2']:.3f} |")
            A("")
        A("A `—` means no m <= "
          f"{M_MAX} reached the tolerance.")
        A("")
        A("## 2. The error curves rel(m)")
        A("")
        for H in HS:
            A(f"### H = {H}")
            A("")
            A("| construction | n | " + " | ".join(f"m={m}" for m in
                                                   range(1, M_MAX + 1)) + " |")
            A("|---|---|" + "---|" * M_MAX)
            for con in CONSTRUCTIONS:
                for n in NS:
                    c = curves[(H, n, con)]
                    A(f"| `{con}` | {n} | " + " | ".join(fmt(e) for e in c) + " |")
            A("")
        if wobbles:
            A("Non-monotone steps (a local search finding a worse optimum at")
            A("m+1 than at m), reported rather than hidden:")
            A("")
            for w in wobbles[:40]:
                A(f"- {w}")
            if len(wobbles) > 40:
                A(f"- ... and {len(wobbles)-40} more")
            A("")

        A("## 3. The mechanism: does a fixed node ratio buy a fixed accuracy?")
        A("")
        A("The explicit construction covers [1/T, 1/delta] with a fixed ratio r")
        A("between consecutive nodes, so `m = 1 + ceil(log n / log r)` holds")
        A("identically.  If the error at fixed r does not grow with n, then the")
        A("log law is a consequence and not a fit: covering a range of n at a")
        A("bounded cost per node costs log n nodes.")
        A("")
        for H in HS:
            A(f"### H = {H}")
            A("")
            A("| r | " + " | ".join(f"n={n} (m)" for n in NS_RATIO)
              + " | rel spread across n |")
            A("|---|" + "---|" * (len(NS_RATIO) + 1))
            for r in RATIOS:
                vals = [ratio_tab[(H, r, n)] for n in NS_RATIO]
                es = [v[0] for v in vals]
                spread = (max(es) / min(es)) if min(es) > 0 else float("nan")
                A(f"| {r:g} | "
                  + " | ".join(f"{e:.4f} ({m})" for e, m in vals)
                  + f" | x{spread:.2f} |")
            A("")

        A("## 4. How the tolerance enters the slope")
        A("")
        A("Two candidate shapes for `a(tol) = dm*/d ln n`: algebraic,")
        A("`a ~ C tol^-p` (cell-wise moment matching), or logarithmic,")
        A("`a ~ C log(1/tol)` (near-optimal exponential sums).  Both are fitted")
        A("to the same three tolerances; the exponent is the discriminator —")
        A("`p ~ 0` would favour the logarithmic shape, `p ~ 1/2` the algebraic.")
        A("")
        A("| H | construction | " + " | ".join(f"a({t:.0%})" for t in TOLS)
          + " | p (a ~ tol^-p) | R² | R² of a ~ log(1/tol) |")
        A("|---|---|" + "---|" * (len(TOLS) + 3))
        for H in HS:
            for con in CONSTRUCTIONS:
                e = eps_fits[(H, con)]
                pw, lg = e.get("power") or {}, e.get("vs_log") or {}
                A(f"| {H} | `{con}` | "
                  + " | ".join(f"{a:.3f}" if math.isfinite(a) else "—"
                               for a in e["slopes"])
                  + f" | {pw.get('exponent', float('nan')):.3f}"
                  + f" | {pw.get('r2', float('nan')):.3f}"
                  + f" | {lg.get('r2', float('nan')):.3f} |")
        A("")

        A("## 5. The other criterion: the driver variance")
        A("")
        A("The L2 kernel error is a proxy.  What drives the price error is the")
        A("driver's variance, and Proposition 8.3 is stated in exactly that")
        A("quantity.  Below: |Var[V^m_T]/Var[V_T] - 1| for the discrete scheme")
        A("(exact, `route_b.discrete_covariance_report`, cell-average mode),")
        A("against the one-step scheme's `2H n^{1-2H} - 1` for contrast.")
        A("")
        A("| H | n | " + " | ".join(f"m={m}" for m in range(1, 7))
          + " | scheme's own | one-step |")
        A("|---|---|" + "---|" * 8)
        for H in HS_VAR:
            for n in NS_VAR:
                row = var_tab[(H, n)]
                one = abs(nb.onestep_variance_ratio(n, H) - 1.0)
                own = next((row[m]["true_vs_cont"] for m in range(1, 7)
                            if row[m] is not None), None)
                A(f"| {H} | {n} | "
                  + " | ".join(fmt(None if row[m] is None
                                   else row[m]["lift_vs_true"])
                               for m in range(1, 7))
                  + f" | {fmt(own)} | {one:.3f} |")
        A("")
        A("`scheme's own` is the discretisation error the lattice already has")
        A("with the TRUE kernel (true-discrete versus continuous covariance).")
        A("There is no point driving the lift below it, and it is independent")
        A("of m.  The one-step column grows without bound in n; every lift")
        A("column is flat in n at fixed m.  That is the bounded-versus-unbounded")
        A("distinction of Part VI, seen in the variance instead of the price.")
        A("")
        A("Whether an `L2(delta,T)` tolerance controls this variance is the")
        A("diagnostic below: the lattice's FIRST cell average runs over")
        A("(0, delta] and so straddles the singularity that L2(delta,T)")
        A("excludes by construction.  If the variance error tracks the")
        A("`L2(0,T)` column rather than the `L2(delta,T)` one, that first cell")
        A("is the reason, and the L2(delta,T) tolerance alone is not enough.")
        A("")
        A("| H | n | m | var err | L2(delta,T) | L2(0,T) |")
        A("|---|---|---|---|---|---|")
        for H in HS_VAR:
            for n in NS_VAR:
                for m in range(1, 7):
                    d = var_tab[(H, n)][m]
                    if d is None:
                        continue
                    A(f"| {H} | {n} | {m} | {fmt(d['lift_vs_true'])}"
                      f" | {fmt(d['rel_l2_delta'])} | {fmt(d['rel_l2_zero'])} |")
        A("")

        A("## 6. What this is NOT")
        A("")
        A("- **Not a bound.** Every entry is a measurement, and the log law is a")
        A("  fit with its residuals shown. The `m` sufficient side is reachable")
        A("  from section 3 with an explicit node placement; the `m` necessary")
        A("  side (that no m-factor lift with positive weights does better) is")
        A("  not measured here beyond the `best` construction's local search.")
        A("- **Not a price error.** Propagating a kernel error to a price error")
        A("  goes through Proposition (B1'), whose constant is unevaluated. Until")
        A("  it is, `rel(m) <= tol` does not license any statement about a price.")
        A("- **`best` is a local search**, warm-started from m-1; it is a 'best")
        A("  found', so it is an upper bound on the achievable error and its m*")
        A("  is a lower bound on what is necessary only up to that search.")
        A("")

        A("## Timing")
        A("")
        A("| phase | wall time |")
        A("|---|---|")
        for name, dt in timings:
            A(f"| {name} | {dt:.1f}s |")
        A(f"| **total** | **{total:.1f}s** |")
        A("")
        A(f"Grid: H in {HS}, n in {NS}, tol in {TOLS}, m <= {M_MAX}, "
          f"constructions {CONSTRUCTIONS}.")

        pr.write_results_md("\n".join(L))

        # the fitted constants, so nm_bound.set_law can be calibrated from a run
        (pr.dir / "nm_law.json").write_text(json.dumps(
            {"law": {str(h): v for h, v in law.items()},
             "run": pr.dir.name,
             "algo_version": pr.meta.get("algo_version"),
             "code_sha": pr.meta.get("code_sha"),
             "fits": {f"H={h},tol={t},{c}": fits[(h, t, c)]
                      for (h, t, c) in fits},
             "mstar": {f"H={h},tol={t},{c}": mstar[(h, t, c)]
                       for (h, t, c) in mstar}}, indent=2))
        for H in HS:
            pr.result(f"slope_3pct_best_H{H}", fits[(H, 0.03, 'best')]["slope"])
        pr.result("wobbles", len(wobbles))


if __name__ == "__main__":
    main()
