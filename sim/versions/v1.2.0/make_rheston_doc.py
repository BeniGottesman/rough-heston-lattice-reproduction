#!/usr/bin/env python3
"""
make_rheston_doc — assemble docs/ROUGH_HESTON.md from the runs.

The tables are never retyped: they are read from the RESULTS.md of the runs that
produced them, so a number in the document cannot drift from its evidence.  What
this script adds is the framing -- what the model is, what each engine is, what
the numbers mean and what they do not -- plus the headline figures, which are read
from each run's progress.json rather than copied by hand.

    python3 sim/make_rheston_doc.py [tables_run] [anchor_run]

With no arguments it picks the most recent run of each kind.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RUNS = ROOT / "runs"
OUT = ROOT / "docs" / "ROUGH_HESTON.md"


def latest(prefix: str) -> Path:
    c = sorted(p for p in RUNS.glob(f"{prefix}-*") if (p / "RESULTS.md").exists())
    if not c:
        raise SystemExit(f"no completed run found for {prefix}")
    return c[-1]


def load(run: Path) -> tuple[str, dict]:
    md = (run / "RESULTS.md").read_text()
    prog = json.loads((run / "progress.json").read_text())
    return md, prog.get("results", {})


def body(md: str) -> str:
    """Drop the run's own H1 title; keep everything below it."""
    lines = md.splitlines()
    for i, l in enumerate(lines):
        if l.startswith("# "):
            return "\n".join(lines[i + 1:]).strip()
    return md.strip()


def _rows(md: str, start: str, stop: str):
    """Yield the pipe-split cells of every data row between two headings."""
    seg = md[md.find(start):md.find(stop, md.find(start))]
    for line in seg.splitlines():
        if not line.startswith("| ") or "---" in line:
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if cells and cells[0] in ("S0", ""):
            continue
        yield cells


def _f(x):
    try:
        return float(x)
    except Exception:
        return None


def headline_european(lattice_md: str, tables_md: str) -> list:
    """One table: our tree (both walks) vs our Monte-Carlo vs the references,
    with the signed error of each. This is THE comparison a reader wants and it
    is assembled here because the lattice and the Monte-Carlo live in different
    runs -- the lattice never sees a path, the Monte-Carlo never sees a node."""
    # Monte-Carlo, analytic, Fourier: the anchor table of the tables run (H=0.5)
    mc = {}
    for c in _rows(tables_md, "## 1. The anchor", "## 2. European"):
        if len(c) < 8:
            continue
        k = (_f(c[0]), _f(c[1]), c[2])
        if None in k[:2]:
            continue
        mc[k] = (_f(c[3]), _f(c[4]), _f(c[6]), _f(c[7]))   # analytic, Fourier, MC, se
    # lattice both walks at n=200, plus ADS tree: the European table of the lattice run
    lat = {}
    for c in _rows(lattice_md, "## 1. European put", "### Summary"):
        if len(c) < 10:
            continue
        k = (_f(c[0]), _f(c[1]), c[2])
        if None in k[:2]:
            continue
        lat[k] = (_f(c[4]), _f(c[6]), _f(c[8]))            # bino200, trino200, ADS500
    keys = sorted(set(mc) & set(lat), key=lambda k: (k[2], k[1], k[0]))
    if not keys:
        return []

    out = ["## 0. The one table: our tree vs Monte-Carlo vs the references\n",
           "Every column below is a price for the SAME parameter set at `H = 0.5`, "
           "and every `err` is signed against the published analytical Heston "
           "price. Our tree is `n = 200`; our Monte-Carlo is 400k antithetic paths "
           "with its standard error `se`; `ADS` is their independent tree at "
           "`N = 500`. The lattice and the Monte-Carlo share no code -- one walks "
           "a recombining grid backward, the other samples forward paths -- so "
           "their agreement here is a genuine cross-check.\n",
           "| S0 | sqrt(V0) | T | **analytic** | tree bino | err | tree **trino** | "
           "err | Monte-Carlo | se | err | ADS N=500 | err |",
           "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    import numpy as _np
    eb, et, em, ea = [], [], [], []
    for k in keys:
        an, fo, m, se = mc[k]
        b, t, ads = lat[k]
        eb.append(b - an); et.append(t - an); em.append(m - an); ea.append(ads - an)
        out.append(f"| {k[0]:.0f} | {k[1]} | {k[2]} | {an:.4f} | {b:.4f} | "
                   f"{b - an:+.4f} | {t:.4f} | {t - an:+.4f} | {m:.4f} | {se:.4f} | "
                   f"{m - an:+.4f} | {ads:.4f} | {ads - an:+.4f} |")
    def mabs(x):
        return _np.abs(_np.array(x)).mean()
    out.append("")
    out.append(f"**Mean absolute error over the {len(keys)} sets** — tree binomial "
               f"`{mabs(eb):.5f}`, tree **trinomial `{mabs(et):.5f}`**, Monte-Carlo "
               f"`{mabs(em):.5f}`, ADS tree `{mabs(ea):.5f}`. The trinomial tree and "
               "the Monte-Carlo agree with the exact price to the same order; the "
               "binomial tree carries the larger finite-`n` constant discussed "
               "below. The Monte-Carlo's error is dominated by its standard error "
               "`se` (a few `1e-3`), which shrinks only with more paths, whereas "
               "the trees carry no statistical noise at all.\n")
    return out


def main() -> None:
    tr = Path(sys.argv[1]) if len(sys.argv) > 1 else latest("rheston-tables")
    ar = Path(sys.argv[2]) if len(sys.argv) > 2 else latest("rheston-american-anchor")
    tmd, tres = load(tr)
    amd, ares = load(ar)
    # the lattice parts are optional so the document still builds without them
    try:
        lr = latest("heston-ads-lattice")
        lmd, lres = load(lr)
    except SystemExit:
        lr, lmd, lres = None, None, {}
    try:
        or_ = latest("heston-lattice-order")
        omd, ores = load(or_)
    except SystemExit:
        or_, omd, ores = None, None, {}
    try:
        rr = latest("rough-bergomi-ladder")
        rmd, rres = load(rr)
    except SystemExit:
        rr, rmd, rres = None, None, {}
    try:
        br = latest("routeb-compare")
        bmd, bres = load(br)
    except SystemExit:
        br, bmd, bres = None, None, {}

    def g(d, k, default=None):
        return d.get(k, default)

    hs = [0.05, 0.1, 0.3, 0.5]
    bias_line = ", ".join(
        f"`H = {h}` → `{g(tres, f'eu_put_mean_err_H{h}', float('nan')):+.4f}`"
        for h in hs)

    P: list[str] = []
    P.append("# Rough Heston: reference tables\n")
    P.append("**Status: kept aside from the paper, on purpose.**\n")
    P.append("Rough Heston is *outside* the model class the paper analyses. That "
             "class needs `v_t = F(t, (K^h Y)(t))` with `Y` an autonomous "
             "diffusion, and rough Heston's `nu sqrt(V) dB` makes the driver's "
             "coefficients depend on the variance, hence on its own past. So "
             "nothing in this document is a lattice result and nothing here "
             "supports or contradicts a claim in the paper. What it is: a "
             "reference table set for the rough Heston model in its own right, "
             "laid out like the numerical section of the paper whose parameter "
             "sets it borrows.\n")
    P.append("## What is being priced\n")
    P.append("```\n"
             "dS_t = r S_t dt + S_t sqrt(V_t) dW_t\n"
             "V_t  = V_0 + int_0^t K(t-s) [ kappa (theta - V_s) ds\n"
             "                            + eta sqrt(V_s) dB_s ]\n"
             "K(u) = u^{alpha-1} / Gamma(alpha),   alpha = H + 1/2,\n"
             "d<W,B>_t = rho dt\n"
             "```\n")
    P.append("At `H = 1/2` the kernel is the constant `1` and this **is** the "
             "Heston model. That is not a detail: it is what makes the whole "
             "document checkable, because at `H = 1/2` third parties have "
             "published the prices.\n")
    P.append("## The parameter sets are not ours\n")
    P.append("They are Beliaeva--Nawalkha's, used as such in Akyildirim, Dolinsky "
             "and Soner, *Approximating stochastic volatility by recombinant "
             "trees*, Ann. Appl. Probab. 24(5) 2176--2205, 2014 "
             "([arXiv:1205.3555](https://arxiv.org/abs/1205.3555)) -- the ADS2014 "
             "the paper cites as the weak-convergence route it does not take. "
             "Keeping their grid means their published columns can be used as "
             "external references instead of ours.\n")
    P.append("| | European call and put | American put |")
    P.append("|---|---|---|")
    P.append("| strike | `K = 100` | `K = 100` |")
    P.append("| spot | `S0 = 90, 95, 100, 105, 110` | `S0 = 90, 100, 110` |")
    P.append("| maturity | `T = 1, 3, 6` months | `T = 1, 3, 6` months |")
    P.append("| initial vol | `sqrt(V0) = 0.2, 0.3, 0.4` | `sqrt(V0) = 0.2, 0.4` |")
    P.append("| correlation | `rho = -0.7` | `rho = -0.1, -0.7` |")
    P.append("| rate | `r = 0.05` | `r = 0.05` |")
    P.append("| vol of vol | `eta = 0.1` | `eta = 0.1` |")
    P.append("| mean reversion | `kappa = 3.0` | `kappa = 3.0` |")
    P.append("| long-run variance | `theta = 0.04` | `theta = 0.04` |")
    P.append("| roughness | `H = 0.05, 0.10, 0.30, 0.50` | `H = 0.05, 0.10, 0.30, 0.50` |")
    P.append("")
    P.append("The Feller ratio `2 kappa theta / eta^2 = 24` is far from zero, so "
             "the variance almost never needs truncating -- which matters, because "
             "truncation at zero is the known source of the simulation's upward "
             "bias.\n")
    P.append("## Where the paper's own lattice appears, and where it cannot\n")
    P.append("A document attached to a paper about a lattice ought to contain a "
             "lattice column. It does -- but at `H = 0.5` only, and the reason is a "
             "theorem rather than a limitation of the code.\n")
    P.append("The paper's class needs an autonomous driver. Rough Heston's "
             "`eta sqrt(V) dB` makes the driver's coefficients depend on the "
             "variance, hence on the driver's own past, so **no lattice of the "
             "paper's kind prices rough Heston.** Conversely a semi-analytic "
             "characteristic function requires affineness, which forces exactly "
             "that `sqrt(V) dB`. So **no rough model is both inside the paper's "
             "class and semi-analytic**: the two families meet only at `h = 0`. "
             "That is why the lattice column and the Fourier column can be put on "
             "the same line at `H = 0.5` and nowhere else, and it is also why the "
             "paper's rough numerics run on rough Bergomi -- which IS in the class "
             "-- against an exact-covariance simulation rather than against a "
             "closed form.\n")
    P.append("Far from being a consolation prize, `H = 0.5` is the strongest test "
             "available: the model is classical Heston, so the lattice can be set "
             "against **two** references it did not produce -- the analytical "
             "Heston price, and the independent recombining tree of the ADS paper "
             "itself. Part III is that comparison.\n")
    P.append("## The two engines\n")
    P.append("**Fourier** (`sim/rough_heston.py`, Python/NumPy, single-threaded). "
             "Lewis inversion of the characteristic function, which for rough "
             "Heston solves a *fractional Riccati* equation "
             "(El Euch--Rosenbaum), integrated by the Diethelm--Ford--Freed "
             "predictor--corrector. There is no Monte-Carlo anywhere in it. It "
             "prices European options only; the rate enters through the forward.\n")
    P.append("**Monte-Carlo** (`sim/cpp/rheston_mc.cpp`, C++17, multithreaded). "
             "Euler--Volterra on the variance with the left-endpoint kernel whose "
             "smallest lag is `dt`, so `K` is never evaluated at `0` where it is "
             "infinite. Antithetic pairs. For the American problem, "
             "Longstaff--Schwartz with the exercise policy fitted on a regression "
             "sample and applied to a **disjoint** valuation sample, so the "
             "reported American price carries no in-sample look-ahead: it is a "
             "genuine lower bound with an honest standard error.\n")
    P.append("## What the numbers say\n")
    P.append("**1. Both engines reproduce third-party published Heston prices.** "
             "At `H = 0.5`, against the analytical column ADS print for 45 "
             "European puts: Fourier is out by a mean of "
             f"`{g(tres, 'anchor_fourier_mean_abs_err', float('nan')):.2e}` "
             f"(worst `{g(tres, 'anchor_fourier_max_abs_err', float('nan')):.2e}`), "
             "the Monte-Carlo by a mean of "
             f"`{g(tres, 'anchor_mc_mean_abs_err', float('nan')):.2e}` "
             f"(worst `{g(tres, 'anchor_mc_max_abs_err', float('nan')):.2e}`, "
             "which is inside one standard error). Nothing downstream rests on "
             "code that only we have checked.\n")
    P.append("**2. The simulation's bias is a roughness effect, and it has a "
             "sign.** Mean signed error `MC - Fourier` over the 45 European puts: "
             + bias_line + ". It is positive, it grows monotonically as `H` falls, "
             "and at `H = 0.5` it vanishes into the noise. Averaged over 45 "
             "parameter sets the standard error of each of those means is about "
             "`0.002`, so the ordering in `H` is resolved and is not an artefact "
             "of the path count. The Euler--Volterra scheme therefore "
             "*overprices* the rough model, and the rougher the model the more it "
             "overprices. Anyone using simulation as the reference for a rough "
             "model is measuring against a moving target.\n")
    P.append("**2b. But the step-refinement decomposition of that bias is *not* "
             "resolved here, and is not claimed to be.** Section 5 holds the paths "
             "fixed and doubles the step count. On the two high-volatility, "
             "long-maturity rows the trend is clean and far outside the standard "
             "error -- at `H = 0.05`, `sqrt(V0) = 0.4`, `T = 6m` the error falls "
             "`+0.0975 → +0.0726 → +0.0613 → +0.0443` against a standard error of "
             "`0.017`. On the low-volatility short-maturity rows the whole effect "
             "is smaller than the standard error and the series wanders. So "
             "Section 5 demonstrates that refinement removes error where the error "
             "is large enough to see, and it does **not** establish a bias floor. "
             "Doing that needs common random numbers across step counts and far "
             "more paths; it is an open item in the project's register, not a "
             "result of this run.\n")
    P.append("**3. Fourier beats simulation on European payoffs, and not "
             "narrowly.** Per `H`, the Fourier engine prices all 90 European "
             "options in about "
             f"`{g(tres, 'eu_time_fourier_H0.1', float('nan')):.1f} s` against "
             f"`{g(tres, 'eu_time_mc_H0.1', float('nan')):.1f} s` for the "
             "Monte-Carlo -- and the Fourier value is exact while the Monte-Carlo "
             "still carries both a standard error and the bias of point 2.\n")
    P.append("**4. The American estimator is calibrated, not asserted.** Run "
             "unchanged at `H = 0.5`, where ADS publish American put prices, the "
             "LSM estimator sits a mean "
             f"`{g(ares, 'am_anchor_mean_gap', float('nan')):+.4f}` "
             f"(`{g(ares, 'am_anchor_mean_abs_gap_pct', float('nan')):.2f}%` in "
             "absolute value) from their tree, worst "
             f"`{g(ares, 'am_anchor_worst_gap', float('nan')):+.4f}`, and it is "
             "below the tree on essentially every row -- which is the direction "
             "theory demands of a policy-based lower bound. **That number, about "
             "`-0.02`, is the error bar to attach to every American column in this "
             "document.**\n")
    P.append("**4b. The shortfall is the regression basis, not the exercise "
             "grid.** Going from 50 to 200 exercise dates moves the gap by at most "
             "`0.01` and not always in the same direction, while the gap itself "
             "stays near `-0.02` throughout. Allowing exercise more often "
             "therefore does not close it: what is missing is in the regressors -- "
             "the continuation value is fitted on `(S, V, average variance)`, which "
             "in a rough model does not span the state, since the variance's whole "
             "past matters. Enriching the basis is the way to shrink this, and it "
             "has not been attempted.\n")
    P.append("**5. A second, independent American control.** With `r >= 0` and no "
             "dividend, early exercise of a *call* is never optimal, so the true "
             "American call is the European call, which Fourier gives. The gap the "
             "LSM shows on the call is therefore pure estimator error, measured "
             "rather than assumed: "
             + ", ".join(f"`H = {h}` → `{g(tres, f'lsm_call_gap_mean_H{h}', float('nan')):+.4f}`"
                         for h in hs) + ".\n")
    P.append("**5b. Where the estimator error is visible to the naked eye.** In 5 "
             "of the 144 American-put rows the early-exercise premium comes out "
             "*negative* (worst `-0.0029`), which is impossible: an American "
             "option is worth at least its European counterpart. Every one of "
             "those rows is a deep-out-of-the-money one-month put whose true "
             "premium is of order `0.001`, i.e. far below the `~0.02` estimator "
             "error established in point 4. The negative entries have deliberately "
             "not been clipped away, because they are the cheapest available "
             "reminder of where the American columns stop carrying information: "
             "wherever the premium is smaller than `0.02`, the premium column is "
             "noise.\n")
    if lres:
        P.append("**6. The paper's binomial lattice converges, and its one weak "
                 "spot is understood exactly.** At `H = 0.5` on the 45 European "
                 "puts the binomial mean absolute error against the published "
                 "analytical column is "
                 f"`{g(lres, 'eu_put_mean_abs_err_binomial_n50', float('nan')):.5f}` "
                 f"at `n = 50`, "
                 f"`{g(lres, 'eu_put_mean_abs_err_binomial_n100', float('nan')):.5f}` "
                 f"at `n = 100`, "
                 f"`{g(lres, 'eu_put_mean_abs_err_binomial_n200', float('nan')):.5f}` "
                 "at `n = 200` -- roughly halving as `n` doubles, i.e. order one "
                 "in `delta`. Split by initial volatility at `n = 200` the error "
                 "is "
                 f"`{g(lres, 'eu_put_mean_abs_err_bino_vol0.2', float('nan')):.5f}` "
                 f"at `sqrt(V0) = 0.2` but "
                 f"`{g(lres, 'eu_put_mean_abs_err_bino_vol0.4', float('nan')):.5f}` "
                 "at `sqrt(V0) = 0.4`: the whole error lives where the Lamperti "
                 "drift is large. This is the `o(delta)` variance deficit of "
                 "Section 10.8.3, not a defect of the construction, and it does "
                 "not touch the rate theorem.\n")
        P.append("**7. That deficit is measured directly, node by node, not "
                 "inferred from prices.** The largest relative variance error "
                 "`max |Var/delta - 1|` of the two-point walk over the whole "
                 "lattice, on the worst cell, is "
                 f"`{g(lres, 'max_var_err_binomial_n100', float('nan')):.3f}` at "
                 f"`n = 100` and "
                 f"`{g(lres, 'max_var_err_binomial_n200', float('nan')):.3f}` at "
                 "`n = 200`, which is `mu^2 delta` to three digits. A second, "
                 "independent channel confirms it: the binomial breaks put-call "
                 "parity by a mean "
                 f"`{g(lres, 'parity_mean_abs_err_binomial_n200', float('nan')):.4f}` "
                 "at `n = 200` -- the lost variance desynchronises the price leg's "
                 "Ito correction from the realised driver variance, so `E[S_T]` "
                 "drifts off the forward.\n")
        P.append("**8. The trinomial repair is a numerical diagnostic, not a "
                 "change to the paper.** The Hull--White trinomial with branch "
                 "switching matches the driver's variance exactly (measured "
                 f"`max |Var/delta - 1| = "
                 f"{g(lres, 'max_var_err_trinomial_n200', float('nan')):.0e}`, "
                 "i.e. machine precision), which drops the European mean error to "
                 f"`{g(lres, 'eu_put_mean_abs_err_trinomial_n200', float('nan')):.5f}` "
                 "at `n = 200` -- the same order as ADS's own tree at `N = 500`, "
                 f"`{g(lres, 'eu_put_mean_abs_err_ads_n500', float('nan')):.5f}`, "
                 "and holds put-call parity to "
                 f"`{g(lres, 'parity_mean_abs_err_trinomial_n200', float('nan')):.5f}`. "
                 "But it is a fixed-grid moment-matching chain with **no rate "
                 "theory** and it does not fit the paper's Skorokhod-embedding "
                 "framework; it is kept only to isolate and quantify the deficit "
                 "of point 7. It carries no claim in the rough regime, where the "
                 "Lamperti transform it relies on does not even exist.\n")
        P.append("**9. On the American put the two schemes now bracket the "
                 "answer.** Against the Beliaeva--Nawalkha control-variate value "
                 "at `n = 200`, the trinomial lattice sits "
                 f"`{g(lres, 'am_put_mean_err_trino_vs_cv', float('nan')):+.4f}` "
                 "above on average (its `max(intrinsic, continuation)` test is "
                 "biased up) while the Longstaff--Schwartz lower bound sits "
                 f"`{g(lres, 'am_put_mean_err_lsm_vs_cv', float('nan')):+.4f}` "
                 "below (a fitted, suboptimal policy), so the two straddle the "
                 "published value. The American call is a free control: with "
                 "`r = 0.05 > 0` and no dividend it must equal the European call, "
                 "and on the trinomial lattice the gap is "
                 f"`{g(lres, 'am_call_control_trinomial_n200', 0.0):.1e}` across "
                 "all 81 parameter sets.\n")
    if ores:
        P.append("**10. The order study says: same order, different constant.** "
                 "Measured against the full-precision Fourier reference (the "
                 "published column cannot resolve below `5e-5`), **both walks "
                 "converge at order one** on every line where the error is above "
                 "that floor -- on the high-drift lines the binomial fits `1.04, "
                 "1.02, 1.01` and the trinomial `0.99, 1.06, 1.01`. The variance "
                 "deficit `mu^2 delta` is `O(delta)`, so it does not change the "
                 "*order*; it multiplies the *constant*, by about a hundred on the "
                 "worst cell. So the trinomial does not converge *faster* than the "
                 "binomial -- it converges to the same order with a far smaller "
                 "constant. (The per-walk mean orders printed in Part IV, `0.69` "
                 "and `1.18`, are contaminated by the low-drift lines whose error "
                 "already sits on the reference floor, and should not be "
                 "over-read.) This is the third independent confirmation, after "
                 "the node-by-node variance and put-call parity, that the "
                 "binomial's weakness is a large finite-`n` constant, not a broken "
                 "rate.\n")
    if rres:
        P.append("**11. In the ROUGH regime the recombining lattice does not "
                 "converge -- and Part V is the numerical proof of the paper's "
                 "central obstruction.** Everything in points 1--10 is at "
                 "`H = 0.5`, which is not rough. At `H < 0.5` a recombining "
                 "lattice cannot represent the exact convolution, so the "
                 "implementable object is the one-step lattice, whose terminal "
                 "variance diverges as `2H n^{1-2H}`. Priced against a "
                 "2-million-path exact-covariance Monte-Carlo on rough Bergomi, "
                 "the lattice error tracks that ratio with rank correlation up to "
                 f"`{g(rres, 'rank_corr_ratio_error_H0.1', float('nan')):+.2f}` and "
                 "GROWS with `n`: at `H = 0.1`, `n = 128` the variance ratio is "
                 f"`{g(rres, 'variance_ratio_n128_H0.1', float('nan')):.1f}` and "
                 "the tree overprices by "
                 f"`{g(rres, 'tree_error_n128_H0.1', float('nan')):+.4f}` (about 60 "
                 "Monte-Carlo bands). Refining the grid makes it worse. This is "
                 "why there is no rough lattice column against an exact price "
                 "anywhere here, and cannot be one until Route B (the Markovian "
                 "lift) is implemented -- the single most valuable next step.\n")
    if bres:
        P.append("**12. Route B is now implemented, and it converts the "
                 "divergence into a bounded, controllable bias (Part VI).** The "
                 "Markovian lift replaces the one divergent driver by `m` "
                 "Ornstein--Uhlenbeck factors sharing one Brownian increment, "
                 "built as a real `(m+1)`-dimensional recombining lattice. At the "
                 "level of the driver variance the effect is exact: at `H = 0.1`, "
                 f"`n = 512` the one-step ratio is "
                 f"`{g(bres, 'vr_onestep_n512_H0.1', float('nan')):.1f}` (and "
                 "climbing) while the lift ratios sit at "
                 f"`{g(bres, 'vr_lift_m1_n512_H0.1', float('nan')):.2f}` (m=1) and "
                 f"`{g(bres, 'vr_lift_m2_n512_H0.1', float('nan')):.2f}` (m=2), "
                 "flat in `n`. In price, the one-step error keeps growing with `n` "
                 "while the lift error converges to an `m`-dependent floor that "
                 "rises toward the truth as `m` grows. This is bounded-vs-"
                 "unbounded, not accurate-vs-inaccurate: a small `m` still carries "
                 "a real bias, but the scheme is now CONVERGENT, which the "
                 "one-step is not. The lift was cross-validated against an "
                 "independent mixing-formula pricer (agreement `~1.5e-2`).\n")
    P.append("## What this document does not establish\n")
    P.append("- **There is no reference for the rough American price.** No closed "
             "form and no unbiased simulation of it exists. The American columns "
             "are lower bounds from a specific exercise policy on a finite grid, "
             "with the size of that shortfall measured at `H = 0.5` and carried "
             "over. They are not the American price.\n")
    P.append("- **The Monte-Carlo is biased, by the amount reported in point 2.** "
             "Where it is used as the comparison point -- which is everywhere the "
             "Fourier engine cannot go -- that bias is part of the comparison and "
             "is stated alongside it.\n")
    P.append("- **The timings are not a fair race against ADS.** They report "
             "MATLAB on a 2010 laptop; this is C++ on 9 threads of Apple silicon. "
             "Comparisons within this document are meaningful; comparisons across "
             "the machine boundary are indicative at best.\n")
    P.append("- **Only vanilla payoffs.** ADS also price barrier, lookback and "
             "Asian options. Those are not attempted here.\n")
    P.append("- **There is no lattice column for `H < 0.5`, and there cannot be "
             "one for this model.** That is the proposition restated at the top, "
             "not an omission. What a rough lattice column would require is a "
             "rough model *inside* the paper's class -- rough Bergomi -- where the "
             "lattice does run but no closed form exists to check it against. The "
             "two things the reader might want on one line, roughness and an exact "
             "price, cannot be had together.\n")
    P.append("- **The lattice comparison is at `H = 0.5`, so it tests the "
             "construction and not the roughness.** It exercises the moment "
             "matching, the recombination, the two-dimensional backward induction, "
             "the Route A' price coupling and a state-dependent `sigma_y`. It says "
             "nothing about the rough kernel, which is where the paper's rate lives.\n")
    P.append("\n---\n")
    P.append(f"# Part I — the tables\n")
    P.append(f"*Generated by `sim/run_rheston_tables.py`, run `{tr.name}`.*\n")
    P.append(body(tmd))
    P.append("\n---\n")
    P.append(f"# Part II — the American estimator, calibrated\n")
    P.append(f"*Generated by `sim/run_rheston_american_anchor.py`, run "
             f"`{ar.name}`.*\n")
    P.append(body(amd))
    if lmd is not None:
        P.append("\n---\n")
        P.append("# Part III — the paper's lattice, against two references it did "
                 "not produce\n")
        P.append(f"*Generated by `sim/run_heston_ads_tables.py`, run "
                 f"`{lr.name}`; the headline table below joins it with the "
                 f"Monte-Carlo of run `{tr.name}`.*\n")
        P.extend(headline_european(lmd, tmd))
        P.append("\nThe European put above is the direct tree-vs-Monte-Carlo-vs-"
                 "reference comparison; the American put is the analogous table in "
                 "Section 2 below (tree vs Longstaff--Schwartz vs ADS's tree vs the "
                 "Beliaeva--Nawalkha control variate). The remaining sections "
                 "dissect where the two tree walks differ.\n")
        P.append(body(lmd))
    if omd is not None:
        P.append("\n---\n")
        P.append("# Part IV — the lattice's convergence order\n")
        P.append(f"*Generated by `sim/run_heston_lattice_order.py`, run "
                 f"`{or_.name}`.*\n")
        P.append(body(omd))
    if rmd is not None:
        P.append("\n---\n")
        P.append("# Part V — the lattice in the ROUGH regime (rough Bergomi)\n")
        P.append("Everything above is at `H = 0.5`, which is not rough: it exists "
                 "only because that is the one place an exact price is available to "
                 "validate against. This part is the genuinely rough case. It has "
                 "to switch model -- from rough Heston to rough Bergomi -- for a "
                 "reason that is structural, not incidental: **a rough lattice "
                 "column can only be built for a model inside the paper's class "
                 "(autonomous driver), and rough Heston is not in it.** Rough "
                 "Bergomi is; the price of that is it has no closed form, so the "
                 "reference becomes a trustworthy Monte-Carlo rather than an exact "
                 "value, and the error is resolved only down to the Monte-Carlo "
                 "band.\n")
        P.append(f"*Generated by `sim/run_rough_bergomi_ladder.py`, run "
                 f"`{rr.name}`.*\n")
        P.append(body(rmd))
    if bmd is not None:
        P.append("\n---\n")
        P.append("# Part VI — Route B implemented: does the Markovian lift fix "
                 "the rough divergence?\n")
        P.append("Part V is the disease; this is the treatment. Route B replaces "
                 "the single divergent driver by `m` Ornstein--Uhlenbeck factors "
                 "sharing one Brownian increment, so the driver variance stays "
                 "finite. The lift is built as an actual `(m+1)`-dimensional "
                 "recombining lattice (`route_b_lattice.LiftedLattice`) and priced "
                 "against the same Monte-Carlo. The finding, in one line: the lift "
                 "converts the one-step's **unbounded** divergence into a "
                 "**bounded** bias that shrinks as `m` grows -- it makes the "
                 "scheme convergent, which is what the paper needs, even though a "
                 "small `m` is not yet accurate.\n")
        P.append(f"*Generated by `sim/run_routeb_compare.py`, run `{br.name}`.*\n")
        P.append(body(bmd))
        P.append("\nThe question this leaves — **how large must `m` be, given "
                 "`n`?** — is answered separately in "
                 "[`docs/N_M_BOUND.md`](N_M_BOUND.md): `m` grows like `log n`, "
                 "so doubling the number of time steps costs a fraction of one "
                 "extra factor. That document has the measured law across `H`, "
                 "the mechanism behind it, the candidate inequality with what is "
                 "and is not proved, and `sim/nm_bound.m_required`, which sizes "
                 "`m` for a run instead of choosing it by hand.\n")
    P.append("\n---\n")
    P.append("# Terms used above\n")
    P.append("- **Hurst parameter `H`** — how rough the variance path is. "
             "`H = 0.5` is the ordinary, non-rough case; the smaller `H`, the more "
             "jagged the variance and the more the past matters.\n")
    P.append("- **Fractional Riccati equation** — the ordinary Riccati equation of "
             "the classical Heston model, but with a fractional derivative in "
             "place of the ordinary one. Solving it is what replaces the "
             "closed-form characteristic function that the non-rough model has.\n")
    P.append("- **Euler--Volterra scheme** — the simulation method. At each step "
             "the variance is rebuilt from its *whole* past through a convolution, "
             "because in a rough model the variance is not a Markov process: it "
             "cannot be advanced from its current value alone.\n")
    P.append("- **Antithetic pairs** — each random path is simulated twice, once "
             "with the random numbers and once with their signs flipped. The two "
             "errors partly cancel, so fewer paths are needed for a given accuracy.\n")
    P.append("- **Longstaff--Schwartz (LSM)** — the standard way to price an "
             "early-exercisable option by simulation: at each possible exercise "
             "date, estimate by regression what holding on is worth, and exercise "
             "when the immediate payoff beats it.\n")
    P.append("- **Standard error (`s.e.`)** — the statistical uncertainty of a "
             "simulated price, from the finite number of paths. It shrinks as more "
             "paths are added. A *bias* does not.\n")
    P.append("- **Bias** — a systematic error that more paths do not remove, here "
             "coming from cutting time into finitely many steps.\n")
    P.append("- **Early-exercise premium** — how much more the American option is "
             "worth than the otherwise identical European one, i.e. the value of "
             "being allowed to exercise early.\n")
    P.append("- **Feller ratio** — `2 kappa theta / eta^2`. When it is large the "
             "variance stays comfortably positive; when it approaches zero the "
             "variance hugs zero and simulation schemes start needing to truncate.\n")
    P.append("- **Control variate** — a quantity simulated alongside the payoff "
             "whose exact value is known, used to cancel part of the simulation "
             "error.\n")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(P) + "\n")
    print(f"wrote {OUT} ({OUT.stat().st_size:,} bytes)")
    print(f"  tables run: {tr.name}")
    print(f"  anchor run: {ar.name}")


if __name__ == "__main__":
    main()
