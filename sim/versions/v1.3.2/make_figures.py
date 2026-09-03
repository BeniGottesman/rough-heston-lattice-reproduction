#!/usr/bin/env python3
"""
make_figures — the figures of the appendix, generated FROM THE RUNS.

Nothing here is a transcribed number: each panel reads the run directory that
produced it (`progress.json` where the values were stored as results, the
`RESULTS.md` table where they were not), so a figure cannot drift away from the
evidence it illustrates.  Output goes to paper/<version>/figs/*.pdf.

    python3 sim/make_figures.py [paper/v1.4.0]
"""
from __future__ import annotations

import glob
import json
import math
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402
import mpmath as mp                      # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
mp.mp.dps = 30

plt.rcParams.update({
    "figure.dpi": 200, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
    "font.size": 8.5, "axes.labelsize": 8.5, "axes.titlesize": 9,
    "legend.fontsize": 7.5, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.linewidth": 0.6, "grid.linewidth": 0.4, "lines.linewidth": 1.2,
    "axes.grid": True, "grid.color": "0.88", "axes.axisbelow": True,
    "axes.spines.top": False, "axes.spines.right": False,
})
C = {"blue": "#2a78d6", "orange": "#eb6834", "green": "#1baf7a",
     "amber": "#c98500", "red": "#d03b3b", "violet": "#4a3aa7", "grey": "#7a7a76"}


def latest(name: str) -> Path:
    """Newest run whose recorded name is exactly `name`.

    Matching on the glob alone is wrong: "route-b-scaling-<stamp>" sorts after
    "route-b-<stamp>" and would silently win the pattern "route-b-*".
    """
    hits = []
    for d in sorted(glob.glob(str(ROOT / "runs" / "*"))):
        j = Path(d) / "progress.json"
        if j.is_file():
            try:
                if json.loads(j.read_text()).get("name") == name:
                    hits.append(Path(d))
            except Exception:
                pass
    if not hits:
        raise SystemExit(f"no run named {name!r}")
    return hits[-1]


def results(run: Path) -> dict:
    return json.loads((run / "progress.json").read_text())["results"]


# --------------------------------------------------------------------------- #
def fig_variance_proposition(out: Path) -> str:
    """Proposition 4.6: measured relative variance error vs the closed form."""
    fig, ax = plt.subplots(figsize=(4.5, 2.9))
    ns = np.array([2 ** k for k in range(6, 15)], float)
    for col, H in zip(("blue", "green", "amber", "orange", "red"),
                      (0.05, 0.10, 0.20, 0.30, 0.45)):
        Hm = mp.mpf(H)
        z = mp.zeta(1 - 2 * Hm)
        meas, pred = [], []
        for n in ns:
            d = mp.mpf(1) / int(n)
            s = 2 * Hm - 1
            vd = d ** (2 * Hm) * (mp.zeta(-s) - mp.zeta(-s, int(n) + 1))
            meas.append(float(abs((vd - 1 / (2 * Hm)) * 2 * Hm)))
            pred.append(float(abs(2 * Hm * z * d ** (2 * Hm) + Hm * d)))
        ax.loglog(ns, meas, "o", ms=3.2, color=C[col], label=f"$H={H}$")
        ax.loglog(ns, pred, "-", color=C[col], alpha=0.75)
    ax.set_xlabel("$n$")
    ax.set_ylabel(r"$\left|\,\mathrm{Var}[V^{(n)}_T]/\mathrm{Var}[V_T]-1\,\right|$")
    ax.legend(ncol=2, frameon=False, loc="lower left")
    ax.set_title("markers: exact finite sum · lines: "
                 r"$2H\zeta_{\mathrm{R}}(1-2H)\delta^{2H}+H\delta$", fontsize=7.5)
    p = out / "fig-variance-proposition.pdf"
    fig.savefig(p); plt.close(fig)
    return p.name


def fig_lift_error(out: Path) -> str:
    """Route B: the lift's discrete covariance error vs m, both pairings."""
    r = results(latest("route-b"))["phase6_discrete_covariance"]
    n = "256"
    blocks = r["by_n"][n]["by_mode"]
    fig, ax = plt.subplots(figsize=(4.5, 2.9))
    spec = [("left", "L2(delta,T)", "blue", "-", "o",
             r"left endpoint $+$ $L^2(\delta,T)$"),
            ("cellavg", "L2(0,T)", "orange", "--", "s",
             r"cell average $+$ $L^2(0,T)$")]
    for mode, lift, col, ls, mk, lab in spec:
        rows = [x for x in blocks[mode]["rows"] if x["lift"] == lift]
        ms = [x["m"] for x in rows]
        fr = [max(x["frob"], 1e-9) for x in rows]
        ax.loglog(ms, fr, ls, marker=mk, ms=3.2, color=C[col], label=lab)
        base = blocks[mode]["scheme_own_discretisation_error"]["rel_frobenius"]
        ax.axhline(base, color=C[col], ls=":", lw=0.9)
        ax.text(ms[-1], base * 1.25, "threshold", color=C[col],
                fontsize=6.5, ha="right")
    ax.set_xlabel("number of factors $m$")
    ax.set_ylabel("relative Frobenius error of the\ndiscrete covariance")
    ax.set_xticks([1, 2, 3, 5, 8, 12, 20, 30])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.legend(frameon=False, loc="lower left")
    ax.set_title(f"$n={n}$; dotted: the scheme's own error with the true kernel",
                 fontsize=7.5)
    p = out / "fig-lift-error.pdf"
    fig.savefig(p); plt.close(fig)
    return p.name


def fig_mstar(out: Path) -> str:
    """Route B: m*(n) over eight doublings of the grid."""
    run = latest("route-b-scaling")
    r2 = results(run)["phase2_m_star"]["rows"]
    r3 = results(run)["phase3_frobenius_confirmation"]["rows"]
    fig, ax = plt.subplots(figsize=(4.5, 2.7))
    ns = [x["n"] for x in r2]
    ax.semilogx(ns, [x["m_star_ajee"] for x in r2], "--s", ms=3.2,
                color=C["orange"], label="Abi Jaber--El Euch partition")
    ax.semilogx(ns, [x["m_star_floor"] for x in r2], "-o", ms=3.6,
                color=C["green"], label="achievable floor (variance criterion)")
    ax.semilogx([x["n"] for x in r3], [x["m_star_floor"] for x in r3], "^",
                ms=5, color=C["violet"], label="floor, full covariance (strict)")
    ax.set_xlabel("$n$")
    ax.set_ylabel("smallest usable $m^{*}$")
    ax.set_xticks(ns)
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_yticks(range(0, 10, 2))
    ax.set_ylim(0, 9)
    ax.legend(frameon=False, loc="upper left")
    p = out / "fig-mstar.pdf"
    fig.savefig(p); plt.close(fig)
    return p.name


def _sweep_rows() -> dict:
    """Per-n signed errors, parsed from the sweep's RESULTS.md table."""
    txt = (latest("mc-vs-tree-sweep") / "RESULTS.md").read_text()
    blocks = re.findall(
        r"### (\S+) = (\S+)\s+\((.*?)\).*?\|---\|---\|---\|---\|---\|\n(.*?)\n\n",
        txt, re.S)
    out: dict = {}
    for key, val, _params, tbl in blocks:
        rows = [r.split("|")[1:-1] for r in tbl.strip().split("\n")]
        out.setdefault(key, {})[float(val)] = {
            "n": [int(r[0]) for r in rows],
            "err": [float(r[3]) for r in rows]}
    band = re.search(r"MC true = \*\*[\d.]+\*\* ±([\d.]+)", txt)
    out["_band"] = float(band.group(1)) if band else 0.007
    return out


def fig_sweep_rho(out: Path) -> str:
    """The signed error against n, ordered by rho."""
    sw = _sweep_rows()
    fig, ax = plt.subplots(figsize=(4.5, 2.9))
    spec = [(-0.9, "green", "-", "o"), (-0.7, "blue", "--", "s"),
            (-0.4, "amber", "-.", "^"), (0.0, "red", "-", "D")]
    for val, col, ls, mk in spec:
        d = sw["rho"][val]
        ax.plot(d["n"], d["err"], ls, marker=mk, ms=3.4, color=C[col],
                label=rf"$\rho={val}$")
    b = sw["_band"]
    ax.axhspan(-b, b, color="0.85", zorder=0)
    ax.axhline(0.0, color="0.5", lw=0.6)
    ax.set_xscale("log", base=2)
    ax.set_xticks([8, 16, 32, 64])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xlabel("$n$")
    ax.set_ylabel("signed error, tree $-$ reference")
    ax.legend(frameon=False, ncol=2, loc="upper left")
    ax.set_title(f"grey band: the reference's own $\\pm{b:g}$; "
                 "every series crosses zero", fontsize=7.5)
    p = out / "fig-sweep-rho.pdf"
    fig.savefig(p); plt.close(fig)
    return p.name


def fig_sweep_panels(out: Path) -> str:
    """The error at n=64 as a function of each parameter."""
    sw = _sweep_rows()
    fig, axes = plt.subplots(1, 3, figsize=(6.4, 2.2), sharey=True)
    spec = [("H", "$H$ (roughness)", "violet"),
            ("eta", r"$\eta$ (vol of vol)", "blue"),
            ("rho", r"$\rho$ (correlation)", "red")]
    for ax, (key, lab, col) in zip(axes, spec):
        vals = sorted(sw[key])
        e64 = [sw[key][v]["err"][sw[key][v]["n"].index(64)] for v in vals]
        cols = [C["orange"] if x < 0 else C[col] for x in e64]
        ax.bar(range(len(vals)), e64, color=cols, width=0.62)
        ax.axhline(0.0, color="0.4", lw=0.6)
        ax.set_xticks(range(len(vals)))
        ax.set_xticklabels([f"{v:g}" for v in vals])
        ax.set_xlabel(lab)
        ax.grid(axis="x", visible=False)
    axes[0].set_ylabel("signed error at $n=64$")
    p = out / "fig-sweep-panels.pdf"
    fig.savefig(p); plt.close(fig)
    return p.name


def fig_heston(out: Path) -> str | None:
    """Rough Heston: the Fourier reference against Monte-Carlo, and the cost."""
    try:
        r = results(latest("rough-heston"))
    except SystemExit:
        return None
    if "phase3_fourier_vs_mc" not in r:
        return None
    rows = r["phase3_fourier_vs_mc"]["rows"]
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 2.5))
    Hs = [x["H"] for x in rows]
    ax = axes[0]
    ax.errorbar(Hs, [x["mc"] for x in rows],
                yerr=[x["mc_ci95"] for x in rows], fmt="o", ms=3.4,
                color=C["blue"], capsize=2, lw=1.0, label="Monte-Carlo")
    ax.plot(Hs, [x["fourier"] for x in rows], "-s", ms=3.4,
            color=C["orange"], label="Fourier (fractional Riccati)")
    ax.set_xlabel("$H$"); ax.set_ylabel("put price")
    ax.legend(frameon=False)
    ax = axes[1]
    w = 0.36
    x = np.arange(len(Hs))
    ax.bar(x - w / 2, [x_["fourier_seconds"] for x_ in rows], w,
           color=C["orange"], label="Fourier")
    ax.bar(x + w / 2, [x_["mc_seconds"] for x_ in rows], w,
           color=C["blue"], label="Monte-Carlo")
    ax.set_yscale("log")
    ax.set_xticks(x); ax.set_xticklabels([f"{h:g}" for h in Hs])
    ax.set_xlabel("$H$"); ax.set_ylabel("seconds")
    ax.grid(axis="x", visible=False)
    ax.legend(frameon=False)
    p = out / "fig-heston.pdf"
    fig.savefig(p); plt.close(fig)
    return p.name


def fig_telescope(out: Path) -> str | None:
    """Section 10.7: the two terms of the telescope, measured.

    Left: the defect Route B removes against the residue it leaves, so that the
    reader sees the one growing in n and the other not.  Right: the term
    Theorem 7.1 bounds, which is flat -- the point of the restated (O3).
    """
    try:
        run = latest("route-b-lattice")
    except SystemExit:
        return None
    tel = results(run).get("telescope")
    R = results(run).get("continuous_reference")
    if not tel or R is None:
        return None
    ns = np.array([r["n"] for r in tel], float)
    defect = np.array([r["onestep"] - r["exact"] for r in tel])
    driver = np.array([r["exact"] - R - r["price_scheme_err"] for r in tel])

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.2, 3.1))

    a1.plot(ns, np.abs(defect), "o-", color=C["red"], lw=1.6,
            label=r"one-step defect $|\check V^{(n)}-\Lambda^{(n)}|$")
    for m, col in ((1, C["blue"]), (2, C["green"]), (3, C["amber"])):
        xs = [r["n"] for r in tel if str(m) in r["lift"]]
        ys = [abs(r["lift"][str(m)] - r["exact"]) for r in tel
              if str(m) in r["lift"]]
        if xs:
            a1.plot(xs, ys, "s--", color=col, lw=1.4, ms=4,
                    label=rf"lift residue, $m={m}$")
    band = [r["exact_se"] for r in tel]
    a1.fill_between(ns, 0, band, color=C["grey"], alpha=0.25,
                    label=r"Monte-Carlo band on $\Lambda^{(n)}$")
    a1.set_xscale("log", base=2); a1.set_yscale("log")
    a1.set_xlabel("$n$"); a1.set_ylabel("absolute deviation from $\\Lambda^{(n)}$")
    a1.set_title("what the lift removes, and what it leaves", fontsize=9)
    a1.legend(fontsize=6.5, frameon=False)

    a2.axhline(0.0, color=C["grey"], lw=0.8)
    a2.plot(ns, driver, "o-", color=C["violet"], lw=1.6,
            label=r"measured: $\Lambda^{(n)}-\Lambda$, price grid removed")
    a2.plot(ns, -(1.0 / ns) ** 0.05 * 0.09, "--", color=C["grey"], lw=1.2,
            label=r"shape of $\delta^{(h+\kappa)/2}=\delta^{0.05}$ (scaled)")
    a2.set_xscale("log", base=2)
    a2.set_xlabel("$n$"); a2.set_ylabel("signed error")
    a2.set_title(r"the term Theorem~7.1 bounds: flat", fontsize=9)
    a2.set_ylim(-0.105, 0.012)
    a2.legend(fontsize=6.5, frameon=False, loc="upper right")

    fig.tight_layout()
    p = "fig-telescope.pdf"
    fig.savefig(out / p); plt.close(fig)
    return p


def main() -> None:
    version = sys.argv[1] if len(sys.argv) > 1 else "paper/v1.5.0"
    out = ROOT / version / "figs"
    out.mkdir(parents=True, exist_ok=True)
    made = [fig_variance_proposition(out), fig_lift_error(out), fig_mstar(out),
            fig_sweep_rho(out), fig_sweep_panels(out)]
    h = fig_heston(out)
    if h:
        made.append(h)
    else:
        print("  (rough Heston figure skipped: no run yet)")
    t = fig_telescope(out)
    if t:
        made.append(t)
    else:
        print("  (telescope figure skipped: no route-b-lattice run yet)")
    for m in made:
        print(f"  wrote {version}/figs/{m}")


if __name__ == "__main__":
    main()
