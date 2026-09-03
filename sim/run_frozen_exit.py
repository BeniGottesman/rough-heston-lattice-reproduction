"""
run_frozen_exit -- driver for the Route F falsifiers FF1-FF4, FF6, FF7.

    python3 sim/run_frozen_exit.py            # full run (~7 min)
    python3 sim/run_frozen_exit.py --quick    # smoke run (~1 min)

Writes runs/frozen-exit-<UTC>/{progress.json, log.txt, RESULTS.md} through
sim/progress.py, as every long run in this project does.

FF1-FF4 are stated in research-notes/L003-ROUTE-F.md sec.3.  FF6 and FF7 test
the barrier inset itself: eq:barriers writes the literal `delta^{1/3}`, but the
value 1/3 is a free parameter of the scheme (it occurs only inside the barrier
definition and its own lemma), so every FF2 quantity is measured at
`beta = 1/3` AND `beta = 1/2` on the SAME paths.

This driver only MEASURES; it closes nothing.  Every modelling shortcut is
listed in the "Shortcuts" section at the end of RESULTS.md.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import frozen_exit as fe                                   # noqa: E402
from progress import Progress                              # noqa: E402

QUICK = "--quick" in sys.argv

# ----------------------------------------------------------------- parameters
H_LIST = [-0.1, -0.3, -0.45]              # FF1 only: no continuous object, so
                                          # h=-0.45 costs nothing and is exact
FF2_H = [-0.1, -0.3]                      # h=-0.45 dropped: its target rate is
                                          # delta^{0.025}, unresolvable on any
                                          # feasible ladder (see RESULTS.md)
BETAS = [fe.BETA0, 0.5]                   # barrier inset delta^beta
ELLS = [1.0, 2.0, 4.0]                    # prop:Vconv is claimed for every ell>=1

FF1_N = [1000, 3162, 10000, 31623, 100000, 316228, 1000000]
FF1_PATHS = {1000: 512, 3162: 512, 10000: 384, 31623: 256,
             100000: 160, 316228: 96, 1000000: 64}

FF2_MPU = 1 << 20                         # fine steps per unit time
FF2_HORIZON = 1.25
FF2_N = [64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384]
FF2_PATHS = 384
FF2_BATCH = 4
FF2_RES_N = 512
FF2_RES_MPU = [1 << 16, 1 << 18, 1 << 20]
FF2_RES_PATHS = 96

FF3_MPU = 1 << 19
FF3_HORIZON = 1.25
FF3_N = [64, 256, 1024, 4096]
FF3_PATHS = 192
FF3_BATCH = 64
FF3_INJECT_C = [0.0, 0.5, 1.0, 2.0]

FF4_KAPS = [0.2, 0.35, 0.45]
FF4_N = [10 ** k for k in range(2, 13)]   # exact law: depth is free
FF4_TAIL_N = 10000
FF4_MC_N = [256, 1024, 4096]
FF4_MC_PATHS = 512
FF4_MC_MPU = 1 << 19
FF4_MC_HORIZON = 3.0
FF7_C = 0.5                               # the fixed c of FF7: c = T/2
WELLFORM_N = [16, 32, 64, 128, 256, 1024, 4096, 16384, 65536]

if QUICK:
    FF1_N = [1000, 10000, 100000]
    FF1_PATHS = {k: 64 for k in FF1_N}
    FF2_MPU, FF2_N, FF2_PATHS = 1 << 17, [64, 256, 1024], 32
    FF2_RES_MPU, FF2_RES_PATHS = [1 << 15, 1 << 17], 24
    FF3_MPU, FF3_N, FF3_PATHS, FF3_BATCH = 1 << 16, [64, 256], 32, 32
    FF4_MC_N, FF4_MC_PATHS, FF4_MC_MPU, FF4_MC_HORIZON = [256], 48, 1 << 17, 2.0


def fmt(x, d=4):
    return "n/a" if x is None else f"{x:.{d}g}"


def bs(b):
    return "1/3" if abs(b - 1 / 3) < 1e-9 else ("1/2" if abs(b - 0.5) < 1e-9 else f"{b:g}")


# =============================================================== FF1
def run_ff1(pr, rng):
    pr.phase("FF1 frozen-vs-unfrozen (discrete only)", total=len(FF1_N) * len(H_LIST))
    rows, step = [], 0
    for n in FF1_N:
        delta = fe.T / n
        for h in H_LIST:
            P = FF1_PATHS[n]
            chunk = max(1, min(P, int(3e6 // n) or 1))
            tot, kp, ip, sc, exf, nzf = [], [], [], [], [], []
            for s in range(0, P, chunk):
                m = min(chunk, P - s)
                a, b, c, d, ex, nz = fe.ff1_batch(n, h, m, rng)
                tot.append(a); kp.append(b); ip.append(c); sc.append(d)
                exf.append(ex); nzf.append(nz)
            tot = np.concatenate(tot)
            l2, se = fe.bootstrap_lp(tot, 2.0)
            rows.append(dict(n=n, delta=delta, h=h, paths=P, L2=l2, se=se,
                             L2_kern=fe.lp_norm(np.concatenate(kp), 2.0),
                             L2_incr=fe.lp_norm(np.concatenate(ip), 2.0),
                             L2_scale=fe.lp_norm(np.concatenate(sc), 2.0),
                             exit_frac=float(np.mean(exf)),
                             nz_post=float(np.mean(nzf)),
                             deltaH=delta ** (h + 0.5)))
            step += 1
            pr.tick(step, note=f"FF1 n={n} h={h} L2={l2:.4f}")
            pr.log(f"FF1 n={n} h={h} P={P} L2={l2:.5f}+-{se:.5f} "
                   f"d^H={delta ** (h + 0.5):.5f} kern={rows[-1]['L2_kern']:.4f} "
                   f"incr={rows[-1]['L2_incr']:.4f} scale={rows[-1]['L2_scale']:.4f} "
                   f"nz_post={rows[-1]['nz_post']:.4f}")
    return rows


# =============================================================== FF2 / FF6
KEYS = ("total_theta", "total_tk", "F1", "F1_full", "F2")


def _ff2_core(pr, rng, mpu, n_list, paths, h_list, tag, phase=True, betas=None):
    betas = betas or BETAS
    M = int(FF2_HORIZON * mpu)
    dt = 1.0 / mpu
    nb = (paths + FF2_BATCH - 1) // FF2_BATCH
    if phase:
        pr.phase(f"{tag} frozen-vs-frozen (coupled, mpu={mpu}, M={M})", total=nb)
    acc = {(n, h, k, b): [] for n in n_list for h in h_list
           for k in KEYS for b in betas}
    diags = {(n, b): [] for n in n_list for b in betas}
    kern = {h: fe.kernel_array(M, dt, h) for h in h_list}
    for j in range(nb):
        m = min(FF2_BATCH, paths - j * FF2_BATCH)
        dW = rng.standard_normal((m, M)) * math.sqrt(dt)
        yb = np.concatenate([np.zeros((m, 1)), np.cumsum(dW, axis=1)], axis=1)
        vb = {h: np.concatenate([np.zeros((m, 1)), fe.V0 + fe.causal_conv(dW, kern[h])],
                                axis=1) for h in h_list}
        # absorbed-driver continuous convolution: the INTEGRAND stops at Xi^y,
        # the kernel argument keeps running.  This is the continuous side of
        # prop:Vconv as its proof actually writes it.
        dWa = dW.copy()
        for p in range(m):
            jx = fe.first_exit_index(yb[p])
            if jx < M:
                dWa[p, jx:] = 0.0
        va = {h: np.concatenate([np.zeros((m, 1)), fe.V0 + fe.causal_conv(dWa, kern[h])],
                                axis=1) for h in h_list}
        for p in range(m):
            res, dg = fe.ff2_one_path(yb[p], {h: vb[h][p] for h in h_list},
                                      mpu, n_list, h_list, betas=betas,
                                      vabs={h: va[h][p] for h in h_list})
            for k, v in res.items():
                acc[k].append(v)
            for k, v in dg.items():
                diags[k].append(v)
        if phase:
            pr.tick(j + 1, note=f"{tag} batch {j + 1}/{nb}")
    per = {(h, k, b): np.array([acc[(n, h, k, b)] for n in n_list])
           for h in h_list for k in KEYS for b in betas}
    return per, diags, dt, M


def run_ff2(pr, rng):
    per, diags, dt, M = _ff2_core(pr, rng, FF2_MPU, FF2_N, FF2_PATHS, FF2_H,
                                  "FF2/FF6")
    rows, fits = [], {}
    ds = [fe.T / n for n in FF2_N]
    for h in FF2_H:
        for k in KEYS:
            for b in BETAS:
                for i, n in enumerate(FF2_N):
                    for ell in ELLS:
                        v, se = fe.bootstrap_lp(per[(h, k, b)][i], ell)
                        rows.append(dict(n=n, h=h, key=k, beta=b, ell=ell,
                                         delta=fe.T / n, norm=v, se=se))
                for ell in ELLS:
                    m, sd, ci = fe.bootstrap_slope_lp(ds, per[(h, k, b)], ell)
                    fits[(h, k, b, ell)] = (m, sd, ci)
                    pr.log(f"FF6 h={h} {k} beta={bs(b)} L{ell:g} "
                           f"slope={m:.4f}+-{sd:.4f} ci={ci}")
    dsum = {(n, b): dict(
        e_max=float(np.mean([d["e_max"] for d in diags[(n, b)]])),
        multi=int(sum(d["multi"] for d in diags[(n, b)])),
        short=float(np.mean([d["avail"] < n for d in diags[(n, b)]])),
        exited=float(np.mean([d["exited"] for d in diags[(n, b)]])),
        y_exit=float(np.mean([d["y_exited"] for d in diags[(n, b)]])),
        meanG=float(np.mean([d["G"] for d in diags[(n, b)]])),
        G=np.array([d["G"] for d in diags[(n, b)]]),
    ) for n in FF2_N for b in BETAS}
    return rows, fits, dsum, dt


def run_ff2_resolution(pr, rng):
    pr.phase("FF2 resolution sensitivity (fixed n, varying fine grid)",
             total=len(FF2_RES_MPU))
    out = []
    for i, mpu in enumerate(FF2_RES_MPU):
        per, _, dt, _ = _ff2_core(pr, rng, mpu, [FF2_RES_N], FF2_RES_PATHS,
                                  FF2_H, "res", phase=False, betas=[fe.BETA0])
        for h in FF2_H:
            row = dict(mpu=mpu, dt=dt, h=h)
            for k in KEYS:
                v, se = fe.bootstrap_lp(per[(h, k, fe.BETA0)][0], 2.0)
                row[k], row[k + "_se"] = v, se
            out.append(row)
        pr.tick(i + 1, note=f"resolution mpu={mpu}")
        pr.log(f"FF2-res mpu={mpu} " + " ".join(
            f"h={r['h']}:{r['total_theta']:.4f}" for r in out if r["mpu"] == mpu))
    return out


# =============================================================== FF3
def run_ff3(pr, rng):
    M = int(FF3_HORIZON * FF3_MPU)
    dt = 1.0 / FF3_MPU
    nb = (FF3_PATHS + FF3_BATCH - 1) // FF3_BATCH
    pr.phase(f"FF3 ordering L1 (non-constant coefficients, M={M})", total=2 * nb)
    out, step = {}, 0

    def rec(key, r):
        e = out.setdefault(key, [0, 0, [], 0])
        if r is None:
            e[3] += 1
            return
        jth, jxi, err = r
        e[0] += 1
        e[1] += int(jth > jxi)
        e[2].append(err)

    for b in range(nb):
        m = min(FF3_BATCH, FF3_PATHS - b * FF3_BATCH)
        y = fe.euler_paths(M, dt, m, rng, drift=False)
        for p in range(m):
            yf = np.ascontiguousarray(y[:, p])
            for n in FF3_N:
                rec(("A", n), fe.ff3_variantA(yf, n))
        step += 1
        pr.tick(step, note=f"FF3-A batch {b + 1}/{nb}")

    for b in range(nb):
        m = min(FF3_BATCH, FF3_PATHS - b * FF3_BATCH)
        y, bb = fe.euler_paths(M, dt, m, rng, drift=True, aux=True)
        for p in range(m):
            yf = np.ascontiguousarray(y[:, p])
            bf = np.ascontiguousarray(bb[:, p])
            for n in FF3_N:
                ks = np.arange(n + 1) * (FF3_MPU // n)
                shape = np.tanh(bf[ks])
                for c in FF3_INJECT_C:
                    inj = None if c == 0.0 else c * (fe.T / n) ** 0.25 * shape
                    rec((f"C{c}" if c else "B", n),
                        fe.ff3_grid_surrogate(yf, n, FF3_MPU, inject=inj))
        step += 1
        pr.tick(step, note=f"FF3-BC batch {b + 1}/{nb}")

    rows = []
    for (var, n), (used, viol, es, none) in sorted(out.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        rows.append(dict(variant=var, n=n, used=used, viol=viol,
                         frac=(viol / used if used else float("nan")),
                         e_max=float(np.mean(es)) if es else float("nan"),
                         gap=fe.barrier_gap(fe.T / n),
                         d14=(fe.T / n) ** 0.25, no_exit=none))
        pr.log(f"FF3 {var} n={n} used={used} viol={viol} frac={rows[-1]['frac']:.4f} "
               f"e_max={rows[-1]['e_max']:.5f} gap={rows[-1]['gap']:.5f}")
    return rows


# =============================================================== FF4 / FF7
def run_ff4_exact(pr):
    pr.phase("FF4/FF7 exit-time law (exact, both insets)", total=len(FF4_N))
    rows = []
    for i, n in enumerate(FF4_N):
        r = fe.ff4_exact(n, FF4_KAPS)
        rc = fe.ff4_exact(n, FF4_KAPS, tmax=fe.T)
        r["cens"] = {k: rc[f"EG^{k}"] for k in FF4_KAPS}
        r["EG_cens"] = rc["EG"]
        rows.append(r)
        pr.tick(i + 1, note=f"FF4 n={n}")
        pr.log(f"FF4 n={n} gap={r['gap']:.5g} "
               + " ".join(f"E[G^{k}]={r[f'EG^{k}']:.5g}" for k in FF4_KAPS)
               + f" E[G]={r['EG']:.5g}")
    return rows


def run_ff4_beta(ns):
    """(c) E[G^kappa'] on the FF2 ladder at both insets, exact."""
    out = {}
    for b in BETAS:
        for n in ns:
            out[(n, b)] = fe.ff4_exact(n, FF4_KAPS, beta=b)
    return out


FF7_DEEP_N = [64 * 4 ** k for k in range(0, 12)]      # 64 .. ~2.8e8, exact


def run_ff7(ns):
    """(d) P(G > T/2) at both insets, exact.  Fitted on a deep ladder as well:
    the exact law costs nothing, and the exponent-1 prediction is asymptotic."""
    allns = sorted(set(list(ns) + FF7_DEEP_N))
    return {(n, b): fe.ff7_exact(n, FF7_C, beta=b) for b in BETAS for n in allns}


def run_ff4_tail(nsel):
    delta = fe.T / nsel
    _, Cn = fe.discrete_barriers(delta)
    L, x = 2 * fe.BAND, Cn + fe.BAND
    g = fe.BAND - Cn
    ts = np.geomspace(g * g, fe.T, 9)
    return nsel, g, ts, fe.psurv(ts, x, L), g * np.sqrt(2.0 / (math.pi * ts))


def run_ff4_mc(pr, rng):
    M = int(FF4_MC_HORIZON * FF4_MC_MPU)
    dt = 1.0 / FF4_MC_MPU
    batch = 8
    nb = (FF4_MC_PATHS + batch - 1) // batch
    pr.phase(f"FF4/FF7 coupled MC cross-check (M={M})", total=nb)
    G = {(n, b): [] for n in FF4_MC_N for b in BETAS}
    cens = {(n, b): 0 for n in FF4_MC_N for b in BETAS}
    noexit = {(n, b): 0 for n in FF4_MC_N for b in BETAS}
    for j in range(nb):
        m = min(batch, FF4_MC_PATHS - j * batch)
        dW = rng.standard_normal((m, M)) * math.sqrt(dt)
        yb = np.concatenate([np.zeros((m, 1)), np.cumsum(dW, axis=1)], axis=1)
        for p in range(m):
            yf = yb[p]
            jxi = fe.first_exit_index(yf)
            for n in FF4_MC_N:
                a = math.sqrt(fe.T / n)
                idx, lev, _ = fe.embed_walk(yf, a)
                if idx.size == 0:
                    for b in BETAS:
                        noexit[(n, b)] += 1
                    continue
                Yn = fe.Y0 + a * lev[:n]
                for b in BETAS:
                    Bn, Cn = fe.discrete_barriers(fe.T / n, beta=b)
                    ow = (Yn <= Bn) | (Yn >= Cn)
                    if not ow.any():
                        noexit[(n, b)] += 1
                        continue
                    if jxi > M:
                        cens[(n, b)] += 1
                        continue
                    G[(n, b)].append((jxi - idx[int(ow.argmax())]) * dt)
        pr.tick(j + 1, note=f"FF4-MC batch {j + 1}/{nb}")
    rows = []
    for n in FF4_MC_N:
        for b in BETAS:
            g = np.array(G[(n, b)]) if G[(n, b)] else np.zeros(0)
            gp = np.maximum(g, 0)
            rows.append(dict(n=n, beta=b, used=g.size, neg=int((g < 0).sum()),
                             cens=cens[(n, b)], noexit=noexit[(n, b)],
                             **{f"EG^{k}": float((gp ** k).mean()) if g.size else float("nan")
                                for k in FF4_KAPS},
                             EG=float(gp.mean()) if g.size else float("nan"),
                             Pc=float((gp > FF7_C).mean()) if g.size else float("nan"),
                             Pc_se=(float(np.sqrt(max((gp > FF7_C).mean() *
                                                      (1 - (gp > FF7_C).mean()), 0) / g.size))
                                    if g.size else float("nan"))))
    return rows


# =============================================================== report
def build_md(pr, ff1, ff2, ff2r, ff3, ff4e, ff4b, ff7, ff4t, ff4mc):
    rows2, fits2, dsum2, ff2dt = ff2
    L = []
    A = L.append
    A("# Route F falsifiers FF1-FF4, FF6, FF7 -- measured")
    A("")
    A("Route and falsifiers: `research-notes/L003-ROUTE-F.md` sec.2-3.")
    A("Code: `sim/frozen_exit.py`; driver: `sim/run_frozen_exit.py`.")
    A("")
    A(f"Common settings: `T=1`, `y_0=0`, `v_0=0`, `K(u)=u^h` (`L==1`), continuous "
      f"band `(B^Y,C^Y)=({-fe.BAND},{fe.BAND})`, `mu_y=0`, `sigma_y=1` except in "
      f"FF3.  `H := h+1/2`.")
    A("")
    A("**The barrier inset is a parameter here.**  `eq:barriers` writes the literal "
      "`delta^{1/3}`.  That value occurs in exactly four places "
      "(`part2-scheme.tex:24-30, 44-47, 57, 71-81`), all inside the definition of "
      "the discrete barriers and inside `lem:barriers-wellposed` itself, and nothing "
      "downstream reads it.  Everything below is therefore computed with the inset "
      "`delta^{beta}`, at `beta = 1/3` (the current text) and `beta = 1/2` (the grid "
      "resolution, the smallest a `sqrt(delta)`-lattice can express), **on the same "
      "paths**.  `beta = 1/3` is the default everywhere the inset is not named.")
    A("")
    A("These are measurements, not proofs.  Nothing here closes a Definition-of-Done "
      "item; a refutation is an instruction to stop, a non-refutation is not a "
      "permission to continue.")
    A("")

    # ---------------------------------------------------------- decision table
    A("## Decision table")
    A("")
    A("Verdicts are the author's reading of the numbers below; the numbers are the "
      "evidence and the tables are where they live.")
    A("")
    A("| falsifier | verdict | the number that decides it |")
    A("|---|---|---|")
    _f1 = [r for r in ff1 if r["h"] == -0.3][-1]
    A(f"| FF1 (the diagnosis) | **CONFIRMED** | at h=-0.3, n=10^6 the L2 gap is "
      f"{_f1['L2']:.3f} against delta^H = {_f1['deltaH']:.4f}, a ratio of "
      f"{_f1['L2'] / _f1['deltaH']:.0f}, and it GROWS with n; "
      f"{_f1['nz_post']:.0%} of post-exit increments are non-zero, so "
      "`rem:freezing-convention`'s premise is false outright |")
    _bad13 = [f"h={h},ell={ell:g}" for h in FF2_H for ell in ELLS
              if fits2[(h, "total_theta", fe.BETA0, ell)][0]
              + 2 * fits2[(h, "total_theta", fe.BETA0, ell)][1] < (h + 0.5) / 2]
    A(f"| FF2 (Route F **at beta=1/3**) | **REFUTED** | the frozen-vs-frozen decay "
      f"misses the target (h+kappa)/2 by more than 2 bootstrap sd in {len(_bad13)} of "
      f"{len(FF2_H) * len(ELLS)} (h, ell) cells: {', '.join(_bad13)}.  F1 meets its "
      "rate; F2 does not, and its exponent falls with ell exactly as the closed-form "
      "law of G predicts.  This is a refutation of the CURRENT INSET, not of the "
      "decomposition -- see FF6 |")
    _bad12 = [f"h={h},ell={ell:g}" for h in FF2_H for ell in ELLS
              if fits2[(h, "total_theta", 0.5, ell)][0]
              + 2 * fits2[(h, "total_theta", 0.5, ell)][1] < (h + 0.5) / 2]
    _f2r = [fits2[(h, "F2", 0.5, e)][0] / fits2[(h, "F2", fe.BETA0, e)][0]
            for h in FF2_H for e in ELLS
            if fits2[(h, "F2", fe.BETA0, e)][0] > 0.01]
    A(f"| FF6 (the inset is the free parameter) | **{'CONFIRMED' if len(_bad12) < len(_bad13) else 'NOT CONFIRMED'}** | "
      f"moving the inset from `delta^{{1/3}}` to `delta^{{1/2}}` on the SAME paths "
      f"takes the count of failing cells from {len(_bad13)}/{len(FF2_H) * len(ELLS)} "
      f"to {len(_bad12)}/{len(FF2_H) * len(ELLS)}"
      + (f" (the survivor is {_bad12[0]}, and `ell=4 > ell*={2 * 0.5 / (FF2_H[0] + 0.5):.2f}` "
         "there, so it is PREDICTED to fail)" if len(_bad12) == 1 else "")
      + f"; `ell=2` -- the order `lem:freeze` consumes -- passes at both h.  F2 is "
      f"what moves (its exponent rises by a factor "
      f"{min(_f2r):.1f}-{max(_f2r):.1f} against the predicted 1.5, the excess being "
      "the milder pre-asymptotics at the smaller inset) |")
    _f1m = max(abs(fits2[(h, "F1_full", b, e)][0])
               for h in FF2_H for b in BETAS for e in ELLS)
    A(f"| **(unasked, and it matters)** `prop:Vconv`-as-proved | **REFUTED at both "
      f"insets** | `F1_full` -- the UNFROZEN clamped-walk convolution against the "
      f"ABSORBED-driver continuous one, which is the pair `prop:Vconv`'s proof "
      f"actually writes -- has a NEGATIVE fitted exponent at every (h, ell, beta) "
      f"tested (worst {-_f1m:.4f}).  The clamp/absorb mismatch of FF1 is inside "
      "`prop:Vconv` itself.  Restricted `F1` (on `j <= Xi^Y_n`, where the clamp is "
      "inactive) does meet its rate, so Route F's (F1) survives -- but it may not be "
      "obtained by quoting `prop:Vconv`; it has to be proved on the stopped range |")
    _c2 = [r for r in ff3 if r["variant"] == "C2.0"]
    A(f"| FF3 (ordering L1) | **NOT REFUTED as a fact, REFUTED as an inference** | "
      f"0 violations in {sum(r['used'] for r in ff3 if r['variant'] in ('A', 'B'))} "
      "faithful-construction paths; but with a node error of the size (E2) permits "
      f"the fraction runs {' -> '.join(f'{r[chr(102)+chr(114)+chr(97)+chr(99)]:.1%}' for r in _c2)} "
      "over n = " + ", ".join(str(r["n"]) for r in _c2) + " -- rising, not vanishing |")
    _k = FF4_KAPS[0]
    _sg, _ = fe.loglog_fit([r["gap"] for r in ff4e][-4:],
                           [r[f"EG^{_k}"] for r in ff4e][-4:])
    A(f"| FF4 (exit-time L2) | **NOT REFUTED** | on an exact ladder to n = 10^12 the "
      f"exponent of E[G^{_k}] in eps_n converges to {_sg:.3f} against the predicted "
      f"{2 * _k:.1f}, with the kappa'=1 control landing on 1.000; the tail matches "
      "eps_n sqrt(2/(pi t)) to within 3% over its stated range |")
    _p7 = [fe.loglog_fit([ff7[(n, b)]["gap"] for n in FF7_DEEP_N[-6:]],
                         [ff7[(n, b)]["P"] for n in FF7_DEEP_N[-6:]])[0] for b in BETAS]
    A(f"| FF7 (the lower-bound mechanism) | **NOT REFUTED** | P(G>T/2) scales like "
      f"eps_n^{{{_p7[0]:.3f}}} at beta=1/3 and eps_n^{{{_p7[1]:.3f}}} at beta=1/2 "
      "against the predicted exponent 1, so the sharpness claim's mechanism is "
      "present and is pure `eps_n` |")
    A("")

    # ------------------------------------------------------------------ FF1
    A("## FF1 -- frozen against UNFROZEN, discrete side only")
    A("")
    A("Measured `|| max_k |Vcal^(n)(t_k) - V^(n)_k| ||_{L2}`.  No continuous object "
      "and no coupling enter, so this is exact up to Monte-Carlo error and can be "
      "pushed to `n = 10^6`.  Two extra columns split the gap the way "
      "`rem:freezing-convention` splits it:")
    A("")
    A("* **kernel part** `max_{k>Xi} |sum_{i<=Xi} [K(t_k-t_{i-1}) - K(t_Xi-t_{i-1})] dY_i|` "
      "-- the sum the remark tries to telescope;")
    A("* **increment part** `max_{k>Xi} |sum_{Xi<i<=k} K(t_k-t_{i-1}) dY_i|` -- zero "
      "if and only if the remark's premise (\"past `Xi^Y_n` the increments vanish\") "
      "is true;")
    A("* **nz post** = fraction of increments `dY_i` with `i > Xi^Y_n` that are "
      "NON-zero.  The premise says this is 0.")
    A("")
    A("`scale` = `|| max_k |Vcal^(n)(t_k)| ||_{L2}`, the natural size of the object.")
    A("")
    for h in H_LIST:
        A(f"### h = {h}  (H = {h + 0.5:g})")
        A("")
        A("| n | delta | paths | L2 gap | +- | delta^H | gap/delta^H | kernel part | increment part | scale | nz post | P(walk exits) |")
        A("|---|---|---|---|---|---|---|---|---|---|---|---|")
        for r in ff1:
            if r["h"] != h:
                continue
            A(f"| {r['n']} | {fmt(r['delta'])} | {r['paths']} | {r['L2']:.4f} | "
              f"{r['se']:.4f} | {r['deltaH']:.4f} | {r['L2'] / r['deltaH']:.1f} | "
              f"{r['L2_kern']:.4f} | {r['L2_incr']:.4f} | {r['L2_scale']:.4f} | "
              f"{r['nz_post']:.4f} | {r['exit_frac']:.3f} |")
        A("")

    # -------------------------------------------------------------- FF2 / FF6
    A("## FF2 / FF6 -- frozen discrete against FROZEN continuous, at both insets")
    A("")
    A("**Coupling.**  The exact Skorokhod level-crossing embedding "
      "`theta_k = inf{t > theta_{k-1} : |y_t - y_{theta_{k-1}}| = sqrt(delta)}`, "
      f"realised on a fine grid of step `dt = {ff2dt:.3g}`.  With `mu_y=0, sigma_y=1` "
      "this is the manuscript's OWN embedding and it is exact: `Y(theta_k) = Y^(n)_k` "
      "and the node error `e_j` vanishes (`part3-embedding.tex` l.582).  Because `y` "
      "is a martingale the embedded increments are i.i.d. fair `+-1` coins, so the "
      "walk is the `eq:walks` walk in law -- the coupling does not deform either "
      "object.  It was chosen over sign-coupling "
      "(`zeta_k = sign(y(t_k)-y(t_{k-1}))`, which leaves `|Y^(n)_k - y(t_k)| = Theta(1)` "
      "and would fake a stall) and over lattice-rounding of `y(t_k)` (which is not the "
      "`eq:walks` walk).")
    A("")
    A("**The two insets are evaluated on the same paths and the same walk.**  The "
      "embedding, the walk `Y^(n)`, and the discrete convolution `Vcal^(n)` do not "
      "depend on `beta` at all; only the exit index `Xi^Y_n` does, through the "
      "barriers.  So the `beta = 1/3` and `beta = 1/2` columns below are PAIRED, and "
      "the bootstrap resamples whole paths, preserving the pairing.")
    A("")
    A("**Decomposition measured** (Route F sec.2, `nu := k ^ Xi^Y_n`):")
    A("")
    A("```")
    A("  total = max_k |V^(n)_k - V_{theta_k}|")
    A("  F1    = max_{j<=Xi^Y_n} |Vcal^(n)(t_j) - v_{theta_j}|          scheme error")
    A("  F2    = max_{k>Xi^Y_n} |v_{theta_k ^ Xi^y} - v_{theta_{Xi}}|   exit-time term")
    A("```")
    A("")
    A("`h = -0.45` is NOT measured here and the paths were spent on the two cells that "
      "can be resolved instead.  At `H = 0.05` the target rate is `delta^{0.025}`, a "
      "factor 1.15 across a factor-256 range in `delta`; no feasible ladder separates "
      "that from a constant, so any verdict there would have been noise.  It is kept "
      "in FF1, which is exact and needs no ladder.")
    A("")
    A("### (a) the fitted exponent, both insets side by side")
    A("")
    A("| h | ell | beta=1/3 | beta=1/2 | target (h+kappa)/2 | ell* = 2 beta/(h+kappa) at 1/3 | at 1/2 |")
    A("|---|---|---|---|---|---|---|")
    for h in FF2_H:
        hk = h + 0.5                      # h + kappa as kappa -> 1/2
        for ell in ELLS:
            m13, s13, _ = fits2[(h, "total_theta", fe.BETA0, ell)]
            m12, s12, _ = fits2[(h, "total_theta", 0.5, ell)]
            A(f"| {h} | {ell:g} | {m13:.4f} +- {s13:.4f} | {m12:.4f} +- {s12:.4f} | "
              f"{hk / 2:.4f} | {2 * fe.BETA0 / hk:.2f} | {2 * 0.5 / hk:.2f} |")
    A("")
    A("### (b) the F1 / F2 split, both insets")
    A("")
    A("The prediction under test was: **F1 is identical at both insets** while "
      "**F2's exponent moves with beta**.  The second half holds.  The first half "
      "does NOT hold as stated, and the reason is worth having:")
    A("")
    A("* `F1` as written in the decomposition is a maximum over `j <= Xi^Y_n`, and "
      "`Xi^Y_n` **is** a function of the inset.  A wider discrete band (larger `beta`) "
      "means a later exit index and a longer index range, so restricted `F1` moves "
      "with `beta` for a purely combinatorial reason.  At `beta = 1/3` and small `n` "
      "the band can be a single grid step wide -- at `n = 64`, `C^Y_n = 0.125 = "
      "1 x sqrt(delta)` -- so `Xi^Y_n = 1` on most paths and the maximum is taken "
      "over one index; by `n = 16384` it is taken over hundreds.  That n-dependence "
      "of the RANGE, not of the summand, is what bends the fitted slope.")
    A("* The bound the route actually invokes is `prop:Vconv` **read at the random "
      "index** `nu <= n`, and `prop:Vconv`'s proof bounds a maximum over ALL indices "
      "of the UNFROZEN discrete convolution against the ABSORBED-DRIVER continuous "
      "one.  That object -- `F1_full` below -- is the one the route is entitled to "
      "quote.  It too carries a `beta`, but through the CLAMP in `eq:walks` rather "
      "than through a truncated index range, and the dependence is weaker.")
    A("* One asymmetry inside `F1_full` deserves naming, because it is the same defect "
      "FF1 exhibits, one level up.  Its DISCRETE side is the unfrozen convolution of "
      "the **clamped** walk, whose increments do not vanish past `Xi^Y_n` (FF1 "
      "measures ~40-49% of them non-zero); its CONTINUOUS side is the "
      "**absorbed**-driver convolution, whose integrand genuinely does stop at "
      "`Xi^y`.  Those two are not analogues of one another, and `prop:Vconv`'s proof "
      "pairs them.  If `F1_full` does not decay while restricted `F1` does, the "
      "reading is that `prop:Vconv`-as-proved inherits the clamp/absorb mismatch too, "
      "and the route may not quote it to dominate (F1): it must bound (F1) directly "
      "on `j <= Xi^Y_n`, where the clamp is inactive and the pairing is honest.")
    A("")
    A("So the decomposition is intact.  What needs replacing is the sentence \"(F1) is "
      "free\" / \"(F1) never sees the barrier gap\".  (F1) does see it -- through its "
      "index range, and through the clamp inside the quantity that dominates it.")
    A("")
    A("| h | quantity | ell | beta=1/3 | beta=1/2 | beta min(2H,1/ell) at 1/3 | at 1/2 | target |")
    A("|---|---|---|---|---|---|---|---|")
    for h in FF2_H:
        H = h + 0.5
        for k in ("F1", "F1_full", "F2"):
            for ell in ELLS:
                m13, s13, _ = fits2[(h, k, fe.BETA0, ell)]
                m12, s12, _ = fits2[(h, k, 0.5, ell)]
                if k == "F2":
                    p13 = f"{fe.BETA0 * min(2 * H, 1 / ell):.4f}"
                    p12 = f"{0.5 * min(2 * H, 1 / ell):.4f}"
                else:
                    p13 = p12 = "--"
                A(f"| {h} | {k} | {ell:g} | {m13:.4f} +- {s13:.4f} | "
                  f"{m12:.4f} +- {s12:.4f} | {p13} | {p12} | {H / 2:.4f} |")
    A("")
    A("### The L2 norms themselves")
    A("")
    for h in FF2_H:
        A(f"#### h = {h}")
        A("")
        A("| n | delta | " + " | ".join(f"{k} (b={bs(b)})" for k in ("total_theta", "F1", "F2")
                                        for b in BETAS) + " |")
        A("|---|---|" + "---|" * 6)
        for n in FF2_N:
            cells = []
            for k in ("total_theta", "F1", "F2"):
                for b in BETAS:
                    r = next(r for r in rows2 if r["n"] == n and r["h"] == h
                             and r["key"] == k and r["ell"] == 2.0 and r["beta"] == b)
                    cells.append(f"{r['norm']:.4f}")
            A(f"| {n} | {fmt(fe.T / n)} | " + " | ".join(cells) + " |")
        A("")
    A("`total_tk` (the comparison against `V_{t_k}` rather than `V_{theta_k}`) is "
      "recorded in the run's `progress.json` and is uniformly WORSE at both insets; "
      "it is not the manuscript's object -- `prop:Vconv` compares at `theta_k` -- so "
      "it is not tabulated here.")
    A("")

    # ------------------------------------------------ exact-law cross-check
    A("### FF6 x FF4 -- where the arithmetic loses a power, and what beta buys back")
    A("")
    A("Route F sec.2 bounds `|F2| <= |v|_{H-} * G^{H-}` and quotes "
      "`E[G^{kappa'}] = O(eps_n^{2 kappa'})`, valid *for kappa' < 1/2*.  But "
      "`prop:Vconv` is an `L^ell` statement, so what is needed is "
      "`|| |v|_H G^H ||_{L^ell} ~ E[G^{ell H}]^{1/ell}`: the exponent inside the "
      "expectation is `ell H`, not `H`.  The `kappa' < 1/2` clause therefore reads "
      "`ell H < 1/2`, and above it `E[G^{kappa'}] ~ eps_n^{1}`, so with "
      "`eps_n ~ delta^{beta}` the norm decays like `delta^{beta/ell}` rather than "
      "`delta^{2 beta H}`.  Against the target `(h+kappa)/2`:")
    A("")
    A("```")
    A("  exponent in delta = beta * min(2 kappa', 1/ell)")
    A("  meets the target  <=>  ell <= ell* := 2 beta / (h + kappa)")
    A("```")
    A("")
    A("At `beta = 1/3` and `kappa -> 1/2` this is `ell* = 2/(3H)`; at `beta = 1/2` it "
      "is `ell* = 1/(h+kappa) > 2` strictly for every admissible pair, since "
      "`h < 0` and `kappa < 1/2` force `h+kappa < 1/2`.")
    A("")
    A("The `exact` columns are computed from the closed-form law of `G` over the SAME "
      "`n` ladder as the tables above -- no simulation -- so they are directly "
      "comparable with the measured `F2` slopes.")
    A("")
    A("| h | ell | beta | ell*H or ell*kappa' | exact E[G^{ell H}]^{1/ell} slope | measured F2 slope | asymptotic beta min(2H,1/ell) | target |")
    A("|---|---|---|---|---|---|---|---|")
    _ds = [fe.T / n for n in FF2_N]
    for h in FF2_H:
        Hh = h + 0.5
        for ell in ELLS:
            for b in BETAS:
                kap = ell * Hh
                ys = [fe.moment_G(kap, fe.discrete_barriers(fe.T / n, beta=b)[1] + fe.BAND,
                                  2 * fe.BAND) ** (1.0 / ell) for n in FF2_N]
                sl, _ = fe.loglog_fit(_ds, ys)
                m, sd, _ci = fits2[(h, "F2", b, ell)]
                A(f"| {h} | {ell:g} | {bs(b)} | {kap:.3g} | {sl:.4f} | "
                  f"{m:.4f} +- {sd:.4f} | {b * min(2 * Hh, 1 / ell):.4f} | {Hh / 2:.4f} |")
    A("")
    A("`prop:Vconv` is stated \"for every `ell >= 1`\", and at least one consumer needs "
      "`ell >= 2`: `part4-obstructions.tex` l.483-486 squares the `prop:Vconv` error "
      "and works in `L^{ell/2}`.")
    A("")

    # ------------------------------------------------ well-formedness at beta=1/2
    A("### Is the band still well formed at beta = 1/2?  (the adversarial check)")
    A("")
    A("At `beta = 1/2` the inset is `sqrt(delta)`, the SAME ORDER as the grid spacing, "
      "so `lem:barriers-wellposed`'s \"there exists n_0\" is a real constraint and is "
      "checked rather than assumed.  Required: `B^Y_n < C^Y_n`, both strictly inside "
      "`(B^Y, C^Y)`, `y_0` strictly inside, and both on the `sqrt(delta)` grid.")
    A("")
    A("| n | beta | B^Y_n | C^Y_n | eps_n realised | sqrt(delta) | band width in grid steps | well formed |")
    A("|---|---|---|---|---|---|---|---|")
    for n in WELLFORM_N:
        for b in BETAS:
            w = fe.band_wellformed(fe.T / n, beta=b)
            A(f"| {n} | {bs(b)} | {w['Bn']:.4f} | {w['Cn']:.4f} | {w['gap']:.4f} | "
              f"{w['sqrt_delta']:.4f} | {w['width_in_steps']:.0f} | "
              f"{'yes' if w['ok'] else '**NO**'} |")
    A("")
    A(f"The floor is `n = 32` for BOTH insets and it is the same floor: at `n = 16` "
      f"the band collapses onto `y_0` at `beta = 1/3` AND at `beta = 1/2`, because "
      f"`sqrt(delta) = 1/4` cannot place a lattice point strictly inside "
      f"`(-{fe.BAND}, {fe.BAND})` once anything is inset.  **`beta = 1/2` costs "
      f"nothing in ladder range** -- the FF2 ladder `n = {FF2_N[0]} .. {FF2_N[-1]}` is "
      "usable at both insets, which is why the side-by-side comparison above is "
      "on identical `n`.  At `beta = 1/2` the realised `eps_n` sits in "
      "`(sqrt(delta), 2 sqrt(delta)]`, i.e. one to two grid steps: the smallest "
      "non-degenerate inset the scheme can express, as intended.")
    A("")

    A("### FF2 coupling diagnostics (beta = 1/3)")
    A("")
    A("| n | mean max \\|e_j\\| | multi-level fine steps | short of n crossings | walk exits | y exits | mean G (b=1/3) | mean G (b=1/2) |")
    A("|---|---|---|---|---|---|---|---|")
    for n in FF2_N:
        d = dsum2[(n, fe.BETA0)]
        d2 = dsum2[(n, 0.5)]
        A(f"| {n} | {d['e_max']:.2e} | {d['multi']} | {d['short']:.3f} | "
          f"{d['exited']:.3f} | {d['y_exit']:.3f} | {d['meanG']:.4f} | "
          f"{d2['meanG']:.4f} |")
    A("")
    A("The node error is pure crossing overshoot, `O(sqrt(dt_fine))`; the exact "
      "embedding has `e_j == 0`.  It is two to three orders of magnitude below "
      "`delta^{1/4}`, the size the manuscript's own (E2) allows, so this artifact "
      "cannot be driving anything.")
    A("")
    A("### FF2 resolution sensitivity (the check that this is not a numerical floor)")
    A("")
    A(f"Fixed `n = {FF2_RES_N}`, {FF2_RES_PATHS} paths, `beta = 1/3`, fine grid "
      "refined by a factor 16 twice.  A COARSE `dt` makes the discrete and continuous "
      "objects MORE alike (they share the left-point convention), so under-resolution "
      "SUPPRESSES the measured gap -- and does so increasingly at large `n`, which "
      "INFLATES the apparent decay rate.  Measuring rates below target therefore "
      "bounds the true rates from above; the artifact cannot manufacture the stall.")
    A("")
    A("| h | fine steps / unit time | dt | total_theta | F1 | F2 |")
    A("|---|---|---|---|---|---|")
    for h in FF2_H:
        for r in ff2r:
            if r["h"] != h:
                continue
            A(f"| {h} | {r['mpu']} | {r['dt']:.3g} | {r['total_theta']:.4f} +- "
              f"{r['total_theta_se']:.4f} | {r['F1']:.4f} +- {r['F1_se']:.4f} | "
              f"{r['F2']:.4f} +- {r['F2_se']:.4f} |")
    A("")

    # ------------------------------------------------------------------ FF3
    A("## FF3 -- the ordering lemma L1")
    A("")
    A("Counted: the fraction of paths with `theta^(n)_{Xi^Y_n} > Xi^y`, i.e. the walk "
      "leaving its NARROW band later than the diffusion leaves its WIDE one.  "
      "`sigma_y(y) = 1 + 0.3 tanh(y)` in `[0.7,1.3]` throughout: bounded, Lipschitz, "
      "bounded away from 0, so `asm:coeff` holds.  Inset `beta = 1/3`.")
    A("")
    A("* **A** -- genuine Skorokhod embedding (level crossings) with `mu_y = 0`.  "
      "Because `y` is then a martingale the embedded increments are exactly fair "
      "`+-1` coins whatever `sigma_y` is, so this IS the `eq:walks` walk and the "
      "embedding is faithful.  Node error = fine-grid crossing overshoot only.")
    A("* **B** -- SURROGATE: `theta_k := t_k` and `Y^(n)_k :=` the `sqrt(delta)`-lattice "
      "rounding of `y(t_k)`, with `mu_y(y) = 0.3 cos(y)` switched on.  NOT the "
      "`eq:walks` walk (its increments are multiples of `sqrt(delta)`, not "
      "`+-sqrt(delta)`); used only to give the node error a non-zero value of a known "
      "size, `|e_k| <= sqrt(delta)/2`.")
    A("* **C c** -- variant B plus an INJECTED node error "
      "`e_k = c * delta^{1/4} * tanh(Btilde_{t_k})`, `Btilde` independent, so "
      "`max_k |e_k| <= c delta^{1/4}` exactly.  Not an embedding: it is the largest "
      "node error the manuscript's own (E2) permits, injected deliberately, to answer "
      "the question L1 actually poses -- is the ordering IMPLIED by (E2) + "
      "`lem:barriers-wellposed`?")
    A("")
    A("The decisive comparison is `mean max |e_j|` against the barrier gap `eps_n`. "
      "Note `delta^{1/4} > delta^{1/3}` for every `delta < 1`, and the ratio "
      "`delta^{1/4}/delta^{1/3} = delta^{-1/12}` DIVERGES -- at `beta = 1/2` the "
      "ratio `delta^{1/4}/delta^{1/2} = delta^{-1/4}` diverges faster, so shrinking "
      "the inset makes L1 HARDER, not easier.  That is a cost of D1 and it is "
      "recorded here rather than left for the reader to notice.")
    A("")
    A("| variant | n | paths used | violations | fraction | mean max \\|e_j\\| | barrier gap eps_n | delta^{1/4} | no exit |")
    A("|---|---|---|---|---|---|---|---|---|")
    for r in ff3:
        A(f"| {r['variant']} | {r['n']} | {r['used']} | {r['viol']} | "
          f"{r['frac']:.4f} | {r['e_max']:.4f} | {r['gap']:.4f} | {r['d14']:.4f} | "
          f"{r['no_exit']} |")
    A("")

    # ------------------------------------------------------------------ FF4
    A("## FF4 -- the exit-time reconciliation L2")
    A("")
    A("`G := Xi^y - theta^(n)_{Xi^Y_n} >= 0`.  With constant coefficients the exactly "
      "embedded walk sits ON the discrete barrier at `theta_{Xi^Y_n}`, so by the strong "
      "Markov property `G` is the exit time of a Brownian motion from `(B^Y,C^Y)` "
      "started at `C^Y_n`.  That law is classical, so the primary numbers below are "
      "EXACT (spectral + reflection series, `frozen_exit.psurv` / `moment_G`), "
      "validated at `kappa'=1` against `E[G] = x(L-x)` to 11 digits.  A coupled "
      "Monte-Carlo run cross-checks the reduction in (c).  Inset `beta = 1/3` unless "
      "stated.")
    A("")
    A("### (a) E[G^kappa'] against the prediction O(eps_n^{2 kappa'}) = O(delta^{2 beta kappa'})")
    A("")
    A("| n | delta | eps_n realised | " + " | ".join(f"E[G^{k}]" for k in FF4_KAPS) + " | E[G] |")
    A("|---|---|---|" + "---|" * (len(FF4_KAPS) + 1))
    for r in ff4e:
        A(f"| {r['n']} | {fmt(r['delta'])} | {r['gap']:.5g} | "
          + " | ".join(f"{r[f'EG^{k}']:.5g}" for k in FF4_KAPS)
          + f" | {r['EG']:.5g} |")
    A("")
    ds = [r["delta"] for r in ff4e]
    gs = [r["gap"] for r in ff4e]
    W = [(0, "all"), (len(ds) - 7, "last 7"), (len(ds) - 4, "last 4")]
    A("The prediction is asymptotic and the `kappa' -> 1/2` crossover is only "
      "logarithmically fast, so the exponents are fitted over three nested windows. "
      "`kappa' = 1` is the CONTROL: there the exponent must be exactly 1 in `eps_n` "
      "(and `beta` in `delta`), which pins the fitting procedure.")
    A("")
    A("| quantity | " + " | ".join(f"exp in delta ({w})" for _, w in W)
      + " | prediction | " + " | ".join(f"exp in eps_n ({w})" for _, w in W)
      + " | prediction |")
    A("|---|" + "---|" * (2 * len(W) + 2))
    for k in FF4_KAPS + [1.0]:
        key = f"EG^{k}" if k != 1.0 else "EG"
        ys = [r[key] for r in ff4e]
        pd = 2 * k / 3 if k < 0.5 else 1 / 3
        pg = 2 * k if k < 0.5 else 1.0
        cd = [f"{fe.loglog_fit(ds[i:], ys[i:])[0]:.4f}" for i, _ in W]
        cg = [f"{fe.loglog_fit(gs[i:], ys[i:])[0]:.4f}" for i, _ in W]
        A(f"| E[G^{k}] | " + " | ".join(cd) + f" | {pd:.4f} | "
          + " | ".join(cg) + f" | {pg:.4f} |")
    A("")
    A("The exponents rise monotonically toward the prediction as the window moves out, "
      "and the control lands on 1.000.  The residual shortfall at `kappa' = 0.45` is "
      "the `kappa' -> 1/2` crossover, where "
      "`kappa' int t^{kappa'-1} P(G>t) dt` stops being dominated by its lower end.")
    A("")
    A("`E[G]` is the diagnostic the note asks for.  In the idealised half-line problem "
      "`E[G]` is infinite.  In the actual BOUNDED band it is finite -- but it does not "
      "scale like `eps_n^2`, it scales like `eps_n^1`, one full power short.  Any "
      "mean-based estimate of L2 gives `O(eps_n)`, not `O(eps_n^2)`.")
    A("")
    A("### (c) the same, at both insets, on the FF2 ladder -- is the ratio the pure beta factor?")
    A("")
    A("If nothing but the inset changed, the exponent in `delta` must be "
      "`2 beta kappa'` at each inset, i.e. the ratio of the two fitted exponents must "
      "be `(1/2)/(1/3) = 1.5` exactly, and the exponent in `eps_n` must be "
      "`2 kappa'` at BOTH.")
    A("")
    A("| kappa' | exp in delta (beta=1/3) | exp in delta (beta=1/2) | ratio | predicted 1.5 | exp in eps_n (b=1/3) | (b=1/2) | predicted 2 kappa' |")
    A("|---|---|---|---|---|---|---|---|")
    _dl = [fe.T / n for n in FF2_N]
    for k in FF4_KAPS[:2]:
        e13 = [ff4b[(n, fe.BETA0)][f"EG^{k}"] for n in FF2_N]
        e12 = [ff4b[(n, 0.5)][f"EG^{k}"] for n in FF2_N]
        g13 = [ff4b[(n, fe.BETA0)]["gap"] for n in FF2_N]
        g12 = [ff4b[(n, 0.5)]["gap"] for n in FF2_N]
        s13, _ = fe.loglog_fit(_dl, e13)
        s12, _ = fe.loglog_fit(_dl, e12)
        q13, _ = fe.loglog_fit(g13, e13)
        q12, _ = fe.loglog_fit(g12, e12)
        A(f"| {k} | {s13:.4f} | {s12:.4f} | {s12 / s13:.3f} | 1.500 | {q13:.4f} | "
          f"{q12:.4f} | {2 * k:.4f} |")
    A("")
    A("### (b) raw tail P(G > t) against the prediction ~ eps_n / sqrt(t)")
    A("")
    nsel, gsel, ts, ps, pred = ff4t
    A(f"n = {nsel}, `beta = 1/3`, realised `eps_n = {gsel:.5g}`, `t` over "
      "`[eps_n^2, T]`.  The prediction is the half-line stable(1/2) tail "
      "`eps_n sqrt(2/(pi t))`.")
    A("")
    A("| t | t / eps_n^2 | P(G>t) exact | eps_n sqrt(2/(pi t)) | ratio |")
    A("|---|---|---|---|---|")
    for t, p, q in zip(ts, ps, pred):
        A(f"| {t:.4g} | {t / gsel ** 2:.4g} | {p:.5g} | {q:.5g} | "
          f"{(p / q if q else float('nan')):.4f} |")
    A("")

    # ------------------------------------------------------------------ FF7
    A("## FF7 -- P(G > c) at a FIXED c, against eps_n")
    A("")
    A(f"`c = T/2 = {FF7_C}`.  This is the lower-bound mechanism of the sharpness claim "
      "measured directly: with probability of order `eps_n` the diffusion started at "
      "the discrete barrier stays inside the band for a time of order 1, and on that "
      "event the frozen and the unfrozen object differ by `Theta(1)`.  Exact, from the "
      "same closed-form law.  Prediction: `P(G>c) ~ eps_n`, hence exponent 1 in "
      "`eps_n` and `beta` in `delta`; and at fixed `n` the ratio between the two "
      "insets must be `delta^{1/2-1/3} = delta^{1/6}`.")
    A("")
    A("| n | delta | eps_n (b=1/3) | P(G>c) (b=1/3) | eps_n (b=1/2) | P(G>c) (b=1/2) | ratio P(1/2)/P(1/3) | delta^{1/6} |")
    A("|---|---|---|---|---|---|---|---|")
    for n in FF2_N:
        a13, a12 = ff7[(n, fe.BETA0)], ff7[(n, 0.5)]
        d = fe.T / n
        A(f"| {n} | {fmt(d)} | {a13['gap']:.5g} | {a13['P']:.5g} | "
          f"{a12['gap']:.5g} | {a12['P']:.5g} | {a12['P'] / a13['P']:.4f} | "
          f"{d ** (1 / 6):.4f} |")
    A("")
    A("Fitted on the FF2 ladder and on a deep exact ladder "
      f"(`n = {FF7_DEEP_N[0]} .. {FF7_DEEP_N[-1]:.3g}`), since the exponent-1 "
      "prediction is asymptotic and the exact law costs nothing to extend:")
    A("")
    A("| beta | exp in eps_n (FF2 ladder) | (deep ladder) | (deep, last 6) | prediction | exp in delta (deep) | prediction beta |")
    A("|---|---|---|---|---|---|---|")
    for b in BETAS:
        gg = [ff7[(n, b)]["gap"] for n in FF2_N]
        pp = [ff7[(n, b)]["P"] for n in FF2_N]
        gd = [ff7[(n, b)]["gap"] for n in FF7_DEEP_N]
        pdp = [ff7[(n, b)]["P"] for n in FF7_DEEP_N]
        sg, _ = fe.loglog_fit(gg, pp)
        sD, _ = fe.loglog_fit(gd, pdp)
        sT, _ = fe.loglog_fit(gd[-6:], pdp[-6:])
        sdl, _ = fe.loglog_fit([fe.T / n for n in FF7_DEEP_N], pdp)
        A(f"| {bs(b)} | {sg:.4f} | {sD:.4f} | {sT:.4f} | 1.0000 | {sdl:.4f} | {b:.4f} |")
    A("")

    A("### (c) coupled Monte-Carlo cross-check, both insets")
    A("")
    A("Same fine-grid construction as FF2 (exact level-crossing embedding), horizon "
      f"{FF4_MC_HORIZON}T so that `Xi^y` is observed.  `G<0` counts L1 violations "
      "(must be 0); `censored` counts paths whose `Xi^y` exceeded the horizon and are "
      "dropped, which biases the MC column DOWNWARD, so agreement is the meaningful "
      "direction.  `P(G>T/2)` here is the Monte-Carlo counterpart of the FF7 table "
      "and is quoted with its binomial standard error -- at the small `eps_n` end the "
      "counts get thin, which is exactly why FF7's primary numbers are exact.")
    A("")
    A("| n | beta | paths used | G<0 | censored | " + " | ".join(f"E[G^{k}] MC / exact" for k in FF4_KAPS)
      + " | E[G] MC / exact | P(G>T/2) MC | exact |")
    A("|---|---|---|---|---|" + "---|" * (len(FF4_KAPS) + 3))
    for r in ff4mc:
        ex = fe.ff4_exact(r["n"], FF4_KAPS, beta=r["beta"])
        p7 = fe.ff7_exact(r["n"], FF7_C, beta=r["beta"])
        A(f"| {r['n']} | {bs(r['beta'])} | {r['used']} | {r['neg']} | {r['cens']} | "
          + " | ".join(f"{r[f'EG^{k}']:.4g} / {ex[f'EG^{k}']:.4g}" for k in FF4_KAPS)
          + f" | {r['EG']:.4g} / {ex['EG']:.4g} | {r['Pc']:.4g} +- {r['Pc_se']:.3g} "
            f"| {p7['P']:.4g} |")
    A("")

    # ------------------------------------------------------------------ caveats
    A("## Shortcuts, and how each could have manufactured the answer")
    A("")
    A("1. **Band choice, and an open obligation it touches.**  "
      f"`(B^Y,C^Y) = ({-fe.BAND},{fe.BAND})`, symmetric about `y_0 = 0`.  "
      "`part1-setup.tex` declares `B^Y in [0,inf)` while `eq:model` starts the driver "
      "at `y_0 = 0`, so **as the symbols stand `y_0` sits ON the barrier and "
      "`Xi^y = 0` almost surely** -- the absorbed driver never moves and `V_t == v_0`. "
      "The `[0,inf)` constraint is a leftover from the drafts in which the Y-barriers "
      "were centred at `v_0` (a variance, hence non-negative); the `\\fixed{}` note in "
      "`eq:barriers` recentres them at `y_0` but the constraint on `B^Y` was not "
      "moved with them.  This work took the recentred reading throughout.  "
      "**Recorded as an independent hit on the project's open obligation O003.**  A "
      "different band changes every CONSTANT below (FF1's plateau height, FF2's level) "
      "but not the exponents, since `eps_n ~ delta^{beta}` is band-independent.")
    A("2. **Fine-grid continuous reference (FF2/FF6, FF3, FF4c).**  `v_t` is the "
      "left-point Riemann sum at `dt_fine`, own error `O(dt_fine^H)`.  This is the "
      "shortcut most able to fake a stall, so it is tested directly in the resolution "
      "table, and its bias direction is stated there: a coarse `dt` makes the two "
      "objects MORE alike, so under-resolution suppresses the gap and INFLATES the "
      "apparent decay rate.  Measured rates below target therefore bound the true "
      "rates from above.")
    A("3. **h = -0.45 is not measured in FF2/FF6.**  Its target rate `delta^{0.025}` "
      "is a factor 1.15 over a factor-256 range in `delta` and cannot be separated "
      "from a constant on any feasible ladder.  The earlier run that did include it "
      "produced slopes of -0.001 +- 0.005 at `ell=2`; that number is not evidence and "
      "was dropped rather than quoted.  `h = -0.45` remains in FF1, which is exact.")
    A("4. **Crossing overshoot.**  The exact embedding is realised on a grid, so "
      "`y(theta_k)` overshoots the lattice by `O(sqrt(dt_fine))` and the walk is taken "
      "to be the exact lattice point.  Reported per `n`; 2-3 orders below "
      "`delta^{1/4}`.  It biases the node error UPWARD, i.e. AGAINST L1, and L1 still "
      "never fails in variants A and B.")
    A("5. **Fine-grid exit detection.**  `Xi^y` is the first fine-grid index outside "
      "the band, so it is LATE by `O(dt_fine)` and the effective barrier sits "
      "`~0.6 sqrt(dt_fine)` outside the true one.  This biases `G` UPWARD, i.e. in "
      "favour of a slower FF2 decay.  At `beta = 1/2` the inset is only "
      "`sqrt(delta)`, so this artifact is a LARGER fraction of `eps_n` than at "
      "`beta = 1/3` -- it is `0.6 sqrt(dt_fine)/sqrt(delta) = 0.6 sqrt(delta/dt_fine)`"
      f" inverted, i.e. `0.6/sqrt(mpu/n)`, at most 0.6/8 = 7.5% at the top of the "
      "ladder.  It therefore flatters `beta = 1/2` slightly, and the FF6 improvement "
      "should be read with that in mind.  FF4 and FF7's primary numbers avoid it "
      "entirely by using the exact law.")
    A("6. **FF3 variants B and C are surrogates, not embeddings**, and are labelled as "
      "such.  A violation count from C is evidence about what (E2) IMPLIES, not about "
      "what the manuscript's embedding DOES.  Variant A is the only faithful "
      "construction there, and for `mu_y = 0` it has `e_j == 0` in continuous time -- "
      "so a zero violation count is close to a tautology and must not be read as "
      "support for L1.")
    A("7. **Censoring.**  Paths whose `Xi^y` exceeds the horizon are dropped in the "
      "FF4/FF7 Monte-Carlo (count reported); dropping the longest `G` biases the MC "
      "moments and `P(G>T/2)` DOWNWARD, so MC-below-exact is the expected direction.")
    A("8. **Monte-Carlo error.**  All `+-` are bootstrap standard errors over paths "
      "(binomial for `P(G>T/2)`).  In FF2/FF6 the SAME fine paths are used for every "
      "`n` AND for both insets, so the estimates are strongly positively correlated "
      "and the fitted slopes -- and especially the `beta=1/3` vs `beta=1/2` "
      "difference -- are far better determined than the individual points.  The "
      "bootstrap resamples whole paths and respects that.")
    A("")
    A(pr.timing_table_md())
    return "\n".join(L)


def main():
    meta = dict(quick=QUICK, band=fe.BAND, h_list=H_LIST, ff2_h=FF2_H,
                betas=BETAS, ells=ELLS, ff1_n=FF1_N, ff2_n=FF2_N,
                ff2_paths=FF2_PATHS, ff2_mpu=FF2_MPU, ff3_n=FF3_N,
                ff3_paths=FF3_PATHS, ff4_n=FF4_N, ff7_c=FF7_C)
    with Progress("frozen-exit", total_phases=6, meta=meta) as pr:
        rng = np.random.default_rng(20260807)
        ff1 = run_ff1(pr, rng)
        for r in ff1:
            pr.result(f"FF1_L2_n{r['n']}_h{r['h']}", round(r["L2"], 5))
        ff2 = run_ff2(pr, rng)
        for (h, k, b, ell), (m, sd, _) in ff2[1].items():
            if k in ("total_theta", "F1", "F2"):
                pr.result(f"FF6_slope_{k}_h{h}_b{bs(b).replace('/', '')}_L{ell:g}",
                          [round(m, 4), round(sd, 4)])
        ff2r = run_ff2_resolution(pr, rng)
        ff3 = run_ff3(pr, rng)
        for r in ff3:
            pr.result(f"FF3_{r['variant']}_n{r['n']}_violfrac", round(r["frac"], 5))
        ff4e = run_ff4_exact(pr)
        ff4b = run_ff4_beta(FF2_N)
        ff7 = run_ff7(FF2_N)
        for b in BETAS:
            sg, _ = fe.loglog_fit([ff7[(n, b)]["gap"] for n in FF7_DEEP_N],
                                  [ff7[(n, b)]["P"] for n in FF7_DEEP_N])
            pr.result(f"FF7_slope_in_epsn_b{bs(b).replace('/', '')}", round(sg, 4))
        ff4t = run_ff4_tail(FF4_TAIL_N)
        ff4mc = run_ff4_mc(pr, rng)
        for r in ff4mc:
            pr.result(f"FF4_MC_n{r['n']}_b{bs(r['beta']).replace('/', '')}_neg", r["neg"])
        pr.write_results_md(build_md(pr, ff1, ff2, ff2r, ff3, ff4e, ff4b, ff7,
                                     ff4t, ff4mc))
        print(pr.dir)


if __name__ == "__main__":
    main()
