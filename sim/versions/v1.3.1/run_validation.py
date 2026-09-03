#!/usr/bin/env python3
"""
run_validation — the validation stack of paper §10.3, on the simplest example
that still exercises the whole pipeline.

EXAMPLE.  Fractional Brownian driver, constant price volatility:

    Y = W  (sigma_y = 1, mu_y = 0),      K(u) = u^h,  h = H - 1/2
    X = log S,  sigma_x = sigma const,   mu_x = -sigma^2/2
    payoff  g(t,x,v) = (Kstrike - e^x)^+     (American put)

Because sigma_x does not depend on V, the American put has an EXACT
one-dimensional reference (a CRR binomial with many steps).  The two-dimensional
tree must reproduce it: that is test (V8), and it exercises the joint transition
kernel, the recombination, the 2-D backward induction and the stopping rule.
The rough component V is still built and is tested separately by (V1)-(V2).

PHASES
  1  (V1,V2)  variance and covariance of the exact convolution V^(n) versus the
              one-step lattice Vcheck^(n) versus the truth.  Closed form, plus a
              Monte-Carlo cross-check.  Predicted: ratio Var[Vcheck]/Var[V]
              diverges with log-log slope 1-2H.
  2  (V3,V8)  2-D lattice American put versus the CRR reference, and the
              admissibility assertion on every transition probability.
  3  (V13)    self-consistent rate |Lambda^(2n) - Lambda^(n)| ~ delta^r for a
              genuinely rough case (eta > 0), which has no reference.

Progress is written continuously to runs/<name>/progress.json; ask at any time
with `bin/simstatus`.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from progress import Progress  # noqa: E402

# ----------------------------------------------------------------- parameters
T = 1.0
S0 = 100.0
KSTRIKE = 100.0
SIGMA = 0.30           # constant price volatility for phases 1-2
RHO = -0.70
HURSTS = (0.10, 0.20, 0.30, 0.40)
N_PHASE2 = (8, 16, 32, 64, 128, 256)
N_PHASE3 = (8, 16, 32, 64, 128)
ETA = 0.30             # rough Bergomi vol-of-vol for phase 3
XI0 = SIGMA ** 2       # flat forward variance


# =========================================================================== #
# phase 1 — variance / covariance of the rough component
# =========================================================================== #
def var_exact_convolution(n: int, H: float) -> float:
    """Var[V^(n)_T] = delta^{2H} sum_{j=1}^{n} j^{2h},  h = H - 1/2."""
    d = T / n
    h = H - 0.5
    j = np.arange(1, n + 1, dtype=float)
    return d ** (2 * H) * np.sum(j ** (2 * h))


def var_one_step(n: int, H: float) -> float:
    """Var[Vcheck^(n)_T] = n * delta^{2H}."""
    return n * (T / n) ** (2 * H)


def var_truth(H: float) -> float:
    """Var[V_T] = T^{2H} / (2H)."""
    return T ** (2 * H) / (2 * H)


def cov_exact_convolution(n: int, H: float, kj: int, kk: int) -> float:
    d = T / n
    h = H - 0.5
    i = np.arange(1, min(kj, kk) + 1, dtype=float)
    return d ** (2 * H) * np.sum((kj - i + 1) ** h * (kk - i + 1) ** h)


def cov_truth(H: float, s: float, t: float) -> float:
    """int_0^{s^t} (s-u)^h (t-u)^h du, computed numerically (h in (-1/2,0))."""
    h = H - 0.5
    m = min(s, t)
    u = (np.arange(200_001) / 200_000.0) * m
    f = np.where(u < m, (s - u) ** h * (t - u) ** h, 0.0)
    return float(np.trapezoid(f, u)) if hasattr(np, "trapezoid") else float(np.trapz(f, u))


def mc_var_check(n: int, H: float, paths: int, rng) -> tuple[float, float]:
    """Monte-Carlo Var of both discretisations, as an independent cross-check."""
    d = T / n
    h = H - 0.5
    z = rng.choice(np.array([-1.0, 1.0]), size=(paths, n))
    w = np.power(np.arange(n, 0, -1, dtype=float), h)      # (n-i+1)^h, i=1..n
    v_exact = (d ** H) * (z @ w)                           # delta^H sum (k-i+1)^h zeta_i
    v_step = (d ** H) * z.sum(axis=1)
    return float(v_exact.var()), float(v_step.var())


def phase1(pr: Progress, rng) -> dict:
    out: dict = {}
    pr.phase("V1/V2 variance+covariance of the rough component",
             total=len(HURSTS) * 2)
    ns = np.array([2 ** p for p in range(3, 13)])
    for H in HURSTS:
        ve = np.array([var_exact_convolution(int(n), H) for n in ns])
        vs = np.array([var_one_step(int(n), H) for n in ns])
        vt = var_truth(H)
        ratio_step = vs / vt
        ratio_exact = ve / vt
        slope = float(np.polyfit(np.log(ns), np.log(ratio_step), 1)[0])
        out[f"H={H:.2f}"] = {
            "slope_ratio_onestep_vs_n": round(slope, 4),
            "predicted_slope_1_minus_2H": round(1 - 2 * H, 4),
            "ratio_exact_at_n4096": round(float(ratio_exact[-1]), 4),
            "ratio_onestep_at_n4096": round(float(ratio_step[-1]), 2),
        }
        pr.tick(note=f"V1 H={H}")
        # covariance: is Cov a function of j^k only?
        n = 512
        c_far = cov_exact_convolution(n, H, n // 4, n)
        c_near = cov_exact_convolution(n, H, n // 4, n // 4 + 1)
        c_true_far = cov_truth(H, T / 4, T)
        c_true_near = cov_truth(H, T / 4, T / 4 + T / n)
        out[f"H={H:.2f}"].update({
            "cov_exact_far_over_true": round(c_far / c_true_far, 3),
            "cov_exact_near_over_true": round(c_near / c_true_near, 3),
            "cov_onestep_far_over_near": 1.0,   # by construction: fn of j^k only
            "cov_true_far_over_near": round(c_true_far / c_true_near, 3),
        })
        pr.tick(note=f"V2 H={H}")
        pr.partial.update({f"slope_H{H:.2f}": round(slope, 3)})
    # Monte-Carlo cross-check at one H
    mc_e, mc_s = mc_var_check(256, 0.10, 40_000, rng)
    out["mc_check_H0.10_n256"] = {
        "mc_var_exact": round(mc_e, 6),
        "closed_form_exact": round(var_exact_convolution(256, 0.10), 6),
        "mc_var_onestep": round(mc_s, 6),
        "closed_form_onestep": round(var_one_step(256, 0.10), 6),
    }
    return out


# =========================================================================== #
# phase 2/3 — the two-dimensional lattice
# =========================================================================== #
def crr_american_put(n: int, sigma: float) -> float:
    """CRR binomial American put — the one-dimensional reference."""
    d = T / n
    u = np.exp(sigma * np.sqrt(d))
    up, dn = u, 1.0 / u   # risk-neutral probability with r = 0
    p = (1.0 - dn) / (up - dn)
    j = np.arange(n + 1)
    s = S0 * up ** j * dn ** (n - j)
    v = np.maximum(KSTRIKE - s, 0.0)
    for k in range(n - 1, -1, -1):
        j = np.arange(k + 1)
        s = S0 * up ** j * dn ** (k - j)
        v = p * v[1:k + 2] + (1.0 - p) * v[0:k + 1]
        v = np.maximum(v, KSTRIKE - s)
    return float(v[0])


def zmax_feasible(H: float, eta: float, rho: float, support: str = "4pt") -> float:
    """Widest |zeta| for which the transition kernel stays a probability.

    With sigma_y = 1 (so beta = 0) and sbar_X = sup sigma_x, one has
    A = sigma_x / sbar_X and, expanding |a| + |c| <= 1/4 for the FOUR-point
    kernel (4.10),

        |a| + |c| <= 1/4   <=>   (1 - A^2) + 2|rho| A <= 1 + A^2   <=>   |rho| <= A.

    For the NINE-point kernel of Theorem 9.1 the condition is instead
    rho^2 <= inf A / sup A = inf A (Corollary 9.3), i.e. exactly twice the
    admissible width in zeta.  Since sigma_x varies by exp(eta sqrt(2H) Z) over
    |zeta| <= Z, we get

        4-point:  Z <=     log(1/|rho|) / (eta sqrt(2H)),
        9-point:  Z <= 2 * log(1/|rho|) / (eta sqrt(2H)).
    """
    if eta == 0.0 or rho == 0.0:
        return np.inf
    base = np.log(1.0 / abs(rho)) / (eta * np.sqrt(2.0 * H))
    return base if support == "4pt" else 2.0 * base


def zmax_admissible(n: int, H: float, eta: float, rho: float) -> float:
    """The finite-delta refinement of `zmax_feasible`.

    Keeping the drift term, |a| + |c| <= 1/4 reads  A^2 - |rho| A - eps >= 0 with
    eps = |mu_x| sqrt(delta) / sbar_X = O(sqrt(delta)), so the requirement is

        A_inf >= ( |rho| + sqrt(rho^2 + 4 eps) ) / 2,

    which is strictly stronger than |rho| <= A_inf and tightens as delta grows.
    Solved by three fixed-point iterations, eps depending on Z through sup v.
    """
    if eta == 0.0 or rho == 0.0:
        return np.inf
    c = eta * np.sqrt(2.0 * H)
    sqd = np.sqrt(T / n)
    Z = zmax_feasible(H, eta, rho, "4pt")
    for _ in range(3):
        sup_v = XI0 * np.exp(c * Z)
        eps = 0.5 * np.sqrt(sup_v) * sqd
        a_inf = 0.5 * (abs(rho) + np.sqrt(rho ** 2 + 4.0 * eps))
        Z = max(0.0, -np.log(min(a_inf, 1.0)) / c)
    return float(Z)


def lattice_american_put(n: int, H: float, eta: float, rho: float,
                         zmax: float = np.inf,
                         pr: Progress | None = None) -> tuple[float, int, float]:
    """Two-dimensional recombining lattice of paper §4-§5.

    Returns (value, n_violations, max_violation) where a violation is a
    transition probability outside [0,1] — test (V3).

    State at step k: (l, m, xi_prev, zeta_prev), l,m in {0..k}.
      X = log S0 + (2l-k)*hx        hx = sbar_x * sqrt(delta)
      Ycoord = (2m-k)*hy            hy = sbar_y * sqrt(delta) = sqrt(delta)
      zeta_coord = (2m-k)*delta^H   (the discrete GFO value)
      v = xi0*exp(eta*sqrt(2H)*zeta_coord - eta^2 t^{2H}/2)   (eta>0)
        = xi0                                                  (eta=0)
    """
    d = T / n
    h = H - 0.5
    sqd = np.sqrt(d)
    dH = d ** H                       # K(delta)*sqrt(delta)
    sbar_y = 1.0

    def variance(k: int, m_idx: np.ndarray) -> np.ndarray:
        if eta == 0.0:
            return np.full(m_idx.shape, XI0)
        zc = (2 * m_idx - k) * dH
        if np.isfinite(zmax):                       # absorbing barrier on the
            zc = np.clip(zc, -zmax, zmax)           # variance coordinate
        t = k * d
        return XI0 * np.exp(eta * np.sqrt(2 * H) * zc - 0.5 * eta ** 2 * t ** (2 * H))

    # sbar_x: normalisation of Remark 4.11 — sup of sigma_x over the tree
    if eta == 0.0:
        sbar_x = np.sqrt(XI0)
    else:
        zcap = min(n * dH, zmax) if np.isfinite(zmax) else n * dH
        sbar_x = np.sqrt(XI0 * np.exp(eta * np.sqrt(2 * H) * zcap))
    hx = sbar_x * sqd

    viol = 0
    maxviol = 0.0

    # terminal layer
    l = np.arange(n + 1)[:, None]
    m = np.arange(n + 1)[None, :]
    logS = np.log(S0) + (2 * l - n) * hx + 0.0 * m       # broadcast to (n+1,n+1)
    val = np.maximum(KSTRIKE - np.exp(logS), 0.0)
    val = np.ascontiguousarray(np.broadcast_to(val, (n + 1, n + 1)))
    val = np.repeat(np.repeat(val[:, :, None, None], 2, axis=2), 2, axis=3)

    for k in range(n - 1, -1, -1):
        li = np.arange(k + 1)[:, None]
        mi = np.arange(k + 1)[None, :]
        v_here = variance(k, mi)                       # (1, k+1)
        sig_x = np.sqrt(v_here)                        # (1, k+1)
        mu_x = -0.5 * v_here
        A = sig_x / sbar_x                             # in (0,1]
        alpha = 0.5 * (A ** 2 - 1.0)                   # eq (4.8)
        beta = 0.0                                     # sigma_y = 1
        # kernel eq (4.10) coefficients, broadcast over (l, m)
        base_a = (mu_x * sqd / sbar_x) / (4.0 * (1.0 + alpha))
        c_coef = (rho * A * 1.0) / (4.0 * (1.0 + alpha) * (1.0 + beta))
        new = np.empty((k + 1, k + 1, 2, 2))
        for ip in (0, 1):
            xi_prev = 2 * ip - 1
            a_c = base_a + (alpha * xi_prev) / (4.0 * (1.0 + alpha))
            for jp in (0, 1):
                zeta_prev = 2 * jp - 1
                b_c = (beta * zeta_prev) / (4.0 * (1.0 + beta))
                cont = np.zeros((k + 1, k + 1))
                for i in (0, 1):
                    si = 2 * i - 1
                    for j in (0, 1):
                        sj = 2 * j - 1
                        p = 0.25 + si * a_c + sj * b_c + si * sj * c_coef
                        p = np.broadcast_to(p, (k + 1, k + 1))
                        bad = (p < -1e-12) | (p > 1.0 + 1e-12)
                        if bad.any():
                            viol += int(bad.sum())
                            maxviol = max(maxviol, float(np.max(
                                np.maximum(-p[bad], p[bad] - 1.0))))
                        cont = cont + p * val[i:i + k + 1, j:j + k + 1, i, j]
                logS_k = np.log(S0) + (2 * li - k) * hx
                exer = np.maximum(KSTRIKE - np.exp(logS_k), 0.0)
                new[:, :, ip, jp] = np.maximum(
                    np.broadcast_to(exer, (k + 1, k + 1)), cont)
        val = new
        if pr is not None and (k % max(1, n // 20) == 0):
            pr.tick(pr.step, inner=f"n={n} k={k}")
    return float(val[0, 0, 0, 0]), viol, maxviol


def phase2(pr: Progress) -> dict:
    pr.phase("V3/V8 2-D lattice American put vs CRR reference",
             total=len(N_PHASE2) + 1)
    ref = crr_american_put(20_000, SIGMA)
    pr.result("crr_reference_n20000", round(ref, 6))
    pr.tick()
    rows = []
    for n in N_PHASE2:
        v, viol, mv = lattice_american_put(n, H=0.10, eta=0.0, rho=RHO, pr=pr)
        rows.append({"n": n, "lattice": round(v, 6),
                     "abs_err": round(abs(v - ref), 6),
                     "prob_violations": viol,
                     "max_violation": round(mv, 6)})
        pr.tick(note=f"n={n}")
        pr.partial.update({f"put_n{n}": round(v, 4),
                           f"err_n{n}": round(abs(v - ref), 5)})
        pr.log(f"n={n} value={v:.6f} err={abs(v-ref):.6f} viol={viol}")
    ns = np.array([r["n"] for r in rows], dtype=float)
    er = np.array([r["abs_err"] for r in rows], dtype=float)
    ok = er > 0
    slope = float(np.polyfit(np.log(T / ns[ok]), np.log(er[ok]), 1)[0]) if ok.sum() > 2 else float("nan")
    return {"reference_crr": round(ref, 6), "rows": rows,
            "observed_error_slope_in_delta": round(slope, 3)}


def phase3(pr: Progress, zmax: float, label: str) -> dict:
    pr.phase(label, total=len(N_PHASE3))
    vals, viols = {}, {}
    for n in N_PHASE3:
        v, viol, mv = lattice_american_put(n, H=0.10, eta=ETA, rho=RHO,
                                          zmax=zmax, pr=pr)
        viols[n] = viol
        vals[n] = v
        pr.tick(note=f"n={n}")
        pr.partial.update({f"put_n{n}": round(v, 4), f"viol_n{n}": viol})
        pr.log(f"eta={ETA} n={n} value={v:.6f} viol={viol} maxviol={mv:.3g}")
    ns = sorted(vals)
    diffs = [(ns[i], abs(vals[ns[i + 1]] - vals[ns[i]])) for i in range(len(ns) - 1)]
    d = np.array([T / n for n, _ in diffs], dtype=float)
    e = np.array([x for _, x in diffs], dtype=float)
    ok = e > 0
    r = float(np.polyfit(np.log(d[ok]), np.log(e[ok]), 1)[0]) if ok.sum() > 2 else float("nan")
    return {"zmax": (None if not np.isfinite(zmax) else round(float(zmax), 3)),
            "values": {str(k): round(v, 6) for k, v in vals.items()},
            "prob_violations": {str(k): int(x) for k, x in viols.items()},
            "successive_diffs": [[n, round(x, 6)] for n, x in diffs],
            "observed_rate_r": round(r, 3),
            "theoretical_H_over_2": 0.05}


# =========================================================================== #
def main() -> None:
    rng = np.random.default_rng(20260804)
    meta = {"example": "fractional driver + constant-vol American put",
            "T": T, "S0": S0, "K": KSTRIKE, "sigma": SIGMA, "rho": RHO,
            "eta_phase3": ETA}
    sd_zeta = float(np.sqrt(T ** (2 * 0.10) / (2 * 0.10)))
    zf = min(zmax_admissible(n, 0.10, ETA, RHO) for n in N_PHASE3)
    meta["zmax_asymptotic_4point"] = round(float(zmax_feasible(0.10, ETA, RHO, "4pt")), 3)
    meta["zmax_asymptotic_9point_Thm9.1"] = round(float(zmax_feasible(0.10, ETA, RHO, "9pt")), 3)
    meta["zmax_used_fixed_across_n"] = round(float(zf), 3)
    meta["zmax_used_in_sd_of_zeta_T"] = round(float(zf) / sd_zeta, 2)
    with Progress("validation", total_phases=4, meta=meta) as pr:
        r1 = phase1(pr, rng)
        pr.result("phase1_V1_V2", r1)
        r2 = phase2(pr)
        pr.result("phase2_V3_V8", r2)
        r3 = phase3(pr, np.inf,
                    "V13a rough Bergomi WITHOUT barriers (expected to diverge)")
        pr.result("phase3a_no_barriers", r3)
        r4 = phase3(pr, zf,
                    "V13b rough Bergomi WITH the admissible barrier")
        pr.result("phase3b_with_barriers", r4)

        lines = ["# Validation run — results", "",
                 f"Example: {meta['example']}", "",
                 "## Phase 1 (V1, V2) — the rough component", ""]
        lines.append("| H | slope of Var[Vcheck]/Var[V] vs n | predicted 1-2H | "
                     "Var[V^(n)]/Var[V] at n=4096 | Var[Vcheck]/Var[V] at n=4096 |")
        lines.append("|---|---|---|---|---|")
        for H in HURSTS:
            d = r1[f"H={H:.2f}"]
            lines.append(f"| {H:.2f} | {d['slope_ratio_onestep_vs_n']} | "
                         f"{d['predicted_slope_1_minus_2H']} | "
                         f"{d['ratio_exact_at_n4096']} | {d['ratio_onestep_at_n4096']} |")
        lines += ["", "## Phase 2 (V3, V8) — lattice vs CRR reference", "",
                  f"CRR reference (n=20000): **{r2['reference_crr']}**", "",
                  "| n | lattice | abs err | prob. violations |", "|---|---|---|---|"]
        for row in r2["rows"]:
            lines.append(f"| {row['n']} | {row['lattice']} | {row['abs_err']} | "
                         f"{row['prob_violations']} |")
        lines += ["", f"observed error slope in delta: "
                      f"{r2['observed_error_slope_in_delta']}", "",
                  "## Phase 3a (V13) — rough Bergomi, NO barrier", "",
                  "| n | value | prob. violations |", "|---|---|---|"]
        for n in N_PHASE3:
            lines.append(f"| {n} | {r3['values'][str(n)]} | "
                         f"{r3['prob_violations'][str(n)]} |")
        lines += ["", f"observed r = {r3['observed_rate_r']} — meaningless: the "
                      "value collapses, the scheme diverges.", "",
                  "## Phase 3b (V13) — rough Bergomi, WITH the admissible barrier",
                  "", f"zeta clamped to |zeta| <= {r4['zmax']} "
                      f"(= {meta['zmax_used_in_sd_of_zeta_T']} standard "
                      f"deviations of zeta_T; the 9-point kernel of Thm 9.1 "
                      f"would allow {meta['zmax_asymptotic_9point_Thm9.1']})", "",
                  "| n | value | prob. violations |", "|---|---|---|"]
        for n in N_PHASE3:
            lines.append(f"| {n} | {r4['values'][str(n)]} | "
                         f"{r4['prob_violations'][str(n)]} |")
        lines += ["", f"observed r = {r4['observed_rate_r']} "
                      f"(theory H/2 = {r4['theoretical_H_over_2']})", ""]
        pr.write_results_md("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
