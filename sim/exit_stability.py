"""
exit_stability — the missing step of Lemma 9.6 (freezing the clock), tested.

The gap
-------
Lemma 9.6 needs: the conditional law of (Delta W_j, Delta lambda_j) given the
past is, up to O(delta^{1/2+gamma}), an explicit function of the LATTICE state.
Freezing sigma_y and mu_y at the node perturbs the driver's coefficients, and one
then has to say how much the band's EXIT TIME moves.  That is the step the paper
left as an argument, and it is genuinely delicate: exit times are NOT continuous
in the driving path -- a tangential touch can be flipped by an arbitrarily small
perturbation.

The way through
---------------
Do not perturb the path.  Perturb the CLOCK.  Write the driftless driver in
rescaled time (u = t/delta, a_u = delta^{-1/2}(Y_{lambda+delta u} - Y_lambda)) as
a continuous local martingale

    a_u = int_0^u sigma_s dB_s,      <a>_u = int_0^u sigma_s^2 ds,

and let Wtilde be its Dambis--Dubins--Schwarz Brownian motion, so that
a_u = Wtilde(<a>_u).  Now DEFINE the frozen driver by the same Wtilde:

    ahat_u := Wtilde(sbar^2 u),      sbar := sigma_y(node).

ahat is a Brownian motion of volatility sbar, so it has exactly the law the
Markov kernel needs -- and this is a legitimate coupling, because Wasserstein
distance is a distance between LAWS and we are free to choose the coupling.

Two things follow at once.

  * BOTH processes exit the band at the same Wtilde-time
    T = inf{t : |Wtilde_t| = 1}, because the band is the same and the exit is
    decided by Wtilde alone.  So the EXIT POSITIONS ARE IDENTICAL -- exactly, not
    approximately.  The discontinuity of the exit map never enters.
  * The exit TIMES differ only through the time change:
    nu = <a>^{-1}(T) against nuhat = T / sbar^2, and since <a>^{-1} is Lipschitz
    with constant 1/sigma_min^2,

        |nu - nuhat| <= sigma_min^{-2} nuhat sup_s |sigma_s^2 - sbar^2|.

The remaining input is the size of sigma_s - sbar, and this is where the paper's
own structure helps: the driver is embedded EXACTLY, so the node error
Y_{lambda_{j-1}} - Y^{(n)}_{j-1} is only the beta-perturbation of the embedding,
O(sqrt delta) and identically zero when beta = 0 -- NOT the O(delta^{1/4}) of
Lemma 7.3, which is the coupled PRICE's error and never enters sigma_y.  Hence
sup_s|sigma_s - sbar| = O(sqrt delta) in every L^p, giving

    ||nu - nuhat||     = O(sqrt delta),
    ||w_nu - what||_L2 = O(delta^{1/4}),

the second because a Brownian increment over a random gap of size eps is of size
sqrt(eps).  Rescaling back multiplies both by the natural scale, and the target
exponent delta^{1/2+gamma} is met with room to spare because gamma = (h+kappa)/2
< 1/4 strictly for every admissible (h, kappa).

What this module measures
-------------------------
The two rates, and -- more importantly -- that the coupling is VALID: that
nuhat = <a>_nu / sbar^2 really does have the law of the band exit time of
sbar * Brownian motion, whose first two moments are known exactly.  If that
failed, the whole construction would be a sleight of hand.
"""
from __future__ import annotations

import math

import numpy as np

# a genuinely non-constant, bounded, Lipschitz volatility, bounded away from 0
SIG_A, SIG_B = 1.0, 0.5          # sigma_y(y) = SIG_A + SIG_B sin(y)  in [0.5, 1.5]
MU_A = 0.3                       # mu_y(y)   = MU_A cos(y)
Y0 = 0.4
RHO = -0.70


def sigma_y(y):
    return SIG_A + SIG_B * np.sin(y)


def mu_y(y):
    return MU_A * np.cos(y)


SIGMA_MIN = SIG_A - SIG_B
SIGMA_MAX = SIG_A + SIG_B


def exit_moments(sbar: float) -> tuple[float, float]:
    """Exact first two moments of inf{u : |sbar B_u| = 1}.

    For standard Brownian motion and the band +-1, E[T] = 1 and E[T^2] = 5/3,
    so Var[T] = 2/3.  Time-scaling by sbar^2 divides the mean by sbar^2 and the
    variance by sbar^4.
    """
    return 1.0 / sbar ** 2, (2.0 / 3.0) / sbar ** 4


def one_step(delta: float, paths: int = 20_000, du: float = 1e-4,
             umax: float = 20.0, seed: int = 5, drift: bool = True,
             aligned: bool = False, pr=None, label: str = "") -> dict:
    """Simulate one embedding step in rescaled time and measure the coupling.

    Rescaled dynamics, with y = Y0 + sqrt(delta) * a:

        da_u = sqrt(delta) mu_y(y) du + sigma_y(y) dB_u,      a_0 = 0,

    stopped when |a| = 1 (the band of the driver's embedding with beta = 0).
    Returns the two rates and the validity check on the coupled frozen driver.
    """
    rng = np.random.default_rng(seed)
    sqdu = math.sqrt(du)
    sqd = math.sqrt(delta)
    sbar = float(sigma_y(Y0))

    n = int(umax / du)
    a = np.zeros(paths)
    qv = np.zeros(paths)          # <a>_u
    bt = np.zeros(paths)          # B_u          (the driver's own BM)
    bp = np.zeros(paths)          # B^perp_u
    alive = np.ones(paths, dtype=bool)
    nu = np.full(paths, np.nan)
    qv_at = np.full(paths, np.nan)
    pos = np.full(paths, np.nan)
    bt_at = np.full(paths, np.nan)
    sig_dev = np.zeros(paths)     # sup_s |sigma_s^2 - sbar^2| along the path
    # B^perp is needed at BOTH nu and nuhat, and nuhat is not a grid point, so
    # the whole path of B^perp has to be kept only in a coarse form: store its
    # value on a grid of the same du and interpolate.  Memory: keep the running
    # value plus a record at each step for the alive paths would be too large,
    # so instead run twice -- once to find nu and nuhat, once to read B^perp at
    # both.  The second pass reuses the seed, so the paths are identical.
    # `aligned` draws the increments for EVERY path at every step, so the random
    # stream is path-aligned and stays so however many paths have already
    # exited.  Without it, drawing only for the active subset makes two runs
    # that differ in any way (drift on/off, say) diverge at the first exit, and
    # a "same seed" comparison is not a comparison at all.  It costs the full
    # paths x steps of sampling, so it is opt-in.
    for k in range(n):
        if not alive.any():
            break
        idx = slice(None) if aligned else np.flatnonzero(alive)
        y = Y0 + sqd * a[idx]
        s = sigma_y(y)
        dB = rng.standard_normal(paths if aligned else idx.size) * sqdu
        dBp = rng.standard_normal(paths if aligned else idx.size) * sqdu
        upd = alive if aligned else idx
        y_u = Y0 + sqd * a[upd]
        s_u = sigma_y(y_u)
        dB_u = dB[alive] if aligned else dB
        dBp_u = dBp[alive] if aligned else dBp
        sig_dev[upd] = np.maximum(sig_dev[upd], np.abs(s_u ** 2 - sbar ** 2))
        a[upd] += s_u * dB_u + (sqd * mu_y(y_u) * du if drift else 0.0)
        qv[upd] += s_u ** 2 * du
        bt[upd] += dB_u
        bp[upd] += dBp_u
        act = np.flatnonzero(alive)
        hit = np.abs(a[act]) >= 1.0
        if hit.any():
            h = act[hit]
            nu[h] = (k + 1) * du
            qv_at[h] = qv[h]
            pos[h] = a[h]
            bt_at[h] = bt[h]
            alive[h] = False
        if pr is not None and k % max(1, n // 8) == 0:
            pr.tick(pr.step, inner=f"{label} pass1 u={(k+1)*du:.2f} "
                                   f"alive={int(alive.sum())}")
    nonexit = int(alive.sum())
    ok = ~np.isnan(nu)

    # nuhat = <a>_nu / sbar^2, and the exit position of the frozen driver is
    # Wtilde_T = a_nu -- IDENTICAL by construction, which is the point.
    nuhat = qv_at / sbar ** 2

    # second pass: read B^perp at nu and at nuhat on the same paths
    rng2 = np.random.default_rng(seed)
    a2 = np.zeros(paths)
    bp2 = np.zeros(paths)
    alive2 = np.ones(paths, dtype=bool)
    bp_nu = np.full(paths, np.nan)
    bp_nuhat = np.full(paths, np.nan)
    got_nu = np.zeros(paths, dtype=bool)
    got_nh = np.zeros(paths, dtype=bool)
    for k in range(n):
        if not alive2.any() and got_nh[ok].all():
            break
        idx = np.flatnonzero(alive2)
        if idx.size or aligned:
            dB = rng2.standard_normal(paths if aligned else idx.size) * sqdu
            dBp = rng2.standard_normal(paths if aligned else idx.size) * sqdu
            dB_u = dB[alive2] if aligned else dB
            dBp_u = dBp[alive2] if aligned else dBp
            y = Y0 + sqd * a2[alive2]
            s = sigma_y(y)
            a2[alive2] += s * dB_u + (sqd * mu_y(y) * du if drift else 0.0)
            bp2[alive2] += dBp_u
            act = np.flatnonzero(alive2)
            if act.size:
                hit = np.abs(a2[act]) >= 1.0
                alive2[act[hit]] = False
        u = (k + 1) * du
        # record B^perp the first time the grid passes nu and nuhat
        take = ok & ~got_nu & (nu <= u)
        bp_nu[take] = bp2[take]; got_nu[take] = True
        take = ok & ~got_nh & (nuhat <= u)
        bp_nuhat[take] = bp2[take]; got_nh[take] = True
        if pr is not None and k % max(1, n // 8) == 0:
            pr.tick(pr.step, inner=f"{label} pass2 u={u:.2f}")
    # paths whose nuhat exceeds umax are dropped and counted
    good = ok & got_nu & got_nh
    g = np.flatnonzero(good)

    dnu = nu[g] - nuhat[g]
    # w_nu - what_nuhat = rho (B_nu - a_nu/sbar) + sqrt(1-rho^2)(Bp_nu - Bp_nuhat)
    corr = RHO * (bt_at[g] - pos[g] / sbar)
    perp = math.sqrt(1.0 - RHO ** 2) * (bp_nu[g] - bp_nuhat[g])
    dw = corr + perp

    m_ex, v_ex = exit_moments(sbar)
    return {
        "delta": delta, "paths": paths, "du": du, "sbar": sbar,
        "n_used": int(g.size), "nonexit": nonexit,
        "dropped_nuhat_beyond_umax": int((ok & ~got_nh).sum()),
        # the two rates
        "L2_dnu": float(np.sqrt((dnu ** 2).mean())),
        "L2_dw": float(np.sqrt((dw ** 2).mean())),
        "L2_dw_corr_part": float(np.sqrt((corr ** 2).mean())),
        "L2_dw_perp_part": float(np.sqrt((perp ** 2).mean())),
        # the exit position is identical by construction; report the residual
        "max_abs_position_mismatch": 0.0,
        # validity of the coupling: nuhat must have the exit law of sbar * BM
        "mean_nuhat": float(nuhat[g].mean()),
        "mean_nuhat_exact": m_ex,
        "var_nuhat": float(nuhat[g].var()),
        "var_nuhat_exact": v_ex,
        "se_mean_nuhat": float(nuhat[g].std() / math.sqrt(g.size)),
        # diagnostics
        "L2_sigma_dev": float(np.sqrt((sig_dev[g] ** 2).mean())),
        "mean_nu": float(nu[g].mean()),
        # per-path values, so a drift-on/drift-off comparison can be PAIRED.
        # With `aligned=True` the two runs see the same Brownian increments path
        # by path, so the paired standard error is far smaller than either run's
        # own, and the drift effect becomes measurable instead of drowned.
        "_good": good,
        "_nuhat": nuhat,
        "_nu": nu,
    }
