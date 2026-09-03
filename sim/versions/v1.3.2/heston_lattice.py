"""
heston_lattice — the paper's lattice at h = 0, so that it can at last be measured
against a reference containing no Monte-Carlo.

Why this file exists.  The independent reference of `rough_heston` prices rough
Heston, and the lattice of `route_aprime` prices rough Bergomi, so the two cannot
be compared directly.  That is not an oversight but a structural fact, and it is
worth stating sharply:

    the paper's model class needs  v_t = F(t, (K^h Y)(t))  with Y an AUTONOMOUS
    diffusion, whereas rough Heston has  nu sqrt(V) dB,  so its driver's
    coefficients depend on the variance -- that is, on the driver's own past.
    Rough Heston is therefore OUTSIDE the model class, and no lattice of the
    paper's kind prices it.

Conversely, the models with a semi-analytic characteristic function are affine,
and in the rough setting affineness forces exactly that sqrt(V) dB structure.  So

    no ROUGH model both fits the paper's class and has a semi-analytic price,

and the lattice can never be checked against a Fourier reference in the rough
regime.  The best available is the boundary case h = 0, where the two classes
finally meet: at h = 0 the operator K^0 is the identity, so v = v_0 + y with

    dy = lambda(theta - v_0 - y) dt + nu sqrt(v_0 + y) dB,

which IS autonomous in y -- this is why Akyildirim--Dolinsky--Soner could build a
lattice for Heston at all -- and the same model is classical Heston, which the
Fourier pricer handles at alpha = 1.

What this test therefore buys, and what it does not.  It measures the lattice
construction -- moment matching, recombination, the two-dimensional backward
induction and the Route A' coupling of the price -- against a reference with no
Monte-Carlo in it, on a model genuinely inside the paper's class, and it exercises
a STATE-DEPENDENT sigma_y, which rough Bergomi (sigma_y = 1) never touches.  It
does not test the rough part: h = 0 is not rough.

Construction.  A unit-volatility driver is needed for a recombining walk, so the
variance is put through the Lamperti transform U = 2 sqrt(v) / nu, giving

    dU = [ (2 lambda theta / nu^2 - 1/2) / U - lambda U / 2 ] dt + dB,

with unit diffusion coefficient.  U lives on u_0 + sqrt(delta) Z and the drift is
carried by the up-probability p = (1 + mu_U sqrt(delta)) / 2.  The price is
coupled to the driver by the randomised rounding of Route A', at the finer scale
a_X, which is what removes the admissibility constraint.
"""
from __future__ import annotations

import math

import numpy as np

from route_aprime import interp_kernel

T = 1.0
S0 = 100.0
KSTRIKE = 100.0


def lamperti(v: np.ndarray, nu: float) -> np.ndarray:
    return 2.0 * np.sqrt(np.maximum(v, 0.0)) / nu


def inv_lamperti(u: np.ndarray, nu: float) -> np.ndarray:
    return 0.25 * nu * nu * u * u


def drift_U(u: np.ndarray, lam: float, theta: float, nu: float) -> np.ndarray:
    """The Lamperti drift, with the 1/U term guarded away from the origin."""
    c = 2.0 * lam * theta / (nu * nu) - 0.5
    return c / np.maximum(u, 1e-12) - 0.5 * lam * u


def cir_moments(V0: float, theta: float, lam: float, nu: float,
                t: float) -> tuple[float, float]:
    """Exact mean and standard deviation of v_t under the CIR dynamics.

        E[v_t]   = theta + (V0 - theta) e^{-lam t}
        Var[v_t] = V0 (nu^2/lam)(e^{-lam t} - e^{-2 lam t})
                   + theta (nu^2/(2 lam))(1 - e^{-lam t})^2

    Used to place the variance barrier in units the parameters themselves set, so
    that a given barrier width means the same thing at every (lam, nu, theta).
    """
    e = math.exp(-lam * t)
    mean = theta + (V0 - theta) * e
    var = (V0 * (nu * nu / lam) * (e - e * e)
           + theta * (nu * nu / (2.0 * lam)) * (1.0 - e) ** 2)
    return mean, math.sqrt(max(var, 0.0))


def barrier_U(V0: float, theta: float, lam: float, nu: float, t: float,
              width_sd: float, delta: float | None = None,
              drift_safety: float = 2.0, nscan: int = 64) -> tuple[float, float]:
    """Absorbing band for the Lamperti driver, at +-width_sd of v_s over s in [0,t].

    The band is taken over the WHOLE horizon and not at the final date alone.  An
    earlier version used the moments of v_t only, which is safe when V0 = theta --
    the case Section 10.8.3 tests, where the mean is constant -- and wrong as soon
    as V0 is far from theta: the mean reverts away from V0, so a band built around
    E[v_t] can fail to contain v_0 itself, the driver is clipped from the first
    step, and the variance is held below its own initial value.  On the
    Beliaeva--Nawalkha grid with sqrt(V0) = 0.4 and theta = 0.04 the old band
    topped out at 0.1244 against V0 = 0.16 at T = 6 months, and the put came out
    0.26 too low.  Scanning s in [0,t] and including s = 0 -- where the mean is V0
    and the standard deviation is zero -- makes v_0 an endpoint of the band by
    construction, so it can never be excluded.

    Two constraints fight over the LOWER end.  The variance cannot be negative, so
    the band would like to reach down to u = 0; but the Lamperti drift carries a
    term c/u with c = 2 lambda theta / nu^2 - 1/2, so at u near zero the drift
    diverges and the binomial up-probability
    p = (1 + mu sqrt(delta))/2 leaves [0,1].  A one-step +-sqrt(delta) walk simply
    cannot represent a drift larger than 1/sqrt(delta).

    We therefore floor the band at u_lo = drift_safety * c * sqrt(delta), the
    smallest level at which |mu sqrt(delta)| <= 1/drift_safety, when delta is
    supplied.  The cost is a variance floor that is higher on coarse grids and
    vanishes as delta -> 0, which is the right direction; the alternative -- what
    an earlier version of this file did -- is to clip the offending probabilities
    and leave a third of the driver states with the wrong dynamics.
    """
    v_lo, v_hi = math.inf, -math.inf
    for i in range(nscan + 1):
        mean, sd = cir_moments(V0, theta, lam, nu, t * i / nscan)
        v_lo = min(v_lo, mean - width_sd * sd)
        v_hi = max(v_hi, mean + width_sd * sd)
    v_lo = max(0.0, v_lo)
    u_lo = float(lamperti(np.array([v_lo]), nu)[0])
    u_hi = float(lamperti(np.array([v_hi]), nu)[0])
    if delta is not None:
        c = 2.0 * lam * theta / (nu * nu) - 0.5
        u_lo = max(u_lo, drift_safety * abs(c) * math.sqrt(delta))
    return max(u_lo, 1e-6), min(u_hi, 1e6)


def heston_put_lattice(n: int, V0: float, theta: float, lam: float, nu: float,
                       rho: float, T_: float = T, s0: float = S0,
                       k: float = KSTRIKE, mref: int = 4,
                       american: bool = False,
                       barrier_sd: float | None = None,
                       drift_floor: bool = False,
                       r: float = 0.0, payoff_type: str = "put") -> dict:
    """Two-dimensional lattice for classical Heston, Route A' coupling.

    The driver walk is +-sqrt(delta) in the Lamperti variable, so it recombines
    and its state after k steps is the number of up-moves.  The price lives on a
    grid of spacing a_X = sqrt(V0) sqrt(delta) / mref and every conditional law is
    placed on it by the interpolation kernel of Route A', which is a probability
    with the exact first moment for free.

    `barrier_sd` absorbs the driver at that many standard deviations of v_{T},
    exactly as Section 4 absorbs Y at B^Y and C^Y.  Without it the walk reaches
    u_0 + sqrt(delta) n, which at n = 64 means a variance of 2.25 -- a volatility
    of 150 per cent -- and the price grid is sized for that absurdity: 30209 nodes
    against 10233 at four standard deviations.  Since the width of the grid enters
    the cost twice, through the number of nodes and through the number of offsets,
    the saving is quadratic.

    `drift_floor` selects between the two ways of regularising the 1/U drift near
    the origin.  With it False the band reaches down to zero and the handful of
    up-probabilities that leave [0,1] are clipped; with it True the band is floored
    at the level where the binomial walk can still represent the drift, which
    removes every violation but imposes a variance floor that bites on coarse
    grids.  The two agree for n >= 32 and differ below, so the rate must be fitted
    on n >= 32 where the answer does not depend on this choice.

    `r` is the risk-free rate.  It enters twice and both are needed: the log-price
    drift gains `r delta` per step, and the continuation value is discounted by
    `exp(-r delta)` per step.  With the default `r = 0` every arithmetic operation
    is unchanged from the version that produced Section 10.8.3, so the numbers
    recorded there are not disturbed.

    `payoff_type` is `"put"` or `"call"`.  For an American call with `r >= 0` and
    no dividend the early-exercise feature is worthless, so `american=True` with
    `payoff_type="call"` must return the European call: a free control on the
    whole backward induction, which the run exercises.
    """
    if payoff_type not in ("put", "call"):
        raise ValueError(f"payoff_type must be 'put' or 'call', got {payoff_type!r}")
    d = T_ / n
    sqd = math.sqrt(d)
    u0 = float(lamperti(np.array([V0]), nu)[0])

    # driver states: u(k, m) = u0 + sqrt(d) (2m - k), m = 0..k
    def u_of(kk: int, m: np.ndarray) -> np.ndarray:
        return u0 + sqd * (2.0 * m - kk)

    # price grid, wide enough for the largest conditional spread
    if barrier_sd is None:
        u_lo, u_hi = 1e-3, u0 + sqd * n
    else:
        u_lo, u_hi = barrier_U(V0, theta, lam, nu, T_, barrier_sd,
                               delta=d if drift_floor else None)
    v_max = float(inv_lamperti(np.array([min(u_hi, u0 + sqd * n)]), nu)[0])
    s_ref = math.sqrt(V0)
    a_X = s_ref * sqd / mref
    sig_max = math.sqrt(max(v_max, V0))
    kmax = int(math.ceil(4.5 * sig_max * math.sqrt(1.0 - rho ** 2) * sqd / a_X
                         + abs(rho) * sig_max * sqd / a_X)) + 1
    offsets = np.arange(-kmax, kmax + 1)
    reach = kmax * n
    nx = 2 * reach + 1
    ix0 = reach
    logS = math.log(s0) + (np.arange(nx) - ix0) * a_X
    payoff = (np.maximum(k - np.exp(logS), 0.0) if payoff_type == "put"
              else np.maximum(np.exp(logS) - k, 0.0))

    val = np.repeat(payoff[:, None], n + 1, axis=1)
    disc = math.exp(-r * d)
    neg = 0
    for kk in range(n - 1, -1, -1):
        m = np.arange(kk + 1, dtype=float)
        # absorption: the effective driver is the clamped one, and it is used for
        # the drift as well as the variance, so the dynamics stay consistent
        u_here = np.clip(u_of(kk, m), u_lo, u_hi)
        v_here = inv_lamperti(u_here, nu)
        sig = np.sqrt(v_here)
        mu = drift_U(u_here, lam, theta, nu)
        p_up = 0.5 * (1.0 + mu * sqd)
        bad = (p_up < 0.0) | (p_up > 1.0)
        neg += int(bad.sum())
        p_up = np.clip(p_up, 0.0, 1.0)
        lo, hi = max(0, ix0 - kmax * kk), min(nx, ix0 + kmax * kk + 1)
        cont = np.zeros((hi - lo, kk + 1))
        for j, zeta in ((1, +1.0), (0, -1.0)):
            w = p_up if zeta > 0 else (1.0 - p_up)
            # compensate the driver drift so that E[dX] = mu_x d exactly
            zc = zeta - mu * sqd
            mu_M = sig * rho * sqd * zc - 0.5 * v_here * d + r * d
            sd_M = sig * math.sqrt(1.0 - rho ** 2) * sqd
            P = interp_kernel(mu_M, sd_M, a_X, offsets)
            nxt = val[:, j:j + kk + 1]
            for oi, o in enumerate(offsets):
                src = np.clip(np.arange(lo, hi) + o, 0, nx - 1)
                cont += w[None, :] * P[oi][None, :] * nxt[src, :]
        nv = np.zeros((nx, kk + 1))
        nv[lo:hi, :] = disc * cont
        nv[:lo, :] = val[:lo, :kk + 1]
        nv[hi:, :] = val[hi:, :kk + 1]
        if american:
            ex = payoff[:, None]
            nv = np.maximum(nv, np.repeat(ex, kk + 1, axis=1))
        val = nv
    return {"value": float(val[ix0, 0]), "n": n, "mref": mref,
            "a_X": a_X, "offsets": len(offsets), "grid": nx,
            "driver_probability_violations": neg,
            "u0": u0, "u_lo": u_lo, "u_hi": u_hi,
            "barrier_sd": barrier_sd,
            "v_at_upper_barrier": float(inv_lamperti(np.array([u_hi]), nu)[0]),
            "walk_reach_u": u0 + sqd * n,
            "feller_ratio": 2.0 * lam * theta / (nu * nu)}
