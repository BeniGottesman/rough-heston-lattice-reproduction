"""
mc_reference — a Monte-Carlo reference for rough Bergomi good enough to arbitrate
the tree, which the previous one was not.

Three defects of `route_aprime.european_put_mc_rough` are removed.

1.  It discretised the fractional driver with the LEFT-ENDPOINT kernel
    ker_j = sqrt(2H)(j df)^h.  Finding F022 of this project measures exactly what
    that costs: the relative variance error is O(df^{2H}), and at H = 0.1 with
    2048 steps it still loses **19%** of Var[V^H_T].  A reference carrying a 19%
    error in the driver cannot adjudicate a 0.15 discrepancy in the price.
    Here the driver is sampled with its EXACT covariance: the coefficients are
    the Cholesky factor of the true covariance matrix
    Cov(V^H_s, V^H_t) = 2H int_0^{s^t}(s-u)^h (t-u)^h du, so
    Var[V^H_t] = t^{2H} and every cross-covariance holds to machine precision,
    at every number of steps.  The factor stays lower-triangular in the SAME
    normals that drive the price, so the correlation rho is exact too.

2.  No control variate.  Setting eta = 0 in the same paths gives a Black--Scholes
    payoff whose expectation is known in closed form, and the two payoffs are
    strongly correlated, so
        estimator = mean(payoff_rough - payoff_flat) + BS_exact
    is unbiased and has a far smaller variance.  It also cancels most of the
    log-Euler error in the price leg, since that error is common to both.
    At eta = 0 the two payoffs coincide pathwise, so the estimator returns
    BS_exact with zero variance — a built-in consistency check.

3.  No antithetic sampling.  Each chunk now uses (Z, Z') and (-Z, -Z').

O018, repaired.  Antithetic pairing makes the paths dependent, so the standard error
must be measured on the PAIR AVERAGES and not on the paths.  Writing `rho` for the
within-pair correlation, `SE_correct = SE_iid * sqrt(1 + rho)`, and the sign of `rho`
is not a detail: on the raw put the pairing works and `rho ~ -0.58`, but the eta = 0
control variate cancels the monotone dependence on the driver -- exactly the
antisymmetric part the pairing exploited -- leaving an even residual with `rho ~ +0.19`.
The two variance-reduction devices are therefore NOT independent, and with both switched
on the old i.i.d. formula understated the band by 5 to 9 per cent.  `stderr` is now the
pair-level figure; the old one is returned beside it as `stderr_iid_SUPERSEDED` so the
correction stays visible.  The PRICE is unchanged, and `verify/verify_mc_band.py` fails
on the old formula.

What remains is the time-stepping of the log price only: v is taken at the LEFT
endpoint of each interval (pairing it with the same interval's dB would be
anticipating and biases the price upward through rho — the 27.24-instead-of-11.9
bug of an earlier session).  That error is O(df) and is what `nfine` controls.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.integrate import quad
from scipy.stats import norm

T = 1.0
S0 = 100.0
KSTRIKE = 100.0
XI0 = 0.30 ** 2


def bs_put(s0: float = S0, k: float = KSTRIKE, sigma: float = math.sqrt(XI0),
           t: float = T) -> float:
    """Black--Scholes put, zero rate — the control variate's exact mean."""
    sd = sigma * math.sqrt(t)
    d1 = (math.log(s0 / k) + 0.5 * sd * sd) / sd
    return float(k * norm.cdf(-(d1 - sd)) - s0 * norm.cdf(-d1))


def cov_VH(H: float, s: float, t: float) -> float:
    """Cov(V^H_s, V^H_t) = 2H int_0^{s^t} (s-u)^h (t-u)^h du, h = H - 1/2.

    Normalised so that Var[V^H_t] = t^{2H}, which is the convention the tree
    and `route_aprime` use.
    """
    h = H - 0.5
    lo = min(s, t)
    if lo <= 0.0:
        return 0.0
    if abs(s - t) < 1e-15:
        return float(s ** (2.0 * H))
    # u = min - x, then x = z^p with p = 1/(1+h) so that x^h dx = p dz:
    # the x^h singularity at the near end is removed exactly and the integrand
    # is bounded, which plain adaptive quadrature on the raw form is not.
    c = abs(s - t)
    p = 1.0 / (1.0 + h)
    zmax = lo ** (1.0 / p)
    val, _ = quad(lambda z: (z ** p + c) ** h * p, 0.0, zmax,
                  limit=400, epsabs=1e-14, epsrel=1e-12)
    return float(2.0 * H * val)


_CHOL: dict = {}


def driver_factor(H: float, nfine: int) -> np.ndarray:
    """Lower-triangular L with L L^T = exact covariance of (V^H_{t_1..t_n}).

    Cached: it depends only on (H, nfine), not on eta or rho.
    """
    key = (round(H, 12), nfine)
    if key in _CHOL:
        return _CHOL[key]
    df = T / nfine
    ts = np.arange(1, nfine + 1) * df
    C = np.empty((nfine, nfine))
    for i in range(nfine):
        for j in range(i, nfine):
            C[i, j] = C[j, i] = cov_VH(H, float(ts[i]), float(ts[j]))
    # tiny jitter only if the exact matrix is numerically indefinite
    jit = 0.0
    scale = float(np.trace(C)) / nfine
    for _ in range(10):
        try:
            L = np.linalg.cholesky(C + jit * scale * np.eye(nfine))
            break
        except np.linalg.LinAlgError:
            jit = 1e-14 if jit == 0.0 else jit * 10.0
    else:
        raise np.linalg.LinAlgError("driver covariance not factorisable")
    _CHOL[key] = L
    return L


def european_put_mc(H: float, eta: float, rho: float, nfine: int = 1024,
                    paths: int = 200_000, chunk: int = 10_000, seed: int = 7,
                    control: bool = True, exact_driver: bool = True,
                    cap: float | None = None, pr=None) -> dict:
    """Rough Bergomi European put by Monte-Carlo.

    exact_driver=False falls back to the old left-endpoint kernel, so the cost of
    that choice can be measured rather than asserted.
    Returns price, standard error, and the diagnostics needed to trust it.
    """
    rng = np.random.default_rng(seed)
    df = T / nfine
    h = H - 0.5
    ts = np.arange(1, nfine + 1) * df
    bs = bs_put()

    if exact_driver:
        L = driver_factor(H, nfine)
    else:
        ker = np.sqrt(2.0 * H) * ts ** h
        Lp = 1
        while Lp < 2 * nfine:
            Lp *= 2
        fker = np.fft.rfft(ker, Lp)

    tot = tot2 = 0.0
    ptot = ptot2 = 0.0
    pn = 0
    done = 0
    var_acc = 0.0                     # to report the driver's realised variance
    while done < paths:
        p = min(chunk, paths - done)
        half = (p + 1) // 2
        Z = rng.standard_normal((half, nfine))
        Zp = rng.standard_normal((half, nfine))
        Z = np.concatenate([Z, -Z], axis=0)[:p]          # antithetic
        Zp = np.concatenate([Zp, -Zp], axis=0)[:p]

        if exact_driver:
            with np.errstate(all="ignore"):    # matmul raises spurious FE flags
                VH = Z @ L.T                              # exact covariance
            assert np.all(np.isfinite(VH)), "non-finite driver sample"
        else:
            dW_ = Z * np.sqrt(df)
            VH = np.fft.irfft(np.fft.rfft(dW_, Lp, axis=1) * fker[None, :],
                              Lp, axis=1)[:, :nfine]
        var_acc += float((VH[:, -1] ** 2).sum())
        if cap is not None:
            # Mirrors the tree's variance barrier.  The tree clips its discrete
            # driver at zmax = cap/sqrt(2H), so that eta*sqrt(2H)*zc lands in
            # [-cap*eta, +cap*eta]; clipping V^H at +-cap puts the log-variance
            # in the same window, which makes tree and reference comparable and
            # isolates the absorption error from the scheme error.
            VH = np.clip(VH, -cap, cap)

        v = XI0 * np.exp(eta * VH - 0.5 * eta ** 2 * ts[None, :] ** (2.0 * H))
        # LEFT endpoint: v on (t_{k-1}, t_k] is v(t_{k-1}).  Pairing v[k] with
        # dB[k] would anticipate and bias the price up through rho.
        v = np.concatenate([np.full((p, 1), XI0), v[:, :-1]], axis=1)

        dW = Z * np.sqrt(df)
        dB = rho * dW + math.sqrt(1.0 - rho ** 2) * Zp * np.sqrt(df)
        sB = dB.sum(axis=1)

        logST = (math.log(S0) - 0.5 * (v * df).sum(axis=1)
                 + (np.sqrt(v) * dB).sum(axis=1))
        pay = np.maximum(KSTRIKE - np.exp(logST), 0.0)

        if control:
            logflat = math.log(S0) - 0.5 * XI0 * T + math.sqrt(XI0) * sB
            payflat = np.maximum(KSTRIKE - np.exp(logflat), 0.0)
            sample = pay - payflat
        else:
            sample = pay
        tot += sample.sum()
        tot2 += (sample ** 2).sum()
        # O018.  The paths are drawn in ANTITHETIC PAIRS -- path i is the mirror of
        # path i + half -- so they are not independent and the path-level i.i.d.
        # standard error below is wrong.  The estimator is the mean of the pair
        # averages, so its dispersion must be measured on those.  Accumulated here
        # rather than derived later, because `sample` does not survive the chunk.
        npair = p - half                      # complete pairs in this chunk
        if npair > 0:
            pa = 0.5 * (sample[:npair] + sample[half:half + npair])
            ptot += float(pa.sum()); ptot2 += float((pa ** 2).sum()); pn += npair
        if half > npair:                      # odd `p`: unpaired primals stand alone
            res = sample[npair:half]
            ptot += float(res.sum()); ptot2 += float((res ** 2).sum()); pn += res.size
        done += p
        if pr is not None:
            pr.tick(note=f"mc {done}/{paths}")

    mean = tot / paths                    # the point estimate is UNCHANGED
    # the naive figure, kept only so the correction is visible and testable
    se_iid = math.sqrt(max(0.0, tot2 / paths - mean ** 2) / paths)
    if pn < 2:
        raise AssertionError("fewer than two variance units: no band can be reported")
    pmean = ptot / pn
    pvar = max(0.0, ptot2 / pn - pmean ** 2)
    se = math.sqrt(pvar / pn)             # of the MEAN OF PAIR AVERAGES
    # rho_pair is recovered from the two variances rather than stored per path:
    # Var(A) = Var(s)(1 + rho)/2, so rho = 2 Var(A)/Var(s) - 1.
    var_s = max(0.0, tot2 / paths - mean ** 2)
    rho_pair = (2.0 * pvar / var_s - 1.0) if var_s > 0.0 else 0.0
    price = mean + bs if control else mean
    return {"price": float(price), "stderr": float(se),
            "ci95": float(1.96 * se),
            "stderr_iid_SUPERSEDED": float(se_iid),
            "rho_pair": float(rho_pair),
            "band_inflation": float(se / se_iid) if se_iid > 0.0 else 1.0,
            "variance_units": int(pn),
            "variance_unit": "antithetic pair",
            "paths": paths, "nfine": nfine,
            "control": control, "exact_driver": exact_driver,
            "bs_control_mean": bs,
            "driver_var_T_realised": var_acc / paths,
            "driver_var_T_target": T ** (2.0 * H),
            "driver_var_ratio": (var_acc / paths) / T ** (2.0 * H)}
