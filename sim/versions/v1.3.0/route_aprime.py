"""
route_aprime — the hybrid scheme of paper §9.2 (Route A'), implemented.

Construction actually used
--------------------------
The driver is embedded exactly: with sigma_y = 1, mu_y = 0 and beta = 0 the
one-dimensional embedding gives, at its own clock,

    Delta N_j = Delta B_j = sqrt(delta) * zeta_j,     zeta_j = +-1 w.p. 1/2,

so the Y / zeta coordinate is a symmetric walk on the sqrt(delta) lattice, and
the rough coordinate is varsigma = (2m - k) * delta^H, exactly as in the
synchronous scheme.

The price is COUPLED, not embedded.  Writing W = rho B + sqrt(1-rho^2) B^perp
with B^perp independent of B (hence of zeta_j and of the clock),

    M_j := Delta M_j + delta * mu_x
         = sigma_x * ( rho sqrt(delta) zeta_j + sqrt(1-rho^2) sqrt(delta) G_j )
           + delta * mu_x,        G_j ~ N(0,1) independent,

and the lattice increment is the randomised rounding of M_j to a_X * Z:

    xi_j = rd_{a_X}(M_j ; U_j),      a_X := s_ref * sqrt(delta) / mref .

`mref` is the refinement factor of the price lattice; the theory takes
mref = delta^{-gamma}, and the accumulated rounding error is
O( sqrt(n a_X^2) ) = O( s_ref sqrt(T) / mref ), i.e. controlled by mref alone.

Because randomised rounding is mean-preserving, the transition law in lattice
units is the *linear-interpolation* kernel

    p_o(state, zeta) = E[ Lambda( M_j / a_X - o ) ],   Lambda(x) = (1-|x|)^+ ,

which satisfies sum_o p_o = 1 and sum_o o * p_o = E[M_j]/a_X EXACTLY.  Two
consequences, both tested in `run_route_aprime.py`:

  * every p_o is a probability BY CONSTRUCTION (Lambda >= 0, partition of
    unity), so the admissibility condition of Lemma 4.8 — the one that caps the
    barrier at 0.84 standard deviations for the four-point kernel — simply does
    not arise;
  * the correlation is carried by the true increment: E[a_X xi_j * a_Y zeta_j]
    = E[M_j Delta N_j] = rho sigma_x sigma_y delta exactly, with no condition on
    the step ratio.

State space: (price node, rough node) only — no (xi_prev, zeta_prev) memory is
needed here since beta = 0 and the price kernel depends on the past through the
state alone.
"""
from __future__ import annotations

import numpy as np

T = 1.0
S0 = 100.0
KSTRIKE = 100.0
XI0 = 0.30 ** 2
GH_NODES = 25                      # Gauss-Hermite nodes for E[Lambda(...)]


def _hermite():
    x, w = np.polynomial.hermite_e.hermegauss(GH_NODES)
    return x, w / w.sum()           # standard normal nodes / weights


_GX, _GW = _hermite()


def interp_kernel(mu: np.ndarray, sd: np.ndarray, a: float,
                  offsets: np.ndarray) -> np.ndarray:
    """p_o = E[Lambda(M/a - o)] for M ~ N(mu, sd^2), vectorised.

    mu, sd : shape (S,)      offsets : shape (K,)
    returns  shape (K, S), columns summing to 1 (up to quadrature error).
    """
    z = (mu[None, :] + sd[None, :] * _GX[:, None]) / a          # (G, S)
    d = z[None, :, :] - offsets[:, None, None]                  # (K, G, S)
    lam = np.maximum(0.0, 1.0 - np.abs(d))
    return np.tensordot(_GW, lam, axes=([0], [1]))              # (K, S)


def route_aprime_american_put(n: int, H: float, eta: float, rho: float,
                              zmax: float = np.inf, mref: int = 4,
                              pr=None) -> dict:
    """American put under Route A'. Returns a dict of diagnostics."""
    d = T / n
    sqd = np.sqrt(d)
    dH = d ** H
    s_ref = np.sqrt(XI0)
    a_X = s_ref * sqd / mref

    def variance(k: int, m_idx: np.ndarray) -> np.ndarray:
        if eta == 0.0:
            return np.full(m_idx.shape, XI0, dtype=float)
        zc = (2.0 * m_idx - k) * dH
        if np.isfinite(zmax):
            zc = np.clip(zc, -zmax, zmax)
        t = k * d
        return XI0 * np.exp(eta * np.sqrt(2.0 * H) * zc
                            - 0.5 * eta ** 2 * t ** (2 * H))

    # widest |offset| ever needed: 5 sd plus the correlated shift plus 1
    sig_max = np.sqrt(variance(n, np.array([n if np.isfinite(zmax) else n])).max()) \
        if eta > 0 else s_ref
    if eta > 0:
        zc_max = min(n * dH, zmax) if np.isfinite(zmax) else n * dH
        sig_max = np.sqrt(XI0 * np.exp(eta * np.sqrt(2.0 * H) * zc_max))
    kmax = int(np.ceil(4.5 * sig_max * np.sqrt(1 - rho ** 2) * sqd / a_X
                       + abs(rho) * sig_max * sqd / a_X
                       + 0.5 * sig_max ** 2 * d / a_X)) + 1
    offsets = np.arange(-kmax, kmax + 1)

    reach = kmax * n                       # max lattice displacement of X
    nx = 2 * reach + 1
    ix0 = reach
    logS = np.log(S0) + (np.arange(nx) - ix0) * a_X
    payoff_all = np.maximum(KSTRIKE - np.exp(logS), 0.0)        # (nx,)

    # terminal layer: (nx, n+1)
    val = np.repeat(payoff_all[:, None], n + 1, axis=1)

    neg = 0
    mass_err = 0.0
    for k in range(n - 1, -1, -1):
        mi = np.arange(k + 1, dtype=float)
        v_here = variance(k, mi)                                # (k+1,)
        sig = np.sqrt(v_here)
        mu_x = -0.5 * v_here
        # only the band reachable in k steps matters
        lo = max(0, ix0 - kmax * k)
        hi = min(nx, ix0 + kmax * k + 1)
        w = hi - lo
        cont = np.zeros((w, k + 1))
        for j, zeta in ((1, +1.0), (0, -1.0)):
            mu_M = sig * rho * sqd * zeta + mu_x * d
            sd_M = sig * np.sqrt(1.0 - rho ** 2) * sqd
            P = interp_kernel(mu_M, sd_M, a_X, offsets)         # (K, k+1)
            neg += int((P < -1e-12).sum())
            mass_err = max(mass_err, float(np.abs(P.sum(axis=0) - 1.0).max()))
            nxt = val[:, j:j + k + 1]                           # (nx, k+1)
            for oi, o in enumerate(offsets):
                src = np.clip(np.arange(lo, hi) + o, 0, nx - 1)
                cont += 0.5 * P[oi][None, :] * nxt[src, :]
        val = np.full((nx, k + 1), 0.0)
        val[lo:hi, :] = np.maximum(payoff_all[lo:hi, None], cont)
        val[:lo, :] = payoff_all[:lo, None]
        val[hi:, :] = payoff_all[hi:, None]
        if pr is not None and (k % max(1, n // 10) == 0):
            pr.tick(pr.step, inner=f"A' n={n} k={k}")
    return {"value": float(val[ix0, 0]),
            "negative_probabilities": neg,
            "max_mass_error": round(mass_err, 12),
            "a_X_over_sigma_sqrt_delta": round(float(a_X / (s_ref * sqd)), 4),
            "kmax": int(kmax), "nx": int(nx), "mref": int(mref)}


# =========================================================================== #
# European variants: the tree without early exercise, and a Monte-Carlo of the
# TRUE rough model.  Their gap measures the covariance discrepancy of §8, which
# Route A' does NOT fix (see Remark on the scope of (A'1)).
# =========================================================================== #
def route_aprime_european_put(n: int, H: float, eta: float, rho: float,
                              zmax: float = np.inf, mref: int = 4) -> float:
    """Same lattice, exercise only at T."""
    d = T / n
    sqd = np.sqrt(d)
    dH = d ** H
    s_ref = np.sqrt(XI0)
    a_X = s_ref * sqd / mref

    def variance(k, m_idx):
        if eta == 0.0:
            return np.full(np.shape(m_idx), XI0, dtype=float)
        zc = (2.0 * np.asarray(m_idx, dtype=float) - k) * dH
        if np.isfinite(zmax):
            zc = np.clip(zc, -zmax, zmax)
        return XI0 * np.exp(eta * np.sqrt(2.0 * H) * zc
                            - 0.5 * eta ** 2 * (k * d) ** (2 * H))

    zc_max = min(n * dH, zmax) if np.isfinite(zmax) else n * dH
    sig_max = s_ref if eta == 0 else np.sqrt(
        XI0 * np.exp(eta * np.sqrt(2.0 * H) * zc_max))
    kmax = int(np.ceil(4.5 * sig_max * np.sqrt(1 - rho ** 2) * sqd / a_X
                       + abs(rho) * sig_max * sqd / a_X)) + 1
    offsets = np.arange(-kmax, kmax + 1)
    reach = kmax * n
    nx = 2 * reach + 1
    ix0 = reach
    logS = np.log(S0) + (np.arange(nx) - ix0) * a_X
    payoff = np.maximum(KSTRIKE - np.exp(logS), 0.0)
    val = np.repeat(payoff[:, None], n + 1, axis=1)
    for k in range(n - 1, -1, -1):
        mi = np.arange(k + 1, dtype=float)
        v_here = variance(k, mi)
        sig = np.sqrt(v_here)
        lo, hi = max(0, ix0 - kmax * k), min(nx, ix0 + kmax * k + 1)
        cont = np.zeros((hi - lo, k + 1))
        for j, zeta in ((1, +1.0), (0, -1.0)):
            mu_M = sig * rho * sqd * zeta - 0.5 * v_here * d
            sd_M = sig * np.sqrt(1.0 - rho ** 2) * sqd
            P = interp_kernel(mu_M, sd_M, a_X, offsets)
            nxt = val[:, j:j + k + 1]
            for oi, o in enumerate(offsets):
                src = np.clip(np.arange(lo, hi) + o, 0, nx - 1)
                cont += 0.5 * P[oi][None, :] * nxt[src, :]
        nv = np.zeros((nx, k + 1))
        nv[lo:hi, :] = cont
        nv[:lo, :] = val[:lo, :k + 1]
        nv[hi:, :] = val[hi:, :k + 1]
        val = nv
    return float(val[ix0, 0])


def european_put_mc_rough(H: float, eta: float, rho: float,
                          nfine: int = 2048, paths: int = 200_000,
                          chunk: int = 20_000, seed: int = 7) -> tuple[float, float]:
    """Monte-Carlo of the TRUE rough Bergomi European put.

    V^H_t = sqrt(2H) * discrete Ito sum with the exact fractional kernel, so the
    covariance structure is the real one.  Returns (price, standard error).
    """
    rng = np.random.default_rng(seed)
    df = T / nfine
    h = H - 0.5
    ker = np.sqrt(2.0 * H) * (np.arange(1, nfine + 1) * df) ** h    # (nfine,)
    L = 1
    while L < 2 * nfine:
        L *= 2
    fker = np.fft.rfft(ker, L)
    tgrid = np.arange(1, nfine + 1) * df
    tot = 0.0
    tot2 = 0.0
    done = 0
    while done < paths:
        p = min(chunk, paths - done)
        dW = rng.standard_normal((p, nfine)) * np.sqrt(df)
        dWp = rng.standard_normal((p, nfine)) * np.sqrt(df)
        # VH[k] = sum_{i<=k} ker[k-i] dW[i]  (causal convolution)
        VH = np.fft.irfft(np.fft.rfft(dW, L, axis=1) * fker[None, :],
                          L, axis=1)[:, :nfine]
        v = XI0 * np.exp(eta * VH - 0.5 * eta ** 2 * tgrid[None, :] ** (2 * H))
        # v must be evaluated at the LEFT end of each interval: VH[k] already
        # contains dW[k], so pairing v[k] with dB[k] would be anticipating and
        # (since dB is correlated with dW) would bias the price upwards.
        v = np.concatenate([np.full((p, 1), XI0), v[:, :-1]], axis=1)
        dB = rho * dW + np.sqrt(1.0 - rho ** 2) * dWp
        logST = (np.log(S0) - 0.5 * (v * df).sum(axis=1)
                 + (np.sqrt(v) * dB).sum(axis=1))
        pay = np.maximum(KSTRIKE - np.exp(logST), 0.0)
        tot += pay.sum()
        tot2 += (pay ** 2).sum()
        done += p
    mean = tot / paths
    var = max(0.0, tot2 / paths - mean ** 2)
    return float(mean), float(np.sqrt(var / paths))
