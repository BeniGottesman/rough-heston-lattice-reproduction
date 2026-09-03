"""
frozen_exit -- Route F (research-notes/L003-ROUTE-F.md), falsifiers FF1-FF4.

What is being separated
-----------------------
Two operations are both called "freezing" in the manuscript and the note
L003-ROUTE-F.md states the split exactly.  With K(u) = u^h, h in (-1/2, 0),
H := h + 1/2, T = 1, delta = T/n, t_k = k delta:

  * absorbed-driver convolution   Vcal(u)  = v0 + int_0^{u ^ Xi^y} K(u-s) dy_s
    -- the integrand stops, the KERNEL ARGUMENT KEEPS RUNNING;
  * frozen convolution            V_t      = v_{t ^ Xi^y} = Vcal(t ^ Xi^y)
    -- both stop.

and on the discrete side, with the CLAMPED walk of eq:walks and its exit index
Xi^Y_n of eq:XiN,

  * Vcal^(n)(t_k) = v0 + sum_{i<=k} K(t_k - t_{i-1}) (Yn_i - Yn_{i-1})   unfrozen
  * V^(n)_k       = Vcal^(n)(t_{k ^ Xi^Y_n})                             frozen (dfn:Vexact)

Note that for t <= Xi^y one has Vcal(t) = v_t = v0 + int_0^t K(t-s) dy_s with
the RAW driver, so V_t = v_{t ^ Xi^y} can be read off the unabsorbed continuous
convolution.  Same on the discrete side for k <= Xi^Y_n: the clamp in eq:walks
is inactive before the exit index.  Both frozen objects are therefore computable
from the two UNFROZEN convolutions plus the two exit times, which is what this
module does.

The four measurements
--------------------
FF1  ||max_k |Vcal^(n)(t_k) - V^(n)_k| ||_{L2}  -- purely discrete, no coupling
     needed, so it can be pushed to n = 10^6.  Predicted by the note to sit at
     a positive constant instead of decaying like delta^H.

FF2  ||max_k |V^(n)_k - V_{theta_k}| ||_{L2}   -- frozen against frozen, on a
     COMMON Brownian path.  Coupling: the exact Skorokhod (level-crossing)
     embedding, theta_k = the successive hitting times of the sqrt(delta)
     lattice.  For mu_y = 0 this is the manuscript's own embedding and it is
     EXACT: Y(theta_k) = Yn_k and the node error e_j vanishes (part3, l.582).
     It is realised on a fine grid, so the only defects are (a) crossing
     overshoot O(sqrt(dt_fine)) and (b) the left-point quadrature error
     O(dt_fine^H) of the continuous convolution.  Both are reported.

FF3  P( theta^(n)_{Xi^Y_n} > Xi^y ), the failure of the ordering lemma L1, under
     three constructions of increasing node error (see run_ff3).

FF4  G := Xi^y - theta^(n)_{Xi^Y_n} >= 0.  For constant coefficients the walk
     sits EXACTLY on the discrete barrier at theta_{Xi^Y_n}, so by the strong
     Markov property G is the exit time of a Brownian motion from (B^Y, C^Y)
     started at the discrete barrier.  That law is known in closed form, so the
     primary FF4 numbers are computed from the spectral / reflection series
     rather than by Monte Carlo, and a coupled MC run cross-checks the
     reduction.

Conventions used everywhere
---------------------------
T = 1, y_0 = 0, v_0 = 0 (it cancels), K(u) = u^h (L == 1),
continuous band (B^Y, C^Y) = (-BAND, +BAND).

NOTE on a manuscript inconsistency: part1-setup.tex declares B^Y in [0, inf)
while y_0 = 0, so y_0 is never strictly inside the Y-band as written.  The
[0, inf) constraint is a leftover from the drafts in which the Y-barriers were
centred at v_0 (a variance, hence positive); the \fixed{} note in eq:barriers
recentres them at y_0.  We take the recentred reading and use a symmetric band
about y_0 = 0.
"""
from __future__ import annotations

import math

import numpy as np
from numpy.fft import irfft, rfft

# ------------------------------------------------------------------ constants
T = 1.0
Y0 = 0.0
V0 = 0.0
BAND = 0.5                     # (B^Y, C^Y) = (-BAND, +BAND)
BY, CY = -BAND, BAND


# ------------------------------------------------------------------- barriers
BETA0 = 1.0 / 3.0        # the inset exponent as eq:barriers writes it today


def discrete_barriers(delta: float, band: float = BAND, y0: float = Y0,
                      beta: float = BETA0):
    """B^Y_n, C^Y_n of eq:barriers, closed form of part2-scheme.tex.

    `beta` generalises the barrier INSET: eq:barriers writes the literal
    `delta^{1/3}`, and the value 1/3 occurs in exactly four places
    (part2-scheme.tex:24-30, 44-47, 57, 71-81), all inside the definition of the
    discrete barriers and inside `lem:barriers-wellposed` itself.  Nothing
    downstream reads the value, so `delta^{beta}` with `beta in (0, 1/2]` is a
    strict generalisation recovering the current text at `beta = 1/3`.  At
    `beta = 1/2` the inset is the grid resolution `sqrt(delta)` itself, the
    smallest a `sqrt(delta)`-lattice can express.
    """
    sq = math.sqrt(delta)
    eps = delta ** beta
    Bn = y0 + sq * (1 + math.floor(((-band) + eps - y0) / sq))
    Cn = y0 + sq * (math.ceil((band - eps - y0) / sq) - 1)
    return Bn, Cn


def barrier_gap(delta: float, band: float = BAND, y0: float = Y0,
                beta: float = BETA0) -> float:
    """eps_n as realised: C^Y - C^Y_n, in (delta^beta, delta^beta + sqrt delta]."""
    _, Cn = discrete_barriers(delta, band, y0, beta)
    return band - Cn


def band_wellformed(delta: float, band: float = BAND, y0: float = Y0,
                    beta: float = BETA0) -> dict:
    """Is the discrete band well formed at this (delta, beta)?

    `lem:barriers-wellposed` asserts, for n large enough, B^Y_n < C^Y_n, both
    strictly inside (B^Y, C^Y), and both on the sqrt(delta) grid.  At beta = 1/2
    the inset is the same order as the grid spacing, so the lemma's "n large
    enough" is a real constraint and has to be checked, not assumed.
    """
    sq = math.sqrt(delta)
    Bn, Cn = discrete_barriers(delta, band, y0, beta)
    return {
        "Bn": Bn, "Cn": Cn, "gap": band - Cn, "sqrt_delta": sq,
        "width_in_steps": (Cn - Bn) / sq,
        "ordered": Bn < Cn,
        "inside": (-band < Bn) and (Cn < band),
        "y0_strictly_inside": (Bn < y0 < Cn),
        "on_grid": abs((Cn - y0) / sq - round((Cn - y0) / sq)) < 1e-9,
        "ok": (Bn < Cn) and (-band < Bn) and (Cn < band) and (Bn < y0 < Cn),
    }


# --------------------------------------------------------------- convolutions
def _pow2(m: int) -> int:
    p = 1
    while p < m:
        p <<= 1
    return p


def causal_conv(dz: np.ndarray, kern: np.ndarray) -> np.ndarray:
    """out[:, k-1] = sum_{i=1..k} kern[k-i] * dz[:, i-1],  k = 1..N.

    With kern[j] = ((j+1) dt)^h this is exactly
    sum_{i<=k} K(t_k - t_{i-1}) dz_i, the left-point convolution used by both
    dfn:Vexact and (in the fine-grid limit) the continuous v_t.
    """
    n = dz.shape[1]
    L = _pow2(2 * n)
    F = rfft(dz, L, axis=1) * rfft(kern, L)
    return irfft(F, L, axis=1)[:, :n]


def kernel_array(n: int, dt: float, h: float) -> np.ndarray:
    return (np.arange(1, n + 1, dtype=np.float64) * dt) ** h


# ---------------------------------------------------------------- FF1: purely
#                                                                   discrete
def ff1_batch(n: int, h: float, paths: int, rng):
    """The FF1 statistic and its two-part decomposition, for `paths` walks.

    Returns per-path
        total   = max_k |Vcal^(n)(t_k) - V^(n)_k|
        kernpt  = max_{k>Xi} |sum_{i<=Xi} [K(t_k-t_{i-1}) - K(t_Xi-t_{i-1})] dY_i|
                  -- the piece rem:freezing-convention tries to telescope,
        incrpt  = max_{k>Xi} |sum_{Xi<i<=k} K(t_k-t_{i-1}) dY_i|
                  -- the piece that the sentence "past Xi^Y_n the increments
                  vanish" claims is zero,
        scale   = max_k |Vcal^(n)(t_k)|          (for scale),
    plus the exit fraction and the fraction of post-exit increments that are
    NON-zero (the direct test of that sentence's premise).
    """
    delta = T / n
    sq = math.sqrt(delta)
    Bn, Cn = discrete_barriers(delta)

    zeta = rng.integers(0, 2, size=(paths, n)) * 2.0 - 1.0
    raw = Y0 + sq * np.cumsum(zeta, axis=1)
    Y = np.clip(raw, Bn, Cn)                       # eq:walks CLAMPS
    out = (Y <= Bn) | (Y >= Cn)                    # exit from the OPEN band
    has = out.any(axis=1)
    Xi = np.where(has, out.argmax(axis=1) + 1, n)  # 1-based index k

    dY = np.empty_like(Y)
    dY[:, 0] = Y[:, 0] - Y0
    dY[:, 1:] = np.diff(Y, axis=1)

    kk = np.arange(1, n + 1)[None, :]
    post = kk > Xi[:, None]
    nz_post = float((np.abs(dY) > 0)[post].mean()) if post.any() else 0.0

    kern = kernel_array(n, delta, h)
    Vcal = V0 + causal_conv(dY, kern)                      # Vcal[:, k-1]
    dYpre = np.where(post, 0.0, dY)
    Vpre = V0 + causal_conv(dYpre, kern)

    idx = (np.arange(paths), Xi - 1)
    at_xi = Vcal[idx]
    total = np.abs(Vcal - at_xi[:, None])
    kernpt = np.abs(Vpre - Vpre[idx][:, None])
    incrpt = np.abs(Vcal - Vpre)
    for A in (total, kernpt, incrpt):
        A[~post] = 0.0
    return (total.max(axis=1), kernpt.max(axis=1), incrpt.max(axis=1),
            np.abs(Vcal).max(axis=1), float(has.mean()), nz_post)


# ------------------------------------------------ Skorokhod level embedding
def embed_walk(yf: np.ndarray, a: float, y0: float = Y0):
    """Exact Skorokhod embedding of the +-a lattice walk into a continuous path.

    theta_k = inf{ t > theta_{k-1} : |y_t - y_{theta_{k-1}}| = a }.  Because the
    path is continuous and starts on the lattice y0 + a Z, theta_k is simply the
    k-th time the path reaches a NEW adjacent lattice level.  Writing
    c_j := floor((y_j - y0)/a), the level crossed between fine indices j and j+1
    is max(c_j, c_{j+1}); the embedded walk is that sequence with consecutive
    duplicates removed (the duplicates are re-crossings of the level the walk is
    already sitting on).

    Returns (fine_index, level) for k = 1, 2, ... and the count of fine steps
    that jumped more than one lattice level (which would break the argument).
    """
    c = np.floor((yf - y0) / a).astype(np.int64)
    d = np.diff(c)
    ev = np.flatnonzero(d)
    if ev.size == 0:
        return np.empty(0, np.int64), np.empty(0, np.int64), 0
    multi = int(np.count_nonzero(np.abs(d[ev]) > 1))
    lev = np.maximum(c[ev], c[ev + 1])
    keep = np.empty(lev.size, dtype=bool)
    keep[0] = lev[0] != 0
    keep[1:] = lev[1:] != lev[:-1]
    return ev[keep] + 1, lev[keep], multi


def first_exit_index(yf: np.ndarray, band: float = BAND) -> int:
    """First fine index at which |y| >= band; len(yf) if it never happens."""
    o = np.abs(yf) >= band
    j = int(o.argmax()) if o.any() else len(yf)
    return j


# ------------------------------------------------------------------- FF2 core
def ff2_one_path(yf: np.ndarray, vfull: dict, mpu: int, n_list, h_list,
                 band: float = BAND, betas=(BETA0,), vabs: dict | None = None):
    """One fine Brownian path -> the FF2/FF6 statistics for every (n, h, beta).

    Route F splits (L003-ROUTE-F.md sec.2), with nu := k ^ Xi^Y_n:

        V^(n)_k - V_{theta_k}
          = [ Vcal^(n)(t_nu) - Vcal(theta_nu) ]                (F1)
          + [ Vcal(theta_nu) - Vcal(theta_k ^ Xi^y) ]          (F2)

    The EMBEDDING, the walk and the discrete convolution do not depend on beta;
    only the exit index Xi^Y_n does, through the discrete barriers.  So all
    betas are evaluated on exactly the same path and the same walk, which is
    what makes the side-by-side comparison paired.

    Returns, per (n, h, key, beta):
        total_theta = max_k |V^(n)_k - V_{theta_k}|      (prop:Vconv's object)
        total_tk    = max_k |V^(n)_k - V_{t_k}|
        F1          = max_{j<=Xi^Y_n} |Vcal^(n)(t_j) - v_{theta_j}|
        F2          = max_{k>Xi^Y_n} |v_{theta_k ^ Xi^y} - v_{theta_{Xi}}|
    and per (n, beta) the exit-time gap G.
    """
    M = yf.size - 1
    jXi_y = first_exit_index(yf, band)
    res, diag = {}, {}
    for n in n_list:
        delta = T / n
        a = math.sqrt(delta)
        idx, lev, multi = embed_walk(yf, a)
        avail = idx.size
        kmax = min(n, avail)
        idx, lev = idx[:kmax], lev[:kmax]
        Yn = Y0 + a * lev                                   # exact lattice walk
        e_max = float(np.max(np.abs(Yn - yf[idx]))) if kmax else 0.0
        dY = np.empty(kmax)
        if kmax:
            dY[0] = Yn[0] - Y0
            dY[1:] = np.diff(Yn)
        theta_idx = np.concatenate(([0], idx))              # theta_0 = 0
        rd_theta = np.minimum(theta_idx, jXi_y)
        step = mpu // n                                     # t_k on the fine grid
        rd_tk = np.minimum(np.minimum(np.arange(kmax + 1) * step, M), jXi_y)
        Vc = {}
        for h in h_list:
            if kmax:
                c = causal_conv(dY[None, :], kernel_array(kmax, delta, h))[0]
                Vc[h] = np.concatenate(([V0], c))           # index k = 0..kmax
            else:
                Vc[h] = np.array([V0])
        for beta in betas:
            Bn, Cn = discrete_barriers(delta, band, Y0, beta)
            Vcl = {}
            if vabs is not None and kmax:
                Ycl = np.clip(Yn, Bn, Cn)          # eq:walks CLAMPS, not absorbs
                dYc = np.empty(kmax)
                dYc[0] = Ycl[0] - Y0
                dYc[1:] = np.diff(Ycl)
                for h in h_list:
                    c = causal_conv(dYc[None, :], kernel_array(kmax, delta, h))[0]
                    Vcl[h] = np.concatenate(([V0], c))
            elif vabs is not None:
                Vcl = {h: np.array([V0]) for h in h_list}
            outw = (Yn <= Bn) | (Yn >= Cn)
            Xi = (int(outw.argmax()) + 1) if outw.any() else kmax
            kfr = np.minimum(np.arange(kmax + 1), Xi)
            G = (min(jXi_y, M) - theta_idx[Xi]) * (1.0 / mpu)
            for h in h_list:
                vf = vfull[h]
                Vn = Vc[h][kfr]
                res[(n, h, "total_theta", beta)] = float(np.max(np.abs(Vn - vf[rd_theta])))
                res[(n, h, "total_tk", beta)] = float(np.max(np.abs(Vn - vf[rd_tk])))
                # (F1) of the route: nu = k ^ Xi^Y_n, so both arguments are in
                # the PRE-absorption regime (theta_nu <= Xi^y by L1) and Vcal = v
                # there.  Correct as written; its index RANGE depends on beta.
                res[(n, h, "F1", beta)] = float(
                    np.max(np.abs(Vc[h][:Xi + 1] - vf[rd_theta[:Xi + 1]])))
                # F1_full: prop:Vconv's OWN object, the quantity the route says
                # dominates (F1) -- max over ALL indices 0..n of the UNFROZEN
                # discrete convolution against the ABSORBED-DRIVER continuous
                # one.  Both sides keep their kernel argument running past the
                # exit; only the integrand stops.  Note the discrete side must
                # use the CLAMPED walk of eq:walks (which is why this term is
                # still beta-dependent -- through the clamp, not through a
                # truncated index range).
                if vabs is not None:
                    res[(n, h, "F1_full", beta)] = float(
                        np.max(np.abs(Vcl[h] - vabs[h][theta_idx])))
                tail = rd_theta[Xi:]
                res[(n, h, "F2", beta)] = float(
                    np.max(np.abs(vf[tail] - vf[rd_theta[Xi]])))
            diag[(n, beta)] = dict(avail=avail, kmax=kmax, Xi=Xi, e_max=e_max,
                                   multi=multi, exited=bool(outw.any()),
                                   y_exited=bool(jXi_y <= M), G=G)
    return res, diag


# -------------------------------------------------------- FF4: the exact law
def _psurv_reflect(t: np.ndarray, x: float, L: float, mmax: int = 8):
    """P(tau > t) for BM in (0, L) from x, reflection series (fast for small t)."""
    from math import erf
    s = np.sqrt(np.maximum(t, 1e-300))
    out = np.zeros_like(t)
    Phi = lambda z: 0.5 * (1.0 + np.vectorize(erf)(z / math.sqrt(2.0)))
    for m in range(-mmax, mmax + 1):
        out += (Phi((L - x - 2 * m * L) / s) - Phi((-x - 2 * m * L) / s)
                - Phi((L + x - 2 * m * L) / s) + Phi((x - 2 * m * L) / s))
    return np.clip(out, 0.0, 1.0)


def _psurv_spectral(t: np.ndarray, x: float, L: float, kmax: int = 4000):
    """P(tau > t) eigenfunction series (fast for large t)."""
    k = np.arange(1, kmax + 1, 2, dtype=np.float64)
    coef = (4.0 / (k * math.pi)) * np.sin(k * math.pi * x / L)
    lam = (k * math.pi / L) ** 2 / 2.0
    return np.clip((coef[None, :] * np.exp(-lam[None, :] * t[:, None])).sum(axis=1),
                   0.0, 1.0)


def psurv(t: np.ndarray, x: float, L: float) -> np.ndarray:
    """P(G > t) for BM in (0, L) started at x, both series glued."""
    t = np.asarray(t, dtype=np.float64)
    out = np.empty_like(t)
    cut = 0.15 * L * L
    lo, hi = t < cut, t >= cut
    if lo.any():
        out[lo] = _psurv_reflect(t[lo], x, L)
    if hi.any():
        out[hi] = _psurv_spectral(t[hi], x, L)
    return out


def moment_G(kap: float, x: float, L: float, tmax: float | None = None,
             lo: float = 1e-14, npts: int = 6000) -> float:
    """E[(G ^ tmax)^kap] = kap int_0^tmax t^{kap-1} P(G>t) dt + tmax^kap P(G>tmax).

    tmax = None gives the uncensored moment (the tail is exponential in a bounded
    band, so the integral converges).  Log-grid trapezoid on u = log t, which is
    exact to machine noise for this smooth integrand; validated at kap = 1 where
    E[G] = x(L-x) exactly.
    """
    hi = tmax if tmax is not None else 50.0 * L * L
    u = np.linspace(math.log(lo), math.log(hi), npts)
    t = np.exp(u)
    f = kap * (t ** kap) * psurv(t, x, L)          # kap t^{kap-1} P dt = kap t^kap P du
    val = float(np.trapezoid(f, u)) if hasattr(np, "trapezoid") else float(np.trapz(f, u))
    if tmax is not None:
        val += (tmax ** kap) * float(psurv(np.array([tmax]), x, L)[0])
    return val


def ff4_exact(n: int, kaps, band: float = BAND, tmax: float | None = None,
              beta: float = BETA0):
    """FF4(a) primary numbers, from the strong-Markov reduction.

    At theta_{Xi^Y_n} the exactly-embedded walk sits ON the discrete barrier
    C^Y_n (or B^Y_n; symmetric).  By strong Markov, G is the exit time of a BM
    from (B^Y, C^Y) started there: a band of width L = 2 band, entered at
    distance g_n = C^Y - C^Y_n from the near edge.
    """
    delta = T / n
    Bn, Cn = discrete_barriers(delta, band, Y0, beta)
    L = 2.0 * band
    x = Cn - (-band)                      # position inside (0, L)
    g = band - Cn
    return {"n": n, "delta": delta, "gap": g, "beta": beta,
            "eps_n": delta ** beta + math.sqrt(delta),
            "x": x, "L": L,
            **{f"EG^{k}": moment_G(k, x, L, tmax) for k in kaps},
            "EG": moment_G(1.0, x, L, tmax),
            "EG_exact_uncensored": x * (L - x)}


# ------------------------------------------------------------------- fitting
def loglog_fit(xs, ys):
    """slope of log y against log x (xs = delta values)."""
    lx, ly = np.log(np.asarray(xs, float)), np.log(np.asarray(ys, float))
    A = np.vstack([lx, np.ones_like(lx)]).T
    coef, *_ = np.linalg.lstsq(A, ly, rcond=None)
    return float(coef[0]), float(coef[1])


def bootstrap_slope(deltas, per_path, B: int = 400, seed: int = 11):
    """per_path: array (len(deltas), P) of per-path maxima, PATHS SHARED across
    deltas.  Resamples whole paths, recomputes the L2 norms, refits the slope."""
    rng = np.random.default_rng(seed)
    P = per_path.shape[1]
    sl = np.empty(B)
    for b in range(B):
        j = rng.integers(0, P, P)
        l2 = np.sqrt((per_path[:, j] ** 2).mean(axis=1))
        sl[b], _ = loglog_fit(deltas, l2)
    return float(sl.mean()), float(sl.std(ddof=1)), (float(np.percentile(sl, 2.5)),
                                                     float(np.percentile(sl, 97.5)))


def bootstrap_l2(x: np.ndarray, B: int = 400, seed: int = 7):
    rng = np.random.default_rng(seed)
    P = x.size
    v = np.empty(B)
    for b in range(B):
        v[b] = math.sqrt((x[rng.integers(0, P, P)] ** 2).mean())
    return math.sqrt((x ** 2).mean()), float(v.std(ddof=1))


# ------------------------------------------------------------- FF3 machinery
SIG3_A, SIG3_B = 1.0, 0.3        # sigma_y(y) = 1 + 0.3 tanh(y) in [0.7, 1.3]
MU3 = 0.3                        # mu_y(y)    = 0.3 cos(y), bounded + Lipschitz


def sigma3(y):
    return SIG3_A + SIG3_B * np.tanh(y)


def mu3(y):
    return MU3 * np.cos(y)


def euler_paths(M: int, dt: float, paths: int, rng, drift: bool = True,
                aux: bool = False):
    """Euler-Maruyama for dy = mu_y(y) dt + sigma_y(y) dW on a fine grid.

    Stored as (M+1, paths) so the time loop writes contiguous rows.  Optionally
    returns an independent auxiliary Brownian path used to shape the injected
    node error of FF3 variant C.
    """
    sq = math.sqrt(dt)
    y = np.zeros((M + 1, paths), dtype=np.float64)
    cur = np.zeros(paths)
    b = np.zeros((M + 1, paths), dtype=np.float64) if aux else None
    curb = np.zeros(paths)
    for j in range(M):
        dW = rng.standard_normal(paths) * sq
        cur = cur + sigma3(cur) * dW + (mu3(cur) * dt if drift else 0.0)
        y[j + 1] = cur
        if aux:
            curb = curb + rng.standard_normal(paths) * sq
            b[j + 1] = curb
    return (y, b) if aux else y


def ff3_variantA(yf: np.ndarray, n: int, band: float = BAND):
    """Genuine Skorokhod embedding (level crossings) into y; mu_y == 0 so the
    embedded increments are fair +-1 coins and the walk IS the eq:walks walk."""
    delta = T / n
    a = math.sqrt(delta)
    Bn, Cn = discrete_barriers(delta, band)
    idx, lev, _ = embed_walk(yf, a)
    if idx.size == 0:
        return None
    lev = lev[:n]
    idx = idx[:n]
    Yn = Y0 + a * lev
    outw = (Yn <= Bn) | (Yn >= Cn)
    if not outw.any():
        return None
    k = int(outw.argmax())
    jXi_y = first_exit_index(yf, band)
    e = float(np.max(np.abs(Yn[:k + 1] - yf[idx[:k + 1]])))
    return idx[k], jXi_y, e


def ff3_grid_surrogate(yf: np.ndarray, n: int, mpu: int, band: float = BAND,
                       inject: np.ndarray | None = None):
    """Surrogate: theta_k := t_k and Y^(n)_k := sqrt(delta)-lattice rounding of
    y(t_k) (+ an injected node error).  NOT the eq:walks walk -- its increments
    are multiples of sqrt(delta), not +-sqrt(delta).  Used only to interrogate
    the MECHANISM of L1: can the node error push y(theta_{Xi^Y_n}) outside the
    continuous band before the walk leaves the narrow one?
    """
    delta = T / n
    a = math.sqrt(delta)
    Bn, Cn = discrete_barriers(delta, band)
    step = mpu // n
    ks = np.arange(n + 1) * step
    ys = yf[ks]
    tgt = ys if inject is None else ys + inject
    Yn = Y0 + a * np.round((tgt - Y0) / a)
    outw = (Yn <= Bn) | (Yn >= Cn)
    if not outw.any():
        return None
    k = int(outw.argmax())
    jXi_y = first_exit_index(yf, band)
    e = float(np.max(np.abs(Yn[:k + 1] - ys[:k + 1])))
    return ks[k], jXi_y, e


def lp_norm(x: np.ndarray, ell: float) -> float:
    return float(np.mean(np.abs(x) ** ell) ** (1.0 / ell))


def bootstrap_lp(x: np.ndarray, ell: float, B: int = 400, seed: int = 7):
    rng = np.random.default_rng(seed)
    P = x.size
    v = np.array([lp_norm(x[rng.integers(0, P, P)], ell) for _ in range(B)])
    return lp_norm(x, ell), float(v.std(ddof=1))


def bootstrap_slope_lp(deltas, per_path: np.ndarray, ell: float,
                       B: int = 400, seed: int = 11):
    """per_path: (len(deltas), P), PATHS SHARED across deltas."""
    rng = np.random.default_rng(seed)
    P = per_path.shape[1]
    sl = np.empty(B)
    for b in range(B):
        j = rng.integers(0, P, P)
        nm = [lp_norm(per_path[i, j], ell) for i in range(per_path.shape[0])]
        sl[b], _ = loglog_fit(deltas, nm)
    return float(sl.mean()), float(sl.std(ddof=1)), (float(np.percentile(sl, 2.5)),
                                                     float(np.percentile(sl, 97.5)))


def ff7_exact(n: int, c: float, band: float = BAND, beta: float = BETA0) -> dict:
    """FF7: P(G > c) at a FIXED c, against eps_n.

    This is the lower-bound mechanism of the sharpness claim measured directly:
    with probability of order eps_n the diffusion started at the discrete barrier
    stays inside the band for a time of order 1, and on that event the frozen and
    unfrozen objects differ by Theta(1).  If P(G > c) does NOT scale like eps_n,
    the sharpness claim is wrong.
    """
    delta = T / n
    _, Cn = discrete_barriers(delta, band, Y0, beta)
    return {"n": n, "delta": delta, "beta": beta, "gap": band - Cn,
            "P": float(psurv(np.array([c]), Cn + band, 2.0 * band)[0])}
