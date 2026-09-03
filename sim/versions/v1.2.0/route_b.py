"""
route_b — the Markovian lift of paper section "Route B", quantified.

The Riemann--Liouville kernel with h = H - 1/2 in (-1/2, 0) has the completely
monotone representation

    K(t) = t^h / Gamma(1+h) = int_0^infty exp(-s t) mu_h(ds),
    mu_h(ds) = c_h s^{-h-1} ds,      c_h = 1 / (Gamma(-h) Gamma(1+h)).

An m-factor lift replaces mu_h by a finite measure sum_i w_i delta_{s_i}, giving
an (m+1)-dimensional MARKOV state and a kernel

    K^m(t) = sum_{i=1}^m w_i exp(-s_i t).

Everything here is analytic: the L^2 distance ||K - K^m||_{L^2(0,T)} and the
covariance surface of the lifted Gaussian component have closed forms, so no
Monte-Carlo and no lattice is needed to answer the question that decides whether
Route B is usable at all --- how large must m be.

Two node/weight constructions are provided.

  ajee      Abi Jaber--El Euch: partition (0, eta_m] into m cells, put on each
            cell the exact mu_h mass and the mu_h-barycentre of s.  Mass and
            first moment of mu_h are matched cell by cell, and every weight is
            positive by construction, which is what the lift needs (each w_i is
            the variance loading of a genuine OU factor).  The partition is a
            two-parameter family, optimised numerically for each m.

  best      the achievable floor: nodes free, weights the non-negative
            least-squares solution.  This is not a construction one would
            implement --- it is an upper bound on how good ANY m-factor lift
            with positive weights can be, so that a negative verdict on m
            cannot be blamed on a poor choice of nodes.

Conventions.  K carries the 1/Gamma(1+h) normalisation of the Laplace
representation above.  The bare kernel u^h used in the handoff differs by the
constant factor Gamma(1+h), which matters when an L^2 norm is compared with a
dimensionless rate; both are reported and the caller must say which it means.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize, nnls
from scipy.special import gamma as Gamma
from scipy.special import gammainc


# --------------------------------------------------------------------- kernel
def c_h(h: float) -> float:
    """Total-mass constant of mu_h.  Reflection: Gamma(-h)Gamma(1+h)=pi/sin(-pi h)."""
    return 1.0 / (Gamma(-h) * Gamma(1.0 + h))


def K_true(t, h: float):
    return np.power(t, h) / Gamma(1.0 + h)


def K_lift(t, w: np.ndarray, s: np.ndarray):
    t = np.atleast_1d(np.asarray(t, float))
    return (w[None, :] * np.exp(-np.outer(t, s))).sum(axis=1)


def norm_K_sq(h: float, T: float, t0: float = 0.0) -> float:
    """||K||^2_{L^2(t0,T)} = (T^{2h+1} - t0^{2h+1}) / ((2h+1) Gamma(1+h)^2)."""
    a = 2.0 * h + 1.0
    return (T ** a - t0 ** a) / (a * Gamma(1.0 + h) ** 2)


def _b_vector(s: np.ndarray, h: float, T: float, t0: float = 0.0) -> np.ndarray:
    """b_i = int_{t0}^T K(t) exp(-s_i t) dt, via the lower incomplete gamma."""
    a = h + 1.0
    # scipy's gammainc is the REGULARISED lower incomplete gamma
    low = (gammainc(a, s * T) - gammainc(a, s * t0)) * Gamma(a)
    return s ** (-a) * low / Gamma(1.0 + h)


def _gram(s: np.ndarray, T: float, t0: float = 0.0) -> np.ndarray:
    """G_ij = int_{t0}^T exp(-(s_i+s_j)t) dt."""
    ss = s[:, None] + s[None, :]
    return (np.exp(-ss * t0) - np.exp(-ss * T)) / ss


def l2_error(w: np.ndarray, s: np.ndarray, h: float, T: float,
             t0: float = 0.0) -> float:
    """||K - K^m||_{L^2(t0,T)}, exact.  Guarded against cancellation noise.

    t0 = 0 is the norm that (B1) of the paper is stated in.  t0 = delta is the
    norm the DISCRETE scheme actually sees, since a lattice with time step delta
    never evaluates the kernel closer to the origin than delta; the singularity
    of K at 0, which no finite sum of exponentials can follow, is then excluded.
    The two give very different answers, so the caller must say which it means.
    """
    A = norm_K_sq(h, T, t0)
    b = _b_vector(s, h, T, t0)
    G = _gram(s, T, t0)
    e2 = A - 2.0 * w @ b + w @ G @ w
    return math.sqrt(max(e2, 0.0))


_GL_X, _GL_W = np.polynomial.legendre.leggauss(20)


def l2_error_quadrature(w: np.ndarray, s: np.ndarray, h: float, T: float,
                        t0: float = 0.0, panels: int = 800) -> float:
    """Same quantity by composite Gauss--Legendre — an independent cross-check.

    Two features must be resolved or the check is worse than useless: the
    t^{2h} singularity of K at the origin, and one boundary layer of width 1/s_i
    per factor, which the optimal partitions push down to 10^-15.  Adaptive
    quadrature steps straight over those layers and returns a number that looks
    plausible and is wrong, so it is not used.

    For t0 = 0 the substitution t = z^p with p = 1/(1+2h) removes the
    singularity exactly: t^{2h} dt = p z^{p(2h+1)-1} dz = p dz.  Panels are then
    geometric in z, which resolves every layer, and the residual piece below the
    smallest panel is bounded by p z_min / Gamma(1+h)^2 and is negligible.
    """
    def integrate(f, a: float, b: float, npan: int) -> float:
        edges = np.geomspace(a, b, npan + 1)
        tot = 0.0
        for lo, hi in zip(edges[:-1], edges[1:]):
            x = 0.5 * (hi - lo) * _GL_X + 0.5 * (hi + lo)
            tot += 0.5 * (hi - lo) * (_GL_W * f(x)).sum()
        return tot

    if t0 > 0.0:
        return math.sqrt(max(integrate(
            lambda t: (K_true(t, h) - K_lift(t, w, s)) ** 2, t0, T, panels), 0.0))

    p = 1.0 / (1.0 + 2.0 * h)
    zmax = T ** (1.0 / p)
    zmin = zmax * 1e-14

    def g(z):
        t = z ** p
        return (K_true(t, h) - K_lift(t, w, s)) ** 2 * p * z ** (p - 1.0)

    tail = p * zmin / Gamma(1.0 + h) ** 2           # bound on (0, zmin)
    return math.sqrt(max(integrate(g, zmin, zmax, panels) + tail, 0.0))


def singularity_share(w: np.ndarray, s: np.ndarray, h: float, T: float,
                      t0: float) -> float:
    """Fraction of ||K-K^m||^2_{L^2(0,T)} contributed by (0, t0).

    If this is close to 1 the lift is being judged almost entirely on a
    neighbourhood of the origin that the discrete scheme never visits.
    """
    tot = l2_error(w, s, h, T, 0.0) ** 2
    if tot <= 0.0:
        return 0.0
    tail = l2_error(w, s, h, T, t0) ** 2
    return max(0.0, min(1.0, (tot - tail) / tot))


# ------------------------------------------------------- AJEE nodes & weights
def ajee_from_partition(eta: np.ndarray, h: float) -> tuple[np.ndarray, np.ndarray]:
    """Cell masses and mu_h-barycentres for the partition 0 = eta_0 < ... < eta_m.

        w_i = c_h (eta_i^{-h} - eta_{i-1}^{-h}) / (-h)
        s_i = (-h)/(1-h) * (eta_i^{1-h} - eta_{i-1}^{1-h})
                          / (eta_i^{-h}   - eta_{i-1}^{-h})

    Both are finite at eta_0 = 0 because -h > 0.  Truncation is at the TOP:
    mu_h is integrable at the origin and infinite at infinity, so the mass
    beyond eta_m is what gets dropped, and that is what limits small m.
    """
    e = np.concatenate(([0.0], np.asarray(eta, float)))
    mh = -h
    d_low = e[1:] ** mh - e[:-1] ** mh
    d_hi = e[1:] ** (1.0 - h) - e[:-1] ** (1.0 - h)
    w = c_h(h) * d_low / mh
    s = (mh / (1.0 - h)) * d_hi / d_low
    return w, s


def _partition_power(m: int, eta_max: float, p: float) -> np.ndarray:
    """eta_i = eta_max (i/m)^p — p = 1 is uniform, p > 1 clusters near zero."""
    i = np.arange(1, m + 1, dtype=float)
    return eta_max * (i / m) ** p


def _partition_geometric(m: int, eta_max: float, r: float) -> np.ndarray:
    """eta_i = eta_max r^{-(m-i)}, r > 1 — geometric, the AJEE default shape."""
    i = np.arange(1, m + 1, dtype=float)
    return eta_max * r ** (i - m)


def ajee_optimised(m: int, h: float, T: float, t0: float = 0.0,
                   seeds: list | None = None) -> dict:
    """Best AJEE lift for this m: grid search over both partition families,
    then a local Nelder--Mead refinement of the two shape parameters.

    The partition is optimised for the norm actually being used, so that a
    comparison between the t0 = 0 and t0 = delta criteria is fair.

    `seeds` is a list of (family, eta_max, shape) triples added as extra
    Nelder--Mead starting points — pass the optimum found at m-1 so that the
    error is monotone in m.  Without it the search occasionally lands in a worse
    local optimum at large m and the table stops being readable as a rate.
    """
    best = None
    seeds = seeds or []
    families = (("power", _partition_power,
                 np.geomspace(1.0, 1e6, 40), np.linspace(1.0, 6.0, 26)),
                ("geometric", _partition_geometric,
                 np.geomspace(1.0, 1e6, 40), np.geomspace(1.05, 20.0, 26)))
    for name, build, grid_max, grid_shape in families:
        for em in grid_max:
            for sh in grid_shape:
                if m == 1 and name == "geometric":
                    continue                      # degenerate: same as power
                try:
                    w, s = ajee_from_partition(build(m, em, sh), h)
                    err = l2_error(w, s, h, T, t0)
                except Exception:
                    continue
                if np.any(~np.isfinite(s)) or np.any(s <= 0):
                    continue
                if best is None or err < best["l2"]:
                    best = {"family": name, "eta_max": em, "shape": sh,
                            "l2": err, "w": w, "s": s}

        def obj(x, build=build, name=name):
            em, sh = math.exp(x[0]), x[1]
            lo = 1.0 if name == "power" else 1.0000001
            if sh < lo or em <= 0:
                return 1e3
            try:
                w, s = ajee_from_partition(build(m, em, sh), h)
                if np.any(~np.isfinite(s)) or np.any(s <= 0):
                    return 1e3
                return l2_error(w, s, h, T, t0)
            except Exception:
                return 1e3

        starts = []
        if best is not None and best["family"] == name:
            starts.append((best["eta_max"], best["shape"]))
        starts += [(em, sh) for fam, em, sh in seeds if fam == name]
        for em0, sh0 in starts:
            r = minimize(obj, [math.log(em0), sh0], method="Nelder-Mead",
                         options={"xatol": 1e-8, "fatol": 1e-14, "maxiter": 4000})
            if best is None or r.fun < best["l2"]:
                em, sh = math.exp(r.x[0]), r.x[1]
                w, s = ajee_from_partition(build(m, em, sh), h)
                best = {"family": name, "eta_max": em, "shape": sh,
                        "l2": float(r.fun), "w": w, "s": s}
    return best


# ---------------------------------------------------- achievable floor (NNLS)
def _nnls_weights(s: np.ndarray, h: float, T: float,
                  t0: float = 0.0) -> tuple[np.ndarray, float]:
    """Weights >= 0 minimising the L^2 error for FIXED nodes.

    min_w  w'Gw - 2w'b + A  with  G = L L'  gives  min_w ||L'w - L^{-1}b||^2,
    a plain non-negative least-squares problem.  G is a Gram matrix of
    exponentials and is severely ill-conditioned for large m, so a jittered
    Cholesky is used and the resulting error is recomputed from the exact
    formula rather than read off the least-squares residual.
    """
    G = _gram(s, T, t0)
    b = _b_vector(s, h, T, t0)
    jit = 0.0
    scale = float(np.trace(G)) / len(s)
    for _ in range(12):
        try:
            L = np.linalg.cholesky(G + jit * scale * np.eye(len(s)))
            break
        except np.linalg.LinAlgError:
            jit = 1e-14 if jit == 0.0 else jit * 10.0
    else:
        return np.zeros(len(s)), float("inf")
    rhs = np.linalg.solve(L, b)
    # nnls raises RuntimeError('Maximum number of iterations reached') on the
    # worst-conditioned systems -- large m, or H near 0, where the exponentials
    # are nearly linearly dependent. Treat that exactly like the Cholesky
    # failure above: an infinite error, so the outer search walks away from that
    # node configuration instead of the whole run dying on one evaluation.
    try:
        w, _ = nnls(L.T, rhs)
    except (RuntimeError, ValueError):
        try:
            w, _ = nnls(L.T, rhs, maxiter=50 * len(s) ** 2)
        except Exception:
            return np.zeros(len(s)), float("inf")
    return w, l2_error(w, s, h, T, t0)


def best_lift(m: int, h: float, T: float, s_init: np.ndarray,
              t0: float = 0.0, extra_inits: list | None = None) -> dict:
    """Nodes free (in log scale), weights NNLS-optimal.  Local search only:
    this is a 'best found', not a certified global optimum, and it is used as an
    upper bound on the achievable accuracy of an m-factor lift.

    `extra_inits` takes further node vectors to start from — pass the m-1
    optimum padded with one node, so the result is monotone in m.  In 20 or 30
    dimensions a single Nelder--Mead start is not reliable, and without this the
    error at m = 30 can come out worse than at m = 20.
    """
    def obj(x):
        s = np.exp(np.sort(x))
        if np.any(~np.isfinite(s)) or np.any(s <= 0):
            return 1e3
        return _nnls_weights(s, h, T, t0)[1]

    cands = [np.asarray(s_init, float)]
    for e in (extra_inits or []):
        e = np.asarray(e, float)
        if len(e) == m:
            cands.append(e)
        elif len(e) == m - 1:                 # pad: one node beyond the largest
            cands.append(np.append(e, e.max() * max(2.0, e.max() / max(e.min(), 1e-30))))
    budget = 400 * max(m, 4) + 2000
    best = None
    for c in cands:
        if len(c) != m or np.any(c <= 0):
            continue
        r = minimize(obj, np.log(np.sort(c)), method="Nelder-Mead",
                     options={"xatol": 1e-9, "fatol": 1e-15,
                              "maxiter": budget, "maxfev": budget})
        s = np.exp(np.sort(r.x))
        w, err = _nnls_weights(s, h, T, t0)
        if best is None or err < best["l2"]:
            best = {"l2": float(err), "w": w, "s": s,
                    "active_factors": int(np.sum(w > 0))}
    return best


# --------------------------------------------------------- covariance surface
def cov_true(u: float, v: float, h: float) -> float:
    """Cov(V_u, V_v) = int_0^{u^v} K(u-r) K(v-r) dr for the true kernel."""
    lo = min(u, v)
    if lo <= 0.0:
        return 0.0
    g2 = Gamma(1.0 + h) ** 2
    if abs(u - v) < 1e-15:
        return u ** (2.0 * h + 1.0) / ((2.0 * h + 1.0) * g2)

    def f(r):
        return (u - r) ** h * (v - r) ** h
    val, _ = quad(f, 0.0, lo, points=[lo], limit=400,
                  epsabs=1e-13, epsrel=1e-11)
    return val / g2


def cov_lift(u: float, v: float, w: np.ndarray, s: np.ndarray) -> float:
    """Same for K^m, in closed form.

    Written in the stable form: with a = min(u,v), b = max(u,v),

        Cov^m = sum_ij w_i w_j [exp(-s_j (b-a)) - exp(-s_i a - s_j b)] / (s_i+s_j),

    every exponent non-positive.  The naive form carries exp(+(s_i+s_j) a),
    which overflows as soon as the largest node is of order 10^5 — and the
    optimal partitions put nodes at 10^5 and beyond.
    """
    a, b = (u, v) if u <= v else (v, u)
    if a <= 0.0:
        return 0.0
    ss = s[:, None] + s[None, :]
    e1 = np.broadcast_to(np.exp(-s * (b - a))[None, :], ss.shape)
    e2 = np.exp(-s[:, None] * a - s[None, :] * b)
    return float((np.outer(w, w) * (e1 - e2) / ss).sum())


def covariance_report(w: np.ndarray, s: np.ndarray, h: float, T: float,
                      grid: int = 12, t_min: float = 0.0) -> dict:
    """Relative error of the covariance surface on a grid of (u,v) in (0,T]^2.

    This is the quantity Route B exists to fix: Propositions 8.3 and 8.4 say the
    one-step lattice gets it wrong by a factor that DIVERGES, and that its
    covariance is a function of u^v alone.  A lift that passes the L^2 test but
    fails here would be useless.
    """
    ts = np.linspace(max(T / grid, t_min), T, grid)
    worst, worst_at, num, den = 0.0, None, 0.0, 0.0
    for u in ts:
        for v in ts:
            ct = cov_true(float(u), float(v), h)
            cm = cov_lift(float(u), float(v), w, s)
            if ct > 0:
                rel = abs(cm - ct) / ct
                if rel > worst:
                    worst, worst_at = rel, (round(float(u), 4), round(float(v), 4))
            num += (cm - ct) ** 2
            den += ct ** 2
    vT_true = cov_true(T, T, h)
    vT_lift = cov_lift(T, T, w, s)
    return {"t_min": t_min,
            "max_rel_error": worst,
            "max_rel_error_at": worst_at,
            "rel_frobenius_error": math.sqrt(num / den),
            "var_T_true": vT_true,
            "var_T_lift": vT_lift,
            "var_T_ratio": vT_lift / vT_true}


_CONT_CACHE: dict = {}


def continuous_cov_matrix(h: float, T: float, n: int) -> np.ndarray:
    """Cov(V_{k delta}, V_{l delta}) for the TRUE kernel, cached.

    Every entry costs an adaptive quadrature, and the matrix depends only on
    (h, T, n) — not on the lift — so it is computed once per grid.
    """
    key = (round(h, 12), round(T, 12), n)
    if key not in _CONT_CACHE:
        d = T / n
        M = np.empty((n, n))
        for i in range(n):
            for j in range(i, n):
                M[i, j] = M[j, i] = cov_true(float((i + 1) * d),
                                             float((j + 1) * d), h)
        _CONT_CACHE[key] = M
    return _CONT_CACHE[key]


def _cell_average_true(lags: np.ndarray, h: float, d: float) -> np.ndarray:
    """(1/d) int_{(l-1)d}^{ld} K(u) du, the L^2 projection of K on the grid.

    I(x) = int_0^x K = x^{h+1} / ((h+1) Gamma(1+h)), so the cell average is
    (I(ld) - I((l-1)d))/d.  This is the discretisation that matches the
    variance to O(d^{...}) instead of losing 40% of it, and it is the one the
    'exact convolution' column of the first validation run measured.
    """
    a = h + 1.0
    I = lambda x: x ** a / (a * Gamma(1.0 + h))          # noqa: E731
    return (I(lags) - I(np.maximum(lags - d, 0.0))) / d


def _cell_average_lift(lags: np.ndarray, w: np.ndarray, s: np.ndarray,
                       d: float) -> np.ndarray:
    """Same cell average for K^m = sum_i w_i exp(-s_i u), in closed form."""
    lo = np.maximum(lags - d, 0.0)
    num = np.exp(-np.outer(lo, s)) - np.exp(-np.outer(lags, s))
    return (num * (w / s)[None, :]).sum(axis=1) / d


def discrete_covariance_report(w: np.ndarray, s: np.ndarray, h: float,
                               T: float, n: int, mode: str = "cellavg") -> dict:
    """The comparison the SCHEME actually depends on.

    A lattice with n steps builds the volatility factor as the Riemann sum
    sqrt(delta) sum_{j<k} G((k-j) delta) zeta_j, so its covariance is

        C^G_{k,l} = delta sum_{j < k^l} G((k-j)delta) G((l-j)delta),

    and the smallest kernel lag it ever evaluates is delta.  Comparing K^m with
    K in L^2(delta,T) is therefore the right norm for the KERNEL but not for the
    COVARIANCE: C_{k,l} at small k sums lags all the way down to delta, and the
    weight of those short lags in the sum is what a lift optimised away from the
    origin gets wrong.  This function measures the discrepancy directly, with no
    norm standing in for it, on three comparisons:

      lift vs true, both discrete   the error the lift itself introduces;
      true discrete vs continuous   the quadrature error the scheme already has,
                                    so that the first can be judged against it;
      lift discrete vs continuous   the total.
    """
    d = T / n
    lags = np.arange(1, n + 1, dtype=float) * d
    if mode == "left":
        kt, km = K_true(lags, h), K_lift(lags, w, s)
    elif mode == "cellavg":
        kt = _cell_average_true(lags, h, d)
        km = _cell_average_lift(lags, w, s, d)
    else:
        raise ValueError(f"unknown mode {mode!r}")
    if not (np.all(np.isfinite(kt)) and np.all(np.isfinite(km))):
        raise FloatingPointError("non-finite discrete kernel")

    def cov_matrix(g: np.ndarray) -> np.ndarray:
        # M[k,j] = g((k-j) delta) for 0 <= j < k, lower-triangular Toeplitz
        idx = np.arange(n)
        lag = idx[:, None] - idx[None, :]
        M = np.where(lag >= 0, g[np.clip(lag, 0, n - 1)], 0.0)
        with np.errstate(all="ignore"):        # matmul raises spurious FE flags
            C = d * (M @ M.T)
        assert np.all(np.isfinite(C)), "non-finite discrete covariance"
        return C

    Ct, Cm = cov_matrix(kt), cov_matrix(km)
    cont = continuous_cov_matrix(h, T, n)

    def rel(A: np.ndarray, B: np.ndarray) -> dict:
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.abs(A - B) / np.abs(B)
        r = np.where(np.isfinite(r), r, 0.0)
        k = int(np.argmax(r))
        i, j = divmod(k, n)
        return {"max_rel_error": float(r.max()),
                "max_rel_error_at": (round((i + 1) * d, 6), round((j + 1) * d, 6)),
                "rel_frobenius": float(np.linalg.norm(A - B) / np.linalg.norm(B)),
                "var_T_ratio": float(A[-1, -1] / B[-1, -1])}

    return {"n": n, "delta": d, "mode": mode,
            "lift_vs_true_discrete": rel(Cm, Ct),
            "true_discrete_vs_continuous": rel(Ct, cont),
            "lift_discrete_vs_continuous": rel(Cm, cont)}


# ------------------------------------------------------------------ the rates
def lattice_error(n: int, h: float, kappa: float, T: float = 1.0) -> float:
    """The rate of Theorem 7.1: delta^{(h+kappa)/2}, delta = T/n."""
    return (T / n) ** (0.5 * (h + kappa))


def onestep_variance_ratio(n: int, H: float) -> float:
    """Proposition 8.3: Var[Vcheck_T]/Var[V_T] = 2H n^{1-2H}. Diverges."""
    return 2.0 * H * n ** (1.0 - 2.0 * H)


def cost_exponent(m: int, H: float, gamma: float | None = None) -> float:
    """Cost in the target accuracy eps.

    Rate delta^{H/2} forces n ~ eps^{-2/H}; an (m+1)-dimensional backward
    induction costs O(n^{m+1}), and Route A' inside it adds gamma.
    Returns q with cost = eps^{-q}.
    """
    p = (m + 1.0) + (0.0 if gamma is None else gamma)
    return p * 2.0 / H
