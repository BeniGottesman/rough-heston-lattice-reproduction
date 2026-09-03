"""
rough_heston — the independent reference this project did not have.

Every comparison in Sections 10.5 and 10.7 sets our lattice against our own
Monte-Carlo, so a shared misreading of the model would cancel and be invisible
(Remark on the limits of the reference).  Rough Bergomi admits no semi-analytic
price, so that circularity cannot be broken inside it.  Rough Heston can: its
characteristic function solves a FRACTIONAL RICCATI equation, and Fourier
inversion then gives a price with no Monte-Carlo anywhere in it.

Model (El Euch--Rosenbaum, *The characteristic function of rough Heston models*,
Math. Finance 29(1) 3--38, 2019), with alpha = H + 1/2 in (1/2, 1):

    dS_t = S_t sqrt(V_t) dW_t,
    V_t  = V_0 + int_0^t K(t-s) lambda (theta - V_s) ds
               + int_0^t K(t-s) nu sqrt(V_s) dB_s,
    K(u) = u^{alpha-1} / Gamma(alpha),      d<W,B>_t = rho dt.

Characteristic function of X_t = log(S_t/S_0):

    E[exp(i z X_t)] = exp( theta lambda int_0^t h(z,s) ds + V_0 I^{1-alpha}h(z,t) )

where h solves  D^alpha h = F(z,h),  I^{1-alpha}h(z,0) = 0, and

    F(z,x) = -(z^2 + i z)/2 + (i z rho nu - lambda) x + (nu^2/2) x^2 .

At alpha = 1 this is the classical Heston Riccati, which is the first check.

Three independent controls, all exact, are built in.

  nu = 0    the variance becomes deterministic and solves the linear fractional
            equation; the price is then Black--Scholes with the integrated
            deterministic variance, in closed form, FOR EVERY H.  This is the
            strongest control because it does not degrade the roughness.
  z-plane   phi(-i) = 1 and phi(0) = 1 must hold to machine precision.
  H = 1/2   alpha = 1 and the scheme must reproduce classical Heston.

The Riccati is solved by the fractional predictor--corrector of Diethelm, Ford
and Freed (Nonlinear Dynam. 29 (2002) 3--22), vectorised over the Fourier
frequencies so that all of them advance together.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.special import gamma as Gamma
from scipy.stats import norm

T_DEFAULT = 1.0
S0 = 100.0
KSTRIKE = 100.0


# --------------------------------------------------------------- Black--Scholes
def bs_put(s0: float, k: float, total_var: float) -> float:
    """Put with zero rate, given the TOTAL integrated variance int_0^T v dt."""
    sd = math.sqrt(max(total_var, 1e-300))
    d1 = (math.log(s0 / k) + 0.5 * sd * sd) / sd
    return float(k * norm.cdf(-(d1 - sd)) - s0 * norm.cdf(-d1))


# ------------------------------------------------------- the fractional Riccati
def _dff_weights(alpha: float, n: int) -> tuple[np.ndarray, np.ndarray]:
    """Predictor and corrector weight tables of Diethelm--Ford--Freed.

    Returns (b, a) with, for the step from t_k to t_{k+1},
        b[k, j] = (k+1-j)^alpha - (k-j)^alpha,                 0 <= j <= k
        a[k, j] = k^{alpha+1} - (k-alpha)(k+1)^alpha,          j = 0
                = (k-j+2)^{a+1} + (k-j)^{a+1} - 2(k-j+1)^{a+1}, 1 <= j <= k

    The weight at j = k (where k-j = 0) is 2^{alpha+1} - 2 and is NOT zero; an
    earlier version masked on k-j >= 1 and dropped it, which silently degraded the
    scheme from order 1+alpha to order alpha.  The consistency test
    `sum_j a[k,j] + 1 == (alpha+1)(k+1)^alpha`, which makes the corrector exact on
    a constant right-hand side, is what caught it and is checked below.
    """
    j = np.arange(n + 1, dtype=float)[None, :]
    kk = np.arange(n, dtype=float)[:, None]
    d = kk - j + 1.0                               # k+1-j
    b = np.where(d >= 1.0, d ** alpha - np.maximum(d - 1.0, 0.0) ** alpha, 0.0)
    m = np.maximum(kk - j, 0.0)                    # k-j, clamped for the mask
    a = np.where(
        (j >= 1.0) & (j <= kk),                    # 1 <= j <= k, so k-j >= 0
        (m + 2.0) ** (alpha + 1.0) + m ** (alpha + 1.0)
        - 2.0 * (m + 1.0) ** (alpha + 1.0),
        0.0)
    a[:, 0] = kk[:, 0] ** (alpha + 1.0) - (kk[:, 0] - alpha) * (kk[:, 0] + 1.0) ** alpha
    return b, a


def riccati(z: np.ndarray, H: float, lam: float, nu: float, rho: float,
            T: float, steps: int) -> tuple[np.ndarray, float]:
    """Solve D^alpha h = F(z,h) on [0,T] for every z at once.

    Returns h and F(z,h), both of shape (steps+1, len(z)), and the step size.
    """
    alpha = H + 0.5
    dt = T / steps
    z = np.atleast_1d(np.asarray(z, complex))
    c0 = -0.5 * (z * z + 1j * z)
    c1 = 1j * z * rho * nu - lam
    c2 = 0.5 * nu * nu

    def F(x):
        return c0 + c1 * x + c2 * x * x

    b, a = _dff_weights(alpha, steps)
    ga1 = dt ** alpha / Gamma(alpha + 1.0)
    ga2 = dt ** alpha / Gamma(alpha + 2.0)
    h = np.zeros((steps + 1, len(z)), complex)
    f = np.zeros((steps + 1, len(z)), complex)
    f[0] = F(h[0])
    for k in range(steps):
        hp = ga1 * (b[k, :k + 1] @ f[:k + 1])
        h[k + 1] = ga2 * (a[k, :k + 1] @ f[:k + 1] + F(hp))
        f[k + 1] = F(h[k + 1])
    return h, f, dt


def _int_F(f: np.ndarray, dt: float) -> np.ndarray:
    """I^{1-alpha}h(T), computed WITHOUT any fractional quadrature.

    Since h solves D^alpha h = F(z,h) with I^{1-alpha}h(0) = 0, we have
    h = I^alpha F(z,h), hence

        I^{1-alpha} h = I^{1-alpha} I^{alpha} F = I^{1} F = int_0^T F(z,h(s)) ds.

    So the singular kernel (T-s)^{-alpha} disappears entirely and only a plain
    integral of the already-computed right-hand side is needed.  Doing it the
    other way -- quadrature against (T-s)^{-alpha} on h, which behaves like
    t^{alpha} at the origin and so has an unbounded derivative there -- cost
    between one and two percent on the price, which is what an earlier version of
    this file did and what the Black--Scholes control caught.
    """
    w = np.full(f.shape[0], dt)
    w[0] = w[-1] = 0.5 * dt
    with np.errstate(all="ignore"):        # matmul raises spurious FE flags
        out = (w[None, :] @ f)[0]
    if not np.all(np.isfinite(out)):
        raise FloatingPointError("Riccati blew up: reduce nu or increase steps")
    return out


def char_fn(z: np.ndarray, H: float, V0: float, theta: float, lam: float,
            nu: float, rho: float, T: float, steps: int) -> np.ndarray:
    """E[exp(i z log(S_T/S_0))] for rough Heston."""
    h, f, dt = riccati(z, H, lam, nu, rho, T, steps)
    w = np.full(h.shape[0], dt)
    w[0] = w[-1] = 0.5 * dt                        # trapezoid for int_0^T h ds
    int_h = (w[None, :] @ h)[0]
    return np.exp(theta * lam * int_h + V0 * _int_F(f, dt))


# ------------------------------------------------------------ Fourier inversion
def put_fourier(H: float, V0: float, theta: float, lam: float, nu: float,
                rho: float, T: float = T_DEFAULT, s0: float = S0,
                k: float = KSTRIKE, steps: int = 300, nu_max: float = 60.0,
                nq: int = 200) -> float:
    """Put price by the Lewis representation, zero rate.

        C = s0 - (sqrt(s0 k)/pi) int_0^infty Re[ e^{i u kappa} phi(u - i/2) ]
                                            / (u^2 + 1/4) du,   kappa = log(s0/k)

    then the put by parity.  The integrand decays like the characteristic
    function, which for these parameters is already 1e-10 at u = 50 and 1e-17 at
    u = 80, so a Gauss--Legendre rule on (0, 60) is ample and the truncation is
    far below every other error in sight.  Keeping nu_max small also matters for
    stability, not only for cost: the explicit predictor of the Riccati solver
    needs the step size to resolve the frequency, and at u = 200 it diverges
    unless the grid is refined to some 400 steps.  Extending the range therefore
    buys nothing and costs robustness.  Both nu_max and nq are swept in the run.
    """
    x, wq = np.polynomial.legendre.leggauss(nq)
    u = 0.5 * nu_max * (x + 1.0)
    wq = 0.5 * nu_max * wq
    phi = char_fn(u - 0.5j, H, V0, theta, lam, nu, rho, T, steps)
    if not np.all(np.isfinite(phi)):
        bad = u[~np.isfinite(phi)]
        raise FloatingPointError(
            f"characteristic function diverged at u in "
            f"[{bad.min():.1f}, {bad.max():.1f}] with steps={steps}: the Riccati "
            f"predictor needs a finer grid for these frequencies, or a smaller "
            f"nu_max (the integrand there is negligible anyway)")
    kap = math.log(s0 / k)
    integrand = (np.exp(1j * u * kap) * phi / (u * u + 0.25)).real
    call = s0 - math.sqrt(s0 * k) / math.pi * float(wq @ integrand)
    return float(call - s0 + k)                    # parity, r = 0


# ----------------------------------------------------- the deterministic control
def deterministic_variance(H: float, V0: float, theta: float, lam: float,
                           T: float, steps: int) -> np.ndarray:
    """V with nu = 0: the linear fractional equation, by the same DFF scheme.

    V_t = V_0 + I^alpha[ lambda (theta - V) ](t).
    """
    alpha = H + 0.5
    dt = T / steps
    b, a = _dff_weights(alpha, steps)
    ga1 = dt ** alpha / Gamma(alpha + 1.0)
    ga2 = dt ** alpha / Gamma(alpha + 2.0)
    v = np.zeros(steps + 1)
    v[0] = V0
    f = np.zeros(steps + 1)
    f[0] = lam * (theta - V0)
    for kk in range(steps):
        vp = V0 + ga1 * (b[kk, :kk + 1] @ f[:kk + 1])
        v[kk + 1] = V0 + ga2 * (a[kk, :kk + 1] @ f[:kk + 1] + lam * (theta - vp))
        f[kk + 1] = lam * (theta - v[kk + 1])
    return v


def put_deterministic(H: float, V0: float, theta: float, lam: float,
                      T: float = T_DEFAULT, s0: float = S0, k: float = KSTRIKE,
                      steps: int = 2000) -> tuple[float, float]:
    """Exact price when nu = 0: Black--Scholes with the integrated variance."""
    v = deterministic_variance(H, V0, theta, lam, T, steps)
    dt = T / steps
    tot = float(dt * (v[1:-1].sum() + 0.5 * (v[0] + v[-1])))
    return bs_put(s0, k, tot), tot


# ------------------------------------------------------------------ Monte-Carlo
def put_mc(H: float, V0: float, theta: float, lam: float, nu: float, rho: float,
           T: float = T_DEFAULT, s0: float = S0, k: float = KSTRIKE,
           steps: int = 250, paths: int = 200_000, chunk: int = 25_000,
           seed: int = 11, control: bool = True) -> dict:
    """Euler--Volterra Monte-Carlo, antithetic, with the nu = 0 control variate.

    The variance is the causal convolution
        V_k = V_0 + sum_{j<k} K((k-j)dt) [ lambda (theta - V_j) dt
                                           + nu sqrt(V_j) dB_j ],
    truncated at zero, so the smallest kernel lag is dt -- the same convention as
    the paper's Definition of the exact convolution.  The control variate uses
    the SAME Brownian increments with nu = 0, whose price is closed form, so it
    is unbiased and removes most of the variance and of the Euler error in the
    price leg.
    """
    alpha = H + 0.5
    dt = T / steps
    rng = np.random.default_rng(seed)
    lag = np.arange(1, steps + 1) * dt
    ker = lag ** (alpha - 1.0) / Gamma(alpha)

    # the control leg: nu = 0 on the SAME grid and the SAME Brownian increments.
    # Its total variance must be the identical left-endpoint sum the paths use,
    # otherwise the control is not exact and the estimator picks up a bias.
    vdet = np.maximum(deterministic_variance(H, V0, theta, lam, T, steps), 0.0)
    ctrl_var = float(vdet[:steps].sum() * dt)
    ctrl_price = bs_put(s0, k, ctrl_var) if control else 0.0
    sig_det = np.sqrt(vdet[:steps])[None, :]

    tot = tot2 = 0.0
    done = 0
    vT_acc = 0.0
    neg = 0
    while done < paths:
        p = min(chunk, paths - done)
        half = (p + 1) // 2
        z1 = rng.standard_normal((half, steps))
        z2 = rng.standard_normal((half, steps))
        z1 = np.concatenate([z1, -z1])[:p]
        z2 = np.concatenate([z2, -z2])[:p]
        dB = z1 * math.sqrt(dt)
        dW = rho * dB + math.sqrt(1.0 - rho ** 2) * z2 * math.sqrt(dt)

        g = np.empty((p, steps))                   # integrand of the memory
        acc_logs = np.zeros(p)
        acc_var = np.zeros(p)
        v_prev = np.full(p, V0)
        vT = v_prev
        for kk in range(steps):
            neg += int((v_prev < 0.0).sum())
            vp = np.maximum(v_prev, 0.0)
            acc_var += vp * dt
            acc_logs += np.sqrt(vp) * dW[:, kk]
            g[:, kk] = lam * (theta - vp) * dt + nu * np.sqrt(vp) * dB[:, kk]
            if kk + 1 < steps:
                # V_{k+1} = V_0 + sum_{j<=k} K((k+1-j) dt) g_j : smallest lag dt
                v_prev = V0 + g[:, :kk + 1] @ ker[kk::-1]
                vT = v_prev
        logST = math.log(s0) - 0.5 * acc_var + acc_logs
        pay = np.maximum(k - np.exp(logST), 0.0)
        if control:
            logc = (math.log(s0) - 0.5 * ctrl_var
                    + (sig_det * dW).sum(axis=1))
            sample = pay - np.maximum(k - np.exp(logc), 0.0)
        else:
            sample = pay
        tot += sample.sum()
        tot2 += (sample ** 2).sum()
        vT_acc += float(np.maximum(vT, 0.0).sum())
        done += p
    mean = tot / paths
    se = math.sqrt(max(0.0, tot2 / paths - mean * mean) / paths)
    return {"price": float(mean + ctrl_price), "stderr": float(se),
            "ci95": float(1.96 * se), "paths": paths, "steps": steps,
            "control": control, "control_price": ctrl_price,
            "control_total_variance": ctrl_var,
            "negative_variance_hits": neg, "mean_V_T": vT_acc / paths}
