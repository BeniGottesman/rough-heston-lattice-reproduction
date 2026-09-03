"""
route_b_lattice — Route B (the Markovian lift) built as an actual LATTICE.

Everything Route B had until now was a kernel-approximation calculation:
||K - K^m||, the covariance surface, m*(n), a cost exponent.  No lattice was
ever constructed, so the one question a lattice exists to answer -- does the
lift repair the price? -- had never been asked.  It also hid a defect.

The defect
----------
The paper's lifted factor is

    Z^i_k = exp(-s_i delta) ( Z^i_{k-1} + sqrt(delta) zeta_k ),     zeta = +-1,

equivalently  Z^i_k = sqrt(delta) sum_{j<=k} exp(-s_i (k-j+1) delta) zeta_j.
Detrending by exp(s_i t_k) turns this into a walk with deterministic but
UNEQUAL step magnitudes,  sqrt(delta) exp(s_i t_{j-1}),  and a walk whose step
magnitudes differ does not recombine: Z^i_k takes 2^k distinct values, not
k+1.  Only s_i = 0 recombines, and s_i = 0 is the constant kernel, i.e. the
one-step scheme the lift is meant to replace.  So Route B as stated does not
make the dynamic programme feasible, and the O(n^m) state count claimed for it
is not attained.  (Verified by exact enumeration; see the run's phase 1.)

The repair
----------
Do to each factor what Route A' already does to the price: randomised rounding
onto a grid, which is mean-exact and therefore costs only a martingale.  Work
with the WEIGHTED factors

    Zhat^i := w_i Z^i,     V^H_k = sum_i Zhat^i_k,
    Zhat^i_k = exp(-s_i delta) ( Zhat^i_{k-1} + w_i sqrt(delta) zeta_k ),

and put every Zhat^i on the SAME spacing a.  Two consequences, both used below:

  * the rounding error of Zhat^i is a martingale with increments <= a/2, so the
    error in V^H is at most m sqrt(n) a / 2, which is O(delta^gamma) for
    a = 2 delta^{gamma + 1/2} / (m sqrt(T)) -- the rate of Proposition 7.5, so
    the repair is free at the order at which the scheme works;
  * V^H = a * (sum_i j_i) is a function of the SUM of the grid indices alone,
    so the price kernel has to be built for O(sum_i N_i) values and not
    O(prod_i N_i).

The node count per factor is N_i ~ w_i sd(Z^i) / a.  The weights blow up like
K(delta) as delta falls -- sum_i w_i = K^m(0) chases the singularity -- but
sd(Z^i) ~ (2 s_i)^{-1/2} shrinks at exactly the compensating rate, so the
products w_i sd(Z^i) stay O(1) and the state space is
O(n^{(m+1)(1/2+gamma)}) with an O(1) constant.  That is the corrected cost.

Normalisation
-------------
The rough Bergomi code of this project uses K(u) = sqrt(2H) u^h, for which
Var[V^H_t] = t^{2H}.  route_b.py optimises lifts for K_true(u) = u^h/Gamma(1+h).
Both K and K^m are linear in the weights, so a lift for one is a lift for the
other after multiplying the weights by sqrt(2H) Gamma(1+h); `lift_for` does it.

Barriers
--------
Each Zhat^i is absorbed at +- zbar_z standard deviations, and log S at
+- zbar_x standard deviations.  Both are in the model (Assumption on the
coefficients puts absorbing barriers on X) and both costs are measured by
widening them, as in the Heston lattice of section 10.8.3.
"""
from __future__ import annotations

import itertools
import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import route_b as rb                                          # noqa: E402
from route_aprime import interp_kernel, T, S0, KSTRIKE, XI0   # noqa: E402


# ------------------------------------------------------------------ the lift
_LIFT_CACHE: dict = {}


def lift_for(n: int, m: int, H: float, T_: float = T) -> tuple[np.ndarray, np.ndarray]:
    """(what, s) for K(u) = sqrt(2H) u^h, optimised in L^2(delta, T).

    Nodes free, weights non-negative least squares -- the 'achievable floor'
    branch of route_b.best_lift, which is the branch m*(n) was measured in.
    Warm-started from m-1 so the error is monotone in m.
    """
    key = (n, m, H, T_)
    if key in _LIFT_CACHE:
        return _LIFT_CACHE[key]
    h = H - 0.5
    d = T_ / n
    scale = math.sqrt(2.0 * H) * float(rb.Gamma(1.0 + h))
    prev = None
    for mm in range(1, m + 1):
        s0 = np.geomspace(1.0, 50.0, mm) if mm > 1 else np.array([2.0])
        b = rb.best_lift(mm, h, T_, s0, t0=d,
                         extra_inits=[prev] if prev is not None else None)
        prev = b["s"]
        _LIFT_CACHE[(n, mm, H, T_)] = (b["w"] * scale, b["s"])
    return _LIFT_CACHE[key]


def factor_sd(s: np.ndarray, n: int, T_: float = T) -> np.ndarray:
    """sd(Z^i_k) maximised over k, for the unit-weight factor."""
    d = T_ / n
    lags = np.arange(1, n + 1) * d
    return np.sqrt(np.array([d * np.sum(np.exp(-2.0 * si * lags)) for si in s]))


def exact_var_VH(n: int, H: float, T_: float = T,
                 w: np.ndarray | None = None,
                 s: np.ndarray | None = None) -> float:
    """Var[V^H_n] for the exact-convolution lattice, true kernel or lift."""
    d = T_ / n
    lags = np.arange(1, n + 1) * d
    if w is None:
        K = math.sqrt(2.0 * H) * lags ** (H - 0.5)
    else:
        K = (w[None, :] * np.exp(-np.outer(lags, s))).sum(axis=1)
    return float(d * np.sum(K ** 2))


# --------------------------------------------------------- state-count checks
def exact_states(k: int, s: float, d: float, tol: float = 1e-12) -> int:
    """Distinct values of the UNROUNDED factor after k steps, by enumeration."""
    vals = np.array([0.0])
    dec, step = math.exp(-s * d), math.sqrt(d)
    for _ in range(k):
        vals = np.sort(np.concatenate([dec * (vals + step), dec * (vals - step)]))
        keep = [vals[0]]
        cut = tol * max(1.0, float(np.abs(vals).max()))
        for v in vals[1:]:
            if v - keep[-1] > cut:
                keep.append(v)
        vals = np.array(keep)
    return len(vals)


# ------------------------------------------------------------- the lattice
class LiftedLattice:
    """The (m+1)-dimensional recombining lattice of Route A' + Route B."""

    def __init__(self, n: int, H: float, eta: float, rho: float, m: int,
                 mref: int = 4, mref_z: int = 4,
                 zbar_z: float = 3.0, zbar_x: float = 4.5,
                 T_: float = T):
        self.n, self.H, self.eta, self.rho, self.m = n, H, eta, rho, m
        self.T = T_
        self.d = d = T_ / n
        self.sqd = sqd = math.sqrt(d)
        self.what, self.s = lift_for(n, m, H, T_)
        self.sd_z = factor_sd(self.s, n, T_)

        # ---- factor grids: one common spacing, mref_z controls the rounding
        # error in V^H units:  m sqrt(n) a / 2 = sqrt(T) / (2 mref_z).
        self.a = sqd / (m * mref_z)
        self.rounding_budget = math.sqrt(T_) / (2.0 * mref_z)
        halfw = zbar_z * self.what * self.sd_z                 # range of Zhat^i
        self.half = np.maximum(1, np.ceil(halfw / self.a).astype(int))
        self.N = 2 * self.half + 1                             # nodes per factor
        self.S = int(np.prod(self.N))
        # flattened index <-> per-factor indices
        self.grids = [np.arange(N) - hf for N, hf in zip(self.N, self.half)]
        mesh = np.meshgrid(*self.grids, indexing="ij")
        self.jsum = sum(g.ravel() for g in mesh)               # (S,) integer sum
        self.strides = np.array(
            [int(np.prod(self.N[i + 1:])) for i in range(m)], dtype=np.int64)

        # ---- price grid
        s_ref = math.sqrt(XI0)
        self.a_X = a_X = s_ref * sqd / mref
        self.mref = mref
        vmax = self._variance_from_jsum(np.array([self.jsum.max()]), n).max()
        sig_max = math.sqrt(float(vmax))
        self.kmax = int(np.ceil(4.5 * sig_max * math.sqrt(1 - rho ** 2) * sqd / a_X
                                + abs(rho) * sig_max * sqd / a_X
                                + 0.5 * sig_max ** 2 * d / a_X)) + 1
        self.offsets = np.arange(-self.kmax, self.kmax + 1)
        reach_x = zbar_x * s_ref * math.sqrt(T_)
        self.ix_half = int(np.ceil(reach_x / a_X))
        self.nx = 2 * self.ix_half + 1
        self.logS = math.log(S0) + (np.arange(self.nx) - self.ix_half) * a_X
        self.payoff = np.maximum(KSTRIKE - np.exp(self.logS), 0.0)

        self._build_transitions()

    # -------------------------------------------------------------- variance
    def _variance_from_jsum(self, jsum, k: int) -> np.ndarray:
        """v_k = xi0 exp(eta V^H_k - 0.5 eta^2 t_k^{2H}),  V^H = a * jsum."""
        if self.eta == 0.0:
            return np.full(np.shape(jsum), XI0, dtype=float)
        VH = self.a * np.asarray(jsum, dtype=float)
        t = k * self.d
        return XI0 * np.exp(self.eta * VH
                            - 0.5 * self.eta ** 2 * t ** (2 * self.H))

    # ----------------------------------------------------------- transitions
    def _build_transitions(self) -> None:
        """For each zeta, per factor: lower successor index and its probability.

        Randomised rounding is mean-exact on the interior; at the barrier the
        target is projected onto the extreme node, which is absorption and is
        the only place the mean is not preserved.  `absorbed_mass` records how
        much probability that projection moves, so the barrier cost is visible.
        """
        self.tr = {}
        self.absorbed_mass = 0.0
        for zeta in (+1.0, -1.0):
            per = []
            for i in range(self.m):
                g = self.grids[i].astype(float)                 # index coords
                z = g * self.a
                tgt = math.exp(-self.s[i] * self.d) * (
                    z + self.what[i] * self.sqd * zeta)
                ti = tgt / self.a
                lim = float(self.half[i])
                clipped = np.clip(ti, -lim, lim)
                self.absorbed_mass = max(
                    self.absorbed_mass,
                    float(np.abs(ti - clipped).max()) * self.a)
                lo = np.floor(clipped).astype(np.int64)
                frac = clipped - lo
                lo = np.clip(lo, -self.half[i], self.half[i])
                hi = np.clip(lo + 1, -self.half[i], self.half[i])
                per.append((lo + self.half[i], hi + self.half[i], frac))
            self.tr[zeta] = per

    def _contract_factors(self, val: np.ndarray, zeta: float) -> np.ndarray:
        """U[ix, J] = sum over successor states J' of p(J -> J') val[ix, J'].

        The successor law factorises over the m factors, so the sum runs over
        the 2^m corners of a box.
        """
        per = self.tr[zeta]
        U = np.zeros_like(val)
        for corner in itertools.product((0, 1), repeat=self.m):
            idx = np.zeros(self.S, dtype=np.int64)
            pr = np.ones(self.S)
            for i, c in enumerate(corner):
                lo, hi, frac = per[i]
                # broadcast the per-factor arrays over the flattened state
                shape = [1] * self.m
                shape[i] = self.N[i]
                sel_lo = np.broadcast_to(lo.reshape(shape), tuple(self.N)).ravel()
                sel_hi = np.broadcast_to(hi.reshape(shape), tuple(self.N)).ravel()
                sel_fr = np.broadcast_to(frac.reshape(shape), tuple(self.N)).ravel()
                idx += self.strides[i] * (sel_hi if c else sel_lo)
                pr = pr * (sel_fr if c else (1.0 - sel_fr))
            U += pr[None, :] * val[:, idx]
        return U

    # ------------------------------------------------------------ the DP
    def price(self, american: bool = False, pr=None, label: str = "") -> dict:
        n, nx, S = self.n, self.nx, self.S
        kmax = self.kmax
        val = np.repeat(self.payoff[:, None], S, axis=1)
        neg = 0
        mass_err = 0.0
        for k in range(n - 1, -1, -1):
            cont = np.zeros((nx, S))
            v_here = self._variance_from_jsum(self.jsum, k)
            sig = np.sqrt(v_here)
            sd_M = sig * math.sqrt(1.0 - self.rho ** 2) * self.sqd
            for zeta in (+1.0, -1.0):
                U = self._contract_factors(val, zeta)
                # pad so the price shift is a slice, and absorb at the barrier
                Up = np.empty((nx + 2 * kmax, S))
                Up[kmax:kmax + nx] = U
                Up[:kmax] = self.payoff[0]
                Up[kmax + nx:] = self.payoff[-1]
                mu_M = sig * self.rho * self.sqd * zeta - 0.5 * v_here * self.d
                P = interp_kernel(mu_M, sd_M, self.a_X, self.offsets)   # (K,S)
                neg += int((P < -1e-12).sum())
                mass_err = max(mass_err,
                               float(np.abs(P.sum(axis=0) - 1.0).max()))
                acc = np.zeros((nx, S))
                for oi, o in enumerate(self.offsets):
                    acc += P[oi][None, :] * Up[kmax + o:kmax + o + nx]
                cont += 0.5 * acc
            val = np.maximum(self.payoff[:, None], cont) if american else cont
            # the price barrier is absorbing: value there is the payoff
            val[0, :] = self.payoff[0]
            val[-1, :] = self.payoff[-1]
            if pr is not None and n >= 8 and (k % max(1, n // 8) == 0):
                pr.tick(pr.step, inner=f"{label} k={k}")
        j0 = int(np.ravel_multi_index(tuple([hf for hf in self.half]),
                                      tuple(self.N)))
        return {"value": float(val[self.ix_half, j0]),
                "negative_probabilities": neg,
                "max_mass_error": round(mass_err, 12),
                "states_per_factor": self.N.tolist(),
                "state_space": self.S,
                "nx": self.nx,
                "kmax": int(self.kmax),
                "grid_a": self.a,
                "rounding_budget_VH": self.rounding_budget,
                "barrier_projection": self.absorbed_mass,
                "w_hat": self.what.tolist(),
                "s": self.s.tolist(),
                "w_sd_products": (self.what * self.sd_z).tolist()}

    # ------------------------------------------------- the telescope, measured
    def price_eta0(self) -> float:
        """The same price scheme with v == xi0: a one-dimensional DP.

        At eta = 0 the driver drops out of the price entirely, so this is the
        exact mean of the Monte-Carlo's eta = 0 payoff under the identical
        rounding and the identical barrier.  It is therefore an exact control
        variate, and it also pins the price scheme's own discretisation error.
        """
        nx, kmax = self.nx, self.kmax
        v = float(XI0)
        sig = math.sqrt(v)
        sd_M = sig * math.sqrt(1.0 - self.rho ** 2) * self.sqd
        val = self.payoff.copy()
        for _ in range(self.n):
            acc = np.zeros(nx)
            for zeta in (+1.0, -1.0):
                mu = np.array([sig * self.rho * self.sqd * zeta - 0.5 * v * self.d])
                P = interp_kernel(mu, np.array([sd_M]), self.a_X,
                                  self.offsets)[:, 0]
                Up = np.empty(nx + 2 * kmax)
                Up[kmax:kmax + nx] = val
                Up[:kmax] = self.payoff[0]
                Up[kmax + nx:] = self.payoff[-1]
                for oi, o in enumerate(self.offsets):
                    acc += 0.5 * P[oi] * Up[kmax + o:kmax + o + nx]
            val = acc
            val[0] = self.payoff[0]
            val[-1] = self.payoff[-1]
        return float(val[self.ix_half])

    def mc_same_price_scheme(self, driver: str, paths: int = 400_000,
                             chunk: int = 50_000, seed: int = 23,
                             control: bool = True) -> dict:
        """Monte-Carlo with the lattice's OWN price scheme and a chosen driver.

        The point is to isolate the two terms of the telescope
            |Lambda - Lambda^{(n,m)}| <= |Lambda - Lambda^{(n)}|
                                       + |Lambda^{(n)} - Lambda^{(n,m)}|.
        Holding the price scheme fixed and swapping only the driver turns each
        term into a difference of two numbers that differ in one ingredient:

          'exact'        V^{(n)}: the exact convolution, true kernel, all n lags
          'lift'         V^{(n,m)}: the lift, unrounded  -> isolates ||K-K^m||_n
          'lift_rounded' the lift as the LATTICE runs it  -> must reproduce the
                         lattice's own DP value to Monte-Carlo error, which is
                         the implementation control
          'onestep'      Vcheck: the one-step scheme of Definition Vcheck
        """
        rng = np.random.default_rng(seed)
        n, d, sqd = self.n, self.d, self.sqd
        lags = np.arange(1, n + 1) * d
        tpow = (np.arange(n) * d) ** (2 * self.H)
        # lower-triangular kernel matrix G[k, j] = K((k-j+1) delta), j <= k
        kk, jj = np.tril_indices(n)
        if driver == "exact":
            Kv = math.sqrt(2.0 * self.H) * lags ** (self.H - 0.5)
        elif driver in ("lift", "lift_rounded"):
            Kv = (self.what[None, :] * np.exp(-np.outer(lags, self.s))).sum(axis=1)
        elif driver == "onestep":
            Kv = np.full(n, math.sqrt(2.0 * self.H) * d ** (self.H - 0.5))
        else:
            raise ValueError(driver)
        G = np.zeros((n, n))
        G[kk, jj] = Kv[kk - jj]
        s_ref = math.sqrt(XI0)
        tot = tot2 = 0.0
        done = 0
        vh_last = []
        while done < paths:
            p = min(chunk, paths - done)
            zeta = np.where(rng.random((p, n)) < 0.5, 1.0, -1.0)
            if driver == "lift_rounded":
                VH = np.zeros((p, n))
                j = np.zeros((p, self.m), dtype=np.int64)
                u = rng.random((p, n, self.m))
                for k in range(n):
                    for zi, zv in enumerate((+1.0, -1.0)):
                        sel = zeta[:, k] == zv
                        if not sel.any():
                            continue
                        per = self.tr[zv]
                        for i in range(self.m):
                            lo, hi, frac = per[i]
                            gi = j[sel, i] + self.half[i]
                            f = frac[gi]
                            nxt = np.where(u[sel, k, i] < f, hi[gi], lo[gi])
                            col = j[:, i].copy()
                            col[sel] = nxt - self.half[i]
                            j[:, i] = col
                    VH[:, k] = self.a * j.sum(axis=1)
                VH = np.concatenate([np.zeros((p, 1)), VH[:, :-1]], axis=1)
            else:
                # Apple's BLAS leaves spurious divide/overflow/invalid flags set
                # on this matmul: the inputs are finite (G is a truncated power
                # kernel, zeta is +-1) and the output agrees with a non-BLAS
                # einsum to 9e-16.  Suppress the flags, assert what matters.
                with np.errstate(divide="ignore", over="ignore",
                                 invalid="ignore"):
                    VH = sqd * (zeta @ G.T)
                if not np.isfinite(VH).all():
                    raise FloatingPointError(
                        f"non-finite driver for {driver!r} at n={n}")
                if driver == "onestep":
                    # the one-step scheme caps its rough coordinate zc at
                    # zmax = 3/sqrt(2H), i.e. it caps V^H = sqrt(2H) zc at 3
                    VH = np.clip(VH, -3.0, 3.0)
                VH = np.concatenate([np.zeros((p, 1)), VH[:, :-1]], axis=1)
            v = XI0 * np.exp(self.eta * VH - 0.5 * self.eta ** 2 * tpow[None, :]) \
                if self.eta > 0 else np.full((p, n), XI0)
            Gn = rng.standard_normal((p, n))
            U = rng.random((p, n))

            def run(vv: np.ndarray) -> np.ndarray:
                """The lattice's price scheme, driven by the variance path vv."""
                sig = np.sqrt(vv)
                ix = np.zeros(p, dtype=np.int64)
                alive = np.ones(p, dtype=bool)
                for k in range(n):
                    mu_M = (sig[:, k] * self.rho * sqd * zeta[:, k]
                            - 0.5 * vv[:, k] * d)
                    sd_M = sig[:, k] * math.sqrt(1.0 - self.rho ** 2) * sqd
                    ti = (mu_M + sd_M * Gn[:, k]) / self.a_X
                    lo = np.floor(ti)
                    step = (lo + (U[:, k] < (ti - lo))).astype(np.int64)
                    ix = np.where(alive, ix + step, ix)
                    hit = np.abs(ix) >= self.ix_half
                    ix = np.clip(ix, -self.ix_half, self.ix_half)
                    alive = alive & ~hit
                return np.maximum(
                    KSTRIKE - np.exp(math.log(S0) + ix * self.a_X), 0.0)

            pay = run(v)
            if control and self.eta > 0.0:
                # same zeta, same G, same U, v == xi0: an exact control whose
                # mean is price_eta0(), so only the eta-effect is sampled
                pay = pay - run(np.full((p, n), float(XI0)))
            tot += pay.sum()
            tot2 += (pay ** 2).sum()
            vh_last.append(VH[:, -1])
            done += p
        mean = tot / paths
        var = max(0.0, tot2 / paths - mean ** 2)
        shift = self.price_eta0() if (control and self.eta > 0.0) else 0.0
        return {"price": float(mean + shift), "stderr": float(math.sqrt(var / paths)),
                "control_shift": float(shift),
                "var_VH_penultimate": float(np.concatenate(vh_last).var())}

    # --------------------------------------------------- forward diagnostics
    def simulate_VH_variance(self, paths: int = 200_000, seed: int = 11) -> dict:
        """Var[V^H_n] realised by the ROUNDED chain, vs the lift's own target.

        Isolates what the rounding and the barrier cost in the driver alone,
        with no option pricing in the way.
        """
        rng = np.random.default_rng(seed)
        j = np.zeros((paths, self.m), dtype=np.int64)
        for _ in range(self.n):
            zeta = np.where(rng.random(paths) < 0.5, 0, 1)      # index into tr
            u = rng.random((paths, self.m))
            for zi, zv in enumerate((+1.0, -1.0)):
                sel = zeta == zi
                if not sel.any():
                    continue
                per = self.tr[zv]
                for i in range(self.m):
                    lo, hi, frac = per[i]
                    gi = j[sel, i] + self.half[i]
                    f = frac[gi]
                    up = u[sel, i] < f
                    nxt = np.where(up, hi[gi], lo[gi]) - self.half[i]
                    col = j[:, i]
                    col[sel] = nxt
                    j[:, i] = col
        VH = self.a * j.sum(axis=1)
        tgt = exact_var_VH(self.n, self.H, self.T, self.what, self.s)
        return {"var_chain": float(VH.var()),
                "var_lift_target": tgt,
                "var_true_kernel": exact_var_VH(self.n, self.H, self.T),
                "var_continuous": self.T ** (2 * self.H),
                "ratio_chain_over_target": float(VH.var() / tgt)}
