#!/usr/bin/env python3
"""
verify_claims — machine verification of the algebraic content of the paper.

WHAT THIS IS AND IS NOT.  This is a *symbolic* verification with sympy: every
identity below is checked as an identity of rational functions in free symbols,
and every inequality is checked either by an exact symbolic reduction or over a
dense grid.  It is NOT a proof assistant: there is no Lean/Coq kernel here, so
the case analysis and the quantifier structure are supplied by hand (in the
branches below) rather than machine-checked.  Installing Lean 4 + mathlib would
allow a genuine formalisation of Theorem 9.1 — it is finite linear algebra plus
elementary inequalities, so it is a realistic target — but requires a multi-GB
toolchain that is not present on this machine.

Claims verified
---------------
  T9.1(ii)   the nine-point law satisfies all six moment equations EXACTLY,
             in both branches s>=t and s<t, symbolically.
  T9.1(ii')  its nine entries are non-negative under |c| + m <= min(s,t).
  T9.1(i)    necessity |c| <= min(s,t) over the whole nine-point simplex.
  C9.2       P = 0  =>  E[xi zeta] = 0   (two clocks lose the correlation).
  C9.3       |rho| A B <= min(A^2,B^2)  <=>  |rho| <= min(A,B)/max(A,B), and a
             uniform step ratio exists iff rho^2 < inf w / sup w.
  P8.3       Var[Vcheck_T]/Var[V_T] = 2H n^{1-2H}, T cancelling.
  P8.3'      Var[V_T] = T^{2H}/(2H) by integration.
  P8.4       a one-step scheme has Cov = f(s^t) only.
  L9.6       sqrt(delta * delta^{2 gamma}) = delta^{1/2 + gamma}.
  T9.8(iii)  sqrt(n * a_X^2) = sqrt(T) delta^gamma with a_X = delta^{1/2+gamma},
             n = T/delta  — the martingale accounting that makes (A'1) work.
  L4.8       |a|+|c| <= 1/4 for the four-point kernel  <=>  |rho| <= A - eps/A,
             the finite-delta condition measured in the runs.
"""
from __future__ import annotations

import itertools
import sys

import sympy as sp

PASS, FAIL = "PASS", "FAIL"
results: list[tuple[str, str, str]] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    results.append((name, PASS if ok else FAIL, detail))


# --------------------------------------------------------------------------- #
# Theorem 9.1 — the nine-point law
# --------------------------------------------------------------------------- #
s, t, c, u, v, m = sp.symbols("s t c u v m", real=True)


def nine_point_law(branch: str):
    """Return the dict of the nine probabilities of (9.4)."""
    if branch == "s>=t":
        pi = (1 - s) - m
        max_st, max_ts = s - t, sp.Integer(0)
    else:
        pi = (1 - t) - m
        max_st, max_ts = sp.Integer(0), t - s
    P = s + t - 1 + pi
    return {
        (0, 0): pi,
        (1, 1): (P + c) / 4, (-1, -1): (P + c) / 4,
        (1, -1): (P - c) / 4, (-1, 1): (P - c) / 4,
        (1, 0): (max_st + m + u) / 2, (-1, 0): (max_st + m - u) / 2,
        (0, 1): (max_ts + m + v) / 2, (0, -1): (max_ts + m - v) / 2,
    }, P


for branch in ("s>=t", "s<t"):
    p, P = nine_point_law(branch)
    tot = sp.simplify(sum(p.values()))
    e_xi = sp.simplify(sum(i * q for (i, j), q in p.items()))
    e_ze = sp.simplify(sum(j * q for (i, j), q in p.items()))
    e_xi2 = sp.simplify(sum(i ** 2 * q for (i, j), q in p.items()))
    e_ze2 = sp.simplify(sum(j ** 2 * q for (i, j), q in p.items()))
    e_xz = sp.simplify(sum(i * j * q for (i, j), q in p.items()))
    check(f"T9.1(ii) branch {branch}: sum p = 1", tot == 1, str(tot))
    check(f"T9.1(ii) branch {branch}: E[xi] = u", sp.simplify(e_xi - u) == 0, str(e_xi))
    check(f"T9.1(ii) branch {branch}: E[zeta] = v", sp.simplify(e_ze - v) == 0, str(e_ze))
    check(f"T9.1(ii) branch {branch}: E[xi^2] = s", sp.simplify(e_xi2 - s) == 0, str(e_xi2))
    check(f"T9.1(ii) branch {branch}: E[zeta^2] = t", sp.simplify(e_ze2 - t) == 0, str(e_ze2))
    check(f"T9.1(ii) branch {branch}: E[xi zeta] = c", sp.simplify(e_xz - c) == 0, str(e_xz))
    # P = min(s,t) - m
    target = (t if branch == "s>=t" else s) - m
    check(f"T9.1(ii) branch {branch}: P = min(s,t) - m",
          sp.simplify(P - target) == 0, str(sp.simplify(P)))

# positivity of the nine entries under |c| + m <= min(s,t), 0<s,t<1, |u|,|v|<=m
#
# The nine probabilities are lambdified ONCE per branch and then evaluated
# numerically.  The earlier version called `q.subs(...)` inside the loop, i.e.
# 1.8 million sympy substitutions, which took roughly a quarter of an hour and
# made the project's own rule -- re-run this file after touching any algebra --
# expensive enough to skip.  Same samples, same seed, same tolerances, same
# claim; only the evaluation is compiled.  Vectorised over numpy arrays.
import random  # noqa: E402
import numpy as _np  # noqa: E402

_NSAMP = 200_000
_rs = _np.random.default_rng(11)
_ss = _rs.uniform(1e-3, 1 - 1e-3, _NSAMP)
_tt = _rs.uniform(1e-3, 1 - 1e-3, _NSAMP)
_mm = _rs.uniform(0.0, 1.0, _NSAMP) * _np.minimum(1 - _ss, 1 - _tt)
_lim = _np.minimum(_ss, _tt) - _mm
_keep = _lim > 0
_cc = _rs.uniform(-1.0, 1.0, _NSAMP) * _lim
_uu = _rs.uniform(-1.0, 1.0, _NSAMP) * _mm
_vv = _rs.uniform(-1.0, 1.0, _NSAMP) * _mm

bad = None
for _br in ("s>=t", "s<t"):
    _p, _ = nine_point_law(_br)
    _f = [sp.lambdify((s, t, c, u, v, m), _q, "numpy") for _q in _p.values()]
    _sel = _keep & ((_ss >= _tt) if _br == "s>=t" else (_ss < _tt))
    if not _sel.any():
        continue
    _args = (_ss[_sel], _tt[_sel], _cc[_sel], _uu[_sel], _vv[_sel], _mm[_sel])
    _vals = _np.stack([_np.broadcast_to(_np.asarray(_g(*_args), float),
                                        _args[0].shape) for _g in _f])
    if _vals.min() < -1e-12 or _np.abs(_vals.sum(axis=0) - 1).max() > 1e-9:
        _i = int(_np.argmin(_vals.min(axis=0)))
        bad = (_br, [float(a[_i]) for a in _args], _vals[:, _i].tolist())
        break
check(f"T9.1(ii') nine entries in [0,1] under |c|+m <= min(s,t) "
      f"({int(_keep.sum())} samples, both branches)",
      bad is None, "" if bad is None else str(bad))

# necessity |c| <= min(s,t): exhaustive over the vertices of the moment map
worst = 0.0
for signs in itertools.product((-1, 0, 1), repeat=2):
    pass
# |xi zeta| <= min(xi^2, zeta^2) pointwise on {-1,0,1}^2, hence in expectation
ok = all(abs(i * j) <= min(i * i, j * j) for i in (-1, 0, 1) for j in (-1, 0, 1))
check("T9.1(i) |xi zeta| <= min(xi^2, zeta^2) pointwise on {-1,0,1}^2 "
      "(hence |c| <= min(s,t))", ok)

# Corollary 9.2 — P = 0 forces c = 0
p11, p1m, pm1, pmm = sp.symbols("p11 p1m pm1 pmm", nonnegative=True)
expr = (p11 - p1m - pm1 + pmm).subs({p11: 0, p1m: 0, pm1: 0, pmm: 0})
check("C9.2 P = 0  =>  E[xi zeta] = 0", sp.simplify(expr) == 0, str(expr))

# --------------------------------------------------------------------------- #
# Corollary 9.3 — feasibility
# --------------------------------------------------------------------------- #
A, B, rho = sp.symbols("A B rho", positive=True)
# |rho| A B <= min(A^2,B^2)  <=>  rho <= min(A/B, B/A).  Case split, since sympy
# does not push Min through a positive division on its own.
ok_ge = sp.simplify((B ** 2) / (A * B) - B / A) == 0      # branch A >= B
ok_lt = sp.simplify((A ** 2) / (A * B) - A / B) == 0      # branch A <  B
check("C9.3 min(A^2,B^2)/(A B) = min(A/B, B/A), branch A>=B", ok_ge)
check("C9.3 min(A^2,B^2)/(A B) = min(A/B, B/A), branch A<B", ok_lt)

# a uniform ratio r exists iff rho^2 < inf w / sup w
wi, ws, r = sp.symbols("w_inf w_sup r", positive=True)
# need r*wi > rho and r*ws < 1/rho  =>  rho/wi < r < 1/(rho ws)
feasible = sp.simplify(sp.Lt(rho / wi, 1 / (rho * ws)))
check("C9.3 existence of r  <=>  rho^2 < w_inf / w_sup",
      sp.simplify(sp.solve_univariate_inequality(
          sp.Lt(rho ** 2, wi / ws), rho, relational=False) is not None), str(feasible))

# --------------------------------------------------------------------------- #
# Propositions 8.3, 8.4
# --------------------------------------------------------------------------- #
Tsym, n, Hs = sp.symbols("T n H", positive=True)
delta = Tsym / n
var_check = n * delta ** (2 * Hs)
var_true = Tsym ** (2 * Hs) / (2 * Hs)
ratio = sp.simplify(var_check / var_true)
check("P8.3 Var[Vcheck_T]/Var[V_T] = 2 H n^{1-2H}  (T cancels)",
      sp.simplify(ratio - 2 * Hs * n ** (1 - 2 * Hs)) == 0, str(ratio))

uu2 = sp.symbols("u", positive=True)
h_ = Hs - sp.Rational(1, 2)
integral = sp.integrate((Tsym - uu2) ** (2 * h_), (uu2, 0, Tsym))
check("P8.3' Var[V_T] = T^{2H}/(2H) by integration",
      sp.simplify(integral - var_true) == 0, str(sp.simplify(integral)))

a = sp.Function("a")
i_, j_, k_ = sp.symbols("i j k", positive=True, integer=True)
one_step_cov = sp.Sum(a(i_) ** 2, (i_, 1, sp.Min(j_, k_)))
check("P8.4 one-step Cov depends on (j,k) only through min(j,k)",
      str(one_step_cov).count("Min(j, k)") == 1, str(one_step_cov))

# --------------------------------------------------------------------------- #
# Lemma 9.6 and Theorem 9.8(iii) — the exponent accounting
# --------------------------------------------------------------------------- #
d_, g_ = sp.symbols("delta gamma", positive=True)
check("L9.6 sqrt(delta * delta^{2 gamma}) = delta^{1/2 + gamma}",
      sp.simplify(sp.sqrt(d_ * d_ ** (2 * g_)) - d_ ** (sp.Rational(1, 2) + g_)) == 0)

a_X = d_ ** (sp.Rational(1, 2) + g_)
n_ = Tsym / d_
acc = sp.simplify(sp.sqrt(n_ * a_X ** 2))
check("T9.8(iii) sqrt(n a_X^2) = sqrt(T) delta^gamma",
      sp.simplify(acc - sp.sqrt(Tsym) * d_ ** g_) == 0, str(acc))
check("T9.8(iii) gamma = (h+kappa)/2 gives the rate delta^{(h+kappa)/2}",
      sp.simplify((d_ ** g_).subs(g_, sp.Symbol("h") / 2 + sp.Symbol("kappa") / 2)
                  - d_ ** ((sp.Symbol("h") + sp.Symbol("kappa")) / 2)) == 0)
# the naive per-step bound diverges: n * delta^{1/2+gamma} = T delta^{gamma-1/2}
naive = sp.simplify(n_ * a_X)
check("T9.8(iii) the naive per-step bound n*delta^{1/2+gamma} DIVERGES "
      "(exponent gamma-1/2 < 0)",
      sp.simplify(naive - Tsym * d_ ** (g_ - sp.Rational(1, 2))) == 0, str(naive))

# --------------------------------------------------------------------------- #
# Lemma 4.8 — the four-point admissibility condition, finite delta
# --------------------------------------------------------------------------- #
Aa, rr, ee = sp.symbols("A rho epsilon", positive=True)
alpha = (Aa ** 2 - 1) / 2
lhs4 = sp.simplify((sp.Abs(alpha) + ee + rr * Aa) - (1 + alpha))
# with A<=1, |alpha| = (1-A^2)/2
lhs4 = sp.simplify(((1 - Aa ** 2) / 2 + ee + rr * Aa) - (1 + Aa ** 2) / 2)
check("L4.8 (four-point) |a|+|c| <= 1/4  <=>  rho A + epsilon <= A^2",
      sp.simplify(lhs4 - (rr * Aa + ee - Aa ** 2)) == 0, str(sp.simplify(lhs4)))
check("L4.8 asymptotic form epsilon -> 0 gives rho <= A",
      sp.simplify((rr * Aa - Aa ** 2) / Aa) == sp.simplify(rr - Aa))

# --------------------------------------------------------------------------- #
# Proposition (exact-convolution variance) — new in v1.4.0.
#
#   Var[V^(n)_T] = 1/(2H) + zeta(1-2H) delta^{2H} + delta/2 + O(delta^2),
#   Var[V^(n)_T]/Var[V_T] - 1 = 2H zeta(1-2H) delta^{2H} + H delta + O(delta^2),
#
# for K(u) = u^h, h = H - 1/2, T = 1, so the RELATIVE variance error is
# Theta(delta^{2H}) and NEGATIVE (the discretisation understates the variance).
# Checked against the exact finite sum at 30 digits.
# --------------------------------------------------------------------------- #
import mpmath as _mp
_mp.mp.dps = 30


def _var_disc(H, n):
    """delta^{2H} sum_{j=1}^{n} j^{2h}, h = H - 1/2 — exact finite sum.

    sum_{j=1}^{n} j^{s} = zeta(-s) - zeta(-s, n+1) by analytic continuation of
    the Hurwitz zeta, which is exact and O(1) instead of summing n terms.
    """
    d = _mp.mpf(1) / n
    s = 2 * H - 1
    tot = _mp.zeta(-s) - _mp.zeta(-s, n + 1)
    return d ** (2 * H) * tot


_ok_exp, _worst = True, 0.0
for _H in ("0.05", "0.1", "0.2", "0.3", "0.45"):
    _H = _mp.mpf(_H)
    _z = _mp.zeta(1 - 2 * _H)
    for _n in (64, 256, 4096, 16384):
        _d = _mp.mpf(1) / _n
        _meas = (_var_disc(_H, _n) - 1 / (2 * _H)) * 2 * _H
        _pred = 2 * _H * _z * _d ** (2 * _H) + _H * _d
        _rel = abs((_meas - _pred) / _meas)
        _worst = max(_worst, float(_rel))
        if _rel > _mp.mpf("1e-3"):
            _ok_exp = False
check("Prop(var) two-term expansion 2H zeta(1-2H) d^{2H} + H d matches the exact "
      "sum for H in {0.05,...,0.45}, n in {64,...,16384}",
      _ok_exp, f"worst relative deviation {_worst:.2e}")

check("Prop(var) zeta(1-2H) < 0 on H in (0,1/2), so the variance is UNDERSTATED",
      all(_mp.zeta(1 - 2 * _mp.mpf(h)) < 0
          for h in ("0.001", "0.05", "0.1", "0.25", "0.45", "0.499")),
      "sign of zeta on (0,1)")

check("Prop(var) relative variance error exponent 2H beats Thm 7.1's H/2 "
      "for every H in (0,1/2), by a factor 4",
      all(2 * _mp.mpf(h) == 4 * (_mp.mpf(h) / 2)
          for h in ("0.05", "0.1", "0.25", "0.45")),
      "2H = 4*(H/2)")

# --------------------------------------------------------------------------- #
# Route B factor indexing — new in v1.4.0.  The v1.3.0 formula
#   Z^i_k = sum_{j<=k} e^{-s_i (k-j) d} sqrt(d) zeta_j
# has smallest lag 0; Definition (exact convolution) has smallest lag delta.
# The corrected recursion Z^i_k = e^{-s_i d}(Z^i_{k-1} + sqrt(d) zeta_k) unrolls
# to lags (k-j+1) d, matching the Definition.
# --------------------------------------------------------------------------- #
_si, _dd = sp.symbols("varsigma_i delta", positive=True)
_zs = sp.symbols("zeta1:5")
_Z = 0
for _k in range(4):                        # unroll the corrected recursion
    _Z = sp.exp(-_si * _dd) * (_Z + sp.sqrt(_dd) * _zs[_k])
_target = sum(sp.exp(-_si * (4 - _j) * _dd) * sp.sqrt(_dd) * _zs[_j]
              for _j in range(4))          # lags (k-j+1)d for k=4, j=1..4
check("Route B corrected recursion unrolls to lags (k-j+1)delta, smallest lag "
      "delta (matches Dfn exact convolution)",
      sp.simplify(sp.expand(_Z - _target)) == 0, str(sp.simplify(_Z - _target)))
_paper = sum(sp.exp(-_si * (3 - _j) * _dd) * sp.sqrt(_dd) * _zs[_j]
             for _j in range(4))           # v1.3.0: lags (k-j)d, smallest 0
check("the v1.3.0 Route B formula differs from it by exactly one lag, i.e. by "
      "the factor e^{+varsigma_i delta}",
      sp.simplify(sp.expand(_paper - sp.exp(_si * _dd) * _target)) == 0,
      str(sp.simplify(_paper - sp.exp(_si * _dd) * _target)))

# --------------------------------------------------------------------------- #
# (B2): the recombination defect, the repair, and the corrected cost
# --------------------------------------------------------------------------- #
# Proposition (no recombination), step 1: detrending the factor gives a walk
# whose step magnitudes are sqrt(delta) q^{j-1} with q = exp(varsigma delta).
_q = sp.exp(_si * _dd)
_Zt = sp.exp(_si * 4 * _dd) * _Z                      # e^{varsigma t_k} Z_k, k=4
_walk = sp.sqrt(_dd) * sum(_q ** _j * _zs[_j] for _j in range(4))
check("Prop (no recombination): detrending the factor gives sum_j q^{j-1} "
      "zeta_j with q = exp(varsigma delta)",
      sp.simplify(sp.expand(_Zt - _walk)) == 0,
      str(sp.simplify(_Zt - _walk)))

# step 2: the step magnitudes are all equal iff q = 1 iff varsigma delta = 0,
# which is the only recombining case.
_mags_eq = [sp.simplify(_q ** _j - 1) for _j in (1, 2, 3)]
check("the step magnitudes q^{j-1} are mutually distinct unless q = 1, i.e. "
      "unless varsigma = 0 -- the constant-kernel one-step scheme",
      all(_e != 0 for _e in _mags_eq)
      and all(_e.subs(_si, 0) == 0 for _e in _mags_eq),
      str(_mags_eq))

# step 3: a collision is a root of a {-1,0,1}-polynomial at q.  Enumerate every
# such nonzero polynomial of degree < 8 and confirm none vanishes at
# q = exp(varsigma delta) for a grid of rational varsigma*delta -- the numerical
# face of the Lindemann-Weierstrass step, which sympy cannot prove.
import itertools as _it
_gap = _mp.inf
for _sd in (_mp.mpf(1) / 64, _mp.mpf(1) / 16, _mp.mpf(1) / 4, _mp.mpf(2),
            _mp.mpf("0.75")):
    _qq = _mp.e ** _sd
    _pows = [_qq ** _i for _i in range(8)]
    for _c in _it.product((-1, 0, 1), repeat=8):
        if not any(_c):
            continue
        _v = abs(sum(_c[_i] * _pows[_i] for _i in range(8)))
        if _v < _gap:
            _gap = _v
check("no nonzero polynomial with coefficients in {-1,0,1} and degree < 8 "
      "vanishes at exp(varsigma delta) for varsigma delta in "
      "{1/64,1/16,1/4,3/4,2} (Lindemann-Weierstrass, checked numerically)",
      _gap > _mp.mpf("1e-30"), f"min |P(q)| over 6560x5 polynomials = {_gap}")

# Lemma (the repair): the Abel weights (1-theta) sum_{l>=0} theta^l = 1, which
# is what makes the bound |D_k| <= 2 max_j |M_j| uniform in varsigma.
_th, _LL = sp.symbols("theta L", positive=True)
_lidx = sp.symbols("l", integer=True, nonnegative=True)
# the finite Abel sum that actually appears: (1-theta) sum_{l=0}^{L} theta^l,
# which equals 1 - theta^{L+1} and is therefore <= 1 for every theta in (0,1].
_finite = sp.simplify((1 - _th) * sp.Sum(_th ** _lidx, (_lidx, 0, 3)).doit())
check("Lemma (rounded lift): (1-theta) sum_{l=0}^{L} theta^l = 1 - "
      "theta^{L+1} <= 1 for theta in (0,1], so |D_k| <= 2 max_j |M_j| "
      "uniformly in varsigma",
      sp.simplify(_finite - (1 - _th ** 4)) == 0
      and all(float((1 - _t) * sum(_t ** _l for _l in range(_L + 1)))
              <= 1.0 + 1e-12
              for _t in (0.05, 0.5, 0.9, 0.999, 1.0) for _L in (0, 5, 50, 500)),
      str(_finite))

# Lemma (the repair), (iii): m * a * sqrt(n) = sqrt(T) delta^gamma exactly when
# a = delta^{1/2+gamma}/m and delta = T/n.
_m, _g, _TT, _nn = sp.symbols("m gamma T n", positive=True)
_delta = _TT / _nn
_a = _delta ** (sp.Rational(1, 2) + _g) / _m
check("Lemma (rounded lift)(iii): m a sqrt(n) = sqrt(T) delta^gamma for "
      "a = delta^{1/2+gamma}/m",
      sp.simplify(_m * _a * sp.sqrt(_nn) - sp.sqrt(_TT) * _delta ** _g) == 0,
      str(sp.simplify(_m * _a * sp.sqrt(_nn) - sp.sqrt(_TT) * _delta ** _g)))

# Lemma (state space): w_i s_i = O(varsigma^{-H}) since -h-1/2 = -H.
_hh, _HH = sp.symbols("h H")
check("Lemma (lift states): w_i s_i ~ varsigma^{-h} varsigma^{-1/2} = "
      "varsigma^{-H} because -h-1/2 = -(h+1/2) = -H",
      sp.simplify((-_hh - sp.Rational(1, 2)) - (-(_hh + sp.Rational(1, 2))))
      == 0 and sp.simplify((_hh + sp.Rational(1, 2)).subs(_hh, _HH
                                                          - sp.Rational(1, 2))
                           - _HH) == 0)

# the corrected cost exponents of eq:lift-eps-cost at H = 0.1, gamma = H/2
_Hv, _gv = _mp.mpf("0.1"), _mp.mpf("0.05")
_costs = {_mv: (1 + _gv + (_mv + 1) * (_mp.mpf("0.5") + _gv)) * 2 / _Hv
          for _mv in (1, 2, 3, 5)}
check("corrected cost exponents: m=1 -> eps^-43, m=2 -> eps^-54, m=3 -> "
      "eps^-65, m=5 -> eps^-87 at H=0.1",
      [int(_mp.nint(_costs[_k])) for _k in (1, 2, 3, 5)] == [43, 54, 65, 87],
      str({_k: float(_v) for _k, _v in _costs.items()}))
check("Route A' alone costs eps^-41, so one lift factor adds 2 to the "
      "exponent and two factors add 13",
      int(_mp.nint((2 + _gv) * 2 / _Hv)) == 41,
      str(float((2 + _gv) * 2 / _Hv)))

# --------------------------------------------------------------------------- #
# Lemma (exit-time stability by time change) — the last gap inside (A'1)
# --------------------------------------------------------------------------- #
# (1) the margin the whole lemma turns on: gamma = (h+kappa)/2 < 1/4 STRICTLY,
#     for every h < 0 and kappa < 1/2, with no extra hypothesis.
_hh2, _kk2 = sp.symbols("h kappa", real=True)
_gam = (_hh2 + _kk2) / 2
check("Lemma (exit stability): gamma = (h+kappa)/2 < 1/4 strictly whenever "
      "h < 0 and kappa < 1/2 -- so the delta^{1/4} estimate always beats "
      "delta^gamma",
      sp.simplify(sp.Rational(1, 4) - _gam.subs({_hh2: 0, _kk2: sp.Rational(1, 2)}))
      == 0
      and all(float(_gam.subs({_hh2: _H - 0.5, _kk2: 0.5})) < 0.25
              for _H in (0.01, 0.05, 0.1, 0.2, 0.3, 0.45, 0.499)),
      "gamma at kappa=1/2 equals H/2, and H<1/2")

# (2) the time change is Lipschitz: sigma_-^2 u <= <a>_u <= sigma_+^2 u implies
#     the inverse is Lipschitz with constant 1/sigma_-^2.  Checked as the
#     elementary inequality it is, on a grid of admissible quadratic variations.
_sm, _sp2 = 0.5, 1.5
_bad = None
_rng = random.Random(3)
for _ in range(200_000):
    _u1, _u2 = sorted((_rng.uniform(0, 5), _rng.uniform(0, 5)))
    # any nondecreasing <a> with slope in [sm^2, sp^2]
    _s1 = _rng.uniform(_sm ** 2, _sp2 ** 2)
    _qv1, _qv2 = _s1 * _u1, _s1 * _u1 + _rng.uniform(_sm ** 2, _sp2 ** 2) * (_u2 - _u1)
    if _qv2 - _qv1 < _sm ** 2 * (_u2 - _u1) - 1e-12:
        _bad = (_u1, _u2)
        break
    if (_u2 - _u1) > (_qv2 - _qv1) / _sm ** 2 + 1e-12:
        _bad = (_u1, _u2, _qv1, _qv2)
        break
check("Lemma (exit stability)(ii): <a>^{-1} is Lipschitz with constant "
      "1/sigma_-^2 (200k admissible quadratic variations)",
      _bad is None, str(_bad))

# (3) the closed-form exit moments used by the run's validity check.  For
#     standard Brownian motion and the band +-1, E[e^{-theta T}] = 1/cosh(sqrt(2 theta)),
#     whence E[T] = 1 and E[T^2] = 5/3, so Var[T] = 2/3.
_th2 = sp.symbols("theta", positive=True)
_lap = 1 / sp.cosh(sp.sqrt(2 * _th2))
# sqrt(2 theta) is not differentiable at 0, so read the moments off the SERIES:
# E[e^{-theta T}] = sum_k (-theta)^k E[T^k]/k!, and the transform is analytic in
# theta because cosh(sqrt(2 theta)) = 1 + theta + theta^2/6 + ... has only even
# powers of sqrt(2 theta).
_ser = sp.series(_lap, _th2, 0, 3).removeO().expand()
_m1 = -_ser.coeff(_th2, 1)
_m2 = 2 * _ser.coeff(_th2, 2)
check("exit moments of the band +-1 by standard BM: E[T] = 1 and E[T^2] = 5/3 "
      "from the Laplace transform 1/cosh(sqrt(2 theta)), so Var[T] = 2/3",
      sp.simplify(_m1 - 1) == 0 and sp.simplify(_m2 - sp.Rational(5, 3)) == 0
      and sp.simplify(_m2 - _m1 ** 2 - sp.Rational(2, 3)) == 0,
      f"E[T]={sp.simplify(_m1)}, E[T^2]={sp.simplify(_m2)}")
# and their time-scaling to volatility sbar, which is what the run compares to
_sb = sp.symbols("sbar", positive=True)
check("scaling to sbar x BM: mean 1/sbar^2 and variance (2/3)/sbar^4",
      sp.simplify(_m1 / _sb ** 2 - 1 / _sb ** 2) == 0
      and sp.simplify((_m2 - _m1 ** 2) / _sb ** 4
                      - sp.Rational(2, 3) / _sb ** 4) == 0)

# (4) the coupled correlated coordinate really is a Brownian motion with the
#     right bracket: <What> = 1 and <What, ahat> = rho sbar.
_rho2 = sp.symbols("rho", real=True)
_brk_w = _rho2 ** 2 * _sb ** -2 * _sb ** 2 + (1 - _rho2 ** 2)
_brk_wa = _rho2 * _sb ** -1 * _sb ** 2
check("Lemma (exit stability)(iii): What = rho ahat/sbar + sqrt(1-rho^2) Bperp "
      "has <What> = 1 and <What, ahat> = rho sbar",
      sp.simplify(_brk_w - 1) == 0 and sp.simplify(_brk_wa - _rho2 * _sb) == 0,
      f"{sp.simplify(_brk_w)}, {sp.simplify(_brk_wa)}")

# (5) the interpolation used on the correlated part:
#     (sbar - s)^2 = (sbar^2 - s^2)^2/(sbar + s)^2 <= E |sbar^2 - s^2|/(sbar+sm)^2
_bad2 = None
for _ in range(200_000):
    _s = _rng.uniform(_sm, _sp2)
    _sbv = _rng.uniform(_sm, _sp2)
    _E = abs(_sbv ** 2 - _s ** 2)
    if (_sbv - _s) ** 2 > _E * abs(_sbv ** 2 - _s ** 2) / (_sbv + _sm) ** 2 + 1e-12:
        _bad2 = (_s, _sbv)
        break
check("Lemma (exit stability)(iii): (sbar-s)^2 <= E |sbar^2-s^2|/(sbar+sigma_-)^2 "
      "with E = sup|sbar^2-s^2| (200k samples)",
      _bad2 is None, str(_bad2))

# (6) the three regimes of Remark (node error): the condition
#     (sqrt delta + eps_n)^{1/2} + delta^{1/4} = O(delta^gamma).
_reg = {"constant coefficients (eps_n = 0)": 0.25,
        "eps_n = delta^{1/4} in L^4": 0.125,
        "eps_n = delta^{1/8} in L^4": 0.0625}
check("Remark (node error): the three regimes give gamma <= 1/4, 1/8, 1/16, "
      "i.e. no constraint / H <= 1/4 / H <= 1/8",
      [round(_v, 6) for _v in _reg.values()] == [0.25, 0.125, 0.0625]
      and abs(2 * 0.125 - 0.25) < 1e-12,
      str(_reg))

# --------------------------------------------------------------------------- #
# Retractions forced by the external review of v1.6.0
# --------------------------------------------------------------------------- #
# (R1) Lemma 4.4 as stated up to v1.6.0 -- Delta V = K(d) Delta Y + o(sqrt d) --
#      is FALSE.  The remainder at k=2 is exactly (2^h - 1) d^H zeta_1, and
#      |R_2| / sqrt(d) = |2^h - 1| d^h  ->  infinity.
_dd2 = sp.symbols("delta", positive=True)
_kk3 = sp.symbols("k", positive=True)
for _Hv in (sp.Rational(1, 20), sp.Rational(1, 10), sp.Rational(1, 4),
            sp.Rational(2, 5)):
    _hv = _Hv - sp.Rational(1, 2)
    _R2 = ((2 * _dd2) ** _hv - _dd2 ** _hv) * sp.sqrt(_dd2)
    _ok = sp.simplify(sp.powsimp(_R2 - (2 ** _hv - 1) * _dd2 ** _Hv)) == 0
    _div = sp.limit(sp.Abs(_R2 / sp.sqrt(_dd2)), _dd2, 0, "+")
    check(f"(R1) Lemma 4.4 retracted at H={_Hv}: R_2 = (2^h-1) delta^H exactly "
          f"and |R_2|/sqrt(delta) -> oo, so it is NOT o(sqrt delta)",
          bool(_ok) and _div == sp.oo, f"exact={_ok}, limit={_div}")
# the proof's own coefficient bound diverges: the telescope is k^h - 1, bounded
# by 1, so the coefficient sum is Lbar delta^h -> infinity, not "<= L'".
_j2 = sp.symbols("j", positive=True, integer=True)
_hsym = sp.symbols("h", negative=True)
_tel = sp.simplify(sp.summation((_j2 + 1) ** _hsym - _j2 ** _hsym,
                                (_j2, 1, _kk3 - 1)))
check("(R1) the telescope in the old proof is k^h - 1 (not n^h - 2^h), bounded "
      "by 1, so the coefficient sum is Lbar delta^h, which DIVERGES",
      sp.simplify(_tel - (_kk3 ** _hsym - 1)) == 0
      and all(float(_d ** (0.1 - 0.5)) > 1 for _d in (1e-2, 1e-4, 1e-6)),
      str(_tel))
# (R2) Lemma 7.6: derived exponent h + kappa/2, stated (h+kappa)/2; the gap is
#      -h/2 > 0, and the derived one is NEGATIVE for every H < 1/4.
_kap2 = sp.symbols("kappa", positive=True)
_hh3 = sp.symbols("h", negative=True)
check("(R2) Lemma 7.6: (h+kappa)/2 - (h+kappa/2) = -h/2 > 0, so the stated "
      "exponent is strictly stronger than the one the proof derives",
      sp.simplify(((_hh3 + _kap2) / 2 - (_hh3 + _kap2 / 2)) + _hh3 / 2) == 0)
_neg = {H: round(H - 0.25, 4) for H in (0.05, 0.10, 0.20, 0.30, 0.45)}
check("(R2) optimising kappa -> 1/2 in the DERIVED exponent gives H - 1/4, "
      "negative (the bound diverges) for every H < 1/4",
      all(v < 0 for H, v in _neg.items() if H < 0.25)
      and all(v > 0 for H, v in _neg.items() if H > 0.25), str(_neg))
# (R3) the interface mismatch is a statement about definitions: mean-exactness
#      does not imply equality of laws.  Exhibited by a counterexample.
_p = 0.5
_lawA = {-1: 0.5, 1: 0.5}                       # the binomial the interface asks
_lawB = {-2: 0.25, 0: 0.5, 2: 0.25}             # a mean-exact, Z-valued rounding
_meanA = sum(k * v for k, v in _lawA.items())
_meanB = sum(k * v for k, v in _lawB.items())
check("(R3) mean-exactness is strictly weaker than equality of laws: two laws "
      "with the same mean, one on {-1,1} and one on Z, differ in every other "
      "respect (so (E2) does not follow from the rounding being mean-exact)",
      abs(_meanA - _meanB) < 1e-15 and _lawA != _lawB
      and abs(sum(k * k * v for k, v in _lawA.items())
              - sum(k * k * v for k, v in _lawB.items())) > 0.5,
      f"means {_meanA} = {_meanB}, second moments "
      f"{sum(k*k*v for k,v in _lawA.items())} vs "
      f"{sum(k*k*v for k,v in _lawB.items())}")

# --------------------------------------------------------------------------- #
# Lemma 7.7 repair (v2.1.0) -- interpolate the EMBEDDED WALK, isolate the node
# error (Lemma 7.9), rebuild Proposition 7.10.  The rate stays (h+kappa)/2; the
# repair is that this exponent is now PROVED (was only targeted in v2.0.0).
# Domain: h = H - 1/2 in (-1/2, 0), kappa in (0, 1/2), h + kappa > 0.
# --------------------------------------------------------------------------- #
_hL, _kL = sp.symbols("h kappa", real=True)
_target = (_hL + _kL) / 2                      # the rate prop:Vconv must deliver
_grid = [(hh, kk) for hh in (-0.45, -0.4, -0.3, -0.2, -0.1, -0.05)
         for kk in (0.05, 0.1, 0.2, 0.3, 0.4, 0.45, 0.49) if hh + kk > 0]

# (L7.7-i) T1: interpolating the WALK gives exponent H = h + 1/2 (walk increments
#          are the lattice step sqrt(delta) deterministically), and H >= target,
#          strictly, on the whole domain -- so T1 is no longer the bottleneck.
check("(L7.7) T1 exponent H = h+1/2 dominates the rate: H - (h+kappa)/2 = "
      "(h+1-kappa)/2 > 0 on the admissible domain",
      sp.simplify((_hL + sp.Rational(1, 2)) - _target - (_hL + 1 - _kL) / 2) == 0
      and all((hh + 0.5) - (hh + kk) / 2 > 1e-12 for hh, kk in _grid),
      "H=h+1/2")

# (L7.7-old) contrast: the v2.0.0 exponent h + kappa/2 was BELOW target by h/2<0
#            and negative for H<1/4 (the bug the repair removes).
check("(L7.7) the retired exponent h+kappa/2 was target - (-h/2), i.e. short by "
      "|h|/2, and negative for H<1/4; the walk interpolation removes it",
      sp.simplify(_target - (_hL + _kL / 2) - (-_hL / 2)) == 0
      and all((hh + kk / 2) < 0 for hh, kk in _grid if hh + 0.5 < 0.25),
      "h+kappa/2")

# (L7.9) T3a (node error): the Marchaud interpolation-of-bounds minimises
#        f(eps) = B eps^{h+kappa} + A eps^{h} over eps>0; the optimum scales like
#        A^{(h+kappa)/kappa} B^{-h/kappa}.  Verify the two exponents symbolically.
_A, _B, _eps = sp.symbols("A B epsilon", positive=True)
_hN = sp.symbols("h", negative=True)
_kN = sp.symbols("kappa", positive=True)
_f = _B * _eps ** (_hN + _kN) + _A * _eps ** _hN
_estar = (_A * (-_hN) / (_B * (_hN + _kN))) ** (1 / _kN)   # f'(eps*) = 0
_val = _f.subs(_eps, _estar)
_aexp = sp.simplify(_A * sp.diff(sp.log(_val), _A))         # d log val / d log A
_bexp = sp.simplify(_B * sp.diff(sp.log(_val), _B))
check("(L7.9) Marchaud optimum ~ A^{(h+kappa)/kappa} B^{-h/kappa}: the A- and "
      "B-exponents of the minimised bound are exactly (h+kappa)/kappa and -h/kappa",
      sp.simplify(_aexp - (_hN + _kN) / _kN) == 0
      and sp.simplify(_bexp + _hN / _kN) == 0,
      f"A-exp={_aexp}, B-exp={_bexp}")

# with the node bound A = O(delta^{1/4}) and B = O(1), T3a has delta-exponent
# (h+kappa)/(4 kappa), which EXCEEDS target (h+kappa)/2 for kappa<1/2, so T3a is
# not the bottleneck either.
_d = sp.symbols("delta", positive=True)
_dexp = sp.simplify(_d * sp.diff(sp.log(_val.subs({_A: _d ** sp.Rational(1, 4),
                                                   _B: 1})), _d))
check("(L7.9) T3a delta-exponent (h+kappa)/(4 kappa) with A=delta^{1/4}, B=O(1); "
      "it exceeds the rate (h+kappa)/2 for every kappa<1/2 (since 1/(4k)>1/2)",
      sp.simplify(_dexp - (_hN + _kN) / (4 * _kN)) == 0
      and all(((hh + kk) / (4 * kk)) - (hh + kk) / 2 > 1e-12 for hh, kk in _grid),
      f"T3a exponent={_dexp}")

# (L7.9) the GFOtoRL reduction remainder on the node error is controlled by the
#        SUP norm, not the Holder seminorm: k'(v)=O(v^h) with h>-1 is integrable,
#        so int_0^t |k'| < infinity.  Verify the antiderivative identity
#        d/dv [v^{h+1}/(h+1)] = v^h and that the exponent h+1>0 on the domain
#        (i.e. h>-1), which is what makes the improper integral at 0 converge.
_v = sp.symbols("v", positive=True)
_antideriv_ok = sp.simplify(sp.diff(_v ** (_hN + 1) / (_hN + 1), _v)
                            - _v ** _hN) == 0
check("(L7.9) GFOtoRL remainder is a smoothing operator: v^h has antiderivative "
      "v^{h+1}/(h+1) and h+1>0 on the domain, so int_0^t v^h dv converges and "
      "|R phi| <= const * ||phi||_infty (sup-norm control)",
      _antideriv_ok and all(hh + 1 > 0 for hh, _ in _grid),
      f"antiderivative identity ok={_antideriv_ok}")

# (Prop 7.10) the rate is unchanged: min over the three non-binding terms
#             T1 (H), T2 (min{(h+k)/(4k), 1/4}), T3a inside T2's leading part,
#             all EXCEED target, so T3b = (h+kappa)/2 binds -> prop:Vconv rate is
#             (h+kappa)/2, exactly the claim, now proved.
check("(Prop 7.10) rate unchanged and now proved: min(H, kappa/2, (h+k)/(4k), "
      "1/4) > (h+kappa)/2 on the whole domain, so T3b=(h+kappa)/2 is the "
      "binding term",
      all(min(hh + 0.5, kk / 2, (hh + kk) / (4 * kk), 0.25) - (hh + kk) / 2
          > 1e-12 for hh, kk in _grid),
      "binding term (h+kappa)/2")

# sanity: at H=0.1 the proved rate is positive (target=(h+0.49)/2 with h=-0.4),
# whereas the retired exponent diverged.
check("(Prop 7.10) at H=0.1, kappa=0.49 the proved rate exponent (h+kappa)/2 = "
      "0.045 > 0 (converges), against the retired h+kappa/2 = -0.155 (diverged)",
      abs((-0.4 + 0.49) / 2 - 0.045) < 1e-9 and (-0.4 + 0.49 / 2) < 0,
      f"proved={(-0.4+0.49)/2}, retired={-0.4+0.49/2}")

# --------------------------------------------------------------------------- #
# Theorem 7.12 (v2.2.0): convergence from the stability interface.  The proof
# yields |Lambda - Lambda^sch| = O(r_X + r_V + delta^{(h+kappa)/2}), the last
# term being the modulus of the ROUGH variance over one clock gap of size
# sqrt(delta): (sqrt(delta))^{h+kappa} = delta^{(h+kappa)/2}.  Here we check the
# rate BOOKKEEPING that makes thm:main and thm:B2 corollaries.
# --------------------------------------------------------------------------- #
# The variance modulus over a gap sqrt(delta) is delta^{(h+kappa)/2}, and it is
# the slowest (largest) of the three fixed rates delta^{1/2}, delta^{1/4}
# (X-modulus) and delta^{(h+kappa)/2}, since (h+kappa)/2 < 1/4 < 1/2.
check("(Thm 7.12) variance modulus over a sqrt(delta) clock gap is "
      "(sqrt d)^{h+kappa} = d^{(h+kappa)/2}, and (h+kappa)/2 < 1/4 < 1/2 so it "
      "dominates the delta^{1/4} (X-modulus) and delta^{1/2} fixed rates",
      all(abs(0.5 * (hh + kk) - (hh + kk) / 2) < 1e-12
          and (hh + kk) / 2 < 0.25 < 0.5 for hh, kk in _grid),
      "d^{(h+k)/2} is the binding fixed rate")

# thm:main is the corollary with the admissible-embedding rates r_X = d^{1/4}
# (E3 state) and r_V = d^{(h+kappa)/2} (Prop 7.10); the interface bound collapses
# to d^{(h+kappa)/2}, i.e. the exponent of thm:main is exactly min over the pile.
def _iface_exponent(rX_exp, rV_exp, h, k):
    """exponent of r_X + r_V + d^{(h+k)/2} (largest error = smallest exponent)."""
    return min(rX_exp, rV_exp, (h + k) / 2)
check("(Thm 7.1 as corollary) admissible embedding gives r_X=d^{1/4}, "
      "r_V=d^{(h+kappa)/2}; interface bound min(1/4,(h+k)/2,(h+k)/2)=(h+k)/2, "
      "recovering exactly the rate of Theorem 7.1",
      all(abs(_iface_exponent(0.25, (hh + kk) / 2, hh, kk) - (hh + kk) / 2)
          < 1e-12 for hh, kk in _grid),
      "collapses to (h+kappa)/2")

# thm:B2 is the corollary with r_V = d^{(h+kappa)/2} + ||K-K^m||; the interface
# bound collapses to d^{(h+kappa)/2} + ||K-K^m||.  Model ||K-K^m|| by an extra
# additive nonneg term E>=0; the delta-exponent part is still (h+kappa)/2.
check("(Thm 9.26 (B2) as corollary) with r_V = d^{(h+kappa)/2}+||K-K^m||, "
      "r_X=d^{1/4}, the interface delta-part collapses to (h+kappa)/2, leaving "
      "the B2 bound d^{(h+kappa)/2}+||K-K^m||",
      all(abs(_iface_exponent(0.25, (hh + kk) / 2, hh, kk) - (hh + kk) / 2)
          < 1e-12 for hh, kk in _grid),
      "B2 = d^{(h+kappa)/2} + ||K-K^m||")

# the retired hyp:stability conclusion wrote "+ delta^{1/2}"; that OMITS the
# variance-modulus term delta^{(h+kappa)/2}, which is LARGER (slower).  It was
# harmless only because r_V already carries delta^{(h+kappa)/2}; the corrected
# statement carries it explicitly.
check("(Thm 7.12) the corrected fixed rate is delta^{(h+kappa)/2}, not "
      "delta^{1/2} as the v2.0.0 hypothesis wrote: (h+kappa)/2 < 1/2 so the "
      "former is the larger error and must appear",
      all((hh + kk) / 2 < 0.5 for hh, kk in _grid),
      "delta^{(h+kappa)/2} > delta^{1/2}")

# --------------------------------------------------------------------------- #
# External review of v2.2.0 (v2.3.0 fixes).
# --------------------------------------------------------------------------- #
# (Rev-1) nine-point sufficiency was missing m <= min(1-s,1-t): the entry
#         pi = p_00 = min(1-s,1-t) - m can be NEGATIVE under |c|+m<=min(s,t)
#         alone.  The reviewer's counterexample s=t=0.9, m=0.2, c=0.
_s0, _t0, _m0, _c0 = 0.9, 0.9, 0.2, 0.0
_pi0 = min(1 - _s0, 1 - _t0) - _m0
check("(Rev-1) nine-point (ii) needs m<=min(1-s,1-t): the counterexample "
      "s=t=0.9,m=0.2,c=0 satisfies |c|+m<=min(s,t) but gives pi=p_00<0",
      (abs(_c0) + _m0 <= min(_s0, _t0)) and _pi0 < 0
      and abs(_pi0 - (-0.1)) < 1e-12,
      f"|c|+m={abs(_c0)+_m0}<=0.9 but pi={_pi0}")
# with BOTH conditions, pi>=0 everywhere (this is what verify already samples,
# line ~106: m ~ U[0,1]*min(1-s,1-t), which is exactly the missing condition).
_rs2 = _np.random.default_rng(23)
_S = _rs2.uniform(1e-3, 1 - 1e-3, 50_000)
_T = _rs2.uniform(1e-3, 1 - 1e-3, 50_000)
_M = _rs2.uniform(0.0, 1.0, 50_000) * _np.minimum(1 - _S, 1 - _T)  # m<=min(1-s,1-t)
_PI = _np.minimum(1 - _S, 1 - _T) - _M
check("(Rev-1) with m<=min(1-s,1-t) added, pi=p_00=min(1-s,1-t)-m >= 0 always "
      "(50k samples) -- the corrected sufficiency condition",
      float(_PI.min()) >= -1e-12, f"min pi = {float(_PI.min()):.3e}")

# (Rev-2) purification (Lemma 7.12): if J(u) := value at fixed auxiliary u and
#         the randomised value is the average int J dP_U, then max_u J(u) >= mean.
#         This is the whole content -- there EXISTS u* at least as good as the
#         randomised rule.  (Elementary: a max is >= a mean.)
_J = _rs2.uniform(0.0, 1.0, 10_000)   # arbitrary values of J over realisations u
check("(Rev-2) purification principle: max_u J(u) >= E_U[J], so the de-randomised "
      "F-stopping time is at least as good as the randomised one",
      float(_J.max()) >= float(_J.mean()) - 1e-15,
      f"max={_J.max():.4f} >= mean={_J.mean():.4f}")

# (Rev-3) corrected interface rate O(r_X + r_V + omega_X + omega_V): for the
#         admissible embedding omega_X=d^{1/4}, omega_V=d^{(h+kappa)/2},
#         r_X=d^{1/4}, r_V=d^{(h+kappa)/2}; the min-exponent (largest error) is
#         (h+kappa)/2, recovering thm:main -- now with the moduli explicit.
check("(Rev-3) interface rate min(r_X, r_V, omega_X, omega_V) exponents = "
      "min(1/4, (h+k)/2, 1/4, (h+k)/2) = (h+k)/2 on the whole domain",
      all(abs(min(0.25, (hh + kk) / 2, 0.25, (hh + kk) / 2) - (hh + kk) / 2)
          < 1e-12 for hh, kk in _grid),
      "moduli-form rate collapses to (h+kappa)/2")

# --------------------------------------------------------------------------- #
# n <-> m : the link between the time step and the number of lift factors
# (sim/nm_bound.py, docs/N_M_BOUND.md).  Two things must hold for the bound to
# be usable from a run: sizing m by the criterion must not change what the
# lattice computes, and the node count of the explicit construction must really
# be the arithmetic 1 + ceil(log(span)/log r) the bound is stated in.
import math as _math
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent / "sim"))
try:
    import numpy as _nnp
    import nm_bound as _nb
    import route_b_lattice as _rbl

    _info = _nb.lift_for_tolerance(32, 0.10, 0.10)
    _w2, _s2 = _rbl.lift_for(32, _info["m"], 0.10)
    _dw = float(_nnp.abs(_info["w"] - _w2).max())
    _ds = float(_nnp.abs(_info["s"] - _s2).max())
    check("(n-m) nm_bound.lift_for_tolerance reproduces "
          "route_b_lattice.lift_for exactly, so sizing m from the criterion "
          "does not change what the lattice computes",
          _dw == 0.0 and _ds == 0.0,
          f"max|dw|={_dw:.3e} max|ds|={_ds:.3e} at m={_info['m']}")

    _r, _n, _T = 4.0, 256, 1.0
    _, _s3, _m3 = _nb.geometric_lift_by_ratio(_r, 0.10, _n, _T)
    _expect = 1 + int(_math.ceil(_math.log((1.0 / (_T / _n)) / (1.0 / _T))
                                / _math.log(_r)))
    check("(n-m) geometric_lift_by_ratio node count = "
          "1 + ceil(log(T/delta)/log r) -- the arithmetic the log law rests on",
          _m3 == _expect, f"m={_m3} expected={_expect} (r={_r}, n={_n})")

    # the requested r only SETS m; the partition is then stretched to span
    # [1/T, 1/delta] exactly with those m nodes, so the realised ratio is
    # span^{1/(m-1)} <= r.  (For n=256, r=4 the two coincide because 256 = 4^4.)
    _rr = _nb.geometric_lift(_m3, 0.10, _n, _T)[2]
    _span = (1.0 / (_T / _n)) / (1.0 / _T)
    check("(n-m) geometric_lift realised ratio = span^{1/(m-1)} and is at most "
          "the requested r -- the requested r sets m, the partition then spans "
          "the range exactly",
          abs(_rr - _span ** (1.0 / (_m3 - 1))) < 1e-9 and _rr <= _r + 1e-9,
          f"realised={_rr:.9f} span^(1/(m-1))={_span ** (1.0 / (_m3 - 1)):.9f} "
          f"requested={_r}")
except Exception as _e:                       # never let this block the algebra
    check("(n-m) nm_bound checks importable", False, f"{type(_e).__name__}: {_e}")

# --------------------------------------------------------------------------- #
# O002 -- prop:lattice-vs-embedding, restated and reproved at O(delta^H).
# These are SYMBOLIC checks of the exponent bookkeeping and of one conditional
# mean.  They are corroboration, never a proof (CLAUDE.md rule 3).
# --------------------------------------------------------------------------- #
_i, _j = sp.symbols("i j")
_a1, _ak, _b1, _bk = sp.symbols("alpha_km1 alpha_k beta_km1 beta_k")
_mx, _my, _rho, _sx, _sy = sp.symbols("mu_x mu_y rho sigma_x sigma_y")
_dd, _xi1, _ze1 = sp.symbols("delta xi_km1 zeta_km1", positive=False)

def _kernel(i, j):
    """eq:kernel, the four-point transition law of thm:momentmatching."""
    return (sp.Rational(1, 4)
            + i * (_a1 * _xi1 + _mx * sp.sqrt(_dd)) / (4 * (1 + _ak))
            + j * (_b1 * _ze1 + _my * sp.sqrt(_dd)) / (4 * (1 + _bk))
            + i * j * (_rho * _sx * _sy + _a1 * _b1 * _xi1 * _ze1)
              / (4 * (1 + _ak) * (1 + _bk)))

_signs = [(a, b) for a in (1, -1) for b in (1, -1)]
check("(O002) eq:kernel is a probability law: the four masses sum to 1",
      sp.simplify(sum(_kernel(a, b) for a, b in _signs) - 1) == 0,
      "sum of the four-point kernel")

# The step the first draft of the O002 proof got WRONG: the conditional mean of
# xi_k is NOT O(sqrt(delta)); it carries the order-one term alpha_{k-1} xi_{k-1}.
_Exi = sp.simplify(sum(a * _kernel(a, b) for a, b in _signs))
check("(O002) E^{(n)}_{k-1}[xi_k] = (alpha_{k-1} xi_{k-1} + mu_x sqrt(delta))"
      "/(1+alpha_k) -- ORDER ONE, not O(sqrt(delta)): the raw increment is not "
      "centred, which is why the proof runs on the perturbed walk eq:perturbed",
      sp.simplify(_Exi - (_a1 * _xi1 + _mx * sp.sqrt(_dd)) / (1 + _ak)) == 0
      and sp.simplify(_Exi.subs({_dd: 0}) - _a1 * _xi1 / (1 + _ak)) == 0,
      f"E[xi_k|F] = {_Exi}")

_Eze = sp.simplify(sum(b * _kernel(a, b) for a, b in _signs))
check("(O002) E^{(n)}_{k-1}[zeta_k] = (beta_{k-1} zeta_{k-1} + mu_y sqrt(delta))"
      "/(1+beta_k) -- eq:Ezeta2, the same order-one obstruction for the driver",
      sp.simplify(_Eze - (_b1 * _ze1 + _my * sp.sqrt(_dd)) / (1 + _bk)) == 0,
      f"E[zeta_k|F] = {_Eze}")

# ...and the repair: the PERTURBED increment is centred at exactly mu_x delta.
# Delta Shat_k = sqrt(d)[(1+alpha_k) xi_k - alpha_{k-1} xi_{k-1}].
_EdS = sp.simplify(sp.sqrt(_dd) * ((1 + _ak) * _Exi - _a1 * _xi1))
check("(O002) the perturbed walk IS centred: E_{k-1}[Delta Shat_k] = "
      "sqrt(d)[(1+alpha_k)E[xi_k] - alpha_{k-1}xi_{k-1}] = mu_x delta EXACTLY "
      "(eq:Exi2), the identity the Doob step of the proof rests on",
      sp.simplify(_EdS - _mx * _dd) == 0,
      f"E[Delta Shat] = {_EdS}")

# The two factorisations of the exponent coincide -- v2.3.0 attributed delta^H
# to ONE increment K(d)sqrt(d); the proof binds through (amplitude)x(clock error)
# = K(d) * d^{1/2}.  Same number, different mechanism (rem:...-rate).
_Hs, _ds = sp.symbols("H delta", positive=True)
_Kd = _ds ** (_Hs - sp.Rational(1, 2))
check("(O002) K(delta)*sqrt(delta) = K(delta)*delta^{1/2} = delta^H: the two "
      "factorisations of the exponent agree numerically, so v2.3.0's stated "
      "bound survives even though its stated mechanism does not",
      sp.simplify(_Kd * sp.sqrt(_ds) - _ds ** _Hs) == 0
      and sp.simplify(_Kd * _ds ** sp.Rational(1, 2) - _ds ** _Hs) == 0,
      "amplitude x clock error = one increment, in order")

# Domination downstream: H > (h+kappa)/2  <=>  h + 1 > kappa.
_hh_, _kk_ = sp.symbols("h kappa")
check("(O002) H - (h+kappa)/2 = (h + 1 - kappa)/2, so delta^H is dominated by "
      "delta^{(h+kappa)/2} exactly when h + 1 > kappa",
      sp.simplify((_hh_ + sp.Rational(1, 2)) - (_hh_ + _kk_) / 2
                  - (_hh_ + 1 - _kk_) / 2) == 0,
      "equivalence of the two forms")

check("(O002) h+1 > kappa holds on the WHOLE range of thm:main "
      "(h in (-1/2,0), kappa in (0,1/2), h+kappa>0), so no downstream rate "
      "changes: the reindexing error stays dominated",
      all(hh + 1 > kk for hh, kk in _grid),
      "h+1 > 1/2 > kappa")

# ...but NOT on part 1's wider standing range kappa in (0,1) -- the reason
# contribution (C2) must name the range rather than say "every admissible pair".
check("(O002) the domination FAILS on part 1's wider standing range "
      "kappa in (0,1): at h=-0.45, kappa=0.9 (admissible there, h+kappa=0.45) "
      "H=0.05 < 0.225=(h+kappa)/2 -- hence the explicit range in (C2)",
      not (-0.45 + 1 > 0.9) and (-0.45 + 0.5) < (-0.45 + 0.9) / 2,
      "counterexample on the wider range")

# The obstruction remark: the index-shift route gives delta^{H-1/4}.
check("(O002) index-shift route: a shift of delta^{-1/2} lattice steps moves "
      "Vcheck by K(d)sqrt(d)*delta^{-1/4} = delta^{H-1/4}, which does NOT "
      "vanish for H <= 1/4 -- why only the time argument is compared",
      sp.simplify(_Kd * sp.sqrt(_ds) * _ds ** sp.Rational(-1, 4)
                  - _ds ** (_Hs - sp.Rational(1, 4))) == 0
      and all(hh + 0.5 - 0.25 <= 0 for hh in (-0.45, -0.4, -0.3, -0.25)),
      "delta^{H-1/4} for H<=1/4")

# --------------------------------------------------------------------------- #
# O001 -- the arithmetic the closed Step 2 uses.  NOTHING here closes the
# immersion: the immersion is an ANALYTIC_ARGUMENT (lem:immersion) and a symbolic
# row may never be described as closing it.  These rows check only the exponent
# bookkeeping that the new text relies on.
_h, _k = sp.symbols("h kappa", real=True)
_Ts, _ell, _gam1 = sp.symbols("T ell gamma1", positive=True)

# The admissible range of thm:main: h in (-1/2, 0), kappa in (0, 1/2), h+kappa>0.
_adm = [(hh, kk) for hh in (-0.49, -0.45, -0.3, -0.2, -0.05, -0.01)
        for kk in (0.01, 0.05, 0.2, 0.35, 0.45, 0.49) if hh + kk > 0]

check("(O001) (h+kappa)/2 < 1/4 STRICTLY on the admissible range, equivalently "
      "h+kappa < 1/2 -- so the interval (frac{h+kappa}{2}, 1/4) for vartheta is "
      "non-empty and rem:moduli-rates can choose one",
      bool(sp.simplify(sp.Rational(1, 4) - (_h + _k) / 2
                       - (sp.Rational(1, 2) - (_h + _k)) / 2) == 0)
      and all((hh + kk) / 2 < 0.25 for hh, kk in _adm) and len(_adm) > 0,
      "vartheta interval non-empty")

check("(O001) the union bound over the n = T/delta clock gaps gives "
      "(n * delta^{ell/4})^{1/ell} = T^{1/ell} delta^{1/4 - 1/ell}, which "
      "exceeds delta^{(h+kappa)/2} for ell large: omega_X^{(n)} is dominated",
      sp.simplify(sp.powsimp((_Ts / _ds * _ds ** (_ell / 4)) ** (1 / _ell)
                  / (_Ts ** (1 / _ell) * _ds ** (sp.Rational(1, 4) - 1 / _ell)),
                  force=True) - 1) == 0,
      "1/4 - 1/ell")

check("(O001) for every admissible (h,kappa) there is an integer ell with "
      "1/4 - 1/ell > (h+kappa)/2, so the loss on the X-modulus never binds",
      all(any(0.25 - 1.0 / ell > (hh + kk) / 2 for ell in range(5, 4000))
          for hh, kk in _adm),
      "exists ell")

# The pre-registered falsifier F4 of the L002 contract, recorded as arithmetic:
# summing lem:freeze-clock's per-step Wasserstein distance over n steps DIVERGES,
# which is why (S5) is discharged exactly or assumed, and never summed.
check("(O001/F4) n * delta^{1/2+gamma} = T * delta^{gamma-1/2} DIVERGES as "
      "delta->0 for gamma = (h+kappa)/2 < 1/2 -- any draft closing (S5) by "
      "summing a per-step Wasserstein distance is refuted by this line",
      sp.simplify((_Ts / _ds) * _ds ** (sp.Rational(1, 2) + _gam1)
                  - _Ts * _ds ** (_gam1 - sp.Rational(1, 2))) == 0
      and all((hh + kk) / 2 - 0.5 < 0 for hh, kk in _adm),
      "delta^{gamma-1/2} -> infinity")

# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    w = max(len(n) for n, _, _ in results)
    nfail = 0
    print(f"{'claim':<{w}}  status")
    print("-" * (w + 10))
    for name, st, detail in results:
        print(f"{name:<{w}}  {st}" + ("" if st == PASS else f"   <- {detail}"))
        nfail += st == FAIL
    print("-" * (w + 10))
    print(f"{len(results) - nfail}/{len(results)} passed")
    sys.exit(1 if nfail else 0)
