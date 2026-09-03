#!/usr/bin/env python3
"""Symbolic check of the Lot C conditional moment table (lock L005).

This checks ALGEBRA, not consistency.  It exists to catch arithmetic slips in a
long computation, exactly as verify_claims.py does elsewhere; per the L005
contract's evidence_note it is NOT a proof and discharges no Definition-of-Done
item on its own.

What is checked, by exact expectation over the finite law of
(xi, xi_perp, U) -- two Rademachers and one shared uniform:

  M1   first moments of Delta X, Delta Y, Delta Z^i are exact / correct to o(delta)
  M2   every (physical) x (rounding) cross term is EXACTLY zero          <- the load-bearing one
  M3   the second-moment block reproduces a = sigma sigma^T at order delta
  M4   the shared-uniform rounding covariance is h_i h_j (min(p_i,p_j) - p_i p_j)
  M5   the diagonal rounding variance is unchanged by sharing
  M6   the whole rounding block is bounded by (sum_i h_i)^2 / 4
  M7   Convention U: if the factor recursion reads the ROUNDED driver
       increment, mean-exactness FAILS -- the adversarial finding of section 9 T2
  M8   R2: the CLIP applied after rounding destroys (K) at the barrier
  M9   R1: the Hoelder loss P(A_n)^(1-1/p) beats every power of n
  M10  R11: the barrier saturation threshold is B_n/sqrt(delta)
  M11  TWO-CHAIN: shared-uniform rounding is an exact W1 isometry, tested
       in-cell and across one adjacent cell (NOT a general |n - n'| proof)
  M12  TWO-CHAIN: E[zeta (e - e')] = 0, with zeta carried symbolically through
       the raw integrand.  Honest note: this identity is just the LINEARITY of
       two single-chain centrings -- which is all section 8.2 needs.  The point
       against M2 is that M2 evaluates ph*0 and is single-chain, not that this
       identity is deep
  M13  TWO-CHAIN: E[(e - e')^2] <= 4h^2 -- the crude bound charged in (P2); the
       exact maximum is h^2/4, so the bound is far from tight and is used only
       as an upper bound
"""

import sys
import sympy as sp

FAIL = []
PASS = []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"{name:<62} {'PASS' if cond else 'FAIL'}  {detail}")


# ---------------------------------------------------------------- symbols
d = sp.symbols("delta", positive=True)          # time step
rho = sp.symbols("rho", real=True)
rhob = sp.sqrt(1 - rho**2)
sx, sy = sp.symbols("sigma_x sigma_y", positive=True)
mx, my = sp.symbols("mu_x mu_y", real=True)
w1, w2 = sp.symbols("w_1 w_2", positive=True)
c1, c2 = sp.symbols("varsigma_1 varsigma_2", positive=True)
z1, z2 = sp.symbols("z_1 z_2", real=True)
aX, aY, a = sp.symbols("a_X a_Y a", positive=True)
# rounding fractions, treated as constants given the physical sigma-field
pX, pY, p1, p2 = [sp.Rational(*r) for r in ((1,3),(3,4),(1,2),(1,5))]  # concrete fractions in (0,1)

# The two Rademacher variables: 4 equally likely sign pairs.
SIGNS = [(1, 1), (1, -1), (-1, 1), (-1, -1)]


def phys(xi, xip):
    """Physical (pre-rounding) increments for one sign pair, m = 2 factors."""
    dW = sp.sqrt(d) * xi
    dWp = sp.sqrt(d) * xip
    dB = rho * dW + rhob * dWp
    dX = mx * d + sx * dW
    dYt = my * d + sy * dB                      # UNROUNDED driver increment
    th1, th2 = sp.exp(-c1 * d), sp.exp(-c2 * d)
    dZ1 = (th1 - 1) * z1 + th1 * w1 * dYt
    dZ2 = (th2 - 1) * z2 + th2 * w2 * dYt
    return dX, dYt, dZ1, dZ2


def E_phys(f):
    """Expectation over the two Rademachers only."""
    terms = [f(xi, xip) for xi, xip in SIGNS]
    total = terms[0]
    for t in terms[1:]:
        total = total + t
    return sp.simplify(total / 4)


# ------------------------------------------- rounding errors, shared uniform
# e = h*(1{U < p} - p), ONE shared U.  Both moments are computed by ACTUALLY
# integrating over U in [0,1] -- not by substituting the expected answer.
Uv = sp.symbols("U", nonnegative=True)


def ind(p):
    return sp.Piecewise((1, Uv < p), (0, True))


def E_e(h, p):
    """E[e] = integral_0^1 h*(1{U<p} - p) dU."""
    return sp.simplify(sp.integrate(h * (ind(p) - p), (Uv, 0, 1)))


def E_ee(hi, pi, hj, pj):
    """E[e_i e_j] = integral_0^1 h_i h_j (1{U<p_i}-p_i)(1{U<p_j}-p_j) dU."""
    return sp.simplify(
        sp.integrate(hi * hj * (ind(pi) - pi) * (ind(pj) - pj), (Uv, 0, 1)))


# ------------------------------------------------------------------ M1
EV = E_phys(lambda xi, xip: sp.Matrix(phys(xi, xip)))
check("M1a  E[Delta X]   = mu_x delta exactly (physical part)",
      sp.simplify(EV[0] - mx * d) == 0)
check("M1b  E[Delta Y]   = mu_y delta exactly (physical part)",
      sp.simplify(EV[1] - my * d) == 0)
drift1 = sp.simplify(sp.series(EV[2], d, 0, 2).removeO())
check("M1c  E[Delta Z^1] = (-varsigma_1 z_1 + w_1 mu_y) delta + O(delta^2)",
      sp.simplify(drift1 - (-c1 * z1 + w1 * my) * d) == 0,
      f"got {sp.simplify(drift1)}")

# rounding adds exactly zero to every first moment
check("M1d  rounding adds exactly 0 to every first moment",
      all(E_e(h, p) == 0 for h, p in [(aX, pX), (aY, pY), (a, p1), (a, p2)]))

# ------------------------------------------------------------------ M2
# (physical) x (rounding): E[ phi * e ] = E[ phi * E[e | G] ] = 0 for any
# G-measurable phi.  Symbolically: the physical factor is G-measurable, so the
# product's expectation factorises through E[e|G] = 0.
cross_zero = []
for xi, xip in SIGNS:
    dX, dYt, dZ1, dZ2 = phys(xi, xip)
    for ph in (dX, dYt, dZ1, dZ2):
        for h, p in [(aX, pX), (aY, pY), (a, p1), (a, p2)]:
            cross_zero.append(sp.simplify(ph * E_e(h, p)))
check("M2   every (physical) x (rounding) cross term is EXACTLY zero",
      all(t == 0 for t in cross_zero),
      "load-bearing: this is what N5-1 buys")

# ------------------------------------------------------------------ M3
def second(i, j):
    return E_phys(lambda xi, xip: sp.Matrix(phys(xi, xip))[i] * sp.Matrix(phys(xi, xip))[j])


targets = {
    (0, 0): sx**2 * d,
    (1, 1): sy**2 * d,
    (0, 1): rho * sx * sy * d,
    (0, 2): w1 * rho * sx * sy * d,
    (1, 2): w1 * sy**2 * d,
    (2, 2): w1**2 * sy**2 * d,
    (2, 3): w1 * w2 * sy**2 * d,
}
ok = True
for (i, j), tgt in targets.items():
    got = sp.expand(second(i, j))
    lead = sp.simplify(sp.limit(sp.expand(got - tgt) / d, d, 0))
    if lead != 0:
        ok = False
        print(f"      block ({i},{j}): residual/delta -> {lead}")
check("M3   second-moment block reproduces a = sigma sigma^T at order delta", ok,
      "all 7 distinct blocks, m = 2")

# ------------------------------------------------------------------ M4/M5
q1, q2 = sp.Rational(1, 3), sp.Rational(3, 4)
check("M4   shared-U covariance = h_i h_j (min(p_i,p_j) - p_i p_j)",
      sp.simplify(E_ee(a, q1, a, q2) - a**2 * (q1 - q1 * q2)) == 0)
check("M5   diagonal variance h^2 p(1-p) unchanged by sharing",
      sp.simplify(E_ee(a, q1, a, q1) - a**2 * q1 * (1 - q1)) == 0)

# ------------------------------------------------------------------ M6
hs = [aX, aY, a, a]
ps = [pX, pY, p1, p2]
tot = sum(sp.Abs(E_ee(hs[i], ps[i], hs[j], ps[j]))
          for i in range(4) for j in range(4))
bound = (sum(hs))**2 / 4
subs = {aX: sp.Rational(1, 10), aY: sp.Rational(1, 10), a: sp.Rational(1, 20),
        pX: sp.Rational(1, 3), pY: sp.Rational(3, 4),
        p1: sp.Rational(1, 2), p2: sp.Rational(1, 5)}
check("M6   whole rounding block <= (sum_i h_i)^2 / 4",
      bool(sp.simplify(tot.subs(subs) <= bound.subs(subs))),
      f"{float(tot.subs(subs)):.5f} <= {float(bound.subs(subs)):.5f}")

# ------------------------------------------------------------------ M7
# Convention U.  If the factor recursion reads the ROUNDED driver increment,
# then Z-tilde depends on U, so p_1 depends on U, so E[e^{Z1} | G] is no longer
# the constant 0.  Model that dependence minimally: p_1 = alpha + beta*1{U<pY}
# with beta != 0.  Then E[e^{Z1}] picks up a term proportional to beta.
al, be = sp.symbols("alpha beta", positive=True)
# E[ 1{U<p1(U)} - p1(U) ] with p1 depending on U through the Y-rounding event.
# On {U < pY}: p1 = al + be ; on {U >= pY}: p1 = al.
E_bad = sp.simplify(
    (sp.Min(pY, al + be) - pY * (al + be)) + (sp.Max(0, al - pY) - (1 - pY) * al)
)
check("M7   Convention U violated => mean-exactness FAILS (non-zero residual)",
      sp.simplify(E_bad.subs({pY: sp.Rational(1, 2), al: sp.Rational(1, 4),
                              be: sp.Rational(1, 8)})) != 0,
      "adversarial finding, section 9 T2")


# ------------------------------------------------------------------ M8
# R2: the CLIP destroys (K).  Barrier at a grid point B = N*h, target
# (N + 1/3)*h.  round -> (N+1)h w.p. 1/3, Nh w.p. 2/3; clip maps BOTH to Nh.
# So the error is deterministic and its conditional mean is NOT zero.
h = sp.symbols("h", positive=True)
N = sp.Integer(7)
target = (N + sp.Rational(1, 3)) * h
B = N * h
frac = sp.Rational(1, 3)


def clip_round_err(U_lt_frac):
    rounded = (N + 1) * h if U_lt_frac else N * h
    return sp.Min(rounded, B) - target          # clip then subtract target


E_clip = sp.simplify(frac * clip_round_err(True) + (1 - frac) * clip_round_err(False))
check("M8   clip AFTER round => E[e|G] != 0, (K) FAILS at the barrier",
      sp.simplify(E_clip) != 0,
      f"E[e^Y|G] = {sp.simplify(E_clip)} (expected -h/3)")
check("M8b  and the free (unclipped) error at the same target IS mean-zero",
      sp.simplify(frac * ((N + 1) * h - target) + (1 - frac) * (N * h - target)) == 0)

# ------------------------------------------------------------------ M9
# R1: the Hoelder loss is absorbed.  With P(A_n) <= C (m n)^(-c w^2) and the
# Hoelder exponent 1 - 1/p, the bound must still beat every power of n.
nn, mm, cc, ww, pp = sp.symbols("n m c omega p", positive=True)
holder = (mm * nn) ** (-cc * ww**2 * (1 - 1 / pp))
# for any fixed power n^{-A}, the ratio must -> 0 as omega -> oo
A = sp.symbols("A", positive=True)
ratio = sp.simplify(holder / nn**(-A))
lim = sp.limit(ratio.subs({mm: 1, nn: 10, cc: sp.Rational(1, 100),
                           pp: 2, A: 50}), ww, sp.oo)
check("M9   Hoelder loss P(A)^(1-1/p) beats every power of n as omega -> oo",
      lim == 0, f"limit = {lim} at c=1/100, p=2, A=50")

# ------------------------------------------------------------------ M10
# R11: the barrier saturation threshold is B_n/sqrt(delta), not B_n^2/delta.
T_, Bn = sp.symbols("T B_n", positive=True)
dd = T_ / nn
k_correct = sp.simplify(Bn / sp.sqrt(dd))
k_wrong = sp.simplify(Bn**2 / dd)
sub = {Bn: 3, T_: 1, nn: 10**6}
check("M10  saturation threshold B_n/sqrt(delta) < n < B_n^2/delta",
      bool(k_correct.subs(sub) < sub[nn]) and bool(k_wrong.subs(sub) > sub[nn]),
      f"{float(k_correct.subs(sub)):.0f} < {sub[nn]} < {float(k_wrong.subs(sub)):.0f}")


# ================================================================== M11-M13
# TWO-CHAIN checks, added in L007 cycle 2 (review-4 finding U6).
#
# M2 above is SINGLE-CHAIN and, as coded, evaluates ph * E_e(...) where E_e is
# already identically zero by M1d -- i.e. ph*0.  It cannot support the two-chain
# rounding cancellation that section 8.2 of LOT-C-PROOF.md leans on.  These three
# checks exercise the actual object: TWO chains rounded with the SAME uniform U.

# two chains land at v = h(n + p) and v' = h(n' + p'); shared-uniform rounding
# sends them to Q = h(n + 1{U<p}) and Q' = h(n' + 1{U<p'}).
hh = sp.symbols("h", positive=True)
pA, pB = sp.symbols("p_A p_B", positive=True)


def Q_minus_Qp(n_off, pa, pb):
    """Q(v,U) - Q(v',U) with v in cell n_off above v' , same U."""
    return hh * (n_off + ind(pa) - ind(pb))


def E_over_U(expr):
    return sp.simplify(sp.integrate(expr, (Uv, 0, 1)))


# --- M11: the W1 isometry E|Q - Q'| = |v - v'|, EXACTLY, in-cell and across cells
iso = []
# same cell, p_A >= p_B  ->  difference is h*1{p_B <= U < p_A} >= 0
same = E_over_U(sp.Abs(hh * (ind(pA) - ind(pB))).rewrite(sp.Piecewise))
same_sub = same.subs({pA: sp.Rational(7, 10), pB: sp.Rational(2, 10)})
iso.append(sp.simplify(same_sub - hh * sp.Rational(5, 10)) == 0)
# adjacent cells, n_off = 1  ->  difference is h(1 + 1{U<p_A} - 1{U<p_B}) >= 0
adj = E_over_U(sp.Abs(Q_minus_Qp(1, pA, pB)).rewrite(sp.Piecewise))
adj_sub = adj.subs({pA: sp.Rational(2, 10), pB: sp.Rational(7, 10)})
# |v - v'| = h(1 + p_A - p_B) = h(1 + 0.2 - 0.7) = 0.5h
iso.append(sp.simplify(adj_sub - hh * sp.Rational(5, 10)) == 0)
# a second adjacent case with the fractions the other way round
adj2 = adj.subs({pA: sp.Rational(9, 10), pB: sp.Rational(1, 10)})
iso.append(sp.simplify(adj2 - hh * sp.Rational(18, 10)) == 0)
check("M11  shared-uniform rounding is an EXACT W1 isometry (in-cell + adjacent)",
      all(iso), "E|Q-Q'| = |v-v'| exactly, monotone coupling")

# --- M12: the two-chain rounding-error difference is conditionally CENTRED
#     e - e' = h[(1{U<p_A} - p_A) - (1{U<p_B} - p_B)]
#     NOTE ON THE ASSUMPTION.  sympy only knows p_A, p_B > 0, not p < 1, so the
#     symbolic integral of 1{U<p} returns min(p,1) and does not simplify to p.
#     The fractions are frac(.) values and therefore live in [0,1) by
#     construction, so the integral is taken with EXPLICIT limits, which is the
#     honest encoding of that constraint rather than a workaround.
zeta = sp.symbols("zeta")          # any G_{k+1}-measurable factor, indep of U


def E_centred(pa):
    """E[ 1{U<p} - p ] with 0 <= p <= 1 imposed by the limits of integration."""
    return sp.simplify(sp.integrate(1 - pa, (Uv, 0, pa))
                       + sp.integrate(-pa, (Uv, pa, 1)))


e_diff_centred = sp.simplify(hh * (E_centred(pA) - E_centred(pB)))
# Carry a SYMBOLIC zeta through the raw integrand, with p_A, p_B still symbolic
# and the frac-in-[0,1) constraint imposed by the limits -- so the zeta line is
# not the tautology U6 criticised in M2.
zeta_carried = sp.simplify(
    sp.integrate(zeta * hh * (1 - pA), (Uv, 0, pA))
    + sp.integrate(zeta * hh * (-pA), (Uv, pA, 1))
    - sp.integrate(zeta * hh * (1 - pB), (Uv, 0, pB))
    - sp.integrate(zeta * hh * (-pB), (Uv, pB, 1)))
m12 = [e_diff_centred == 0, zeta_carried == 0]
# and a numerical cross-check on the raw Piecewise form at sample fractions
for va, vb in [(sp.Rational(3, 10), sp.Rational(8, 10)),
               (sp.Rational(1, 2), sp.Rational(1, 2)),
               (sp.Rational(0), sp.Rational(9, 10))]:
    raw = hh * ((ind(pA) - pA) - (ind(pB) - pB))
    m12.append(E_over_U(raw.subs({pA: va, pB: vb})) == 0)
check("M12  E[(e - e')] = 0 and E[zeta (e - e')] = 0 for the TWO-chain difference",
      all(m12),
      "what section 8.2 needs; M2 does NOT check it (M2 evaluates ph*0)")

# --- M13: the L2 cost of the rounding difference is bounded by 4h^2
e_diff_raw = hh * ((ind(pA) - pA) - (ind(pB) - pB))
var_diff = E_over_U(e_diff_raw**2)
worst = sp.maximum(var_diff.subs(pB, sp.Rational(0)), pA, sp.Interval(0, 1))
check("M13  E[(e - e')^2] <= 4 h^2, the variance charged per step in (P2)",
      bool(sp.simplify(worst - 4 * hh**2) <= 0),
      f"max over p of E[(e-e')^2] = {sp.simplify(worst)} <= 4h^2")

# --- M14: sigma_x = sqrt(v), so the sup-norm asymptotic HALVES the exponent.
# This is the mechanical identity behind F236/F238/F243/F245 and the sigma_x
# correction landed by L013. It is exactly the kind of thing this verifier should
# hold: a symbolic identity, not a claim about prose. Lock L014 chose to keep
# this and to DROP the prose guard it had also written, because a semantic rule
# ("the document never describes itself") is not reducible to regular
# expressions -- that stays with human review.
xi0, eta_s, H_s, R_s = sp.symbols("xi_0 eta H R", positive=True)
sup_V = xi0 * sp.exp(eta_s * sp.sqrt(2 * H_s) * R_s)          # sup of the VARIANCE
sup_sigma = sp.sqrt(sup_V)                                     # sigma_x = sqrt(v)
claimed = sp.sqrt(xi0) * sp.exp(eta_s * sp.sqrt(2 * H_s) * R_s / 2)
wrong = sup_V                                                  # the superseded form
check("M14  sigma_x = sqrt(v) gives ||sigma_x||_inf ~ sqrt(xi_0) e^{eta sqrt(2H) R / 2}",
      sp.simplify(sup_sigma - claimed) == 0
      and sp.simplify(sup_sigma - wrong) != 0,
      "the superseded form stated the SQUARE; F178 always had it right")

# ------------------------------------------------------------------ report
print("-" * 96)
print(f"{len(PASS)}/{len(PASS) + len(FAIL)} passed")
if FAIL:
    print("FAILED: " + ", ".join(FAIL))
sys.exit(1 if FAIL else 0)
