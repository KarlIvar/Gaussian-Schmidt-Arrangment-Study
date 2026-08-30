"""The phase atlas: Schmidt circles of the ideal triangle colored by the phase.

Companion to phase-atlas.md (and outlook.md 1.6).  For each odd level n the
level-n Schmidt circles inside the ideal triangle T(0, 1, oo) are enumerated
with the alpha_circles machinery, mapped to their form classes of discriminant
D = 1 - n^2, and each primitive class is given its phase unit

    u_f = eps * Theta_f,     eps = n + sqrt(n^2 - 1),

computed by TWO independent routes (any disagreement is a bug):

  route A (canonical matrix; moduli-invariants.md 5.4-5.5): the closed form
      u = -eps mu^-2 h2(b1)/h2(b2), evaluated as eps * j'(m1) X'(m2) /
      conj(j'(conj m2)) on the Lemma-A matrix X_f of class-formula-proof.md
      (this is scripts/moduli_invariants.theta_integral);
  route B (derivative of the modular correspondence; first-power-descent.md):
      u = -r0 h2(b1)/h2(r^-1 b1), r0 = (n-1)/2, with the second lattice
      produced by exact integer HNF arithmetic and both h2-values evaluated
      after exact SL2(Z)-reduction.  By Proposition 1.3 + Theorem 3.3 there,
      this IS the value Phi_y/Phi_x(beta1, beta2) of the modular polynomial
      Phi_{r0} -- evaluated through its uniformization instead of its integer
      coefficients.

j' = -2 pi i E4^2 E6 / Delta is computed via theta constants at full working
precision (never a library default), and mp.dps is set inside functions, after
imports (CLAUDE.md guard rails).  Exactness anchors:

  * every u_f at odd n <= 17 is checked to be a root of the PUBLISHED integer
    level polynomial Q_n of moduli-invariants.md 5.9 (whose entries for
    n <= 13 were derived by the exact rational arithmetic of
    first-power-descent.md Theorem 4.2, and certified at n = 15, 17);
  * at n = 5 the phase is additionally evaluated as Phi_y/Phi_x(beta1, beta2)
    from the CLASSICAL integer modular polynomial Phi_2 (itself re-certified
    here against the functional identity Phi_2(j(2 tau), j(tau)) = 0).

Circles whose form is imprimitive get the same treatment when their core
discriminant is < -4 (Theta is still a well-defined class invariant); the
elliptic-core classes g*(1,1,1) and g*(1,0,1) are excluded (j' vanishes at
rho and i -- the 0/0 already met at alpha = 2) and drawn neutrally.

Usage (from the repo root):
    python3 scripts/make_phase_atlas.py --selftest          # laws + anchors
    python3 scripts/make_phase_atlas.py --figures           # all PNGs
    python3 scripts/make_phase_atlas.py --signs             # sign-law data
    python3 scripts/make_phase_atlas.py --figures --levels 7 15
"""
import sys
import math
import argparse
from fractions import Fraction
from math import gcd

sys.path.insert(0, 'scripts')

from mpmath import (mp, mpf, mpc, fabs, nstr, pi as mppi, nint, atan2, log,
                    sqrt as msqrt, arg as marg)

from involution_classmap import (classes_of_disc, reduce_form, compose,
                                 inverse, is_primitive)
from involution_experiments import inv_sl2
from proof_check import build_P
import moduli_invariants as MI
from alpha_circles import alpha_circles, hurwitz, phi as disk_phi, phi_circle

# --------------------------------------------------------------------------
# published exact level polynomials Q_n (moduli-invariants.md 5.9) --
# the anchors of the cross-check; coefficients highest degree first.
# --------------------------------------------------------------------------

PUBLISHED_Q = {
    3: [1, 1],
    5: [6647, 30594194, 6647],
    7: [11891, 80674200806, 11891],
    9: [10565574794063311, 73919532109765731422845124,
        -118807282021266004510100774, 73919532109765731422845124,
        10565574794063311],
    11: [76575720951, 466015525084217238173676,
         -216521978405797871634733654, 466015525084217238173676,
         76575720951],
    13: [722610532225, 3464286958371072692766958316,
         4603575719671472165025576604518, 3464286958371072692766958316,
         722610532225],
    15: [231902488879724417597324208272447,
         820767611794540060586926641691129674477616253164712,
         -3736910623206271622009575975230452428835195604215370332,
         4648900280945215092152336462310642642094781730724532480920,
         -1793168822031308965307730451798651893626330340481570873697350,
         4648900280945215092152336462310642642094781730724532480920,
         -3736910623206271622009575975230452428835195604215370332,
         820767611794540060586926641691129674477616253164712,
         231902488879724417597324208272447],
    17: [819697933195874886721, 2045693252535068591803195236166141118259260,
         4207718806754010047859719149317312347901942534,
         2045693252535068591803195236166141118259260,
         819697933195874886721],
}

# classical modular polynomial Phi_2 (re-certified in selftest before use)
PHI2 = {(3, 0): 1, (0, 3): 1, (2, 2): -1, (2, 1): 1488, (1, 2): 1488,
        (2, 0): -162000, (0, 2): -162000, (1, 1): 40773375,
        (1, 0): 8748000000, (0, 1): 8748000000, (0, 0): -157464000000000}

LEVELS = list(range(3, 42, 2))          # the atlas levels
SPOT = [101]                            # spot levels (dataset + contact sheet)


# --------------------------------------------------------------------------
# modular forms: F = E4^2 E6 / Delta (j' = -2 pi i F), evaluated after exact
# SL2(Z) reduction (weight-2 cocycle) so theta series always converge fast
# --------------------------------------------------------------------------

def sl2_reduce(tau):
    """(tau_red, (c, d)) with tau_red = gamma tau in the fundamental domain,
    (c, d) the bottom row of gamma."""
    a_, b_, c, d = 1, 0, 0, 1
    t = tau
    for _ in range(100000):
        k = int(nint(t.real))
        if k:
            t -= k
            a_, b_ = a_ - k * c, b_ - k * d
        if abs(t) < 1 - mpf(10) ** (-30):
            t = -1 / t
            a_, b_, c, d = -c, -d, a_, b_
        else:
            break
    return t, (c, d)


def F_at(tau):
    """F(tau) = E4^2 E6 / Delta via theta constants (weight 2, reduced)."""
    tr, (c, d) = sl2_reduce(tau)
    E4, E6, Dl = MI._E4E6D(tr)
    return (E4 * E4 * E6 / Dl) / (c * tau + d) ** 2


def J_at(tau):
    tr, _ = sl2_reduce(tau)
    E4, _, Dl = MI._E4E6D(tr)
    return E4 ** 3 / Dl


# --------------------------------------------------------------------------
# the two phase routes
# --------------------------------------------------------------------------

def u_route_A(n, f):
    """u = eps * Theta(X_f): canonical-matrix route (= -eps mu^-2 h2/h2)."""
    eps = n + msqrt(mpf(n * n - 1))
    th, _br, _m1, _m2 = MI.theta_integral(inv_sl2(build_P(n, f)[0]))
    return eps * th


def _hnf2(rows):
    """HNF basis (A,0),(B,C) of the integer row lattice, A > 0, C > 0."""
    rows = [list(r) for r in rows if any(r)]
    while True:
        nz = [r for r in rows if r[1] != 0]
        if len(nz) <= 1:
            break
        nz.sort(key=lambda r: abs(r[1]))
        p = nz[0]
        for r in nz[1:]:
            q = r[1] // p[1]
            r[0] -= q * p[0]
            r[1] -= q * p[1]
        rows = [r for r in rows if any(r)]
    sec = [r for r in rows if r[1] != 0]
    assert len(sec) == 1
    B, C = sec[0]
    A = 0
    for r in rows:
        if r[1] == 0:
            A = gcd(A, r[0])
    assert A > 0
    if C < 0:
        B, C = -B, -C
    B %= A
    return A, B, C


def routeB_lattice(n, f):
    """Exact HNF data (A, B, C) of (2a) * r * b1 in the basis (1, omega0):
    r^-1 b1 = (A/(2 a r0)) [1, tau2], tau2 = (B + C omega0)/A."""
    a, b, _c = f
    r0 = (n - 1) // 2
    return _hnf2([(2 * a * r0, 0), (-r0 * b, 2 * r0), (0, 2 * a),
                  ((1 - n * n) // 2, -b)])


def u_route_B(n, f):
    """u = -r0 h2(b1)/h2(r^-1 b1): the derivative-of-Phi_{r0} route."""
    a, b, _c = f
    r0 = (n - 1) // 2
    A, B, C = routeB_lattice(n, f)
    w0 = mpc(0, msqrt(mpf(n * n - 1)) / 2)          # omega0 = sqrt(D)/2
    m1 = (mpf(-b) + 2 * w0) / (2 * a)
    tau2 = (B + C * w0) / A
    return -r0 * F_at(m1) * (mpf(A) / (2 * a * r0)) ** 2 / F_at(tau2)


def core_data(f):
    """(content g, core form f/g, core discriminant)."""
    g = gcd(gcd(f[0], f[1]), f[2])
    core = (f[0] // g, f[1] // g, f[2] // g)
    return g, core, core[1] ** 2 - 4 * core[0] * core[2]


def elliptic_core(f):
    """True iff the phase is 0/0 at this class: core disc -3 or -4."""
    return core_data(f)[2] in (-3, -4)


# --------------------------------------------------------------------------
# per-level dataset
# --------------------------------------------------------------------------

class LevelData:
    """All phase data of one level: classes, u-values, laws, circles."""

    def __init__(self, n, dps=100, verify_dps=140, want_circles=True):
        self.n = n
        self.D = 1 - n * n
        self.dps = verify_dps
        self.r = reduce_form((n - 1) // 2, 0, (n + 1) // 2)
        self.classes = classes_of_disc(self.D)
        self.prim = [f for f in self.classes if is_primitive(f)]
        self.h = len(self.prim)
        self.u = {}                 # class -> mpc (at verify_dps)
        self.route_err = {}         # class -> |uA - uB| relative (primitives)
        self.stab_err = {}          # class -> two-precision drift
        self.excluded = set()       # elliptic-core classes (no phase)
        self.imprimitive_u = {}     # imprimitive class -> mpc (route A only)

        for f in self.classes:
            if elliptic_core(f):
                self.excluded.add(f)
                continue
            mp.dps = dps
            uA1 = u_route_A(n, f)
            mp.dps = verify_dps
            uA2 = u_route_A(n, f)
            scale = max(fabs(uA2), mpf(10) ** (-dps // 2))
            self.stab_err[f] = float(fabs(uA1 - uA2) / scale)
            if is_primitive(f):
                uB2 = u_route_B(n, f)
                self.route_err[f] = float(fabs(uA2 - uB2) / scale)
                self.u[f] = uA2
            else:
                self.imprimitive_u[f] = uA2
        mp.dps = verify_dps

        # class-group maps on primitive classes
        self.twin = {f: compose(self.r, f, self.D) for f in self.prim}
        self.inv = {f: reduce_form(*inverse(f)) for f in self.prim}
        self.ambiguous = [f for f in self.prim if self.inv[f] == f]
        self.sign = {f: (1 if self.u[f].real > 0 else -1)
                     for f in self.ambiguous}

        self.circles = alpha_circles(n) if want_circles else []
        for c in self.circles:
            c['class'] = reduce_form(c['q'], -c['x'], c['m'])

    def log_abs_norm(self, f):
        """log|u_f| / log|u_principal| (the principal class is the max)."""
        top = max(float(log(fabs(self.u[g]))) for g in self.prim)
        top = max(top, 1e-9)
        return float(log(fabs(self.u[f]))) / top

    def check_laws(self, tol=None):
        """Pixel-independent re-check of the proved laws on the computed
        data.  Returns a dict of worst-case errors."""
        tol = tol or mpf(10) ** (-(self.dps - 45))
        worst = {'conj': mpf(0), 'twin': mpf(0), 'ambig': mpf(0),
                 'routes': max(self.route_err.values()) if self.route_err else 0.0,
                 'stability': max(self.stab_err.values()) if self.stab_err else 0.0}
        for f in self.prim:
            u = self.u[f]
            # law 1 (conjugation): u_{f^-1} = conj(u_f)
            worst['conj'] = max(worst['conj'],
                                fabs(self.u[self.inv[f]] - u.conjugate())
                                / max(fabs(u), mpf(1)))
            # law 2 (twin): u_{r f} u_f = 1
            worst['twin'] = max(worst['twin'],
                                fabs(self.u[self.twin[f]] * u - 1))
        for f in self.ambiguous:
            worst['ambig'] = max(worst['ambig'],
                                 fabs(self.u[f].imag) / max(fabs(self.u[f]), mpf(1)))
        # principal class is the max of |u|
        principal = reduce_form(1, 0, (self.n * self.n - 1) // 4)
        assert principal in self.prim
        worst['principal_max'] = all(
            fabs(self.u[principal]) >= fabs(self.u[f]) - mpf(10) ** (-20)
            for f in self.prim)
        ok = (worst['conj'] < tol and worst['twin'] < tol
              and worst['ambig'] < tol and worst['routes'] < float(tol)
              and worst['stability'] < float(tol) and worst['principal_max'])
        return ok, worst

    def check_published_Q(self):
        """max_f |Q_n(u_f)| / scale against the published integer polynomial
        (None when the level has no published polynomial)."""
        if self.n not in PUBLISHED_Q:
            return None
        co = PUBLISHED_Q[self.n]
        worst = mpf(0)
        for f in self.prim:
            v, scale, x = mpc(0), mpf(0), self.u[f]
            for cf in co:
                v = v * x + cf
                scale = scale * fabs(x) + fabs(mpf(cf))
            worst = max(worst, fabs(v) / max(scale, mpf(1)))
        return worst


def phi2_check_and_anchor():
    """(certification residual of Phi_2, worst n=5 anchor error).

    Certifies the hard-coded Phi_2 by the functional identity
    Phi_2(j(2 tau), j(tau)) = 0 and by symmetry, then evaluates
    u = Phi_y/Phi_x(beta1, beta2) at both primitive classes of n = 5
    (beta2 = j(r^-1 b1) from the exact route-B lattice) and compares with
    route A.  This is the exact Phi_y/Phi_x route of first-power-descent.md,
    run through the INTEGER coefficients of Phi_2."""
    mp.dps = 140
    assert all(PHI2[(i, j)] == PHI2[(j, i)] for (i, j) in PHI2)
    worst_id = mpf(0)
    for tau in (mpc(0.31, 1.22), mpc(-0.42, 0.91), mpc(0.11, 2.03)):
        x, y = J_at(2 * tau), J_at(tau)
        v = sum(c * x ** i * y ** j for (i, j), c in PHI2.items())
        scale = sum(fabs(mpf(c)) * fabs(x) ** i * fabs(y) ** j
                    for (i, j), c in PHI2.items())
        worst_id = max(worst_id, fabs(v) / scale)

    n = 5
    worst_anchor = mpf(0)
    for f in [g for g in classes_of_disc(1 - n * n) if is_primitive(g)]:
        a, b, _c = f
        A, B, C = routeB_lattice(n, f)
        w0 = mpc(0, msqrt(mpf(n * n - 1)) / 2)
        m1 = (mpf(-b) + 2 * w0) / (2 * a)
        beta1, beta2 = J_at(m1), J_at((B + C * w0) / A)
        phx = sum(c * i * beta1 ** (i - 1) * beta2 ** j
                  for (i, j), c in PHI2.items() if i)
        phy = sum(c * j * beta1 ** i * beta2 ** (j - 1)
                  for (i, j), c in PHI2.items() if j)
        u_phi = phy / phx
        uA = u_route_A(n, f)
        worst_anchor = max(worst_anchor, fabs(u_phi - uA) / fabs(uA))
    return worst_id, worst_anchor


# --------------------------------------------------------------------------
# sign law: ambiguous-class signs against conductor-aware genus characters
# --------------------------------------------------------------------------

def fundamental_discs_dividing(D):
    """All fundamental discriminants d1 (including 1) with d1 | D and
    D/d1 = 0 or 1 mod 4 (the splittings supporting a genus character)."""
    absD = abs(D)
    divs = [d for d in range(1, absD + 1) if absD % d == 0]
    out = []
    for d in divs:
        for d1 in (d, -d):
            if d1 == 1 or _is_fund_disc(d1):
                d2, rem = divmod(D, d1)
                if rem == 0 and d2 % 4 in (0, 1):
                    out.append(d1)
    return sorted(set(out), key=abs)


def _is_fund_disc(d):
    if d == 1:
        return True
    if d % 4 == 1:
        return _squarefree(abs(d))
    if d % 4 == 0:
        m = d // 4
        return m % 4 in (2, 3) and _squarefree(abs(m))
    return False


def _squarefree(m):
    k = 2
    while k * k <= m:
        if m % (k * k) == 0:
            return False
        k += 1
    return True


def kronecker(a, b):
    """Kronecker symbol (a/b)."""
    if b == 0:
        return 1 if abs(a) == 1 else 0
    if a % 2 == 0 and b % 2 == 0:
        return 0
    v = 0
    while b % 2 == 0:
        b //= 2
        v += 1
    k = 1
    if v % 2 == 1 and a % 8 in (3, 5):
        k = -1
    if b < 0:
        b = -b
        if a < 0:
            k = -k
    while a != 0:
        v = 0
        while a % 2 == 0:
            a //= 2
            v += 1
        if v % 2 == 1 and b % 8 in (3, 5):
            k = -k
        if a % 4 == 3 and b % 4 == 3:
            k = -k
        a, b = b % abs(a), abs(a)
    return k if b == 1 else 0


def represented_coprime(f, m, count=6):
    """Values of f coprime to m (small search)."""
    a, b, c = f
    out = []
    for x in range(-12, 13):
        for y in range(-12, 13):
            v = a * x * x + b * x * y + c * y * y
            if v and gcd(v, m) == 1:
                out.append(v)
                if len(out) >= count:
                    return out
    return out


def genus_character(d1, f, D):
    """chi_{d1}([f]) via Kronecker(d1, represented value coprime to D), or
    None if the candidate is NOT well-defined on this class (conductor
    trouble) -- the conductor-aware filter."""
    vals = represented_coprime(f, 2 * abs(D) * abs(d1) if d1 != 1 else 2 * abs(D))
    if not vals:
        return None
    chars = {kronecker(d1, v) for v in vals}
    return chars.pop() if len(chars) == 1 else None


def multiplicativity_test(L):
    """Is psi := -sign a homomorphism Cl[2] -> {+-1}?  Tests every product of
    ambiguous classes with the real Gauss composition.  Returns (bool, list of
    violating triples (f, g, fg))."""
    bad = []
    for i, f in enumerate(L.ambiguous):
        for g in L.ambiguous[i:]:
            fg = compose(f, g, L.D)
            assert fg in L.sign, (L.n, f, g, fg)
            if (-L.sign[f]) * (-L.sign[g]) != -L.sign[fg]:
                bad.append((f, g, fg))
    return not bad, bad


def _real_sign(v, floor=None):
    """sign of a certified-real value; 0 when v is genuinely complex."""
    floor = floor or mpf(10) ** (-mp.dps // 2)
    if fabs(v.imag) > floor * (1 + fabs(v)):
        return 0
    assert fabs(v.real) > floor, v
    return 1 if v.real > 0 else -1


def sign_decomposition(L):
    """Where the sign of u_f on an ambiguous class comes from.  On such a
    class  u = -r0 (A/(2 a r0))^2 F(m1)/F(tau2)  with every factor real, so

        sign(u_f) = -sgn F(m1) * sgn F(tau2),

    and each kernel value splits as F(tau) = F(tau_red) / (c tau + d)^2 into
    an E6-region sign at the reduced point and the ORIENTATION sign of the
    SL2(Z)-reduction cocycle (negative exactly when c tau + d is purely
    imaginary, i.e. when the reduction sends the rational foot Re(tau) to
    the cusp).  Returns a per-class table
        (f, sign u, e1, c1, e2, c2, foot of tau2)
    with e = E6-region sign of the reduced point, c = cocycle sign (0 when
    the value is not individually real: the |tau| = 1 arc classes)."""
    n = L.n
    table = []
    for f in L.ambiguous:
        a, b, _c = f
        A, B, C = routeB_lattice(n, f)
        w0 = mpc(0, msqrt(mpf(n * n - 1)) / 2)
        m1 = (mpf(-b) + 2 * w0) / (2 * a)
        tau2 = (B + C * w0) / A
        row = [f, L.sign[f]]
        prod = 1
        for tau in (m1, tau2):
            tr, (cc, dd) = sl2_reduce(tau)
            E4, E6, Dl = MI._E4E6D(tr)
            e = _real_sign(E4 * E4 * E6 / Dl)
            coc = _real_sign((cc * tau + dd) ** 2)
            row += [e, coc]
            prod *= (e * coc)
        # the product identity is unconditional whenever the factors are real
        if 0 not in row[2:]:
            assert -prod == L.sign[f], (n, f, row)
        table.append(tuple(row) + (Fraction(B, A),))
    return table


def b0_sign_rule(n, f):
    """The closed-form sign law for divisor-type (b = 0) ambiguous classes.

    For f = (a, 0, c) with ac = r0 s0 and gcd(a, c) = 1, split a = a_r a_s,
    c = c_r c_s along r0 = a_r c_r, s0 = a_s c_s (a_r = gcd(a, r0) etc.).
    The exact route-B lattice is tau2 = i y2 with y2 = sqrt(a_r c_s/(a_s c_r))
    and the archimedean mechanism gives

        sign(u_f) = +1  iff  a_r c_s < a_s c_r  (iff y2 < 1),

    i.e. -1 iff the r-twisted divisor a_r c_s exceeds sqrt(N)/2.  Returns the
    predicted sign, or None for the 2-adic (b != 0) ambiguous classes."""
    a, b, c = f
    if b != 0:
        return None
    r0, s0 = (n - 1) // 2, (n + 1) // 2
    a_r, a_s = gcd(a, r0), gcd(a, s0)
    c_r, c_s = gcd(c, r0), gcd(c, s0)
    assert a_r * a_s == a and c_r * c_s == c and a_r * c_r == r0 \
        and a_s * c_s == s0, (n, f)
    # the exact HNF of the twisted lattice in closed form:
    # tau2 = i y2 with y2 = sqrt(a_r c_s / (a_s c_r))
    A, B, C = routeB_lattice(n, f)
    assert (A, B, C) == (2 * r0 * a_s, 0, 2 * a_r), (n, f, (A, B, C))
    return 1 if a_r * c_s < a_s * c_r else -1


def sign_analysis(level_list, out=sys.stdout):
    """The ambiguous-class sign table for all levels, and the character hunt:
    which genus characters chi satisfy sign(u_f) = -chi(f) on every ambiguous
    class of the level."""
    w = out.write
    w("sign law data (outlook.md 1.1): sign(u_f) on ambiguous classes\n")
    w("level | ambiguous classes (reduced form: sign) | matching d1 with "
      "sign = -chi_{d1}\n")
    verdict = {}
    for L in level_list:
        n, D = L.n, L.D
        pattern = "".join('+' if L.sign[f] > 0 else '-' for f in L.ambiguous)
        matches = []
        for d1 in fundamental_discs_dividing(D):
            ok = True
            for f in L.ambiguous:
                ch = genus_character(d1, f, D)
                if ch is None or -ch != L.sign[f]:
                    ok = False
                    break
            if ok:
                matches.append(d1)
        mult_ok, bad = multiplicativity_test(L)
        dec = sign_decomposition(L)
        verdict[n] = (pattern, matches, mult_ok, dec)
        amb = ", ".join(f"{f}:{'+' if L.sign[f] > 0 else '-'}"
                        for f in L.ambiguous)
        w(f"n={n:3d} D={D:6d} h={L.h:2d}  pattern {pattern:10s} [{amb}]\n")
        w(f"      splitting genus characters chi_d1 with sign = -chi: "
          f"{matches if matches else 'NONE'}\n")
        w(f"      -sign multiplicative on Cl[2]: "
          f"{'yes' if mult_ok else 'NO: ' + str(bad[:3])}\n")

        def pm(s):
            return {1: '+', -1: '-', 0: 'o'}[s]
        w("      decomposition -sgnF(m1) sgnF(tau2) [E6-region/cocycle]: "
          + "  ".join(f"{r[0][0]},{r[0][1]}:{pm(r[2])}{pm(r[3])}|"
                      f"{pm(r[4])}{pm(r[5])}@{r[6]}" for r in dec) + "\n")
        # closed-form b = 0 rule
        rule_ok, rule_n = True, 0
        for f in L.ambiguous:
            pred = b0_sign_rule(n, f)
            if pred is not None:
                rule_n += 1
                rule_ok &= (pred == L.sign[f])
        verdict[n] += (rule_ok, rule_n)
        w(f"      closed-form b=0 rule sign = -sgn(a_r c_s - a_s c_r): "
          f"{'MATCHES' if rule_ok else 'FAILS'} ({rule_n} classes)\n")
    chars_fail = [n for n, v in verdict.items() if not v[2]]
    rule_fail = [n for n, v in verdict.items() if not v[4]]
    w("\nVERDICT:\n")
    w(f"  levels where -sign is NOT multiplicative on Cl[2], hence not ANY "
      f"character\n  (genus hypothesis refuted): "
      f"{chars_fail if chars_fail else 'none'}\n")
    w(f"  closed-form b=0 divisor rule fails at: "
      f"{rule_fail if rule_fail else 'no level -- matches every b=0 class'}\n")
    w("  on every b = 0 class both E6-region signs are +, so the sign is "
      "carried by the\n  orientation cocycle of the reduction of the twisted "
      "CM point tau2 (foot printed after @).\n")
    return verdict


# --------------------------------------------------------------------------
# figures
# --------------------------------------------------------------------------

def _hue_color(argu, cmap):
    """color for arg u in (-pi, pi] on the cyclic map: -pi/pi -> red,
    0 -> cyan; twins (u -> 1/u) get mirrored hues, conjugates the same."""
    t = (argu / (2 * math.pi) + 0.5) % 1.0
    return cmap(t)


MIRROR_KW = dict(color='0.35', linewidth=0.6, linestyle=(0, (1.5, 1.5)),
                 zorder=4)


def _draw_mirrors(ax, model, ymax):
    """the mirror lines Re z in {0, 1/2, 1} and |z| = 1, |z-1| = 1
    (2-torsion geometry, circle-composition.md 2).  In the disk model the
    vertical lines need log-spaced heights: linear sampling turns the mapped
    geodesic into isolated dash fragments near the boundary."""
    segs = []
    if model == 'half':
        ts = [ymax * k / 400 for k in range(1, 401)]
    else:
        ts = [math.exp(v / 25) for v in range(-200, 201)]
    for x0 in (0.0, 0.5, 1.0):
        segs.append([complex(x0, t) for t in ts])
    for c0 in (0.0, 1.0):
        segs.append([c0 + math.e ** complex(0, math.pi * t)
                     for t in [k / 600 for k in range(1, 600)]])
    for s in segs:
        if model == 'half':
            pts = [z for z in s if -0.2 <= z.real <= 1.2 and 0 < z.imag <= ymax]
        else:
            pts = [disk_phi(z) for z in s if z.imag > 1e-9]
        ax.plot([z.real for z in pts], [z.imag for z in pts], **MIRROR_KW)


def _draw_triangle(ax, model, ymax):
    tri_kw = dict(color='0.4', linewidth=1.2, zorder=1)
    if model == 'half':
        ax.plot([0, 0], [0, ymax], **tri_kw)
        ax.plot([1, 1], [0, ymax], **tri_kw)
        th = [math.pi * t / 400 for t in range(401)]
        ax.plot([0.5 + 0.5 * math.cos(t) for t in th],
                [0.5 * math.sin(t) for t in th], **tri_kw)
        ax.plot([-0.2, 1.2], [0, 0], color='0.6', linewidth=0.8, zorder=1)
        ax.set_xlim(-0.12, 1.12)
        ax.set_ylim(-0.03, ymax)
    else:
        th = [2 * math.pi * t / 600 for t in range(601)]
        ax.plot([math.cos(t) for t in th], [math.sin(t) for t in th],
                color='0.6', linewidth=0.8, zorder=1)
        ts = [math.exp(v / 25) for v in range(-200, 201)]
        edges = ([complex(0, t) for t in ts], [complex(1, t) for t in ts],
                 [0.5 + 0.5 * math.e ** complex(0, math.pi * k / 600)
                  for k in range(1, 600)])
        for e in edges:
            wpts = [disk_phi(z) for z in e]
            ax.plot([z.real for z in wpts], [z.imag for z in wpts], **tri_kw)
        ax.set_xlim(-1.03, 1.03)
        ax.set_ylim(-1.03, 1.03)
        ax.set_axis_off()
    ax.set_aspect('equal')


def _panel(ax, L, model, mode, cm, twins=True, signs=True, mirrors=True,
           label=True, ycap=None):
    """one atlas panel: mode 'arg' (hue by arg u) or 'logabs'."""
    import matplotlib
    from matplotlib.patches import Circle as MplCircle
    n = L.n
    ymax = n / 2 + 0.6
    if ycap:
        ymax = min(ymax, ycap)
    _draw_triangle(ax, model, ymax)
    if mirrors:
        _draw_mirrors(ax, model, ymax)

    hsv = matplotlib.colormaps['hsv']
    cw = matplotlib.colormaps['coolwarm']

    def face(f):
        if f in L.excluded:
            return (0.75, 0.75, 0.75, 0.9), ':'
        if f in L.imprimitive_u:
            u = L.imprimitive_u[f]
            edge = ':'
        else:
            u = L.u[f]
            edge = '-'
        if mode == 'arg':
            col = _hue_color(float(marg(u)), hsv)
            return (col[0], col[1], col[2], 0.85), edge
        t = float(log(fabs(u)))
        top = max(max(abs(float(log(fabs(L.u[g])))) for g in L.prim), 1e-9)
        col = cw(0.5 + 0.5 * max(-1.0, min(1.0, t / top)))
        return (col[0], col[1], col[2], 0.9), edge

    reps = {}
    for c in L.circles:
        fc, ls = face(c['class'])
        if model == 'half':
            z0, rr = c['center'], c['radius']
        else:
            z0, rr = phi_circle(c['center'], c['radius'])
        ax.add_patch(MplCircle((z0.real, z0.imag), rr, facecolor=fc,
                               edgecolor=(0, 0, 0, 0.5), linewidth=0.4,
                               linestyle=ls, zorder=2))
        f = c['class']
        if f not in reps or c['q'] < reps[f][0]:
            reps[f] = (c['q'], z0, rr)
        if signs and f in L.sign and 300 * rr >= 4.5:
            ax.text(z0.real, z0.imag, '+' if L.sign[f] > 0 else chr(0x2212),
                    ha='center', va='center', zorder=6, clip_on=True,
                    fontsize=min(11.0, 300 * rr), color='black')
    if twins:
        done = set()
        for f in L.u:
            g = L.twin.get(f)
            if g is None or g == f or (g, f) in done or f not in reps \
                    or g not in reps:
                continue
            done.add((f, g))
            z1, z2 = reps[f][1], reps[g][1]
            ax.plot([z1.real, z2.real], [z1.imag, z2.imag], color='0.15',
                    linewidth=0.6, alpha=0.55, zorder=5)
    if label:
        ax.set_title(f"n = {n}:  h = {L.h}, {len(L.circles)} circles, "
                     f"3H = {3 * hurwitz(n * n - 1)}", fontsize=9)


def figure_level(L, save_dir):
    """the per-level figure: two large disk panels (arg u, log|u|) over two
    half-plane strips cropped to the low region (the disk view is faithful
    globally; the strip shows the crowded bottom in the natural chart)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(13.5, 11))
    gs = fig.add_gridspec(2, 2, height_ratios=[2.4, 1.0], hspace=0.16)
    axd1 = fig.add_subplot(gs[0, 0])
    axd2 = fig.add_subplot(gs[0, 1])
    axh1 = fig.add_subplot(gs[1, 0])
    axh2 = fig.add_subplot(gs[1, 1])
    _panel(axd1, L, 'disk', 'arg', None)
    _panel(axd2, L, 'disk', 'logabs', None)
    _panel(axh1, L, 'half', 'arg', None, ycap=2.4)
    _panel(axh2, L, 'half', 'logabs', None, ycap=2.4)
    axd1.set_title(f"disk model, hue = arg u_f   (n = {L.n})", fontsize=10)
    axd2.set_title(f"disk model, shade = log|u_f| / log|u_1|   (n = {L.n})",
                   fontsize=10)
    axh1.set_title("half-plane detail (y < 2.4), hue = arg u_f", fontsize=9)
    axh2.set_title("half-plane detail, shade = log|u_f|/log|u_1|", fontsize=9)
    fig.suptitle(
        f"Phase atlas, level n = {L.n} (D = {L.D}, h = {L.h}): "
        "hue = arg u (red = u < 0, cyan = u > 0); +/- = sign on ambiguous "
        "classes;\nthin segments join r-twin pairs (u -> 1/u: mirrored hue, "
        "inverted shade); dotted curves = mirror geodesics; dotted edge = "
        "imprimitive, gray = elliptic-core", fontsize=10)
    import os
    os.makedirs(save_dir, exist_ok=True)
    path = os.path.join(save_dir, f"phase-atlas-n{L.n}.png")
    fig.savefig(path, dpi=140, bbox_inches='tight')
    plt.close(fig)
    return path


def figure_contact_sheet(levels, save_dir, mode='arg', tag=''):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cols = 5
    rows = (len(levels) + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(3.1 * cols, 3.1 * rows))
    axf = [axes[i // cols][i % cols] if rows > 1 else axes[i % cols]
           for i in range(rows * cols)]
    for ax in axf:
        ax.set_axis_off()
    for ax, L in zip(axf, levels):
        _panel(ax, L, 'disk', mode, None, twins=False, signs=True,
               mirrors=False)
        ax.set_title(f"n = {L.n}  (h = {L.h})", fontsize=8)
    fig.suptitle("The Gaussian Schmidt arrangement, phase portrait: "
                 f"level-n circles of the ideal triangle, hue = arg u_f"
                 if mode == 'arg' else
                 "shade = log|u_f|/log|u_1|", fontsize=12)
    import os
    os.makedirs(save_dir, exist_ok=True)
    path = os.path.join(save_dir, f"phase-atlas-contact{tag}.png")
    fig.savefig(path, dpi=140, bbox_inches='tight')
    plt.close(fig)
    return path


# ------------------------- Euclidean companion (euclidean 6.5) ------------

def euclidean_level(n):
    """[(zeta, u, form, ambiguous, RI)] for the 2 N_e(n) disks of curvature
    2n in the unit square; RI in {'R','I',None} by the PROVED center
    criterion of euclidean-moduli-invariants.md 5.2."""
    import euclidean_moduli_invariants as EU
    Om = EU.Omega_period()
    out = []
    seen = set()
    for Lat in EU.sublattices(n):
        if not EU.is_primitive(*Lat):
            continue
        for Lat2 in (Lat, EU.i_times(*Lat, n)):
            if Lat2 in seen:
                continue
            seen.add(Lat2)
            X = EU.build_X(*Lat2)
            zx, zy = EU.zeta_of(X)
            zx, zy = zx % (2 * n), zy % (2 * n)
            Xc = ((EU.cval(X[0][0]), EU.cval(X[0][1])),
                  (EU.cval(X[1][0]), EU.cval(X[1][1])))
            u = EU.theta_E(Xc) / Om
            form = EU.form_of_lattice(*Lat2)
            fi = EU.reduce_form(form[0], -form[1], form[2])
            ambiguous = (fi == form)
            RI = None
            if ambiguous:
                RI = 'R' if zy % n == 0 else ('I' if zx % n == 0 else None)
            out.append(dict(zeta=(zx, zy), u=u, form=form,
                            ambiguous=ambiguous, RI=RI))
    return out


def euclidean_check(n):
    """assert the R/I criterion against the computed u (pipeline check);
    returns worst deviation from reality/imaginarity on ambiguous disks."""
    worst = mpf(0)
    for d in euclidean_level(n):
        if not d['ambiguous']:
            continue
        u = d['u']
        assert d['RI'] in ('R', 'I'), (n, d)
        dev = fabs(u.imag) if d['RI'] == 'R' else fabs(u.real)
        worst = max(worst, dev / fabs(u))
    return worst


def figure_euclidean(ns, save_dir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle as MplCircle
    hsv = matplotlib.colormaps['hsv']
    fig, axes = plt.subplots(1, len(ns), figsize=(6.4 * len(ns), 6.8))
    if len(ns) == 1:
        axes = [axes]
    for ax, n in zip(axes, ns):
        mp.dps = 60
        data = euclidean_level(n)
        ax.set_aspect('equal')
        for k in range(2):
            ax.plot([0, 1], [k, k], color='0.5', lw=0.8)
            ax.plot([k, k], [0, 1], color='0.5', lw=0.8)
        for x0 in (0.0, 0.5, 1.0):
            ax.plot([x0, x0], [0, 1], **MIRROR_KW)
            ax.plot([0, 1], [x0, x0], **MIRROR_KW)
        r = 1.0 / (2 * n)
        for d in data:
            zx, zy = d['zeta']
            for tx in (0, 2 * n):
                for ty in (0, 2 * n):
                    cx, cy = (zx + tx) / (2 * n), (zy + ty) / (2 * n)
                    if -r < cx < 1 + r and -r < cy < 1 + r:
                        col = _hue_color(float(marg(d['u'])), hsv)
                        ax.add_patch(MplCircle(
                            (cx, cy), r, facecolor=(col[0], col[1], col[2], 0.85),
                            edgecolor=(0, 0, 0, 0.5), linewidth=0.4, zorder=2))
                        if d['RI']:
                            ax.text(cx, cy, d['RI'], ha='center', va='center',
                                    fontsize=min(11.0, 260 * r), zorder=5)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.set_title(f"curvature 2n = {2 * n}:  N_e = {len(data)} disks, "
                     f"hue = arg u; R/I = proved center criterion", fontsize=9)
    fig.suptitle("Euclidean companion: Schmidt disks of the unit square, "
                 "hue = arg u (lemniscatic phase); the R/I overlay is the "
                 "proved criterion of euclidean-moduli-invariants.md 5.2",
                 fontsize=11)
    import os
    os.makedirs(save_dir, exist_ok=True)
    path = os.path.join(save_dir, "phase-atlas-euclidean.png")
    fig.savefig(path, dpi=140, bbox_inches='tight')
    plt.close(fig)
    return path


# --------------------------------------------------------------------------
# selftest and main
# --------------------------------------------------------------------------

def build_levels(levels, verbose=True):
    out = []
    for n in levels:
        L = LevelData(n)
        out.append(L)
        if verbose:
            ok, worst = L.check_laws()
            q = L.check_published_Q()
            qs = f"  Q_n root check {nstr(q, 3)}" if q is not None else ""
            print(f"n={n:3d}  h={L.h:2d}  laws {'OK' if ok else 'FAIL'} "
                  f"(conj {nstr(worst['conj'], 3)}, twin {nstr(worst['twin'], 3)}, "
                  f"routes {worst['routes']:.1e}){qs}")
            assert ok, (n, worst)
            if q is not None:
                assert q < mpf(10) ** (-60), (n, q)
    return out


def selftest():
    print("=" * 78)
    print("phase atlas selftest: laws, route agreement, exact anchors")
    print("=" * 78)
    lv = build_levels(LEVELS + SPOT)
    print("\nPhi_2 anchor (exact modular polynomial route at n = 5):")
    rid, ranc = phi2_check_and_anchor()
    print(f"  Phi_2(j(2t), j(t)) identity residual: {nstr(rid, 3)}")
    print(f"  |Phi_y/Phi_x(beta1, beta2) - u_f| / |u_f| at n=5: {nstr(ranc, 3)}")
    assert rid < mpf(10) ** (-100) and ranc < mpf(10) ** (-100)
    print("\nEuclidean R/I criterion check (proved law vs computed u):")
    mp.dps = 60
    for n in (5, 7, 9, 12, 13):
        dev = euclidean_check(n)
        print(f"  n={n:2d}: worst deviation {nstr(dev, 3)}")
        assert dev < mpf(10) ** (-40)
    print("\nweighted circle count = 3 H(n^2-1) on every atlas level:")
    for L in lv:
        W = sum(c['weight'] for c in L.circles)
        assert W == 3 * hurwitz(L.n ** 2 - 1), L.n
    print("  OK")
    print("\nALL SELFTESTS PASSED")
    return lv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--figures', action='store_true')
    ap.add_argument('--signs', action='store_true')
    ap.add_argument('--levels', type=int, nargs='*', default=None)
    ap.add_argument('--outdir', default='figures')
    args = ap.parse_args()

    if args.selftest:
        selftest()
        return
    levels = args.levels or LEVELS
    if args.signs:
        lv = build_levels(sorted(set(levels)))
        sign_analysis(lv)
        return
    if args.figures:
        lv = build_levels(sorted(set(levels)))
        for L in lv:
            print("wrote", figure_level(L, args.outdir))
        if args.levels is None:
            spot = [LevelData(n) for n in SPOT]
            print("wrote", figure_contact_sheet(lv + spot, args.outdir, 'arg'))
            print("wrote", figure_contact_sheet(lv + spot, args.outdir,
                                                'logabs', tag='-logabs'))
            print("wrote", figure_euclidean([7, 12], args.outdir))
        return
    ap.print_help()


if __name__ == '__main__':
    main()
