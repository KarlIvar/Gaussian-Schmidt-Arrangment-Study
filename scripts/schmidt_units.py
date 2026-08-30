"""Schmidt units: the unit theorem for the r-twisted Delta-ratios R_f, the
first-power object w_f, the per-class valuation law, and the Robert index.

Companion to schmidt-units.md.  Four phases:

Phase 1 (data).  For odd n, disc D = 1 - n^2, the r-twisted Delta-ratio
  R_f = u_f^6 * b2^4(b2-1728)^3 / (b1^4(b1-1728)^3)   (b1 = j(f), b2 = j(rf))
  satisfies the LATTICE LEMMA
      R_f = r0^6 * Delta(b_1) / Delta(r^{-1} b_1),      r0 = (n-1)/2,
  (b_1 = [1, m_1] the proper O-ideal of the class, r = [r0, sqrt(D)/2] the
  ambiguous twist ideal) -- verified classwise to full precision.  The level
  polynomials prod_f (x - R_f) are certified to be palindromic INTEGER
  polynomials with constant term 1 for all odd n <= 21: the R_f are algebraic
  units (Theorem 1 of schmidt-units.md proves this for every n).  On
  imprimitive strata the same lattice formula gives units as well, with
  R = 1 identically exactly when the induced twist ideal r*O' is principal.

Phase 2 (per-class valuations).  The Euclidean Delta-data
  G_c = n^12 Delta(Lambda_c)/Delta(Z[i]) obey the per-class law: for every
  prime P over p, p^k || n,
      v_P(G_c) = 0 (p split),  12(p^k-1)/((p-1)p^(k-1)(p+1)) (p inert),
                 6(2^k-1)/2^k (p = 2),
  independent of the class c.  Verified here through Newton polygons of the
  certified integer polynomials D_n(x) = prod_c (x - G_c): a single slope of
  the predicted value at every p | n, for every computed level (n <= 27, 49).

Phase 3 (first power).  w_f := u_f * (g2^2 g3)(tau_{rf}) / (g2^2 g3)(tau_f)
  with gamma_2 = E4/eta^8, gamma_3 = E6/eta^12 at the reduced CM points.
  Identities w_f^6 = R_f, w_f w_{rf} = 1, conj(w_f) = w_{f^{-1}} (the last
  away from boundary reduced forms) are verified exactly.  The minimal
  m(n) | 6 with {w_f^m} Galois-stable (certified integer polynomial) is
  tabulated for odd 3 <= n <= 35.

Phase 4 (Robert index).  (a) The chi-eigenprojections of the unit systems
  against fundamental real quadratic units, from the proved KLF genus closed
  forms.  (b) The cubic layer at Euclidean n = 9, 11, 13: the Delta-coset
  units theta_u = theta / M(n)^{1/3} in the real cubic subfield L_3 of the
  ring class field; certified k-th-root descent (rigorous via Friedman's
  regulator bound R > 0.2052) finds the fundamental unit, whence
      R_L, h_L = L'(0,chi_3)/R_L (certified integer), and
      [O_{L_3}^x : <-1, theta_u>] = 8 h_L.

Certification policy (CLAUDE.md guard rails): precision is set AFTER
imports; integers/rationals accepted only with >= max(20, dps/5) spare
digits (absolute error); no PSLQ anywhere in this script.

Usage:
    python3 scripts/schmidt_units.py                 # phases 1-4, defaults
    python3 scripts/schmidt_units.py --selftest      # full re-verification
    python3 scripts/schmidt_units.py phase1 [n ...]  # single phase
    python3 scripts/schmidt_units.py phase3 --nmax 35
Requires mpmath; sympy used only for irreducibility cross-checks.
"""
import sys
import os
import time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction
from math import gcd

from mpmath import (mp, mpf, mpc, fabs, nstr, exp, pi, log, arg,
                    sqrt as msqrt)

import moduli_invariants as MI
import euclidean_moduli_invariants as EU
import phase_klf as PK
from involution_classmap import (classes_of_disc, is_primitive, compose,
                                 reduce_form)
from involution_experiments import inv_sl2
from proof_check import build_P

# precision is set in main()/selftest(), never at import (guard rail 2)

FRIEDMAN = mpf("0.2052")     # unconditional lower bound for regulators
                             # (E. Friedman, Invent. math. 98 (1989))


# ======================================================================
# A.  elementary helpers
# ======================================================================

def factorint(m):
    m = abs(m)
    out = {}
    d = 2
    while d * d <= m:
        while m % d == 0:
            out[d] = out.get(d, 0) + 1
            m //= d
        d += 1 if d == 2 else 2
    if m > 1:
        out[m] = out.get(m, 0) + 1
    return out


def vp(m, p):
    v = 0
    while m and m % p == 0:
        m //= p
        v += 1
    return v


def Ne(n):
    out = n
    for q in factorint(n):
        if q % 2 == 1:
            out = out * (q - (1 if q % 4 == 1 else -1)) // q
    return out


def content(f):
    return gcd(gcd(f[0], abs(f[1])), f[2])


def pell_unit(d):
    u = 1
    while True:
        for pm in (-4, 4):
            v = d * u * u + pm
            if v >= 0:
                t = int(v ** 0.5)
                for tt in (t - 1, t, t + 1, t + 2):
                    if tt >= 0 and tt * tt == v:
                        return tt, u
        u += 1


def newton_polygon(coeffs, p):
    """coeffs highest-degree first, integers; returns [(v_p(root), mult)]"""
    h = len(coeffs) - 1
    pts = [(h - i, Fraction(vp(c, p))) for i, c in enumerate(coeffs) if c]
    pts.sort()
    hull = []
    for pt in pts:
        while len(hull) >= 2:
            (x1, y1), (x2, y2) = hull[-2], hull[-1]
            if (y2 - y1) * (pt[0] - x1) >= (pt[1] - y1) * (x2 - x1):
                hull.pop()
            else:
                break
        hull.append(pt)
    return [(-(hull[k + 1][1] - hull[k][1]) / (hull[k + 1][0] - hull[k][0]),
             hull[k + 1][0] - hull[k][0]) for k in range(len(hull) - 1)]


def cert_int(x, what):
    v, sp = EU.cert_integer(x)
    assert v is not None, f"{what}: not integer-certified ({sp})"
    return v, sp


def poly_int_certify(roots, what):
    co = EU.poly_from_roots(roots)
    out, msp = [], mp.inf
    for c in co:
        v, sp = cert_int(c, what)
        out.append(v)
        msp = min(msp, sp)
    return out, msp


def poly_int_try(roots):
    co = EU.poly_from_roots(roots)
    out = []
    for c in co:
        v, sp = EU.cert_integer(c)
        if v is None:
            return None
        out.append(v)
    return out


def spare(sp):
    return "inf" if sp == mp.inf else f"{float(sp):.0f}"


def irreducible_ZZ(coeffs):
    """irreducibility over Q of an integer polynomial, via sympy."""
    try:
        import sympy
        x = sympy.symbols('x')
        poly = sum(c * x ** (len(coeffs) - 1 - i)
                   for i, c in enumerate(coeffs))
        fl = sympy.factor_list(poly)[1]
        return len(fl) == 1 and fl[0][1] == 1
    except ImportError:
        return None


# ======================================================================
# B.  hyperbolic machinery: u_f, R_f, and the lattice lemma
# ======================================================================

def hyper_data(n):
    """primitive classes with u_f, j-values, twist map."""
    D = 1 - n * n
    prim = [f for f in classes_of_disc(D) if is_primitive(f)]
    rn = reduce_form((n - 1) // 2, 0, (n + 1) // 2)
    eps = n + msqrt(mpf(n * n - 1))
    u = {f: eps * MI.theta_integral(inv_sl2(build_P(n, f)[0]))[0]
         for f in prim}
    jv = {f: MI.J(MI.cm_point(f)) for f in prim}
    rmap = {f: compose(rn, f, D) for f in prim}
    R = {}
    for f in prim:
        b1, b2 = jv[f], jv[rmap[f]]
        R[f] = u[f] ** 6 * (b2 / b1) ** 4 * ((b2 - 1728) / (b1 - 1728)) ** 3
    return prim, u, jv, rmap, R


def hnf2(rows):
    """2-column integer HNF; returns ((g,0), (x,y))."""
    rows = [list(r) for r in rows if tuple(r) != (0, 0)]
    while True:
        nz = sorted([r for r in rows if r[1] != 0], key=lambda r: abs(r[1]))
        if len(nz) <= 1:
            break
        a0 = nz[0]
        for r in rows:
            if r is a0 or r[1] == 0:
                continue
            q = r[1] // a0[1]
            r[0] -= q * a0[0]
            r[1] -= q * a0[1]
        if len([r for r in rows if r[1] != 0]) <= 1:
            break
    v2 = next(r for r in rows if r[1] != 0)
    g = 0
    for r in rows:
        if r[1] == 0:
            g = gcd(g, r[0])
    if g:
        v2[0] %= g
    return (g, 0), tuple(v2)


def twist_lattice(f, n):
    """exact HNF basis (pairs of Fractions: x + y sqrt(D)) of r^{-1}[1, m1],
    m1 = (-b + sqrt D)/(2a); works for imprimitive classes too."""
    D = 1 - n * n
    r0 = (n - 1) // 2
    a, b, c = f
    mx = (Fraction(-b, 2 * a), Fraction(1, 2 * a))

    def mulK(u_, v_):
        return (u_[0] * v_[0] + u_[1] * v_[1] * D,
                u_[0] * v_[1] + u_[1] * v_[0])
    one = (Fraction(1), Fraction(0))
    omF = (Fraction(0), Fraction(1, 2))
    gens = []
    for g1 in ((Fraction(r0), Fraction(0)), omF):
        for g2 in (one, mx):
            p = mulK(g1, g2)
            gens.append((p[0] / r0, p[1] / r0))
    den = 1
    for (x, y) in gens:
        den = den * x.denominator // gcd(den, x.denominator)
        den = den * y.denominator // gcd(den, y.denominator)
    rows = [(int(x * den), int(y * den)) for (x, y) in gens]
    v1, v2 = hnf2(rows)
    return ((Fraction(v1[0], den), Fraction(v1[1], den)),
            (Fraction(v2[0], den), Fraction(v2[1], den)))


def R_lattice(f, n):
    """R via the lattice lemma: r0^6 Delta([1,m1]) / Delta(r^{-1}[1,m1])."""
    r0 = (n - 1) // 2
    a, b, c = f
    sq = mpc(0, msqrt(mpf(n * n - 1)))
    m1 = mpc(mpf(-b) / (2 * a), msqrt(mpf(n * n - 1)) / (2 * a))
    w1, w2 = twist_lattice(f, n)

    def emb(v):
        return v[0] + v[1] * sq

    def Dlat(za, zb):
        tau = za / zb
        if tau.imag < 0:
            za, zb = zb, za
            tau = za / zb
        return EU.Dq_at(tau) / zb ** 12
    return mpf(r0) ** 6 * Dlat(m1, mpc(1)) / Dlat(emb(w1), emb(w2))


# ======================================================================
# C.  Phase 1: the unit polynomials, and the imprimitive strata
# ======================================================================

# per-level working precision for the R-polynomials
R_DPS = {3: 250, 5: 250, 7: 250, 9: 250, 11: 250, 13: 250, 15: 300,
         17: 300, 19: 300, 21: 420}

# certified record (regression targets; phase-kronecker-limit Thm 5 for
# n = 9..15, new at n = 5, 7, 17, 19, 21)
R_POLYS = {
    3: [1, -1],
    5: [1, -34, 1],
    7: [1, -2702, 1],
    9: [1, -339524, -95354, -339524, 1],
    11: [1, -56529284, 1538876166, -56529284, 1],
    13: [1, -11382984004, 885435408006, -11382984004, 1],
    15: [1, -2628641876392, -21595933374628, -1373071731101336,
         9740462908109254, -1373071731101336, -21595933374628,
         -2628641876392, 1],
    17: [1, -673122277718404, 8553847041196806, -673122277718404, 1],
    19: [1, -186863535844922888, 44665915402536486036508,
         131633377547326082495944, 179879784238619113420870,
         131633377547326082495944, 44665915402536486036508,
         -186863535844922888, 1],
    21: [1, -55348592922774901452, 87282798056992201360611266,
         107123155246958016867315700580, 268188996917013884893702692723695,
         -1954092767897474560785802455843992,
         4260659151306307321316379019419804,
         -1954092767897474560785802455843992,
         268188996917013884893702692723695, 107123155246958016867315700580,
         87282798056992201360611266, -55348592922774901452, 1],
}


def phase1_level(n, verbose=True):
    mp.dps = R_DPS.get(n, 420)
    say = print if verbose else (lambda *a, **k: None)
    D = 1 - n * n
    prim, u, jv, rmap, R = hyper_data(n)
    tol = mpf(10) ** (-(mp.dps * 3) // 5)
    # lattice lemma, classwise
    devmax = mpf(0)
    for f in prim:
        dev = fabs(R_lattice(f, n) - R[f]) / max(1, fabs(R[f]))
        devmax = max(devmax, dev)
    assert devmax < tol, (n, nstr(devmax, 3))
    ints, sp = poly_int_certify([R[f] for f in prim], f"R-poly n={n}")
    pal = ints == ints[::-1] or ints == [-c for c in ints[::-1]]
    assert pal and ints[0] == 1 and abs(ints[-1]) == 1, (n, ints)
    if n in R_POLYS:
        assert ints == R_POLYS[n], (n, "R-poly regression mismatch")
    irr = irreducible_ZZ(ints) if n > 3 else True
    say(f"n={n:2d} h={len(prim)}: lattice lemma dev {nstr(devmax, 3)};  "
        f"R-poly integer, palindromic, const {ints[-1]} "
        f"(spare {spare(sp)}); irreducible: {irr}")
    if verbose and len(str(ints)) < 200:
        say(f"      {ints}")
    # imprimitive strata (lattice-lemma definition)
    from collections import defaultdict
    strata = defaultdict(list)
    for f in classes_of_disc(D):
        g = content(f)
        if g > 1:
            strata[g].append(f)
    stratum_polys = {}
    for g, fs in sorted(strata.items()):
        Rs = [R_lattice(f, n) for f in fs]
        co, sps = poly_int_certify(Rs, f"n={n} content-{g} stratum")
        pal = co == co[::-1]
        assert co[0] == 1 and abs(co[-1]) == 1, (n, g, co)
        # principality of the induced twist ideal <=> R = 1 identically
        allone = all(fabs(Rr - 1) < tol for Rr in Rs)
        say(f"      content-{g} stratum ({len(fs)} classes, disc "
            f"{D // (g * g)}): poly {co}  (unit: const {co[-1]}, "
            f"palindromic {pal});  R = 1 identically: {allone}")
        stratum_polys[g] = (co, allone)
    return ints, stratum_polys


def phase1(levels=None, verbose=True):
    print("=" * 78)
    print("PHASE 1: the r-twisted Delta-ratios are units "
          "(lattice lemma + level polynomials)")
    print("=" * 78)
    out = {}
    for n in levels or sorted(R_DPS):
        out[n] = phase1_level(n, verbose)
    return out


# ======================================================================
# D.  Phase 2: Euclidean per-class valuation law (Newton polygons)
# ======================================================================

def predicted_w(p, k):
    """the per-class valuation v_P(G_c) at p^k || n (v(p) = 1)."""
    if p == 2:
        return Fraction(6 * (2 ** k - 1), 2 ** k)
    if p % 4 == 1:
        return Fraction(0)
    return Fraction(12 * (p ** k - 1), (p - 1) * p ** (k - 1) * (p + 1))


def euclid_D_poly(n):
    forms = EU.reduced_forms(-4 * n * n)
    cv = EU.class_values(n)
    Gs = []
    for f in forms:
        (a, b, dd) = cv[f][2]
        (_, _), (c_, dd_) = EU.build_X(a, b, dd)
        cc = EU.cval(c_)
        Gs.append(mpf(n) ** 12 * EU.Dq_at(EU.cval(dd_) / cc)
                  / (cc ** 12 * EU.Dq_at(mpc(0, 1))))
    ints, sp = poly_int_certify(Gs, f"D_{n}")
    return ints, sp, len(forms)


def phase2(levels=None, verbose=True):
    print("=" * 78)
    print("PHASE 2: per-class valuations v_P(G_c) -- Newton polygons of D_n")
    print("=" * 78)
    levels = levels or [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
                        18, 21, 25, 27, 49]
    for n in levels:
        mp.dps = 460 if n >= 40 else 250
        ints, sp, h = euclid_D_poly(n)
        msgs = []
        for p, k in sorted(factorint(n).items()):
            npg = newton_polygon(ints, p)
            w = predicted_w(p, k)
            ok = len(npg) == 1 and npg[0][0] == w and npg[0][1] == h
            assert ok, (n, p, npg, w)
            msgs.append(f"p={p}^{k}: single slope {w} OK")
        if verbose:
            print(f"n={n:2d} h={h:2d} (spare {spare(sp)}): " + ";  ".join(msgs))
    # the split-ladder witness: v_5(H(0)) = v_5(H(1728)) = 0 at n = 15
    mp.dps = 250
    forms = EU.reduced_forms(-900)
    cv = EU.class_values(15)
    H, spH = poly_int_certify([cv[f][1] for f in forms], "H_-900")
    H0 = H[-1]
    H1728 = sum(H[k] * 1728 ** (len(H) - 1 - k) for k in range(len(H)))
    assert vp(H0, 5) == 0 and vp(H1728, 5) == 0
    if verbose:
        print("split-ladder witness at n=15: v_5(H(0)) = v_5(H(1728)) = 0, "
              "so v_5(u_c^2) = -4 for every class (the 5^{4k} ladder)")


# ======================================================================
# E.  Phase 3: the first-power object w_f and its coherence exponent m(n)
# ======================================================================

def eta(tau):
    q = exp(2 * pi * mpc(0, 1) * tau)
    out = exp(pi * mpc(0, 1) * tau / 12)
    qn = q
    while fabs(qn) > mpf(10) ** (-mp.dps - 30):
        out *= (1 - qn)
        qn *= q
    return out


def g23(tau):
    """gamma_2^2 gamma_3 = (E4/eta^8)^2 (E6/eta^12), canonical branch."""
    E4, E6, _ = MI._E4E6D(tau)
    e = eta(tau)
    return (E4 / e ** 8) ** 2 * (E6 / e ** 12)


def w_values(n):
    D = 1 - n * n
    prim, u, jv, rmap, R = hyper_data(n)
    w = {}
    for f in prim:
        t1, t2 = MI.cm_point(f), MI.cm_point(rmap[f])
        w[f] = u[f] * g23(t2) / g23(t1)
    return prim, w, rmap, R


W_DPS = {3: 250, 5: 250, 7: 250, 9: 250, 11: 250, 13: 250, 15: 300, 17: 300,
         19: 300, 21: 420, 23: 350, 25: 350, 27: 400, 29: 420, 31: 420,
         33: 480, 35: 600}

M_TABLE = {3: 1, 5: 3, 7: 1, 9: 3, 11: 6, 13: 1, 15: 1, 17: 3, 19: 2,
           21: 3, 23: 6, 25: 1, 27: 6, 29: 6, 31: 6, 33: 6, 35: 6}

W_POLYS = {3: [1, 1],
           7: [1, 4, 1],
           13: [1, 50, 123, 50, 1],
           15: [1, 112, -630, 1568, -2109, 1568, -630, 112, 1],
           25: [1, 13346, 24067, -4876, 11653, -4876, 24067, 13346, 1]}


def phase3(nmax=35, verbose=True):
    print("=" * 78)
    print("PHASE 3: first-power Schmidt units w_f = u_f * "
          "(g2^2 g3)(tau_rf)/(g2^2 g3)(tau_f)")
    print("=" * 78)
    table = {}
    for n in range(3, nmax + 1, 2):
        mp.dps = W_DPS.get(n, 500)
        tol = mpf(10) ** (-(mp.dps * 3) // 5)
        prim, w, rmap, R = w_values(n)
        # exact first-power laws
        boundary = {f for f in prim if f[1] != 0 and
                    (f[1] == f[0] or f[0] == f[2])}
        for f in prim:
            assert fabs(w[f] ** 6 - R[f]) < tol * (1 + fabs(R[f])), (n, f)
            assert fabs(w[f] * w[rmap[f]] - 1) < tol, (n, f)
            finv = reduce_form(f[0], -f[1], f[2])
            if f not in boundary and finv not in boundary:
                assert fabs(w[finv] - w[f].conjugate()) < tol * \
                    (1 + fabs(w[f])), (n, f, "mirror")
        mfound = None
        for m in (1, 2, 3, 6):
            co = poly_int_try([w[f] ** m for f in prim])
            if co is not None:
                mfound = m
                break
        assert mfound == M_TABLE[n], (n, mfound)
        if mfound == 1 and n in W_POLYS:
            assert co == W_POLYS[n], (n, co)
        irr = irreducible_ZZ(co) if len(co) < 12 else '-'
        if verbose:
            extra = f"  poly {co}" if len(str(co)) < 90 else ""
            print(f"n={n:2d} h={len(prim):2d}: laws OK "
                  f"(boundary classes: {len(boundary)});  minimal m = "
                  f"{mfound}{extra}"
                  + (f"  irreducible: {irr}" if mfound == 1 else ""))
        table[n] = mfound
    if verbose:
        print(f"m(n) table: {table}")
    return table


# ======================================================================
# F.  Phase 4: the Robert index
# ======================================================================

# proved genus closed forms (phase-kronecker-limit.md sections 4 and 6):
# hyperbolic odd real characters: (n, ks) -> (d2, m_chi) with
# L'(0,chi) = m_chi log eps_{d2};  Euclidean: analogous.
HYP_GENUS = {(11, (0, 1)): (8, Fraction(2)),
             (11, (1, 0)): (40, Fraction(2, 3)),
             (13, (0, 1)): (24, Fraction(1)),
             (13, (1, 0)): (21, Fraction(1)),
             (15, (0, 1)): (28, Fraction(1)),
             (15, (2, 1)): (56, Fraction(1, 2))}
EU_GENUS = {(3,): (12, Fraction(1, 3)), (5,): (5, Fraction(2)),
            (7,): (28, Fraction(1)), (9,): (12, Fraction(4, 3)),
            (11,): (44, Fraction(1)), (13,): (13, Fraction(2))}


def phase4a(verbose=True):
    print("=" * 78)
    print("PHASE 4a: quadratic layer -- eigenprojections against log eps_d")
    print("=" * 78)
    mp.dps = 250
    tol = mpf(10) ** (-(mp.dps * 3) // 5)
    # hyperbolic: sum_f chi(f) log|R_f| = -24 m_chi log eps_d
    for n in (11, 13, 15):
        D = 1 - n * n
        prim, u, jv, rmap, R = hyper_data(n)
        coords, orders = PK.group_data(prim, D)
        chis = PK.characters(coords, orders)
        rn = reduce_form((n - 1) // 2, 0, (n + 1) // 2)
        for (ks, cord, chi, isreal) in chis:
            if cord != 2 or (n, ks) not in HYP_GENUS:
                continue
            if fabs(chi[rn] + 1) > 0.1:
                continue          # even characters project to 0
            d2, mchi = HYP_GENUS[(n, ks)]
            SR = sum(chi[f] * log(fabs(R[f])) for f in prim).real
            t, uu = pell_unit(d2)
            le = log((t + uu * msqrt(mpf(d2))) / 2)
            I = SR / le
            Iex = -24 * mchi
            assert fabs(I - mpf(Iex.numerator) / Iex.denominator) < tol * 30, \
                (n, ks, nstr(I, 20))
            if verbose:
                print(f"hyp n={n} chi{ks}: sum chi log|R| = {Iex} "
                      f"log eps_{d2}   (residual "
                      f"{nstr(fabs(I - mpf(Iex.numerator)/Iex.denominator), 3)})")
    # Euclidean: sum_c chi(c) log|G_c| = -12 m_chi log eps_d
    for n in (3, 5, 7, 9, 11, 13):
        Dd = -4 * n * n
        forms = EU.reduced_forms(Dd)
        cv = EU.class_values(n)
        coords, orders = PK.group_data(forms, Dd)
        chis = PK.characters(coords, orders)
        Gs = {}
        for f in forms:
            (a, b, dd) = cv[f][2]
            (_, _), (c_, dd_) = EU.build_X(a, b, dd)
            cc = EU.cval(c_)
            Gs[f] = mpf(n) ** 12 * EU.Dq_at(EU.cval(dd_) / cc) \
                / (cc ** 12 * EU.Dq_at(mpc(0, 1)))
        for (ks, cord, chi, isreal) in chis:
            if cord != 2:
                continue
            d2, mchi = EU_GENUS[(n,)]
            SG = sum(chi[f] * log(fabs(Gs[f])) for f in forms).real
            t, uu = pell_unit(d2)
            le = log((t + uu * msqrt(mpf(d2))) / 2)
            I = SG / le
            Iex = -12 * mchi
            assert fabs(I - mpf(Iex.numerator) / Iex.denominator) < tol * 30, \
                (n, ks, nstr(I, 20))
            if verbose:
                print(f"euc n={n} chi{ks}: sum chi log|G| = {Iex} "
                      f"log eps_{d2}   (residual "
                      f"{nstr(fabs(I - mpf(Iex.numerator)/Iex.denominator), 3)})")
            break


def kth_root_descent(triple, kmax):
    """triple = [real, z, conj z], product +-1; descend through integral
    k-th roots (certified minimal polynomials); returns (triple, index)."""
    cur = triple[:]
    index = 1
    changed = True
    while changed:
        changed = False
        for k in range(kmax, 1, -1):
            r, z = cur[0], cur[1]
            cands = []
            if r > 0:
                cands = [fabs(r) ** (mpf(1) / k)]
                if k % 2 == 0:
                    cands.append(-cands[0])
            elif k % 2 == 1:
                cands = [-fabs(r) ** (mpf(1) / k)]
            for rr in cands:
                for j in range(k):
                    zz = fabs(z) ** (mpf(1) / k) * \
                        exp(mpc(0, 1) * (arg(z) + 2 * pi * j) / k)
                    co = poly_int_try([rr, zz, zz.conjugate()])
                    if co is not None and abs(co[-1]) == 1:
                        cur = [rr, zz, zz.conjugate()]
                        index *= k
                        changed = True
                        break
                if changed:
                    break
            if changed:
                break
    return cur, index


CUBIC_RECORD = {
    9: ([1, -11708931, 115597311109635, -1], [1, 15, 57, -1], 1, 8),
    11: ([1, 2297078781, 2651044211389651971, -1], [1, -25, 201, -1], 1, 8),
    13: ([1, 446643445245, 61048319249786206560771, -1],
         [1, -1, 9, -1], 3, 24),
}


def phase4b(verbose=True):
    print("=" * 78)
    print("PHASE 4b: cubic layer -- fundamental units, h_L, "
          "and the index [O_L^x : <±1, theta_u>] = 8 h_L")
    print("=" * 78)
    mp.dps = 250
    for n in (9, 11, 13):
        t0 = time.time()
        Dd = -4 * n * n
        forms = EU.reduced_forms(Dd)
        cv = EU.class_values(n)
        coords, orders = PK.group_data(forms, Dd)
        chis = PK.characters(coords, orders)
        chi3 = next(c for (ks, o, c, r) in chis if o == 3)
        prods = [None, None, None]
        for f in forms:
            wv = chi3[f]
            i = 0 if fabs(wv - 1) < 0.1 else (1 if wv.imag > 0 else 2)
            (a, b, dd) = cv[f][2]
            (_, _), (c_, dd_) = EU.build_X(a, b, dd)
            cc = EU.cval(c_)
            x = mpf(n) ** 12 * EU.Dq_at(EU.cval(dd_) / cc) \
                / (cc ** 12 * EU.Dq_at(mpc(0, 1)))
            prods[i] = x if prods[i] is None else prods[i] * x
        # unit-normalize by M(n)^{1/3} (the per-class law makes this integral)
        Mn = 1
        for p, k in factorint(n).items():
            if p == 2:
                Mn *= p ** (3 * (2 ** k - 1) * Ne(n // p ** k))
            elif p % 4 == 3:
                Mn *= p ** (6 * (p ** k - 1) // (p - 1) * Ne(n // p ** k))
        M13 = round(abs(Mn) ** (1 / 3))
        assert M13 ** 3 == abs(Mn)
        tru = [mpc(p) / M13 for p in prods]
        co_u, sp_u = poly_int_certify(tru, f"n={n} unit-normalized cubic")
        assert abs(co_u[-1]) == 1, (n, co_u)
        tru.sort(key=lambda z: fabs(z.imag))
        assert fabs(tru[0].imag) < mpf(10) ** (-mp.dps // 3)
        triple = [tru[0].real, tru[1], tru[2]]
        theta_log = fabs(log(fabs(triple[0])))
        fund, idx = kth_root_descent(triple, 40)
        fund_co, sp_f = poly_int_certify(fund, f"n={n} fundamental unit")
        RL = fabs(log(fabs(fund[0])))
        # rigorous fundamentality: any residual index J satisfies
        # R_L = J * R_true >= J * 0.2052 (Friedman), and every k <= 40 was
        # tested without a further integral root
        assert RL / FRIEDMAN < 40, (n, nstr(RL, 10))
        # L'(0, chi3) by the independent Epstein evaluation
        Lp, reps, M = PK.epstein_Lprime0(forms, Dd, chis)
        i3 = next(i for i, (ks, o, c, r) in enumerate(chis) if o == 3)
        Lval = Lp[i3].real
        hL = Lval / RL
        v, sp = EU.cert_integer(hL)
        assert v is not None and v >= 1, (n, nstr(hL, 20))
        index = theta_log / RL
        vi, spi = EU.cert_integer(index)
        assert vi == 8 * v, (n, vi, v)
        # Stark relation log|theta_u| = -8 L'(0,chi3)
        assert fabs(theta_log - 8 * Lval) < mpf(10) ** (-(mp.dps * 3) // 5)
        rec = CUBIC_RECORD[n]
        assert (co_u, fund_co, v, vi) == rec, (n, co_u, fund_co, v, vi)
        if verbose:
            print(f"n={n}: theta_u cubic {co_u}")
            print(f"      fundamental unit {fund_co}  R_L = {nstr(RL, 20)}")
            print(f"      h_L = L'(0,chi3)/R_L = {v} (spare {spare(sp)});  "
                  f"index [O_L^x : <±1, theta_u>] = {vi} = 8 h_L   "
                  f"({time.time() - t0:.0f}s)")


# ======================================================================
# G.  driver
# ======================================================================

def main(argv):
    mode = 'all'
    levels = []
    nmax = 35
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == '--selftest':
            mode = 'all'
        elif a in ('phase1', 'phase2', 'phase3', 'phase4', 'all'):
            mode = a
        elif a == '--nmax':
            i += 1
            nmax = int(argv[i])
        else:
            levels.append(int(a))
        i += 1
    mp.dps = 250
    if mode in ('phase1', 'all'):
        phase1(levels or None)
    if mode in ('phase2', 'all'):
        phase2(levels or None)
    if mode in ('phase3', 'all'):
        phase3(nmax)
    if mode in ('phase4', 'all'):
        phase4a()
        phase4b()
    print("ALL CHECKS PASSED")


if __name__ == '__main__':
    main(sys.argv[1:])
