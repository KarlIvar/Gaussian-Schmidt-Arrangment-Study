"""The Schmidt line over other imaginary quadratic fields.

Companion to other-fields.md (item 3 of PROGRAM.md).  For K = Q(sqrt d_K),
O_K = Z[w], with d_K in {-8, -3, -7, -11} (class number one; -4 is the
Gaussian case of Papers I-II and is re-run as an anchor) and d_K = -5 (class
number two), the script verifies, in exact O_K-lattice arithmetic (Hermite
normal forms with respect to (1, w)) plus certified Delta-values:

 (A) the classification of S_K = PSL_2(O_K).Rhat by the congruence
       A = sqrt|d| q,  B = i beta,  C = sqrt|d| m,  beta = 1 mod sqrt(d) O_K,
       N(beta) = 1 + |d| q m                        (orbit BFS = congruence set),
     the level set {alpha = -sgn(q) Re beta} and the orientation sign of each
     level, and the descent (an explicit X in SL_2(O_K) for every circle);
 (B) the hyperbolic dictionary: level-alpha circles <-> positive definite forms
     of discriminant D_K(alpha) = 4(alpha^2 - 1)/d_K, the census 3H(|D|), and
     the CLASS FORMULA of the involution sigma(X) = conj(X)^{-1}:
       sigma[f] = [r_alpha] [f]^{-s},  s = orientation sign of the level,
     r_alpha the ambiguous ideal of norm r_0 = (2alpha-2)/gcd(2alpha-2,|d_K|)
     (checked at every class of every listed level, through the descent),
     together with the three lemmas of its proof (A': the explicit unimodular
     P; B': the Gram form g_s(u) = alpha N(u) + s Re(u^2); C': iota(K) = t a_f);
 (C) the unit theorem: R_f = r_0^6 Delta(b_f)/Delta(r^{-1} b_f) has a
     palindromic integer level polynomial with constant term 1 at the first
     levels of each field (imprimitive strata: units too), and the KLF
     sum chi log|R_f| = -24 L'(0, chi) on odd characters, 0 on even ones;
 (D) the Euclidean dictionary N_e(n) = (w_K/2) h(O_n): the Delta-data
     G_c = n^12 Delta(Lambda_c)/Delta(O_K) with D_n(x) = prod (x - G_c) in Z[x],
     the Delta-mass law |M(n)| = prod_{p^k || n, p not split}
     p^{(24/(w_K e_p)) (p^k-1)/(p-1) N_e(n/p^k)} with its sign, and the
     per-class valuations w_p(k) = 12 (p^k - 1)/(e_p (p-1) N_e(p^k)) through
     Newton polygons (a single slope at every p | n);
 (E) the Kronecker limit formula sum chi log|G_c| = -12 L'(0,chi) against the
     independent incomplete-gamma evaluation, and the genus characters:
     L'_prim(0,chi) = (2h(d_1)/w(d_1)) h(d_2) log eps_{d_2} with the real field
     Q(sqrt p) (p = 1 mod 4) or Q(sqrt(p|d_K|)) (p = 3 mod 4) at prime levels;
 (F) the Euler system: the fiber lemma over O_K and the norm relations
       N(G_{c'}) = (-1)^[l=2] G_c^{l+1}/(G_{[l]c} G_{[lbar]c})   (l split, l not | n)
                 = (-1)^[l=2] l^12 G_c^{l+1}                     (l inert)
                 = (-1)^[l=2] pi^12 G_c^{l+1}/G_{[p]c}           (l = p^2 ramified, (pi) = p)
                 = (-1)^[l=2] G_c^{l+1}/G^{(n/l)}_{c^+}          (l | n)
     (the sign is that of A(2) = -2^{-24}; over Q(i) it is visible only at the ramified 2)
     at >= 20 pairs (n, l) per field, and the Hecke recursion;
 (G) the full Robert index with PARI/GP:
       [O_{H_n}^x : mu(H_n) V_n] = 24^(h-1) prod|L'(0,chi)| / R_{H_n}
                                 = (w_K 24^(h-1) / w_{H_n}) h_{H_n} prod C_chi(0),
     mu(H_n) from the conductor lemma, exact exponent matrices (bnfisunit),
     Smith forms, the 24th-root saturation, the cubic layer 8 h_{L_3} C(0) at one
     cubic level per field (bnfcertify at degree <= 16: unconditional);
 (H) class number two, K = Q(sqrt -5), n = 2, 3: the single cusp sees only
     ker(Pic(O_n) -> Pic(O_K)); the two-cusp Delta-data G_c = n^12
     Delta(Lambda)/Delta(O_K Lambda) indexed by all of Pic(O_n); the corrected
     limit formula on the characters pulled back from Pic(O_K); the index
     with the factor w_K/h_K (experimental).

Certification policy (CLAUDE.md guard rails): precision is set in main(),
never at import; every integer/rational read off a real number passes the
absolute-error criterion with >= max(20, dps/5) spare digits; lattices, class
groups, characters, projections and ideal twists are exact; PARI's class
numbers and units are GRH-conditional unless bnfcertify returned 1; no PSLQ.

Usage:
    python3 scripts/other_fields.py --selftest         # everything (~10 min; needs gp)
    python3 scripts/other_fields.py geometry -8        # phases A-B for one field
    python3 scripts/other_fields.py units -3           # phase C
    python3 scripts/other_fields.py euclid -7          # phases D-E
    python3 scripts/other_fields.py euler -11          # phase F
    python3 scripts/other_fields.py index -8 5         # phase G at one level
    python3 scripts/other_fields.py h2                 # phase H
Requires mpmath, sympy (Smith forms) and PARI/GP (`gp`) for phases G-H.
"""
import sys
import os
import re
import time
import shutil
import subprocess
import tempfile
from math import gcd
from fractions import Fraction
from collections import deque
from itertools import product as iproduct

from mpmath import (mp, mpf, mpc, exp, pi, log, fabs, nstr, sqrt, nint, e1,
                    matrix, det)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from involution_classmap import compose, classes_of_disc, is_primitive as form_is_primitive, reduce_form

# precision is set in main(), never here (guard rail 2)

GP = shutil.which("gp")

FIELDS = [-8, -3, -7, -11]          # class number one, transported
ANCHOR = -4                          # the Gaussian case (Papers I-II)


# ==========================================================================
# A.  exact arithmetic in O_K = Z[w] and in O_K-lattices (HNF)
# ==========================================================================

def ext_gcd(a, b):
    if b == 0:
        return (abs(a), (1 if a >= 0 else -1), 0)
    g, s, t = ext_gcd(b, a % b)
    return (g, t, s - (a // b) * t)


def hnf(vecs):
    """HNF (d, b, a) of the Z-span of integer 2-vectors (x, y) <-> x + y w:
    the lattice Z d + Z (b + a w), d > 0, a > 0, 0 <= b < d."""
    vecs = [(x, y) for (x, y) in vecs if (x, y) != (0, 0)]
    xs, g = 0, 0
    for (x, y) in vecs:
        if y == 0:
            continue
        gg, s, t = ext_gcd(g, y)
        xs, g = s * xs + t * x, gg
    assert g > 0, "rank < 2"
    d = 0
    for (x, y) in vecs:
        d = gcd(d, abs(x - (y // g) * xs))
    assert d > 0, "rank < 2"
    return (d, xs % d, g)


def factor(n):
    f = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def vp(m, p):
    k = 0
    while m and m % p == 0:
        m //= p
        k += 1
    return k


class Field:
    """K = Q(sqrt dK), O_K = Z[w], w^2 = t w - n0:  w = sqrt(dK)/2 (dK = 0 mod 4),
    w = (1 + sqrt dK)/2 (dK = 1 mod 4).  Elements are pairs (a, b) = a + b w."""

    def __init__(self, dK):
        assert dK < 0 and dK % 4 in (0, 1)
        self.dK = dK
        if dK % 4 == 0:
            self.t, self.n0 = 0, -dK // 4
        else:
            self.t, self.n0 = 1, (1 - dK) // 4
        self.wK = 4 if dK == -4 else (6 if dK == -3 else 2)
        self.hK = {-3: 1, -4: 1, -7: 1, -8: 1, -11: 1, -19: 1, -20: 2, -24: 2, -15: 2}[dK]
        self.one = (1, 0)
        self.w = (0, 1)
        self.name = {-4: "Q(i)", -8: "Q(sqrt -2)", -20: "Q(sqrt -5)"}.get(dK, f"Q(sqrt {dK})")

    # ---- elements ----
    def mul(self, x, y):
        a, b = x
        c, e = y
        return (a * c - b * e * self.n0, a * e + b * c + b * e * self.t)

    def add(self, x, y):
        return (x[0] + y[0], x[1] + y[1])

    def sub(self, x, y):
        return (x[0] - y[0], x[1] - y[1])

    def neg(self, x):
        return (-x[0], -x[1])

    def smul(self, k, x):
        return (k * x[0], k * x[1])

    def conj(self, x):
        return (x[0] + x[1] * self.t, -x[1])

    def norm(self, x):
        a, b = x
        return a * a + a * b * self.t + b * b * self.n0

    def re(self, x):
        """real part (Fraction)"""
        return Fraction(x[0]) + Fraction(x[1] * self.t, 2)

    def imK(self, x):
        """2 Im(x)/sqrt|dK| = the w-coordinate of x"""
        return x[1]

    def sqrtd(self):
        """sqrt(dK) = 2w - t as an element"""
        return (-self.t, 2)

    def units(self):
        if self.dK == -4:
            return [(1, 0), (-1, 0), (0, 1), (0, -1)]
        if self.dK == -3:
            return [(1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)]
        return [(1, 0), (-1, 0)]

    def chi(self, p):
        """Kronecker symbol (dK/p) for a prime p"""
        d = self.dK
        if p == 2:
            if d % 2 == 0:
                return 0
            return 1 if d % 8 == 1 else -1
        r = d % p
        if r == 0:
            return 0
        return 1 if pow(r, (p - 1) // 2, p) == 1 else -1

    def Ne(self, n):
        """n prod_{p | n} (1 - chi(p)/p) = number of primitive index-n sublattices of O_K"""
        r = n
        for p in factor(n):
            r = r * (p - self.chi(p)) // p
        return r

    def ramified(self):
        return [p for p in factor(-self.dK)]

    def prime_above(self, p):
        """an element pi with N(pi) = p (split or ramified p; h_K = 1), else None"""
        if self.chi(p) == -1:
            return None
        for b in range(0, 2 * p + 2):
            for a in range(-2 * p - 2, 2 * p + 3):
                if self.norm((a, b)) == p:
                    return (a, b)
        raise ValueError((self.dK, p))

    def to_complex(self, x):
        return x[0] + x[1] * self.w_complex()

    def w_complex(self):
        if self.t == 0:
            return mpc(0, sqrt(mpf(self.n0)))
        return mpc(mpf(1) / 2, sqrt(mpf(-self.dK)) / 2)

    # ---- lattices: HNF triples (d, b, a) = Z d + Z (b + a w) ----
    def vecs(self, L):
        d, b, a = L
        return [(d, 0), (b, a)]

    def index(self, L):
        return L[0] * L[2]

    def scale(self, L, lam):
        return hnf([self.mul(v, lam) for v in self.vecs(L)])

    def span_OK(self, L):
        return hnf(self.vecs(L) + [self.mul(v, self.w) for v in self.vecs(L)])

    def is_primitive(self, L):
        return self.index(self.span_OK(L)) == 1

    def lmul(self, L1, L2):
        return hnf([self.mul(u, v) for u in self.vecs(L1) for v in self.vecs(L2)])

    def conj_lat(self, L):
        return hnf([self.conj(v) for v in self.vecs(L)])

    def extend(self, L, n):
        """O_n L, O_n = Z + n O_K"""
        return hnf(self.vecs(L) + [self.mul(v, (0, n)) for v in self.vecs(L)])

    def canon(self, L):
        """canonical representative of the homothety class {u L : u unit}"""
        return min(self.scale(L, u) for u in self.units())

    def intersect_with_ideal(self, L, lam):
        """(L cap lam O_K)/lam as an HNF lattice"""
        d, b, a = L
        N = self.norm(lam)
        lc = self.conj(lam)
        out = []
        for u in range(N):
            for v in range(N):
                x = (u * d + v * b, v * a)
                y = self.mul(x, lc)
                if y[0] % N == 0 and y[1] % N == 0:
                    out.append((y[0] // N, y[1] // N))
        out += [self.mul(v, lc) for v in self.vecs(L)]      # N L / lam = L conj(lam)
        return hnf(out)

    def sublattices(self, L, ell):
        """the ell + 1 sublattices of index ell of the lattice L"""
        d, b, a = L
        v1, v2 = (d, 0), (b, a)
        subs = []
        for j in range(ell):
            subs.append(hnf([self.smul(ell, v1), self.add(v2, self.smul(j, v1))]))
        subs.append(hnf([v1, self.smul(ell, v2)]))
        assert len(set(subs)) == ell + 1
        return subs


class ClassGroup:
    """Pic(O_n) for h_K = 1: canonical primitive index-n sublattices of O_K
    (each class represented by w_K/2 lattices, one canonical)."""

    def __init__(self, K, n):
        self.K, self.n = K, n
        reps = set()
        count = 0
        for a in range(1, n + 1):
            if n % a:
                continue
            d = n // a
            for b in range(d):
                L = (d, b, a)
                if K.is_primitive(L):
                    count += 1
                    reps.add(K.canon(L))
        assert count == K.Ne(n), (K.dK, n, count, K.Ne(n))
        self.reps = sorted(reps)
        self.h = len(self.reps)
        if K.hK == 1:
            assert self.h * K.wK // 2 == K.Ne(n) or n == 1, (K.dK, n, self.h)
        self.idx = {c: k for k, c in enumerate(self.reps)}
        self.one = K.canon((1, 0, n))
        assert self.one in self.idx

    def cls(self, L):
        c = self.K.canon(L)
        assert c in self.idx, (self.K.dK, self.n, L)
        return c

    def prod(self, c1, c2):
        P = self.K.lmul(c1, c2)
        assert self.K.index(P) == self.n and self.K.is_primitive(P)
        return self.cls(P)

    def inv(self, c):
        return self.cls(self.K.conj_lat(c))

    def power(self, c, k):
        r = self.one
        for _ in range(k):
            r = self.prod(r, c)
        return r

    def order(self, c):
        k, r = 1, c
        while r != self.one:
            r = self.prod(r, c)
            k += 1
        return k

    def project(self, c, n2):
        """image in Pic(O_{n2}), n2 | n"""
        L = self.K.extend(c, n2)
        assert self.K.index(L) == n2 and self.K.is_primitive(L), (self.n, n2, c, L)
        return self.K.canon(L)

    def ideal_class(self, lam):
        return self.cls(self.K.intersect_with_ideal(self.one, lam))

    def twist(self, c, lam):
        """[lam] c via Lambda cap lam O_K"""
        L = self.K.intersect_with_ideal(c, lam)
        assert self.K.index(L) == self.n and self.K.is_primitive(L), (self.n, c, lam, L)
        return self.cls(L)


def abelian_structure(cg):
    """(coords, orders): coordinates in a direct-product decomposition"""
    h = cg.h
    ords = {c: cg.order(c) for c in cg.reps}
    known = {cg.one: ()}
    gens, orders = [], []
    while len(known) < h:
        found = False
        for g in sorted(cg.reps, key=lambda c: (-ords[c], c)):
            if g in known:
                continue
            o = ords[g]
            new, ok, p = {}, True, cg.one
            for k in range(o):
                for c, co in known.items():
                    x = cg.prod(c, p) if k else c
                    if x in new:
                        ok = False
                        break
                    new[x] = co + (k,)
                if not ok:
                    break
                p = cg.prod(p, g)
            if ok and len(new) == len(known) * o:
                gens.append(g)
                orders.append(o)
                known = new
                found = True
                break
        assert found, "no direct generator found"
    ng = len(gens)
    coords = {c: tuple(list(co) + [0] * (ng - len(co))) for c, co in known.items()}
    return coords, orders


def all_characters(coords, orders):
    """[(ks, order, phase dict class -> Fraction in [0,1))]"""
    from math import lcm
    out = []
    for ks in iproduct(*[range(o) for o in orders]):
        cord = 1
        for k, o in zip(ks, orders):
            if k:
                cord = lcm(cord, o // gcd(k, o))
        ph = {}
        for c, co in coords.items():
            p = Fraction(0)
            for k, o, e in zip(ks, orders, co):
                p += Fraction(k * e, o)
            ph[c] = p % 1
        out.append((ks, cord, ph))
    return out


def cval(phase):
    if phase == 0:
        return mpc(1)
    if phase == Fraction(1, 2):
        return mpc(-1)
    return exp(2 * pi * mpc(0, 1) * mpf(phase.numerator) / phase.denominator)


# ==========================================================================
# certification helpers
# ==========================================================================

def need_digits():
    return max(20, mp.dps // 5)


def spare_of(err):
    err = fabs(err)
    return mp.dps if err == 0 else int(-log(err, 10))


def cert_int(x, what, need=None):
    need = need_digits() if need is None else need
    if isinstance(x, mpc):
        assert fabs(x.imag) < mpf(10) ** (-need), (what, "imaginary part", nstr(x.imag, 5))
        x = x.real
    k = int(nint(x))
    err = fabs(x - k)
    sp = spare_of(err)
    assert sp >= need, (what, nstr(x, 30), "spare", sp, "need", need)
    return k, sp


def poly_from_roots(roots):
    co = [mpc(1)]
    for r in roots:
        new = [mpc(0)] * (len(co) + 1)
        for i, c in enumerate(co):
            new[i] += c
            new[i + 1] -= c * r
        co = new
    return co


def cert_int_poly(roots, what):
    co = poly_from_roots(roots)
    out, worst = [], mp.dps
    for c in co:
        k, sp = cert_int(c, what)
        out.append(k)
        worst = min(worst, sp)
    return out, worst


# ==========================================================================
# Delta-values
# ==========================================================================

_DELTA_CACHE = {}


def delta_q(tau):
    """Delta_q(tau) = q prod (1 - q^k)^24, q = e^{2 pi i tau}  (without (2pi)^12)"""
    key = (mp.dps, nstr(tau.real, 40), nstr(tau.imag, 40))
    if key in _DELTA_CACHE:
        return _DELTA_CACHE[key]
    q = exp(2 * pi * mpc(0, 1) * tau)
    p = mpf(1)
    aq = fabs(q)
    k = 1
    while True:
        p *= (1 - q ** k) ** 24
        if aq ** k < mpf(10) ** (-mp.dps - 8):
            break
        k += 1
    val = q * p
    _DELTA_CACHE[key] = val
    return val


def delta_lattice(z1, z2):
    """Delta(Z z1 + Z z2)/(2 pi)^12 for complex z1, z2 (any orientation)"""
    tau = z1 / z2
    if tau.imag < 0:
        z1, z2 = z2, z1
        tau = z1 / z2
    return delta_q(tau) / z2 ** 12


def G_value(K, L, n, ref=None):
    """G = n^12 Delta(L)/Delta(ref), ref = O_K by default (L, ref HNF triples)"""
    w = K.w_complex()
    d, b, a = L
    num = delta_lattice(b + a * w, mpf(d))
    if ref is None:
        den = delta_lattice(w, mpf(1))
    else:
        d2, b2, a2 = ref
        den = delta_lattice(b2 + a2 * w, mpf(d2))
    return mpf(n) ** 12 * num / den


# ==========================================================================
# A/B.  geometry: circles (q, beta, m), descent, the level, the class formula
# ==========================================================================

def circle_of(K, X):
    """(q, beta, m) of X(Rhat), X = ((a,b),(c,d)) in SL_2(O_K):
    A = 2 Im(c conj d) = sqrt|d| q, B = i(a conj d - b conj c) = i beta, C = 2 Im(a conj b)"""
    (a, b), (c, d) = X
    q = K.imK(K.mul(c, K.conj(d)))
    beta = K.sub(K.mul(a, K.conj(d)), K.mul(b, K.conj(c)))
    m = K.imK(K.mul(a, K.conj(b)))
    assert K.norm(beta) == 1 + (-K.dK) * q * m, (X, q, beta, m)
    return (q, beta, m)


def mat_mul(K, X, Y):
    (a, b), (c, d) = X
    (e, f), (g, h) = Y
    return ((K.add(K.mul(a, e), K.mul(b, g)), K.add(K.mul(a, f), K.mul(b, h))),
            (K.add(K.mul(c, e), K.mul(d, g)), K.add(K.mul(c, f), K.mul(d, h))))


def mat_inv(K, X):
    (a, b), (c, d) = X
    return ((d, K.neg(b)), (K.neg(c), a))


def mat_conj(K, X):
    (a, b), (c, d) = X
    return ((K.conj(a), K.conj(b)), (K.conj(c), K.conj(d)))


def mat_det(K, X):
    (a, b), (c, d) = X
    return K.sub(K.mul(a, d), K.mul(b, c))


def T_mat(K, lam):
    return (((1, 0), lam), ((0, 0), (1, 0)))


S_MAT = (((0, 0), (-1, 0)), ((1, 0), (0, 0)))
I_MAT = (((1, 0), (0, 0)), ((0, 0), (1, 0)))


def is_congruent_one(K, beta):
    """beta = 1 mod sqrt(dK) O_K"""
    x = K.sub(beta, K.one)
    y = K.mul(x, K.conj(K.sqrtd()))
    return y[0] % (-K.dK) == 0 and y[1] % (-K.dK) == 0


def nearest(K, num, den):
    """a lattice point of O_K nearest to num/den (den > 0)"""
    a, b = num
    v = Fraction(b, den)
    u = Fraction(a, den)
    best, bestn = None, None
    for vv in (int(v) - 1, int(v), int(v) + 1, int(v) + 2):
        for uu in (int(u) - 1, int(u), int(u) + 1, int(u) + 2):
            lam = (uu, vv)
            nn = K.norm(K.sub(num, K.smul(den, lam)))
            if bestn is None or nn < bestn:
                best, bestn = lam, nn
    return best


def translate(K, circ, lam):
    """the circle (q, beta, m) translated by lam: beta' = beta + sqrt(d) q lam,
    m' = m + imK(conj(lam) beta) + q N(lam)"""
    q, beta, m = circ
    beta2 = K.add(beta, K.smul(q, K.mul(K.sqrtd(), lam)))
    m2 = m + K.imK(K.mul(K.conj(lam), beta)) + q * K.norm(lam)
    return (q, beta2, m2)


def invert(K, circ):
    q, beta, m = circ
    return (m, K.conj(beta), q)


def descend(K, circ):
    """an X in SL_2(O_K) with circle_of(X) == circ (Theorem A, sufficiency)"""
    d = -K.dK
    g = I_MAT
    cur = circ
    steps = 0
    while cur[0] != 0:
        q, beta, m = cur
        num = K.mul(beta, K.conj(K.sqrtd()))          # beta conj(sqrt d): lam ~ -num/(|d| q)
        lam = nearest(K, K.neg(num), d * q) if q > 0 else nearest(K, num, d * (-q))
        cur = translate(K, cur, lam)
        g = mat_mul(K, T_mat(K, lam), g)
        cur = invert(K, cur)
        g = mat_mul(K, S_MAT, g)
        steps += 1
        assert steps < 500, "descent does not terminate"
    q0, beta0, m0 = cur
    us = [u for u in K.units() if K.mul(K.conj(u), K.conj(u)) == beta0]
    assert us, ("terminal beta not a unit square", K.dK, cur)
    u = us[0]
    b = K.mul(K.smul(-m0, K.w), K.conj(u))
    Y = ((K.conj(u), b), ((0, 0), u))
    assert circle_of(K, Y) == cur
    X = mat_mul(K, mat_inv(K, g), Y)
    assert mat_det(K, X) == K.one
    assert circle_of(K, X) == circ, (circ, circle_of(K, X))
    return X


def orbit_bfs(K, Q, NB):
    d = -K.dK
    start = (0, (1, 0), 0)
    seen = {start}
    dq = deque([start])
    gens = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    while dq:
        c = dq.popleft()
        nbrs = [translate(K, c, lam) for lam in gens] + [invert(K, c)]
        for nb in nbrs:
            q2, b2, m2 = nb
            if abs(q2) > Q or abs(m2) > Q or K.norm(b2) > NB:
                continue
            if nb not in seen:
                seen.add(nb)
                dq.append(nb)
    return seen


def hnf2_general(v1, v2):
    """HNF of the lattice spanned by two integer 2-vectors: (g0, (p, r)) = Z (g0, 0) + Z (p, r)"""
    rows = [list(v1), list(v2)]
    while rows[0][1] != 0 and rows[1][1] != 0:
        if abs(rows[0][1]) < abs(rows[1][1]):
            rows = [rows[1], rows[0]]
        k = rows[0][1] // rows[1][1]
        rows[0] = [rows[0][0] - k * rows[1][0], rows[0][1] - k * rows[1][1]]
    if rows[1][1] == 0:
        rows = [rows[1], rows[0]]
    g0 = abs(rows[0][0])
    v = rows[1]
    if v[1] < 0:
        v = [-v[0], -v[1]]
    v[0] %= g0
    return (g0, tuple(v))


def residues_mod_qsqrtd(K, q):
    b1 = K.smul(q, K.sqrtd())
    b2 = K.smul(q, K.mul(K.sqrtd(), K.w))
    return hnf2_general(b1, b2)


def reduce_mod(x, H):
    g0, (p, r) = H
    a, b = x
    k = b // r
    a, b = a - k * p, b - k * r
    return (a % g0, b)


def classification_check(K, Qmax=5, say=print):
    """orbit BFS residues (beta mod q sqrt(d) O_K) == congruence residues, |q| <= Qmax;
    returns the number of circles of curvature sqrt|d| q per translation class."""
    d = -K.dK
    orb = orbit_bfs(K, 2 * Qmax + 4, 3000 * d)
    counts = {}
    for q in list(range(1, Qmax + 1)) + list(range(-Qmax, 0)):
        H = residues_mod_qsqrtd(K, abs(q))
        got = set(reduce_mod(beta, H) for (qq, beta, m) in orb if qq == q)
        pred = set()
        g0, (p, r) = H
        for a in range(g0):
            for b in range(r):
                beta = (a, b)
                if is_congruent_one(K, beta) and (K.norm(beta) - 1) % (d * abs(q)) == 0:
                    pred.add(reduce_mod(beta, H))
        assert got == pred, (K.dK, q, len(got), len(pred))
        counts[q] = len(got)
        if K.hK == 1:
            assert len(got) == K.Ne(abs(q)), (K.dK, q, len(got), K.Ne(abs(q)))
    say(f"  (A1) orbit of M_0 (BFS, {len(orb)} circles) = congruence set for every |q| <= {Qmax}: "
        f"circles per translation class N_e(q) = {[counts[q] for q in range(1, Qmax + 1)]}")
    # levels and orientation signs from the orbit
    levels = {}
    for (q, beta, m) in orb:
        if q == 0:
            continue
        alpha = -K.re(beta) * (1 if q > 0 else -1)
        if alpha > 1:
            levels.setdefault(alpha, set()).add(1 if q > 0 else -1)
    return levels


def level_data(K, k):
    """for 2 alpha = k: (D, s) if the level occurs, else None.
    s = orientation sign: the circles of level alpha have sgn(q) = s
    (both signs for K = Q(i))."""
    d = -K.dK
    if (k * k - 4) % d:
        return None
    D = -(k * k - 4) // d
    signs = []
    for s in (1, -1):
        # beta = -s(alpha + b sqrt d/2) must be = 1 mod sqrt(d) for some b of the right parity:
        # test with b = 0 (k even) or b = 1 (k odd)
        b = k % 2
        if K.t == 0:
            if k % 2:
                continue
            beta = (-s * k // 2, -s * b)
        else:
            num = -s * k + s * b       # -s alpha - s b/2 = (-s k + ... )/2 ; beta = (num)/2 + (-s b) w with w = (1+sqrt d)/2
            # -s(alpha + b sqrt d /2) = -s alpha - s b (2w - 1)/2 = (-s k + s b)/2 - s b w
            if num % 2:
                continue
            beta = (num // 2, -s * b)
        if is_congruent_one(K, beta):
            signs.append(s)
    if not signs:
        return None
    return (D, signs)


def circles_of_level(K, k, s):
    """[(form f, circle (q, beta, m))] for all reduced forms of disc D_K(alpha), 2 alpha = k,
    orientation sign s"""
    D, signs = level_data(K, k)
    assert s in signs
    out = []
    for f in classes_of_disc(D):
        a, b, c = f
        x = -s * b
        if K.t == 0:
            beta = (-s * k // 2, x)
        else:
            beta = ((-s * k - x) // 2, x)
        assert K.norm(beta) == 1 + (-K.dK) * s * a * s * c, (f, s, beta)
        assert is_congruent_one(K, beta)
        out.append((f, (s * a, beta, s * c)))
    return out


def form_of_circle(K, circ):
    q, beta, m = circ
    s = 1 if q > 0 else -1
    x = K.imK(beta)
    alpha = -s * K.re(beta)
    return (abs(q), -s * x, abs(m)), alpha


def sigma_image(K, X):
    """the class read from sigma(X)(Rhat), reflected into H when it lies below"""
    Xs = mat_inv(K, mat_conj(K, X))
    circ = circle_of(K, Xs)
    q, beta, m = circ
    f1, alpha1 = form_of_circle(K, circ)
    reflected = False
    if alpha1 < 0:
        circ = (q, K.neg(K.conj(beta)), m)
        f1, alpha1 = form_of_circle(K, circ)
        reflected = True
    return reduce_form(*f1), alpha1, (1 if q > 0 else -1), reflected


def twist_ideal_data(K, k, s):
    """(r0, s0, raw form (r0, r0 or 0, .) of the ideal r_alpha, reduced class of r_alpha, reduced class of s_alpha):
    t0 = 2(alpha + s)/|d| is the norm of the ideal carrying the Gram form g_s,
    r0 = (k - 2)/gcd(k - 2, |d|), s0 = (k + 2)/gcd(k + 2, |d|)   [dK odd]
    r0 = (k - 2)/gcd(k - 2, 2|d|) ..., i.e. r0 = (alpha - 1)/gcd(alpha - 1, |d|/2)-type [dK even]"""
    d = -K.dK
    D = -(k * k - 4) // d
    if K.t == 1:
        r0 = (k - 2) // gcd(k - 2, d)
        s0 = (k + 2) // gcd(k + 2, d)
        assert r0 * s0 == -D
        fr = (r0, r0, (r0 + s0) // 4)
        fs = (s0, s0, (r0 + s0) // 4)
        assert (r0 + s0) % 4 == 0 and fr[1] ** 2 - 4 * fr[0] * fr[2] == D
    else:
        alpha = k // 2
        r0 = (alpha - 1) // gcd(alpha - 1, d // 2)
        s0 = (alpha + 1) // gcd(alpha + 1, d // 2)
        assert 4 * r0 * s0 == -D, (K.dK, k, r0, s0, D)
        fr = (r0, 0, s0)
        fs = (s0, 0, r0)
    return r0, s0, fr, reduce_form(*fr), reduce_form(*fs)


def gram_gs(K, alpha, s, u):
    """g_s(u) = alpha N(u) + s Re(u^2)  (Fraction)"""
    u2 = K.mul(u, u)
    return alpha * K.norm(u) + s * K.re(u2)


def gram_gs_pol(K, alpha, s, u, v):
    """polarization: g_s(u, v) = (g_s(u+v) - g_s(u) - g_s(v))/2"""
    return (gram_gs(K, alpha, s, K.add(u, v)) - gram_gs(K, alpha, s, u) - gram_gs(K, alpha, s, v)) / 2


def lattice_K(K, f, s, k):
    """the lattice K_f = {u in O_K : conj(u) - beta u in a sqrt(d) O_K} as an HNF triple,
    with beta = -s(alpha + b sqrt d / 2)"""
    a, b, c = f
    d = -K.dK
    x = -s * b
    beta = (-s * k // 2, x) if K.t == 0 else ((-s * k - x) // 2, x)
    # sublattice of O_K = Z^2: brute force over a fundamental domain of a*sqrt(d)*O_K
    mod = a * d
    gens = []
    for p_ in range(mod):
        for q_ in range(mod):
            u = (p_, q_)
            y = K.sub(K.conj(u), K.mul(beta, u))
            z = K.mul(y, K.conj(K.sqrtd()))      # y/(a sqrt d) integral iff y conj(sqrt d)/(a |d|) integral
            if z[0] % (a * d) == 0 and z[1] % (a * d) == 0:
                gens.append(u)
    gens += [(mod, 0), (0, mod)]
    return hnf(gens), beta


def sqrtD_ideal_basis(D, a, b):
    """Z-basis of the lattice [a, (-b + sqrt D)/2] in K' = Q(sqrt D): pairs (x, y) = x + y sqrt D (Fractions).
    With b -> -b this is Paper I's ideal a_f = [a, (b + sqrt D)/2] of the form (a, b, c)."""
    return [(Fraction(a), Fraction(0)), (Fraction(-b, 2), Fraction(1, 2))]


def Kp_mul(D, u, v):
    return (u[0] * v[0] + u[1] * v[1] * D, u[0] * v[1] + u[1] * v[0])


def Kp_hnf(D, gens):
    """HNF Z-basis of the Z-span of elements x + y sqrt D (Fractions): returns (den, (g0, 0), (p, r))
    meaning the lattice (1/den)(Z (g0, 0) + Z (p, r))"""
    den = 1
    for (x, y) in gens:
        den = den * x.denominator // gcd(den, x.denominator)
        den = den * y.denominator // gcd(den, y.denominator)
    rows = [(int(x * den), int(y * den)) for (x, y) in gens]
    # HNF via repeated hnf2_general
    cur = None
    for r in rows:
        if r == (0, 0):
            continue
        if cur is None:
            cur = [r]
            continue
        cur.append(r)
    # reduce a list of vectors to a basis
    basis = [list(v) for v in cur]
    # column 2 gcd
    while True:
        nz = [v for v in basis if v[1] != 0]
        if len(nz) <= 1:
            break
        nz.sort(key=lambda v: abs(v[1]))
        v0 = nz[0]
        for v in basis:
            if v is v0 or v[1] == 0:
                continue
            kq = v[1] // v0[1]
            v[0] -= kq * v0[0]
            v[1] -= kq * v0[1]
    v2 = next(v for v in basis if v[1] != 0)
    g0 = 0
    for v in basis:
        if v[1] == 0:
            g0 = gcd(g0, abs(v[0]))
    assert g0 > 0
    if v2[1] < 0:
        v2 = [-v2[0], -v2[1]]
    v2[0] %= g0
    return (den, (g0, 0), (v2[0], v2[1]))


def Kp_lattice_eq(D, gens1, gens2):
    return Kp_hnf(D, gens1) == Kp_hnf(D, gens2)


def class_formula_level(K, k, s, say=print, check_lemmas=True):
    """verifies sigma[f] = [r_alpha][f]^{-s} at every primitive class of the level 2 alpha = k,
    and the lemmas A', B', C' of the proof.  Returns (h, number of classes checked)."""
    d = -K.dK
    D, signs = level_data(K, k)
    alpha = Fraction(k, 2)
    r0, s0, fr_raw, fr, fs = twist_ideal_data(K, k, s)
    prim = [f for f in classes_of_disc(D) if form_is_primitive(f)]
    assert form_is_primitive(fr_raw), (K.dK, k, fr_raw)
    if fr not in prim:
        fr = (fr[0], -fr[1], fr[2])
    assert fr in prim, (K.dK, k, fr)
    unit_form = next(g for g in prim if g[0] == 1)
    assert compose(fr, fr, D) == unit_form, "r_alpha not 2-torsion"
    # r and s represent the same class
    fs2 = fs if fs in prim else (fs[0], -fs[1], fs[2])
    assert fs2 == fr, (K.dK, k, fr, fs2)
    checked = 0
    for (f, circ) in circles_of_level(K, k, s):
        X = descend(K, circ)
        f2, alpha2, s2, refl = sigma_image(K, X)
        assert alpha2 == alpha, (K.dK, k, f, alpha2)
        assert s2 == -1, "sigma-image always negatively oriented"
        assert refl == (s == 1), "reflection needed iff the level is positively oriented"
        if not form_is_primitive(f):
            continue
        finv = reduce_form(f[0], -f[1], f[2])
        pred = compose(fr, finv if s == 1 else f, D)
        assert f2 == pred, (K.dK, k, s, f, f2, pred)
        checked += 1
        if check_lemmas:
            # Lemma A': K_f has index a, P = [[u1, v1],[u2, v2]] with v = s(conj u - beta u)/(a sqrt d)
            a, b, c = f
            Kf, beta = lattice_K(K, f, s, k)
            assert K.index(Kf) == a, (K.dK, k, f, Kf)
            u1, u2 = (Kf[0], 0), (Kf[1], Kf[2])
            if K.imK(K.mul(u1, K.conj(u2))) != s * a:
                u1, u2 = u2, u1
            assert K.imK(K.mul(u1, K.conj(u2))) == s * a
            vs = []
            for u in (u1, u2):
                y = K.smul(s, K.sub(K.conj(u), K.mul(beta, u)))
                z = K.mul(y, K.conj(K.sqrtd()))
                assert z[0] % (a * d) == 0 and z[1] % (a * d) == 0
                vs.append((z[0] // (a * d), z[1] // (a * d)))
            P = ((u1, vs[0]), (u2, vs[1]))
            assert mat_det(K, P) == K.one, (K.dK, k, f, P)
            Xp = mat_inv(K, P)
            assert circle_of(K, Xp) == circ, (circle_of(K, Xp), circ)
            # Lemma B': the sigma-image is the form (2/(a|d|)) g_s in the basis (u2, -u1)
            g11 = gram_gs(K, alpha, s, u2)
            g22 = gram_gs(K, alpha, s, u1)
            g12 = gram_gs_pol(K, alpha, s, u2, K.neg(u1))
            fB = (Fraction(2, a * d) * g11, Fraction(2, a * d) * 2 * g12, Fraction(2, a * d) * g22)
            assert all(v.denominator == 1 for v in fB), (K.dK, k, f, fB)
            fB = tuple(int(v) for v in fB)
            f2p, alpha2p, s2p, reflp = sigma_image(K, Xp)
            assert reduce_form(*fB) == f2p, (K.dK, k, f, fB, f2p)
            # Lemma C': iota(K_f) = t * a_f in K' = Q(sqrt D), iota(x + y sqrt d/2) = t0 x + (y/2) sqrt D,
            # t0 = 2(alpha + s)/|d| (the norm of the ideal carrying g_s); a_f = [a, (b + sqrt D)/2] is the
            # ideal of f in Paper I's convention (positively oriented basis, norm form = f)
            t0 = Fraction(2 * (alpha + s), d)
            assert t0.denominator == 1
            t0 = int(t0)

            def iota(u):
                x = K.re(u)
                y = Fraction(K.imK(u))
                return (t0 * x, y / 2)
            iotaK = [iota(u1), iota(u2)]
            # the ideal t = iota(O_K): t0 x + y/2 sqrt D over O_K
            t_ideal = [iota(K.one), iota(K.w)]
            # a_f = [a, (b + sqrt D)/2]; product t * a_f
            af = sqrtD_ideal_basis(D, a, -b)
            prod_gens = [Kp_mul(D, u, v) for u in t_ideal for v in af]
            assert Kp_lattice_eq(D, iotaK, prod_gens), (K.dK, k, f, "iota(K) != t a_f")
            # norm check: N(iota(u)) = t0 * (2/|d|) g_s(u)
            for u in (u1, u2, K.add(u1, u2)):
                z = iota(u)
                Nz = z[0] * z[0] - D * z[1] * z[1]
                assert Nz == t0 * Fraction(2, d) * gram_gs(K, alpha, s, u)
    say(f"  (B) 2alpha = {k:2d}, s = {s:+d}, D = {D:6d}, h = {len(prim):2d} (all forms {len(classes_of_disc(D))}): "
        f"sigma[f] = [r][f]^{'-1' if s == 1 else '+1'} at every primitive class, r = {fr} "
        f"(norm {r0}; s-ideal norm {s0}), lemmas A'-C' verified" + (" (checked)" if check_lemmas else ""))
    return len(prim), checked


def hurwitz_H(N):
    """Hurwitz class number H(N) (forms of disc -N, weights 1/2 at (a,0,a), 1/3 at (a,a,a))"""
    tot = Fraction(0)
    for (a, b, c) in classes_of_disc(-N):
        if a == b == c:
            tot += Fraction(1, 3)
        elif b == 0 and a == c:
            tot += Fraction(1, 2)
        else:
            tot += 1
    return tot


def geometry_phase(dK, kmax=None, say=print):
    K = Field(dK)
    t0 = time.time()
    say("=" * 78)
    say(f"PHASE A/B  K = {K.name}, d_K = {dK}, O_K = Z[w], w^2 = {K.t} w - {K.n0}, w_K = {K.wK}")
    say("=" * 78)
    levels = classification_check(K, say=say)
    # level set: 2 alpha = k with k = -2 s mod |d| (odd d), etc.
    d = -dK
    got = sorted((int(2 * a), tuple(sorted(sg))) for a, sg in levels.items())
    pred = []
    for k in range(3, 400):
        ld = level_data(K, k)
        if ld is not None:
            pred.append((k, tuple(sorted(ld[1]))))
        if len(pred) >= 8:
            break
    assert got[:8] == pred[:8], (dK, got[:8], pred[:8])
    assert all(level_data(K, k) is not None and tuple(sorted(level_data(K, k)[1])) == sg for (k, sg) in got), (dK, got)
    say(f"  (A2) level set (2 alpha, orientation signs): {got[:10]} ...   = predicted from beta = 1 mod sqrt(d)")
    # the hyperbolic dictionary and the class formula at the first levels
    ks = [k for k in range(3, 200) if level_data(K, k) is not None]
    if kmax is not None:
        ks = [k for k in ks if k <= kmax]
    else:
        ks = ks[:10]
    tot = 0
    for k in ks:
        D, signs = level_data(K, k)
        # census: number of level-alpha circles in the ideal triangle = 3 H(|D|): the forms of disc D
        for s in signs:
            circ = circles_of_level(K, k, s)
            assert len(circ) == len(classes_of_disc(D))
            h, ch = class_formula_level(K, k, s, say=say, check_lemmas=True)
            tot += ch
    say(f"  phase A/B done: class formula checked at {tot} primitive classes over {len(ks)} levels "
        f"({time.time() - t0:.1f} s)")
    return tot


# ==========================================================================
# C.  the hyperbolic units R_f = r0^6 Delta(b)/Delta(r^{-1} b)
# ==========================================================================

def R_unit(K, k, s, f, D, r0, fr):
    """R_f by the lattice lemma (Paper II Lemma 4.1 with the twist ideal r_alpha):
    b = [a, (-b + sqrt D)/2], r = ideal of the form fr = (r0, r0 or 0, .), r^{-1} b = (1/r0) r b."""
    a, b, c = f
    bf = sqrtD_ideal_basis(D, a, b)
    r_id = sqrtD_ideal_basis(D, fr[0], fr[1])       # [r0, (-b_r + sqrt D)/2]
    prod_gens = [Kp_mul(D, u, v) for u in r_id for v in bf]
    den, g1, g2 = Kp_hnf(D, prod_gens)
    sD = mpc(0, sqrt(mpf(-D)))

    def emb(x, y):
        return x + y * sD
    w1 = emb(Fraction(g1[0], den), Fraction(g1[1], den)) / r0
    w2 = emb(Fraction(g2[0], den), Fraction(g2[1], den)) / r0
    z1 = emb(Fraction(-b, 2), Fraction(1, 2))
    z2 = mpc(a, 0)
    return mpf(r0) ** 6 * delta_lattice(z1, z2) / delta_lattice(w1, w2)


def rep_numbers_form(Q, M):
    A, B, C = Q
    D = B * B - 4 * A * C
    r = [0] * (M + 1)
    ymax = int((4 * A * M / (-D)) ** 0.5) + 2
    for y in range(-ymax, ymax + 1):
        disc = B * B * y * y - 4 * A * (C * y * y - M)
        if disc < 0:
            continue
        sq = int(disc ** 0.5) + 2
        for x in range((-B * y - sq) // (2 * A) - 2, (-B * y + sq) // (2 * A) + 3):
            if x == 0 and y == 0:
                continue
            m = A * x * x + B * x * y + C * y * y
            if 1 <= m <= M:
                r[m] += 1
    return r


def epstein_Lprime0_forms(forms, D, chis, w=2):
    """[L'(0,chi)] for characters given as dicts form -> complex value, by the independent
    incomplete-gamma evaluation (Paper II Lemma 2.3)"""
    alpha = 2 * pi / sqrt(mpf(-D))
    M = int((mp.dps + 15) * log(mpf(10)) / alpha) + 10
    reps = {f: rep_numbers_form(f, M) for f in forms}
    Iv = [None] * (M + 1)
    for m in range(1, M + 1):
        am = alpha * m
        Iv[m] = exp(-am) / am + e1(am)
    out = []
    for chi in chis:
        tot = mpc(0)
        for f in forms:
            s_ = mpf(0)
            rf = reps[f]
            for m in range(1, M + 1):
                if rf[m]:
                    s_ += rf[m] * Iv[m]
            tot += chi[f] * s_
        out.append(tot / w)
    return out, M


def form_group_data(forms, D):
    """abelian coordinates of the primitive form classes under composition"""
    unit = next(f for f in forms if f[0] == 1)
    h = len(forms)

    def elem_order(f):
        g, k = f, 1
        while g != unit:
            g = compose(g, f, D)
            k += 1
        return k
    ords = {f: elem_order(f) for f in forms}
    known = {unit: ()}
    gens, orders = [], []
    while len(known) < h:
        for g in sorted(forms, key=lambda f: (-ords[f], f)):
            if g in known:
                continue
            o = ords[g]
            new, p, ok = {}, unit, True
            for k in range(o):
                for f, co in known.items():
                    x = compose(f, p, D) if k else f
                    if x in new:
                        ok = False
                        break
                    new[x] = co + (k,)
                if not ok:
                    break
                p = compose(p, g, D)
            if ok and len(new) == len(known) * o:
                gens.append(g)
                orders.append(o)
                known = new
                break
        else:
            raise RuntimeError("no direct generator")
    ng = len(gens)
    coords = {f: tuple(list(c) + [0] * (ng - len(c))) for f, c in known.items()}
    return coords, orders


# per-field working precision for the unit polynomials
UNIT_DPS = 300


def units_phase(dK, kmax=None, say=print, with_klf=True):
    K = Field(dK)
    t0 = time.time()
    say("=" * 78)
    say(f"PHASE C  hyperbolic units R_f over K = {K.name}   (dps = {mp.dps})")
    say("=" * 78)
    need = need_digits()
    ks = [k for k in range(3, 200) if level_data(K, k) is not None]
    ks = [k for k in ks if k <= (kmax if kmax else 30)]
    results = {}
    for k in ks:
        D, signs = level_data(K, k)
        s = signs[0]
        r0, s0, fr_raw, fr, fs = twist_ideal_data(K, k, s)
        forms = classes_of_disc(D)
        prim = [f for f in forms if form_is_primitive(f)]
        if fr not in prim:
            fr = (fr[0], -fr[1], fr[2])
        R = {f: R_unit(K, k, s, f, D, r0, fr_raw) for f in forms}
        co, sp = cert_int_poly([R[f] for f in prim], f"R-polynomial 2alpha={k}")
        assert abs(co[-1]) == 1, (dK, k, co)
        assert co == co[::-1] or [(-1) ** i * c for i, c in enumerate(co)] == co[::-1] or len(prim) == 1, (dK, k, co)
        # twist law R_{r f} = 1/R_f, conjugation R_{f^-1} = conj(R_f)
        worst = mp.dps
        for f in prim:
            g = compose(fr, f, D)
            worst = min(worst, spare_of(R[g] * R[f] - 1))
            finv = reduce_form(f[0], -f[1], f[2])
            worst = min(worst, spare_of(R[finv] - R[f].conjugate()))
        assert worst >= need, (dK, k, worst)
        # imprimitive strata: units as well (integer polynomial with constant term +-1)
        impr = [f for f in forms if not form_is_primitive(f)]
        impr_note = ""
        if impr:
            co2, sp2 = cert_int_poly([R[f] for f in impr], f"imprimitive R 2alpha={k}")
            assert abs(co2[-1]) == 1
            impr_note = f"; imprimitive stratum ({len(impr)} classes): unit polynomial too (spare {sp2})"
        # KLF on the R_f: sum chi log|R| = -24 L' (odd chi), 0 (even chi)
        klf_note = ""
        if with_klf and len(prim) > 1:
            coords, orders = form_group_data(prim, D)
            chars = all_characters(coords, orders)
            byco = {v: kk for kk, v in coords.items()}
            rn_co = coords[fr]
            chis = []
            for (ks_, cord, ph) in chars:
                if cord == 1:
                    continue
                chis.append({f: cval(ph[f]) for f in prim})
            Lp, M = epstein_Lprime0_forms(prim, D, chis)
            wk = mp.dps
            nodd = 0
            for chi, l in zip(chis, Lp):
                S = sum(chi[f] * log(fabs(R[f])) for f in prim)
                if fabs(chi[fr] + 1) < mpf("0.5"):
                    nodd += 1
                    wk = min(wk, spare_of((S + 24 * l) / max(1, fabs(l))))
                else:
                    wk = min(wk, spare_of(S))
            assert wk >= need, (dK, k, wk)
            klf_note = f"; KLF: {nodd} odd characters with sum chi log|R| = -24 L'(0,chi), even sums 0 (spare {wk}, {M} terms)"
        say(f"  2alpha = {k:2d} (s = {s:+d}, D = {D:5d}, h = {len(prim):2d}, r_alpha = {fr}, r0 = {r0}): "
            f"prod_f (x - R_f) = {co if len(co) <= 5 else str(co[:3])[:-1] + ', ..., ' + str(co[-1]) + ']'} "
            f"(integer, palindromic, constant term {co[-1]}; spare {sp}); twist and conjugation laws (spare {worst})"
            f"{impr_note}{klf_note}")
        results[k] = co
    say(f"  phase C done ({time.time() - t0:.1f} s)")
    return results


# ==========================================================================
# D.  the Euclidean Delta-data, the mass law, the per-class valuations
# ==========================================================================

def mass_exponents(K, n):
    """{p: v_p(|M(n)|)} from the transported mass law: (24/(w_K e_p)) (p^k - 1)/(p - 1) N_e(n/p^k), p not split"""
    out = {}
    for p, k in factor(n).items():
        if K.chi(p) == 1:
            continue
        e = 2 if K.chi(p) == 0 else 1
        val = Fraction(24, K.wK * e) * Fraction(p ** k - 1, p - 1) * K.Ne(n // p ** k)
        assert val.denominator == 1, (K.dK, n, p, val)
        out[p] = int(val)
    return out


def mass_abs(K, n):
    M = 1
    for p, e in mass_exponents(K, n).items():
        M *= p ** e
    return M


def mass_log(K, n):
    s = mpf(0)
    for p, e in mass_exponents(K, n).items():
        s += e * log(mpf(p))
    return s


def per_class_valuation(K, p, k):
    """w_p(k) = 12 (p^k - 1)/(e_p (p - 1) N_e(p^k)) for a non-split p (Fraction)"""
    e = 2 if K.chi(p) == 0 else 1
    return Fraction(12 * (p ** k - 1), e * (p - 1) * K.Ne(p ** k))


def newton_polygon(coeffs, p):
    """slopes of the Newton polygon of sum coeffs[i] x^(deg - i) at p: list of (slope, length)"""
    deg = len(coeffs) - 1
    pts = [(deg - i, vp(c, p)) for i, c in enumerate(coeffs) if c != 0]
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
    out = []
    for (x1, y1), (x2, y2) in zip(hull, hull[1:]):
        out.append((Fraction(y1 - y2, x2 - x1), x2 - x1))
    return out


def euclid_phase(dK, nmax=None, say=print, with_klf=True):
    K = Field(dK)
    t0 = time.time()
    say("=" * 78)
    say(f"PHASE D/E  Euclidean Delta-data over K = {K.name}   (dps = {mp.dps})")
    say("=" * 78)
    need = need_digits()
    ns = list(range(2, (nmax if nmax else 12) + 1))
    sgnOK = 1 if K.t == 0 else -1          # sign of Delta(O_K): q = e^{2 pi i w} < 0 iff Re w = 1/2
    pK = K.ramified()[0] if K.ramified() else None
    rec = {}
    for n in ns:
        cg = ClassGroup(K, n)
        h = cg.h
        assert K.wK // 2 * h == K.Ne(n)
        G = {c: G_value(K, c, n) for c in cg.reps}
        co, sp = cert_int_poly([G[c] for c in cg.reps], f"D_{n}")
        Mn = abs(co[-1])
        assert Mn == mass_abs(K, n), (dK, n, Mn, mass_abs(K, n))
        sgn = 1 if co[-1] * (-1) ** h > 0 else -1
        # certified sign law
        is_ram_power = pK is not None and n > 1 and set(factor(n)) == {pK} and (dK != -4 or n >= 4)
        pred_sign = (sgnOK ** h) * (-1 if is_ram_power else 1)
        assert sgn == pred_sign, (dK, n, sgn, pred_sign)
        # Newton polygons: one slope w_p(k) at every non-split p | n; slope 0 at split p
        slopes = []
        for p, k in factor(n).items():
            npg = newton_polygon(co, p)
            if K.chi(p) == 1:
                assert npg == [] or all(sl == 0 for sl, ln in npg), (dK, n, p, npg)
                slopes.append(f"{p}: 0")
            else:
                assert len(npg) == 1 and npg[0][0] == per_class_valuation(K, p, k) and npg[0][1] == h, (dK, n, p, npg)
                slopes.append(f"{p}^{k}: {per_class_valuation(K, p, k)}")
        # Galois/conjugation law on the data: conj(G_c) = G_{c^-1}
        wc = mp.dps
        for c in cg.reps:
            wc = min(wc, spare_of(G[cg.inv(c)] - G[c].conjugate()))
        assert wc >= need
        fac = factor(Mn) if Mn > 1 else {}
        say(f"  n = {n:2d}: h = {h:2d} (N_e = {K.Ne(n):2d}), D_n in Z[x] (spare {sp}), D_n(0) = (-1)^h M(n), "
            f"|M(n)| = {' '.join(f'{p}^{e}' for p, e in sorted(fac.items())) if fac else '1'} = mass law, "
            f"sign {'+' if sgn > 0 else '-'} = law; Newton slopes {', '.join(slopes) if slopes else '-'}; "
            f"conj(G_c) = G_(c^-1) (spare {wc})")
        rec[n] = (co, h)
        if with_klf and h > 1:
            klf_level(K, cg, G, say=say)
    say(f"  phase D/E done ({time.time() - t0:.1f} s)")
    return rec


# ==========================================================================
# E.  the Kronecker limit formula and the genus characters
# ==========================================================================

def rep_numbers_lattice(K, L, M):
    """r[m] = #{x in Lambda : N(x) = m}, Lambda = Z d + Z (b + a w), m <= M  (norm form of O_K)"""
    d, b, a = L
    r = [0] * (M + 1)
    # x = u d + v (b + a w) = (u d + v b) + v a w ; N = X^2 + t X Y + n0 Y^2 with X = ud + vb, Y = va
    vmax = int((4 * M / (-K.dK)) ** 0.5 / a) + 2
    for v in range(-vmax, vmax + 1):
        Y = v * a
        # X^2 + t X Y + n0 Y^2 <= M
        disc = K.t * K.t * Y * Y - 4 * (K.n0 * Y * Y - M)
        if disc < 0:
            continue
        sq = int(disc ** 0.5) + 2
        for X in range((-K.t * Y - sq) // 2 - 1, (-K.t * Y + sq) // 2 + 2):
            m = X * X + K.t * X * Y + K.n0 * Y * Y
            if 1 <= m <= M and (X - v * b) % d == 0:
                r[m] += 1
    return r


def epstein_Lprime_lattices(K, cg, chars):
    """{ks: L'(0,chi)} by the independent evaluation, discriminant n^2 d_K, w = 2"""
    n = cg.n
    Dn = n * n * K.dK
    alpha = 2 * pi / sqrt(mpf(-Dn))
    M = int((mp.dps + 15) * log(mpf(10)) / alpha) + 10
    reps = {c: rep_numbers_lattice(K, c, M) for c in cg.reps}
    Iv = [None] * (M + 1)
    for m in range(1, M + 1):
        am = alpha * m
        Iv[m] = exp(-am) / am + e1(am)
    out = {}
    for (ks, cord, ph) in chars:
        if cord == 1:
            continue
        tot = mpc(0)
        for c in cg.reps:
            rc = reps[c]
            s = mpf(0)
            for m in range(1, M + 1):
                if rc[m]:
                    s += rc[m] * Iv[m]
            tot += cval(ph[c]) * s
        out[ks] = tot / 2
    return out, M


def gp_run(script, timeout=3600):
    assert GP is not None, "PARI/GP (`gp`) not found on the PATH"
    with tempfile.NamedTemporaryFile("w", suffix=".gp", delete=False) as fh:
        fh.write("default(parisize, 1000000000);\n" + script + "\nquit;\n")
        path = fh.name
    try:
        out = subprocess.run([GP, "-q", "-f", path], capture_output=True, text=True,
                             timeout=timeout, stdin=subprocess.DEVNULL)
    finally:
        os.unlink(path)
    errs = [ln for ln in (out.stdout + out.stderr).splitlines() if "***" in ln and "Warning" not in ln]
    if out.returncode != 0 or errs:
        raise RuntimeError("gp failed:\n" + out.stdout + out.stderr)
    return out.stdout


def parse_tagged(out):
    d = {}
    for line in out.splitlines():
        line = line.strip()
        m = re.match(r"^([A-Z0-9]+) (.*)$", line)
        if m:
            d[m.group(1)] = m.group(2).strip()
    return d


_QUAD_CACHE = {}


def quad_data(disc):
    """(h, w, log eps or None) of the quadratic order of discriminant disc, from PARI
    (fundamental discriminants: unconditional; quadclassunit is GRH-based for h, we
    use qfbclassno for h and quadunit for the fundamental unit)"""
    if disc in _QUAD_CACHE:
        return _QUAD_CACHE[disc]
    script = f"""
default(realprecision, {mp.dps + 20});
d = {disc};
print("H ", qfbclassno(d));
if(d < 0, w = if(d == -3, 6, if(d == -4, 4, 2)); le = 0, w = 2; le = log(abs(quadunit(d))));
print("W ", w);
print("LE ", le);
"""
    dd = parse_tagged(gp_run(script))
    res = (int(dd["H"]), int(dd["W"]), mpf(dd["LE"]))
    _QUAD_CACHE[disc] = res
    return res


def is_disc(x):
    return x % 4 in (0, 1) and x != 0 and x != 1


def fundamental_disc(D):
    """largest fundamental discriminant d with D = f^2 d"""
    f = 1
    D0 = D
    k = 2
    while k * k <= abs(D0):
        while D0 % (k * k) == 0 and is_disc(D0 // (k * k)):
            D0 //= k * k
            f *= k
        k += 1
    return D0, f


def genus_prediction(K, n):
    """for a prime level n = p not dividing d_K: the real field of the genus character,
    d_2 = p (p = 1 mod 4) or p|d_K| (p = 3 mod 4), and d_1 = n^2 d_K / d_2 (fundamental)"""
    p = n
    pstar = p if p % 4 == 1 else -p
    if pstar > 0:
        d2, d1 = pstar, K.dK * pstar
    else:
        d1, d2 = pstar, K.dK * pstar
    assert d1 * d2 == n * n * K.dK
    return d1, d2


def klf_level(K, cg, G, say=print):
    """S_chi = sum chi log|G_c| = -12 L'(0,chi) for every nontrivial chi (independent evaluation);
    genus characters: closed form at prime levels"""
    need = need_digits()
    n = cg.n
    coords, orders = abelian_structure(cg)
    chars = all_characters(coords, orders)
    f = {c: log(fabs(G[c])) for c in cg.reps}
    Lp, M = epstein_Lprime_lattices(K, cg, chars)
    worst = mp.dps
    S = {}
    for (ks, cord, ph) in chars:
        if cord == 1:
            continue
        s = sum(cval(ph[c]) * f[c] for c in cg.reps)
        assert fabs(s.imag) < mpf(10) ** (-need)
        S[ks] = s.real
        l = Lp[ks]
        assert fabs(l.imag) < mpf(10) ** (-need)
        worst = min(worst, spare_of((S[ks] + 12 * l.real) / max(1, fabs(l.real))))
    assert worst >= need, (K.dK, n, worst)
    note = f"      (E) sum chi log|G| = -12 L'(0,chi) for all {len(S)} nontrivial characters (independent Epstein evaluation, {M} terms; spare >= {worst})"
    # genus character at a prime level p not dividing d_K, p >= 3: L'(0,chi_2) = (2h(d1)/w(d1)) h(d2) log eps_{d2}
    if n >= 3 and list(factor(n).items()) == [(n, 1)] and K.chi(n) != 0 and GP is not None:
        d1, d2 = genus_prediction(K, n)
        real = [(ks, cord, ph) for (ks, cord, ph) in chars if cord == 2]
        assert len(real) == 1, (K.dK, n, len(real))
        ks = real[0][0]
        h1, w1, _ = quad_data(d1)
        h2, w2, le2 = quad_data(d2)
        closed = Fraction(2 * h1, w1) * h2 * le2
        l = Lp[ks].real
        spg = spare_of((l - closed) / max(1, fabs(l)))
        assert spg >= need, (K.dK, n, nstr(l, 25), nstr(closed, 25), d1, d2)
        note += (f"\n      genus character: L'(0,chi_2) = (2h({d1})/w) h({d2}) log eps_{d2} = "
                 f"(2*{h1}/{w1})*{h2}*log eps: real field Q(sqrt {fundamental_disc(d2)[0]}) (spare {spg})")
    say(note)
    return S, Lp


# ==========================================================================
# main (phases F-H follow in the second half of the file)
# ==========================================================================


# ==========================================================================
# F.  the Euler system: fiber lemma and norm relations over O_K
# ==========================================================================

def pi12(K, pi):
    """pi^12 as a rational integer (a unit-independent quantity for the ramified prime)"""
    x = K.one
    for _ in range(12):
        x = K.mul(x, pi)
    assert x[1] == 0, (K.dK, pi, x)
    return x[0]


def norm_relation(K, n, l, say=print):
    """the fiber lemma and the norm relation at (n, l); returns the minimal spare digits"""
    m = n * l
    big, small = ClassGroup(K, m), ClassGroup(K, n)
    fiber = {c: [] for c in small.reps}
    for c2 in big.reps:
        fiber[big.project(c2, n)].append(c2)
    chi_l = K.chi(l)
    expected = l if n % l == 0 else l - chi_l
    for c in small.reps:
        assert len(fiber[c]) == expected, (K.dK, n, l, c, len(fiber[c]))
    # the O_n-stable sublattices: 1 + chi(l) of them for l not | n, one for l | n (Lemma 1.2)
    for c in small.reps[:3]:
        subs = K.sublattices(c, l)
        stable = [L for L in subs if K.extend(L, n) == L]
        assert len(stable) == (1 if n % l == 0 else 1 + chi_l), (K.dK, n, l, c, len(stable))
        nonst = [L for L in subs if L not in stable]
        assert all(K.is_primitive(L) and K.extend(L, m) == L for L in nonst)
        assert set(K.canon(L) for L in nonst) == set(fiber[c]), (K.dK, n, l, c)
    worst = mp.dps
    if n % l:
        if chi_l == 1:
            lam = K.prime_above(l)
            lams = [lam, K.conj(lam)]
            case = "split"
        elif chi_l == -1:
            lams = []
            case = "inert"
        else:
            lam = K.prime_above(l)
            lams = [lam]
            case = "ramified"
    else:
        lams = None
        case = "l | n"
    for c in small.reps:
        P = mpf(1)
        for c2 in fiber[c]:
            P *= G_value(K, c2, m)
        Gc = G_value(K, c, n)
        if lams is None:
            n2 = n // l
            Gplus = mpf(1) if n2 == 1 else G_value(K, small.project(c, n2), n2)
            rhs = Gc ** (l + 1) / Gplus
            if l == 2:
                rhs = -rhs
        else:
            rhs = Gc ** (l + 1)
            for lam in lams:
                rhs /= G_value(K, small.twist(c, lam), n)
            if case == "inert":
                rhs *= mpf(l) ** 12
            elif case == "ramified":
                rhs *= pi12(K, lams[0])
            if l == 2:
                rhs = -rhs          # the sign of A(2) = -2^{-24} (all cases)
        worst = min(worst, spare_of(P / rhs - 1))
    need = need_digits()
    assert worst >= need, (K.dK, n, l, worst)
    say(f"  n = {n:2d}, l = {l:2d} -> {m:3d} [{case:8s}] h = {small.h:2d} -> {big.h:3d}, fiber {expected}: "
        f"norm relation at every class (spare {worst})"
        + (f"; pi^12 = {pi12(K, lams[0])}" if case == "ramified" else "") + ("; sign (-1)^[l=2] = -1" if l == 2 else ""))
    return worst


def euler_pairs(K, top=60):
    """>= 20 pairs (n, l) per field covering the four cases, n l <= top"""
    pairs = []
    for l in (2, 3, 5, 7):
        for n in range(2, 13):
            if n * l <= top and n % l:
                pairs.append((n, l))
    for l in (2, 3, 5):
        for n in (l, 2 * l, 3 * l, l * l):
            if n * l <= top and n >= 2:
                pairs.append((n, l))
    pairs = sorted(set(pairs))
    # balance: keep at most 7 per case, all four cases present
    cases = {"split": [], "inert": [], "ramified": [], "l | n": []}
    for (n, l) in pairs:
        if n % l == 0:
            cases["l | n"].append((n, l))
        elif K.chi(l) == 1:
            cases["split"].append((n, l))
        elif K.chi(l) == -1:
            cases["inert"].append((n, l))
        else:
            cases["ramified"].append((n, l))
    out = []
    for k, v in cases.items():
        out += v[:7]
    return sorted(out)


def hecke_recursion(K, n, l, say=print):
    """S_chi(nl) = P_l(chi) S_chi(n) for every nontrivial chi of Pic(O_n), l not | n"""
    need = need_digits()
    cg = ClassGroup(K, n)
    big = ClassGroup(K, n * l)
    coords, orders = abelian_structure(cg)
    chars = all_characters(coords, orders)
    f = {c: log(fabs(G_value(K, c, n))) for c in cg.reps}
    fb = {c: log(fabs(G_value(K, c, n * l))) for c in big.reps}
    proj = {c: big.project(c, n) for c in big.reps}
    worst = mp.dps
    mults = []
    for (ks, cord, ph) in chars:
        if cord == 1:
            continue
        S = sum(cval(ph[c]) * f[c] for c in cg.reps)
        Sb = sum(cval(ph[proj[c]]) * fb[c] for c in big.reps)
        chi_l = K.chi(l)
        if chi_l == -1:
            P = mpc(l + 1)
        else:
            lam = K.prime_above(l)
            z = cval(ph[cg.ideal_class(lam)])
            P = (l + 1 - z - 1 / z) if chi_l == 1 else (l + 1 - z)
        worst = min(worst, spare_of((Sb - P * S) / max(1, fabs(S))))
        Pint, _ = cert_int(P, "P_l")
        mults.append(Pint)
    assert worst >= need, (K.dK, n, l, worst)
    say(f"  Hecke recursion {n} -> {n * l} (l = {l}, {'split' if K.chi(l) == 1 else 'inert' if K.chi(l) == -1 else 'ramified'}): "
        f"S_chi(nl) = P_l(chi) S_chi(n) for all {len(mults)} characters, P_l = {sorted(set(mults))} (spare {worst})")
    return mults


def euler_phase(dK, say=print, top=60):
    K = Field(dK)
    t0 = time.time()
    say("=" * 78)
    say(f"PHASE F  Euler system over K = {K.name}   (dps = {mp.dps})")
    say("=" * 78)
    pairs = euler_pairs(K, top)
    worst = mp.dps
    for (n, l) in pairs:
        worst = min(worst, norm_relation(K, n, l, say=say))
    say(f"  all {len(pairs)} pairs: minimal spare digits {worst}")
    # Hecke recursion at the first level with h >= 2, three primes
    n0 = next(n for n in range(2, 20) if ClassGroup(K, n).h >= 2)
    for l in (2, 3, 5):
        if n0 % l == 0 or n0 * l > top:
            continue
        hecke_recursion(K, n0, l, say=say)
    say(f"  phase F done ({time.time() - t0:.1f} s)")
    return len(pairs), worst


# ==========================================================================
# G.  the full Robert index with PARI/GP
# ==========================================================================

def torsion_prediction(K, n):
    """w(H_n) = |mu(H_n)| from the conductor lemma (Lemma of other-fields.md)"""
    z3 = (K.dK == -3) or (n % 3 == 0)
    z4 = (K.dK == -4) or (K.dK == -8 and n % 2 == 0) or (n % 4 == 0)
    z8 = (K.dK == -4 and n % 4 == 0) or (K.dK == -8 and n % 2 == 0) or (n % 8 == 0)
    w = 2
    if z4:
        w = 4
    if z8:
        w = 8
    if z3:
        w *= 3
    return w, z4, z8, z3


def parse_int_vector(s):
    s = s.strip()
    assert s.startswith("[") and s.endswith("]"), s
    body = s[1:-1].strip()
    if not body:
        return []
    return [int(t) for t in body.replace(",", " ").split()]


def parse_int_matrix(s):
    s = s.strip()
    m = re.match(r"^Mat\((.*)\)$", s)
    if m:
        inner = m.group(1).strip()
        if inner.startswith("["):
            return [parse_int_vector(inner)]
        return [[int(inner)]]
    assert s.startswith("[") and s.endswith("]"), s
    body = s[1:-1].strip()
    if not body:
        return []
    return [[int(t) for t in row.replace(",", " ").split()] for row in body.split(";")]


def smith_invariants(rows):
    from sympy import Matrix, ZZ
    from sympy.matrices.normalforms import smith_normal_form
    S = smith_normal_form(Matrix(rows), domain=ZZ)
    inv = [abs(int(S[i, i])) for i in range(min(S.shape))]
    return [x for x in inv if x != 1]


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def factors_through(cg, ph, m):
    if m == 1:
        return all(p == 0 for p in ph.values())
    img_one = cg.project(cg.one, m)
    for c in cg.reps:
        if cg.project(c, m) == img_one and ph[c] != 0:
            return False
    return True


def primitive_level(cg, ph):
    levels = [m for m in divisors(cg.n) if factors_through(cg, ph, m)]
    for m1 in levels:
        for m2 in levels:
            assert gcd(m1, m2) in levels, ("levels not gcd-closed", cg.n, levels)
    return min(levels), levels


def character_on_level(cg, ph, m):
    K = cg.K
    cgm = ClassGroup(K, m) if m >= 2 else None
    if cgm is None:
        return None, None
    out = {}
    for c in cg.reps:
        cm = cg.project(c, m)
        if cm in out:
            assert out[cm] == ph[c]
        else:
            out[cm] = ph[c]
    assert set(out) == set(cgm.reps)
    return cgm, out


def multiplier_chain(cg, ph, m):
    """C_chi(0) = prod_l a_{k_l}: a_1 = P_l(chi) [l not | m] or l + 1 [l | m];
    a_{k+1} = (l+1) a_k - m_k a_{k-1}, m_1 = l - chi_K(l) [l not | m] or l [l | m], m_k = l (k >= 2);
    P_l = l + 1 - chi(l) - chi(l)^{-1} (split), l + 1 (inert), l + 1 - chi(p) (ramified)"""
    K = cg.K
    n = cg.n
    fac = factor(n // m)
    C = mpf(1)
    desc = []
    for l, k in sorted(fac.items()):
        mprime = n // l ** k
        if mprime % l == 0:
            a1 = mpf(l + 1)
            m1 = l
            case = f"l={l} | m: a_1 = {l + 1}"
        else:
            chi_l = K.chi(l)
            if chi_l == -1:
                a1 = mpf(l + 1)
                case = f"l={l} inert: a_1 = {l + 1}"
            else:
                cgp, php = character_on_level(cg, ph, mprime)
                assert cgp is not None
                lam = K.prime_above(l)
                z = cval(php[cgp.ideal_class(lam)])
                if chi_l == 0:
                    a1 = (l + 1 - z).real
                    case = f"l={l} ramified: a_1 = l+1-chi(p) = {nstr(a1, 8)}"
                else:
                    a1 = (l + 1 - z - 1 / z).real
                    case = f"l={l} split: a_1 = l+1-2Re chi(l) = {nstr(a1, 8)}"
            m1 = l - chi_l
        a_prev, a_cur = mpf(1), a1
        for j in range(2, k + 1):
            mk = m1 if j == 2 else l
            a_prev, a_cur = a_cur, (l + 1) * a_cur - mk * a_prev
        C *= a_cur
        desc.append(case + (f", a_{k} = {nstr(a_cur, 8)}" if k > 1 else ""))
    return C, "; ".join(desc) if desc else "primitive"


def regulator_matrix(cg, f):
    others = [c for c in cg.reps if c != cg.one]
    Mx = matrix(len(others), len(others))
    for i, a in enumerate(others):
        ainv = cg.inv(a)
        base = f[ainv]
        for j, c in enumerate(others):
            Mx[i, j] = 2 * (f[cg.prod(ainv, c)] - base)
    return Mx


def gp_field_script(K, n, prec, certdeg, Dn):
    return f"""default(realprecision, {prec});
dK = {K.dK}; n = {n};
t0 = getabstime();
P = polclass(n^2*dK, 0, 'y);
Q = polredbest(polcompositum(y^2 - dK, P)[1]);
print("TPOL ", (getabstime()-t0)/1000.);
t0 = getabstime();
bnf = bnfinit(Q, 1);
print("TBNF ", (getabstime()-t0)/1000.);
print("POL ", Q);
print("DEG ", poldegree(Q));
print("H ", bnf.no);
print("CYC ", bnf.cyc);
print("W ", bnf.tu[1]);
print("Z4 ", #nfroots(bnf.nf, x^2+1));
print("Z3 ", #nfroots(bnf.nf, x^2+3));
print("S2 ", #nfroots(bnf.nf, x^2-2));
t0 = getabstime();
cert = if(poldegree(Q) <= {certdeg}, bnfcertify(bnf), -1);
print("CERT ", cert);
print("TCERT ", (getabstime()-t0)/1000.);
Dn = Pol({Dn}, x);
r = nfroots(bnf.nf, Dn);
print("NROOTS ", #r);
E = matrix(#r-1, #bnf.fu);
T = vector(#r-1);
for(i=2, #r, e = bnfisunit(bnf, r[i]/r[1]); for(j=1, #bnf.fu, E[i-1,j] = e[j]); T[i-1] = lift(e[#e]));
print("EXP ", E);
print("TORS ", T);
print("EDET ", matdet(E));
print("REG ", bnf.reg);
"""


def gp_cubic_script(K, prec, cubic_theta_u):
    return f"""default(realprecision, {prec});
cu = Pol({cubic_theta_u}, y);
L3 = polredbest(cu);
print("L3POL ", L3);
bL = bnfinit(L3, 1);
print("L3H ", bL.no);
print("L3REG ", bL.reg);
print("L3FU ", lift(bL.fu[1]));
print("L3CERT ", bnfcertify(bL));
rt = nfroots(bL.nf, subst(cu, y, x));
print("L3NROOTS ", #rt);
e = bnfisunit(bL, rt[1]);
print("L3EXP ", e[1], " ", lift(e[2]));
F = polredbest(polcompositum(y^2 - ({K.dK}), cu)[1]);
print("FPOL ", F);
bF = bnfinit(F, 1);
print("FH ", bF.no);
print("FREG ", bF.reg);
print("FW ", bF.tu[1]);
print("FCERT ", bnfcertify(bF));
rF = nfroots(bF.nf, subst(cu, y, x));
print("FNROOTS ", #rF);
for(i=1, #rF, e = bnfisunit(bF, rF[i]); print("FEXP", i, " ", e[1], " ", e[2], " ", lift(e[3])));
"""


def gp_quad_script(prec, quad):
    return f"""default(realprecision, {prec});
qu = Pol({quad}, y);
K2 = polredbest(qu);
print("K2POL ", K2);
bQ = bnfinit(K2, 1);
print("K2FU ", lift(bQ.fu[1]));
print("K2DISC ", nfdisc(K2));
print("K2REG ", bQ.reg);
rq = nfroots(bQ.nf, subst(qu, y, x));
print("K2NROOTS ", #rq);
e = bnfisunit(bQ, rq[1]);
print("K2EXP ", e[1], " ", lift(e[2]));
"""


def iroot(M, k):
    r = int(round(M ** (1.0 / k)))
    while r ** k < M:
        r += 1
    while r ** k > M:
        r -= 1
    return r if r ** k == M else None


def index_level(dK, n, prec_pari=120, certdeg=16, say=print, layers=True):
    K = Field(dK)
    t0 = time.time()
    cg = ClassGroup(K, n)
    h = cg.h
    coords, orders = abelian_structure(cg)
    chars = all_characters(coords, orders)
    need = need_digits()
    say("=" * 78)
    say(f"PHASE G  K = {K.name}, level n = {n}: Pic(O_n) = {' x '.join('Z/%d' % o for o in orders)}, h = {h}, "
        f"[H_n : Q] = {2 * h}   (dps = {mp.dps})")
    say("=" * 78)
    rec = {"dK": dK, "n": n, "h": h, "orders": orders}
    G = {c: G_value(K, c, n) for c in cg.reps}
    f = {c: log(fabs(G[c])) for c in cg.reps}
    Dn, spD = cert_int_poly([G[c] for c in cg.reps], f"D_{n}")
    assert abs(Dn[-1]) == mass_abs(K, n)
    # (V1) Dedekind determinant
    S = {}
    for (ks, cord, ph) in chars:
        if cord == 1:
            continue
        s = sum(cval(ph[c]) * f[c] for c in cg.reps)
        assert fabs(s.imag) < mpf(10) ** (-need)
        S[ks] = s.real
        assert fabs(S[ks]) > mpf("0.01"), (dK, n, ks, "S_chi vanishes?")
    prodS = mpf(1)
    for ks in S:
        prodS *= S[ks]
    Mx = regulator_matrix(cg, f)
    dM = det(Mx) if h > 1 else mpf(1)
    rhs = 2 ** (h - 1) * prodS
    spV1 = spare_of((dM - rhs) / max(1, fabs(rhs)))
    assert spV1 >= need, (dK, n, nstr(dM, 20), nstr(rhs, 20))
    RV = fabs(dM)
    say(f"  (V1) det(2(f(a^-1 c) - f(a^-1))) = 2^(h-1) prod S_chi (signed): spare {spV1};  R(V_n) = {nstr(RV, 20)}")
    # (V2)/(V3) independent L' and multipliers
    Lp, M = epstein_Lprime_lattices(K, cg, chars)
    worst2 = mp.dps
    prodC = mpf(1)
    prodL = mpf(1)
    Cinfo = {}
    for (ks, cord, ph) in chars:
        if cord == 1:
            continue
        l = Lp[ks].real
        worst2 = min(worst2, spare_of((S[ks] + 12 * l) / max(1, fabs(l))))
        prodL *= fabs(l)
        m, levels = primitive_level(cg, ph)
        C, desc = multiplier_chain(cg, ph, m)
        if m < n:
            cgm, phm = character_on_level(cg, ph, m)
            Lm, _ = epstein_Lprime_lattices(K, cgm, [((), 2, phm)])
            spC = spare_of((l - C * Lm[()].real) / max(1, fabs(l)))
            assert spC >= need, (dK, n, ks, m, nstr(l / Lm[()].real, 20), nstr(C, 20))
        else:
            spC = None
        Cinfo[ks] = (cord, m, C, desc, spC)
        prodC *= C
    assert worst2 >= need
    prodC_int, spPC = cert_int(prodC, "prod C")
    say(f"  (V2) S_chi = -12 L'(0,chi) (independent evaluation, {M} terms): spare >= {worst2}")
    say(f"  (V3) primitive levels / multipliers: " + "; ".join(
        f"chi{ks} (ord {cord}) from level {m}, C = {nstr(C, 6)}" + (f" [{desc}]" if m < n else "")
        for ks, (cord, m, C, desc, spC) in sorted(Cinfo.items())) + f";  prod C = {prodC_int}")
    # (V4) PARI
    out = gp_run(gp_field_script(K, n, prec_pari, certdeg, Dn))
    d = parse_tagged(out)
    hH, w, cyc, cert, Rh, deg = int(d["H"]), int(d["W"]), parse_int_vector(d["CYC"]), int(d["CERT"]), mpf(d["REG"]), int(d["DEG"])
    assert deg == 2 * h and int(d["NROOTS"]) == h, (dK, n, deg, d["NROOTS"])
    w_pred, z4, z8, z3 = torsion_prediction(K, n)
    assert w == w_pred, (dK, n, w, w_pred)
    assert (int(d["Z4"]) > 0) == z4 and (int(d["Z3"]) > 0) == z3 and (int(d["S2"]) > 0 and z4) == z8, (dK, n, d["Z4"], d["Z3"], d["S2"], z4, z8, z3)
    grh = "unconditional (bnfcertify = 1)" if cert == 1 else "GRH-conditional"
    index_real = RV / Rh
    index_int, spI = cert_int(index_real, f"index {dK} n={n}", need=max(40, need))
    pred = Fraction(K.wK * 24 ** (h - 1), w) * hH * prodC_int
    assert pred == index_int, (dK, n, index_int, pred)
    say(f"  (V4) PARI: H_n = Q[y]/({d['POL'][:60]}{'...' if len(d['POL']) > 60 else ''}); h(H_n) = {hH}, Cl = {cyc or '[]'}, "
        f"w = {w} (= torsion lemma: zeta_4 {z4}, zeta_8 {z8}, zeta_3 {z3}), R = {nstr(Rh, 20)} -- {grh}")
    say(f"       index 24^(h-1) prod|L'| / R = {index_int} (certified, spare {spI}) = (w_K 24^(h-1)/w) h(H_n) prod C "
        f"= ({K.wK} * 24^{h - 1} / {w}) * {hH} * {prodC_int}: OK")
    # (V5) exact index
    E = parse_int_matrix(d["EXP"]) if h > 1 else []
    edet = int(d["EDET"]) if h > 1 else 1
    assert abs(edet) == index_int, (dK, n, edet, index_int)
    snf = smith_invariants(E) if h > 1 else []
    sat = 1
    quo = []
    for dd in snf:
        sat *= gcd(dd, 24)
        if dd // gcd(dd, 24) > 1:
            quo.append(dd // gcd(dd, 24))
    sat_index = index_int // sat
    say(f"  (V5) exact (bnfisunit): |det E| = {abs(edet)} = index; O^x/mu V = {' x '.join('Z/%d' % x for x in snf) or 'trivial'}; "
        f"24th-root saturation: [O^x : W_n] = {sat_index} ({'= h prod C' if sat_index == hH * prodC_int else 'DIFFERENT from h prod C'}), "
        f"O^x/W_n = {' x '.join('Z/%d' % x for x in quo) or 'trivial'}")
    rec.update({"hH": hH, "w": w, "cyc": cyc, "cert": cert, "index": index_int, "prodC": prodC_int, "snf": snf,
                "sat_index": sat_index, "spI": spI, "pol": d["POL"], "D": Dn})
    # (V6) layers
    if layers and h > 1:
        Mn = mass_abs(K, n)
        r2 = iroot(Mn, 2)
        if r2 is not None:
            for (ks, cord, ph) in chars:
                if cord != 2:
                    continue
                ker = [c for c in cg.reps if ph[c] == 0]
                oth = [c for c in cg.reps if ph[c] != 0]
                th = mpc(1)
                for c in ker:
                    th *= G[c]
                th2 = mpc(1)
                for c in oth:
                    th2 *= G[c]
                co, sp = cert_int_poly([th / r2, th2 / r2], f"quadratic layer chi{ks}")
                assert abs(co[-1]) == 1
                dq = parse_tagged(gp_run(gp_quad_script(prec_pari, co)))
                e2 = int(dq["K2EXP"].split()[0])
                Rq = mpf(dq["K2REG"])
                e_pred, spe = cert_int(6 * Lp[ks].real / Rq, "6 L'/log eps")
                assert e_pred == abs(e2), (dK, n, ks, e2, e_pred)
                say(f"  (V6) quadratic layer chi{ks}: theta^(2) = prod_ker G/|M|^(1/2) = +-eps^{e2} in Q(sqrt {dq['K2DISC']}) "
                    f"(minimal polynomial {co}, eps = PARI's fundamental unit): index |e| = {abs(e2)} = 6 L'(0,chi)/log eps (spare {spe})")
                rec.setdefault("quad", []).append((ks, co, int(dq["K2DISC"]), e2))
        cub = [(ks, ph) for (ks, cord, ph) in chars if cord == 3]
        M13 = iroot(Mn, 3)
        if cub and M13 is not None:
            ks, ph = cub[0]
            cosets = {}
            for c in cg.reps:
                cosets.setdefault(ph[c], []).append(c)
            assert len(cosets) == 3
            prods = []
            for key in sorted(cosets):
                p = mpc(1)
                for c in cosets[key]:
                    p *= G[c]
                prods.append(p)
            co_u, sp_u = cert_int_poly([p / M13 for p in prods], "theta_u cubic")
            assert abs(co_u[-1]) == 1, (dK, n, co_u)
            dc = parse_tagged(gp_run(gp_cubic_script(K, prec_pari, co_u)))
            eL, hL, hF, wF = int(dc["L3EXP"].split()[0]), int(dc["L3H"]), int(dc["FH"]), int(dc["FW"])
            certL = int(dc["L3CERT"]) == 1 and int(dc["FCERT"]) == 1
            rows = [[int(t) for t in dc[f"FEXP{i}"].split()] for i in (1, 2, 3)]
            assert all(rows[0][j] + rows[1][j] + rows[2][j] == 0 for j in (0, 1)), rows
            d1 = [rows[0][j] - rows[1][j] for j in (0, 1)]
            d2 = [rows[1][j] - rows[2][j] for j in (0, 1)]
            idxF = abs(d1[0] * d2[1] - d1[1] * d2[0])
            C3i, _ = cert_int(Cinfo[ks][2], "C_chi3(0)")
            assert abs(eL) == 8 * hL * C3i, (dK, n, eL, hL, C3i)
            assert wF == K.wK, (dK, n, wF)
            assert idxF == Fraction(K.wK * 24 ** 2, wF) * hF * C3i ** 2 == 576 * hF * C3i ** 2, (dK, n, idxF, hF, C3i)
            say(f"  (V6) cubic layer chi{ks}: theta_u = theta/{M13}, minimal polynomial {co_u}; L_3 = Q[y]/({dc['L3POL']}): "
                f"h = {hL}, theta_u = +-eta^{eL} -> [O_L3^x : <-1, theta_u>] = {abs(eL)} = 8 h_L3 C(0); "
                f"sextic F = K(theta_u): h_F = {hF}, w_F = {wF}, [O_F^x : mu V^A] = {idxF} = 576 h_F C^2 "
                f"({'unconditional' if certL else 'GRH'})")
            rec["cubic"] = (co_u, hL, eL, hF, idxF, dc["L3POL"])
    say(f"  level done in {time.time() - t0:.1f} s (PARI: polynomial {float(d['TPOL']):.2f} s, bnfinit {float(d['TBNF']):.2f} s, bnfcertify {float(d['TCERT']):.2f} s)")
    return rec


# ==========================================================================
# H.  class number two: K = Q(sqrt -5)
# ==========================================================================

class ClassGroupH2:
    """Pic(O_n) for K = Q(sqrt -5): classes of proper O_n-lattices Lambda with O_K Lambda in
    {O_K, L}, L = p_2 = (2, 1 + sqrt -5), represented by (span type, HNF sublattice), modulo +-1."""

    def __init__(self, K, n):
        assert K.dK == -20
        self.K, self.n = K, n
        self.L = hnf([(2, 0), (1, 1)])          # p_2 = Z 2 + Z (1 + w)
        assert K.span_OK(self.L) == self.L and K.index(self.L) == 2
        reps = []
        for typ, base in ((0, (1, 0, 1)), (1, self.L)):
            v1, v2 = K.vecs(base)
            for a in divisors(n):
                d = n // a
                for b in range(d):
                    Lam = hnf([K.smul(d, v1), K.add(K.smul(b, v1), K.smul(a, v2))])
                    if K.span_OK(Lam) == base:
                        reps.append((typ, Lam))
        self.reps = sorted(set(reps))
        self.h = len(self.reps)
        assert self.h == 2 * K.Ne(n), (n, self.h, K.Ne(n))
        self.one = (0, (1, 0, n))
        assert self.one in self.reps
        self.idx = {c: k for k, c in enumerate(self.reps)}

    def normalize(self, Lam):
        """(span type, lattice) of a proper O_n-lattice Lam (any O_K-span), scaled into O_K or L"""
        K = self.K
        span = K.span_OK(Lam)
        # span is an O_K-ideal; principal iff its class is trivial: test by finding a generator g with
        # (g) = span, i.e. span = g O_K: try elements of small norm equal to index(span)
        idx = K.index(span)
        for typ, base in ((0, (1, 0, 1)), (1, self.L)):
            # is span = g * base for some g?  N(g) = idx / index(base)
            if idx % K.index(base):
                continue
            Ng = idx // K.index(base)
            for b in range(0, int(Ng ** 0.5) + 2):
                for a in range(-int(Ng ** 0.5) - 2, int(Ng ** 0.5) + 3):
                    g = (a, b)
                    if K.norm(g) == Ng and K.scale(base, g) == span:
                        # divide Lam by g: Lam = g * Lam'  <=>  Lam' = conj(g) Lam / N(g)
                        Lp = hnf([K.mul(v, K.conj(g)) for v in K.vecs(Lam)])
                        dd, bb, aa = Lp
                        assert dd % Ng == 0 and bb % Ng == 0 and aa % Ng == 0
                        Lp = (dd // Ng, bb // Ng, aa // Ng)
                        assert K.span_OK(Lp) == base
                        return (typ, Lp)
        raise ValueError("no generator found")

    def cls(self, Lam):
        c = self.normalize(Lam)
        assert c in self.idx, (self.n, Lam, c)
        return c

    def prod(self, c1, c2):
        return self.cls(self.K.lmul(c1[1], c2[1]))

    def inv(self, c):
        return self.cls(self.K.conj_lat(c[1]))

    def order(self, c):
        k, r = 1, c
        while r != self.one:
            r = self.prod(r, c)
            k += 1
        return k

    def project(self, c, n2):
        assert n2 == 1
        return c[0]


def h2_phase(say=print, prec_pari=120, levels=(2, 3)):
    K = Field(-20)
    t0 = time.time()
    say("=" * 78)
    say(f"PHASE H  class number two: K = Q(sqrt -5), d_K = -20, h_K = 2, w_K = 2   (dps = {mp.dps})   [experimental]")
    say("=" * 78)
    need = need_digits()
    recs = {}
    # level 1: the two classes of O_K, the Epstein L'(0, chi_K) of discriminant -20 (w = 2)
    forms1 = [(1, 0, 5), (2, 2, 3)]
    chiK = {(1, 0, 5): mpc(1), (2, 2, 3): mpc(-1)}
    LpK, _ = epstein_Lprime0_forms(forms1, -20, [chiK])
    LpK = LpK[0].real
    for n in levels:
        cg = ClassGroupH2(K, n)
        h = cg.h
        # (H1) Stange: the single cusp (lattices in O_K) sees exactly ker(Pic(O_n) -> Pic(O_K)),
        #      a subgroup of index h_K = 2 with N_e(n) = h/2 elements
        ker = [c for c in cg.reps if c[0] == 0]
        assert len(ker) == K.Ne(n) == h // 2
        for c1 in ker:
            for c2 in ker:
                assert cg.prod(c1, c2)[0] == 0
        coords, orders = abelian_structure(cg)
        chars = all_characters(coords, orders)
        G = {c: G_value(K, c[1], n, ref=None if c[0] == 0 else cg.L) for c in cg.reps}
        f = {c: log(fabs(G[c])) for c in cg.reps}
        Dn, spD = cert_int_poly([G[c] for c in cg.reps], f"D_{n} (h_K = 2)")
        wc = min(spare_of(G[cg.inv(c)] - G[c].conjugate()) for c in cg.reps)
        assert wc >= need
        say(f"  n = {n}: Pic(O_n) = {' x '.join('Z/%d' % o for o in orders)}, h = {h}; single cusp sees the kernel "
            f"({len(ker)} classes = N_e(n) = h/h_K); two-cusp Delta-data G_c = n^12 Delta(Lambda)/Delta(O_K Lambda): "
            f"D_n = {Dn} in Z[x] (spare {spD}), conj(G_c) = G_(c^-1) (spare {wc})")
        # (H2) KLF: S_chi = -12 L'(0, chi) unless chi is pulled back from Pic(O_K); then
        #      S_chi = -12 L'(0,chi) + 12 |ker| L'(0, chi_K)
        # independent L'(0,chi): the forms of the classes (norm forms of the lattices / N(span))
        forms = {}
        for c in cg.reps:
            d, b, a = c[1]
            typ = c[0]
            X1, Y1 = d, 0
            X2, Y2 = b, a
            A = K.norm((X1, Y1))
            Cc = K.norm((X2, Y2))
            B = K.norm((X1 + X2, Y1 + Y2)) - A - Cc
            div = 1 if typ == 0 else 2
            assert A % div == 0 and B % div == 0 and Cc % div == 0
            fm = (A // div, B // div, Cc // div)
            assert fm[1] ** 2 - 4 * fm[0] * fm[2] == n * n * K.dK, (c, fm)
            forms[c] = fm
        chis = []
        keys = []
        for (ks, cord, ph) in chars:
            if cord == 1:
                continue
            chis.append({forms[c]: cval(ph[c]) for c in cg.reps})
            keys.append((ks, cord, ph))
        # representation numbers of the forms need distinct forms per class: check
        assert len(set(forms.values())) == h, "non-distinct forms"
        Lp, M = epstein_Lprime0_forms(list(forms.values()), n * n * K.dK,
                                      [{fm: chi[fm] for fm in forms.values()} for chi in chis])
        kerpi = len(ker)
        worst = mp.dps
        base_info = []
        for (ks, cord, ph), l in zip(keys, Lp):
            S = sum(cval(ph[c]) * f[c] for c in cg.reps)
            assert fabs(S.imag) < mpf(10) ** (-need)
            S = S.real
            l = l.real
            from_base = all(ph[c] == ph[cg.one] for c in ker)
            if from_base:
                # chi = chi_K o pi: correction 12 |ker| L'(0, chi_K); and L'(0,chi) = C L'(0,chi_K)
                sp = spare_of((S + 12 * l - 12 * kerpi * LpK) / max(1, fabs(l)))
                Cb = l / LpK
                Cb_int, _ = cert_int(Cb, "base multiplier")
                base_info.append((ks, Cb_int, sp))
            else:
                sp = spare_of((S + 12 * l) / max(1, fabs(l)))
            worst = min(worst, sp)
        assert worst >= need, (n, worst)
        say(f"      (H2) sum chi log|G| = -12 L'(0,chi) for the {len(keys) - len(base_info)} characters not from Pic(O_K); "
            f"for the pulled-back character: sum chi log|G| = -12 L'(0,chi) + 12 |ker pi| L'(0,chi_K) with "
            f"|ker pi| = {kerpi}, L'(0,chi) = C L'(0,chi_K), C = {[b[1] for b in base_info]} (spare >= {worst}, {M} terms)")
        # (H3) the index with PARI: R(V_n) = |det| of the Dedekind matrix; H_n from polclass(n^2 dK)
        S = {}
        for (ks, cord, ph) in keys:
            S[ks] = sum(cval(ph[c]) * f[c] for c in cg.reps).real
        Mx = regulator_matrix(cg, f)
        dM = det(Mx)
        prodS = mpf(1)
        for ks in S:
            prodS *= S[ks]
        assert spare_of((dM - 2 ** (h - 1) * prodS) / max(1, fabs(prodS))) >= need
        RV = fabs(dM)
        out = gp_run(gp_field_script(K, n, prec_pari, 16, Dn))
        d = parse_tagged(out)
        hH, w, cyc, cert, Rh, deg = int(d["H"]), int(d["W"]), parse_int_vector(d["CYC"]), int(d["CERT"]), mpf(d["REG"]), int(d["DEG"])
        assert deg == 2 * h and int(d["NROOTS"]) == h
        index_real = RV / Rh
        index_int, spI = cert_int(index_real, f"index h_K=2 n={n}", need=max(40, need))
        edet = int(d["EDET"])
        assert abs(edet) == index_int
        # prediction: 24^(h-1) (w_K/h_K) (h_H/w_H) prod_chi C_chi(0) prod_base |1 - |ker|/C|
        prodC = 1
        corr = Fraction(1)
        for (ks, Cb, sp) in base_info:
            prodC *= Cb
            corr *= Fraction(abs(Cb - kerpi), Cb)
        pred = Fraction(24 ** (h - 1) * K.wK, K.hK * w) * hH * prodC * corr
        say(f"      (H3) PARI: H_n of degree {deg}, h(H_n) = {hH}, Cl = {cyc or '[]'}, w = {w}, R = {nstr(Rh, 15)} "
            f"({'bnfcertify = 1' if cert == 1 else 'GRH'}); exact index [O^x : mu V_n] = {index_int} (|det E| = {abs(edet)}, spare {spI}) "
            f"= 24^(h-1) (w_K/h_K)(h_H/w_H) prod C |1 - |ker|/C| = {pred}: {'OK' if pred == index_int else 'MISMATCH'}")
        assert pred == index_int, (n, index_int, pred)
        recs[n] = {"h": h, "orders": orders, "hH": hH, "w": w, "index": index_int, "C": [b[1] for b in base_info],
                   "D": Dn, "cert": cert, "snf": smith_invariants(parse_int_matrix(d["EXP"]))}
    say(f"  phase H done ({time.time() - t0:.1f} s)")
    return recs


# ==========================================================================
# selftest, records, main
# ==========================================================================

# regression records: (dK, n) -> (h, w(H_n), prod C, index, Smith invariants)
INDEX_LEVELS = {-8: [2, 3, 4, 5, 6], -3: [4, 5, 6, 7, 8], -7: [3, 4, 5, 6], -11: [2, 3, 4, 5]}
INDEX_RECORD = {
    (-8, 2): (2, 8, 1, 6, [6]),
    (-8, 3): (2, 6, 1, 8, [8]),
    (-8, 4): (4, 8, 3, 10368, [6, 24, 72]),
    (-8, 5): (6, 2, 1, 7962624, [24, 24, 24, 24, 24]),
    (-8, 6): (4, 24, 24, 27648, [2, 24, 576]),
    (-3, 4): (2, 12, 1, 12, [12]),
    (-3, 5): (2, 6, 1, 24, [24]),
    (-3, 6): (3, 6, 1, 576, [24, 24]),
    (-3, 7): (2, 6, 1, 24, [24]),
    (-3, 8): (4, 24, 3, 10368, [6, 24, 72]),
    (-7, 3): (4, 6, 1, 4608, [8, 24, 24]),
    (-7, 4): (2, 4, 1, 12, [12]),
    (-7, 5): (6, 2, 1, 7962624, [24, 24, 24, 24, 24]),
    (-7, 6): (4, 6, 45, 207360, [8, 72, 360]),
    (-11, 2): (3, 2, 1, 576, [24, 24]),
    (-11, 3): (2, 6, 1, 8, [8]),
    (-11, 4): (6, 4, 9, 35831808, [12, 24, 24, 72, 72]),
    (-11, 5): (4, 2, 1, 27648, [24, 24, 48]),
    (-4, 5): (2, 4, 1, 24, [24]),
}


def selftest(quick=False):
    T0 = time.time()
    print("=" * 78)
    print("other_fields.py --selftest   (mpmath; PARI/GP %s)" % ("found" if GP else "MISSING"))
    print("=" * 78)
    times = {}
    for dK in FIELDS + [ANCHOR]:
        t = time.time()
        mp.dps = 60
        geometry_phase(dK)
        times[(dK, "geometry")] = time.time() - t
    for dK in FIELDS:
        t = time.time()
        mp.dps = 200
        units_phase(dK, kmax=14 if quick else 18)
        times[(dK, "units")] = time.time() - t
    for dK in FIELDS:
        t = time.time()
        mp.dps = 150
        euclid_phase(dK, nmax=7 if quick else 9)
        times[(dK, "euclid")] = time.time() - t
    for dK in FIELDS:
        t = time.time()
        mp.dps = 80
        euler_phase(dK, top=40 if quick else 60)
        times[(dK, "euler")] = time.time() - t
    recs = []
    if GP is not None:
        for dK in FIELDS + [ANCHOR]:
            t = time.time()
            mp.dps = 150
            for n in (INDEX_LEVELS.get(dK, [5]) if not quick else INDEX_LEVELS.get(dK, [5])[:2]):
                r = index_level(dK, n)
                recs.append(r)
                if (dK, n) in INDEX_RECORD:
                    h, w, pC, idx, snf = INDEX_RECORD[(dK, n)]
                    assert (r["h"], r["w"], r["prodC"], r["index"]) == (h, w, pC, idx), (dK, n, r["h"], r["w"], r["prodC"], r["index"])
                    assert sorted(r["snf"]) == sorted(snf), (dK, n, r["snf"])
            times[(dK, "index")] = time.time() - t
        t = time.time()
        mp.dps = 150
        h2_phase()
        times[(-5, "h2")] = time.time() - t
    print("=" * 78)
    print("SUMMARY  [O_{H_n}^x : mu(H_n) V_n] = (w_K 24^(h-1)/w) h(H_n) prod C_chi(0)")
    print("=" * 78)
    print(f"{'dK':>4} {'n':>3} {'Pic':>10} {'h':>3} {'w':>3} {'h(H_n)':>7} {'Cl(H_n)':>9} {'prodC':>6} {'index':>22} {'sat':>8} status")
    for r in recs:
        st = "certified" if r["cert"] == 1 else "GRH"
        print(f"{r['dK']:>4} {r['n']:>3} {'x'.join('Z/%d' % o for o in r['orders']):>10} {r['h']:>3} {r['w']:>3} "
              f"{r['hH']:>7} {str(r['cyc'] or '[]'):>9} {r['prodC']:>6} {r['index']:>22} {r['sat_index']:>8} {st}")
    print("timings (s): " + ", ".join(f"{k[1]}[{k[0]}]={v:.1f}" for k, v in times.items()))
    print(f"ALL CHECKS PASSED   ({time.time() - T0:.0f}s)")


def main(argv):
    args = [a for a in argv[1:] if not a.startswith("--")]
    flags = [a for a in argv[1:] if a.startswith("--")]
    if "--selftest" in flags or not args:
        selftest(quick="--quick" in flags)
        return
    what = args[0]
    if what == "geometry":
        mp.dps = 60
        geometry_phase(int(args[1]), kmax=int(args[2]) if len(args) > 2 else None)
    elif what == "units":
        mp.dps = 200
        units_phase(int(args[1]), kmax=int(args[2]) if len(args) > 2 else 18)
    elif what == "euclid":
        mp.dps = 150
        euclid_phase(int(args[1]), nmax=int(args[2]) if len(args) > 2 else 9)
    elif what == "euler":
        mp.dps = 80
        euler_phase(int(args[1]))
    elif what == "index":
        mp.dps = 150
        for n in args[2:]:
            index_level(int(args[1]), int(n))
    elif what == "h2":
        mp.dps = 150
        h2_phase()
    else:
        print(__doc__)


if __name__ == "__main__":
    main(sys.argv)
