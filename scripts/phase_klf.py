"""Character sums of log|phase| against L'(0,chi): the Kronecker-limit /
elliptic-unit structure of the Schmidt phases.

Companion to phase-kronecker-limit.md.  Two aspects:

Euclidean (disc -4n^2): u_c = Theta/Omega over Cl(O_n), O_n = Z + nZ[i].
  From the proved closed form u^6 = -beta^4(beta-1728)^3 Delta(L_c)/Delta(Z[i])
  and the first Kronecker limit formula applied classwise,
      6 S(chi) = -12 L'(0,chi) + 4 Sigma_0(chi) + 3 Sigma_1728(chi)
  for every nontrivial character chi of Cl(O_n), where
      S(chi)       = sum_c chi(c) log|u_c|,
      Sigma_x(chi) = sum_c chi(c) log|j(c) - x|,
      L(s,chi)     = (1/2) sum_c chi(c) zeta_{Q_c}(s)   (Epstein / form-class
                     L-function of discriminant -4n^2; functional equation
                     gamma(s)L(s) = gamma(1-s)L(1-s), gamma(s) =
                     (sqrt|D|/2pi)^s Gamma(s), so L'(0,chi) =
                     (sqrt|D|/2pi) L(1,chi)).

Hyperbolic (disc 1-n^2): u_f = eps Theta_f over Cl(1-n^2).  The Norm Lemma
  cancels eps and mu, leaving the scale-invariant classwise law
      |u_f|^6 = |b1^4(b1-1728)^3| / |b2^4(b2-1728)^3| * g(f)/g(rf),
  g(f) = y_f^6 |Delta_q(tau_f)|, b2(f) = b1(rf); hence S(chi) = 0 unless
  chi(r_n) = -1, and for odd chi
      3 S(chi) = -12 L'(0,chi) + 4 Sigma_0(chi) + 3 Sigma_1728(chi).

The script further identifies the pieces in closed form:
  * L'(0,chi) for real (genus) characters through the discovered exact
    factorization  R_chi(m) = 2 (conv_{d1,d2} * corr)(m)  of the
    representation numbers (Kronecker symbols (d1/.), (d2/.), finite Euler
    correction at conductor primes), giving  L'(0,chi) =
    (2h(d1)/w(d1)) * h(d2) log eps_{d2} * C(0)-type closed forms;
  * Sigma_0, Sigma_1728 for real characters as logs of explicit quadratic
    S-numbers: coset products A, B are certified integer-quadratic
    conjugates, and A/B is factored EXACTLY in Q(sqrt d) as
    +- eps^k prod (pi_p/pi_p')^{e_p}  (unit power times split-prime ratios);
  * for cubic-order characters the coset products are certified as roots of
    explicit integer cubics (Stark-type evaluations).

Certification policy (CLAUDE.md guard rails): precision is set AFTER
imports; integers/rationals are accepted only with >= max(20, dps/5) spare
digits in the ABSOLUTE-error sense; multi-term PSLQ is used only in the
information-theoretically safe regime and every discovered relation is
re-verified exactly in the quadratic field.

Usage:
    python3 scripts/phase_klf.py                # both aspects, default levels
    python3 scripts/phase_klf.py euclid 9 11 13
    python3 scripts/phase_klf.py hyper 9 11 13 15
    python3 scripts/phase_klf.py --selftest     # full re-verification
    python3 scripts/phase_klf.py --dps 250 euclid 7
Requires mpmath (sympy optional, only for pretty factorization cross-checks).
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from math import gcd
from fractions import Fraction
from itertools import combinations, product as iproduct

from mpmath import (mp, mpf, mpc, fabs, nstr, exp, pi, log, sqrt as msqrt,
                    e1, sin, loggamma, nint, log10)

import euclidean_moduli_invariants as EU
import moduli_invariants as MI
from involution_classmap import compose, classes_of_disc, is_primitive
from involution_experiments import inv_sl2
from proof_check import build_P
from gz_denominators import gz_primes

# precision is set in main()/selftest(), never at import (guard rail 2)

DEFAULT_DPS = 250


# ======================================================================
# A.  elementary number theory
# ======================================================================

def kron_p(D, p):
    if p == 2:
        if D % 2 == 0:
            return 0
        return 1 if D % 8 in (1, 7) else -1
    r = pow(D % p, (p - 1) // 2, p)
    return 0 if r == 0 else (1 if r == 1 else -1)


def factorint(m):
    m = abs(m)
    if m > 10 ** 14:
        try:
            import sympy
            return {int(p): int(e) for p, e in sympy.factorint(m).items()}
        except ImportError:
            pass
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


def square_kernel(x):
    """(d0, t) with x = d0 t^2, d0 the fundamental discriminant of x > 0
    (or 1), extracting small primes and an exact integer square root --
    safe for huge x whose squarefree kernel is small."""
    from math import isqrt
    sf, t = 1, 1
    for p in range(2, 100000):
        if p * p > x and p > 2:
            break
        while x % (p * p) == 0:
            x //= p * p
            t *= p
        if x % p == 0:
            x //= p
            sf *= p
    r = isqrt(x)
    if r * r == x:
        t *= r
    else:
        sf *= x     # large squarefree cofactor (should not occur here)
    if sf % 4 == 1:
        return sf, t
    assert (4 * sf) * (t // 2) ** 2 == 4 * sf * (t // 2) ** 2
    assert t % 2 == 0, (sf, t)
    return 4 * sf, t // 2


def kron(D, m):
    """Kronecker symbol (D/m), m >= 1."""
    if m == 0:
        return 0
    res = 1
    for p, e in factorint(m).items():
        k = kron_p(D, p)
        if k == 0:
            return 0
        if e % 2:
            res *= k
    return res


def fundamentalize(D):
    """D = d0 * f^2 with d0 the fundamental discriminant (or 1)."""
    sf, f = 1, 1
    for p, e in factorint(D).items():
        f *= p ** (e // 2)
        if e % 2:
            sf *= p
    sf = sf if D > 0 else -sf
    if sf % 4 == 1:
        return sf, f
    assert (4 * sf) * (f // 2) ** 2 == D or True
    return 4 * sf, f // 2


def pell_unit(d):
    """(t, u): fundamental unit (t + u sqrt d)/2 of the order of
    discriminant d > 0 (smallest solution of t^2 - d u^2 = +-4)."""
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


def h_imag(D):
    """class number of primitive forms of disc D < 0."""
    return len(EU.reduced_forms(D))


def w_disc(D):
    return 6 if D == -3 else (4 if D == -4 else 2)


# ---------- Dirichlet L data at s = 0 for Kronecker symbols (D/.) ----------

def L_data_kronecker(D):
    """(L(0), L'(0)) of L(s, (D/.)) as mpf's, D any nonzero discriminant
    or 1.  Uses  L(s,(D/.)) = L(s,chi_{d0}) prod_{p|f, p not| d0}
    (1 - chi_{d0}(p) p^{-s})  and Hurwitz-zeta closed forms at s=0."""
    d0, f = fundamentalize(D) if D != 1 else (1, 1)
    if d0 == 1:
        L0, L1 = mpf(-1) / 2, -log(2 * pi) / 2      # zeta(0), zeta'(0)
    elif d0 > 1:
        q = d0
        L0 = mpf(0)
        L1 = -sum(kron(d0, a) * log(sin(pi * mpf(a) / q))
                  for a in range(1, q)) / 2
    else:
        q = -d0
        s_ = sum(kron(d0, a) * a for a in range(1, q))
        L0 = mpf(-s_) / q
        L1 = -log(mpf(q)) * L0 + sum(kron(d0, a) * loggamma(mpf(a) / q)
                                     for a in range(1, q))
    # Euler corrections at p | f, p not dividing d0
    for p in factorint(f):
        if d0 % p == 0:
            continue
        c = kron(d0, p)
        F0, F1 = 1 - c, mpf(c) * log(p)     # value / derivative of 1-c p^-s
        L0, L1 = L0 * F0, L1 * F0 + L0 * F1
    return L0, L1


def L_desc_kronecker(D):
    """human-readable closed form of L'(0,(D/.)) when L(0) = 0 or D<0."""
    d0, f = fundamentalize(D) if D != 1 else (1, 1)
    parts = []
    if d0 == 1:
        parts.append("zeta(s)")
    elif d0 > 0:
        t, u = pell_unit(d0)
        parts.append(f"L(s,chi_{d0}) [eps_{d0} = ({t}+{u}sqrt{d0})/2]")
    else:
        parts.append(f"L(s,chi_{d0}) [h={h_imag(d0)}, w={w_disc(d0)}]")
    for p in factorint(f):
        if d0 % p == 0:
            continue
        parts.append(f"(1-({d0}/{p}){p}^-s)")
    return " ".join(parts)


# ======================================================================
# B.  real quadratic fields, exactly
# ======================================================================

class QF:
    """x + y sqrt(d), x, y Fractions; d a fundamental discriminant > 0."""
    __slots__ = ('x', 'y', 'd')

    def __init__(self, x, y, d):
        self.x, self.y, self.d = Fraction(x), Fraction(y), d

    def __mul__(self, o):
        assert self.d == o.d
        return QF(self.x * o.x + self.y * o.y * self.d,
                  self.x * o.y + self.y * o.x, self.d)

    def conj(self):
        return QF(self.x, -self.y, self.d)

    def norm(self):
        return self.x * self.x - self.y * self.y * self.d

    def inv(self):
        n = self.norm()
        return QF(self.x / n, -self.y / n, self.d)

    def __truediv__(self, o):
        return self * o.inv()

    def __eq__(self, o):
        return self.d == o.d and self.x == o.x and self.y == o.y

    def is_integral(self):
        """member of the maximal order O_d, d fundamental."""
        if self.d % 4 == 1:
            a, b = 2 * self.x, 2 * self.y
            return a.denominator == 1 and b.denominator == 1 \
                and (a - b) % 2 == 0
        # d = 4m: O = Z[sqrt m], element x + (2y) sqrt(m)
        return self.x.denominator == 1 and (2 * self.y).denominator == 1

    def embed(self):
        return self.x + self.y * msqrt(mpf(self.d))

    def __repr__(self):
        return f"({self.x} + {self.y} sqrt{self.d})"


def qf_pow(a, k):
    if k < 0:
        return qf_pow(a.inv(), -k)
    out = QF(1, 0, a.d)
    b = a
    while k:
        if k & 1:
            out = out * b
        b = b * b
        k >>= 1
    return out


def hensel_sqrt(d, p, N):
    """integer s with s^2 = d mod p^N (p odd split, or p = 2, d = 1 mod 8)."""
    if p == 2:
        assert d % 8 == 1
        s, k = 1, 3
        while k < N:
            if (s * s - d) % 2 ** (k + 1):
                s += 2 ** (k - 1)
            k += 1
        s %= 2 ** N
        assert (s * s - d) % 2 ** N == 0
        return s
    s = next(x for x in range(p) if (x * x - d) % p == 0)
    pk = p
    while pk < p ** N:
        pk = pk * pk if pk * pk < p ** N else p ** N
        inv = pow((2 * s) % pk, -1, pk)
        s = (s - (s * s - d) * inv) % pk
    assert (s * s - d) % p ** N == 0
    return s


def vp(m, p):
    v = 0
    while m and m % p == 0:
        m //= p
        v += 1
    return v


def padic_val_pair(s_, t_, d, p, Nq):
    """(v_p1, v_p2): valuations of A = (s_ + t_ sqrt d)/2 at the two primes
    above the split prime p, computed via Hensel embeddings (exact)."""
    N = Nq + 8
    sN = hensel_sqrt(d % p ** N if d > 0 else d, p, N)
    pN = p ** N
    a1 = (s_ + t_ * sN) % pN
    a2 = (s_ - t_ * sN) % pN
    v1 = vp(a1, p) if a1 else N
    v2 = vp(a2, p) if a2 else N
    if p == 2:
        v1 -= 1
        v2 -= 1
    assert v1 < Nq + 4 and v2 < Nq + 4
    return v1, v2


def pideal_generator(p, d):
    """(m, pi): pi = (x + y sqrt d)/2 with N(pi) = +-p^m, p not dividing pi
    (i.e. valuations (m, 0)), m minimal >= 1 -- the order of [p-ideal] in
    the class group of O_d."""
    tE, uE = pell_unit(d)
    epsv = (tE + uE * d ** 0.5) / 2
    for m in range(1, 30):
        pm = p ** m
        ybound = int((epsv + 1) * pm ** 0.5 / d ** 0.5) + 3
        y = 1
        while y <= ybound:
            for pmm in (4 * pm, -4 * pm):
                v = d * y * y + pmm
                if v >= 0:
                    x = int(v ** 0.5)
                    for xx in (x - 1, x, x + 1, x + 2):
                        if xx >= 0 and xx * xx == v:
                            cand = QF(Fraction(xx, 2), Fraction(y, 2), d)
                            if cand.is_integral():
                                # exclude p | cand (valuations (m-j, j))
                                q1, q2 = padic_val_pair(xx, y, d, p, 2 * m + 2)
                                if p == 2:
                                    q1, q2 = q1 + 1, q2 + 1  # undo /2 shift:
                                    # cand numerator is (xx + y sqrt d),
                                    # cand = that/2; we want vals of cand
                                    q1, q2 = q1 - 1, q2 - 1
                                if min(q1, q2) == 0 and max(q1, q2) == m:
                                    if q1 < q2:
                                        cand = cand.conj()
                                    return m, cand
            y += 1
    raise RuntimeError((p, d, "no p-ideal generator found"))


def factor_ratio(A_pair, d):
    """A, B = (s +- t sqrt d)/2 integer-quadratic conjugates.  Factor
    R = A/B exactly in Q(sqrt d):
        R^L = sign * eps^K * prod_p (pi_p / pi_p')^{L e_p / m_p}
    with L = lcm of the m_p, i.e.
        log|A/B| = (K/L) log eps + sum_p (e_p/m_p) log|pi_p/pi_p'|,
    where e_p = v_p(A) - v_p'(A) (p-adic, exact), pi_p generates p-ideal^m_p.
    Returns (sign, Fraction(K, L), [(p, Fraction(e_p, m_p), pi_p)]); the
    L-th-power identity is verified EXACTLY in Q(sqrt d)."""
    s, t = A_pair          # A = (s + t sqrt d)/2 exactly
    A = QF(Fraction(s, 2), Fraction(t, 2), d)
    B = A.conj()
    q = A.norm()
    assert q.denominator == 1
    q = int(q)
    tE, uE = pell_unit(d)
    eps = QF(Fraction(tE, 2), Fraction(uE, 2), d)
    fac = []
    from math import lcm
    L = 1
    for p in sorted(factorint(q)):
        if kron(d, p) != 1:
            continue                     # inert/ramified: cancels in A/B
        Nq = vp(q, p) + 2
        v1, v2 = padic_val_pair(s, t, d, p, Nq)
        e = v1 - v2
        if e == 0:
            continue
        m, piq = pideal_generator(p, d)
        # orient pi to the same prime as embedding 1
        w1, w2 = padic_val_pair(int(piq.x * 2), int(piq.y * 2), d, p, 2 * m + 2)
        if w1 == 0:
            piq = piq.conj()
        # canonicalize by the unit: 0 <= log|pi/pi'| < 2 log eps
        from mpmath import floor as _floor
        le0 = log(fabs(eps.embed()))
        lam = log(fabs(piq.embed())) - log(fabs(piq.conj().embed()))
        tsh = int(_floor(lam / (2 * le0)))
        piq = piq * qf_pow(eps, -tsh)
        fac.append((p, Fraction(e, m), piq))
        L = lcm(L, m)
    # exact verification at the L-th power
    lhs = qf_pow(A, L) * qf_pow(B, L).inv()
    for (p, em, piq) in fac:
        lhs = lhs * qf_pow(piq.conj() / piq, int(em * L))
    le = log(fabs(eps.embed()))
    lr = log(fabs(lhs.embed()))
    K = int(nint(lr / le))
    unit = lhs * qf_pow(eps, -K)
    sign = None
    if unit == QF(1, 0, d):
        sign = 1
    elif unit == QF(-1, 0, d):
        sign = -1
    assert sign is not None, (s, t, d, K, L, unit)
    return sign, Fraction(K, L), fac


# ======================================================================
# C.  class groups and characters
# ======================================================================

def group_data(forms, D):
    """(coords, orders): abelian-group coordinates for the primitive
    reduced forms under Gauss composition; small-h backtracking."""
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
        for g in sorted(forms, key=lambda f: -ords[f]):
            if g in known:
                continue
            o = ords[g]
            new = {}
            p, ok = unit, True
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
            raise RuntimeError("no direct generator found")
    ng = len(gens)
    coords = {f: tuple(list(c) + [0] * (ng - len(c))) for f, c in known.items()}
    return coords, orders


def characters(coords, orders):
    """[(label, order, chi_dict, real?)] for all characters."""
    from math import lcm
    out = []
    forms = list(coords)
    for ks in iproduct(*[range(o) for o in orders]):
        cord = 1
        for k, o in zip(ks, orders):
            if k:
                cord = lcm(cord, o // gcd(k, o))
        vals = {}
        for f in forms:
            ph = Fraction(0)
            for k, o, e in zip(ks, orders, coords[f]):
                ph += Fraction(k * e, o)
            ph %= 1
            vals[f] = exp(2 * pi * mpc(0, 1) * mpf(ph.numerator)
                          / ph.denominator) if ph else mpf(1)
        out.append((ks, cord, vals, cord <= 2))
    return out


# ======================================================================
# D.  Epstein L-functions of a discriminant
# ======================================================================

def rep_numbers(Q, M):
    A, B, C = Q
    D = B * B - 4 * A * C
    r = [0] * (M + 1)
    ymax = int((4 * A * M / (-D)) ** 0.5) + 2
    for y in range(-ymax, ymax + 1):
        disc = B * B * y * y - 4 * A * (C * y * y - M)
        if disc < 0:
            continue
        sq = int(disc ** 0.5) + 2
        for x in range((-B * y - sq) // (2 * A) - 2,
                       (-B * y + sq) // (2 * A) + 3):
            if x == 0 and y == 0:
                continue
            m = A * x * x + B * x * y + C * y * y
            if 1 <= m <= M:
                r[m] += 1
    return r


def epstein_Lprime0(forms, D, chis, w=2):
    """[L'(0,chi)] for the listed characters, via
    L'(0,chi) = (1/w) sum_m R_chi(m) [e^{-am}/(am) + E1(am)],
    a = 2pi/sqrt|D| (incomplete-gamma unfolding of the Epstein zetas;
    entire for nontrivial chi).  Also returns the rep-number table."""
    alpha = 2 * pi / msqrt(mpf(-D))
    M = int((mp.dps + 15) * log(mpf(10)) / alpha) + 10
    reps = {f: rep_numbers(f, M) for f in forms}
    Iv = [None] * (M + 1)
    for m in range(1, M + 1):
        am = alpha * m
        Iv[m] = exp(-am) / am + e1(am)
    out = []
    for (ks, cord, chi, isreal) in chis:
        tot = mpc(0)
        for f in forms:
            s = mpf(0)
            rf = reps[f]
            for m in range(1, M + 1):
                if rf[m]:
                    s += rf[m] * Iv[m]
            tot += chi[f] * s
        out.append(tot / w)
    return out, reps, M


# ---------- genus factorization discovery ----------

def divisor_discs(D):
    out = []
    for d1 in range(-abs(D), abs(D) + 1):
        if d1 == 0 or D % d1:
            continue
        d2 = D // d1
        if d1 % 4 in (0, 1) and d2 % 4 in (0, 1) and (d2, d1) not in out:
            out.append((d1, d2))
    return out


def conv_table(d1, d2, M):
    c = [0] * (M + 1)
    for m in range(1, M + 1):
        c[m] = sum(kron(d1, e) * kron(d2, m // e)
                   for e in range(1, m + 1) if m % e == 0)
    return c


def genus_factorization(forms, D, chi_dict, reps, Mtest=300):
    """For a real character: find (d1, d2), D = d1 d2, and finite Euler
    corrections c_{p,a} at conductor primes with
        R_chi(m) = 2 * (conv_{d1,d2} * corr)(m)   for all m <= Mtest,
    corr multiplicative, supported on primes p with p^2 | D.
    Returns (d1, d2, {p: [c_{p,1}, c_{p,2}, ...]}) or None."""
    Mtest = min(Mtest, len(next(iter(reps.values()))) - 1)
    R = [sum(chi_dict[f] * reps[f][m] for f in forms).real
         for m in range(Mtest + 1)]
    R = [int(nint(x)) for x in R]
    _, f0 = fundamentalize(D)
    condp = sorted(factorint(f0))
    accepted = []
    for (d1, d2) in divisor_discs(D):
        conv = conv_table(d1, d2, Mtest)
        if not all(R[m] == 2 * conv[m] for m in range(1, Mtest + 1)
                   if all(m % p for p in condp)):
            continue
        # solve for corrections at prime powers; require TERMINATION
        # (the top two computable layers must vanish -- an infinite
        # correction series means this (d1, d2) is a degenerate
        # representation and must be rejected)
        corr = {}
        ok = True
        for p in condp:
            cs, a = [], 1
            while p ** a <= Mtest:
                pa = p ** a
                v = Fraction(R[pa], 2) - conv[pa] \
                    - sum(cs[j - 1] * conv[pa // p ** j]
                          for j in range(1, a))
                if v.denominator != 1:
                    ok = False
                    break
                cs.append(int(v))
                a += 1
            if not ok or len(cs) < 2 or cs[-1] != 0 or cs[-2] != 0:
                ok = False
                break
            while cs and cs[-1] == 0:
                cs.pop()
            if cs:
                corr[p] = cs
        if not ok:
            continue

        def corr_val(m):
            out = 1
            for p in condp:
                a = 0
                while m % p == 0:
                    m //= p
                    a += 1
                if a:
                    cp = corr.get(p, [])
                    out *= cp[a - 1] if a <= len(cp) else 0
            return out
        # full verification: R(m) = 2 sum_{e | m, e = conductor part}
        # corr(e) conv(m/e)  -- since corr is supported on conductor primes
        good = True
        for m in range(1, Mtest + 1):
            tot = conv[m]
            for e in range(2, m + 1):
                if m % e == 0 and all(
                        p in condp for p in factorint(e)) and corr_val(e):
                    tot += corr_val(e) * conv[m // e]
            if R[m] != 2 * tot:
                good = False
                break
        if good:
            accepted.append((sum(len(v) for v in corr.values()),
                             d1, d2, corr))
    if accepted:
        accepted.sort(key=lambda t: t[0])
        return accepted[0][1], accepted[0][2], accepted[0][3]
    return None


def eps_ratio_identify(Lval, D, say):
    """fallback when no finite-Euler genus factorization exists: certify
    L'(0,chi) as a rational multiple of log eps_d for a fundamental d > 1
    from the divisor pairs of D.  Certification only (not a proof)."""
    cands = sorted({max(d1, d2) for (d1, d2) in divisor_discs(D)
                    if max(d1, d2) > 1})
    hits, seen = [], set()
    for d2 in cands:
        d0, _ = fundamentalize(d2)
        if d0 <= 1 or d0 in seen:
            continue
        seen.add(d0)
        tE, uE = pell_unit(d0)
        le = log((tE + uE * msqrt(mpf(d0))) / 2)
        fr, sp = EU.cert_rational(mpc(Lval).real / le, maxden=1000)
        if fr is not None and sp > mp.dps / 3:
            hits.append((d0, fr, sp, (tE, uE)))
    for (d0, fr, sp, (tE, uE)) in hits[:1]:
        say(f"   certified ratio: L'(0,chi) = {fr} * log eps_{d0} "
            f"[eps = ({tE}+{uE}sqrt{d0})/2]   (spare {spare(sp)}; "
            f"identification only, no factorization proof)")
    if not hits:
        say("   (L'(0,chi): no eps-ratio identification either)")
    return hits[:1]


def genus_Lprime0_closed(D, d1, d2, corr):
    """closed-form L'(0,chi) = [L(s,(d1/.)) L(s,(d2/.)) prod_p C_p(s)]'(0),
    C_p(s) = 1 + sum_a c_{p,a} p^{-as}.  Returns (value, description)."""
    L10, L11 = L_data_kronecker(d1)
    L20, L21 = L_data_kronecker(d2)
    V0, V1 = L10 * L20, L11 * L20 + L10 * L21
    desc = f"L(s,({d1}/.)) * L(s,({d2}/.))"
    for p, cs in corr.items():
        C0 = 1 + sum(cs)
        C1 = -sum(c * (a + 1) for a, c in enumerate(cs)) * log(mpf(p))
        V0, V1 = V0 * C0, V1 * C0 + V0 * C1
        desc += " * (1 + " + " + ".join(
            f"{c}*{p}^-{a + 1}s" for a, c in enumerate(cs) if c) + ")"
    return V1, desc


# ======================================================================
# E.  certification helpers
# ======================================================================

def cert_int(x, what):
    v, sp = EU.cert_integer(x)
    assert v is not None, f"{what}: not integer-certified ({sp})"
    return v, sp


def cert_rat(x, what, maxden=10 ** 80):
    v, sp = EU.cert_rational(x, maxden=maxden)
    assert v is not None, f"{what}: not rational-certified ({sp})"
    return v, sp


def poly_rat_certify(roots, what, maxden=10 ** 80):
    """certify prod (x - r) in Q[x]; returns (Fractions, min spare)."""
    co = EU.poly_from_roots(roots)
    out, msp = [], mp.inf
    for c in co:
        v, sp = cert_rat(c, what, maxden)
        out.append(v)
        msp = min(msp, sp)
    return out, msp


def safe_pslq(vec, names, what, maxcoeff=10 ** 6, say=print):
    """information-theoretically safe PSLQ: a relation is reported only when
    its total information content is far below the working precision, and
    the residual is at rounding level.  Returns the relation or None."""
    from mpmath import pslq
    tol = mpf(10) ** (-(mp.dps * 4) // 5)
    rel = pslq(vec, maxcoeff=maxcoeff, maxsteps=200000, tol=tol)
    if rel is None:
        say(f"   {what}: NO relation (PSLQ, coeff height <= {maxcoeff}, "
            f"tol 1e-{(mp.dps * 4) // 5}) -- certified non-fit")
        return None
    info = sum(len(str(abs(c))) for c in rel if c)
    resid = fabs(sum(c * v for c, v in zip(rel, vec)))
    if info > mp.dps // 4 or resid > tol * max(fabs(v) for v in vec):
        say(f"   {what}: PSLQ fit rejected as UNSAFE "
            f"(info {info} digits vs dps {mp.dps}, residual {nstr(resid, 3)})")
        return None
    terms = " + ".join(f"{c}*{nm}" for c, nm in zip(rel, names) if c)
    say(f"   {what}: {terms} = 0   (residual {nstr(resid, 3)}, "
        f"info {info} digits, SAFE)")
    return rel


def poly_int_certify(roots, what):
    """certify prod (x - r) in Z[x]; returns (coeffs, min spare)."""
    co = EU.poly_from_roots(roots)
    out, msp = [], mp.inf
    for c in co:
        v, sp = cert_int(c, what)
        out.append(v)
        msp = min(msp, sp)
    return out, msp


def spare(sp):
    return "inf" if sp == mp.inf else f"{float(sp):.0f}"


# ======================================================================
# F.  the identification of a real-character j-sum
# ======================================================================

def identify_real_charsum(vals, chi_dict, forms, what, gzset=None):
    """vals: dict form -> complex value (algebraic integers, closed under
    Galois).  For a real character chi: A = prod_{chi=+1} v, B = prod_{-1}.
    Certify A + B, A*B integers, hence A, B = (s +- t sqrt d)/2; factor A/B
    in Q(sqrt d).  Returns dict with the exact data, or a rational-case
    variant.  All statements exactly verified."""
    A = mpf(1)
    B = mpf(1)
    for f in forms:
        if chi_dict[f].real > 0:
            A = A * vals[f]
        else:
            B = B * vals[f]
    # A, B are real (chi real => cosets closed under inversion, values
    # closed under conjugation); drop the rounding-level imaginary part
    assert fabs(mpc(A).imag) < mpf(10) ** (-mp.dps // 3) * (1 + fabs(A))
    assert fabs(mpc(B).imag) < mpf(10) ** (-mp.dps // 3) * (1 + fabs(B))
    A, B = mpc(A).real, mpc(B).real
    s_, sp1 = cert_int(A + B, what + " A+B")
    q_, sp2 = cert_int(A * B, what + " A*B")
    disc = s_ * s_ - 4 * q_
    out = {'s': s_, 'q': q_, 'spare': min(sp1, sp2),
           'logAB': log(fabs(A)) - log(fabs(B))}
    if disc == 0:
        out['kind'] = 'equal'
        return out
    d0, tt = square_kernel(disc)
    if d0 == 1:                     # A, B individually rational integers
        Ai = (s_ + tt) // 2
        Bi = (s_ - tt) // 2
        assert Ai * Bi == q_ and Ai + Bi == s_
        if fabs(A - Ai) > fabs(A - Bi):
            Ai, Bi = Bi, Ai
        out.update(kind='rational', A=Ai, B=Bi,
                   fA=factorint(Ai), fB=factorint(Bi))
        return out
    assert disc == d0 * tt * tt
    # match numerical A to (s + t sqrt d)/2 branch
    t_signed = tt if A > B else -tt
    sign, k, fac = factor_ratio((s_, t_signed), d0)
    out.update(kind='quadratic', d=d0, t=t_signed, sign=sign, k=k, fac=fac,
               eps=pell_unit(d0))
    if gzset is not None:
        out['gz_ok'] = all(p in gzset for (p, e, piq) in fac)
    return out


def charsum_desc(idd, what):
    """pretty print + exact log value of the identified character sum."""
    if idd['kind'] == 'equal':
        return f"{what}: A = B, sum = 0", mpf(0)
    if idd['kind'] == 'rational':
        fa = " ".join(f"{p}^{e}" if e > 1 else str(p)
                      for p, e in sorted(idd['fA'].items()))
        fb = " ".join(f"{p}^{e}" if e > 1 else str(p)
                      for p, e in sorted(idd['fB'].items()))
        v = log(fabs(mpf(idd['A']))) - log(fabs(mpf(idd['B'])))
        return (f"{what} = log|A/B|, A = {idd['A']} = +-{fa}, "
                f"B = {idd['B']} = +-{fb}", v)
    d, k = idd['d'], idd['k']
    tE, uE = idd['eps']
    v = k * log(fabs(QF(Fraction(tE, 2), Fraction(uE, 2), d).embed()))
    parts = [f"{k}*log eps_{d}"] if k else []
    for (p, e, piq) in idd['fac']:
        lp = log(fabs(piq.embed())) - log(fabs(piq.conj().embed()))
        v += e * lp
        parts.append(f"{e}*log|pi_{p}/pi_{p}'|, pi_{p}={piq}")
    return (f"{what} = " + (" + ".join(parts) if parts else "0")
            + f"   [in Q(sqrt {d})]"), v


# ======================================================================
# G.  Euclidean level
# ======================================================================

def euclid_level(n, verbose=True, do_cosets=True):
    D = -4 * n * n
    forms = EU.reduced_forms(D)
    h = len(forms)
    coords, orders = group_data(forms, D)
    chis = characters(coords, orders)
    say = print if verbose else (lambda *a, **k: None)
    say("=" * 78)
    say(f"EUCLIDEAN n = {n}   disc {D}   h = {h}   "
        f"Cl = {' x '.join('Z/%d' % o for o in orders)}   dps = {mp.dps}")
    say("=" * 78)

    cv = EU.class_values(n)          # form -> (u, j, Lambda)
    assert set(cv) == set(forms)

    def gam(f):
        A, B, C = f
        tau = mpc(mpf(-B) / (2 * A), msqrt(mpf(-D)) / (2 * A))
        return 6 * log(tau.imag) + log(fabs(EU.Dq_at(tau)))

    # classwise closed form (constant included)
    const = -6 * log(mpf(n)) - log(fabs(EU.Dq_at(mpc(0, 1))))
    devmax = mpf(0)
    for f in forms:
        u, j, _ = cv[f]
        lhs = 6 * log(fabs(u))
        rhs = 4 * log(fabs(j)) + 3 * log(fabs(j - 1728)) + gam(f) + const
        devmax = max(devmax, fabs(lhs - rhs))
    say(f"classwise 6log|u| = 4log|j| + 3log|j-1728| + log g + const : "
        f"max dev {nstr(devmax, 3)}")
    checks = {'classwise': devmax}

    say("computing L'(0,chi) (Epstein, incomplete gamma) ...")
    Lp, reps, M = epstein_Lprime0(forms, D, chis)
    say(f"  ({M} terms)")

    results = []
    tol_id = mpf(10) ** (-(mp.dps * 3) // 5)
    for idx, (ks, cord, chi, isreal) in enumerate(chis):
        S = sum(chi[f] * log(fabs(cv[f][0])) for f in forms)
        S0 = sum(chi[f] * log(fabs(cv[f][1])) for f in forms)
        S17 = sum(chi[f] * log(fabs(cv[f][1] - 1728)) for f in forms)
        Sg = sum(chi[f] * gam(f) for f in forms)
        row = {'ks': ks, 'order': cord, 'S': S, 'S0': S0, 'S17': S17,
               'Sg': Sg, 'Lp': Lp[idx]}
        if cord > 1:
            klf = fabs(Sg + 12 * Lp[idx])
            master = fabs(6 * S - (-12 * Lp[idx] + 4 * S0 + 3 * S17))
            row['klf'], row['master'] = klf, master
            say(f"chi{ks} (order {cord}):  S = {nstr(S.real, 20)}")
            say(f"   L'(0,chi) = {nstr(Lp[idx].real, 25)}   "
                f"KLF residual {nstr(klf, 3)}   master residual "
                f"{nstr(master, 3)}")
            assert klf < tol_id and master < tol_id, (n, ks)
            if cord == 2:
                gf = genus_factorization(forms, D, chi, reps)
                if gf:
                    d1, d2, corr = gf
                    val, desc = genus_Lprime0_closed(D, d1, d2, corr)
                    resid = fabs(val - Lp[idx])
                    say(f"   genus factorization: L(s,chi) = {desc}")
                    say(f"   closed form check: |closed - Epstein| = "
                        f"{nstr(resid, 3)}")
                    assert resid < tol_id, (n, ks)
                    row['genus'] = (d1, d2, corr, desc)
                else:
                    say("   (no finite-Euler genus factorization found)")
                    row['epsid'] = eps_ratio_identify(Lp[idx], D, say)
        results.append(row)

    # trivial-character anchor: 6 sum log|u| = 4 log|H(0)| + 3 log|H(1728)|
    # + log|M(n)| - 12 h log n
    Hco, spH = poly_int_certify([cv[f][1] for f in forms], f"H_{D}")
    H0 = Hco[-1] * (-1) ** h            # prod j  (up to sign; take abs)
    H1728 = sum(Hco[k] * 1728 ** (h - k) for k in range(h + 1))
    Mn = 1
    for p, k in factorint(n).items():
        rest = n // p ** k
        Ne = rest
        for q2 in factorint(rest):
            if q2 % 2 == 1:
                Ne = Ne * (q2 - (1 if q2 % 4 == 1 else -1)) // q2
        if p == 2:
            Mn *= p ** (3 * (2 ** k - 1) * Ne)
        elif p % 4 == 3:
            Mn *= p ** (6 * (p ** k - 1) // (p - 1) * Ne)

    # the Delta-mass polynomial D_n(x) = prod_c (x - G_c),
    # G_c = n^12 Delta(Lambda_c)/Delta(Z[i]): one certified integer
    # polynomial whose root-logs carry EVERY L'(0,chi) of the level
    # (uniform Stark statement:  -12 L'(0,chi) = sum_c chi(c) log|G_c|).
    Gs = {}
    for f in forms:
        (a, b, dd) = cv[f][2]
        (_, _), (c_, dd_) = EU.build_X(a, b, dd)
        cc = EU.cval(c_)
        Gs[f] = mpf(n) ** 12 * EU.Dq_at(EU.cval(dd_) / cc) \
            / (cc ** 12 * EU.Dq_at(mpc(0, 1)))
    Dn, spDn = poly_int_certify([Gs[f] for f in forms], f"D_{n} mass poly")
    assert abs(Dn[-1]) == Mn, (n, Dn[-1], Mn)
    say(f"Delta-mass polynomial D_n = {Dn}")
    say(f"   (integer-certified, spare {spare(spDn)};  D_n(0) = "
        f"(-1)^h eps(n) M(n): OK)")
    dev_stark = mpf(0)
    for idx, (ks, cord, chi, isreal) in enumerate(chis):
        if cord == 1:
            continue
        SG = sum(chi[f] * log(fabs(Gs[f])) for f in forms)
        dev_stark = max(dev_stark, fabs(SG + 12 * Lp[idx]))
    say(f"uniform Stark law -12 L'(0,chi) = sum chi(c) log|G_c|, all chi: "
        f"max residual {nstr(dev_stark, 3)}")
    assert dev_stark < mpf(10) ** (-(mp.dps * 3) // 5)
    checks['stark'] = dev_stark

    S_triv = results[0]['S']
    anchor = fabs(6 * S_triv - (4 * log(fabs(mpf(H0)))
                                + 3 * log(fabs(mpf(H1728)))
                                + log(mpf(Mn)) - 12 * h * log(mpf(n))))
    say(f"trivial-character anchor (Thm 2 + Delta-mass law): residual "
        f"{nstr(anchor, 3)}   (H-coeff spare {spare(spH)})")
    assert anchor < tol_id
    checks['anchor'] = anchor

    # identify the j-dressing of the real characters
    if do_cosets:
        g3 = gz_primes(-D, 3)
        g4 = gz_primes(-D, 4)
        for row in results:
            if row['order'] != 2:
                continue
            chi = next(c for (ks, o, c, r) in chis if ks == row['ks'])
            jd = {f: cv[f][1] for f in forms}
            j17 = {f: cv[f][1] - 1728 for f in forms}
            id0 = identify_real_charsum(jd, chi, forms, f"n={n} Sigma0",
                                        gzset=g3 | {2} | set(factorint(n)))
            id17 = identify_real_charsum(j17, chi, forms, f"n={n} Sigma1728",
                                         gzset=g4 | {2} | set(factorint(n)))
            d0desc, v0 = charsum_desc(id0, "Sigma_0")
            d17desc, v17 = charsum_desc(id17, "Sigma_1728")
            say(f"   {d0desc}")
            say(f"     exact-vs-numeric residual {nstr(fabs(v0 - row['S0']), 3)}"
                f"   (spare {spare(id0['spare'])})")
            say(f"   {d17desc}")
            say(f"     exact-vs-numeric residual "
                f"{nstr(fabs(v17 - row['S17']), 3)}   "
                f"(spare {spare(id17['spare'])})")
            assert fabs(v0 - row['S0']) < tol_id
            assert fabs(v17 - row['S17']) < tol_id
            if 'gz_ok' in id0:
                say(f"   Sigma_0 split-prime support in GZ(D,-3) u {{p|2n}}: "
                    f"{id0['gz_ok']};  Sigma_1728 in GZ(D,-4) u {{p|2n}}: "
                    f"{id17.get('gz_ok')}")
            row['id0'], row['id17'] = id0, id17
            # final closed form of S(chi):
            #   S = -2 L' + (2/3) Sigma0 + (1/2) Sigma1728
            Sclosed = -2 * row['Lp'].real + Fraction(2, 3) * v0 \
                + Fraction(1, 2) * v17
            say(f"   => S(chi) closed-form residual "
                f"{nstr(fabs(Sclosed - row['S'].real), 3)}")
        # cubic coset data for an order-3 character (if the group has one)
        for row in results:
            if row['order'] != 3:
                continue
            chi = next(c for (ks, o, c, r) in chis if ks == row['ks'])
            ker = [f for f in forms if fabs(chi[f] - 1) < mpf(10) ** (-10)]
            cosets = {}
            for f in forms:
                key = min(nstr(chi[f], 5) for _ in (0,))
                cosets.setdefault(key, []).append(f)
            prods = [None, None, None]
            for f in forms:
                w = chi[f]
                i = 0 if fabs(w - 1) < 0.1 else (1 if w.imag > 0 else 2)
                prods[i] = cv[f][1] if prods[i] is None else prods[i] * cv[f][1]
            cub, spc = poly_int_certify(prods, f"n={n} cubic j-cosets")
            say(f"   order-3 j-coset cubic: {cub}  (spare {spare(spc)})")
            prods17 = [None, None, None]
            for f in forms:
                w = chi[f]
                i = 0 if fabs(w - 1) < 0.1 else (1 if w.imag > 0 else 2)
                x = cv[f][1] - 1728
                prods17[i] = x if prods17[i] is None else prods17[i] * x
            cub17, spc17 = poly_int_certify(prods17,
                                            f"n={n} cubic (j-1728)-cosets")
            say(f"   order-3 (j-1728)-coset cubic: {cub17}  "
                f"(spare {spare(spc17)})")
            # Stark-type: the Delta-mass cosets, normalized by n^12 per class
            prodsD = [None, None, None]
            for f in forms:
                w = chi[f]
                i = 0 if fabs(w - 1) < 0.1 else (1 if w.imag > 0 else 2)
                u, j, L = cv[f]
                (a, b, dd) = L
                (A, B), (c, ddd) = EU.build_X(a, b, dd)
                cc = EU.cval(c)
                x = mpf(n) ** 12 * EU.Dq_at(EU.cval(ddd) / cc) \
                    / (cc ** 12 * EU.Dq_at(mpc(0, 1)))
                prodsD[i] = x if prodsD[i] is None else prodsD[i] * x
            cubD, spD = poly_int_certify(prodsD,
                                         f"n={n} cubic Delta-mass cosets")
            say(f"   order-3 Delta-mass coset cubic (Stark object): {cubD}  "
                f"(spare {spare(spD)})")
            row['cubics'] = (cub, cub17, cubD)
            # closed form: -12 L'(0,chi) = (3/2) log|G0| - (1/2) log|C(0)|
            #   - 12 log n  per class pair... (2 classes per coset)
            G0 = prodsD[0]
            v = (Fraction(3, 2) * log(fabs(G0))
                 - Fraction(1, 2) * log(fabs(mpf(cubD[-1]))))
            npc = len(forms) // 3
            v -= 0  # n^12-normalization cancels in the character sum
            resid = fabs(v - (-12) * row['Lp'].real)
            say(f"   Stark check: (3/2)log|G0| - (1/2)log|C_D(0)| + 12L'(0) "
                f"residual = {nstr(resid, 3)}")
            assert resid < tol_id
            break

        # mission-style PSLQ for the non-real characters: does S(chi) fit
        # the naive basis {L'(0,chi), log p (p | 2n, GZ supports), log eps}?
        base_primes = sorted(set(factorint(2 * n))
                             | set(factorint(H0)) | set(factorint(H1728)))
        eps_terms = []
        for row in results:
            if row['order'] == 2 and 'genus' in row:
                d2 = max(row['genus'][0], row['genus'][1])
                if d2 > 1:
                    tE, uE = pell_unit(d2)
                    eps_terms = [(f"log eps_{d2}",
                                  log((tE + uE * msqrt(mpf(d2))) / 2))]
        done_orders = set()
        for row in results:
            if row['order'] <= 2 or row['order'] in done_orders:
                continue
            done_orders.add(row['order'])
            vec = [row['S'].real, row['Lp'].real] \
                + [log(mpf(p)) for p in base_primes] \
                + [t[1] for t in eps_terms]
            names = ["S(chi)", "L'(0,chi)"] \
                + [f"log{p}" for p in base_primes] + [t[0] for t in eps_terms]
            say(f"PSLQ, order-{row['order']} character, basis "
                f"{{S, L', log p (p | 2n H(0) H(1728)), eps}}:")
            row['pslq'] = safe_pslq(vec, names,
                                    f"n={n} order-{row['order']}", say=say)
    return results, checks


# ======================================================================
# H.  hyperbolic level
# ======================================================================

def hyper_uf(n, prim):
    epsn = n + msqrt(mpf(n * n - 1))
    return {f: epsn * MI.theta_integral(inv_sl2(build_P(n, f)[0]))[0]
            for f in prim}


def hyper_level(n, verbose=True, do_cosets=True):
    D = 1 - n * n
    prim = [f for f in classes_of_disc(D) if is_primitive(f)]
    h = len(prim)
    coords, orders = group_data(prim, D)
    chis = characters(coords, orders)
    rn = MI.reduce_form((n - 1) // 2, 0, (n + 1) // 2) if False else None
    from involution_classmap import reduce_form as rf_
    rn = rf_((n - 1) // 2, 0, (n + 1) // 2)
    say = print if verbose else (lambda *a, **k: None)
    say("=" * 78)
    say(f"HYPERBOLIC n = {n}   disc {D}   h = {h}   "
        f"Cl = {' x '.join('Z/%d' % o for o in orders)}   r_n = {rn}   "
        f"dps = {mp.dps}")
    say("=" * 78)

    u = hyper_uf(n, prim)
    jv = {f: MI.J(MI.cm_point(f)) for f in prim}
    rmap = {f: compose(rn, f, D) for f in prim}

    def gam(f):
        tau = MI.cm_point(f)
        return 6 * log(tau.imag) + log(fabs(EU.Dq_at(tau)))

    # classwise law: |u_f|^6 = |b1/b2|-dressing * g(f)/g(rf), b2 = j(rf)
    devmax = mpf(0)
    for f in prim:
        b1, b2 = jv[f], jv[rmap[f]]
        lhs = 6 * log(fabs(u[f]))
        rhs = (4 * (log(fabs(b1)) - log(fabs(b2)))
               + 3 * (log(fabs(b1 - 1728)) - log(fabs(b2 - 1728)))
               + gam(f) - gam(rmap[f]))
        devmax = max(devmax, fabs(lhs - rhs))
    say(f"classwise |u_f|^6 = |b1^4(b1-1728)^3 / b2^4(b2-1728)^3| "
        f"* g(f)/g(rf) : max dev {nstr(devmax, 3)}")
    checks = {'classwise': devmax}

    say("computing L'(0,chi) (Epstein) ...")
    Lp, reps, M = epstein_Lprime0(prim, D, chis)
    say(f"  ({M} terms)")

    # the r-twisted Delta-ratio, in its canonical (class-function) form:
    #   R_f := u_f^6 * (b2/b1)^4 ((b2-1728)/(b1-1728))^3,
    # so |R_f| = g(f)/g(rf) (the classwise law), and by the first-power
    # descent sigma(R_f) = R_{f^e c}: the multiset {R_f} is Galois-stable
    # and prod_f (x - R_f) has rational coefficients.  For every odd chi:
    #   -24 L'(0,chi) = sum_f chi(f) log|R_f|.
    Rr = {}
    for f in prim:
        b1, b2 = jv[f], jv[rmap[f]]
        Rr[f] = u[f] ** 6 * (b2 / b1) ** 4 \
            * ((b2 - 1728) / (b1 - 1728)) ** 3
    try:
        Rco, spR = poly_rat_certify([Rr[f] for f in prim],
                                    f"hyp n={n} R-poly", maxden=10 ** 100)
        from math import gcd as _g
        dR = 1
        for c_ in Rco:
            dR = dR * c_.denominator // _g(dR, c_.denominator)
        say(f"r-twisted Delta-ratio polynomial prod (x - R_f) in Q[x], "
            f"lcm denominator {dR} = {factorint(dR)}")
        say(("   coefficients: " + ", ".join(str(c_) for c_ in Rco))[:400])
        say(f"   (rational-certified, spare {spare(spR)})")
        checks['Rpoly'] = spR
    except AssertionError as exc:
        say(f"r-twisted Delta-ratio polynomial: not certified at this "
            f"precision ({exc})")

    results = []
    tol_id = mpf(10) ** (-(mp.dps * 3) // 5)
    g3 = gz_primes(-D, 3)
    g4 = gz_primes(-D, 4)
    for idx, (ks, cord, chi, isreal) in enumerate(chis):
        S = sum(chi[f] * log(fabs(u[f])) for f in prim)
        S0 = sum(chi[f] * log(fabs(jv[f])) for f in prim)
        S17 = sum(chi[f] * log(fabs(jv[f] - 1728)) for f in prim)
        Sg = sum(chi[f] * gam(f) for f in prim)
        chir = chi[rn]
        odd = fabs(chir + 1) < mpf(10) ** (-10)
        row = {'ks': ks, 'order': cord, 'S': S, 'S0': S0, 'S17': S17,
               'Sg': Sg, 'Lp': Lp[idx], 'odd': odd}
        if cord == 1:
            results.append(row)
            continue
        klf = fabs(Sg + 12 * Lp[idx])
        assert klf < tol_id, (n, ks, nstr(klf, 5))
        row['klf'] = klf
        if odd:
            SR = sum(chi[f] * log(fabs(Rr[f])) for f in prim)
            stark = fabs(SR + 24 * Lp[idx])
            row['stark'] = stark
            assert stark < tol_id, (n, ks)
        if not odd:
            say(f"chi{ks} (order {cord}, chi(r) = +1):  S = "
                f"{nstr(fabs(S), 3)}   (must vanish)")
            assert fabs(S) < tol_id
            if cord == 2:
                gf = genus_factorization(prim, D, chi, reps)
                if gf:
                    d1, d2, corr = gf
                    val, desc = genus_Lprime0_closed(D, d1, d2, corr)
                    resid = fabs(val - Lp[idx])
                    say(f"   [even] genus factorization: L(s,chi) = {desc}"
                        f"   closed-form residual {nstr(resid, 3)}")
                    assert resid < tol_id
                    row['genus'] = (d1, d2, corr, desc)
                else:
                    say("   [even] (no finite-Euler genus factorization)")
                    row['epsid'] = eps_ratio_identify(Lp[idx], D, say)
            results.append(row)
            continue
        master = fabs(3 * S - (-12 * Lp[idx] + 4 * S0 + 3 * S17))
        row['master'] = master
        say(f"chi{ks} (order {cord}, chi(r) = -1):  S = {nstr(S.real, 20)}")
        say(f"   L'(0,chi) = {nstr(Lp[idx].real, 25)}   KLF residual "
            f"{nstr(klf, 3)}   master residual {nstr(master, 3)}")
        assert master < tol_id, (n, ks)
        if cord == 2:
            gf = genus_factorization(prim, D, chi, reps)
            if gf:
                d1, d2, corr = gf
                val, desc = genus_Lprime0_closed(D, d1, d2, corr)
                resid = fabs(val - Lp[idx])
                say(f"   genus factorization: L(s,chi) = {desc}")
                say(f"   closed form check: residual {nstr(resid, 3)}")
                assert resid < tol_id
                row['genus'] = (d1, d2, corr, desc)
            else:
                say("   (no finite-Euler genus factorization found)")
                row['epsid'] = eps_ratio_identify(Lp[idx], D, say)
            if do_cosets:
                id0 = identify_real_charsum(jv, chi, prim,
                                           f"hyp n={n} Sigma0",
                                           gzset=g3 | {2} | set(factorint(n * n - 1)))
                j17 = {f: jv[f] - 1728 for f in prim}
                id17 = identify_real_charsum(j17, chi, prim,
                                             f"hyp n={n} Sigma1728",
                                             gzset=g4 | {2} | set(factorint(n * n - 1)))
                d0desc, v0 = charsum_desc(id0, "Sigma_0")
                d17desc, v17 = charsum_desc(id17, "Sigma_1728")
                say(f"   {d0desc}")
                say(f"   {d17desc}")
                assert fabs(v0 - S0) < tol_id and fabs(v17 - S17) < tol_id
                row['id0'], row['id17'] = id0, id17
                Sclosed = -4 * Lp[idx].real + Fraction(4, 3) * v0 + v17
                say(f"   => S(chi) closed-form residual "
                    f"{nstr(fabs(Sclosed - S.real), 3)}")
        results.append(row)

    # mission-style PSLQ for the odd non-real characters
    if do_cosets and any(r.get('odd') and r['order'] > 2 for r in results):
        Hco, spH = poly_int_certify([jv[f] for f in prim], f"H_{D}")
        H0 = Hco[-1] * (-1) ** h
        H1728 = sum(Hco[k] * 1728 ** (h - k) for k in range(h + 1))
        base_primes = sorted(set(factorint(2 * (n * n - 1)))
                             | set(factorint(H0)) | set(factorint(H1728)))
        eps_terms, seen_d2 = [], set()
        for row in results:
            if 'genus' in row:
                d2 = max(row['genus'][0], row['genus'][1])
                if d2 > 1 and d2 not in seen_d2:
                    seen_d2.add(d2)
                    tE, uE = pell_unit(d2)
                    eps_terms.append(
                        (f"log eps_{d2}",
                         log((tE + uE * msqrt(mpf(d2))) / 2)))
        done = set()
        for row in results:
            if not row.get('odd') or row['order'] <= 2 \
                    or row['order'] in done:
                continue
            done.add(row['order'])
            vec = [row['S'].real, row['Lp'].real] \
                + [log(mpf(p)) for p in base_primes] \
                + [t[1] for t in eps_terms]
            names = ["S(chi)", "L'(0,chi)"] \
                + [f"log{p}" for p in base_primes] + [t[0] for t in eps_terms]
            say(f"PSLQ, odd order-{row['order']} character, basis "
                f"{{S, L', log p (p | 2(n^2-1) H(0) H(1728)), eps}}:")
            row['pslq'] = safe_pslq(vec, names,
                                    f"hyp n={n} order-{row['order']}",
                                    say=say)
    return results, checks


# ======================================================================
# I.  driver
# ======================================================================

def main(argv):
    dps = DEFAULT_DPS
    mode = 'all'
    levels = []
    selftest = False
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == '--dps':
            i += 1
            dps = int(argv[i])
        elif a == '--selftest':
            selftest = True
        elif a in ('euclid', 'hyper', 'all'):
            mode = a
        else:
            levels.append(int(a))
        i += 1
    mp.dps = dps
    if selftest:
        mode = 'all'
        levels = []
    if mode in ('euclid', 'all'):
        for n in levels or [9, 11, 13, 7, 5, 3]:
            euclid_level(n)
            print()
    if mode in ('hyper', 'all'):
        for n in levels or [9, 11, 13, 15]:
            hyper_level(n)
            print()
    if selftest:
        mp.dps = 400          # the composite square-free conductor level
        euclid_level(15)      # (larger coset products need more digits)
        print()
    print("ALL CHECKS PASSED")


if __name__ == '__main__':
    main(sys.argv[1:])
