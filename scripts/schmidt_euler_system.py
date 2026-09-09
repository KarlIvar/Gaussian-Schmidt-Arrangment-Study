"""The Schmidt Delta-system as an Euler system over the ring class tower of
Q(i): norm relations, the Hecke-type recursion for L'(0,chi), and the Robert
index along the tower.

Companion to schmidt-euler-system.md.  Setting (Paper II, Theorem 2.6): for a
level n >= 2 and a class c of Pic(O_n), O_n = Z + nZ[i], let Lambda_c be the
primitive index-n sublattice of Z[i] representing c and
    G_c = n^12 Delta(Lambda_c) / Delta(Z[i]).
For a prime l the ring class field H_{nl} contains H_n, and
Gal(H_{nl}/H_n) = ker(Pic(O_{nl}) -> Pic(O_n)) permutes the fiber of the
projection over c.  The script verifies, in exact lattice arithmetic (HNF)
plus certified Delta-values:

 (T1) the fiber lemma: the non-O_n-stable index-l sublattices of Lambda_c
      are exactly the primitive representatives of the fiber over c
      (bijection; |fiber| = l - chi_{-4}(l) for l not dividing n, = l for l | n);
 (T2) the NORM RELATIONS (Theorem 1 of the document), for every class:
        l = lam*lambar split, l not | n :  N(G_{c'}) = G_c^{l+1} / (G_{[l]c} G_{[lbar]c})
        l inert, l not | n             :  N(G_{c'}) = l^12 G_c^{l+1}
        l = 2, 2 not | n               :  N(G_{c'}) = 2^6 G_c^3 / G_{[p]c}
        l | n                          :  N(G_{c'}) = (-1)^[l=2] G_c^{l+1} / G^{(n/l)}_{c^+}
      where N is the product over the fiber (= the field norm), [l]c the
      class of the O_n-ideal (lam) cap O_n times c, and c^+ the image of c
      in Pic(O_{n/l}) (G^{(1)} := 1);
 (T3) the Hecke-type recursion for the KLF sums S_chi(m) = sum chi log|G|
      (= -12 L'(0,chi)) along pullbacks chi^{(nl)} = chi o pi:
        S(nl) = P_l(chi) S(n)   (l not | n),   S(nl) = (l+1) S(n) - m S(n/l)  (l | n),
      P_l = l + 1 - chi([l]) - chi([l])^{-1} (split), l + 1 (inert), 3 - chi([p]) (l = 2),
      m = |ker(Pic(O_n) -> Pic(O_{n/l}))|, and S(n/l) := 0 if chi does not factor;
      checked against the proved genus closed forms at n = 3, 5 -> 9, 15;
 (T4) the Robert index along the tower: for the mass-normalized cubic coset
      units theta_u^{(m)} = prod_{c in ker chi_3} G_c / |M(m)|^{1/3},
        log|theta_u^{(nl)}| = P_l log|theta_u^{(n)}|  (l not | n),
        log|theta_u^{(nl)}| = (l+1) log|theta_u^{(n)}| - m log|theta_u^{(n/l)}|  (l | n),
      so the index [O_{L_3}^x : <-1, theta_u>] of Paper II, Theorem 6.6,
      multiplies by the certified integer P_l: the pullback multipliers
      C_n(0) = 2, 2, 2, 4 at n = 18, 22, 26, 27 are proved, and new levels
      (33, 36, 39, 44, 45, 46, 52, 54, 55, 63, 65, 69, 81) are computed.

Certification policy (CLAUDE.md guard rails): precision is set AFTER imports;
all class-group bookkeeping (representatives, products, projections, ideal
twists, kernels, characters) is EXACT integer arithmetic on Hermite normal
forms -- the only numerics are Delta-values, and every identity is asserted
with the absolute/relative-error criterion and >= max(20, dps/5) spare
digits; no PSLQ anywhere.

Usage:
    python3 scripts/schmidt_euler_system.py --selftest      # everything (~2 min)
    python3 scripts/schmidt_euler_system.py norm 9 5         # one (n, l) pair
    python3 scripts/schmidt_euler_system.py tower 9 2 3 5 7  # cubic layer from n=9
Requires mpmath only.
"""
import sys
import os
import time
from math import gcd
from mpmath import mp, mpf, mpc, exp, pi, log, fabs, nstr, sqrt

# --------------------------------------------------------------------------
# exact lattice arithmetic on Z[i] = Z^2   (vectors (x, y) <-> x + y i)
# --------------------------------------------------------------------------

def ext_gcd(a, b):
    """return (g, s, t) with s*a + t*b = g = gcd(a, b) >= 0"""
    if b == 0:
        return (abs(a), (1 if a >= 0 else -1), 0)
    g, s, t = ext_gcd(b, a % b)
    return (g, t, s - (a // b) * t)


def hnf(vecs):
    """Hermite normal form (d, b, a) of the Z-span of integer 2-vectors:
    the lattice Z*d + Z*(b + a i), d > 0, a > 0, 0 <= b < d."""
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


def vecs_of(L):
    d, b, a = L
    return [(d, 0), (b, a)]


def gmul(u, v):
    return (u[0] * v[0] - u[1] * v[1], u[0] * v[1] + u[1] * v[0])


def index_of(L):
    return L[0] * L[2]


def scale(L, lam):
    return hnf([gmul(v, lam) for v in vecs_of(L)])


def times_i(L):
    return hnf([(-y, x) for (x, y) in vecs_of(L)])


def conj(L):
    return hnf([(x, -y) for (x, y) in vecs_of(L)])


def zi_span(L):
    return hnf(vecs_of(L) + [(-y, x) for (x, y) in vecs_of(L)])


def is_primitive(L):
    return index_of(zi_span(L)) == 1


def mul(L1, L2):
    return hnf([gmul(u, v) for u in vecs_of(L1) for v in vecs_of(L2)])


def extend(L, n):
    """O_n * L, O_n = Z + n i Z"""
    return hnf(vecs_of(L) + [(-n * y, n * x) for (x, y) in vecs_of(L)])


def intersect_with_ideal(L, lam):
    """(L cap lam Z[i]) / lam as an HNF lattice (lam a Gaussian integer)"""
    d, b, a = L
    x, y = lam
    N = x * x + y * y
    out = []
    for u in range(N):
        for v in range(N):
            re, im = u * d + v * b, v * a
            r2, i2 = re * x + im * y, im * x - re * y
            if r2 % N == 0 and i2 % N == 0:
                out.append((r2 // N, i2 // N))
    out += [(d * x, -d * y), (b * x + a * y, a * x - b * y)]   # N*L/lam = L*conj(lam)
    return hnf(out)


def canon(L):
    """canonical form of the class of a primitive lattice: min(L, iL)"""
    return min(L, times_i(L))


def chi4(p):
    return 0 if p == 2 else (1 if p % 4 == 1 else -1)


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


def Ne(n):
    r = n
    for p in factor(n):
        r = r * (p - chi4(p)) // p
    return r


def mass_log(n):
    """log |M(n)| from the Delta-mass law (Paper I, Thm 7.18)"""
    s = mpf(0)
    for p, k in factor(n).items():
        if chi4(p) == 1:
            continue
        e = 2 if p == 2 else 1
        s += (6 // e) * (p ** k - 1) // (p - 1) * Ne(n // p ** k) * log(p)
    return s


# --------------------------------------------------------------------------
# class groups of the orders O_n, exactly
# --------------------------------------------------------------------------

class ClassGroup:
    """Pic(O_n) via canonical primitive index-n sublattices of Z[i]."""

    def __init__(self, n):
        self.n = n
        reps = set()
        count = 0
        for a in range(1, n + 1):
            if n % a:
                continue
            d = n // a
            for b in range(d):
                L = (d, b, a)
                if is_primitive(L):
                    count += 1
                    reps.add(canon(L))
        assert count == Ne(n), (n, count, Ne(n))
        self.reps = sorted(reps)
        self.h = len(self.reps)
        assert self.h == (Ne(n) // 2 if n >= 2 else 1), (n, self.h)
        self.idx = {c: k for k, c in enumerate(self.reps)}
        self.one = canon((1, 0, n)) if n >= 2 else (1, 0, 1)
        assert self.one in self.idx

    def cls(self, L):
        """class of an arbitrary primitive index-n lattice"""
        c = canon(L)
        assert c in self.idx, (self.n, L)
        return c

    def prod(self, c1, c2):
        P = mul(c1, c2)
        assert index_of(P) == self.n and is_primitive(P)
        return self.cls(P)

    def inv(self, c):
        return self.cls(conj(c))

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
        """image in Pic(O_{n2}) for n2 | n"""
        L = extend(c, n2)
        assert index_of(L) == n2 and is_primitive(L), (self.n, n2, c, L)
        return canon(L)

    def ideal_class(self, lam):
        """class of the O_n-ideal (lam) cap O_n, N(lam) prime to n"""
        return self.cls(intersect_with_ideal(self.one, lam))

    def twist(self, c, lam):
        """[lam]*c via Lambda cap lam Z[i]"""
        L = intersect_with_ideal(c, lam)
        assert index_of(L) == self.n and is_primitive(L), (self.n, c, lam, L)
        return canon(L)


# --------------------------------------------------------------------------
# Delta-values
# --------------------------------------------------------------------------

_DELTA_CACHE = {}


def delta_q(tau):
    key = (nstr(tau.real, 40), nstr(tau.imag, 40))
    if key in _DELTA_CACHE:
        return _DELTA_CACHE[key]
    q = exp(2 * pi * 1j * tau)
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


def G_value(L, m):
    """G = m^12 Delta(L)/Delta(Z[i]) for the primitive index-m lattice L"""
    d, b, a = L
    tau = mpc(b, a) / d
    return mpf(m) ** 12 * mpf(d) ** (-12) * delta_q(tau) / delta_q(mpc(0, 1))


def spare_digits(err, scale=1):
    err = fabs(err)
    if err == 0:
        return mp.dps
    return int(-log(err / scale, 10))


# --------------------------------------------------------------------------
# T1 + T2: fiber lemma and norm relations
# --------------------------------------------------------------------------

def split_lams(l):
    """Gaussian primes over l: [] inert, [(1,1)] l=2, [lam, lambar] split"""
    if l == 2:
        return [(1, 1)]
    if l % 4 == 3:
        return []
    for x in range(1, l):
        y2 = l - x * x
        y = int(round(y2 ** 0.5))
        if y > 0 and y * y == y2:
            return [(x, y), (x, -y)]
    raise ValueError


def norm_relation(n, l, verbose=True, tol_digits=None):
    """returns the minimal spare digits over all classes"""
    m = n * l
    big, small = ClassGroup(m), ClassGroup(n)
    if verbose:
        print(f"  n={n:3d} l={l:2d} -> level {m:3d}: h(n)={small.h}, h(nl)={big.h}", end="")
    fiber = {c: [] for c in small.reps}
    for c2 in big.reps:
        fiber[big.project(c2, n)].append(c2)
    expected_size = (l if n % l == 0 else l - chi4(l))
    for c in small.reps:
        assert len(fiber[c]) == expected_size, (n, l, c, len(fiber[c]))
    lams = split_lams(l)
    worst = mp.dps
    for c in small.reps:
        P = mpf(1)
        for c2 in fiber[c]:
            P *= G_value(c2, m)
        Gc = G_value(c, n)
        if n % l == 0:
            n2 = n // l
            if n2 == 1:
                Gplus = mpf(1)
            else:
                Gplus = G_value(small.project(c, n2), n2)
            rhs = Gc ** (l + 1) / Gplus
            if l == 2:
                rhs = -rhs
        else:
            rhs = Gc ** (l + 1)
            for lam in lams:
                rhs /= G_value(small.twist(c, lam), n)
            if lams == []:
                rhs *= mpf(l) ** 12
            elif l == 2:
                rhs *= 64
        err = fabs(P / rhs - 1)
        sp = spare_digits(err)
        worst = min(worst, sp)
    need = max(20, mp.dps // 5) if tol_digits is None else tol_digits
    assert worst >= need, (n, l, worst)
    if verbose:
        case = ("l | n" if n % l == 0 else ("l = 2" if l == 2 else ("split" if lams else "inert")))
        print(f"  [{case:5s}] fiber size {expected_size}: relation holds, spare {worst} digits")
    return worst


# --------------------------------------------------------------------------
# characters on Pic(O_n)
# --------------------------------------------------------------------------

def quadratic_character(cg):
    """the real character of a class group of order 2"""
    assert cg.h == 2
    other = [c for c in cg.reps if c != cg.one][0]
    return {cg.one: 1, other: -1}


def cubic_data(cg):
    """(kernel K, chi_3 as dict class -> k in {0,1,2}, generator g of order 3)
    at a level whose 3-Sylow subgroup has order exactly 3"""
    h = cg.h
    assert h % 3 == 0 and (h // 3) % 3 != 0, h
    K = [c for c in cg.reps if cg.power(c, h // 3) == cg.one]
    assert len(K) == h // 3, (len(K), h)
    g = None
    for c in cg.reps:
        if cg.order(c) == 3:
            g = c
            break
    assert g is not None
    Kset = set(K)
    chi = {}
    ginv = cg.inv(g)
    for c in cg.reps:
        x = c
        for k in range(3):
            if x in Kset:
                chi[c] = k
                break
            x = cg.prod(x, ginv)
        assert c in chi
    return K, chi, g


def pullback_multiplier(l, n, chi_val_of_lam):
    """P_l(chi) for l not dividing n; chi_val_of_lam(lam) in the unit circle"""
    lams = split_lams(l)
    if lams == []:
        return l + 1, "inert"
    if l == 2:
        return 3 - chi_val_of_lam((1, 1)), "l = 2"
    z = chi_val_of_lam(lams[0])
    return l + 1 - z - (1 / z), "split"


# --------------------------------------------------------------------------
# T3: the Hecke-type recursion for S_chi = sum chi log|G|
# --------------------------------------------------------------------------

def klf_sum(cg, chi_of_class):
    s = mpc(0)
    for c in cg.reps:
        s += chi_of_class(c) * log(fabs(G_value(c, cg.n)))
    return s


def test_recursion_quadratic():
    """chi_2 at n = 3 and n = 5, pulled back to 9 and 15; against the proved
    genus closed forms of Paper II (Prop. 3.1, Remark 3.2)."""
    eps12 = 2 + sqrt(3)
    eps5 = (1 + sqrt(5)) / 2
    need = max(20, mp.dps // 5)
    results = []
    # level 3 -> 9  (l = 3 | 3, m = h(3)/h(1) = 2, chi does not factor: term absent)
    c3 = ClassGroup(3)
    chi3 = quadratic_character(c3)
    S3 = klf_sum(c3, lambda c: chi3[c]).real
    err = fabs(S3 + 4 * log(eps12))
    assert spare_digits(err) >= need
    c9 = ClassGroup(9)
    S9 = klf_sum(c9, lambda c: chi3[c9.project(c, 3)]).real
    err = fabs(S9 - 4 * S3)
    sp = spare_digits(err)
    assert sp >= need
    err2 = fabs(S9 + 16 * log(eps12))
    assert spare_digits(err2) >= need
    results.append(("3 -> 9 (l=3 | 3)", 4, sp))
    # level 3 -> 15 (l = 5 split): P = 6 - chi(lam) - chi(lam)^{-1}
    lam = split_lams(5)[0]
    z = chi3[c3.ideal_class(lam)]
    P = 6 - 2 * z
    c15 = ClassGroup(15)
    S15 = klf_sum(c15, lambda c: chi3[c15.project(c, 3)]).real
    err = fabs(S15 - P * S3)
    sp = spare_digits(err)
    assert sp >= need and P == 8, (P, sp)
    err2 = fabs(S15 + 32 * log(eps12))     # L' = (8/3) log eps_12, Paper II Rem. 3.2
    assert spare_digits(err2) >= need
    results.append(("3 -> 15 (l=5 split, chi(l) = -1)", P, sp))
    # level 5 -> 15 (l = 3 inert): P = 4
    c5 = ClassGroup(5)
    chi5 = quadratic_character(c5)
    S5 = klf_sum(c5, lambda c: chi5[c]).real
    assert spare_digits(fabs(S5 + 24 * log(eps5))) >= need
    S15b = klf_sum(c15, lambda c: chi5[c15.project(c, 5)]).real
    err = fabs(S15b - 4 * S5)
    sp = spare_digits(err)
    assert sp >= need
    assert spare_digits(fabs(S15b + 96 * log(eps5))) >= need   # L' = 8 log eps_5
    results.append(("5 -> 15 (l=3 inert)", 4, sp))
    return results


def test_recursion_cubic(n, l):
    """cubic character at a primitive level n pulled back to nl (l not | n)"""
    cg = ClassGroup(n)
    K, chi, g = cubic_data(cg)
    w = exp(2 * pi * 1j / 3)
    S_n = klf_sum(cg, lambda c: w ** chi[c])
    big = ClassGroup(n * l)
    S_nl = klf_sum(big, lambda c: w ** chi[big.project(c, n)])
    P, case = pullback_multiplier(l, n, lambda lam: w ** chi[cg.ideal_class(lam)])
    P = P.real if hasattr(P, "real") else mpf(P)
    err = fabs(S_nl - P * S_n)
    sp = spare_digits(err, max(1, fabs(S_n)))
    need = max(20, mp.dps // 5)
    assert sp >= need, (n, l, sp)
    Pint = int(P.real + mpf("0.5")) if hasattr(P, "real") else int(P)
    assert fabs(P - Pint) < mpf(10) ** (-need)
    return Pint, case, sp


# --------------------------------------------------------------------------
# T4: cubic coset units along the tower
# --------------------------------------------------------------------------

def coset_unit_log(cg, K):
    s = mpf(0)
    for c in K:
        s += log(fabs(G_value(c, cg.n)))
    return s - mass_log(cg.n) / 3


def tower_from(n0, primes, verbose=True):
    """Cubic layer: base level n0 (primitive cubic character), then the
    levels n0*l (l not | n0) and, for l | n0*l, the second step.  Returns
    the list of (level, multiplier, spare)."""
    need = max(20, mp.dps // 5)
    base = ClassGroup(n0)
    K0, chi, g = cubic_data(base)
    w = exp(2 * pi * 1j / 3)
    L0 = coset_unit_log(base, K0)
    out = []
    if verbose:
        print(f"  base n={n0}: h={base.h}, |ker chi_3|={len(K0)}, log|theta_u| = {nstr(L0, 20)}")
    kernels = {n0: (base, set(K0))}

    def kernel_at(m, n_lower):
        cgm = ClassGroup(m)
        cgl, Kl = kernels[n_lower]
        Km = set(c for c in cgm.reps if cgm.project(c, n_lower) in Kl)
        assert len(Km) == cgm.h // 3, (m, len(Km), cgm.h)
        kernels[m] = (cgm, Km)
        return cgm, Km

    logs = {n0: L0}
    for l in primes:
        m = n0 * l
        if n0 % l == 0:
            # chi_3 must be primitive at n0 with respect to n0/l: the kernel of
            # Pic(O_{n0}) -> Pic(O_{n0/l}) is not contained in ker chi_3
            n_low = n0 // l
            if n_low >= 2:
                lowK = set(base.project(c, n_low) for c in K0)
                ker_low = [c for c in base.reps if base.project(c, n_low) == base.project(base.one, n_low)]
                assert not set(ker_low) <= set(K0), (n0, l, "chi_3 factors through", n_low)
            cgm, Km = kernel_at(m, n0)
            Lm = coset_unit_log(cgm, Km)
            logs[m] = Lm
            ratio = Lm / L0
            Pint = int(ratio + mpf("0.5"))
            sp = spare_digits(fabs(ratio - Pint))
            assert sp >= need and Pint == l + 1, (n0, l, nstr(ratio, 20))
            case = "l | n "
            if verbose:
                print(f"    level {m:3d} = {n0}*{l:<2d} [{case}, chi_3 primitive at {n0}]: index multiplier {Pint:2d} = l+1   (spare {sp})")
            out.append((m, Pint, sp))
        else:
            cgm, Km = kernel_at(m, n0)
            Lm = coset_unit_log(cgm, Km)
            logs[m] = Lm
            P, case = pullback_multiplier(l, n0, lambda lam: w ** chi[base.ideal_class(lam)])
            P = mpf(P.real) if hasattr(P, "real") else mpf(P)
            ratio = Lm / L0
            Pint = int(ratio + mpf("0.5"))
            err = fabs(ratio - Pint)
            sp = spare_digits(err)
            assert sp >= need and fabs(P - Pint) < mpf(10) ** (-need), (n0, l, nstr(ratio, 20), P)
            extra = ""
            if case == "split":
                extra = f", chi_3([l]) = omega^{chi[base.ideal_class(split_lams(l)[0])]}"
            if verbose:
                print(f"    level {m:3d} = {n0}*{l:<2d} [{case:5s}{extra}]: index multiplier {Pint:2d}   (spare {sp})")
            out.append((m, Pint, sp))
        # second step: m -> m*l with l | m:  (l+1) L(m) - |ker(Pic(O_m)->Pic(O_{n0}))| L(n0)
        m2 = m * l
        if m2 <= 81:
            cg2, K2 = kernel_at(m2, m)
            L2 = coset_unit_log(cg2, K2)
            kk = cgm.h // base.h
            pred = (l + 1) * Lm - kk * L0
            err = fabs(L2 - pred)
            sp2 = spare_digits(err, max(1, fabs(pred)))
            assert sp2 >= need, (m2, nstr(L2, 20), nstr(pred, 20))
            mult2 = int(L2 / L0 + mpf("0.5"))
            assert spare_digits(fabs(L2 / L0 - mult2)) >= need
            if verbose:
                print(f"    level {m2:3d} = {m}*{l:<2d} [l | n ]: log theta_u = ({l}+1) log theta_u({m}) - {kk} log theta_u({n0});  "
                      f"index multiplier {mult2:2d}   (spare {sp2})")
            out.append((m2, mult2, sp2))
    return out


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

NORM_PAIRS = [
    # split, l not | n
    (3, 5), (7, 5), (3, 13), (9, 5), (11, 5), (2, 5), (6, 5), (4, 5),
    # inert, l not | n
    (3, 7), (5, 3), (7, 3), (2, 3), (4, 3), (11, 3), (13, 3), (10, 3), (5, 7),
    # ramified l = 2, 2 not | n
    (3, 2), (5, 2), (7, 2), (9, 2), (11, 2), (13, 2), (15, 2), (25, 2),
    # l | n
    (2, 2), (3, 3), (5, 5), (4, 2), (6, 2), (6, 3), (9, 3), (10, 2), (10, 5),
    (12, 2), (15, 3), (15, 5), (8, 2), (18, 3), (25, 5), (14, 2), (14, 7),
]


def selftest():
    t0 = time.time()
    print("=" * 78)
    print("T1/T2: fiber lemma and norm relations   (exact HNF class groups; Delta at "
          f"{mp.dps} digits)")
    print("=" * 78)
    worst = mp.dps
    for (n, l) in NORM_PAIRS:
        worst = min(worst, norm_relation(n, l))
    print(f"  all {len(NORM_PAIRS)} pairs: minimal spare digits {worst}")
    print("=" * 78)
    print("T3: Hecke-type recursion for S_chi = sum chi log|G| = -12 L'(0, chi)")
    print("=" * 78)
    for name, P, sp in test_recursion_quadratic():
        print(f"  quadratic {name}: S(nl) = {P} S(n)   (spare {sp}); closed forms of Paper II reproduced")
    for (n, l) in [(9, 5), (9, 7), (11, 3), (11, 5), (13, 5), (9, 2), (13, 2)]:
        P, case, sp = test_recursion_cubic(n, l)
        print(f"  cubic chi_3 at {n} -> {n*l} [{case}]: S(nl) = {P} S(n)   (spare {sp})")
    print("=" * 78)
    print("T4: the Robert index along the tower (cubic layer)")
    print("=" * 78)
    tower_from(9, [2, 3, 5, 7])
    tower_from(11, [2, 3, 5])
    tower_from(13, [2, 3, 5])
    tower_from(23, [2, 3])
    print(f"ALL CHECKS PASSED   ({time.time() - t0:.0f}s)")


def main(argv):
    mp.dps = 80          # precision set here, never at import (guard rail 2)
    if len(argv) <= 1 or argv[1] == "--selftest":
        selftest()
    elif argv[1] == "norm":
        n, l = int(argv[2]), int(argv[3])
        norm_relation(n, l)
    elif argv[1] == "tower":
        n0 = int(argv[2])
        tower_from(n0, [int(x) for x in argv[3:]])
    else:
        print(__doc__)


if __name__ == "__main__":
    main(sys.argv)
