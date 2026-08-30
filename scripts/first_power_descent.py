"""Machine verification for first-power-descent.md.

Four parts (run all with `python3 scripts/first_power_descent.py all` or
`--selftest`):

(a) `phi`   -- exact construction of the classical modular polynomials
              Phi_m(x, y), m <= MMAX (default 10), from integer q-expansions:
              the b-summed power sums of j((a*tau+b)/d) collapse to integer
              q-series (the cyclotomic sums are Ramanujan-type and evaluate
              in Z), each power sum is peeled against powers of j -- a
              holomorphic SL_2(Z)-invariant function that is O(q) vanishes,
              so the peeling is an identity, verified with slack -- and
              Newton's identities assemble Phi_m.  Validations: integrality,
              symmetry Phi_m(x,y) = Phi_m(y,x), the classical Phi_2, and
              Phi_m(j(m tau), j(tau)) = 0 numerically at a random point.

(b) `omega` -- exact verification (rational arithmetic in B = Q(i, sqrt N),
              no rounding anywhere) that omega_f = eps nu^2/(r0 mu^2) = 1
              at every primitive class of every odd n <= 21, together with
              the three identities of Lemma 3.2 (N o varsigma = (s0/2) g;
              m2 = conj(sigma1)/conj(sigma2) matches the sigma-circle Gram;
              nu = a omega0 / conj(sigma2) satisfies the lattice equation
              b2 = nu r^{-1} b1).

(c) `uf`    -- numerical verification of Theorem 3.3's closed form
              u_f = Phi_y/Phi_x(beta1, beta2), Phi = Phi_{(n-1)/2},
              beta1 = j(b1), beta2 = j(r^{-1} b1), against the canonical-
              matrix route u = eps Theta(X) of moduli-invariants.md, to
              >= 65 digits at every primitive class, odd n <= 13.

(d) `pin`   -- the exact computation of Pi_n = prod_f (x - u_f) of Theorem
              4.2: H_D by certified integer rounding, the pairing root
              beta2(t) as the degree-1 gcd of (H_D(y), Phi_{r0}(t, y)) over
              F1 = Q[t]/(H_D), u(t) = Phi_y/Phi_x in F1, and Pi_n as the
              exact characteristic polynomial of multiplication by u(t).
              The primitive integer form is matched against the published
              level polynomials Q_n (moduli-invariants.md 5.9), proved
              irreducible and squarefree by exact factorization.

Guard rails (CLAUDE.md): mp.dps is set inside functions, never at import
time; integer identifications use the absolute-error criterion with
explicit spare digits; everything finite is exact integer/Fraction
arithmetic.
"""
import sys
import time
from fractions import Fraction
from math import gcd, isqrt

sys.path.insert(0, 'scripts')


# ====================================================================
# (a)  Exact modular polynomials Phi_m from integer q-expansions
# ====================================================================

def j_q_coefficients(nmax):
    """Integer coefficients c(-1..nmax) of j = 1/q + 744 + 196884 q + ...

    j = E4^3 / Delta with E4 = 1 + 240 sum sigma_3(k) q^k and
    Delta/q = prod (1-q^k)^24, all in exact integer arithmetic.
    Returns the list [c(-1), c(0), ..., c(nmax)].
    """
    N = nmax + 2
    # E4
    e4 = [0] * N
    e4[0] = 1
    for k in range(1, N):
        s3 = sum(d ** 3 for d in range(1, k + 1) if k % d == 0)
        e4[k] = 240 * s3
    # eta-product prod(1 - q^k) by the pentagonal number theorem
    eta = [0] * N
    eta[0] = 1
    k = 1
    while True:
        g1 = k * (3 * k - 1) // 2
        g2 = k * (3 * k + 1) // 2
        if g1 >= N and g2 >= N:
            break
        s = -1 if k % 2 else 1
        if g1 < N:
            eta[g1] += s
        if g2 < N:
            eta[g2] += s
        k += 1

    def mul(a, b):
        c = [0] * N
        for i, ai in enumerate(a):
            if ai:
                for jj, bj in enumerate(b):
                    if i + jj >= N:
                        break
                    if bj:
                        c[i + jj] += ai * bj
        return c

    p2 = mul(eta, eta)
    p4 = mul(p2, p2)
    p8 = mul(p4, p4)
    p16 = mul(p8, p8)
    dq = mul(p8, p16)             # Delta/q = prod(1-q^k)^24
    # invert Delta/q (leading coefficient 1)
    inv = [0] * N
    inv[0] = 1
    for nn in range(1, N):
        inv[nn] = -sum(dq[k] * inv[nn - k] for k in range(1, nn + 1))
    e43 = mul(mul(e4, e4), e4)
    jq = mul(e43, inv)            # j * q as a series in q
    # c(m) = coefficient of q^{m}, m = -1..nmax  <->  jq[m+1]
    return jq[:nmax + 2]


def _restricted_b_sum(a, d, midx):
    """S = sum_{b mod d, gcd(b, gcd(a,d)) = 1} zeta_d^{b*midx}, an integer."""
    g = gcd(a, d)
    total = 0
    e = 1
    while e * e <= g or e <= g:
        if e > g:
            break
        if g % e == 0:
            # mu(e)
            mu, x, sq = 1, e, False
            f = 2
            xx = e
            while f * f <= xx:
                if xx % f == 0:
                    xx //= f
                    if xx % f == 0:
                        sq = True
                        break
                    mu = -mu
                f += 1
            if not sq:
                if xx > 1:
                    mu = -mu
                de = d // e
                if midx % de == 0:
                    total += mu * de
        e += 1
    return total


def build_phi(m, slack=12, verbose=False):
    """Exact Phi_m(x, y) as {(i, jexp): int} for x^i y^jexp; x-monic deg psi.

    Method: for each divisor block (a, d), ad = m, the b-summed power sums
    P_k^{(a,d)}(tau) = sum_b j((a tau + b)/d)^k have integer q-expansions
    after the Ramanujan collapse; summing blocks gives the full power sum
    p_k, a polynomial in j found by exact peeling (tail asserted zero
    through q^slack); Newton's identities give the elementary symmetric
    functions e_k(y) in Z[y] and Phi_m(x,y) = x^psi + sum (-1)^k e_k x^(psi-k).
    """
    blocks = []
    psi = 0
    for a in range(1, m + 1):
        if m % a == 0:
            d = m // a
            g = gcd(a, d)
            nb = sum(1 for b in range(d) if gcd(b, g) == 1)
            blocks.append((a, d, nb))
            psi += nb
    kmax = psi
    emax = kmax * m
    # base series lengths per block: exponent a*midx/d up to slack, with
    # buffer kmax for the iterated-power truncation boundary; the j-powers
    # for the peeling need coefficients up to slack + emax (a factor at
    # exponent -(e-1) pairs with one at slack + e - 1).
    need_c = max([kmax + (slack * d) // a + kmax + 2 for a, d, _ in blocks]
                 + [slack + emax + 2])
    c = j_q_coefficients(need_c)          # c[t] = coefficient of q^{t-1}

    # j series (as dict exponent -> int) for the peeling powers
    jser = {t - 1: c[t] for t in range(len(c))
            if t - 1 <= slack + emax and c[t]}

    def ser_mul(s1, s2, lo, hi):
        out = {}
        for e1, v1 in s1.items():
            for e2, v2 in s2.items():
                e = e1 + e2
                if lo <= e <= hi:
                    out[e] = out.get(e, 0) + v1 * v2
        return {e: v for e, v in out.items() if v}

    # jpow[e] exact for exponents <= slack + (emax - e); this cascades so
    # that every jpow[e] is exact through exponent slack.
    jpow = [None] * (emax + 1)
    jpow[0] = {0: 1}
    for e in range(1, emax + 1):
        jpow[e] = ser_mul(jpow[e - 1], jser, -e, slack + emax - e)

    # power sums p_k as q-series
    pk_series = [None] * (kmax + 1)
    for a, d, _nb in blocks:
        L = kmax + (slack * d) // a + 2
        base = c[:L + 1]                  # index t <-> midx = t - 1
        base_d = {t - 1: base[t] for t in range(len(base)) if base[t]}
        pw_d = dict(base_d)
        for k in range(1, kmax + 1):
            if k > 1:
                hi = (slack * d) // a + 1 + (kmax - k)
                pw_d = ser_mul(pw_d, base_d, -k, hi)
            # collapse the b-sum: exponent a*midx/d, weight S(midx)
            tgt = pk_series[k] if pk_series[k] is not None else {}
            for midx, v in pw_d.items():
                S = _restricted_b_sum(a, d, midx)
                if S:
                    num = a * midx
                    assert num % d == 0, (m, a, d, midx)
                    e = num // d
                    if e <= slack:
                        tgt[e] = tgt.get(e, 0) + v * S
            pk_series[k] = tgt

    # peel each p_k against powers of j
    def peel(ser, deg):
        ser = dict(ser)
        coeffs = [0] * (deg + 1)          # coeffs[e] * j^e
        for e in range(deg, -1, -1):
            al = ser.get(-e, 0)
            coeffs[e] = al
            if al:
                for ee, v in jpow[e].items():
                    ser[ee] = ser.get(ee, 0) - al * v
                ser = {k2: v2 for k2, v2 in ser.items() if v2}
        tail = {e: v for e, v in ser.items() if v and e <= slack}
        assert not tail, (m, 'nonzero peeling tail', sorted(tail)[:3])
        return coeffs

    p_poly = [None] + [peel(pk_series[k], k * m) for k in range(1, kmax + 1)]

    # Newton's identities over Q[j]:  k e_k = sum_{i=1}^k (-1)^{i-1} e_{k-i} p_i
    def poly_mul(p, q):
        out = [Fraction(0)] * (len(p) + len(q) - 1)
        for i, pi in enumerate(p):
            if pi:
                for jj, qj in enumerate(q):
                    if qj:
                        out[i + jj] += Fraction(pi) * qj
        return out

    def poly_add(p, q, scal=1):
        n = max(len(p), len(q))
        out = [Fraction(0)] * n
        for i in range(n):
            if i < len(p):
                out[i] += p[i]
            if i < len(q):
                out[i] += Fraction(scal) * q[i]
        return out

    e_poly = [[Fraction(1)]]
    for k in range(1, kmax + 1):
        acc = [Fraction(0)]
        for i in range(1, k + 1):
            term = poly_mul(e_poly[k - i], [Fraction(x) for x in p_poly[i]])
            acc = poly_add(acc, term, (-1) ** (i - 1))
        ek = [x / k for x in acc]
        for x in ek:
            assert x.denominator == 1, (m, k, 'non-integral e_k')
        e_poly.append(ek)

    phi = {(psi, 0): 1}
    for k in range(1, kmax + 1):
        for jexp, v in enumerate(e_poly[k]):
            v = int(v)
            if v:
                phi[(psi - k, jexp)] = phi.get((psi - k, jexp), 0) + ((-1) ** k) * v
    phi = {k2: v for k2, v in phi.items() if v}
    if verbose:
        print(f"    Phi_{m}: degree {psi}, {len(phi)} monomials")
    return phi, psi


PHI2_CLASSICAL = {
    (3, 0): 1, (0, 3): 1, (2, 2): -1,
    (2, 1): 1488, (1, 2): 1488,
    (2, 0): -162000, (0, 2): -162000,
    (1, 1): 40773375,
    (1, 0): 8748000000, (0, 1): 8748000000,
    (0, 0): -157464000000000,
}


def phi_eval(phi, x, y):
    return sum(v * x ** i * y ** jj for (i, jj), v in phi.items())


def phi_partial(phi, which):
    out = {}
    for (i, jj), v in phi.items():
        if which == 'x' and i > 0:
            out[(i - 1, jj)] = out.get((i - 1, jj), 0) + i * v
        if which == 'y' and jj > 0:
            out[(i, jj - 1)] = out.get((i, jj - 1), 0) + jj * v
    return {k: v for k, v in out.items() if v}


def run_phi(mmax=10, verbose=True):
    """Part (a): build and validate Phi_m for m <= mmax."""
    from mpmath import mp, mpc, fabs
    t0 = time.time()
    phis = {}
    for m in range(1, mmax + 1):
        phi, psi = build_phi(m, verbose=verbose)
        # symmetry (m > 1)
        if m > 1:
            assert all(phi.get((jj, i), 0) == v for (i, jj), v in phi.items()), \
                (m, 'Phi_m not symmetric')
        phis[m] = phi
    assert phis[1] == {(1, 0): 1, (0, 1): -1}, 'Phi_1 != x - y'
    assert phis.get(2) == PHI2_CLASSICAL, 'Phi_2 != classical polynomial'
    # numeric root check Phi_m(j(m tau), j(tau)) = 0 at a generic point
    mp.dps = 120
    import moduli_invariants as MI
    tau = mpc('0.31377', '1.07231')
    for m in range(1, mmax + 1):
        jm, j1 = MI.J(m * tau), MI.J(tau)
        val = phi_eval(phis[m], jm, j1)
        scale = max(fabs(v * jm ** i * j1 ** jj)
                    for (i, jj), v in phis[m].items())
        rel = fabs(val) / scale
        assert rel < mp.mpf(10) ** (-(mp.dps - 30)), (m, 'root check', rel)
    if verbose:
        print(f"  (a) Phi_m exact for m <= {mmax}: integrality, symmetry, "
              f"classical Phi_2, numeric root checks at 120 digits "
              f"[{time.time()-t0:.1f} s]")
    return phis


# ====================================================================
# Exact arithmetic in K = Q(sqrt D) and B = Q(i, sqrt N)
# ====================================================================

class K2:
    """x + y sqrt(D) with x, y in Q (D fixed per instance)."""
    __slots__ = ('x', 'y', 'D')

    def __init__(self, x, y, D):
        self.x, self.y, self.D = Fraction(x), Fraction(y), D

    def __add__(s, o):
        return K2(s.x + o.x, s.y + o.y, s.D)

    def __sub__(s, o):
        return K2(s.x - o.x, s.y - o.y, s.D)

    def __mul__(s, o):
        if isinstance(o, (int, Fraction)):
            return K2(s.x * o, s.y * o, s.D)
        return K2(s.x * o.x + s.D * s.y * o.y, s.x * o.y + s.y * o.x, s.D)

    __rmul__ = __mul__

    def conj(s):                      # sqrt(D) -> -sqrt(D)
        return K2(s.x, -s.y, s.D)

    def inv(s):
        nrm = s.x * s.x - s.D * s.y * s.y
        assert nrm != 0
        return K2(s.x / nrm, -s.y / nrm, s.D)

    def __eq__(s, o):
        return s.x == o.x and s.y == o.y

    def __repr__(s):
        return f"({s.x}+{s.y}*sqrtD)"


class B4:
    """z + w*i with z, w in K = Q(sqrt D); D = -N < 0, sqrt D = i sqrt N."""
    __slots__ = ('z', 'w')

    def __init__(self, z, w):
        self.z, self.w = z, w

    def __add__(s, o):
        return B4(s.z + o.z, s.w + o.w)

    def __sub__(s, o):
        return B4(s.z - o.z, s.w - o.w)

    def __mul__(s, o):
        return B4(s.z * o.z - s.w * o.w, s.z * o.w + s.w * o.z)

    def inv(s):
        # (z + wi)^{-1} = (z - wi)/(z^2 + w^2)
        den = s.z * s.z + s.w * s.w    # in K
        di = den.inv()
        return B4(s.z * di, K2(0, 0, s.z.D) - s.w * di)

    def __eq__(s, o):
        return s.z == o.z and s.w == o.w

    def __repr__(s):
        return f"[{s.z} + ({s.w})*i]"


def lattice_hnf_q(gens):
    """Canonical basis of the Z-lattice spanned by pairs of Fractions."""
    from fractions import Fraction as F
    den = 1
    for a, b in gens:
        den = den * F(a).denominator // gcd(den, F(a).denominator)
        den = den * F(b).denominator // gcd(den, F(b).denominator)
    rows = [[int(F(a) * den), int(F(b) * den)] for a, b in gens]
    rows = [r for r in rows if r != [0, 0]]
    while True:
        nz = [r for r in rows if r[0] != 0]
        if len(nz) <= 1:
            break
        nz.sort(key=lambda r: abs(r[0]))
        p = nz[0]
        for r in nz[1:]:
            qq = r[0] // p[0]
            r[0] -= qq * p[0]
            r[1] -= qq * p[1]
        rows = [r for r in rows if r != [0, 0]]
    top = [r for r in rows if r[0] != 0]
    rest = [r for r in rows if r[0] == 0 and r[1] != 0]
    while len(rest) > 1:
        rest.sort(key=lambda r: abs(r[1]))
        p = rest[0]
        for r in rest[1:]:
            qq = r[1] // p[1]
            r[1] -= qq * p[1]
        rest = [r for r in rest if r[1] != 0]
    t = top[0] if top else [0, 0]
    s = rest[0] if rest else [0, 0]
    if t[0] < 0:
        t = [-t[0], -t[1]]
    if s[1] < 0:
        s = [0, -s[1]]
    if s[1]:
        t[1] %= s[1]
    return (F(t[0], den), F(t[1], den)), (F(s[0], den), F(s[1], den))


def k2_lattice(gens):
    """HNF key for the Z-lattice spanned by K2 elements (coords x, y)."""
    return lattice_hnf_q([(g.x, g.y) for g in gens])


def run_omega(nmax=21, verbose=True):
    """Part (b): omega_f = 1 exactly at every primitive class, odd n <= nmax."""
    from involution_classmap import classes_of_disc, is_primitive
    from involution_experiments import inv_sl2, M_of_X
    from proof_check import build_P
    import moduli_invariants as MI
    t0 = time.time()
    total = 0
    for n in range(3, nmax + 1, 2):
        N = n * n - 1
        D = 1 - n * n
        r0, s0 = (n - 1) // 2, (n + 1) // 2
        for f in classes_of_disc(D):
            if not is_primitive(f):
                continue
            a, b, cc = f
            P, (w1, w2) = build_P(n, f)
            (s1, t1), (s2, t2) = w1, w2
            assert s2 == 0, 'K_basis not in HNF normal form'
            X = inv_sl2(P)
            # sigma_k = s_k s0 + t_k omega0, omega0 = sqrt(D)/2
            sg1 = K2(Fraction(s1 * s0), Fraction(t1, 2), D)
            sg2 = K2(Fraction(s2 * s0), Fraction(t2, 2), D)

            # Lemma 3.2(1): N(varsigma(w)) = (s0/2) g(w), polarized
            def gform(u, v):
                return (n + 1) * u[0] * v[0] + (n - 1) * u[1] * v[1]
            for u, v, su, sv in ((w1, w1, sg1, sg1), (w2, w2, sg2, sg2),
                                 (w1, w2, sg1, sg2)):
                nrm = su * sv.conj()
                half = nrm + nrm.conj()          # trace = 2*Re
                assert half.y == 0 and half.x == Fraction(s0 * gform(u, v)), \
                    (n, f, 'Lemma 3.2(1)')

            # m2 = conj(sg1)/conj(sg2), and it matches the sigma-circle
            m2 = sg1.conj() * sg2.conj().inv()
            Msig = M_of_X(MI.sigma_mat(X))
            (A2, B2), (_, C2) = Msig
            A2v, bx, by = A2[0], B2[0], B2[1]
            if A2v < 0:
                A2v, bx, by, C2v = -A2v, -bx, -by, -C2[0]
            else:
                C2v = C2[0]
            x2, y2 = -bx, -by                  # zeta' = x' + y' i, y' = -n
            assert y2 == -n, (n, f, 'sigma-circle not at level -n')
            # lower hyperbolic center = (x' - i sqrt N)/A' = (x' - sqrtD)/A2v
            m2_direct = K2(Fraction(x2, A2v), Fraction(-1, A2v), D)
            assert m2 == m2_direct, (n, f, 'Lemma 3.2(2)')

            # nu = a*omega0/conj(sg2) satisfies b2 = nu * r^{-1} b1
            omega0 = K2(0, Fraction(1, 2), D)
            nu = omega0 * a * sg2.conj().inv()
            one = K2(1, 0, D)
            m1 = K2(Fraction(-b, 2 * a), Fraction(1, 2 * a), D)
            b1 = [one, m1]
            rid = [K2(r0, 0, D), omega0]
            rinv_b1 = [g1 * g2 * Fraction(1, r0) for g1 in rid for g2 in b1]
            b2 = [one, K2(0, 0, D) - m2]
            lhs = k2_lattice(b2)
            rhs = k2_lattice([nu * g for g in rinv_b1])
            assert lhs == rhs, (n, f, 'Lemma 3.2(3)')

            # omega_f = eps nu^2 / (r0 mu^2) = 1 in B
            KD = lambda x, y: K2(x, y, D)
            m2B = B4(m2, KD(0, 0))
            crow, drow = X[1]
            cB = B4(KD(crow[0], 0), KD(crow[1], 0))
            dB = B4(KD(drow[0], 0), KD(drow[1], 0))
            mu = cB * m2B + dB
            # eps = n + sqrt N = n + (-sqrtD) * i
            eps = B4(KD(n, 0), KD(0, -1))
            nuB = B4(nu, KD(0, 0))
            om = eps * nuB * nuB * (mu * mu * B4(KD(r0, 0), KD(0, 0))).inv()
            assert om == B4(KD(1, 0), KD(0, 0)), (n, f, 'omega_f != 1')
            total += 1
    if verbose:
        print(f"  (b) omega_f = 1 exact (with Lemma 3.2) on {total} primitive "
              f"classes, odd n <= {nmax} [{time.time()-t0:.1f} s]")
    return total


# ====================================================================
# (c)  u_f = Phi_y/Phi_x(beta1, beta2) numerically
# ====================================================================

def twisted_tau(n, f):
    """Exact basis ratio tau2 (as K2 pair) of r^{-1} b1 for class f, plus
    the K2 basis (al, be) with tau2 = be/al, Im tau2 > 0."""
    D = 1 - n * n
    r0 = (n - 1) // 2
    a, b, _c = f
    one = K2(1, 0, D)
    omega0 = K2(0, Fraction(1, 2), D)
    m1 = K2(Fraction(-b, 2 * a), Fraction(1, 2 * a), D)
    gens = [g1 * g2 * Fraction(1, r0) for g1 in (K2(r0, 0, D), omega0)
            for g2 in (one, m1)]
    (t0x, t0y), (s0x, s0y) = k2_lattice(gens)
    al = K2(t0x, t0y, D)
    be = K2(s0x, s0y, D)
    # tau = be/al; ensure Im > 0 (y-coordinate of be/al > 0 since sqrtD = i sqrtN)
    tau = be * al.inv()
    if tau.y < 0:
        be = K2(0, 0, D) - be
        tau = K2(tau.x, -tau.y, D)
    return tau, al, be


def run_uf(nmax=13, dps=110, target=65, phis=None, verbose=True):
    """Part (c): route A (eps*Theta on the canonical matrix) vs the
    correspondence-derivative closed form, >= `target` digits."""
    from mpmath import mp, mpc, mpf, fabs
    from mpmath import sqrt as msqrt
    mp.dps = dps
    from involution_classmap import classes_of_disc, is_primitive
    from involution_experiments import inv_sl2
    from proof_check import build_P
    import moduli_invariants as MI
    if phis is None:
        phis = {}
        for m in set((n - 1) // 2 for n in range(3, nmax + 1, 2)):
            phis[m], _ = build_phi(m)
    t0 = time.time()
    worst = mp.inf
    total = 0
    for n in range(3, nmax + 1, 2):
        D = 1 - n * n
        r0 = (n - 1) // 2
        eps = n + msqrt(mpf(n * n - 1))
        for f in classes_of_disc(D):
            if not is_primitive(f):
                continue
            X = inv_sl2(build_P(n, f)[0])
            uA = eps * MI.theta_integral(X)[0]
            # correspondence route
            a, b, _c = f
            m1 = mpc(mpf(-b) / (2 * a), msqrt(mpf(n * n - 1)) / (2 * a))
            beta1 = MI.J(m1)
            tau2, _, _ = twisted_tau(n, f)

            def mpq(fr):
                return mpf(fr.numerator) / mpf(fr.denominator)

            t2 = mpc(mpq(tau2.x), mpq(tau2.y) * msqrt(mpf(n * n - 1)))
            beta2 = MI.J(t2)
            phi = phis[r0]
            num = phi_eval(phi_partial(phi, 'y'), beta1, beta2)
            den = phi_eval(phi_partial(phi, 'x'), beta1, beta2)
            uB = num / den
            err = fabs(uA - uB) / max(mpf(1), fabs(uA))
            digs = -mp.log10(err) if err > 0 else mp.inf
            worst = min(worst, digs)
            assert digs >= target, (n, f, float(digs))
            total += 1
    if verbose:
        w = 'inf' if worst == mp.inf else f"{float(worst):.0f}"
        print(f"  (c) u_f = Phi_y/Phi_x(beta1, beta2) vs eps*Theta route: "
              f"{total} classes, odd n <= {nmax}, agreement >= {w} digits "
              f"(target {target}) [{time.time()-t0:.1f} s]")
    return total


# ====================================================================
# (d)  Exact Pi_n via F1 = Q[t]/(H_D)
# ====================================================================

# published level polynomials Q_n (moduli-invariants.md 5.9), leading first
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
}


def hilbert_class_poly(D, dps=120, max_err=None):
    """H_D in Z[t] by certified rounding of prod (t - j(tau_f))."""
    from mpmath import mp, mpc, mpf, fabs, nint
    from mpmath import sqrt as msqrt
    mp.dps = dps
    import moduli_invariants as MI
    max_err = mpf(10) ** (-40) if max_err is None else max_err
    forms = []
    b = D % 2
    while b * b <= -D // 3 + 1:
        num = b * b - D
        if num % 4 == 0:
            ac = num // 4
            for a in range(max(b, 1), isqrt(ac) + 1):
                if ac % a == 0:
                    c = ac // a
                    if a > c or (b > a):
                        continue
                    if gcd(gcd(a, b), c) != 1:
                        continue
                    forms.append((a, b, c))
                    if 0 < b < a < c:
                        forms.append((a, -b, c))
        b += 2
    js = [MI.J(mpc(mpf(-bb) / (2 * aa), msqrt(mpf(4 * aa * cc - bb * bb)) / (2 * aa)))
          for aa, bb, cc in forms]
    co = [mpc(1)]
    for r in js:
        new = [mpc(0)] * (len(co) + 1)
        for i, ci in enumerate(co):
            new[i] += ci * (-r)
            new[i + 1] += ci
        co = new
    ints = []
    for ci in co:
        ii = int(nint(ci.real))
        assert fabs(ci - ii) < max_err, ('H_D rounding', D, fabs(ci - ii))
        ints.append(ii)
    return list(reversed(ints)), len(js)     # highest-degree first, degree


class F1:
    """Q[t]/(H) arithmetic; H monic integer, irreducible (verified)."""

    def __init__(self, H):
        self.H = [Fraction(x) for x in H]      # highest first, monic
        self.h = len(H) - 1

    def elem(self, coeffs_low):
        v = [Fraction(x) for x in coeffs_low] + \
            [Fraction(0)] * (self.h - len(coeffs_low))
        return tuple(v[:self.h])

    def zero(self):
        return self.elem([])

    def one(self):
        return self.elem([1])

    def t(self):
        if self.h >= 2:
            return self.elem([0, 1])
        # h = 1: t = -c mod (t + c)
        return self.elem([-self.H[1]])

    def add(self, a, b):
        return tuple(x + y for x, y in zip(a, b))

    def sub(self, a, b):
        return tuple(x - y for x, y in zip(a, b))

    def scal(self, a, s):
        return tuple(x * s for x in a)

    def mul(self, a, b):
        prod = [Fraction(0)] * (2 * self.h - 1)
        for i, ai in enumerate(a):
            if ai:
                for j, bj in enumerate(b):
                    if bj:
                        prod[i + j] += ai * bj
        # reduce mod H (monic, highest first)
        Hlow = list(reversed(self.H))          # lowest first
        for d in range(len(prod) - 1, self.h - 1, -1):
            c = prod[d]
            if c:
                for k in range(self.h + 1):
                    prod[d - k] -= c * Hlow[self.h - k]
        return tuple(prod[:self.h])

    def is_zero(self, a):
        return all(x == 0 for x in a)

    def inv(self, a):
        # extended Euclid between lift(a) and H in Q[t]
        def pdeg(p):
            d = len(p) - 1
            while d >= 0 and p[d] == 0:
                d -= 1
            return d

        def pdivmod(num, den):
            num = list(num)
            dd, dn = pdeg(den), pdeg(num)
            q = [Fraction(0)] * (max(dn - dd + 1, 1))
            while dn >= dd >= 0:
                c = num[dn] / den[dd]
                q[dn - dd] = c
                for k in range(dd + 1):
                    num[dn - dd + k] -= c * den[k]
                dn = pdeg(num)
            return q, num

        r0 = list(reversed(self.H))
        r1 = list(a)
        s0, s1 = [Fraction(0)], [Fraction(1)]
        while pdeg(r1) > 0:
            q, r = pdivmod(r0, r1)
            r0, r1 = r1, r
            # s = s0 - q*s1
            qs = [Fraction(0)] * (pdeg(q) + pdeg(s1) + 2)
            for i in range(len(q)):
                if q[i]:
                    for j in range(len(s1)):
                        if s1[j]:
                            qs[i + j] += q[i] * s1[j]
            ns = [ (s0[i] if i < len(s0) else Fraction(0)) -
                   (qs[i] if i < len(qs) else Fraction(0))
                   for i in range(max(len(s0), len(qs))) ]
            s0, s1 = s1, ns
        c = r1[pdeg(r1)]
        assert pdeg(r1) == 0 and c != 0, 'element not invertible'
        inv = [x / c for x in s1]
        return self.elem(inv)


def polyF_gcd(F, A, B):
    """Monic gcd of polynomials over the field F1 (lists of F1-elems, low first)."""
    def deg(p):
        d = len(p) - 1
        while d >= 0 and F.is_zero(p[d]):
            d -= 1
        return d

    def make_monic(p):
        d = deg(p)
        li = F.inv(p[d])
        return [F.mul(x, li) for x in p[:d + 1]]

    A, B = list(A), list(B)
    while deg(B) >= 0:
        B = make_monic(B)
        db, da = deg(B), deg(A)
        while da >= db:
            lead = A[da]
            if not F.is_zero(lead):
                for k in range(db + 1):
                    A[da - db + k] = F.sub(A[da - db + k], F.mul(lead, B[k]))
            da = deg(A)
        A, B = B, A[:max(da + 1, 0)]
    return make_monic(A)


def run_pin(levels=(3, 5, 7, 9, 11, 13), dps=140, phis=None, verbose=True):
    """Part (d): exact Pi_n, table match, irreducibility, squarefreeness."""
    import sympy as sp
    t00 = time.time()
    if phis is None:
        phis = {}
        for m in set((n - 1) // 2 for n in levels):
            phis[m], _ = build_phi(m)
    x = sp.symbols('x')
    for n in levels:
        t0 = time.time()
        D = 1 - n * n
        r0 = (n - 1) // 2
        H, h = hilbert_class_poly(D, dps=dps)
        HP = sp.Poly(H, x)
        assert HP.is_irreducible, (n, 'H_D reducible?!')
        F = F1(H)
        phi = phis[r0]
        # Phi_{r0}(t, y) as polynomial in y over F1
        degy = max(jj for (_i, jj) in phi)
        phi_y = [F.zero() for _ in range(degy + 1)]
        tpowers = [F.one()]
        degx = max(i for (i, _jj) in phi)
        for _ in range(degx):
            tpowers.append(F.mul(tpowers[-1], F.t()))
        for (i, jj), v in phi.items():
            phi_y[jj] = F.add(phi_y[jj], F.scal(tpowers[i], Fraction(v)))
        # H_D(y) over F1
        Hlow = list(reversed(H))
        H_y = [F.elem([Fraction(c)]) for c in Hlow]
        g = polyF_gcd(F, H_y, phi_y)
        assert len(g) == 2, (n, 'pairing gcd degree != 1', len(g))
        beta2 = F.sub(F.zero(), g[0])          # monic y + g0 -> root -g0
        # u(t) = Phi_y / Phi_x at (t, beta2)
        dphix, dphiy = phi_partial(phi, 'x'), phi_partial(phi, 'y')

        def eval_F(p2):
            b2pow = [F.one()]
            degj = max(jj for (_i, jj) in p2) if p2 else 0
            for _ in range(degj):
                b2pow.append(F.mul(b2pow[-1], beta2))
            acc = F.zero()
            for (i, jj), v in p2.items():
                acc = F.add(acc, F.scal(F.mul(tpowers[i], b2pow[jj]),
                                        Fraction(v)))
            return acc

        u = F.mul(eval_F(dphiy), F.inv(eval_F(dphix)))
        # characteristic polynomial of multiplication by u on F1
        cols = []
        e = F.one()
        for k in range(h):
            if k:
                e = F.mul(e, F.t())
            cols.append(F.mul(u, e))
        M = sp.Matrix(h, h, lambda i, j: sp.Rational(cols[j][i]))
        Pi = M.charpoly(x).as_expr()
        Pip = sp.Poly(Pi, x)
        # primitive integer form
        co = Pip.all_coeffs()
        den = 1
        for c in co:
            den = sp.ilcm(den, sp.fraction(sp.Rational(c))[1])
        ico = [int(c * den) for c in co]
        cont = 0
        for c in ico:
            cont = gcd(cont, c)
        ico = [c // cont for c in ico]
        if ico[0] < 0:
            ico = [-c for c in ico]
        assert ico == PUBLISHED_Q[n], (n, 'mismatch with published Q_n', ico)
        QP = sp.Poly(ico, x)
        assert QP.is_irreducible, (n, 'Q_n reducible')
        assert sp.gcd(QP, QP.diff(x)).total_degree() == 0, (n, 'not squarefree')
        if verbose:
            print(f"  (d) n={n}: H_{D} certified (deg {h}), pairing gcd "
                  f"degree 1, exact Pi_{n} == published Q_{n}, irreducible, "
                  f"squarefree [{time.time()-t0:.1f} s]")
    if verbose:
        print(f"  (d) total [{time.time()-t00:.1f} s]")


# ====================================================================

def main(argv):
    mode = argv[1] if len(argv) > 1 else 'all'
    mmax = 10
    if '--mmax' in argv:
        mmax = int(argv[argv.index('--mmax') + 1])
    t0 = time.time()
    if mode in ('all', '--selftest', 'selftest'):
        print("first_power_descent.py selftest")
        phis = run_phi(mmax=mmax, verbose=True)
        run_omega(nmax=21)
        run_uf(nmax=13, phis={m: phis[m] for m in phis if m <= 6})
        run_pin(phis={m: phis[m] for m in phis if m <= 6})
        print(f"ALL PASS [{time.time()-t0:.1f} s total]")
    elif mode == 'phi':
        run_phi(mmax=mmax)
    elif mode == 'omega':
        run_omega(nmax=21)
    elif mode == 'uf':
        run_uf(nmax=13)
    elif mode == 'pin':
        run_pin()
    else:
        print(__doc__)


if __name__ == '__main__':
    main(sys.argv)
