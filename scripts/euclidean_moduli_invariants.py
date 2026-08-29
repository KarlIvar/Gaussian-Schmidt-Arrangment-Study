"""Euclidean moduli invariants of Schmidt disks: numerical laboratory for
euclidean-moduli-invariants.md.

Setting: X in SL_2(Z[i]) with X(H) a *bounded* disk D (a Schmidt disk).
With bottom row (c, d):
    curvature  kappa(D) = 2 Im(conj(c) d) = 2n > 0,
    w  = d/c = -X^{-1}(infinity)  in  H      (the reflected pole),
    Lambda = Z c + Z d   -- a primitive index-n sublattice of Z[i],
    zeta = i (a conj(d) - b conj(c)),  center(D) = zeta / 2n.

Invariants of the double coset  N(Z[i]) X SL_2(Z),  N(Z[i]) = translations:
    five:  the circle mod translation (3 real) and beta = j(w) (2 real);
    sixth: Theta(X) = j'(w) / c^2  =  -Res_{z = X^{-1}(inf)} X(z) * j'(w),
    normalized  u = Theta / Omega,  Omega = Gamma(1/4)^4 / (8 pi^2)
    (= pomega^2/pi, pomega the lemniscate constant).

Experiments:
  A. exact structure: disks <-> primitive sublattices <-> ring classes of
     O_n = Z + nZ[i]; N_e(n) = 2 h(-4n^2); conductor exactly n; {L, iL} pairs.
  B. invariance of Theta; the fiber (rotations about the pole): |Theta|
     constant, arg linear at rate 2; period identities; u^6 closed form.
  C. ring class polynomials: prod over disks (x - j(w_D)) = H_{-4n^2}(x)^2,
     integer-certified; trace slicing  sum_{g|n} Tr_E(n/g) = t(4n^2) (Zagier).
  D. the phase units: laws, reality pattern, the integer polynomials
     P2_n(y) = prod_c (y - (lambda_n v_c)^2)  and  P6_n(y) (sixth powers),
     minimal lambda_n, factorizations, Delta-mass, Galois groups.

Usage:
    python3 scripts/euclidean_moduli_invariants.py            # A + B + C
    python3 scripts/euclidean_moduli_invariants.py phase      # D (slower)
    python3 scripts/euclidean_moduli_invariants.py all
Requires mpmath; experiment D also uses sympy (factorization, Galois groups).
"""
import sys
from math import gcd
from fractions import Fraction

from mpmath import (mp, mpf, mpc, fabs, nstr, jtheta, exp, pi, sqrt as msqrt,
                    arg, nint, gamma, quad, log10, expm, matrix, floor)

MAXN = 13          # phase-experiment levels 2..MAXN
MAXN_STRUCT = 24   # exact-structure levels

# ======================= Gaussian integers (exact) =======================
# represented as pairs (x, y) = x + y i

def gadd(u, v): return (u[0] + v[0], u[1] + v[1])
def gsub(u, v): return (u[0] - v[0], u[1] - v[1])
def gmul(u, v): return (u[0]*v[0] - u[1]*v[1], u[0]*v[1] + u[1]*v[0])
def gconj(u): return (u[0], -u[1])
def gneg(u): return (-u[0], -u[1])
def gnorm(u): return u[0]*u[0] + u[1]*u[1]

UNITS = [(1, 0), (-1, 0), (0, 1), (0, -1)]

def gdivmod(u, v):
    n = gnorm(v)
    p = gmul(u, gconj(v))
    q = ((2*p[0] + n) // (2*n), (2*p[1] + n) // (2*n))
    return q, gsub(u, gmul(q, v))

def gxgcd(u, v):
    r0, r1, s0, s1, t0, t1 = u, v, (1, 0), (0, 0), (0, 0), (1, 0)
    while r1 != (0, 0):
        q, r = gdivmod(r0, r1)
        r0, r1 = r1, r
        s0, s1 = s1, gsub(s0, gmul(q, s1))
        t0, t1 = t1, gsub(t0, gmul(q, t1))
    return r0, s0, t0

# ================= primitive index-n sublattices of Z[i] =================
# Hermite data (a, b, d):  Lambda = Z(b + a i) + Z d,  a d = n,  0 <= b < d.

def sublattices(n):
    for a in range(1, n + 1):
        if n % a == 0:
            d = n // a
            for b in range(d):
                yield (a, b, d)

def hnf_index(vs):
    """index in Z^2 of the span of integer vectors vs (0 if rank < 2)."""
    rows = [list(v) for v in vs if any(v)]
    while True:
        nz = [r for r in rows if r[0] != 0]
        if len(nz) <= 1:
            break
        nz.sort(key=lambda r: abs(r[0]))
        p = nz[0]
        for r in nz[1:]:
            q = r[0] // p[0]
            r[0] -= q * p[0]
            r[1] -= q * p[1]
        rows = [r for r in rows if any(r)]
    first = [r for r in rows if r[0] != 0]
    if not first:
        return 0
    g = 0
    for r in rows:
        if r[0] == 0:
            g = gcd(g, r[1])
    return abs(first[0][0] * g) if g else 0

def is_primitive(a, b, d):
    """Lambda Z[i] = Z[i]: the Z-span of {d, di, b+ai, i(b+ai)} is Z[i]."""
    return hnf_index([(d, 0), (0, d), (b, a), (-a, b)]) == 1

def in_lattice(v, a, b, d):
    x, y = v
    if y % a:
        return False
    return (x - (y // a) * b) % d == 0

def conductor(a, b, d, n):
    """smallest y | n with (x + y i) Lambda in Lambda for some x."""
    c, dd = (d, 0), (b, a)
    for y in range(1, n + 1):
        if n % y:
            continue
        for x in range(n):
            if in_lattice(gmul((x, y), c), a, b, d) and \
               in_lattice(gmul((x, y), dd), a, b, d):
                return y
    return None

def i_times(a, b, d, n):
    """Hermite data of i * Lambda (generators (-a, b), (0, d))."""
    (x1, y1), (x2, y2) = (-a, b), (0, d)
    aa = gcd(y1, y2)
    g, u, v = ext_gcd(y1, y2)
    b0 = u * x1 + v * x2
    dre = abs((-y2 // g) * x1 + (y1 // g) * x2)
    return (aa, b0 % dre, dre)

def ext_gcd(p, q):
    if q == 0:
        return (abs(p), 1 if p >= 0 else -1, 0)
    g, u, v = ext_gcd(q, p % q)
    return g, v, u - (p // q) * v

def build_X(a, b, d):
    """X in SL_2(Z[i]) with bottom row (c, dd) = (d, b + a i);
    Im(conj(c) dd) = a d = n > 0, so X(H) is the bounded disk."""
    c, dd = (d, 0), (b, a)
    g, s, t = gxgcd(dd, c)
    assert g in UNITS, (a, b, d)
    ginv = {(1, 0): (1, 0), (-1, 0): (-1, 0),
            (0, 1): (0, -1), (0, -1): (0, 1)}[g]
    A, B = gmul(ginv, s), gneg(gmul(ginv, t))
    assert gsub(gmul(A, dd), gmul(B, c)) == (1, 0)
    return (A, B), (c, dd)

def zeta_of(X):
    (A, B), (c, dd) = X
    return gmul((0, 1), gsub(gmul(A, gconj(dd)), gmul(B, gconj(c))))

def classification_classes(n):
    """zeta mod 2n with x even, y odd, x^2 + y^2 = 1 mod 4n
    (circle-classification.md)."""
    return {(x, y) for x in range(0, 2 * n, 2) for y in range(1, 2 * n, 2)
            if (x * x + y * y) % (4 * n) == 1}

# ==================== binary quadratic forms, disc -4n^2 ====================
def reduced_forms(D, primitive=True):
    out = []
    b0 = 1
    while 3 * b0 * b0 <= -D:
        b0 += 1
    for A in range(1, b0 + 1):
        for B in range(-A + 1, A + 1):
            num = B * B - D
            if num % (4 * A):
                continue
            C = num // (4 * A)
            if C < A or (A == C and B < 0):
                continue
            if primitive and gcd(gcd(A, B), C) != 1:
                continue
            out.append((A, B, C))
    return out

def reduce_form(A, B, C):
    while True:
        Bm = B % (2 * A)
        if Bm > A:
            Bm -= 2 * A
        C += (Bm * Bm - B * B) // (4 * A)
        B = Bm
        if A > C or (A == C and B < 0):
            A, B, C = C, -B, A
            continue
        if B == -A:
            B = A
            continue
        return (A, B, C)

def form_of_lattice(a, b, d):
    """class of Lambda in Cl(O_n): reduced form of (d^2, -2bd, a^2+b^2)."""
    return reduce_form(d * d, -2 * b * d, a * a + b * b)

# ==================== modular forms (theta constants) ====================
def E4E6D(tau):
    q = exp(mpc(0, 1) * pi * tau)
    t2, t3, t4 = jtheta(2, 0, q), jtheta(3, 0, q), jtheta(4, 0, q)
    a, b, c = t2 ** 4, t3 ** 4, t4 ** 4
    return ((a*a + b*b + c*c) / 2, (b + c) * (a + b) * (c - a) / 2,
            (t2 * t3 * t4) ** 8 / 256)

def sl2_reduce(tau):
    """(tau_red, (c, d)): tau_red = gamma tau in the fundamental domain,
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

def J_at(tau):
    tr, _ = sl2_reduce(tau)
    E4, _, Dl = E4E6D(tr)
    return E4 ** 3 / Dl

def Jp_at(tau):
    """j'(tau), computed at the reduced point (weight 2)."""
    tr, (c, d) = sl2_reduce(tau)
    E4, E6, Dl = E4E6D(tr)
    return (-2 * pi * mpc(0, 1) * E4 * E4 * E6 / Dl) / (c * tau + d) ** 2

def Dq_at(tau):
    """Delta_q = eta^24 (weight 12), at the reduced point."""
    tr, (c, d) = sl2_reduce(tau)
    return E4E6D(tr)[2] / (c * tau + d) ** 12

def Omega_period():
    return gamma(mpf(1) / 4) ** 4 / (8 * pi ** 2)

def cval(z):
    return mpc(z[0], z[1])

def theta_E(Xc):
    c, d = Xc[1][0], Xc[1][1]
    w = d / c
    assert w.imag > 0
    return Jp_at(w) / c ** 2

# ==================== certified recognition ====================
def cert_integer(x):
    """(int, spare-digits) or (None, reason).

    The criterion is the ABSOLUTE error |x - nint(x)|: a value of digit-size
    S computed at dps carries an error floor ~ 10^(S - dps), so a genuine
    integer certifies with spare ~ dps - S digits, while a huge non-integer
    (fractional part O(1)) is rejected no matter how large it is.  (A
    relative-error criterion silently accepts every number of size > dps/5 -
    the same class of trap as moduli-invariants.md section 5.6.)"""
    if fabs(x.imag) > mpf(10) ** (-mp.dps // 3) * (1 + fabs(x)):
        return None, "not real"
    xi = nint(x.real)
    err = fabs(x.real - xi)
    if err == 0:
        return int(xi), mp.inf
    spare = float(-log10(err))
    if spare > max(20, mp.dps / 5.0):
        return int(xi), spare
    return None, f"absolute error only 1e-{spare:.0f}"

def cert_rational(x, maxden=10 ** 60):
    if fabs(x.imag) > mpf(10) ** (-mp.dps // 3) * (1 + fabs(x)):
        return None, "not real"
    a = x.real
    p0, q0, p1, q1 = 1, 0, 0, 1
    t = a
    for _ in range(4 * mp.dps):
        cfl = int(floor(t))
        p0, p1 = cfl * p0 + p1, p0
        q0, q1 = cfl * q0 + q1, q0
        if q0 == 0:
            break
        err = fabs(a - mpf(p0) / mpf(q0))
        if err == 0:
            return Fraction(p0, q0), mp.inf
        acc = float(-log10(err / (1 + fabs(a))))
        info = len(str(abs(p0))) + len(str(abs(q0)))
        if acc >= mp.dps - max(10, mp.dps // 15):
            if acc - info < mp.dps / 5:
                return None, f"unsafe fit ({acc - info:.0f} spare)"
            return Fraction(p0, q0), acc - info
        frac = t - cfl
        if frac == 0:
            return Fraction(p0, q0), mp.inf
        if q0 > maxden:
            return None, "denominator blowup"
        t = 1 / frac
    return None, "no fit"

def poly_from_roots(roots):
    coeffs = [mpc(1)]
    for r in roots:
        new = [mpc(0)] * (len(coeffs) + 1)
        for k, cf in enumerate(coeffs):
            new[k] += cf
            new[k + 1] -= r * cf
        coeffs = new
    return coeffs

# ==================== per-level data ====================
def class_reps(n):
    """one lattice per ring class (partner is i*Lambda)."""
    seen = {}
    for L in sublattices(n):
        if not is_primitive(*L):
            continue
        f = form_of_lattice(*L)
        if f not in seen:
            seen[f] = L
    return seen

def class_values(n):
    """form -> (u, j, L) at current precision."""
    Om = Omega_period()
    out = {}
    for f, L in class_reps(n).items():
        (A, B), (c, dd) = build_X(*L)
        w = cval(dd) / cval(c)
        out[f] = (Jp_at(w) / cval(c) ** 2 / Om, J_at(w), L)
    return out

DPS = {2: 140, 3: 180, 4: 220, 5: 260, 6: 320, 7: 360, 8: 420, 9: 500,
       10: 500, 11: 540, 12: 600, 13: 600, 14: 680, 15: 760, 16: 840}

# ==================== A. exact structure ====================
def experiment_A():
    print("=" * 76)
    print("A. structure: disks <-> primitive sublattices <-> ring classes")
    print("=" * 76)
    print(" n   N_e  #prim  2h(-4n^2)  bijection  classmap  conductor  i-pair")
    for n in range(1, MAXN_STRUCT + 1):
        classes = classification_classes(n)
        prims = [L for L in sublattices(n) if is_primitive(*L)]
        zetas = []
        for L in prims:
            z = zeta_of(build_X(*L))
            zetas.append((z[0] % (2 * n), z[1] % (2 * n)))
        bij = sorted(zetas) == sorted(classes) and len(set(zetas)) == len(zetas)
        forms = reduced_forms(-4 * n * n)
        from collections import Counter
        cnt = Counter(form_of_lattice(*L) for L in prims)
        want = 1 if n == 1 else 2
        cm = set(cnt) == set(forms) and all(v == want for v in cnt.values())
        co = all(conductor(*L, n) == n for L in prims)
        ip = True
        if n > 1:
            for L in prims:
                Li = i_times(*L, n)
                ip &= is_primitive(*Li) and Li != L and \
                    form_of_lattice(*Li) == form_of_lattice(*L)
        tick = lambda b: "OK " if b else "FAIL"
        print(f"{n:3d} {len(classes):5d} {len(prims):5d} "
              f"{(2 if n > 1 else 1) * len(forms):8d}"
              f"      {tick(bij)}      {tick(cm)}      {tick(co)}      {tick(ip)}")
        assert bij and cm and co and ip, n
    print("(bijection: multiset of zeta mod 2n classes = classification classes;")
    print(" classmap: every class of disc -4n^2 hit exactly twice (once, n=1);")
    print(" conductor: O(Lambda) = Z + nZ[i] exactly; i-pair: partner = i*Lambda)")

# ==================== B. invariance, fiber, period ====================
def experiment_B():
    mp.dps = 120
    print("=" * 76)
    print("B. invariance of Theta, the fiber, the period")
    print("=" * 76)
    import random
    rng = random.Random(11)

    def matc(P, Q):
        return [[P[0][0]*Q[0][0] + P[0][1]*Q[1][0],
                 P[0][0]*Q[0][1] + P[0][1]*Q[1][1]],
                [P[1][0]*Q[0][0] + P[1][1]*Q[1][0],
                 P[1][0]*Q[0][1] + P[1][1]*Q[1][1]]]

    def Xc_of(X):
        (A, B), (c, d) = X
        return [[cval(A), cval(B)], [cval(c), cval(d)]]

    TOL = mpf(10) ** (-90)
    for n in (5, 6, 7):
        for L in list(class_reps(n).values())[:2]:
            Xc = Xc_of(build_X(*L))
            th0 = theta_E(Xc)
            for _ in range(4):
                T = [[mpf(1), mpc(rng.randint(-9, 9), rng.randint(-9, 9))],
                     [mpf(0), mpf(1)]]
                G = [[mpf(1), mpf(0)], [mpf(0), mpf(1)]]
                for _ in range(6):
                    g = rng.choice([[[1, rng.randint(-3, 3)], [0, 1]],
                                    [[0, -1], [1, 0]]])
                    G = matc(G, [[mpf(g[0][0]), mpf(g[0][1])],
                                 [mpf(g[1][0]), mpf(g[1][1])]])
                th1 = theta_E(matc(T, matc(Xc, G)))
                assert fabs(th1 - th0) < TOL * (1 + fabs(th0)), (n, L)
            c, d = Xc[1][0], Xc[1][1]
            w = d / c
            assert fabs(fabs(th0) - w.imag * fabs(Jp_at(w)) / n) \
                < TOL * (1 + fabs(th0))
    print("two-sided invariance (T_alpha X gamma, random) and")
    print("|Theta| = Im(w)|j'(w)|/n : OK to 90+ digits, n = 5, 6, 7")

    n = 5
    L = list(class_reps(n).values())[1]
    Xc = Xc_of(build_X(*L))
    w = Xc[1][1] / Xc[1][0]
    p = -w.conjugate()      # the fiber rotates about the pole -w
    x, y = p.real, p.imag
    K = matrix([[-x / y, (x*x + y*y) / y], [-1 / y, x / y]])
    print()
    print("fiber = real rotations fixing the pole X^{-1}(inf):")
    for t in (mpf(0), mpf(1) / 5, mpf(2) / 5):
        H = expm(K * t)
        Y = matc(Xc, [[H[0, 0], H[0, 1]], [H[1, 0], H[1, 1]]])
        th = theta_E(Y)
        print(f"  t = {nstr(t, 3):>5}:  |Theta| = {nstr(fabs(th), 25)}   "
              f"arg Theta = {nstr(arg(th), 25)}")
    print("  (|Theta| constant to all digits; arg moves linearly at rate 2)")

    print()
    Om6 = ((2 * pi) ** 6 * Dq_at(mpc(0, 1))).real
    Om = Omega_period()
    pom = gamma(mpf(1) / 4) ** 2 / (2 * msqrt(2 * pi))
    # pomega = 2 int_0^1 dt/sqrt(1-t^4) = B(1/4,1/2)/2 (substitute s = t^4)
    lem_exact = gamma(mpf(1) / 4) * gamma(mpf(1) / 2) / (2 * gamma(mpf(3) / 4))
    lem_quad = 2 * quad(lambda s: 1 / msqrt(1 - s ** 4), [0, 1])
    assert fabs(Om ** 6 - Om6) < mpf(10) ** (-100) * Om6
    assert fabs(pom ** 2 / pi - Om) < mpf(10) ** (-100)
    assert fabs(lem_exact - pom) < mpf(10) ** (-100)
    assert fabs(lem_quad - pom) < mpf(10) ** (-25)
    print(f"Omega = Gamma(1/4)^4/(8 pi^2) = pomega^2/pi = ((2pi)^6 eta(i)^24)^(1/6)")
    print(f"      = {nstr(Om, 40)}   (pomega = {nstr(pom, 30)}, lemniscate)")

    # u^6 closed form on a few disks
    for n in (3, 5, 8):
        for f, (u, jw, L) in list(class_values(n).items())[:2]:
            (A, B), (c, dd) = build_X(*L)
            w = cval(dd) / cval(c)
            rhs = -jw ** 4 * (jw - 1728) ** 3 * Dq_at(w) / \
                (cval(c) ** 12 * Dq_at(mpc(0, 1)))
            assert fabs(u ** 6 - rhs) < mpf(10) ** (-80) * (1 + fabs(rhs))
    print("u^6 = -j(w)^4 (j(w)-1728)^3 Delta(Lambda)/Delta(Z[i]) : OK (n=3,5,8)")

# ==================== C. ring class polynomials, traces ====================
def hurwitz_trace(D):
    """Zagier's t(|D|): all reduced forms of disc D (imprimitive included),
    weights 1/2, 1/3 at the forms equivalent to a(x^2+y^2), a(x^2+xy+y^2)."""
    tot = mpf(0)
    for (A, B, C) in reduced_forms(D, primitive=False):
        g = gcd(gcd(A, B), C)
        f0 = (A // g, B // g, C // g)
        w = 2 if f0 == (1, 0, 1) else (3 if f0 == (1, 1, 1) else 1)
        tau = mpc(mpf(-B) / (2 * A), msqrt(mpf(-D)) / (2 * A))
        tot += (J_at(tau).real - 744) / w
    return tot

def experiment_C():
    print("=" * 76)
    print("C. singular moduli: prod over disks (x - j) = H_{-4n^2}^2; traces")
    print("=" * 76)
    Hs = {}
    for n in range(2, MAXN + 1):
        mp.dps = max(160, DPS.get(n, 400) // 2)
        cv = class_values(n)
        js = [cv[f][1] for f in sorted(cv)]
        Hc = poly_from_roots(js)
        Hint, sp_min = [], mp.inf
        for cf in Hc:
            v, sp = cert_integer(cf)
            assert v is not None, (n, sp)
            Hint.append(v)
            sp_min = min(sp_min, sp)
        Hs[n] = Hint
        # disk side: all 2h disks
        js2 = []
        for L in sublattices(n):
            if not is_primitive(*L):
                continue
            (A, B), (c, dd) = build_X(*L)
            js2.append(J_at(cval(dd) / cval(c)))
        P2h = poly_from_roots(js2)
        # H^2 exactly
        sq = [0] * (2 * len(Hint) - 1)
        for i1, c1 in enumerate(Hint):
            for i2, c2 in enumerate(Hint):
                sq[i1 + i2] += c1 * c2
        ok = all(cert_integer(cf)[0] == sq[k] for k, cf in enumerate(P2h))
        shown = Hint if len(str(Hint)) < 100 else \
            f"[deg {len(Hint)-1}, coeffs up to {max(len(str(abs(x))) for x in Hint)} digits]"
        print(f"n={n:2d}: H_{{-4n^2}} = {shown}")
        print(f"       integer-certified (>= {nstr(mpf(sp_min), 4)} spare digits); "
              f"disk product = H^2: {'OK' if ok else 'FAIL'}")
        assert ok
    print()
    print("traces: Tr_E(n) := sum over Cl(O_n) of (j - 744);  slicing check")
    print("        sum_{g | n} Tr_E(n/g) = t(4 n^2)   (Zagier's trace)")
    mp.dps = 100
    TrE = {1: mpf(492)}
    for n in range(2, MAXN + 1):
        h = len(Hs[n]) - 1
        TrE[n] = mpf(-Hs[n][1] - 744 * h)   # sum of roots - 744h
    for n in range(1, MAXN + 1):
        divs = [g for g in range(1, n + 1) if n % g == 0]
        t = sum(TrE[n // g] for g in divs)
        tz = hurwitz_trace(-4 * n * n)
        ok = fabs(t - tz) < mpf(10) ** (-40) * (1 + fabs(tz))
        print(f"  n={n:2d}: Tr_E = {int(nint(TrE[n])):>28d}   slice = t(4n^2): "
              f"{'OK' if ok else 'FAIL'}")
        assert ok

# ==================== D. the phase units ====================
def exact_square_poly(us, scale=1):
    P = poly_from_roots([(mpf(scale) * u) ** 2 for u in us])
    out = []
    for cf in P:
        v, sp = cert_integer(cf)
        if v is not None:
            out.append((Fraction(v), sp))
            continue
        vr, spr = cert_rational(cf)
        if vr is None:
            return None, spr
        out.append((vr, spr))
    return out, min(sp for _, sp in out)

def min_lambda(frs):
    import sympy
    L = 1
    for k, cf in enumerate(frs):
        den = cf.denominator
        if den == 1:
            continue
        for p, e in sympy.factorint(den).items():
            need = -(-e // (2 * k))
            have = 0
            m = L
            while m % p == 0:
                have += 1
                m //= p
            if have < need:
                L *= p ** (need - have)
    return L

def experiment_D(levels=None):
    import sympy
    print("=" * 76)
    print("D. the phase units u = Theta/Omega on the Bianchi group")
    print("=" * 76)
    levels = levels or list(range(2, MAXN + 1))
    rows = []
    for n in levels:
        mp.dps = DPS.get(n, 600)
        cv = class_values(n)
        fs = sorted(cv)
        us = [cv[f][0] for f in fs]
        js = [cv[f][1] for f in fs]
        h = len(fs)

        # laws: u(iL) = -u(L); mirror u(conj disk) = conj u  (lattice form:
        # u(conj Lambda) = -conj u(Lambda)); reality on ambiguous classes
        for f in fs[:3]:
            (a, b, d) = cv[f][2]
            u0 = cv[f][0]
            tol = mpf(10) ** (-mp.dps // 2) * (1 + fabs(u0))
            Li = i_times(a, b, d, n)
            Om = Omega_period()
            (A, B), (c, dd) = build_X(*Li)
            ui = Jp_at(cval(dd) / cval(c)) / cval(c) ** 2 / Om
            assert fabs(ui + u0) < tol, (n, f)
            Lc = (a, (-b) % d, d)
            (A, B), (c, dd) = build_X(*Lc)
            uc = Jp_at(cval(dd) / cval(c)) / cval(c) ** 2 / Om
            assert fabs(uc + u0.conjugate()) < tol
        # reality pattern + criterion (real <=> n | y(zeta); imag <=> n | x)
        pat = []
        for f in fs:
            u, _, L = cv[f]
            zs = []
            for LL in sublattices(n):
                if is_primitive(*LL) and form_of_lattice(*LL) == f:
                    z = zeta_of(build_X(*LL))
                    zs.append((z[0] % (2 * n), z[1] % (2 * n)))
            t = fabs(u)
            isR = fabs(u.imag) < mpf(10) ** (-mp.dps // 3) * t
            isI = fabs(u.real) < mpf(10) ** (-mp.dps // 3) * t
            assert isR == any(zy % n == 0 for _, zy in zs), (n, f)
            assert isI == any(zx % n == 0 for zx, _ in zs), (n, f)
            pat.append('R' if isR else ('I' if isI else 'C'))

        # P2 with minimal lambda
        fr0, sp0 = exact_square_poly(us)
        assert fr0 is not None, (n, sp0)
        lam = min_lambda([cf for cf, _ in fr0])
        fr, spmin = exact_square_poly(us, scale=lam)
        coeffs = [int(cf) for cf, _ in fr]
        assert all(cf.denominator == 1 for cf, _ in fr)
        # P6 exactly from P2: roots (lam v)^6 = ((lam v)^2)^3 * lam^4-correction
        # (lam v)^6 = ((lam v)^2)^3 / lam^0 ... use resultant on scaled poly
        y, z = sympy.symbols('y z')
        P2 = sum(sympy.Integer(cf) * z ** (h - k) for k, cf in enumerate(coeffs))
        P6 = sympy.Poly(sympy.resultant(P2, y - z ** 3, z), y)
        if P6.LC() < 0:      # sympy's Res carries (-1)^(deg f * deg g)
            P6 = sympy.Poly(-P6.as_expr(), y)
        # numerical cross-check of P6: recompute the sixth powers at elevated
        # precision and certify each coefficient of prod(y - (lam v)^6)
        # against the exact resultant coefficient (absolute agreement).
        c6 = [int(cc) for cc in P6.all_coeffs()]
        size6 = max(len(str(abs(cc))) for cc in c6)
        mp.dps = size6 + 200
        cv6 = class_values(n)
        P6n = poly_from_roots([(mpf(lam) * cv6[f][0]) ** 6 for f in fs])
        ok6 = all(fabs(cf - c6[k]) < mpf(10) ** (-100) for k, cf in enumerate(P6n))
        assert ok6, n
        mp.dps = DPS.get(n, 600)
        # factorization / Galois
        poly2 = sympy.Poly(P2, z)
        fl = sympy.factor_list(poly2.as_expr())[1]
        degs = sorted(sympy.Poly(g, z).degree() for g, m in fl for _ in range(m))
        gal = "-"
        if h <= 6 and degs == [h]:
            try:
                from sympy.polys.numberfields.galoisgroups import galois_group
                G = galois_group(poly2)[0]
                names = {(2, 2): "Z/2", (4, 8): "D4", (4, 4): "V4/C4",
                         (6, 12): "D6", (6, 6): "C6", (3, 6): "S3"}
                gal = names.get((h, G.order()), f"order {G.order()}")
            except Exception:
                gal = "?"
        elif degs == [h]:
            disc = sympy.discriminant(poly2)
            sq = sympy.sqrt(disc).is_integer if disc > 0 else False
            gal = f"deg {h}; disc {'is' if sq else 'not'} a square"
        spf = float(spmin)
        rows.append((n, h, lam, pat, coeffs, degs, gal, spf, ok6))
        pretty = coeffs if h <= 2 else f"[degree {h}]"
        print(f"n={n:2d} (h={h}, lambda_n={lam}): pattern {''.join(pat)},  "
              f"P2 in Z[y]: {pretty}")
        print(f"       spare digits >= {spf:.0f}; P6 (exact cube-resultant) "
              f"numerically confirmed: {'OK' if ok6 else 'FAIL'};  "
              f"factors {degs};  Galois {gal}")
    print()
    print("summary table:")
    print("  n   h  lambda_n  reality   irreducible  Galois")
    for (n, h, lam, pat, coeffs, degs, gal, spm, ok6) in rows:
        print(f"{n:3d} {h:3d}  {lam:7d}  {''.join(pat):9s} "
              f"{'yes' if degs == [h] else str(degs):>11s}  {gal}")
    return rows

def experiment_D_tagging(levels=None):
    import sympy
    print("=" * 76)
    print("D'. prime support of P2(0), and the Delta-mass")
    print("=" * 76)
    levels = levels or list(range(2, MAXN + 1))
    for n in levels:
        mp.dps = DPS.get(n, 600)
        cv = class_values(n)
        fs = sorted(cv)
        us, js = [cv[f][0] for f in fs], [cv[f][1] for f in fs]
        h = len(fs)
        fr, _ = exact_square_poly(us, scale=n * n)
        P20 = int(fr[-1][0])
        Hc = poly_from_roots(js)
        Hint = [cert_integer(cf)[0] for cf in Hc]
        H0 = Hint[-1]
        H1728 = sum(Hint[k] * 1728 ** (h - k) for k in range(h + 1))
        prodb = Fraction((-1) ** h * H0)
        prodb17 = Fraction((-1) ** h * H1728)
        prodv6 = Fraction((-1) ** h * P20, n ** (4 * h)) ** 3
        mass = prodv6 / (Fraction((-1) ** h) * prodb ** 4 * prodb17 ** 3)
        Mn = mass * Fraction(n ** (12 * h))     # prod_c n^12 R_c
        fac = lambda m: ('-' if m < 0 else '') + \
            ('1' if abs(m) == 1 else str(sympy.factorint(abs(m))))
        print(f"n={n:2d}: P2(0) {fac(P20)}")
        print(f"       H(0) {fac(H0)}   H(1728) {fac(H1728)}")
        assert Mn.denominator == 1
        print(f"       mass prod_c n^12 Delta(L_c)/Delta(Z[i]) = {fac(int(Mn))}")
        # support check: primes of P2(0) inside n, H(0), H(1728)
        sup = set(sympy.factorint(abs(P20)))
        cover = set(sympy.factorint(abs(H0 * H1728 * n * 2)))
        assert sup <= cover, (n, sup - cover)
    print("(support of P2(0) always inside primes of 2n, H(0), H(1728))")

def experiment_mass(levels=None):
    """direct verification of the Delta-mass law, independent of P2:
       M(n) = prod_c n^12 Delta(Lambda_c)/Delta(Z[i])
            = +- prod_{p^k || n, p ne 1 mod 4} p^( (6/e_p)(p^k-1)/(p-1) N_e(n/p^k) )
       with e_2 = 2, e_p = 1 for p = 3 mod 4; split primes contribute 1;
       observed sign: -1 exactly at n = 4, 8, 16 (n = 2^k, k >= 2)."""
    import sympy
    from mpmath import log as mlog
    print("=" * 76)
    print("D''. the Delta-mass, directly")
    print("=" * 76)
    mp.dps = 120
    levels = levels or (list(range(2, 17)) + [18, 21, 25, 27, 49])
    for n in levels:
        Di = Dq_at(mpc(0, 1))
        tot = mpc(1)
        reps = class_reps(n)
        for f, L in reps.items():
            (A, B), (c, dd) = build_X(*L)
            cc = cval(c)
            tot *= mpf(n) ** 12 * Dq_at(cval(dd) / cc) / (cc ** 12 * Di)
        pred = {}
        for p, k in sympy.factorint(n).items():
            rest = n // p ** k
            Ne = rest
            for q in sympy.factorint(rest):
                if q % 2 == 1:
                    Ne = Ne * (q - (1 if q % 4 == 1 else -1)) // q
            if p == 2:
                pred[p] = 3 * (2 ** k - 1) * Ne
            elif p % 4 == 3:
                pred[p] = 6 * (p ** k - 1) // (p - 1) * Ne
            else:
                pred[p] = 0
        predlog = sum(e * mlog(mpf(p)) for p, e in pred.items())
        ok = fabs(mlog(fabs(tot)) - predlog) < mpf(10) ** (-60) * (1 + predlog)
        assert fabs(tot.imag) < mpf(10) ** (-60) * fabs(tot)
        sgn = '+' if tot.real > 0 else '-'
        pretty = " * ".join(f"{p}^{e}" for p, e in sorted(pred.items()) if e) or "1"
        print(f"  n={n:3d} (h={len(reps):3d}):  M = {sgn}{pretty}   "
              f"{'OK' if ok else 'FAIL'}")
        assert ok, n
    print("(law: v_p(M) = (6/e_p) (p^k-1)/(p-1) N_e(n/p^k), e_2 = 2,")
    print(" inert e_p = 1, split primes silent; sign -1 at n = 4, 8, 16)")

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "base"
    if mode in ("base", "all"):
        experiment_A()
        print()
        experiment_B()
        print()
        experiment_C()
    if mode in ("phase", "all"):
        levels = [int(a) for a in sys.argv[2:]] or None
        experiment_D(levels)
        print()
        experiment_D_tagging(levels)
        print()
        experiment_mass()
