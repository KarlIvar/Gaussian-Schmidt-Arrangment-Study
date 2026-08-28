#!/usr/bin/env python3
"""orthospectrum_verify.py — machine checks for spectral-geometry/orthogeodesics.md

Setting.  P_R = the totally geodesic plane of H^3 above the real line
({(x, 0, t)}), whose Gamma = PSL_2(Z[i])-orbit has the Schmidt circles as
boundary shadows.  For a level-n circle omega (curvature 2q, centre
(x + ni)/(2q), x^2 + n^2 - 1 = 4qm, N = n^2 - 1) the checks certify:

 1. [exact]  With X in SL_2(Z[i]) realising omega = X(Rhat) (built by the
    descent of circle-classification.md) and Z = X conj(X)^{-1} the Cartan
    image (= the composition of the inversions in omega and in Rhat):
      - tr Z = +-2n  (the trace identity of involution.md);
      - Z (m1, 1)^T = lam (m1, 1)^T  with  m1 = (x + i sqrt(N))/(2q)  the
        hyperbolic centre / CM point, and  lam = +-(n +- sqrt(N))  a unit
        of Z[sqrt(N)] — verified in exact arithmetic in the biquadratic
        field Q(i, sqrt(N)).  Hence the axis of Z is the geodesic with
        ideal endpoints {m1, conj(m1)} and translation length
        2 log(n + sqrt(N)) = 2 arccosh(n).
 2. [exact + numeric]  Orthogeodesic geometry.  The geodesic with endpoints
    m1, conj(m1) meets P_R orthogonally at F1 = (x/(2q), 0, sqrt(N)/(2q))
    ("the CM point is the foot") and the hemisphere over omega orthogonally
    at F2 = (x/(2q), N/(2qn), sqrt(N)/(2qn)), and
        cosh d(F1, F2) = n   (exactly),
    matching the inversive product <M_0, M_omega> = alpha = n: the invariant
    alpha of hyperbolic-counting.md is cosh of the ortho-distance.
 3. [exact]  Counts.  The weighted number of level-n circles with hyperbolic
    centre in the ideal triangle equals 3 H(n^2-1) (hyperbolic-counting.md),
    i.e. H(n^2-1) per PSL_2(Z)-fundamental-domain: the ortho-length spectrum
    of the immersed modular surface is {arccosh n} with Hurwitz class number
    multiplicities.
 4. [numeric]  Growth.  S_H(X) = sum_{n <= X} H(n^2-1) grows like a constant
    times X^2, so the sigma-family of closed geodesics (lengths
    2 arccosh n = 2 log eps_n) has #{length <= L} ~ e^L: exactly HALF the
    topological entropy 2 of the geodesic flow on the Bianchi orbifold.
    (A bulk Hurwitz sieve is cross-checked against the hurwitz() routine of
    alpha_circles.py.)
 5. [exact + numeric]  The atom/gasket stratum.  The atoms of the half-plane
    monoid (half-plane-monoid.md) are the disks of the strip Apollonian
    gasket; the per-period count N(K) of atoms of curvature <= K is computed
    by exact Descartes recursion on Hermitian triples and its growth exponent
    is fitted: N(K) ~ K^delta with delta = 1.30568... (the Hausdorff
    dimension of the gasket = Patterson-Sullivan exponent of the infinite
    covolume Apollonian group, base eigenvalue lambda_0 = delta(2-delta)).

Usage:  python3 scripts/orthospectrum_verify.py          (a few seconds)
"""
import math
import sys
import os
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import mpmath as mp

from omega import G, Mat, Disk, disk_of, matrix_for_disk, inversive, H as HPLANE
from alpha_circles import alpha_circles, hurwitz

mp.mp.dps = 40
PASS, FAIL = "PASS", "FAIL"
results = []


def report(name, ok, detail=""):
    results.append(ok)
    print(f"[{PASS if ok else FAIL}] {name}" + (f"  {detail}" if detail else ""))


# ------------------------------------------------------------ exact arithmetic
class CF:
    """Complex rational a + b i, a, b in Q."""
    __slots__ = ("re", "im")

    def __init__(self, re=0, im=0):
        self.re, self.im = Fraction(re), Fraction(im)

    @staticmethod
    def of_g(g):
        return CF(g.re, g.im)

    def __add__(s, o): return CF(s.re + o.re, s.im + o.im)
    def __sub__(s, o): return CF(s.re - o.re, s.im - o.im)

    def __mul__(s, o):
        return CF(s.re * o.re - s.im * o.im, s.re * o.im + s.im * o.re)

    def __eq__(s, o): return s.re == o.re and s.im == o.im
    def __repr__(s): return f"({s.re}+{s.im}i)"

    def is_zero(s): return s.re == 0 and s.im == 0


class BQ:
    """Element u + v sqrt(N) of the biquadratic field Q(i, sqrt(N)),
    u, v in Q(i);  N > 0 a fixed non-square."""
    __slots__ = ("u", "v", "N")

    def __init__(self, u, v, N):
        self.u, self.v, self.N = u, v, N

    def __add__(s, o): return BQ(s.u + o.u, s.v + o.v, s.N)
    def __sub__(s, o): return BQ(s.u - o.u, s.v - o.v, s.N)

    def __mul__(s, o):
        NN = CF(s.N)
        return BQ(s.u * o.u + NN * (s.v * o.v), s.u * o.v + s.v * o.u, s.N)

    def __eq__(s, o): return s.u == o.u and s.v == o.v
    def is_zero(s): return s.u.is_zero() and s.v.is_zero()
    def __repr__(s): return f"{s.u} + {s.v} sqrt({s.N})"


def bq_of_g(g, N):
    return BQ(CF.of_g(g), CF(0), N)


# ------------------------------------------------------------ 1. exact axis
def conj_mat(X):
    return Mat(X.a.conj(), X.b.conj(), X.c.conj(), X.d.conj())


def check_axis_exact():
    total = bad = 0
    lam_kinds = set()
    for n in (3, 5, 7, 9, 11, 13, 15):
        N = n * n - 1
        for c in alpha_circles(n):
            q, x, m = c['q'], c['x'], c['m']
            D = Disk(-2 * q, G(x, n), -2 * m)
            assert D.det() == -1
            X = matrix_for_disk(D)
            assert disk_of(X) == D
            Z = X * conj_mat(X).inv()
            tr = Z.a + Z.d
            total += 1
            if not (tr.im == 0 and abs(tr.re) == 2 * n):
                bad += 1
                continue
            # exact eigenvector (m1, 1),  m1 = x/(2q) + (i/(2q)) sqrt(N)
            m1 = BQ(CF(Fraction(x, 2 * q)), CF(0, Fraction(1, 2 * q)), N)
            one = BQ(CF(1), CF(0), N)
            top = bq_of_g(Z.a, N) * m1 + bq_of_g(Z.b, N) * one
            lam = bq_of_g(Z.c, N) * m1 + bq_of_g(Z.d, N) * one
            # eigen condition: top = lam * m1
            if not (top - lam * m1).is_zero():
                bad += 1
                continue
            # lam = +-n +- sqrt(N):  u-part rational +-n, v-part +-1
            u, v = lam.u, lam.v
            if not (u.im == 0 and v.im == 0 and abs(u.re) == n and abs(v.re) == 1):
                bad += 1
                continue
            lam_kinds.add((u.re > 0, v.re > 0, (u.re > 0) == (v.re > 0)))
    ok = bad == 0
    eps_units = sorted({k[2] for k in lam_kinds})
    report("1. exact axis: tr(X conj(X)^-1) = +-2n and (m1,1) is an eigenvector "
           "with eigenvalue +-(n +- sqrt(n^2-1))", ok,
           f"{total} circles across odd n <= 15, all exact"
           + (f", |multiplier| = eps^2 in every case" if ok else f", {bad} failures"))
    # translation length: |lam/lam'| = eps^2, i.e. length 2 log eps = 2 arccosh n
    n = 7
    eps = n + math.sqrt(48)
    report("   translation length 2 log(n + sqrt(n^2-1)) = 2 arccosh n",
           abs(2 * math.log(eps) - 2 * math.acosh(n)) < 1e-14,
           f"(n = 7: {2*math.log(eps):.12f})")


# ------------------------------------------------------------ 2. geometry
def check_orthogeodesic_geometry():
    bad_details = []
    ncase = 0
    for n in (3, 4, 5, 7, 9, 12, 15):
        N = n * n - 1
        for c in alpha_circles(n)[:6]:
            q, x = c['q'], c['x']
            ncase += 1
            qq, xx, nn, NN = (mp.mpf(q), mp.mpf(x), mp.mpf(n), mp.mpf(N))
            sN = mp.sqrt(NN)
            a = xx / (2 * qq)
            # geodesic with endpoints m1 = a + i sN/(2q), conj(m1):
            # semicircle {(a + iy, t): y^2 + t^2 = R^2},  R = sN/(2q)
            R = sN / (2 * qq)
            F1 = (a, mp.mpf(0), R)
            F2 = (a, NN / (2 * qq * nn), sN / (2 * qq * nn))
            z0y = nn / (2 * qq)          # centre of omega = a + i n/(2q)
            r = 1 / (2 * qq)
            # F2 on the hemisphere over omega?
            on = (F2[1] - z0y) ** 2 + F2[2] ** 2 - r ** 2
            # cosh distance between the feet:
            coshd = 1 + ((F2[1] - F1[1]) ** 2 + (F2[2] - F1[2]) ** 2) / (2 * F1[2] * F2[2])
            # orthogonality at F2: tangent of the semicircle (0, t, -y) parallel
            # to the Euclidean radial direction of the hemisphere (0, y - z0y, t)
            cross = F2[2] * F2[2] - (-F2[1]) * (F2[1] - z0y)   # t*t + y*(y-z0y)
            # parallel  <=>  (0,t,-y) x (0, y-z0y, t) = 0  <=>  t^2 + y(y-z0y) = 0? no:
            # cross product z-component: t*(y - z0y) - (-y)*t ... do it explicitly:
            v1 = (mp.mpf(0), F2[2], -F2[1])          # tangent to geodesic at F2
            v2 = (mp.mpf(0), F2[1] - z0y, F2[2])     # radial (normal) direction
            crossmag = abs(v1[1] * v2[2] - v1[2] * v2[1])
            n1 = mp.sqrt(v1[1] ** 2 + v1[2] ** 2)
            n2 = mp.sqrt(v2[1] ** 2 + v2[2] ** 2)
            sin_angle = crossmag / (n1 * n2)
            # inversive product with the plane over Rhat:
            D = Disk(-2 * q, G(x, n), -2 * (x * x + N) // (4 * q)) \
                if (x * x + N) % (4 * q) == 0 else None
            inv_ok = (D is not None and abs(inversive(HPLANE, D)) == n)
            checks = [abs(on) < mp.mpf('1e-30'),
                      abs(coshd - n) < mp.mpf('1e-30'),
                      sin_angle < mp.mpf('1e-30'),
                      inv_ok]
            if not all(checks):
                bad_details.append((n, q, x, checks))
    report("2. orthogeodesic geometry: foot F2 on the hemisphere, "
           "cosh d(F1, F2) = n exactly, geodesic orthogonal to the hemisphere, "
           "|<M_0, M_omega>| = n", not bad_details,
           f"{ncase} circles across n in {{3,4,5,7,9,12,15}} at 40 digits"
           + ("" if not bad_details else f"; failures: {bad_details[:3]}"))


# ------------------------------------------------------------ 3. counts
def check_counts():
    ok = True
    rows = []
    for n in range(2, 21):
        w = sum(c['weight'] for c in alpha_circles(n))
        h3 = 3 * hurwitz(n * n - 1)
        ok &= (w == h3)
        rows.append((n, w))
    report("3. weighted level-n count in the ideal triangle = 3 H(n^2-1), "
           "i.e. H(n^2-1) ortho-arcs per PSL_2(Z)-fundamental domain "
           "(n = 2..20)", ok)


# ------------------------------------------------------------ 4. bulk Hurwitz
def hurwitz_sieve6(Y):
    """6*H(N) for all 0 <= N <= Y, by enumerating (a, b, c), 0 <= b <= a <= c,
    N = 4ac - b^2, with the weights of alpha_circles.hurwitz."""
    H6 = np.zeros(Y + 1, dtype=np.int64)
    bmax = int(math.isqrt(Y // 3)) + 1
    for b in range(0, bmax + 1):
        a = max(b, 1)
        while 4 * a * a - b * b <= Y:
            cmax = (Y + b * b) // (4 * a)
            cs = np.arange(a, cmax + 1, dtype=np.int64)
            Ns = 4 * a * cs - b * b
            if b == 0:
                w = np.full(len(cs), 6, dtype=np.int64)
                w[0] = 3                          # (a, 0, a)
            elif b == a:
                w = np.full(len(cs), 6, dtype=np.int64)
                w[0] = 2                          # (a, a, a)
            else:
                w = np.full(len(cs), 12, dtype=np.int64)
                w[0] = 6                          # (a, b, a)
            np.add.at(H6, Ns, w)
            a += 1
    return H6


def check_hurwitz_growth():
    Y = 10 ** 6
    print("  (sieving Hurwitz class numbers up to 10^6 ...)")
    H6 = hurwitz_sieve6(Y)
    # cross-check against the exact routine
    rng = np.random.default_rng(20260828)
    ok = True
    for N in list(rng.integers(3, Y, 25)) + [n * n - 1 for n in range(2, 31)]:
        N = int(N)
        if N % 4 in (1, 2):
            ok &= (H6[N] == 0)
        else:
            ok &= (Fraction(int(H6[N]), 6) == hurwitz(N))
    report("4. bulk Hurwitz sieve H(N), N <= 10^6, agrees with "
           "alpha_circles.hurwitz (25 random N + all N = n^2-1, n <= 30)", ok)

    ns = np.arange(2, 1001)
    Hdiag = H6[ns * ns - 1] / 6.0
    S = np.cumsum(Hdiag)

    def S_at(X):
        return S[X - 2]

    vals = [(X, S_at(X) / X ** 2) for X in (125, 250, 500, 1000)]
    stable = abs(vals[-1][1] / vals[-2][1] - 1) < 0.05
    report("4. S_H(X) = sum_(n<=X) H(n^2-1) ~ C X^2 (multiplicities of the "
           "ortho-length spectrum)", stable,
           "S/X^2 = " + ", ".join(f"{v:.4f} (X={X})" for X, v in vals)
           + f";  naive random-discriminant heuristic pi/(12 zeta(3)) = "
             f"{math.pi/(12*float(mp.zeta(3))):.4f}")

    # entropy: closed sigma-geodesics have length l_n = 2 arccosh n;
    # #{l <= L} = S_H(cosh(L/2)) ~ (C/4) e^L : slope 1 = half of h_top = 2.
    slopes = []
    Ls = [10.0, 12.0, 14.0]
    for L1, L2 in zip(Ls, Ls[1:]):
        P1 = S_at(int(math.cosh(L1 / 2)))
        P2 = S_at(int(math.cosh(L2 / 2)))
        slopes.append((math.log(P2) - math.log(P1)) / (L2 - L1))
    ok = all(0.9 < sl < 1.1 for sl in slopes[-1:])
    report("4. entropy of the sigma-family: d log #/dL -> 1 "
           "(half the topological entropy 2 of the geodesic flow)", ok,
           f"slopes at L = 10->12->14: " + ", ".join(f"{sl:.3f}" for sl in slopes))


# ------------------------------------------------------------ 5. gasket atoms
ROOT = ((0, (0, -1), 0),        # {Im z < 0}
        (0, (0, 1), -2),        # {Im z > 1}
        (-2, (0, 1), 0),        # |z - i/2| < 1/2
        (-2, (2, 1), -2))       # |z - 1 - i/2| < 1/2


def d_add(D1, D2):
    return (D1[0] + D2[0], (D1[1][0] + D2[1][0], D1[1][1] + D2[1][1]), D1[2] + D2[2])


def d_scale(k, D):
    return (k * D[0], (k * D[1][0], k * D[1][1]), k * D[2])


def d_key(D):
    A, (bre, bim), C = D
    return (A, bre % (-A), bim)          # translation z -> z+1 shifts Re B by A


def translate_quad(quad, k):
    """Apply z -> z + k to a quadruple:  B_re += A * k, C adjusts too, but we
    only need a canonical *key*, so track (A, B, C) exactly:
    pull-back by T(k) sends (A, B, C) -> (A, B + A k, C + A k^2 + 2 k Re B)."""
    out = []
    for (A, (bre, bim), C) in quad:
        out.append((A, (bre + A * k, bim), C + A * k * k + 2 * k * bre))
    return tuple(out)


def quad_key(quad):
    """Canonical form of a quadruple modulo the translation z -> z + 1:
    normalise by the youngest (largest-curvature) bounded member."""
    j = min(range(4), key=lambda i: quad[i][0])          # most negative A
    A, (bre, _), _ = quad[j]
    a = -A
    k = bre // a                        # B += A k  brings bre into [0, a)
    return tuple(sorted(translate_quad(quad, k)))


def gasket_count(Kmax, sample_cap=60):
    """All bounded strip-gasket disks with curvature |A| <= Kmax, one per
    period, by exact Descartes recursion  D_j' = 2(D_k + D_l + D_m) - D_j.
    Quadruples are deduplicated modulo the translation z -> z+1 (otherwise
    the walk drifts along the strip forever).  Returns (sorted curvature
    list, sample of small disks as (A, B, C) triples)."""
    seen = {}
    samples = []
    for D in ROOT[2:]:
        k = d_key(D)
        if k not in seen:
            seen[k] = -D[0]
            samples.append(D)
    seen_quads = {quad_key(ROOT)}
    stack = [(ROOT, -1)]
    while stack:
        quad, last = stack.pop()
        tot = quad[0]
        for j in range(1, 4):
            tot = d_add(tot, quad[j])
        for j in range(4):
            if j == last:
                continue
            newD = d_add(d_scale(2, d_add(tot, d_scale(-1, quad[j]))), d_scale(-1, quad[j]))
            A = newD[0]
            if A >= 0 or -A > Kmax:
                continue
            k = d_key(newD)
            if k not in seen:
                seen[k] = -A
                if -A <= 60 and len(samples) < sample_cap:
                    samples.append(newD)
            newq = tuple(newD if i == j else quad[i] for i in range(4))
            qk = quad_key(newq)
            if qk not in seen_quads:
                seen_quads.add(qk)
                stack.append((newq, j))
    return sorted(seen.values()), samples


def check_gasket():
    from omega import is_schmidt
    Kmax = 30000
    curv, samples = gasket_count(Kmax)
    curv = np.array(curv)
    # every generated disk is a Schmidt disk contained in H (alpha >= 1)
    small_ok = all(
        is_schmidt(Disk(A, G(*B), C)) and B[1] >= 1
        for (A, B, C) in samples)
    counts = {K: int(np.sum(curv <= K)) for K in (2, 8, 20, 100, 1000, 10000, 30000)}
    ok_small = counts[2] == 1 and counts[8] == 3
    Ks = np.array([100, 300, 1000, 3000, 10000, 30000], dtype=float)
    Ns = np.array([np.sum(curv <= K) for K in Ks], dtype=float)
    A = np.vstack([np.log(Ks), np.ones(len(Ks))]).T
    delta, _ = np.linalg.lstsq(A, np.log(Ns), rcond=None)[0]
    ok = ok_small and abs(delta - 1.30568) < 0.03 and small_ok
    report("5. atoms = strip-gasket disks: per-period count N(K) ~ K^delta, "
           "delta = 1.30568... (Apollonian dimension; lambda_0 = delta(2-delta))",
           ok, f"N(K) = {counts}, fitted delta = {delta:.4f}")


# ------------------------------------------------------------ main
def main():
    print("orthospectrum_verify.py — checks for spectral-geometry/orthogeodesics.md")
    print("=" * 74)
    check_axis_exact()
    check_orthogeodesic_geometry()
    check_counts()
    check_hurwitz_growth()
    check_gasket()
    print("=" * 74)
    n_ok = sum(results)
    print(f"{n_ok}/{len(results)} checks pass")
    return 0 if n_ok == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
