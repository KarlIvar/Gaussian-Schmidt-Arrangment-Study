"""Modular invariants of Schmidt circles: numerical laboratory for
moduli-invariants.md.

For X in SL_2(C) with circle X(Rhat) in the upper half-plane:
  m1 = hyperbolic center of X(Rhat),
  m2 = hyperbolic center of conj(X)^{-1}(Rhat)  (lower half-plane for
       integral Schmidt X; conj(m2) is the CM point of the sigma-hat class),
  alpha = coth(hyp radius), beta1 = j(m1), beta2 = j(conj(m2)).

Experiments:
  A. Intertwining: X maps {m2, conj(m2)} onto {m1, conj(m1)}
     (fixed points of Z' = conj(X)^{-1} X and Z = X conj(X)^{-1}).
  B. The sixth invariant Theta(X):
       branch '+' (X(m2) = m1):        Theta = j'(m1) X'(m2) / conj(j'(conj m2))
       branch '-' (X(m2) = conj(m1)):  Theta = j'(m1) conj(X'(m2)) / j'(conj m2)
     Invariance under gammaL X gammaR; fiber behaviour; arithmetic nature.
  C. Hilbert class polynomials from beta1-values.
  D. Simultaneous (n-1)/2- and (n+1)/2-isogenies between beta2 and conj(beta1).
  E. Trace table along D = n^2 - 1.
"""
import sys
import random
from math import gcd
sys.path.insert(0, 'scripts')
from mpmath import mp, mpc, mpf, kleinj, diff, fabs, expm, matrix, nstr, arg, pi, nint
from mpmath import sqrt as msqrt

from involution_experiments import matmul, inv_sl2, M_of_X, gc as gconj
from involution_classmap import (classes_of_disc, reduce_form, compose,
                                 is_primitive)
from proof_check import build_P

if mp.dps < 60:
    mp.dps = 60
TOL = mpf(10) ** (-20)


def _E4E6D(tau):
    """E4, E6, Delta via theta constants -- full working precision
    (numerical differentiation of kleinj caps accuracy near 50 digits,
    which silently corrupted early high-precision experiments)."""
    from mpmath import jtheta, exp, pi
    q = exp(mpc(0, 1) * pi * tau)
    t2, t3, t4 = jtheta(2, 0, q), jtheta(3, 0, q), jtheta(4, 0, q)
    a, b, c = t2 ** 4, t3 ** 4, t4 ** 4
    E4 = (a * a + b * b + c * c) / 2
    E6 = (b + c) * (a + b) * (c - a) / 2
    Dl = (t2 * t3 * t4) ** 8 / 256
    return E4, E6, Dl


def J(tau):
    E4, _, Dl = _E4E6D(tau)
    return E4 ** 3 / Dl


def Jp(tau):
    from mpmath import pi
    E4, E6, Dl = _E4E6D(tau)
    return -2 * pi * mpc(0, 1) * E4 * E4 * E6 / Dl


def inv_red(f):
    return reduce_form(f[0], -f[1], f[2])


def cm_point(f):
    a, b, c = f
    return mpc(mpf(-b) / (2 * a), msqrt(mpf(4 * a * c - b * b)) / (2 * a))


def to_c(g):
    return [[mpc(*g[i][j]) for j in range(2)] for i in range(2)]


def moebius(Mc, z):
    return (Mc[0][0] * z + Mc[0][1]) / (Mc[1][0] * z + Mc[1][1])


def hyp_center_of_herm(Mherm):
    (A, B), (_, C) = Mherm
    a, bx, by = A[0], B[0], B[1]
    if a < 0:
        a, bx, by = -a, -bx, -by
    ctr = mpc(mpf(-bx) / a, mpf(-by) / a)
    r = mpf(1) / a
    s = 1 if ctr.imag > 0 else -1
    return mpc(ctr.real, s * msqrt(ctr.imag ** 2 - r ** 2))


def sigma_mat(X):
    return inv_sl2(tuple(tuple(gconj(e) for e in row) for row in X))


def theta_of_complex(Xc, m1, m2):
    """Theta given the two hyperbolic centers and complex matrix Xc."""
    im = moebius(Xc, m2)
    der = 1 / (Xc[1][0] * m2 + Xc[1][1]) ** 2
    if fabs(im - m1) <= fabs(im - m1.conjugate()):
        assert fabs(im - m1) < TOL, fabs(im - m1)
        return Jp(m1) * der / Jp(m2.conjugate()).conjugate(), '+'
    assert fabs(im - m1.conjugate()) < TOL, fabs(im - m1.conjugate())
    return Jp(m1) * der.conjugate() / Jp(m2.conjugate()), '-'


def theta_integral(X):
    m1 = hyp_center_of_herm(M_of_X(X))
    m2 = hyp_center_of_herm(M_of_X(sigma_mat(X)))
    th, br = theta_of_complex(to_c(X), m1, m2)
    return th, br, m1, m2


def circ3(p, q, r):
    ax, ay, bx, by, cx, cy = p.real, p.imag, q.real, q.imag, r.real, r.imag
    Dd = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    ux = ((ax ** 2 + ay ** 2) * (by - cy) + (bx ** 2 + by ** 2) * (cy - ay)
          + (cx ** 2 + cy ** 2) * (ay - by)) / Dd
    uy = ((ax ** 2 + ay ** 2) * (cx - bx) + (bx ** 2 + by ** 2) * (ax - cx)
          + (cx ** 2 + cy ** 2) * (bx - ax)) / Dd
    ctr = mpc(ux, uy)
    return ctr, fabs(p - ctr)


def theta_floating(Xc):
    """Theta for arbitrary complex 2x2 (det 1): circles found numerically."""
    def hypc(ctr, r):
        s = 1 if ctr.imag > 0 else -1
        return mpc(ctr.real, s * msqrt(ctr.imag ** 2 - r ** 2))
    pts = [mpf(1) / 3, mpf(3) / 2, mpf(-7) / 5]
    m1 = hypc(*circ3(*[moebius(Xc, t) for t in pts]))
    det = Xc[0][0] * Xc[1][1] - Xc[0][1] * Xc[1][0]
    Xci = [[Xc[1][1] / det, -Xc[0][1] / det], [-Xc[1][0] / det, Xc[0][0] / det]]
    Xsig = [[e.conjugate() for e in row] for row in Xci]
    m2 = hypc(*circ3(*[moebius(Xsig, t) for t in pts]))
    return theta_of_complex(Xc, m1, m2)[0]


GENS = [((1, 1), (0, 1)), ((1, -1), (0, 1)), ((0, -1), (1, 0))]


def rand_gamma(rng, k=4):
    def gz(m):
        return (((m[0][0], 0), (m[0][1], 0)), ((m[1][0], 0), (m[1][1], 0)))
    out = gz(((1, 0), (0, 1)))
    for _ in range(k):
        out = matmul(out, gz(rng.choice(GENS)))
    return out


def hilbert_poly_coeffs(js):
    coeffs = [mpc(1)]                       # highest degree first
    for r in js:
        coeffs = [coeffs[0]] + [coeffs[k + 1] - r * coeffs[k]
                                for k in range(len(coeffs) - 1)] + [-r * coeffs[-1]]
        # rebuild properly: multiply polynomial by (x - r)
    return coeffs


def poly_from_roots(js):
    coeffs = [mpc(1)]
    for r in js:
        new = [mpc(0)] * (len(coeffs) + 1)
        for k, c in enumerate(coeffs):
            new[k] += c                     # x * c x^k
            new[k + 1] -= r * c
        coeffs = new
    return coeffs                            # highest degree first


def run():
    rng = random.Random(11)

    print("=" * 74)
    print("A/B. intertwining, and the sixth invariant Theta")
    print("=" * 74)
    for n in (3, 5, 7, 11):
        D = 1 - n * n
        for f in classes_of_disc(D):
            if not is_primitive(f):
                continue
            X = inv_sl2(build_P(n, f)[0])
            th, br, m1, m2 = theta_integral(X)
            for _ in range(2):
                Y = matmul(rand_gamma(rng), matmul(X, rand_gamma(rng)))
                th2 = theta_integral(Y)[0]
                assert fabs(th2 - th) < TOL * (1 + fabs(th)), (n, f)
            print(f"n={n:2d} f={str(f):13s} branch {br}  "
                  f"|Theta| = {nstr(fabs(th), 15):>20s}   "
                  f"arg/pi = {nstr(arg(th) / pi, 15):>20s}")
    print("(two-sided invariance verified to ~20 digits for every line)")

    print()
    print("fiber test at n=5, f=(1,0,6): X h_t, h_t = exp(t W') real rotations")
    print("fixing the sigma-circle -- (alpha, beta1, beta2) constant along t:")
    X = inv_sl2(build_P(5, (1, 0, 6))[0])
    Xc0 = to_c(X)
    Ms = M_of_X(sigma_mat(X))
    (A, B), (_, C) = Ms
    a0, bx0, c0 = A[0], B[0], C[0]
    if a0 < 0:
        a0, bx0, c0 = -a0, -bx0, -c0
    Wg = [mpf(bx0) / 2, mpf(c0) / 2, mpf(-a0) / 2, mpf(-bx0) / 2]
    for t in (mpf(0), mpf(1) / 7, mpf(2) / 7, mpf(1) / 2):
        Ht = expm(matrix([[t * Wg[0], t * Wg[1]], [t * Wg[2], t * Wg[3]]]))
        XH = [[Xc0[i][0] * Ht[0, j] + Xc0[i][1] * Ht[1, j] for j in range(2)]
              for i in range(2)]
        th = theta_floating(XH)
        print(f"  t = {nstr(t, 4):>7s}:  |Theta| = {nstr(fabs(th), 18)}   "
              f"arg/pi = {nstr(arg(th) / pi, 18)}")

    print()
    print("=" * 74)
    print("C. Hilbert class polynomials from the beta1-values")
    print("=" * 74)
    for n in (3, 5, 7, 9, 11):
        D = 1 - n * n
        prim = [f for f in classes_of_disc(D) if is_primitive(f)]
        js = [J(cm_point(f)) for f in prim]
        coeffs = poly_from_roots(js)
        out, okint = [], True
        for c in coeffs:
            ci = nint(c.real)
            okint &= (fabs(c - ci) < mpf(10) ** (-12) * (1 + fabs(c)))
            out.append(int(ci))
        print(f"n={n:2d}  disc {D:5d}  H(x) = {out}   "
              f"{'integer OK' if okint else 'NOT integral!'}")
        if D == -8:
            assert out == [1, -8000]
    print("(n=3 anchor: H_{-8}(x) = x - 8000 confirmed)")

    print()
    print("=" * 74)
    print("D. simultaneous modular equations for the sigma-hat pairs")
    print("=" * 74)
    for n in (3, 5, 7, 11, 13):
        D = 1 - n * n
        rn = reduce_form((n - 1) // 2, 0, (n + 1) // 2)
        for f in classes_of_disc(D):
            if not is_primitive(f):
                continue
            f2 = compose(rn, inv_red(f), D)
            tau2, tau1c = cm_point(f2), cm_point(inv_red(f))
            j2 = J(tau2)
            hits = []
            for M in ((n - 1) // 2, (n + 1) // 2):
                found = None
                for a in [d for d in range(1, M + 1) if M % d == 0]:
                    dd = M // a
                    for b in range(dd):
                        if gcd(gcd(a, b), dd) != 1:
                            continue
                        if fabs(J((a * tau1c + b) / dd) - j2) < \
                           mpf(10) ** (-20) * (1 + fabs(j2)):
                            found = (a, b, dd)
                            break
                    if found:
                        break
                hits.append((M, found))
            assert all(h for _, h in hits), (n, f, hits)
            print(f"n={n:2d} f={str(f):13s} -> f2={str(f2):13s}  "
                  + "  ".join(f"deg {M}: tau2 ~ ({h[0]}t+{h[1]})/{h[2]}"
                              for M, h in hits))

    print()
    print("=" * 74)
    print("E. traces of singular moduli along the discriminants n^2 - 1")
    print("=" * 74)
    print("  n     D     Tr_n  (sum of (j-744)/w over primitive classes)")
    for n in range(2, 14):
        D = 1 - n * n
        tot = mpf(0)
        for f in classes_of_disc(D):
            if not is_primitive(f):
                continue
            a, b, c = f
            w = 3 if a == b == c else (2 if (b == 0 and a == c) else 1)
            tot += (J(cm_point(f)) - 744).real / w
        ti = nint(tot)
        assert fabs(tot - ti) < mpf(10) ** (-8) * (1 + fabs(tot))
        print(f" {n:3d} {D:6d}   {int(ti)}")


def deep(levels=(3, 5, 7, 9, 11, 13, 15, 17)):
    """The unit-cocycle laws for u = eps * Theta on the class groups:
      (i)  u(f^{-1}) = conj(u(f))
      (ii) u(r_n f) * u(f) = 1        (=> |u| = 1 on sigma-hat-fixed classes,
                                          prod_f u(f) = +-1)
      (iii) u + 1/u on each r-pair lies in the real genus field
            Q(sqrt(d) : d | n^2-1)   (PSLQ identification)."""
    mp.dps = 80

    def squarefree_divisors(N):
        out = []
        for d in range(2, N + 1):
            if N % d == 0 and all(d % (p * p) for p in range(2, int(d ** .5) + 1)):
                out.append(d)
        return out

    for n in levels:
        D, N = 1 - n * n, n * n - 1
        rn = reduce_form((n - 1) // 2, 0, (n + 1) // 2)
        eps = n + msqrt(mpf(N))
        prim = [f for f in classes_of_disc(D) if is_primitive(f)]
        u = {}
        for f in prim:
            u[f] = eps * theta_integral(inv_sl2(build_P(n, f)[0]))[0]
        ok_i = all(fabs(u[inv_red(f)] - u[f].conjugate())
                   < mpf(10) ** (-40) * (1 + fabs(u[f])) for f in prim)
        ok_ii = all(fabs(u[compose(rn, f, D)] * u[f] - 1) < mpf(10) ** (-40)
                    for f in prim)
        tot = mpf(1)
        for f in prim:
            tot *= u[f]
        print(f"n={n:2d} h={len(prim)}: (i) {'OK' if ok_i else 'FAIL'}  "
              f"(ii) {'OK' if ok_ii else 'FAIL'}  prod u = {nstr(tot, 10)}")
        from mpmath import pslq
        done = set()
        for f in prim:
            fr = compose(rn, f, D)
            key = tuple(sorted([f, fr]))
            if key in done:
                continue
            done.add(key)
            S = u[f] + u[fr]
            if fabs(S.imag) < mpf(10) ** (-40):
                sq = squarefree_divisors(N)
                rel = pslq([S.real, mpf(1)] + [msqrt(mpf(d)) for d in sq],
                           maxcoeff=10 ** 20, maxsteps=200000)
                if rel and rel[0] != 0:
                    terms = " + ".join(
                        f"({-mpf(r) / rel[0]}){b}" for r, b in
                        zip(rel[1:], ['1'] + [f"sqrt({d})" for d in sq]) if r)
                    print(f"   pair {f},{fr}: u + 1/u = {terms}")
                else:
                    print(f"   pair {f},{fr}: u + 1/u = {nstr(S.real, 30)} "
                          f"(coefficients exceed PSLQ bound)")
            else:
                print(f"   pair {f},{fr}: u + 1/u = {nstr(S, 25)} (complex pair)")


def laws(levels=(3, 5, 7, 9, 11, 13)):
    """Representative-level verification of every lemma in the proof of the
    phase-geodesic functional equations (moduli-invariants.md section 5):
      Lemma A: Z (m1,1)^T = -eps (m1,1)^T  and the sigma-side analogue;
      Lemma B: branches X(m2)=m1, conj(X)^{-1}(m1)=m2, X^{-1}(conj m1)=conj m2;
      Theorem (b): Theta(X) conj(Theta(X^{-1})) = eps^{-2}   (rep level);
      Theorem (a): Theta(R conj(X) R) = conj(Theta(X)), R = diag(-1,1);
      kernel independence: same laws for g = j'/j, with change factor
      conj(beta2)/beta1."""
    TOLL = mpf(10) ** (-25)

    def conj_mat(X):
        return tuple(tuple(gconj(e) for e in row) for row in X)

    def theta_kernel(X, g):
        m1 = hyp_center_of_herm(M_of_X(X))
        m2 = hyp_center_of_herm(M_of_X(sigma_mat(X)))
        Xc = to_c(X)
        im = moebius(Xc, m2)
        der = 1 / (Xc[1][0] * m2 + Xc[1][1]) ** 2
        if fabs(im - m1) <= fabs(im - m1.conjugate()):
            return g(m1) * der / g(m2.conjugate()).conjugate()
        return g(m1) * der.conjugate() / g(m2.conjugate())

    for n in levels:
        D = 1 - n * n
        eps = n + msqrt(mpf(n * n - 1))
        for f in classes_of_disc(D):
            if not is_primitive(f):
                continue
            X = inv_sl2(build_P(n, f)[0])
            m1 = hyp_center_of_herm(M_of_X(X))
            m2 = hyp_center_of_herm(M_of_X(sigma_mat(X)))
            Xc = to_c(X)
            Zc = to_c(matmul(X, inv_sl2(conj_mat(X))))
            assert fabs(Zc[1][0] * m1 + Zc[1][1] + eps) < TOLL          # Lemma A
            assert fabs(moebius(Xc, m2) - m1) < TOLL                     # Lemma B
            assert fabs(moebius(to_c(sigma_mat(X)), m1) - m2) < TOLL
            assert fabs(moebius(to_c(inv_sl2(X)), m1.conjugate())
                        - m2.conjugate()) < TOLL
            thX = theta_integral(X)[0]
            thY = theta_integral(inv_sl2(X))[0]
            assert fabs(thX * thY.conjugate() - 1 / eps ** 2) < TOLL     # (b)
            Rm = (((-1, 0), (0, 0)), ((0, 0), (1, 0)))
            thS = theta_integral(matmul(Rm, matmul(conj_mat(X), Rm)))[0]
            assert fabs(thS - thX.conjugate()) < TOLL * (1 + fabs(thX))  # (a)
            g2 = lambda t: Jp(t) / J(t)
            t2X = theta_kernel(X, g2)
            t2Y = theta_kernel(inv_sl2(X), g2)
            assert fabs(t2X * t2Y.conjugate() - 1 / eps ** 2) < TOLL
            beta1, beta2 = J(m1), J(m2.conjugate())
            assert fabs(t2X / thX - beta2.conjugate() / beta1) < \
                TOLL * (1 + fabs(beta2 / beta1))
        print(f"n={n:2d}: Lemmas A, B, Theorems (a), (b), kernel independence: OK")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == 'deep':
        deep()
    elif len(sys.argv) > 1 and sys.argv[1] == 'laws':
        laws()
    else:
        run()
