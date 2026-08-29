"""Machine verification for middle-kernel.md: the sums S_f(X,Y) =
sum_{gamma in Gamma} f(X gamma Y) for Gamma-bi-invariant f.

Main object: for radial f = k(alpha), the MASTER FORMULA

    alpha(X gamma Y) = ab + sqrt((a^2-1)(b^2-1)) * cosh d_H(z_X, gamma z_Y),

a = alpha(X), b = alpha(Y), z_X = m_1(sigma(X)) (the CM point carrying
beta_2(X)), z_Y = m_1(Y) (the CM point carrying beta_1(Y)), d_H = hyperbolic
distance, gamma acting by Moebius on H.  Hence S_f is the automorphic kernel
of the MODULAR SURFACE evaluated at the two Heegner points.

Checks:
  1. [master]   the formula, 45 digits, on random Omega-pairs and random gamma;
                and that the produced pairing values are odd integers.
  2. [eps]      eps_{ab+P} = eps_a * eps_b (exact algebra, numeric).
  3. [layer]    the layer integral: int_0^inf k_s(ab + P cosh r) sinh r dr
                = eps_a^{-s} eps_b^{-s} / (s P) for k_s(n) = eps_n^{-s}/sqrt(n^2-1).
  4. [count]    #{gamma in PSL(2,Z) : alpha(X gamma Y) <= T} vs 6T/P
                (complete orbit enumeration, growing T).
  5. [rank1]    direct S_{f_s}(X,Y) vs the rank-one main term (6/s) f_s f_s;
                and the class-DEPENDENCE of S (no-go witness): same levels,
                different classes give different S.

Run: venv/bin/python scripts/middle_kernel.py
"""

import os
import sys
from math import gcd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from omega import (G, Mat, Disk, ID, S, T, disk_of, in_omega, inversive,
                   matrix_for_disk)
from mpmath import mp, mpc, mpf, sqrt as msqrt, fabs, quad, cosh, sinh, exp

mp.dps = 60
TOL = mpf(10) ** (-45)

CHECKS = {"pass": 0, "fail": 0}


def check(label, ok, detail=""):
    CHECKS["pass" if ok else "fail"] += 1
    if not ok:
        print(f"  FAIL [{label}] {detail}")


def sigma(X):
    return Mat(X.a.conj(), X.b.conj(), X.c.conj(), X.d.conj()).inv()


def z_point(V):
    """Hyperbolic center of the disk V(H) (bounded, level >= 2)."""
    D = disk_of(V)
    q, x, y = D.n, D.B.re, D.B.im
    assert q > 0 and y >= 2
    return mpc(mpf(x) / (2 * q), msqrt(mpf(y * y - 1)) / (2 * q))


def coshd(z, w):
    return 1 + fabs(z - w) ** 2 / (2 * z.imag * w.imag)


def eps(n):
    return n + msqrt(mpf(n) * n - 1)


# ------------------------------------------------------------------ sampling
FORM_POOL = [(3, (1, 0, 2)), (5, (1, 0, 6)), (5, (2, 0, 3)),
             (7, (1, 0, 12)), (7, (3, 0, 4)), (9, (3, 2, 7)),
             (9, (1, 0, 20)), (11, (2, 0, 15))]
UNITS = [ID, S, T(1), T(-1), S * T(2), T(3) * S]


def rep_of(n, f):
    a, b, c = f
    return matrix_for_disk(Disk(-2 * a, G(-b, n), -2 * c))


def rand_rep(rng):
    n, f = FORM_POOL[rng.randrange(len(FORM_POOL))]
    return UNITS[rng.randrange(len(UNITS))] * rep_of(n, f) \
        * UNITS[rng.randrange(len(UNITS))]


def rand_gamma(rng):
    g = ID
    for _ in range(rng.randrange(1, 5)):
        g = g * (S if rng.random() < 0.4 else T(rng.randrange(-3, 4)))
    return g


# ------------------------------------------------------------------ 1. master
def master_checks(trials=120, seed=3):
    import random
    rng = random.Random(seed)
    for _ in range(trials):
        X, Y, g = rand_rep(rng), rand_rep(rng), rand_gamma(rng)
        a, b = disk_of(X).alpha, disk_of(Y).alpha
        n = disk_of(X * g * Y).alpha
        z = z_point(sigma(X))
        w = z_point(Y)
        gw = (mpc(g.a.re, g.a.im) * w + mpc(g.b.re, g.b.im)) \
            / (mpc(g.c.re, g.c.im) * w + mpc(g.d.re, g.d.im))
        P = msqrt(mpf(a * a - 1) * (b * b - 1))
        pred = a * b + P * coshd(z, gw)
        check("master", fabs(pred - n) < TOL, f"{pred} vs {n}")
        check("odd", n % 2 == 1 and n >= 1, f"n={n}")


# ------------------------------------------------------------------ 2-3. eps & layer
def eps_and_layer_checks():
    for (a, b) in [(3, 3), (3, 7), (5, 9), (7, 11)]:
        P = msqrt(mpf(a * a - 1) * (b * b - 1))
        check("eps", fabs(eps(a * b + P) - eps(a) * eps(b)) < TOL,
              f"a={a} b={b}")
        for s in (mpf(2), mpf(3.5)):
            def integrand(r):
                n = a * b + P * cosh(r)
                return eps(n) ** (-s) / msqrt(n * n - 1) * sinh(r)
            got = quad(integrand, [0, 5, 15, 40])
            want = eps(a) ** (-s) * eps(b) ** (-s) / (s * P)
            check("layer", fabs(got - want) < mpf(10) ** (-30) * want,
                  f"a={a} b={b} s={s}: {got} vs {want}")


# ------------------------------------------------------------------ orbit enumeration
def orbit_pairings(z, w, Cmax):
    """All gamma in PSL(2,Z) with cosh d(z, gamma w) <= Cmax, returned as the
    exact multiset of cosh-values (complete enumeration)."""
    out = []
    lo = float(w.imag / (z.imag * (Cmax + msqrt(Cmax * Cmax - 1))))
    hi = float(w.imag * (Cmax + msqrt(Cmax * Cmax - 1)) / z.imag)
    rmax = int(float(msqrt(mpf(hi)) / w.imag)) + 2
    for r in range(0, rmax + 1):
        smin = -int(abs(r * float(w.real)) + float(msqrt(mpf(hi))) + 2)
        for s in range(smin, -smin + 1):
            if r == 0 and s != 1:
                continue
            if r > 0 and gcd(r, s) != 1:
                continue
            den2 = fabs(r * w + s) ** 2
            if den2 < lo * (1 - 1e-12) or den2 > hi * (1 + 1e-12):
                continue
            # solve p*s - q*r = 1
            if r == 0:
                p, q = 1, 0
            else:
                # extended gcd for p*s - q*r = 1
                x0, y0, gg = _xgcd(s, -r)
                assert gg == 1
                p, q = x0, y0
            g0w = (p * w + q) / (r * w + s)
            im = g0w.imag
            rad2 = 2 * z.imag * im * Cmax - (z.imag - im) ** 2 \
                - 0 * im  # |Re window|^2 bound from cosh d <= Cmax
            # cosh d <= C  <=>  (Re diff)^2 <= 2 Im z Im gw (C-1) - (Im z - Im gw)^2 + ...
            b2 = 2 * z.imag * im * (Cmax - 1) - (z.imag - im) ** 2
            if b2 < 0:
                continue
            halfw = msqrt(b2)
            t0 = float(z.real - g0w.real)
            for t in range(int(t0 - float(halfw)) - 1,
                           int(t0 + float(halfw)) + 2):
                c = coshd(z, g0w + t)
                if c <= Cmax * (1 + mpf(10) ** (-40)):
                    out.append(c)
    return out


def _xgcd(a, b):
    if b == 0:
        return (1 if a >= 0 else -1), 0, abs(a)
    x0, x1, y0, y1 = 1, 0, 0, 1
    while b:
        qq, a, b = a // b, b, a % b
        x0, x1 = x1, x0 - qq * x1
        y0, y1 = y1, y0 - qq * y1
    if a < 0:
        a, x0, y0 = -a, -x0, -y0
    return x0, y0, a


# ------------------------------------------------------------------ 4-5. counting and rank-one
def kernel_checks():
    # two class-pairs at the SAME levels (a, b) = (9, 9), different classes
    pairs = {}
    for tagX, fX in (("g", (3, 2, 7)), ("e", (1, 0, 20))):
        X = rep_of(9, fX)
        Y = rep_of(9, (1, 0, 20))
        a, b = 9, 9
        P = msqrt(mpf(a * a - 1) * (b * b - 1))
        z, w = z_point(sigma(X)), z_point(Y)
        Cmax = mpf(4000) / P
        vals = orbit_pairings(z, w, Cmax)
        ns = sorted(float(a * b + P * c) for c in vals)
        # integrality of the pairing values
        for nn in ns[:200]:
            check("int-pairing", abs(nn - round(nn)) < 1e-30
                  and round(nn) % 2 == 1, f"{nn}")
        # counting law: N(T) vs 6T/P
        print(f"  levels (9,9), classes ({fX}, (1,0,20)):")
        for Tcut in (500, 1500, 3000, 4000):
            N = sum(1 for nn in ns if nn <= Tcut)
            print(f"    N({Tcut}) = {N:5d}   6T/P = {float(6*Tcut/P):8.1f}"
                  f"   ratio {float(N * P / (6 * Tcut)):.4f}")
        # rank-one comparison at s = 2 and 3
        for s in (mpf(2), mpf(3)):
            Ssum = sum(eps(a * b + P * c) ** (-s)
                       / msqrt((a * b + P * c) ** 2 - 1) for c in vals)
            # tail estimate beyond the enumeration window (density 6/P)
            nmax = a * b + P * Cmax
            tail = 6 / P * eps(nmax) ** (-s) / s
            main = (6 / s) * eps(a) ** (-s) / msqrt(mpf(a * a - 1)) \
                * eps(b) ** (-s) / msqrt(mpf(b * b - 1))
            print(f"    s={float(s)}: S = {float(Ssum):.6e} (+tail<{float(tail):.1e})"
                  f"   rank-one main = {float(main):.6e}"
                  f"   S/main = {float(Ssum/main):.4f}")
            pairs[(tagX, float(s))] = Ssum
    # class dependence (no-go witness): same (a,b), different S
    for s in (2.0, 3.0):
        d = fabs(pairs[("g", s)] - pairs[("e", s)]) / pairs[("e", s)]
        print(f"  class-dependence at s={s}: |S_g - S_e|/S_e = {float(d):.3%}")
        check("class-dep", d > mpf(10) ** (-6), "S should differ across classes")


if __name__ == "__main__":
    master_checks()
    eps_and_layer_checks()
    kernel_checks()
    print(f"\nchecks: {CHECKS['pass']} passed, {CHECKS['fail']} failed")
    sys.exit(1 if CHECKS["fail"] else 0)
