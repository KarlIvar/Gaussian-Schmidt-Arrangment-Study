"""Machine verification for product-cocycle.md: invariant-level properties of
the double-coset product (X, Y) -> X Gamma Y.

Three properties are checked:

  R1 (beta_1 rigidity, exact):  <omega_1(X), omega_1(X g Y)> = alpha(Y)
      for every g in Gamma: the first circles of the whole product family lie
      at constant inversive distance alpha(Y) from omega_1(X).
  R2 (beta_2 rigidity, exact):  <omega_2(X g Y), omega_2(Y)> = alpha(X).
  R1' (CM-distance form, 50 digits): cosh d_H(m_1(X), m_1(W)) =
      (alpha_X alpha_W - alpha_Y) / sqrt((alpha_X^2-1)(alpha_W^2-1)).
  C (phase cocycle, 50 digits): with W = X.Y (representative level, all
      branches '+', junction point p = Y(m_2(W)) = X^{-1}(m_1(W))):

      W'(m2W) X'(m2X) Y'(m2Y)
          = [ (m1W - m1X)(p - m1Y) / ((p - m2X)(m2W - m2Y)) ]^2 ,

      which upon multiplying by kernel ratios gives
      Theta(W) Theta(X) Theta(Y) = K * Q^2 with K the product of kernel
      values at the six hyperbolic centers.  Also verified with the j'-kernel.

Run: python3 scripts/product_cocycle.py [trials]
"""

import os
import sys
import random
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from omega import (G, Mat, Disk, ID, S, T, disk_of, in_omega, inversive,
                   matrix_for_disk)
from mpmath import mp, mpc, mpf, sqrt as msqrt, fabs

mp.dps = 70
TOL = mpf(10) ** (-50)

CHECKS = {"pass": 0, "fail": 0}


def check(label, ok, detail=""):
    CHECKS["pass" if ok else "fail"] += 1
    if not ok:
        print(f"  FAIL [{label}] {detail}")


# ------------------------------------------------------------------ helpers
def sigma(X):
    return Mat(X.a.conj(), X.b.conj(), X.c.conj(), X.d.conj()).inv()


def to_c(X):
    return [[mpc(X.a.re, X.a.im), mpc(X.b.re, X.b.im)],
            [mpc(X.c.re, X.c.im), mpc(X.d.re, X.d.im)]]


def moebius(Mc, z):
    return (Mc[0][0] * z + Mc[0][1]) / (Mc[1][0] * z + Mc[1][1])


def der(Mc, z):
    return 1 / (Mc[1][0] * z + Mc[1][1]) ** 2


def m1_of(V):
    """Hyperbolic center of the disk V(H) (level >= 3 assumed)."""
    D = disk_of(V)
    q, x, y = D.n, D.B.re, D.B.im
    assert q > 0 and y >= 3
    return mpc(mpf(x) / (2 * q), msqrt(mpf(y * y - 1)) / (2 * q))


def m2_of(V, m1):
    """The fixed point m of Z' = conj(V)^{-1} V with V(m) = m1 (branch +)."""
    Vc = to_c(V)
    Vbar = [[e.conjugate() for e in row] for row in Vc]
    d0 = Vbar[0][0] * Vbar[1][1] - Vbar[0][1] * Vbar[1][0]
    Vbi = [[Vbar[1][1] / d0, -Vbar[0][1] / d0],
           [-Vbar[1][0] / d0, Vbar[0][0] / d0]]
    Z = [[Vbi[0][0] * Vc[0][0] + Vbi[0][1] * Vc[1][0],
          Vbi[0][0] * Vc[0][1] + Vbi[0][1] * Vc[1][1]],
         [Vbi[1][0] * Vc[0][0] + Vbi[1][1] * Vc[1][0],
          Vbi[1][0] * Vc[0][1] + Vbi[1][1] * Vc[1][1]]]
    a, b, c, d = Z[0][0], Z[0][1], Z[1][0], Z[1][1]
    disc = (a + d) ** 2 - 4
    for sgn in (1, -1):
        root = ((a - d) + sgn * msqrt(disc)) / (2 * c)
        if fabs(moebius(Vc, root) - m1) < TOL:
            return root
    return None                      # '-' branch (does not occur here)


def cosh_dist(p, q):
    return 1 + fabs(p - q) ** 2 / (2 * p.imag * q.imag)


# ------------------------------------------------------------------ sampling
FORM_POOL = [(3, (1, 0, 2)), (5, (1, 0, 6)), (5, (2, 0, 3)),
             (7, (1, 0, 12)), (7, (3, 0, 4)), (7, (4, 4, 4)),
             (9, (3, 2, 7)), (9, (4, 0, 5)), (11, (2, 0, 15))]
UNITS = [ID, S, T(1), T(-1), T(2), S * T(1), T(-1) * S, T(3) * S * T(-2)]


def rand_rep(rng):
    """Random representative of a random Schmidt class of level >= 3."""
    n, f = FORM_POOL[rng.randrange(len(FORM_POOL))]
    a, b, c = f
    X = matrix_for_disk(Disk(-2 * a, G(-b, n), -2 * c))
    return UNITS[rng.randrange(len(UNITS))] * X * UNITS[rng.randrange(len(UNITS))]


def rand_gamma(rng):
    g = ID
    for _ in range(rng.randrange(1, 5)):
        g = g * (S if rng.random() < 0.4 else T(rng.randrange(-3, 4)))
    return g


# ------------------------------------------------------------------ checks
def run(trials=200, seed=7):
    rng = random.Random(seed)
    done = 0
    while done < trials:
        X, Y, g = rand_rep(rng), rand_rep(rng), rand_gamma(rng)
        Xh = X * g                                  # X-hat absorbs gamma
        W = Xh * Y
        DX, DY, DW = disk_of(X), disk_of(Y), disk_of(W)
        aX, aY, aW = DX.alpha, DY.alpha, DW.alpha
        if aW < 3:
            continue                                 # degenerate level

        # R1 / R2 (exact integer arithmetic)
        check("R1", inversive(DX, DW) == aY,
              f"{inversive(DX, DW)} vs {aY}")
        check("R2", inversive(disk_of(sigma(W)), disk_of(sigma(Y))) == aX,
              f"n={aW}")

        # R1' CM-distance form
        m1X, m1Y, m1W = m1_of(Xh), m1_of(Y), m1_of(W)
        lhs = cosh_dist(m1X, m1W)
        rhs = (mpf(aX * aW - aY)
               / msqrt(mpf(aX * aX - 1) * mpf(aW * aW - 1)))
        check("R1'-dist", fabs(lhs - rhs) < TOL, f"{lhs} vs {rhs}")

        # C: the derivative triple identity
        m2X, m2Y, m2W = m2_of(Xh, m1X), m2_of(Y, m1Y), m2_of(W, m1W)
        check("branch+", None not in (m2X, m2Y, m2W), "minus branch met")
        if None in (m2X, m2Y, m2W):
            continue
        Xc, Yc, Wc = to_c(Xh), to_c(Y), to_c(W)
        p = moebius(Yc, m2W)                        # junction point
        check("junction", fabs(moebius(Xc, p) - m1W) < TOL, "X(p) != m1W")
        lhs = der(Wc, m2W) * der(Xc, m2X) * der(Yc, m2Y)
        Q = (m1W - m1X) * (p - m1Y) / ((p - m2X) * (m2W - m2Y))
        check("cocycle", fabs(lhs - Q * Q) < TOL * max(1, fabs(lhs)),
              f"|lhs - Q^2| = {fabs(lhs - Q*Q)}")
        done += 1

    # full Theta-level identity with the j'-kernel on a few samples
    from moduli_invariants import Jp
    rng = random.Random(11)
    done = 0
    while done < 5:
        X, Y, g = rand_rep(rng), rand_rep(rng), rand_gamma(rng)
        Xh = X * g
        W = Xh * Y
        if disk_of(W).alpha < 3:
            continue
        ms = {}
        for name, V in (("X", Xh), ("Y", Y), ("W", W)):
            m1 = m1_of(V)
            m2 = m2_of(V, m1)
            ms[name] = (V, m1, m2)
        def theta(V, m1, m2):
            h = m2 if m2.imag > 0 else m2.conjugate()   # H-member of pair
            return Jp(m1) * der(to_c(V), m2) / Jp(h).conjugate()
        lhs = theta(*ms["W"]) * theta(*ms["X"]) * theta(*ms["Y"])
        Kk = mpf(1)
        for name in ("W", "X", "Y"):
            V, m1, m2 = ms[name]
            h = m2 if m2.imag > 0 else m2.conjugate()
            Kk = Kk * Jp(m1) / Jp(h).conjugate()
        p = moebius(to_c(Y), ms["W"][2])
        Q = ((ms["W"][1] - ms["X"][1]) * (p - ms["Y"][1])
             / ((p - ms["X"][2]) * (ms["W"][2] - ms["Y"][2])))
        check("theta-cocycle", fabs(lhs - Kk * Q * Q) < TOL * fabs(lhs),
              f"rel err {fabs(lhs - Kk*Q*Q)/fabs(lhs)}")
        done += 1

    print(f"checks: {CHECKS['pass']} passed, {CHECKS['fail']} failed")
    return CHECKS["fail"]


if __name__ == "__main__":
    trials = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    sys.exit(1 if run(trials) else 0)
