"""Verification of the elementary circle-composition and circle-inverse
recipes stated in circle-composition.md.

Level n (odd), D = 1 - n^2.  A circle at level n is (q, x): radius 1/(2q),
center (x + ni)/(2q), with x^2 = D mod 4q; its form is f = (q, -x, m).

Checked here, for every ordered pair of primitive classes at every odd
n <= 25 (and samples up to 41):

  (1) INVERSE = MIRROR: the mirror circle (q, -x) lies in the inverse class.
      Also: unit-circle inversion (q, x, m) -> (m, x mod 2m, q) inverts the
      class.
  (2) COMPOSITION = CRT/MAGNIFICATION: if gcd(q1, q2) = 1, the unique
      circle (q1*q2, x3) with x3 = x1 mod 2q1 and x3 = x2 mod 2q2 lies in
      the Gauss composition of the two classes.  (Slide-to-coprime is done
      by finding an equivalent form with leading coefficient coprime to q1.)
  (3) AMBIGUOUS = ON A MIRROR: a class is 2-torsion iff its reduced form
      has b = 0, b = a, or a = c (circle centered on Re = 0, Re = 1/2, or
      the unit semicircle respectively).
"""
import sys
from math import gcd
sys.path.insert(0, 'scripts')
from involution_classmap import (classes_of_disc, reduce_form, compose,
                                 inverse, is_primitive)


def principal_form(D):
    return reduce_form(1, 0, -D // 4)


def inv_red(f):
    return reduce_form(f[0], -f[1], f[2])


def rep_with_coprime_lead(f, M, D, bound=25):
    """An SL2(Z)-equivalent form whose leading coefficient is coprime to M."""
    a, b, c = f
    for s in range(-bound, bound + 1):
        for t in range(-bound, bound + 1):
            if gcd(s, t) != 1:
                continue
            val = a * s * s + b * s * t + c * t * t
            if val > 0 and gcd(val, M) == 1:
                # complete (s, t) to an SL2 matrix [[s, u], [t, v]], sv - tu = 1
                # extended euclid on (s, t)
                def egcd(p, q):
                    if q == 0:
                        return (p, 1, 0)
                    g, x, y = egcd(q, p % q)
                    return (g, y, x - (p // q) * y)
                g, v, mu = egcd(s, t)     # s*v + t*mu = g = +-1
                if g == -1:
                    v, mu = -v, -mu
                u, vv = -mu, v            # s*vv - t*u = 1
                assert s * vv - t * u == 1
                a2 = val
                b2 = 2 * a * s * u + b * (s * vv + t * u) + 2 * c * t * vv
                c2 = a * u * u + b * u * vv + c * vv * vv
                assert b2 * b2 - 4 * a2 * c2 == D
                return (a2, b2, c2)
    raise RuntimeError((f, M))


def crt_compose_circles(n, q1, x1, q2, x2):
    """The CRT circle (q1 q2, x3), x3 = x1 (2q1), x2 (2q2); gcd(q1,q2)=1."""
    assert gcd(q1, q2) == 1 and (x1 - x2) % 2 == 0
    k = ((x2 - x1) // 2 * pow(q1, -1, q2)) % q2
    x3 = (x1 + 2 * q1 * k) % (2 * q1 * q2)
    N = n * n - 1
    assert (x3 * x3 + N) % (4 * q1 * q2) == 0
    m3 = (x3 * x3 + N) // (4 * q1 * q2)
    return q1 * q2, x3, m3


def run(nmax=25):
    pairs = ambigs = 0
    for n in range(3, nmax + 1, 2):
        D = 1 - n * n
        classes = [f for f in classes_of_disc(D) if is_primitive(f)]
        one = principal_form(D)
        for F1 in classes:
            a1, b1, c1 = F1
            q1, x1 = a1, (-b1) % (2 * a1)

            # (1) inverse = mirror
            m1 = (x1 * x1 + n * n - 1) // (4 * q1)
            assert reduce_form(q1, x1, m1) == inv_red(F1), (n, F1)
            # unit-circle inversion (q, x, m) -> (m, x, q) also inverts the class
            assert reduce_form(m1, -x1, q1) == inv_red(F1), (n, F1)

            # (3) ambiguous <=> reduced shape (mirror-centered circle)
            two_torsion = (compose(F1, F1, D) == one)
            shape = (b1 == 0 or b1 == a1 or a1 == c1)
            assert two_torsion == shape, (n, F1)
            ambigs += shape

            # (2) composition = CRT
            for F2 in classes:
                a2p, b2p, c2p = rep_with_coprime_lead(F2, q1, D)
                q2, x2 = a2p, (-b2p) % (2 * a2p)
                q3, x3, m3 = crt_compose_circles(n, q1, x1, q2, x2)
                assert reduce_form(q3, -x3, m3) == compose(F1, F2, D), \
                    (n, F1, F2)
                pairs += 1
    print(f"inverse=mirror, ambiguous=mirror-shape ({ambigs} ambiguous), "
          f"and CRT composition verified on {pairs} ordered pairs "
          f"(odd n <= {nmax})")


if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 25)
