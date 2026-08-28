"""Matrix-level composition of Schmidt circles: verification for
circle-composition.md section 6.

For X in SL_2(Z[i]) with circle at level n, oriented so that
tr(X conj(X)^{-1}) = -2n, the fun-fact matrix Z_X = X conj(X)^{-1} gives

    W_X := (i/2)(Z_X + n I) = [[-x/2, m], [-q, x/2]]      (integral!)

i.e. the circle data (q, x, m) -- hence the form (q, -x, m) -- can be read
off Z_X affinely, with no centers, radii, or Bezout computations.

Composition at matrix level: slide X, Y by translations (left mult by T^k,
i.e. conjugation of Z) until both numerators equal the CRT value x3; then
(for gcd(q_X, q_Y) = 1)

    W_3 = diag(q_Y, 1)^{-1} W_X diag(q_Y, 1) = diag(q_X, 1)^{-1} W_Y diag(q_X, 1)

is integral, EQUAL from both sides, and reads off the composed circle
(q_X q_Y, x3).  A plain matrix product does none of this: circle(XY) is the
Moebius transport X(circle of Y) (wrong level, not class-well-defined),
while tr(Z_X Z_Y^{-1}) = 2n^2 + 2 x1 x2 - 4(q1 m2 + q2 m1) is a geometric
pair invariant.
"""
import sys
from math import gcd
import itertools
sys.path.insert(0, 'scripts')
from involution_experiments import matmul, inv_sl2, M_of_X, ga, gc
from involution_classmap import (classes_of_disc, reduce_form, compose,
                                 is_primitive, oriented_data)
from proof_check import build_P
from composition_check import crt_compose_circles


def conj_mat(X):
    return tuple(tuple(gc(e) for e in row) for row in X)


def cartan(X):
    return matmul(X, inv_sl2(conj_mat(X)))


def W_of(Z, n):
    """W = i(Z + nI)/2; returns (q, x, m, W) with W = [[-x/2, m], [-q, x/2]]."""
    ZnI = ((ga(Z[0][0], (n, 0)), Z[0][1]), (Z[1][0], ga(Z[1][1], (n, 0))))
    W = tuple(tuple((-e[1], e[0]) for e in row) for row in ZnI)   # times i
    for row in W:
        for e in row:
            assert e[1] == 0 and e[0] % 2 == 0, (Z, n)
    W = tuple(tuple(e[0] // 2 for e in row) for row in W)
    assert W[0][0] == -W[1][1]
    return -W[1][0], -2 * W[0][0], W[0][1], W


def Tk(k):
    return (((1, 0), (k, 0)), ((0, 0), (1, 0)))


def small_coprime_rep(f, M, D):
    """Equivalent form with SMALL leading coefficient coprime to M."""
    a, b, c = f
    best = None
    for bound in (3, 6, 12, 25, 50):
        for s in range(-bound, bound + 1):
            for t in range(-bound, bound + 1):
                if gcd(s, t) != 1:
                    continue
                val = a * s * s + b * s * t + c * t * t
                if val > 0 and gcd(val, M) == 1 and (best is None or val < best[0]):
                    best = (val, s, t)
        if best:
            break
    val, s, t = best

    def egcd(p, q):
        if q == 0:
            return (p, 1, 0)
        g, x, y = egcd(q, p % q)
        return (g, y, x - (p // q) * y)
    g, v, mu = egcd(s, t)
    if g == -1:
        v, mu = -v, -mu
    u, vv = -mu, v
    assert s * vv - t * u == 1
    b2 = 2 * a * s * u + b * (s * vv + t * u) + 2 * c * t * vv
    c2 = a * u * u + b * u * vv + c * vv * vv
    assert b2 * b2 - 4 * val * c2 == D
    return (val, b2, c2)


def run(levels=(3, 5, 9, 11, 15, 21)):
    # 1. W read-off, robust under slides
    ok = 0
    for n in levels:
        D = 1 - n * n
        for f in classes_of_disc(D):
            if not is_primitive(f):
                continue
            X = inv_sl2(build_P(n, f)[0])
            for k in (0, 1, -3):
                Z = cartan(matmul(Tk(k), X))
                assert ga(Z[0][0], Z[1][1]) == (-2 * n, 0)
                q, x, m, _ = W_of(Z, n)
                assert reduce_form(q, -x, m) == reduce_form(*f)
                ok += 1
    print(f"W = i(Z + nI)/2 read-off: OK on {ok} matrices")

    # 2. aligned diag-conjugation composes
    pairs = 0
    for n in levels:
        D = 1 - n * n
        classes = [f for f in classes_of_disc(D) if is_primitive(f)]
        for F1, F2 in itertools.product(classes, classes):
            q1, x1 = F1[0], (-F1[1]) % (2 * F1[0])
            g2 = small_coprime_rep(F2, q1, D)
            q2, x2 = g2[0], (-g2[1]) % (2 * g2[0])
            q3, x3, m3 = crt_compose_circles(n, q1, x1, q2, x2)
            X = inv_sl2(build_P(n, (q1, -x1, (x1 * x1 + n * n - 1) // (4 * q1)))[0])
            Y = inv_sl2(build_P(n, (q2, -x2, (x2 * x2 + n * n - 1) // (4 * q2)))[0])
            X = matmul(Tk((x3 - x1) // (2 * q1)), X)
            Y = matmul(Tk((x3 - x2) // (2 * q2)), Y)
            qX, xX, mX, WX = W_of(cartan(X), n)
            qY, xY, mY, WY = W_of(cartan(Y), n)
            assert xX == xY == x3 and (qX, qY) == (q1, q2)
            assert mX % q2 == 0 and mY % q1 == 0
            W3a = ((WX[0][0], WX[0][1] // q2), (WX[1][0] * q2, WX[1][1]))
            W3b = ((WY[0][0], WY[0][1] // q1), (WY[1][0] * q1, WY[1][1]))
            assert W3a == W3b, (n, F1, F2)
            assert (-W3a[1][0], -2 * W3a[0][0], W3a[0][1]) == (q3, x3, m3)
            assert reduce_form(q3, -x3, m3) == compose(F1, F2, D)
            pairs += 1
    print(f"aligned diag(q,1)-conjugation composes: OK on {pairs} ordered pairs")

    # 3. a plain product is transport, not composition
    n = 9
    X = inv_sl2(build_P(n, (3, 2, 7))[0])
    Y = inv_sl2(build_P(n, (4, 0, 5))[0])
    S = (((0, 0), (-1, 0)), ((1, 0), (0, 0)))
    d1 = oriented_data(M_of_X(matmul(X, Y)))
    d2 = oriented_data(M_of_X(matmul(X, matmul(S, Y))))
    print(f"circle(XY):   (eps, level, class) = {d1}")
    print(f"circle(XSY):  (eps, level, class) = {d2}   <- same circles, "
          f"different middle rep")

    # 4. trace pairing identity
    for n in (5, 9):
        D = 1 - n * n
        cl = [f for f in classes_of_disc(D) if is_primitive(f)]
        for F1, F2 in itertools.product(cl, cl):
            Z1 = cartan(inv_sl2(build_P(n, F1)[0]))
            Z2 = cartan(inv_sl2(build_P(n, F2)[0]))
            q1, x1, m1, _ = W_of(Z1, n)
            q2, x2, m2, _ = W_of(Z2, n)
            tr = ga(*[matmul(Z1, inv_sl2(Z2))[i][i] for i in range(2)])
            assert tr == (2 * n * n + 2 * x1 * x2 - 4 * (q1 * m2 + q2 * m1), 0)
    print("trace pairing tr(Z_X Z_Y^{-1}) = 2n^2 + 2x1x2 - 4(q1m2 + q2m1): OK")


if __name__ == "__main__":
    run()
