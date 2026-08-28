"""Determine the involution induced by sigma(X) = conj(X)^{-1} on form classes.

For each odd n, classes of positive definite forms of disc D = 1 - n^2
correspond to SL_2(Z)-classes of Schmidt circles with alpha = n. For every
class we build an explicit X in SL_2(Z[i]) whose circle realizes it, compute
the circle of sigma(X), and compare with the conjecture

    sigma[f] = [r_n] * [f]^{-1},   r_n = ((n-1)/2, 0, (n+1)/2).

Composition is implemented via concordant forms (primitive classes only;
imprimitive classes are reported separately).
"""
import sys, math
sys.path.insert(0, 'scripts')
from involution_experiments import (I2, M0, herm_from_triple, herm_triple,
                                    matrix_for_circle, M_of_X, sigma,
                                    to_form_y, reduce_form, opp)
from math import gcd


def classes_of_disc(D):
    """All reduced positive definite forms of discriminant D < 0
    (imprimitive included)."""
    out = []
    b = D % 2
    while 3 * b * b <= -D:
        M4 = b * b - D
        a = max(b, 1)
        while a * a <= M4 // 4:
            if (M4 // 4) % a == 0 and M4 % 4 == 0:
                c = M4 // 4 // a
                out.append((a, b, c))
                if 0 < b < a and a < c:
                    out.append((a, -b, c))
            a += 1
        b += 2
    return sorted(out)


def transform(f, u, v):
    """Equivalent form via unimodular matrix with first column (u, v)."""
    a, b, c = f
    g, r, s = extended(u, v)          # u*s - v*r = 1 form
    assert g == 1
    # matrix [[u, -r], [v, s']] with det 1: find r0, s0 with u*s0 - v*r0 = 1
    # extended(u, v) returns (g, x, y) with u*x + v*y = g
    _, xx, yy = extended(u, v)
    r0, s0 = -yy, xx                  # u*s0 - v*r0 = u*xx + v*yy = 1
    A = a*u*u + b*u*v + c*v*v
    B = 2*a*u*r0 + b*(u*s0 + v*r0) + 2*c*v*s0
    C = a*r0*r0 + b*r0*s0 + c*s0*s0
    return (A, B, C)


def extended(x, y):
    if y == 0:
        return (abs(x), 1 if x >= 0 else -1, 0)
    g, p, q = extended(y, x % y)
    return (g, q, p - (x // y) * q)


def is_primitive(f):
    return gcd(gcd(f[0], f[1]), f[2]) == 1


def compose(f1, f2, D):
    """Gauss composition of primitive forms of discriminant D."""
    a1, b1, c1 = f1
    # find a2' represented by f2, coprime to 2*a1
    found = None
    for u in range(-12, 13):
        for v in range(-12, 13):
            if gcd(u, v) != 1:
                continue
            val = f2[0]*u*u + f2[1]*u*v + f2[2]*v*v
            if val != 0 and gcd(val, 2 * a1) == 1:
                found = (u, v)
                break
        if found:
            break
    assert found, (f1, f2)
    a2, b2, c2 = transform(f2, *found)
    assert gcd(a2, 2 * a1) == 1
    # concordant middle coefficient: B = b1 mod 2a1, B = b2 mod 2a2
    g, p, q = extended(2 * a1, 2 * a2)     # 2a1*p + 2a2*q = g = 2
    assert (b2 - b1) % g == 0
    B = b1 + 2 * a1 * p * ((b2 - b1) // g)
    mod = 2 * a1 * a2
    B %= mod
    assert (B - b1) % (2 * a1) == 0 and (B - b2) % (2 * a2) == 0
    assert (B * B - D) % (4 * a1 * a2) == 0
    return reduce_form(a1 * a2, B, (B * B - D) // (4 * a1 * a2))


def inverse(f):
    return (f[0], -f[1], f[2])


def circle_of_class(f, n):
    """Oriented Hermitian matrix (A>0, y_M = n) of the circle of class f."""
    a, b, c = f
    return herm_from_triple(2 * a, (b, -n), 2 * c)   # B = -zeta, zeta = -b + n i


def oriented_data(M):
    A, B, C = herm_triple(M)
    eps = 1 if A > 0 else -1
    yM = -B[1]
    q, x, m, y = to_form_y(M)
    return eps, yM, reduce_form(q, -x, m)


def run(nmax=41):
    print("n : class  --sigma-->  image   (conjecture r*f^{-1})   eps' yM'")
    all_ok = True
    for n in range(3, nmax + 1, 2):
        D = 1 - n * n
        r = reduce_form((n - 1) // 2, 0, (n + 1) // 2)
        for f in classes_of_disc(D):
            M = circle_of_class(f, n)
            X = matrix_for_circle(M)
            assert M_of_X(X) == M
            eps2, yM2, f2 = oriented_data(M_of_X(sigma(X)))
            assert yM2 == n, "trace/alpha preservation failed?!"
            if is_primitive(f):
                pred = compose(r, inverse(f), D)
                ok = (f2 == pred)
                all_ok &= ok
                tag = "" if ok else "   <-- MISMATCH"
                print(f"{n:2d}: {f} -> {f2}  (pred {pred})  eps'={eps2:+d}{tag}")
            else:
                print(f"{n:2d}: {f} -> {f2}   [imprimitive, content "
                      f"{gcd(gcd(f[0], f[1]), f[2])}]  eps'={eps2:+d}")
    print("\nconjecture sigma[f] = [r_n][f]^{-1} on all primitive classes:",
          "CONFIRMED" if all_ok else "FAILED")


if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 25
    run(nmax)
