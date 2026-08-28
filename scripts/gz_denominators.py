"""Denominators of the certified Schmidt-phase fractions are Gross-Zagier
primes (moduli-invariants.md section 5.7).

For each certified denominator at level n (disc D = 1 - n^2), every prime
factor lies in GZ(D,-3) u GZ(D,-4), where GZ(D,-e) is the set of primes
dividing (e|D| - x^2)/4 for admissible x -- the primes of Deuring-collision
of the level-n CM points with j = 0 (e = 3) resp. j = 1728 (e = 4).
Observed: GZ(D,-3)-primes enter squared, GZ(D,-4)-primes to first power.
"""
from math import isqrt

CERTIFIED = {
    5:  [6647],
    7:  [11891],
    9:  [36559082332399],
    11: [8508413439, 76575720951],
    13: [38032133275, 722610532225],
}


def factor(n):
    f, d = {}, 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def gz_primes(absD, e):
    out, M = set(), e * absD
    for x in range(M % 2, isqrt(M) + 1, 2):
        v = M - x * x
        if v > 0 and v % 4 == 0:
            out |= set(factor(v // 4))
    return out


def run():
    for n, dens in CERTIFIED.items():
        absD = n * n - 1
        g3, g4 = gz_primes(absD, 3), gz_primes(absD, 4)
        for d in dens:
            fac = factor(d)
            assert all(p in g3 or p in g4 for p in fac), (n, d)
            desc = " * ".join(f"{p}^{e}" if e > 1 else str(p)
                              for p, e in sorted(fac.items()))
            tags = {p: ("GZ-3" if p in g3 else "") + ("/GZ-4" if p in g4 else "")
                    for p in fac}
            print(f"n={n:2d}: {d} = {desc}   {tags}")
    print("all certified denominator primes are Gross-Zagier primes: OK")


if __name__ == "__main__":
    run()
