"""Machine verification for atomic-census.md.

The Gamma-double cosets of the monoid Omega (Gamma = SL(2,Z)) are the
SL(2,Z)-classes of Schmidt disks in H, i.e. (level n, form class of disc 1-n^2).
Unique factorisation in Omega assigns to each class a *word* in the atom
alphabet (its Apollonian address).  This script computes the full address
census for all odd levels 3 <= n <= NMAX and verifies, in exact arithmetic:

  1. [cells]        the census total at level n is h_+(1-n^2), the number of
                    reduced positive definite forms of disc 1-n^2 (all forms,
                    imprimitive included), and the word map is well defined
                    (re-multiplying the factors reproduces the class).
  2. [atoms]        depth-1 classes = maximal disks; their levels/counts match
                    the atom table of half-plane-monoid.md.
  3. [ford-first]   #{classes at level n whose first letter is [T_i]}
                        = sum_{q=1}^{(n-1)/2} #{x mod 2q : x even,
                                                x^2 + n^2 = 1 mod 4q}.
  4. [ford-ford]    #{classes with word exactly ([T_i],[T_i])} = phi(c) if
                    n = 2c^2+1, else 0.
  5. [reversal]     word(sigma(X)) = reverse of the letterwise sigma-images of
                    word(X); the induced map on level-n classes is an
                    involution (the class map hat-sigma of involution.md).
  6. [superadd]     eps_{n_i} >= eps_{a_i} * eps_{n_{i+1}} along every chain
                    (eps_n = n + sqrt(n^2-1)); checked exactly.
  7. [pairing]      alpha(XY) = <M_X, M_{Y^{-1}}> on random products in Omega.
  8. [asym]         (separate mode) ford-first density -> 3/pi.

Run:  python3 scripts/atomic_census.py [NMAX]        (default 41)
      python3 scripts/atomic_census.py asym [NMAX]   (density experiment)
"""

import os
import sys
from fractions import Fraction
from math import isqrt, gcd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from omega import (G, Mat, Disk, ID, H, T, S, disk_of, in_omega, inversive,
                   matrix_for_disk, factor, is_maximal)

CHECKS = {"pass": 0, "fail": 0}


def check(label, ok, detail=""):
    CHECKS["pass" if ok else "fail"] += 1
    if not ok:
        print(f"  FAIL [{label}] {detail}")


# ------------------------------------------------------------------ forms
def reduce_form(a, b, c):
    """Canonical Gauss-reduced representative of the class of (a,b,c),
    positive definite."""
    assert a > 0 and b * b - 4 * a * c < 0
    while True:
        if c < a or (c == a and b < 0):
            a, b, c = c, -b, a
            continue
        # normalise b into (-a, a]
        if b > a or b <= -a:
            k = (a - b) // (2 * a)          # b + 2ak in (-a, a]
            b2 = b + 2 * a * k
            c2 = a * k * k + b * k + c
            b, c = b2, c2
            continue
        if a == c and b < 0:
            b = -b
            continue
        return (a, b, c)


def reduced_forms(D):
    """All reduced positive definite integral forms of discriminant D < 0
    (imprimitive included).  |list| = h_+(D)."""
    out = []
    amax = isqrt(-D // 3) + 1
    for a in range(1, amax + 1):
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a):
                continue
            c = (b * b - D) // (4 * a)
            if c < a:
                continue
            if a == c and b < 0:
                continue
            out.append((a, b, c))
    return out


def disk_from_form(n, f):
    """Level-n Schmidt disk with f_D = (q,-x,m) = f, i.e. x = -b."""
    a, b, c = f
    return Disk(-2 * a, G(-b, n), -2 * c)


def form_of_disk(D):
    """(q, -x, m), canonically reduced -> the class of the disk."""
    return reduce_form(D.n, -D.B.re, D.m)


# ------------------------------------------------------------------ words
def atom_key(A):
    """Class key of an atom A: 'T' for the Ford stratum (alpha = 1), else
    (alpha, reduced form)."""
    D = disk_of(A)
    al = D.alpha
    if al == 1:
        return "T"
    return (al, form_of_disk(D))


def key_str(k):
    if k == "T":
        return "T"
    al, f = k
    return f"A{al}"


def sigma(X):
    """sigma(X) = conj(X)^{-1}."""
    Xc = Mat(X.a.conj(), X.b.conj(), X.c.conj(), X.d.conj())
    return Xc.inv()


def word_and_suffix_levels(X, nmax=800):
    """Atomic factorisation of X: returns (atom matrices, unit, keys,
    suffix levels [alpha(A_1...A_k U), alpha(A_2...A_k U), ..., alpha(U)=1])."""
    atoms, U = factor(X, nmax)
    keys = [atom_key(A) for A in atoms]
    # suffix levels from the right
    levels = [1]                      # alpha of the unit's disk H
    P = U
    for A in reversed(atoms):
        P = A * P
        levels.append(disk_of(P).alpha)
    levels.reverse()                  # levels[i] = alpha(A_{i+1} ... A_k U)
    return atoms, U, keys, levels


# ------------------------------------------------------------------ closed forms
def ford_first_count(n):
    """sum_{q=1}^{(n-1)/2} #{x mod 2q : x even, x^2 + n^2 = 1 mod 4q}."""
    total = 0
    for q in range(1, (n - 1) // 2 + 1):
        total += sum(1 for x in range(0, 2 * q, 2)
                     if (x * x + n * n - 1) % (4 * q) == 0)
    return total


def euler_phi(c):
    out, m, p = 1, c, 2
    while p * p <= m:
        if m % p == 0:
            e = 0
            while m % p == 0:
                m //= p
                e += 1
            out *= (p - 1) * p ** (e - 1)
        p += 1
    if m > 1:
        out *= m - 1
    return out


def eps_super_ok(nbig, a, nsmall):
    """Exact check of eps_{nbig} >= eps_a * eps_{nsmall}, i.e.
    nbig - a*nsmall >= 0 and (nbig - a*nsmall)^2 >= (a^2-1)(nsmall^2-1)."""
    t = nbig - a * nsmall
    return t >= 0 and t * t >= (a * a - 1) * (nsmall * nsmall - 1)


# ------------------------------------------------------------------ census
def census(NMAX=41, verbose=True):
    grand = {}
    for n in range(3, NMAX + 1, 2):
        Dn = 1 - n * n
        forms = reduced_forms(Dn)
        rows = {}                     # reduced form -> (keys, levels)
        for f in forms:
            D = disk_from_form(n, f)
            X = matrix_for_disk(D)
            check("rep", disk_of(X) == D, f"n={n} f={f}")
            atoms, U, keys, levels = word_and_suffix_levels(X)
            # 1. multiply back
            P = U
            for A in reversed(atoms):
                P = A * P
            check("refactor", disk_of(P) == D, f"n={n} f={f}")
            check("level", levels[0] == n, f"n={n} f={f} levels={levels}")
            rows[f] = (atoms, U, keys, levels)

        # 2. atoms = depth 1 = maximal
        for f, (atoms, U, keys, levels) in rows.items():
            mx = is_maximal(disk_from_form(n, f), nmax=800)
            check("maximal", (len(keys) == 1) == mx, f"n={n} f={f}")
        a_n = sum(1 for r in rows.values() if len(r[2]) == 1)

        # 3. ford-first closed form
        bT = sum(1 for r in rows.values() if len(r[2]) >= 1 and r[2][0] == "T")
        check("ford-first", bT == ford_first_count(n),
              f"n={n}: census {bT} vs formula {ford_first_count(n)}")

        # 4. ford-ford = phi(c) at n = 2c^2+1
        nTT = sum(1 for r in rows.values() if r[2] == ["T", "T"])
        c2 = (n - 1) // 2
        expect = euler_phi(isqrt(c2)) if isqrt(c2) ** 2 == c2 else 0
        check("ford-ford", nTT == expect, f"n={n}: {nTT} vs {expect}")

        # 5. reversal under sigma
        perm = {}
        for f, (atoms, U, keys, levels) in rows.items():
            X = matrix_for_disk(disk_from_form(n, f))
            Y = sigma(X)
            check("sigma-in-omega", in_omega(Y), f"n={n} f={f}")
            DY = disk_of(Y)
            check("sigma-level", DY.alpha == n, f"n={n} f={f}")
            fY = form_of_disk(DY)
            perm[f] = fY
            _, _, keysY, _ = word_and_suffix_levels(Y)
            pred = [atom_key(sigma(A)) for A in reversed(atoms)]
            check("reversal", keysY == pred,
                  f"n={n} f={f}: {keysY} vs {pred}")
        for f, fY in perm.items():
            check("sigma-involution", perm.get(fY) == f, f"n={n} f={f}")

        # 6. superadditivity along every chain
        for f, (atoms, U, keys, levels) in rows.items():
            for i, k in enumerate(keys):
                ai = 1 if k == "T" else k[0]
                check("superadd", eps_super_ok(levels[i], ai, levels[i + 1]),
                      f"n={n} f={f} step {i}: {levels[i]} vs "
                      f"({ai},{levels[i+1]})")

        grand[n] = (len(forms), a_n, bT, rows)
        if verbose:
            words = {}
            for _, _, keys, _ in rows.values():
                w = "·".join(key_str(k) for k in keys)
                words[w] = words.get(w, 0) + 1
            R = len(forms) - a_n - bT
            ws = ", ".join(f"{w}:{c}" for w, c in
                           sorted(words.items(), key=lambda t: (len(t[0]), t[0])))
            print(f"n={n:3d}  h+={len(forms):3d}  atoms={a_n}  "
                  f"ford-first={bT:3d}  deep={R:2d}   [{ws}]")
    return grand


# ------------------------------------------------------------------ ideal composition
def ideal_hnf(gens):
    """HNF basis of the Z-module in Z^2 spanned by gens = [(u,v)]:
    returns (d, e, g) with module = Z(d,0) + Z(e,g), d,g > 0, 0 <= e < d."""
    e, g = 0, 0
    resid = []
    for (u, v) in gens:
        if v == 0:
            resid.append(u)
            continue
        if g == 0:
            e, g = u, v
            continue
        # extended gcd of g and v
        a0, b0, gg = _xgcd(g, v)
        # new second vector: a0*(e,g) + b0*(u,v); residual: combination with v-part 0
        e2 = a0 * e + b0 * u
        resid.append(e * (v // gg) - u * (g // gg))
        e, g = e2, gg
    if g < 0:
        e, g = -e, -g
    d = 0
    for u in resid:
        d = gcd(d, u)
    d = abs(d)
    assert d > 0 and g > 0
    e %= d
    return d, e, g


def _xgcd(a, b):
    if b == 0:
        return (1 if a >= 0 else -1), 0, abs(a)
    x0, x1, y0, y1 = 1, 0, 0, 1
    while b:
        q, a, b = a // b, b, a % b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    if a < 0:
        a, x0, y0 = -a, -x0, -y0
    return x0, y0, a


def compose_with_r(n, f):
    """The class [r_n * a_f] as a reduced form; f = (a,b,c) primitive of
    disc D = 1-n^2, a_f = Z a + Z(b/2 + omega) (convention of
    class-formula-proof.md section 0), r_n = [(n-1)/2, omega], omega =
    sqrt(D)/2.  Elements are (u,v) = u + v*omega, omega^2 = D/4."""
    a, b, c = f
    D = 1 - n * n
    r0 = (n - 1) // 2
    # generators of r_n * a_f
    g1 = (a * r0, 0)                       # a * r0
    g2 = (0, a)                            # a * omega
    g3 = ((b // 2) * r0, r0)               # (b/2 + omega) * r0
    g4 = (D // 4, b // 2)                  # (b/2 + omega) * omega
    d, e, g = ideal_hnf([g1, g2, g3, g4])
    assert d % g == 0 and e % g == 0
    A = d // g
    B = 2 * e // g
    Cnum = e * e - g * g * (D // 4)
    assert Cnum % (d * g) == 0
    C = Cnum // (d * g)
    assert B * B - 4 * A * C == D
    return reduce_form(A, B, C)


def content(f):
    return gcd(gcd(f[0], f[1]), f[2])


def pure_twist_checks(grand):
    """[f_{D(sigma X)}] = [r_n][f] on primitive classes (pure twist theorem)."""
    for n, (_, _, _, rows) in grand.items():
        for f in rows:
            if content(f) != 1:
                continue
            X = matrix_for_disk(disk_from_form(n, f))
            got = form_of_disk(disk_of(sigma(X)))
            want = compose_with_r(n, f)
            check("pure-twist", got == want, f"n={n} f={f}: {got} vs {want}")


# ------------------------------------------------------------------ extremal depth
def depth_checks(grand):
    """1 <= depth <= (n+1)/2, attained exactly by the principal class and by
    r_n = ((n-1)/2, 0, (n+1)/2)."""
    for n, (_, _, _, rows) in grand.items():
        kmax = (n + 1) // 2
        attain = {f for f, r in rows.items() if len(r[2]) == kmax}
        princ = (1, 0, (n * n - 1) // 4)
        rn = reduce_form((n - 1) // 2, 0, (n + 1) // 2)
        check("depth-bound", all(len(r[2]) <= kmax for r in rows.values()),
              f"n={n}")
        check("depth-attain", attain == {princ, rn},
              f"n={n}: {attain} vs {{{princ}, {rn}}}")


# ------------------------------------------------------------------ second-letter kernel
def second_letter_checks(grand):
    """N_{(T,[A])}(n) = #{x mod 2q' : the level-a disk (q', x) lies in the
    class [A]}, q' = (n - a)/2  (for [A] = [T_i] this is phi(c), n = 2c^2+1)."""
    for n, (_, _, _, rows) in grand.items():
        # collect the atom letters that occur anywhere up to this level
        letters = {"T"}
        for r in rows.values():
            letters.update(k for k in r[2] if k != "T")
        for k in letters:
            a = 1 if k == "T" else k[0]
            if (n - a) % 2 or n <= a:
                continue
            q = (n - a) // 2
            pred = 0
            for x in range(0, 2 * q, 2):
                if (x * x + a * a - 1) % (4 * q):
                    continue
                D2 = Disk(-2 * q, G(x, a), -2 * ((x * x + a * a - 1) // (4 * q)))
                if k == "T":
                    ok = is_maximal(D2, nmax=800)
                else:
                    ok = form_of_disk(D2) == k[1] and is_maximal(D2, nmax=800)
                if ok:
                    pred += 1
            got = sum(1 for r in rows.values()
                      if len(r[2]) == 2 and r[2][0] == "T" and r[2][1] == k)
            check("second-letter", got == pred,
                  f"n={n} letter {key_str(k)}: census {got} vs kernel {pred}")


# ------------------------------------------------------------------ hard-coded sigma checks
def sigma_class_checks(grand):
    """In the inner-disk normalisation of Omega the involution acts as the
    PURE TWIST [f] -> [r_n][f] (no inversion; the inversion in involution.md's
    hat-sigma = [r_n][f]^{-1} comes from its reflection convention).
    Decisive test: n = 9, Cl(-80) = Z/4 generated by g = (3,2,7), r_9 = g^2 =
    (4,0,5): pure twist swaps (3,2,7) <-> (3,-2,7); twist-inverse would fix
    them.  At n = 11 (2-torsion class group) the two formulas agree."""
    want = {
        11: [((1, 0, 30), (5, 0, 6)), ((2, 0, 15), (3, 0, 10))],
        9: [((1, 0, 20), (4, 0, 5)), ((4, 0, 5), (1, 0, 20)),
            ((3, 2, 7), (3, -2, 7)), ((3, -2, 7), (3, 2, 7))],
    }
    for n, pairs in want.items():
        if n not in grand:
            continue
        rows = grand[n][3]
        perm = {}
        for f in rows:
            X = matrix_for_disk(disk_from_form(n, f))
            perm[f] = form_of_disk(disk_of(sigma(X)))
        for f, g in pairs:
            check("sigma-classmap", perm.get(f) == g,
                  f"n={n}: sigma{f} = {perm.get(f)} expected {g}")


# ------------------------------------------------------------------ pairing formula
def pairing_checks(trials=200, seed=1):
    import random
    random.seed(seed)
    # a pool of elements of Omega: atoms at levels 1, 7, 17 and translates
    pool = [matrix_for_disk(Disk(0, G(0, 1), -2))]                   # T_i
    pool.append(matrix_for_disk(Disk(-2, G(0, 1), 0)))               # Ford at 0
    pool.append(matrix_for_disk(Disk(-8, G(4, 7), -8)))              # atom a=7
    pool.append(matrix_for_disk(Disk(-18, G(-6, 17), -18)))          # atom a=17
    units = [ID, S, T(1), T(-1), S * T(2), T(3) * S]
    import itertools
    for _ in range(trials):
        X = random_omega(pool, units)
        Y = random_omega(pool, units)
        lhs = disk_of(X * Y).alpha
        rhs = inversive(disk_of(X.inv()), disk_of(Y))
        check("pairing", lhs == rhs, f"{lhs} vs {rhs}")
        # and the mirrored unfolding  alpha(YX) = <M_{Y^{-1}}, M_X>
        check("pairing2", disk_of(Y * X).alpha
              == inversive(disk_of(Y.inv()), disk_of(X)), "mirror")
        # superadditivity of the product itself
        a, b = disk_of(X).alpha, disk_of(Y).alpha
        check("superadd-prod", eps_super_ok(lhs, a, b), f"{lhs},{a},{b}")


def random_omega(pool, units):
    import random
    X = units[random.randrange(len(units))]
    for _ in range(random.randrange(1, 4)):
        X = X * pool[random.randrange(len(pool))]
        X = X * units[random.randrange(len(units))]
    assert in_omega(X)
    return X


# ------------------------------------------------------------------ density
def h_plus(D):
    return len(reduced_forms(D))


def asym(NMAX=801):
    from math import pi
    sb = sh = 0
    print("  X    sum bT    sum h+    ratio    3/pi = %.6f" % (3 / pi))
    marks = [101, 201, 401, 601, 801, 1201, 1601]
    for n in range(3, NMAX + 1, 2):
        sb += ford_first_count(n)
        sh += h_plus(1 - n * n)
        if n in marks or n == NMAX:
            print(f"{n:5d}  {sb:8d}  {sh:8d}   {sb/sh:.6f}")


# ------------------------------------------------------------------ main
if __name__ == "__main__":
    args = sys.argv[1:]
    if args and args[0] == "asym":
        asym(int(args[1]) if len(args) > 1 else 801)
        sys.exit(0)
    NMAX = int(args[0]) if args else 41
    grand = census(NMAX)
    sigma_class_checks(grand)
    pure_twist_checks(grand)
    depth_checks(grand)
    second_letter_checks(grand)
    pairing_checks()
    print(f"\nchecks: {CHECKS['pass']} passed, {CHECKS['fail']} failed")
    sys.exit(1 if CHECKS["fail"] else 0)
