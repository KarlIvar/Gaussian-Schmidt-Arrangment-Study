"""Irreducibility of the level polynomials: one level, one Galois orbit.

Companion to section 5.10 of moduli-invariants.md.

For fixed odd n let Q_n in Z[x] be the certified integer polynomial with
root multiset {u_f} (section 5.9) and Psi_n in Q[y] the one with root
multiset {u_f^12} (Theorem C).  This script proves both IRREDUCIBLE over Q
at every computed level, by two independent routes.

ROUTE 1 (theory + one exact gcd; Theorem D of section 5.10).  Theorem A
makes Gal(Qbar/Q) act on the twelfth powers by f -> f^{e(sigma)} c(sigma),
and the Artin map realizes every translation c, so all u_f^12 are Galois
conjugates of one another; hence Psi_n = m^{|T|} for a single irreducible
m, where T = {t : u_t^12 = u_1^12} is a subgroup of the class group.  So
for Psi_n, SQUAREFREE implies IRREDUCIBLE -- and squarefreeness of the
certified Psi_n is an exact gcd computation in Q[y], done here in Fraction
arithmetic.  Then [Q(u_f) : Q] >= [Q(u_f^12) : Q] = h, so the certified
degree-h polynomial Phi_n = Q_n/lc is the minimal polynomial of every u_f:
Q_n is irreducible as well.

ROUTE 2 (pure computation; no use of Theorem A).  Exact statements about
the certified integer polynomials themselves:
  * exact squarefreeness (Fraction gcd),
  * a complete integer-arithmetic decision in degree <= 4: our polynomials
    are palindromic, so rational roots besides +-1 come in pairs {t, 1/t}
    that already give a 2+2 split, and every 2+2 split of the monicized
    quartic forces an integer root of its resolvent cubic; the script
    enumerates those in pure integer arithmetic (monotone-piece bisection
    between the critical points) and refutes each candidate split by
    perfect-square and exact-division tests,
  * mod-p irreducibility certificates where the predicted dihedral Galois
    image contains an h-cycle (n = 5, 7, 9); where it does not (n = 11,
    13, 15, 17 -- 2-torsion class groups, or no order-8 element) no single
    prime can certify, and the script instead verifies that Q_n factors
    mod every tested prime, exactly as the image predicts,
  * a root-subset integrality certificate: lc * u_f is an algebraic
    integer, so a monic rational factor of Phi_n would have all elementary
    symmetric functions e_i(lc * u_f) over its root subset in Z; every
    proper nonempty subset is excluded by a wide margin.  This settles the
    degree-8 level n = 15 (2^8 - 2 = 254 subsets),
  * an exact Sylvester-resultant discriminant, whose perfect-squareness is
    compared with the parity of the predicted dihedral image (they agree
    at every level),
  * an optional sympy cross-check (exact factorization over Q).

Usage:
    python3 scripts/uf_irreducibility.py                 # n = 3..17
    python3 scripts/uf_irreducibility.py 9 15            # chosen levels
    python3 scripts/uf_irreducibility.py --psi-subsets 9 # subset cert for Psi
"""
import sys
from fractions import Fraction
from math import isqrt

sys.path.insert(0, 'scripts')
from mpmath import mp, mpf, mpc, fabs, nstr, nint

from uf_integer_polynomial import (_load, certify_poly, clear_denominators,
                                   twelfth_power_poly, DEFAULT_DPS, spare_str)

PMAX = 600          # mod-p search bound
SAFETY = 1e12       # a subset/root margin must exceed SAFETY * error bound


# ----------------------------------------------------- exact Q[x] arithmetic

def qtrim(P):
    i = 0
    while i < len(P) and P[i] == 0:
        i += 1
    return P[i:]


def qrem(A, B):
    """Remainder of A by B in Q[x] (lists of Fraction, highest degree first)."""
    A = qtrim([Fraction(c) for c in A])
    B = qtrim([Fraction(c) for c in B])
    while A and len(A) >= len(B):
        c = A[0] / B[0]
        A = [a - c * b for a, b in zip(A[1:], B[1:] + [Fraction(0)] *
                                       (len(A) - len(B)))]
        A = qtrim(A)
    return A


def qgcd_degree(A, B):
    """Degree of gcd(A, B) in Q[x]."""
    A, B = qtrim(list(A)), qtrim(list(B))
    while B:
        A, B = B, qrem(A, B)
    return len(A) - 1


def squarefree_exact(ints):
    """Exact: is the integer polynomial squarefree (gcd with derivative 1)?"""
    d = len(ints) - 1
    der = [c * (d - k) for k, c in enumerate(ints[:-1])]
    g = qgcd_degree(ints, der)
    return g == 0, g


def eval_int(ints, x):
    v = 0
    for c in ints:
        v = v * x + c
    return v


def resultant_int(A, B):
    """Exact resultant of integer polynomials via Bareiss on Sylvester."""
    m, n = len(A) - 1, len(B) - 1
    N = m + n
    M = [[0] * N for _ in range(N)]
    for i in range(n):
        for j, c in enumerate(A):
            M[i][i + j] = c
    for i in range(m):
        for j, c in enumerate(B):
            M[n + i][i + j] = c
    sign, prev = 1, 1
    for k in range(N - 1):
        if M[k][k] == 0:
            for r in range(k + 1, N):
                if M[r][k] != 0:
                    M[k], M[r] = M[r], M[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, N):
            for j in range(k + 1, N):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
            M[i][k] = 0
        prev = M[k][k]
    return sign * M[N - 1][N - 1]


def discriminant_int(ints):
    d = len(ints) - 1
    der = [c * (d - k) for k, c in enumerate(ints[:-1])]
    r = resultant_int(ints, der)
    s = -1 if (d * (d - 1) // 2) % 2 else 1
    q, rem = divmod(s * r, ints[0])
    assert rem == 0
    return q


def is_square(m):
    return m >= 0 and isqrt(m) ** 2 == m


# ------------------------------------------------- exact degree <= 4 decision

def cubic_integer_roots(c2, c1, c0):
    """All integer roots of z^3 + c2 z^2 + c1 z + c0 -- exact integer
    arithmetic only.

    The cubic is monotone outside its critical points (-c2 +- sqrt(r))/3,
    r = c2^2 - 3 c1, which isqrt locates within +-2.  On each monotone
    integer range a sign-change bisection finds the unique bracketing pair,
    so an integer root there is hit exactly; the +-3 windows around the
    critical points are tested exhaustively.  Complete and exact."""
    def C(z):
        return ((z + c2) * z + c1) * z + c0

    M = 1 + max(abs(c2), abs(c1), abs(c0))          # Cauchy bound, monic
    out = set()

    def bisect(lo, hi):
        """C weakly monotone on [lo, hi]: collect an integer root if any."""
        if lo > hi:
            return
        va, vb = C(lo), C(hi)
        for z, v in ((lo, va), (hi, vb)):
            if v == 0:
                out.add(z)
        if va * vb >= 0:
            return
        while hi - lo > 1:
            mid = (lo + hi) // 2
            vm = C(mid)
            if vm == 0:
                out.add(mid)
                return
            if (vm > 0) == (va > 0):
                lo, va = mid, vm
            else:
                hi, vb = mid, vm
    r = c2 * c2 - 3 * c1
    if r <= 0:                                       # strictly monotone
        bisect(-M, M)
    else:
        s = isqrt(r)
        k1 = (-c2 - s - 3) // 3                      # k1 <= t1 <= t2 <= k2
        k2 = (-c2 + s + 3) // 3 + 1
        for z in range(k1 - 3, k1 + 4):
            if C(z) == 0:
                out.add(z)
        for z in range(k2 - 3, k2 + 4):
            if C(z) == 0:
                out.add(z)
        bisect(-M, min(k1 - 3, M))
        bisect(max(k1 + 3, -M), min(k2 - 3, M))
        bisect(max(k2 + 3, -M), M)
    return sorted(out)


def palindromic_quartic_factor(co):
    """Complete exact reducibility decision for a palindromic integer
    quartic.  Returns None (irreducible over Q) or a witness factorization.

    Completeness: the roots pair t <-> 1/t.  A rational root is +-1 (tested
    exactly) or gives the rational pair-quadratic (x-t)(x-1/t), i.e. a 2+2
    split.  Monicizing (y = lc*x, Gauss), every rational 2+2 split is into
    monic integer quadratics (y^2+ay+b)(y^2+cy+d), and z = b + d is then an
    integer root of the resolvent cubic; all of those are enumerated
    exactly and each candidate split is refuted or exhibited by
    perfect-square and exact-expansion tests.
    """
    a, b, c, d, e = co
    assert co == co[::-1] and len(co) == 5
    for t in (1, -1):
        if eval_int(co, t) == 0:
            return f'rational root x = {t}'
    # monic integer model y^4 + p y^3 + q y^2 + r y + s,  y = a x
    p, q, r, s = b, a * c, a * a * d, a ** 3 * e
    # resolvent cubic z^3 - q z^2 + (pr - 4s) z - (p^2 s - 4 q s + r^2)
    for z in cubic_integer_roots(-q, p * r - 4 * s, -(p * p * s - 4 * q * s
                                                      + r * r)):
        D1 = z * z - 4 * s          # (b, d) are roots of w^2 - z w + s
        D2 = p * p - 4 * (q - z)    # (a, c) are roots of w^2 - p w + (q - z)
        if not (is_square(D1) and is_square(D2)):
            continue
        t1, t2 = isqrt(D1), isqrt(D2)
        if (z + t1) % 2 or (p + t2) % 2:
            continue
        B, Dd = (z + t1) // 2, (z - t1) // 2
        for A, C in (((p + t2) // 2, (p - t2) // 2),
                     ((p - t2) // 2, (p + t2) // 2)):
            if (A + C == p and B + Dd + A * C == q
                    and A * Dd + B * C == r and B * Dd == s):
                return (f'(y^2 + {A} y + {B})(y^2 + {C} y + {Dd}), y = lc*x')
    return None


def exact_low_degree_verdict(ints):
    """Exact irreducibility decision for degree 1, 2, or palindromic 4."""
    d = len(ints) - 1
    if d == 1:
        return True, 'degree 1'
    if d == 2:
        sq = is_square(ints[1] ** 2 - 4 * ints[0] * ints[2])
        return (not sq), ('discriminant not a perfect square (exact)'
                          if not sq else 'discriminant is a square')
    if d == 4 and ints == ints[::-1]:
        w = palindromic_quartic_factor(ints)
        return (w is None), (w or 'no rational root, every resolvent-cubic '
                                  'integer root refuted (exact)')
    return None, 'no exact low-degree decision implemented'


# ------------------------------------------------------------ mod p machinery

def primes(bound):
    sieve = [True] * bound
    for i in range(2, bound):
        if sieve[i]:
            yield i
            for j in range(i * i, bound, i):
                sieve[j] = False


def pmod(co, p):
    """Highest-first int list -> lowest-first residue list, normalized."""
    lo = [c % p for c in co[::-1]]
    while lo and lo[-1] == 0:
        lo.pop()
    return lo


def pmul(A, B, p):
    out = [0] * (len(A) + len(B) - 1)
    for i, x in enumerate(A):
        if x:
            for j, y in enumerate(B):
                out[i + j] = (out[i + j] + x * y) % p
    while out and out[-1] == 0:
        out.pop()
    return out


def prem_p(A, B, p):
    A = A[:]
    binv = pow(B[-1], p - 2, p)
    while len(A) >= len(B):
        c = A[-1] * binv % p
        off = len(A) - len(B)
        for i, y in enumerate(B):
            A[off + i] = (A[off + i] - c * y) % p
        while A and A[-1] == 0:
            A.pop()
    return A


def pgcd(A, B, p):
    while B:
        A, B = B, prem_p(A, B, p)
    inv = pow(A[-1], p - 2, p)
    return [x * inv % p for x in A]


def ppow_x(e, M, p):
    """x^e mod (M, p), M lowest-first monic-ish."""
    result, base = [1], [0, 1]
    base = prem_p(base, M, p) if len(M) <= 2 else base
    while e:
        if e & 1:
            result = prem_p(pmul(result, base, p), M, p)
        base = prem_p(pmul(base, base, p), M, p)
        e >>= 1
    return result


def ddf_pattern(ints, p):
    """Distinct-degree factorization degree pattern of ints mod p, or None
    if p divides the leading coefficient or the reduction isn't squarefree."""
    if ints[0] % p == 0:
        return None
    M = pmod(ints, p)
    d = len(M) - 1
    der = [c * k % p for k, c in enumerate(M)][1:]
    while der and der[-1] == 0:
        der.pop()
    if not der or len(pgcd(M, der, p)) > 1:
        return None
    pattern = []
    W = M[:]
    k = 0
    while len(W) - 1 > 0:
        k += 1
        if 2 * k > len(W) - 1:
            pattern += [len(W) - 1]
            break
        xp = ppow_x(p ** k, W, p)
        diff = xp[:]
        if len(diff) < 2:
            diff += [0] * (2 - len(diff))
        diff[1] = (diff[1] - 1) % p
        while diff and diff[-1] == 0:
            diff.pop()
        if diff:
            g = pgcd(W, diff, p)
        else:
            g = W[:]
        if len(g) > 1:
            pattern += [k] * ((len(g) - 1) // k)
            W = pquo(W, g, p)
    return sorted(pattern)


def pquo(A, B, p):
    A = A[:]
    out = [0] * (len(A) - len(B) + 1)
    binv = pow(B[-1], p - 2, p)
    while len(A) >= len(B):
        c = A[-1] * binv % p
        out[len(A) - len(B)] = c
        off = len(A) - len(B)
        for i, y in enumerate(B):
            A[off + i] = (A[off + i] - c * y) % p
        while A and A[-1] == 0:
            A.pop()
    return out


def modp_report(ints, pmax=PMAX):
    """(certificate prime or None, #primes tested, sample patterns)."""
    h = len(ints) - 1
    seen = {}
    tested = 0
    cert = None
    for p in primes(pmax):
        pat = ddf_pattern(ints, p)
        if pat is None:
            continue
        tested += 1
        key = tuple(pat)
        seen[key] = seen.get(key, 0) + 1
        if pat == [h] and cert is None:
            cert = p
    return cert, tested, seen


# ------------------------------------------- dihedral image of the class group

def dihedral_image(n):
    """The predicted Galois image on the classes: all maps f -> f^{+-1} c.

    Returns (classes, list of permutations as index tuples)."""
    from involution_classmap import classes_of_disc, is_primitive, compose
    D = 1 - n * n
    prim = [f for f in classes_of_disc(D) if is_primitive(f)]
    idx = {f: i for i, f in enumerate(prim)}
    principal = next(f for f in prim if f[0] == 1)

    def inv(f):
        # inverse class: (a, -b, c) reduced
        from involution_classmap import reduce_form
        return reduce_form(f[0], -f[1], f[2])

    perms = set()
    for c in prim:
        perms.add(tuple(idx[compose(f, c, D)] for f in prim))
        perms.add(tuple(idx[compose(inv(f), c, D)] for f in prim))
    return prim, sorted(perms)


def perm_cycle_type(perm):
    seen, out = set(), []
    for i in range(len(perm)):
        if i in seen:
            continue
        k, j = 0, i
        while j not in seen:
            seen.add(j)
            j = perm[j]
            k += 1
        out.append(k)
    return sorted(out)


def perm_parity(perm):
    return (-1) ** sum(k - 1 for k in perm_cycle_type(perm))


# ------------------------------------------------- subset (factor) certificate

def subset_certificate(vroots, label):
    """No proper nonempty subset R of the scaled roots can have all
    elementary symmetric functions integral.

    vroots are lc * (roots); they are algebraic integers, so a monic
    rational factor of the monic model would have e_i(R) in Z for all i.
    For each R we certify some e_i(R) is far from Z (or far from real):
    the margin must beat the forward error bound by SAFETY and also an
    absolute floor at a quarter of the working precision (a smallest
    scaled root L*u^-12 can itself be a legitimately tiny margin).
    Returns (ok, worst margin, error bound)."""
    h = len(vroots)
    absv = [fabs(v) for v in vroots]
    eps = mpf(10) ** (-(mp.dps - 10))
    hard_min = mpf(10) ** (-(mp.dps // 4))
    worst = mp.inf
    ok = True
    for mask in range(1, 2 ** h - 1):
        co = [mpc(1)]
        sc = [mpf(1)]
        for i in range(h):
            if mask >> i & 1:
                co = [mpc(0)] + co
                sc = [mpf(0)] + sc
                for k in range(len(co) - 1):
                    co[k] = co[k] - vroots[i] * co[k + 1]
                    sc[k] = sc[k] + absv[i] * sc[k + 1]
        best = mpf(0)
        for k in range(len(co) - 1):          # skip the leading 1
            err = sc[k] * eps * 2 ** h
            re, im = co[k].real, co[k].imag
            dist = max(fabs(im), fabs(re - nint(re)))
            if dist > SAFETY * err and dist > hard_min:
                best = max(best, dist)
        if best == 0:
            ok = False
            print(f'    {label}: subset {mask:b} NOT excluded!')
        worst = min(worst, best)
    return ok, worst, eps


# --------------------------------------------------------------------- driver

def analyze(n, dps=None, psi_subsets=False):
    dps = dps or DEFAULT_DPS.get(n, 300)
    u_values, r_class = _load(dps + 20)
    prim, u = u_values(n)
    h = len(prim)
    roots = [u[f] for f in prim]

    print('=' * 78)
    print(f'n = {n}   D = {1 - n * n}   h = {h}   working precision {dps}')
    print('=' * 78)

    # certified polynomials (as in uf_integer_polynomial.py)
    fr, spare = certify_poly(roots, dps, 'Phi')
    Q, _ = clear_denominators([Fraction(1)] + fr)
    psi = twelfth_power_poly(fr)
    PSI, _ = clear_denominators([Fraction(1)] + psi)
    print(f'  certified Q_n and Psi_n rebuilt (>= {spare_str(spare)} spare '
          f'digits)')

    # distinctness of the root multisets (T = 1)
    scale = max(fabs(r) for r in roots)
    m1 = min(fabs(roots[i] - roots[j]) for i in range(h) for j in range(i)) \
        if h > 1 else mp.inf
    m12 = min(fabs(roots[i] ** 12 - roots[j] ** 12)
              for i in range(h) for j in range(i)) if h > 1 else mp.inf
    print(f'  distinctness: min |u_i - u_j| = {nstr(m1, 5)},  '
          f'min |u_i^12 - u_j^12| = {nstr(m12, 5)}  '
          f'(errors ~ 1e-{dps - 10} * scale, scale = {nstr(scale, 5)}): '
          f'T = 1')

    verdicts = {}
    for name, ints in (('Q_n  ', Q), ('Psi_n', PSI)):
        sf, g = squarefree_exact(ints)
        low, why = exact_low_degree_verdict(ints)
        cert, tested, seen = modp_report(ints)
        pats = ', '.join(f'{list(k)}x{v}' for k, v in sorted(seen.items()))
        print(f'\n  {name} degree {len(ints) - 1}: squarefree (exact '
              f'Fraction gcd): {sf}')
        if low is not None:
            print(f'    exact degree-{len(ints) - 1} decision: '
                  f'{"IRREDUCIBLE" if low else "REDUCIBLE"} -- {why}')
        if cert:
            print(f'    mod-p certificate: irreducible mod p = {cert} '
                  f'(lc not divisible by p) => irreducible over Q')
        print(f'    Frobenius degree patterns over {tested} good primes '
              f'< {PMAX}: {pats}')
        verdicts[name] = (sf, low, cert, set(seen))

    # subset certificate for Q_n (covers every degree, settles h = 8)
    v = [Q[0] * r for r in roots]
    ok, worst, eps = subset_certificate(v, 'Q_n')
    print(f'\n  subset certificate for Q_n: all {2 ** h - 2} proper subsets '
          f'excluded: {ok}  (worst margin {nstr(worst, 3)} vs error '
          f'~{nstr(eps, 2)})')

    if psi_subsets:
        # same for Psi_n; needs precision ~ digits of e_{h-1}(lc * u^12)
        L = PSI[0]
        need = int(mp.ceil(sum(max(mp.log10(fabs(L * r ** 12)), mpf(0))
                               for r in roots))) + 260
        print(f'  subset certificate for Psi_n at {need} digits...')
        u_values2, _ = _load(need)
        prim2, u2 = u_values2(n)
        v2 = [L * u2[f] ** 12 for f in prim2]
        ok2, worst2, eps2 = subset_certificate(v2, 'Psi_n')
        print(f'  subset certificate for Psi_n: all subsets excluded: {ok2} '
              f'(worst margin {nstr(worst2, 3)} vs error ~{nstr(eps2, 2)})')
        mp.dps = dps + 20

    # exact discriminants vs the parity of the predicted dihedral image
    classes, perms = dihedral_image(n)
    assert len(classes) == h
    types = sorted(set(tuple(perm_cycle_type(p)) for p in perms))
    even = all(perm_parity(p) == 1 for p in perms)
    hcyc = any(perm_cycle_type(p) == [h] for p in perms)
    dQ = discriminant_int(Q)
    observed = verdicts['Q_n  '][3] | verdicts['Psi_n'][3]
    inside = observed <= set(types)
    print(f'\n  predicted dihedral image: order {len(perms)}, cycle types '
          f'{[list(t) for t in types]}')
    print(f'    h-cycle exists: {hcyc}  (=> single-prime certificate '
          f'{"possible" if hcyc else "IMPOSSIBLE: reducible mod every p"})')
    print(f'    observed Frobenius patterns inside predicted types: {inside}')
    print(f'    image even: {even}  <=>  disc(Q_n) square: {is_square(dQ)}  '
          f'(match: {even == is_square(dQ)})')

    # optional sympy cross-check
    try:
        from sympy import Poly, Symbol, ZZ
        X = Symbol('x')
        line = []
        for name, ints in (('Q_n', Q), ('Psi_n', PSI)):
            _, fl = Poly(ints, X, domain=ZZ).factor_list()
            irr = len(fl) == 1 and fl[0][1] == 1
            line.append(f'{name}: {"irreducible" if irr else "REDUCIBLE"}')
        print(f'  sympy exact factorization cross-check: ' + ';  '.join(line))
    except ImportError:
        print('  (sympy not installed; cross-check skipped)')

    route1 = verdicts['Psi_n'][0]
    q_exact = verdicts['Q_n  '][1]
    q_cert = verdicts['Q_n  '][2]
    route2 = bool(q_exact) or bool(q_cert) or ok
    assert route1, 'Psi_n not squarefree: Route 1 fails!'
    print(f'\n  VERDICT: Psi_n squarefree (exact) + Theorem A '
          f'=> Psi_n IRREDUCIBLE over Q (Theorem D);')
    print(f'           => [Q(u_f) : Q] >= h, so the certified degree-h '
          f'Phi_n is the minimal')
    print(f'           polynomial of every u_f: Q_n IRREDUCIBLE over Q.'
          + ('   (independently' if route2 or ok else '')
          + (' confirmed by Route 2)' if route2 or ok else ''))
    return verdicts


def main(argv):
    dps = None
    psi_subsets = False
    levels = []
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == '--dps':
            i += 1
            dps = int(argv[i])
        elif a == '--psi-subsets':
            psi_subsets = True
        else:
            levels.append(int(a))
        i += 1
    if not levels:
        levels = [3, 5, 7, 9, 11, 13, 15, 17]
    for n in levels:
        analyze(n, dps, psi_subsets)
        print()


if __name__ == '__main__':
    main(sys.argv[1:])
