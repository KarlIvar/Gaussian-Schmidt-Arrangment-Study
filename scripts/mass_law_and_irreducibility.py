"""Proof verification for the Delta-mass law, and irreducibility of the
level polynomials.  Companion to euclidean-moduli-invariants.md sections
5.4-5.5 (Theorem 4 and the irreducibility block).

The Delta-mass law (now a theorem):

    M(n) := prod_{c in Cl(O_n)} n^12 Delta(Lambda_c)/Delta(Z[i])
          = eps(n) * prod_{p^k || n, p not split} p^((6/e_p)(p^k-1)/(p-1) N_e(n/p^k)),

eps(n) = -1 iff n = 2^k (k >= 2).  Proof ingredients checked here:

  M1 (constancy):  A(n)(tau) := prod over ALL index-n sublattices of
      Delta(Lambda)/Delta(L_tau) is constant = (-1)^t(n) prod_{d|n} d^(-12d),
      t(n) = #even divisors.  [checked numerically at random tau]
  M2 (stratification):  Lambda = (delta) Lambda'' with Lambda'' primitive
      gives the exact recursion
      A(n) = prod_{m|n} gamma(m)^(-12 N_e(n/m)) Q(n/m)^r(m),
      gamma(m) = generator of the product of all ideals of norm m,
      Q(k) = product over primitive index-k sublattices.
      [recursion solved in exact rationals for n <= 60 and matched
       against the closed form Q(n) = C(n)^2 n^(-24h)]
  M3 (generating functions):  the three per-prime identities
      alpha + G V = R S  are checked symbolically in sympy.
  M4 (sign): per ambiguous class, Delta(Lambda_c) is real with sign +
      except for the x = n classes (center on Re = 1/2), whose conjugation-
      stable lattices are the Hermite b = d/2 ones; R-classes are positive
      via the (1-i)-rotation trick; the x=n class count is 2^omega(m) for
      4 | n (m = odd part), else 0.  [all checked per class, n <= 20]

Irreducibility:

  H_{-4n^2}        irreducible (classical: ring class field theory; checked)
  H_{-4n^2}^2      reducible by construction (a square)
  P2_n             irreducible at every computed level (2 <= n <= 16),
                   with exact squarefreeness certificates
  P6_n             factored exactly; irreducibility decided per level
  P2_n(x^2)        (the first-power disk polynomial) decided per level

Usage:
    python3 scripts/mass_law_and_irreducibility.py           # mass proof
    python3 scripts/mass_law_and_irreducibility.py irred     # polynomials
    python3 scripts/mass_law_and_irreducibility.py all
"""
import sys
from math import gcd
from fractions import Fraction

sys.path.insert(0, 'scripts')
from euclidean_moduli_invariants import (
    gadd, gsub, gmul, gconj, gnorm, gneg, UNITS, sublattices, is_primitive,
    build_X, zeta_of, reduced_forms, form_of_lattice, E4E6D, sl2_reduce,
    Dq_at, cval, class_reps, class_values, exact_square_poly, min_lambda,
    poly_from_roots, cert_integer, DPS, Omega_period, i_times)
from mpmath import mp, mpf, mpc, fabs, nstr, pi, mpmathify, rand, log
import sympy


def Ne(n):
    out = n
    for p, _ in sympy.factorint(n).items():
        if p % 2 == 1:
            out = out * (p - (1 if p % 4 == 1 else -1)) // p
    return out


def mass_law_value(n):
    """(|M(n)| as exact integer, eps(n))"""
    val = 1
    for p, k in sympy.factorint(n).items():
        rest = n // p ** k
        if p == 2:
            val *= 2 ** (3 * (2 ** k - 1) * Ne(rest))
        elif p % 4 == 3:
            val *= p ** (6 * (p ** k - 1) // (p - 1) * Ne(rest))
    eps = -1 if (n & (n - 1)) == 0 and n % 4 == 0 else 1
    return val, eps


# ------------------------- M1: the full mass is constant -------------------
def check_M1(levels=(2, 3, 4, 5, 6, 8, 9, 12)):
    print("M1: A(n)(tau) constant = (-1)^t(n) prod d^(-12d)   (random tau)")
    mp.dps = 60
    for n in levels:
        t_even = sum(1 for d in range(1, n + 1) if n % d == 0 and d % 2 == 0)
        target = mpf(1)
        for d in range(1, n + 1):
            if n % d == 0:
                target /= mpf(d) ** (12 * d)
        target *= (-1) ** t_even
        for k in range(3):
            tau = mpc(mpf(rand()) * 2 - 1, mpf('0.6') + mpf(rand()))
            A = mpc(1)
            for a in range(1, n + 1):
                if n % a:
                    continue
                d = n // a
                for b in range(d):
                    w = (a * tau + b) / d
                    A *= Dq_at(w) / (mpf(d) ** 12 * Dq_at(tau))
            assert fabs(A - target) < mpf(10) ** (-30) * fabs(target), (n, k)
        print(f"  n={n:2d}: A(n) = {'-' if t_even % 2 else '+'}prod d^-12d  OK "
              f"(3 random tau, 30+ digits)")


# --------------------- M2: exact recursion in rationals --------------------
def gaussian_prime_over(p):
    """pi = x+yi with x^2+y^2 = p (p = 1 mod 4)."""
    for x in range(1, p):
        y2 = p - x * x
        if y2 <= 0:
            break
        y = sympy.integer_nthroot(y2, 2)[0]
        if y * y == y2:
            return (x, y)
    raise RuntimeError(p)


def ideals_of_norm(m):
    """list of generators (Gaussian pairs) of all ideals of norm m."""
    gens = [(1, 0)]
    for p, a in sympy.factorint(m).items():
        new = []
        if p == 2:
            pw = (1, 0)
            for _ in range(a):
                pw = gmul(pw, (1, 1))
            for g in gens:
                new.append(gmul(g, pw))
        elif p % 4 == 3:
            if a % 2:
                return []
            pw = (p ** (a // 2), 0)
            for g in gens:
                new.append(gmul(g, pw))
        else:
            piv = gaussian_prime_over(p)
            pivc = gconj(piv)
            for u in range(a + 1):
                pw = (1, 0)
                for _ in range(u):
                    pw = gmul(pw, piv)
                for _ in range(a - u):
                    pw = gmul(pw, pivc)
                for g in gens:
                    new.append(gmul(g, pw))
        gens = new
    return gens


def gamma12(m):
    """gamma(m)^12 as an exact rational integer (sign included)."""
    prod = (1, 0)
    for g in ideals_of_norm(m):
        prod = gmul(prod, g)
    out = (1, 0)
    for _ in range(12):
        out = gmul(out, prod)
    assert out[1] == 0, m
    return out[0]


def check_M2_M3_exact(N=60):
    print(f"M2/M3: exact recursion  ->  Q(n) = C(n)^2 n^(-24h)   (n <= {N})")
    A = {}
    for n in range(1, N + 1):
        t_even = sum(1 for d in sympy.divisors(n) if d % 2 == 0)
        val = Fraction(1)
        for d in sympy.divisors(n):
            val /= Fraction(d) ** (12 * d)
        A[n] = val * (-1) ** t_even
    Q = {1: Fraction(1)}
    for n in range(2, N + 1):
        rest = A[n]
        for m in sympy.divisors(n):
            if m == 1:
                continue
            g12 = gamma12(m)
            r = len(ideals_of_norm(m))
            if r == 0:
                continue
            rest = rest * Fraction(g12) ** Ne(n // m) / Q[n // m] ** r
        Q[n] = rest   # the m = 1 stratum: gamma(1) = 1, r(1) = 1
        h = len(reduced_forms(-4 * n * n))
        C, _ = mass_law_value(n)
        assert Q[n] == Fraction(C) ** 2 / Fraction(n) ** (24 * h), n
    print(f"  Q(n) from the recursion equals C(n)^2 n^(-24h) for all n <= {N}: OK")

    # sign bookkeeping of the recursion vs M1
    for n in range(2, N + 1):
        s = 1
        for m in sympy.divisors(n):
            r = len(ideals_of_norm(m))
            if r and gamma12(m) < 0:
                s *= (-1) ** Ne(n // m)
        t_even = sum(1 for d in sympy.divisors(n) if d % 2 == 0)
        assert s == (-1) ** t_even, n
    print(f"  sign consistency  prod sign(gamma^12)^Ne = (-1)^t(n): OK (n <= {N})")

    # M3: the generating-function identities, symbolically
    p, T = sympy.symbols('p T', positive=True)
    cases = [
        # (chi, R, G, S)  with V = (1 - chi T)/(1 - p T),
        # alpha = -12 p T /((1-T)(1-pT)^2), S = target Sum q_j T^j
        (1, 1 / (1 - T) ** 2, 12 * T / (1 - T) ** 3,
         -12 * T * (p - 1) / (1 - p * T) ** 2),
        (-1, 1 / (1 - T ** 2), 12 * T ** 2 / (1 - T ** 2) ** 2,
         12 * T / ((1 - T) * (1 - p * T)) - 12 * T * (1 + p) / (1 - p * T) ** 2),
        (0, 1 / (1 - T), 6 * T / (1 - T) ** 2,
         3 * 2 * T / ((1 - T) * (1 - p * T)) - 12 * T * p / (1 - p * T) ** 2),
    ]
    # chi = 0 is p = 2: V = 1/(1-2T); we keep symbolic p and substitute at the
    # end for the ramified case (S above uses nu_j = p^j, s_j = 3(p^j - 1)).
    for chi, R, G, S in cases:
        V = (1 - chi * T) / (1 - p * T)
        alpha = -12 * p * T / ((1 - T) * (1 - p * T) ** 2)
        expr = sympy.simplify(alpha + G * V - R * S)
        if chi == 0:
            expr = sympy.simplify(expr.subs(p, 2))
        assert expr == 0, chi
    print("  M3 generating-function identities alpha + G V = R S: OK "
          "(split, inert, ramified; symbolic)")


# ------------------------------ M4: the sign -------------------------------
def hermite_of_lattice(gens):
    """Hermite (a, b, d): lattice = Z(b + a i) + Z d from generator pairs."""
    (x1, y1), (x2, y2) = gens
    def ext(pq, q):
        if q == 0:
            return (abs(pq), 1 if pq >= 0 else -1, 0)
        g, u, v = ext(q, pq % q)
        return g, v, u - (pq // q) * v
    g, u, v = ext(y1, y2)
    a = abs(g)
    b0 = u * x1 + v * x2
    dre = abs((-y2 // g) * x1 + (y1 // g) * x2)
    return (a, b0 % dre, dre)


def check_M4(N=20):
    print("M4: the sign, class by class")
    for n in range(2, N + 1):
        mp.dps = max(80, 30 + 8 * n)
        reps = {}
        for L in sublattices(n):
            if is_primitive(*L):
                reps.setdefault(form_of_lattice(*L), []).append(L)
        m_odd = n
        while m_odd % 2 == 0:
            m_odd //= 2
        pred_xn = 2 ** len(sympy.factorint(m_odd)) if n % 4 == 0 else 0
        found_xn = 0
        sign_M = 1
        for f, Ls in reps.items():
            amb = (f[1] == 0 or f[0] == f[1] or f[0] == f[2])
            (a, b, d) = Ls[0]
            (A, B), (c, dd) = build_X(a, b, d)
            cc = cval(c)
            R = mpf(n) ** 12 * Dq_at(cval(dd) / cc) / (cc ** 12 * Dq_at(mpc(0, 1)))
            if not amb:
                continue
            # ambiguous: R real; classify I (x = 0 or n) vs R-class (y = n)
            assert fabs(R.imag) < mpf(10) ** (-20) * fabs(R), (n, f)
            zs = [zeta_of(build_X(*L)) for L in Ls]
            xs = {z[0] % (2 * n) for z in zs}
            ys = {z[1] % (2 * n) for z in zs}
            if any(y % n == 0 for y in ys):          # R-class
                # (1-i)Lambda is conjugation-stable and rhombic; Delta > 0
                g1 = gmul((1, -1), (d, 0))
                g2 = gmul((1, -1), (b, a))
                ha, hb, hd = hermite_of_lattice([g1, g2])
                assert 2 * hb == hd, (n, f, "rotated lattice not rhombic")
                assert R.real > 0, (n, f)
            else:                                     # I-class
                ha, hb, hd = a, b, d                  # conj-stable itself
                assert (2 * hb) % hd == 0, (n, f, "not Hermite-symmetric")
                if hb == 0:
                    assert 0 in xs and R.real > 0, (n, f)
                else:
                    assert all(x % n == 0 and x % (2 * n) != 0 for x in xs), (n, f)
                    assert R.real < 0, (n, f)
                    found_xn += 1
                    sign_M *= -1
        assert found_xn == pred_xn, (n, found_xn, pred_xn)
        _, eps = mass_law_value(n)
        assert sign_M == eps, (n, sign_M, eps)
        print(f"  n={n:2d}: ambiguous signs, Hermite flavors, x=n count "
              f"{found_xn} = {pred_xn}, sign(M) = {'+' if eps > 0 else '-'}: OK")


# --------------------------- irreducibility -------------------------------
LAMBDA = {2: 1, 3: 1, 4: 4, 5: 25, 6: 1, 7: 49, 8: 32, 9: 27, 10: 25,
          11: 121, 12: 4, 13: 169}


def irreducibility(levels=None):
    print("irreducibility of the level polynomials")
    print(" n   h | H irr | P2 irr sqfree | P6 irr sqfree | P2(x^2) irr")
    x, y, z = sympy.symbols('x y z')
    levels = levels or list(range(2, 17))
    for n in levels:
        mp.dps = DPS.get(n, 600)
        cv = class_values(n)
        fs = sorted(cv)
        h = len(fs)
        # H
        Hc = poly_from_roots([cv[f][1] for f in fs])
        Hint = [cert_integer(cf)[0] for cf in Hc]
        assert all(v is not None for v in Hint), n
        Hpoly = sympy.Poly(sum(sympy.Integer(cf) * x ** (h - k)
                               for k, cf in enumerate(Hint)), x)
        H_irr = len(sympy.factor_list(Hpoly.as_expr())[1]) == 1 and \
            sympy.factor_list(Hpoly.as_expr())[1][0][1] == 1
        # P2 (minimal lambda)
        us = [cv[f][0] for f in fs]
        lam = LAMBDA.get(n)
        if lam is None:
            fr0, _ = exact_square_poly(us)
            lam = min_lambda([cf for cf, _ in fr0])
        fr, _ = exact_square_poly(us, scale=lam)
        c2 = [int(cf) for cf, _ in fr]
        P2 = sympy.Poly(sum(sympy.Integer(cf) * z ** (h - k)
                            for k, cf in enumerate(c2)), z)
        fl2 = sympy.factor_list(P2.as_expr())[1]
        P2_irr = len(fl2) == 1 and fl2[0][1] == 1 and \
            sympy.Poly(fl2[0][0], z).degree() == h
        P2_sqf = sympy.gcd(P2, P2.diff(z)).degree() == 0
        # P6 with roots (n^2 v)^6 = (n^2/lam)^6 z^3 at z = (lam v)^2:
        # clear the rational scale exactly
        s6 = sympy.Rational(n * n, lam) ** 6
        P6e = sympy.resultant(P2.as_expr(), y - s6 * z ** 3, z)
        P6 = sympy.Poly(sympy.expand(P6e), y)
        if P6.LC() < 0:
            P6 = sympy.Poly(-P6.as_expr(), y)
        P6 = sympy.Poly(P6.as_expr() * sympy.lcm([sympy.fraction(cc)[1]
                        for cc in P6.all_coeffs()]), y)
        assert P6.LC() == 1 and all(cc.is_integer for cc in P6.all_coeffs()), n
        fl6 = sympy.factor_list(P6.as_expr())[1]
        P6_irr = len(fl6) == 1 and fl6[0][1] == 1 and \
            sympy.Poly(fl6[0][0], y).degree() == h
        P6_sqf = sympy.gcd(P6, P6.diff(y)).degree() == 0
        # P2(x^2)
        P2x2 = sympy.Poly(P2.as_expr().subs(z, x ** 2), x)
        flx = sympy.factor_list(P2x2.as_expr())[1]
        Px_irr = len(flx) == 1 and flx[0][1] == 1 and \
            sympy.Poly(flx[0][0], x).degree() == 2 * h
        tick = lambda b: 'yes' if b else 'NO '
        print(f"{n:2d}  {h:2d} |  {tick(H_irr)}  |  {tick(P2_irr)}  {tick(P2_sqf)}"
              f"   |  {tick(P6_irr)}  {tick(P6_sqf)}   |   {tick(Px_irr)}")
        if not (H_irr and P2_irr and P2_sqf and P6_irr and P6_sqf):
            degs = lambda fl, v: sorted(sympy.Poly(g, v).degree()
                                        for g, mm in fl for _ in range(mm))
            print(f"     factor degrees: H {degs(sympy.factor_list(Hpoly.as_expr())[1], x)}, "
                  f"P2 {degs(fl2, z)}, P6 {degs(fl6, y)}, P2(x^2) {degs(flx, x)}")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "mass"
    if mode in ("mass", "all"):
        check_M1()
        print()
        check_M2_M3_exact()
        print()
        check_M4()
    if mode in ("irred", "all"):
        levels = [int(a) for a in sys.argv[2:]] or None
        irreducibility(levels)
