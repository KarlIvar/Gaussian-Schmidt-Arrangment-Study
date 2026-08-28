"""The phase units u_f of one level are the roots of ONE integer polynomial.

Companion to section 5.9 of moduli-invariants.md.

For fixed odd n let D = 1 - n^2, h = h(D), and let u_f = eps * Theta_f
(f running over the h primitive classes) be the phase units of
moduli-invariants.md section 4.  Theorem A (section 5.8) says the multiset
{u_f^12} is permuted by Gal(Qbar/Q), hence

    Psi_n(y) = prod_f (y - u_f^12)  in  Q[y]        (monic, self-reciprocal)

and P_n(x) = den(Psi_n) * Psi_n(x^12) in Z[x] kills every u_f.  Conjecturally
(certified here) already the first-power polynomial

    Phi_n(x) = prod_f (x - u_f)     in  Q[x]

is rational, and Q_n(x) = den(Phi_n) * Phi_n(x) in Z[x] is the small answer.

What this script does:
  * computes u_f at high precision,
  * certifies the coefficients of Phi_n as rationals (safe two-term PSLQ /
    continued-fraction fits: the fit is accepted only when it is
    overdetermined by a large number of spare digits),
  * clears denominators -> Q_n in Z[x], factors the leading coefficient
    (Gross-Zagier primes, section 5.7), checks Q_n(u_f) = 0 numerically,
  * derives Psi_n exactly from Phi_n by Newton's identities (rational
    arithmetic) and, independently, certifies Psi_n directly from the
    numerical u_f^12 -- the unconditional statement, which does not use the
    rationality of Phi_n.

Usage:
    python3 scripts/uf_integer_polynomial.py                 # n = 3..15
    python3 scripts/uf_integer_polynomial.py --dps 500 7 9   # chosen levels
    python3 scripts/uf_integer_polynomial.py --direct 11     # + direct Psi
"""
import sys
from fractions import Fraction
from math import gcd

sys.path.insert(0, 'scripts')
from mpmath import mp, mpf, mpc, fabs, nstr, mpmathify

DEFAULT_DPS = {3: 120, 5: 160, 7: 200, 9: 260, 11: 300, 13: 340, 15: 460,
               17: 420}
# extra precision for the direct (unconditional) certification of Psi_n:
# its coefficients are 12th powers, hence ~12x as long.
DIRECT_DPS = {3: 120, 5: 400, 7: 600, 9: 900, 11: 1100, 13: 1300,
               15: 2600, 17: 1800}


def _load(dps):
    """(Re)import the numerics at working precision dps."""
    mp.dps = dps
    import moduli_invariants as MI
    from involution_classmap import classes_of_disc, reduce_form, is_primitive
    from involution_experiments import inv_sl2
    from proof_check import build_P
    from mpmath import sqrt as msqrt

    def u_values(n):
        D = 1 - n * n
        prim = [f for f in classes_of_disc(D) if is_primitive(f)]
        eps = n + msqrt(mpf(n * n - 1))
        return prim, {f: eps * MI.theta_integral(inv_sl2(build_P(n, f)[0]))[0]
                      for f in prim}

    def r_class(n):
        return reduce_form((n - 1) // 2, 0, (n + 1) // 2)

    return u_values, r_class


# ---------------------------------------------------------------- rationals

def rat_recognize(x, dps, guard=None):
    """Recognize the real mpf x as a rational p/q by continued fractions.

    Returns (Fraction(p, q), spare) where `spare` is the number of decimal
    digits by which the fit is overdetermined: the rational carries
    len(p)+len(q) digits of information and reproduces x to `acc` digits, so
    spare = acc - (digits(p) + digits(q)).  A fit is trustworthy only when
    spare is large (we demand >= dps/5).  Returns (None, reason) on failure.
    """
    guard = guard if guard is not None else max(10, dps // 20)
    acc_target = dps - guard
    a = mpf(x)
    # continued fraction with exact integer convergents
    p0, q0, p1, q1 = 1, 0, 0, 1        # h_{-1}/k_{-1}, h_{-2}/k_{-2}
    t = a
    best = None
    for _ in range(4 * dps):
        c = int(mp.floor(t))
        p0, p1 = c * p0 + p1, p0
        q0, q1 = c * q0 + q1, q0
        if q0 == 0:
            break
        err = fabs(a - mpf(p0) / mpf(q0))
        if err == 0:
            best = (Fraction(p0, q0), mp.inf)
            break
        acc = -mp.log10(err) + (mp.log10(fabs(a)) if a != 0 else 0)
        info = len(str(abs(p0))) + len(str(abs(q0)))
        if acc >= acc_target:
            best = (Fraction(p0, q0), float(acc - info))
            break
        frac = t - c
        if frac == 0:
            best = (Fraction(p0, q0), mp.inf)
            break
        t = 1 / frac
    if best is None:
        return None, 'no convergent reached target accuracy'
    fr, spare = best
    if spare < dps / 5.0:
        return None, f'unsafe fit (only {spare:.0f} spare digits)'
    return fr, spare


def spare_str(spare):
    """Format a spare-digit count that may be mp.inf (an exact fit)."""
    return 'inf' if spare == mp.inf else f'{float(spare):.0f}'


def poly_from_roots(roots):
    """Coefficients of prod (x - r), highest degree first."""
    co = [mpc(1)]
    for r in roots:
        new = [mpc(0)] * (len(co) + 1)
        for k, c in enumerate(co):
            new[k] += c
            new[k + 1] -= r * c
        co = new
    return co


def clear_denominators(fracs):
    """[Fraction] -> (integer coefficients, common denominator)."""
    den = 1
    for f in fracs:
        den = den * f.denominator // gcd(den, f.denominator)
    ints = [int(f * den) for f in fracs]
    g = 0
    for c in ints:
        g = gcd(g, c)
    if g > 1:
        ints = [c // g for c in ints]
        den //= gcd(den, g)
    return ints, den


def factor_int(m):
    m = abs(m)
    out = []
    d = 2
    while d * d <= m:
        e = 0
        while m % d == 0:
            m //= d
            e += 1
        if e:
            out.append((d, e))
        d += 1 if d == 2 else 2
    if m > 1:
        out.append((m, 1))
    return out


def show_factor(m):
    return " * ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in factor_int(m))


def gz_tags(m, n):
    """Tag the primes of m by the Gross-Zagier sets of section 5.7."""
    from gz_denominators import gz_primes
    absD = n * n - 1
    g3, g4 = gz_primes(absD, 3), gz_primes(absD, 4)
    tags, ok = [], True
    for p, e in factor_int(m):
        where = [s for s, S in (('-3', g3), ('-4', g4)) if p in S]
        ok &= bool(where)
        tags.append(f"{p}^{e}" if e > 1 else f"{p}")
        tags[-1] += "[" + ",".join(where) + "]" if where else "[NOT GZ]"
    return " ".join(tags), ok


def poly_str(co, var='x'):
    deg = len(co) - 1
    parts = []
    for k, c in enumerate(co):
        if c == 0:
            continue
        d = deg - k
        mono = '' if d == 0 else (var if d == 1 else f'{var}^{d}')
        parts.append(('' if not parts and c > 0 else (' + ' if c > 0 else ' - '))
                     + (str(abs(c)) if abs(c) != 1 or d == 0 else '') + mono)
    return ''.join(parts) if parts else '0'


# ------------------------------------------------- symmetric-function algebra

def power_sums_from_monic(a, m_max):
    """Newton: monic x^h + a[0]x^{h-1} + ... + a[h-1] -> power sums s_1..s_max."""
    h = len(a)
    s = []
    for m in range(1, m_max + 1):
        if m <= h:
            acc = -m * a[m - 1]
            for i in range(1, m):
                acc -= a[i - 1] * s[m - i - 1]
            s.append(acc)
        else:
            acc = Fraction(0)
            for i in range(1, h + 1):
                acc -= a[i - 1] * s[m - i - 1]
            s.append(acc)
    return s


def monic_from_power_sums(s, h):
    """Newton the other way: power sums s_1..s_h -> monic coefficients."""
    a = []
    for m in range(1, h + 1):
        acc = s[m - 1]
        for i in range(1, m):
            acc += a[i - 1] * s[m - i - 1]
        a.append(-acc / m)
    return a


def twelfth_power_poly(phi_monic):
    """Phi (monic, rational, degree h) -> Psi(y) = prod (y - root^12)."""
    h = len(phi_monic)
    s = power_sums_from_monic(phi_monic, 12 * h)
    return monic_from_power_sums([s[12 * k - 1] for k in range(1, h + 1)], h)


# ----------------------------------------------------------------- pipeline

def certify_poly(roots, dps, label, var='x'):
    """prod (var - root) -> monic rational coefficients (certified)."""
    co = poly_from_roots(roots)
    assert abs(co[0] - 1) < mpf(10) ** (-dps // 2)
    fracs, spares = [], []
    for c in co[1:]:
        scale = max(fabs(c), mpf(1))
        if fabs(c.imag) > mpf(10) ** (-dps // 2) * scale:
            raise ValueError(f'{label}: non-real coefficient {nstr(c, 10)}')
        fr, sp = rat_recognize(c.real, dps)
        if fr is None:
            raise ValueError(f'{label}: coefficient not recognized ({sp})')
        fracs.append(fr)
        spares.append(sp)
    return fracs, min(spares)


def residual(int_co, roots, sub=1):
    """max relative |P(u)| over the roots, P given by integer coefficients."""
    worst = mpf(0)
    for r in roots:
        v = mpc(0)
        scale = mpf(0)
        x = r ** sub
        for c in int_co:
            v = v * x + c
            scale = scale * fabs(x) + fabs(mpf(c))
        worst = max(worst, fabs(v) / max(scale, mpf(1)))
    return worst


def level(n, dps=None, direct=False):
    dps = dps or DEFAULT_DPS.get(n, 300)
    u_values, r_class = _load(dps + 20)
    prim, u = u_values(n)
    h = len(prim)
    roots = [u[f] for f in prim]

    print("=" * 78)
    print(f"n = {n}   D = {1 - n * n}   h = {h}   r_n = {r_class(n)}   "
          f"working precision {dps} digits")
    print("=" * 78)

    # -- sanity: the two proved laws that make the root multiset symmetric
    prod = mpc(1)
    for r in roots:
        prod *= r
    print(f"  prod_f u_f = {nstr(prod, 12)}   (must be +-1: law 2)")

    # -- first-power polynomial (sharp form; rationality certified, not proved)
    fr, spare = certify_poly(roots, dps, 'Phi')
    ints, den = clear_denominators([Fraction(1)] + fr)
    res = residual(ints, roots)
    print(f"\n  Phi_n(x) = prod_f (x - u_f) has rational coefficients "
          f"(>= {spare_str(spare)} spare digits)")
    print(f"  Q_n(x) = {poly_str(ints)}")
    tags, gzok = gz_tags(ints[0], n)
    print(f"    leading coefficient {ints[0]} = {show_factor(ints[0])}")
    print(f"    Gross-Zagier tags (section 5.7): {tags}"
          f"   {'all GZ primes: OK' if gzok else 'NON-GZ PRIME PRESENT'}")
    print(f"    palindromic: {ints == ints[::-1] or ints == [-c for c in ints[::-1]]}"
          f";  max_f |Q_n(u_f)|/scale = {nstr(res, 3)}")

    # -- unconditional polynomial, degree 12h
    phi_monic = fr
    psi = twelfth_power_poly(phi_monic)
    pints, pden = clear_denominators([Fraction(1)] + psi)
    big = [0] * (12 * h + 1)
    for k, c in enumerate(pints):
        big[12 * k] = c
    res12 = residual(big, roots)
    print(f"\n  Psi_n(y) = prod_f (y - u_f^12) from Phi_n by Newton's identities:")
    print(f"    integer form: {poly_str(pints, 'y')[:300]}"
          f"{' ...' if len(poly_str(pints, 'y')) > 300 else ''}")
    print(f"    P_n(x) = Psi-integer-form at y = x^12: degree {12 * h}, "
          f"max_f |P_n(u_f)|/scale = {nstr(res12, 3)}")

    if direct:
        d2 = DIRECT_DPS.get(n, 12 * dps)
        u_values2, _ = _load(d2 + 20)
        prim2, u2 = u_values2(n)
        roots12 = [u2[f] ** 12 for f in prim2]
        fr2, spare2 = certify_poly(roots12, d2, 'Psi', 'y')
        pints2, _ = clear_denominators([Fraction(1)] + fr2)
        agree = (pints2 == pints)
        print(f"\n  DIRECT certification of Psi_n at {d2} digits "
              f"(no use of Phi_n's rationality): coefficients rational, "
              f">= {spare_str(spare2)} spare digits;  agrees with the Newton "
              f"form: {agree}")
        mp.dps = dps + 20
    return ints, pints


def pair_sums(n, dps=None):
    """The r-pair sums S = u_f + u_{r f} and their arithmetic (sections 5.6,
    5.8, 5.9): rational, or conjugate over an explicitly identified real
    quadratic field."""
    dps = dps or DEFAULT_DPS.get(n, 300)
    u_values, r_class = _load(dps + 20)
    from involution_classmap import compose
    from mpmath import sqrt as msqrt
    prim, u = u_values(n)
    D, rn = 1 - n * n, r_class(n)
    seen, S = set(), []
    for f in prim:
        fr = compose(rn, f, D)
        key = tuple(sorted([f, fr]))
        if key in seen:
            continue
        seen.add(key)
        S.append((key, u[f] + u[fr]))
    print(f"  pair sums at n = {n} ({len(S)} pairs of the {len(prim)} classes):")
    rats = []
    for key, s_ in S:
        fr, sp = (rat_recognize(s_.real, dps) if fabs(s_.imag) <
                  mpf(10) ** (-dps // 2) * max(fabs(s_), mpf(1)) else (None, 'complex'))
        rats.append(fr)
        print(f"    {key[0]},{key[1]}: S = {nstr(s_, 20)}"
              + (f"  = {fr}  ({spare_str(sp)} spare digits)" if fr else
                 f"  (not rational: {sp})"))
    if len(S) == 2 and not all(rats):
        e1, e2 = S[0][1] + S[1][1], S[0][1] * S[1][1]
        f1, s1 = rat_recognize(e1.real, dps)
        f2, s2 = rat_recognize(e2.real, dps)
        if f1 and f2:
            print(f"    the two are quadratic conjugates over Q:")
            print(f"      S_A + S_B = {f1}   ({spare_str(s1)} spare digits)")
            print(f"      S_A S_B   = {f2}   ({spare_str(s2)} spare digits)")
            disc = f1 * f1 - 4 * f2
            r = msqrt(mpf(disc.numerator) / mpf(disc.denominator))
            for m in [d for d in range(2, 2 * (n * n - 1) + 1)
                      if (2 * (n * n - 1)) % d == 0
                      and all(d % (q * q) for q in range(2, int(d ** .5) + 1))]:
                fm, sm = rat_recognize(r / msqrt(mpf(m)), dps)
                if fm:
                    print(f"      S_A - S_B = ({fm}) sqrt({m})"
                          f"   ({spare_str(sm)} spare digits)  ->  the "
                          f"pair-sums lie in Q(sqrt({m}))")
                    break


def main(argv):
    dps = None
    direct = False
    pairs_only = False
    levels = []
    it = iter(range(len(argv)))
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == '--dps':
            i += 1
            dps = int(argv[i])
        elif a == '--direct':
            direct = True
        elif a == '--pairs':
            pairs_only = True
        else:
            levels.append(int(a))
        i += 1
    if not levels:
        levels = [3, 5, 7, 9, 11, 13, 15]
    for n in levels:
        if pairs_only:
            print("=" * 78)
            print(f"n = {n}")
            print("=" * 78)
            pair_sums(n, dps)
        else:
            level(n, dps, direct)
        print()


if __name__ == '__main__':
    main(sys.argv[1:])
