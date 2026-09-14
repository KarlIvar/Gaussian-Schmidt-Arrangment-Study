"""Horizontal families: Atkin-Lehner-coupled Heegner points and the slices
(D, [r]) of the Schmidt unit systems.

Companion to horizontal-families.md (item 4 of PROGRAM.md).  A SLICE is a pair
(D, [r]): a negative discriminant D and the class of a primitive ambiguous
invertible ideal r of the order O_D (conj(r) = r, r^2 = (N r)); its UNIT SYSTEM
is {R_b = N(r)^6 Delta(b)/Delta(r^{-1} b) : b in Pic(O_D)} (units by Paper II
Thm 4.2 / other-fields.md Thm 4).  The Gaussian arrangement realizes the slices
D = 1 - n^2 with [r] = [((n-1)/2, 0, (n+1)/2)] (odd n, S) and
[(n-1, n-1, n/2)] (even n, iS); Stange's arrangements S_K realize
D = -(k^2-4)/|d_K| at the level 2 alpha = k.  The script verifies:

 (S) SLICES: the primitive ambiguous forms of disc D are (r0, 0, s0)
     [4 r0 s0 = |D|] and (r0, r0, (r0+s0)/4) [r0 s0 = |D|], gcd = 1; their
     ideals are exactly the primitive ambiguous invertible ideals; they come in
     pairs {r, s} with r s = (theta), theta purely imaginary, N(theta) = r0 s0;
     pairs <-> Cl(D)[2] (each 2-torsion class contains exactly one pair),
     #forms = 2^mu (Gauss's mu), [r] = 1 iff min(r0, s0) = 1;
 (R) REALIZATION at primitive strata: S_K (nine class-number-one fields, with
     iS at Q(i)) has a level of discriminant D iff |d_K||D| + 4 = k^2; the
     twist of the level k is (r0, r0, (r0+s0)/4) or (r0, 0, s0) with
     r0 = (k-2)/g_-, s0 = (k+2)/g_+, g_-+ = gcd(k -+ 2, |d_K|), g_+ s0 - g_- r0 = 4,
     g_- g_+ in {|d_K|, 4|d_K|}; the coincidence tables (equal D at several
     fields: same class when |Cl[2]| = 2, different classes otherwise) and the
     slices with nontrivial class that no primitive stratum reaches;
 (P) IMPRIMITIVE STRATA: the stratum of conductor f at level k carries the
     slice (D/f^2, [r_alpha O_{D/f^2}]) with the class computed exactly (HNF);
     the levels carrying D' are the solutions of k^2 - |d_K||D'| f^2 = 4 (Pell),
     and the stratum twist class along the powers of the fundamental solution;
 (U) UNITS on slices by the lattice formula (200 digits): integer palindromic
     polynomials with constant term 1, the laws, the KLF
     sum chi log|R| = -24 L'(0, chi) on odd characters (chi(r) = -1) against the
     independent incomplete-gamma evaluation, equal slices -> equal
     polynomials (D = -15 thrice, -20, -24, -36, the imprimitive carriers),
     different classes -> different systems (D = -84 thrice, D = -160 twice);
 (E) the EVEN LEVELS of iS = PSL_2(Z[i]).(i Rhat): classification of the
     second-kind orbit S_K^perp = PSL_2(O_K).(i Rhat) by BFS, the hyperbolic
     dictionary (level y, discriminant 4 - y^2|d_K|), descent, the class
     formula sigma[f] = [(n-1, n-1, n/2)][f]^{-s} at n = 4, ..., 16 with the
     lemmas A''-C'' of its proof (odd-D Gram form n N(u) - s Im(u^2)),
     the census 3H(n^2-1), and the same orbit at the other fields;
 (G) GENUS structure: the assigned characters of r_alpha are (g_+/p) at p | r0
     and (-g_-/p) at p | s0, the odd genus characters and their real fields,
     the closed forms (2h(d1)/w(d1)) h(d2) log eps_{d2} C(0) with PARI;
 (I) the ODD INDEX on slices with PARI (robert-index-full.md Thm 3 with an
     arbitrary twist): [E^- : <R_b>] = 24^{h/2} (2^{h/2-1}/Q^-)
     (h_H w_{H+})/(h_{H+} w_H) prod_{chi odd} C_chi(0);
 (F) the PHASE on slices: u_b = Phi_y/Phi_x(beta1, beta2), Phi = Phi_{r0}
     (smaller-norm twist), the j-dressing u^6 = R beta1^4 (beta1-1728)^3 /
     (beta2^4 (beta2-1728)^3), the partner law (s gives -u), the exact slice
     polynomial via Q[t]/(H_D) and its irreducibility;
 (X) OTHER CARRIERS: Lemma A' at the principal cusp of class-number-two fields
     (explicit P, no descent), and S_K^perp at the other fields.

Certification policy (CLAUDE.md guard rails): precision is set in main(),
never at import; integers/rationals are accepted only with
>= max(20, dps/5) spare digits in absolute error; forms, ideals, class groups,
characters, stratum twists and Pell solutions are exact; PARI results are
GRH-conditional unless bnfcertify returned 1; no PSLQ anywhere.

Usage:
    python3 scripts/horizontal_families.py --selftest      # everything (needs gp)
    python3 scripts/horizontal_families.py slices|realize|pell|units|even|genus|index|phase|other
Requires mpmath, sympy and PARI/GP (`gp`) for the index and genus phases.
"""
import sys
import os
import time
from math import gcd, isqrt
from fractions import Fraction
from collections import deque
from itertools import product as iproduct

from mpmath import (mp, mpf, mpc, exp, pi, log, fabs, nstr, sqrt, nint, e1,
                    matrix, det, jtheta)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import other_fields as OF
from other_fields import (Field, hnf, factor, vp, cert_int, cert_int_poly, spare_of,
                          need_digits, delta_lattice, poly_from_roots, sqrtD_ideal_basis,
                          Kp_mul, Kp_hnf, Kp_lattice_eq, epstein_Lprime0_forms,
                          form_group_data, all_characters, cval, circle_of, mat_mul,
                          mat_inv, mat_conj, mat_det, is_congruent_one, gp_run,
                          parse_tagged, parse_int_vector, parse_int_matrix,
                          smith_invariants, hurwitz_H, nearest, T_mat, S_MAT, I_MAT)
from involution_classmap import (classes_of_disc, compose, reduce_form,
                                 is_primitive as form_is_primitive)

# precision is set in main(), never here (guard rail 2)

CLASS_NUMBER_ONE = [-3, -4, -7, -8, -11, -19, -43, -67, -163]
EUCLIDEAN = [-3, -4, -7, -8, -11]
HK = {-3: 1, -4: 1, -7: 1, -8: 1, -11: 1, -19: 1, -43: 1, -67: 1, -163: 1,
      -15: 2, -20: 2, -24: 2, -35: 2, -40: 2, -51: 2, -52: 2, -88: 2, -91: 2}


class FieldX(Field):
    """other_fields.Field for any fundamental discriminant in HK."""

    def __init__(self, dK):
        assert dK < 0 and dK % 4 in (0, 1)
        self.dK = dK
        if dK % 4 == 0:
            self.t, self.n0 = 0, -dK // 4
        else:
            self.t, self.n0 = 1, (1 - dK) // 4
        self.wK = 4 if dK == -4 else (6 if dK == -3 else 2)
        self.hK = HK[dK]
        self.one = (1, 0)
        self.w = (0, 1)
        self.name = {-4: "Q(i)", -8: "Q(sqrt -2)"}.get(dK, f"Q(sqrt {dK})" if dK % 4 or dK == -4
                                                        else f"Q(sqrt {dK // 4})")


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def is_disc(D):
    return D < 0 and D % 4 in (0, 1)


def snf2(C):
    """Smith form of a 2x2 integer matrix: unimodular U, V with U C V = diag(s1, s2),
    s1 | s2, s_i >= 0.  Returns (U, (s1, s2), V)."""
    A = [list(C[0]), list(C[1])]
    U = [[1, 0], [0, 1]]
    V = [[1, 0], [0, 1]]

    def rowop(i, j, p, q, r, s):       # rows i, j <- p*row_i + q*row_j, r*row_i + s*row_j
        for M in (A, U):
            ri, rj = M[i][:], M[j][:]
            M[i] = [p * ri[k] + q * rj[k] for k in range(2)]
            M[j] = [r * ri[k] + s * rj[k] for k in range(2)]

    def colop(i, j, p, q, r, s):
        for M in (A, V):
            ci = [M[0][i], M[1][i]]
            cj = [M[0][j], M[1][j]]
            for k in range(2):
                M[k][i] = p * ci[k] + q * cj[k]
                M[k][j] = r * ci[k] + s * cj[k]

    def bezout(x, y):
        # returns (g, p, q) with p x + q y = g >= 0
        g, p, q = OF.ext_gcd(x, y)
        if g < 0:
            g, p, q = -g, -p, -q
        return g, p, q

    for _ in range(100):
        if A[0][0] == 0:
            # bring a nonzero entry to (0,0)
            if A[1][0] != 0:
                rowop(0, 1, 0, 1, 1, 0)
            elif A[0][1] != 0:
                colop(0, 1, 0, 1, 1, 0)
            elif A[1][1] != 0:
                rowop(0, 1, 0, 1, 1, 0)
                colop(0, 1, 0, 1, 1, 0)
            else:
                break
        # clear column 0
        if A[1][0] != 0:
            g, p, q = bezout(A[0][0], A[1][0])
            rowop(0, 1, p, q, -A[1][0] // g, A[0][0] // g)
            continue
        if A[0][1] != 0:
            g, p, q = bezout(A[0][0], A[0][1])
            colop(0, 1, p, q, -A[0][1] // g, A[0][0] // g)
            continue
        if A[1][1] % A[0][0] != 0:
            rowop(0, 1, 1, 1, 0, 1)
            continue
        break
    if A[0][0] < 0:
        rowop(0, 1, -1, 0, 0, 1)
    if A[1][1] < 0:
        rowop(0, 1, 1, 0, 0, -1)
    assert A[0][1] == 0 and A[1][0] == 0 and (A[0][0] == 0 or A[1][1] % A[0][0] == 0)
    # check
    for i in range(2):
        for j in range(2):
            assert sum(U[i][k] * C[k][l] * V[l][j] for k in range(2) for l in range(2)) == A[i][j]
    assert abs(U[0][0] * U[1][1] - U[0][1] * U[1][0]) == 1
    assert abs(V[0][0] * V[1][1] - V[0][1] * V[1][0]) == 1
    return U, (A[0][0], A[1][1]), V


def preimage_lattice(C, delta):
    """{u in Z^2 : C u in delta Z^2} as a list of two basis vectors (C integer 2x2)."""
    U, (s1, s2), V = snf2(C)
    m1 = delta // gcd(s1, delta) if s1 else 1
    m2 = delta // gcd(s2, delta) if s2 else 1
    b1 = (V[0][0] * m1, V[1][0] * m1)
    b2 = (V[0][1] * m2, V[1][1] * m2)
    return [b1, b2]


# ==========================================================================
# S.  slices: ambiguous forms, ideals, pairs, Cl(D)[2]
# ==========================================================================

def ambiguous_forms(D):
    """primitive ambiguous forms of disc D with their (r0, s0):
    (a, 0, c), 4ac = |D|, gcd(a, c) = 1 -> (r0, s0) = (a, c);
    (a, a, c), a(4c - a) = |D|, gcd(a, c) = 1 -> (r0, s0) = (a, 4c - a)."""
    m = -D
    out = []
    if m % 4 == 0:
        for a in divisors(m // 4):
            c = m // 4 // a
            if gcd(a, c) == 1:
                out.append(((a, 0, c), a, c))
    for a in divisors(m):
        if (m // a + a) % 4 == 0:
            c = (m // a + a) // 4
            if c > 0 and gcd(a, c) == 1:
                out.append(((a, a, c), a, m // a))
    for f, r0, s0 in out:
        assert f[1] ** 2 - 4 * f[0] * f[2] == D
    return out


def ideal_of_form(D, f):
    """Z-basis [a, (-b + sqrt D)/2] of the ideal of f = (a, b, c) (pairs x + y sqrt D)"""
    return sqrtD_ideal_basis(D, f[0], f[1])


def order_basis(D):
    """[1, (D + sqrt D)/2] (D odd) or [1, sqrt D / 2] (D even)"""
    one = (Fraction(1), Fraction(0))
    om = (Fraction(D, 2), Fraction(1, 2)) if D % 2 else (Fraction(0), Fraction(1, 2))
    return [one, om]


def lat_mul(D, L1, L2):
    return [Kp_mul(D, u, v) for u in L1 for v in L2]


def lat_conj(L):
    return [(u[0], -u[1]) for u in L]


def lat_scale(L, q):
    return [(u[0] * q, u[1] * q) for u in L]


def lat_eq(D, L1, L2):
    return Kp_hnf(D, L1) == Kp_hnf(D, L2)


def hnf_basis(D, L):
    """oriented Z-basis (w1, w2) of the lattice L, w1 rational > 0, Im(w2/w1) > 0"""
    den, g1, g2 = Kp_hnf(D, L)
    w1 = (Fraction(g1[0], den), Fraction(0))
    w2 = (Fraction(g2[0], den), Fraction(g2[1], den))
    assert w1[0] > 0 and w2[1] > 0
    return w1, w2


def lat_norm(D, L):
    """index-norm N(L) = covol(L)/covol(O_D) (Fraction)"""
    w1, w2 = hnf_basis(D, L)
    return w1[0] * w2[1] * 2        # covol(O_D) = sqrt|D|/2; covol(L) = w1 * Im(w2) = w1 * y2 sqrt|D|


def form_of_lattice(D, L):
    """the form (a, b, c) with ideal_of_form(D, (a,b,c)) = L up to the HNF basis choice
    (norm form of the oriented basis divided by N(L), with the sign convention of ideal_of_form)"""
    w1, w2 = hnf_basis(D, L)
    N = lat_norm(D, L)
    A = (w1[0] * w1[0]) / N
    B = -(2 * w1[0] * w2[0]) / N
    Cc = (w2[0] * w2[0] - D * w2[1] * w2[1]) / N
    assert A.denominator == 1 and B.denominator == 1 and Cc.denominator == 1, (D, L, A, B, Cc)
    A, B, Cc = int(A), int(B), int(Cc)
    assert B * B - 4 * A * Cc == D, (D, A, B, Cc)
    return (A, B, Cc)


def is_ideal(D, L):
    """L is an O_D-module (fractional ideal)"""
    om = order_basis(D)[1]
    return lat_eq(D, L + [Kp_mul(D, u, om) for u in L], L)


def is_ambiguous_ideal(D, L):
    return lat_eq(D, L, lat_conj(L))


def class_of_lattice(D, L):
    """reduced form of the class of the fractional ideal L (must be proper over O_D)"""
    f = form_of_lattice(D, L)
    return reduce_form(*f)


def content_in_order(D, L):
    """largest positive integer c with L subset c O_D (for an integral ideal L), and L/c"""
    w1, w2 = hnf_basis(D, L)
    one, om = order_basis(D)
    # coordinates of w1, w2 in the basis (1, om): w = x*1 + y*om  ->  y = w_y / om_y, x = w_x - y om_x
    cs = []
    for w in (w1, w2):
        y = w[1] / om[1]
        x = w[0] - y * om[0]
        assert x.denominator == 1 and y.denominator == 1, (D, L, w)
        cs += [int(x), int(y)]
    c = 0
    for v in cs:
        c = gcd(c, abs(v))
    return c, lat_scale(L, Fraction(1, c))


def gauss_mu(D):
    """number of assigned characters (Cox, Thm 3.15 for D = -4n; the odd primes for odd D)"""
    m = -D
    if D % 4 == 1:
        return len(factor(m))
    n = m // 4
    r = len([p for p in factor(n) if p != 2])
    if n % 4 == 3:
        return r
    if n % 4 in (1, 2):
        return r + 1
    if n % 8 == 4:
        return r + 1
    return r + 2


def class_group(D):
    """(prim, coords, orders, unit, mul, inv, two_torsion) for the primitive classes of disc D"""
    prim = [f for f in classes_of_disc(D) if form_is_primitive(f)]
    unit = next(f for f in prim if f[0] == 1)
    if len(prim) == 1:
        coords, orders = {unit: ()}, []
    else:
        coords, orders = form_group_data(prim, D)
    byco = {v: k for k, v in coords.items()}

    def mul(f1, f2):
        return byco[tuple((a + b) % o for a, b, o in zip(coords[f1], coords[f2], orders))]

    def inv(f):
        return byco[tuple((-a) % o for a, o in zip(coords[f], orders))]
    tt = [f for f in prim if mul(f, f) == unit]
    return prim, coords, orders, unit, mul, inv, tt


def slices_of_disc(D, check=True):
    """[(pair (r0, s0) with r0 <= s0, forms (fr, fs), class, theta-scale)] one per ambiguous class;
    verifies the slice theorem at D (D < -4)."""
    amb = ambiguous_forms(D)
    prim, coords, orders, unit, mul, inv, tt = class_group(D)
    byr0 = {}
    for f, r0, s0 in amb:
        byr0.setdefault(r0, []).append((f, s0))
    pairs = {}
    seen = set()
    for f, r0, s0 in amb:
        key = (min(r0, s0), max(r0, s0))
        if key in seen:
            continue
        seen.add(key)
        partner = [(g, t0) for (g, t0) in byr0.get(s0, []) if t0 == r0]
        assert len(partner) == 1, (D, f, r0, s0, partner)
        g = partner[0][0]
        fr, fs = (f, g) if r0 <= s0 else (g, f)
        cl = reduce_form(*fr)
        if cl not in coords:
            cl = (cl[0], -cl[1], cl[2])
        assert cl in coords, (D, fr, cl)
        cl2 = reduce_form(*fs)
        if cl2 not in coords:
            cl2 = (cl2[0], -cl2[1], cl2[2])
        assert cl2 == cl, (D, fr, fs, cl, cl2)
        if check:
            Lr, Ls = ideal_of_form(D, fr), ideal_of_form(D, fs)
            assert is_ideal(D, Lr) and is_ambiguous_ideal(D, Lr), (D, fr)
            assert is_ideal(D, Ls) and is_ambiguous_ideal(D, Ls), (D, fs)
            assert lat_norm(D, Lr) == fr[0] and lat_norm(D, Ls) == fs[0]
            # r^2 = (r0)
            assert lat_eq(D, lat_mul(D, Lr, Lr), lat_scale(order_basis(D), Fraction(fr[0]))), (D, fr)
            # r s = (theta), theta = q sqrt D with q^2 |D| = r0 s0
            q2 = Fraction(fr[0] * fs[0], -D)
            q = Fraction(isqrt(q2.numerator), isqrt(q2.denominator))
            assert q * q == q2, (D, fr, fs, q2)
            theta = (Fraction(0), q)
            assert lat_eq(D, lat_mul(D, Lr, Ls), lat_mul(D, [theta], order_basis(D))), (D, fr, fs)
            # class from the lattice agrees with the form reduction
            assert class_of_lattice(D, Lr) in (cl, (cl[0], -cl[1], cl[2])), (D, fr)
        pairs[key] = (fr, fs, cl, q if check else None)
    if check:
        # pairs <-> Cl(D)[2]
        classes = [v[2] for v in pairs.values()]
        assert len(set(classes)) == len(classes) == len(tt), (D, classes, tt)
        assert set(classes) == set(tt), (D, set(classes), set(tt))
        assert len(amb) == 2 ** gauss_mu(D), (D, len(amb), gauss_mu(D))
        assert len(pairs) == 2 ** (gauss_mu(D) - 1)
        for (r0, s0), (fr, fs, cl, q) in pairs.items():
            assert (cl == unit) == (min(r0, s0) == 1), (D, r0, s0, cl)
        # the (a, b, a)-type ambiguous forms are NOT ambiguous ideals (for D < -4)
        for f in prim:
            a, b, c = f
            if a == c and b != 0 and f not in [v[0] for v in amb] + [v[1] for v in amb]:
                L = ideal_of_form(D, f)
                assert not is_ambiguous_ideal(D, L), (D, f)
    return pairs, (prim, coords, orders, unit, mul, inv, tt)


def slices_phase(Dmax=400, say=print):
    t0 = time.time()
    say("=" * 78)
    say(f"PHASE S  slices (D, [r]) of every discriminant -{Dmax} <= D < -4")
    say("=" * 78)
    nD = 0
    nslices = 0
    nontriv = 0
    examples = {}
    for D in range(-5, -Dmax - 1, -1):
        if not is_disc(D):
            continue
        pairs, cg = slices_of_disc(D)
        nD += 1
        nslices += len(pairs)
        nontriv += sum(1 for v in pairs.values() if v[2] != cg[3])
        if D in (-15, -20, -24, -56, -84, -160, -168, -480):
            examples[D] = sorted((k, v[0], v[2]) for k, v in pairs.items())
    say(f"  {nD} discriminants, {nslices} slices ({nontriv} with nontrivial class): at every D the primitive")
    say(f"  ambiguous forms are (r0,0,s0)/(r0,r0,(r0+s0)/4), number 2^mu, their ideals are ambiguous with")
    say(f"  r^2 = (r0), the pairs {{r, s}} satisfy r s = (theta), N(theta) = r0 s0, pairs <-> Cl(D)[2],")
    say(f"  [r] = 1 iff min(r0, s0) = 1; the (a,b,a)-forms are not ambiguous ideals")
    for D, ex in examples.items():
        say(f"  D = {D}: pairs " + ", ".join(f"{k} -> class {cl}" for k, fr, cl in ex))
    say(f"  phase S done ({time.time() - t0:.1f} s)")
    return nD, nslices


# ==========================================================================
# R.  realization at primitive strata
# ==========================================================================

def slice_of_level(dK, k):
    """the slice of the level 2 alpha = k of S_K (k = 0 mod 4 at Q(i): iS):
    (D, r0, s0, g_-, g_+, fr, fs, kind) or None if (k^2 - 4) is not 0 mod |dK|"""
    d = -dK
    if (k * k - 4) % d:
        return None
    D = -(k * k - 4) // d
    gm, gp_ = gcd(k - 2, d), gcd(k + 2, d)
    r0, s0 = (k - 2) // gm, (k + 2) // gp_
    assert gp_ * s0 - gm * r0 == 4
    if gm * gp_ == d:
        assert r0 * s0 == -D and (r0 + s0) % 4 == 0, (dK, k, r0, s0, D)
        fr, fs = (r0, r0, (r0 + s0) // 4), (s0, s0, (r0 + s0) // 4)
        kind = "odd"
    else:
        assert gm * gp_ == 4 * d and 4 * r0 * s0 == -D, (dK, k, gm, gp_, r0, s0, D)
        fr, fs = (r0, 0, s0), (s0, 0, r0)
        kind = "even"
    assert fr[1] ** 2 - 4 * fr[0] * fr[2] == D and fs[1] ** 2 - 4 * fs[0] * fs[2] == D
    assert form_is_primitive(fr) and form_is_primitive(fs), (dK, k, fr, fs)
    return D, r0, s0, gm, gp_, fr, fs, kind


def twist_class(D, fr):
    """reduced class of the ambiguous form fr in the labels of classes_of_disc"""
    cl = reduce_form(*fr)
    prim = [f for f in classes_of_disc(D) if form_is_primitive(f)]
    if cl not in prim:
        cl = (cl[0], -cl[1], cl[2])
    assert cl in prim, (D, fr, cl)
    return cl


def levels_of_field(dK, kmax):
    out = []
    for k in range(3, kmax + 1):
        sl = slice_of_level(dK, k)
        if sl is None:
            continue
        D = sl[0]
        if D >= -3:
            continue          # D = -3: k = 5 at Q(sqrt -7); D = -4: the alpha = 2 levels; excluded (units)
        out.append((k, sl))
    return out


def realized_slices(fields, Dmin):
    """{D: [(dK, k, class, r0, s0)]} over all levels with |D| <= |Dmin|"""
    real = {}
    for dK in fields:
        d = -dK
        kmax = isqrt(d * (-Dmin) + 4) + 1
        for k, sl in levels_of_field(dK, kmax):
            D, r0, s0, gm, gp_, fr, fs, kind = sl
            if D < Dmin:
                continue
            cl = twist_class(D, fr)
            real.setdefault(D, []).append((dK, k, cl, r0, s0, kind))
    return real


def realize_phase(say=print, Dmin=-300):
    t0 = time.time()
    say("=" * 78)
    say(f"PHASE R  realization at primitive strata, nine class-number-one fields, |D| <= {-Dmin}")
    say("=" * 78)
    # the criterion: |dK||D| + 4 a square, and k = +-2 mod |dK| automatic; iS at Q(i) for k = 0 mod 4
    for dK in CLASS_NUMBER_ONE:
        d = -dK
        for D in range(-4, Dmin - 1, -1):
            if not is_disc(D):
                continue
            k2 = d * (-D) + 4
            k = isqrt(k2)
            has = (k * k == k2)
            sl = slice_of_level(dK, k) if has else None
            assert has == (sl is not None and sl[0] == D), (dK, D, k)
            if has:
                assert (k % d in (2, d - 2)) or (d == 4 and k % 4 == 0) or (d == 3 and k % 3 in (1, 2)), (dK, D, k)
    say("  (R1) S_K has a level of discriminant D iff |dK||D| + 4 is a perfect square k^2; k = +-2 mod |dK|")
    say("       is then automatic (at Q(i): k = 2 mod 4 gives S, k = 0 mod 4 gives iS)")
    # first slices per field (prompt's table)
    for dK in CLASS_NUMBER_ONE:
        rows = []
        for k, sl in levels_of_field(dK, 400)[:10]:
            D, r0, s0, gm, gp_, fr, fs, kind = sl
            cl = twist_class(D, fr)
            unit = (cl[0] == 1)
            rows.append(f"{k}: D={D} {fr}{'~' + str(cl) if cl != fr else ''}{' (trivial)' if unit else ''}")
        say(f"  {FieldX(dK).name:12s} " + "; ".join(rows))
    real = realized_slices(CLASS_NUMBER_ONE, Dmin)
    # the twist is trivial exactly at k = |dK| +- 2 (r0 = 1 or s0 = 1)
    for dK in CLASS_NUMBER_ONE:
        for k, sl in levels_of_field(dK, 2 * (-dK) + 10):
            D, r0, s0 = sl[:3]
            if D >= -4:
                continue
            triv = (twist_class(D, sl[5])[0] == 1)
            assert triv == (min(r0, s0) == 1) == (k in (-dK - 2, -dK + 2) or (dK == -4 and k == 6) or (dK == -3 and k in (4, 5))), (dK, k, r0, s0, triv)
    say("  (R2) the twist class is trivial exactly when min(r0, s0) = 1, i.e. at k = |dK| +- 2")
    # coincidences
    coinc = {D: v for D, v in real.items() if len(v) >= 2}
    say(f"  (R3) {len(real)} discriminants realized at primitive strata (|D| <= {-Dmin}), {len(coinc)} at >= 2 fields:")
    checks = {
        -15: {(-3, 7): (2, 1, 2), (-4, 8): (2, 1, 2), (-11, 13): (1, 1, 4), (-19, 17): (1, 1, 4)},
        -20: {(-3, 8): (2, 2, 3), (-7, 12): (2, 2, 3)},
        -24: {(-4, 10): (2, 0, 3), (-8, 14): (2, 0, 3)},
        -36: {(-7, 16): (2, 2, 5), (-11, 20): (2, 2, 5)},
        -84: {(-3, 16): (5, 4, 5), (-8, 26): (3, 0, 7), (-19, 40): (2, 2, 11)},
        -160: {(-3, 22): (7, 6, 7), (-11, 42): (4, 4, 11)},
        -180: {(-8, 38): (5, 0, 9), (-43, 88): (2, 2, 23)},
        -195: {(-4, 28): (7, 1, 7), (-7, 37): (5, 5, 11)},
        -224: {(-3, 26): (8, 8, 9), (-4, 30): (7, 0, 8)},
    }
    for D in sorted(coinc, reverse=True):
        v = coinc[D]
        _, cg = slices_of_disc(D, check=False)
        prim, coords, orders, unit, mul, inv, tt = cg
        desc = ", ".join(f"{FieldX(dK).name} k={k}: {cl}{' (trivial)' if cl == unit else ''}" for dK, k, cl, r0, s0, kind in v)
        say(f"      D = {D} (Cl = {' x '.join('Z/%d' % o for o in orders) or '1'}, |Cl[2]| = {len(tt)}): {desc}")
        if D in checks:
            for dK, k, cl, r0, s0, kind in v:
                if (dK, k) in checks[D]:
                    exp_ = reduce_form(*checks[D][(dK, k)])
                    if exp_ not in prim:
                        exp_ = (exp_[0], -exp_[1], exp_[2])
                    assert cl == exp_, (D, dK, k, cl, exp_)
        if len(tt) == 2:
            cls = set(cl for dK, k, cl, r0, s0, kind in v if cl != unit)
            assert len(cls) <= 1, (D, cls)
    say("      (the predicted classes of the prompt at D = -15, -20, -24, -36, -84, -160, -180, -195, -224 all confirmed;")
    say("       with |Cl[2]| = 2 all nontrivial realizations of a D agree)")
    # unrealized nontrivial slices |D| <= 120
    unreal = []
    for D in range(-5, -121, -1):
        if not is_disc(D):
            continue
        pairs, cg = slices_of_disc(D, check=False)
        unit = cg[3]
        got = set(cl for dK, k, cl, r0, s0, kind in real.get(D, []))
        for key, (fr, fs, cl, q) in pairs.items():
            if cl != unit and cl not in got:
                unreal.append((D, cl, key))
    say(f"  (R4) slices with nontrivial class realized by no primitive stratum, |D| <= 120 ({len(unreal)} slices): " +
        ", ".join(f"{D} {cl}" for D, cl, key in unreal))
    none_at_all = sorted(set(D for D, cl, key in unreal if D not in real), reverse=True)
    say(f"       discriminants with a nontrivial class and no realization at all: {none_at_all} (the prompt's list);"
        f" D = -96 and -120 are realized at one class only")
    assert none_at_all == [-56, -72, -88, -91, -100, -104, -115, -116], none_at_all
    assert sorted(set(D for D, cl, key in unreal), reverse=True) == [-56, -72, -88, -91, -96, -100, -104, -115, -116, -120]
    # (R5) the t-parametrization: k = t|dK| +- 2 (odd dK) has the pair {t, t|dK| +- 4}; k = 4mt +- 2 (dK = -4m) has {t, mt +- 1};
    #      conversely the fields realizing a pair {r0, s0} are read off from (s0 -+ 4)/r0 resp. (s0 -+ 1)/r0
    for dK in CLASS_NUMBER_ONE:
        d = -dK
        for k, sl in levels_of_field(dK, 40 * d):
            D, r0, s0, gm, gp_, fr, fs, kind = sl
            if D >= -4:
                continue
            cands = []
            if d % 2:
                if (k - 2) % d == 0:
                    t = (k - 2) // d
                    cands.append({t, t * d + 4})
                if (k + 2) % d == 0:
                    t = (k + 2) // d
                    cands.append({t, t * d - 4})
            elif d == 4 and k % 4 == 0:
                t = k // 4                          # iS: k = 4t, pair {2t-1, 2t+1}
                cands.append({2 * t - 1, 2 * t + 1})
            else:
                m = d // 4
                if (k - 2) % (4 * m) == 0:
                    t = (k - 2) // (4 * m)
                    cands.append({t, m * t + 1})
                if (k + 2) % (4 * m) == 0:
                    t = (k + 2) // (4 * m)
                    cands.append({t, m * t - 1})
            assert {r0, s0} in cands, (dK, k, r0, s0, cands)
    # every pair {r0, s0} with |D| <= |Dmin|: the realizing fields are exactly those predicted
    nchk = 0
    for D in range(-5, Dmin - 1, -1):
        if not is_disc(D):
            continue
        pairs, cg = slices_of_disc(D, check=False)
        got = {}
        for dK, k, cl, r0, s0, kind in real.get(D, []):
            got.setdefault(cl, set()).add(dK)
        for (r0, s0), (fr, fs, cl, q) in pairs.items():
            predf = set()
            if r0 * s0 == -D:                       # odd type D = -r0 s0: odd d_K, and iS
                for (a, b) in ((r0, s0), (s0, r0)):
                    for sg in (1, -1):
                        qq = b - sg * 4
                        if qq > 0 and qq % a == 0 and (qq // a) % 2 == 1 and -(qq // a) in CLASS_NUMBER_ONE:
                            predf.add(-(qq // a))
                if abs(r0 - s0) == 2 and r0 % 2 == 1:
                    predf.add(-4)
            else:                                   # even type D = -4 r0 s0: Q(i) (S) and Q(sqrt -2)
                for (a, b) in ((r0, s0), (s0, r0)):
                    for sg in (1, -1):
                        qq = b - sg
                        if qq > 0 and qq % a == 0 and qq // a in (1, 2):
                            predf.add(-4 * (qq // a))
            assert predf == got.get(cl, set()), (D, r0, s0, cl, predf, got.get(cl))
            nchk += 1
    say(f"  (R5) k = t|dK| +- 2 carries the pair {{t, t|dK| +- 4}} of type D = -r0 s0 (odd dK); k = 4mt +- 2 the pair {{t, mt +- 1}} of type")
    say(f"       D = -4 r0 s0 (dK = -4m); k = 4t (iS) the pair {{2t-1, 2t+1}}: twists of norm t.  Conversely a slice with pair {{r0, s0}}")
    say(f"       is realized iff |dK| = (s0 -+ 4)/r0 or (r0 -+ 4)/s0 (odd type, odd dK), |r0 - s0| = 2 (odd type, iS),")
    say(f"       m = (s0 -+ 1)/r0 or (r0 -+ 1)/s0 in {{1, 2}} (even type, Q(i) and Q(sqrt -2)): {nchk} pairs checked, |D| <= {-Dmin}")
    say(f"  phase R done ({time.time() - t0:.1f} s)")
    return real, unreal


# ==========================================================================
# P.  imprimitive strata: the Pell equation and the stratum twist class
# ==========================================================================

def extend_twist(D, fr, f):
    """r_alpha O_{D'} for D = f^2 D': (content c, primitive form fr' of the ambiguous ideal r' = r O'/c,
    its class in Cl(D')) -- exact HNF arithmetic in Q(sqrt D')."""
    assert D % (f * f) == 0
    Dp = D // (f * f)
    assert is_disc(Dp)
    # r as a lattice in Q(sqrt D'): sqrt D = f sqrt D'
    r0, br = fr[0], fr[1]
    L = [(Fraction(r0), Fraction(0)), (Fraction(-br, 2), Fraction(f, 2))]
    Op = order_basis(Dp)
    rOp = lat_mul(Dp, L, Op)
    assert is_ideal(Dp, rOp) and is_ambiguous_ideal(Dp, rOp), (D, fr, f)
    assert lat_norm(Dp, rOp) == r0, (D, fr, f, lat_norm(Dp, rOp))
    c, rp = content_in_order(Dp, rOp)
    assert lat_norm(Dp, rp) * c * c == r0
    frp = form_of_lattice(Dp, rp)
    assert form_is_primitive(frp) and frp[0] * c * c == r0, (D, fr, f, frp, c)
    assert frp[1] % frp[0] == 0, ("primitive part not an ambiguous form", D, fr, f, frp)
    # (r')^2 = (r'_0)
    assert lat_eq(Dp, lat_mul(Dp, rp, rp), lat_scale(Op, Fraction(frp[0]))), (D, fr, f)
    cl = twist_class(Dp, frp)
    return c, frp, cl


def fundamental_pell(m):
    """the fundamental solution (k, f) of k^2 - m f^2 = 4, k > 2, f > 0 (m > 0 not a square),
    through PARI's quadunit of the order of discriminant m (m = 0, 1 mod 4) or 4m."""
    assert m > 0 and isqrt(m) ** 2 != m
    if OF.GP is None:
        f = 1
        while True:
            k2 = m * f * f + 4
            k = isqrt(k2)
            if k * k == k2:
                return k, f
            f += 1
            assert f < 10 ** 6, m
    disc = m if m % 4 in (0, 1) else 4 * m
    out = gp_run(f"""u = quadunit({disc}); print("A ", component(u, 2)); print("B ", component(u, 3)); print("N ", norm(u));""")
    d = parse_tagged(out)
    a, b, nrm = int(d["A"]), int(d["B"]), int(d["N"])
    # u = a + b w, w = (1 + sqrt disc)/2 (disc odd) or sqrt(disc)/2 = sqrt(disc/4) (disc even)
    if disc % 2:
        x, y = Fraction(2 * a + b, 2), Fraction(b, 2)          # u = x + y sqrt(disc)
    else:
        x, y = Fraction(a), Fraction(b, 2)
    # convert to sqrt m: sqrt(disc) = sqrt(m) (disc = m) or 2 sqrt(m) (disc = 4m)
    if disc == 4 * m:
        y = 2 * y
    assert x * x - m * y * y == nrm
    if nrm == -1:
        x, y = x * x + m * y * y, 2 * x * y
    assert x * x - m * y * y == 1
    k, f = 2 * x, 2 * y
    assert k.denominator == 1 and f.denominator == 1
    k, f = int(k), int(f)
    if k < 0:
        k, f = -k, -f
    if f < 0:
        f = -f
    assert k * k - m * f * f == 4 and k > 2 and f > 0
    return k, f


def pell_powers(m, k1, f1, J):
    """(k_j, f_j) for ((k1 + f1 sqrt m)/2)^j, j = 1..J"""
    out = [(k1, f1)]
    k, f = k1, f1
    for _ in range(J - 1):
        k, f = (k * k1 + m * f * f1) // 2, (k * f1 + f * k1) // 2
        assert k * k - m * f * f == 4
        out.append((k, f))
    return out


def stratum_law(dK, Dp, J=4):
    """the stratum twist classes of D' along the powers of the fundamental solution of
    k^2 - |dK||D'| f^2 = 4: [(j, k, f, content, r'_0, class)]"""
    m = (-dK) * (-Dp)
    if isqrt(m) ** 2 == m:
        return None
    k1, f1 = fundamental_pell(m)
    out = []
    for j, (k, f) in enumerate(pell_powers(m, k1, f1, J), start=1):
        sl = slice_of_level(dK, k)
        assert sl is not None
        D, r0, s0, gm, gp_, fr, fs, kind = sl
        assert D == f * f * Dp, (dK, Dp, j, k, f, D)
        c, frp, cl = extend_twist(D, fr, f)
        # the partner s gives the same class (rational multiple of the same pair)
        c2, fsp, cl2 = extend_twist(D, fs, f)
        assert cl2 == cl, (dK, Dp, j, cl, cl2)
        out.append((j, k, f, c, frp[0], cl, frp))
    return out


def pell_phase(say=print, J=4, Dlist=None, fields=None):
    t0 = time.time()
    say("=" * 78)
    say("PHASE P  imprimitive strata: the Pell equation k^2 - |dK||D'| f^2 = 4 and the stratum twist class")
    say("=" * 78)
    # (P1) the instances of the prompt, exact
    inst = [(-3, 26, 2, -56, (2, 0, 7)), (-8, 22, 2, -15, (2, 1, 2)), (-4, 18, 2, -20, (1, 0, 5)),
            (-4, 34, 2, -72, (2, 0, 9)), (-19, 74, 2, -72, (1, 0, 18)), (-7, 54, 2, -104, (2, 0, 13)),
            (-4, 102, 5, -104, (1, 0, 26)), (-3, 52, 3, -100, (2, 2, 13)), (-8, 198, 7, -100, (1, 0, 25)),
            (-4, 30, 2, -56, (2, 0, 7))]
    for dK, k, f, Dp, exp_cl in inst:
        sl = slice_of_level(dK, k)
        D, r0, s0, gm, gp_, fr, fs, kind = sl
        assert D == f * f * Dp, (dK, k, D, f, Dp)
        c, frp, cl = extend_twist(D, fr, f)
        exp_cl = twist_class(Dp, exp_cl)
        assert cl == exp_cl, (dK, k, f, Dp, cl, exp_cl)
        say(f"  (P1) {FieldX(dK).name} k = {k}: D = {D} = {f}^2 ({Dp}), twist {fr} (norm {r0}) extends to "
            f"{c} * {frp}, class {cl}{' (trivial)' if cl[0] == 1 else ''}")
    # (P2) the law along the Pell solutions
    fields = CLASS_NUMBER_ONE if fields is None else fields
    Dlist = [-15, -20, -24, -36, -56, -72, -84, -88, -91, -100, -104, -115, -116, -160] if Dlist is None else Dlist
    table = {}
    for Dp in Dlist:
        pairs, cg = slices_of_disc(Dp, check=False)
        unit = cg[3]
        for dK in fields:
            law = stratum_law(dK, Dp, J=J)
            if law is None:
                say(f"  (P2) {FieldX(dK).name}, D' = {Dp}: |dK||D'| is a square (K = Q(sqrt D')): no level, the vertical (Euclidean) case")
                continue
            table[(dK, Dp)] = law
            desc = "; ".join(f"j={j}: k={k if k < 10 ** 12 else str(k)[:8] + '...(' + str(len(str(k))) + ' digits)'}, f={f if f < 10 ** 12 else str(f)[:6] + '...'}, "
                             f"c={c}, r'0={rp0}, class {cl}{'=1' if cl == unit else ''}" for j, k, f, c, rp0, cl, frp in law)
            say(f"  (P2) {FieldX(dK).name:12s} D' = {Dp:5d} (|Cl[2]| = {len(cg[6])}): {desc}")
    # (P3) the law: class at j depends only on the parity of j, and the odd-j class is the 'genus class'
    ok_parity = True
    for (dK, Dp), law in table.items():
        cls_odd = set(cl for j, k, f, c, rp0, cl, frp in law if j % 2 == 1)
        cls_even = set(cl for j, k, f, c, rp0, cl, frp in law if j % 2 == 0)
        if len(cls_odd) > 1 or len(cls_even) > 1:
            ok_parity = False
            say(f"      !! {FieldX(dK).name} D' = {Dp}: the class is NOT a function of the parity of j: odd {cls_odd}, even {cls_even}")
    say(f"  (P3) the stratum twist class is a function of the parity of j at every pair: {ok_parity}")
    # (P4) the law: the odd part of r'_0 is prod_{p odd, v_p(k-2) > v_p(k+2)} p^{v_p(D')}; class(j) = c^j
    #      (trivial for even j) with c the class at the fundamental solution; for the pairs of type
    #      (r'_0, 0, s'_0) the same rule holds at p = 2 with the exponent v_2(|D'|/4)
    ok_law = True
    ok_two = True
    for (dK, Dp), law in table.items():
        cls_odd = None
        for j, k, f, c, rp0, cl, frp in law:
            pred = 1
            for p, e in factor(-Dp).items():
                if p != 2 and vp(k - 2, p) > vp(k + 2, p):
                    pred *= p ** e
            odd_part = rp0
            while odd_part % 2 == 0:
                odd_part //= 2
            if pred != odd_part:
                ok_law = False
                say(f"      !! odd-part law fails at {FieldX(dK).name}, D' = {Dp}, j = {j}: predicted {pred}, got {odd_part}")
            if frp[1] == 0:
                e2 = vp(-Dp // 4, 2)
                pred2 = 2 ** e2 if vp(k - 2, 2) > vp(k + 2, 2) else 1
                if pred2 != rp0 // odd_part:
                    ok_two = False
                    say(f"      !! 2-part law fails at {FieldX(dK).name}, D' = {Dp}, j = {j}: predicted {pred2}, got {rp0 // odd_part}")
            if j % 2 == 0 and cl[0] != 1:
                ok_law = False
            if j % 2 == 1:
                cls_odd = cl if cls_odd is None else cls_odd
                if cl != cls_odd:
                    ok_law = False
    say(f"  (P4) odd part of r'_0 = prod_{{p odd : v_p(k-2) > v_p(k+2)}} p^{{v_p(D')}} and class(j) = c^j (c = the class at j = 1, trivial at even j)"
        f" at every pair: {ok_law}; the 2-part obeys the same rule (exponent v_2(|D'|/4)) on the pairs of type (r'_0, 0, s'_0): {ok_two}")
    # (P5) which nontrivial slices |D'| <= 120 (and -160) are realized by the strata of the nine fields
    realized_strata = {}
    for Dp in sorted(set(Dlist), reverse=True):
        pairs, cg = slices_of_disc(Dp, check=False)
        unit = cg[3]
        got = {}
        for dK in fields:
            law = table.get((dK, Dp))
            if law:
                c1 = law[0][5]
                if c1 != unit:
                    got.setdefault(c1, []).append((FieldX(dK).name, law[0][1], law[0][2]))
        missing = [cl for key, (fr, fs, cl, q) in pairs.items() if cl != unit and cl not in got]
        realized_strata[Dp] = (got, missing)
        say(f"  (P5) D' = {Dp}: " + "; ".join(f"{cl} at {', '.join(f'{nm} k={k} (f={f})' for nm, k, f in v)}" for cl, v in got.items())
            + (f"; NOT realized by any stratum of the nine fields: {missing}" if missing else "; every nontrivial class realized"))
    say(f"  phase P done ({time.time() - t0:.1f} s)")
    return table, realized_strata


# ==========================================================================
# U.  units on slices
# ==========================================================================

def R_slice(D, fr, f):
    """R_f = r0^6 Delta(b_f)/Delta(r^{-1} b_f) for the slice (D, [fr]), fr an ambiguous form of norm r0
    (f primitive or not: on an imprimitive class the twist is r O_{D'})"""
    return OF.R_unit(None, None, None, f, D, fr[0], fr)


def slice_units(D, fr, say=None, with_klf=True, forms=None, tag=""):
    """unit polynomial of the slice (D, [fr]) over the primitive classes, the laws, the KLF on odd characters;
    returns dict(poly, R, odd, rn, prim, coords, orders, Lp)"""
    need = need_digits()
    prim, coords, orders, unit, mul, inv, tt = class_group(D)
    rn = twist_class(D, fr)
    R = {f: R_slice(D, fr, f) for f in prim}
    co, sp = cert_int_poly([R[f] for f in prim], f"R-polynomial {tag} D={D} r={fr}")
    assert abs(co[-1]) == 1, (D, fr, co)
    assert co == co[::-1] or len(prim) == 1, (D, fr, co)
    worst = mp.dps
    for f in prim:
        worst = min(worst, spare_of(R[mul(rn, f)] * R[f] - 1))
        worst = min(worst, spare_of(R[inv(f)] - R[f].conjugate()))
    assert worst >= need, (D, fr, worst)
    out = {"poly": co, "R": R, "rn": rn, "prim": prim, "coords": coords, "orders": orders, "spare": sp,
           "h": len(prim), "unit": unit, "mul": mul, "inv": inv, "tt": tt}
    if with_klf and len(prim) > 1:
        chars = all_characters(coords, orders)
        chis, labels = [], []
        for (ks, cord, ph) in chars:
            if cord == 1:
                continue
            chis.append({f: cval(ph[f]) for f in prim})
            labels.append((ks, cord, ph))
        Lp, M = epstein_Lprime0_forms(prim, D, chis)
        odd, even = [], []
        wk = mp.dps
        for (ks, cord, ph), chi, l in zip(labels, chis, Lp):
            S = sum(chi[f] * log(fabs(R[f])) for f in prim)
            if ph[rn] == Fraction(1, 2):
                odd.append((ks, cord, l))
                wk = min(wk, spare_of((S + 24 * l) / max(1, fabs(l))))
            else:
                even.append((ks, cord))
                wk = min(wk, spare_of(S))
        assert wk >= need, (D, fr, wk)
        assert len(odd) == (len(prim) // 2 if rn != unit else 0), (D, fr, len(odd))
        out.update({"odd": odd, "even": even, "Lp": Lp, "klf_spare": wk, "M": M})
    if say:
        say(f"  {tag} D = {D}, [r] = {rn} (r0 = {fr[0]}), h = {len(prim)}: prod (x - R) = "
            f"{co if len(co) <= 5 else str(co[:3])[:-1] + ', ..., ' + str(co[-1]) + ']'} (spare {sp}); laws (spare {worst})"
            + (f"; KLF: {len(out['odd'])} odd characters, sum chi log|R| = -24 L'(0,chi), even sums 0 (spare {out['klf_spare']})" if "odd" in out else ""))
    return out


def stratum_values(D, fr, f):
    """{primitive form f' of D' = D/f^2 : R on the class f * f' of D with the twist r_alpha}"""
    Dp = D // (f * f)
    out = {}
    for fp in classes_of_disc(Dp):
        if not form_is_primitive(fp):
            continue
        g = (f * fp[0], f * fp[1], f * fp[2])
        out[fp] = R_slice(D, fr, g)
    return out


UNIT_TARGETS = {
    # (dK, k): slices at the primitive strata, first levels of each field (D >= -200)
    -3: [7, 8, 10, 11, 13, 14, 16, 17, 19, 20, 22, 23],
    -4: [8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28],
    -7: [12, 16, 19, 23, 26, 30, 33, 37],
    -8: [14, 18, 22, 26, 30, 34, 38],
    -11: [20, 24, 31, 35, 42, 46],
    -19: [36, 40, 55, 59],
    -43: [84, 88],
    -67: [132, 136],
    -163: [324, 328],
}


# regression record: (D, twist class) -> unit polynomial (leading first)
UNIT_RECORD = {
    (-15, (2, 1, 2)): [1, -7, 1], (-20, (2, 2, 3)): [1, 18, 1], (-24, (2, 0, 3)): [1, -34, 1], (-35, (3, 1, 3)): [1, -322, 1],
    (-36, (2, 2, 5)): [1, 194, 1], (-63, (4, 1, 4)): [1, -31279, 45681, -31279, 1],
    (-84, (5, 4, 5)): [1, 686308, 7365318, 686308, 1], (-84, (3, 0, 7)): [1, -297220, 7312902, -297220, 1],
    (-84, (2, 2, 11)): [1, 27940, 76614, 27940, 1], (-99, (5, 1, 5)): [1, -4468994, 1],
    (-160, (7, 6, 7)): [1, 5677418956, 13748395686, 5677418956, 1], (-160, (4, 4, 11)): [1, 2146623116, 4535235366, 2146623116, 1],
    (-195, (7, 1, 7)): [1, -174683640964, 2849307736326, -174683640964, 1], (-195, (5, 5, 11)): [1, -112143790724, 2842988279046, -112143790724, 1],
    (-56, (2, 0, 7)): [1, -1988, -3194, -1988, 1], (-72, (2, 0, 9)): [1, -9602, 1], (-100, (2, 2, 13)): [1, 103682, 1],
}


def units_phase(say=print, targets=None, extra=True):
    t0 = time.time()
    say("=" * 78)
    say(f"PHASE U  units on slices by the lattice formula (dps = {mp.dps})")
    say("=" * 78)
    targets = UNIT_TARGETS if targets is None else targets
    systems = {}        # (D, class) -> (poly, realization list, odd labels)
    for dK in targets:
        for k in targets[dK]:
            sl = slice_of_level(dK, k)
            D, r0, s0, gm, gp_, fr, fs, kind = sl
            tag = f"{FieldX(dK).name} k={k}{' (iS)' if dK == -4 and k % 4 == 0 else ''}:"
            res = slice_units(D, fr, say=say, tag=tag)
            # independence of the partner s (same values)
            Rs = {f: R_slice(D, fs, f) for f in res["prim"]}
            assert min(spare_of(Rs[f] - res["R"][f]) for f in res["prim"]) >= need_digits()
            key = (D, res["rn"])
            if key in UNIT_RECORD:
                assert res["poly"] == UNIT_RECORD[key], (key, res["poly"])
            systems.setdefault(key, []).append((dK, k, res["poly"], sorted(res.get("odd", []))))
    # (U1) equal slices -> equal polynomials and odd characters; different classes -> different
    byD = {}
    for (D, cl), lst in systems.items():
        byD.setdefault(D, {})[cl] = lst
    for D in sorted(byD, reverse=True):
        for cl, lst in byD[D].items():
            polys = set(tuple(p) for dK, k, p, odd in lst)
            assert len(polys) == 1, (D, cl, polys)
            odds = set(tuple(odd) for dK, k, p, odd in lst)
            assert len(odds) == 1, (D, cl, odds)
        if len(byD[D]) > 1:
            polys = [tuple(lst[0][2]) for lst in byD[D].values()]
            nontriv = [tuple(lst[0][2]) for cl, lst in byD[D].items() if cl[0] != 1]
            assert len(set(nontriv)) == len(nontriv), (D, "different classes with equal polynomials", nontriv)
            odds = [tuple(lst[0][3]) for cl, lst in byD[D].items() if cl[0] != 1]
            assert len(set(odds)) == len(odds), (D, "different classes with equal odd characters")
            say(f"  (U1) D = {D}: {len(byD[D])} classes realized -> " + "; ".join(
                f"{cl}: {[(FieldX(dK).name, k) for dK, k, p, odd in lst]}" for cl, lst in byD[D].items())
                + " -- distinct unit polynomials and distinct odd character sets")
        elif len(list(byD[D].values())[0]) > 1:
            lst = list(byD[D].values())[0]
            say(f"  (U1) D = {D}, class {list(byD[D])[0]}: same polynomial at {[(FieldX(dK).name, k) for dK, k, p, odd in lst]}")
    if extra:
        # (U2) imprimitive carriers: stratum values = slice values classwise
        carriers = [(-8, 22, 2, (-3, 7)), (-3, 26, 2, (-4, 30)), (-4, 30, 2, None), (-4, 18, 2, None),
                    (-4, 34, 2, None), (-7, 54, 2, None), (-3, 52, 3, None)]
        for dK, k, f, other in carriers:
            sl = slice_of_level(dK, k)
            D, r0, s0, gm, gp_, fr, fs, kind = sl
            Dp = D // (f * f)
            c, frp, cl = extend_twist(D, fr, f)
            strat = stratum_values(D, fr, f)
            direct = {fp: R_slice(Dp, frp, fp) for fp in strat}
            worst = min(spare_of(strat[fp] - direct[fp]) for fp in strat)
            assert worst >= need_digits(), (dK, k, f, worst)
            co, sp = cert_int_poly(list(strat.values()), f"stratum {dK} {k} {f}")
            assert abs(co[-1]) == 1
            if (Dp, cl) in UNIT_RECORD:
                assert co == UNIT_RECORD[(Dp, cl)], (Dp, cl, co)
            note = ""
            if other is not None:
                sl2 = slice_of_level(other[0], other[1])
                D2, fr2 = sl2[0], sl2[5]
                if D2 == Dp:
                    res2 = slice_units(Dp, fr2, with_klf=False)
                    assert res2["poly"] == co, (dK, k, f, co, res2["poly"])
                    note = f"; = the primitive-stratum polynomial of {FieldX(other[0]).name} k={other[1]}"
                else:
                    c2, frp2, cl2 = extend_twist(D2, fr2, isqrt(D2 // Dp))
                    strat2 = stratum_values(D2, fr2, isqrt(D2 // Dp))
                    co2, sp2 = cert_int_poly(list(strat2.values()), "stratum 2")
                    assert co2 == co and cl2 == cl, (dK, k, f, other, co, co2, cl, cl2)
                    note = f"; = the conductor-{isqrt(D2 // Dp)} stratum of {FieldX(other[0]).name} k={other[1]}"
            say(f"  (U2) {FieldX(dK).name} k={k}, conductor-{f} stratum (D = {D} = {f}^2 ({Dp})): twist extends to {c} * {frp}, "
                f"class {cl}{' (trivial)' if cl[0] == 1 else ''}; R on the stratum = R of the slice ({Dp}, {cl}) classwise (spare {worst}); "
                f"polynomial {co if len(co) <= 5 else str(co[:2])[:-1] + ', ...]'} (spare {sp}){note}")
    say(f"  phase U done ({time.time() - t0:.1f} s)")
    return systems


# ==========================================================================
# E/X.  the second-kind orbit S_K^perp = PSL_2(O_K).(i Rhat)  (= iS at Q(i))
# ==========================================================================
# Circles X(i Rhat), X in SL_2(O_K), as Hermitian matrices [[A, B],[conj B, C]] (det -1):
#   A = -Tr(c conj d), B = a conj d + b conj c, C = -Tr(a conj b),  N(B) = 1 + A C.
# Translation by lam:  (A, B, C) -> (A, B - lam A, C - Tr(conj(lam) B) + N(lam) A);
# inversion z -> -1/z: (A, B, C) -> (C, -conj B, A).
# In the normalization of S_K ([[sqrt|d| q, i beta], ...]) this is beta = i B, q = -A/sqrt|d|:
# S_K^perp is the orbit of the second unimodular Hermitian form h_1 = Tr(x conj y) instead of
# h_0 = 2 Im(x conj y); at Q(i) the two forms are isometric (diag(i,1)) and S^perp = iS.

def tr(K, x):
    return 2 * x[0] + x[1] * K.t


def circle_perp(K, X):
    (a, b), (c, d) = X
    A = -tr(K, K.mul(c, K.conj(d)))
    B = K.add(K.mul(a, K.conj(d)), K.mul(b, K.conj(c)))
    C = -tr(K, K.mul(a, K.conj(b)))
    assert K.norm(B) == 1 + A * C, (X, A, B, C)
    return (A, B, C)


def translate_perp(K, circ, lam):
    A, B, C = circ
    return (A, K.sub(B, K.smul(A, lam)), C - tr(K, K.mul(K.conj(lam), B)) + K.norm(lam) * A)


def invert_perp(K, circ):
    A, B, C = circ
    return (C, K.neg(K.conj(B)), A)


def descend_perp(K, circ):
    """an X in SL_2(O_K) with circle_perp(K, X) == circ (Euclidean fields).  Terminal circles:
    the lines (0, u, C), u a unit, and the unit circle +-(1, 0, -1) (odd d_K)."""
    g = I_MAT
    cur = circ
    steps = 0
    while cur[0] != 0 and not (abs(cur[0]) == 1 and cur[1] == (0, 0)):
        A, B, C = cur
        lam = nearest(K, B, A) if A > 0 else nearest(K, K.neg(B), -A)
        cur = translate_perp(K, cur, lam)
        g = mat_mul(K, T_mat(K, lam), g)
        cur = invert_perp(K, cur)
        g = mat_mul(K, S_MAT, g)
        steps += 1
        assert steps < 500, "descent does not terminate"
    A0, B0, C0 = cur
    Y = None
    if A0 == 0:
        assert K.norm(B0) == 1
        us = [u for u in K.units() if K.mul(u, u) == B0]
        us2 = [u for u in K.units() if K.mul(u, u) == K.neg(B0)]
        assert us or us2, ("terminal B not +- a unit square", K.dK, cur)
        # Y = [[u, b], [0, conj u]] (B = u^2, C = -Tr(u conj b)) or [[b, -u], [conj u, 0]] (B = -u^2, C = Tr(b conj u))
        for x in range(-abs(C0) - 2, abs(C0) + 3):
            for yy in range(-abs(C0) - 2, abs(C0) + 3):
                b = (x, yy)
                if us and -tr(K, K.mul(us[0], K.conj(b))) == C0:
                    Y = ((us[0], b), ((0, 0), K.conj(us[0])))
                    break
                if us2 and tr(K, K.mul(b, K.conj(us2[0]))) == C0:
                    Y = ((b, K.neg(us2[0])), (K.conj(us2[0]), (0, 0)))
                    break
            if Y:
                break
    else:
        # the unit circle (1, 0, -1): X0 = [[-1, -conj w], [1, -w]] (Tr w = 1); (-1, 0, 1) = S X0
        assert K.t == 1, ("unit circle at even d_K", K.dK, cur)
        X0 = ((( -1, 0), K.neg(K.conj(K.w))), ((1, 0), K.neg(K.w)))
        assert circle_perp(K, X0) == (1, (0, 0), -1)
        Y = X0 if A0 == 1 else mat_mul(K, S_MAT, X0)
    assert Y is not None, ("terminal circle not realized", K.dK, cur)
    assert circle_perp(K, Y) == cur
    X = mat_mul(K, mat_inv(K, g), Y)
    assert mat_det(K, X) == K.one
    assert circle_perp(K, X) == circ, (circ, circle_perp(K, X))
    return X


def orbit_bfs_perp(K, Q, NB):
    start = (0, (1, 0), 0)
    seen = {start}
    dq = deque([start])
    gens = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    while dq:
        c = dq.popleft()
        nbrs = [translate_perp(K, c, lam) for lam in gens] + [invert_perp(K, c)]
        for nb in nbrs:
            A2, B2, C2 = nb
            if abs(A2) > Q or abs(C2) > Q or K.norm(B2) > NB:
                continue
            if nb not in seen:
                seen.add(nb)
                dq.append(nb)
    return seen


def residues_mod_A(K, A):
    return OF.hnf2_general((A, 0), (0, A))


def perp_congruence(K, A, B):
    """membership of (A, B, (N(B)-1)/A) in S_K^perp (Theorem on S_K^perp): A, C in Tr(O_K),
    N(B) = 1 + A C, and B = 1 mod 2 at K = Q(i) (the orbit of R hat has B = i mod 2 there)"""
    two = 2 if K.t == 0 else 1
    if A % two:
        return False
    if (K.norm(B) - 1) % A:
        return False
    C = (K.norm(B) - 1) // A
    if C % two:
        return False
    if K.dK == -4 and not (B[0] % 2 == 1 and B[1] % 2 == 0):
        return False
    return True


def classification_perp(K, Qmax=6, say=print):
    """BFS orbit of i Rhat; for every |A| <= Qmax the residues B mod A O_K that occur are exactly the
    predicted ones; returns the level set {y: signs} (y = -sgn(A) im_K(B) > 0, alpha = y sqrt|d_K|/2)"""
    d = -K.dK
    orb = orbit_bfs_perp(K, 2 * Qmax + 6, 4000 * d)
    counts = {}
    for A in [a for a in range(1, Qmax + 1)] + [-a for a in range(1, Qmax + 1)]:
        H = residues_mod_A(K, abs(A))
        got = set(OF.reduce_mod(B, H) for (AA, B, C) in orb if AA == A)
        g0, (p, r) = H
        pred = set()
        for x in range(g0):
            for y in range(r):
                B = (x, y)
                if perp_congruence(K, A, B):
                    pred.add(OF.reduce_mod(B, H))
        assert got == pred, (K.dK, A, sorted(got - pred)[:5], sorted(pred - got)[:5], len(got), len(pred))
        counts[A] = len(got)
    levels = {}
    for (A, B, C) in orb:
        if A == 0:
            continue
        y = -K.imK(B) * (1 if A > 0 else -1)
        if y > 0 and y * y * d > 4:
            levels.setdefault(y, set()).add(1 if A > 0 else -1)
    say(f"  (X1) {K.name}: orbit of i Rhat (BFS, {len(orb)} circles) = {{A, C in Tr(O_K), N(B) = 1 + AC"
        f"{', B = 1 mod 2' if K.dK == -4 else ''}}} for every |A| <= {Qmax}: residues per A = 1..{Qmax}: "
        f"{[counts[a] for a in range(1, Qmax + 1)]} (same for A < 0); levels y with both orientations: "
        f"{sorted(y for y, sg in levels.items() if sg == {1, -1})[:8]}")
    return levels


def disc_perp(K, y):
    """discriminant of the level-y forms of S_K^perp, and the scale (1 or 2) dividing (A, 2 Re B, C)"""
    d = -K.dK
    if K.t == 1:
        return 4 - y * y * d, 1
    return (4 - y * y * d) // 4, 2


def form_perp(K, circ):
    """(form (a, b, c), level y, orientation sign) of a circle of S_K^perp in H"""
    A, B, C = circ
    s = 1 if A > 0 else -1
    X = 2 * B[0] + B[1] * K.t          # 2 Re B
    y = -s * K.imK(B)
    assert y > 0
    sc = 1 if K.t == 1 else 2
    assert s * A % sc == 0 and s * X % sc == 0 and s * C % sc == 0
    f = (s * A // sc, s * X // sc, s * C // sc)
    assert f[1] ** 2 - 4 * f[0] * f[2] == disc_perp(K, y)[0], (circ, f, disc_perp(K, y))
    return f, y, s


def circle_perp_of_form(K, f, y, s):
    """the level-y circle of the form f with orientation sign s (in the R hat frame)"""
    a, b, c = f
    sc = 1 if K.t == 1 else 2
    A, C = s * sc * a, s * sc * c
    X = s * sc * b                    # = 2 Re B
    B2 = -s * y                       # im_K(B)
    assert (X - B2 * K.t) % 2 == 0
    B = ((X - B2 * K.t) // 2, B2)
    circ = (A, B, C)
    assert K.norm(B) == 1 + A * C, (f, y, s, circ)
    return circ


def perp_level_of_tau(K, X):
    """the R hat level (y, or None for a line / a circle through R hat) of tau(X)(i Rhat),
    tau(X) = E conj(X)^{-1} E the adjoint involution of h_1"""
    Y = mat_inv(K, mat_conj(K, X))
    A, B, C = circle_perp(K, Y)     # tau(X)(i Rhat) = -(conj(X)^{-1}(i Rhat)); the sign does not change y
    if A == 0:
        return None
    return abs(K.imK(B))


def perp_dictionary(K, ymax, say=print):
    """every form of disc D_perp(y) (imprimitive included), both orientations, is realized at level y
    (constructive: descent); the class number census; the involution tau does not preserve the level."""
    rows = []
    tot = 0
    tau_changes = 0
    for y in range(1, ymax + 1):
        D, sc = disc_perp(K, y)
        if K.t == 0 and y % 2:
            # even d_K: only even y (1 - m y^2 must be a discriminant)
            assert not is_disc(D) or D >= -4 or True
            continue
        if D >= -3:
            continue
        forms = classes_of_disc(D)
        levels_tau = set()
        for f in forms:
            for s in (1, -1):
                circ = circle_perp_of_form(K, f, y, s)
                X = descend_perp(K, circ)
                tot += 1
                yt = perp_level_of_tau(K, X)
                levels_tau.add(yt)
        if levels_tau != {y}:
            tau_changes += 1
        rows.append((y, D, len(forms), 3 * hurwitz_H(-D), sorted(str(v) for v in levels_tau)[:4]))
    return rows, tot, tau_changes


# ---- iS at Q(i) through the coset det X = i (Paper I's normalization) ----

def iS_class_formula(n, check_lemmas=True):
    """Level n even of iS: for every primitive form f of disc 1 - n^2 and both orientations s,
    the explicit P in GL_2(Z[i]) with det P = -i (Lemma A'': K_f = {u : -i conj(u) - beta u in 2a Z[i]},
    beta = -s(n + b i), v_k = s(-i conj(u_k) - beta u_k)/(2 i a)), X = P^{-1} of determinant i realizing
    the circle [[2sa, i beta], [-i conj(beta), 2sc]], and the class of sigma(X) = conj(X)^{-1}:
    sigma[f] = [r_n][f]^{-s}, r_n = (n-1, n-1, n/2) ~ (n/2, 1, n/2), reflection into H iff s = +1.
    Lemma B'': the Gram form of the image is g_s(u) = n N(u) - s Im(u^2) = 2 (n/2, -s, n/2)(x, y);
    Lemma C'': iota(x + i y) = x n/2 + y (-s + sqrt D)/2 maps K_f onto t a_f, t the ideal of (n/2, s, n/2)."""
    K = FieldX(-4)
    D = 1 - n * n
    prim = [f for f in classes_of_disc(D) if form_is_primitive(f)]
    rn = twist_class(D, (n - 1, n - 1, n // 2))
    assert rn == twist_class(D, (n // 2, 1, n // 2))
    delta = (0, -1)                     # det P = -i, det X = i
    M0 = (((0, 0), (0, 1)), ((0, -1), (0, 0)))
    count = 0
    for s in (1, -1):
        for f in classes_of_disc(D):
            a, b, c = f
            beta = (-s * n, -s * b)
            circ = (s * a, beta, s * c)
            assert K.norm(beta) == 1 + 4 * a * c

            def phi(u):
                return K.sub(K.mul(delta, K.conj(u)), K.mul(beta, u))
            e1, e2 = phi((1, 0)), phi((0, 1))
            Cm = [[e1[0], e2[0]], [e1[1], e2[1]]]
            Kf = hnf(preimage_lattice(Cm, 2 * a))
            assert K.index(Kf) == a, (n, f, Kf)
            u1, u2 = (Kf[0], 0), (Kf[1], Kf[2])
            if K.imK(K.mul(u1, K.conj(u2))) != s * a:
                u1, u2 = u2, u1
            assert K.imK(K.mul(u1, K.conj(u2))) == s * a
            vs = []
            for u in (u1, u2):
                z = phi(u)
                assert z[0] % (2 * a) == 0 and z[1] % (2 * a) == 0
                w = K.mul((0, -1), (z[0] // (2 * a), z[1] // (2 * a)))     # z/(2ia) = -i z/(2a)
                vs.append(K.smul(s, w))
            P = ((u1, vs[0]), (u2, vs[1]))
            assert mat_det(K, P) == delta, (n, f, mat_det(K, P))
            Pd = mat_conj(K, ((P[0][0], P[1][0]), (P[0][1], P[1][1])))
            M = mat_mul(K, mat_mul(K, Pd, M0), P)
            assert M == (((2 * s * a, 0), K.mul((0, 1), beta)), (K.mul((0, -1), K.conj(beta)), (2 * s * c, 0))), (n, f, M)
            X = mat_mul(K, mat_inv(K, P), (((0, 1), (0, 0)), ((0, 0), (0, 1))))   # adj(P)/delta = i adj(P)
            assert mat_det(K, X) == (0, 1)
            assert OF.circle_of(K, X) == circ, (n, f, OF.circle_of(K, X), circ)
            f2, alpha2, s2, refl = OF.sigma_image(K, X)
            assert alpha2 == n and s2 == -1 and refl == (s == 1), (n, f, alpha2, s2, refl)
            if not form_is_primitive(f):
                continue
            finv = reduce_form(a, -b, c)
            pred = compose(rn, finv if s == 1 else f, D)
            assert f2 == pred, (n, s, f, f2, pred)
            count += 1
            if not check_lemmas:
                continue
            # Lemma B'': N = X^dagger M_0 X has N_11 = g(u2)/a, N_22 = g(u1)/a, g(u) = n N(u) - s Im(u^2)

            def g(u):
                return n * K.norm(u) - s * K.mul(u, u)[1]

            def gpol(u, v):
                return (g(K.add(u, v)) - g(u) - g(v)) // 2
            Xd = mat_conj(K, ((X[0][0], X[1][0]), (X[0][1], X[1][1])))
            N = mat_mul(K, mat_mul(K, Xd, M0), X)
            assert N[0][0] == (g(u2) // a, 0) and N[1][1] == (g(u1) // a, 0) and g(u2) % a == 0, (n, f, N)
            fB = (g(u2) // a, -2 * gpol(u1, u2) // a, g(u1) // a)
            assert (2 * gpol(u1, u2)) % a == 0 and all(v % 2 == 0 for v in fB)
            fB = tuple(v // 2 for v in fB)
            assert fB[1] ** 2 - 4 * fB[0] * fB[2] == D and reduce_form(*fB) == f2, (n, f, fB, f2)
            # Lemma C'': iota(x + i y) = x n/2 + y (-s + sqrt D)/2; N(iota(u)) = (n/4) g(u); iota(K_f) = t a_f
            t0 = n // 2

            def iota(u):
                x, y = u
                return (Fraction(x * t0) - Fraction(y * s, 2), Fraction(y, 2))
            for u in (u1, u2, K.add(u1, u2)):
                z = iota(u)
                assert z[0] * z[0] - D * z[1] * z[1] == Fraction(t0 * g(u), 2), (n, f, u)
            t_ideal = [iota((1, 0)), iota((0, 1))]
            tform = (t0, s, t0)
            assert lat_eq(D, t_ideal, ideal_of_form(D, tform)), (n, s, t_ideal)
            af = sqrtD_ideal_basis(D, a, -b)
            assert lat_eq(D, [iota(u1), iota(u2)], lat_mul(D, t_ideal, af)), (n, f, "iota(K_f) != t a_f")
            assert twist_class(D, tform) == rn
    return count


def even_levels_phase(say=print, nmax=16, fields_perp=None, ymax=None):
    t0 = time.time()
    say("=" * 78)
    say("PHASE E  the even levels of iS = PSL_2(Z[i]).(i Rhat), and the second-kind orbits S_K^perp")
    say("=" * 78)
    K = FieldX(-4)
    levels = classification_perp(K, Qmax=6, say=say)
    assert sorted(levels)[:6] == [2, 4, 6, 8, 10, 12] and all(levels[y] == {1, -1} for y in levels)
    say("  (E1) iS = {curvature 2a, centre (x + n i)/(2a): x odd, n even, x^2 + n^2 = 1 mod 4a}; levels: the even n,")
    say("       both orientations; level-n circles = all forms of disc 1 - n^2 (census 3H(n^2-1))")
    tot = 0
    rows = []
    for n in range(4, nmax + 1, 2):
        D = 1 - n * n
        cnt = iS_class_formula(n, check_lemmas=True)
        tot += cnt
        prim = [f for f in classes_of_disc(D) if form_is_primitive(f)]
        rn = twist_class(D, (n - 1, n - 1, n // 2))
        rows.append(f"n={n}: D={D}, h={len(prim)}, r_n={rn}, 3H={3 * hurwitz_H(n * n - 1)}")
    say("  (E2) " + "; ".join(rows))
    say(f"  (E2) sigma[f] = [r_n][f]^(-s) (r_n = (n-1, n-1, n/2) ~ (n/2, 1, n/2); reflection iff s = +1) at "
        f"{tot} (class, orientation) pairs, even n <= {nmax}, with Lemmas A''-C'' (det P = -i, Gram form "
        f"n N(u) - s Im(u^2), iota(K_f) = t a_f) at every one")
    fields_perp = [-3, -7, -8, -11] if fields_perp is None else fields_perp
    out = {}
    for dK in fields_perp:
        K2 = FieldX(dK)
        classification_perp(K2, Qmax=5, say=say)
        ym = ymax if ymax else (6 if K2.t == 1 else 8)
        rows, cnt, tch = perp_dictionary(K2, ym, say=say)
        out[dK] = rows
        say(f"  (X2) {K2.name}: level y <-> all forms of disc {'4 - y^2|d_K|' if K2.t == 1 else '1 - m y^2 (y even)'}"
            f" (both orientations; {cnt} circles descended): " +
            "; ".join(f"y={y}: D={D}, {nf} forms, 3H={H}, tau-levels {lt}" for y, D, nf, H, lt in rows))
        say(f"  (X2) {K2.name}: the involution tau(X) = E conj(X)^{{-1}} E changes the R-hat level at {tch}/{len(rows)} levels"
            f" -- no level-preserving involution, hence no twist: S_K^perp realizes discriminants, not slices")
    say(f"  phase E/X done ({time.time() - t0:.1f} s)")
    return out


# ==========================================================================
# X'.  other carriers: the principal cusp at h_K > 1 and the non-Euclidean fields (explicit P, no descent)
# ==========================================================================

def lattice_K_fast(K, f, s, k):
    """K_f = {u in O_K : conj(u) - beta u in a sqrt(d_K) O_K} through a Smith form (no brute force)"""
    a, b, c = f
    d = -K.dK
    x = -s * b
    beta = (-s * k // 2, x) if K.t == 0 else ((-s * k - x) // 2, x)
    sdc = K.conj(K.sqrtd())

    def psi(u):
        y = K.sub(K.conj(u), K.mul(beta, u))
        return K.mul(y, sdc)                       # in a |d| O_K iff u in K_f
    e1, e2 = psi((1, 0)), psi((0, 1))
    Cm = [[e1[0], e2[0]], [e1[1], e2[1]]]
    return hnf(preimage_lattice(Cm, a * d)), beta


def class_formula_explicit(K, k, s, say=print, check_lemmas=True):
    """Theorem 3 of other-fields.md at the level 2 alpha = k of S_K through Lemma A' alone (no descent):
    valid for every K (h_K > 1 and the non-Euclidean fields included)."""
    d = -K.dK
    sl = slice_of_level(K.dK, k)
    D, r0, s0, gm, gp_, fr, fs, kind = sl
    alpha = Fraction(k, 2)
    prim = [f for f in classes_of_disc(D) if form_is_primitive(f)]
    frc = twist_class(D, fr)
    checked = 0
    for f in classes_of_disc(D):
        a, b, c = f
        Kf, beta = lattice_K_fast(K, f, s, k)
        assert K.index(Kf) == a, (K.dK, k, f, Kf)
        assert is_congruent_one(K, beta) and K.norm(beta) == 1 + d * a * c
        u1, u2 = (Kf[0], 0), (Kf[1], Kf[2])
        if K.imK(K.mul(u1, K.conj(u2))) != s * a:
            u1, u2 = u2, u1
        assert K.imK(K.mul(u1, K.conj(u2))) == s * a
        vs = []
        for u in (u1, u2):
            y = K.smul(s, K.sub(K.conj(u), K.mul(beta, u)))
            z = K.mul(y, K.conj(K.sqrtd()))
            assert z[0] % (a * d) == 0 and z[1] % (a * d) == 0
            vs.append((z[0] // (a * d), z[1] // (a * d)))
        P = ((u1, vs[0]), (u2, vs[1]))
        assert mat_det(K, P) == K.one, (K.dK, k, f, P)
        X = mat_inv(K, P)
        circ = (s * a, beta, s * c)
        assert OF.circle_of(K, X) == circ, (OF.circle_of(K, X), circ)
        f2, alpha2, s2, refl = OF.sigma_image(K, X)
        assert alpha2 == alpha and s2 == -1 and refl == (s == 1), (K.dK, k, f, alpha2, s2, refl)
        if not form_is_primitive(f):
            continue
        finv = reduce_form(a, -b, c)
        pred = compose(frc, finv if s == 1 else f, D)
        assert f2 == pred, (K.dK, k, s, f, f2, pred)
        checked += 1
        if check_lemmas:
            g11 = OF.gram_gs(K, alpha, s, u2)
            g22 = OF.gram_gs(K, alpha, s, u1)
            g12 = OF.gram_gs_pol(K, alpha, s, u2, K.neg(u1))
            fB = (Fraction(2, a * d) * g11, Fraction(2, a * d) * 2 * g12, Fraction(2, a * d) * g22)
            assert all(v.denominator == 1 for v in fB)
            assert reduce_form(*[int(v) for v in fB]) == f2, (K.dK, k, f, fB, f2)
            t0 = Fraction(2 * (alpha + s), d)
            assert t0.denominator == 1
            t0 = int(t0)

            def iota(u):
                return (t0 * K.re(u), Fraction(K.imK(u)) / 2)
            t_ideal = [iota(K.one), iota(K.w)]
            af = sqrtD_ideal_basis(D, a, -b)
            assert lat_eq(D, [iota(u1), iota(u2)], lat_mul(D, t_ideal, af)), (K.dK, k, f, "iota(K_f) != t a_f")
            for u in (u1, u2, K.add(u1, u2)):
                z = iota(u)
                assert z[0] * z[0] - D * z[1] * z[1] == t0 * Fraction(2, d) * OF.gram_gs(K, alpha, s, u)
    return D, frc, len(prim), checked


def other_carriers_phase(say=print):
    t0 = time.time()
    say("=" * 78)
    say("PHASE X'  other carriers: Lemma A' without descent at the non-Euclidean fields and at the principal cusp of h_K = 2")
    say("=" * 78)
    # (X3) regression: the explicit route agrees with other_fields' descent route at the Euclidean fields
    for dK, ks in ((-3, [7, 8, 11]), (-8, [14, 22]), (-7, [12, 16]), (-11, [20, 24]), (-4, [10, 14, 18])):
        K = FieldX(dK)
        for k in ks:
            D, signs = OF.level_data(K, k)
            for s in signs:
                D2, frc, h, ch = class_formula_explicit(K, k, s, say=say)
                Kf_fast, _ = lattice_K_fast(K, (2, 1, (1 - D) // 8) if False else classes_of_disc(D)[-1], s, k)
                Kf_slow, _ = OF.lattice_K(K, classes_of_disc(D)[-1], s, k)
                assert Kf_fast == Kf_slow
    say("  (X3) explicit-P route = descent route (other_fields.py) at 14 Euclidean levels; K_f by Smith form = brute force")
    # (X4) the non-Euclidean class-number-one fields: the first levels
    for dK, ks in ((-19, [17, 21, 36, 40]), (-43, [41, 45, 84, 88]), (-67, [65, 69, 132, 136]), (-163, [161, 165, 324, 328])):
        K = FieldX(dK)
        rows = []
        for k in ks:
            sl = slice_of_level(dK, k)
            s = -1 if (k - 2) % (-dK) == 0 else 1      # 2 alpha = -2s mod |d_K|
            assert (k + 2 * s) % (-dK) == 0
            D, frc, h, ch = class_formula_explicit(K, k, s, say=say)
            rows.append(f"k={k} (s={s:+d}): D={D}, h={h}, r={frc}{' (trivial)' if frc[0] == 1 else ''}, {ch} classes")
        say(f"  (X4) {K.name}: sigma[f] = [r_alpha][f]^(-s) with Lemmas A'-C' at " + "; ".join(rows))
    # (X5) the principal cusp at class number two: the levels k = +-2 mod |d_K| (Theorem 1(4) of other-fields.md)
    for dK in (-20, -24, -40, -52):
        K = FieldX(dK)
        d = -dK
        rows = []
        for k in range(3, 3 * d + 3):
            if (k * k - 4) % d or k % d not in (2, d - 2):
                continue
            sl = slice_of_level(dK, k)
            D, r0, s0, gm, gp_, fr, fs, kind = sl
            if D >= -4:
                continue
            ok = []
            for s in (1, -1):
                b0 = k % 2
                beta = (-s * k // 2, -s * b0) if K.t == 0 else ((-s * k + s * b0) // 2, -s * b0)
                if is_congruent_one(K, beta):
                    ok.append(s)
            assert len(ok) == 1, (dK, k, ok)
            s = ok[0]
            D2, frc, h, ch = class_formula_explicit(K, k, s, say=say)
            rows.append(f"k={k} (s={s:+d}): D={D}, h={h}, r={frc}{' (trivial)' if frc[0] == 1 else ''}, {ch} classes")
            if len(rows) >= 4:
                break
        others = [k for k in range(3, 2 * d + 3) if (k * k - 4) % d == 0 and k % d not in (2, d - 2)]
        say(f"  (X5) {K.name} (h_K = {K.hK}): sigma[f] = [r_alpha][f]^(-s) with Lemmas A'-C' at the principal-cusp levels "
            + "; ".join(rows))
        say(f"       residues k^2 = 4 mod {d} that are not levels of the principal cusp (k != +-2 mod {d}): k = {others}"
            f" -> D = {[-(k * k - 4) // d for k in others]} (candidates for the non-principal cusp; not realized here)")
    say(f"  phase X' done ({time.time() - t0:.1f} s)")


# ==========================================================================
# G.  genus structure of a slice
# ==========================================================================

def legendre(m, p):
    r = m % p
    if r == 0:
        return 0
    return 1 if pow(r, (p - 1) // 2, p) == 1 else -1


def assigned_characters(D):
    """[(name, chi)] the assigned characters of disc D (Cox, Thm 3.15 for D = -4n; odd primes for odd D)"""
    m = -D
    out = []

    def delta(x):
        return -1 if x % 4 == 3 else 1

    def eps(x):
        return -1 if x % 8 in (3, 5) else 1
    if D % 4 == 1:
        for p in factor(m):
            out.append((f"({p})", lambda x, p=p: legendre(x, p)))
        return out
    n = m // 4
    for p in factor(n):
        if p != 2:
            out.append((f"({p})", lambda x, p=p: legendre(x, p)))
    if n % 4 == 3:
        pass
    elif n % 4 == 1 or n % 8 == 4:
        out.append(("delta", delta))
    elif n % 8 == 2:
        out.append(("delta*eps", lambda x: delta(x) * eps(x)))
    elif n % 8 == 6:
        out.append(("eps", eps))
    else:
        out.append(("delta", delta))
        out.append(("eps", eps))
    assert len(out) == gauss_mu(D)
    return out


def represented_coprime(f, D, bound=30):
    a, b, c = f
    for x in range(0, bound):
        for y in range(-bound, bound):
            m = a * x * x + b * x * y + c * y * y
            if m > 0 and gcd(m, D) == 1:
                return m
    raise ValueError((f, D))


def genus_vector(D, f, chars=None):
    chars = assigned_characters(D) if chars is None else chars
    m = represented_coprime(f, D)
    return tuple(chi(m) for name, chi in chars)


def genus_phase(say=print, slices=None):
    t0 = time.time()
    say("=" * 78)
    say("PHASE G  genus structure of the slices: the assigned characters of r_alpha, odd genus characters, real fields")
    say("=" * 78)
    # (G1) at every level of the nine fields with |D| <= 300: genus vector of r_alpha = ((g_+/p) at p | r0, (-g_-/p) at p | s0)
    nlev = 0
    real_fields = {}
    two_adic = {}
    for dK in CLASS_NUMBER_ONE:
        for k, sl in levels_of_field(dK, isqrt(300 * (-dK) + 4) + 1):
            D, r0, s0, gm, gp_, fr, fs, kind = sl
            if D < -300 or D >= -4:
                continue
            chars = assigned_characters(D)
            gv = genus_vector(D, fr)
            # genus theory: the genus vector is a class function, a homomorphism with kernel Cl^2
            prim, coords, orders, unit, mul, inv, tt = class_group(D)
            G = {f: genus_vector(D, f, chars) for f in prim}
            assert G[reduce_form(*fr)] == gv or G[twist_class(D, fr)] == gv
            for f in prim:
                for g in prim:
                    assert tuple(x * y for x, y in zip(G[f], G[g])) == G[mul(f, g)]
            sq = set(mul(f, f) for f in prim)
            assert set(f for f in prim if all(x == 1 for x in G[f])) == sq
            assert len(set(G.values())) == 2 ** (gauss_mu(D) - 1)
            # the odd-prime law
            for (name, chi), val in zip(chars, gv):
                if name.startswith("("):
                    p = int(name[1:-1])
                    if r0 % p == 0:
                        assert s0 % p and val == legendre(gp_, p), (dK, k, p, val, gp_)
                    else:
                        assert s0 % p == 0 and val == legendre(-gm, p), (dK, k, p, val, gm)
                else:
                    two_adic.setdefault((dK, name), set()).add((k % (16 * (-dK)) if dK % 2 else k % 32, val))
            nlev += 1
            # the odd quadratic characters and their real fields: chi = chi_{d1,d2}, D = d1 d2, d2 > 0
            quad = [(ks, cord, ph) for (ks, cord, ph) in all_characters(coords, orders) if cord == 2]
            rn = twist_class(D, fr)
            odd_real = [(ks, ph) for ks, cord, ph in quad if ph[rn] == Fraction(1, 2)]
            fields = []
            for ks, ph in odd_real:
                d2s = []
                for d1 in range(-abs(D), 0):
                    if D % d1 == 0 and is_disc(d1) and (D // d1) > 0 and (D // d1) % 4 in (0, 1):
                        d2 = D // d1
                        okc = all(((1 if ph[f] == 0 else -1) == legendre_kron(d1, represented_coprime(f, D))) for f in prim)
                        if okc:
                            d2s.append(d2)
                assert d2s, (dK, k, ks)
                fields.append(min(d2s))
            real_fields[(dK, k)] = (D, rn, fields)
    say(f"  (G1) at {nlev} levels (nine fields, |D| <= 300): the genus vector of r_alpha is (g_+/p) at odd p | r0 and (-g_-/p) at odd p | s0"
        f" (g_-+ = gcd(k -+ 2, |d_K|)); the assigned characters are class functions, multiplicative, kernel Cl(D)^2")
    for (dK, name), vals in sorted(two_adic.items()):
        by = {}
        for kk, v in vals:
            by.setdefault(kk, set()).add(v)
        assert all(len(v) == 1 for v in by.values()), (dK, name, by)
        say(f"       2-adic character {name} of r_alpha at {FieldX(dK).name}: determined by k mod {16 * (-dK) if dK % 2 else 32}: "
            + ", ".join(f"{kk}:{'+' if list(v)[0] == 1 else '-'}" for kk, v in sorted(by.items())))
    say("  (G2) odd genus characters (chi(r) = -1) and the real quadratic fields Q(sqrt d2) of their closed forms:")
    for (dK, k), (D, rn, fields) in sorted(real_fields.items()):
        if fields:
            say(f"       {FieldX(dK).name} k={k}: D={D}, r={rn}: d2 in {fields}")
    # (G3) the closed form L'(0,chi) = (2h(d1)/w(d1)) h(d2) log eps_{d2} C(0) on odd real characters of some slices, with PARI
    if OF.GP is not None and slices is not None:
        for D, fr in slices:
            res = slice_units(D, fr, with_klf=True)
            prim, coords, orders = res["prim"], res["coords"], res["orders"]
            rn = res["rn"]
            for (ks, cord, l) in res["odd"]:
                if cord != 2:
                    continue
                ph = next(p for (kk, co, p) in all_characters(coords, orders) if kk == ks)
                found = None
                for d1 in range(-abs(D), 0):
                    if D % d1 == 0 and is_disc(d1) and (D // d1) % 4 in (0, 1):
                        d2 = D // d1
                        if all(((1 if ph[f] == 0 else -1) == legendre_kron(d1, represented_coprime(f, D))) for f in prim):
                            found = (d1, d2)
                            break
                assert found, (D, fr, ks)
                d1, d2 = found
                d1s, d2s = OF.fundamental_disc(d1)[0], OF.fundamental_disc(d2)[0]
                out = gp_run(f"""print("H1 ", qfbclassno({d1s})); print("W1 ", if({d1s} == -3, 6, if({d1s} == -4, 4, 2)));
print("H2 ", qfbclassno({d2s})); u = quadunit({d2s}); print("LE ", log(abs(u)));""")
                dd = parse_tagged(out)
                closed = mpf(2 * int(dd["H1"])) / int(dd["W1"]) * int(dd["H2"]) * mpf(dd["LE"])
                ratio = l.real / closed
                # C(0) is a rational with small denominator: certify
                C0 = None
                for den in (1, 2, 3, 4, 6, 8, 12, 24):
                    if spare_of(ratio * den - nint(ratio * den)) >= need_digits():
                        C0 = Fraction(int(nint(ratio * den)), den)
                        break
                assert C0 is not None, (D, fr, ks, nstr(ratio, 30))
                say(f"  (G3) D = {D}, r = {rn}, odd genus character chi_{{{d1},{d2}}} (fundamental {d1s}, {d2s}): "
                    f"L'(0,chi) = (2h({d1s})/w) h({d2s}) log eps_{d2s} * C(0), C(0) = {C0}; real field Q(sqrt {d2s})")
    say(f"  phase G done ({time.time() - t0:.1f} s)")
    return real_fields


def legendre_kron(d, m):
    """Kronecker symbol (d/m) for m > 0 coprime to d"""
    return OF.Field.chi(type("F", (), {"dK": d})(), 2) if False else kron(d, m)


def kron(d, m):
    """Kronecker symbol (d/m), m > 0"""
    res = 1
    for p, e in factor(m).items():
        if p == 2:
            if d % 2 == 0:
                return 0
            v = 1 if d % 8 in (1, 7) else -1
        else:
            v = legendre(d, p)
        if v == 0:
            return 0
        res *= v ** e
    return res


# ==========================================================================
# I.  the odd index on a slice (robert-index-full.md Thm 3 with an arbitrary twist), PARI
# ==========================================================================

GP_HEADER = """default(parisize, 2000000000);
default(realprecision, {prec0});
"""


def odd_index_slice(D, fr, res=None, prec_pari=120, certdeg=16, say=print, tag=""):
    """[E^- : <R_b>] = 24^{h/2} (2^{h/2-1}/Q^-) (h_H w_H+)/(h_H+ w_H) prod_{chi odd} C_chi(0) on the slice (D, [fr]):
    the twisted Dedekind determinant (I1), the multipliers through exact projections (I2), PARI (I3)."""
    from robert_index_full import project_form, fundamental_part
    need = need_digits()
    t0 = time.time()
    res = slice_units(D, fr, with_klf=True) if res is None else res
    prim, coords, orders, R, rn, mul, inv, unit = (res["prim"], res["coords"], res["orders"], res["R"], res["rn"],
                                                   res["mul"], res["inv"], res["unit"])
    h = len(prim)
    assert rn != unit, "trivial slice: no odd units"
    co = res["poly"]
    dK, fc = fundamental_part(D)
    chars = [(ks, cord, ph) for (ks, cord, ph) in all_characters(coords, orders) if cord > 1]
    chis = [{f: cval(ph[f]) for f in prim} for (ks, cord, ph) in chars]
    Lp, M = epstein_Lprime0_forms(prim, D, chis)
    fl = {f: log(fabs(R[f])) for f in prim}
    odd, S = [], {}
    for (ks, cord, ph), chi, l in zip(chars, chis, Lp):
        s_ = sum(chi[f] * fl[f] for f in prim)
        S[ks] = s_
        if ph[rn] == Fraction(1, 2):
            odd.append((ks, cord, ph, chi, l))
            assert spare_of((s_ + 24 * l) / max(1, fabs(l))) >= need
        else:
            assert fabs(s_) < mpf(10) ** (-need)
    assert len(odd) == h // 2
    # (I1) twisted group determinant over a transversal of <r>
    T, seen = [], set()
    for f in prim:
        if f not in seen:
            T.append(f)
            seen.add(f)
            seen.add(mul(f, rn))
    N_ = matrix(len(T), len(T))
    for i, a in enumerate(T):
        ai = inv(a)
        for j, f in enumerate(T):
            N_[i, j] = 2 * fl[mul(ai, f)]
    dN = det(N_)
    prodS, prodL = mpc(1), mpf(1)
    for ks, cord, ph, chi, l in odd:
        prodS *= S[ks]
        prodL *= fabs(l)
    assert spare_of((dN - prodS) / max(1, fabs(prodS))) >= need, (D, fr, nstr(dN, 20), nstr(prodS, 20))
    # (I2) multipliers through exact projections
    cond_divs = [d_ for d_ in range(1, fc + 1) if fc % d_ == 0]
    proj = {}
    for fp in cond_divs:
        if fp == fc:
            continue
        Dp = fp * fp * dK
        forms_p = [g for g in classes_of_disc(Dp) if form_is_primitive(g)]
        pm = {}
        for f in prim:
            g = project_form(f, D, dK, fp)
            if g not in forms_p:
                g = (g[0], -g[1], g[2])
            assert g in forms_p, (D, fp, f, g)
            pm[f] = g
        assert set(pm.values()) == set(forms_p), (D, fp)
        proj[fp] = (Dp, forms_p, pm)
    prodC = mpf(1)
    Cdesc = []
    for ks, cord, ph, chi, l in odd:
        levels = [fc]
        for fp in cond_divs:
            if fp == fc:
                continue
            Dp, forms_p, pm = proj[fp]
            val, ok = {}, True
            for f in prim:
                g = pm[f]
                if g in val and val[g] != ph[f]:
                    ok = False
                    break
                val[g] = ph[f]
            if ok:
                levels.append(fp)
        fmin = min(levels)
        if fmin == fc:
            C = mpf(1)
            Cdesc.append(f"chi{ks} (order {cord}) primitive")
        else:
            Dp, forms_p, pm = proj[fmin]
            val = {pm[f]: ph[f] for f in prim}
            Lpp, _ = epstein_Lprime0_forms(forms_p, Dp, [{g: cval(val[g]) for g in forms_p}])
            C = (l / Lpp[0]).real
            Cdesc.append(f"chi{ks} (order {cord}) from conductor {fmin} (disc {Dp}), C = {nstr(C, 8)}")
        prodC *= C
    prodC_int, spC = cert_int(prodC, "prod C (slice)")
    # (I3) PARI
    Rpol = "Pol(" + str(co) + ", x)"
    script = GP_HEADER.format(prec0=prec_pari) + f"""
P = polclass({D}, 0, 'y);
H = polredbest(polcompositum(y^2 - ({D}), P)[1]);
print("HPOL ", H);
bnf = bnfinit(H, 1);
print("H ", bnf.no); print("CYC ", bnf.cyc); print("W ", bnf.tu[1]); print("REG ", bnf.reg);
print("CERT ", if(poldegree(H) <= {certdeg}, bnfcertify(bnf), -1));
r = nfroots(bnf.nf, {Rpol});
print("NROOTS ", #r);
nu = #bnf.fu;
E = matrix(#r, nu); TT = vector(#r);
for(i=1, #r, e = bnfisunit(bnf, r[i]); for(j=1, nu, E[i,j] = e[j]); TT[i] = lift(e[#e]));
print("TORS ", TT);
G = nfgaloisconj(bnf);
sd = nfroots(bnf.nf, x^2 - ({D}))[1];
tau = 0; ntau = 0;
for(k=1, #G, ok = 1; for(i=1, #r, if(nfgaloisapply(bnf, G[k], r[i]) != 1/r[i], ok = 0)); if(ok && nfgaloisapply(bnf, G[k], sd) == sd, tau = G[k]; ntau++));
print("NTAU ", ntau);
Tm = matrix(nu, nu);
for(j=1, nu, e = bnfisunit(bnf, nfgaloisapply(bnf, tau, bnf.fu[j])); for(i=1, nu, Tm[i,j] = e[i]));
Km = matkerint(matid(nu) + Tm); Kp = matkerint(matid(nu) - Tm);
print("KMRANK ", #Km); print("KPRANK ", #Kp);
X = matrix(#Km, #r);
for(i=1, #r, x = matinverseimage(Km, E[i,]~); for(k=1, #Km, X[k,i] = x[k]));
sn = matsnf(X); idx = 1; for(k=1, #sn, if(sn[k] != 0, idx *= abs(sn[k]))); print("INDEX ", idx);
print("XSNF ", sn);
w = bnf.tu[1]; zt = bnf.tu[2]; et = bnfisunit(bnf, nfgaloisapply(bnf, tau, zt)); t = lift(et[#et]);
g0 = gcd(w, t - 1); gk = g0;
for(k=1, #Kp, u = factorback(bnf.fu, Kp[,k]); e = bnfisunit(bnf, nfgaloisapply(bnf, tau, u) / u); gk = gcd(gk, lift(e[#e])));
ee = g0 / gk;
Qm = ee * abs(matdet(concat(Kp, Km))); print("QMINUS ", Qm); print("EPLUS ", ee); print("TAUMU ", t);
z = sd; c = 1;
while(poldegree(minpoly(z)) < {h} && c < 40, z = sd + sum(i=1, #r, c^i * (r[i] + 1/r[i])); c++);
print("HPDEG ", poldegree(minpoly(z)));
Hp = polredbest(minpoly(z));
print("HPPOL ", Hp);
bp = bnfinit(Hp, 1);
print("HPH ", bp.no); print("HPCYC ", bp.cyc); print("HPW ", bp.tu[1]); print("HPREG ", bp.reg);
print("HPCERT ", if(poldegree(Hp) <= {certdeg}, bnfcertify(bp), -1));
"""
    d = parse_tagged(gp_run(script))
    assert int(d["NROOTS"]) == h and int(d["NTAU"]) == 1, (D, fr, d.get("NROOTS"), d.get("NTAU"))
    assert int(d["HPDEG"]) == h and int(d["KMRANK"]) == h // 2 and int(d["KPRANK"]) == h // 2 - 1, (D, fr, d["HPDEG"], d["KMRANK"], d["KPRANK"])
    idx, Qm, eplus = int(d["INDEX"]), int(d["QMINUS"]), int(d["EPLUS"])
    hH, wH, RH = int(d["H"]), int(d["W"]), mpf(d["REG"])
    hHp, wHp, RHp = int(d["HPH"]), int(d["HPW"]), mpf(d["HPREG"])
    cert = int(d["CERT"]) == 1 and int(d["HPCERT"]) == 1
    Rminus = fabs(dN) / idx
    spB = spare_of(((RH / RHp) - 2 ** (h // 2 - 1) * Rminus / Qm) / (RH / RHp))
    lhs = (hH * RH / wH) / (hHp * RHp / wHp)
    spA = spare_of((lhs * prodC_int - prodL) / prodL)
    pred = Fraction(24 ** (h // 2) * 2 ** (h // 2 - 1), Qm) * Fraction(hH * wHp, hHp * wH) * prodC_int
    snf = [abs(int(x)) for x in parse_int_vector(d["XSNF"]) if int(x) != 0]
    sat = 1
    for dd in snf:
        sat *= gcd(dd, 24)
    say(f"  (I) {tag}D = {D}, [r] = {rn}, h = {h}: R-polynomial {co if len(co) <= 5 else str(co[:2])[:-1] + ', ...]'}; "
        f"odd characters {len(odd)}; det(2 log|R_(a^-1 f)|)_(T x T) = prod_odd S_chi (spare >= {need}); {'; '.join(Cdesc)}; prod C = {prodC_int}")
    say(f"      PARI: H = ring class field of disc {D} (degree {2 * h}): h_H = {hH}, Cl(H) = {parse_int_vector(d['CYC'])}, w_H = {wH}; "
        f"H^+ = H^tau: h = {hHp}, Cl = {parse_int_vector(d['HPCYC'])}, w = {wHp}; Q^- = {Qm}; "
        f"{'bnfcertify = 1 for H and H^+ (unconditional)' if cert else 'GRH-conditional'}")
    say(f"      R_H/R_H+ = 2^(h/2-1) R^-(E^-)/Q^- (spare {spB}); (hR/w)_H/(hR/w)_H+ = prod L'/prod C (spare {spA}); "
        f"exact index [E^- : <R_b>] = {idx} = 24^{h // 2} * {Fraction(idx, 24 ** (h // 2))} = formula {pred}: {'OK' if pred == idx else 'MISMATCH'}; "
        f"Smith form {snf}, 24-saturated index {idx // sat}   ({time.time() - t0:.1f} s)")
    assert spB >= need and spA >= need and pred == idx, (D, fr, idx, pred, spB, spA)
    return {"D": D, "rn": rn, "h": h, "idx": idx, "Qm": Qm, "eplus": eplus, "hH": hH, "hHp": hHp, "wH": wH, "wHp": wHp,
            "cyc": parse_int_vector(d["CYC"]), "cycp": parse_int_vector(d["HPCYC"]), "prodC": prodC_int, "cert": cert,
            "snf": snf, "sat": idx // sat, "poly": co}


INDEX_TARGETS = [
    # (D, twist form, tag)
    (-15, (3, 3, 2), "iS n=4: "), (-35, (5, 5, 3), "iS n=6: "), (-63, (7, 7, 4), "iS n=8: "), (-99, (9, 9, 5), "iS n=10: "),
    (-195, (13, 13, 7), "iS n=14: "),
    (-84, (14, 14, 5), "Q(sqrt-3) k=16: "), (-84, (3, 0, 7), "Q(sqrt-2) alpha=13: "), (-84, (2, 2, 11), "Q(sqrt-19) alpha=20: "),
    (-160, (8, 8, 7), "Q(sqrt-3) k=22: "), (-160, (4, 4, 11), "Q(sqrt-11) k=42: "), (-160, (5, 0, 8), "no arrangement: "),
    (-168, (3, 0, 14), "Q(sqrt-5) k=58 (principal cusp): "),
]

INDEX_RECORD = {}     # filled by the selftest run; regression targets recorded in horizontal-families.md


def index_phase(say=print, targets=None, with_480=False):
    t0 = time.time()
    say("=" * 78)
    say(f"PHASE I  the odd index on slices with PARI (dps = {mp.dps})")
    say("=" * 78)
    assert OF.GP is not None, "PARI/GP needed"
    targets = INDEX_TARGETS if targets is None else targets
    if with_480:
        targets = targets + [(-480, (38 - 2, 38 - 2, (36 + 40) // 4), "Q(sqrt-3) k=38: "), (-480, (56 // 7 * 0 + 8, 8, (8 + 60) // 4), "Q(sqrt-7) k=58: "),
                             (-480, (15, 0, 8), "Q(sqrt-2) alpha=31: ")]
    recs = {}
    for D, fr, tag in targets:
        rec = odd_index_slice(D, fr, say=say, tag=tag)
        recs[(D, rec["rn"])] = rec
        if D in INDEX_RECORD_EXPECTED:
            exp_ = INDEX_RECORD_EXPECTED[D].get(rec["rn"])
            if exp_ is not None:
                assert (rec["idx"], rec["Qm"], rec["hH"], rec["hHp"], rec["wH"], rec["wHp"]) == exp_, (D, rec["rn"], rec["idx"], rec["Qm"], exp_)
    say(f"  phase I done ({time.time() - t0:.1f} s)")
    return recs


# regression record: (D, twist class) -> ([E^- : <R_b>], Q^-, h_H, h_{H+}, w_H, w_{H+}); all bnfcertify = 1
INDEX_RECORD_EXPECTED = {
    -15: {(2, 1, 2): (4, 1, 1, 2, 6, 2)},
    -35: {(3, 1, 3): (12, 1, 1, 2, 2, 2)},
    -63: {(4, 1, 4): (576, 2, 1, 1, 6, 6)},
    -99: {(5, 1, 5): (8, 1, 1, 1, 6, 2)},
    -195: {(7, 1, 7): (576, 2, 4, 4, 6, 6)},
    -84: {(5, 4, 5): (96, 2, 1, 2, 12, 4), (3, 0, 7): (144, 2, 1, 2, 12, 6), (2, 2, 11): (48, 2, 1, 2, 12, 2)},
    -160: {(7, 6, 7): (288, 2, 1, 2, 8, 2), (4, 4, 11): (144, 2, 1, 1, 8, 2), (5, 0, 8): (576, 2, 1, 2, 8, 4)},
    -168: {(3, 0, 14): (96, 2, 1, 2, 6, 2)},
}


# ==========================================================================
# F.  the phase u = Phi_y/Phi_x on a slice
# ==========================================================================

PHI_CACHE = {}


def is_prime(m):
    return m > 1 and all(m % p for p in range(2, isqrt(m) + 1))


def phi_from_pari(m):
    """exact Phi_m(x, y) from PARI's polmodular (prime m), as {(i, j): coeff}"""
    import sympy as sp
    out = gp_run(f"""default(parisize, 400000000);
P = polmodular({m}, 0, 'x, 'y);
print("PHI ", P);""")
    d = parse_tagged(out)
    x, y = sp.symbols("x y")
    expr = sp.sympify(d["PHI"].replace("^", "**"))
    P = sp.Poly(expr, x, y)
    phi = {}
    for (i, j), c in P.terms():
        phi[(i, j)] = int(c)
    # symmetry check
    for (i, j), c in phi.items():
        assert phi.get((j, i)) == c, (m, i, j)
    return phi


def phi_exact(m):
    from first_power_descent import build_phi
    if m in PHI_CACHE:
        return PHI_CACHE[m]
    if m <= 12 or not is_prime(m):
        phi, psi = build_phi(m)
    else:
        phi = phi_from_pari(m)
    PHI_CACHE[m] = phi
    return phi


def j_theta(tau):
    """j(tau) via theta constants (moduli_invariants._E4E6D)"""
    q = exp(mpc(0, 1) * pi * tau)
    t2, t3, t4 = jtheta(2, 0, q), jtheta(3, 0, q), jtheta(4, 0, q)
    a, b, c = t2 ** 4, t3 ** 4, t4 ** 4
    E4 = (a * a + b * b + c * c) / 2
    Dl = (t2 * t3 * t4) ** 8 / 256
    return E4 ** 3 / Dl


def slice_phase(D, fr, fs, res=None, say=print, tag="", exact=True, target=100):
    """u_b = Phi_y/Phi_x(beta1, beta2), Phi = Phi_{r0} (r the twist of smaller norm), at every primitive class:
    the j-dressing u^6 = R beta1^4 (beta1-1728)^3/(beta2^4 (beta2-1728)^3), the partner law (Phi_{s0} gives -u),
    and (exact) the slice polynomial through Q[t]/(H_D) with its irreducibility."""
    from first_power_descent import phi_eval, phi_partial, hilbert_class_poly, F1, polyF_gcd
    import sympy as sp
    need = need_digits()
    t0 = time.time()
    res = slice_units(D, fr, with_klf=False) if res is None else res
    prim, R = res["prim"], res["R"]
    if fr[0] > fs[0]:
        fr, fs = fs, fr
    r0, s0 = fr[0], fs[0]
    assert r0 < s0
    phi = phi_exact(r0)
    phx, phy = phi_partial(phi, "x"), phi_partial(phi, "y")
    phi_s = phi_exact(s0) if (s0 <= 12 or is_prime(s0)) and s0 <= 40 else None
    sD = mpc(0, sqrt(mpf(-D)))
    u = {}
    worst = mp.dps
    worst_zero = mp.dps
    worst_partner = mp.dps
    for f in prim:
        a, b, c = f
        tau1 = (mpf(-b) + sD) / (2 * a)
        beta1 = j_theta(tau1)
        L = lat_mul(D, ideal_of_form(D, fr), ideal_of_form(D, f))
        w1, w2 = hnf_basis(D, L)
        tau2 = (mpf(w2[0].numerator) / w2[0].denominator + (mpf(w2[1].numerator) / w2[1].denominator) * sD) / (mpf(w1[0].numerator) / w1[0].denominator)
        beta2 = j_theta(tau2)
        val = phi_eval(phi, beta1, beta2)
        scale = max(fabs(v * beta1 ** i * beta2 ** jj) for (i, jj), v in phi.items())
        worst_zero = min(worst_zero, spare_of(fabs(val) / scale))
        num, den = phi_eval(phy, beta1, beta2), phi_eval(phx, beta1, beta2)
        u[f] = num / den
        lhs = u[f] ** 6
        rhs = R[f] * beta1 ** 4 * (beta1 - 1728) ** 3 / (beta2 ** 4 * (beta2 - 1728) ** 3)
        worst = min(worst, spare_of((lhs - rhs) / max(1, fabs(rhs))))
        if phi_s is not None:
            us = phi_eval(phi_partial(phi_s, "y"), beta1, beta2) / phi_eval(phi_partial(phi_s, "x"), beta1, beta2)
            worst_partner = min(worst_partner, spare_of((us + u[f]) / max(1, fabs(u[f]))))
    assert worst >= target and worst_zero >= target, (D, fr, worst, worst_zero)
    if phi_s is not None:
        assert worst_partner >= target, (D, fr, worst_partner)
    out = {"u": u, "spare": worst, "zero": worst_zero, "partner": worst_partner if phi_s is not None else None}
    msg = (f"  (F) {tag}D = {D}, r = {fr} (norm {r0}), s = {fs} (norm {s0}), h = {len(prim)}: Phi_{r0}(beta1, beta2) = 0 (spare {worst_zero}); "
           f"u^6 = R beta1^4(beta1-1728)^3/(beta2^4(beta2-1728)^3) at every class (spare {worst})"
           + (f"; Phi_{s0} gives -u (spare {worst_partner})" if phi_s is not None else f"; partner Phi_{s0} not computed"))
    if exact:
        dps0 = mp.dps
        H, hh = hilbert_class_poly(D, dps=max(dps0, 160))
        mp.dps = dps0
        assert hh == len(prim)
        x = sp.symbols("x")
        assert sp.Poly(H, x).is_irreducible
        F = F1(H)
        degy = max(jj for (_i, jj) in phi)
        degx = max(i for (i, _jj) in phi)
        tpowers = [F.one()]
        for _ in range(degx):
            tpowers.append(F.mul(tpowers[-1], F.t()))
        phi_y = [F.zero() for _ in range(degy + 1)]
        for (i, jj), v in phi.items():
            phi_y[jj] = F.add(phi_y[jj], F.scal(tpowers[i], Fraction(v)))
        H_y = [F.elem([Fraction(cc)]) for cc in reversed(H)]
        g = polyF_gcd(F, H_y, phi_y)
        branch_note = ""
        if len(g) > 2:
            # several primitive invertible ideals of norm r0: isolate the Atkin-Lehner branch on the fiber
            # product X_0(r0) x X_0(s0) (Phi_{s0}), or divide out the branches of the other ambiguous ideals of
            # norm r0 through their own partners
            def phi_in_y(ph):
                dy = max(jj for (_i, jj) in ph)
                dx = max(i for (i, _jj) in ph)
                tp = [F.one()]
                for _ in range(dx):
                    tp.append(F.mul(tp[-1], F.t()))
                out_ = [F.zero() for _ in range(dy + 1)]
                for (i, jj), v in ph.items():
                    out_[jj] = F.add(out_[jj], F.scal(tp[i], Fraction(v)))
                return out_
            e0 = len(g) - 1
            if phi_s is not None:
                g = polyF_gcd(F, g, phi_in_y(phi_s))
                branch_note = f"; {e0} ideals of norm {r0}: the branch isolated on X_0({r0}) x X_0({s0})"
            else:
                others = [(f2, s2) for f2, r2, s2 in ambiguous_forms(D) if r2 == r0 and reduce_form(*f2) != reduce_form(*fr)
                          and twist_class(D, f2) != twist_class(D, fr)]
                for f2, s2 in others:
                    if s2 <= 12 or (is_prime(s2) and s2 <= 40):
                        hq = polyF_gcd(F, g, phi_in_y(phi_exact(s2)))
                        if len(hq) == 2:
                            g = polyF_divexact(F, g, hq)
                branch_note = f"; {e0} ideals of norm {r0}: the other ambiguous branches divided out through their partners"
        assert len(g) == 2, (D, fr, "pairing gcd degree != 1: the branch could not be isolated exactly", len(g))
        beta2t = F.sub(F.zero(), g[0])

        def eval_F(p2):
            b2pow = [F.one()]
            degj = max(jj for (_i, jj) in p2) if p2 else 0
            for _ in range(degj):
                b2pow.append(F.mul(b2pow[-1], beta2t))
            acc = F.zero()
            for (i, jj), v in p2.items():
                acc = F.add(acc, F.scal(F.mul(tpowers[i], b2pow[jj]), Fraction(v)))
            return acc
        ut = F.mul(eval_F(phy), F.inv(eval_F(phx)))
        cols = []
        e = F.one()
        for k in range(hh):
            if k:
                e = F.mul(e, F.t())
            cols.append(F.mul(ut, e))
        Mx = sp.Matrix(hh, hh, lambda i, j: sp.Rational(cols[j][i]))
        Pi = sp.Poly(Mx.charpoly(x).as_expr(), x)
        co_ = Pi.all_coeffs()
        den = 1
        for cc in co_:
            den = sp.ilcm(den, sp.fraction(sp.Rational(cc))[1])
        ico = [int(cc * den) for cc in co_]
        cont = 0
        for cc in ico:
            cont = gcd(cont, cc)
        ico = [cc // cont for cc in ico]
        if ico[0] < 0:
            ico = [-cc for cc in ico]
        QP = sp.Poly(ico, x)
        irred = QP.is_irreducible
        sqf = sp.gcd(QP, QP.diff(x)).total_degree() == 0
        # the numerical roots agree with the exact polynomial
        num_co = poly_from_roots([u[f] for f in prim])
        spn = min(spare_of(num_co[i] * ico[0] - ico[i]) for i in range(len(ico)))
        assert spn >= need, (D, fr, spn)
        out.update({"Pi": ico, "irreducible": irred, "squarefree": sqf})
        msg += (f"; exact slice polynomial (via Q[t]/(H_D), pairing gcd of degree 1{branch_note}): "
                f"{ico if len(str(ico)) < 90 else str(ico[:2])[:-1] + ', ... (' + str(len(ico) - 1) + ' more)]'} "
                f"{'irreducible' if irred else 'REDUCIBLE'}, {'squarefree' if sqf else 'NOT squarefree'}; numerics match (spare {spn})")
    say(msg + f"   ({time.time() - t0:.1f} s)")
    return out


def polyF_divexact(F, A, B):
    """exact quotient A / B of polynomials over F1 (lists of F1 elements, low first; B monic)"""
    A = list(A)
    db = len(B) - 1
    assert F.is_zero(F.sub(B[db], F.one()))
    q = [F.zero() for _ in range(len(A) - db)]
    for k in range(len(A) - 1, db - 1, -1):
        c = A[k]
        q[k - db] = c
        for i in range(db + 1):
            A[k - db + i] = F.sub(A[k - db + i], F.mul(c, B[i]))
    assert all(F.is_zero(a) for a in A), "division not exact"
    return q


PHASE_TARGETS = [
    (-15, (3, 3, 2), (5, 5, 2), "iS n=4 / Q(sqrt-3) k=7: "), (-35, (5, 5, 3), (7, 7, 3), "iS n=6: "),
    (-63, (7, 7, 4), (9, 9, 4), "iS n=8: "), (-99, (9, 9, 5), (11, 11, 5), "iS n=10: "), (-195, (13, 13, 7), (15, 15, 7), "iS n=14: "),
    (-24, (2, 0, 3), (3, 0, 2), "Q(i) n=5: "), (-20, (2, 2, 3), (10, 10, 3), "Q(sqrt-3) k=8: "), (-36, (2, 2, 5), (18, 18, 5), "Q(sqrt-7) k=16: "),
    (-84, (6, 6, 5), (14, 14, 5), "Q(sqrt-3) k=16: "), (-84, (3, 0, 7), (7, 0, 3), "Q(sqrt-2) alpha=13: "), (-84, (2, 2, 11), (42, 42, 11), "Q(sqrt-19) alpha=20: "),
    (-160, (8, 8, 7), (20, 20, 7), "Q(sqrt-3) k=22: "), (-160, (4, 4, 11), (40, 40, 11), "Q(sqrt-11) k=42: "), (-160, (5, 0, 8), (8, 0, 5), "no arrangement: "),
    (-168, (3, 0, 14), (14, 0, 3), "Q(sqrt-5) k=58: "), (-56, (2, 0, 7), (7, 0, 2), "strata (Q(sqrt-3) k=26, Q(i) n=15): "),
]


def phase_phase(say=print, targets=None):
    t0 = time.time()
    say("=" * 78)
    say(f"PHASE F  the phase u = Phi_y/Phi_x(beta1, beta2) on slices (dps = {mp.dps})")
    say("=" * 78)
    targets = PHASE_TARGETS if targets is None else targets
    out = {}
    for D, fr, fs, tag in targets:
        out[(D, fr)] = slice_phase(D, fr, fs, say=say, tag=tag)
    say(f"  phase F done ({time.time() - t0:.1f} s)")
    return out


# ==========================================================================
# driver
# ==========================================================================

def selftest(quick=False):
    t00 = time.time()
    print("horizontal_families.py selftest" + (" (quick)" if quick else ""))
    timings = {}
    mp.dps = 30
    t = time.time(); slices_phase(Dmax=400); timings["S"] = time.time() - t
    t = time.time(); realize_phase(); timings["R"] = time.time() - t
    t = time.time(); pell_phase(J=4); timings["P"] = time.time() - t
    mp.dps = 200
    t = time.time(); units_phase(); timings["U"] = time.time() - t
    mp.dps = 30
    t = time.time(); even_levels_phase(nmax=16); timings["E"] = time.time() - t
    t = time.time(); other_carriers_phase(); timings["X"] = time.time() - t
    mp.dps = 120
    t = time.time(); genus_phase(slices=[(-84, (14, 14, 5)), (-84, (3, 0, 7)), (-84, (2, 2, 11)), (-120, (5, 0, 6)),
                                        (-168, (6, 0, 7)), (-168, (3, 0, 14)), (-160, (8, 8, 7)), (-160, (4, 4, 11))]); timings["G"] = time.time() - t
    mp.dps = 150
    t = time.time(); index_phase(); timings["I"] = time.time() - t
    mp.dps = 150
    t = time.time(); phase_phase(); timings["F"] = time.time() - t
    print("timings (s): " + ", ".join(f"{k} {v:.1f}" for k, v in timings.items()) + f"; total {time.time() - t00:.1f}")
    print("ALL PASS")


def main(argv):
    mode = argv[1] if len(argv) > 1 else "--selftest"
    if mode in ("--selftest", "selftest", "all"):
        selftest()
        return
    mp.dps = 30
    if mode == "slices":
        slices_phase(Dmax=int(argv[2]) if len(argv) > 2 else 400)
    elif mode == "realize":
        realize_phase()
    elif mode == "pell":
        pell_phase(J=int(argv[2]) if len(argv) > 2 else 4)
    elif mode == "units":
        mp.dps = 200
        units_phase()
    elif mode == "even":
        even_levels_phase(nmax=int(argv[2]) if len(argv) > 2 else 16)
    elif mode == "other":
        other_carriers_phase()
    elif mode == "genus":
        mp.dps = 120
        genus_phase(slices=[(-84, (14, 14, 5)), (-84, (3, 0, 7)), (-84, (2, 2, 11)), (-120, (5, 0, 6)),
                            (-168, (6, 0, 7)), (-168, (3, 0, 14)), (-160, (8, 8, 7)), (-160, (4, 4, 11))])
    elif mode == "index":
        mp.dps = 150
        index_phase(with_480="--with-480" in argv)
    elif mode == "phase":
        mp.dps = 150
        phase_phase()
    else:
        print(__doc__)


if __name__ == "__main__":
    main(sys.argv)
