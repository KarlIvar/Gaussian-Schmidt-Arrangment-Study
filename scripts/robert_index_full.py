"""The full Robert index of the Schmidt Delta-units.

Companion to robert-index-full.md.  Setting (Paper II, Thm 2.6; schmidt-euler-
system.md, Cor. 0.1): for a Euclidean level n >= 2 the Delta-data
    G_c = n^12 Delta(Lambda_c)/Delta(Z[i]),   c in Pic(O_n),  O_n = Z + nZ[i],
are algebraic integers of the ring class field H_n with
sigma_a(G_c) = G_{a^{-1}c} and conj(G_c) = G_{c^{-1}}, the same ideal (G_c) for
every class, and sum_c chi(c) log|G_c| = -12 L'(0,chi) for every nontrivial
character chi (Epstein L-function of the level).  The group
    V_n = < G_c/G_1 : c != 1 >  subset  O_{H_n}^x
is Galois-stable of rank h - 1 = unit rank of H_n, and the theorem of
robert-index-full.md says

    [O_{H_n}^x : mu(H_n) V_n] = 24^(h-1) prod_{chi != 1} |L'(0,chi)| / R_{H_n}
                              = (4 * 24^(h-1) / w_{H_n}) h_{H_n} prod_{chi != 1} C_chi(0),

C_chi(0) = L'(0,chi)/L'_prim(0,chi) the (integer) imprimitivity multiplier of
Paper II Lemma 6.3 / schmidt-euler-system.md Thm 2 (= 1 when chi factors
through no proper divisor level).

The script verifies, level by level (n = 3, 5, 7, 9, 11, 13, 15 and, with
--with-23, n = 23):

 (V1) the Dedekind-determinant identity behind the regulator of V_n:
      det( 2(f(a^{-1}c) - f(a^{-1})) )_{a,c != 1} = 2^(h-1) prod_{chi != 1} S_chi,
      f(c) = log|G_c|, S_chi = sum_c chi(c) f(c)  (signed identity; V_n has
      full rank because no S_chi vanishes);
 (V2) S_chi = -12 L'(0,chi) against the INDEPENDENT incomplete-gamma evaluation
      of the Epstein L'(0,chi) (Paper II Lemma 2.3; no modular quantity);
 (V3) the imprimitivity multipliers: for every chi the primitive level m | n is
      found exactly (kernels of the projections Pic(O_n) -> Pic(O_m)), C_chi(0)
      is computed from the recursion of schmidt-euler-system.md Thm 2 (exact
      character values), and checked against L'(0,chi)/L'(0,chi^{(m)}) with
      both L' from the independent evaluation;
 (V4) PARI/GP (bnfinit at 120 digits, bnfcertify where it finishes): h_{H_n},
      R_{H_n}, w_{H_n}; the index 24^(h-1) prod|L'| / R_{H_n} is certified as an
      integer with >= 40 spare digits (absolute-error criterion), and equals
      (4 * 24^(h-1)/w) h_{H_n} prod C_chi(0);
 (V5) the EXACT index: the roots of the certified integer polynomial
      D_n(x) = prod_c (x - G_c) are found inside H_n (nfroots), the units
      G_c/G_1 are written in PARI's fundamental-unit basis (bnfisunit), and
      the determinant of the integer exponent matrix reproduces the index;
      the Smith normal form of O_{H_n}^x / mu V_n is recorded;
 (V6) the layers: for every real character the quadratic-layer unit
      theta^(2) = prod_{ker} G / |M|^(1/2) has index 6 m_chi in Z[eps_{d_2}]^x
      (exact, against PARI's fundamental unit); at the cubic levels the coset
      unit theta_u (Paper II Thm 6.6) is re-expressed in the fundamental unit
      of the real cubic field L_3 (exponent +-8 h_{L_3} C_n(0), exact via
      bnfisunit -- an independent confirmation of Paper II's table), and in the
      sextic F = K(theta_u) (totally complex, NOT CM: unit rank 2) the ratio
      group N_{H_n/F}(V_n) has index (4 24^2/w_F) h_F C_chi3(0)^2 = 576 h_F C^2
      (Theorem 2 of the document, the layer theorem);
 (V7) the 24th-root saturation W_n = {u : u^24 in mu V_n}: from the Smith
      invariants d_i of the exponent matrix, [W_n : mu V_n] = prod gcd(d_i, 24)
      and [O_{H_n}^x : W_n] = h_{H_n} prod C_chi(0) at every level (§5.5 of the
      document; conjectured in general); the regression record below.

Certification policy (CLAUDE.md guard rails): precision is set in main(),
never at import; every integer read off a real number passes the absolute-
error criterion with >= max(20, dps/5) spare digits (>= 40 demanded for the
index); class groups, characters, kernels and projections are exact (HNF
arithmetic of schmidt_euler_system.py); PARI's class numbers and units are
labelled GRH-conditional unless bnfcertify returned 1; no PSLQ anywhere.

 (H)  the hyperbolic side (Theorem 3 of the document), odd n = 5..21: the units
      R_f of Paper II Thm 4.2 lie in the odd units E^- = ker(1 + sigma_r) of the
      ring class field H of disc 1 - n^2, and
        [E^- : <R_f>] = 24^(h/2) (2^(h/2-1)/Q^-) (h_H w_H+)/(h_H+ w_H) prod_{chi odd} C_chi(0),
      H^+ = H^{sigma_r}, Q^- = [E : E^+ E^-]; twisted group determinant, exact
      projections for the multipliers, PARI (tau from nfgaloisconj, kernels on
      the unit lattice, exact index, H^+), the regulator relation and the
      relative class number formula; the sextic layer at n = 21 (index 27648).

Usage:
    python3 scripts/robert_index_full.py --selftest            # Euclidean n <= 15 + hyperbolic (~2 min)
    python3 scripts/robert_index_full.py --selftest --with-23  # + Euclidean degree 24 (~+1 min)
    python3 scripts/robert_index_full.py 9                     # one Euclidean level
    python3 scripts/robert_index_full.py hyp 21                # one hyperbolic level
Requires mpmath, sympy (Smith form) and PARI/GP (`gp` on the PATH).
"""
import sys
import os
import re
import time
import shutil
import subprocess
import tempfile
from math import gcd
from fractions import Fraction
from itertools import product as iproduct

from mpmath import (mp, mpf, mpc, log, exp, pi, fabs, nstr, sqrt, e1, cos,
                    matrix, det, nint, floor)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from schmidt_euler_system import (ClassGroup, G_value, mass_log, factor, chi4,
                                  split_lams, Ne, extend, canon, index_of,
                                  is_primitive, vecs_of, hnf)

# precision is set in main(), never here (guard rail 2)

GP = shutil.which("gp")


# --------------------------------------------------------------------------
# certification helpers
# --------------------------------------------------------------------------

def need_digits():
    return max(20, mp.dps // 5)


def cert_int(x, what, need=None):
    """x (mpf/mpc) is certified as an integer: absolute error <= 10^-need."""
    need = need_digits() if need is None else need
    if isinstance(x, mpc):
        assert fabs(x.imag) < mpf(10) ** (-need), (what, "imaginary part", nstr(x.imag, 5))
        x = x.real
    k = int(nint(x))
    err = fabs(x - k)
    sp = int(-log(err, 10)) if err != 0 else mp.dps
    assert sp >= need, (what, nstr(x, 30), "spare", sp, "need", need)
    return k, sp


def spare_of(err):
    err = fabs(err)
    return mp.dps if err == 0 else int(-log(err, 10))


def poly_from_roots(roots):
    co = [mpc(1)]
    for r in roots:
        new = [mpc(0)] * (len(co) + 1)
        for i, c in enumerate(co):
            new[i] += c
            new[i + 1] -= c * r
        co = new
    return co


def cert_int_poly(roots, what):
    co = poly_from_roots(roots)
    out, worst = [], mp.dps
    for c in co:
        k, sp = cert_int(c, what)
        out.append(k)
        worst = min(worst, sp)
    return out, worst


# --------------------------------------------------------------------------
# abelian structure and characters of Pic(O_n)   (exact)
# --------------------------------------------------------------------------

def abelian_structure(cg):
    """(coords, orders): coordinates of every class in a direct-product
    decomposition Pic(O_n) = prod Z/o_i (greedy with direct-product test)."""
    h = cg.h
    ords = {c: cg.order(c) for c in cg.reps}
    known = {cg.one: ()}
    gens, orders = [], []
    while len(known) < h:
        found = False
        for g in sorted(cg.reps, key=lambda c: -ords[c]):
            if g in known:
                continue
            o = ords[g]
            new, ok, p = {}, True, cg.one
            for k in range(o):
                for c, co in known.items():
                    x = cg.prod(c, p) if k else c
                    if x in new:
                        ok = False
                        break
                    new[x] = co + (k,)
                if not ok:
                    break
                p = cg.prod(p, g)
            if ok and len(new) == len(known) * o:
                gens.append(g)
                orders.append(o)
                known = new
                found = True
                break
        assert found, "no direct generator found"
    ng = len(gens)
    coords = {c: tuple(list(co) + [0] * (ng - len(co))) for c, co in known.items()}
    return coords, orders


def all_characters(coords, orders):
    """list of (ks, order, phase dict class -> Fraction in [0,1))"""
    from math import lcm
    out = []
    for ks in iproduct(*[range(o) for o in orders]):
        cord = 1
        for k, o in zip(ks, orders):
            if k:
                cord = lcm(cord, o // gcd(k, o))
        ph = {}
        for c, co in coords.items():
            p = Fraction(0)
            for k, o, e in zip(ks, orders, co):
                p += Fraction(k * e, o)
            ph[c] = p % 1
        out.append((ks, cord, ph))
    return out


def cval(phase):
    """exp(2 pi i phase) at working precision (exact for phase 0, 1/2)"""
    if phase == 0:
        return mpc(1)
    if phase == Fraction(1, 2):
        return mpc(-1)
    return exp(2 * pi * mpc(0, 1) * mpf(phase.numerator) / phase.denominator)


# --------------------------------------------------------------------------
# the Delta-data and the Epstein L'(0,chi)
# --------------------------------------------------------------------------

def rep_numbers(L, M):
    """r[m] = #{x in Lambda : N(x) = m}, Lambda = Z d + Z(b + a i), m <= M"""
    d, b, a = L
    r = [0] * (M + 1)
    vmax = int((M ** 0.5) / a) + 2
    for v in range(-vmax, vmax + 1):
        im2 = (v * a) ** 2
        if im2 > M:
            continue
        rem = M - im2
        s = int(rem ** 0.5) + 2
        # |u d + v b| <= s
        for u in range((-v * b - s) // d - 1, (-v * b + s) // d + 2):
            re = u * d + v * b
            m = re * re + im2
            if 1 <= m <= M:
                r[m] += 1
    return r


def epstein_Lprime(cg, chars, verbose=False):
    """{ks: L'(0,chi)} for the nontrivial characters, by the independent
    incomplete-gamma evaluation (Paper II Lemma 2.3):
    L'(0,chi) = (1/2) sum_m R_chi(m) [e^{-alpha m}/(alpha m) + E_1(alpha m)],
    alpha = 2 pi / sqrt|D| = pi/n.  Inputs: exact representation numbers."""
    n = cg.n
    alpha = pi / n
    M = int((mp.dps + 15) * log(mpf(10)) / alpha) + 10
    reps = {c: rep_numbers(c, M) for c in cg.reps}
    Iv = [None] * (M + 1)
    for m in range(1, M + 1):
        am = alpha * m
        Iv[m] = exp(-am) / am + e1(am)
    out = {}
    for (ks, cord, ph) in chars:
        if cord == 1:
            continue
        tot = mpc(0)
        for c in cg.reps:
            rc = reps[c]
            s = mpf(0)
            for m in range(1, M + 1):
                if rc[m]:
                    s += rc[m] * Iv[m]
            tot += cval(ph[c]) * s
        out[ks] = tot / 2
    if verbose:
        print(f"    (independent Epstein evaluation: {M} terms)")
    return out


def delta_data(cg):
    """{c: G_c} (mpc) and {c: log|G_c|}"""
    G = {c: G_value(c, cg.n) for c in cg.reps}
    f = {c: log(fabs(G[c])) for c in cg.reps}
    return G, f


def mass_int(n):
    """|M(n)| as an exact integer (Paper I, Thm 7.18)"""
    M = 1
    for p, k in factor(n).items():
        if chi4(p) == 1:
            continue
        e = 2 if p == 2 else 1
        M *= p ** ((6 // e) * (p ** k - 1) // (p - 1) * Ne(n // p ** k))
    return M


# --------------------------------------------------------------------------
# V1: the regulator of V_n as a Dedekind determinant
# --------------------------------------------------------------------------

def regulator_matrix(cg, f):
    """M_{a,c} = 2 (f(a^{-1} c) - f(a^{-1})),  a, c != 1  (rows: the complex
    place of the embedding sigma_a; columns: the unit G_c/G_1)"""
    others = [c for c in cg.reps if c != cg.one]
    Mx = matrix(len(others), len(others))
    for i, a in enumerate(others):
        ainv = cg.inv(a)
        base = f[ainv]
        for j, c in enumerate(others):
            Mx[i, j] = 2 * (f[cg.prod(ainv, c)] - base)
    return Mx


# --------------------------------------------------------------------------
# V3: primitive levels and the multipliers C_chi(0)
# --------------------------------------------------------------------------

def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def factors_through(cg, ph, m):
    """does the character with phases ph factor through Pic(O_m)?  (m | n)"""
    if m == 1:
        return all(p == 0 for p in ph.values())
    img_one = cg.project(cg.one, m)
    for c in cg.reps:
        if cg.project(c, m) == img_one and ph[c] != 0:
            return False
    return True


def primitive_level(cg, ph):
    """minimal m | n through which chi factors; asserts gcd-closure"""
    levels = [m for m in divisors(cg.n) if factors_through(cg, ph, m)]
    for m1 in levels:
        for m2 in levels:
            assert gcd(m1, m2) in levels, ("levels not gcd-closed", cg.n, levels)
    return min(levels), levels


def character_on_level(cg, ph, m):
    """chi as a phase dict on Pic(O_m) (chi factors through m); returns
    (ClassGroup(m), phases)"""
    cgm = ClassGroup(m) if m >= 2 else None
    if cgm is None:
        return None, None
    out = {}
    for c in cg.reps:
        cm = cg.project(c, m)
        if cm in out:
            assert out[cm] == ph[c]
        else:
            out[cm] = ph[c]
    assert set(out) == set(cgm.reps)
    return cgm, out


def multiplier_chain(cg, ph, m):
    """C_chi(0) = prod_l a_{k_l}(chi) for the pullback from the primitive
    level m to n (schmidt-euler-system.md Thm 2 iterated):
       a_0 = 1,  a_1 = P_l(chi) [l not | m] or l + 1 [l | m],
       a_{k+1} = (l+1) a_k - m_k a_{k-1},  m_1 = l - chi_{-4}(l) [l not | m] or l [l | m],
       m_k = l (k >= 2).
    Returns (value as mpf, description)."""
    n = cg.n
    fac = factor(n // m)
    C = mpf(1)
    desc = []
    for l, k in sorted(fac.items()):
        mprime = n // l ** k          # the level before l is adjoined
        if l % mprime == 0 and mprime % l == 0:
            pass
        if mprime % l == 0:           # l | m: chi primitive at m w.r.t. l
            a1 = mpf(l + 1)
            m1 = l
            case = f"l={l} | m: a_1 = l+1 = {l+1}"
        else:
            lams = split_lams(l)
            if lams == []:
                a1 = mpf(l + 1)
                case = f"l={l} inert: a_1 = {l+1}"
            else:
                cgp, php = character_on_level(cg, ph, mprime)
                if cgp is None:
                    # mprime = 1: chi trivial -- cannot happen for chi != 1
                    raise AssertionError
                z = cval(php[cgp.ideal_class(lams[0])])
                if l == 2:
                    a1 = (3 - z).real
                    case = f"l=2: a_1 = 3 - chi(p) = {nstr(a1, 8)}"
                else:
                    a1 = (l + 1 - z - 1 / z).real
                    case = f"l={l} split: a_1 = l+1-2Re chi(l) = {nstr(a1, 8)}"
            m1 = l - chi4(l)
        a_prev, a_cur = mpf(1), a1
        for j in range(2, k + 1):
            mk = m1 if j == 2 else l
            a_prev, a_cur = a_cur, (l + 1) * a_cur - mk * a_prev
        C *= a_cur
        desc.append(case + (f", a_{k} = {nstr(a_cur, 8)}" if k > 1 else ""))
    return C, "; ".join(desc) if desc else "primitive"


# --------------------------------------------------------------------------
# PARI/GP
# --------------------------------------------------------------------------

def gp_version():
    if GP is None:
        return None
    out = subprocess.run([GP, "--version-short"], capture_output=True, text=True)
    return (out.stdout + out.stderr).strip()


def run_gp(script, timeout=3600):
    assert GP is not None, "PARI/GP (`gp`) not found on the PATH"
    with tempfile.NamedTemporaryFile("w", suffix=".gp", delete=False) as fh:
        fh.write(script + "\nquit;\n")
        path = fh.name
    try:
        out = subprocess.run([GP, "-q", "-f", path], capture_output=True, text=True,
                             timeout=timeout, stdin=subprocess.DEVNULL)
    finally:
        os.unlink(path)
    errs = [ln for ln in (out.stdout + out.stderr).splitlines()
            if "***" in ln and "Warning" not in ln]
    if out.returncode != 0 or errs:
        raise RuntimeError("gp failed:\n" + out.stdout + out.stderr)
    return out.stdout


def parse_int_vector(s):
    s = s.strip()
    assert s.startswith("[") and s.endswith("]"), s
    body = s[1:-1].strip()
    if not body:
        return []
    return [int(t) for t in body.replace(",", " ").split()]


def parse_int_matrix(s):
    """PARI matrix printed as [a, b; c, d] -> list of rows (also Mat(a), Mat([a, b]))"""
    s = s.strip()
    m = re.match(r"^Mat\((.*)\)$", s)
    if m:
        inner = m.group(1).strip()
        if inner.startswith("["):
            return [parse_int_vector(inner)]
        return [[int(inner)]]
    assert s.startswith("[") and s.endswith("]"), s
    body = s[1:-1].strip()
    if not body:
        return []
    return [[int(t) for t in row.replace(",", " ").split()] for row in body.split(";")]


GP_HEADER = """default(parisize, 2000000000);
default(realprecision, {prec0});
"""


def gp_field_script(n, prec, certdeg, Dn, extra=""):
    """gp script: H_n = polredbest(polcompositum(y^2+1, polclass(-4n^2))),
    bnfinit (with units) at `prec` digits; class number, cyclic structure,
    torsion, regulator, bnfcertify (degree <= certdeg), the roots of D_n in
    H_n and the bnfisunit exponent matrix of the units G_c/G_1."""
    return GP_HEADER.format(prec0=prec) + f"""
n = {n};
t0 = getabstime();
P = polclass(-4*n^2, 0, 'y);
Q = polredbest(polcompositum(y^2+1, P)[1]);
print("TPOL ", (getabstime()-t0)/1000.);
t0 = getabstime();
bnf = bnfinit(Q, 1);
print("TBNF ", (getabstime()-t0)/1000.);
print("POL ", Q);
print("DEG ", poldegree(Q));
print("H ", bnf.no);
print("CYC ", bnf.cyc);
print("W ", bnf.tu[1]);
print("SQRT2 ", #nfroots(bnf.nf, x^2-2));
print("SQRT3 ", #nfroots(bnf.nf, x^2-3));
t0 = getabstime();
cert = if(poldegree(Q) <= {certdeg}, bnfcertify(bnf), -1);
print("CERT ", cert);
print("TCERT ", (getabstime()-t0)/1000.);
Dn = Pol({Dn}, x);
r = nfroots(bnf.nf, Dn);
print("NROOTS ", #r);
E = matrix(#r-1, #bnf.fu);
T = vector(#r-1);
for(i=2, #r, e = bnfisunit(bnf, r[i]/r[1]); for(j=1, #bnf.fu, E[i-1,j] = e[j]); T[i-1] = lift(e[#e]));
print("EXP ", E);
print("TORS ", T);
print("EDET ", matdet(E));
print("SNF ", matsnf(E));
{extra}
print("REG ", bnf.reg);
"""


def gp_layers_script(n, prec, cubic_theta, M13, cubic_theta_u, quad):
    """the cubic layer (L_3 = Q(theta_u), F = K(theta_u)) and the quadratic
    layer (Q(theta^(2))) built directly from the certified coset polynomials.
    cubic_theta: prod (x - theta_i) (theta_i the coset products);
    cubic_theta_u: prod (x - theta_i/M^(1/3)); quad: x^2 - t x +- 1 of the
    quadratic-layer unit."""
    s = GP_HEADER.format(prec0=prec)
    if cubic_theta_u is not None:
        s += f"""
cu = Pol({cubic_theta_u}, y);
L3 = polredbest(cu);
print("L3POL ", L3);
bL = bnfinit(L3, 1);
print("L3H ", bL.no);
print("L3REG ", bL.reg);
print("L3FU ", lift(bL.fu[1]));
print("L3CERT ", bnfcertify(bL));
rt = nfroots(bL.nf, subst(cu, y, x));
print("L3NROOTS ", #rt);
e = bnfisunit(bL, rt[1]);
print("L3EXP ", e[1], " ", lift(e[2]));
F = polredbest(polcompositum(y^2+1, cu)[1]);
print("FPOL ", F);
bF = bnfinit(F, 1);
print("FH ", bF.no);
print("FREG ", bF.reg);
print("FW ", bF.tu[1]);
print("FCERT ", bnfcertify(bF));
rF = nfroots(bF.nf, subst(cu, y, x));
print("FNROOTS ", #rF);
for(i=1, #rF, e = bnfisunit(bF, rF[i]); print("FEXP", i, " ", e[1], " ", e[2], " ", lift(e[3])));
"""
    if quad is not None:
        s += f"""
qu = Pol({quad}, y);
K2 = polredbest(qu);
print("K2POL ", K2);
bQ = bnfinit(K2, 1);
print("K2FU ", lift(bQ.fu[1]));
print("K2REG ", bQ.reg);
print("K2DISC ", nfdisc(K2));
rq = nfroots(bQ.nf, subst(qu, y, x));
print("K2NROOTS ", #rq);
e = bnfisunit(bQ, rq[1]);
print("K2EXP ", e[1], " ", lift(e[2]));
"""
    return s


def parse_tagged(out):
    d = {}
    for line in out.splitlines():
        line = line.strip()
        m = re.match(r"^([A-Z0-9]+) (.*)$", line)
        if m:
            d[m.group(1)] = m.group(2).strip()
    return d


# --------------------------------------------------------------------------
# the level computation
# --------------------------------------------------------------------------

def smith_invariants(rows):
    """Smith normal form invariants of an integer matrix (sympy)"""
    from sympy import Matrix, ZZ
    from sympy.matrices.normalforms import smith_normal_form
    S = smith_normal_form(Matrix(rows), domain=ZZ)
    inv = [abs(int(S[i, i])) for i in range(min(S.shape))]
    return [x for x in inv if x != 1]


def level(n, prec_pari=120, certdeg=16, verbose=True, layers=True):
    say = print if verbose else (lambda *a, **k: None)
    t0 = time.time()
    cg = ClassGroup(n)
    h = cg.h
    coords, orders = abelian_structure(cg)
    chars = all_characters(coords, orders)
    say("=" * 78)
    say(f"LEVEL n = {n}:  Pic(O_n) = {' x '.join('Z/%d' % o for o in orders)}, h = {h}, "
        f"[H_n : Q] = {2 * h}, unit rank {h - 1}   (dps = {mp.dps})")
    say("=" * 78)
    need = need_digits()
    rec = {"n": n, "h": h, "orders": orders}

    # ---- the Delta-data, D_n, the regulator determinant (V1) ----
    G, f = delta_data(cg)
    Dn, spD = cert_int_poly([G[c] for c in cg.reps], f"D_{n}")
    assert abs(Dn[-1]) == mass_int(n), (n, Dn[-1], mass_int(n))
    say(f"D_n certified (integer coefficients, spare {spD}; D_n(0) = (-1)^h M(n), "
        f"M(n) = {mass_int(n)})")
    S = {}
    for (ks, cord, ph) in chars:
        if cord == 1:
            continue
        s = sum(cval(ph[c]) * f[c] for c in cg.reps)
        assert fabs(s.imag) < mpf(10) ** (-need), (n, ks, "S_chi not real")
        S[ks] = s.real
        assert fabs(S[ks]) > 1, (n, ks, "S_chi vanishes?")
    prodS = mpf(1)
    for ks in S:
        prodS *= S[ks]
    Mx = regulator_matrix(cg, f)
    dM = det(Mx) if h > 1 else mpf(1)
    rhs = 2 ** (h - 1) * prodS
    spV1 = spare_of((dM - rhs) / max(1, fabs(rhs)))
    assert spV1 >= need, (n, nstr(dM, 20), nstr(rhs, 20))
    RV = fabs(dM)
    say(f"(V1) det(2(f(a^-1 c) - f(a^-1)))_(a,c != 1) = 2^(h-1) prod_chi S_chi: "
        f"relative spare {spV1};  R(V_n) = {nstr(RV, 25)}")

    # ---- the independent L'(0,chi) (V2) and the multipliers (V3) ----
    Lp = epstein_Lprime(cg, chars)
    worst2 = mp.dps
    prodC = mpf(1)
    Cinfo = {}
    prodL = mpf(1)
    for (ks, cord, ph) in chars:
        if cord == 1:
            continue
        l = Lp[ks]
        assert fabs(l.imag) < mpf(10) ** (-need), (n, ks)
        l = l.real
        sp = spare_of((S[ks] + 12 * l) / max(1, fabs(l)))
        worst2 = min(worst2, sp)
        prodL *= fabs(l)
        m, levels = primitive_level(cg, ph)
        C, desc = multiplier_chain(cg, ph, m)
        # check against the independent evaluation at the primitive level
        if m < n:
            cgm, phm = character_on_level(cg, ph, m)
            coords_m, orders_m = abelian_structure(cgm)
            Lm = epstein_Lprime(cgm, [((), 2, phm)])[()]
            spC = spare_of((l - C * Lm.real) / max(1, fabs(l)))
            assert spC >= need, (n, ks, m, nstr(l / Lm.real, 20), nstr(C, 20))
        else:
            spC = None
        Cinfo[ks] = (cord, m, C, desc, spC)
        prodC *= C
    assert worst2 >= need, (n, worst2)
    say(f"(V2) S_chi = -12 L'(0,chi) against the independent Epstein evaluation, all "
        f"{len(S)} characters: relative spare >= {worst2}")
    prodC_int, spPC = cert_int(prodC, "prod C_chi(0)")
    say(f"(V3) primitive levels and multipliers C_chi(0):")
    for ks, (cord, m, C, desc, spC) in sorted(Cinfo.items()):
        tail = "" if spC is None else f"   [L'(0,chi)/L'(0,chi^(m)) certified, spare {spC}]"
        say(f"      chi{ks} (order {cord}): primitive level {m}, C = {nstr(C, 10)}  {desc}{tail}")
    say(f"      prod_chi C_chi(0) = {prodC_int}   (spare {spPC})")
    rec.update({"D": Dn, "RV": RV, "prodL": prodL, "prodC": prodC_int, "S": S, "Lp": Lp})

    # ---- PARI: the field, h, R, w; the exact exponent matrix (V4, V5) ----
    out = run_gp(gp_field_script(n, prec_pari, certdeg, Dn))
    d = parse_tagged(out)
    hH = int(d["H"])
    w = int(d["W"])
    cyc = parse_int_vector(d["CYC"])
    cert = int(d["CERT"])
    Rh = mpf(d["REG"])
    deg = int(d["DEG"])
    assert deg == 2 * h, (n, deg, h)
    assert int(d["NROOTS"]) == h, (n, d["NROOTS"])
    assert w % 4 == 0
    sq2, sq3 = int(d["SQRT2"]) > 0, int(d["SQRT3"]) > 0
    w_pred = 4 * (2 if sq2 else 1) * (3 if sq3 else 1)
    assert w == w_pred, (n, w, sq2, sq3)
    grh = "unconditional (bnfcertify = 1)" if cert == 1 else "GRH-conditional (bnfcertify not run)"
    say(f"(V4) PARI/GP: H_n = Q[y]/({d['POL']})")
    say(f"      h(H_n) = {hH}, Cl = {cyc if cyc else '[trivial]'}, w(H_n) = {w} "
        f"(sqrt2 in H_n: {sq2}, sqrt3 in H_n: {sq3}), R(H_n) = {nstr(Rh, 30)}   -- {grh}")
    index_real = RV / Rh
    spare_need = max(40, need)
    index_int, spI = cert_int(index_real, f"index at n={n}", need=spare_need)
    say(f"      24^(h-1) prod|L'(0,chi)| / R(H_n) = {index_int}   (certified integer, spare {spI} >= {spare_need})")
    pred = Fraction(4 * 24 ** (h - 1), w) * hH * prodC_int
    assert pred.denominator == 1 and pred.numerator == index_int, (n, index_int, pred)
    say(f"      = (4 * 24^(h-1) / w) h(H_n) prod C_chi(0) = (4 * 24^{h-1} / {w}) * {hH} * {prodC_int}: OK")
    # exact index from the exponent matrix
    E = parse_int_matrix(d["EXP"]) if h > 1 else []
    edet = int(d["EDET"]) if h > 1 else 1
    assert abs(edet) == index_int, (n, edet, index_int)
    snf = smith_invariants(E) if h > 1 else []
    tors = parse_int_vector(d["TORS"]) if h > 1 else []
    say(f"(V5) exact: units G_c/G_1 in the fundamental-unit basis (bnfisunit): |det| = {abs(edet)} "
        f"= index;  O_H^x / mu V_n = {' x '.join('Z/%d' % x for x in snf) if snf else 'trivial'}")
    say(f"      torsion components of G_c/G_1 (exponent of the generator of mu, mod w): {tors}")
    # 24th-root saturation W_n = {u : u^24 in mu V_n}: [W_n : mu V_n] = prod gcd(d_i, 24)
    if h > 1:
        gE = 0
        for row in E:
            for x in row:
                gE = gcd(gE, abs(x))
        sat = 1
        quo = []
        for dd in snf:
            sat *= gcd(dd, 24)
            if dd // gcd(dd, 24) > 1:
                quo.append(dd // gcd(dd, 24))
        # invariants equal to 1 were dropped from snf; they contribute gcd 1
        sat_index = index_int // sat
        assert index_int % sat == 0
        say(f"      saturation: every G_c/G_1 lies in mu (O_H^x)^{gE}; [W_n : mu V_n] = prod gcd(d_i, 24) = {sat}; "
            f"[O_H^x : W_n] = {sat_index}, O_H^x/W_n = {' x '.join('Z/%d' % x for x in quo) if quo else 'trivial'}"
            f"   (h(H_n) prod C = {hH * prodC_int}: {'EQUAL' if sat_index == hH * prodC_int else 'DIFFERENT'})")
    else:
        gE, sat, sat_index, quo = 1, 1, 1, []
    say(f"      PARI timings: polynomial {d['TPOL']} s, bnfinit {d['TBNF']} s, bnfcertify {d['TCERT']} s")
    rec.update({"hH": hH, "w": w, "cyc": cyc, "cert": cert, "RH": Rh, "index": index_int,
                "snf": snf, "pol": d["POL"], "E": E, "tors": tors, "spI": spI,
                "gE": gE, "sat": sat, "sat_index": sat_index, "quo": quo})

    # ---- the layers (V6) ----
    if layers and h > 1:
        layer_results = layers_of(cg, chars, G, f, S, n, prec_pari, say, Cinfo)
        rec["layers"] = layer_results
    say(f"level {n} done in {time.time() - t0:.1f} s")
    return rec


def factorize_int(m):
    return factor(m)


def layers_of(cg, chars, G, f, S, n, prec, say, Cinfo):
    """quadratic layer (all real characters) and cubic layer (if 3 | h)."""
    out = {}
    need = need_digits()
    # quadratic layer: for each real chi_2, theta^(2) = prod_{ker} G / |M|^(1/2)
    Mn = mass_int(n)
    r = int(round(Mn ** 0.5))
    while r * r < Mn:
        r += 1
    while r * r > Mn:
        r -= 1
    if r * r == Mn:
        for (ks, cord, ph) in chars:
            if cord != 2:
                continue
            ker = [c for c in cg.reps if ph[c] == 0]
            oth = [c for c in cg.reps if ph[c] != 0]
            th = mpc(1)
            for c in ker:
                th *= G[c]
            th2 = mpc(1)
            for c in oth:
                th2 *= G[c]
            u1, u2 = th / r, th2 / r
            co, sp = cert_int_poly([u1, u2], f"quadratic layer chi{ks}")
            assert abs(co[-1]) == 1, (n, ks, co)
            # PARI: Q(theta^(2)), fundamental unit, exponent
            d = parse_tagged(run_gp(gp_layers_script(n, prec, None, None, None, co)))
            e2 = int(d["K2EXP"].split()[0])
            disc2 = int(d["K2DISC"])
            say(f"(V6) quadratic layer chi{ks}: theta^(2) = prod_ker G / |M|^(1/2) is a unit of "
                f"Q(sqrt {disc2}), minimal polynomial {co}; theta^(2) = +-eps^{e2} "
                f"(eps = {d['K2FU']}, fundamental unit certified by PARI): index {abs(e2)}")
            # consistency with the character sum: log|theta^(2)| = S_chi / 2
            sp2 = spare_of((log(fabs(u1)) - S[ks] / 2) / max(1, fabs(S[ks])))
            assert sp2 >= need
            out[("quad", ks)] = (co, disc2, e2)
    # cubic layer
    cub = [(ks, ph) for (ks, cord, ph) in chars if cord == 3]
    if cub:
        ks, ph = cub[0]
        cosets = {}
        for c in cg.reps:
            cosets.setdefault(ph[c], []).append(c)
        assert len(cosets) == 3
        prods = []
        for key in sorted(cosets):
            p = mpc(1)
            for c in cosets[key]:
                p *= G[c]
            prods.append(p)
        co_theta, sp = cert_int_poly(prods, "coset cubic")
        M13 = int(round(Mn ** (1.0 / 3)))
        while M13 ** 3 < Mn:
            M13 += 1
        while M13 ** 3 > Mn:
            M13 -= 1
        assert M13 ** 3 == Mn, (n, Mn)
        co_u, sp_u = cert_int_poly([p / M13 for p in prods], "theta_u cubic")
        assert abs(co_u[-1]) == 1, (n, co_u)
        d = parse_tagged(run_gp(gp_layers_script(n, prec, co_theta, M13, co_u, None)))
        eL = int(d["L3EXP"].split()[0])
        hL = int(d["L3H"])
        RL = mpf(d["L3REG"])
        RF = mpf(d["FREG"])
        hF = int(d["FH"])
        wF = int(d["FW"])
        certL = int(d["L3CERT"]) == 1 and int(d["FCERT"]) == 1
        assert int(d["FNROOTS"]) == 3
        rows = [[int(t) for t in d[f"FEXP{i}"].split()] for i in (1, 2, 3)]
        assert all(rows[0][j] + rows[1][j] + rows[2][j] == 0 for j in (0, 1)), rows
        idxF0 = abs(rows[0][0] * rows[1][1] - rows[0][1] * rows[1][0])
        # the ratio group V^A = N_{H_n/F}(V_n) = <theta_u^(b)/theta_u^(b')> has index 3
        # in mu<theta_u, theta_u'> (the three conjugates multiply to +-1)
        d1 = [rows[0][j] - rows[1][j] for j in (0, 1)]
        d2 = [rows[1][j] - rows[2][j] for j in (0, 1)]
        idxF = abs(d1[0] * d2[1] - d1[1] * d2[0])
        assert idxF == 3 * idxF0
        C3 = Cinfo[ks][2]
        C3i, _ = cert_int(C3, "C_chi3(0)")
        assert wF == 4, wF
        predF = Fraction(4 * 24 ** 2, wF) * hF * C3i ** 2
        assert predF == idxF, (n, idxF, predF)
        # the same theorem for the real cubic layer: 8 h_L3 C_n(0) (Paper II, Thm 6.6)
        say(f"(V6) cubic layer chi{ks}: coset cubic {co_theta}; theta_u = theta/{M13}: minimal polynomial {co_u}")
        say(f"      L_3 = Q(theta_u) = Q[y]/({d['L3POL']}): h = {hL}, R = {nstr(RL, 20)}, fundamental unit "
            f"{d['L3FU']}; theta_u = +-eta^{eL}  ->  [O_L3^x : <-1, theta_u>] = {abs(eL)} = 8 h_L3 C_n(0)"
            f"   ({'unconditional' if certL else 'GRH'})")
        say(f"      F = K(theta_u) = Q[y]/({d['FPOL']}) (totally complex, NOT CM: unit rank 2): h = {hF}, w = {wF}, "
            f"R = {nstr(RF, 20)}; exponent vectors of the three conjugates of theta_u: {rows};")
        say(f"      [O_F^x : mu_F <theta_u, theta_u'>] = {idxF0};  ratio group V^A = N_(H_n/F)(V_n): "
            f"[O_F^x : mu_F V^A] = 3 * {idxF0} = {idxF} = (4 24^2/w_F) h_F C_n(0)^2 = 576 * {hF} * {C3i}^2: OK")
        out["cubic"] = (co_theta, co_u, hL, eL, hF, wF, idxF, d["L3POL"], d["FPOL"])
    return out


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

SELFTEST_LEVELS = [3, 5, 7, 9, 11, 13, 15]
HYP_LEVELS = [5, 7, 9, 11, 13, 15, 17, 19, 21]

# hyperbolic regression record: n -> (h_H, h_H+, w_H, w_H+, Q^-, prod C, index)
HYP_RECORD = {
    5: (1, 2, 6, 2, 1, 1, 4),
    7: (1, 1, 12, 6, 1, 1, 12),
    9: (1, 1, 4, 4, 2, 1, 576),
    11: (2, 4, 6, 2, 2, 1, 96),
    13: (1, 2, 6, 6, 2, 1, 288),
    15: (1, 2, 8, 2, 8, 16, 663552),
    17: (1, 1, 24, 8, 2, 4, 768),
    19: (1, 1, 6, 6, 4, 1, 663552),
    21: (4, 4, 2, 2, 32, 1, 191102976),
}

# regression record: n -> (h(H_n), w, prod C, index, Smith invariants of E_n)
RECORD = {
    3: (1, 12, 1, 8, [8]),
    5: (1, 4, 1, 24, [24]),
    7: (1, 4, 1, 13824, [24, 24, 24]),
    9: (1, 12, 4, 10616832, [96, 24, 24, 24, 8]),
    11: (1, 4, 1, 7962624, [24, 24, 24, 24, 24]),
    13: (3, 4, 1, 23887872, [72, 24, 24, 24, 24]),
    15: (2, 12, 32, 97844723712, [384, 48, 48, 24, 24, 24, 8]),
    23: (12, 4, 1, 18260173718028288, [144, 48] + [24] * 9),
}


def selftest(with_23=False):
    t0 = time.time()
    print("=" * 78)
    print("robert_index_full.py --selftest   (mpmath dps = %d; PARI/GP %s)" % (mp.dps, gp_version()))
    print("=" * 78)
    levels = SELFTEST_LEVELS + ([23] if with_23 else [])
    recs = []
    for n in levels:
        r = level(n)
        hH, w, pC, idx, snf = RECORD[n]
        assert (r["hH"], r["w"], r["prodC"], r["index"]) == (hH, w, pC, idx), (n, r["hH"], r["w"], r["prodC"], r["index"])
        assert sorted(r["snf"]) == sorted(snf), (n, r["snf"])
        assert r["sat_index"] == hH * pC, (n, r["sat_index"])
        recs.append(r)
    hrecs = []
    for n in HYP_LEVELS:
        r = hyperbolic_odd_index(n)
        hH, hHp, wH, wHp, Qm, pC, idx = HYP_RECORD[n]
        assert (r["hH"], r["hHp"], r["wH"], r["wHp"], r["Qm"], r["prodC"], r["idx"]) == (hH, hHp, wH, wHp, Qm, pC, idx), (n, r)
        hrecs.append(r)
    if 21 in HYP_LEVELS:
        r21 = hrecs[HYP_LEVELS.index(21)]
        assert r21["sextic"][6] == 27648 and r21["sextic"][1] == 8 and r21["sextic"][2] == 4, r21["sextic"]
    print("=" * 78)
    print("SUMMARY   index = [O_{H_n}^x : mu(H_n) V_n] = 24^(h-1) prod|L'| / R_H = (4 24^(h-1)/w) h(H_n) prod C")
    print("=" * 78)
    print(f"{'n':>3} {'Pic':>10} {'h':>3} {'w':>3} {'h(H_n)':>7} {'Cl(H_n)':>9} {'prodC':>6} "
          f"{'index':>28} {'factored':>22} {'spare':>6} status")
    for r in recs:
        fac = " ".join(f"{p}^{k}" if k > 1 else f"{p}" for p, k in sorted(factor(r['index']).items())) if r['index'] > 1 else "1"
        st = "certified" if r["cert"] == 1 else "GRH"
        print(f"{r['n']:>3} {'x'.join('Z/%d' % o for o in r['orders']):>10} {r['h']:>3} {r['w']:>3} "
              f"{r['hH']:>7} {str(r['cyc'] or '[]'):>9} {r['prodC']:>6} {r['index']:>28} {fac:>22} {r['spI']:>6} {st}")
    print("=" * 78)
    print("HYPERBOLIC   [E^- : <R_f>] = 24^(h/2) (2^(h/2-1)/Q^-) (h_H w_H+)/(h_H+ w_H) prod_(chi odd) C   (Theorem 3)")
    print("=" * 78)
    print(f"{'n':>3} {'D':>6} {'Cl':>10} {'h_H':>4} {'w_H':>4} {'h_H+':>5} {'w_H+':>5} {'Q^-':>4} {'prodC':>6} {'index':>12} {'index/24^(h/2)':>15} status")
    for r in hrecs:
        st = "certified" if r["cert"] else "GRH"
        print(f"{r['n']:>3} {1 - r['n'] ** 2:>6} {'x'.join('Z/%d' % o for o in r['orders']):>10} {r['hH']:>4} {r['wH']:>4} "
              f"{r['hHp']:>5} {r['wHp']:>5} {r['Qm']:>4} {r['prodC']:>6} {r['idx']:>12} {str(Fraction(r['idx'], 24 ** (r['h'] // 2))):>15} {st}")
    print(f"ALL CHECKS PASSED   ({time.time() - t0:.0f}s)")
    return recs


def main(argv):
    mp.dps = 150          # precision set here, never at import (guard rail 2)
    args = [a for a in argv[1:] if not a.startswith("--")]
    flags = [a for a in argv[1:] if a.startswith("--")]
    if GP is None:
        print("PARI/GP (`gp`) is required: apt-get install -y pari-gp")
        sys.exit(1)
    if "--selftest" in flags or not args:
        selftest(with_23="--with-23" in flags)
    elif args[0] == "hyp":
        for a in args[1:]:
            hyperbolic_odd_index(int(a))
    else:
        for a in args:
            level(int(a))



# --------------------------------------------------------------------------
# the hyperbolic side: the odd-unit index of the R_f  (Theorem 3 of the document)
# --------------------------------------------------------------------------

def fundamental_part(D):
    """D = f^2 dK with dK a fundamental discriminant; returns (dK, f)"""
    f = 1
    m = -D
    # largest square dividing D such that the cofactor is a discriminant
    k = 2
    while k * k <= m:
        while D % (k * k) == 0 and ((D // (k * k)) % 4 in (0, 1)):
            D //= k * k
            f *= k
        k += 1
    # D may still be 4 * (1 mod 4) = a non-fundamental with f 2
    if D % 4 == 0 and (D // 4) % 4 == 1:
        D //= 4
        f *= 2
    return D, f


def project_form(f, D, dK, fp):
    """the class of the primitive form f of discriminant D = fc^2 dK in the
    class group of the order of conductor fp (fp | fc): the reduced form of the
    ideal O_{fp} * [a, (-b + sqrt D)/2].  Elements of K' = Q(sqrt dK) are
    pairs (x, y) = x + y sqrt(dK) of Fractions."""
    from euclidean_moduli_invariants import reduce_form
    a, b, c = f
    fc2 = D // dK
    fc = int(round(fc2 ** 0.5))
    assert fc * fc == fc2 and fc % fp == 0
    alpha = (Fraction(a), Fraction(0))
    beta = (Fraction(-b, 2), Fraction(fc, 2))
    om = (Fraction(1, 2), Fraction(1, 2)) if dK % 4 == 1 else (Fraction(0), Fraction(1, 2))
    omf = (om[0] * fp, om[1] * fp)

    def mulK(u, v):
        return (u[0] * v[0] + u[1] * v[1] * dK, u[0] * v[1] + u[1] * v[0])

    gens = [alpha, beta, mulK(omf, alpha), mulK(omf, beta)]
    den = 1
    for (x, y) in gens:
        den = den * x.denominator // gcd(den, x.denominator)
        den = den * y.denominator // gcd(den, y.denominator)
    rows = [(int(x * den), int(y * den)) for (x, y) in gens]
    d_, b_, a_ = hnf(rows)               # lattice Z (d_, 0) + Z (b_, a_)
    g1 = (Fraction(d_, den), Fraction(0))
    g2 = (Fraction(b_, den), Fraction(a_, den))
    # orientation: Im(g2/g1) > 0  <=>  x1 y2 - x2 y1 > 0
    if g1[0] * g2[1] - g2[0] * g1[1] < 0:
        g2 = (-g2[0], -g2[1])

    def N(u):
        return u[0] * u[0] - dK * u[1] * u[1]
    A = N(g1)
    C = N(g2)
    B = 2 * (g1[0] * g2[0] - dK * g1[1] * g2[1])
    covol = abs(g1[0] * g2[1] - g2[0] * g1[1])
    covol_O = abs(Fraction(1) * omf[1])       # det of (1, 0), (omf)
    Nideal = covol / covol_O
    A, B, C = A / Nideal, B / Nideal, C / Nideal
    assert A.denominator == 1 and B.denominator == 1 and C.denominator == 1, (f, A, B, C)
    A, B, C = int(A), int(B), int(C)
    assert B * B - 4 * A * C == fp * fp * dK, (f, A, B, C)
    g = gcd(gcd(A, B), C)
    assert g == 1, (f, A, B, C)
    return reduce_form(A, B, C)


def hyperbolic_odd_index(n, prec_pari=120, verbose=True, layer=True):
    """Level n odd, D = 1 - n^2, K' = Q(sqrt D), H the ring class field of the
    order of discriminant D, r the twist class, tau = sigma_r, H^+ = H^tau.
    The units R_f of Paper II Thm 4.2 satisfy R_{rf} = 1/R_f, so they lie in
    the odd units E^- = {u : tau(u) u in mu} of rank h/2.  Verified:
      (H1) the twisted Dedekind determinant det(2 log|R_{a^-1 f}|)_{a,f in T}
           (T a transversal of <r>) equals prod_{chi odd} S_chi, S_chi =
           sum chi(f) log|R_f| = -24 L'(0,chi) (independent evaluation);
      (H2) the multipliers: for every odd chi the primitive conductor f' | f
           (D = f^2 dK) is found exactly through the projections
           Cl(D) -> Cl(f'^2 dK) (ideal extension, in exact arithmetic), and
           C_chi(0) = L'(0,chi) / L'(0,chi_{f'}) with both L' from the
           independent evaluation; prod C certified as an integer;
      (H3) PARI: H, h_H, w_H, R_H; the h units R_f (nfroots of the certified
           R-polynomial, bnfisunit); tau as the unique automorphism of H/K'
           inverting every R_f; E^- = ker(1 + tau) on the unit lattice; the
           EXACT index [E^- : mu <R_f>]; H^+ (degree h over Q) with h_{H^+},
           R_{H^+}, w_{H^+}; Q^- = [E : E^+ E^-] with E^+ the image of
           O_{H^+}^x (the tau-fixed units modulo the ones fixed only up to a
           root of unity); the regulator relation
           R_H / R_{H^+} = 2^{h/2-1} R^-(E^-) / Q^-; and the identity
             index = 24^{h/2} (2^{h/2-1}/Q^-) (h_H w_{H^+})/(h_{H^+} w_H) prod_{chi odd} C_chi(0);
      (H4) the sextic layer (odd chi of order 6 with kernel A): the coset units
           u_b = N_{H/F}(R_f), F = H^A, and the same identity in F.
    Class numbers and units are GRH-conditional unless bnfcertify succeeded
    (attempted at degree <= 16)."""
    from involution_classmap import classes_of_disc, is_primitive as prim_form
    from euclidean_moduli_invariants import reduce_form
    import schmidt_units as SU
    import phase_klf as PK
    say = print if verbose else (lambda *a, **k: None)
    need = need_digits()
    t0 = time.time()
    D = 1 - n * n
    dK, fc = fundamental_part(D)
    r0, s0 = (n - 1) // 2, (n + 1) // 2
    prim = [f for f in classes_of_disc(D) if prim_form(f)]
    h = len(prim)
    say("=" * 78)
    say(f"HYPERBOLIC LEVEL n = {n}:  D = {D} = {fc}^2 ({dK}), h = {h}, [H : Q] = {2 * h}   (dps = {mp.dps})")
    say("=" * 78)
    R = {f: SU.R_lattice(f, n) for f in prim}
    co, spR = cert_int_poly([R[f] for f in prim], f"R_{n}")
    if n in SU.R_POLYS:
        assert co == SU.R_POLYS[n], (n, co)
    assert abs(co[-1]) == 1
    coords, orders = PK.group_data(prim, D)
    chis = PK.characters(coords, orders)
    rn = reduce_form(r0, 0, s0)
    if rn not in coords:
        rn = (rn[0], -rn[1], rn[2])
    assert rn in coords
    byco = {v: k for k, v in coords.items()}

    def cmul(f1, f2):
        return byco[tuple((a + b) % o for a, b, o in zip(coords[f1], coords[f2], orders))]

    def cinv(f):
        return byco[tuple((-a) % o for a, o in zip(coords[f], orders))]

    fl = {f: log(fabs(R[f])) for f in prim}
    for f in prim:
        assert spare_of(R[cmul(f, rn)] * R[f] - 1) >= need, (n, f)
    Lp, reps, M = PK.epstein_Lprime0(prim, D, chis)
    odd, S = [], {}
    worst = mp.dps
    for idx_, (ks, cord, chi, isreal) in enumerate(chis):
        s_ = sum(chi[f] * fl[f] for f in prim)
        S[ks] = s_
        if cord == 1:
            continue
        if fabs(chi[rn] + 1) < 0.1:
            odd.append((ks, cord, chi, Lp[idx_]))
            worst = min(worst, spare_of((s_ + 24 * Lp[idx_]) / max(1, fabs(Lp[idx_]))))
        else:
            assert fabs(s_) < mpf(10) ** (-need), (n, ks, "even sum does not vanish")
    assert worst >= need, (n, worst)
    assert len(odd) == h // 2
    say(f"R-polynomial certified (spare {spR}); Pic = {' x '.join('Z/%d' % o for o in orders)}, twist class "
        f"{rn}; {len(odd)} odd characters; S_chi = -24 L'(0,chi) (spare >= {worst}); even sums vanish")
    # (H1) twisted group determinant over a transversal
    T, seen = [], set()
    for f in prim:
        if f not in seen:
            T.append(f)
            seen.add(f)
            seen.add(cmul(f, rn))
    N_ = matrix(len(T), len(T))
    for i, a in enumerate(T):
        ai = cinv(a)
        for j, f in enumerate(T):
            N_[i, j] = 2 * fl[cmul(ai, f)]
    dN = det(N_)
    prodS = mpc(1)
    prodL = mpf(1)
    for ks, cord, chi, l in odd:
        prodS *= S[ks]
        prodL *= fabs(l)
    sp1 = spare_of((dN - prodS) / max(1, fabs(prodS)))
    assert sp1 >= need, (n, nstr(dN, 20), nstr(prodS, 20))
    say(f"(H1) det(2 log|R_(a^-1 f)|)_(T x T) = prod_(chi odd) S_chi = 24^(h/2) prod|L'|: relative spare {sp1};  "
        f"R^-(<R_f>) = {nstr(fabs(dN), 25)}")
    # (H2) multipliers through exact projections to the smaller orders
    cond_divs = [d for d in range(1, fc + 1) if fc % d == 0]
    proj = {}
    for fp in cond_divs:
        if fp == fc:
            continue
        Dp = fp * fp * dK
        forms_p = [g for g in classes_of_disc(Dp) if prim_form(g)]
        pm = {}
        for f in prim:
            g = project_form(f, D, dK, fp)
            if g not in forms_p:
                g = (g[0], -g[1], g[2])
            assert g in forms_p, (n, fp, f, g)
            pm[f] = g
        assert set(pm.values()) == set(forms_p), (n, fp)
        proj[fp] = (Dp, forms_p, pm)
    Cvals = {}
    prodC = mpf(1)
    for ks, cord, chi, l in odd:
        levels = [fc]
        for fp in cond_divs:
            if fp == fc:
                continue
            Dp, forms_p, pm = proj[fp]
            ok = True
            val = {}
            for f in prim:
                g = pm[f]
                if g in val:
                    if fabs(val[g] - chi[f]) > mpf(10) ** (-need):
                        ok = False
                        break
                else:
                    val[g] = chi[f]
            if ok:
                levels.append(fp)
        for m1 in levels:
            for m2 in levels:
                assert gcd(m1, m2) in levels, (n, ks, levels)
        fmin = min(levels)
        if fmin == fc:
            C = mpf(1)
            desc = "primitive"
        else:
            Dp, forms_p, pm = proj[fmin]
            assert len(forms_p) > 1, (n, ks, fmin)
            val = {pm[f]: chi[f] for f in prim}
            Lp_p, _, _ = PK.epstein_Lprime0(forms_p, Dp, [(ks, cord, val, cord <= 2)])
            C = (l / Lp_p[0]).real
            desc = f"from conductor {fmin} (disc {Dp})"
        Cvals[ks] = (cord, fmin, C, desc)
        prodC *= C
    prodC_int, spC = cert_int(prodC, "prod C (hyperbolic)")
    say(f"(H2) odd characters: " + "; ".join(f"chi{ks} (order {cord}) {desc}, C = {nstr(C, 8)}"
                                            for ks, (cord, fmin, C, desc) in sorted(Cvals.items())))
    say(f"      prod_(chi odd) C_chi(0) = {prodC_int}   (spare {spC})")
    # (H3) PARI
    Rpol = "Pol(" + str(co) + ", x)"
    certdeg = 16
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
    d = parse_tagged(run_gp(script))
    assert int(d["NROOTS"]) == h and int(d["NTAU"]) == 1, (n, d.get("NROOTS"), d.get("NTAU"))
    assert int(d["HPDEG"]) == h, (n, "H^+ not generated", d["HPDEG"])
    assert int(d["KMRANK"]) == h // 2 and int(d["KPRANK"]) == h // 2 - 1, (n, d["KMRANK"], d["KPRANK"])
    idx = int(d["INDEX"])
    Qm = int(d["QMINUS"])
    eplus = int(d["EPLUS"])
    hH, wH, RH = int(d["H"]), int(d["W"]), mpf(d["REG"])
    hHp, wHp, RHp = int(d["HPH"]), int(d["HPW"]), mpf(d["HPREG"])
    cert = int(d["CERT"]) == 1 and int(d["HPCERT"]) == 1
    Rminus = fabs(dN) / idx
    spB = spare_of(((RH / RHp) - 2 ** (h // 2 - 1) * Rminus / Qm) / (RH / RHp))
    lhs = (hH * RH / wH) / (hHp * RHp / wHp)
    spA = spare_of((lhs * prodC_int - prodL) / prodL)
    pred = Fraction(24 ** (h // 2) * 2 ** (h // 2 - 1), Qm) * Fraction(hH * wHp, hHp * wH) * prodC_int
    say(f"(H3) PARI: H = Q[y]/({d['HPOL'][:50]}...): h_H = {hH}, Cl = {parse_int_vector(d['CYC'])}, w_H = {wH}, "
        f"R_H = {nstr(RH, 20)}; tau = sigma_r (unique automorphism of H/K' inverting all R_f); rank E^- = {h // 2}")
    say(f"      H^+ = H^tau: h = {hHp}, Cl = {parse_int_vector(d['HPCYC'])}, w = {wHp}, R = {nstr(RHp, 20)};  "
        f"Q^- = [E : E^+ E^-] = {Qm}  ([ker(1-tau) : E^+] = {eplus}, tau on mu: zeta -> zeta^{d['TAUMU']});  "
        f"{'bnfcertify = 1 for H and H^+' if cert else 'GRH-conditional'}")
    say(f"      R_H / R_H+ = 2^(h/2-1) R^-(E^-) / Q^-  with R^-(E^-) = det N / index: spare {spB}")
    say(f"      (hR/w)_H / (hR/w)_H+ = prod_(chi odd) L'_prim(0,chi) = prod L' / prod C: spare {spA}")
    say(f"      exact index [E^- : mu <R_f>] = {idx} = 24^{h // 2} * {Fraction(idx, 24 ** (h // 2))}"
        f"  =  24^(h/2) (2^(h/2-1)/Q^-) (h_H w_H+)/(h_H+ w_H) prod C = {pred}: {'OK' if pred == idx else 'MISMATCH'}")
    assert spB >= need, (n, nstr(RH / RHp, 20), nstr(2 ** (h // 2 - 1) * Rminus / Qm, 20))
    assert spA >= need, (n, nstr(lhs * prodC_int, 20), nstr(prodL, 20))
    assert pred == idx, (n, idx, pred)
    rec = {"n": n, "h": h, "orders": orders, "idx": idx, "Qm": Qm, "eplus": eplus, "hH": hH, "hHp": hHp,
           "wH": wH, "wHp": wHp, "RH": RH, "RHp": RHp, "cyc": parse_int_vector(d["CYC"]),
           "cycp": parse_int_vector(d["HPCYC"]), "prodC": prodC_int, "cert": cert, "pol": d["HPOL"],
           "polp": d["HPPOL"], "dK": dK, "fc": fc}
    # (H4) the sextic layer: an odd character of order 6 (if any)
    six = [(ks, cord, chi, l) for (ks, cord, chi, l) in odd if cord == 6]
    if layer and six:
        ks, cord, chi, l = six[0]
        A = [f for f in prim if fabs(chi[f] - 1) < 0.1]
        assert len(A) == h // 6
        cos = {}
        for f in prim:
            key = min(cmul(f, a) for a in A)
            cos.setdefault(key, []).append(f)
        assert len(cos) == 6
        U = []
        for key in sorted(cos):
            p_ = mpc(1)
            for f in cos[key]:
                p_ *= R[f]
            U.append(p_)
        cu, spU = cert_int_poly(U, "sextic coset units")
        assert abs(cu[-1]) == 1
        oddB = [(k2, c2, ch2, l2) for (k2, c2, ch2, l2) in odd if all(fabs(ch2[a] - 1) < 0.1 for a in A)]
        assert len(oddB) == 3
        prodLB, prodCB = mpf(1), mpf(1)
        for k2, c2, ch2, l2 in oddB:
            prodLB *= fabs(l2)
            prodCB *= Cvals[k2][2]
        prodCB_int, _ = cert_int(prodCB, "prod C (sextic)")
        script2 = GP_HEADER.format(prec0=prec_pari) + f"""
cu = Pol({cu}, y);
F = polredbest(polcompositum(y^2 - ({D}), cu)[1]);
print("FPOL ", F);
bF = bnfinit(F, 1);
print("FH ", bF.no); print("FCYC ", bF.cyc); print("FW ", bF.tu[1]); print("FREG ", bF.reg);
u = nfroots(bF.nf, subst(cu, y, x));
print("UNROOTS ", #u);
nu = #bF.fu;
E = matrix(#u, nu);
for(i=1, #u, e = bnfisunit(bF, u[i]); for(j=1, nu, E[i,j] = e[j]));
G = nfgaloisconj(bF);
sd = nfroots(bF.nf, x^2 - ({D}))[1];
tau = 0; ntau = 0;
for(k=1, #G, ok = 1; for(i=1, #u, if(nfgaloisapply(bF, G[k], u[i]) != 1/u[i], ok = 0)); if(ok && nfgaloisapply(bF, G[k], sd) == sd, tau = G[k]; ntau++));
print("NTAU ", ntau);
Tm = matrix(nu, nu);
for(j=1, nu, e = bnfisunit(bF, nfgaloisapply(bF, tau, bF.fu[j])); for(i=1, nu, Tm[i,j] = e[i]));
Km = matkerint(matid(nu) + Tm); Kp = matkerint(matid(nu) - Tm);
print("KMRANK ", #Km); print("KPRANK ", #Kp);
X = matrix(#Km, #u);
for(i=1, #u, x = matinverseimage(Km, E[i,]~); for(k=1, #Km, X[k,i] = x[k]));
sn = matsnf(X); idx = 1; for(k=1, #sn, if(sn[k] != 0, idx *= abs(sn[k]))); print("INDEX ", idx);
w = bF.tu[1]; zt = bF.tu[2]; et = bnfisunit(bF, nfgaloisapply(bF, tau, zt)); t = lift(et[#et]);
g0 = gcd(w, t - 1); gk = g0;
for(k=1, #Kp, v = factorback(bF.fu, Kp[,k]); e = bnfisunit(bF, nfgaloisapply(bF, tau, v) / v); gk = gcd(gk, lift(e[#e])));
ee = g0 / gk;
Qm = ee * abs(matdet(concat(Kp, Km))); print("QMINUS ", Qm); print("EPLUS ", ee); print("TAUMU ", t);
z = sd; c = 1;
while(poldegree(minpoly(z)) < 6 && c < 40, z = sd + sum(i=1, #u, c^i * (u[i] + 1/u[i])); c++);
print("FPDEG ", poldegree(minpoly(z)));
Fp = polredbest(minpoly(z));
print("FPPOL ", Fp);
bp = bnfinit(Fp, 1);
print("FPH ", bp.no); print("FPW ", bp.tu[1]); print("FPREG ", bp.reg);
print("FCERT ", bnfcertify(bF)); print("FPCERT ", bnfcertify(bp));
"""
        d2 = parse_tagged(run_gp(script2))
        assert int(d2["FPDEG"]) == 6, (n, "F^+ not generated", d2["FPDEG"])
        assert int(d2["UNROOTS"]) == 6 and int(d2["NTAU"]) == 1 and int(d2["KMRANK"]) == 3 and int(d2["KPRANK"]) == 2
        idxF, QmF = int(d2["INDEX"]), int(d2["QMINUS"])
        hF, wF, RF = int(d2["FH"]), int(d2["FW"]), mpf(d2["FREG"])
        hFp, wFp, RFp = int(d2["FPH"]), int(d2["FPW"]), mpf(d2["FPREG"])
        lhsF = (hF * RF / wF) / (hFp * RFp / wFp)
        spF = spare_of((lhsF * prodCB_int - prodLB) / prodLB)
        predF = Fraction(24 ** 3 * 2 ** 2, QmF) * Fraction(hF * wFp, hFp * wF) * prodCB_int
        certF = int(d2["FCERT"]) == 1 and int(d2["FPCERT"]) == 1
        say(f"(H4) sextic layer chi{ks} (odd, order 6): coset units u_b = N_(H/F)(R_f), minimal polynomial {cu}")
        say(f"      F = K'(u) = Q[y]/({d2['FPOL']}): h_F = {hF}, Cl = {parse_int_vector(d2['FCYC'])}, w = {wF}, R = {nstr(RF, 15)}; "
            f"F^+: h = {hFp}, w = {wFp}, R = {nstr(RFp, 15)}; Q^- = {QmF}   ({'bnfcertify = 1 for both' if certF else 'GRH'})")
        say(f"      exact index [E^-(F) : mu <u_b>] = {idxF}  =  24^3 (2^2/Q^-) (h_F w_F+)/(h_F+ w_F) prod C = {predF}: "
            f"{'OK' if predF == idxF else 'MISMATCH'};  class number formula check spare {spF}")
        assert predF == idxF, (n, idxF, predF)
        assert spF >= need, (n, spF)
        rec["sextic"] = (cu, hF, hFp, wF, wFp, QmF, idxF, certF, d2["FPOL"], d2["FPPOL"], prodCB_int)
    say(f"hyperbolic level {n} done in {time.time() - t0:.1f} s")
    return rec


if __name__ == "__main__":
    main(sys.argv)
