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

Usage:
    python3 scripts/robert_index_full.py --selftest            # n <= 15 (~3 min)
    python3 scripts/robert_index_full.py --selftest --with-23  # + degree 24 (~+2 min)
    python3 scripts/robert_index_full.py 9                     # one level
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
                                  is_primitive, vecs_of)

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
    else:
        for a in args:
            level(int(a))


if __name__ == "__main__":
    main(sys.argv)
