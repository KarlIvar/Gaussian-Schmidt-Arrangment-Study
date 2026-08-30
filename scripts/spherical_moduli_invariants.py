"""Spherical moduli invariants of Schmidt circles: numerical laboratory for
spherical-moduli-invariants.md.

Setting: the Riemann sphere of diameter 1 tangent to C at 0 (north pole
(0,0,1)); stereographic projection from the north pole.  A Schmidt circle
with Hermitian data M = (2q, -zeta; -conj zeta, 2m) (curvature 2q, center
zeta/2q, co-curvature 2m, det M = -1) projects to a circle of angular
radius theta on the sphere with

    cot theta  =  q + m  =:  c        (the spherical level).

Moduli problem: functions of X in SL_2(C) invariant under X -> u X gamma',
u in Gamma_sph = SU(2) cap SL_2(Z[i]) (the finite polar group, V_4 in PSL),
gamma' in SL_2(Z).  Invariants of X = (a b; c d):

    S  = Re(X^dagger X) = [A, 2Br, C],  A = |a|^2+|c|^2, C = |b|^2+|d|^2,
         B = conj(a) b + conj(c) d,  Br = Re B     -- the SHAPE form,
         disc S = -4(1 + c^2), c = -Im B = level;
    z_S = (-Br + i sqrt(1+c^2))/A  in H            -- the shape point;
    T  = X^T X  (complex symmetric, det 1, Z[i]-entries),
    Theta(X) = T(z_S, 1) * j'(z_S)                 -- the PHASE,
    u^2 := Theta^2 / ((2 pi)^2 eta(tau_c)^8),  tau_c = i sqrt(c^2+1).

Experiments:
  A. exact structure: level-c circles; count = r3(c^2+1)/3 = 4 H(4(c^2+1));
     free V4-action; descent X per circle; the shape bijection
     {level-c double cosets} <-> {ALL reduced forms of disc -4(c^2+1)};
     mirror/antipodal laws; the resultant identity Res(S,T) = -4c^2 (exact).
  B. invariance of Theta (sign character on Gamma_sph); the fiber
     (rotations about the circle's axis): |Theta| constant, arg linear at
     rate -1; anchors T(z_S,1) = -2c/eps_c on the polar lines,
     eps_c = c + sqrt(c^2+1) (the negative-Pell unit); |T(conj z_S,1)| =
     eps_c^2 |T(z_S,1)| per class; spherical radius by direct projection.
  C. shape polynomials: prod over level-c double cosets (x - j(z_S)) =
     prod_{f^2 | c^2+1} H_{-4(c^2+1)/f^2}(x), integer-certified; the trace
     slice t(4(c^2+1)) of Zagier's weight-3/2 form.
  D. the phase: laws (orientation twist u^2(-) = eps^4 u^2(+); mirror
     conj u^2[S] = u^2[S^-1]; the pair product), exact Z[sqrt d]-
     coefficients of P_c(y) = prod_{2H}(y - u^2) at the 2-torsion levels
     (d = the real genus field), the completed integer polynomial
     N_{F/Q} P_c in Z[y], the level norm m_c, factorizations.
  E. the iS-companion (PGL_2) family: per-class rational ratios.

Usage:
    python3 scripts/spherical_moduli_invariants.py            # A + B + C
    python3 scripts/spherical_moduli_invariants.py phase      # D + E
    python3 scripts/spherical_moduli_invariants.py all
Requires mpmath; D also uses sympy.
"""
import sys
import random
from math import isqrt, gcd
from fractions import Fraction

from mpmath import (mp, mpf, mpc, jtheta, exp, pi, sqrt as msqrt, fabs, arg,
                    nstr, nint, log10, floor, im, re, pslq, matrix, eig, norm,
                    gamma as mgamma)

CMAX_STRUCT = 20      # exact structure levels
CMAX_SHAPE = 8        # shape-polynomial levels
PHASE_LEVELS = (1, 2, 3, 4, 6)   # phase levels (2-torsion: 1,2,3,6; test: 4)
DPS_SHAPE = {1: 60, 2: 60, 3: 80, 4: 120, 5: 160, 6: 200, 7: 260, 8: 300}
DPS_PHASE = {1: 160, 2: 220, 3: 260, 4: 340, 6: 560}

# ======================= Gaussian integers (exact) =======================
def gadd(u, v): return (u[0]+v[0], u[1]+v[1])
def gsub(u, v): return (u[0]-v[0], u[1]-v[1])
def gmul(u, v): return (u[0]*v[0]-u[1]*v[1], u[0]*v[1]+u[1]*v[0])
def gconj(u): return (u[0], -u[1])
def gneg(u): return (-u[0], -u[1])
def gnorm(u): return u[0]*u[0]+u[1]*u[1]

def mmul(P, Q):
    (a,b),(c,d) = P
    (e,f),(g,h) = Q
    return ((gadd(gmul(a,e), gmul(b,g)), gadd(gmul(a,f), gmul(b,h))),
            (gadd(gmul(c,e), gmul(d,g)), gadd(gmul(c,f), gmul(d,h))))

def minv(P):
    (a,b),(c,d) = P
    return ((d, gneg(b)), (gneg(c), a))

def mdet(P):
    (a,b),(c,d) = P
    return gsub(gmul(a,d), gmul(b,c))

def mconj(P):
    (a,b),(c,d) = P
    return ((gconj(a), gconj(b)), (gconj(c), gconj(d)))

ONE = (((1,0),(0,0)), ((0,0),(1,0)))
IOTA = (((0,1),(0,0)), ((0,0),(0,-1)))     # diag(i,-i): z -> -z
S0 = (((0,0),(-1,0)), ((1,0),(0,0)))       # z -> -1/z
ANTI = S0                                   # antipodal = ANTI o conj (as maps)
GSPH = [ONE, IOTA, S0, mmul(IOTA, S0)]     # Gamma_sph mod +-1

def Tmat(lam):
    return (((1,0), lam), ((0,0), (1,0)))

# ============== circles as Hermitian data (A, B, C), B in Z[i] ==============
def herm_of_X(X):
    """M_X = (X^-1)^dagger M0 X^-1: (2 Im(c conj d), i(a conj d - b conj c);
    ..., 2 Im(a conj b)) -- circle-classification.md conventions."""
    (a,b),(c,d) = X
    A = 2*(gmul(c, gconj(d))[1])
    C = 2*(gmul(a, gconj(b))[1])
    B = gmul((0,1), gsub(gmul(a, gconj(d)), gmul(b, gconj(c))))
    return (A, B, C)

def herm_transform(M, g):
    """data of g(circle): (g^-1)^dagger M g^-1, exact integer arithmetic."""
    A, B, C = M
    gi = minv(g)
    (p,q_),(r,s) = gi
    Ag = (A,0); Cg = (C,0); Bc = gconj(B)
    m11 = gadd(gmul(gconj(p), Ag), gmul(gconj(r), Bc))
    m12 = gadd(gmul(gconj(p), B), gmul(gconj(r), Cg))
    m21 = gadd(gmul(gconj(q_), Ag), gmul(gconj(s), Bc))
    m22 = gadd(gmul(gconj(q_), B), gmul(gconj(s), Cg))
    n11 = gadd(gmul(m11, p), gmul(m12, r))
    n12 = gadd(gmul(m11, q_), gmul(m12, s))
    n22 = gadd(gmul(m21, q_), gmul(m22, s))
    assert n11[1] == 0 and n22[1] == 0
    return (n11[0], n12, n22[0])

# =================== enumeration: level-c circles, counts ===================
def circles_of_level(c, parity=1):
    """oriented circles with q+m = c >= 1, q,m >= 0: (q, m, zeta);
    parity=1: zeta == i (mod 2), the PSL2-family S; parity=0: zeta == 1,
    the PGL2-companion iS."""
    out = []
    for q in range(0, c+1):
        m = c - q
        N = 4*q*m + 1
        r = isqrt(N)
        for x in range(-r, r+1):
            y2 = N - x*x
            y = isqrt(y2)
            if y*y != y2 or x % 2 != (1 - parity):
                continue
            for yy in ({y, -y} if y else {0}):
                if yy % 2 == parity:
                    out.append((q, m, (x, yy)))
    return sorted(set(out))

def r3(n):
    cnt = 0
    r = isqrt(n)
    for x in range(-r, r+1):
        for y in range(-r, r+1):
            z2 = n - x*x - y*y
            if z2 < 0:
                continue
            z = isqrt(z2)
            if z*z == z2:
                cnt += 2 if z else 1
    return cnt

def reduced_forms(D):
    """ALL reduced positive forms (a,b,c) of discriminant D < 0 (including
    imprimitive) -- their number is the Hurwitz class number H(|D|) (no
    weight issues for |D| > 4)."""
    out = []
    b = D % 2
    while 3*b*b <= -D:
        q4 = b*b - D
        q = q4 // 4
        a = max(b, 1)
        while a*a <= q:
            if a and q % a == 0:
                cc = q // a
                if -a < b <= a <= cc:
                    out.append((a, b, cc))
                    if 0 < b < a < cc:
                        out.append((a, -b, cc))
            a += 1
        b += 2
    return sorted(out)

def reduce_form(a, b, c):
    assert a > 0 and b*b - 4*a*c < 0
    while True:
        if -a < b <= a <= c:
            if a == c and b < 0:
                b = -b
            return (a, b, c)
        if not (-a < b <= a):
            k = (a - b) // (2*a)
            b2 = b + 2*a*k
            c = a*k*k + b*k + c
            b = b2
        if a > c:
            a, b, c = c, -b, a
        elif a == c and b < 0:
            b = -b
        elif -a < b <= a <= c:
            return (a, b, c)

# V4 = Gamma_sph action on circle data
def v4_orbits(datas):
    acts = [lambda d: d,
            lambda d: (d[0], d[1], gneg(d[2])),
            lambda d: (d[1], d[0], gconj(d[2])),
            lambda d: (d[1], d[0], gneg(gconj(d[2])))]
    datas = set(datas)
    orbits, seen = [], set()
    for d in sorted(datas):
        if d in seen:
            continue
        orb = sorted({f(d) for f in acts})
        for e in orb:
            assert e in datas, (d, e)
            seen.add(e)
        orbits.append(orb)
    return orbits

def M_of_data(data):
    q, m, z = data
    return (2*q, gneg(z), 2*m)

# ============================= descent =============================
def descent_X(M, base='S'):
    """X in SL2(Z[i]) with herm_of_X-style data of X(base-line) equal to M.
    base 'S': the real line (terminal B = +-i); base 'iS': the imaginary
    axis (terminal B = +-1)."""
    term = (0,1) if base == 'S' else (1,0)
    G = ONE
    sign = 1
    cur = M
    for _ in range(10000):
        A, B, C = cur
        if A == 0 and C == 0:
            break
        if A == 0:
            assert B == term or B == gneg(term), cur
            g = Tmat((B[0]*C//2, 0)) if B[0] else Tmat((0, B[1]*C//2))
            cur = herm_transform(cur, g); G = mmul(g, G)
            continue
        if A < 0:
            cur = (-A, gneg(B), -C); sign = -sign
            continue
        lam = (round(B[0]/A), round(B[1]/A))
        if lam != (0,0):
            g = Tmat(lam)
            cur = herm_transform(cur, g); G = mmul(g, G)
            A, B, C = cur
        assert abs(C) < abs(A) or C == 0, cur
        cur = herm_transform(cur, S0); G = mmul(S0, G)
    A, B, C = cur
    assert A == 0 and C == 0 and (B == term or B == gneg(term)), cur
    if B == gneg(term):
        sign = -sign
    X = minv(G) if sign == 1 else mmul(minv(G), IOTA)
    assert mdet(X) == (1,0)
    if base == 'S':
        assert herm_of_X(X) == M, (herm_of_X(X), M)
    else:
        assert herm_transform((0, (1,0), 0), X) == M
    return X

def random_SL2Z(rng, k=6):
    T = (((1,0),(1,0)), ((0,0),(1,0)))
    Ti = (((1,0),(-1,0)), ((0,0),(1,0)))
    S = (((0,0),(-1,0)), ((1,0),(0,0)))
    g = ONE
    for _ in range(k):
        g = mmul(g, rng.choice([T, Ti, S]))
    return g

# ================== invariants: S, T, z_S, Theta, u^2 ==================
def shape_T_of(X, family='S'):
    """((A, 2Br, C), level, (T11, T12, T22)); family 'iS' applies the
    z -> iz conjugation R = diag(zeta8, 1/zeta8): S-form uses Im H12,
    T-form gets (i T11, T12, -i T22)."""
    (a,b),(c_,d) = X
    A = gnorm(a) + gnorm(c_)
    C = gnorm(b) + gnorm(d)
    Bf = gadd(gmul(gconj(a), b), gmul(gconj(c_), d))
    T11 = gadd(gmul(a,a), gmul(c_,c_))
    T12 = gadd(gmul(a,b), gmul(c_,d))
    T22 = gadd(gmul(b,b), gmul(d,d))
    if family == 'S':
        return (A, 2*Bf[0], C), -Bf[1], (T11, T12, T22)
    Tt11 = gmul((0,1), T11)
    Tt22 = gmul((0,-1), T22)
    return (A, 2*Bf[1], C), Bf[0], (Tt11, T12, Tt22)

def res_ST(Sf, Ts):
    """resultant of S(z,1) and T(z,1), exact in Z[i]:
    Res(az^2+bz+c, a'z^2+b'z+c') = (ac'-a'c)^2 - (ab'-a'b)(bc'-b'c)."""
    a, b, c = (Sf[0],0), (Sf[1],0), (Sf[2],0)
    ap, bp, cp = Ts[0], gmul((2,0), Ts[1]), Ts[2]
    t1 = gsub(gmul(a, cp), gmul(ap, c))
    t2 = gsub(gmul(a, bp), gmul(ap, b))
    t3 = gsub(gmul(b, cp), gmul(bp, c))
    return gsub(gmul(t1, t1), gmul(t2, t3))

# ==================== modular forms (theta constants) ====================
def E4E6D(tau):
    q = exp(mpc(0, 1) * pi * tau)
    t2, t3, t4 = jtheta(2, 0, q), jtheta(3, 0, q), jtheta(4, 0, q)
    a, b, c = t2 ** 4, t3 ** 4, t4 ** 4
    return ((a*a + b*b + c*c) / 2, (b + c) * (a + b) * (c - a) / 2,
            (t2 * t3 * t4) ** 8 / 256)

def sl2_reduce(tau):
    a_, b_, c, d = 1, 0, 0, 1
    t = tau
    for _ in range(100000):
        k = int(nint(t.real))
        if k:
            t -= k
            a_, b_ = a_ - k * c, b_ - k * d
        if abs(t) < 1 - mpf(10) ** (-30):
            t = -1 / t
            a_, b_, c, d = -c, -d, a_, b_
        else:
            break
    return t, (c, d)

def J_at(tau):
    tr, _ = sl2_reduce(tau)
    E4, _, Dl = E4E6D(tr)
    return E4 ** 3 / Dl

def Jp_at(tau):
    tr, (c, d) = sl2_reduce(tau)
    E4, E6, Dl = E4E6D(tr)
    return (-2 * pi * mpc(0, 1) * E4 * E4 * E6 / Dl) / (c * tau + d) ** 2

def Dq_at(tau):
    tr, (c, d) = sl2_reduce(tau)
    return E4E6D(tr)[2] / (c * tau + d) ** 12

def cval(z):
    return mpc(z[0], z[1])

def zS_of(Sf):
    A, twoBr, C = Sf
    Br = twoBr / 2
    return (-Br + mpc(0,1)*msqrt(mpf(A)*C - mpf(Br)**2)) / A

def Theta_of(X, family='S'):
    Sf, lev, Ts = shape_T_of(X, family)
    zS = zS_of(Sf)
    Tz = cval(Ts[0])*zS**2 + 2*cval(Ts[1])*zS + cval(Ts[2])
    return Tz * Jp_at(zS), lev, zS, Tz, Sf, Ts

def u2_of(X, norm2, family='S'):
    th, lev, zS, Tz, Sf, Ts = Theta_of(X, family)
    return th*th/norm2, lev, Sf

# ==================== certified recognition ====================
def cert_int(x, spare=20):
    """integer with absolute-error certification: |x - n| < 10^-spare and
    imaginary part below the same threshold."""
    if abs(im(x)) > mpf(10)**(-spare):
        return None
    n = nint(re(x))
    if abs(x - n) < mpf(10)**(-spare):
        return int(n)
    return None

def cert_quad(x, d, spare=20, maxc=None):
    """certified p + q sqrt(d), p, q integers, via PSLQ + residual check."""
    if abs(im(x)) > mpf(10)**(-spare):
        return None
    digs = int(floor(log10(max(mpf(1), abs(x))))) + 1
    r = pslq([re(x), mpf(1), msqrt(d)], maxcoeff=maxc or 10**(digs+18),
             maxsteps=3000000)
    if r is None or r[0] == 0:
        return None
    a, b, cc = r
    p, q = Fraction(-b, a), Fraction(-cc, a)
    err = abs(x - (mpf(p.numerator)/p.denominator
                   + (mpf(q.numerator)/q.denominator)*msqrt(d)))
    if err < mpf(10)**(-spare):
        return (p, q)
    return None

def esym(vals):
    esf = [mpc(1)]
    for u in vals:
        new = [esf[0]]
        for k in range(1, len(esf)):
            new.append(esf[k] + u*esf[k-1])
        new.append(u*esf[-1])
        esf = new
    return esf[1:]

# ============================ experiment A ============================
def experiment_A():
    print("=" * 72)
    print("A. Exact structure: counts, V4, shape bijection, laws")
    print("=" * 72)
    for c in range(1, CMAX_STRUCT + 1):
        datas = circles_of_level(c)
        n1 = len(datas)
        n2 = r3(c*c + 1)
        forms = reduced_forms(-4*(c*c+1))
        H = len(forms)
        orbs = v4_orbits(datas)
        free = all(len(o) == 4 for o in orbs)
        ok = (3*n1 == n2) and (n1 == 4*H) and free and (len(orbs) == H)
        assert ok, (c, n1, n2, H)
        print(f"  c={c:3d}: circles {n1:4d} = r3({c*c+1})/3 = 4 H({4*(c*c+1)})"
              f" = 4*{H};  V4 free, {H} orbits  OK")
    print("\n  shape bijection + resultant identity + laws:")
    for c in range(1, CMAX_STRUCT + 1):
        forms = set(reduced_forms(-4*(c*c+1)))
        got = []
        for orb in v4_orbits(circles_of_level(c)):
            reds = set()
            for d in orb:
                q, m, (x, yy) = d
                X = descent_X(M_of_data(d))
                Sf, lev, Ts = shape_T_of(X)
                assert lev == c
                disc = Sf[1]**2 - 4*Sf[0]*Sf[2]
                assert disc == -4*(c*c+1)
                # resultant identity, exact in Z[i]:
                #   Res(S, T) = 4(y^2 - (c^2+1)) = -4((q-m)^2 + x^2),
                # y = Im zeta; so T(z_S,1) T(conj z_S,1) = Res/A^2, and
                # Theta = 0 would force q = m, x = 0: only the c = 0 line.
                assert res_ST(Sf, Ts) == (4*(yy*yy - (c*c+1)), 0), (c, d)
                assert yy*yy - (c*c+1) == -((q-m)**2 + x*x)
                reds.add(reduce_form(*Sf))
            assert len(reds) == 1, (c, orb, reds)
            got.append(reds.pop())
        assert len(got) == len(set(got)) == len(forms) and set(got) == forms
        print(f"  c={c:3d}: {len(forms):3d} double cosets <-> all reduced forms"
              f" of disc {-4*(c*c+1)} (bijection, exact);"
              f"  Res(S,T) = 4(y^2 - c^2 - 1)  OK")
    # mirror & antipodal (exact, spot levels)
    for c in (5, 7):
        for d in circles_of_level(c)[:6]:
            q, m, z = d
            X = descent_X(M_of_data(d))
            red = reduce_form(*shape_T_of(X)[0])
            inv = reduce_form(red[0], -red[1], red[2])
            Xm = descent_X(M_of_data((q, m, gconj(z))))
            assert reduce_form(*shape_T_of(Xm)[0]) == inv
            Xa = mmul(ANTI, mconj(X))         # antipodal partner
            Sfa, leva, _ = shape_T_of(Xa)
            assert leva == -c and reduce_form(*Sfa) == red
        print(f"  c={c}: mirror = inverse class; antipodal = same class at"
              f" level -c  OK")

# ============================ experiment B ============================
def experiment_B():
    print("=" * 72)
    print("B. Theta: invariance, sign character, fiber, anchors, radius")
    print("=" * 72)
    mp.dps = 60
    rng = random.Random(20260830)
    # spherical radius by direct stereographic projection
    for c in (1, 3, 7):
        q, m, z = next(d for d in circles_of_level(c) if d[0] > 0)
        cen = cval(z)/(2*q); r = mpf(1)/(2*q)
        pts = [cen + r*exp(mpc(0,1)*mpf(t)) for t in (0.3, 1.7, 4.1)]
        def proj(w):
            den = 1 + abs(w)**2
            return matrix([w.real/den, w.imag/den, abs(w)**2/den])
        vs = [(proj(w) - matrix([0,0,mpf(1)/2]))/mpf('0.5') for w in pts]
        u_ = vs[1] - vs[0]; v_ = vs[2] - vs[0]
        nrm = matrix([u_[1]*v_[2]-u_[2]*v_[1], u_[2]*v_[0]-u_[0]*v_[2],
                      u_[0]*v_[1]-u_[1]*v_[0]])
        nrm = nrm/norm(nrm)
        costh = sum(nrm[i]*vs[0][i] for i in range(3))
        cot = fabs(costh/msqrt(1-costh**2))
        assert fabs(cot - (q+m)) < mpf(10)**-40
        print(f"  c={c}: stereographic cot(theta) = {nstr(cot, 15)} = q+m  OK")
    # invariance + sign character
    for d in rng.sample(circles_of_level(5), 3):
        X = descent_X(M_of_data(d))
        th0 = Theta_of(X)[0]
        worst = mpf(0)
        for _ in range(6):
            Y = mmul(rng.choice(GSPH), mmul(X, random_SL2Z(rng)))
            th1 = Theta_of(Y)[0]
            worst = max(worst, min(abs(th1-th0), abs(th1+th0)))
        for u, want in ((ONE, 1), (S0, 1), (IOTA, -1), (mmul(IOTA, S0), -1)):
            th1 = Theta_of(mmul(u, X))[0]
            got = 1 if abs(th1-th0) < abs(th1+th0) else -1
            assert got == want
        print(f"  d={d}: two-sided invariance to {nstr(worst, 3)};"
              f" sign character (+,+,-,-) on (1, S0, iota, iota S0)  OK")
    # fiber: rotations about the circle's axis
    for c in (1, 2, 5):
        d = circles_of_level(c)[1]
        X = descent_X(M_of_data(d))
        A, B, C = M_of_data(d)
        Mm = matrix([[mpf(A), cval(B)], [cval(gconj(B)), mpf(C)]])
        E, V = eig(Mm)
        v1 = matrix([V[0,0], V[1,0]]); v1 = v1/norm(v1)
        v2 = matrix([V[0,1], V[1,1]]); v2 = v2/norm(v2)
        U = matrix([[v1[0], v2[0]], [v1[1], v2[1]]])
        th0 = Theta_of(X)[0]
        Xc = [[cval(X[0][0]), cval(X[0][1])], [cval(X[1][0]), cval(X[1][1])]]
        for tt in (mpf('0.1'), mpf('0.35')):
            D = matrix([[exp(mpc(0,1)*tt/2), 0], [0, exp(-mpc(0,1)*tt/2)]])
            R = U * D * U.transpose_conj()
            Y = [[R[0,0]*Xc[0][0]+R[0,1]*Xc[1][0], R[0,0]*Xc[0][1]+R[0,1]*Xc[1][1]],
                 [R[1,0]*Xc[0][0]+R[1,1]*Xc[1][0], R[1,0]*Xc[0][1]+R[1,1]*Xc[1][1]]]
            a, b = Y[0]; c_, dd = Y[1]
            Af = abs(a)**2 + abs(c_)**2
            Bf = a.conjugate()*b + c_.conjugate()*dd
            zS = (-Bf.real + mpc(0,1)*msqrt(Af*(abs(b)**2+abs(dd)**2) - Bf.real**2))/Af
            Tz = (a*a+c_*c_)*zS**2 + 2*(a*b+c_*dd)*zS + (b*b+dd*dd)
            th = Tz * Jp_at(zS)
            assert fabs(abs(th)/abs(th0) - 1) < mpf(10)**-45
            assert fabs(arg(th/th0) + tt) < mpf(10)**-45
        print(f"  c={c}: fiber rotation by t: |Theta| constant,"
              f" arg Theta drops by exactly t (rate 1)  OK")
    # anchors: polar lines; the T-norm lemma: eps T(z_S,1) A/2 in K = Q(i sqrt N)
    # with K-norm (q-m)^2 + x^2 -- the Pell unit is exactly the archimedean
    # discrepancy between the two roots of S
    print("  anchors and the T-norm lemma:")
    for c in (1, 2, 3, 5, 6):
        N = c*c + 1
        eps = c + msqrt(N)
        X = descent_X(M_of_data((0, c, (0, 1))))
        th, lev, zS, Tz, Sf, Ts = Theta_of(X)
        assert Sf == (1, 0, N) and lev == c
        assert fabs(Tz - (-2*c/eps)) < mpf(10)**-45
        for d in circles_of_level(c):
            q, m, (x, yy) = d
            Y = descent_X(M_of_data(d))
            _, _, zS2, Tz2, Sf2, Ts2 = Theta_of(Y)
            v = eps * Tz2 * Sf2[0] / 2
            # v in K = Q(i sqrt N): rational real part, im/sqrt(N) rational
            for comp in (re(v), im(v)/msqrt(N)):
                if fabs(comp) < mpf(10)**-40:
                    continue
                pr = pslq([comp, mpf(1)], maxcoeff=10**8, maxsteps=200000)
                assert pr is not None, (c, d, nstr(v, 20))
            assert fabs(abs(v)**2 - ((q-m)**2 + x*x)) < mpf(10)**-40
            Tzbar = (cval(Ts2[0])*zS2.conjugate()**2
                     + 2*cval(Ts2[1])*zS2.conjugate() + cval(Ts2[2]))
            assert fabs(abs(Tzbar)/abs(Tz2) - eps**2) < mpf(10)**-40
        print(f"  c={c}: T(z_S,1) = -2c/eps on the lines; at every class"
              f" eps T(z_S,1) A/2 in K with norm (q-m)^2 + x^2,"
              f" and |T(conj z_S,1)|/|T(z_S,1)| = eps^2  OK")

# ============================ experiment C ============================
def hilbert_like_poly(D, dps):
    """prod over ALL reduced forms of disc D of (x - j(form)), certified
    integer coefficients (the product over f^2 | |D|/4 of the ring class
    polynomials H_{D/f^2})."""
    mp.dps = dps
    js = []
    for (a, b, cc) in reduced_forms(D):
        tau = (-b + mpc(0,1)*msqrt(-D)) / (2*a)
        js.append(J_at(tau))
    ek = esym(js)
    coeffs = []
    for k, v in enumerate(ek, 1):
        n = cert_int(v, spare=max(20, dps//5))
        assert n is not None, (D, k, nstr(v, 30))
        coeffs.append(n)
    return coeffs, js

def experiment_C():
    print("=" * 72)
    print("C. Shape polynomials and the Zagier trace slice t(4(c^2+1))")
    print("=" * 72)
    for c in range(1, CMAX_SHAPE + 1):
        N = c*c + 1
        dps = DPS_SHAPE[c]
        mp.dps = dps
        # j-values over the level's double cosets
        js = []
        for orb in v4_orbits(circles_of_level(c)):
            X = descent_X(M_of_data(orb[0]))
            Sf, lev, Ts = shape_T_of(X)
            js.append(J_at(zS_of(Sf)))
        ekj = esym(js)
        coeffs, js2 = hilbert_like_poly(-4*N, dps)
        okpoly = True
        for k, v in enumerate(ekj, 1):
            n = cert_int(v, spare=max(20, dps//5))
            if n is None or n != coeffs[k-1]:
                okpoly = False
        # trace slice
        tr = sum(js) - 744*len(js)
        trn = cert_int(tr, spare=max(20, dps//5))
        # strata
        fs = [f for f in range(1, isqrt(N)+1) if N % (f*f) == 0]
        strata = " * ".join(f"H_{-4*N//(f*f)}" for f in fs)
        print(f"  c={c}: prod(x - j(z_S)) = {strata}(x)"
              f" (degree {len(js)}) certified integer: {okpoly};"
              f"  trace t({4*N}) = {trn}")
        assert okpoly and trn is not None
    print("  (c=1 cross-check: t(8) = 7256, the hyperbolic n=3 slice value)")

# ============================ experiment D ============================
GENUS_D = {1: None, 2: 5, 3: 5, 6: 37}   # real genus field Q(sqrt d) or Q

def phase_values(c, dps, family='S'):
    mp.dps = dps
    N = c*c + 1
    tau_c = mpc(0,1)*msqrt(N)
    norm2 = (2*pi)**2 * Dq_at(tau_c) ** (mpf(1)/3)
    vals = []
    base = 'S' if family == 'S' else 'iS'
    for orb in v4_orbits(circles_of_level(c, parity=1 if family=='S' else 0)):
        X = descent_X(M_of_data(orb[0]), base=base)
        for Y in (X, mmul(X, IOTA)):
            u2, lev, Sf = u2_of(Y, norm2, family)
            assert abs(lev) == c
            red = reduce_form(*Sf)
            vals.append((red, '+' if lev > 0 else '-', u2))
    return vals

def experiment_D():
    import sympy
    print("=" * 72)
    print("D. The phase: laws and level polynomials")
    print("=" * 72)
    y = sympy.symbols('y')
    for c in PHASE_LEVELS:
        N = c*c + 1
        dps = DPS_PHASE[c]
        vals = phase_values(c, dps)
        eps = c + msqrt(N)
        byco = {(r, o): u for r, o, u in vals}
        H = len(vals)//2
        print(f"--- c={c} (N={N}, h={H}, dps {dps}) ---")
        # laws
        for (r, o), u in byco.items():
            if o == '+':
                rat = byco[(r, '-')]/u
                assert fabs(rat - eps**4) < mpf(10)**(-dps//2), (c, r)
            inv = reduce_form(r[0], -r[1], r[2])
            assert abs(byco[(inv, o)].conjugate() - u) < mpf(10)**(-dps//3)
        print(f"  orientation twist u2(-) = eps^4 u2(+) and mirror law: OK"
              f" (eps = {c} + sqrt({N}), the negative-Pell unit)")
        us = [u for _, _, u in vals]
        ek = esym(us)
        d = GENUS_D[c] if c in GENUS_D else None
        if c not in GENUS_D:
            # non-2-torsion class group: certified negative -- e_1 fits no
            # rational or quadratic value at this precision
            e1 = ek[0]
            digs = int(floor(log10(abs(e1)))) + 1
            hits = []
            for dd in (0, 2, 3, 5, 13, 17, 26, 34):
                if dd == 0:
                    r = pslq([re(e1), mpf(1)], maxcoeff=10**(digs+14),
                             maxsteps=1000000)
                else:
                    r = pslq([re(e1), mpf(1), msqrt(dd)],
                             maxcoeff=10**(digs+14), maxsteps=1000000)
                if r is not None and r[0] != 0:
                    hits.append(dd)
            print(f"  class group not 2-torsion: e_1 fits none of Q,"
                  f" Q(sqrt d) for d in (2,3,5,13,17,26,34): "
                  f"{'confirmed' if not hits else f'UNEXPECTED {hits}'}")
        if c in GENUS_D:
            if d is None:
                # rational level: P_c in Z[y] directly
                coeffs = [cert_int(v, spare=dps//5) for v in ek]
                assert all(n is not None for n in coeffs)
                P = y**len(us) + sum((-1)**k * n * y**(len(us)-k)
                                     for k, n in enumerate(coeffs, 1))
                print(f"  P_{c}(y) = {sympy.expand(P)}  in Z[y] (degree {len(us)})")
                Pp = sympy.Poly(P, y)
                print(f"    irreducible: {len(sympy.factor_list(Pp)[1]) == 1 and sympy.factor_list(Pp)[1][0][1] == 1}")
                print(f"    constant {sympy.factorint(coeffs[-1])}")
            else:
                fits = []
                okf = True
                for k, v in enumerate(ek, 1):
                    f = cert_quad(v, d, spare=dps//6)
                    fits.append(f)
                    if f is None:
                        okf = False
                print(f"  e_k of the 2H-multiset in Z[sqrt{d}]"
                      f" (real genus field): {'all certified' if okf else 'FAIL'}")
                assert okf
                sq = sympy.sqrt(d)
                P = y**len(us)
                Ps = y**len(us)
                for k, (p, q) in enumerate(fits, 1):
                    coef = sympy.Rational(p.numerator, p.denominator) \
                           + sympy.Rational(q.numerator, q.denominator)*sq
                    P += (-1)**k * coef * y**(len(us)-k)
                    coefs = sympy.Rational(p.numerator, p.denominator) \
                            - sympy.Rational(q.numerator, q.denominator)*sq
                    Ps += (-1)**k * coefs * y**(len(us)-k)
                    print(f"    e_{k} = {p} + {q} sqrt{d}")
                Pfull = sympy.expand(P*Ps)
                Pp = sympy.Poly(Pfull, y)
                assert all(co.is_integer for co in Pp.all_coeffs())
                print(f"  completed N_F/Q P_{c}(y) in Z[y], degree {Pp.degree()}:")
                print(f"    {Pp.as_expr()}")
                fl = sympy.factor_list(Pp)
                print(f"    factor degrees: {[(f.degree(), m) for f, m in fl[1]]}")
                cst = Pp.all_coeffs()[-1]
                print(f"    constant = {sympy.factorint(cst)}")
                # verify the certified polynomial kills the numeric values
                worst = mpf(0)
                for u in us:
                    tot = mpc(0)
                    for kk, co in enumerate(Pp.all_coeffs()):
                        tot = tot*u + int(co)
                    worst = max(worst, abs(tot)/abs(u)**Pp.degree())
                print(f"    residual |P(u2)|/|u2|^deg at all 2H phases:"
                      f" < {nstr(worst, 3)}")
        # the level norm m_c via PSLQ on power sums (works when the values
        # are quadratic over Q, i.e. the real genus field = coefficient field)
        p1 = sum(us); pm1 = sum(1/u for u in us)
        digs = int(floor(log10(abs(p1)))) + 1
        r = pslq([re(p1), re(pm1), mpf(1)], maxcoeff=10**(digs+16),
                 maxsteps=3000000)
        if r and r[0]:
            mv = Fraction(r[1], r[0])
            ok2 = True
            for rr in (2, 3):
                t = sum(u**rr for u in us) + (mpf(mv.numerator)/mv.denominator)**rr \
                    * sum(u**(-rr) for u in us)
                if cert_int(t, spare=dps//5) is None:
                    ok2 = False
            if ok2 and mv != 0:
                import sympy as sp
                print(f"  level norm m_{c} = u2 * sigma(u2) = {mv}"
                      f" = {sp.factorint(mv.numerator) if mv.denominator == 1 else mv}")
        print()

def experiment_E():
    print("=" * 72)
    print("E. The iS companion family: rational per-class ratios")
    print("=" * 72)
    for c in (1, 2, 3, 4):
        dps = DPS_PHASE.get(c, 300)
        vs = {(r, o): u for r, o, u in phase_values(c, dps, 'S')}
        vi = {(r, o): u for r, o, u in phase_values(c, dps, 'iS')}
        print(f"--- c={c}: u2_iS / u2_S per class ---")
        for (r, o), u in sorted(vi.items(), key=str):
            if (r, o) not in vs:
                continue
            rat = u / vs[(r, o)]
            fr = None
            if abs(im(rat)) < mpf(10)**(-dps//3):
                rr = pslq([re(rat), mpf(1)], maxcoeff=10**12, maxsteps=200000)
                if rr and rr[0]:
                    fr = Fraction(-rr[1], rr[0])
            print(f"  [{r}] {o}: ratio = {fr if fr is not None else nstr(rat, 18)}")

def experiment_symbolic():
    """Symbolic proof of the resultant identity
        Res(S, T) = 4(y^2 - (q+m)^2 - 1)
    for X in SL2(C): holomorphic and antiholomorphic entries are treated as
    independent variables, det X = 1 and its conjugate are imposed exactly,
    and the difference is simplified to zero."""
    import sympy as sp
    a, b, c, ab, bb, cb = sp.symbols('a b c ab bb cb')
    d = (1 + b*c)/a
    db = (1 + bb*cb)/ab
    q2 = c*db - cb*d                     # 2i Im(c conj d)
    m2 = a*bb - ab*b                     # 2i Im(a conj b)
    zeta = sp.I*(a*db - b*cb)
    zetab = -sp.I*(ab*d - bb*c)
    y2i = zeta - zetab                   # 2i y
    A = a*ab + c*cb
    C = b*bb + d*db
    Br2 = (ab*b + cb*d) + (a*bb + c*db)  # 2 Re B
    T11 = a*a + c*c; T12 = a*b + c*d; T22 = b*b + d*d
    f1 = A*T22 - T11*C
    f2 = A*2*T12 - T11*Br2
    f3 = Br2*T22 - 2*T12*C
    Res = f1*f1 - f2*f3
    y_sq = -(y2i**2)/4
    qm_sq = -((q2 + m2)**2)/4
    claim = sp.together(Res - 4*(y_sq - qm_sq - 1))
    zero = sp.simplify(sp.numer(claim)) == 0
    print("symbolic identity Res(S,T) = 4(y^2 - (q+m)^2 - 1) on SL2(C):",
          "PROVED" if zero else "FAILED")
    assert zero

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "struct"
    if which in ("struct", "all"):
        experiment_A()
        experiment_B()
        experiment_C()
    if which in ("phase", "all"):
        experiment_D()
        experiment_E()
    if which in ("symbolic", "all"):
        experiment_symbolic()
