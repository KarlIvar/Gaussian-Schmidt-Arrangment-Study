#!/usr/bin/env python3
"""spectral_verify.py — machine checks for spectral-geometry/eisenstein-catalan.md

Conventions.  H^3 = {(z,t) : z in C, t > 0}, metric (|dz|^2 + dt^2)/t^2,
volume dV = dx dy dt / t^3, Laplacian eigenfunction t^{1+s} <-> lambda = 1 - s^2.
Gamma = PSL_2(Z[i]);  Gamma_oo = stabiliser of oo  (translations by Z[i] and
the rotation z -> -z);  the Eisenstein series is

    E(P, s) = sum_{Gamma_oo \ Gamma} t(gamma P)^{1+s}
            = (1 / (4 zeta_K(1+s))) * sum_{(c,d) != (0,0)}
                  ( t / (|cz+d|^2 + |c|^2 t^2) )^{1+s} ,

where zeta_K(s) = zeta(s) L(s, chi_{-4}) is the Dedekind zeta function of
K = Q(i) (the second identity removes the coprimality condition and is what
the brute-force summation below uses).

Checks (all numerical statements verified here; the exact statements are
proved in the markdown):

 1. Humbert's volume:  vol(PSL_2(Z[i]) \ H^3) = G/3,  G = Catalan's constant
    — by direct integration of the Elstrodt–Grunewald–Mennicke fundamental
    domain  {|Re z| <= 1/2, 0 <= Im z <= 1/2, |z|^2 + t^2 >= 1}.
 2. zeta_K(s) = zeta(s) L(s, chi_{-4})  vs. the raw lattice sum
    (1/4) sum_{(m,n) != 0} (m^2+n^2)^{-s}.
 3. The scattering term:  the constant term of E(P, s) is
        t^{1+s} + phi(s) t^{1-s},   phi(s) = (pi/s) zeta_K(s) / zeta_K(s+1),
    by brute-force lattice summation of E and Fourier averaging.
    Residue check:  Res_{s=1} phi(s) = 3/(2G) = 1/(2 vol).
 4. Dirichlet series of the Euclidean count (euclidean-counting.md):
        sum_n N_e(n) n^{-s} = zeta(s-1) zeta(s) / zeta_K(s),
    partial sums of the sieved N_e against mpmath values.
 5. The kinematic lemma:  with  J_s(v) = int_v^oo (1+u^2)^{-(2+s)/2} du
    and  C_s = int_0^oo cosh(u)^{-s} du = (sqrt(pi)/2) G(s/2)/G((s+1)/2):
        int_0^oo J_s(x - 1/(4x)) dx = C_s ,
    and its Mellin refinement
        K_s(w) = int_0^oo J_s(x - 1/(4x)) x^{w-1} dx
               = 2^{s-w} G((1+s+w)/2) G((1+s-w)/2) / (w G(1+s)),  0 < w < 1+s,
    with K_s(1) = C_s.
 6. The horosphere limit.  For the plane-orbit kernel
        Theta_s(P) = sum_{planes of Gamma . P_R} cosh(dist(P, plane))^{-(2+s)},
    the constant term  CT(t) = int_{[0,1]^2} Theta_s((z,t)) dz  equals exactly
        CT(t) = sqrt(pi) G((1+s)/2)/G((2+s)/2) * t                 [lines]
              + pi t sum_{n>=1} (N_e(n)/n) J_s(nt - 1/(4nt))       [hemispheres]
    and CT(t) -> pi C_s / G  as t -> 0+  (equidistribution of expanding
    horospheres, with the limit = area(Y)-vs-vol(M) ratio; the numerics
    verify the limit against the exact series).
 7. The Mellin/period identity:
        int_0^oo Hem_s(t) t^{w-2} dt
            = pi K_s(w) zeta(w) zeta(w+1) / zeta_K(w+1),   1 < w < 1+s,
    by direct numerical t-integration of the exact series.
 8. Residue bookkeeping at w = 1:
        pi K_s(1) zeta(2)/zeta_K(2) = pi C_s / G
            = (3/(2G)) * (pi/3) C_s = Res_{w=1} E * <Theta_s>-unfolding,
    i.e. the Maass–Selberg residue and the counting residue agree.

Usage:  python3 scripts/spectral_verify.py         (~1-2 minutes)
        python3 scripts/spectral_verify.py --fast  (skips the slow lattice sum)
"""
import sys
import math
import numpy as np
import mpmath as mp

mp.mp.dps = 30

GCAT = mp.catalan
PASS, FAIL = "PASS", "FAIL"
results = []


def report(name, ok, detail=""):
    results.append(ok)
    print(f"[{PASS if ok else FAIL}] {name}" + (f"  {detail}" if detail else ""))


def rel_err(a, b):
    a, b = mp.re(mp.mpc(a)), mp.re(mp.mpc(b))
    return abs(a - b) / max(abs(b), mp.mpf('1e-30'))


# ---------------------------------------------------------------- L-functions
def Lchi4(s):
    """Dirichlet beta function L(s, chi_{-4}) via Hurwitz zeta."""
    return mp.mpf(4) ** (-s) * (mp.zeta(s, mp.mpf(1) / 4) - mp.zeta(s, mp.mpf(3) / 4))


def zetaK(s):
    """Dedekind zeta of Q(i)."""
    return mp.zeta(s) * Lchi4(s)


def phi_scatter(s):
    """Scattering term phi(s) = (pi/s) zeta_K(s)/zeta_K(s+1)."""
    return mp.pi / s * zetaK(s) / zetaK(s + 1)


# ---------------------------------------------------------------- 1. volume
def check_volume():
    inner = lambda x: mp.quad(lambda y: 1 / (2 * (1 - x * x - y * y)), [0, mp.mpf(1) / 2])
    vol = mp.quad(inner, [-mp.mpf(1) / 2, mp.mpf(1) / 2])
    target = GCAT / 3
    ok = rel_err(vol, target) < mp.mpf('1e-12')
    report("1. Humbert volume: integral over EGM fundamental domain = G/3",
           ok, f"integral = {mp.nstr(vol, 15)}, G/3 = {mp.nstr(target, 15)}")
    # Humbert's formula |d|^{3/2} zeta_K(2) / (4 pi^2) for d = -4:
    humbert = 8 * zetaK(2) / (4 * mp.pi ** 2)
    report("   Humbert formula 8*zeta_K(2)/(4 pi^2) = G/3",
           rel_err(humbert, target) < mp.mpf('1e-25'),
           f"= {mp.nstr(humbert, 15)}")
    return vol


# ---------------------------------------------------------------- 2. zeta_K
def check_zetaK():
    s = 3
    R = 220
    xs = np.arange(-R, R + 1)
    X, Y = np.meshgrid(xs, xs)
    q = X * X + Y * Y
    q = q[(q > 0) & (q <= R * R)].astype(float)
    lattice = 0.25 * np.sum(q ** (-s))
    # tail beyond radius R: ~ (1/4) * 2 pi int_R^oo r^{1-2s} dr
    tail = 0.25 * 2 * math.pi * R ** (2 - 2 * s) / (2 * s - 2)
    val = mp.mpf(lattice + tail)
    ok = rel_err(val, zetaK(3)) < mp.mpf('1e-6')
    report("2. zeta_K(3): lattice sum = zeta(3) L(3, chi_-4)", ok,
           f"lattice = {mp.nstr(val, 10)}, product = {mp.nstr(zetaK(3), 10)}")


# ---------------------------------------------------------------- 3. scattering
def E_bruteforce(zs, t, s, Rc=12, wrel=8.0):
    """E(P, s) at the points P = (z, t), z in the array zs, via the all-pairs
    lattice sum divided by 4 zeta_K(1+s); c- and d-tails are corrected by
    their integral approximations."""
    t = float(t)
    s = float(s)
    zs = np.asarray(zs, dtype=complex)
    total = np.zeros(len(zs))
    # c = 0 part, exactly: sum_{d != 0} (t/|d|^2)^{1+s} = 4 zeta_K(1+s) t^{1+s}
    c0 = float(4 * zetaK(1 + s)) * t ** (1 + s)
    total += c0
    cs = []
    for cre in range(-Rc, Rc + 1):
        for cim in range(-Rc, Rc + 1):
            if cre == 0 and cim == 0:
                continue
            if cre * cre + cim * cim <= Rc * Rc:
                cs.append(complex(cre, cim))
    zmid = np.mean(zs)
    for c in cs:
        ac2 = abs(c) ** 2
        W = wrel * abs(c) * max(t, 1.0)          # keep |cz + d| <= W
        centre = -c * zmid
        half = int(math.ceil(W + abs(c) * 1.5)) + 1
        dre = np.arange(math.floor(centre.real) - half, math.floor(centre.real) + half + 1)
        dim = np.arange(math.floor(centre.imag) - half, math.floor(centre.imag) + half + 1)
        DR, DI = np.meshgrid(dre, dim)
        D = (DR + 1j * DI).ravel()
        # broadcast over the z-points
        q = np.abs(c * zs[:, None] + D[None, :]) ** 2       # (nz, nd)
        mask = q <= W * W
        vals = np.where(mask, (t / (q + ac2 * t * t)) ** (1 + s), 0.0)
        total += vals.sum(axis=1)
        # d-tail (density-1 integral over |cz+d| > W):
        total += (math.pi / s) * t ** (1 + s) * (W * W + ac2 * t * t) ** (-s)
    # c-tail: for fixed z the full d-sum over one c is ~ (pi/s) t^{1-s} |c|^{-2s},
    # so the missing mass is (pi/s) t^{1-s} sum_{|c| > Rc} |c|^{-2s}, and the
    # lattice sum is evaluated exactly through 4 zeta_K(s):
    lattice_head = sum(abs(c) ** (-2 * s) for c in cs)
    total += (math.pi / s) * t ** (1 - s) * (float(4 * zetaK(s)) - lattice_head)
    return total / float(4 * zetaK(1 + s))


def check_scattering(fast=False):
    for (s, t) in [(2.0, 2.0)] if fast else [(2.0, 2.0), (1.5, 2.5)]:
        ngrid = 4
        xs = (np.arange(ngrid) + 0.5) / ngrid
        zs = np.array([complex(x, y) for x in xs for y in xs])
        E = E_bruteforce(zs, t, s)
        a0 = E.mean()
        phi_num = (a0 - t ** (1 + s)) / t ** (1 - s)
        phi_exact = phi_scatter(s)
        ok = rel_err(phi_num, phi_exact) < mp.mpf('3e-3')
        report(f"3. scattering phi({s}) from lattice sum of E at t = {t}", ok,
               f"numeric = {mp.nstr(mp.mpf(phi_num), 8)}, "
               f"(pi/s) zK(s)/zK(s+1) = {mp.nstr(phi_exact, 8)}, "
               f"rel err = {mp.nstr(rel_err(phi_num, phi_exact), 3)}")
    # residue of phi at s = 1: (s-1) phi(s) -> pi * Res zeta_K / zeta_K(2)
    res = mp.limit(lambda s: (s - 1) * phi_scatter(s), 1)
    target = 3 / (2 * GCAT)
    ok = rel_err(res, target) < mp.mpf('1e-15')
    report("   Res_{s=1} phi = 3/(2G) = 1/(2 vol)", ok,
           f"residue = {mp.nstr(res, 12)}, 3/(2G) = {mp.nstr(target, 12)}")
    reszK = mp.limit(lambda s: (s - 1) * zetaK(s), 1)
    report("   Res_{s=1} zeta_K = pi/4 (class number formula)",
           rel_err(reszK, mp.pi / 4) < mp.mpf('1e-15'),
           f"= {mp.nstr(reszK, 12)}")


# ---------------------------------------------------------------- 4. N_e series
def ne_sieve(X):
    """N_e(n) for n <= X:  N_e(2^a) = 2^a,  N_e(p^a) = p^{a-1}(p - chi_{-4}(p))."""
    spf = np.zeros(X + 1, dtype=np.int64)
    for p in range(2, X + 1):
        if spf[p] == 0:
            spf[p::p] = np.where(spf[p::p] == 0, p, spf[p::p])
    ne = np.zeros(X + 1, dtype=np.int64)
    ne[1] = 1
    for n in range(2, X + 1):
        p = spf[n]
        m, a = n, 0
        while m % p == 0:
            m //= p
            a += 1
        if p == 2:
            loc = 1 << a
        else:
            chi = 1 if p % 4 == 1 else -1
            loc = p ** (a - 1) * (p - chi)
        ne[n] = ne[m] * loc
    return ne


def check_dirichlet(ne):
    X = len(ne) - 1
    for s in (3, 4):
        n = np.arange(1, X + 1, dtype=float)
        partial = float(np.sum(ne[1:] / n ** s))
        # tail: N_e(n) ~ n/G on average => tail ~ (1/G) X^{2-s}/(s-2)
        tail = float(1 / GCAT) * X ** (2 - s) / (s - 2)
        val = mp.mpf(partial + tail)
        target = mp.zeta(s - 1) * mp.zeta(s) / zetaK(s)
        also = mp.zeta(s - 1) / Lchi4(s)
        ok = (rel_err(val, target) < mp.mpf('1e-5')
              and rel_err(target, also) < mp.mpf('1e-25'))
        report(f"4. sum N_e(n)/n^{s} = zeta(s-1) zeta(s)/zeta_K(s) "
               f"( = zeta(s-1)/L(s,chi) )", ok,
               f"sieved = {mp.nstr(val, 10)}, zeta ratio = {mp.nstr(target, 10)}")
    # mean value: sum_{n<=X} N_e(n) * 2G/X^2 -> 1
    S = float(ne.sum())
    ratio = mp.mpf(S) * 2 * GCAT / mp.mpf(X) ** 2
    report(f"   mean count: (2G/X^2) sum_(n<=X) N_e(n) at X = {X}",
           abs(ratio - 1) < mp.mpf('2e-3'), f"= {mp.nstr(ratio, 8)}")


# ---------------------------------------------------------------- 5. kinematic
def C_s(s):
    return mp.sqrt(mp.pi) / 2 * mp.gamma(mp.mpf(s) / 2) / mp.gamma((mp.mpf(s) + 1) / 2)


def J_num(s, v):
    """J_s(v) = int_v^oo (1+u^2)^{-(2+s)/2} du = int_{atan v}^{pi/2} cos^s."""
    return mp.re(mp.quad(lambda ph: mp.cos(ph) ** s, [mp.atan(v), mp.pi / 2]))


def K_closed(s, w):
    s, w = mp.mpf(s), mp.mpf(w)
    return (mp.mpf(2) ** (s - w) * mp.gamma((1 + s + w) / 2)
            * mp.gamma((1 + s - w) / 2) / (w * mp.gamma(1 + s)))


def check_kinematic():
    for s in (2, mp.mpf('3.7')):
        # substitute x = e^theta / 2:  int J_s(sinh th) e^th/2 dth
        val = mp.quad(lambda th: J_num(s, mp.sinh(th)) * mp.e ** th / 2,
                      [-mp.inf, mp.inf])
        ok = rel_err(val, C_s(s)) < mp.mpf('1e-10')
        report(f"5. kinematic lemma at s = {s}: int J_s(x - 1/4x) dx = C_s", ok,
               f"integral = {mp.nstr(val, 12)}, C_s = {mp.nstr(C_s(s), 12)}")
    for (s, w) in [(2, mp.mpf('1.3')), (mp.mpf('3.7'), mp.mpf('0.9')),
                   (mp.mpf('2.2'), mp.mpf('2.7'))]:
        val = mp.quad(lambda th: J_num(s, mp.sinh(th))
                      * mp.e ** (w * th) / mp.mpf(2) ** w, [-mp.inf, mp.inf])
        ok = rel_err(val, K_closed(s, w)) < mp.mpf('1e-9')
        report(f"   Mellin K_s(w) at (s,w) = ({s},{w}): Gamma closed form", ok,
               f"num = {mp.nstr(val, 10)}, closed = {mp.nstr(K_closed(s, w), 10)}")
    report("   K_s(1) = C_s (duplication)",
           all(rel_err(K_closed(s, 1), C_s(s)) < mp.mpf('1e-25')
               for s in (2, 3, mp.mpf('4.4'))))


# ---------------------------------------------------------------- 6. CT limit
def J2_np(v):
    """J_2(v) = (pi/2 - arctan v)/2 - v/(2(1+v^2)), vectorised;
    arctan2(1, v) = pi/2 - arctan(v) avoids cancellation, and the
    asymptotic branch 1/(3v^3) avoids float overflow for huge v."""
    v = np.asarray(v, dtype=float)
    small = v < 1e6
    vs = np.where(small, v, 1.0)
    exact = np.arctan2(1.0, vs) / 2 - vs / (2 * (1 + vs * vs))
    with np.errstate(over='ignore'):
        asym = 1.0 / (3.0 * np.where(small, 1.0, v) ** 3)
    return np.where(small, exact, asym)


def J1_np(v):
    """J_1(v) = 1 - v/sqrt(1+v^2), with asymptotic branch 1/(2v^2)."""
    v = np.asarray(v, dtype=float)
    small = v < 1e6
    vs = np.where(small, v, 1.0)
    exact = 1 - vs / np.sqrt(1 + vs * vs)
    with np.errstate(over='ignore'):
        asym = 1.0 / (2.0 * np.where(small, 1.0, v) ** 2)
    return np.where(small, exact, asym)


def hem(t, ne, s=2, vmax=250.0):
    """pi t sum_n (N_e(n)/n) J_s(nt - 1/(4nt)) with tail-integral correction,
    for s = 1, 2 (elementary J).  The tail uses the closed antiderivatives
        int J_1 = v - sqrt(1+v^2),   int J_2 = (1/2) v (pi/2 - arctan v)."""
    Jf = J2_np if s == 2 else J1_np
    nmax = min(len(ne) - 1, max(100, int(vmax / t)))
    n = np.arange(1, nmax + 1, dtype=float)
    v = n * t - 1 / (4 * n * t)
    series = math.pi * t * float(np.sum(ne[1:nmax + 1] / n * Jf(v)))
    # tail n > nmax with density N_e(n)/n ~ 1/G:
    # pi t (1/G) int_{nmax}^oo J_s(xt - 1/(4xt)) dx ~ (pi/G) int_a^oo J_s(v) dv
    a = nmax * t - 1 / (4 * nmax * t)
    if s == 2:
        # int_a^oo J_2 = (1/2)(1 - a arctan(1/a)) = 1/(6a^2) - 1/(10a^4) + ...
        tailint = (1 / (6 * a * a) - 1 / (10 * a ** 4) if a > 1e4
                   else 0.5 * (1 - a * math.atan2(1.0, a)))
    else:
        # int_a^oo J_1 = sqrt(1+a^2) - a = 1/(2a) - 1/(8a^3) + ...
        tailint = (1 / (2 * a) - 1 / (8 * a ** 3) if a > 1e4
                   else math.sqrt(1 + a * a) - a)
    return series + (math.pi / float(GCAT)) * tailint


def line_term(t, s):
    return float(mp.sqrt(mp.pi) * mp.gamma((1 + mp.mpf(s)) / 2)
                 / mp.gamma((2 + mp.mpf(s)) / 2)) * t


def check_ct_limit(ne):
    for s in (2, 1):
        target = float(mp.pi * C_s(s) / GCAT)
        vals = []
        for t in (0.2, 0.1, 0.05, 0.02):
            ct = line_term(t, s) + hem(t, ne, s=s)
            vals.append((t, ct))
        errs = [abs(ct - target) for _, ct in vals]
        ok = errs[-1] < 0.02 * target and errs[-1] < errs[0]
        detail = ", ".join(f"CT({t}) = {ct:.5f}" for t, ct in vals)
        report(f"6. horosphere limit (s = {s}): CT(t) -> pi C_s/G = {target:.6f}",
               ok, detail)


# ---------------------------------------------------------------- 7. Mellin
def check_mellin(ne):
    s = 2
    for w in (2.5, 1.5):
        target = mp.pi * K_closed(s, w) * mp.zeta(w) * mp.zeta(w + 1) / zetaK(w + 1)
        if w > 2:
            val = mp.quad(lambda t: hem(float(t), ne, s=s) * mp.mpf(t) ** (w - 2),
                          [0, 0.5, 2, 10, mp.inf])
        else:
            # split off [0, t0] with CT-limit head:  Hem ~ pi C_s/G  as t -> 0
            t0 = 0.05
            head = mp.quad(lambda t: (mp.mpf(hem(float(t), ne, s=s))
                                      - mp.pi * C_s(s) / GCAT) * mp.mpf(t) ** (w - 2),
                           [0.003, t0])
            head += mp.pi * C_s(s) / GCAT * mp.mpf(t0) ** (w - 1) / (w - 1)
            val = head + mp.quad(lambda t: hem(float(t), ne, s=s) * mp.mpf(t) ** (w - 2),
                                 [t0, 0.5, 2, 10, mp.inf])
        ok = rel_err(val, target) < mp.mpf('5e-3')
        report(f"7. Mellin identity at (s,w) = (2,{w}): "
               f"int Hem t^(w-2) dt = pi K_s(w) z(w)z(w+1)/zK(w+1)", ok,
               f"integral = {mp.nstr(val, 8)}, closed form = {mp.nstr(target, 8)}")


# ---------------------------------------------------------------- 8. residues
def check_residues():
    # <Theta_s, E(., w)> = (1/2) int_0^oo CT(t) t^{w-2} dt   (half square:
    # the cusp cross-section is C/Z[i] modulo the rotation z -> -z), so
    #   Res_{w=1} <Theta_s, E> = (1/2) pi K_s(1) zeta(2)/zeta_K(2)
    # must equal  Res_{w=1} E * int Theta_s = (3/(2G)) * (pi/3) C_s.
    for s in (2, mp.mpf('3.3')):
        lhs = mp.pi * K_closed(s, 1) * mp.zeta(2) / zetaK(2) / 2
        rhs = (3 / (2 * GCAT)) * (mp.pi / 3) * C_s(s)
        ok = rel_err(lhs, rhs) < mp.mpf('1e-25')
        report(f"8. residue bookkeeping at w=1 (s = {s}): "
               f"(1/2) pi K_s(1) zeta(2)/zeta_K(2) = Res E * int Theta_s", ok,
               f"both = {mp.nstr(lhs, 12)} = pi C_s/(2G)")
    # and the counting constant:  X^2/(2G) = (X^2/pi) * area(Y)/vol(M)
    lhs = 1 / (2 * GCAT)
    rhs = (1 / mp.pi) * (mp.pi / 6) / (GCAT / 3)
    report("   counting constant 1/(2G) = (1/pi) area(Y)/vol(M), "
           "area(Y) = pi/6, vol(M) = G/3",
           rel_err(lhs, rhs) < mp.mpf('1e-25'), f"both = {mp.nstr(lhs, 12)}")


# ---------------------------------------------------------------- main
def main():
    fast = "--fast" in sys.argv
    print("spectral_verify.py — checks for spectral-geometry/eisenstein-catalan.md")
    print("=" * 74)
    check_volume()
    check_zetaK()
    check_scattering(fast=fast)
    print("  (sieving N_e up to 10^6 ...)")
    ne = ne_sieve(10 ** 6)
    check_dirichlet(ne)
    check_kinematic()
    check_ct_limit(ne)
    check_mellin(ne)
    check_residues()
    print("=" * 74)
    n_ok = sum(results)
    print(f"{n_ok}/{len(results)} checks pass")
    return 0 if n_ok == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
