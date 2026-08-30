"""The Duke-Imamoglu-Toth comparison for the Schmidt phase (outlook.md 2.6).

The level-n circle data has a real-quadratic shadow: the Cartan geodesic has
trace -2n, length 2 log eps_n, eps_n = n + sqrt(n^2-1) -- discriminant
d = 4(n^2 - 1) = -4D.  This script computes, from the definitions, the
DIT-type invariants of that discriminant and tests whether the phase u_f is
their linear shadow.

Real-quadratic side (per SL2(Z)-class A of forms of disc d = 4(n^2-1)):
  * cycle integrals  C_A(f) = int_{Gamma_Q \\ S_Q} f(z) sqrt(d) dz / Q(z,1),
    computed as int_1^{eps^2} f(sigma(iy)) dy/y with sigma = [[w+, w-],[1,1]]
    (so C_A(1) = 2 log eps exactly), for f = j - 744 and for the Gamma-
    invariant f = log(y^6 |Delta(z)|) (the Eisenstein / Kronecker-limit
    kernel);  a second, independent parametrization (dtheta/sin theta along
    the semicircle) and random SL2(Z)-transforms of Q validate every value;
  * the Rademacher symbol Psi(gamma_A) of the automorph (the linking number
    of the modular knot with the trefoil, Ghys), computed twice: from the
    Dedekind-sum formula Phi = (a+d)/c - 12 sign(c) s(d,|c|),
    Psi = Phi - 3 sign(c(a+d)), and from the minus-continued-fraction cycle
    Psi = sum (b_i - 3); the two must agree.

Phase side (imported from make_phase_atlas): u_f for the primitive classes
of disc 1 - n^2.

Comparisons (safe PSLQ, guard rails of CLAUDE.md):
  * log|u_1| and log|u_f| against [1, log eps, Tr_d(j - 744), sum E_A, pi^2]
    and against the per-class cycle-integral basis;
  * log|u_f| against [1, log eps] alone (S-unit shadow);
  * arg u_f against [pi, log eps], and the rationality of arg u_f / pi
    (the cusp-degeneration / Dedekind-sum probe of outlook.md 2.4);
  * every found relation must survive a +80-digit recomputation and a level
    not used in the fit; the log prints margins for every non-match.

Usage:
    python3 scripts/dit_comparison.py --selftest          # anchors
    python3 scripts/dit_comparison.py                     # full run, n<=17
    python3 scripts/dit_comparison.py --cusp              # arg u / pi probe
"""
import sys
import argparse
from fractions import Fraction
from math import gcd, isqrt

sys.path.insert(0, 'scripts')

from mpmath import (mp, mpf, mpc, fabs, nstr, sqrt as msqrt, log as mlog,
                    pi as mppi, quad, floor as mfloor, ceil as mceil, arg as marg,
                    pslq, exp as mexp)

import moduli_invariants as MI
from make_phase_atlas import (LevelData, sl2_reduce, J_at, LEVELS)


# ==========================================================================
# indefinite binary quadratic forms of discriminant d > 0
# ==========================================================================

def is_reduced_indef(f, d):
    a, b, c = f
    s = isqrt(d)
    # sqrt(d) - 2|a| < b < sqrt(d), 0 < b  (standard Gauss condition)
    if b <= 0 or b * b >= d:
        return False
    return (b + 2 * abs(a) > s or (b + 2 * abs(a)) ** 2 > d) and \
        (b - 2 * abs(a) < s and (2 * abs(a) - b) ** 2 < d or 2 * abs(a) <= b)


def reduced_indef_forms(d):
    """all Gauss-reduced indefinite forms of disc d (non-square d > 0):
    |sqrt(d) - 2|a|| < b < sqrt(d)."""
    out = []
    s = isqrt(d)
    assert s * s != d, "square discriminant not supported"
    for b in range(1, s + 1):
        if (d - b * b) % 4:
            continue
        m = (d - b * b) // 4          # = -ac > 0
        for a in range(1, m + 1):
            if m % a:
                continue
            c = -m // a
            for (aa, cc) in ((a, c), (c, a), (-a, -c), (-c, -a)):
                # b in (sqrt d - 2|aa|, sqrt d)
                t = 2 * abs(aa)
                lo_ok = (b + t > s) if (b + t) ** 2 > d else False
                if (b + t) ** 2 > d and b * b < d:
                    if (t - b) ** 2 < d or t <= b:
                        out.append((aa, b, cc))
    return sorted(set(out))


def rho_step(f, d):
    """the reduction-cycle neighbor of a reduced indefinite form."""
    a, b, c = f
    s = isqrt(d)
    # b' = -b mod 2|c|, with sqrt(d) - 2|c| < b' < sqrt(d)
    t = 2 * abs(c)
    b2 = (-b) % t
    # lift to the window (s - t, s]... choose b2 with (s - t) < b2 <= s
    while b2 * b2 < d and (b2 + t) * (b2 + t) <= d:
        b2 += t
    while b2 > 0 and b2 * b2 > d:
        b2 -= t
    # now b2 is the largest value < sqrt d in its class
    c2 = (b2 * b2 - d) // (4 * c)
    return (c, b2, c2)


def form_classes(d):
    """[(cycle of reduced forms)] -- one cycle per SL2(Z)-class."""
    forms = set(reduced_indef_forms(d))
    cycles = []
    seen = set()
    for f0 in sorted(forms):
        if f0 in seen:
            continue
        cyc = [f0]
        seen.add(f0)
        f = rho_step(f0, d)
        guard = 0
        while f != f0:
            assert f in forms, (d, f0, f)
            cyc.append(f)
            seen.add(f)
            f = rho_step(f, d)
            guard += 1
            assert guard < 10000
        cycles.append(cyc)
    return cycles


def content(f):
    return gcd(gcd(abs(f[0]), abs(f[1])), abs(f[2]))


def pell_4(d):
    """minimal (t, u), t, u > 0, with t^2 - d u^2 = 4."""
    u = 1
    while True:
        t2 = d * u * u + 4
        t = isqrt(t2)
        if t * t == t2:
            return t, u
        u += 1


def automorph(f, d, tu=None):
    """generator of the (positive) automorph group of f in SL2(Z)."""
    a, b, c = f
    t, u = tu or pell_4(d)
    g = ((t - b * u) // 2, -c * u, a * u, (t + b * u) // 2)
    assert g[0] * g[3] - g[1] * g[2] == 1
    # check the form is preserved: Q(gx) = Q(x)
    A, B, C, D = g
    a2 = a * A * A + b * A * C + c * C * C
    b2 = 2 * a * A * B + b * (A * D + B * C) + 2 * c * C * D
    c2 = a * B * B + b * B * D + c * D * D
    assert (a2, b2, c2) == (a, b, c), (f, g)
    return g


# ==========================================================================
# cycle integrals
# ==========================================================================

def roots_of(f, d):
    """(w_plus, w_minus) = ((-b +- sqrt d)/(2a))."""
    a, b, _c = f
    sd = msqrt(mpf(d))
    return (-b + sd) / (2 * a), (-b - sd) / (2 * a)


def cycle_integral(f, d, fun, eps=None):
    """int_1^{eps^2} fun(sigma(iy)) dy/y,  sigma = [[w+, w-], [1, 1]].
    fun(1) --> exactly 2 log eps."""
    wp, wm = roots_of(f, d)
    eps = eps or (lambda tu: (tu[0] + pell_4(d)[1] * msqrt(mpf(d))) / 2)(pell_4(d))
    L = 2 * mlog(eps)

    def g(s):
        y = mexp(s)
        z = (wp * mpc(0, y) + wm) / (mpc(0, y) + 1)
        if z.imag < 0:
            z = z.conjugate()          # a < 0 puts sigma(iy) in -H; fold back
        return fun(z)
    # smooth periodic-ish integrand; subdivide for safety
    pieces = max(4, int(float(L)) + 1)
    pts = [L * k / pieces for k in range(pieces + 1)]
    return quad(g, pts)


def cycle_integral_theta(f, d, fun, eps=None):
    """independent parametrization: dtheta/sin(theta) along the semicircle
    |z - W| = R, endpoints matched to the y-parametrization window."""
    wp, wm = roots_of(f, d)
    W, R = (wp + wm) / 2, fabs(wp - wm) / 2
    if eps is None:
        t, u = pell_4(d)
        eps = (t + u * msqrt(mpf(d))) / 2

    def z_of_y(y):
        z = (wp * mpc(0, y) + wm) / (mpc(0, y) + 1)
        return mpc(z.real, fabs(z.imag))
    th1 = marg(z_of_y(mpf(1)) - W)
    th2 = marg(z_of_y(eps ** 2) - W)
    lo, hi = (th1, th2) if th1 < th2 else (th2, th1)

    def g(th):
        from mpmath import sin
        z = W + R * mexp(mpc(0, th))
        return fun(z) / sin(th)
    pieces = 6
    pts = [lo + (hi - lo) * k / pieces for k in range(pieces + 1)]
    return quad(g, pts)


def kernel_j744(z):
    return J_at(z) - 744


def kernel_logdisc(z):
    """log(Im(z)^6 |Delta(z)|): Gamma-invariant, real."""
    zr, _ = sl2_reduce(z)
    _E4, _E6, Dl = MI._E4E6D(zr)
    return 6 * mlog(zr.imag) + mlog(fabs(Dl))


def transform_form(f, g):
    """f composed with g = (p, q, r, s) in SL2(Z): Q(px + qy, rx + sy)."""
    a, b, c = f
    p, q, r, s = g
    return (a * p * p + b * p * r + c * r * r,
            2 * a * p * q + b * (p * s + q * r) + 2 * c * r * s,
            a * q * q + b * q * s + c * s * s)


# ==========================================================================
# Dedekind sums, Rademacher symbol, minus continued fractions
# ==========================================================================

def dedekind_sum(h, k):
    """s(h, k) as an exact Fraction (k > 0)."""
    assert k > 0
    h %= k
    # reciprocity-based fast evaluation
    if k == 1:
        return Fraction(0)
    if h == 0:
        return Fraction(0)
    # s(h,k) + s(k,h) = -1/4 + (h/k + k/h + 1/(hk))/12
    s = Fraction(0)
    sign = 1
    stack = []
    while h != 0:
        stack.append((h, k))
        h, k = k % h, h
    # rebuild via reciprocity backwards
    val = Fraction(0)
    for (hh, kk) in reversed(stack):
        # val currently = s(kk mod hh, hh) = s(kk, hh)
        val = -val - Fraction(1, 4) + \
            (Fraction(hh, kk) + Fraction(kk, hh) + Fraction(1, hh * kk)) / 12
    return val


def dedekind_sum_direct(h, k):
    """brute-force s(h,k) for validation."""
    def saw(x):
        num, den = x.numerator, x.denominator
        if num % den == 0:
            return Fraction(0)
        return x - Fraction(num // den) - Fraction(1, 2)
    return sum(saw(Fraction(j, k)) * saw(Fraction(h * j, k))
               for j in range(1, k))


def rademacher_phi(g):
    """Rademacher's Phi function on SL2(Z) (c != 0 branch)."""
    a, b, c, d = g
    if c == 0:
        assert a == d and abs(a) == 1
        return Fraction(b, d)
    sgn = 1 if c > 0 else -1
    return Fraction(a + d, c) - 12 * sgn * dedekind_sum(d, abs(c))


def rademacher_psi(g):
    """Psi = Phi - 3 sign(c(a+d)): the conjugation-invariant Rademacher
    symbol (= linking number with the trefoil for hyperbolic g)."""
    a, b, c, d = g
    val = rademacher_phi(g) - 3 * (1 if c * (a + d) > 0 else -1 if c * (a + d) < 0 else 0)
    assert val.denominator == 1
    return int(val)


def minus_cf_cycle(f, d):
    """the minus-continued-fraction period of the first root of f: the
    cycle (b_1, ..., b_k) of w = b - 1/w', all b_i >= 2 on the cycle.

    Pure iteration of the step  w -> 1/(ceil(w) - w)  on EXACT forms (the
    root is recomputed from the integer form each step, so nothing
    accumulates); the first repeated form closes the primitive cycle."""
    sd = msqrt(mpf(d))

    def root(a, b, c):
        return (-b + sd) / (2 * a)

    cur, hist, ms = f, {}, []
    for step in range(4000):
        if cur in hist:
            cyc = ms[hist[cur]:]
            assert all(m >= 2 for m in cyc), (f, d, cyc)
            return cyc
        hist[cur] = step
        w = root(*cur)
        m = int(mceil(w))
        assert fabs(w - m) > mpf(10) ** (-mp.dps // 2), (f, d, cur)
        ms.append(m)
        cur = transform_form(cur, (m, -1, 1, 0))
    raise RuntimeError(("minus CF did not cycle", f, d))


def psi_from_cf(f, d):
    """Psi(gamma_Q) = sum (b_i - 3) over the minus-CF cycle (Rademacher =
    Zagier-reduction data; the trefoil-linking form of the symbol)."""
    bs = minus_cf_cycle(f, d)
    return sum(bi - 3 for bi in bs), bs


# ==========================================================================
# safe PSLQ
# ==========================================================================

def safe_pslq(vec, names, dps, maxcoeff=10 ** 8, label=''):
    """PSLQ with the certification guard rails: a relation is reported only
    if (terms) x (coefficient digits) is far below the working precision and
    the residual is at the noise floor; otherwise the margin is logged.
    Returns (relation or None, message)."""
    tol = mpf(10) ** (-(dps - 20))
    rel = pslq([mpf(v) if not isinstance(v, mpf) else v for v in vec],
               maxcoeff=maxcoeff, maxsteps=200000)
    if rel is None:
        return None, (f"    {label}: NO integer relation "
                      f"(coeffs up to {maxcoeff}, {dps} digits)")
    resid = fabs(sum(r * v for r, v in zip(rel, vec)))
    scale = max(fabs(mpf(v)) for v in vec)
    info = sum(len(str(abs(r))) for r in rel if r)
    nterms = sum(1 for r in rel if r)
    spare = -mlog(resid / scale, 10) - info if resid > 0 else mp.inf
    if resid > tol * scale:
        return None, (f"    {label}: pslq candidate {rel} REJECTED "
                      f"(residual {nstr(resid / scale, 3)} above noise floor)")
    if info * nterms > dps / 3:
        return None, (f"    {label}: pslq candidate {rel} UNSAFE "
                      f"((terms x digits) = {info * nterms} vs dps = {dps})")
    msg = (f"    {label}: RELATION {rel} on {names} "
           f"(residual {nstr(resid / scale, 3)}, spare {nstr(spare, 4)})")
    return rel, msg


def pslq_target(target, tname, basis, bnames, dps, maxcoeff, label):
    """PSLQ of a target against a basis that may carry internal rational
    dependencies (the E_A are logs of algebraic invariants, so cross-stratum
    relations among THEM are genuine mathematics, not noise).  Internal
    relations (zero coefficient on the target) are logged, one participant is
    removed, and the test repeats -- so the verdict is about the target's
    membership in the rational span, not about the basis's redundancy."""
    lines = []
    basis, bnames = list(basis), list(bnames)
    for _ in range(len(basis) + 1):
        rel, msg = safe_pslq([target] + basis, [tname] + bnames, dps,
                             maxcoeff=maxcoeff, label=label)
        if rel is None:
            lines.append(msg)
            return None, lines
        if rel[0] != 0:
            lines.append(msg)
            return rel, lines
        lines.append(msg.replace('RELATION',
                                 'internal basis relation (target absent)'))
        idx = max(i for i in range(1, len(rel)) if rel[i] != 0)
        basis.pop(idx - 1)
        bnames.pop(idx - 1)
    return None, lines


# ==========================================================================
# the per-level comparison
# ==========================================================================

class DITLevel:
    def __init__(self, n, dps=150):
        self.n = n
        self.N = n * n - 1
        self.d = 4 * self.N
        # the j-kernel reaches e^{pi sqrt(d)} at the top of the tallest
        # geodesic, costing ~ pi sqrt(d)/ln 10 digits of the budget; work at
        # an elevated precision so the certified accuracy stays at `dps`
        import math
        self.dps_eff = dps + int(math.pi * math.sqrt(self.d) / math.log(10)) + 15
        mp.dps = self.dps_eff
        self.dps = dps
        t, u = pell_4(self.d)
        assert (t, u) == (2 * n, 1), (n, t, u)   # eps_d = eps_n on the nose
        self.eps = n + msqrt(mpf(self.N))
        self.cycles = form_classes(self.d)
        self.reps = [cyc[0] for cyc in self.cycles]
        self.Cj = {}
        self.E = {}
        self.psi = {}
        self.eps_of = {}      # PRIMITIVE stabilizer eigenvalue of the class:
        self.wrap = {}        # eps_n = eps_of^wrap (imprimitive classes with
        # a smaller core discriminant ride a shorter primitive geodesic that
        # the level geodesic wraps; the cycle integral runs over ONE
        # primitive period, as in DIT)
        for f in self.reps:
            g = content(f)
            core = (f[0] // g, f[1] // g, f[2] // g)
            dc = self.d // (g * g)
            tc, uc = pell_4(dc)
            epsQ = (tc + uc * msqrt(mpf(dc))) / 2
            k = int(round(float(mlog(self.eps) / mlog(epsQ))))
            assert fabs(epsQ ** k - self.eps) < mpf(10) ** (-dps + 30), (n, f)
            self.eps_of[f], self.wrap[f] = epsQ, k
            self.Cj[f] = cycle_integral(f, self.d, kernel_j744, eps=epsQ)
            self.E[f] = cycle_integral(f, self.d, kernel_logdisc, eps=epsQ)
            gam = automorph(core, dc, (tc, uc))
            p1 = rademacher_psi(gam)
            p2, _bs = psi_from_cf(f, self.d)
            assert p1 == p2, (n, f, p1, p2)
            self.psi[f] = p1
        mp.dps = dps

    def validate(self, f, trials=2, seed=11):
        """SL2-invariance + second parametrization, worst relative error."""
        import random
        rng = random.Random(seed)
        mp.dps = self.dps_eff
        worst = mpf(0)
        v0 = self.Cj[f]
        th = cycle_integral_theta(f, self.d, kernel_j744, eps=self.eps_of[f])
        worst = max(worst, fabs(th - v0) / max(fabs(v0), mpf(1)))
        for _ in range(trials):
            g = (1, 0, 0, 1)
            for _k in range(4):
                g2 = rng.choice([(1, 1, 0, 1), (1, -1, 0, 1), (0, -1, 1, 0)])
                g = (g[0] * g2[0] + g[1] * g2[2], g[0] * g2[1] + g[1] * g2[3],
                     g[2] * g2[0] + g[3] * g2[2], g[2] * g2[1] + g[3] * g2[3])
            f2 = transform_form(f, g)
            v2 = cycle_integral(f2, self.d, kernel_j744, eps=self.eps_of[f])
            worst = max(worst, fabs(v2 - v0) / max(fabs(v0), mpf(1)))
        mp.dps = self.dps
        return worst


def selftest():
    print("=" * 78)
    print("dit_comparison selftest")
    print("=" * 78)
    mp.dps = 60
    # Dedekind sums: reciprocity routine vs brute force
    for (h, k) in [(1, 1), (1, 3), (2, 5), (5, 7), (7, 12), (13, 30), (25, 43)]:
        assert dedekind_sum(h, k) == dedekind_sum_direct(h, k), (h, k)
    print("  Dedekind sums: reciprocity = brute force  OK")
    # Rademacher anchors: d = 5 golden class has Psi = 0, both routes
    g5 = automorph((1, 1, -1), 5)
    assert rademacher_psi(g5) == 0
    p, bs = psi_from_cf((1, 1, -1), 5)
    assert p == 0 and all(b == 3 for b in bs), (p, bs)
    print("  Rademacher anchor d=5 (golden geodesic): Psi = 0, minus-CF "
          f"cycle {bs}  OK")
    # nonzero anchor, hand-provable: (1,-4,1), d = 12, automorph [[4,-1],[1,0]]:
    # Phi = 4/1 - 12 s(0,1) = 4, Psi = 4 - 3 = 1; minus-CF cycle (4)
    g12 = automorph((1, -4, 1), 12)
    assert g12 == (4, -1, 1, 0), g12
    assert rademacher_psi(g12) == 1
    p12, bs12 = psi_from_cf((1, -4, 1), 12)
    assert p12 == 1 and bs12 == [4], (p12, bs12)
    print("  Rademacher anchor d=12: gamma = [[4,-1],[1,0]], Psi = 1 both "
          "routes  OK")
    # conjugation invariance of Psi
    import random
    rng = random.Random(3)
    g = automorph((1, 4, -4), 32)
    base = rademacher_psi(g)
    for _ in range(6):
        h = (1, 0, 0, 1)
        for _k in range(5):
            h2 = rng.choice([(1, 1, 0, 1), (1, -1, 0, 1), (0, -1, 1, 0)])
            h = (h[0] * h2[0] + h[1] * h2[2], h[0] * h2[1] + h[1] * h2[3],
                 h[2] * h2[0] + h[3] * h2[2], h[2] * h2[1] + h[3] * h2[3])
        det = h[0] * h[3] - h[1] * h[2]
        assert det == 1
        hi = (h[3], -h[1], -h[2], h[0])
        gc = (0, 0, 0, 0)
        # gc = h g h^-1
        m1 = (h[0] * g[0] + h[1] * g[2], h[0] * g[1] + h[1] * g[3],
              h[2] * g[0] + h[3] * g[2], h[2] * g[1] + h[3] * g[3])
        gc = (m1[0] * hi[0] + m1[1] * hi[2], m1[0] * hi[1] + m1[1] * hi[3],
              m1[2] * hi[0] + m1[3] * hi[2], m1[2] * hi[1] + m1[3] * hi[3])
        assert rademacher_psi(gc) == base
    print(f"  Psi conjugation-invariant (d=32 automorph, Psi = {base})  OK")
    # cycle integrals at n = 3 (d = 32): C(1) = 2 log eps, invariance,
    # second parametrization
    mp.dps = 90
    L3 = DITLevel(3, dps=90)
    for f in L3.reps:
        c1 = cycle_integral(f, 32, lambda z: mpf(1))
        assert fabs(c1 - 2 * mlog(L3.eps)) < mpf(10) ** (-70), f
        dev = L3.validate(f)
        print(f"  d=32 class {f}: C(1) = 2 log eps OK; invariance + "
              f"theta-parametrization dev {nstr(dev, 3)}")
        assert dev < mpf(10) ** (-60)
    print("ALL SELFTESTS PASSED")


def compare_level(n, dps=150, phase_dps=170):
    """the full comparison at one level; returns log lines."""
    out = []
    D = DITLevel(n, dps=dps)
    mp.dps = phase_dps
    L = LevelData(n, dps=phase_dps - 40, verify_dps=phase_dps,
                  want_circles=False)
    mp.dps = dps
    logeps = mlog(D.eps)
    out.append(f"n = {n}: d = {D.d}, {len(D.reps)} real-quadratic classes "
               f"(primitive {sum(1 for f in D.reps if content(f) == 1)}), "
               f"h(1-n^2) = {L.h} phase classes; 2 log eps = "
               f"{nstr(2 * logeps, 20)}")
    Tr = sum(D.Cj[f] for f in D.reps)
    Etot = sum(D.E[f] for f in D.reps)
    out.append(f"  Tr_d(j - 744) = {nstr(Tr, 25)}")
    out.append(f"  sum_A E_A     = {nstr(Etot, 25)}")
    out.append(f"  Psi(gamma_A)  = {[D.psi[f] for f in D.reps]}  "
               f"(minus-CF = Dedekind route, both)")
    wraps = {f: k for f, k in D.wrap.items() if k > 1}
    if wraps:
        out.append(f"  wrapped classes (level geodesic = primitive^k): "
                   f"{[(f, k) for f, k in wraps.items()]}")
    for f in D.reps[:2]:
        dev = D.validate(f, trials=1)
        out.append(f"  validation {f}: {nstr(dev, 3)}")
        assert dev < mpf(10) ** (-dps + 25)

    # ---- the PSLQ battery ------------------------------------------------
    # representatives: one log|u_f| per r-twin pair (the pair mate is -log)
    seen = set()
    reps = []
    for f in L.prim:
        key = tuple(sorted([f, L.twin[f]]))
        if key not in seen:
            seen.add(key)
            reps.append(f)
    principal = min(L.prim, key=lambda f: (f[0], abs(f[1])))
    assert principal == (1, 0, (n * n - 1) // 4)

    basis_names = ['1', 'log_eps', 'Tr_d(j-744)', 'sum_E', 'pi^2']
    basis = [mpf(1), logeps, Tr.real, Etot, mppi ** 2]
    # one representative per conjugate pair of d-classes (A and A^-1 share
    # Re Cj and E exactly -- feeding both hands PSLQ a degenerate relation)
    pair_reps = []
    for f in D.reps:
        dup = any(fabs(D.E[f] - D.E[g]) < mpf(10) ** (-dps + 30) and
                  fabs(D.Cj[f].real - D.Cj[g].real) < mpf(10) ** (-dps + 30)
                  * max(1, fabs(D.Cj[g].real)) for g in pair_reps)
        if not dup:
            pair_reps.append(f)
    out.append(f"  conjugate-pair representatives: {pair_reps}")
    percls = [D.Cj[f].real for f in pair_reps] + [D.E[f] for f in pair_reps]
    percls_names = [f'Cj{f}' for f in pair_reps] + [f'E{f}' for f in pair_reps]
    imcls = [D.Cj[f].imag for f in pair_reps
             if fabs(D.Cj[f].imag) > mpf(10) ** (-dps + 30)]
    imcls_names = [f'ImCj{f}' for f in pair_reps
                   if fabs(D.Cj[f].imag) > mpf(10) ** (-dps + 30)]

    for f in reps:
        u = L.u[f]
        lu = mlog(fabs(u))
        if fabs(lu) < mpf(10) ** (-40):
            out.append(f"  class {f}: |u| = 1 (skip log tests)")
            continue
        tag = 'principal' if f == principal else str(f)
        _r1, m1 = pslq_target(lu, 'log|u|', basis, basis_names, dps, 10 ** 8,
                              f'{tag}: log|u| vs aggregate basis')
        out.extend(m1)
        _r2, m2 = pslq_target(lu, 'log|u|', [mpf(1), logeps], ['1', 'log_eps'],
                              dps, 10 ** 10, f'{tag}: log|u| vs [1, log eps]')
        out.extend(m2)
        if len(percls) + 1 <= 13:
            _r3, m3 = pslq_target(lu, 'log|u|', percls, percls_names, dps,
                                  10 ** 6,
                                  f'{tag}: log|u| vs per-class integrals')
            out.extend(m3)
        au = marg(u)
        if fabs(au) > mpf(10) ** (-40) and fabs(fabs(au) - mppi) > mpf(10) ** (-40):
            _r4, m4 = pslq_target(au, 'arg u', [mppi, logeps],
                                  ['pi', 'log_eps'], dps, 10 ** 10,
                                  f'{tag}: arg u vs [pi, log eps]')
            out.extend(m4)
            if imcls and len(imcls) + 2 <= 10:
                _r6, m6 = pslq_target(au, 'arg u', [mppi] + imcls,
                                      ['pi'] + imcls_names, dps, 10 ** 6,
                                      f'{tag}: arg u vs Im cycle integrals')
                out.extend(m6)
    # pair-sums against the aggregate basis (algebraic vs transcendental)
    for f in reps[:2]:
        S = L.u[f] + 1 / L.u[f]
        if fabs(S.imag) < mpf(10) ** (-40):
            _r5, m5 = pslq_target(S.real, 'S_x', basis, basis_names, dps,
                                  10 ** 6, f'S_x{f} vs aggregate basis')
            out.extend(m5)
    return out


def cusp_probe(levels=None, dps=110):
    """outlook.md 2.4 at the cheap end: is arg u_f / pi rational anywhere
    (Dedekind-sum-type quantization), and do the Ford levels n = 2c^2 + 1
    behave specially?  Prints the best small-denominator approximations."""
    mp.dps = dps
    levels = levels or LEVELS
    print("cusp-degeneration probe: rationality of arg u_f / pi "
          "(non-ambiguous classes)")
    ford = {2 * c * c + 1 for c in range(1, 6)}
    worst_close = []
    for n in levels:
        L = LevelData(n, dps=dps - 30, verify_dps=dps, want_circles=False)
        mp.dps = dps
        for f in L.prim:
            u = L.u[f]
            au = marg(u)
            if fabs(u.imag) < mpf(10) ** (-40):
                continue                      # ambiguous: exactly 0 or pi
            x = au / mppi
            # best rational approximation with denominator <= 10^4
            from fractions import Fraction as Fr
            best = None
            p0, q0, p1, q1 = 1, 0, 0, 1
            t = x
            for _ in range(200):
                cfl = int(mfloor(t))
                p0, p1 = cfl * p0 + p1, p0
                q0, q1 = cfl * q0 + q1, q0
                if q0 == 0 or q0 > 10 ** 4:
                    break
                err = fabs(x - mpf(p0) / q0)
                best = (p0, q0, err)
                fr = t - cfl
                if fr == 0:
                    break
                t = 1 / fr
            tagf = ' [Ford level]' if n in ford else ''
            if best and best[2] < mpf(10) ** (-30):
                print(f"  n={n} {f}: arg u/pi = {best[0]}/{best[1]} "
                      f"EXACTLY (err {nstr(best[2], 3)}){tagf}  <-- quantized!")
            else:
                worst_close.append((float(best[2] * best[1] ** 2), n, f, best))
    if worst_close:
        worst_close.sort()
        print("  no non-ambiguous class has arg u/pi rational "
              f"(denominators <= 10^4, {dps} digits).  closest approaches "
              "(err * q^2, the CF quality measure):")
        for q, n, f, best in worst_close[:5]:
            print(f"    n={n} {f}: arg u/pi ~ {best[0]}/{best[1]}, "
                  f"err {nstr(best[2], 3)}, quality {q:.3g}")
        print(f"  (a Dedekind-sum-type quantization would give exact small "
              f"fractions; {len(worst_close)} classes tested)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--cusp', action='store_true')
    ap.add_argument('--levels', type=int, nargs='*', default=None)
    ap.add_argument('--dps', type=int, default=150)
    args = ap.parse_args()
    if args.selftest:
        selftest()
        return
    if args.cusp:
        cusp_probe(args.levels)
        return
    levels = args.levels or [3, 5, 7, 9, 11, 13, 15, 17]
    print("=" * 78)
    print("DIT comparison: real-quadratic cycle-integral invariants of "
          "d = 4(n^2-1)\nagainst the Schmidt phase u_f of disc 1 - n^2")
    print("=" * 78)
    for n in levels:
        for line in compare_level(n, dps=args.dps):
            print(line)
        print()


if __name__ == '__main__':
    main()
