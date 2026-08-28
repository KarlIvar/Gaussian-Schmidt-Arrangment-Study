#!/usr/bin/env python3
"""alpha_circles.py — circles of the Gaussian Schmidt arrangement with a given
value of the hyperbolic invariant alpha, inside the ideal triangle (0, 1, oo).

For a circle omega contained in the upper half-plane, with Euclidean center
x + yi and radius r, set

    alpha(omega) = y / r = y * curvature = coth(hyperbolic radius),

which is invariant under the isometries PSL_2(R) of the hyperbolic plane.
For circles of the Schmidt arrangement, alpha is a positive integer n, and
(writing the curvature as 2q) the circles with alpha(omega) = n are exactly

    center (x + n i)/(2q), radius 1/(2q),  where  x^2 + n^2 - 1 = 4 q m

for integers q >= 1, m, x  (n odd: x even, circles of S = PSL_2(Z[i]).Rhat;
n even: x odd, circles of the companion family i.S).  The hyperbolic center
is (x + i*sqrt(n^2-1))/(2q) — the root of the positive definite binary
quadratic form (q, -x, m) of discriminant 1 - n^2.

The hyperbolic center lies in the closed ideal triangle T with vertices
0, 1, oo  iff  0 <= x <= 2q  and  x(2q - x) <= n^2 - 1;  it lies on an edge
of T iff one of these holds with equality (such circles get weight 1/2).
The weighted number of circles equals 3 * H(n^2 - 1), where H is the Hurwitz
class number.  See hyperbolic-counting.md for the derivations.

Usage:
    python3 alpha_circles.py 7                # window; press 'm' to switch model
    python3 alpha_circles.py 7 --model disk   # start in the Poincare disk model
    python3 alpha_circles.py 7 --save figures --no-show
    python3 alpha_circles.py 7 --list
    python3 alpha_circles.py --selftest
"""
import argparse
import math
import sys
from array import array
from fractions import Fraction

# --------------------------------------------------------------------------
# square roots modulo prime powers, and modulo 4q by CRT
# --------------------------------------------------------------------------

def tonelli(a, p):
    """A square root of a mod odd prime p, or None."""
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    # Tonelli-Shanks
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        t2, i = t, 0
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c = i, b * b % p
        t, r = t * c % p, r * b % p
    return r


def _sqrts_unit_pp(a, p, e):
    """All u in [0, p^e) with u^2 = a (mod p^e), for a coprime to p."""
    pe = p ** e
    a %= pe
    if p == 2:
        if e == 1:
            return [1]
        if e == 2:
            return [1, 3] if a % 4 == 1 else []
        if a % 8 != 1:
            return []
        r = 1
        for k in range(3, e):
            if (r * r - a) % (1 << (k + 1)):
                r += 1 << (k - 1)
        r %= pe
        half = pe >> 1
        return sorted({r, pe - r, (r + half) % pe, (pe - r + half) % pe})
    r = tonelli(a % p, p)
    if r is None:
        return []
    pk = p
    while pk < pe:
        pk *= p
        r = (r - (r * r - a) * pow(2 * r % pk, -1, pk)) % pk
    return sorted({r % pe, (pe - r) % pe})


def sqrts_mod_pp(a, p, e):
    """All x in [0, p^e) with x^2 = a (mod p^e)."""
    pe = p ** e
    a %= pe
    if a == 0:
        step = p ** ((e + 1) // 2)
        return list(range(0, pe, step))
    v = 0
    while a % p == 0:
        a //= p
        v += 1
    if v % 2:
        return []
    e1 = e - v
    us = _sqrts_unit_pp(a % (p ** e1), p, e1)
    h = p ** (v // 2)
    pe1 = p ** e1
    return sorted({(h * (u + t * pe1)) % pe for u in us for t in range(h)})


def sqrts_mod(a, factors):
    """All x modulo prod(p^e) with x^2 = a, factors = [(p, e), ...]."""
    roots, mod = [0], 1
    for p, e in factors:
        pe = p ** e
        rs = sqrts_mod_pp(a, p, e)
        if not rs:
            return []
        inv = pow(mod % pe, -1, pe)
        roots = [r0 + mod * (((r1 - r0) * inv) % pe) for r0 in roots for r1 in rs]
        mod *= pe
    return sorted(set(roots))


def spf_sieve(limit):
    """smallest prime factor for 0..limit"""
    spf = array('i', range(limit + 1))
    for i in range(2, int(limit ** 0.5) + 1):
        if spf[i] == i:
            for j in range(i * i, limit + 1, i):
                if spf[j] == j:
                    spf[j] = i
    return spf


def factor_with_spf(n, spf):
    f = {}
    while n > 1:
        p = spf[n]
        f[p] = f.get(p, 0) + 1
        n //= p
    return f

# --------------------------------------------------------------------------
# the circles with alpha = n whose hyperbolic center lies in the triangle T
# --------------------------------------------------------------------------

def alpha_circles(n):
    """All Schmidt circles omega with alpha(omega) = n and hyperbolic center in
    the closed ideal triangle T(0, 1, oo).  Returns a list of dicts with keys
    q, x, m, weight, center (complex, Euclidean), radius, hyp_center."""
    if n < 2:
        raise ValueError("need n >= 2 (n = 1 gives horocycles tangent to the boundary)")
    N = n * n - 1
    qmax = N // 4 + 1 if n % 2 else (N + 1) // 2
    spf = spf_sieve(4 * qmax)
    ynum = math.sqrt(N)
    out = []
    for q in range(1, qmax + 1):
        f = factor_with_spf(q, spf)
        f[2] = f.get(2, 0) + 2                      # modulus 4q
        for x in sqrts_mod(-N, sorted(f.items())):
            if x > 2 * q or x * (2 * q - x) > N:
                continue
            m = (x * x + N) // (4 * q)
            w = Fraction(1, 2) if (x == 0 or x == 2 * q or x * (2 * q - x) == N) else Fraction(1)
            out.append(dict(q=q, x=x, m=m, weight=w,
                            center=complex(x / (2 * q), n / (2 * q)),
                            radius=1 / (2 * q),
                            hyp_center=complex(x / (2 * q), ynum / (2 * q))))
    out.sort(key=lambda c: (c['q'], c['x']))
    return out

# --------------------------------------------------------------------------
# Hurwitz class numbers (independent check)
# --------------------------------------------------------------------------

def hurwitz(N):
    """Hurwitz class number H(N): weighted count of positive definite integer
    binary quadratic forms of discriminant -N (imprimitive included)."""
    if N == 0:
        return Fraction(-1, 12)
    if N < 0 or N % 4 in (1, 2):
        return Fraction(0)
    total = Fraction(0)
    b = N % 2
    while 3 * b * b <= N:
        M = (b * b + N) // 4
        a = max(b, 1)
        while a * a <= M:
            if M % a == 0:
                c = M // a
                mult = 1 if (b == 0 or b == a or a == c) else 2
                if a == b == c:
                    w = Fraction(1, 3)
                elif b == 0 and a == c:
                    w = Fraction(1, 2)
                else:
                    w = Fraction(1)
                total += mult * w
            a += 1
        b += 2
    return total

# --------------------------------------------------------------------------
# geometry helpers (Poincare disk)
# --------------------------------------------------------------------------

ZETA = complex(0.5, math.sqrt(3) / 2)   # centroid of T; phi(ZETA) = 0

def phi(z):
    """Half-plane -> unit disk, sending T to an ideal triangle with vertices
    at the cube roots of unity: phi(z) = (z - zeta)/(z - conj(zeta))."""
    return (z - ZETA) / (z - ZETA.conjugate())


def mobius_circle(a, b, c, d, z0, r):
    """Image of the circle |z - z0| = r under z -> (az+b)/(cz+d)."""
    den = abs(c * z0 + d) ** 2 - abs(c) ** 2 * r * r
    center = ((a * z0 + b) * (c * z0 + d).conjugate() - a * c.conjugate() * r * r) / den
    radius = abs(a * d - b * c) * r / abs(den)
    return center, radius


def phi_circle(z0, r):
    return mobius_circle(1, -ZETA, 1, -ZETA.conjugate(), z0, r)

# --------------------------------------------------------------------------
# plotting
# --------------------------------------------------------------------------

def _draw(ax, n, circles, model):
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors
    from matplotlib.patches import Circle as MplCircle

    ax.clear()
    ax.set_aspect('equal')
    qs = [c['q'] for c in circles]
    norm = mcolors.LogNorm(vmin=1, vmax=max(max(qs), 2))
    cmap = cm.viridis

    def style(c):
        fc = cmap(norm(c['q']))
        return dict(facecolor=(fc[0], fc[1], fc[2], 0.35), edgecolor=fc,
                    linewidth=1.1, linestyle='--' if c['weight'] != 1 else '-')

    tri_kw = dict(color='0.35', linewidth=1.4, zorder=1)

    if model == 'half':
        ymax = n / 2 + 0.75
        ax.plot([0, 0], [0, ymax], **tri_kw)
        ax.plot([1, 1], [0, ymax], **tri_kw)
        th = [math.pi * t / 400 for t in range(401)]
        ax.plot([0.5 + 0.5 * math.cos(t) for t in th],
                [0.5 * math.sin(t) for t in th], **tri_kw)
        ax.plot([-0.25, 1.25], [0, 0], color='0.6', linewidth=1.0, zorder=1)
        for c in circles:
            ax.add_patch(MplCircle((c['center'].real, c['center'].imag),
                                   c['radius'], zorder=2, **style(c)))
            ax.plot([c['hyp_center'].real], [c['hyp_center'].imag], '.',
                    color='black', markersize=3, zorder=3)
        ax.set_xlim(-0.15, 1.15)
        ax.set_ylim(-0.05, ymax)
        ax.set_title(f"upper half-plane:  ideal triangle (0, 1, ∞),  α = {n}")
    else:
        th = [2 * math.pi * t / 600 for t in range(601)]
        ax.plot([math.cos(t) for t in th], [math.sin(t) for t in th],
                color='0.6', linewidth=1.0, zorder=1)
        # triangle edges, mapped pointwise
        ts = [math.exp(u / 25) for u in range(-200, 201)]         # 1e-3.5 .. 1e3.5
        edges = ([complex(0, t) for t in ts],
                 [complex(1, t) for t in ts],
                 [0.5 + 0.5 * complex(math.cos(a), math.sin(a))
                  for a in [math.pi * k / 600 for k in range(1, 600)]])
        for e in edges:
            w = [phi(z) for z in e]
            ax.plot([z.real for z in w], [z.imag for z in w], **tri_kw)
        for c in circles:
            z0, r = phi_circle(c['center'], c['radius'])
            ax.add_patch(MplCircle((z0.real, z0.imag), r, zorder=2, **style(c)))
            hz = phi(c['hyp_center'])
            ax.plot([hz.real], [hz.imag], '.', color='black', markersize=3, zorder=3)
        ax.set_xlim(-1.05, 1.05)
        ax.set_ylim(-1.05, 1.05)
        ax.set_title(f"Poincaré disk:  ideal triangle at cube roots of unity,  α = {n}")
        ax.set_axis_off()

    W = sum(c['weight'] for c in circles)
    ax.figure.suptitle(
        f"Schmidt circles with α = {n}:  {len(circles)} circles, "
        f"weighted count {W} = 3·H({n * n - 1})   "
        f"(dashed = on an edge, weight ½;  dot = hyperbolic center;  color = curvature)",
        fontsize=10)


def show_or_save(n, circles, model, save_dir, show):
    import matplotlib
    if not show:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 8))
    state = {'model': model}
    _draw(ax, n, circles, state['model'])

    if save_dir:
        import os
        os.makedirs(save_dir, exist_ok=True)
        for mdl in ('half', 'disk'):
            _draw(ax, n, circles, mdl)
            path = os.path.join(save_dir, f"alpha{n}-{mdl}.png")
            fig.savefig(path, dpi=150, bbox_inches='tight')
            print(f"saved {path}")
        _draw(ax, n, circles, state['model'])

    if show:
        def toggle(_event=None):
            state['model'] = 'disk' if state['model'] == 'half' else 'half'
            _draw(ax, n, circles, state['model'])
            fig.canvas.draw_idle()

        def on_key(event):
            if event.key == 'm':
                toggle()

        fig.canvas.mpl_connect('key_press_event', on_key)
        try:
            from matplotlib.widgets import Button
            bax = fig.add_axes([0.82, 0.015, 0.16, 0.05])
            btn = Button(bax, 'switch model (m)')
            btn.on_clicked(toggle)
            fig._model_button = btn        # keep a reference alive
        except Exception:
            pass
        print("press 'm' (or use the button) to switch between the models")
        plt.show()

# --------------------------------------------------------------------------
# self-test
# --------------------------------------------------------------------------

def selftest():
    import random
    random.seed(2)

    # 1. modular square roots against brute force
    for _ in range(400):
        mod = random.randint(1, 600)
        a = random.randint(0, mod - 1)
        f = {}
        mm = mod
        d = 2
        while d * d <= mm:
            while mm % d == 0:
                f[d] = f.get(d, 0) + 1
                mm //= d
            d += 1
        if mm > 1:
            f[mm] = f.get(mm, 0) + 1
        got = sqrts_mod(a, sorted(f.items()))
        want = [x for x in range(mod) if (x * x - a) % mod == 0]
        assert got == want, (a, mod, got, want)
    print("modular square roots vs brute force: OK")

    # 2. Hurwitz class numbers against a table
    table = {3: Fraction(1, 3), 4: Fraction(1, 2), 7: 1, 8: 1, 11: 1,
             12: Fraction(4, 3), 15: 2, 16: Fraction(3, 2), 19: 1, 20: 2,
             23: 3, 24: 2, 27: Fraction(4, 3), 28: 2, 31: 3, 32: 3,
             35: 2, 36: Fraction(5, 2), 39: 4, 40: 2, 47: 5, 48: Fraction(10, 3)}
    for N, v in table.items():
        assert hurwitz(N) == Fraction(v), (N, hurwitz(N), v)
    print("Hurwitz class numbers vs table: OK")

    # 3. weighted circle count == 3 H(n^2 - 1)
    for n in range(2, 41):
        cs = alpha_circles(n)
        W = sum(c['weight'] for c in cs)
        assert W == 3 * hurwitz(n * n - 1), (n, W, 3 * hurwitz(n * n - 1))
        for c in cs:  # every circle satisfies the defining congruence & parity
            assert (c['x'] ** 2 + n * n - 1) == 4 * c['q'] * c['m']
            assert c['x'] % 2 != n % 2
    print("weighted count == 3*H(n^2-1) for 2 <= n <= 40: OK")

    # 4. the configuration is closed under the order-3 rotation
    #    rho: z -> 1/(1-z) of T and under the mirror z -> 1 - conj(z)
    for n in (2, 3, 4, 5, 7, 9, 12, 15):
        cs = alpha_circles(n)
        key = {(c['q'], c['x']) for c in cs}
        for c in cs:
            q, x, m = c['q'], c['x'], c['m']
            # Hermitian matrix (A, B, C) = (2q, -(x + n i), 2m); rho^{-1} = [[1,-1],[1,0]]
            # M' = (rho^{-1})^T M rho^{-1} for real matrices
            p11, p12, p21, p22 = 1, -1, 1, 0
            A, Bx, By, C = 2 * q, -x, -n, 2 * m
            A2 = p11 * (A * p11 + Bx * p21) + p21 * (Bx * p11 + C * p21)
            C2 = p12 * (A * p12 + Bx * p22) + p22 * (Bx * p12 + C * p22)
            Bx2 = p11 * (A * p12 + Bx * p22) + p21 * (Bx * p12 + C * p22)
            By2 = By * (p11 * p22 - p12 * p21)
            if A2 < 0:
                A2, Bx2, By2, C2 = -A2, -Bx2, -By2, -C2
            assert By2 == -n and A2 % 2 == 0 and C2 % 2 == 0
            assert (A2 // 2, -Bx2) in key, (n, q, x, (A2 // 2, -Bx2))
            assert (q, 2 * q - x) in key                      # mirror symmetry
    print("closure under the rotation of T and the mirror symmetry: OK")

    # 5. Moebius circle-image formula against sampled points
    for _ in range(50):
        a, b, c, d = (complex(random.uniform(-2, 2), random.uniform(-2, 2)) for _ in range(4))
        if abs(a * d - b * c) < 1e-3:
            continue
        z0 = complex(random.uniform(-2, 2), random.uniform(0.5, 3))
        r = random.uniform(0.05, 0.4)
        w0, rw = mobius_circle(a, b, c, d, z0, r)
        for k in range(12):
            z = z0 + r * complex(math.cos(k), math.sin(k))
            w = (a * z + b) / (c * z + d)
            assert abs(abs(w - w0) - rw) < 1e-8
    print("Moebius circle-image formula: OK")

    print("\nall self-tests passed")

# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('n', type=int, nargs='?', help='the value of alpha (integer >= 2)')
    ap.add_argument('--model', choices=['half', 'disk'], default='half',
                    help='initial model (default: half-plane)')
    ap.add_argument('--save', metavar='DIR', help='save PNGs of both models to DIR')
    ap.add_argument('--no-show', action='store_true', help='do not open a window')
    ap.add_argument('--list', action='store_true', help='print the circles as a table')
    ap.add_argument('--no-plot', action='store_true', help='skip plotting entirely')
    ap.add_argument('--selftest', action='store_true', help='run the self-tests and exit')
    args = ap.parse_args()

    if args.selftest:
        selftest()
        return
    if args.n is None:
        ap.error('give a value of alpha (or --selftest)')

    n = args.n
    circles = alpha_circles(n)
    W = sum(c['weight'] for c in circles)
    H = hurwitz(n * n - 1)
    print(f"alpha = {n}:  {len(circles)} circles in the ideal triangle, "
          f"weighted count {W}")
    print(f"3 * H({n * n - 1}) = {3 * H}   "
          f"({'match' if W == 3 * H else 'MISMATCH!'})")

    if args.list:
        print(f"\n{'q':>6} {'x':>7} {'m':>7} {'weight':>6}   center (Euclidean)      radius")
        for c in circles:
            print(f"{c['q']:>6} {c['x']:>7} {c['m']:>7} {str(c['weight']):>6}   "
                  f"({c['x']}/{2*c['q']}) + ({n}/{2*c['q']})i {'':6} 1/{2*c['q']}")

    if not args.no_plot:
        show_or_save(n, circles, args.model, args.save, show=not args.no_show)


if __name__ == '__main__':
    main()
