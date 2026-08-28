"""Numerical verification of the classification of circles in the Gaussian
Schmidt arrangement (orbit of Rhat under PSL_2(Z[i])).

We represent an oriented circle/line by an integral Hermitian matrix
    M = [[A, B], [conj(B), C]],  A, C in Z, B in Z[i],  det M = AC - |B|^2 = -1,
with incidence relation (conj(z), 1) M (z, 1)^T = 0, i.e.
    A|z|^2 + conj(B) z + B conj(z) + C = 0.
For A != 0 this is the circle with center -B/A and radius 1/|A|.

The extended real line Rhat corresponds to M0 = [[0, i], [-i, 0]].
If w = g z with g in SL_2(Z[i]), the image circle has matrix
    M' = (g^{-1})^dagger M (g^{-1}).
For orbit enumeration we may equivalently close under M -> h^dagger M h
for h ranging over generators and their inverses.

Claims verified here (with curvature normalized positive, A = 2n > 0,
zeta := 2n * center = -B = x + iy):
  (1) x is even, y is odd (i.e. zeta = i mod 2), and x^2 + y^2 = 1 mod 4n;
  (2) conversely EVERY residue class (x mod 2n, y mod 2n) satisfying these
      conditions occurs (so the classification is exact);
  (3) the lines of the arrangement are exactly Im z = k, k in Z
      (no vertical lines in the PSL_2 orbit);
  (4) N_e(n) := #{circles of curvature 2n with center in [0,1)^2}
             = n * prod_{p | n, p odd} (1 - chi_{-4}(p)/p);
  (5) sum_{n <= X} N_e(n) ~ X^2 / (2G),  G = Catalan's constant.
"""
from collections import deque, defaultdict

# ---------- Gaussian integer / matrix helpers (tuples (re, im)) ----------
def gm(p, q):
    a, b = p; c, d = q
    return (a * c - b * d, a * d + b * c)

def ga(p, q):
    return (p[0] + q[0], p[1] + q[1])

def gc(p):
    return (p[0], -p[1])

def matmul(X, Y):
    return tuple(
        tuple(ga(gm(X[i][0], Y[0][j]), gm(X[i][1], Y[1][j])) for j in range(2))
        for i in range(2)
    )

def dagger(X):
    return ((gc(X[0][0]), gc(X[1][0])), (gc(X[0][1]), gc(X[1][1])))

def det(X):
    return ga(gm(X[0][0], X[1][1]), (lambda t: (-t[0], -t[1]))(gm(X[0][1], X[1][0])))

# ---------- generators of SL_2(Z[i]): T, T^-1, T_i, T_i^-1, S, S^-1 ----------
GENS = []
for lam in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
    GENS.append((((1, 0), lam), ((0, 0), (1, 0))))          # translations
GENS.append((((0, 0), (-1, 0)), ((1, 0), (0, 0))))           # S : z -> -1/z
GENS.append((((0, 0), (1, 0)), ((-1, 0), (0, 0))))           # S^-1

M0 = (((0, 0), (0, 1)), ((0, -1), (0, 0)))                   # Rhat

def pack(M):
    (A, B), (Bc, C) = M
    assert A[1] == 0 and C[1] == 0, "diagonal must be real"
    return (A[0], B[0], B[1], C[0])

def unpack(t):
    A, bx, by, C = t
    return (((A, 0), (bx, by)), ((bx, -by), (C, 0)))

def orbit(LIM):
    """BFS closure of M0 under M -> h^dagger M h, keeping |A|,|C| <= LIM.
    Since det = -1 forces |B|^2 = 1 + AC, the state space is finite."""
    start = pack(M0)
    visited = {start}
    dq = deque([start])
    while dq:
        M = unpack(dq.popleft())
        for h in GENS:
            t2 = pack(matmul(dagger(h), matmul(M, h)))
            A, bx, by, C = t2
            if abs(A) <= LIM and abs(C) <= LIM and t2 not in visited:
                assert A * C - bx * bx - by * by == -1
                visited.add(t2)
                dq.append(t2)
    return visited

def collect(visited):
    """Group circles by half-curvature n; store residue class of zeta mod 2n.
    zeta := (2n) * center = -B for the representative with A = 2n > 0."""
    classes = defaultdict(set)
    lines = set()
    for (A, bx, by, C) in visited:
        if A == 0:
            lines.add((bx, by))
        else:
            if A < 0:
                A, bx, by, C = -A, -bx, -by, -C
            x, y = (-bx) % A, (-by) % A
            assert A % 2 == 0
            classes[A // 2].add((x, y))
    return classes, lines

def predicted(n):
    k = 2 * n
    return {(x, y) for x in range(0, k, 2) for y in range(1, k, 2)
            if (x * x + y * y) % (4 * n) == 1}

def factorize(n):
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f

def Ne_formula(n):
    res = n
    for p in factorize(n):
        if p > 2:
            chi = 1 if p % 4 == 1 else -1
            res = res // p * (p - chi)
    return res

# ---------- sanity check of the curvature/center formulas via a numeric fit ----------
def numeric_crosscheck(trials=200, seed=1):
    import random, cmath
    random.seed(seed)
    for _ in range(trials):
        g = (((1, 0), (0, 0)), ((0, 0), (1, 0)))
        for _ in range(random.randint(3, 8)):
            g = matmul(g, random.choice(GENS))
        (a, b), (c, d) = g
        ac, bc_, cc, dc = (complex(*a), complex(*b), complex(*c), complex(*d))
        Im_cdbar = (cc * dc.conjugate()).imag
        if Im_cdbar == 0:
            continue  # image is a line
        # image of three real points
        pts = []
        for t in (0.37, 1.91, -2.34):
            pts.append((ac * t + bc_) / (cc * t + dc))
        z1, z2, z3 = pts
        # circumcenter
        ax_, ay = z1.real, z1.imag
        bx_, by_ = z2.real, z2.imag
        cx, cy = z3.real, z3.imag
        D = 2 * (ax_ * (by_ - cy) + bx_ * (cy - ay) + cx * (ay - by_))
        ux = ((ax_**2 + ay**2) * (by_ - cy) + (bx_**2 + by_**2) * (cy - ay)
              + (cx**2 + cy**2) * (ay - by_)) / D
        uy = ((ax_**2 + ay**2) * (cx - bx_) + (bx_**2 + by_**2) * (ax_ - cx)
              + (cx**2 + cy**2) * (bx_ - ax_)) / D
        r = abs(z1 - complex(ux, uy))
        # predicted from the Hermitian-matrix formulas
        A_pred = 2 * Im_cdbar
        B_pred = 1j * (ac * dc.conjugate() - bc_ * cc.conjugate())
        ctr_pred = -B_pred / A_pred
        assert abs(abs(1 / A_pred) - r) < 1e-9, (g, r, A_pred)
        assert abs(ctr_pred - complex(ux, uy)) < 1e-9
    print("numeric cross-check of curvature/center formulas: OK")

def main():
    numeric_crosscheck()

    classes_small, lines_small = collect(orbit(60))
    classes, lines = collect(orbit(120))

    print("lines found (values of B up to sign):", sorted(lines))
    assert lines == {(0, 1), (0, -1)}, "only horizontal lines expected"

    NMAX = 15
    print(f"\n n | #orbit(60) #orbit(120) #congruence  formula   match")
    all_ok = True
    for n in range(1, NMAX + 1):
        pred = predicted(n)
        got60 = classes_small.get(n, set())
        got120 = classes.get(n, set())
        ok = (got120 == pred) and (got60 == pred) and (len(pred) == Ne_formula(n))
        all_ok &= ok
        print(f"{n:2d} | {len(got60):8d} {len(got120):10d} {len(pred):11d} "
              f"{Ne_formula(n):8d}   {'OK' if ok else 'MISMATCH'}")
    assert all_ok
    print("\nexact set equality orbit == congruence classes for n <= 15: OK")

    # parity claim across the whole computed orbit
    for n, cl in classes.items():
        for (x, y) in cl:
            assert x % 2 == 0 and y % 2 == 1
    print("parity zeta = i (mod 2) holds across entire computed orbit: OK")

    # local factor at 2: r(2^k) = 2^(k+1) for k >= 2
    for k in range(2, 12):
        q = 2 ** k
        cnt = sum(1 for x in range(q) for y in range(q) if (x * x + y * y) % q == 1)
        assert cnt == 2 ** (k + 1), (k, cnt)
    print("r(2^k) = 2^(k+1) for 2 <= k <= 11: OK")

    # congruence count == formula for larger n
    for n in range(1, 201):
        assert len(predicted(n)) == Ne_formula(n), n
    print("congruence count == n * prod(1 - chi(p)/p) for n <= 200: OK")

    # asymptotics: sum_{n<=X} N_e(n) * 2G / X^2 -> 1
    G = 0.9159655941772190
    X = 10 ** 6
    spf = list(range(X + 1))
    for i in range(2, int(X ** 0.5) + 1):
        if spf[i] == i:
            for j in range(i * i, X + 1, i):
                if spf[j] == j:
                    spf[j] = i
    Ne = [0] * (X + 1)
    Ne[1] = 1
    for n in range(2, X + 1):
        p = spf[n]
        m = n
        a = 0
        while m % p == 0:
            m //= p
            a += 1
        # multiplicative: Ne(p^a * m) = Ne(p^a) * Ne(m)
        if p == 2:
            loc = p ** a
        else:
            chi = 1 if p % 4 == 1 else -1
            loc = p ** (a - 1) * (p - chi)
        Ne[n] = loc * Ne[m]
    S = 0
    print("\n      X    sum_{n<=X} N_e(n) * 2G/X^2")
    targets = {10 ** k for k in range(2, 7)}
    for n in range(1, X + 1):
        S += Ne[n]
        if n in targets:
            print(f"{n:8d}   {S * 2 * G / n**2:.6f}")

if __name__ == "__main__":
    main()
