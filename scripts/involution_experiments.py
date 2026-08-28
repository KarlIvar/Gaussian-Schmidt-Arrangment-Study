"""Experiments for the involution sigma(X) = conj(X)^{-1} on SL_2(Z[i]).

Verifies (exactly, in integer arithmetic):
  (1) sigma(X) = J^{-1} X^dagger J  (adjoint anti-involution of M0 = iJ);
  (2) X * conj(X)^{-1} = -i J conj(M_X)   where M_X is the Hermitian matrix
      of the circle X(Rhat)  (the "fun fact" / Cartan embedding);
  (3) fixed points of sigma are exactly X = M0 M with M integral Hermitian
      of determinant -1 (matrices of circle inversions);
  (4) M_{sigma(X)} = -conj(M_{X^{-1}});
  (5) the induced map on SL_2(Z)-classes of circles: with a circle recorded
      as (form f = (q,-x,m), height y), sigma sends (f, y) -> (?, y) and the
      question is whether ? = f or f^op (class-group inversion) — tested on
      random X and on a targeted order-4 class of discriminant -80;
  (6) classification of ALL integral binary Hermitian forms of det -1 into
      SL_2(Z[i])-unitary classes (= twisted Galois classes): BFS from seeds.
"""
import random
from collections import deque

# ---------- Gaussian integers as (re, im) tuples ----------
def gm(p, q): a, b = p; c, d = q; return (a*c - b*d, a*d + b*c)
def ga(p, q): return (p[0] + q[0], p[1] + q[1])
def gneg(p): return (-p[0], -p[1])
def gc(p): return (p[0], -p[1])

def matmul(X, Y):
    return tuple(tuple(ga(gm(X[i][0], Y[0][j]), gm(X[i][1], Y[1][j]))
                       for j in range(2)) for i in range(2))
def dagger(X):
    return ((gc(X[0][0]), gc(X[1][0])), (gc(X[0][1]), gc(X[1][1])))
def conjm(X):
    return tuple(tuple(gc(X[i][j]) for j in range(2)) for i in range(2))
def det(X):
    return ga(gm(X[0][0], X[1][1]), gneg(gm(X[0][1], X[1][0])))
def inv_sl2(X):
    (a, b), (c, d) = X
    assert det(X) == (1, 0)
    return ((d, gneg(b)), (gneg(c), a))
def scal(s, X):  # multiply matrix by Gaussian scalar s
    return tuple(tuple(gm(s, X[i][j]) for j in range(2)) for i in range(2))

I2 = (((1, 0), (0, 0)), ((0, 0), (1, 0)))
J = (((0, 0), (1, 0)), (((-1), 0), (0, 0)))
M0 = (((0, 0), (0, 1)), ((0, -1), (0, 0)))          # iJ
S = (((0, 0), (-1, 0)), ((1, 0), (0, 0)))
G0 = (((0, 1), (0, 0)), ((0, 0), (0, -1)))          # diag(i, -i): z -> -z

def T(mu):  # translation z -> z + mu, mu = (re, im)
    return (((1, 0), mu), ((0, 0), (1, 0)))

GENS = [T((1, 0)), T((-1, 0)), T((0, 1)), T((0, -1)), S, inv_sl2(S)]

def sigma(X):
    return conjm(inv_sl2(X))

def M_of_X(X):
    """Hermitian matrix of the circle X(Rhat): (X^{-1})^dagger M0 X^{-1}."""
    Xi = inv_sl2(X)
    return matmul(dagger(Xi), matmul(M0, Xi))

def herm_triple(M):
    """(A, B, C) with A, C integers, B Gaussian."""
    (A, B), (Bc, C) = M
    assert A[1] == 0 and C[1] == 0 and Bc == gc(B)
    return A[0], B, C[0]

def herm_from_triple(A, B, C):
    return (((A, 0), B), ((gc(B)), (C, 0)))

def random_X(rng, length=None):
    X = I2
    for _ in range(length or rng.randint(4, 12)):
        X = matmul(X, rng.choice(GENS))
    return X

# ---------- circle -> (form, y), Gauss reduction ----------
def to_form_y(M):
    """Unoriented circle data: normalize +-M so the curvature is positive.
    Returns (q, x, m, y): radius 1/(2q), center (x + yi)/(2q), form (q,-x,m).
    None for lines (A = 0)."""
    A, B, C = herm_triple(M)
    if A == 0:
        return None
    if A < 0:
        A, B, C = -A, gneg(B), -C
    x, y = -B[0], -B[1]
    assert A % 2 == 0 and C % 2 == 0
    return A // 2, x, C // 2, y

def reduce_form(a, b, c):
    """Gauss reduction of a positive definite integral form."""
    assert a > 0 and b*b - 4*a*c < 0
    while True:
        # normalize: -a < b <= a
        if not (-a < b <= a):
            k = (a - b) // (2 * a)        # b + 2ak in (-a, a]
            b, c = b + 2*a*k, c + b*k + a*k*k
        if a > c:
            a, b, c = c, -b, a
            continue
        break
    if b < 0 and (a == -b or a == c):
        b = -b
    return (a, b, c)

def opp(f):
    a, b, c = f
    return reduce_form(a, -b, c)

# ---------- build X realizing a prescribed circle (descent) ----------
def matrix_for_circle(Mtarget):
    """X in SL_2(Z[i]) with M_of_X(X) == Mtarget, for Mtarget in the Schmidt
    family (integral Hermitian, det -1, off-diagonal = i mod 2)."""
    M = Mtarget
    X = I2
    def apply(h):
        nonlocal M, X
        M = matmul(dagger(h), matmul(M, h))
        X = matmul(X, h)
    for _ in range(10_000):
        A, B, C = herm_triple(M)
        if A != 0:
            # translate: B -> B + A*mu, choose mu = -round(B/A)
            mu = (-round(B[0] / A), -round(B[1] / A))
            if mu != (0, 0):
                apply(T(mu))
                continue
            A, B, C = herm_triple(M)
            if C == 0 and abs(A) > 1:   # |B| = 1 here; invert to reach a line
                apply(S)
                continue
            if abs(C) < abs(A) or C == 0:
                apply(S)
                continue
            apply(S)
        else:
            if B not in ((0, 1), (0, -1)):
                raise ValueError("not in the Schmidt family")
            if B == (0, -1):
                apply(G0)      # flips B
                continue
            if C != 0:
                # line (0, i, C): C -> C + 2 Im(mu); pick mu = -C/2 * i
                apply(T((0, -C // 2)))
                continue
            break
    assert M == M0
    # M0 = (X^{-1} ... ), with our step order Mtarget = (X^{-1})^† M0 X^{-1}
    Xr = X
    assert M_of_X(Xr) == Mtarget
    return Xr

# ---------- experiments ----------
def exp_identities(trials=500, seed=5):
    rng = random.Random(seed)
    minus_i = (0, -1)
    for _ in range(trials):
        X = random_X(rng)
        # (1) sigma(X) = J^{-1} X^dagger J
        Jinv = scal((-1, 0), J)   # J^{-1} = -J
        assert sigma(X) == matmul(Jinv, matmul(dagger(X), J))
        # (2) X conj(X)^{-1} = -i J conj(M_X)
        Y = matmul(X, conjm(inv_sl2(X)))
        assert Y == scal(minus_i, matmul(J, conjm(M_of_X(X))))
        # (4) M_{sigma X} = -conj(M_{X^{-1}})
        lhs = M_of_X(sigma(X))
        rhs = scal((-1, 0), conjm(M_of_X(inv_sl2(X))))
        assert lhs == rhs
        # (3) X0 = M0 * M is sigma-fixed for any Hermitian det -1 M
        M = M_of_X(X)
        X0 = matmul(M0, M)
        assert det(X0) == (1, 0) and sigma(X0) == X0
    print(f"identities (1)-(4) on {trials} random X: OK")

def exp_class_map(trials=4000, seed=7):
    rng = random.Random(seed)
    stats = dict(total=0, y_eq=0, f_eq=0, f_op=0, distinguishing=0,
                 dist_f_eq=0, dist_f_op=0)
    for _ in range(trials):
        X = random_X(rng)
        d = to_form_y(M_of_X(X))
        if d is None:
            continue
        q, x, m, y = d
        if y * y <= 1:
            continue                      # degenerate / indefinite
        d2 = to_form_y(M_of_X(sigma(X)))
        q2, x2, m2, y2 = d2
        f = reduce_form(q, -x, m)
        f2 = reduce_form(q2, -x2, m2)
        stats['total'] += 1
        stats['y_eq'] += (y2 == y)
        stats['f_eq'] += (f2 == f)
        stats['f_op'] += (f2 == opp(f))
        if opp(f) != f:                   # a class of order > 2: distinguishes
            stats['distinguishing'] += 1
            stats['dist_f_eq'] += (f2 == f)
            stats['dist_f_op'] += (f2 == opp(f))
    print("class-map statistics over random X:", stats)
    return stats

def exp_targeted():
    # discriminant -80 (alpha = 9): class group Z/4, classes of (3, +-2, 7)
    # are inverse to each other and NOT ambiguous.
    for f in [(3, 2, 7), (3, -2, 7)]:
        a, b, c = f
        x, y = -b, 9
        M = herm_from_triple(2*a, (-x, -y), 2*c)
        X = matrix_for_circle(M)
        q2, x2, m2, y2 = to_form_y(M_of_X(sigma(X)))
        f2 = reduce_form(q2, -x2, m2)
        print(f"  form {f} (y=9)  --sigma-->  form {f2} (y={y2});   "
              f"f reduced = {reduce_form(*f)}, opposite = {opp(f)}")

def exp_twisted_classes(LIM=40):
    """Classify integral Hermitian det -1 forms under M -> h^dagger M h,
    h in SL_2(Z[i]) (the twisted / unitary equivalence)."""
    def orbit(seed):
        seen = {seed}
        dq = deque([seed])
        while dq:
            M = dq.popleft()
            for h in GENS + [G0]:
                M2 = matmul(dagger(h), matmul(M, h))
                A, B, C = herm_triple(M2)
                if abs(A) <= LIM and abs(C) <= LIM and M2 not in seen:
                    seen.add(M2)
                    dq.append(M2)
        return seen

    # all integral Hermitian det -1 in the box
    allM = set()
    for A in range(-LIM, LIM + 1):
        for C in range(-LIM, LIM + 1):
            n2 = 1 + A * C
            if n2 < 0:
                continue
            for bx in range(-LIM - 1, LIM + 2):
                r = n2 - bx * bx
                if r < 0:
                    continue
                by = int(round(r ** 0.5))
                for byy in {by, -by}:
                    if bx*bx + byy*byy == n2:
                        allM.add(herm_from_triple(A, (bx, byy), C))
    seeds = [("S (real line, B=i)", M0),
             ("iS (vertical lines, B=1)", herm_from_triple(0, (1, 0), 0)),
             ("iS reversed (B=-1)", herm_from_triple(0, (-1, 0), 0)),
             ("unit circle (1,0,-1)", herm_from_triple(1, (0, 0), -1)),
             ("unit circle reversed", herm_from_triple(-1, (0, 0), 1)),
             ("odd (1,1+i,1)", herm_from_triple(1, (1, 1), 1)),
             ("odd (-1,1+i,-1)", herm_from_triple(-1, (1, 1), -1))]
    orbits = []
    covered = set()
    for name, seedM in seeds:
        if seedM in covered:
            print(f"  seed {name}: already covered by a previous class")
            continue
        o = orbit(seedM)
        orbits.append((name, o))
        covered |= o
    # restrict comparison to a safe inner box (BFS is complete well inside LIM)
    inner = LIM // 2
    def in_inner(M):
        A, B, C = herm_triple(M)
        return abs(A) <= inner and abs(C) <= inner
    left = [M for M in allM if in_inner(M) and M not in covered]
    print(f"  det -1 Hermitian forms with |A|,|C| <= {inner}: "
          f"{sum(1 for M in allM if in_inner(M))}; "
          f"classes found: {len(orbits)}; unassigned: {len(left)}")
    for name, o in orbits:
        cnt = sum(1 for M in o if in_inner(M))
        print(f"    class {name}: {cnt} forms in the inner box")
    if left:
        print("    UNASSIGNED examples:", [herm_triple(M) for M in left[:5]])

if __name__ == "__main__":
    exp_identities()
    exp_class_map()
    print("targeted test (discriminant -80, class group Z/4):")
    exp_targeted()
    print("twisted (unitary) classes of det -1 integral Hermitian forms:")
    exp_twisted_classes()
