"""Machine verification of the claims in half-plane-monoid.md."""
import sys, random, itertools
from fractions import Fraction
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from omega import *
from omega import _pull_back

random.seed(20260827)
OK = []


def check(name, cond, extra=""):
    OK.append(cond)
    print(("  ok  " if cond else " FAIL ") + name + (("   " + extra) if extra else ""))


# ---------------------------------------------------------------- random group elements
GENS = [T(1), T(I), T(-1), T(-I), S, S * S * S]


def random_element(k=8):
    X = ID
    for _ in range(k):
        X = X * random.choice(GENS)
    return X


def maps_H_into_H_numeric(X, trials=400):
    """Sample points of H and test Im(Xz) > 0 (with tangency tolerance)."""
    for _ in range(trials):
        z = complex(random.uniform(-6, 6), random.expovariate(0.5) + 1e-9)
        w = X.act(z)
        if w == "inf" or w.imag <= -1e-12:
            return False
    return True


print("1. the criterion  X in Omega  <=>  A <= 0, C <= 0, Im B >= 1")
bad = 0
for _ in range(400):
    X = random_element(random.randint(1, 10))
    if in_omega(X) != maps_H_into_H_numeric(X):
        bad += 1
check("criterion agrees with numerical sampling on 400 random elements", bad == 0, f"{bad} mismatches")
check("criterion = Re(a conj(d) - b conj(c)) >= 1 form",
      all((lambda D: (D.A <= 0 and D.C <= 0 and D.B.im >= 1) ==
                     (D.A <= 0 and D.C <= 0 and (X.a * X.d.conj() - X.b * X.c.conj()).re >= 1))(disk_of(X))
          for X in (random_element(6) for _ in range(300))))
def cartan(X):
    """X * conj(X)^{-1}."""
    return X * Mat(X.a.conj(), X.b.conj(), X.c.conj(), X.d.conj()).inv()


check("Cartan embedding: X conj(X)^-1 = [[-iB,-iC],[iA,iB^-]], so tr = 2 Im B = 2 alpha",
      all((lambda X, D, Y: Y == Mat(-I * D.B, -I * D.C, I * D.A, I * D.B.conj())
                           and (Y.a + Y.d) == G(2 * D.B.im, 0))(X, disk_of(X), cartan(X))
          for X in (random_element(6) for _ in range(200))))

print("\n2. Omega is a monoid; its units are exactly SL(2,Z)")
els = [random_element(random.randint(1, 9)) for _ in range(200)]
om = [X for X in els if in_omega(X)]
check(f"closed under multiplication ({len(om)} elements of Omega sampled)",
      all(in_omega(X * Y) for X in om for Y in om))
check("identity in Omega", in_omega(ID))
units = [X for X in om if in_omega(X.inv())]
check("X and X^-1 both in Omega  <=>  entries are rational integers",
      all((X.a.im == X.b.im == X.c.im == X.d.im == 0) for X in units))
check("every SL(2,Z) element is a unit of Omega",
      all(in_omega(U) and in_omega(U.inv())
          for U in [Mat(1, 1, 0, 1), Mat(0, -1, 1, 0), Mat(2, 1, 1, 1), Mat(-3, 5, 4, -7), Mat(1, 0, 5, 1)]))
check("i*GL(2,Z)^{det=-1} is NOT in Omega (it maps H to the lower half plane)",
      not any(in_omega(Mat(I * a, I * b, I * c, I * d))
              for (a, b, c, d) in [(1, 0, 0, -1), (0, 1, 1, 0), (2, 1, 1, 0)]))

print("\n3. the Schmidt arrangement is a packing: <M,M'> is always odd")
disks = [disk_of(X) for X in els]
prods = [inversive(D1, D2) for D1, D2 in itertools.combinations(disks, 2)]
check(f"inversive product odd for all {len(prods)} pairs", all(p % 2 == 1 for p in prods))
check("hence |<M,M'>| >= 1: no two Schmidt circles cross transversally",
      all(abs(p) >= 1 for p in prods))
check("all sampled circles are Schmidt-classified", all(is_schmidt(D) for D in disks))


print("\n4. maximal disks = the Apollonian strip gasket")
N = 14
gk = {D for D in gasket(N)}
check("every gasket circle is in the Schmidt arrangement", all(is_schmidt(D) for D in gk),
      f"{len(gk)} gasket disks of curvature <= {2*N} inside H")
check("gasket disks are pairwise non-overlapping (inversive product <= -1)",
      all(inversive(D1, D2) <= -1 for D1, D2 in itertools.combinations(list(gk), 2)))
check("every gasket disk is a maximal element of D", all(is_maximal(D, N) for D in gk))

strip = all_disks_in_H(N, lambda n: range(0, 2 * n + 1, 2))
maxima = [D for D in strip if is_maximal(D, N)]
gk_win = {D for D in gk if D.A != 0 and 0 <= D.B.re <= -D.A}
mx_win = {D for D in maxima if D.A != 0 and 0 <= D.B.re <= -D.A}
check(f"brute-force maximal disks = gasket disks (curvature <= {2*N}, one period 0 <= x <= 2n)",
      gk_win == mx_win,
      f"{len(mx_win)} maximal vs {len(gk_win)} gasket; "
      f"extra={sorted((D.n, D.B.re, D.B.im) for D in mx_win - gk_win)[:4]} "
      f"missing={sorted((D.n, D.B.re, D.B.im) for D in gk_win - mx_win)[:4]}")
print(f"     ({len(strip)} Schmidt disks in H with curvature <= {2*N}, 0 <= x <= 2n; "
      f"{len(maxima)} maximal)")


print("\n5. atomic factorisation: existence, uniqueness up to associates, additive length")
sample = [X for X in (random_element(random.randint(1, 9)) for _ in range(400)) if in_omega(X)]
facs = [factor(X) for X in sample]
def prod(l, U=ID):
    P = ID
    for A in l: P = P * A
    return P * U
check(f"X = A_1...A_k U reproduces X for all {len(sample)} sampled X in Omega",
      all(prod(a, u) == X for X, (a, u) in zip(sample, facs)))
check("every factor returned is an atom", all(is_atom(A) for a, _ in facs for A in a))
check("the cofactor is a unit", all(in_omega(u) and in_omega(u.inv()) for _, u in facs))
check("X is a unit  <=>  k = 0", all((len(a) == 0) == in_omega(X.inv()) for X, (a, _) in zip(sample, facs)))
check("X is an atom <=>  k = 1", all((len(a) == 1) == is_atom(X) for X, (a, _) in zip(sample, facs)))

pairs = [(X, Y) for X in sample[:40] for Y in sample[:40]]
check(f"length is additive: L(XY) = L(X) + L(Y) on {len(pairs)} pairs",
      all(len(factor(X * Y)[0]) == len(factor(X)[0]) + len(factor(Y)[0]) for X, Y in pairs))

# uniqueness up to associates: build a product of atoms with random units interleaved
atoms_pool = [matrix_for_disk(D) for D in list(gasket(6))[:12]]
units_pool = [Mat(1, 1, 0, 1), Mat(1, 0, 1, 1), Mat(0, -1, 1, 0), ID, Mat(1, -1, 0, 1)]
uniq = True
for _ in range(200):
    k = random.randint(1, 4)
    chosen = [random.choice(atoms_pool) * random.choice(units_pool) for _ in range(k)]
    X = prod(chosen)
    a, u = factor(X)
    if len(a) != k:
        uniq = False; break
    # partial products must define the same nested disks
    P1 = P2 = ID
    for j in range(k):
        P1, P2 = P1 * chosen[j], P2 * a[j]
        if disk_of(P1) != disk_of(P2):
            uniq = False; break
check("any product of k atoms factors back into exactly k atoms, with the same "
      "chain of disks (uniqueness up to associates)", uniq)

print("\n6. the algorithm: atoms through a point, and nested chains")
pts = [complex(random.uniform(-2, 2), random.uniform(0.02, 3)) for _ in range(60)]
chains = [apollonian_chain(z, 6) for z in pts]
check("every disk of every chain contains z",
      all(all(D.contains(z) for _, D in ch) for z, (ch, _) in zip(pts, chains)))
check("the chains are strictly nested (inversive product >= 1, all distinct)",
      all(all(inversive(ch[j + 1][1], ch[j][1]) >= 1 and ch[j + 1][1] != ch[j][1]
              for j in range(len(ch) - 1)) for ch, _ in chains))
check("every disk of every chain is a Schmidt disk inside H",
      all(is_schmidt(D) and D.A <= 0 and D.C <= 0 and D.B.im >= 1 for ch, _ in chains for _, D in ch))
check("each step is an atom", all(is_atom(X) for ch, _ in chains for X, _ in ch))
check("the intermediate points X_j^-1 ... X_1^-1 z stay in H",
      all(w.imag > 0 for _, w in chains))
lens = [len(ch) for ch, _ in chains]
print(f"     (chain lengths reached with curvature bound 400: min {min(lens)}, max {max(lens)})")



print("\n7. invariants of atoms: alpha, quadratic forms, the Ford stratum")
def gcd(a, b):
    import math
    return math.gcd(a, b)

atoms7 = [matrix_for_disk(D) for D in list(gasket(20))]
us = [Mat(1, 1, 0, 1), Mat(0, -1, 1, 0), Mat(1, 0, -1, 1), Mat(2, 1, 1, 1), ID]
check("alpha(X) = Im B(X) is a two-sided SL(2,Z)-class invariant",
      all(disk_of(u * X * v).B.im == disk_of(X).B.im for X in atoms7 for u in us for v in us))
check("alpha(X) = tr(X conj(X)^-1)/2, and X in Omega forces alpha >= 1",
      all(disk_of(X).B.im >= 1 for X in atoms7))

fords = [D for D in (disk_of(X) for X in atoms7) if D.A != 0 and D.B.im == 1]
def is_ford(D):
    from math import isqrt
    n, x = D.n, D.B.re
    q = isqrt(n)
    return q * q == n and x % (2 * q) == 0 and gcd(abs(x // (2 * q)), q) == 1
check(f"every atom disk with alpha = 1 is a Ford disk ({len(fords)} tested)",
      all(is_ford(D) for D in fords))
def sl2_to(p, q):
    """U in SL(2,Z) with U(infinity) = p/q."""
    from math import gcd as _g
    assert _g(abs(p), abs(q)) == 1
    # solve p*s - q*r = 1
    old_r, r = p, q
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r:
        quo = old_r // r
        old_r, r = r, old_r - quo * r
        old_s, s = s, old_s - quo * s
        old_t, t = t, old_t - quo * t
    if old_r == -1:
        old_s, old_t = -old_s, -old_t
    return Mat(p, -old_t, q, old_s)          # det = p*old_s + q*old_t = 1

def ford_pq(D):
    from math import isqrt
    q = isqrt(D.n); return (D.B.re // (2 * q), q)
check("every Ford disk equals U({Im z > 1}) for an explicit U in SL(2,Z); so its matrix "
      "lies in SL(2,Z) T_i SL(2,Z)",
      all(_pull_back(Disk(0, I, -2), sl2_to(*ford_pq(D)).inv()) == D for D in fords))

# containment in a Ford disk <=> f_D(p,q) <= (alpha-1)/2, tangency <=> = (alpha+1)/2
def ford(p, q):
    return Disk(-2 * q * q, G(2 * p * q, 1), -2 * p * p)
tests = []
for D in [d for d in all_disks_in_H(9, lambda n: range(-2 * n, 2 * n + 1, 2))][:400]:
    for (p, q) in [(0, 1), (1, 1), (1, 2), (2, 3), (-1, 3), (1, 0), (3, 4)]:
        if gcd(abs(p), abs(q)) != 1:
            continue
        f = D.n * p * p - D.B.re * p * q + D.m * q * q
        tests.append((inversive(D, ford(p, q)) >= 1) == (2 * f <= D.B.im - 1))
        tests.append((inversive(D, ford(p, q)) == -1) == (2 * f == D.B.im + 1))
check(f"D subset Ford(p/q) <=> f_D(p,q) <= (alpha-1)/2, and tangency <=> f_D(p,q) = (alpha+1)/2"
      f"  [{len(tests)} instances]", all(tests))

def fD(D, p, q):
    return D.n * p * p - D.B.re * p * q + D.m * q * q

eq = True
for D in [d for d in all_disks_in_H(7, lambda n: range(-2 * n, 2 * n + 1, 2))][:60]:
    for U in [Mat(1, 1, 0, 1), Mat(1, 0, 1, 1), Mat(0, -1, 1, 0), Mat(2, 1, 1, 1), Mat(3, 2, 4, 3)]:
        E = _pull_back(D, U.inv())                      # E = U(D)
        p_, r_, q_, s_ = U.a.re, U.b.re, U.c.re, U.d.re
        for (v, w) in [(1, 0), (0, 1), (1, 1), (2, -3)]:
            if fD(E, v, w) != fD(D, s_ * v - r_ * w, -q_ * v + p_ * w):
                eq = False
check("D -> f_D is SL(2,Z)-equivariant: f_{U(D)}(v) = f_D(U^{-1}v); since (alpha, f_D) "
      "determines D, associate classes of atoms <-> classes of forms of disc 1-alpha^2", eq)

check("Gauss reduction of f_D realises the least curvature in the SL(2,Z)-orbit of D",
      all((lambda D: (lambda U, E: all(-_pull_back(D, (V).inv()).A >= -E.A
                                       for V in [u1 * u2 for u1 in us for u2 in us]))(*reduce_disk(D)))(D)
          for D in [d for d in all_disks_in_H(6, lambda n: range(0, 2 * n + 1, 2))][:120]))

D_bad = Disk(-32, G(-16, 31), -38)           # radius 1/32, tangent to Im z = 1 at -1/2
check("the disk of radius 1/32 tangent to Im z = 1 at -1/2 is a Schmidt disk in H",
      is_schmidt(D_bad) and D_bad.B.im >= 1 and D_bad.A <= 0 and D_bad.C <= 0)
check("...it lies in no Ford disk (min f = 16 > (alpha-1)/2 = 15) yet is NOT an atom: "
      "it sits inside the atom of curvature 8",
      (not any(inversive(D_bad, ford(p, q)) >= 1
               for q in range(0, 40) for p in range(-40, 41) if gcd(abs(p), abs(q)) == 1))
      and maximal_disk_containing(D_bad) == Disk(-8, G(-4, 7), -8))

print("\n8. the involution sigma(X) = conj(X)^-1 is an anti-automorphism of Omega")
def sigma(X):
    return Mat(X.a.conj(), X.b.conj(), X.c.conj(), X.d.conj()).inv()

om8 = [X for X in (random_element(random.randint(1, 9)) for _ in range(300)) if in_omega(X)]
check(f"sigma(Omega) = Omega ({len(om8)} elements tested)",
      all(in_omega(sigma(X)) for X in om8) and all(sigma(sigma(X)) == X for X in om8))
check("sigma reverses products: sigma(XY) = sigma(Y) sigma(X)",
      all(sigma(X * Y) == sigma(Y) * sigma(X) for X in om8[:30] for Y in om8[:30]))
check("sigma preserves alpha and the factorisation length",
      all(disk_of(sigma(X)).B.im == disk_of(X).B.im and
          len(factor(sigma(X))[0]) == len(factor(X)[0]) for X in om8))
check("sigma permutes the atoms", all(is_atom(sigma(matrix_for_disk(D)))
                                      for D in list(gasket(12))))

print("\n9. the alpha-spectrum of the atoms is an Apollonian orbit")
seen9, frontier, vals = set(), [(-1, 1, 1, 1)], {-1, 1}
for _ in range(10):
    nxt = []
    for q in frontier:
        for j in range(4):
            new = 2 * (sum(q) - q[j]) - q[j]
            if abs(new) > 600:
                continue
            nq = tuple(new if k == j else q[k] for k in range(4))
            if nq in seen9:
                continue
            seen9.add(nq); vals.add(new); nxt.append(nq)
    frontier = nxt
gk_alpha = sorted({D.B.im for D in gasket(250) if D.A != 0})
check("alpha-values of the atoms = positive entries of the orbit of (-1,1,1,1) under "
      "alpha_j -> 2(sum of the others) - alpha_j",
      sorted(v for v in vals if 0 < v <= 500) == [a for a in gk_alpha if a <= 500],
      f"{[a for a in gk_alpha if a <= 500]}")

print("\n" + ("ALL %d CHECKS PASSED" % len(OK) if all(OK) else "SOME CHECKS FAILED"))
