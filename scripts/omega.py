"""
The monoid Omega = { X in SL(2,Z[i]) : X(H) subset H }, H = upper half plane.

Conventions follow circle-classification.md:  a generalised circle is encoded by
a Hermitian matrix  M = [[A, B], [conj B, C]]  (A, C real, B complex),
      z on circle  <=>  A|z|^2 + conj(B) z + B conj(z) + C = 0,
and a *generalised disk* is the positive side  {z : A|z|^2 + 2 Re(conj(B) z) + C > 0}.

For X = [[a,b],[c,d]] in SL(2,Z[i]) the disk X(H) has matrix

      M(X) = (X^-1)^dagger M_0 X^-1 = [[ 2 Im(c conj d),  i(a conj d - b conj c) ],
                                       [        conj(.) ,  2 Im(a conj b)        ]]

with M_0 = [[0,i],[-i,0]] the matrix of H itself.  We write (A,B,C) = M(X) and
      n = -A/2,   zeta = B = x + iy,   m = -C/2 .
Then  X in Omega  <=>  A <= 0, C <= 0, Im B >= 1  (Proposition 1 of half-plane-monoid.md).

All Gaussian-integer arithmetic is exact.
"""

from fractions import Fraction
from math import isqrt


# ----------------------------------------------------------------- Gaussian integers
class G:
    """Exact Gaussian integer a + b i."""
    __slots__ = ("re", "im")

    def __init__(self, re=0, im=0):
        self.re = int(re)
        self.im = int(im)

    def __add__(s, o): o = _g(o); return G(s.re + o.re, s.im + o.im)
    __radd__ = __add__

    def __sub__(s, o): o = _g(o); return G(s.re - o.re, s.im - o.im)
    def __rsub__(s, o): return _g(o) - s

    def __mul__(s, o):
        o = _g(o)
        return G(s.re * o.re - s.im * o.im, s.re * o.im + s.im * o.re)
    __rmul__ = __mul__

    def __neg__(s): return G(-s.re, -s.im)
    def conj(s): return G(s.re, -s.im)
    def norm(s): return s.re * s.re + s.im * s.im          # |z|^2
    def __eq__(s, o): o = _g(o); return s.re == o.re and s.im == o.im
    def __hash__(s): return hash((s.re, s.im))
    def __complex__(s): return complex(s.re, s.im)

    def __repr__(s):
        if s.im == 0: return str(s.re)
        if s.re == 0: return ("i" if s.im == 1 else "-i" if s.im == -1 else f"{s.im}i")
        return f"{s.re}{'+' if s.im > 0 else '-'}{'' if abs(s.im)==1 else abs(s.im)}i"

    def divmod_g(s, o):
        """Gaussian division with remainder: s = q*o + r, |r| <= |o|/sqrt(2)."""
        o = _g(o)
        d = o.norm()
        p = s * o.conj()
        q = G(_round_half(p.re, d), _round_half(p.im, d))
        return q, s - q * o


def _g(o):
    return o if isinstance(o, G) else G(o, 0)


def _round_half(num, den):
    """Nearest integer to num/den (den > 0), ties away from zero."""
    q, r = divmod(num, den)
    return q + 1 if 2 * r >= den else q


ZERO, ONE, I = G(0, 0), G(1, 0), G(0, 1)


# ----------------------------------------------------------------- 2x2 matrices
class Mat:
    """2x2 matrix over Z[i] (or over Z[i]-Hermitian entries), det not enforced."""
    __slots__ = ("a", "b", "c", "d")

    def __init__(self, a, b, c, d):
        self.a, self.b, self.c, self.d = _g(a), _g(b), _g(c), _g(d)

    def __mul__(s, o):
        return Mat(s.a * o.a + s.b * o.c, s.a * o.b + s.b * o.d,
                   s.c * o.a + s.d * o.c, s.c * o.b + s.d * o.d)

    def det(s): return s.a * s.d - s.b * s.c
    def dagger(s): return Mat(s.a.conj(), s.c.conj(), s.b.conj(), s.d.conj())

    def inv(s):
        assert s.det() == ONE, "inverse only for det 1"
        return Mat(s.d, -s.b, -s.c, s.a)

    def __eq__(s, o): return (s.a, s.b, s.c, s.d) == (o.a, o.b, o.c, o.d)
    def __hash__(s): return hash((s.a, s.b, s.c, s.d))
    def __repr__(s): return f"[{s.a}, {s.b}; {s.c}, {s.d}]"

    def act(s, z):
        """Moebius action on a complex number (float) or on 'inf'."""
        a, b, c, d = (complex(s.a), complex(s.b), complex(s.c), complex(s.d))
        if z == "inf":
            return "inf" if c == 0 else a / c
        den = c * z + d
        if den == 0:
            return "inf"
        return (a * z + b) / den


ID = Mat(1, 0, 0, 1)
S = Mat(0, -1, 1, 0)                 # z -> -1/z
W = Mat(I, 0, 0, -I)                 # z -> -z   (maps H to the lower half plane)


def T(lam):
    """z -> z + lam."""
    return Mat(1, _g(lam), 0, 1)


# ----------------------------------------------------------------- circles / disks
class Disk:
    """The generalised disk {z : A|z|^2 + 2 Re(conj(B) z) + C > 0}, det = AC-|B|^2 = -1.

    Stored as (A, B, C) with A, C integers and B in Z[i].  For A < 0 this is the
    open disk of centre B/(-A) and radius 1/(-A); for A = 0 a half plane.
    """
    __slots__ = ("A", "B", "C")

    def __init__(self, A, B, C):
        self.A, self.B, self.C = int(A), _g(B), int(C)

    # normalised data of circle-classification.md
    @property
    def n(self): return -self.A // 2          # curvature/2   (A = -2n)
    @property
    def m(self): return -self.C // 2          # co-curvature/2
    @property
    def zeta(self): return self.B             # = 2n * centre
    @property
    def alpha(self): return self.B.im         # the invariant alpha = Im B = y

    def det(self): return self.A * self.C - self.B.norm()

    def is_halfplane(self): return self.A == 0

    def centre(self):
        assert self.A != 0
        return complex(Fraction(self.B.re, -self.A), Fraction(self.B.im, -self.A))

    def radius(self): return Fraction(1, -self.A) if self.A else None

    def height(self):
        """t, for the half plane {Im z > t}."""
        assert self.A == 0 and self.B == I
        return Fraction(-self.C, 2)

    def contains(self, z):
        """Is z (a float complex or an exact QC point) strictly inside?"""
        if self.A == 0:                       # B = i: 2 Im z + C > 0
            return 2 * z.imag + self.C > 0
        n = self.n                            # |2n z - zeta| < 1
        if hasattr(z, "abs2"):
            w = z * (2 * n) - QC(self.B.re, self.B.im)
            return w.abs2() < 1
        return abs(2 * n * z - complex(self.B)) < 1

    def value(self, z):
        return self.A * abs(z) ** 2 + 2 * (complex(self.B).conjugate() * z).real + self.C

    def __eq__(s, o):
        return isinstance(o, Disk) and (s.A, s.B, s.C) == (o.A, o.B, o.C)
    def __hash__(s): return hash((s.A, s.B, s.C))

    def __repr__(s):
        if s.A == 0:
            if s.B == I:
                return f"halfplane Im z > {s.height()}"
            return f"halfplane [A=0, B={s.B}, C={s.C}]"
        c, r = s.centre(), s.radius()
        return f"disk centre {c.real}{'+' if c.imag>=0 else '-'}{abs(c.imag)}i radius {r}"

    def pretty(s):
        if s.A == 0:
            return f"{{Im z > {s.height()}}}"
        c, r = s.centre(), s.radius()
        return f"|z - ({c.real}{'+' if c.imag>=0 else '-'}{abs(c.imag)}i)| < {r}"


H = Disk(0, I, 0)                     # the upper half plane itself
M0 = H


def disk_of(X):
    """M(X) = matrix of the disk X(H)."""
    a, b, c, d = X.a, X.b, X.c, X.d
    A = 2 * (c * d.conj()).im
    B = I * (a * d.conj() - b * c.conj())
    C = 2 * (a * b.conj()).im
    return Disk(A, B, C)


def in_omega(X):
    """X(H) subset H ?"""
    D = disk_of(X)
    return D.A <= 0 and D.C <= 0 and D.B.im >= 1


def inversive(D1, D2):
    """<M1,M2> = Re(B1 conj B2) - (A1 C2 + A2 C1)/2.  = -1 * inversive distance.

    Two disks are nested iff this is >= 1, have disjoint interiors iff <= -1,
    and their circles cross iff |.| < 1.
    """
    return (D1.B * D2.B.conj()).re - (D1.A * D2.C + D2.A * D1.C) // 2


def is_schmidt(D):
    """Membership test for the underlying (unoriented) circle, from
    circle-classification.md (Theorem):  curvature 2n, zeta = 2n*centre with
    zeta = i (mod 2) and |zeta|^2 = 1 (mod 4n)."""
    if D.det() != -1:
        return False
    A, B, C = (D.A, D.B, D.C) if D.A > 0 or (D.A == 0 and D.C >= 0) else (-D.A, -D.B, -D.C)
    if A % 2 or C % 2:
        return False
    if B.re % 2 != 0 or B.im % 2 != 1:          # zeta = i mod 2
        return False
    n = A // 2
    if n == 0:
        return B.norm() == 1                    # the lines Im z = k
    return (B.norm() - 1) % (4 * n) == 0


# ----------------------------------------------------------------- descent: disk -> matrix
def _pull_back(D, h):
    """M -> h^dagger M h : the matrix of the disk h^{-1}(D)."""
    MM = h.dagger() * Mat(D.A, D.B, D.B.conj(), D.C) * h
    assert MM.a.im == 0 and MM.d.im == 0
    return Disk(MM.a.re, MM.b, MM.d.re)


def matrix_for_disk(D, trace=False):
    """Return X in SL(2,Z[i]) with X(H) = D, by the reduction of
    circle-classification.md.  Raises ValueError if D is not a Schmidt disk."""
    if not is_schmidt(D):
        raise ValueError(f"{D} is not in the Schmidt arrangement")
    M, X, steps = D, ID, []
    while True:
        if M.A == 0:
            break
        q, _ = (-M.B).divmod_g(G(M.A))          # lam ~ -B/A  : shrink |B|
        lam = q
        if lam != ZERO:
            M, X = _pull_back(M, T(lam)), X * T(lam)
            steps.append(("T", lam))
        if M.C == 0:                            # tangent to R at 0: one more S
            M, X = _pull_back(M, S), X * S
            steps.append(("S", None))
            break
        M, X = _pull_back(M, S), X * S
        steps.append(("S", None))
    # now A = 0, B = +- i, C = -2t : the line Im z = t
    assert M.A == 0 and M.B.norm() == 1, M
    t = (-M.C // 2) if M.B == I else (M.C // 2)     # the line Im z = t
    if t != 0:
        lam = G(0, t)
        M, X = _pull_back(M, T(lam)), X * T(lam)
        steps.append(("T", lam))
    assert M.A == 0 and M.C == 0, M
    if M.B == -I:                               # wrong orientation: post-compose z -> -z
        M, X = _pull_back(M, W), X * W
        steps.append(("W", None))
    assert M == H, M
    if trace:
        return X, steps
    return X


# ----------------------------------------------------------------- atoms
def schmidt_disks_through(z, n):
    """All Schmidt disks of curvature 2n containing the point z (float complex)."""
    out = []
    if n == 0:
        t = 1
        while t <= z.imag:                      # {Im z > t} for t = 1 .. floor(Im z)
            out.append(Disk(0, I, -2 * t))
            t += 1
        return out
    from math import floor
    wx, wy = 2 * n * z.real, 2 * n * z.imag     # need |zeta - w| < 1
    for x in range(floor(wx) - 1, floor(wx) + 3):
        if x % 2:
            continue
        for y in range(floor(wy) - 1, floor(wy) + 3):
            if y % 2 == 0 or y < 1:
                continue
            if (x - wx) ** 2 + (y - wy) ** 2 >= 1:
                continue
            if (x * x + y * y - 1) % (4 * n):
                continue
            m = (x * x + y * y - 1) // (4 * n)
            out.append(Disk(-2 * n, G(x, y), -2 * m))
    return out


def maximal_disk_through(z, nmax=200):
    """The unique largest Schmidt disk inside H containing z -- i.e. the atom disk of z.

    Returns None if none is found with curvature <= 2*nmax, which happens exactly when
    z lies in (or extremely close to) the residual set of the gasket.
    """
    if z.imag > 1:
        return Disk(0, I, -2)                   # {Im z > 1}
    for n in range(1, nmax + 1):
        found = schmidt_disks_through(z, n)
        if found:
            assert len(found) == 1, (z, found)
            return found[0]
    return None


def atom_through(z, nmax=200, reduce_first=True):
    """An atom X of Omega with z in X(H) (None if the search fails).

    Because D(XU) = D(X) for units U, we may first move z into the standard
    fundamental domain of SL(2,Z) by a unit: there the atom has small curvature,
    so the search is O(1) instead of O(1/Im z).
    """
    U = ID
    if reduce_first:
        U, z = reduce_to_fundamental_domain(z)
    D = maximal_disk_through(z, nmax)
    if D is None:
        return None
    return U.inv() * matrix_for_disk(D)


def apollonian_chain(z, length, nmax=200):
    """Nested generalised disks D_1 > D_2 > ... > D_length, all containing z.

    D_k = X_1 X_2 ... X_k (H) with every X_j an atom of Omega, so that
    X_1 X_2 ... X_k is the (unique) atomic factorisation of the corresponding
    element of Omega.  Returns (list of (X_j, D_k), last point) and stops early
    if no atom is found (z in the residual set of the gasket).
    """
    out, Y, w = [], ID, z
    for _ in range(length):
        Xj = atom_through(w, nmax)
        if Xj is None:
            break
        Y = Y * Xj
        out.append((Xj, disk_of(Y)))
        w = act(Xj.inv(), w)
        if w is None or w.imag <= 0:
            break
    return out, w


# ----------------------------------------------------------------- exact search / maximality
def disks_through_exact(zx, zy, n):
    """All Schmidt disks of curvature 2n containing the exact rational point zx+i*zy."""
    out = []
    if n == 0:
        return [Disk(0, I, -2 * t) for t in range(1, int(zy) + 1) if zy > t]
    from math import floor
    wx, wy = 2 * n * zx, 2 * n * zy
    for x in range(floor(wx) - 1, floor(wx) + 3):
        if x % 2:
            continue
        for y in range(floor(wy) - 1, floor(wy) + 3):
            if y % 2 == 0 or y < 1:
                continue
            if (x - wx) ** 2 + (y - wy) ** 2 >= 1:
                continue
            if (x * x + y * y - 1) % (4 * n):
                continue
            out.append(Disk(-2 * n, G(x, y), -2 * ((x * x + y * y - 1) // (4 * n))))
    return out


def biggest_through_exact(zx, zy, nmax=400):
    if zy > 1:
        return Disk(0, I, -2)
    for n in range(1, nmax + 1):
        got = disks_through_exact(zx, zy, n)
        if got:
            assert len(got) == 1, (zx, zy, got)
            return got[0]
    return None


def is_maximal(D, nmax=400):
    """D subset D' forces centre(D) in D'; so D is maximal iff the largest Schmidt
    disk through its centre is D itself."""
    if D.A == 0:
        return D.C == -2                                   # {Im z > 1}
    zx, zy = Fraction(D.B.re, -D.A), Fraction(D.B.im, -D.A)
    return biggest_through_exact(zx, zy, nmax) == D


def all_disks_in_H(N, xrange):
    """All Schmidt disks contained in H with curvature 2n <= 2N, centre in a window."""
    out = [Disk(0, I, -2 * t) for t in range(1, N + 1)]
    for n in range(1, N + 1):
        for y in range(1, 4 * n + 2, 2):                   # y >= 1 <=> disk inside H
            for x in xrange(n):
                if (x * x + y * y - 1) % (4 * n) == 0:
                    out.append(Disk(-2 * n, G(x, y), -2 * ((x * x + y * y - 1) // (4 * n))))
    return out


# ---------------------------------------------------------------- the strip gasket
def gasket(N):
    """Apollonian gasket generated by the Descartes quadruple
       {Im z < 0}, {Im z > 1}, |z - i/2| = 1/2, |z - (1+i/2)| = 1/2,
    via the Lagarias-Mallows-Wilks move  M4' = 2(M1+M2+M3) - M4.
    Returns the disks of curvature <= 2N that lie inside H."""
    def move(q, j):
        s = [q[k] for k in range(4) if k != j]
        A = 2 * sum(d.A for d in s) - q[j].A
        B = 2 * (s[0].B + s[1].B + s[2].B) - q[j].B
        C = 2 * sum(d.C for d in s) - q[j].C
        return Disk(A, B, C)

    start = (Disk(0, -I, 0), Disk(0, I, -2), Disk(-2, I, 0), Disk(-2, G(2, 1), -2))
    seen, quads, found = set(), [start], set()
    while quads:
        nxt = []
        for q in quads:
            for j in range(4):
                D = move(q, j)
                if D.A >= 0 or -D.A > 2 * N:
                    continue
                if not (-2 * (-D.A) <= D.B.re <= 3 * (-D.A)):   # centre in [-2, 3]
                    continue
                key = (D.A, D.B.re, D.B.im, D.C)
                if key in seen:
                    continue
                seen.add(key)
                found.add(D)
                nq = tuple(D if k == j else q[k] for k in range(4))
                nxt.append(nq)
        quads = nxt
    for D in start:
        if D.A <= 0 and (D.A != 0 or D.C < 0):
            found.add(D)
    return found


def reduce_disk(D):
    """U in SL(2,Z) with U(D) of least curvature in its SL(2,Z)-orbit.

    This is Gauss reduction of the associated binary quadratic form
    f_D = (n, -x, m) (of discriminant 1 - alpha^2); the moves are z -> z + k and
    z -> -1/z, and the least curvature reached is 2*min(f_D).  Returns (U, U(D)).
    """
    if D.A == 0:
        return ID, D
    U, E = ID, D
    for _ in range(10000):
        if E.A == 0:                                   # curvature 0: already minimal
            return U, E
        if E.m < E.n:                                  # S : n <-> m, x -> -x
            U, E = S * U, _pull_back(E, S)
            continue
        k = -_nearest_int(Fraction(E.B.re, 2 * E.n))   # translate x into (-n, n]
        if k:
            U, E = T(k) * U, _pull_back(E, T(-k))
            continue
        return U, E
    raise RuntimeError("disk reduction did not terminate")


def maximal_disk_containing(D, nmax=400):
    """The unique atom disk containing the Schmidt disk D (D != H).

    Uses: D subset D' forces centre(D) in D', and two Schmidt disks never cross.
    """
    if D == H:
        return None
    if D.A == 0:                                   # {Im z > t}, t >= 1
        return Disk(0, I, -2)
    U, E = reduce_disk(D)                          # cost drops to O(min f_D)
    if E.A == 0:                                   # a half plane {Im z > t}, t >= 1
        G = Disk(0, I, -2)
    else:
        G = biggest_through_exact(Fraction(E.B.re, -E.A), Fraction(E.B.im, -E.A),
                                  max(nmax, E.n))
    if G is None:
        return None
    G = _pull_back(G, U)                           # U^{-1}(G)
    assert inversive(D, G) >= 1
    return G


def factor(X, nmax=400):
    """Atomic factorisation  X = A_1 A_2 ... A_k U  (A_j atoms of Omega, U a unit).

    Returns (list_of_atoms, unit).  The factorisation is unique up to associates.
    """
    assert in_omega(X), "X is not in Omega"
    atoms, Y = [], X
    while disk_of(Y) != H:
        G = maximal_disk_containing(disk_of(Y), nmax)
        if G is None:
            raise RuntimeError("no atom found within the curvature bound")
        A = matrix_for_disk(G)
        atoms.append(A)
        Y = A.inv() * Y
    return atoms, Y


def is_atom(X):
    """X is an atom of Omega  <=>  X(H) is a maximal proper Schmidt disk."""
    D = disk_of(X)
    return in_omega(X) and D != H and maximal_disk_containing(D) == D


# ----------------------------------------------------------------- exact rational points
class QC:
    """A point of Q(i) (exact rational complex number), or of any field via Fraction."""
    __slots__ = ("re", "im")

    def __init__(self, re=0, im=0):
        self.re, self.im = Fraction(re), Fraction(im)

    @staticmethod
    def of(z):
        if isinstance(z, QC): return z
        if isinstance(z, G): return QC(z.re, z.im)
        if isinstance(z, complex): raise TypeError("use floats or QC, not complex, here")
        return QC(z, 0)

    def __add__(s, o): o = QC.of(o); return QC(s.re + o.re, s.im + o.im)
    __radd__ = __add__
    def __sub__(s, o): o = QC.of(o); return QC(s.re - o.re, s.im - o.im)
    def __neg__(s): return QC(-s.re, -s.im)
    def __mul__(s, o):
        o = QC.of(o)
        return QC(s.re * o.re - s.im * o.im, s.re * o.im + s.im * o.re)
    __rmul__ = __mul__
    def conj(s): return QC(s.re, -s.im)
    def abs2(s): return s.re * s.re + s.im * s.im
    def __truediv__(s, o):
        o = QC.of(o); d = o.abs2()
        if d == 0: raise ZeroDivisionError
        p = s * o.conj()
        return QC(p.re / d, p.im / d)
    def __eq__(s, o): o = QC.of(o); return s.re == o.re and s.im == o.im
    def __hash__(s): return hash((s.re, s.im))
    def __complex__(s): return complex(float(s.re), float(s.im))
    def __repr__(s): return f"{s.re}{'+' if s.im >= 0 else '-'}{abs(s.im)}i"

    @property
    def imag(self): return self.im
    @property
    def real(self): return self.re
    def __abs__(s): return float(s.abs2()) ** 0.5


def act(X, z):
    """Moebius action, exact on QC points, floating point on complex ones."""
    if isinstance(z, QC):
        a, b, c, d = QC.of(X.a), QC.of(X.b), QC.of(X.c), QC.of(X.d)
        den = c * z + d
        if den == QC(0, 0):
            return None                      # infinity
        return (a * z + b) / den
    return X.act(z)


def reduce_to_fundamental_domain(z, maxit=10000):
    """U in SL(2,Z) with U z in the standard fundamental domain of SL(2,Z)
    (|Re| <= 1/2, |z| >= 1).  Exact for QC input.  Returns (U, U z)."""
    U, w = ID, z
    for _ in range(maxit):
        k = -_nearest_int(w.real)
        if k:
            U, w = T(k) * U, act(T(k), w)
        if (w.re * w.re + w.im * w.im if isinstance(w, QC) else abs(w) ** 2) < 1:
            U, w = S * U, act(S, w)
        else:
            return U, w
    raise RuntimeError("reduction did not terminate")


def _nearest_int(x):
    from math import floor
    return int(floor(x + Fraction(1, 2) if isinstance(x, Fraction) else x + 0.5))
