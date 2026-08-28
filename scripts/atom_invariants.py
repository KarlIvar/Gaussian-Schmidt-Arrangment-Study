"""Invariants of the atoms of Omega.

Each atom disk D carries
  * alpha(D) = Im B = tr(X conj(X)^-1)/2  (the invariant alpha of hyperbolic-counting.md),
  * the positive semidefinite binary quadratic form f_D = (n, -x, m) of discriminant
    1 - alpha^2 (see involution.md),
and the pair (alpha, [f_D]) is a complete invariant of the SL(2,Z)-orbit (= of the
associate class SL(2,Z) X SL(2,Z) of the atom).  This script tabulates which pairs
actually occur among the atoms, i.e. among the disks of the Apollonian gasket.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from omega import gasket


def reduce_form(a, b, c):
    """Gauss reduction of a positive (semi)definite form a x^2 + b xy + c y^2."""
    if b * b == 4 * a * c:
        return ("degenerate",)                       # alpha = 1: the Ford class
    while True:
        if a > c:
            a, b, c = c, -b, a; continue
        if b <= -a:
            k = (a - b) // (2 * a); b, c = b + 2 * k * a, a * k * k + b * k + c; continue
        if b > a:
            k = -((b + a - 1) // (2 * a)); b, c = b + 2 * k * a, a * k * k + b * k + c; continue
        break
    if b < 0 and (a == c or -b == a):
        b = -b
    return (a, b, c)


def main(curv=500):
    disks = [D for D in gasket(curv // 2) if D.A != 0]
    table = {}
    for D in disks:
        table.setdefault(D.B.im, set()).add(reduce_form(D.n, -D.B.re, D.m))
    print(f"atoms of curvature <= {curv} with centre in [-2,3]:  {len(disks)} disks\n")
    print(f"{'alpha':>6} {'(alpha-1)/2':>12} {'#classes':>9}  {'min f':>7}   reduced forms")
    for y in sorted(table):
        fs = sorted(table[y])
        mins = sorted({1 if f[0] == "degenerate" else f[0] for f in fs})
        show = ", ".join("Ford" if f[0] == "degenerate" else str(f) for f in fs[:4])
        if len(fs) > 4:
            show += ", ..."
        print(f"{y:>6} {(y-1)//2:>12} {len(fs):>9}  {str(mins)[1:-1]:>7}   {show}")
    print("\n(alpha = 2q^2-1 comes from the atoms tangent to the line Im z = 1 at p/q;")
    print(" the remaining alpha are carried by atoms tangent to neither boundary line.)")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 500)
