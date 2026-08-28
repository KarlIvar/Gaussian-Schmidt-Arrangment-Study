#!/usr/bin/env python3
"""Nested Schmidt disks around a point of the upper half plane.

  python3 scripts/apollonian_chain_demo.py "31/101 + 7/53 i" 12
  python3 scripts/apollonian_chain_demo.py "0.3141 + 0.2718i" 15 --float

Given z in H the program repeatedly finds an *atom* X of the monoid
    Omega = { X in SL(2,Z[i]) : X(H) subset H }
with z in X(H), and returns the nested generalised disks

    H  >  D_1 = X_1(H)  >  D_2 = X_1X_2(H)  >  ...  >  D_n,

all of which contain z.  X_1 X_2 ... X_n is then the atomic factorisation of an
element of Omega whose disk contains z.  Rational input is handled in exact
arithmetic; the search fails (and the chain stops) exactly when the current
point lies in the residual set of the Apollonian gasket.
"""
import sys, os, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction
from omega import (QC, ID, H, act, atom_through, disk_of, is_atom, in_omega,
                   inversive, factor, is_schmidt)


def parse_point(text, use_float=False):
    """'31/101 + 7/53 i', '0.3+0.2i', 'i/2' ..."""
    t = text.replace(" ", "")
    m = re.fullmatch(r"([+-]?[^+-]*)([+-][^+-]*i)", t) or re.fullmatch(r"()([+-]?[^+-]*i)", t)
    if not m:
        raise SystemExit(f"cannot parse point {text!r}; use e.g. '1/3+4/7i' or '0.3+0.2i'")
    re_part, im_part = m.group(1) or "0", m.group(2)[:-1]
    if im_part in ("", "+"): im_part = "1"
    if im_part == "-": im_part = "-1"
    if use_float:
        return complex(float(Fraction(re_part)), float(Fraction(im_part)))
    return QC(Fraction(re_part), Fraction(im_part))


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 1
    use_float = "--float" in argv
    argv = [a for a in argv if a != "--float"]
    z = parse_point(argv[1], use_float)
    n = int(argv[2]) if len(argv) > 2 else 10
    nmax = int(argv[3]) if len(argv) > 3 else 4000
    if z.imag <= 0:
        raise SystemExit("z must lie in the upper half plane")

    print(f"z = {z}    (Im z = {float(z.imag):.6f})\n")
    print(f"{'k':>3}  {'curvature':>10}  {'radius':>12}  centre / half plane")
    print("  " + "-" * 78)
    Y, w, disks, atoms = ID, z, [], []
    for k in range(1, n + 1):
        X = atom_through(w, nmax)
        if X is None:
            print(f"\n  no atom of curvature <= {2*nmax} contains the reduced point "
                  f"{w} -- it lies in (or very near) the residual set of the gasket; stopping.")
            break
        assert is_atom(X) and in_omega(X)
        Y = Y * X
        D = disk_of(Y)
        assert D.contains(z), "internal error: z left the disk"
        disks.append(D); atoms.append(X)
        if D.A == 0:
            print(f"{k:>3}  {0:>10}  {'infinite':>12}  {{Im z > {D.height()}}}")
        else:
            c = D.centre()
            print(f"{k:>3}  {-D.A:>10}  {'1/'+str(-D.A):>12}  "
                  f"{float(c.real):+.12f} {float(c.imag):+.12f} i")
        w = act(X.inv(), w)

    print()
    for k, X in enumerate(atoms, 1):
        print(f"  X_{k} = {X}")
    print(f"\n  product X_1...X_{len(atoms)} = {Y}")
    if disks:
        print(f"  its disk is D_{len(disks)}; z inside all D_k: "
              f"{all(D.contains(z) for D in disks)}; "
              f"strictly nested: {all(inversive(disks[j+1], disks[j]) >= 1 and disks[j+1] != disks[j] for j in range(len(disks)-1))}; "
              f"all in the Schmidt arrangement: {all(is_schmidt(D) for D in disks)}")
        a, u = factor(Y)
        print(f"  re-factoring the product returns {len(a)} atoms "
              f"(unique factorisation): {len(a) == len(atoms)}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
