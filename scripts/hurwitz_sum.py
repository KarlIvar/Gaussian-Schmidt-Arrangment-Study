"""Numerical experiment: growth of S(X) = sum_{n<=X} 3*H(n^2-1).

H(N) = Hurwitz class number, computed by enumerating reduced positive
definite forms (a,b,c), b^2-4ac = -N, -a < b <= a <= c, b >= 0 if a == c,
with weights 1/2 for k*(1,0,1) and 1/3 for k*(1,1,1).

Also: volume of the Picard orbifold = G/3 (Humbert), G = Catalan.
We fit S(X) ~ c3 * X^2 and compare c3 against candidate constants.
"""
import math
from math import isqrt

def hurwitz(N):
    # N > 0, N ≡ 0 or 3 (mod 4); returns float (rational with denom 6)
    if N % 4 not in (0, 3):
        return 0.0
    total = 0.0
    amax = isqrt(N // 3) + 1
    for a in range(1, amax + 1):
        # b ≡ N mod 2 (since b^2 ≡ -N ≡ b^2 mod 4 -> b parity fixed by N mod 4)
        # b^2 + N ≡ 0 mod 4 requires b even if N ≡ 0 (4), b odd if N ≡ 3 (4)
        bstart = 0 if N % 4 == 0 else 1
        for b in range(bstart, a + 1, 2):
            t = b * b + N
            if t % (4 * a):
                continue
            c = t // (4 * a)
            if c < a:
                continue
            # reduced domain: -a < b <= a <= c ; if a == c then b >= 0
            # enumerate b >= 0 and count b and -b separately when both reduced
            if b == 0:
                mult = 1
            elif b == a:
                mult = 1          # (a,-a,c) ~ (a,a,c): only b=+a counted
            elif a == c:
                mult = 1          # b>=0 required when a==c
            else:
                mult = 2          # +-b both reduced and inequivalent
            # weights
            if a == b == c:
                w = 1.0 / 3.0
            elif b == 0 and a == c:
                w = 0.5
            else:
                w = 1.0
            total += mult * w
    return total

def main(X=800):
    checkpoints = [100, 200, 400, 600, 800, 1000, 1200]
    checkpoints = [x for x in checkpoints if x <= X] + ([X] if X not in checkpoints else [])
    S = 0.0
    results = {}
    for n in range(2, X + 1):
        S += 3.0 * hurwitz(n * n - 1)
        if n in checkpoints:
            results[n] = S
    G = 0.915965594177219015054603514932384110774  # Catalan
    volM = G / 3.0
    areaT = math.pi / 3.0   # area of modular surface PSL2(Z)\H^2
    print("vol(Picard orbifold) = G/3 =", volM)
    print()
    print(" X      S(X)=sum 3H(n^2-1)    S/X^2     S/(X^2) Richardson")
    prev = None
    for x in sorted(results):
        c = results[x] / x**2
        rich = ""
        if prev is not None:
            xp, cp = prev
            # if S = cX^2 + bX, then c(X) = c + b/X; Richardson: (X*c(X) - Xp*c(Xp))/(X-Xp)
            rich = (x * c - xp * cp) / (x - xp)
            rich = f"{rich:.6f}"
        print(f"{x:5d}  {results[x]:16.2f}  {c:.6f}   {rich}")
        prev = (x, results[x] / x**2)
    cfit = prev[1]
    print()
    print("candidate ratios (using last checkpoint c3 ~= S/X^2):")
    cands = {
        "pi^2/(3G) (=(area)^2/vol /pi^0...)": math.pi**2 / (3 * G),
        "pi/(3G)": math.pi / (3 * G),
        "1/(2G) (Euclid count const)": 1 / (2 * G),
        "pi/6": math.pi / 6,
        "pi^2/16": math.pi**2 / 16,
        "3/pi": 3 / math.pi,
        "2/pi": 2 / math.pi,
        "areaT^2/volM = (pi/3)^2/(G/3)": areaT**2 / volM,
        "areaT^2/(2 volM)": areaT**2 / (2 * volM),
        "areaT^2/(4 volM)": areaT**2 / (4 * volM),
        "areaT^2/(8 volM)": areaT**2 / (8 * volM),
        "pi^2/(18 zeta(3))": math.pi**2 / (18 * 1.2020569031595942854),
        "7 zeta(3)/(2 pi^2) *3": 3*7*1.2020569031595942854/(2*math.pi**2),
    }
    for name, v in cands.items():
        print(f"  c3 / [{name} = {v:.6f}] = {cfit / v:.4f}")

if __name__ == "__main__":
    import sys
    X = int(sys.argv[1]) if len(sys.argv) > 1 else 800
    main(X)
