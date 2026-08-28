# Euclidean counting in the Gaussian Schmidt arrangement

This document builds on the classification in [circle-classification.md](circle-classification.md): a circle of Euclidean curvature \(2n\) (\(n \ge 1\)) with center \(z_0\) lies in \(\mathcal{S} = \mathrm{PSL}_2(\mathbb{Z}[i]) \cdot \hat{\mathbb{R}}\) iff
$$
\zeta := 2n z_0 = x + iy \in \mathbb{Z}[i],
\qquad x \equiv 0,\ y \equiv 1 \ (\mathrm{mod}\ 2),
\qquad x^2 + y^2 \equiv 1 \pmod{4n}.
$$

## 1. The counting function

Since \(\mathcal{S}\) is invariant under translation by \(\mathbb{Z}[i]\), the natural density of curvature-\(2n\) circles is the number of them per unit area, i.e. with center in the fundamental square:
$$
N_e(n) \;:=\; \#\{\,\mathcal{C} \in \mathcal{S} \;:\; \text{curvature}(\mathcal{C}) = 2n,\ \ \text{center}(\mathcal{C}) \in [0,1) + [0,1)i \,\}.
$$

**Reduction to a congruence count.** Centers in \([0,1)^2\) correspond bijectively to residues \(\zeta \bmod 2n\) (in each coordinate). The defining conditions descend to \((\mathbb{Z}/2n)^2\):

- \((x + 2n)^2 \equiv x^2 \pmod{4n}\), so the congruence \(x^2+y^2 \equiv 1 \pmod{4n}\) only depends on \(x, y \bmod 2n\);
- \(2n\) is even, so the parities of \(x, y\) are also well defined mod \(2n\).

Hence
$$
N_e(n) \;=\; \#\bigl\{ (x,y) \in (\mathbb{Z}/2n\mathbb{Z})^2 \;:\; x \text{ even},\ y \text{ odd},\ x^2+y^2 \equiv 1 \ (\mathrm{mod}\ 4n) \bigr\}.
$$
(Each residue class is exactly one circle: center and radius determine the circle, and the co-curvature is determined by \((n,\zeta)\).)

## 2. Closed form

Let
$$
r(q) := \#\{(x,y) \bmod q : x^2 + y^2 \equiv 1 \pmod q\}.
$$

**Step 1: \(N_e(n) = \tfrac18 r(4n)\).** Lift the count to \((\mathbb{Z}/4n)^2\): each class mod \(2n\) has \(4\) lifts, and the congruence is invariant under lifting, so the number of solutions mod \(4n\) *without* parity constraints is \(4 \cdot \#\{\text{solutions mod } 2n\}\). Moreover any solution of \(x^2 + y^2 \equiv 1 \pmod 4\) has exactly one of \(x, y\) odd, and the swap \((x,y) \mapsto (y,x)\) is a parity-reversing involution on solutions; so the parity condition (\(x\) even, \(y\) odd) cuts the count exactly in half:
$$
N_e(n) = \frac12 \cdot \frac14\, r(4n) = \frac{r(4n)}{8}.
$$

**Step 2: local factors.** \(r\) is multiplicative in \(q\) by CRT. Write \(n = 2^e n'\) with \(n'\) odd, so \(4n = 2^{e+2} n'\).

*Odd \(p\), \(p^a \| n'\):* the conic \(x^2 + y^2 = 1\) over \(\mathbb{F}_p\) is smooth with
$$
r(p) = p - \chi_{-4}(p), \qquad \chi_{-4}(p) = (-1)^{(p-1)/2}
$$
points (the standard count: sweeping lines through the rational point \((1,0)\) parametrizes all but the points at infinity, of which there are \(1 + \chi_{-4}(p)\)). Smoothness (the gradient \((2x, 2y)\) never vanishes on the conic mod odd \(p\)) gives, by Hensel lifting,
$$
r(p^a) = p^{a-1}\bigl(p - \chi_{-4}(p)\bigr).
$$

*\(p = 2\):* for \(k \ge 2\) one checks (directly for \(k = 2, 3\); by the ramified analogue of Hensel from there, verified numerically for \(k \le 11\))
$$
r(2^k) = 2^{k+1}.
$$

**Step 3: assemble.**
$$
N_e(n) = \frac18 \cdot 2^{e+3} \prod_{p^a \| n'} p^{a-1}(p - \chi_{-4}(p))
= 2^e \prod_{p^a \| n'} p^a \Bigl(1 - \frac{\chi_{-4}(p)}{p}\Bigr),
$$
that is:

> **Theorem (exact count).**
> $$
> \boxed{\;N_e(n) \;=\; n \prod_{\substack{p \mid n \\ p \text{ odd}}} \Bigl(1 - \frac{\chi_{-4}(p)}{p}\Bigr)\;}
> \qquad
> \chi_{-4}(p) = \begin{cases} +1 & p \equiv 1 \ (4) \\ -1 & p \equiv 3 \ (4). \end{cases}
> $$
> In particular \(N_e\) is multiplicative, with \(N_e(2^a) = 2^a\), \(N_e(p^a) = p^{a-1}(p - \chi_{-4}(p))\).

**Small values** (each verified against a brute-force enumeration of the group orbit, exact set equality of residue classes for \(n \le 15\); see [scripts/verify_classification.py](scripts/verify_classification.py)):

| \(n\) | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| \(N_e(n)\) | 1 | 2 | 4 | 4 | 4 | 8 | 8 | 8 | 12 | 8 | 12 | 16 | 12 | 16 | 16 |

So per unit area there is \(1\) circle of curvature \(2\) (the Ford circle), \(2\) of curvature \(4\), \(4\) of curvature \(6\), and so on.

## 3. Asymptotics: enter Catalan

### Dirichlet series

Since \(N_e\) is multiplicative, its Dirichlet series has an Euler product. For odd \(p\):
$$
\sum_{a \ge 0} \frac{N_e(p^a)}{p^{as}}
= 1 + \Bigl(1 - \tfrac{\chi_{-4}(p)}{p}\Bigr) \sum_{a\ge1} p^{a(1-s)}
= 1 + \Bigl(1 - \tfrac{\chi_{-4}(p)}{p}\Bigr) \frac{p^{1-s}}{1 - p^{1-s}}
= \frac{1 - \chi_{-4}(p)\, p^{-s}}{1 - p^{1-s}},
$$
and at \(p = 2\) the factor is \((1 - 2^{1-s})^{-1}\). Since \(\chi_{-4}(2) = 0\),
$$
\sum_{n \ge 1} \frac{N_e(n)}{n^s} \;=\; \frac{\zeta(s-1)}{L(s, \chi_{-4})}.
$$

The pole at \(s = 2\) has residue \(1/L(2, \chi_{-4})\) — and
$$
L(2, \chi_{-4}) \;=\; 1 - \frac{1}{3^2} + \frac{1}{5^2} - \frac{1}{7^2} + \cdots \;=\; G,
$$
**Catalan's constant**, \(G = 0.9159655941\ldots\)

### Elementary derivation of the partial sums

No Tauberian theorem is needed. Expanding the Euler product of \(1/L(s,\chi_{-4})\) termwise is the Möbius identity
$$
N_e(n) = \sum_{d \mid n} \mu(d)\, \chi_{-4}(d)\, \frac{n}{d}
$$
(both sides are multiplicative; check on prime powers, using \(\chi_{-4}(2) = 0\) to kill even \(d\)). Then
$$
\sum_{n \le X} N_e(n)
= \sum_{d \le X} \mu(d)\chi_{-4}(d) \sum_{q \le X/d} q
= \frac{X^2}{2} \sum_{d=1}^{\infty} \frac{\mu(d)\chi_{-4}(d)}{d^2} + O\!\Bigl(X \log X\Bigr),
$$
and \(\sum_d \mu(d)\chi_{-4}(d) d^{-2} = \prod_p (1 - \chi_{-4}(p) p^{-2}) = L(2,\chi_{-4})^{-1}\). Hence:

> **Theorem (mean count).**
> $$
> \sum_{n \le X} N_e(n) \;=\; \frac{X^2}{2G} \;+\; O(X \log X),
> \qquad G = L(2, \chi_{-4}) = \text{Catalan's constant}.
> $$

Numerically (sieve up to \(10^6\)): \(\dfrac{2G}{X^2}\sum_{n\le X} N_e(n) = 1.000782,\ 1.000089,\ 1.000006,\ 1.000001\) at \(X = 10^3, 10^4, 10^5, 10^6\) — the convergence is visibly consistent with the stated error term.

### Reformulations

- **By curvature.** The number of circles of \(\mathcal{S}\) per unit area with curvature \(\le \kappa\) is
  $$
  \sim \frac{\kappa^2}{8G} \;\approx\; 0.13646\, \kappa^2 .
  $$
- **By radius.** The number per unit area with radius \(\ge \rho\) is \(\sim \dfrac{1}{8 G \rho^2}\).
- **On average**, there are \(\sim \dfrac{2n}{G}\) circles of curvature exactly \(2n\) per unit area — i.e. the arithmetic factor \(\prod_{p \mid n}(1 - \chi_{-4}(p)/p)\) has mean value \(1/G \approx 1.0917\). (Its own Dirichlet series is \(\zeta(s)/L(s+1,\chi_{-4})\), giving \(\sum_{n \le X} \prod_{p\mid n, \text{odd}}(1 - \chi_{-4}(p)/p) \sim X/G\).)

### Fluctuations around the mean

\(N_e(n)/n = \prod_{p \mid n \text{ odd}}(1 - \chi_{-4}(p)/p)\) oscillates: it exceeds \(1\) when \(n\) favors primes \(\equiv 3 \pmod 4\) and dips below when \(n\) favors \(p \equiv 1 \pmod 4\). By Mertens' theorem in arithmetic progressions (\(\sum_{p \le Y,\, p \equiv a (4)} 1/p = \tfrac12 \log\log Y + O(1)\)), the extreme orders are
$$
\limsup_{n} \frac{N_e(n)}{n (\log\log n)^{1/2}} \in (0,\infty),
\qquad
\liminf_{n} \frac{N_e(n) (\log\log n)^{1/2}}{n} \in (0,\infty),
$$
so the normalized count fluctuates within \((\log\log n)^{\pm 1/2}\) — mild, as expected for a multiplicative weight.

### Geometric corollaries

- **Total circumference per curvature level.** The circles of curvature \(2n\) in the unit square have total length \(\frac{2\pi}{2n} N_e(n) = \pi \prod_{p \mid n \text{ odd}}(1 - \chi_{-4}(p)/p)\), which is \(\asymp \pi\) for every \(n\) and \(\to \pi/G \approx 3.4298\) on average. Summed over \(n \le X\) the total length of the arrangement grows linearly, \(\sim (\pi/G) X\).
- **Total disk area per unit area** up to curvature \(2X\):
  $$
  \sum_{n \le X} N_e(n)\, \frac{\pi}{4n^2} \;\sim\; \frac{\pi}{4G} \log X \;\to\; \infty,
  $$
  a logarithmic divergence quantifying the (well-known) density of the arrangement: the disks bounded by Schmidt circles cover every point infinitely often, at logarithmic rate in the curvature cutoff.

## 4. Remarks

1. **Where Catalan comes from.** Here \(G\) enters purely arithmetically, as the value \(L(2, \chi_{-4})\) produced by the local solution densities of \(x^2 + y^2 \equiv 1 \pmod{4n}\) — equivalently through \(\zeta_{\mathbb{Q}(i)}(2) = \zeta(2) L(2,\chi_{-4})\). (The same constant is famously the hyperbolic-volume quantum of the Bianchi orbifold \(\mathrm{PSL}_2(\mathbb{Z}[i]) \backslash \mathbb{H}^3\); no hyperbolic geometry was needed above, but the agreement is of course not a coincidence and is worth revisiting when we study intrinsic/hyperbolic counting.)
2. **\(\mathrm{PGL}_2\) normalization.** For the arrangement \(\mathcal{S} \sqcup i\mathcal{S}\) of the full unit-determinant group, every count doubles: \(N_e^{\mathrm{PGL}}(n) = 2n \prod (1 - \chi_{-4}(p)/p)\) and \(\sum_{n \le X} N_e^{\mathrm{PGL}}(n) \sim X^2/G\). Catalan's constant appears either way.
3. **Contrast with Apollonian counting.** This is *Euclidean* counting: all circles of a given curvature, everywhere, weighted per unit area — a completely explicit multiplicative function. Counting circles inside a *fixed bounded region* that is not translation-full (or inside a single tangency component, as in Apollonian circle counting) is a genuinely harder equidistribution problem; the exact formula for \(N_e\) is the natural baseline for those questions.

## 5. Verification

All statements above with finite content are machine-checked in [scripts/verify_classification.py](scripts/verify_classification.py):
orbit enumeration \(=\) congruence classes (sets, \(n \le 15\)), \(N_e\) congruence count \(=\) closed form (\(n \le 200\)), \(r(2^k) = 2^{k+1}\) (\(k \le 11\)), and the \(X^2/(2G)\) asymptotic (up to \(X = 10^6\)).
