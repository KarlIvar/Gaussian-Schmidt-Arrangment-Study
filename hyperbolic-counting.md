# Schmidt circles in the ideal triangle: the invariant \(\alpha\), an algorithm, and Hurwitz class numbers

This document treats the circles of the Gaussian Schmidt arrangement as objects of the *hyperbolic* plane \(\mathbb{H}^2 = \{\operatorname{Im} z > 0\}\). It defines the isometry invariant \(\alpha\), describes the algorithm implemented in [scripts/alpha_circles.py](scripts/alpha_circles.py) for listing all circles with \(\alpha(\omega) = n\) inside the ideal triangle \(T\) with vertices \(0, 1, \infty\), and proves the counting formula
$$
\#\{\text{circles } \omega \text{ in } T \text{ with } \alpha(\omega) = n\}_{\text{weighted}} \;=\; 3\,H(n^2 - 1),
$$
where \(H\) is the Hurwitz class number and circles whose (hyperbolic) center lies on an edge of \(T\) count half.

Throughout, [circle-classification.md](circle-classification.md) is used: a circle of curvature \(2q\) with center \(\tfrac{x + yi}{2q}\) lies in \(\mathcal{S} = \mathrm{PSL}_2(\mathbb{Z}[i]) \cdot \hat{\mathbb{R}}\) iff \(x\) is even, \(y\) is odd and \(x^2 + y^2 \equiv 1 \pmod{4q}\) — and in the rotated companion \(i\mathcal{S}\) (the rest of the \(\mathrm{PGL}_2(\mathbb{Z}[i])\)-arrangement) iff the same congruence holds with \(x\) odd, \(y\) even.

## 1. The invariant \(\alpha\)

For a Euclidean circle \(\omega \subset \mathbb{H}^2\) with center \(x + yi\) and radius \(r\) (so \(y > r\)), define
$$
\alpha(\omega) \;=\; \frac{y}{r} \;=\; y \cdot \operatorname{curv}(\omega).
$$

**\(\alpha\) is invariant under the isometries \(\mathrm{PSL}_2(\mathbb{R})\) of \(\mathbb{H}^2\).** Indeed, a Euclidean circle inside \(\mathbb{H}^2\) is the same thing as a *hyperbolic* circle: the hyperbolic circle with hyperbolic center \(a + bi\) and hyperbolic radius \(\rho\) meets the vertical geodesic through its center at \(a + ibe^{\pm\rho}\), and is the Euclidean circle with
$$
\text{center } a + ib\cosh\rho, \qquad \text{radius } b\sinh\rho .
$$
Therefore
$$
\alpha(\omega) = \frac{b \cosh \rho}{b \sinh\rho} = \coth\bigl(\rho_{\mathrm{hyp}}(\omega)\bigr),
$$
a function of the hyperbolic radius alone — manifestly invariant. (Conversely \(\rho = \operatorname{arcoth}\alpha = \tfrac12\log\tfrac{\alpha+1}{\alpha-1}\).) Since the Schmidt arrangement is invariant under \(\mathrm{PSL}_2(\mathbb{Z}) = \mathrm{PSL}_2(\mathbb{Z}[i]) \cap \mathrm{PSL}_2(\mathbb{R})\), the isometries in that subgroup permute Schmidt circles with a fixed value of \(\alpha\). (So do the reflections \(z \mapsto -\bar z\), \(z \mapsto 1 - \bar z\), etc., which also preserve each family \(\mathcal{S}\), \(i\mathcal{S}\).)

**On Schmidt circles \(\alpha\) is a positive integer.** For a circle of curvature \(2q\) and center \(\tfrac{x+yi}{2q}\),
$$
\alpha = \frac{y/(2q)}{1/(2q)} = y ,
$$
the integer numerator of the height. Circles contained in \(\mathbb{H}^2\) have \(\alpha > 1\); the value \(\alpha = 1\) corresponds to circles *tangent* to \(\hat{\mathbb{R}}\) (horocycle-like, e.g. the Ford circles), which we exclude. So the possible values are the integers \(n \ge 2\), and every one occurs:

- **\(n\) odd** — circles of \(\mathcal{S}\) itself (\(x\) even);
- **\(n\) even** — circles of the companion \(i\mathcal{S}\) (\(x\) odd).

The algorithm below treats both uniformly: the congruence \(x^2 + n^2 \equiv 1 \pmod 4\) *forces* the correct parity of \(x\), so no case split is ever needed. (If one insists on the strict \(\mathrm{PSL}_2\) convention, the counts below apply verbatim to odd \(n\) and are \(0\) for even \(n\); in the \(\mathrm{PGL}_2\) arrangement they apply to all \(n \ge 2\).)

## 2. Circles with \(\alpha = n\) are quadratic forms; centers are CM points

Fix \(n \ge 2\) and write \(N = n^2 - 1\). By the classification, circles with \(\alpha(\omega) = n\) are exactly the
$$
\omega_{q,x}: \quad \text{radius } \frac{1}{2q}, \quad \text{center } \frac{x + ni}{2q},
\qquad
x^2 + n^2 - 1 = 4qm \ \ \text{for some } m \in \mathbb{Z},
$$
with \(q \ge 1\) (and \(2m\) = the co-curvature). Encode \(\omega_{q,x}\) as the integral binary quadratic form
$$
f_\omega = (q, -x, m), \qquad \operatorname{disc} f_\omega = x^2 - 4qm = 1 - n^2 = -N < 0,
$$
positive definite since \(q > 0\). This is a **bijection** between circles with \(\alpha = n\) and positive definite integral forms of discriminant \(1 - n^2\) (imprimitive forms included — primitivity is *not* automatic here, unlike the determinant \(-1\) Hermitian matrices), and it is **equivariant**: the \(\mathrm{PSL}_2(\mathbb{Z})\)-action on circles (restriction of the Möbius action, i.e. \(M \mapsto (g^{-1})^{\mathsf T} M g^{-1}\) on Hermitian matrices, whose real part is the form) matches \(\mathrm{SL}_2(\mathbb{Z})\)-equivalence of forms.

Two structural facts fall out:

1. **The hyperbolic center of \(\omega_{q,x}\) is the root of its form.** From §1, the hyperbolic center is \(a + bi\) with \(a = \tfrac{x}{2q}\) and \(b = \sqrt{(n/2q)^2 - (1/2q)^2} = \tfrac{\sqrt{N}}{2q}\), i.e.
$$
z_{\mathrm{hyp}}(\omega_{q,x}) \;=\; \frac{x + i\sqrt{n^2 - 1}}{2q}
\;=\; \text{root in } \mathbb{H}^2 \text{ of } q z^2 - x z + m ,
$$
a **CM point of discriminant \(1 - n^2\)**.
2. **The hyperbolic radius depends only on \(n\):** \(\rho = \operatorname{arcoth} n\).

> **Reformulation.** The circles of the (\(\mathrm{PGL}_2\)-)Schmidt arrangement contained in \(\mathbb{H}^2\) are precisely the hyperbolic circles of radius \(\operatorname{arcoth} n\) centered at the CM points of discriminant \(1 - n^2\), for \(n = 2, 3, 4, \dots\)

## 3. Membership in the ideal triangle

Let \(T\) be the closed ideal triangle with vertices \(0, 1, \infty\): the region \(0 \le \operatorname{Re} z \le 1\) outside the open disk \(|z - \tfrac12| < \tfrac12\). We locate a circle by its hyperbolic center (this is the right notion of "where the circle is"; the user-facing convention "circles on the edge count \(\tfrac12\)" refers to the hyperbolic center lying on an edge of \(T\)). With \(z_{\mathrm{hyp}} = \tfrac{x + i\sqrt N}{2q}\):

- \(0 \le \operatorname{Re} z_{\mathrm{hyp}} \le 1 \iff 0 \le x \le 2q\);
- \(\bigl|z_{\mathrm{hyp}} - \tfrac12\bigr|^2 \ge \tfrac14 \iff \left(\tfrac{x}{2q}\right)^2 - \tfrac{x}{2q} + \tfrac{N}{4q^2} \ge 0 \iff x(2q - x) \le N\).

So:
$$
\boxed{\;\omega_{q,x} \text{ has hyperbolic center in } T \iff 0 \le x \le 2q \ \text{ and } \ x(2q-x) \le n^2 - 1\;}
$$
with the center on an **edge** of \(T\) exactly when one of the three constraints is an equality: \(x = 0\) (left edge), \(x = 2q\) (right edge), \(x(2q - x) = N\) (bottom edge). Since \(x(2q-x) = N \iff 2qx = x^2 + N = 4qm \iff x = 2m\), the three edge conditions are simply
$$
x = 0, \qquad x = 2q, \qquad x = 2m .
$$
In form language, the triangle condition is \(0 \le -b \le 2\min(a, c)\) for \(f = (a,b,c)\) — a "triangle reduction" which (as §5 shows) selects each \(\mathrm{SL}_2(\mathbb{Z})\)-class exactly three times with the edge/automorphism weights, in contrast to Gauss reduction \(|b| \le a \le c\) which selects it once.

**Consequences of parity.** For even \(n\), \(x\) is odd, so no edge condition can ever hold (\(0\), \(2q\), \(2m\) are all even): *every* circle is interior and the weighted count is an honest count. For odd \(n\), edge circles are ubiquitous (for \(n = 3, 5\) *all* circles lie on edges).

## 4. The algorithm

Input \(n \ge 2\); write \(N = n^2 - 1\).

1. **Bound \(q\).** If \(0 < x < 2q\) then \(x(2q - x) \ge 2q - 1\), so \(q \le \tfrac{N+1}{2}\); if moreover \(x\) and \(2q - x\) are both even (the case \(n\) odd) then \(x(2q-x) \ge 2(2q-2)\), giving \(q \le \tfrac{N}{4} + 1\). If \(x \in \{0, 2q\}\), integrality of \(m\) forces \(4q \mid N\), so \(q \le \tfrac N4\). Hence
$$
q_{\max} = \tfrac{N}{4} + 1 \ (n \text{ odd}), \qquad q_{\max} = \tfrac{N+1}{2} \ (n \text{ even}).
$$
2. **Sweep the curvatures.** For each \(q = 1, \dots, q_{\max}\), solve the quadratic congruence
$$
x^2 \equiv 1 - n^2 \pmod{4q}
$$
by factoring \(4q\) (smallest-prime-factor sieve up to \(4q_{\max}\), computed once), extracting square roots modulo each prime power (Tonelli–Shanks + Hensel for odd \(p\); the standard \(2\)-adic case analysis at \(p = 2\); valuation handling when \(p \mid 1 - n^2\)), and combining by CRT.
3. **Filter and weight.** Keep the roots \(x \in [0, 2q]\) with \(x(2q - x) \le N\); set \(m = \tfrac{x^2 + N}{4q}\); assign weight \(\tfrac12\) if \(x \in \{0,\, 2q,\, 2m\}\), else \(1\). Output the circle: radius \(\tfrac1{2q}\), center \(\tfrac{x + ni}{2q}\), hyperbolic center \(\tfrac{x + i\sqrt N}{2q}\).

**Complexity.** The sieve and the sweep cost \(O(n^2 \operatorname{polylog} n)\) bit operations altogether (each modulus \(4q \le 2n^2 + 2\) has \(O(\log)\) prime-power factors; most \(q\) are rejected at the first non-residue). The output has size \(3H(n^2-1) = n^{1 + o(1)}\). In practice: \(n = 301\) in \(0.05\) s, \(n = 1000\) (3648 circles) in \(2\) s. An *output-sensitive* \(n^{1+\varepsilon}\) alternative exists — enumerate Gauss-reduced forms of discriminant \(1 - n^2\) (cost \(O(\sqrt N \cdot d)\)) and move each class into \(T\) by explicit coset representatives — but the congruence sweep is simpler, exact, and fast enough.

**Display.** The script draws the configuration in the upper half-plane (triangle \(0, 1, \infty\)) and, toggled by the key `m` or a button (or `--model disk`), in the Poincaré disk via the Möbius map
$$
\varphi(z) = \frac{z - \zeta}{z - \bar\zeta}, \qquad \zeta = \tfrac{1 + i\sqrt3}{2} \ (\text{the centroid of } T),
$$
which sends \(T\) to the ideal triangle with vertices at the cube roots of unity and the centroid to \(0\). Circles are mapped exactly (image center/radius in closed form, self-tested against sampled points). Since the rotation \(\rho: z \mapsto 1/(1-z)\) of order \(3\) lies in \(\mathrm{PSL}_2(\mathbb{Z})\), preserves \(T\), the arrangement, and \(\alpha\), the configuration is invariant under it — so the disk picture has exact threefold symmetry (in fact full \(D_3\) symmetry, including the reflections). See [figures/](figures/) for \(\alpha = 7\) and \(\alpha = 15\).

## 5. The count is \(3\,H(n^2 - 1)\)

Recall the **Hurwitz class number** \(H(N)\) (\(N \equiv 0, 3 \bmod 4\)): the number of \(\mathrm{SL}_2(\mathbb{Z})\)-classes of positive definite integral binary quadratic forms of discriminant \(-N\), *imprimitive forms included*, each class weighted by \(\tfrac{2}{|\mathrm{Aut}|}\) — i.e. weight \(\tfrac13\) for the classes of multiples of \(u^2 + uv + v^2\), weight \(\tfrac12\) for multiples of \(u^2 + v^2\), weight \(1\) otherwise. Equivalently: weight \(1/s\), where \(s \in \{1,2,3\}\) is the order of the stabilizer in \(\mathrm{PSL}_2(\mathbb{Z})\). Note \(N = n^2 - 1 \equiv 0 \bmod 4\) for odd \(n\) and \(\equiv 3 \bmod 4\) for even \(n\), so \(H(n^2-1)\) is always of the admissible type.

> **Theorem.** For every integer \(n \ge 2\), the number of Schmidt circles \(\omega\) with \(\alpha(\omega) = n\) and hyperbolic center in \(T\), counted with weight \(\tfrac12\) when the center lies on an edge of \(T\), equals \(3\,H(n^2 - 1)\).

*Proof.* By §2 the circles with \(\alpha = n\) form a single \(\mathrm{PSL}_2(\mathbb{Z})\)-set isomorphic (equivariantly) to the positive definite forms of discriminant \(1 - n^2\); classes correspond, and the stabilizer of a circle equals the stabilizer of its hyperbolic center \(z_0\) (an isometry fixing \(z_0\) and preserving the arrangement fixes the circle, and conversely). It therefore suffices to prove the

**Incidence lemma.** *Let \(z_0 \in \mathbb{H}^2\) with \(s = |\mathrm{Stab}_{\mathrm{PSL}_2(\mathbb{Z})}(z_0)|\). Then the orbit of \(z_0\), counted in \(T\) with weight \(\tfrac12\) on the edges, has total weight \(3/s\).*

The triangle \(T\) is a tile of the Farey tessellation, on which \(\Gamma = \mathrm{PSL}_2(\mathbb{Z})\) acts transitively with tile stabilizer of order \(3\) (generated by \(\rho: z \mapsto 1/(1-z)\)). Let \(t(z) \in \{1, 2\}\) be the number of closed tiles containing \(z\) (\(1\) interior, \(2\) on an edge; tessellation vertices are ideal, so never in \(\mathbb{H}^2\)); \(t\) is constant on \(\Gamma\)-orbits, and the weight of \(z \in \bar T\) is exactly \(1/t(z)\). Count the set \(\{\gamma \in \Gamma : \gamma z_0 \in \bar T\}\) in two ways: fibering over the point \(\gamma z_0\) gives \(s \cdot \#(\Gamma z_0 \cap \bar T)\); fibering over the tile \(\gamma^{-1} T \ni z_0\) gives \(3 \cdot t(z_0)\). Hence the weighted count is \(\#(\Gamma z_0 \cap \bar T)/t(z_0) = 3/s\). \(\square\)

Summing \(3/s\) over classes gives \(3 \sum_{\text{classes}} 1/s = 3H(n^2-1)\). \(\blacksquare\)

**Which stabilizers actually occur.**

- \(s = 2\) would require the class of \(f\cdot(1,0,1)\), i.e. \(n^2 - 4f^2 = 1\): impossible for \(n \ge 2\). So no circle is ever centered at an order-\(2\) elliptic point — the weights \(1\) and \(\tfrac12\) seen by the algorithm are the whole story, matching the theorem's bookkeeping (\(s=2\) never contributes).
- \(s = 3\) requires the class of \(f\cdot(1,1,1)\), i.e. the **Pell equation** \(n^2 - 3f^2 = 1\): solutions \(n = 2, 7, 26, 97, 362, \dots\) (\(n_{k+1} = 4n_k - n_{k-1}\)). For exactly these \(n\), exactly one class has \(s = 3\), and its unique representative in \(T\) is the circle centered at the **centroid** \(\zeta = \tfrac{1+i\sqrt3}{2}\) (weight \(1\), contributing \(1 = 3 \cdot \tfrac13\)). For \(n = 7\): the form \(4(u^2+uv+v^2) = (4,4,4)\), circle of curvature \(8\) centered at \(\tfrac{4 + 7i}{8}\).

**First examples.**

- \(n = 2\): \(3H(3) = 1\). The unique circle is the **inscribed circle of the ideal triangle** (center \(\tfrac12 + i\), radius \(\tfrac12\); hyperbolic center = centroid) — Schmidt's original "Farey circle". It lies in \(i\mathcal{S}\).
- \(n = 3\): \(3H(8) = 3\), realized as \(6\) edge circles of curvatures \(2, 4, 6\).
- \(n = 7\): \(3H(48) = 10 = 18 \cdot \tfrac12 + 1\) (\(18\) edge circles + the centroid circle).

**Data** (all verified by the script; "circles" is the raw count, "weighted" the edge-weighted one):

| \(n\) | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 | 21 | 22 | 23 | 24 | 25 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| circles | 1 | 6 | 6 | 12 | 6 | 19 | 15 | 27 | 9 | 24 | 30 | 24 | 12 | 48 | 36 | 45 | 12 | 48 | 48 | 48 | 12 | 54 | 63 | 66 |
| weighted | 1 | 3 | 6 | 6 | 6 | 10 | 15 | 18 | 9 | 12 | 30 | 12 | 12 | 36 | 36 | 27 | 12 | 30 | 48 | 36 | 12 | 36 | 63 | 48 |
| \(H(n^2{-}1)\) | \(\tfrac13\) | 1 | 2 | 2 | 2 | \(\tfrac{10}3\) | 5 | 6 | 3 | 4 | 10 | 4 | 4 | 12 | 12 | 9 | 4 | 10 | 16 | 12 | 4 | 12 | 21 | 16 |

By the class number formula, \(H(n^2-1)\) has average order a constant times \(n\) (and \(n^{1+o(1)}\) pointwise), so the triangle holds roughly \(\asymp n\) circles of each \(\alpha\)-level.

## 6. Verification

`python3 scripts/alpha_circles.py --selftest` checks:

1. the modular square-root machinery against brute force (400 random moduli);
2. the Hurwitz class number routine against a hand-checked table (including \(H(32) = 3\), confirmed independently by the Kronecker–Hurwitz relation \(\sum_{t} H(32 - t^2) = \sum_{d \mid 8} \max(d, 8/d) = 24\));
3. **weighted count \(= 3H(n^2-1)\) for all \(2 \le n \le 40\)** (and spot-checked at \(n = 301\), \(1000\));
4. closure of the computed configuration under the order-\(3\) rotation \(\rho\) of \(T\) (exact integer arithmetic on the Hermitian matrices) and under the mirror \(x \mapsto 2q - x\);
5. the Möbius circle-image formula used for the disk view, against sampled points.

Usage: `python3 scripts/alpha_circles.py 7` (interactive; press `m` to switch half-plane \(\leftrightarrow\) disk), `--list` for the table of circles, `--save DIR --no-show` for PNGs of both models.

## 7. Outlook

The identity gives Hurwitz class numbers a direct circle-geometric meaning inside the Schmidt arrangement. Two directions worth pursuing later:

- **Class number relations as circle geometry.** The Kronecker–Hurwitz relations \(\sum_{t \in \mathbb{Z}} H(4M - t^2) = \sum_{d \mid M} \max(d, M/d)\) (and the Eichler–Selberg trace formula behind them) now assert identities between counts of Schmidt circles of *different* hyperbolic radii \(\operatorname{arcoth} n\) in the same triangle. A geometric proof via the arrangement — e.g. tangency or incidence arguments between the \(\alpha\)-levels — could be a genuinely new way to see these relations.
- **Level structure.** Replacing \(T\) by other tiles/regions, or intersecting with congruence conditions on \(q\), should produce class numbers for non-maximal orders and Hurwitz-type counts with level, again with explicit circle configurations attached.
