# The circles of the Gaussian Schmidt arrangement: a Diophantine classification

## Setup and conventions

Let \(\Gamma = \mathrm{PSL}_2(\mathbb{Z}[i])\) act on \(\hat{\mathbb{C}} = \mathbb{C} \cup \{\infty\}\) by Möbius transformations
$$
\begin{pmatrix} a & b \\ c & d \end{pmatrix} \cdot z = \frac{az+b}{cz+d},
\qquad a,b,c,d \in \mathbb{Z}[i],\ ad-bc = 1 .
$$
The **Gaussian Schmidt arrangement** is the orbit of the extended real line,
$$
\mathcal{S} \;=\; \Gamma \cdot \hat{\mathbb{R}} .
$$
Each element of \(\mathcal{S}\) is a circle or a line in \(\mathbb{C}\). The goal of this document is a complete Diophantine description of which circles occur — their radii and centers — proved in both directions (necessity and sufficiency) and verified by machine enumeration of the orbit ([scripts/verify_classification.py](scripts/verify_classification.py)).

> **Convention warning.** We take the group to be \(\mathrm{PSL}_2(\mathbb{Z}[i])\), i.e. determinant \(1\) up to sign. The larger group \(\mathrm{PGL}_2(\mathbb{Z}[i])\) (all unit determinants) gives a strictly larger arrangement — see the Remark at the end. The distinction is real: for example, the imaginary axis is **not** in \(\mathcal{S}\), while the real axis and all its integer translates \(\operatorname{Im} z = k\) are.

## Circles as Hermitian matrices

A circle or line in \(\mathbb{C}\) is the zero set of
$$
A\,z\bar z + \bar B z + B \bar z + C = 0,
\qquad A, C \in \mathbb{R},\; B \in \mathbb{C},\; |B|^2 - AC > 0,
$$
which we encode by the Hermitian matrix
$$
M = \begin{pmatrix} A & B \\ \bar B & C \end{pmatrix},
\qquad
z \in \text{circle} \iff (\bar z, 1)\, M \begin{pmatrix} z \\ 1\end{pmatrix} = 0 .
$$
Completing the square, for \(A \neq 0\) this is the circle with
$$
\text{center } z_0 = -\frac{B}{A}, \qquad \text{radius } r = \frac{\sqrt{|B|^2 - AC}}{|A|} = \frac{\sqrt{-\det M}}{|A|} .
$$
For \(A = 0\) it is a line. \(M\) determines the circle, and conversely the circle determines \(M\) up to a real scalar; oriented circles correspond to \(M\) up to a *positive* scalar, and \(M \mapsto -M\) reverses orientation.

The extended real line \(\operatorname{Im} z = 0\) is \(-iz + i\bar z = 0\), i.e.
$$
M_0 = \begin{pmatrix} 0 & i \\ -i & 0 \end{pmatrix}, \qquad \det M_0 = -1 .
$$

**Equivariance.** If \(w = gz\) then \((z, 1)^{\mathsf T} \propto g^{-1}(w,1)^{\mathsf T}\), so the image circle \(g \cdot \mathcal{C}\) has matrix
$$
M_{g\mathcal C} = (g^{-1})^{\dagger}\, M_{\mathcal C}\, g^{-1} .
$$
For \(g \in \mathrm{SL}_2\) this preserves the determinant. Hence:

> \(\mathcal{S}\) is in bijection with the \(\Gamma\)-orbit of \(M_0\) inside Hermitian matrices of determinant \(-1\) (up to sign, for unoriented circles). Every \(M\) in the orbit has \(\det M = -1\), so every circle of \(\mathcal{S}\) has radius \(1/|A|\).

Since \(z \mapsto -z\) (the class of \(\operatorname{diag}(i,-i)\)) lies in \(\Gamma\) and reverses the orientation of \(\hat{\mathbb{R}}\) (a direct computation gives \(\operatorname{diag}(i,-i)^\dagger M_0 \operatorname{diag}(i,-i)^{-1}\)-conjugation \(= -M_0\)), each circle of \(\mathcal{S}\) occurs with both orientations; we may therefore always normalize \(A \geq 0\).

## The image of \(\hat{\mathbb{R}}\): explicit formulas

Let \(g = \begin{pmatrix} a & b \\ c & d\end{pmatrix} \in \mathrm{SL}_2(\mathbb{Z}[i])\), so \(g^{-1} = \begin{pmatrix} d & -b \\ -c & a \end{pmatrix}\). A two-line computation gives
$$
(g^{-1})^{\dagger} M_0\, g^{-1}
= \begin{pmatrix}
2\operatorname{Im}(c\bar d) & i(a\bar d - b\bar c) \\
\overline{i(a\bar d - b\bar c)} & 2\operatorname{Im}(a\bar b)
\end{pmatrix} .
$$
(These formulas were cross-checked numerically against circumcircles of image points for hundreds of random group elements.) So the image \(g\cdot\hat{\mathbb{R}}\) has

- **curvature** \(\;A = 2\operatorname{Im}(c\bar d)\) — always an *even rational integer*;
- **"co-curvature"** \(\;C = 2\operatorname{Im}(a\bar b)\) — also even;
- \(B = i(a\bar d - b\bar c) \in \mathbb{Z}[i]\).

Write the curvature as \(2n\) (\(n \geq 1\); the case \(n = 0\) gives the lines, treated below) and define the **curvature–center**
$$
\zeta \;:=\; 2n\,z_0 \;\in\; \mathbb{Z}[i], \qquad \zeta = x + iy ,
$$
so the circle has radius \(\tfrac{1}{2n}\) and center \(\tfrac{\zeta}{2n}\). With the normalization \(A = 2n > 0\) we have \(\zeta = -B\).

## Necessity: three constraints

**(i) Determinant.** \(\det M = -1\) reads \(4nm - |\zeta|^2 = -1\) where \(C = 2m\), i.e.
$$
x^2 + y^2 \;=\; 1 + 4nm \quad\text{for some } m \in \mathbb{Z},
\qquad\text{equivalently}\qquad
x^2 + y^2 \equiv 1 \pmod{4n}.
$$
Note the auxiliary \(m\) is not free data: \(2m\) is the co-curvature, determined by \((n, \zeta)\).

**(ii) Parity (the mod-2 invariant).** Reduce modulo \(2\) in \(\mathbb{Z}[i]\). Since \(\bar\alpha \equiv \alpha \pmod 2\) for every \(\alpha \in \mathbb{Z}[i]\) (as \(\alpha - \bar\alpha = 2i\operatorname{Im}\alpha\)),
$$
a\bar d - b \bar c \;\equiv\; ad - bc \;=\; 1 \pmod 2,
\qquad\text{so}\qquad
B = i(a\bar d - b\bar c) \equiv i \pmod 2 .
$$
This is invariant under \(B \mapsto -B\) (because \(-i \equiv i \bmod 2\)), so it holds for both orientations, and \(\zeta = -B\) satisfies
$$
\boxed{\;\zeta \equiv i \pmod 2\;}
\qquad\text{i.e.}\qquad x \text{ even},\; y \text{ odd}.
$$
This single congruence is what distinguishes the \(\mathrm{PSL}_2\)-arrangement from its \(90°\) rotation: it is *not* symmetric under \(x \leftrightarrow y\).

**(iii) Primitivity** is automatic: any common divisor of \(2n\), \(\zeta\), \(2m\) divides \(|\zeta|^2 - 4nm = 1\).

**Lines.** If \(A = 2\operatorname{Im}(c\bar d) = 0\) with \(\gcd(c,d) = 1\), then \(d/c \in \mathbb{Q}\) forces \(c = qu\), \(d = pu\) with \(p, q \in \mathbb{Z}\) coprime and \(u\) a unit; then \(ad - bc = u(ap - bq) = 1\) gives \(ap - bq = \bar u\), hence \(B = i(a\bar d - b\bar c) = i\,\bar u(ap-bq) = i\,\bar u^2 = \pm i\). So the only lines in \(\mathcal{S}\) are
$$
\pm i z \mp i \bar z + C = 0 \iff \operatorname{Im} z = k, \quad k \in \mathbb{Z} \ \ (C = \mp 2k \text{ even}),
$$
the horizontal integer lines. **Vertical lines do not belong to \(\mathcal{S}\)** (they belong to \(i\mathcal{S}\), the \(\mathrm{PGL}_2\)-companion; see the Remark).

## Sufficiency: descent

**Claim.** Every integral Hermitian \(M = \begin{pmatrix} 2n & -\bar\zeta \\ -\zeta & 2m \end{pmatrix}\) with \(\det M = -1\), \(n \geq 1\) and \(\zeta \equiv i \pmod 2\) lies in the \(\Gamma\)-orbit of \(M_0\) (up to sign).

*Proof (reduction algorithm).* Two moves, both realized by elements of \(\Gamma\) and both preserving the constraint set (integrality, even diagonal, \(\det = -1\), off-diagonal \(\equiv i \bmod 2\)):

- **Translate** (\(z \mapsto z + \lambda\), \(\lambda \in \mathbb{Z}[i]\)): replaces \(\zeta \mapsto \zeta + 2n\lambda\) (curvature unchanged). The off-diagonal changes by \(2n\lambda\), even, so the parity class is preserved. Choose \(\lambda\) so that \(\operatorname{Re}\zeta, \operatorname{Im}\zeta \in (-n, n]\); then \(|\zeta|^2 \le 2n^2\).
- **Invert** (\(z \mapsto -1/z\), the element \(S\)): swaps \(A \leftrightarrow C\) and sends \(B \mapsto -\bar B\) (which preserves \(B \equiv i \bmod 2\)).

Start with curvature \(2n > 0\). After translating, the co-curvature satisfies
$$
|2m| = \frac{\bigl| |\zeta|^2 - 1 \bigr|}{2n} \le \frac{2n^2 - 1}{2n} < n < 2n .
$$
If \(m \neq 0\), invert (and negate \(M\) if needed to restore positive curvature): the curvature strictly drops from \(2n\) to \(|2m| < n\). A strictly decreasing sequence of positive even integers must terminate, i.e. we reach \(m = 0\). Then \(|\zeta|^2 = 1\) and \(\zeta \equiv i \pmod 2\) force \(\zeta = \pm i\): the circle of radius \(\tfrac{1}{2n'}\) tangent to \(\hat{\mathbb{R}}\) at \(0\) (for the current curvature \(2n'\)). One more inversion sends it to the line \(\operatorname{Im} z = \mp n'\), a translate of \(\hat{\mathbb{R}}\), which is \(\Gamma\)-equivalent to \(M_0\). Reversing the walk expresses \(M\) as an element of the orbit. \(\blacksquare\)

## Main theorem

> **Theorem (classification).** A circle of radius \(r\) centered at \(z_0\) belongs to the Gaussian Schmidt arrangement \(\mathcal{S} = \mathrm{PSL}_2(\mathbb{Z}[i])\cdot\hat{\mathbb{R}}\) **iff**
> $$
> r = \frac{1}{2n} \ \ (n \in \mathbb{Z}_{\ge 1}), \qquad \zeta := 2n\,z_0 = x + iy \in \mathbb{Z}[i],
> $$
> with
> $$
> x \equiv 0,\ y \equiv 1 \pmod 2, \qquad x^2 + y^2 \equiv 1 \pmod{4n} .
> $$
> The lines of \(\mathcal{S}\) are exactly \(\operatorname{Im} z = k\), \(k \in \mathbb{Z}\).

**Single Diophantine equation.** Substituting \(x = 2u\), \(y = 2v+1\) absorbs the parity conditions, and the congruence becomes (\(4u^2 + 4v^2 + 4v + 1 = 1 + 4nm\)):
$$
\boxed{\;u^2 + v(v+1) \;=\; nm\;}
$$
So the circles of \(\mathcal{S}\) are exactly
$$
\left\{\; \left|\,z - \frac{2u + (2v+1)i}{2n}\,\right| = \frac{1}{2n}
\;:\; u, v, n, m \in \mathbb{Z},\ n \ge 1,\ u^2 + v(v+1) = nm \;\right\},
$$
one circle per solution \((n, u \bmod{\text{—}}, v)\); distinct \((n, u, v)\) give distinct circles, and \(m = (u^2 + v^2 + v)/n\) is determined (its double, \(2m\), is the co-curvature of the circle).

## First examples

- \(n = 1\) (radius \(\tfrac12\)): condition \(x \equiv 0 \bmod 2\), \(y \equiv 1 \bmod 2\), \(x^2 + y^2 \equiv 1 \bmod 4\) — one class, \(\zeta \equiv i \pmod 2\). These are the **Ford circles** \(|z - (k + \tfrac i2)| = \tfrac12\), \(k \in \mathbb{Z}\), tangent to \(\hat{\mathbb{R}}\) at the integers (together with their reflections at height \(-\tfrac12\), which have \(\zeta = -i \equiv i\)). Notably the "sideways" circle \(|z - \tfrac12| = \tfrac12\) is **not** in \(\mathcal{S}\).
- \(n = 2\) (radius \(\tfrac14\)): classes \(\zeta \equiv i,\, 3i \pmod 4\) — circles at heights \(\tfrac14, \tfrac34\) above/below each Gaussian integer.
- \(n = 3\) (radius \(\tfrac16\)): four classes mod \(6\): \(\zeta \equiv i,\ 5i,\ 2+3i,\ 4+3i\).

## Remark: the \(\mathrm{PGL}_2\) variant

\(\mathrm{PGL}_2(\mathbb{Z}[i]) / \mathrm{PSL}_2(\mathbb{Z}[i]) \cong \mathbb{Z}/2\), generated by the class of \(\operatorname{diag}(i, 1)\), i.e. \(z \mapsto iz\). Repeating the parity computation with \(\det g = \pm i\) gives \(\zeta \equiv \pm 1 \pmod 2\) (\(x\) odd, \(y\) even). Hence
$$
\mathrm{PGL}_2(\mathbb{Z}[i]) \cdot \hat{\mathbb{R}} \;=\; \mathcal{S} \,\sqcup\, i\,\mathcal{S},
$$
a **disjoint** union: the two families share no circle and no line (\(i\mathcal S\) has the vertical lines and the \(x\)-odd-\(y\)-even circles). The full \(\mathrm{PGL}_2\)-arrangement is the one with \(4\)-fold rotational symmetry; \(\mathcal{S}\) itself is invariant only under \(z \mapsto z + \mathbb{Z}[i]\), \(z \mapsto -z\), \(z \mapsto \bar z\) (and all of \(\Gamma\), of course). In the \(\mathrm{PGL}_2\) arrangement the classification reads: \(x + y\) odd, \(x^2 + y^2 \equiv 1 \pmod{4n}\), and all counting functions double.

## Verification

[scripts/verify_classification.py](scripts/verify_classification.py) enumerates the actual orbit of \(M_0\) under \(\langle T, T_i, S \rangle\) by BFS over Hermitian matrices with bounded entries (finite state space, since \(\det = -1\) bounds \(|B|\) in terms of \(A, C\)) and confirms:

1. the curvature/center formulas against numerically fitted circumcircles of image points;
2. **exact set equality** between the orbit's residue classes \((\zeta \bmod 2n)\) and the congruence classes of the Theorem, for all \(n \le 15\) (stable across cutoffs \(|A|,|C| \le 60\) and \(\le 120\));
3. the parity invariant \(\zeta \equiv i \pmod 2\) across the entire computed orbit;
4. that the only lines produced are \(\operatorname{Im} z = k\).

See [euclidean-counting.md](euclidean-counting.md) for the counting function these classes define and its asymptotics (where Catalan's constant appears).
