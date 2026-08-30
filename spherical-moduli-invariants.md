# Spherical moduli invariants of Schmidt circles: integer-cotangent radii, Pell units, and the third slice \(d = 4(\ell^2+1)\)

This is the spherical companion of [moduli-invariants.md](moduli-invariants.md) (the modular surface: left and right \(\mathrm{SL}_2(\mathbb{Z})\)) and [euclidean-moduli-invariants.md](euclidean-moduli-invariants.md) (the torus: left translations \(\mathbb{Z}[i]\)). Here the Gaussian Schmidt circles are projected onto the **Riemann sphere** — throughout, the sphere of diameter \(1\) tangent to \(\mathbb{C}\) at the origin, north pole \((0,0,1)\), with stereographic projection from the north pole — and the left factor of the moduli problem becomes the **spherical symmetries inside the Bianchi group**,
$$
\Gamma_{\mathrm{sph}} \;=\; \mathrm{SU}(2)\,\cap\,\mathrm{SL}_2(\mathbb{Z}[i])
\;=\; \Bigl\{ \pm\begin{pmatrix}1&0\\0&1\end{pmatrix},\ \pm\begin{pmatrix}i&0\\0&-i\end{pmatrix},\ \pm\begin{pmatrix}0&1\\-1&0\end{pmatrix},\ \pm\begin{pmatrix}0&i\\i&0\end{pmatrix} \Bigr\},
$$
a **finite** group (the condition \(|a|^2+|b|^2 = 1\) on Gaussian integers forces one entry to be a unit): in \(\mathrm{PSL}_2\) it is the Klein group \(V_4 = \{z,\ -z,\ 1/z,\ -1/z\}\), the rotations by \(\pi\) about the three coordinate axes of the sphere. We study functions of \(X \in \mathrm{SL}_2(\mathbb{C})\) invariant under
$$
X \;\longmapsto\; u\,X\,\gamma', \qquad u \in \Gamma_{\mathrm{sph}},\ \gamma' \in \mathrm{SL}_2(\mathbb{Z}),
$$
evaluated on the Bianchi group. The compactness of the sphere and the finiteness of \(\Gamma_{\mathrm{sph}}\) change the flavor of every answer: counts are finite with **no quotient taken at all**, the class field theory runs in a **third discriminant family** \(-4(\ell^2+1)\), and the role of the hyperbolic \(\varepsilon = n + \sqrt{n^2-1}\) and the lemniscatic \(\Omega\) is taken by the **negative-Pell unit** \(\varepsilon_\ell = \ell + \sqrt{\ell^2+1}\). This document: (1) computes the number of Schmidt circles of a given spherical radius (the answer is \(4H(4(\ell^2+1))\)); (2) constructs the spherical position/shape/rotation invariants; (3) answers the level-polynomial question for both the shape and the rotation invariant — the shape answer is a clean *yes* (a product of class polynomials, each exactly once), the rotation answer is a *yes with a genuinely new caveat*: the phases form only **half** of the natural Galois root system, the completion being by Pell-unit conjugates that no circle realizes. All numerics: [scripts/spherical_moduli_invariants.py](scripts/spherical_moduli_invariants.py) (mpmath, 60–560 digits, absolute-error certification policy of the euclidean study; the structural statements in exact integer arithmetic).

Throughout, \(K_\ell = \mathbb{Q}(\sqrt{-(\ell^2+1)})\), \(H(N)\) is the Hurwitz class number (discriminant \(-N\)), \(H_D(x)\) the (ring) class polynomial, and a Schmidt circle is written in the Hermitian normalization of [circle-classification.md](circle-classification.md): curvature \(2q\), curvature-center \(\zeta = x+iy \equiv i \pmod 2\), co-curvature \(2m\), with \(|\zeta|^2 = 4qm+1\).

## 1. The projection, the radius, and the level

**Proposition 1 (spherical radius).** Let \(\omega\) be a circle with Hermitian data \(M = \begin{pmatrix} 2q & -\zeta \\ -\bar\zeta & 2m\end{pmatrix}\), \(\det M = -1\). Its stereographic image on the sphere is a circle of angular radius \(\theta\) (the angle subtended at the sphere's center; the geodesic radius is \(\theta/2\)) with
$$
\boxed{\;\cot\theta \;=\; q + m \;=:\; \ell\;}
$$
— *half the sum of curvature and co-curvature*, an integer for every Schmidt circle. Lines \(\operatorname{Im} z = k\) have \(\ell = |k|\) (circles through the north pole); the real line itself is the single great circle (\(\ell = 0\)).

*Proof.* A point \((x_0, y_0, t)\) of the sphere satisfies \(x_0^2+y_0^2+t^2 = t\) and projects from \(z\) with \(|z|^2 = t/(1-t)\), \(z = (x_0+iy_0)/(1-t)\). Substituting into \(2q|z|^2 - \zeta\bar z - \bar\zeta z + 2m = 0\) and clearing \((1-t)\) shows the image lies in the plane \(x\,x_0 + y\,y_0 + (m-q)\,t = m\). Its distance to the center \((0,0,\tfrac12)\) is \(\tfrac{|m+q|}{2}\bigl(|\zeta|^2 + (m-q)^2\bigr)^{-1/2}\), and \(|\zeta|^2 + (m-q)^2 = 4qm + 1 + (m-q)^2 = 1 + (m+q)^2\); dividing by the radius \(\tfrac12\) gives \(\cos\theta = |q+m|/\sqrt{1+(q+m)^2}\). \(\square\)

Since \(M\sim -M\), an unoriented circle has \(|\ell| = \cot\theta\) with \(\theta \le \pi/2\); we normalize oriented circles to \(\ell > 0\) and call \(\ell\) the **spherical level**. Two immediate consequences set the scene:

- **Quantization.** The spherical radii of Schmidt circles are quantized to the integer-cotangent angles \(\theta_\ell = \operatorname{arccot} \ell\); no other radius occurs. (Determinant \(-1\) forces \(qm \ge 0\), so \(q, m \ge 0\) after orientation and \(\ell \ge 0\).)
- **The Pell unit.** The eigenvalues of \(M\) are \(\ell \pm \sqrt{\ell^2+1}\): the circle's **poles** (the axis endpoints on the sphere, i.e. the fixed points of the rotations preserving it) are the eigenvectors of \(M\), explicitly \(p_\pm = \zeta\bigl((q-m) \mp \sqrt{\ell^2+1}\bigr)^{-1}\), and the eigenvalue \(\varepsilon_\ell := \ell + \sqrt{\ell^2+1}\) is the fundamental solution of the **negative Pell equation** \(a^2 - (\ell^2+1)b^2 = -1\) (the classical family \(d = \ell^2+1\), \((a,b) = (\ell,1)\)). Moreover \(W := iM \in \mathrm{SL}_2(\mathbb{Z}[i])\) — Gaussian entries, determinant \(+1\) — is the Möbius map **antipodal map \(\circ\) inversion in \(\omega\)** (stereographically, central symmetry composed with the reflection in the plane of \(\omega\)); it fixes the poles, and
$$
\operatorname{tr} W = 2\ell\, i, \qquad \text{multiplier } = -\varepsilon_\ell^{\,2},
$$
the exact spherical analogue of the hyperbolic trace identity \(\operatorname{tr}(X\bar X^{-1}) = -2n\) with its multiplier \(\varepsilon^2 = e^{\ell_{\mathrm{geod}}}\). The real quadratic order \(\mathbb{Z}[\sqrt{\ell^2+1}]\) and the imaginary \(K_\ell\) assemble into the CM field \(B_\ell := \mathbb{Q}(i, \sqrt{\ell^2+1}) = \mathbb{Q}(i, \sqrt{-(\ell^2+1)})\) — the counterpart of the hyperbolic \(B = \mathbb{Q}(i,\sqrt{n^2-1})\).

Machine verification (experiment B): the radius formula against fitted circumcircles of projected points, to 40+ digits.

## 2. Question 1: how many Schmidt circles of a given spherical radius?

The sphere is compact, so — unlike the modular surface and the torus — the census needs no fundamental domain: we count **all** Schmidt circles (and lines; on the sphere they are circles too) of level \(\ell\).

> **Theorem 1 (spherical census).** For \(\ell \ge 1\), the number of Schmidt circles of spherical radius \(\operatorname{arccot}\ell\) on the Riemann sphere is
> $$
> \boxed{\;N_{\mathrm{sph}}(\ell) \;=\; \tfrac13\,r_3(\ell^2+1) \;=\; 4\,H\bigl(4(\ell^2+1)\bigr),\;}
> $$
> where \(r_3\) counts lattice points on the sphere of radius \(\sqrt{\ell^2+1}\) in \(\mathbb{Z}^3\) and \(H\) is the Hurwitz class number. At \(\ell = 0\) there is the single great circle \(\hat{\mathbb{R}}\).

*Proof.* By the classification theorem, level-\(\ell\) circles are exactly the data \((q, m, \zeta)\) with \(q+m = \ell\), \(q, m \ge 0\), \(\zeta = x+iy\), \(x\) even, \(y\) odd, \(x^2+y^2 = 4qm+1\) (the congruence \(x^2+y^2 \equiv 1 \bmod 4q\) is automatic, and every admissible datum occurs; \(q = 0\) gives the lines \(\operatorname{Im} z = \pm\ell\)). Substitute \(t := q - m\) (so \(t \equiv \ell \bmod 2\)): the constraint becomes
$$
x^2 + y^2 + t^2 \;=\; \ell^2 + 1, \qquad x \equiv 0,\ y \equiv 1,\ t \equiv \ell \pmod 2,
$$
and this is a bijection (\(y\) odd excludes \(x = y = 0\), so \(|t| \le \ell\) is automatic). For \(\ell\) even, \(\ell^2+1 \equiv 1 \pmod 4\): every representation has exactly one odd coordinate, and permuting coordinates shows the three positions of the odd one are equinumerous — our parity pattern selects "odd coordinate in position \(y\)", one third of all. For \(\ell\) odd, \(\ell^2+1 \equiv 2 \pmod 4\): exactly two odd coordinates, and the same argument applies to the position of the even one. Hence \(N_{\mathrm{sph}}(\ell) = r_3(\ell^2+1)/3\). Finally \(\ell^2+1 \equiv 1, 2 \pmod 4\) always (and never \(\equiv 7 \bmod 8\)), so Gauss's theorem \(r_3(n) = 12\,H(4n)\) applies. \(\square\)

| \(\ell\) | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| \(H(4(\ell^2+1))\) | 1 | 2 | 2 | 4 | 6 | 2 | 7 | 8 | 4 | 14 |
| \(N_{\mathrm{sph}}(\ell)\) | 4 | 8 | 8 | 16 | 24 | 8 | 28 | 32 | 16 | 56 |

(Machine-verified for all \(\ell \le 20\) by exact enumeration, experiment A.) Two refinements:

- **The finite symmetry acts freely.** \(\Gamma_{\mathrm{sph}}\) acts on level-\(\ell\) circles by \((q,m,\zeta) \mapsto (q,m,-\zeta),\ (m,q,\bar\zeta),\ (m,q,-\bar\zeta)\); a fixed point would force \(\zeta\) real (impossible: \(y\) odd) or \(\zeta = yi\), \(q=m\) with \(y^2 = \ell^2+1\) (impossible for \(\ell \ge 1\)). So the action is free and the count on the quotient orbifold (the sphere modulo its three \(\pi\)-rotations) is **exactly \(H(4(\ell^2+1))\)** — the spherical analogue of "one class number per level".
- **Antipodal symmetry.** The antipodal map \(z \mapsto -1/\bar z\) preserves the arrangement (it is \((z\mapsto -1/z) \circ \mathrm{conj}\), both symmetries of \(\mathcal S\)) and pairs the level-\(\ell\) circles freely.

**The weight-\(3/2\) three-slice.** The hyperbolic census reads Zagier's trace form along \(d = n^2-1\), the euclidean along \(d = 4n^2\); the spherical census reads it along
$$
d \;=\; 4(\ell^2+1) \;=\; 4\ell^2 + 4,
$$
with weights \(3 : 2 : 4\) on the three slices \((3H(n^2-1),\ 2h(-4n^2),\ 4H(4\ell^2+4))\). The trace itself appears through the shape invariant of §3: summing \(j - 744\) over the level's classes gives **Zagier's \(t(4(\ell^2+1))\)** (experiment C, certified):

| \(\ell\) | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| \(t(4(\ell^2+1))\) | 7256 | 1262512 | 425691312 | 178211037024 | 82028232170000 | 39660183801070512 |

matching the published values (\(t(8) = 7256\), \(t(20) = 1262512\)). The three slices meet in exactly one point: \(n^2 - 1 = 4(\ell^2+1)\) forces \((n-2\ell)(n+2\ell) = 5\), i.e. \((n, \ell) = (3, 1)\) — **the discriminant \(-8\) is seen twice by the arrangement**, once hyperbolically at level \(3\), once spherically at the \(45°\)-circles \(\ell = 1\); the hyperbolic phase there was \(u = -1\) exactly, the spherical one is the integer of §5.2(4).

## 3. The moduli problem and the six invariants

\(\Gamma_{\mathrm{sph}}\) is finite, so the double coset space \(\Gamma_{\mathrm{sph}}\backslash \mathrm{SL}_2(\mathbb{C})/\mathrm{SL}_2(\mathbb{Z})\) is six-dimensional, like its two companions. The right-invariant content of \(X = \begin{pmatrix} a & b\\ c & d\end{pmatrix}\) is the pair of columns as a \(\mathbb{Z}\)-lattice in \(\mathbb{C}^2\); identifying \(\mathbb{C}^2\) with the quaternions carries it to a rank-2 sublattice of the **Lipschitz order** \(\mathbb{Z}\langle 1, i, j, k\rangle\), with left multiplication by unit quaternions as the \(\mathrm{SU}(2)\)-action. The complete left-\(\mathrm{SU}(2)\)-invariant is the Gram datum
$$
H \;:=\; X^\dagger X \;=\; \begin{pmatrix} A & B \\ \bar B & C\end{pmatrix},
\qquad A = |a|^2 + |c|^2,\quad C = |b|^2+|d|^2,\quad B = \bar a b + \bar c d,
$$
the **pullback of the round (Fubini–Study) form of the sphere**. Under right multiplication \(H \mapsto \gamma'^{\mathsf T} H \gamma'\); decomposing \(H = S + i(\operatorname{Im}B)\,J_0\) into its real symmetric and antisymmetric parts, the scalar \(\operatorname{Im} B\) is right-invariant and the form transforms by integral congruence. The dictionary with §1 is an identity of Gaussian integers:
$$
q = \operatorname{Im}(c\bar d), \quad m = \operatorname{Im}(a \bar b)
\qquad\Longrightarrow\qquad
\ell \;=\; q + m \;=\; -\operatorname{Im} B .
$$

> **Definition (shape).** The **shape form** of \(X\) is the integral binary quadratic form
> $$
> S \;=\; \operatorname{Re}(X^\dagger X) \;=\; [\,A,\ 2\operatorname{Re}B,\ C\,],
> \qquad \operatorname{disc} S = -4\bigl(1 + \ell^2\bigr),
> $$
> positive definite (always: \(A > 0\)); its root \(z_S = \bigl(-\operatorname{Re}B + i\sqrt{\ell^2+1}\bigr)/A \in \mathbb{H}\) is the **shape point** and \(\beta := j(z_S)\) the shape modulus. In quaternion terms \(S\) is the restriction of the reduced norm to the column lattice.

So a level-\(\ell\) Bianchi point carries a CM point of discriminant \(-4(\ell^2+1)\): **the third discriminant family**. The five "visible" invariants are the circle on the sphere (radius \(\theta\) + the polar axis: 3 real parameters; the poles \(p_\pm\) are \(B_\ell\)-rational — the *position*) and \(\beta\) (2 parameters — the *shape*); as in both companion studies, one dimension is missing, and it is a phase — the *rotation*.

> **Theorem 2 (structure of a level).** For \(\ell \ge 1\):
> 1. Each unoriented level-\(\ell\) circle underlies exactly two right-\(\mathrm{SL}_2(\mathbb{Z})\)-cosets of Bianchi matrices (the two orientations = the two caps bounded by the circle, exchanged by \(X \mapsto X\operatorname{diag}(i,-i)\)); \(\Gamma_{\mathrm{sph}}\) acts freely, and the level-\(\ell\) double cosets at orientation \(+\ell\) number \(H(4(\ell^2+1))\).
> 2. The shape class \([S] \in\) {form classes of discriminant \(-4(\ell^2+1)\)} is a well-defined double-coset invariant, and the induced map
> $$
> \{\text{level-}\ell \text{ double cosets}\} \;\longrightarrow\; \{\text{ALL reduced forms of disc } -4(\ell^2+1)\}
> $$
> is a **bijection** — every class, including the imprimitive strata (the suborder classes \(f^2 \mid \ell^2+1\)), occurs **exactly once**.
> 3. Complex conjugation of circles (\(D \mapsto \bar D\)) inverts the shape class; the antipodal map preserves it while reversing orientation. The lines \(\operatorname{Im} z = \pm\ell\) carry the principal class \([1, 0, \ell^2+1]\), with shape point \(\tau_\ell := i\sqrt{\ell^2+1}\).

Parts 1 and 3 and the count in 2 are proved (free action as in §2; the two-coset claim from \(\mathrm{Stab}(\hat{\mathbb{R}}) = \mathrm{SL}_2(\mathbb{Z}) \sqcup \operatorname{diag}(i,-i)\mathrm{SL}_2(\mathbb{Z})\); the counting shell "\(H\) cosets vs \(H\) classes" is Theorem 1); that the class map is a bijection — rather than merely a map between equinumerous sets — is machine-verified **exactly** for every level \(\ell \le 20\) (experiment A, integer arithmetic end to end: descent word for each circle, exact form reduction, exact set equality). A full proof should follow Gauss's quaternionic proof of \(r_3 = 12H\), which is precisely the statement that primitive Lipschitz pairs with prescribed unimodular pairing sweep each form class equally often. The bijection is the spherical counterpart of "every ring class exactly twice" (euclidean Theorem 1) and of the Hilbert-class-field orbit of the hyperbolic \(\beta_1\) — with the new feature that **the imprimitive strata are on an exactly equal footing** (at \(\ell = 7\): \(\ell^2+1 = 50\), and one of the seven classes is the conductor-\(5\) stratum \([5,0,10]\), carrying \(j = 8000\) of discriminant \(-8\)).

## 4. The sixth invariant: the phase

The euclidean phase glued the weight-2 kernel \(j'\) to the residue datum \(1/c^2\); the spherical phase glues it to the **second Gram form**, the complex symmetric
$$
T \;:=\; X^{\mathsf T} X \;=\; \begin{pmatrix} a^2+c^2 & ab+cd \\ ab+cd & b^2+d^2 \end{pmatrix} \in \mathrm{Sym}_2(\mathbb{Z}[i]),\qquad \det T = 1,
$$
the pullback under \(X\) of the quadratic \(z^2+1\) whose zeros \(\{\pm i\}\) are the poles of the base great circle \(\hat{\mathbb{R}}\).

> **Definition (phase).** \(\displaystyle \Theta(X) \;:=\; T(z_S, 1)\cdot j'(z_S)\), where \(T(z,1) = T_{11}z^2 + 2T_{12}z + T_{22}\) and \(z_S\) is the shape point of \(X\).

**Proposition 2 (invariance).** \(\Theta(X\gamma') = \Theta(X)\) for all \(\gamma' \in \mathrm{SL}_2(\mathbb{Z})\), and \(\Theta(uX) = \chi(u)\,\Theta(X)\) for \(u \in \Gamma_{\mathrm{sph}}\), where \(\chi\) is the character with \(\chi = +1\) on \(\{1, \binom{0\ 1}{-1\ 0}\}\) and \(\chi = -1\) on \(\{\operatorname{diag}(i,-i), \binom{0\ i}{i\ 0}\}\). In particular \(\Theta^2\) is a two-sided invariant.

*Proof.* Right: \(S\) and \(T\) transform by the same congruence \(\gamma'^{\mathsf T}(\cdot)\gamma'\), so with \(z' = \gamma'^{-1}z_S\) (Möbius), \(T'(z',1) = (rz'+s)^2\,T(z_S,1)\) while \(j'(z') = (rz'+s)^{-2} j'(z_S)\): the weights cancel exactly, with no cocycle. Left: \(H\) is unitarily invariant, so \(S, z_S\) are fixed; and \(u^{\mathsf T}u = \pm 1\) on \(\Gamma_{\mathrm{sph}}\) (computed case by case), so \(T \mapsto \pm T\). \(\square\)

The sign \(\chi\) is the exact analogue of the euclidean Borel sign: \(\operatorname{diag}(i,-i)\) is \(z \mapsto -z\), so **the two circles \(\pm\omega\) carry opposite phases**. The remaining structure (all machine-certified, experiment B; the resultant identity is *proved*):

1. **Fiber.** Fix the five invariants; the fiber is the circle of left rotations about \(\omega\)'s own axis. Along a rotation of the sphere by angle \(t\) about the axis, \(|\Theta|\) is constant to all digits and \(\arg\Theta\) moves **linearly at rate exactly \(1\)**: the phase *is* the rotation angle of the packing about the circle's axis. (Euclidean rate: 2 in the matrix parameter, i.e. 1 in the geometric rotation angle — the same statement; the hyperbolic rate \(\sqrt{n^2-1}\) is the loxodromic analogue.) So \((\theta,\ \text{axis},\ \beta,\ \arg\Theta)\) — six real dimensions — specify the double coset up to finite ambiguity.
2. **The resultant identity (proved).** For every \(X \in \mathrm{SL}_2(\mathbb{C})\), as a polynomial identity in the entries (verified symbolically in exact arithmetic on \(\det X = 1\)):
$$
T(z_S,1)\,T(\bar z_S, 1) \;=\; \frac{\operatorname{Res}(S, T)}{A^2}, \qquad
\boxed{\;\operatorname{Res}(S, T) \;=\; 4\bigl(y^2 - \ell^2 - 1\bigr) \;=\; -4\bigl((q-m)^2 + x^2\bigr)\;}
$$
(for Bianchi points, with \(\zeta = x+iy\)). In particular \(\Theta = 0\) would force \(q = m\), \(x = 0\), hence \(y^2 = \ell^2+1\): **the phase never vanishes on a level \(\ell \ge 1\)**.
3. **The \(T\)-norm lemma (certified at every class, \(\ell \le 6\)).** The Pell unit is exactly the archimedean discrepancy between the two roots of \(S\):
$$
\varepsilon_\ell\,\cdot\,\tfrac{A}{2}\,T(z_S, 1) \;\in\; K_\ell = \mathbb{Q}\bigl(\sqrt{-(\ell^2+1)}\bigr),
\qquad
N_{K_\ell/\mathbb{Q}}\bigl(\varepsilon_\ell \tfrac{A}{2} T(z_S,1)\bigr) = (q-m)^2 + x^2,
$$
so \(|T(\bar z_S,1)| = \varepsilon_\ell^2\,|T(z_S,1)|\) — the counterpart of the hyperbolic Norm Lemma \(|\mu|^2 = \varepsilon\,q_1/q_2\). On the polar lines this is exact algebra: the level-\(\ell\) representative \(X = \binom{i\ \ \ell}{0\ -i}\) of the line \(\operatorname{Im} z = \ell\) has
$$
T(z_S, 1) \;=\; -\,2\ell\,\varepsilon_\ell^{-1}, \qquad T(\bar z_S, 1) \;=\; 2\ell\,\varepsilon_\ell\,;
$$
**the principal anchor of the theory, with the negative-Pell unit in the role the real-quadratic \(\varepsilon\) played hyperbolically.**

## 5. Question 3: the level polynomials

### 5.1 The shape invariant: yes, and cleanly

> **Theorem 3 (shape polynomial).** For \(\ell \ge 1\), over the level-\(\ell\) double cosets (equivalently, over the level-\(\ell\) circles on the sphere modulo the free \(V_4\)),
> $$
> \prod_{\text{level-}\ell} \bigl(x - j(z_S)\bigr) \;=\; \prod_{f^2 \mid \ell^2+1} H_{-4(\ell^2+1)/f^2}(x) \;\in\; \mathbb{Z}[x],
> $$
> each ring class polynomial of the family **to the first power**. Over all \(4H\) circles of the sphere the product is the fourth power. In particular the shape values are algebraic integers forming full Galois orbits, and the answer to "are they the roots of a rational polynomial?" is **yes — an integer polynomial, one root per class.**

This is Theorem 2 plus the classical integrality of class polynomials. Certified numerically for \(\ell \le 8\) (experiment C; each coefficient with 20+ spare digits in the absolute-error sense), including the first imprimitive case \(\ell = 7\): \(\prod = H_{-200}(x)\,H_{-8}(x)\), degree \(6+1 = 7\). The trace corollary is the \(t(4(\ell^2+1))\)-slice of §2. Compared with the companions: hyperbolic — Hilbert class polynomial \(H_{1-n^2}\) (fields varying with the level, fixed maximal order); euclidean — \(H_{-4n^2}^2\) (fixed field \(\mathbb{Q}(i)\), conductor = level); spherical — \(\prod_f H_{-4(\ell^2+1)/f^2}\): **fields varying along the shifted squares \(\ell^2+1\), all strata exactly once**. The Schmidt arrangement realizes all three classical aspects of CM theory on its three background geometries.

### 5.2 The rotation invariant: normalization and laws

Normalize by the period of the level's own principal class — the shape point \(\tau_\ell\) of the polar lines (the analogue of the euclidean \(\Omega = 2\pi\eta(i)^4\), whose normalizing point \(i\) was the shape point of the Ford stratum):
$$
u^2(X) \;:=\; \frac{\Theta(X)^2}{(2\pi)^2\,\eta(\tau_\ell)^8}, \qquad \tau_\ell = i\sqrt{\ell^2+1}
$$
(\(\eta(\tau_\ell)^8 = \Delta_q(\tau_\ell)^{1/3} > 0\); \(\Theta^2\) is the honest two-sided invariant, and \(u^2\) is indexed by the \(H\) classes and the orientation \(\pm\ell\)). From \(j'^2 = -(2\pi)^2 j(j-1728)E_4\),
$$
u^2 \;=\; -\,T(z_S,1)^2\;\beta(\beta - 1728)\;\frac{E_4(z_S)}{\eta(\tau_\ell)^8}\,,
$$
so \(u^2\) is **algebraic** (Shimura, weight-4 ratio at two points of \(K_\ell\)), and — certified at every computed class — an **algebraic integer**. Note \(3 \nmid \ell^2+1\) for every \(\ell\) (squares are \(0,1 \bmod 3\)), so the \(\gamma_2\)-obstruction of the euclidean study *never occurs spherically*. The certified laws (experiment D, all levels computed; \(\varepsilon = \varepsilon_\ell\)):

1. *(orientation / cap-swap twist)* \(\;u^2(\text{level } -\ell) = \varepsilon^4\, u^2(\text{level } +\ell)\) **at every class**: the two caps bounded by one circle carry phases differing by exactly \(\varepsilon^2\) — the spherical twist law, in the role of the hyperbolic \(u_f u_{\mathfrak{r}f} = 1\) (there the twist cost \(e^{-\ell_{\mathrm{geod}}} = \varepsilon^{-2}\); here the cap-swap costs \(\varepsilon^{+2}\), and nothing cancels it).
2. *(mirror)* \(\;\overline{u^2[S]} = u^2[S^{-1}]\) within each orientation: complex conjugation of circles is inversion in the class group, as in both companions. On ambiguous classes \(u^2\) is real.
3. *(companion family)* On the \(\mathrm{PGL}_2\)-companion \(i\mathcal{S}\) (which has its own level-\(\ell\) census, in bijection with the same classes), the phase differs by a **ratio in \(K_\ell\)**, orientation-independent, equal to \(-(\ell^2+1)/\ell^2\) at the principal class and rational at every ambiguous class (e.g. \(\ell=3\): \(-10/9\) and \(-10\); \(\ell = 4\): \(-17/16\), \(1/16\), and \(-(4 \pm 3\sqrt{-17})/8\) at the complex pair — norm \((13/8)^2\)). The quarter-turn does *not* leave the field of the level.
4. *(anchor, exact)* At \(\ell = 1\) (the \(45°\) circles, discriminant \(-8\), \(\beta = j(i\sqrt2) = 8000\)):
$$
u^2(\pm) \;=\; -\,2^{17}\,5^4\,7^2\;\varepsilon_1^{\mp 2}, \qquad \varepsilon_1 = 1+\sqrt2,
$$
from \(u^2 = -4\ell^2\varepsilon^{\mp2}\,\beta^{4/3}(\beta-1728)\) with \(\beta^{1/3} = 20\), \(\beta - 1728 = 2^7 7^2\). The support \(\{2, 5, 7\}\) is: the level prime, the primes of \(H_{-8}(0) = -2^6 5^3\), and of \(H_{-8}(1728) = -2^7 7^2\) — Gross–Zagier collisions of the \(45°\) stratum with \(j = 0\) and \(j = 1728\).

### 5.3 The level polynomials: certified, with a genuinely new phenomenon

Write \(P_\ell(y) := \prod (y - u^2)\) over the \(2H\) phases of the level (both orientations). The certified record:

| \(\ell\) | \(h\) | \(\mathrm{Cl}(-4(\ell^2+1))\) | coefficient field of \(P_\ell\) | completed polynomial |
|---|---|---|---|---|
| 1 | 1 | trivial | \(\mathbb{Q}\) | \(P_1 \in \mathbb{Z}[y]\), irreducible |
| 2 | 2 | \(\mathbb{Z}/2\) | \(\mathbb{Q}(\sqrt5)\) | \(N_{F/\mathbb{Q}}P_2 \in \mathbb{Z}[y]\), deg 8 = four quadratics |
| 3 | 2 | \(\mathbb{Z}/2\) | \(\mathbb{Q}(\sqrt5)\) | \(N_{F/\mathbb{Q}}P_3 \in \mathbb{Z}[y]\), deg 8 = two quartics |
| 4 | 4 | \(\mathbb{Z}/4\) | beyond quadratic (certified) | — |
| 6 | 2 | \(\mathbb{Z}/2\) | \(\mathbb{Q}(\sqrt{37})\) | \(N_{F/\mathbb{Q}}P_6 \in \mathbb{Z}[y]\), deg 8 = four quadratics |

- **\(\ell = 1\).** The two phases are a full Galois orbit and
$$
P_1(y) \;=\; y^2 + 24084480000\,y + 16112838246400000000 \;=\; y^2 + 6\cdot 2^{17}5^47^2\,y + \bigl(2^{17}5^47^2\bigr)^2,
$$
monic, integral, irreducible — one level, one integer polynomial, exactly as in the companions.
- **\(\ell = 2\)** (\(F = \mathbb{Q}(\sqrt5)\), which is simultaneously the real genus field of \(-20\) *and* the Pell field \(\mathbb{Q}(\sqrt{\ell^2+1})\)): the phases are exact \(\mathbb{Z}[\sqrt5]\)-integers,
$$
u^2\bigl([1,0,5], +\bigr) = -76984628019200 - 34429848780800\sqrt5,\qquad
u^2\bigl([2,2,3], +\bigr) = 9468658830540800 - 4234512970547200\sqrt5,
$$
the \(e_k\) of the \(2H\)-multiset lie in \(\mathbb{Z}[\sqrt5]\) (all four certified), and the completed rational polynomial is
$$
N_{F/\mathbb{Q}}P_2(y) = y^8 + 30736388299161600\,y^7 - 935977682193755759027976929280000\,y^6 - \cdots + 2^{160}5^{16}11^{24}19^{8},
$$
which factors over \(\mathbb{Q}\) into **four quadratics with one and the same constant term**
$$
m_2 \;=\; u^2\cdot\sigma_{\sqrt5}(u^2) \;=\; -\,2^{40}\,5^4\,11^6\,19^2 \qquad\text{for every class and orientation:}
$$
**the level norm** — the spherical analogue of the hyperbolic norm-one law \(u_f u_{\mathfrak{r}f} = 1\), with the fixed \(S\)-integer \(m_\ell\) in place of \(1\).
- **\(\ell = 3, 6\)**: same picture (coefficients certified in \(\mathbb{Z}[\sqrt5]\), resp. \(\mathbb{Z}[\sqrt{37}]\) — the real genus fields of \(-40\), \(-148\); completed degree-8 integer polynomials; at \(\ell = 6\) again four rational quadratics sharing the constant \(m_6 = -2^{40}3^{26}5^{8}7^{4}11^{6}47^{4}67^{2}107^{4}139^{2}\); at \(\ell = 3\) two rational **quartics** — see below). Residuals of the certified polynomials at the numerical phases: below \(10^{-185}\), \(10^{-210}\), \(10^{-480}\).

**The valuation law of the level norm.** With \(H_0 := \prod_{\mathfrak c} j(\mathfrak c)\) and \(H_{1728} := \prod_{\mathfrak c}(j(\mathfrak c) - 1728)\) (over the \(\sigma\)-orbit of the class pair), the certified norms obey, at every prime,
$$
v_p(m_\ell) \;=\; 2\,v_p(4\ell^2)\;+\;\tfrac43\,v_p(H_0)\;+\;v_p(H_{1728}),
$$
the exact transplant of the euclidean numerator law (§5.5 there): the phase-norm is supported at the level primes and at the **Gross–Zagier collisions of the level's CM points with \(j = 0\) (exponent \(4/3\)) and \(j = 1728\) (exponent \(1\))**. Verified at \(\ell = 1, 2, 6\) prime by prime (e.g. \(\ell = 6\): \(v_3 = 2\cdot 2 + \tfrac43\cdot 6 + 14 = 26\), \(v_2 = 2\cdot4 + \tfrac43\cdot 12 + 16 = 40\), \(v_5 = \tfrac43\cdot 6 = 8\), and \(7^4, 11^6, 47^4, 67^2, 107^4, 139^2\) all on the nose).

**The new phenomenon: the level is half an orbit.** In both companion studies the level's phases formed a Galois-stable multiset — one level, one integer polynomial, roots exactly the phases. **Spherically this fails for \(\ell \ge 2\), provably and exactly.** At \(\ell = 2\) the Galois conjugate \(\sigma_{\sqrt5}(u^2) = m_2/u^2\) of each phase is *positive*, while all four phases of the level are negative real numbers: the conjugates are not phases of any circle. The completed root system is
$$
\{u^2_i\} \;\sqcup\; \{m_\ell/u^2_i\}
$$
— each phase paired with a **Pell-twisted virtual partner** — and only the union is a root system of an integer polynomial. Moreover the Galois action on the phases carries a **unit-valued cocycle**: exactly, at \(\ell = 2\),
$$
\sigma_{\sqrt5}\Bigl(u^2\bigl([1,0,5], +\bigr)\Bigr) \;=\; -\,\varphi^{-2}\; u^2\bigl([2,2,3], -\bigr),
\qquad \varphi = \tfrac{1+\sqrt5}{2},\ \ \varepsilon_2 = \varphi^3
$$
(golden-ratio cocycle!) — translation by the nontrivial class, orientation flip, and multiplication by \(-\varphi^{-2}\). A unit of infinite order cannot be removed by passing to powers of \(u\) (the hyperbolic \(\xi\)-torsion argument fails: its Kronecker step needed unitarity at all archimedean places, and here the cocycle lives in a *real* field). This is the structural price of the compact geometry: the two archimedean places of the Pell field see the two caps of each circle, and only one cap is "the disk".

So the answer to Question 3 for the rotation invariant is: **yes — the normalized squared phases of a level are roots of one explicit monic integer polynomial (certified through the levels above, of degree \(2H\cdot[F_\ell:\mathbb{Q}]\)) — but, unlike the hyperbolic and euclidean settings, they constitute only half of its root system, the other half being their Pell-unit conjugates \(m_\ell/u^2\), which are not phases of circles.** The factorization type over \(\mathbb{Q}\) is itself structured: the irreducible factors are the Galois orbits of the individual phases — quadratics \(y^2 - \mathrm{Tr}(u^2)\,y + m_\ell\) when the Pell field coincides with the coefficient field (\(\ell = 2, 6\): \(\sqrt{\ell^2+1}\) lands in the genus field iff \(\ell\) is even), quartics pairing the two orientations when it does not (\(\ell = 3\): \(\sqrt{10} \notin \mathbb{Q}(\sqrt5)\), and \(\sigma_{\sqrt2}\) realizes the cap swap \(u^2 \mapsto \varepsilon^4 u^2\) inside one orbit). At \(\ell = 4\) (class group \(\mathbb{Z}/4\), the first beyond 2-torsion) the coefficient field of \(P_4\) is certified to exceed every quadratic field — the growth of the coefficient field with the class group mirrors the genus-character fields of the hyperbolic §5.8, and pinning it exactly is the first open problem below.

### 5.4 What changed relative to the hyperbolic and euclidean phases

| | hyperbolic | euclidean | **spherical** |
|---|---|---|---|
| left group | \(\mathrm{SL}_2(\mathbb{Z})\) | translations \(\mathbb{Z}[i]\) | \(\Gamma_{\mathrm{sph}}\) finite (\(V_4\)) |
| geometry | modular surface | torus | sphere (no quotient) |
| radius invariant | \(\coth\rho = n\) | curvature \(2n\) | \(\cot\theta = \ell\) |
| census | \(3H(n^2-1)\) | \(2h(-4n^2)\) | \(4H(4(\ell^2+1))\) |
| CM data of a level | disc \(1-n^2\) | conductor \(n\) over \(\mathbb{Q}(i)\) | disc \(-4(\ell^2+1)\), all strata once |
| trace slice | \(t(n^2-1)\) | \(t(4n^2)\) | \(t(4\ell^2+4)\) |
| kernel gluing | two kernels + derivative | one kernel + residue | one kernel + the form \(T = X^{\mathsf T}X\) |
| fiber rate | \(\sqrt{n^2-1}\) | 2 (weight) | 1 (the rotation itself) |
| normalizing constant | \(\varepsilon = n+\sqrt{n^2-1}\) (Pell \(+1\)) | \(\Omega = \varpi^2/\pi\) | \(\eta(\tau_\ell)\)-period and \(\varepsilon_\ell = \ell+\sqrt{\ell^2+1}\) (Pell \(-1\)) |
| twist law | \(u_f u_{\mathfrak{r}f} = 1\) | \(\Delta\)-mass law | cap swap \(u^2_- = \varepsilon^4 u^2_+\); norm \(u^2\sigma(u^2) = m_\ell\) |
| level polynomial | integer, non-monic, roots = phases | monic integer, roots = phases | monic integer, **phases = half the roots** |
| collision primes | GZ\((1-n^2; -3), (-4)\) in denominators | GZ in numerators over \(n\) | GZ in \(m_\ell\): \(2v(4\ell^2) + \tfrac43 v(H_0) + v(H_{1728})\) |

## 6. Research outlook

**Small (days to weeks).**
1. *The bijection, unconditionally.* Upgrade Theorem 2(2) from "verified \(\ell \le 20\)" to a theorem by the quaternionic counting argument (primitive Lipschitz pairs with \(\bar v_1 v_2 = \bar B - j\), Gauss-style); Theorem 1 already provides the mass check.
2. *The cap-swap law and the \(T\)-norm lemma.* Prove \(u^2_- = \varepsilon^4 u^2_+\) and \(\varepsilon_\ell \tfrac A2 T(z_S,1) \in K_\ell\) exactly (the resultant identity is proved; what remains is the assignment of \(\varepsilon^{-1}\) to the *upper* root, an archimedean statement in the style of the euclidean Lemma D).
3. *The \(\ell = 4\) coefficient field.* Identify \(F_4\) (a subfield of the real part of the ring class field of \(-68\)); conjecturally the fixed field of the multiset stabilizer under the dihedral action with the \(\varepsilon\)-cocycle — the spherical version of the genus-character fields of moduli-invariants.md §5.8.
4. *The companion-family ratio.* Close the certified table of §5.2(3) into a formula (the \(K_\ell\)-element \(-(\text{norm-form data})/\ell^2\); its norm was \((13/8)^2\) at the complex pair of \(\ell = 4\)).

**Medium (weeks to months).**
5. *The cocycle calculus.* The exact law \(\sigma(u^2) = (\text{unit cocycle})\cdot u^2(\text{translate, flip})\) begs for a spherical Shimura-reciprocity bookkeeping: transport of the pair \((S, T)\) under the idelic action, with the \(T\)-norm lemma controlling the \(\varepsilon\)-content (the analogue of the \(\mu\)-cocycle-ideal computation of the hyperbolic §5.5, whose Kronecker step is exactly what fails — understanding *what replaces it* is the core question). The golden-ratio cocycle \(-\varphi^{-2}\) at \(\ell = 2\) is the first data point.
6. *The virtual partners.* The completed root system pairs each circle-phase with \(m_\ell/u^2\). Find the geometric carrier of the missing half: candidates are the transpose moduli problem (rows vs columns, \(XX^\dagger\) — the \(\sigma\)-side of euclidean §6.8), or genuinely non-geometric conjugates (the two archimedean places of \(\mathbb{Q}(\sqrt{\ell^2+1})\) as the two caps).
7. *Kronecker limit formula in the third aspect.* \(\log|u^2|\) against \(L'(0,\chi)\) for class characters of disc \(-4(\ell^2+1)\), with the \(\varepsilon^4\)-twist as the archimedean correction; the negative-Pell family \(d = \ell^2+1\) is classical territory (continued fraction of \(\sqrt{d}\) of period... the \((\ell, 1)\)-solution), and the phase gives it a CM shadow.
8. *The three-slice identity.* The censuses slice \(\sum H(d)q^d\) along \(n^2-1\), \(4n^2\), \(4\ell^2+4\). The Hurwitz–Kronecker relations \(\sum_t H(4m - t^2)\) mix exactly these progressions; a bijective proof **on the arrangement** now has all three families as concrete circles (outlook 6.14 of the euclidean study, with the third slice supplied).

**Large (a paper or program each).**
9. *Equidistribution on the sphere.* As \(\ell \to \infty\) the level-\(\ell\) circles equidistribute (the shape points are CM points of disc \(-4(\ell^2+1)\); Duke-type results apply); the phase adds the rotation angle above each circle. Formulate and prove joint equidistribution of (axis, \(\arg u\)) — the compact companion of the hyperbolic outlook 3.4, where the harmonic analysis is cleanest of all three.
10. *Other Bianchi groups.* Over \(\mathcal{O}_{\sqrt{-d}}\) the spherical unit group \(\mathrm{SU}(2)\cap\mathrm{SL}_2(\mathcal{O})\) is again finite; the level becomes \(\cot\theta \in \tfrac{1}{\sqrt{d}}\mathbb{Z}\)-type data and the discriminant family \(-4d(\ell^2+d)\)-flavored. One packing, three geometries, all fields: the full "Schmidt arrangements model CM theory" program of euclidean 6.13 with its third column filled in.

## 7. Files

- [scripts/spherical_moduli_invariants.py](scripts/spherical_moduli_invariants.py) — all experiments (the level is called `c` in the code). `python3 scripts/spherical_moduli_invariants.py` runs A (exact structure: census, free \(V_4\)-action, descent, shape bijection, mirror/antipodal laws, the resultant identity — all \(\ell \le 20\), exact), B (radius by projection; invariance and sign character; fiber; anchors; the \(T\)-norm lemma), C (shape polynomials and the trace slice, \(\ell \le 8\)); `... phase` runs D (phase laws, the certified \(\mathbb{Z}[\sqrt d]\)-coefficients, completed integer polynomials, level norms with factorizations, the \(\ell = 4\) negative) and E (the \(i\mathcal{S}\)-companion ratios). `... all` runs everything. Requires mpmath; D/E also sympy.
- Certification policy as in euclidean-moduli-invariants.md: integers and quadratic coordinates accepted only with large absolute-error margins (spare digits \(\ge \mathrm{dps}/6\)); PSLQ used only where (terms) × (coefficient digits) is far below working precision; the completed polynomials are re-verified against the numerical phases (residuals \(< 10^{-185}\)).
- The symbolic proof of the resultant identity (holomorphic and antiholomorphic entries as independent variables, \(\det X = 1\) and its conjugate imposed exactly, the difference simplified to \(0\) in sympy): `python3 scripts/spherical_moduli_invariants.py symbolic`.
