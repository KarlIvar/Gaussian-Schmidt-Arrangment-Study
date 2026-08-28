# The involution \(X \mapsto \bar X^{-1}\): the adjoint involution, the Cartan embedding, and a twisted inversion of class groups

We investigate the map
$$
\sigma : \mathrm{SL}_2(\mathbb{Z}[i]) \to \mathrm{SL}_2(\mathbb{Z}[i]), \qquad \sigma(X) = \bar X^{-1}
$$
(entrywise complex conjugation followed by matrix inverse). All identities below are proved and additionally machine-verified in exact arithmetic ([scripts/involution_experiments.py](scripts/involution_experiments.py), [scripts/involution_classmap.py](scripts/involution_classmap.py)). Throughout, \(M_\omega = \begin{pmatrix} A & B \\ \bar B & C\end{pmatrix}\) denotes the integral Hermitian matrix of a circle \(\omega\) as in [circle-classification.md](circle-classification.md), \(M_X := M_{X(\hat{\mathbb{R}})}\), \(J = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}\), and \(M_0 = iJ\) is the matrix of \(\hat{\mathbb{R}}\).

## 1. What \(\sigma\) is: the adjoint involution of the Hermitian form of \(\hat{\mathbb{R}}\)

\(\sigma\) is an **anti-automorphism of order two**: \(\sigma(XY) = \sigma(Y)\sigma(X)\), \(\sigma^2 = \mathrm{id}\). Its structure is best seen through the identity (valid for \(X \in \mathrm{SL}_2\), where \(X^{-1} = J^{-1}X^{\mathsf T}J\)):
$$
\sigma(X) \;=\; J^{-1} X^{\dagger} J \;=\; M_0^{-1} X^{\dagger} M_0 .
$$
That is: **\(\sigma\) is the adjoint anti-involution of the Hermitian form** \(h(u, v) = u^\dagger M_0 v = -2\operatorname{Im}(\bar u_1 v_2)\)-type pairing attached to the real line: \(h(Xu, v) = h(u, \sigma(X)v)\). In the language of central simple algebras, \(\sigma\) is a **unitary involution — an involution of the second kind** — on \(M_2(\mathbb{Q}(i))\): it is \(\mathbb{Q}\)-linear, reverses products, and restricts to Galois conjugation on the center \(\mathbb{Q}(i)\). (This is the standard terminology of Knus–Merkurjev–Rost–Tignol, *The Book of Involutions*.)

Two companion objects:

- The **automorphism** \(\theta(X) = \bar X = \sigma(X)^{-1}\) is the Galois involution — the *real structure* on \(\mathrm{SL}_2(\mathbb{C})\) whose fixed group is the split real form \(\mathrm{SL}_2(\mathbb{R})\). Integrally, \(\theta\) exhibits \((\mathrm{SL}_2(\mathbb{Z}[i]),\, \mathrm{SL}_2(\mathbb{Z}))\) as a **symmetric pair**.
- The unitary group of \(h\): a one-line computation (\(g^\dagger J g = J\) and \(g^{\mathsf T} J g = J\) force \(\bar g = g\)) gives
$$
\mathrm{SU}(M_0;\, \mathbb{Z}[i]) \;=\; \mathrm{SL}_2(\mathbb{Z}) .
$$
So the stabilizer of the *oriented* real line, the fixed group of \(\theta\), and the unitary group of \(M_0\) are one and the same group \(\Gamma := \mathrm{SL}_2(\mathbb{Z})\); circles of the Schmidt arrangement correspond to cosets \(X\Gamma\), and \(\Gamma\)-classes of circles to double cosets \(\Gamma X \Gamma\).

## 2. The fun fact, explained: \(X\bar X^{-1}\) is the inversion in the circle \(X(\hat{\mathbb{R}})\)

The map \(X \mapsto Y := X\,\bar X^{-1} = X\,\theta(X)^{-1}\) is the classical **Cartan embedding** of the symmetric space \(G/G^\theta\) into \(G\): its fibers are exactly the cosets \(X\,\mathrm{SL}_2(\mathbb{Z})\), i.e. **\(Y\) determines, and is determined by, the oriented circle \(X(\hat{\mathbb{R}})\)** — which is the conceptual form of the observed fun fact. The precise dictionary is an exact identity between the Cartan image and the Hermitian matrix of the circle:
$$
X \bar X^{-1} \;=\; -\,i\,J\,\overline{M_X}
\;=\;
\begin{pmatrix} -iB & -iC \\ iA & i\bar B \end{pmatrix},
\qquad
M_X = \begin{pmatrix} A & B \\ \bar B & C\end{pmatrix},
$$
so the entries of \(X\bar X^{-1}\) carry precisely the curvature \(A\), the curvature\(\times\)center \(-B\), and the co-curvature \(C\) — the center, radius, orientation, and nothing more. Rescaling by \(1/(iA)\) puts \(Y\) in the classical shape \(\begin{pmatrix} z_0 & r^2 - |z_0|^2 \\ 1 & -\bar z_0\end{pmatrix}\): the matrix of the **anti-holomorphic inversion in the circle** \(\omega = X(\hat{\mathbb{R}})\),
$$
\tau_\omega(z) = Y(\bar z), \qquad \tau_\omega = X \circ (\text{conj}) \circ X^{-1},
$$
which is geometrically obvious: \(\hat{\mathbb{R}}\) is the fixed locus of \(z \mapsto \bar z\), so \(X(\hat{\mathbb{R}})\) is the fixed locus of the conjugated reflection. Note \(Y\bar Y = I\): the Cartan image consists of *twisted involutions* (cocycles for the Galois action).

**A trace identity.** From the formula, \(\operatorname{tr}(X\bar X^{-1}) = -iB + i\bar B = 2\operatorname{Im} B = -2\,y\), where \(y = -\operatorname{Im}B\) is the (oriented) numerator of the center height — i.e. the invariant \(\alpha\) of [hyperbolic-counting.md](hyperbolic-counting.md), with sign:
$$
\boxed{\;\operatorname{tr}\bigl(X \bar X^{-1}\bigr) = -2\,\alpha_{\pm}\bigl(X(\hat{\mathbb{R}})\bigr).\;}
$$
So \(\alpha\) — previously seen as \(\coth\) of the hyperbolic radius — is *the trace on the Cartan embedding*. This has a clean inversive/hyperbolic meaning: \(\alpha(\omega)\) equals the **inversive product** of \(\omega\) with \(\hat{\mathbb{R}}\), and \(-Y = \tau_\omega \circ \tau_{\hat{\mathbb{R}}}\)-normalized is the product of the reflections of \(\mathbb{H}^3\) in the hemispheres over \(\omega\) and over \(\hat{\mathbb{R}}\): a **hyperbolic translation of length \(2\operatorname{arccosh}(n)\)** along their common perpendicular, with characteristic polynomial \(\lambda^2 + 2n\lambda + 1\) and real quadratic eigenvalue \(-(n + \sqrt{n^2-1})\), for \(\alpha = n\). Every Schmidt circle with \(\alpha = n\) thus determines a closed geodesic in the Bianchi orbifold \(\mathrm{PSL}_2(\mathbb{Z}[i]) \backslash \mathbb{H}^3\) with real trace field \(\mathbb{Q}(\sqrt{n^2-1})\).

## 3. Fixed points of \(\sigma\): inversions, Hilbert 90, and exactly three twisted classes

Since \(\sigma(X) = M_0^{-1}X^\dagger M_0\),
$$
\sigma(X) = X \iff M_0 X \text{ is Hermitian} \iff X = M_0 M \ \text{ with } M \text{ integral Hermitian},\ \det M = -1
$$
(\(\det\) is automatic from \(\det X = 1\)). So the fixed points of \(\sigma\) in \(\mathrm{SL}_2(\mathbb{Z}[i])\) are precisely the matrices \(M_0M\) of circle inversions \(z \mapsto M_0M(\bar z)\) — one for each integral Hermitian form of determinant \(-1\). The natural action preserving the fixed set is the **twisted conjugation** \(X \mapsto \bar g X g^{-1}\), and under \(X = M_0M\) it matches the usual action \(M \mapsto (g^{-1})^\dagger M g^{-1}\) on circles. Classifying fixed points up to twisted conjugacy is a Galois-cohomology problem (\(X\bar X = 1\) is the cocycle condition; over the field, Hilbert 90 makes it trivial; over \(\mathbb{Z}[i]\) it is a genus question). The answer, computed by exhaustive orbit enumeration:

> **The integral binary Hermitian forms of determinant \(-1\) fall into exactly three \(\mathrm{SL}_2(\mathbb{Z}[i])\)-classes** (verified completely for \(|A|, |C| \le 20\)):
> 1. the class of \(M_0\): the Schmidt arrangement \(\mathcal{S}\) (\(B \equiv i \bmod 2\); horizontal lines, even curvatures);
> 2. the class of \(\begin{pmatrix} 0 & 1 \\ 1 & 0\end{pmatrix}\): the rotated family \(i\mathcal{S}\) (\(B \equiv 1 \bmod 2\); vertical lines);
> 3. the class of \(\begin{pmatrix} 1 & 0 \\ 0 & -1\end{pmatrix}\) (the unit circle): all circles of **odd** curvature — a family that is *not* part of the Schmidt arrangement at all.
>
> Orientation-reversed seeds land in the same three classes.

So "the Schmidt arrangement" is exactly one twisted Galois class of real structures on \(\mathbb{P}^1\); the involution \(\sigma\) is the algebraic shadow of this descent picture.

## 4. Main theorem: \(\sigma\) acts on class groups by a twisted inversion

\(\sigma\) maps \(\Gamma X \Gamma \mapsto \Gamma\, \sigma(X)\, \Gamma\) (well-defined because \(\sigma(\Gamma) = \Gamma\)), hence induces an involution \(\hat\sigma\) on \(\Gamma\)-classes of oriented Schmidt circles. By [hyperbolic-counting.md](hyperbolic-counting.md), classes of circles with \(|\alpha| = n\) (\(n\) odd — only these lie in the \(\mathrm{SL}_2(\mathbb{Z}[i])\)-orbit of \(\hat{\mathbb{R}}\)) correspond to classes of positive definite integral binary quadratic forms of discriminant \(1 - n^2\), via \(\omega \leftrightarrow f_\omega = (q, -x, m)\). This is the involution of quadratic forms produced by \(\sigma\). What it does:

> **Theorem** (proved in [class-formula-proof.md](class-formula-proof.md); every step additionally machine-verified for every class of every odd \(n \le 41\)).
> 1. \(\hat\sigma\) preserves \(\alpha\) exactly (proof: \(Y_{\sigma X} = \bar X^{-1} Y_X \bar X\) is conjugate to \(Y_X\), and \(\operatorname{tr} Y = -2\alpha_\pm\)); it reverses orientation, so it sends upper-half-plane circle classes to lower ones.
> 2. Compose with the reflection \(z \mapsto \bar z\) (which returns to the upper half-plane and preserves each form class). On **primitive** form classes of discriminant \(1 - n^2\) the resulting involution is
> $$
> \hat\sigma\colon\; [f] \;\longmapsto\; [\mathfrak{r}_n] \cdot [f]^{-1},
> \qquad
> \mathfrak{r}_n = \Bigl[\Bigl(\tfrac{n-1}{2},\; 0,\; \tfrac{n+1}{2}\Bigr)\Bigr],
> $$
> composition and inversion in the form class group.
> 3. On the imprimitive classes of content \(g\), \(\hat\sigma\) is likewise \(t_0 \cdot [f_0]^{-1}\) on the underlying primitive classes of discriminant \((1-n^2)/g^2\), where the twist \(t_0\) is an ambiguous class — in all computed cases the natural image of \(\mathfrak{r}_n\) with the \(g^2\)-part of its norm removed (e.g. \(n = 15\), \(g = 2\): \(t_0 = [(2,0,7)]\)); consistency across each stratum verified for all \(n \le 41\).

**Discussion of \(\mathfrak{r}_n\).** The twist is a canonical *ambiguous* (2-torsion) class: as an ideal, \(\mathfrak{r} = \bigl(\tfrac{n-1}{2},\, \omega\bigr)\) with \(\omega = \sqrt{(1-n^2)/4}\), and
$$
\mathfrak{r}^2 = \Bigl(\tfrac{n-1}{2}\Bigr), \qquad \mathfrak{r}\,\mathfrak{s} = (\omega) \ \text{ with } \ \mathfrak{s} = \Bigl(\tfrac{n+1}{2},\, \omega\Bigr),\ [\mathfrak{s}] = [\mathfrak{r}].
$$
Since \([f]^{-1} = [\bar{\mathfrak{a}}]\) for the ideal class \([\mathfrak{a}]\) of \(f\), the theorem reads
$$
\hat\sigma\colon\ [\mathfrak{a}] \longmapsto [\mathfrak{r}_n\, \bar{\mathfrak{a}}] :
$$
**Galois conjugation of the class group, twisted by the distinguished ambiguous class \(\mathfrak{r}_n\).** Geometrically, \(\mathfrak{r}_n\) is the class of the circles *tangent to the lines \(\operatorname{Im} z = k\) of the arrangement*: its representative circle has curvature \(n - 1\), center \(\tfrac{ni}{n-1}\), and touches \(\operatorname{Im} z = 1\).

The decisive numerical test is \(n = 11\) (discriminant \(-120\), class group \((\mathbb{Z}/2)^2\), three distinct nontrivial ambiguous classes): \(\hat\sigma\) sends the principal class to \([(5,0,6)] = \mathfrak{r}_{11}\) — not to the norm-2 class \([(2,0,15)]\) nor \([(3,0,10)]\) — and swaps the other two, exactly as \([\mathfrak{r}][f]^{-1}\) predicts.

**Corollaries.**
- \(\hat\sigma\) has a fixed class iff \(\mathfrak{r}_n\) is a square in the class group (iff \(\mathfrak{r}_n\) lies in the principal genus); the fixed classes then form a coset of the 2-torsion. For \(n = 11\) there are none: \(\hat\sigma\) is a fixed-point-free involution of the four classes. For \(n = 9\) (class group \(\mathbb{Z}/4\), \(\mathfrak{r}_9 = [(4,0,5)]\) = the square of the generator) the two order-4 classes \([(3,\pm2,7)]\) are fixed.
- Since \(\hat\sigma\) preserves each ordinary \(\mathrm{SL}_2(\mathbb{Z}[i])\)-conjugacy class of Cartan images \(Y\), the pair \(\{[f],\, [\mathfrak{r}_n][f]^{-1}\}\) is an invariant of the *closed geodesic* attached to the circle in §2: the circle-to-geodesic map identifies exactly the \(\hat\sigma\)-orbits.

## 5. Does it have a name? Summary of identifications

| Aspect | Name |
|---|---|
| \(\sigma(X) = \bar X^{-1} = M_0^{-1}X^\dagger M_0\) on the algebra \(M_2(\mathbb{Q}(i))\) | the **unitary involution** (involution of the **second kind**) adjoint to the Hermitian form \(M_0\) |
| \(\theta(X) = \bar X\), fixed group \(\mathrm{SL}_2(\mathbb{Z}) = \mathrm{SU}(M_0)\) | the **Galois involution / real structure**; \((\mathrm{SL}_2(\mathbb{C}), \mathrm{SL}_2(\mathbb{R}))\) as a **symmetric pair** |
| \(X \mapsto X\bar X^{-1}\) | the **Cartan embedding** of the symmetric space; image = twisted involutions (Galois cocycles); geometrically, **circle inversions** |
| \(\hat\sigma\) on \(\Gamma\backslash G/\Gamma\) | the **Gelfand-trick anti-involution** of the symmetric pair (the standard device proving \((\mathrm{SL}_2(\mathbb{C}), \mathrm{SL}_2(\mathbb{R}))\) is a Gelfand pair) |
| \(\hat\sigma\) on form classes | Galois conjugation **twisted by the ambiguous class \(\mathfrak{r}_n\)** — formally identical to the **Atkin–Lehner/Fricke action on Heegner points** (\(w_N: [\mathfrak{a}] \mapsto [\mathfrak{n}\bar{\mathfrak{a}}]\), as in Gross's *Heegner points on \(X_0(N)\)*), with \(\mathfrak{r}_n\) playing the role of \(\mathfrak{n}\) |

The last row is worth emphasizing: at the archimedean/group level \(\sigma\) is exactly the involution that makes the pair "Gelfand" (all double cosets self-inverse *over \(\mathbb{R}\) or \(\mathbb{C}\)*); arithmetically, on \(\mathrm{SL}_2(\mathbb{Z})\backslash \mathrm{SL}_2(\mathbb{Z}[i])/\mathrm{SL}_2(\mathbb{Z})\), it acts **nontrivially**, and the failure is measured precisely by class groups — the obstruction being the twisted inversion \([\mathfrak{r}_n][f]^{-1}\). The involution "sees" genus theory: it is trivial on classes exactly when every class is ambiguous and \(\mathfrak{r}_n\) is principal.

## 6. Verification

- [scripts/involution_experiments.py](scripts/involution_experiments.py): exact checks of \(\sigma(X) = J^{-1}X^\dagger J\), the Cartan/fun-fact identity \(X\bar X^{-1} = -iJ\overline{M_X}\), \(M_{\sigma X} = -\overline{M_{X^{-1}}}\), the fixed-point characterization, random-sample class statistics, and the three-class classification of det \(-1\) Hermitian forms.
- [scripts/involution_classmap.py](scripts/involution_classmap.py): for every class of every odd \(n \le 41\), builds an explicit \(X\) realizing the class (by the descent algorithm), computes \(\hat\sigma\), and confirms \([\mathfrak{r}_n][f]^{-1}\) via Gauss composition (concordant forms); also confirms \(\alpha\)-preservation, the uniform orientation reversal, and the stratum-consistent ambiguous twists on imprimitive classes.
- [scripts/proof_check.py](scripts/proof_check.py): machine verification of every lemma of the proof in [class-formula-proof.md](class-formula-proof.md) (169 primitive classes, odd \(n \le 41\)).

## 7. Outlook

1. **The class formula is now proved** ([class-formula-proof.md](class-formula-proof.md)): an explicit unitary basis reduces \(\hat\sigma\) to the ideal identity \(\mathfrak{c} = \mathfrak{s}\,\mathfrak{a}_f\) with an orientation flip, exactly the anticipated \(\mathfrak{a} \mapsto \mathfrak{r}\bar{\mathfrak{a}}\) shape. Remaining at this node: a clean statement for the imprimitive strata at the ramified prime \(2\) (where \(\mathfrak{a}_f\) is not proper and \(\mathfrak{s}\mathfrak{a}_f \subsetneq \mathfrak{c}\) can occur).
2. **Fixed points vs. geometry.** \(\hat\sigma\)-fixed classes (\([f]^2 = [\mathfrak{r}_n]\)) should correspond to circles with a distinguished extra symmetry (a \(\sigma\)-fixed representative \(X = M_0M\)); combined with the ideal-triangle count \(3H(n^2-1)\), this suggests refined class-number-relation geometry (ambiguous classes \(\leftrightarrow\) edge circles was already visible in the \(D_3\)-symmetric pictures).
3. **Circles \(\leftrightarrow\) geodesics.** The trace identity makes \(X \mapsto X\bar X^{-1}\) a bridge from Euclidean circle counting (\(N_e\), Catalan) and hyperbolic-plane counting (\(3H(n^2-1)\)) to **real-quadratic closed geodesics in the Bianchi orbifold** of length \(2\operatorname{arccosh} n\), with the \(\hat\sigma\)-pairs as fibers. Comparing the three counting theories along this bridge is a promising route to genuinely new results.
