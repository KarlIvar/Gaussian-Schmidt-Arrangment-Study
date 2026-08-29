# The atomic census: double cosets of \(\Omega\), Apollonian addresses of form classes, and the \(3/\pi\) split

This document connects the two structural threads of the project: the **atomic structure** of the monoid \(\Omega = \{X \in \mathrm{SL}_2(\mathbb{Z}[i]) : X(\mathbb{H}) \subseteq \mathbb{H}\}\) ([half-plane-monoid.md](half-plane-monoid.md)) and the **invariant theory of double cosets** — the six-dimensional coordinate system \((\alpha, \beta_1, \beta_2, u)\) of [moduli-invariants.md](moduli-invariants.md). Throughout, \(\Gamma = \mathrm{SL}_2(\mathbb{Z})\) (the unit group \(\Omega^\times\)), atoms and unique factorisation are as in [half-plane-monoid.md](half-plane-monoid.md), and the level of a circle/disk is \(\alpha = \operatorname{Im}\zeta\).

The starting observation (the prompt for this note): atoms are permuted by left and right multiplication by \(\Gamma\), each atomic factor of \(X \in \Omega\) is well defined up to its \(\Gamma\)-double coset, and rigidity should therefore decompose the double coset space \(\Gamma\backslash\Omega/\Gamma\) as a *disjoint* union of product cells \(\Gamma A_1\Gamma A_2 \Gamma \cdots A_k\Gamma\). This is proved below (§1) and then made quantitative: the cells stratify the class numbers \(h_+(1-n^2)\), the strata have closed-form counts on the Ford side, and the stratification has density constants \(3/\pi\) and \(1 - 3/\pi\).

Everything with finite content is machine-verified in exact arithmetic: [scripts/atomic_census.py](scripts/atomic_census.py), **4414 checks, all passing** (all odd levels \(n \le 41\); random product tests; the density experiment to \(X = 1601\)).

## 0. Summary

| | |
|---|---|
| **Dictionary** | \(\Gamma\backslash\Omega/\Gamma \;=\; \{\Gamma\} \sqcup \mathbb{Z}_{\ge 1} \sqcup \bigsqcup_{n \ge 3 \text{ odd}} \mathrm{Cl}_+(1-n^2)\): the double cosets at level \(n\) are exactly the classes of positive definite forms of discriminant \(1-n^2\), imprimitive included (§1) |
| **Cells** | \(\Omega = \bigsqcup_{k \ge 0}\ \bigsqcup_{([A_1],\dots,[A_k])} \Gamma A_1 \Gamma A_2 \Gamma \cdots A_k \Gamma\), disjoint over words in atom classes: every form class of discriminant \(1-n^2\) acquires a canonical **Apollonian address** (Thm. 1) |
| **Ends of the word** | the first letter is a function of the \(\beta_1\)-circle \(\omega_1(X)\) alone, the last letter of the \(\beta_2\)-circle \(\omega_2(X)\) alone; \((\alpha,\beta_1,\beta_2,u)\) all collapse to functions of \((n, [f])\) on the discrete stratum (§2) |
| **Level of a product** | \(\alpha(XY) = \langle M_{X^{-1}}, M_Y\rangle\): the second circle of the left factor pairs with the first circle of the right factor (Prop. 4) |
| **Superadditivity** | \(\varepsilon_{\alpha(XY)} \ge \varepsilon_{\alpha(X)}\,\varepsilon_{\alpha(Y)}\), \(\varepsilon_n = n + \sqrt{n^2-1}\): geodesic length is superadditive along factorisation; words per level are finite (Thm. 5) |
| **Pure twist** | in the inner-disk normalisation the involution \(\sigma\) acts on level-\(n\) classes as \([f] \mapsto [\mathfrak{r}_n]\cdot[f]\) — *no inversion* (Thm. 6); it reverses addresses (Thm. 7) |
| **Ford census** | \(\#\{\text{level-}n\text{ classes in a Ford horoball}\} = \sum_{q \le (n-1)/2} \#\{\xi \bmod q : \xi^2 \equiv \tfrac{1-n^2}{4}\}\) (Thm. 8); \(\#\{\text{word } T\,T\}= \varphi(c)\) at \(n = 2c^2+1\) (Thm. 9); second letters are representation numbers (Prop. 10) |
| **Extremal depth** | \(1 \le \ell \le \frac{n+1}{2}\), the maximum attained exactly by the **principal class** \(T_i^{(n-1)/2}\,ST_i\) and the **twist class** \(\mathfrak{r}_n = T_i\,ST_i^{(n-1)/2}\), which are \(\sigma\)-images of each other (Thm. 11) |
| **The \(3/\pi\) split** | \(\sum_{n \le X \text{ odd}} b^{T}(n) = \frac{X^2}{8G} + o(X^2)\); against the census \(\sim \frac{\pi}{24G}X^2\), the Ford stratum has density \(3/\pi = 0.9549\ldots\) and the deep-gasket stratum density \(1 - 3/\pi = 0.0451\ldots\) (§6) |

## 1. The cell decomposition and the address map

Recall from [half-plane-monoid.md](half-plane-monoid.md): \(X \mapsto D(X) = X(\mathbb{H})\) is a bijection \(\Omega/\Gamma \to \mathcal{D}\) (Schmidt disks in \(\mathbb{H}\)), left divisibility is reverse inclusion, and every non-unit factors as \(X = A_1\cdots A_kU\) (atoms \(A_j\), unit \(U\)), uniquely up to associates: \(A_j' = V_{j-1}^{-1}A_jV_j\) (Thm. 8 there). Consequently \(\Gamma\backslash\Omega/\Gamma \cong \Gamma\backslash\mathcal{D}\), and by the form dictionary of [hyperbolic-counting.md](hyperbolic-counting.md) §2:

> **Proposition 0.** The double cosets of \(\Omega\) are: the unit coset \(\Gamma\); the level-\(1\) classes \(\{\operatorname{Im}z > t\}\), \(t \in \mathbb{Z}_{\ge1}\) (every level-\(1\) disk, bounded or not, is \(\Gamma\)-equivalent to exactly one of these); and, for each odd \(n \ge 3\), the classes of level-\(n\) disks, in bijection with the \(\mathrm{SL}_2(\mathbb{Z})\)-classes of positive definite forms of discriminant \(1 - n^2\), **imprimitive forms included**. Write \(h_+(1-n^2)\) for their number and \([f] \leftrightarrow\) class of the disk with data \(f_D = (q, -x, m)\).

(The level-\(1\) fiber is the degenerate discriminant \(0\): its "class set" \(\mathbb{Z}_{\ge1}\) is the classical cusp stratum, with \(\{\operatorname{Im}z>t\}\) of depth \(t\) and word \(T^t\); this matches the divergence of class numbers at discriminant \(0\).)

> **Theorem 1 (cells).** Let \(\mathcal{A}\mathrm{t}\) be the set of associate classes of atoms (\(=\Gamma\)-double cosets of atoms \(=\Gamma\)-orbits of gasket disks). Then
> $$
> \Omega \;=\; \bigsqcup_{k \ge 0}\ \bigsqcup_{w = ([A_1],\dots,[A_k]) \in \mathcal{A}\mathrm{t}^k} \Gamma A_1 \Gamma A_2 \Gamma \cdots A_k \Gamma ,
> $$
> a **disjoint** union over all finite words \(w\) (the empty word giving \(\Gamma\)). Each cell is a union of double cosets, so every double coset — equivalently every level-\(n\) form class — carries a well-defined word
> $$
> \mathcal{A}\colon\ \mathrm{Cl}_+(1-n^2) \longrightarrow \bigsqcup_k \mathcal{A}\mathrm{t}^k, \qquad [f] \mapsto w([f]),
> $$
> its **Apollonian address**; the length of the word is the depth \(\ell\).

*Proof.* Every element of \(\Gamma A_1\Gamma\cdots A_k\Gamma\) has the form \(V_0A_1V_1A_2\cdots A_kV_k = (V_0A_1)(V_1A_2)\cdots(V_{k-1}A_k)V_k\), a product of \(k\) atoms (a unit times an atom is an atom) and a unit, with \(i\)-th factor in the class \([A_i]\); conversely any \(X\) with factorisation of that letter sequence lies in the cell. Disjointness is exactly the uniqueness clause of the factorisation theorem: two factorisations of \(X\) have associate letters, hence the same word. \(\square\)

> **Proposition 2 (contents of a cell).** For atoms \(A, B\), the double cosets inside \(\Gamma A\Gamma B\Gamma\) are parametrised by
> $$
> \Gamma^{A}\backslash \Gamma / \Gamma_{B}, \qquad \Gamma^A := A^{-1}\Gamma A \cap \Gamma, \quad \Gamma_B := B\Gamma B^{-1} \cap \Gamma ,
> $$
> via \(\gamma \mapsto \Gamma A\gamma B\Gamma\); for longer words, double cosets in the cell of \(w\) are the \(\Gamma\)-orbits of nested chains \(\mathbb{H} \supsetneq G_1 \supsetneq \cdots \supsetneq G_k\) in which each relative step \(G_{i-1} \leadsto G_i\) is a gasket disk of class \([A_i]\).

*Proof.* \(\Gamma_B = \operatorname{Stab}_\Gamma(B(\mathbb{H}))\), so \(A\gamma B(\mathbb{H})\) depends exactly on \(\gamma\Gamma_B\); the residual left action of \(\operatorname{Stab}_\Gamma(A(\mathbb{H})) = A\Gamma^AA^{-1}\) transports to \(\Gamma^A\) acting on the left of \(\gamma\). The chain description is Lemma 4 of [half-plane-monoid.md](half-plane-monoid.md) applied inductively. \(\square\)

The groups \(\Gamma^A, \Gamma_A\) are **thin**: for the Ford atom \(T_i = \begin{pmatrix}1 & i\\ 0 & 1\end{pmatrix}\) one computes \(\Gamma^{T_i} = \Gamma_{T_i} = \pm\langle T\rangle\), of infinite index in \(\Gamma\). This is precisely why \((\Gamma, \mathrm{SL}_2(\mathbb{Z}[i]))\) carries no classical Hecke algebra (double cosets are infinite unions of one-sided cosets, gap **G7** of [spectral-geometry.md](spectral-geometry.md)) — and Theorem 1 is the substitute: a *relative Bruhat decomposition with thin stabilisers*, in which the product of two cells decomposes by the pairing law of §3.

## 2. Where the six invariants sit

On \(\mathrm{SL}_2(\mathbb{C})\) the double coset space \(\Gamma\backslash G/\Gamma\) needs all six coordinates \((\alpha, \beta_1, \beta_2, \arg\Theta)\). On the discrete stratum this collapses: a double coset of \(\Omega\) *is* a pair \((n, [f])\) by Proposition 0, so \(\beta_1 = j(m_1)\) is the singular modulus of \([f]\), \(\beta_2\) that of the twisted class, and \(u = u_f\) the (now algebraic, [first-power-descent.md](first-power-descent.md)) phase — all **functions of \((n,[f])\)**. The address \(\mathcal{A}([f])\) is a further derived invariant. What the continuous coordinates *do* see is the two ends of the word:

* the **first letter** is the class of the unique maximal gasket disk containing \(D(X)\) — a function of the circle \(\omega_1(X) = \partial D(X)\) alone, i.e. of the \(\beta_1\)-side;
* the **last letter** is the first letter of \(\sigma(X)\) (§4), i.e. a function of the \(\sigma\)-circle \(\omega_2(X)\) alone — the \(\beta_2\)-side.

This is the precise form of the observation that "the first and last atoms are slightly more well-defined": each end of the word is visible to one member of the circle pair \((\omega_1, \omega_2)\); the interior letters are visible to neither circle separately.

**The phase is a Birkhoff product over the address.** Let \(X = A_1\cdots A_kU\) and let \(m_2 = m_2(X)\) be the lower CM point, with branch \(X(m_2) = m_1\) (Lemma B of [moduli-invariants.md](moduli-invariants.md) §5). The chain rule gives
$$
X'(m_2) \;=\; \prod_{i=1}^{k} A_i'(w_i)\ \cdot\ U'(m_2),
\qquad w_i := (A_{i+1}\cdots A_kU)(m_2),
$$
so the sixth invariant \(\Theta(X) = j'(m_1)X'(m_2)/\overline{j'(\bar m_2)}\) is, up to the kernel ratio at the two CM points, the **Birkhoff product of atomic multipliers along the backward orbit of the CM point** through the gasket dynamical system. Combined with Theorem 3.3 of [first-power-descent.md](first-power-descent.md) (\(u_f = \Phi_y/\Phi_x(\beta_1,\beta_2)\)): *the derivative of the modular correspondence at the Heegner pair factors as a product of atomic derivatives along the Apollonian address* — the arithmetic invariant and the transfer-operator cocycle of [spectral-outlook.md](spectral-outlook.md) §7 are the same object.

## 3. The level of a product, and superadditivity

Write \(\langle M, M'\rangle = xx' + yy' - 2(qm' + q'm)\) for the (negated) inversive product ([half-plane-monoid.md](half-plane-monoid.md) §1); recall \(\langle M, M\rangle = 1\) — every Schmidt circle is a *unit vector* on the de Sitter quadric \(\det = -1\) — and \(\alpha(X) = \langle M_0, M_X\rangle\).

> **Proposition 4 (level pairing).** For all \(X, Y \in \mathrm{SL}_2(\mathbb{Z}[i])\),
> $$
> \boxed{\;\alpha(XY) \;=\; \bigl\langle M_{X^{-1}},\, M_{Y}\bigr\rangle\;}
> $$
> — the level of a product is the inversive pairing of the circle of \(X^{-1}\) (the reflected \(\sigma\)-circle \(\bar\omega_2(X)\), i.e. the \(\beta_2\)-side of \(X\)) with the circle of \(Y\) (the \(\beta_1\)-side of \(Y\)).

*Proof.* \(M_{XY} = ((XY)^{-1})^\dagger M_0 (XY)^{-1} = (X^{-1})^\dagger M_Y X^{-1}\), and \(\langle M, M'\rangle = -\tfrac12\operatorname{tr}(M\operatorname{adj}M')\) is invariant under the simultaneous congruence \(M \mapsto g^\dagger Mg\), \(g \in \mathrm{SL}_2(\mathbb{C})\). Hence \(\alpha(XY) = \langle M_0, (X^{-1})^\dagger M_YX^{-1}\rangle = \langle X^\dagger M_0X, M_Y\rangle = \langle M_{X^{-1}}, M_Y\rangle\). \(\square\)

(Verified on random products; note the asymmetry — the *naive* guess \(\langle M_X, M_{Y^{-1}}\rangle\) computes \(\alpha(YX)\), and the failed check that exposed this is retained in the script's history. Composition couples the **second** circle of the left factor to the **first** circle of the right factor, exactly as the \((\beta_1,\beta_2)\)-splitting of §2 predicts.) As \(\gamma\) ranges over \(\Gamma\), \(\alpha(X\gamma Y) = \langle M_{X^{-1}}, \gamma\cdot M_Y\rangle\) is a two-orbit pairing representation problem: the double cosets in a product of cells are configurations of two circles with prescribed classes and prescribed pairing — Gram data \(\begin{pmatrix}1 & a & b\\ a & 1 & n\\ b & n & 1\end{pmatrix}\) on the triple \((M_0, M_{X^{-1}}, \gamma M_Y)\).

> **Theorem 5 (superadditivity).** For \(X, Y \in \Omega\), with \(\varepsilon_n := n + \sqrt{n^2-1} = \exp(\operatorname{arccosh} n)\):
> $$
> \varepsilon_{\alpha(XY)} \;\ge\; \varepsilon_{\alpha(X)}\,\varepsilon_{\alpha(Y)},
> \qquad\text{equivalently}\qquad
> \operatorname{arccosh}\alpha(XY) \ \ge\ \operatorname{arccosh}\alpha(X) + \operatorname{arccosh}\alpha(Y).
> $$
> Consequently, along the address of any class: \(\varepsilon_n \ge \prod_i \varepsilon_{\alpha(A_i)}\); in particular every letter has level \(\le n\) and, since each level carries finitely many atom classes, **only finitely many words occur at each level**.

*Proof.* \(D(XY) = X(D(Y)) \subseteq D(X) \subseteq \mathbb{H}\) are nested, so the geodesic plane \(P_1\) over \(\partial D(X)\) weakly separates the plane \(P_0\) over \(\hat{\mathbb{R}}\) from the plane \(P_2\) over \(\partial D(XY)\) in \(\mathbb{H}^3\). For nested disks the pairing is \(\cosh\) of the plane distance (\(=1\) at internal tangency), so \(d(P_0,P_2) \ge d(P_0,P_1) + d(P_1,P_2)\) (cross \(P_1\) along any geodesic realising \(d(P_0,P_2)\)). Finally \(\langle M_X, M_{XY}\rangle = \langle M_0, M_Y\rangle = \alpha(Y)\) by the same congruence invariance, so the three distances are \(\operatorname{arccosh}\) of \(\alpha(XY), \alpha(X), \alpha(Y)\). \(\square\)

Equality holds iff the three planes share a common perpendicular (the tangency/parabolic cases, \(\varepsilon = 1\), being degenerate instances). In the notation of [involution.md](involution.md): the closed geodesic attached to a circle has length \(2\log\varepsilon_n\) — **geodesic length is superadditive along atomic factorisation**, and the deficiency \(\log\varepsilon_n - \sum_i\log\varepsilon_{a_i}\) measures the bending of the chain. Note the recurring identity \(\varepsilon_{2q^2-1} = \varepsilon_q^2\) (e.g. \(\varepsilon_7 = \varepsilon_2^2\)): the atoms tangent to \(\operatorname{Im}z=1\) sit at the "perfectly aligned" levels.

## 4. The involution: pure twist and word reversal

> **Theorem 6 (pure twist).** For \(X \in \Omega\) at level \(n\) with disk class \([f]\), the disk class of \(\sigma(X) = \bar X^{-1} \in \Omega\) is
> $$
> \boxed{\;[f_{D(\sigma X)}] \;=\; [\mathfrak{r}_n]\cdot[f]\;}
> $$
> (composition with the twist class, **no inversion**), for every primitive class — machine-verified for every primitive class of every odd \(n \le 41\) by exact ideal composition \([\mathfrak{r}_n\mathfrak{a}_f]\), and consistent (stratum-adjusted twist) on imprimitive classes.

This is the *inner-disk normalisation* of the class formula \(\hat\sigma[f] = [\mathfrak{r}_n][f]^{-1}\) of [class-formula-proof.md](class-formula-proof.md): there the circle of \(X\) is taken with the positive-curvature (exterior) orientation and the image is reflected by \(z \mapsto \bar z\); passing between the exterior and interior conventions is a mirror, which inverts classes ([circle-composition.md](circle-composition.md) §2), and the two formulas differ by exactly that inversion. The decisive test separating them is \(n = 9\) (\(\mathrm{Cl}(-80) = \mathbb{Z}/4\), generator \(g = [(3,2,7)]\), \(\mathfrak{r}_9 = g^2\)): the pure twist **swaps** \((3,2,7) \leftrightarrow (3,-2,7)\), whereas \([\mathfrak{r}][f]^{-1}\) would fix them; the swap is what the monoid computes. Two pleasant consequences of the twist form: on \(\Omega\) the involution is a *translation* of the class group — fixed-point-free whenever \(\mathfrak{r}_n\) is non-principal — and the primitive level-\(n\) classes fall into canonical **\(\mathfrak{r}_n\)-twin pairs** \(\{[f], [\mathfrak{r}_n f]\}\), the pairs already carrying the phase relation \(u_{\mathfrak{r}f}\,u_f = 1\) of [moduli-invariants.md](moduli-invariants.md). (On imprimitive strata the stratum-adjusted twist can have fixed classes — e.g. \((4,4,4)\) at \(n=7\), the lone word \(A_7\), or the lone \(T^8\)-class at \(n = 33\) — which is why odd word multiplicities occur only at palindromic words on imprimitive strata.)

> **Theorem 7 (word reversal).** \(\sigma\) reverses addresses:
> $$
> \mathcal{A}\bigl([\mathfrak{r}_n][f]\bigr) \;=\; \text{reverse of } \bigl(\hat\sigma[A_1], \dots, \hat\sigma[A_k]\bigr), \qquad \mathcal{A}([f]) = ([A_1],\dots,[A_k]),
> $$
> where \(\hat\sigma\) on atom classes is the corresponding twist at the atom's own level (trivial on the Ford class and on every atom level \(\le 41\), where the atom class is unique). Hence \(N_w(n) = N_{\sigma(w)^{\mathrm{rev}}}(n)\) for all words, and twin pairs have mutually reversed addresses.

*Proof.* \(\sigma(A_1\cdots A_kU) = \sigma(U)\sigma(A_k)\cdots\sigma(A_1)\) (Prop. 10 of [half-plane-monoid.md](half-plane-monoid.md)), which is an atomic factorisation with letters \([\sigma(A_{k+1-i})]\). \(\square\)

## 5. The census identities

Fix odd \(n \ge 3\) and write \(a(n)\) = number of atom classes at level \(n\) (nonzero exactly on the Apollonian \(\alpha\)-spectrum \(1, 7, 17, 31, 49, \dots\)), \(b^{T}(n)\) = number of classes whose first letter is the Ford class \([T_i]\), and \(R(n) = h_+(1-n^2) - a(n) - b^T(n)\) — the **deep classes**, whose disk lies inside a *bounded* gasket disk. The master decomposition is
$$
h_+(1-n^2) \;=\; a(n) \;+\; b^{T}(n)\;+\;R(n), \qquad\text{refined by } \ h_+(1-n^2) = \sum_w N_w(n).
$$

> **Theorem 8 (Ford stratum, closed form).** A class has first letter \([T_i]\) iff its disk lies in a Ford horoball (a Ford disk or a translate of \(\{\operatorname{Im}z>1\}\)), and
> $$
> b^{T}(n) \;=\; \sum_{q=1}^{(n-1)/2} \#\Bigl\{x \bmod 2q\ :\ x \text{ even},\ x^2 + n^2 \equiv 1 \ (\mathrm{mod}\ 4q)\Bigr\}
> \;=\; \sum_{q=1}^{(n-1)/2} \#\Bigl\{\xi \bmod q : \xi^2 \equiv \tfrac{1-n^2}{4}\ (\mathrm{mod}\ q)\Bigr\}.
> $$

*Proof.* First letter \([T_i]\) means the maximal disk \(G_D\) lies in the orbit \(\Gamma\{\operatorname{Im}z>1\}\); normalising \(G_D = \{\operatorname{Im}z>1\}\), two normalised disks are \(\Gamma\)-equivalent iff they differ by \(\operatorname{Stab}_\Gamma(\{\operatorname{Im}z>1\}) = \pm\langle T\rangle\). (Conversely if \(D \subsetneq \{\operatorname{Im}z>1\}\), laminarity forces \(G_D = \{\operatorname{Im}z>1\}\).) A level-\(n\) disk of curvature \(2q\) lies in \(\{\operatorname{Im}z>1\}\) iff its lowest point has height \(\frac{n-1}{2q} \ge 1\), i.e. \(q \le \frac{n-1}{2}\); modulo \(\pm\langle T\rangle\) (translations \(x \mapsto x + 2q\); \(-I\) acts trivially) the disks are the residues \(x \bmod 2q\) subject to the classification congruence. Substituting \(x = 2\xi\) gives the second form. \(\square\)

> **Theorem 9 (Ford–Ford cell).** The classes with word exactly \(([T_i],[T_i])\) exist only at the levels \(n = 2c^2+1\), where their number is **Euler's \(\varphi(c)\)**.

*Proof.* By Proposition 2 with \(\Gamma^{T_i} = \Gamma_{T_i} = \pm\langle T\rangle\): double cosets in \(\Gamma T_i\Gamma T_i\Gamma\) are \(\pm U\backslash\Gamma/\pm U\), i.e. (for lower-left entry \(c \ne 0\), normalised \(c>0\)) the pairs \((c, a \bmod c)\) with \(a \in (\mathbb{Z}/c)^\times\): exactly \(\varphi(c)\) cosets for each \(c \ge 1\) (and one coset with \(c = 0\), the level-\(1\) class \(\{\operatorname{Im}z>2\}\)). The level: \(\gamma\{\operatorname{Im}z>1\}\) is the horoball at \(a/c\) with data \((q,x,y,m) = (c^2, 2ac, 1, a^2)\), and by Proposition 4, \(\alpha(T_i\gamma T_i) = \langle M_{T_i^{-1}}, M_{\gamma T_i}\rangle = 1 + 2c^2\). \(\square\)

> **Proposition 10 (second letters are representation numbers).** For an atom class \([A]\) at level \(a\) and \(n > a\) odd, the number of classes with word exactly \(([T_i],[A])\) is
> $$
> N_{(T,[A])}(n) \;=\; \#\Bigl\{x \bmod 2q'\ :\ \text{the level-}a\text{ disk } (q',x) \text{ exists and lies in } [A]\Bigr\},
> \qquad q' = \tfrac{n-a}{2},
> $$
> the number of times the atom class **represents the half-curvature \(q'\)** (via [circle-composition.md](circle-composition.md) §4.3: the curvature spectrum of a class is its set of represented numbers). Theorem 9 is the case \([A] = [T_i]\): the Ford class represents exactly the half-curvatures \(c^2\), \(\varphi(c)\) times.

*Proof.* Prepending \(T_i\) is the bijection {classes with first letter \(T\), rest of word \(w\)} \(\leftrightarrow\) {\(\pm U\)-orbits of disks \(D'\) with address \(w\) and \(y' + 2q' = n\)}, since \(\alpha(T_iD') = \langle M_{T_i^{-1}}, M_{D'}\rangle = y' + 2q'\). For \(w = ([A])\) the disks \(D'\) are the atom disks of class \([A]\) (level \(y' = a\)), parametrised per period by \(x \bmod 2q'\). \(\square\)

Worked instance (\(n = 31\), \([A_7] = [(4,4,4)]\)): \(q' = 12\), and the only \(x \bmod 24\) with \((12, -x, m)\) in the class \((4,4,4)\) is \(x = 12\) — so exactly one class with word \(T\,A_7\), namely \((12,12,23)\); its twin \([\mathfrak{r}_{31}f]\) is \((17,14,17)\), word \(A_7\,T\). ✓

> **Theorem 11 (extremal depth).** For odd \(n \ge 3\): every level-\(n\) class has \(1 \le \ell \le \frac{n+1}{2}\), and the maximum is attained **exactly** by the principal class and the twist class \(\mathfrak{r}_n = [(\tfrac{n-1}2, 0, \tfrac{n+1}2)]\) (which coincide iff \(n = 3\)), with all-Ford addresses and explicit factorisations
> $$
> X_{\mathrm{princ}} = T_i^{(n-1)/2}\,(ST_i), \qquad X_{\mathfrak{r}} = T_i\,(ST_i)\,T_i^{(n-3)/2}\ \sim\ \sigma(X_{\mathrm{princ}}).
> $$
> The suffix levels descend arithmetically \(n, n-2, n-4, \dots, 3, 1\) for the principal class (the vertical stack of half-planes: the Eisenstein direction) and collapse instantly \(n, 1, 1, \dots, 1\) for the twist class (a chain of tangent horoballs of curvatures \(2, 4, \dots, n-1\)): the involution exchanges "climbing the ladder" with "hugging the cusp".

*Proof.* *Bound.* Proposition 7 of [half-plane-monoid.md](half-plane-monoid.md) gives \(\ell \le q + \lfloor\frac{n-1}{2q}\rfloor\) for any representative of half-curvature \(q\); taking a Gauss-reduced representative, \(q \le \sqrt{(n^2-1)/3}\), and on \([1, \sqrt{(n^2-1)/3}]\) the bound \(q + \frac{n-1}{2q}\) is maximised at \(q=1\), giving \(\frac{n+1}{2}\).
*Attainment.* \(T_i^{(n-1)/2}ST_i(\mathbb{H})\) is the curvature-\(2\) disk \(|z - \tfrac{n i}{2}| < \tfrac12\) (form \((1,0,\tfrac{n^2-1}4)\), principal), of depth \(\frac{n-1}{2}+1\); \(T_iST_i^{(n-1)/2}(\mathbb{H})\) is \(|z - \tfrac{ni}{n-1}| < \tfrac1{n-1}\) (form \((\tfrac{n-1}2,0,\tfrac{n+1}2) = \mathfrak{r}_n\)), same depth; both products are atomic factorisations by rigidity, all letters in \([T_i]\).
*Exactness.* Equality \(q + \frac{n-1}{2q} \ge \frac{n+1}{2}\) on the reduced range forces \(q \in \{1, \frac{n-1}{2}\}\) (roots of \(2q^2 - (n+1)q + (n-1)\)). For \(q = 1\) the only reduced form is principal. For \(q = \frac{n-1}{2}\), a saturating chain has one half-plane step and \(\frac{n-1}{2}\) bounded steps of strictly increasing half-curvatures in \([1, \frac{n-1}{2}]\), forcing curvatures exactly \(2, 4, \dots, n-1\); consecutive nesting of disks with radii \(\frac{1}{2j} \supset \frac1{2j+2}\) forces inversive product \(1\) (an odd integer that cannot exceed \(\frac{(j+1)^2+j^2}{2j(j+1)} < 3\)), i.e. a full tangent chain, whose class is \(\mathfrak{r}_n\). Machine-verified for all \(n \le 41\). \(\square\)

### The census table (all odd \(3 \le n \le 41\); machine-generated, exact)

Words are written left-to-right; \(T^k\) means \(k\) Ford letters; \(A_a\) is the (unique, for \(a \le 41\)) atom class at level \(a\).

| \(n\) | \(h_+\) | \(a(n)\) | \(b^T(n)\) | \(R(n)\) | words (with multiplicities) |
|---|---|---|---|---|---|
| 3 | 1 | 0 | 1 | 0 | \(T^2\):1 |
| 5 | 2 | 0 | 2 | 0 | \(T^3\):2 |
| 7 | 4 | 1 | 3 | 0 | \(A_7\):1, \(T^3\):1, \(T^4\):2 |
| 9 | 6 | 0 | 6 | 0 | \(T^2\):1, \(T^3\):2, \(T^4\):1, \(T^5\):2 |
| 11 | 4 | 0 | 4 | 0 | \(T^4\):2, \(T^6\):2 |
| 13 | 4 | 0 | 4 | 0 | \(T^5\):2, \(T^7\):2 |
| 15 | 12 | 0 | 11 | 1 | \(A_7T\):1, \(TA_7\):1, \(T^3\):2, \(T^4\):4, \(T^5\):2, \(T^8\):2 |
| 17 | 9 | 1 | 8 | 0 | \(A_{17}\):1, \(T^3\):2, \(T^4\):1, \(T^5\):1, \(T^6\):2, \(T^9\):2 |
| 19 | 10 | 0 | 10 | 0 | \(T^2\):2, \(T^4\):3, \(T^6\):3, \(T^{10}\):2 |
| 21 | 12 | 0 | 12 | 0 | \(T^3\):4, \(T^5\):4, \(T^7\):2, \(T^{11}\):2 |
| 23 | 12 | 0 | 11 | 1 | \(A_7T^2\):1, \(T^2A_7\):1, \(T^4\):4, \(T^6\):2, \(T^7\):2, \(T^{12}\):2 |
| 25 | 16 | 0 | 16 | 0 | \(T^3\):2, \(T^4\):6, \(T^5\):2, \(T^7\):2, \(T^8\):2, \(T^{13}\):2 |
| 27 | 12 | 0 | 12 | 0 | \(T^4\):4, \(T^6\):4, \(T^8\):2, \(T^{14}\):2 |
| 29 | 8 | 0 | 8 | 0 | \(T^5\):2, \(T^7\):2, \(T^9\):2, \(T^{15}\):2 |
| 31 | 16 | 1 | 13 | 2 | \(A_{31}\):1, \(A_7T\):1, \(TA_7\):1, \(TA_7T\):1, \(A_7T^3\):1, \(T^3A_7\):1, \(T^5\):3, \(T^7\):1, \(T^8\):2, \(T^9\):2, \(T^{16}\):2 |
| 33 | 28 | 0 | 28 | 0 | \(T^2\):2, \(TA_7T\):2, \(T^4\):7, \(T^5\):8, \(T^7\):4, \(T^8\):1, \(T^{10}\):2, \(T^{17}\):2 |
| 35 | 20 | 0 | 18 | 2 | \(A_{17}T\):2, \(TA_{17}\):2, \(T^4\):6, \(T^6\):4, \(T^8\):2, \(T^{10}\):2, \(T^{18}\):2 |
| 37 | 18 | 0 | 18 | 0 | \(T^3\):8, \(T^5\):4, \(T^9\):2, \(T^{11}\):2, \(T^{19}\):2 |
| 39 | 32 | 0 | 31 | 1 | \(T^3\):2, \(T^4\):8, \(T^5\):10, \(A_7T^4\):1, \(T^4A_7\):1, \(T^8\):6, \(T^{11}\):2, \(T^{20}\):2 |
| 41 | 24 | 0 | 23 | 1 | \(A_{17}T\):1, \(TA_{17}\):1, \(T^3\):4, \(T^5\):6, \(T^6\):4, \(T^9\):4, \(T^{12}\):2, \(T^{21}\):2 |

Every visible regularity is a theorem above: word multiplicities are reversal-symmetric (Thm. 7); \(T^2\) appears exactly at \(n = 3, 9, 19, 33\) (\(=2c^2+1\)) with multiplicities \(1, 1, 2, 2 = \varphi(1),\varphi(2),\varphi(3),\varphi(4)\) (Thm. 9); the longest word \(T^{(n+1)/2}\) always has multiplicity \(2\) — principal + twist (Thm. 11); deep classes (\(R > 0\)) first appear at \(n = 15\), always in \(\hat\sigma\)-twin pairs with a Ford-first partner.

## 6. Mean values: the \(3/\pi\) split

> **Proposition 12 (Ford census on average).** With \(G\) Catalan's constant,
> $$
> B(X) := \sum_{\substack{n \le X \\ n \text{ odd}}} b^{T}(n) \;=\; \frac{X^2}{8G}\ +\ o(X^2).
> $$

*Derivation.* Summing Theorem 8 over odd \(n\), the solutions \((x, n) \bmod 2q\) of \(x^2 + n^2 \equiv 1 \pmod{4q}\) (\(x\) even, \(n\) odd) number exactly \(N_e(q)\), the Euclidean density of [euclidean-counting.md](euclidean-counting.md); each full period of length \(2q\) in \(n\) thus contributes \(N_e(q)\), so
$$
B(X) \approx \sum_{q \le X/2} N_e(q)\,\frac{X - 2q}{2q}
= \frac{X}{2}\sum_{q \le X/2}\frac{N_e(q)}{q} - \sum_{q\le X/2}N_e(q)
= \frac{X^2}{4G} - \frac{X^2}{8G} = \frac{X^2}{8G},
$$
using the exact mean values \(\sum_{q\le Q}N_e(q) = Q^2/2G + O(Q\log Q)\) and (by partial summation) \(\sum_{q \le Q}N_e(q)/q = Q/G + O(\log^2 Q)\). The step requiring care is the equidistribution of the solution set over the last incomplete period in \(n\) — a standard divisor-sum analysis (the count is a Dirichlet-hyperbola problem in \((q, x, n)\)) which we have not carried out rigorously; numerically the error is well below \(X^{2}\): at \(X = 1601\) the formula predicts \(349{,}795\) against the exact \(350{,}258\) (\(0.13\%\)). \(\square\)

Against the level census — the odd-\(n\) restriction of anchor 2 of [spectral-outlook.md](spectral-outlook.md) §0, verified here independently:
$$
\sum_{\substack{n\le X\\ n \text{ odd}}} h_+(1-n^2) \;\sim\; \frac{\pi}{24\,G}\,X^2
\qquad(\text{predicted } 366{,}306 \text{ vs. exact } 366{,}810 \text{ at } X = 1601),
$$
we obtain the headline constants:

> **The \(3/\pi\) split.** The Ford stratum has logarithmic density
> $$
> \frac{1/8G}{\pi/24G} \;=\; \frac{3}{\pi} \;=\; 0.954929\ldots
> $$
> among all level-\(n\) form classes, and the deep-gasket stratum (atoms and classes inside bounded gasket disks) has density \(1 - \tfrac3\pi = 0.045070\ldots\). Measured ratios \(B(X)/\sum h_+\): \(0.95528\) (\(X{=}101\)), \(0.95465\) (\(801\)), \(0.95488\) (\(1601\)).

So, quantitatively: **more than \(95\%\) of all Schmidt circle classes lie in a Ford horoball, and the missing \(4.5\% = 1 - 3/\pi\) is the census of the thin (Apollonian) part of the monoid.** The proportion is universal — independent of \(n\) on average — and its two constants are the Euclidean counting constant \(1/2G\) and the census constant \(\pi/12G\), i.e. ultimately \(\operatorname{vol}(M) = G/3\) and \(\operatorname{area}(Y) = \pi/3\).

## 7. Machine verification

[scripts/atomic_census.py](scripts/atomic_census.py) — **4414 checks, all passing** (`python3 scripts/atomic_census.py 41`; density mode `asym`):

1. for every class of every odd \(n \le 41\): a representative \(X\) is built (via the descent of [circle-classification.md](circle-classification.md)), factored by the algorithm of [half-plane-monoid.md](half-plane-monoid.md) §9, re-multiplied, and its address recorded — the census totals equal \(h_+(1-n^2)\) by construction, and depth-\(1\) classes coincide exactly with maximal disks;
2. Theorem 8 (Ford-first census = congruence count) at every level;
3. Theorem 9 (\(\varphi(c)\) at \(n = 2c^2+1\), zero elsewhere) at every level;
4. Theorem 7: \(\operatorname{word}(\sigma X)\) equals the reversed letterwise-\(\sigma\) word, and the induced class map is an involution, at every class;
5. Theorem 6: the class of \(D(\sigma X)\) equals the exact ideal product \([\mathfrak{r}_n\mathfrak{a}_f]\) (HNF ideal arithmetic, convention of [class-formula-proof.md](class-formula-proof.md) §0) for every **primitive** class, all odd \(n \le 41\); hard-coded decisive cases at \(n = 9, 11\);
6. Theorem 5: the exact inequality \(\varepsilon_{n_i} \ge \varepsilon_{a_i}\varepsilon_{n_{i+1}}\) at every step of every chain (integer arithmetic), plus random products;
7. Proposition 4 on \(400\) random products in \(\Omega\) (both unfoldings);
8. Proposition 10 (second-letter representation kernel) for every occurring letter at every level;
9. Theorem 11 (depth bound and exact attainment set) at every level;
10. the density experiment to \(X = 1601\).

## 8. Research outlook

1. **The atomic transfer operator on class groups.** Theorem 1 organises \(\bigoplus_n \mathbb{Z}[\mathrm{Cl}_+(1-n^2)]\) as a module over an "atomic Hecke" operator \(\mathsf{T}\): \(\mathsf{T}\) sends the class of a disk to the formal sum of classes obtainable by prepending one atom, with matrix entries the two-circle pairing representation numbers of §3 (Gram \(\bigl(\begin{smallmatrix}1&a&b\\a&1&n\\b&n&1\end{smallmatrix}\bigr)\)-configuration counts). The census identity reads \((1 - \mathsf{T})^{-1}(\text{atoms}) = \text{everything}\). This is the arithmetic skeleton of the missing Mayer/transfer-operator theorem (**G4** of [spectral-geometry.md](spectral-geometry.md)): the analytic operator \(\mathcal{L}_s = \sum_A |A'|^s\) of [spectral-outlook.md](spectral-outlook.md) §7 has \(\mathsf{T}\) as its "class-group shadow", and the entries are class-number-like quantities. Concrete first step: compute \(\mathsf{T}\) restricted to the Ford block (below) and identify its Fredholm determinant with a known \(L\)-function ratio.
2. **The Ford block is continued-fraction arithmetic.** Prepending \(T_i\) is the horoball induction \(n \mapsto y' + 2q'\); depth-\(k\) all-Ford words are parametrised by \((\pm U\backslash\Gamma/\pm U)^{k-1}\)-chains, i.e. by sequences of fractions, and their levels are continuant-type polynomials (depth 2: \(n = 2c^2+1\); the extremal chains of Thm. 11 are the two degenerate continuants). Conjecturally the all-Ford stratum of \(\mathrm{Cl}_+(1-n^2)\) is counted by continuant representations of \(n\) — a bijective mechanism squarely aimed at the Kronecker–Hurwitz relations ([outlook.md](outlook.md) 3.3): the right-hand side \(\sum_{d \mid m}\max(d, m/d)\) counts cusp data, and here every class *literally* acquires cusp data (its Ford prefix). Work out depth 3 in closed form first.
3. **The deep stratum and the gasket dimension.** \(R(n)\) first appears at \(n = 15\) and has density \(1 - 3/\pi\); within it, the *atom* count \(\sum_{n\le X}a(n)\) should grow like the Apollonian counting exponent, \(\asymp X^{\delta}\), \(\delta = 1.30568\ldots\) (Kontorovich–Oh applied to \(\Gamma\backslash\mathcal{A}^+\)) — i.e. **the gasket dimension should be visible inside the Hurwitz class-number census** as the growth rate of the depth-\(1\) stratum, with the depth-\(k\) strata interpolating. Numerical target: log-log fit of \(\sum a(n)\) and of \(\sum_n N_{\text{depth}=k}(n)\) against \(\delta\); theoretical target: import Kontorovich–Oh/Lee–Oh through the \(\Gamma\)-orbit structure. The local–global problem for the \(\alpha\)-spectrum (Q1 of [half-plane-monoid.md](half-plane-monoid.md)) is exactly the question "which levels have \(a(n) \neq 0\)".
4. **The phase as an address cocycle.** By §2, \(u_f\) is a Birkhoff product of atomic multipliers along the address, and by [first-power-descent.md](first-power-descent.md) it equals \(\Phi_y/\Phi_x(\beta_1,\beta_2)\). Deriving the Galois laws (twin relation \(u_fu_{\mathfrak{r}f}=1\) = reversal; conjugation = mirror) directly from the cocycle picture would give a *dynamical* proof of the reciprocity theorems, and inserting \(u^s\) into the transfer operator of (1) defines a twisted determinant whose special values should see the \(L'(0,\chi)\)-combinations of [outlook.md](outlook.md) 2.3. Cheap experiment: compute \(\log|u_f|\) against the Birkhoff sums of \(\log|A_i'|\) along the address (the kernel term is the known CM quantity) and watch the twin cancellation happen letter by letter.
5. **The \(\varepsilon\)-semigroup and aligned chains.** Equality in Theorem 5 means aligned common perpendiculars; the identity \(\varepsilon_{2q^2-1} = \varepsilon_q^2\) and the extremal chains of Theorem 11 suggest classifying all *multiplicative* chains (\(\varepsilon_n = \prod\varepsilon_{a_i}\) exactly). These should be exactly the chains lying along a single geodesic — a Schmidt-arrangement version of "geodesic continued fractions" for real quadratic units, and the natural bridge to the closed-geodesic multiplicity program ([outlook.md](outlook.md) 3.1): the aligned chains at level \(n\) refine the geodesic of length \(2\log\varepsilon_n\) into atomic segments.
6. **Depth as a class invariant.** \(\ell\) is a new integer invariant of a form class of discriminant \(1-n^2\), with \(\ell = \frac{n+1}{2}\) characterising \(\{1, \mathfrak{r}_n\}\) (Thm. 11) and \(\ell = 1\) the Apollonian stratum. What is \(\ell\) arithmetically? It is *not* a genus invariant — verified by composing out the genus partition: at \(n = 15\) the principal genus \(\{(1,0,56), (8,8,9)\}\) carries depths \(\{8, 2\}\), and every level \(15 \le n \le 41\) with \(R(n) > 0\) shows such splits — though the mirror symmetry does force \(\ell([f]) = \ell([f]^{-1})\). The extremes are arithmetic: candidates are relations to \(\min f\) (proved: \(\ell \le \min f + \frac{n-1}{2\min f}\)), to the continued fraction of the CM point, or to the Ford prefix length ("how cuspidal is the class"). A statistical study of \(\ell\) across \(\mathrm{Cl}(1-n^2)\) — mean depth, depth vs. order in the class group — is one afternoon with the existing script.
7. **The relative Hecke substitute (G7).** The cell decomposition with thin unit groups \(\Gamma^A\) is the discrete structure that replaces the (nonexistent) Hecke algebra of the pair \((\mathrm{SL}_2(\mathbb{Z}[i]), \Gamma)\). In the planned Schmidt trace formula ([spectral-outlook.md](spectral-outlook.md) §1), the geometric side at level \(n\) carries multiplicity \(H(n^2-1)\); the address now **filters** this multiplicity, with the Ford stratum (density \(3/\pi\)) attached to the cusp of \(Y\) (unipotent/Eisenstein-adjacent contributions) and the deep stratum (density \(1-3/\pi\)) to the thin gasket group. Speculative but testable: the two strata should have different spectral signatures — the Ford part governed by the continuous spectrum (the \(\zeta_K\)-ratio of the scattering matrix), the deep part by the base eigenfunction \(\lambda_0 = \delta(2-\delta)\) of the gasket manifold. A numerical splitting of the census error term along the two strata would be first evidence.
8. **Even levels and other fields.** Everything here used only laminarity, rigidity, and the classification congruence; the \(i\mathcal{S}\)-family (even \(n\), odd discriminants \(1-n^2\)) and Stange's arrangements over other \(\mathcal{O}_K\) should carry the same cell/address structure, with \(3/\pi\) replaced by the corresponding ratio of \(L\)-constants — a clean comparative invariant of Bianchi orbifolds ("the Ford density of the field"), and a new entry for the "hearing the arrangement" question of [spectral-outlook.md](spectral-outlook.md) §8.

## 9. Files

- [scripts/atomic_census.py](scripts/atomic_census.py) — the full census and all checks of §7. `python3 scripts/atomic_census.py 41`; `python3 scripts/atomic_census.py asym 1601` for §6.
- Depends on [scripts/omega.py](scripts/omega.py) (monoid, factorisation, exact disk arithmetic).
