# The double-coset product \((X, Y) \mapsto X\Gamma Y\) in the \((\alpha, \beta_1, \beta_2, \arg u)\)-coordinates

Companion to [atomic-census.md](atomic-census.md). There the product of two \(\Gamma\)-double cosets (\(\Gamma = \mathrm{SL}_2(\mathbb{Z})\)) — the \(\Gamma\)-indexed family \(\{\Gamma X\gamma Y\Gamma\}_{\gamma \in \Gamma}\) — was analysed in the \(\alpha\)-coordinate: \(\alpha(X\gamma Y) = \langle M_{X^{-1}}, \gamma M_Y\rangle\). This note answers the natural next question: *what does the operation do to the remaining coordinates \(\beta_1, \beta_2, \arg u\) of [moduli-invariants.md](moduli-invariants.md)?*

The answer is that the operation has exactly **three** invariant-level laws, one per coordinate block — and the count is forced: for fixed \(X, Y\) the product family \(\{XgY : g \in \mathrm{SL}_2(\mathbb{R})\}\) is (generically) \(3\)-dimensional inside the \(6\)-dimensional double-coset space, so its image must be cut out by three relations. They are:

| coordinate | law |
|---|---|
| \(\beta_1\) | the \(\beta_1\)-circle of **every** member stays at inversive distance \(\alpha(Y)\) from \(\omega_1(X)\): the product's \(\beta_1\) remembers only \(X\) (Thm. 1) |
| \(\beta_2\) | dually, the \(\beta_2\)-circle stays at inversive distance \(\alpha(X)\) from \(\omega_2(Y)\): the product's \(\beta_2\) remembers only \(Y\) (Thm. 1) |
| \(\arg u\) | the phase is **anti-additive up to an explicit geometric coboundary**: \(\Theta(X\gamma Y)\,\Theta(X)\,\Theta(Y) = K\,Q^2\) with \(K\) a product of kernel values at the six hyperbolic centers and \(Q\) a cross-ratio through the junction point (Thm. 3) |

All verified in [scripts/product_cocycle.py](scripts/product_cocycle.py): **1205 checks passing** (relations R1/R2 in exact integer arithmetic on 200 random products; the distance form and the cocycle to 50 digits; the full \(\Theta\)-identity with the \(j'\)-kernel).

## 1. Marginal rigidity: \(\beta_1\) sees only \(X\), \(\beta_2\) sees only \(Y\)

Write \(W = X\gamma Y\), \(a = \alpha(X)\), \(b = \alpha(Y)\), \(n = \alpha(W)\), and \(\omega_1, \omega_2\) for the two circles of a double coset (\(\omega_1(V) = V(\hat{\mathbb{R}})\), \(\omega_2(V) = \sigma(V)(\hat{\mathbb{R}})\)).

> **Theorem 1 (Gram rigidity).** For every \(\gamma \in \Gamma\) (indeed every \(g \in \mathrm{SL}_2(\mathbb{R})\)):
> $$
> \bigl\langle \omega_1(X\gamma Y),\ \omega_1(X)\bigr\rangle = \alpha(Y),
> \qquad
> \bigl\langle \omega_2(X\gamma Y),\ \omega_2(Y)\bigr\rangle = \alpha(X).
> $$
> Equivalently, the triple \((\hat{\mathbb{R}}, \omega_1(X), \omega_1(W))\) has Gram matrix
> \(\begin{pmatrix} 1 & a & n\\ a & 1 & b\\ n & b & 1\end{pmatrix}\)
> — the *same* Gram as the junction triple \((\hat{\mathbb{R}}, \text{circle}(X^{-1}), \gamma\,\omega_1(Y))\), to which it is congruent by \((X\gamma)^{-1}\) — and dually on the \(\omega_2\)-side with \(a \leftrightarrow b\).

*Proof.* \(\omega_1(W) = X\gamma(\omega_1(Y))\) while \(\omega_1(X) = X(\hat{\mathbb{R}}) = X\gamma(\hat{\mathbb{R}})\); Möbius invariance of the pairing gives \(\langle\omega_1 W, \omega_1 X\rangle = \langle\omega_1(Y), \hat{\mathbb{R}}\rangle = \alpha(Y)\). Dually \(\omega_2(W) = \sigma(Y)\gamma^{-1}(\omega_2(X))\) and \(\omega_2(Y) = \sigma(Y)\gamma^{-1}(\hat{\mathbb{R}})\). \(\square\)

To convert to the CM points \(m_1, m_2\) (whose \(j\)-values are \(\beta_1, \beta_2\)), one small lemma of independent use — the pairing of two Schmidt circles as hyperbolic circles ([hyperbolic-counting.md](hyperbolic-counting.md) §2):

> **Lemma 2.** For circles \(\omega, \omega'\) in \(\mathbb{H}\) with hyperbolic centers \(p, p'\) and levels \(\alpha, \alpha'\):
> $$
> \langle \omega, \omega'\rangle \;=\; \alpha\alpha' \;-\; \cosh d_{\mathbb{H}}(p, p')\,\sqrt{(\alpha^2-1)(\alpha'^2-1)} .
> $$

(Direct computation in the normalisation \(p = i\); note the check \(\omega = \omega' \Rightarrow \langle,\rangle = 1\).) Hence:

> **Corollary (CM-circle law).** For every member \(W\) of the family,
> $$
> \cosh d_{\mathbb{H}}\bigl(m_1(X),\, m_1(W)\bigr) \;=\; \frac{a\,n - b}{\sqrt{(a^2-1)(n^2-1)}},
> \qquad
> \cosh d_{\mathbb{H}}\bigl(m_2(Y),\, m_2(W)\bigr) \;=\; \frac{b\,n - a}{\sqrt{(b^2-1)(n^2-1)}} .
> $$
> So at each level \(n\) attained by the family, the \(\beta_1\)-values of the products are singular moduli \(j(m)\) with \(m\) confined to **one hyperbolic circle centered at the CM point of \(X\)**, of radius determined by the three levels alone; the \(\beta_2\)-values likewise orbit the CM point of \(Y\). In the extremal (aligned) case \(\varepsilon_n = \varepsilon_a\varepsilon_b\) of the superadditivity theorem, both distances are \(0\):
> $$
> \varepsilon_{\alpha(XY)} = \varepsilon_{\alpha(X)}\varepsilon_{\alpha(Y)}
> \iff m_1(XY) = m_1(X) \ \text{ and } \ m_2(XY) = m_2(Y),
> $$
> i.e. **on the aligned stratum the operation literally splices the invariants: \(\beta_1\) from the left factor, \(\beta_2\) from the right factor** — the sharpest possible form of the heuristic that composition couples the \(\beta_2\)-side of \(X\) to the \(\beta_1\)-side of \(Y\) and leaves the outer sides untouched.

(The forward implication of the displayed equivalence is the computation \(an - b = \sqrt{a^2-1}\,\bigl(b\sqrt{a^2-1} + a\sqrt{b^2-1}\bigr)\) at \(n = ab + \sqrt{(a^2-1)(b^2-1)}\), which makes \(\cosh d = 1\) exactly.)

**Symmetry of the whole package.** Since \(\sigma(X\gamma Y) = \sigma(Y)\gamma^{-1}\sigma(X)\), the involution maps the family \(X\Gamma Y\) *member-wise* onto \(\sigma(Y)\Gamma\sigma(X)\): the operation is commutative up to \(\sigma\). On invariants (pure twist, [atomic-census.md](atomic-census.md) Thm. 6, plus the laws \(u_{\mathfrak{r}f}u_f = 1\)): reversing the factors and twisting each by its \(\mathfrak{r}\)-class reproduces the same family with \(u \mapsto 1/u\). This is the "microscopic" form of the address-reversal theorem, and it swaps the two marginal laws of Theorem 1 into each other.

## 2. The phase cocycle

The third relation is an exact multiplicative law for the sixth coordinate. Recall the branch convention \(V(m_2(V)) = m_1(V)\) (Lemma B of [moduli-invariants.md](moduli-invariants.md); it holds for all representatives used here, and the script confirms the \(+\) branch throughout). The only analytic input is the two-point identity for a Möbius map \(T\):
\((Tp - Tq)^2 = T'(p)\,T'(q)\,(p-q)^2\).

> **Theorem 3 (triple-product law).** Let \(X, Y \in \mathrm{SL}_2(\mathbb{C})\) be such that \(X\), \(Y\), \(W = XY\) all carry the \(+\) branch, and let
> $$
> p \;=\; Y(m_2(W)) \;=\; X^{-1}(m_1(W))
> $$
> be the **junction point**. Then, exactly,
> $$
> \boxed{\;W'(m_2W)\; X'(m_2X)\; Y'(m_2Y)
> \;=\;
> \left[\frac{(m_1W - m_1X)\,(p - m_1Y)}{(p - m_2X)\,(m_2W - m_2Y)}\right]^{2}.\;}
> $$

*Proof.* Chain rule: \(W'(m_2W) = X'(p)\,Y'(m_2W)\). Multiply by \(X'(m_2X)Y'(m_2Y)\) and apply the two-point identity twice: \(X'(p)X'(m_2X) = (Xp - Xm_2X)^2/(p - m_2X)^2 = (m_1W - m_1X)^2/(p-m_2X)^2\), and \(Y'(m_2W)Y'(m_2Y) = (Ym_2W - Ym_2Y)^2/(m_2W - m_2Y)^2 = (p - m_1Y)^2/(m_2W-m_2Y)^2\). \(\square\)

> **Corollary 4 (phase cocycle).** For any admissible real kernel \(g\) (e.g. \(g = j'\)) and \(W = X\gamma Y\):
> $$
> \Theta(W)\,\Theta(X)\,\Theta(Y) \;=\; K_\gamma\, Q_\gamma^{\,2},
> \qquad
> K_\gamma = \prod_{V \in \{W, X\gamma, Y\}} \frac{g(m_1(V))}{\overline{g(m_2(V)^{\mathbb{H}})}},
> \quad
> Q_\gamma = \frac{(m_1W - m_1X)(p - m_1Y)}{(p - m_2(X\gamma))(m_2W - m_2Y)},
> $$
> (\(m_2^{\mathbb{H}}\) the upper-half-plane member of the \(m_2\)-pair; \(m_1(X\gamma) = m_1(X)\), \(m_2(X\gamma) = \gamma^{-1}m_2(X)\)). Taking arguments, with \(u = \varepsilon\Theta\) and \(\varepsilon > 0\):
> $$
> \arg u(X\gamma Y) \;\equiv\; -\arg u(X) - \arg u(Y) \;+\; \arg K_\gamma + 2\arg Q_\gamma \pmod{2\pi}.
> $$
> **The operation is anti-additive in the phase, with an explicitly geometric coboundary**: the correction term involves only kernel values at the six hyperbolic centers and one cross-ratio through the junction point — no new modular or transcendental data. The \(\gamma\)-dependence of the family enters only through \(Q_\gamma, K_\gamma\).

Remarks.

1. **This is the missing third relation.** \(|\Theta|\) is a function of \((\alpha, \beta_1, \beta_2)\) ([moduli-invariants.md](moduli-invariants.md) §1), so the modulus part of Corollary 4 is consistent with Theorem 1; the argument part pins \(\arg u(W)\) given the rest — completing the count \(6 - 3 = 3\) of constraints on the product family.
2. **Anti-additivity is forced by the symmetric-pair structure.** The Gelfand anti-involution makes the natural product law contravariant; the sign matches the proved twin law \(u_f\,u_{\mathfrak{r}f} = 1\) (products invert phases), whose derivation (law (b) of [moduli-invariants.md](moduli-invariants.md) §5, \(\Theta(X)\overline{\Theta(X^{-1})} = \varepsilon^{-2}\)) is exactly the degenerate limit \(Y \to X^{-1}\) of Theorem 3, where all kernel factors cancelled for the same telescoping reason.
3. **Iterating along an address** recovers the Birkhoff-product expression of the phase ([atomic-census.md](atomic-census.md) §2): the cocycle is its two-letter germ.
4. **Arithmetic content at Schmidt triples.** For Schmidt classes all the ingredients are algebraic: the \(m\)-points and \(p\) are quadratic surds, and \(K\) is a ratio of \(j'\)-CM-values, algebraic by Shimura. So the identity states: the product of the three phases \(u\) — elements of the ring class fields of the **three different discriminants** \(1-a^2\), \(1-b^2\), \(1-n^2\) — is an explicit product of elementary algebraic factors. Since each factor separately obeys its own dihedral Galois law ([first-power-descent.md](first-power-descent.md)), the identity forces a matching equivariance on \(K_\gamma Q_\gamma^2\): a **reciprocity constraint coupling CM points of three distinct imaginary quadratic fields**. Nothing in the classical single-discriminant theory sees this coupling.

## 3. Verification

[scripts/product_cocycle.py](scripts/product_cocycle.py) (mpmath at 70 digits, exact integer arithmetic where stated; run inside a venv with mpmath):

1. R1 and R2 of Theorem 1 for 200 random products \(X\gamma Y\) of random Schmidt representatives (levels \(3\)–\(11\)) — exact;
2. the CM-distance form (Corollary of Lemma 2) to \(50\) digits;
3. the \(+\)-branch and junction-point consistency \(X(p) = m_1(W)\) on every sample;
4. Theorem 3 to \(50\) digits on 200 samples;
5. Corollary 4 with the \(j'\)-kernel on 5 samples (relative error \(< 10^{-50}\)).

**1205 checks, all passing.**

## 4. Outlook

1. **The cocycle class.** \(c(X, Y) := \arg K + 2\arg Q\) is a \(2\)-cocycle for the (partially defined) product on double cosets. Is it a coboundary — i.e. is there a canonical "phase potential" \(\phi\) on double cosets with \(c = \phi(W) + \phi(X) + \phi(Y)\)? On the Ford stratum, where the phase theory degenerates to the cusp ([outlook.md](outlook.md) 2.4), \(c\) should reduce to a classical Dedekind–Rademacher cocycle; identifying it would tie the phase of a Schmidt circle to \(\eta\)-multiplier arithmetic by a second, product-theoretic route.
2. **Inter-level reciprocity.** Make Remark 4 precise: apply \(\sigma \in \mathrm{Gal}(\bar{\mathbb{Q}}/\mathbb{Q})\) to the triple identity at a Schmidt triple and match the dihedral translation laws of the three factors against the Galois action on \(K Q^2\). The expected output is a trilinear symbol on \(\mathrm{Cl}(1-a^2) \times \mathrm{Cl}(1-b^2) \times \mathrm{Cl}(1-n^2)\)-triples occurring in products — a new structure with no single-level counterpart.
3. **The aligned stratum as a functor.** On aligned products the operation is \((\beta_1, \beta_2)\)-splicing (Corollary to Lemma 2) and the junction cross-ratio degenerates; the phase law should become exactly multiplicative after the \(\varepsilon\)-normalisation. Classify aligned Schmidt products (the \(\varepsilon\)-multiplicative chains of [atomic-census.md](atomic-census.md) §8.5) and check whether \(u\) restricted to them is a homomorphism into the unit group.
4. **Effective multiplication table.** Theorem 1 + the pairing law reduce the multiplication table of double cosets to: choose \(n\) in the pairing spectrum, then a point on each of the two CM-circles, then the phase from Corollary 4. This is a complete, finite recipe for the "relative Hecke product" of gap G7 — worth writing as an algorithm and testing against the depth-2 census of [atomic-census.md](atomic-census.md).
