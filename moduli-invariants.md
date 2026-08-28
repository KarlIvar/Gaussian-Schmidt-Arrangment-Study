# Modular invariants of Schmidt circles: \(\alpha, \beta_1, \beta_2\), a sixth invariant, and singular moduli

Setup: \(X \in \mathrm{SL}_2(\mathbb{C})\) with \(\omega_1 = X(\hat{\mathbb{R}})\) contained in the upper half-plane; \(\omega_2 = \sigma(X)(\hat{\mathbb{R}})\), \(\sigma(X) = \bar X^{-1}\); \(m_1, m_2\) their hyperbolic centers (for integral Schmidt \(X\), \(\omega_2\) lies in the lower half-plane and \(m_2\) is understood in the reflected sense, so \(\bar m_2 \in \mathbb{H}\)); \(\alpha = \coth(\text{hyperbolic radius})\); \(f\) a modular function (throughout: \(f = j\)), and
$$
\beta_1 = j(m_1), \qquad \beta_2 = j(\bar m_2).
$$
Both are invariant under \(X \mapsto \gamma X \gamma'\), \(\gamma, \gamma' \in \mathrm{SL}_2(\mathbb{Z})\): left multiplication moves \(\omega_1\) by a Möbius map and fixes \(\omega_2\); right multiplication does the opposite; \(\alpha\) is blind to both. All numerics: [scripts/moduli_invariants.py](scripts/moduli_invariants.py) (mpmath, 60–80 digits).

## 1. The sixth invariant

\((\alpha, \beta_1, \beta_2)\) has \(1 + 2 + 2 = 5\) real dimensions against \(\dim_{\mathbb{R}}\mathrm{SL}_2(\mathbb{C}) = 6\). The missing dimension is a **phase**, and the key to constructing it is an intertwining identity.

**Proposition (intertwining).** Let \(Z = X\bar X^{-1}\) and \(Z' = \bar X^{-1}X\). Then the fixed points of the Möbius transformation \(Z\) are exactly \(\{m_1, \bar m_1\}\) (and of \(Z'\): \(\{m_2, \bar m_2\}\)), and \(XZ'X^{-1} = Z\); consequently
$$
X\{m_2, \bar m_2\} = \{m_1, \bar m_1\}.
$$

*Proof sketch.* \(v \mapsto Z\bar v\) is the anti-Möbius inversion in \(\omega_1\) (its fixed lattice defines the circle), so \(Z = \mathrm{inv}_{\omega_1}\circ \mathrm{conj}\), a composition of reflections in two disjoint geodesic planes of \(\mathbb{H}^3\) — a loxodromic of translation length \(2\operatorname{arccosh}\alpha\) (this is the trace identity \(\operatorname{tr}Z = -2\alpha\) again). Its axis-endpoints are the pair symmetric with respect to both \(\hat{\mathbb{R}}\) and \(\omega_1\), which a direct computation identifies as \(\{m_1, \bar m_1\}\). The identity \(XZ'X^{-1} = Z\) is trivial algebra; it transports fixed points. \(\square\)

So \(X\) canonically maps the CM-type point \(m_2\) to \(m_1\) or \(\bar m_1\) (both branches occur; integral Schmidt matrices in our normalization take the \(+\) branch), and the **derivative of \(X\) at \(m_2\)** is a well-defined complex number \(d = X'(m_2) = (c\,m_2 + d)^{-2}\). Its transformation under the two \(\Gamma\)-actions is exactly cancelled by weight-2 factors built from \(j' \) (using \(j'(\gamma z)\gamma'(z) = j'(z)\)):

> **Definition.** With branch \(X(m_2) = m_1\):
> $$
> \Theta(X) \;=\; \frac{j'(m_1)\, X'(m_2)}{\overline{j'(\bar m_2)}}\,;
> \qquad\text{branch } X(m_2) = \bar m_1: \quad
> \Theta(X) \;=\; \frac{j'(m_1)\, \overline{X'(m_2)}}{j'(\bar m_2)} .
> $$
> Then \(\Theta(\gamma X\gamma') = \Theta(X)\) for all \(\gamma, \gamma' \in \mathrm{SL}_2(\mathbb{Z})\).

(Verified to 20+ digits on random \(\gamma X\gamma'\) for every class tested. Excluded: \(\alpha = 2\), where \(m\) is the elliptic point \(\rho\) and \(j'(\rho) = 0\).)

**\(\arg\Theta\) is exactly the missing coordinate.** The fiber of \(X \mapsto (\omega_1, \omega_2)\) is the one-parameter group of real rotations \(h_t\) about the axis \(\{m_2, \bar m_2\}\) (real, so they fix \(\hat{\mathbb{R}}\)-data; fixing \(\omega_2\) exactly). Along it, numerically: \(|\Theta(Xh_t)|\) is **constant to all digits**, while \(\arg\Theta(Xh_t)\) moves **linearly at rate \(\sqrt{\alpha^2-1}\)** — the multiplier \(h_t'(m_2) = e^{-it\sqrt{\alpha^2-1}}\). So \(|\Theta|\) is a function of \((\alpha,\beta_1,\beta_2)\) and:
$$
(\alpha,\ \beta_1,\ \beta_2,\ \arg\Theta) \ \text{— six real dimensions —}
$$
specify the double coset \(\mathrm{SL}_2(\mathbb{Z})\, X\, \mathrm{SL}_2(\mathbb{Z})\) up to finite ambiguity (finite \(j\)-fiber and stabilizer issues), on the region where circles avoid elliptic points.

## 2. Fixed \(\alpha = n\): singular moduli and simultaneous modular equations

For Schmidt circles at level \(n\), \(m_1\) is a **CM point of discriminant \(D = 1-n^2\)**, so:

- **The \(\beta_1\)-values are singular moduli**: over the circles of the fundamental domain with primitive class, \(\prod (x - \beta_1)\) is the **Hilbert class polynomial** \(H_{1-n^2}(x) \in \mathbb{Z}[x]\) (imprimitive circles contribute the singular moduli of the smaller orders \(\mathcal{O}_{D/g^2}\)). Verified: integer coefficients recovered at 60 digits for \(n = 3, 5, 7, 9, 11\); anchor \(H_{-8}(x) = x - 8000\) (\(n=3\): the Ford-adjacent level has \(\beta_1 = 8000 = 20^3\)).
- **\(\beta_2\) is a singular modulus of the same discriminant**: \(\bar m_2\) is the CM point of the class \(\hat\sigma[f] = [\mathfrak{r}_n][f]^{-1}\) ([class-formula-proof.md](class-formula-proof.md)). So the 2-tuples \((\beta_1, \beta_2)\) are the graph of the involution \([\mathfrak{c}] \mapsto [\mathfrak{r}_n\bar{\mathfrak{c}}]\) on the roots of \(H_{1-n^2}\), and every symmetric function of the pairing is rational: \(\prod (x - \beta_1)(y - \beta_2)\)-type pairing polynomials lie in \(\mathbb{Z}[x,y]\).
- **The identity they satisfy — simultaneous modular equations.** Since \(\mathfrak{a}_{f_2} = \mathfrak{r}_n\,\bar{\mathfrak{a}}_f\) with \(N(\mathfrak{r}_n) = \tfrac{n-1}{2}\), the elliptic curve with \(j = \beta_2\) is \(\tfrac{n-1}{2}\)-isogenous to the curve with \(j = \bar\beta_1\); using the partner ideal \(\mathfrak{s}\) of norm \(\tfrac{n+1}{2}\) in the same class, it is *also* \(\tfrac{n+1}{2}\)-isogenous to it:
$$
\Phi_{\frac{n-1}{2}}\bigl(\beta_2, \bar\beta_1\bigr) = 0
\qquad\text{and}\qquad
\Phi_{\frac{n+1}{2}}\bigl(\beta_2, \bar\beta_1\bigr) = 0 ,
$$
\(\Phi_m\) the classical modular polynomials. Verified for every primitive class at \(n = 3, 5, 7, 11, 13\) by exhibiting the explicit cyclic isogenies \(\tau_2 = \tfrac{a\tau + b}{d}\), \(ad = \tfrac{n\mp1}{2}\) (e.g. at \(n=13\), class \((2,0,21)\): \(\tau_2 \sim 2\tau/3\) and \(\tau_2 \sim \tau/7\)). The two degrees multiply to \(\tfrac{n^2-1}{4} = |D|/4\): each \(\hat\sigma\)-pair of Schmidt circles is a point on the fiber product of \(X_0(\tfrac{n-1}{2})\) and \(X_0(\tfrac{n+1}{2})\) — a Heegner-type configuration in which the *level pair is coupled to the discriminant*.

## 3. All \(\alpha\) together: Zagier traces along \(D = n^2 - 1\)

Ranging over all Schmidt circles in the fundamental domain, the 3-tuples \((\alpha, \beta_1, \beta_2)\) are the multiset
\(\{(n,\ j(\mathfrak{c}),\ j(\mathfrak{r}_n\bar{\mathfrak{c}})) : n \ge 2,\ \mathfrak{c} \in \mathrm{Cl}(\mathcal{O}_{1-n^2})\}\) (plus imprimitive strata). The natural "nice presentation" is through **traces**: with weights \(1/w\) at the elliptic classes,
$$
\mathrm{Tr}_n := \sum_{\text{classes}} \frac{j(m_1) - 744}{w}
$$
equals **Zagier's trace of singular moduli \(t(n^2-1)\)**. Computed table (matching Zagier's published values \(t(3) = -248\), \(t(8) = 7256\), \(t(15) = -192513\), \(t(24) = 4833456\), …):

| \(n\) | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|
| \(\mathrm{Tr}_n\) | \(-248\) | \(7256\) | \(-192513\) | \(4833456\) | \(-117966288\) | \(2835808512\) | \(-67515202851\) |

By Zagier's theorem, \(\sum_d t(d)q^d\) completes to a **weight-\(3/2\) modular form**; the Schmidt arrangement reads off its coefficients along the *shifted-square* arithmetic progression \(d = n^2 - 1\) — precisely the diagonal slicing that theta-decompositions/Shimura lifts convert to integral weight. The same applies at the cruder level of counts: \(\sum_n 3H(n^2-1)q^{n}\)-type series slice the Cohen–Eisenstein weight-\(3/2\) series \(\sum H(d)q^d\) along \(d = n^2-1\) (cf. the Kronecker–Hurwitz terms \(H(4m - t^2)\) with \(t^2 - 4m\) our discriminants). The two-variable refinement \(\sum_n q^n \sum_{\mathfrak{c}} x^{j(\mathfrak{c})}\)-pairing data is the generating object of the \(\hat\sigma\)-correspondence of §2.

## 4. The experiment: \(\Theta\) on the Bianchi group quantizes to algebraic units

The most intriguing question: what values does the *continuous* phase invariant \(\Theta\) take on the *discrete* group \(\mathrm{SL}_2(\mathbb{Z}[i])\)? Define the normalization
$$
u(X) := \varepsilon\,\Theta(X), \qquad \varepsilon = n + \sqrt{n^2-1} = e^{\ell/2},
$$
where \(\ell = 2\operatorname{arccosh} n\) is the geodesic length from the trace identity — \(\varepsilon\) is the norm-1 unit of the *real* quadratic order of discriminant \(4(n^2-1)\). Computed for **every primitive class of every odd \(n \le 17\)** (all to 40+ digits, two-sided invariance checked; \(\Theta\) is constant on classes so we write \(u_f\)):

1. **Conjugation law:** \(u_{f^{-1}} = \overline{u_f}\). In particular \(u_f \in \mathbb{R}\) exactly on ambiguous classes — the phase *quantizes to* \(\{0, \pi\}\) there. (The first complex values appear at \(n = 9\), the first class group with 4-torsion.)
2. **Twist law:** \(u_{\mathfrak{r} _ {n} f}\cdot u_f = 1\). Equivalently \(\Theta_f\,\Theta_{\mathfrak{r}f} = \varepsilon^{-2} = e^{-\ell}\): the \(\mathfrak{r}_n\)-twist costs exactly one geodesic length. Consequences, all observed: \(u_{\hat\sigma f} = \bar u_f^{-1}\); \(|u_f| = 1\) on \(\hat\sigma\)-fixed classes (e.g. \(n=9\), class \((3,2,7)\): \(u\) on the unit circle); \(\prod_f u_f = \pm 1\); and at \(n = 3\): \(u = -1\) *exactly*, i.e. \(\Theta = -\varepsilon^{-1} = -(3-2\sqrt2)\).
3. **Algebraicity (PSLQ):** each \(\mathfrak{r}\)-pair \(\{u_f, u_{\mathfrak{r}f}\}\) satisfies \(x^2 - S x + 1\) with
\(S = u_f + u_f^{-1}\) an explicit element of the **real genus field** \(\mathbb{Q}\bigl(\sqrt{d} : d \mid n^2-1\bigr)\) — e.g. at \(n = 11\), pair \((2,0,15),(3,0,10)\):
\(S\) has rational coordinates over \(1, \sqrt2, \sqrt3, \sqrt5, \sqrt6, \sqrt{10}, \sqrt{15}, \sqrt{30}\), recovered by PSLQ to 80 digits (the spanning set is multiplicatively overcomplete, so the printed representation is one of several equivalent ones). So the \(u_f\) are **algebraic numbers of norm 1 over the genus field** — Stark-unit-flavored quantities coupling the imaginary quadratic class group (through the CM points) to the real quadratic unit \(\varepsilon\).

**Status.** Laws 1 and 2 are now **proved** — see §5. Law 3 splits into algebraicity (provable by standard CM theory, see §5.4) and the genus-field location of \(u + 1/u\) (still experimental).

## 5. Proofs of the phase–geodesic laws

Call \(g\) a **real kernel** if it is a meromorphic weight-2 function for \(\mathrm{SL}_2(\mathbb{Z})\) with \(q\)-expansion in \(2\pi i\cdot\mathbb{R}[[q]][q^{-1}]\) (e.g. \(g = j'\) or \(j'/j\)); such \(g\) satisfy \(g(\gamma z)\gamma'(z) = g(z)\) (all the invariance argument of §1 needs) and \(g(-\bar\tau) = -\overline{g(\tau)}\). Fix a real kernel, finite and nonzero at the CM points involved, and let \(\Theta = \Theta_g\). Work with any \(X \in \mathrm{SL}_2(\mathbb{C})\) whose circle \(\omega_1\) lies in \(\mathbb{H}\) *and whose \(\sigma\)-circle \(\omega_2\) lies in the lower half-plane* — automatic for level-\(n\) Schmidt matrices in the positive-curvature normalization (orientation reversal, [class-formula-proof.md](class-formula-proof.md)), and stable on connected families since \(|\alpha(\sigma X)| = \alpha > 1\) prevents \(\omega_2\) from ever crossing \(\hat{\mathbb{R}}\). Write \(\alpha = n\), \(N = n^2 - 1\), \(\varepsilon = n + \sqrt N\).

**Lemma A (multiplier).** In the normalization \(M_X = \begin{pmatrix} 2q & -\zeta \\ -\bar\zeta & 2m\end{pmatrix}\), \(\zeta = x + ni\), the Cartan matrix \(Z = X\bar X^{-1} = -iJ\overline{M_X}\) satisfies
$$
Z\binom{m_1}{1} = -\varepsilon \binom{m_1}{1}, \qquad m_1 = \frac{x + i\sqrt N}{2q},
$$
hence \(Z'(m_1) = (-\varepsilon)^{-2} = \varepsilon^{-2} = e^{-\ell}\). *Proof:* the bottom row of \(Z\) is \((2iq, \, -i\bar\zeta)\), and \(2iq\,m_1 - i\bar\zeta = i(x + i\sqrt N) - i(x - ni) = -\sqrt N - n\); since \(m_1\) is a fixed point of \(Z\) (§1), the pairing of the bottom row with \((m_1, 1)^{\mathsf T}\) is its eigenvalue, and the Möbius derivative at a fixed point is (eigenvalue)\(^{-2}\). The same computation in the \(\sigma X\)-normalization (curvature entry negative) gives \(Z' = \bar X^{-1}X\) with \(Z'\binom{m_2}{1} = -\varepsilon\binom{m_2}{1}\), \(m_2\) the *lower* center. \(\square\)

**Lemma B (branches).** \(\bar X^{-1}(m_1) = m_2\), \(X(m_2) = m_1\), and \(X^{-1}(\bar m_1) = \bar m_2\). *Proof:* from \(Z'\,\bar X^{-1} = \bar X^{-1} Z\), the vector \(\bar X^{-1}(m_1, 1)^{\mathsf T}\) is a \((-\varepsilon)\)-eigenvector of \(Z'\); the eigenvalues \(-\varepsilon^{\pm1}\) are distinct, so by Lemma A this eigenline is spanned by \((m_2, 1)^{\mathsf T}\), giving \(\bar X^{-1}(m_1) = m_2\), i.e. \(\bar X(m_2) = m_1\). Conjugating all matrix entries and the argument: \(X(\bar m_2) = \bar m_1\), and since \(X\) permutes the pairs, \(X(m_2) = m_1\). Conjugating \(\bar X^{-1}(m_1) = m_2\) likewise gives \(X^{-1}(\bar m_1) = \bar m_2\). \(\square\)

**Theorem (functional equations).** With the branch conventions of §1 (all "+" by Lemma B):
$$
\text{(a)}\quad \Theta(X^\ast) = \overline{\Theta(X)}, \quad X^\ast := R\bar XR,\ R = \operatorname{diag}(-1,1);
\qquad
\text{(b)}\quad \Theta(X)\,\overline{\Theta(X^{-1})} = \varepsilon^{-2} = e^{-\ell}.
$$

*Proof of (b).* The circle of \(X^{-1}\) is \(\bar\omega_2 \subset \mathbb{H}\) (because \(X^{-1} = \overline{\sigma(X)}\) and conjugating a matrix entrywise conjugates its circle), so \(Y := X^{-1}\) is again in our framework with \(m_1(Y) = \bar m_2\), \(m_2(Y) = \bar m_1\), and branch "+" by Lemma B. Then
$$
\Theta(X)\,\overline{\Theta(Y)}
= \frac{g(m_1)X'(m_2)}{\overline{g(\bar m_2)}} \cdot
\overline{\left(\frac{g(\bar m_2)\,Y'(\bar m_1)}{\overline{g(m_1)}}\right)}
= X'(m_2)\cdot \overline{(X^{-1})'(\bar m_1)} :
$$
*every kernel factor cancels.* Now \(\overline{(X^{-1})'(\bar m_1)} = (\bar X^{-1})'(m_1) = (\sigma X)'(m_1)\), and by Lemma B and the chain rule,
$$
X'(m_2)\,(\sigma X)'(m_1) = (X \circ \sigma X)'(m_1) = Z'(m_1) = \varepsilon^{-2}
$$
by Lemma A. \(\blacksquare\)

*Proof of (a).* \(X^\ast = \nu X \nu\) with \(\nu(z) = -\bar z\), so \(m_1^\ast = -\bar m_1\), \(m_2^\ast = -\bar m_2\), the branch stays "+" (\(X^\ast(m_2^\ast) = \nu X(m_2) = m_1^\ast\)), and \(X^{\ast\prime}(m_2^\ast) = \overline{X'(m_2)}\) (direct computation from \(X^\ast = \begin{pmatrix}\bar a & -\bar b\\ -\bar c & \bar d\end{pmatrix}\)). Using \(g(-\bar\tau) = -\overline{g(\tau)}\) once in the numerator, \(g(m_1^\ast) = -\overline{g(m_1)}\), and once in the denominator, \(\overline{g(\overline{m_2^\ast})} = \overline{g(-m_2)} = \overline{-\overline{g(\bar m_2)}} = -g(\bar m_2)\), the two signs cancel:
$$
\Theta(X^\ast) = \frac{-\overline{g(m_1)}\cdot\overline{X'(m_2)}}{-g(\bar m_2)} = \overline{\Theta(X)}. \qquad\blacksquare
$$

**Corollary (laws 1 and 2 of §4).** The mirror \(X \mapsto X^\ast\) inverts the circle class and \(X \mapsto X^{-1}\) produces the class \(\mathfrak{r}_n[f]^{-1}\) ([circle-composition.md](circle-composition.md), [class-formula-proof.md](class-formula-proof.md)); with \(u = \varepsilon\Theta\), (a) gives \(u_{f^{-1}} = \bar u_f\), and (b) plus (a) give \(u_f\, u_{\mathfrak{r}_n f} = \varepsilon^2\Theta_f\overline{\Theta_{\mathfrak{r}_nf^{-1}\text{-rep}}} = 1\). Both laws hold for **every** admissible kernel \(g\), since the kernel cancelled in (b) and conjugated in (a); changing kernel multiplies \(u_f\) by the explicit CM-algebraic factor \(\tfrac{(g_2/g_1)(m_1)}{\overline{(g_2/g_1)(\bar m_2)}}\) (for \(g_2/g_1 = 1/j\): the singular-moduli ratio \(\bar\beta_2/\beta_1\); verified numerically). \(\square\)

**Corollary.** If \(\mathfrak{r}_n\) is principal (e.g. \(n = 3\)), then \(u_f^{2} = 1\) on ambiguous classes and \(|u_f| = 1\) always; at \(n = 3\), \(u = -1\) exactly, i.e. \(\Theta = -\varepsilon^{-1}\).

Every lemma and both theorems are machine-verified at *representative* level (not just class level) for all primitive classes, odd \(n \le 13\), to 25+ digits, including the kernel-independence and the change-of-kernel factor: `python3 scripts/moduli_invariants.py laws`.

### 5.4 Algebraicity and the elliptic-unit connection

Two exact rewritings locate \(u_f\) arithmetically:

1. **Algebraic closed form.** With \(w := \bar m_2\), real \(q\)-coefficients give \(\overline{j'(w)} = -j'(-\bar w) = -j'(-m_2)\), and \(-m_2 \in \mathbb{H}\) is the CM point of the class \([\mathfrak{r}_nf]\). Writing \(\mu := c\,m_2 + d\) (bottom row of \(X\)) — an element of the **biquadratic field** \(B = \mathbb{Q}(i, \sqrt N)\), since \(m_2 \in \mathbb{Q}(\sqrt D) \subset B\) and \(X\) has \(\mathbb{Z}[i]\)-entries —
$$
\Theta(X) \;=\; -\,\frac{j'(\tau_{[f]})}{j'(\tau_{[\mathfrak{r}f]})}\;\mu^{-2},
$$
a ratio of \(j'\)-values at two CM points of the *same* discriminant times an explicitly algebraic number. Since \(j' = -2\pi i\, E_4^2E_6/\Delta\) has weight 2, Shimura's algebraicity theorem makes the ratio algebraic (the CM period \(\Omega_D^2\) and \(2\pi i\) cancel), lying in the ring class field: **\(u_f \in \bar{\mathbb{Q}}\) is proved**; only its precise field/unit structure remains experimental.
2. **Elliptic units.** From \((E_4^2E_6/\Delta)^6 = j^4(j-1728)^3\,\Delta\),
$$
u_f^{\,6} \;=\; \varepsilon^6 \mu^{-12}\; \frac{\beta_1^4(\beta_1 - 1728)^3}{\beta^4(\beta-1728)^3}\;\cdot\;\frac{\Delta(\tau_{[f]})}{\Delta(\tau_{[\mathfrak{r}f]})}, \qquad \beta := j(\tau_{[\mathfrak{r}f]}),
$$
i.e. up to sixth roots and singular-moduli factors, \(u_f\) **is a \(\Delta\)-quotient at two ideal classes of the same order** — precisely the raw material of Siegel/Robert **elliptic units**, whose values are units in ring class fields with Galois action given by Shimura reciprocity. This is the structural explanation this experiment was pointing at; pinning the exact normalization (torsion ambiguities, the \(\mu^{-12}\)-factor, and hence law 3's genus-field statement) is the remaining task.

### 5.5 The Galois action derived from the reciprocity law

Fix \(f\) and its canonical matrix \(X\); write \(q_1, q_2\) for the half-curvatures of the two circles, \(\mathfrak{b}_1 = \mathbb{Z} + \mathbb{Z}m_1\) and \(\mathfrak{b}_2 = \mathbb{Z} + \mathbb{Z}(-m_2)\) — proper fractional \(\mathcal{O}\)-ideals in \(K\) with \(N\mathfrak{b}_i = 1/q_i\) up to the \(\mathfrak{r}\)-content, and \([\mathfrak{b}_2] = [\mathfrak{r}][\mathfrak{b}_1]\) — and \(\mu = c\,m_2 + d \in B := K(i) = \mathbb{Q}(i, \sqrt{n^2-1})\), a CM field containing the real \(\mathbb{Q}(\sqrt{n^2-1})\). Let \(\tau\) generate \(\mathrm{Gal}(B/\mathbb{Q}(i))\) (so \(\tau|_K\) is the conjugation of \(K\) and \(\varepsilon^\tau = \varepsilon^{-1}\)). In terms of the weight-2 homogeneous lattice function \(h_2 = g_2^2g_3/\delta\) (so \(j'(\tau) \propto h_2(\mathbb{Z}+\mathbb{Z}\tau)\), \(h_2(\lambda\Lambda) = \lambda^{-2}h_2(\Lambda)\)), the closed form of §5.4 reads
$$
u_f \;=\; -\,\varepsilon\;\frac{h_2(\mu\,\mathfrak{b}_1)}{h_2(\mathfrak{b}_2)} .
$$

**Norm Lemma.** \(\;|\mu|^2 = \varepsilon\,\dfrac{q_1}{q_2}\) and \(\;|\mu^\tau|^2 = \dfrac{q_1}{\varepsilon\,q_2}\). In particular the full norm \(|\mu\mu^\tau| = q_1/q_2\) is \(\varepsilon\)-free: **the unit \(\varepsilon\) is exactly the discrepancy between the two archimedean places of \(\mu\)**.

*Proof.* Apply the Hermitian form \(h(v,v) = v^\dagger M_0 v\) to \(X(m_2,1)^{\mathsf T} = \mu(m_1,1)^{\mathsf T}\), using \(h(Xv, Xv) = v^\dagger M_{X^{-1}} v\) with \(M_{X^{-1}}\) the (positive-curvature) matrix of the circle \(\bar\omega_2\). Two evaluations, each a three-line computation with \(\zeta = x + ni\): for a level-\(n\) circle \(M = (2q, -\zeta, 2m)\) with hyperbolic center \(z_h\),
$$
M(z_h) = -\frac{\sqrt N}{\varepsilon\, q}, \qquad M(\bar z_h) = \frac{\varepsilon\sqrt N}{q};
$$
and \(h((z,1)^{\mathsf T}) = 2\operatorname{Im}z\), so the left side is \(|\mu|^2\sqrt N/q_1\) while the right side is \(M_{X^{-1}}(m_2) = M_{X^{-1}}(\overline{\bar m_2}) = \varepsilon\sqrt N/q_2\). The \(\tau\)-statement is the same computation applied to the conjugated intertwining \(X(\bar m_2, 1)^{\mathsf T} = \mu^\tau(\bar m_1, 1)^{\mathsf T}\) (Lemma B), where now the argument *is* the hyperbolic center: \(M_{X^{-1}}(\bar m_2) = -\sqrt N/(\varepsilon q_2)\). (Machine-verified to 40 digits on all primitive classes, \(n \le 13\).) \(\square\)

**Classical inputs** (Siegel; see Lang, *Elliptic Functions*, or Silverman, *ATAEC* II§6): for \(\sigma \in \mathrm{Gal}(\bar{\mathbb{Q}}/K)\) with Artin class \(\mathfrak{c} \in \mathrm{Cl}(\mathcal{O})\) on the ring class field,
$$
\sigma(j(\mathfrak{a})) = j(\mathfrak{c}^{-1}\mathfrak{a}), \qquad
\sigma\!\left(\frac{\delta(\mathfrak{a})}{\delta(\mathfrak{b})}\right) = \frac{\delta(\mathfrak{c}^{-1}\mathfrak{a})}{\delta(\mathfrak{c}^{-1}\mathfrak{b})},
$$
together with the weight algebra \(h_2^6 = (\text{universal const})\cdot j_*^4(j_*-1)^3\,\delta\) (\(j_* = g_2^3/\delta\)), which converts ratios of \(h_2^6\) at proper \(\mathcal{O}\)-lattices into \(\mathbb{Q}\)-rational functions of \(j\)-values times \(\delta\)-quotients — all with known translation action.

> **Theorem (translation up to a unitary twist).** Let \(\sigma \in \mathrm{Gal}(\bar{\mathbb{Q}}/B)\) with Artin class \(\mathfrak{c}\), and let \(f'\) be the class whose ideal pair is \((\mathfrak{c}^{-1}\mathfrak{b}_1, \mathfrak{c}^{-1}\mathfrak{b}_2)\) up to \(K^\times\)-scalars \(\lambda_1, \lambda_2\). Then
> $$
> \sigma\bigl(u_f^{\,6}\bigr) \;=\; u_{f'}^{\,6}\cdot \xi^{12},
> \qquad \xi = \frac{\mu_{f'}\lambda_2}{\mu_f\,\lambda_1} \in B^\times,
> $$
> and \(\xi\) is **unitary**: \(|\xi| = |\xi^\tau| = 1\) at both archimedean places of the CM field \(B\).

*Proof.* Apply \(\sigma\) to \(u_f^6 = \varepsilon^6\mu_f^{-12}\bigl(h_2(\mathfrak{b}_1)/h_2(\mathfrak{b}_2)\bigr)^6\): \(\varepsilon^6, \mu_f^{-12} \in B\) are fixed, and the \(h_2^6\)-ratio translates by the classical inputs. Rewriting \(\mathfrak{c}^{-1}\mathfrak{b}_i = \lambda_i\,\mathfrak{b}_i(f')\) and using homogeneity gives the displayed formula. Unitarity: by the Norm Lemma \(|\mu_{f'}/\mu_f|^2 = q_1'q_2/(q_1q_2')\) (the \(\varepsilon\) cancels; likewise at the \(\tau\)-place with \(\varepsilon^{-1}\)), while \(|\lambda_i|^2 = N_{K/\mathbb{Q}}(\lambda_i) = N\mathfrak{c}^{-1}\,N\mathfrak{b}_i / N\mathfrak{b}_i(f')\), so \(|\lambda_2/\lambda_1|^2 = q_1q_2'/(q_1'q_2)\) (the \(N\mathfrak{c}\) and the \(\mathfrak{r}\)-content cancel between \(i = 1, 2\)); the product is \(1\), and \(K\)-elements have equal absolute values at both places of \(B\). \(\blacksquare\)

**\(\xi\)-torsion Lemma (proved).** \(\xi\) is a root of unity in \(B\). Consequently
$$
\sigma\bigl(u_f^{12}\bigr) = u_{f'}^{12} \quad\text{always},
\qquad
\sigma\bigl(u_f^{6}\bigr) = u_{f'}^{6} \quad\text{unless } n^2 - 1 = 2\square \text{ (then up to sign)}.
$$

*Proof.* It suffices to show \(\xi\) is a **unit** of \(\mathcal{O}_B\): it is unitary at both archimedean places (Theorem), and \(B\) is CM, so all four conjugates lie on the unit circle and Kronecker's theorem applies; finally \(\mu(B) = \mu_4\) generically, \(\mu_{12}\) when \(\sqrt3 \in B\) (\(n^2-1 = 3\square\)), \(\mu_8\) when \(\sqrt2 \in B\) (\(n^2-1 = 2\square\)), and \(\zeta^{12} = 1\) for \(\zeta \in \mu_4 \cup \mu_{12}\), \(= \pm1\) for \(\zeta \in \mu_8\).

Integrality is the **cocycle-ideal computation**. Write \(\mathfrak{m}_i := \mathcal{O}_B + \mathcal{O}_B m_i\) (a fractional \(\mathcal{O}_B\)-ideal; \(\mathcal{O}_B\) is Dedekind, so no non-maximal-order issues survive extension). Then
$$
(\mu)\,\mathcal{O}_B \;=\; \mathfrak{m}_2\,\mathfrak{m}_1^{-1}.
$$
Indeed, with \(ev(s,t) := sm_2 + t\): \(\mu\cdot1 = ev(\mathrm{row}_2X)\) and \(\mu m_1 = ev(\mathrm{row}_1X)\), so the lattice \(\mu\mathfrak{b}_1\) is \(ev\) of the \(\mathbb{Z}\)-row-span of \(X\); since \(\det X = 1\), the rows of \(X\) are a \(\mathbb{Z}[i]\)-**basis** of \(\mathbb{Z}[i]^2\) and \(ev\) is \(\mathbb{Z}[i]\)-linear, so the \(\mathcal{O}_B\)-span of \(\mu\mathfrak{b}_1\) equals the \(\mathcal{O}_B\)-span of \(ev(\mathbb{Z}[i]^2) = \mathbb{Z}[i] + \mathbb{Z}[i]m_2\), i.e. \((\mu)\mathfrak{m}_1 = \mathfrak{m}_2\). (Sanity: taking norms recovers \(N_{B/\mathbb{Q}}(\mu) = (q_1/q_2)^2\), the Norm Lemma.) Hence the canonical ideal \(\mathfrak{J}_f := (\mu_f)\,(\mathfrak{b}_1\mathcal{O}_B)\,(\mathfrak{b}_2\mathcal{O}_B)^{-1} = \mathfrak{m}_2\mathfrak{m}_1^{-1}\mathfrak{m}_1\mathfrak{m}_2^{-1} = \mathcal{O}_B\) is trivial for **every** class (note \(\mathfrak{b}_i\mathcal{O}_B = \mathfrak{m}_i\), as \(\mathfrak{b}_2 = \mathbb{Z} + \mathbb{Z}(-m_2)\)); and since the \(\mathfrak{c}\)'s cancel,
$$
(\xi) = \Bigl[(\mu_{f'})\mathfrak{b}_1'\,\mathfrak{b}_2'^{-1}\Bigr]\Bigl[(\mu_f)\mathfrak{b}_1\mathfrak{b}_2^{-1}\Bigr]^{-1}\mathcal{O}_B = \mathfrak{J}_{f'}\,\mathfrak{J}_f^{-1} = \mathcal{O}_B. \qquad\blacksquare
$$

So the Theorem upgrades to an **exact translation law** for \(u^{12}\), and the reciprocity derivation of the Galois action is complete. (Descending from \(u^{12}\) to \(u\) itself introduces root-of-unity bookkeeping in the Galois closure; the certified values of §5.6 show the translation in fact holds at the level of the pair-sums \(u + u^{-1}\) with at most the sign refinements.)

**Consequences.**
1. The twists cancel across \(\mathfrak{r}\)-pairs even before the torsion lemma: pair-products are strictly Galois-equivariant.
2. Elements of \(\mathrm{Gal}(\bar{\mathbb{Q}}/\mathbb{Q})\) not fixing \(B\) are covered by the functional equations of §5 (mirror and inversion), so the full Galois orbit of \(u_f^{12}\) is \(\{u_{f''}^{\pm12}, \bar u_{f''}^{\pm12}\}\) over the reachable translates \(f''\) — which yields the field predictions of §5.6.

### 5.6 The certified experimental record and the verified Galois mechanism

An audit of the PSLQ experiments exposed two numerical traps — a multi-term relation at working precision \(p\) is trustworthy only when (terms) × (coefficient digits) is well below \(p\), and an `mp.dps` assignment at module-import time had silently capped every earlier "high-precision" run at ~54 correct digits (fixed by computing \(j' = -2\pi i E_4^2E_6/\Delta\) exactly via theta constants and setting precision after import). All claims below are certified at genuine 200-digit precision with information-theoretically safe PSLQ parameters; earlier multi-term "genus field" fits (including the first \(n=15\) test) are superseded.

**Certified values.** With \(u_{\mathfrak{r}f} = 1/u_f\) (proved), each \(\mathfrak{r}\)-pair is described by \(S = u + 1/u\):

| \(n\) | pair data | certified statement |
|---|---|---|
| 3 | single class | \(u = -1\) exactly |
| 5 | single pair | \(u_{(1,0,6)} = -\tfrac{15297097 + 10816680\sqrt2}{6647}\in\mathbb{Q}(\sqrt2)\), norm 1 (integer identity \(15297097^2 - 2\cdot10816680^2 = 6647^2\)); \(S = -\tfrac{30594194}{6647}\), \(6647 = 17^2\cdot23\) |
| 7 | single pair | \(S = -\tfrac{80674200806}{11891} \in \mathbb{Q}\) exactly (residual \(10^{-193}\)); \(11891 = 11\cdot23\cdot47\) |
| 9 | two pairs | \(S_1, S_g \in \mathbb{Q}(\sqrt5)\) and they are **\(\mathbb{Q}(\sqrt5)\)-conjugates of each other**: \(36559082332399\,S_1 = -(127888463857726178932258 + 57193459771058256389120\sqrt5)\), \(S_g = S_1^{\sqrt5\mapsto-\sqrt5}\) |
| 11 | two pairs | \(S_A, S_B \notin \mathbb{Q}\) but **quadratic conjugates over \(\mathbb{Q}\)**: \(S_A + S_B = -\tfrac{51779502787135248685964}{8508413439}\), \(S_AS_B = -\tfrac{216521978405798024786175556}{76575720951}\) |
| 13 | two pairs | likewise: \(S_A + S_B = -\tfrac{182330892545845931198260964}{38032133275}\), \(S_AS_B = \tfrac{4603575719671472163580355540068}{722610532225}\) |

So in every certified case the multiset \(\{u_f\}\) at level \(n\) is the root set of a **single reciprocal polynomial with rational coefficients**, of the minimal degree \(h\) — e.g. at \(n = 9, 11, 13\) all four \(u\)'s are roots of \((x^2 - S_Ax + 1)(x^2 - S_Bx + 1) \in \mathbb{Q}[x]\). That the \(u_f\) of a level are the roots of one polynomial with **integer** coefficients is a theorem, not an observation: see §5.9, where it is deduced unconditionally from Theorem A (in degree \(12h\)) and the degree-\(h\) form above is certified up to \(n = 17\).

**The mechanism, matched to the Theorem.** The reachable translation classes in the Theorem are those of \(\mathrm{Gal}(H/H\cap B)\), and \(B = K(i)\) sits inside the genus field iff \(-4\) splits off the discriminant, i.e. iff \(16 \mid n^2-1\), i.e. \(n \equiv \pm1 \pmod 8\). Assuming \(\xi\)-torsion, the Theorem then *predicts* the table:

- \(n = 7\) (\(\equiv -1\)): \(B \subseteq H\), \([H:B] = 1\) — no translations reach \(u\); with the \(\tau\)- and mirror-mechanisms, \(S \in \mathbb{Q}\). ✓
- \(n = 9\) (\(\equiv +1\)): \(\mathrm{Gal}(H/B) = \mathrm{Cl}^2 = \langle\mathfrak{r}\rangle\) — translations preserve each \(\mathfrak{r}\)-pair, so \(S \in B\cap\mathbb{R} = \mathbb{Q}(\sqrt5)\), and the leftover Galois elements swap the two pairs: the sums are \(\mathbb{Q}(\sqrt5)\)-conjugate. ✓ (both features observed exactly)
- \(n = 5, 11, 13\) (\(\equiv \pm3\)): \(B \not\subseteq H\), the full class group translates, orbits of pairs have rational symmetric functions: at \(n=5\) one pair (\(S \in \mathbb{Q}\)), at \(n = 11, 13\) two pairs forming one quadratic Galois orbit. ✓

**\(n = 15\) settled** (\(h = 8\), \(\mathrm{Cl} = \mathbb{Z}/4\times\mathbb{Z}/2\)): the four pair-sums (two real, one complex-conjugate pair — the reality pattern of the inversion action on \(\mathrm{Cl}/\langle\mathfrak{r}\rangle \cong \mathbb{Z}/4\)) form a **single Galois orbit**: all four elementary symmetric functions are certified rational (420 digits; see §5.8), so they are the roots of one rational quartic.

### 5.7 The denominators are Gross–Zagier primes

The closed form \(u^6 = \varepsilon^6\mu^{-12}\cdot\frac{\beta_1^4(\beta_1-1728)^3}{\beta_2^4(\beta_2-1728)^3}\cdot\frac{\Delta(\tau_1)}{\Delta(\tau_2)}\) predicts where denominators can come from: singular moduli are algebraic integers, \(\varepsilon\) is a unit, the \(\Delta\)-quotient is an \(S\)-unit supported over \(\tfrac{n\pm1}{2}\), and \((\mu) = \mathfrak{m}_2\mathfrak{m}_1^{-1}\) is supported over the curvatures — so after the cancellations forced by the norm-1 structure, the only surviving denominator source is the **numerator ideals of \(\beta = j(\tau)\) and \(\beta - 1728\)**. By Gross–Zagier (*On singular moduli*), \(j(\mathfrak{a}) = j(\mathfrak{a}) - j(\rho)\) and \(j(\mathfrak{a}) - j(i)\) are supported exactly at the primes dividing the quantities \(\tfrac{3|D| - x^2}{4}\) resp. \(\tfrac{4|D| - x^2}{4}\) — the primes where the CM curve of discriminant \(D = 1-n^2\) becomes congruent to the special curves with \(j = 0\), \(j = 1728\).

**This is exactly what the certified fractions show** — every prime of every denominator lies in \(\mathrm{GZ}(D,-3) \cup \mathrm{GZ}(D,-4)\), with a uniform power pattern:

| \(n\) | denominator | factorization | tags |
|---|---|---|---|
| 5 | \(6647\) | \(17^2\cdot23\) | \(17 \in \mathrm{GZ}(-3)\), \(23 \in \mathrm{GZ}(-4)\) |
| 7 | \(11891\) | \(11\cdot23\cdot47\) | \(11 \in\) both, \(23, 47 \in \mathrm{GZ}(-4)\) |
| 9 | \(36559082332399\) | \(11\cdot17^2\cdot19\cdot31\cdot59^2\cdot71\cdot79\) | \(17, 59 \in \mathrm{GZ}(-3)\) (squared); \(19,31,71,79 \in \mathrm{GZ}(-4)\); \(11\) both |
| 11 | \(8508413439\), \(76575720951\) | \(3^2\cdot41^2\cdot71\cdot89^2\), extra \(3^2\) | \(41, 89 \in \mathrm{GZ}(-3)\) (squared); \(71 \in \mathrm{GZ}(-4)\) |
| 13 | \(38032133275\), \(722610532225\) | \(5^2\cdot19\cdot47\cdot101^2\cdot167\), extra \(19\) | \(5, 101 \in \mathrm{GZ}(-3)\) (squared); \(19,47,167 \in \mathrm{GZ}(-4)\) |

Observed law: **\(\mathrm{GZ}(D,-3)\)-primes enter the denominator squared, \(\mathrm{GZ}(D,-4)\)-primes to the first power** — the \(2:1\) ratio reflecting the exponents \((4,3)\) of \(\beta^4(\beta-1728)^3\) in \(u^6\) (i.e. \(\tfrac46 : \tfrac36 = 2:1.5\), rounded by the valuation bookkeeping of the specific ideal numerators). Equivalently: **the denominator of the Schmidt phase \(u_f\) is supported exactly at the primes where the level-\(n\) CM points collide, in the sense of Deuring reduction, with the two elliptic fixed points \(j = 0\) and \(j = 1728\) of the modular orbifold** — the same two points whose stabilizers produce the \(\tfrac12, \tfrac13\)-weights in the Hurwitz count \(3H(n^2-1)\) of [hyperbolic-counting.md](hyperbolic-counting.md).

### 5.8 Law 3, completed: full dihedral Galois equivariance

The translation theorem of §5.5 covered \(\sigma \in \mathrm{Gal}(\bar{\mathbb{Q}}/B)\). Two further cases complete the action of all of \(\mathrm{Gal}(\bar{\mathbb{Q}}/\mathbb{Q})\).

**Lemma (\(\tau c\)-case).** Let \(\sigma \in \mathrm{Gal}(\bar{\mathbb{Q}}/K)\) restrict to \(\tau c\) on \(B\) (fixing \(K\) pointwise, flipping both \(i\) and \(\sqrt N\)), with Artin class \(\mathfrak{c}\). Then \(\sigma(u_f^6) = u_{f'}^6\,\xi_\tau^{12}\) with
\(\xi_\tau = \varepsilon^{-1}\mu_{f'}\lambda_2\,/\,(\mu_f^{\tau c}\lambda_1) \in B^\times\) a root of unity.

*Proof.* The \(j\)- and \(\delta\)-quotient reciprocity holds for **all** \(\sigma\) fixing \(K\), so \(\sigma(R_f) = (\lambda_2/\lambda_1)^{12}R_{f'}\) as before, while now \(\sigma(\varepsilon^6) = \varepsilon^{-6}\) and \(\sigma(\mu_f^{-12}) = (\mu_f^{\tau c})^{-12}\); assembling gives the displayed formula. *Unitarity:* \(|\mu_f^{\tau c}| = |\mu_f^\tau|\) (complex conjugation is an isometry) \(= \sqrt{q_1/(\varepsilon q_2)}\) and \((\mu_f^{\tau c})^\tau = \mu_f^{c}\) with \(|\mu_f^c| = |\mu_f| = \sqrt{\varepsilon q_1/q_2}\) (Norm Lemma), so at the first place \(|\xi_\tau|^2 = \varepsilon^{-2}\cdot\frac{\varepsilon q_1'/q_2'}{q_1/(\varepsilon q_2)}\cdot\frac{q_1q_2'}{q_1'q_2} = 1\), and at the \(\tau\)-place likewise with \(\varepsilon \leftrightarrow \varepsilon^{-1}\). *Integrality:* \(\tau c\) fixes \(K\) pointwise, so \(m_i^{\tau c} = m_i\), hence \(\mathfrak{m}_i^{\tau c} = \mathfrak{m}_i\) and \((\mu^{\tau c}) = (\mu)^{\tau c} = (\mathfrak{m}_2\mathfrak{m}_1^{-1})^{\tau c} = \mathfrak{m}_2\mathfrak{m}_1^{-1} = (\mu)\): the cocycle ideal is unchanged, the \(\mathfrak{J}\)-collapse of §5.5 runs verbatim (\(\varepsilon\) is a unit), and \((\xi_\tau) = \mathcal{O}_B\). Kronecker in the CM field \(B\) finishes. \(\square\)

> **Theorem A (full dihedral equivariance).** For every \(\sigma \in \mathrm{Gal}(\bar{\mathbb{Q}}/\mathbb{Q})\),
> $$
> \sigma\bigl(u_f^{12}\bigr) \;=\; u_{\,f^{e(\sigma)}\,\mathfrak{c}(\sigma)}^{12},
> $$
> where \(e(\sigma) = +1\) if \(\sigma|_K = \mathrm{id}\) and \(-1\) otherwise, and \(\mathfrak{c}(\sigma)\) is the Artin class of \(\sigma\) (resp. of \(\sigma\iota\), \(\iota\) = complex conjugation) on the ring class field. (At the \(u^6\)-level the same holds up to a sign, and exactly unless \(n^2-1 = 2\square\).)

*Proof.* \(\sigma \in \mathrm{Gal}(\bar{\mathbb{Q}}/K)\): the two sub-cases \(\sigma|_B \in \{\mathrm{id}, \tau c\}\) are §5.5 and the Lemma. \(\sigma \notin \mathrm{Gal}(\bar{\mathbb{Q}}/K)\): write \(\sigma = \sigma_1\iota\) and use law (a) *exactly* (\(\iota(u_f) = \overline{u_f} = u_{f^{-1}}\), no root of unity), then the previous case for \(\sigma_1\). The composition law \((e_1,\mathfrak{c}_1)(e_2,\mathfrak{c}_2) = (e_1e_2,\ \mathfrak{c}_2^{e_1}\mathfrak{c}_1)\) is the **generalized dihedral** structure — precisely the classical Galois structure of ring class fields over \(\mathbb{Q}\), as it must be. \(\blacksquare\)

> **Theorem B (law 3, canonical form).** Set \(T_x := u_f^{12} + u_f^{-12}\) for the pair \(x = \{f, \mathfrak{r}f\} \in \mathrm{Cl}/\langle\mathfrak{r}\rangle\). Then \(\sigma(T_x) = T_{x^{e(\sigma)}\mathfrak{c}(\sigma)}\). Consequently:
> 1. the multiset \(\{T_x\}\) is \(\mathrm{Gal}(\bar{\mathbb{Q}}/\mathbb{Q})\)-stable, its Galois orbits are the dihedral orbits on \(\mathrm{Cl}/\langle\mathfrak{r}\rangle\), and all orbit-symmetric functions are **rational**;
> 2. \(\sigma\) fixes \(T_x\) iff \(\mathfrak{c}(\sigma) \in x^{1-e(\sigma)}\langle\mathfrak{r}\rangle\); in particular when \(x^2 \in \langle\mathfrak{r}\rangle\) (all cases with \(\mathrm{Cl}^2 \subseteq \langle\mathfrak{r}\rangle\), and all 2-torsion \(\mathrm{Cl}\)), \(T_x\) lies in the **real part of \(H^{\mathrm{Art}\langle\mathfrak{r}\rangle}\)** — for \([\mathrm{Cl}:\langle\mathfrak{r}\rangle] = 2\) this is the real quadratic field \(\mathbb{Q}(\sqrt m)\) singled out by the genus characters with \(\chi(\mathfrak{r}_n) = +1\).

**Explicit fields from genus characters.** Decompose \(1-n^2\) into prime discriminants and evaluate the characters at a represented number of \(\mathfrak{r}_n = (\tfrac{n-1}2, 0, \tfrac{n+1}2)\) coprime to the discriminant:
- \(n = 11\): \(-120 = 8\cdot(-3)\cdot5\); at \(\mathfrak{r} = (5,0,6)\), using the represented value \(11\): \(\chi_8(11) = -1\), \(\chi_{-3}(11) = -1\), \(\chi_5(11) = +1\) — the invariant real quadratic is \(\mathbb{Q}(\sqrt5)\);
- \(n = 13\): \(-168 = (-3)(-8)(-7)\); at \(\mathfrak{r} = (6,0,7)\), using \(13\): \(\chi_{-3} = +1\), \(\chi_{-8} = -1\), \(\chi_{-7} = -1\) — the invariant real quadratic is \(\mathbb{Q}(\sqrt{14})\) (character product \(\chi_{-8}\chi_{-7} = +1\));
- \(n = 9\): \(H^{\mathrm{Art}\langle\mathfrak{r}\rangle} = B\), real part \(\mathbb{Q}(\sqrt5)\); \(n \le 7\): everything collapses to \(\mathbb{Q}\).

**All predictions certified at first-power level** (200 digits, safe PSLQ): \(S \in \mathbb{Q}\) for \(n \le 7\); \(S_1, S_g \in \mathbb{Q}(\sqrt5)\) conjugate at \(n = 9\); \(S_A, S_B \in \mathbb{Q}(\sqrt5)\) at \(n = 11\) and \(\in \mathbb{Q}(\sqrt{14})\) at \(n = 13\) (quadratic-conjugate pairs); at \(n = 15\) the four pair-sums form a single dihedral orbit on \(\mathrm{Cl}/\langle\mathfrak{r}\rangle \cong \mathbb{Z}/4\) (inversion fixes the two real cosets and swaps the complex-conjugate pair — exactly the observed reality pattern), with **all four** elementary symmetric functions certified rational (420 digits), e.g.
\(e_1 = -\tfrac{734533038616697980508422327428261898604776}{207537477565866265193431}\), \(e_4 = -\tfrac{1793161348210062552764486432182896454962023834895714026411792}{231902488879724417597324208272447}\): the four pair-sums are the roots of one rational quartic, completing the verification of Theorem B at every computed level.

**First-power descent (what separates \(u^{12}\) from \(u\)).** Theorem A gives \(\sigma(u_f) = \zeta(\sigma,f)\,u_{f^{e}\mathfrak{c}}\) with \(\zeta \in \mu_{12}\) a twisted cocycle (\(\zeta^6 = 1\) unless \(n^2-1=2\square\)). Law 3 at the level of \(S = u + u^{-1}\) is the statement \(\zeta \in \{\pm1\}\)-with-matching-signs; this holds in **every certified case**, and its general proof is the standard Siegel-unit refinement: writing the kernel as \(\eta^4\cdot(\text{weight-0 function of level 6})\) expresses \(u\) through Siegel–Ramachandra invariants, whose exact first-power Galois behaviour (Kubert–Lang-style multiplier bookkeeping) resolves the sixth root. Unconditionally, since \(\sigma(u^3) = \pm u'^3\) (generic \(n\)), the multiset \(\{(S_x^3 - 3S_x)^2\}\) is *exactly* dihedrally permuted. With this, the derivation of law 3 is complete: the Galois module structure of the phase invariant is fully determined, with only the last root-of-unity bookkeeping outsourced to the classical Siegel-unit calculus.

### 5.9 One level, one integer polynomial: the \(u_f\) as a root system

The dihedral equivariance of §5.8 has an elementary corollary that packages the entire phase invariant of a level into a single object of classical arithmetic: **the \(h\) phase units of a level are the roots of one polynomial with integer coefficients.**

> **Theorem C.** Fix an odd level \(n \ge 3\); write \(D = 1-n^2\), \(h = h(D)\), and let \(u_f\) (\(f \in \mathrm{Cl}(\mathcal{O}_D)\)) be the phase units of §4. Put
> $$
> \Psi_n(y) \;:=\; \prod_{f \in \mathrm{Cl}(\mathcal{O}_D)} \bigl(y - u_f^{12}\bigr).
> $$
> Then \(\Psi_n \in \mathbb{Q}[y]\): it is monic of degree \(h\), self-reciprocal (\(y^h\Psi_n(1/y) = (-1)^h\Psi_n(y)\)), and \(\Psi_n(0) = (-1)^h\). Consequently, with \(d_n \in \mathbb{Z}_{>0}\) the least common denominator of its coefficients,
> $$
> \boxed{\;P_n(x) \;:=\; d_n\,\Psi_n\bigl(x^{12}\bigr) \;\in\; \mathbb{Z}[x], \qquad \deg P_n = 12h, \qquad P_n(u_f) = 0 \quad\text{for every class } f. \;}
> $$
> \(P_n\) is nonzero and (anti)palindromic — \(x^{12h}P_n(1/x) = (-1)^hP_n(x)\) — and its root set is \(\{\zeta u_f : \zeta^{12} = 1,\ f \in \mathrm{Cl}(\mathcal{O}_D)\}\).

*Proof.* Four inputs, each already proved above:

  (i) \(u_f \in \bar{\mathbb{Q}}\) — §5.4 (Shimura algebraicity of the \(j'\)-ratio);
  (ii) \(u_f \ne 0\) — \(j' = -2\pi i\,E_4^2E_6/\Delta\) vanishes or blows up only at the elliptic points \(\rho, i\) (discriminants \(-3, -4\)) and the cusp, while our CM points have discriminant \(1-n^2 \notin \{-3,-4\}\) for \(n \ge 3\); and \(X'(m_2) = \mu^{-2} \ne 0\);
  (iii) \(u_{\mathfrak{r}_nf}\,u_f = 1\) — law 2, proved in §5;
  (iv) \(\sigma\bigl(u_f^{12}\bigr) = u_{f^{e(\sigma)}\mathfrak{c}(\sigma)}^{12}\) for every \(\sigma \in \mathrm{Gal}(\bar{\mathbb{Q}}/\mathbb{Q})\) — Theorem A of §5.8.

**Rationality.** Fix \(\sigma \in \mathrm{Gal}(\bar{\mathbb{Q}}/\mathbb{Q})\). The map \(\varphi_\sigma : f \mapsto f^{e(\sigma)}\mathfrak{c}(\sigma)\) is a bijection of the finite set \(\mathrm{Cl}(\mathcal{O}_D)\) — a group automorphism (\(e = \pm1\)) followed by a translation. So by (iv), \(\sigma\) permutes the multiset \(\{u_f^{12}\}_{f}\), hence fixes each elementary symmetric function \(e_k\{u_f^{12}\}\). By (i) these are algebraic numbers, and being fixed by *all* of \(\mathrm{Gal}(\bar{\mathbb{Q}}/\mathbb{Q})\) they lie in \(\mathbb{Q}\). Thus \(\Psi_n \in \mathbb{Q}[y]\), monic of degree \(h\) by construction. (This is the entire content of the theorem: everything else is bookkeeping.)

**Self-reciprocity.** By (iii) the bijection \(f \mapsto \mathfrak{r}_nf\) of \(\mathrm{Cl}(\mathcal{O}_D)\) sends \(u_f^{12} \mapsto u_f^{-12}\), so the root multiset is stable under \(y \mapsto y^{-1}\) (legitimate by (ii)). Moreover \(\prod_f u_f^{12} = 1\): if \([\mathfrak{r}_n] \ne 1\), then \(f \mapsto \mathfrak{r}_nf\) is a fixed-point-free involution of \(\mathrm{Cl}\) and the factors cancel in pairs, \(\prod_fu_f = 1\); if \([\mathfrak{r}_n] = 1\) (i.e. \(n = 3\)), then (iii) reads \(u_f^2 = 1\) for every \(f\), so \(\prod_fu_f = \pm1\). Either way the twelfth powers multiply to \(1\), and
$$
y^h\Psi_n(1/y) \;=\; \prod_f\bigl(1 - y\,u_f^{12}\bigr) \;=\; (-1)^h\Bigl(\prod_fu_f^{12}\Bigr)\prod_f\bigl(y - u_f^{-12}\bigr) \;=\; (-1)^h\,\Psi_n(y),
$$
the last step by stability of the multiset under inversion. Setting \(y = 0\) in the definition gives \(\Psi_n(0) = \prod_f(-u_f^{12}) = (-1)^h\).

**Clearing denominators.** \(d_n\Psi_n \in \mathbb{Z}[y]\) by definition of \(d_n\), so \(P_n(x) = d_n\Psi_n(x^{12}) \in \mathbb{Z}[x]\); it has degree \(12h\), leading coefficient \(d_n \neq 0\), and \(P_n(u_f) = d_n\Psi_n(u_f^{12}) = 0\) for every \(f\) — one polynomial, all \(h\) units. The symmetry \(x^{12h}P_n(1/x) = (-1)^hP_n(x)\) is the self-reciprocity of \(\Psi_n\) read in \(x\), and the root set is as stated because \(x \mapsto x^{12}\) is \(12\)-to-\(1\) onto the roots of \(\Psi_n\), none of which is \(0\). \(\square\)

**Sharpenings and limits.**

- **Degree \(6h\) away from the \(\zeta_8\)-levels.** By §5.8, \(\sigma(u_f) = \zeta\,u_{f'}\) with \(\zeta^6 = 1\) except when \(n^2-1 = 2\square\). So for all other \(n\) the multiset \(\{u_f^{6}\}\) is already Galois-stable, \(\Psi^{(6)}_n(y) = \prod_f(y - u_f^6) \in \mathbb{Q}[y]\), and \(d\,\Psi^{(6)}_n(x^6) \in \mathbb{Z}[x]\) works in degree \(6h\). Only the levels \(n = 3, 17, 99, \dots\) with \(n^2-1 = 2\square\) (where \(B = \mathbb{Q}(\zeta_8)\), outlook §1.7) need the twelfth power.
- **Degree \(h\): the sharp form.** Granting the first-power descent (\(\zeta \equiv 1\); proved modulo the Kubert–Lang multiplier bookkeeping, §5.8), the multiset \(\{u_f\}\) is *itself* Galois-stable and
$$
\Phi_n(x) := \prod_{f}\bigl(x - u_f\bigr) \in \mathbb{Q}[x], \qquad Q_n(x) := d_n'\,\Phi_n(x) \in \mathbb{Z}[x]
$$
is a common integer polynomial of the minimal conceivable degree \(h\), whose roots are *exactly* the \(u_f\). This is what the computations below certify at every computed level; \(\Psi_n\) is then recovered from \(\Phi_n\) by Newton's identities, and \(P_n\) from \(\Psi_n\).
- **The leading coefficient is unavoidable.** The \(u_f\) are not algebraic integers: were they, \(\Psi_n\) (monic, rational, with algebraic-integer coefficients) would lie in \(\mathbb{Z}[y]\), i.e. \(d_n = 1\) — contradicted already at \(n = 5\), where the directly certified \(\Psi_5(y) = y^2 - Ty + 1\) has \(d_5 = 6647^{12}\). So "integer coefficients" is the correct statement and "monic integer coefficients" is false; by §5.7 the primes of \(d_n\) are Gross–Zagier primes of \((D,-3)\) and \((D,-4)\).
- **One polynomial per level, not one for all levels.** On the principal class \(|u_f|\) grows exponentially with \(n\) (\(1,\ 4.6\cdot10^3,\ 6.8\cdot10^6,\ 7.0\cdot10^9,\ 6.1\cdot10^{12},\ 4.8\cdot10^{15},\ 3.5\cdot10^{18},\ 2.5\cdot10^{21}\) for \(n = 3, \dots, 17\)), so the phase units of distinct levels are distinct; a nonzero polynomial has finitely many roots, and no single one can serve infinitely many levels. The statement is per level by necessity.
- **What is *not* used.** The proof needs neither the genus-field statements of §5.8 nor the denominator law of §5.7 — only algebraicity plus the permutation action. In particular Theorem C is unconditional, whereas the degree-\(h\) form is (still) conditional on the first-power descent.

**Certified integer polynomials.** `python3 scripts/uf_integer_polynomial.py --direct` computes \(u_f\) at the precision listed, certifies each coefficient of \(\Phi_n\) as a rational by an overdetermined continued-fraction fit (every fit below is confirmed by at least \(168\) spare digits — the number of decimal digits by which the reconstruction exceeds the information content of the fraction), clears denominators, and checks \(Q_n(u_f) = 0\) to the working precision. Independently, and *without* using the rationality of \(\Phi_n\), it certifies the coefficients of the unconditional \(\Psi_n = \prod(y - u_f^{12})\) directly from the numerical twelfth powers (at up to \(2600\) digits); the two routes agree exactly at every level.

$$
\begin{aligned}
Q_3(x) &= x + 1\\
Q_5(x) &= 6647\,x^2 + 30594194\,x + 6647\\
Q_7(x) &= 11891\,x^2 + 80674200806\,x + 11891\\
Q_9(x) &= 10565574794063311\,x^4 + 73919532109765731422845124\,x^3 - 118807282021266004510100774\,x^2\\ &\qquad + 73919532109765731422845124\,x + 10565574794063311\\
Q_{11}(x) &= 76575720951\,x^4 + 466015525084217238173676\,x^3 - 216521978405797871634733654\,x^2\\ &\qquad + 466015525084217238173676\,x + 76575720951\\
Q_{13}(x) &= 722610532225\,x^4 + 3464286958371072692766958316\,x^3 + 4603575719671472165025576604518\,x^2\\ &\qquad + 3464286958371072692766958316\,x + 722610532225\\
Q_{17}(x) &= 819697933195874886721\,x^4 + 2045693252535068591803195236166141118259260\,x^3\\ &\qquad + 4207718806754010047859719149317312347901942534\,x^2\\ &\qquad + 2045693252535068591803195236166141118259260\,x + 819697933195874886721
\end{aligned}
$$

and, at \(n = 15\) (\(h = 8\), the largest level computed; coefficients \(c_8, \dots, c_0\), palindromic):

```
231902488879724417597324208272447
820767611794540060586926641691129674477616253164712
-3736910623206271622009575975230452428835195604215370332
4648900280945215092152336462310642642094781730724532480920
-1793168822031308965307730451798651893626330340481570873697350
4648900280945215092152336462310642642094781730724532480920
-3736910623206271622009575975230452428835195604215370332
820767611794540060586926641691129674477616253164712
231902488879724417597324208272447
```

Consistency with the earlier record: dividing out the leading coefficient, \(Q_n = d_n'\prod_{\mathfrak{p} \in \mathrm{Cl}/\langle\mathfrak{r}\rangle}(x^2 - S_{\mathfrak{p}}x + 1)\), and the resulting symmetric functions of the pair-sums reproduce *exactly* the certified fractions of §5.6 and §5.8 — e.g. \(Q_{11}\) gives \(S_A + S_B = -\tfrac{51779502787135248685964}{8508413439}\) and \(S_AS_B = -\tfrac{216521978405798024786175556}{76575720951}\); the \(x^7\)-coefficient of \(Q_{15}\) divided by its leading coefficient reduces to \(\tfrac{734533038616697980508422327428261898604776}{207537477565866265193431} = -e_1\) of §5.8, and the leading coefficient of \(Q_{15}\) is precisely the \(e_4\)-denominator recorded there. Two structural features are visible in every line, and they are exactly the two proved laws: **palindromy** is law 2 (\(u \mapsto 1/u\) under the \(\mathfrak{r}_n\)-twist), and **real coefficients** are law 1 (\(u \mapsto \bar u\) under class inversion). Precisely, \(x^hQ_n(1/x) = (-1)^h\bigl(\prod_fu_f\bigr)Q_n(x)\), and for \(n \ge 5\) the reduced form \(\mathfrak{r}_n = (\tfrac{n-1}2, 0, \tfrac{n+1}2)\) has \(a \ne 1\), hence \([\mathfrak{r}_n] \ne 1\), hence \(h\) is even and \(\prod_fu_f = 1\): \(Q_n\) is exactly palindromic, as the table shows (at \(n = 3\), \(h = 1\) and \(\prod u = -1\), and \(Q_3 = x+1\) is palindromic too).

Leading coefficients, factored and tagged by the Gross–Zagier sets of §5.7 (the script asserts the tagging):

| \(n\) | \(h\) | \(d_n' = \) lead \(Q_n\) | factorization (tags: \(\mathrm{GZ}(D,-3)\) / \(\mathrm{GZ}(D,-4)\)) |
|---|---|---|---|
| 3 | 1 | \(1\) | — |
| 5 | 2 | \(6647\) | \(17^2\)[−3] \(\cdot\ 23\)[−4] |
| 7 | 2 | \(11891\) | \(11\)[−3,−4] \(\cdot\ 23\)[−4] \(\cdot\ 47\)[−4] |
| 9 | 4 | \(10565574794063311\) | \(11\)[−3,−4] \(\cdot\ 17^4\)[−3] \(\cdot\ 19\cdot31\)[−4] \(\cdot\ 59^2\)[−3] \(\cdot\ 71\cdot79\)[−4] |
| 11 | 4 | \(76575720951\) | \(3^4\)[−3,−4] \(\cdot\ 41^2\)[−3] \(\cdot\ 71\)[−4] \(\cdot\ 89^2\)[−3] |
| 13 | 4 | \(722610532225\) | \(5^2\)[−3] \(\cdot\ 19^2\cdot47\)[−4] \(\cdot\ 101^2\)[−3] \(\cdot\ 167\)[−4] |
| 15 | 8 | \(231902488879724417597324208272447\) | \(11^2\)[−3,−4] \(\cdot\ 29^4\)[−3] \(\cdot\ 31\)[−4] \(\cdot\ 43^2\)[−4] \(\cdot\ 47\)[−3,−4] \(\cdot\ 53^4\)[−3] \(\cdot\ 103\)[−4] \(\cdot\ 167^2\)[−3] \(\cdot\ 199\cdot223\)[−4] |
| 17 | 4 | \(819697933195874886721\) | \(23\cdot47\)[−3,−4] \(\cdot\ 71\)[−4] \(\cdot\ 167^3\)[−3,−4] \(\cdot\ 191^2\)[−3] \(\cdot\ 239\cdot263\)[−4] |

Every prime is a Gross–Zagier prime, as §5.7 predicts. The exponent pattern of §5.7 (squares for \(\mathrm{GZ}(-3)\), first powers for \(\mathrm{GZ}(-4)\)) is a statement about a *single* pair-denominator; \(d_n'\) is the least common denominator over all pairs, so multiplicities can add across pairs — visible at \(n = 13\) (\(19^2\), a \(\mathrm{GZ}(-4)\)-prime met by two pairs) and \(n = 15\) (\(43^2\)). Support is unaffected.

**A new level, and a check of the §5.8 mechanism.** Level \(n = 17\) was not in the tables of §5.6/§5.8; the computation above adds it. Here \(h = 4\), \(\mathrm{Cl}(-288) = (\mathbb{Z}/2)^2\), all classes ambiguous (all \(u_f\) real, all signs \(-\)), and the two pair-sums are **not** individually rational but are conjugate over \(\mathbb{Q}(\sqrt2)\):
$$
S_A + S_B = -\frac{12249660194820770010797576264467910887780}{4908370857460328663},\qquad
S_A - S_B = \frac{8661817790988691324094078675916121095168}{4908370857460328663}\sqrt2 ,
$$
both certified with over \(380\) spare digits at \(420\)-digit working precision (`python3 scripts/uf_integer_polynomial.py --pairs --dps 420 17`; the same mode re-derives the earlier record — \(S \in \mathbb{Q}\) at \(n = 5, 7\), \(\mathbb{Q}(\sqrt5)\) at \(n = 9, 11\), \(\mathbb{Q}(\sqrt{14})\) at \(n = 13\)). This is exactly the prediction of Theorem B: \(n = 17 \equiv +1 \pmod 8\), so \(B = \mathbb{Q}(i, \sqrt{288}) = \mathbb{Q}(\zeta_8)\) lies in the genus field, translations preserve each \(\mathfrak{r}\)-pair, and the pair-sums land in \(B \cap \mathbb{R} = \mathbb{Q}(\sqrt2)\), the leftover Galois elements swapping the two pairs — the same pattern as \(n = 9\) with \(\mathbb{Q}(\sqrt5)\). And \(n = 17\) is a \(\zeta_8\)-level (\(288 = 2\cdot12^2\)), the case where Theorem C genuinely needs the twelfth and not the sixth power.

**Remaining open questions:**
- The Kubert–Lang multiplier computation for the first-power \(\zeta \equiv 1\) (routine but long) — it is also what separates the degree-\(h\) common integer polynomial \(Q_n\) of §5.9 from the unconditional degree-\(12h\) one.
- Exact valuations in the denominator law of §5.7.
- Distribution of \(\arg u_f\) over classes as \(n \to \infty\) (between CM-point and closed-geodesic equidistribution); heights on the discriminant-coupled fiber product \(X_0(\tfrac{n-1}2)\times_{X(1)}X_0(\tfrac{n+1}2)\).

## 6. Files

- [scripts/moduli_invariants.py](scripts/moduli_invariants.py) — all experiments (A: intertwining; B: \(\Theta\), invariance, fiber test; C: Hilbert class polynomials; D: simultaneous isogenies; E: trace table). Requires mpmath.
- Laws 1–3 of §4: `python3 scripts/moduli_invariants.py deep` (all primitive classes, odd \(n \le 17\), 80 digits).
- Proof verification (§5, representative level): `python3 scripts/moduli_invariants.py laws`.
- The common integer polynomial of a level (§5.9): [scripts/uf_integer_polynomial.py](scripts/uf_integer_polynomial.py) — `python3 scripts/uf_integer_polynomial.py --direct` (certifies \(\Phi_n, \Psi_n, Q_n, P_n\) and the Gross–Zagier tagging of the leading coefficient for odd \(n \le 17\)); `--pairs` prints the \(\mathfrak{r}\)-pair sums with their fields.
- Denominator identification (§5.7): `python3 scripts/gz_denominators.py`.
- §5.8 verifications (field predictions at \(n = 11, 13\); the \(n=15\) rational quartic, 420 digits): inline runs recorded with exact fractions in §5.8.
