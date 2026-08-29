# First-power descent: the exact Galois law for \(u_f\), and irreducibility of the level polynomials

This document closes the first open problem of [moduli-invariants.md](moduli-invariants.md) §5.8: the root-of-unity cocycle in the Galois action on the phase invariant is trivial, \(\zeta \equiv 1\), so the translation law of Theorem A holds **at first power**,
$$
\sigma(u_f) \;=\; u_{f^{e(\sigma)}\,\mathfrak{c}(\sigma)} \qquad \text{for every } \sigma \in \mathrm{Gal}(\bar{\mathbb{Q}}/\mathbb{Q}).
$$
The mechanism is the isogeny-differential route recommended in [outlook.md](outlook.md) §2.1: \(u_f\) is, on the nose, the **derivative of the modular correspondence** \(X_0(r_0)\) evaluated at the pair of singular moduli it links — a value of a rational function with \(\mathbb{Q}\)-coefficients at CM points, on which Galois acts with no transcendental (and no root-of-unity) ambiguity. Consequences: \(u_f\) lies in the ring class field, law 3 of §4 becomes a theorem, the certified table of §5.6 becomes a theorem (re-derived by exact arithmetic with no PSLQ input), and the level polynomials
$$
\Pi_n(x) \;=\; \prod_{f \in \mathrm{Cl}(1-n^2)} \bigl(x - u_f\bigr)
$$
are **irreducible over \(\mathbb{Q}\)** — proved below for every computed level, and for general \(n\) modulo only the pairwise distinctness of the \(u_f\).

Notation as in moduli-invariants.md: \(n \ge 3\) odd, \(N = n^2-1\), \(D = 1-n^2\), \(K = \mathbb{Q}(\sqrt D)\), \(\mathcal{O}\) the order of discriminant \(D\) with class group \(\mathrm{Cl}\), \(h = |\mathrm{Cl}|\), ring class field \(H\), \(\varepsilon = n + \sqrt N\), \(\omega_0 = \sqrt D/2\), and
$$
r_0 = \tfrac{n-1}{2}, \qquad s_0 = \tfrac{n+1}{2}, \qquad
\mathfrak{r} = [r_0, \omega_0], \qquad \mathfrak{s} = [s_0, \omega_0],
$$
the ambiguous twist ideals (\(N\mathfrak{r} = r_0\), \(N\mathfrak{s} = s_0\), \([\mathfrak{s}] = [\mathfrak{r}]\), \(\mathfrak{r}^2 = (r_0)\), \(\mathfrak{s}^2 = (s_0)\), and \(\mathfrak{r}\mathfrak{s} = (\omega_0)\), since \([r_0 s_0, \omega_0] = \omega_0[1,\omega_0]\) from \(\omega_0^2 = -r_0s_0\)). Per class \(f = (a,b,c)\): \(m_1 = \frac{-b+\sqrt D}{2a}\), \(\mathfrak{b}_1 = [1, m_1]\), \(\mathfrak{b}_2 = [1, -m_2]\), \(\mu = c\,m_2 + d \in B = \mathbb{Q}(i,\sqrt N)\), and the proven closed form (§5.4–5.5)
$$
u_f \;=\; -\,\varepsilon\;\mu^{-2}\,\frac{h_2(\mathfrak{b}_1)}{h_2(\mathfrak{b}_2)},
\qquad h_2 = g_2^2g_3/\delta \ \ (\text{weight } 2,\ h_2(\lambda\Lambda) = \lambda^{-2}h_2(\Lambda)).
$$
All numerics: [scripts/first_power_descent.py](scripts/first_power_descent.py).

## 1. The derived correspondence function

Let \(m \ge 1\) and let \(\Phi_m(x,y) \in \mathbb{Z}[x,y]\) be the classical modular polynomial of cyclic \(m\)-isogenies: \(\Phi_m(j(\tau'), j(\tau)) = 0\) iff the corresponding lattices are cyclically \(m\)-isogenous; \(\Phi_m\) is irreducible, symmetric for \(m>1\), and \(\mathbb{Q}(X_0(m)) = \mathbb{Q}(j(\tau), j(m\tau))\) with \(\Phi_m(j(\tau), j(m\tau)) = 0\) identically.

**Definition.** For a proper \(\mathcal{O}\)-ideal \(\mathfrak{a}\), set
$$
V(\mathfrak{a}) \;:=\; \frac{h_2(\mathfrak{a})}{h_2(\mathfrak{r}^{-1}\mathfrak{a})}.
$$
Homogeneity of weight 2 kills any common scalar: \(V(\lambda\mathfrak{a}) = V(\mathfrak{a})\), so \(V\) descends to classes and we write \(V_f := V(\mathfrak{b}_1)\); the second lattice has class \([\mathfrak{r} f]\).

**Lemma 1.1 (uniformization).** There is \(\tau_0 \in \mathbb{H}\) with \(\mathfrak{a} = \beta[\tau_0, 1]\) and \(\mathfrak{r}^{-1}\mathfrak{a} = \beta[\tau_0, \tfrac1{r_0}]\) for one common \(\beta\); consequently, with \(F := E_4^2E_6/\Delta\) (so \(h_2([\tau,1]) = C\,F(\tau)\), one universal constant),
$$
V(\mathfrak{a}) \;=\; \frac{1}{r_0^{2}}\,\frac{F(\tau_0)}{F(r_0\tau_0)}
\;=\; \frac{1}{r_0^{2}}\,\frac{j'(\tau_0)}{j'(r_0\tau_0)} ,
\qquad
\bigl(j(\tau_0),\, j(r_0\tau_0)\bigr) = (\beta_1, \beta_2) := \bigl(j(\mathfrak{a}),\, j(\mathfrak{r}^{-1}\mathfrak{a})\bigr).
$$

*Proof.* \(\mathcal{O}/\mathfrak{r} \cong \mathbb{Z}/r_0\) (as \(\mathfrak{r} \cap \mathbb{Z} = r_0\mathbb{Z}\)), so \(\mathfrak{r}^{-1}\mathfrak{a}/\mathfrak{a}\) is cyclic of order \(r_0\); pick a \(\mathbb{Z}\)-basis adapted to the elementary divisors \((1, r_0)\) and orient. Then \([\tau_0, \tfrac1{r_0}] = \tfrac1{r_0}[r_0\tau_0, 1]\) gives the displayed value; the \(2\pi i\)'s and \(C\)'s cancel in the ratio. \(\square\)

**Lemma 1.2 (the pair is a smooth point).** At every primitive class and every odd \(n \ge 3\), the point \((\beta_1, \beta_2)\) is a *smooth* point of the affine curve \(\Phi_{r_0} = 0\), and \(\Phi_x(\beta_1,\beta_2) \ne 0 \ne \Phi_y(\beta_1,\beta_2)\).

*Proof.* Branches of \(\{\Phi_{r_0} = 0\}\) through \((\beta_1,\beta_2)\) correspond to cyclic \(r_0\)-isogenies \(E_{\mathfrak{a}} \to E'\) with \(j(E') = \beta_2\) modulo automorphisms. Since \(j\) determines the lattice up to homothety, and singular moduli of distinct orders are distinct, every such isogeny has kernel \(E[\mathfrak{c}']\) for an **invertible** \(\mathcal{O}\)-ideal \(\mathfrak{c}'\) of norm \(r_0\) with \([\mathfrak{c}'] = [\mathfrak{r}]\). Such ideals correspond to primitive representations of \(r_0\) by the form \((r_0, 0, s_0)\) up to automorphs; since \(s_0 > r_0\), the only representation is \(r_0 = r_0\cdot 1^2 + s_0 \cdot 0^2\): the branch is **unique**. It is smooth and transverse to neither axis because its parametrization \(\tau \mapsto (j(\tau), j(r_0\tau))\) near \(\tau_0\) has \(j'(\tau_0) \ne 0 \ne j'(r_0\tau_0)\): both points are CM of discriminant \(D \notin \{-3, -4\}\) (both lattices are proper \(\mathcal{O}\)-ideals, as \(\mathfrak{r}\) is invertible), and \(j'\) vanishes only at the orbits of \(i\) and \(\rho\). A unique smooth branch means a smooth point of the (reduced, irreducible) curve, so \((\Phi_x, \Phi_y) \ne (0,0)\) there; finally \(\Phi_x\,j'(\tau_0) + \Phi_y\, r_0\,j'(r_0\tau_0) = 0\) (differentiate \(\Phi_{r_0}(j(\tau), j(r_0\tau)) = 0\)) shows that if one partial vanished, so would the other. \(\square\)

**Proposition 1.3 (correspondence-derivative formula).**
$$
V(\mathfrak{a}) \;=\; -\,\frac{1}{r_0}\cdot
\frac{\Phi_y(\beta_1, \beta_2)}{\Phi_x(\beta_1, \beta_2)}, \qquad \Phi = \Phi_{r_0}.
$$
In particular \(V(\mathfrak{a}) \in \mathbb{Q}(\beta_1, \beta_2) \subseteq H\): **the weight-2 ratio is a rational expression in the two singular moduli.**

*Proof.* Differentiate \(\Phi(j(\tau), j(r_0\tau)) \equiv 0\) at \(\tau_0\) and use Lemma 1.1 and \(\Phi_x \neq 0\). \(\square\)

The transcendence is now gone: the period, the \(2\pi i\), and — decisively — the sixth-root ambiguity of the \(\Delta\)-quotient route all cancel *before* any Galois argument is made.

## 2. Reciprocity at first power

**Proposition 2.1.** For every \(\sigma \in \mathrm{Gal}(\bar{\mathbb{Q}}/K)\) with Artin class \(\mathfrak{c} = \mathfrak{c}(\sigma) \in \mathrm{Cl}\) on \(H\):
\(\;\sigma(V(\mathfrak{a})) = V(\mathfrak{c}^{-1}\mathfrak{a})\). Complex conjugation gives \(\iota(V(\mathfrak{a})) = V(\bar{\mathfrak{a}})\). Hence, on classes, the full dihedral law \(\sigma(V_f) = V_{f^{e(\sigma)}\mathfrak{c}(\sigma)}\) for every \(\sigma \in \mathrm{Gal}(\bar{\mathbb{Q}}/\mathbb{Q})\), with \(e, \mathfrak{c}\) as in Theorem A.

*Proof.* Apply \(\sigma\) to the algebraic identity of Proposition 1.3. Since \(\Phi \in \mathbb{Z}[x,y]\),
\(\sigma(V(\mathfrak{a})) = -\tfrac1{r_0}\Phi_y(\sigma\beta_1, \sigma\beta_2)/\Phi_x(\sigma\beta_1, \sigma\beta_2)\), and by the classical main theorem for ring class fields of orders (the "classical inputs" of §5.5; Cox, *Primes of the form* \(x^2+ny^2\), Thm 11.36),
\(\sigma\beta_1 = j(\mathfrak{c}^{-1}\mathfrak{a})\) and \(\sigma\beta_2 = j(\mathfrak{c}^{-1}\mathfrak{r}^{-1}\mathfrak{a}) = j(\mathfrak{r}^{-1}(\mathfrak{c}^{-1}\mathfrak{a}))\): the image pair is the canonical pair of the translated class, which by Lemma 1.2 is again a smooth point evaluated by the same rational formula — i.e. it computes \(V(\mathfrak{c}^{-1}\mathfrak{a})\). For \(\iota\): \(\Phi\) has real (integer) coefficients, \(\iota(j(\mathfrak{a})) = j(\bar{\mathfrak{a}})\), and \(\bar{\mathfrak{r}} = \mathfrak{r}\) (ambiguous), so the image pair is the canonical pair of \(\bar{\mathfrak{a}}\). The dihedral composition is then formal, as in Theorem A. \(\square\)

Note what did *not* enter: no Kronecker limit formula, no elliptic units, no Siegel functions, no level structure beyond the \(j\)-line — only \(\sigma(j(\mathfrak{a})) = j(\mathfrak{c}^{-1}\mathfrak{a})\).

## 3. The unit \(\omega_f\), and its collapse to \(1\)

Since \([\mathfrak{b}_2] = [\mathfrak{r}][\mathfrak{b}_1] = [\mathfrak{r}^{-1}\mathfrak{b}_1]\), there is \(\nu = \nu_f \in K^\times\), unique up to sign, with
\(\mathfrak{b}_2 = \nu\,\mathfrak{r}^{-1}\mathfrak{b}_1\). Homogeneity turns the closed form into
$$
u_f \;=\; -\,\varepsilon\,\mu^{-2}\nu^2\; V_f \;=\; -\,r_0\,\omega_f\,V_f,
\qquad
\omega_f := \frac{\varepsilon\,\nu_f^{\,2}}{r_0\,\mu_f^{\,2}} .
$$
\(\omega_f\) is well-defined (only \(\nu^2\) enters) and canonical — \(u_f\) and \(V_f\) are class functions, so \(\omega_f = -u_f/(r_0V_f)\) does not depend on the representative used to compute \(\mu, \nu\).

**Lemma 3.1 (unit rigidity).** \(\omega_f \in \mu(B)\), the roots of unity of the CM field \(B\).

*Proof.* Exactly the \(\xi\)-torsion argument of §5.5. Archimedean: the Norm Lemma gives \(|\mu|^2 = \varepsilon q_1/q_2\), \(|\mu^\tau|^2 = q_1/(\varepsilon q_2)\), while \(N(\nu) = N\mathfrak{b}_2 \cdot N\mathfrak{r} / N\mathfrak{b}_1 = r_0q_1/q_2\) and \(|\nu^\tau| = |\bar\nu| = |\nu|\); so \(|\omega_f| = |\omega_f^\tau| = 1\) at both archimedean places. Ideal: \((\mu)\mathcal{O}_B = \mathfrak{m}_2\mathfrak{m}_1^{-1}\) (cocycle-ideal computation of §5.5) and \((\nu)\mathcal{O}_B = \mathfrak{m}_2\mathfrak{m}_1^{-1}(\mathfrak{r}\mathcal{O}_B)\) from the definition of \(\nu\), while \(\mathfrak{r}^2 = (r_0)\) and \(\varepsilon\) is a unit, so \((\omega_f) = (\mathfrak{r}^2/r_0)\mathcal{O}_B = \mathcal{O}_B\). Kronecker's theorem in the CM field \(B\) finishes. \(\square\)

By the Norm Lemma, \(\varpi := \mu/\nu\) satisfies \(\varpi\bar\varpi = \varepsilon/r_0\); so \(\omega_f = \varpi^{-2}\cdot(\varepsilon/r_0) = \bar\varpi/\varpi\), and
$$
\omega_f = 1 \iff \varpi = \mu/\nu \text{ is \emph{real}} \iff \mu\bar\nu \in \mathbb{R}.
$$

**Lemma 3.2 (coordinates).** Let \((w_1, w_2) = ((s_1,t_1),(s_2,t_2))\) be the oriented basis (\(t_1s_2 - s_1t_2 = a\)) of the congruence lattice \(\mathcal{K}\) from Lemma A of [class-formula-proof.md](class-formula-proof.md), so that the canonical \(X\) has bottom row \((c, d) = (-u_2, u_1)\), \(u_k = s_k + it_k\). Write \(\varsigma(s,t) := s\,s_0 + t\,\omega_0\) and \(\sigma_k := \varsigma(w_k)\), so that (Lemma C of class-formula-proof.md) \(\varsigma(\mathcal{K}) = \mathfrak{s}\,\mathfrak{a}_f\) with \(\mathfrak{a}_f = [a, b/2 + \omega_0]\). Then:

1. \(N \circ \varsigma = \tfrac{s_0}{2}\,g\) where \(g(s,t) = (n+1)s^2 + (n-1)t^2\) is the Gram form of Lemma B — one line: \(N(s\,s_0 + t\,\omega_0) = s_0^2s^2 + r_0s_0t^2 = \tfrac{s_0}{2}g(s,t)\).
2. Consequently the \(\sigma\)-circle Gram of Lemma B is the normalized norm form of the lattice \([\sigma_2, \sigma_1] = \mathfrak{s}\mathfrak{a}_f\), and its hyperbolic center is
   \(\;m_2 = \bar\sigma_1/\bar\sigma_2\) (the lower root: \(\operatorname{Im}(\sigma_1/\sigma_2) = \tfrac{a s_0\sqrt N/2}{N(\sigma_2)} > 0\) by the orientation, so the conjugate is the lower one).
3. \(\mathfrak{b}_2 = [1, -m_2] = \overline{[\sigma_2, \sigma_1]}/\bar\sigma_2 = \mathfrak{s}\,\bar{\mathfrak{a}}_f/\bar\sigma_2\), and with \(\mathfrak{b}_1 = \bar{\mathfrak{a}}_f/a\) and \(\mathfrak{s} = (\omega_0)\mathfrak{r}^{-1}\):
   \(\;(\nu) = \mathfrak{b}_2(\mathfrak{r}^{-1}\mathfrak{b}_1)^{-1} = (a\,\omega_0/\bar\sigma_2)\), i.e. \(\nu = \pm\,a\,\omega_0/\bar\sigma_2\).

All three identities are machine-verified exactly (rational arithmetic, no rounding) at every primitive class for all odd \(n \le 21\).

**Theorem 3.3 (\(\omega \equiv 1\)).** \(\omega_f = 1\) for every primitive class \(f\) and every odd \(n \ge 3\). Equivalently:
$$
\boxed{\,u_f \;=\; -\,r_0\,V_f \;=\; \frac{\Phi_y(\beta_1, \beta_2)}{\Phi_x(\beta_1, \beta_2)}\,,\qquad
\Phi = \Phi_{r_0},\ \ \beta_1 = j(\mathfrak{b}_1),\ \beta_2 = j(\mathfrak{r}^{-1}\mathfrak{b}_1). }
$$

*Proof.* \(\omega_f\) is basis-independent, so compute in the normal form produced by Lemma A's HNF basis, which always has \(s_2 = 0\) (the second basis vector spans \(\mathcal{K} \cap (0 \times \mathbb{Z})\)). Then \(u_2 = it_2\) and \(\sigma_2 = t_2\omega_0\) is purely imaginary, so \(\bar\sigma_2 = -\sigma_2\) and by Lemma 3.2.3, \(\nu = \pm a\omega_0/(t_2\omega_0) = \pm a/t_2 \in \mathbb{Q}\). For \(\mu\), Lemma 3.2.2 gives
$$
u_2 m_2 = it_2\cdot\frac{\bar\sigma_1}{-t_2\omega_0} = -\,\frac{i\bar\sigma_1}{\omega_0}
= -\,\frac{2\bar\sigma_1}{\sqrt N}
\quad\Longrightarrow\quad
\mu = u_1 - u_2m_2 = s_1 + it_1 + \frac{2(s_1s_0 - it_1\sqrt N/2)}{\sqrt N}
= s_1\,\frac{\sqrt N + 2s_0}{\sqrt N}:
$$
**the imaginary part cancels identically**, \(\mu \in \mathbb{Q}(\sqrt N)\) is real, and hence \(\varpi = \mu/\nu\) is real and \(\omega_f = 1\) — with the value confirmed by
\((\sqrt N + 2s_0)^2 = N + (n+1)^2 + 2(n+1)\sqrt N = 2(n+1)(n + \sqrt N) = 4s_0\varepsilon\), so
\(\mu^2 = \varepsilon s_1^2/r_0\) and
\(\omega_f = \varepsilon\nu^2/(r_0\mu^2) = a^2/(s_1t_2)^2 = 1\), since the orientation \(t_1s_2 - s_1t_2 = a\) with \(s_2 = 0\) forces \(s_1t_2 = -a\). \(\blacksquare\)

(Machine check: \(\omega_f = 1\) exactly — as an element of \(B\), by rational arithmetic — at every primitive class, all odd \(n \le 21\); and \(u_f = \Phi_y/\Phi_x(\beta_1,\beta_2)\) numerically to \(\ge 65\) digits at every class, \(n \le 13\).)

**Theorem 3.4 (first-power descent; \(\zeta \equiv 1\)).** For every \(\sigma \in \mathrm{Gal}(\bar{\mathbb{Q}}/\mathbb{Q})\) and every primitive class \(f\):
$$
\sigma(u_f) = u_{f^{e(\sigma)}\mathfrak{c}(\sigma)} .
$$

*Proof.* \(u_f = -r_0V_f\) (Theorem 3.3) with \(r_0 \in \mathbb{Q}\), and \(V\) obeys the dihedral law (Proposition 2.1). \(\blacksquare\)

**Corollaries.**
1. \(u_f \in H\), with \(\mathbb{Q}(u_f) \subseteq H\) cut out by the stabilizer \(\{\sigma : f^{e(\sigma)}\mathfrak{c}(\sigma) = f\}\).
2. Law 3 of §4 in full: \(S_x = u_f + u_f^{-1}\) satisfies \(\sigma(S_x) = S_{x^{e(\sigma)}\mathfrak{c}(\sigma)}\); the field statements of §5.8 (genus characters, \(\mathbb{Q}(\sqrt5)\) at \(n=9,11\), \(\mathbb{Q}(\sqrt{14})\) at \(n=13\), the single quartic orbit at \(n=15\)) are theorems.
3. Theorems A and B of §5.5/§5.8 are subsumed (take 12th powers), including at the \(\mu_8\)-levels \(n^2 - 1 = 2\square\): the sign allowance at \(u^6\) is never exercised.
4. The multiset \(\{u_f\}\) is stable under \(\mathrm{Gal}(\bar{\mathbb{Q}}/\mathbb{Q})\) and a **single Galois orbit**: translations \(\mathfrak{c}\) alone act transitively on classes, and every \(\mathfrak{c}\) is an Artin class (surjectivity of the Artin map).

## 4. Irreducibility of the level polynomials

**Theorem 4.1 (structure).** \(\Pi_n(x) = \prod_f (x - u_f) \in \mathbb{Q}[x]\), and \(\Pi_n = m^k\) where \(m\) is irreducible over \(\mathbb{Q}\) and \(k = h/\deg m\). Moreover \(k = 1\) — i.e. **\(\Pi_n\) is irreducible** — if and only if the \(u_f\) are pairwise distinct.

*Proof.* Rationality: the coefficients are algebraic and fixed by every \(\sigma\) (Corollary 4). Transitivity makes all irreducible factors equal (each root of each factor is conjugate to each root of any other), so \(\Pi_n = m^k\); the coincidence relation \(u_f = u_g\) is translation-invariant by Theorem 3.4, so it is the coset relation of a subgroup \(J \le \mathrm{Cl}\) with \(k = |J|\), and \(k = 1\) iff all values are distinct. \(\square\)

**Theorem 4.2 (the computed levels).** For \(n = 3, 5, 7, 9, 11, 13\) the polynomial \(\Pi_n\) equals (up to the stated integer normalization) the level polynomial of the §5.6 table, and it is **irreducible over \(\mathbb{Q}\)**, with pairwise distinct roots:

| \(n\) | \(\tilde P_n\) (primitive integer form) | degree | irreducible |
|---|---|---|---|
| 3 | \(x + 1\) | 1 | ✓ |
| 5 | \(6647x^2 + 30594194x + 6647\) | 2 | ✓ |
| 7 | \(11891x^2 + 80674200806x + 11891\) | 2 | ✓ |
| 9 | \(10565574794063311x^4 + 73919532109765731422845124x^3 - 118807282021266004510100774x^2 + (\mathrm{sym})\) | 4 | ✓ |
| 11 | \(76575720951x^4 + 466015525084217238173676x^3 - 216521978405797871634733654x^2 + (\mathrm{sym})\) | 4 | ✓ |
| 13 | \(722610532225x^4 + 3464286958371072692766958316x^3 + 4603575719671472165025576604518x^2 + (\mathrm{sym})\) | 4 | ✓ |

*Proof.* By Theorem 3.3, \(u_f = \Phi_y/\Phi_x(\beta_1, \beta_2)\) with \(\Phi = \Phi_{r_0}\) and \((\beta_1, \beta_2)\) the canonical singular-moduli pair. The computation is then exact rational arithmetic, with **no numerical identification anywhere**:
(i) \(H_D \in \mathbb{Z}[t]\) is the Hilbert class polynomial (integer coefficients recovered with error \(< 10^{-40}\), i.e. rigorously; irreducible over \(\mathbb{Q}\), classical);
(ii) in \(F_1 = \mathbb{Q}[t]/(H_D)\), the pairing root \(\beta_2(t)\) is the root of \(\gcd\bigl(H_D(y), \Phi_{r_0}(t, y)\bigr)\), which has degree 1 — the exact-arithmetic reflection of the uniqueness in Lemma 1.2 (at these levels every prime of \(r_0\) is ramified, so \(\mathfrak{r}\) is the *only* invertible ideal of norm \(r_0\) and no other \(H_D\)-root can pair with \(\beta_1\));
(iii) \(u(t) = \Phi_y(t, \beta_2(t))/\Phi_x(t, \beta_2(t)) \in F_1\), and \(\Pi_n = \) the characteristic polynomial of multiplication by \(u(t)\) on \(F_1\) — this is \(\prod_{H_D(\theta)=0}(x - u(\theta)) = \prod_f (x - u_f)\), since \(\theta \mapsto\) class is the standard bijection and \(u(\theta)\) is the class's invariant by Theorem 3.3;
(iv) the resulting integer polynomials are exactly the table's, are irreducible by exact factorization, and are squarefree (\(\gcd(\Pi_n, \Pi_n') = 1\)), which *proves* the pairwise distinctness required by Theorem 4.1. The modular polynomials \(\Phi_2, \dots, \Phi_6\) were themselves constructed exactly from integer \(q\)-expansions with cyclotomic coefficients (a holomorphic \(\mathrm{SL}_2(\mathbb{Z})\)-invariant function that is \(O(q)\) vanishes, so the peeling against powers of \(j\) is an identity, verified with slack), validated against the classical \(\Phi_2\) and by the symmetry \(\Phi_m(x,y) = \Phi_m(y,x)\). \(\blacksquare\)

**Remark (general \(n\)).** The only unproved hypothesis for general \(n\) is distinctness (\(J = 1\)). Two general facts: \(\mathfrak{r} \notin J\) for \(n \ge 5\) (else \(u_f^2 = 1\) for all \(f\) by law 2, contradicting \(|u_{(1,0,\ast)}| = r_0^2|F(\tau_1)/F(r_0\tau_1)| > 1\), which follows from crude \(q\)-expansion bounds since \(\operatorname{Im}\tau_1 = \tfrac{\sqrt N}{2} > r_0\operatorname{Im}\tau_1 / \ldots\) — the principal-class value is exponentially large); and \(J\) is constrained by laws 1–2 to be \(\hat\sigma\)-stable. Every computed level has \(J = 1\).

**Remark (the certified record is now a theorem).** The §5.6 table was *certified* by safe-parameter PSLQ; the computation above re-derives every entry from \(H_D\) and \(\Phi_{r_0}\) by exact arithmetic. The two independent routes agree verbatim, which simultaneously proves the table and validates the certification methodology.

## 5. What this closes, and what it opens

- **Closed:** open problem 1 of §5.8 (\(\zeta \equiv 1\), "Kubert–Lang multiplier computation") — by the isogeny-differential route of outlook §2.1, in the derivative-of-\(\Phi\) form. The Siegel-unit calculus was never needed: the sixth root that plagued the \(\Delta\)-quotient expression of \(u^6\) is an artifact of factoring \(h_2^6\); the ratio \(h_2/h_2\) itself is already a rational function on \(X_0(r_0)\).
- **New closed form.** \(u_f = \Phi_y/\Phi_x(\beta_1, \beta_2)\) gives a purely algebraic, PSLQ-free definition of the phase, and an exact algorithm for \(\Pi_n\) at any level (cost: \(H_D\) plus \(\Phi_{r_0}\)).
- **Moduli-theoretic meaning.** \(-1/(r_0 V_f) = -\Phi_x/(r_0\Phi_y)\) is \(dj_2/dj_1\) along \(X_0(r_0)\): the phase invariant of a Schmidt circle *is* the derivative of the modular correspondence at the Heegner point it defines — e.g. law 2, \(u_{\mathfrak{r}f}u_f = 1\), is the chain rule for the two branches through the point, and \(V_fV'_f\)-type products over the partner ideal \(\mathfrak{s}\) recover \(V_f/V'_f = -s_0/r_0\) from \(\mathfrak{r}\mathfrak{s} = (\omega_0)\).
- **Open (unchanged):** exact denominator valuations (§5.7); the distinctness subgroup \(J = 1\) for general \(n\); the sign law of outlook §1.1 — note the sign data now lives entirely in \(V_f\), i.e. in the geometry of \(\Phi_{r_0}\) at real points.

## 6. Files

- [scripts/first_power_descent.py](scripts/first_power_descent.py) — (a) exact construction of \(\Phi_m\), \(m \le 10\), with the \(\Phi_2\)/symmetry validations; (b) exact verification \(\omega_f = 1\) at every primitive class, odd \(n \le 21\); (c) numerical verification of \(u_f = \Phi_y/\Phi_x(\beta_1,\beta_2)\) to 65+ digits, \(n \le 13\); (d) the exact computation of \(\Pi_n\) (Theorem 4.2) with irreducibility and table match. Run: `python3 scripts/first_power_descent.py all`.
