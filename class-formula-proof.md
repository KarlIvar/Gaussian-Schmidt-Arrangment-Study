# Proof of the class formula \(\hat\sigma[f] = [\mathfrak{r}_n]\,[f]^{-1}\)

This document proves the main theorem of [involution.md](involution.md) §4. Throughout, \(n \ge 3\) is odd, \(D = 1 - n^2 \equiv 0 \pmod 8\), and

$$
\sigma(X) = \bar X^{-1}, \qquad
M_0 = \begin{pmatrix} 0 & i \\ -i & 0\end{pmatrix}, \qquad
h(u, v) = u^\dagger M_0 v = i(\bar u_1 v_2 - \bar u_2 v_1),
$$

so that \(h(u,u) = 2\operatorname{Im}(u_1 \bar u_2)\), and the circle \(V(\hat{\mathbb{R}})\) has Hermitian matrix \((V^{-1})^\dagger M_0 V^{-1}\). A positive definite form \(f = (a, b, c)\) of discriminant \(D = b^2 - 4ac\) (\(b\) is automatically even) corresponds to the circle
$$
M_f = \begin{pmatrix} 2a & b - ni \\ b + ni & 2c \end{pmatrix}
\qquad \text{(curvature } 2a, \text{ center } \tfrac{-b + ni}{2a}\text{)} .
$$

> **Theorem.** Let \(f\) be a primitive positive definite form of discriminant \(1 - n^2\) and \(X \in \mathrm{SL}_2(\mathbb{Z}[i])\) any matrix whose circle \(X(\hat{\mathbb{R}})\) lies in the \(\mathrm{SL}_2(\mathbb{Z})\)-class of \(M_f\). Then the circle of \(\sigma(X)\) lies in the lower half-plane, with reversed orientation, and its reflection \(z \mapsto \bar z\) back into the upper half-plane lies in the class of
> $$ [\mathfrak{r}_n]\cdot[f]^{-1}, \qquad \mathfrak{r}_n = \bigl[\bigl(\tfrac{n-1}{2},\, 0,\, \tfrac{n+1}{2}\bigr)\bigr], $$
> product and inverse taken in the form class group of discriminant \(1-n^2\).

Every identity and sign below is additionally machine-verified, class by class, for all odd \(n \le 41\) ([scripts/proof_check.py](scripts/proof_check.py)).

## 0. Conventions: the ideal dictionary

Let \(K = \mathbb{Q}(\sqrt D)\) with \(\sqrt D := i\sqrt{n^2-1}\), \(\omega = \tfrac{\sqrt D}{2}\), \(\mathcal{O} = \mathbb{Z}[\omega]\) the order of discriminant \(D\), and \(N(x + y\omega) = x^2 + \tfrac{n^2-1}{4}y^2\). For a proper (invertible) \(\mathcal{O}\)-ideal \(\mathfrak{b}\) with \(\mathbb{Z}\)-basis \((\alpha, \beta)\), call the basis *positively oriented* if \(\operatorname{Im}(\bar\alpha\beta) > 0\), and set \(f_{\mathfrak{b}}(s, t) = N(s\alpha + t\beta)/N(\mathfrak{b})\). Classically (Gauss; see Cox, *Primes of the form \(x^2+ny^2\)*, §7):

- \([\mathfrak{b}] \mapsto [f_\mathfrak{b}]\) is an isomorphism from the proper ideal class group of \(\mathcal{O}\) onto the primitive form class group under composition; a negatively oriented basis produces the opposite (= inverse) class;
- norms are multiplicative on proper ideals: \([\mathfrak{b} : \mathfrak{b}\mathfrak{a}] = N(\mathfrak{a})\);
- in this convention \(\mathfrak{a}_f := \mathbb{Z}a + \mathbb{Z}\bigl(\tfrac b2 + \omega\bigr)\) satisfies \(f_{\mathfrak{a}_f} = f\): indeed \(N\bigl(sa + t(\tfrac b2 + \omega)\bigr) = (sa + \tfrac b2 t)^2 + \tfrac{n^2-1}{4}t^2 = a\,f(s,t)\), using \(\tfrac{b^2 + n^2 - 1}{4} = ac\); the basis is positive since \(\operatorname{Im}(\bar a(\tfrac b2+\omega)) = a\operatorname{Im}\omega > 0\);
- likewise \(\mathfrak{s} := \mathbb{Z}\tfrac{n+1}{2} + \mathbb{Z}\omega\) has \(f_\mathfrak{s} = \bigl(\tfrac{n+1}{2}, 0, \tfrac{n-1}{2}\bigr)\), so \([\mathfrak{s}] = \mathfrak{r}_n\); and \(\bar{\mathfrak{s}} = \mathfrak{s}\), so \([\mathfrak{s}]^2 = [\mathfrak{s}\bar{\mathfrak{s}}] = [(N\mathfrak{s})] = 1\): \(\mathfrak{r}_n\) is \(2\)-torsion (as it must be, for \(\hat\sigma\) to be an involution).

Since \(\hat\sigma\) is well defined on double cosets \(\Gamma X \Gamma\), \(\Gamma = \mathrm{SL}_2(\mathbb{Z})\) ([involution.md](involution.md) §4), it suffices to prove the theorem for **one** \(X\) per class. We construct one explicitly.

## 1. An explicit unitary basis for each class

**Lemma A.** Let \(f = (a,b,c)\), \(D = 1-n^2\), and let
$$
\mathcal{K} := \Bigl\{ (s,t) \in \mathbb{Z}^2 \;:\; (n+1)s \equiv bt, \quad bs \equiv -(n-1)t \pmod{2a} \Bigr\}.
$$
Then \(\mathcal{K}\) has index \(a\) in \(\mathbb{Z}^2\). For any basis \(w_1 = (s_1,t_1)\), \(w_2 = (s_2,t_2)\) of \(\mathcal{K}\) with \(t_1s_2 - s_1t_2 = a\), put \(u_k = s_k + it_k\) and
$$
v_k = \frac{(n + bi)\,u_k + \bar u_k}{2ia} \qquad (k = 1,2).
$$
Then \(P = \begin{pmatrix} u_1 & v_1 \\ u_2 & v_2 \end{pmatrix} \in \mathrm{SL}_2(\mathbb{Z}[i])\) and \(P^\dagger M_0 P = M_f\); equivalently, \(X := P^{-1}\) satisfies \(X(\hat{\mathbb{R}}) = \) the circle \(M_f\).

*Proof.* **Index.** \(\mathcal{K} = \ker\bigl(R \bmod 2a\bigr)\) for \(R = \begin{pmatrix} n+1 & -b \\ b & n-1 \end{pmatrix}\), \(\det R = n^2 - 1 + b^2 = 4ac\). The Smith invariants of \(R\) are \((2, 2ac)\): all entries are even and \(\gcd\bigl(\tfrac{n+1}{2}, \tfrac{n-1}{2}\bigr) = 1\), so the first invariant is exactly \(2\). Hence \(\#\ker(R \bmod 2a) = \gcd(2, 2a)\gcd(2ac, 2a) = 2 \cdot 2a\), and \([\mathbb{Z}^2 : \mathcal{K}] = (2a)^2/4a = a\).

**Integrality of \(v_k\).** \((n+bi)(s+it) + (s-it) = \bigl[(n+1)s - bt\bigr] + i\bigl[bs + (n-1)t\bigr]\), so \(2ia \mid (n+bi)u_k + \bar u_k\) iff \((s_k, t_k) \in \mathcal{K}\).

**Gram identities.** Write \(z := u_1\bar u_2\); then \(\operatorname{Im} z = t_1s_2 - s_1t_2 = a\). Now:
- \(h(p_1, p_1) = 2\operatorname{Im}(u_1\bar u_2) = 2a\).
- \(\det P = u_1v_2 - u_2v_1 = \dfrac{u_1\bar u_2 - u_2 \bar u_1}{2ia} = \dfrac{2i\operatorname{Im}z}{2ia} = 1\).
- \(\bar u_1 v_2 - \bar u_2 v_1 = \dfrac{(n+bi)(\bar u_1 u_2 - \bar u_2 u_1)}{2ia} = \dfrac{(n+bi)(-2ia)}{2ia} = -(n+bi)\), so \(h(p_1,p_2) = i(\bar u_1v_2 - \bar u_2 v_1) = b - ni\). ✓
- \(v_1 \bar v_2 = \dfrac{\bigl[(n+bi)u_1 + \bar u_1\bigr]\bigl[(n - bi)\bar u_2 + u_2\bigr]}{4a^2}\); the cross terms \((n+bi)u_1u_2 + (n-bi)\bar u_1\bar u_2\) are conjugate, hence real, so
$$\operatorname{Im}(v_1\bar v_2) = \frac{(n^2+b^2)\operatorname{Im}z + \operatorname{Im}\bar z}{4a^2} = \frac{(n^2 + b^2 - 1)a}{4a^2} = \frac{4ac}{4a} = c,$$
i.e. \(h(p_2,p_2) = 2c\). \(\square\)

(For example \(n = 3\), \(f = (1,0,2)\): \(\mathcal{K} = \mathbb{Z}^2\), \(w_1 = (0,1)\), \(w_2 = (1,0)\) gives \(P = \begin{pmatrix} i & 1 \\ 1 & -2i \end{pmatrix}\).)

## 2. The Gram matrix of the \(\sigma\)-image

Let \(X = P^{-1}\), so \(\sigma(X) = \bar P\), and (as in [involution.md](involution.md) §4) the image circle is
\(M_{\sigma X} = -\overline{N}\) with \(N := (P^{-1})^\dagger M_0 P^{-1}\), the Gram matrix of \(h\) on the columns \((v_2, -u_2)^{\mathsf T}\), \((-v_1, u_1)^{\mathsf T}\) of \(P^{-1}\).

**Lemma B.** With \(g(s,t) := (n+1)s^2 + (n-1)t^2\) and polarization \(g(w, w') = (n+1)ss' + (n-1)tt'\):
$$
N = \begin{pmatrix} \dfrac{g(w_2)}{a} & -\dfrac{g(w_1,w_2)}{a} - ni \\[2mm] -\dfrac{g(w_1,w_2)}{a} + ni & \dfrac{g(w_1)}{a} \end{pmatrix},
$$
and all displayed rational numbers are integers, the diagonal ones even.

*Proof.* Using \(\tfrac{1}{2ia} = \tfrac{-i}{2a}\) and \(\operatorname{Im}(-i\xi) = -\operatorname{Re}\xi\):
$$
N_{11} = -2\operatorname{Im}(v_2\bar u_2) = \frac{\operatorname{Re}\bigl[(n+bi)|u_2|^2 + \bar u_2^{\,2}\bigr]}{a} = \frac{n|u_2|^2 + \operatorname{Re}(u_2^2)}{a} = \frac{(n+1)s_2^2 + (n-1)t_2^2}{a},
$$
and symmetrically \(N_{22} = g(w_1)/a\). For the off-diagonal entry,
$$
u_1\bar v_2 - \bar u_2 v_1 = \frac{i}{2a}\Bigl[(n-bi)u_1\bar u_2 + (n+bi)u_1\bar u_2 + u_1u_2 + \bar u_1\bar u_2\Bigr] = \frac{i}{a}\bigl[n\,u_1\bar u_2 + \operatorname{Re}(u_1u_2)\bigr],
$$
and with \(u_1\bar u_2 = (s_1s_2 + t_1t_2) + ia\), \(\operatorname{Re}(u_1u_2) = s_1s_2 - t_1t_2\):
$$
N_{12} = i\,(u_1\bar v_2 - \bar u_2 v_1) = -\frac{(n+1)s_1s_2 + (n-1)t_1t_2}{a} - ni .
$$
*Integrality:* multiplying the two congruences defining \(\mathcal{K}\) by \(s\) resp. \(t\) gives \((n+1)s^2 \equiv bst\) and \((n-1)t^2 \equiv -bst \pmod{2a}\), so \(g(w) \equiv 0 \pmod{2a}\) on \(\mathcal{K}\); polarizing, \(2g(w_1,w_2) = g(w_1{+}w_2) - g(w_1) - g(w_2) \equiv 0 \pmod{2a}\). \(\square\)

**Consequences.** The oriented image \(M_{\sigma X} = -\bar N\) has curvature entry \(-g(w_2)/a < 0\) (\(g\) is positive definite): **\(\hat\sigma\) always reverses orientation**, and the unoriented image lies in the lower half-plane with \(|y| = n\): **\(\alpha\) is preserved** — both observed facts are now proved. Reflecting by \(z \mapsto \bar z\) returns the form read from \(\operatorname{Re}N\): in the basis \((w_2, -w_1)\) of \(\mathcal{K}\) (which has the same orientation as \((w_1, w_2)\)),
$$
f' := \Bigl(\frac{g(w_2)}{2a},\; -\frac{g(w_1,w_2)}{a},\; \frac{g(w_1)}{2a}\Bigr) = \frac{1}{2a}\, g\big|_{\mathcal{K}} ,
$$
a form of discriminant \(\operatorname{disc}(g)\cdot a^2/(2a)^2 = 4D\cdot a^2/4a^2 = D\). **So the image class is the \(\mathfrak{r}_n\)-form \(\tfrac12 g = \bigl(\tfrac{n+1}{2}, 0, \tfrac{n-1}{2}\bigr)\), restricted to the index-\(a\) lattice \(\mathcal{K}\) and divided by \(a\).** It remains to recognize this as composition.

## 3. \(\mathcal{K}\) is an ideal, and it factors

Embed \(\iota : \mathbb{Z}^2 \hookrightarrow K\), \(\iota(s,t) = s\,\tfrac{n+1}{2} + t\,\omega\), so \(\iota(\mathbb{Z}^2) = \mathfrak{s}\) and
$$
N(\iota(w)) = \tfrac{(n+1)^2}{4}s^2 + \tfrac{n^2-1}{4}t^2 = \tfrac{n+1}{2}\cdot\tfrac{1}{2}g(w),
\qquad\text{hence}\qquad
f' = \frac{N \circ \iota}{a \cdot \frac{n+1}{2}}\Big|_{\mathcal{K}} .
$$

**Lemma C.** \(\mathfrak{c} := \iota(\mathcal{K})\) is an \(\mathcal{O}\)-ideal, and for \(f\) primitive, \(\mathfrak{c} = \mathfrak{s}\,\mathfrak{a}_f\).

*Proof.* **\(\omega\)-stability.** \(\omega\,\iota(s,t) = \tfrac{D}{4}t + \tfrac{n+1}{2}s\,\omega = \iota\bigl(\tfrac{1-n}{2}t,\ \tfrac{n+1}{2}s\bigr)\). The image pair satisfies the \(\mathcal{K}\)-congruences: substituting,
$$
(n+1)\tfrac{1-n}{2}t - b\tfrac{n+1}{2}s = -\tfrac{n+1}{2}\bigl[bs + (n-1)t\bigr], \qquad
b\tfrac{1-n}{2}t + (n-1)\tfrac{n+1}{2}s = \tfrac{n-1}{2}\bigl[(n+1)s - bt\bigr],
$$
both multiples of the congruences defining \(\mathcal{K}\), hence \(\equiv 0 \pmod{2a}\). So \(\mathcal{K}\) is stable under the (integral) action of \(\omega\), and \(\mathfrak{c}\) is an ideal.

**Containment \(\mathfrak{s}\mathfrak{a}_f \subseteq \mathfrak{c}\).** The four generators of \(\mathfrak{s}\mathfrak{a}_f\) have \(\iota\)-coordinates
$$
\tfrac{n+1}{2}\cdot a \leftrightarrow (a, 0), \quad
\tfrac{n+1}{2}\bigl(\tfrac b2 + \omega\bigr) \leftrightarrow \bigl(\tfrac b2, \tfrac{n+1}{2}\bigr), \quad
a\,\omega \leftrightarrow (0, a), \quad
\bigl(\tfrac b2+\omega\bigr)\omega \leftrightarrow \bigl(\tfrac{1-n}{2}, \tfrac b2\bigr),
$$
and each lies in \(\mathcal{K}\): for the first and third use that \(b\) and \(n \pm 1\) are even; for the second, \((n+1)\tfrac b2 - b\tfrac{n+1}{2} = 0\) and \(b\tfrac b2 + (n-1)\tfrac{n+1}{2} = \tfrac{b^2 + n^2 - 1}{2} = 2ac \equiv 0 \pmod {2a}\); for the fourth, \((n+1)\tfrac{1-n}{2} - b\tfrac b2 = -2ac\) and \(b\tfrac{1-n}{2} + (n-1)\tfrac b2 = 0\).

**Equality.** Both sides sit inside \(\mathfrak{s}\) with the same index: \([\mathfrak{s} : \mathfrak{c}] = [\mathbb{Z}^2 : \mathcal{K}] = a\) (Lemma A), and — here primitivity of \(f\) enters — \(\mathfrak{a}_f\) is proper, so \([\mathfrak{s} : \mathfrak{s}\mathfrak{a}_f] = N(\mathfrak{a}_f) = a\). A containment of equal finite index is an equality. \(\square\)

## 4. Conclusion: orientation and the twist

By Lemma C and the displayed norm identity, \(f'\) is the norm form of the proper ideal \(\mathfrak{c} = \mathfrak{s}\mathfrak{a}_f\), of norm \(N(\mathfrak{c}) = \tfrac{n+1}{2}a\), computed in the basis \(\bigl(\iota(w_2), \iota(-w_1)\bigr)\). This basis is **negatively** oriented: for \(w = (s,t)\), \(w' = (s',t')\),
$$
\operatorname{Im}\bigl(\overline{\iota(w)}\,\iota(w')\bigr) = \tfrac{n+1}{2}\cdot\tfrac{\sqrt{n^2-1}}{2}\,(st' - ts'),
$$
and for the pair \((w_2, -w_1)\): \(s_2(-t_1) - t_2(-s_1) = -(t_1s_2 - s_1t_2) = -a < 0\). A negatively oriented basis yields the opposite class, so by the dictionary of §0:
$$
[f'] \;=\; [f_{\mathfrak{c}}]^{-1} \;=\; \bigl([\mathfrak{s}][\mathfrak{a}_f]\bigr)^{-1} \;=\; [\mathfrak{r}_n]^{-1}\,[f]^{-1} \;=\; [\mathfrak{r}_n]\,[f]^{-1},
$$
the last step because \(\mathfrak{r}_n\) is \(2\)-torsion (\(\bar{\mathfrak{s}} = \mathfrak{s}\)). \(\blacksquare\)

Note where each ingredient of the formula comes from:
- the **inversion** \([f]^{-1}\) comes from the orientation flip of complex conjugation inside \(\sigma\) (the \(-a < 0\) computation);
- the **twist** \(\mathfrak{r}_n\) comes from the quadratic form \(g(s,t) = (n+1)s^2 + (n-1)t^2 = n|u|^2 + \operatorname{Re}(u^2)\) — the norm form twisted by Galois conjugation \(u \mapsto \bar u\), which is exactly the fingerprint \(\sigma\) leaves on the Gram matrix (Lemma B); its class is forced to be \(\bigl(\tfrac{n+1}{2}, 0, \tfrac{n-1}{2}\bigr)\), i.e. the split of \(n^2 - 1\) into the two even parts \(\tfrac{n\pm1}{2}\cdot 2\). (Heuristically: \(g\) is \(2N_{B/K'}\)-shaped for the biquadratic field \(B = \mathbb{Q}(i, \sqrt{n^2-1})\), where \(\sqrt{\varepsilon} = \sqrt{\tfrac{n+1}{2}} + \sqrt{\tfrac{n-1}{2}}\) for the unit \(\varepsilon = n + \sqrt{n^2-1}\) — the same numbers \(\tfrac{n\pm1}{2}\) that define \(\mathfrak{r}_n\).)

## 5. Byproducts and remarks

1. **Intrinsic form of the trace identity.** On the lattice \(\Lambda = P\mathbb{Z}^2\) (Gram \(M_f\)), the Hermitian form decomposes as \(h|_\Lambda = 2f - ni\,\langle\cdot,\cdot\rangle\) (\(\langle u,v\rangle = u^{\mathsf T}Jv\) restricts to the standard unimodular symplectic form on \(\Lambda\)). Writing \(Z = P\bar P^{-1}\) (so \(\bar\lambda = Z^{-1}\lambda\) on \(\Lambda\)) one finds
$$
f(\lambda, \mu) = \bigl\langle \lambda,\; W\mu \bigr\rangle, \qquad W = \tfrac{i}{4}\,(Z - Z^{-1}),
$$
with \(\operatorname{tr}W = 0\); since \(Z + Z^{-1} = (\operatorname{tr}Z)I\), Cayley–Hamilton gives \(W^2 = \tfrac{4 - (\operatorname{tr}Z)^2}{16}I\), while \(\det\)-considerations force \(W^2 = \tfrac{D}{4}I\); equating, \(\operatorname{tr}(P\bar P^{-1}) = \pm 2n\) — an intrinsic re-proof of the trace/\(\alpha\) identity, exhibiting \(\Lambda\) as a module over \(\mathcal{O} = \mathbb{Z}[W]\) (the parity \(b \equiv 0 \bmod 2\), i.e. the Schmidt congruence \(\zeta \equiv i \bmod 2\), is exactly the integrality of \(W\)).
2. **Imprimitive classes.** The only step using primitivity is the index count \([\mathfrak{s} : \mathfrak{s}\mathfrak{a}_f] = N(\mathfrak{a}_f)\) in Lemma C. For \(f\) of content \(g\), \(\mathfrak{a}_f\) is not proper, \(\mathfrak{s}\mathfrak{a}_f \subsetneq \mathfrak{c}\) is possible, and the image twist becomes the appropriate image of \(\mathfrak{r}_n\) on the content stratum — consistent with the data in [involution.md](involution.md) (stratum-consistent ambiguous twists, verified for \(n \le 41\)); a clean statement at the ramified prime \(2\) is left open.
3. **A closed-form section.** Lemma A is of independent computational value: it writes down, in closed form, an element of \(\mathrm{SL}_2(\mathbb{Z}[i])\) realizing any prescribed Schmidt circle (previously done by the descent algorithm).
4. **Verification.** [scripts/proof_check.py](scripts/proof_check.py) checks, for every primitive class of every odd \(n \le 41\) (169 classes): Lemma A (\(P \in \mathrm{SL}_2(\mathbb{Z}[i])\), Gram \(= M_f\)); Lemma B (entrywise, including the mod-\(2a\) integrality); Lemma C (\(\omega\)-stability of \(\mathcal{K}\) and the lattice equality \(\mathfrak{c} = \mathfrak{s}\mathfrak{a}_f\) via Hermite normal forms); and the final class identity \([f'] = [\mathfrak{r}_n][f]^{-1}\) via Gauss composition.
