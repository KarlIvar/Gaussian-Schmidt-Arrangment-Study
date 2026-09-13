# The full Robert index of the Schmidt \(\Delta\)-units

This document carries out item 2 of [PROGRAM.md](PROGRAM.md): it determines the
index of the group of Schmidt \(\Delta\)-units in the full unit group of the ring
class field, for **every** Euclidean level \(n \ge 2\), and thereby replaces the
conjecture of Paper II (`papers/2-schmidt-elliptic-units`, Conjecture 6.7) by a
theorem. The objects are those of Paper II, Theorem 2.6, and of
[schmidt-euler-system.md](schmidt-euler-system.md), Corollary 0.1: the
\(\Delta\)-data
$$
G_\mathfrak{c} \;=\; n^{12}\,\frac{\Delta(\Lambda_\mathfrak{c})}{\Delta(\mathbb{Z}[i])},
\qquad \mathfrak{c} \in \mathrm{Pic}(\mathcal{O}_n),\quad \mathcal{O}_n = \mathbb{Z} + n\mathbb{Z}[i],
$$
algebraic integers of the ring class field \(H_n\) of \(K = \mathbb{Q}(i)\), and the group
of their ratios
$$
\mathcal{V}_n \;:=\; \bigl\langle G_\mathfrak{c}/G_{\mathfrak{c}'}\bigr\rangle
\;=\; \bigl\langle G_\mathfrak{c}/G_1 : \mathfrak{c} \ne 1\bigr\rangle \;\subset\; \mathcal{O}_{H_n}^\times ,
$$
a Galois-stable group of units (units by the per-class valuation law, Paper II,
Theorem 4.11). The result:

> **Theorem 1 (the full Robert index).** Let \(n \ge 2\), \(h = |\mathrm{Pic}(\mathcal{O}_n)| = h(-4n^2)\),
> \(w_{H_n} = |\mu(H_n)|\), \(h_{H_n}\) and \(R_{H_n}\) the class number and regulator of \(H_n\),
> and for a nontrivial character \(\chi\) of \(\mathrm{Pic}(\mathcal{O}_n)\) let \(L(s,\chi)\) be the
> Epstein \(L\)-function of the level (Paper II, Definition 2.1) and
> \(C_\chi(0) = L'(0,\chi)/L'_{\mathrm{prim}}(0,\chi)\) its imprimitivity multiplier
> (\(=1\) when \(\chi\) factors through no proper divisor level). Then \(\mathcal{V}_n\) is free
> abelian of rank \(h-1 = \operatorname{rank}\mathcal{O}_{H_n}^\times\), meets \(\mu(H_n)\) trivially, and
> $$
> \boxed{\;\bigl[\mathcal{O}_{H_n}^\times : \mu(H_n)\,\mathcal{V}_n\bigr]
> \;=\; \frac{24^{\,h-1}\prod_{\chi\ne1}|L'(0,\chi)|}{R_{H_n}}
> \;=\; \frac{4\cdot 24^{\,h-1}}{w_{H_n}}\; h_{H_n}\prod_{\chi\ne1}C_\chi(0)\;}
> $$
> with \(\prod_{\chi\ne1}C_\chi(0)\) a positive integer, and \(w_{H_n} = 4\cdot 2^{[4\mid n]}\cdot 3^{[3\mid n]}\).
> Equivalently \([\mathcal{O}_{H_n}^\times : \mathcal{V}_n] = 4\cdot 24^{h-1}\,h_{H_n}\prod_\chi C_\chi(0)\).

The Robert index of the Schmidt \(\Delta\)-units is the class number of the ring class
field, times an explicit power of \(2\) and \(3\), times the Euler multipliers of the
imprimitive characters — the Kubert–Lang shape that Paper II conjectured, now for
every level. The same argument gives the index in every subfield of \(H_n\) cut out
by a subgroup of \(\mathrm{Pic}(\mathcal{O}_n)\) (Theorem 2), which contains the cubic layer of
Paper II and settles the sextic layer; the hyperbolic units \(R_f\) of Paper II have
the analogous index in the odd units of their class field, with a genuine \(2\)-adic
invariant \(Q^-\) (Theorem 3, §5.6; \(= 24^6\) exactly at \(n = 21\)); and the exponent matrices of the units in
PARI's fundamental-unit bases show, at every computed level, that modulo the
\(24\)-th-root saturation of \(\mathcal{V}_n\) the index is *exactly*
\(h_{H_n}\prod C_\chi(0)\) (§5, certified; conjectured in general).

Everything finite is verified by
[scripts/robert_index_full.py](scripts/robert_index_full.py) (`--selftest`), which
computes the Dedekind determinant from the \(h\) numbers \(\log|G_\mathfrak{c}|\), the
character sums, the independent incomplete-gamma evaluation of \(L'(0,\chi)\), the
class groups and multipliers exactly, and the class numbers, regulators, torsion
and fundamental units of \(H_n\) with PARI/GP 2.15.4 — `bnfcertify` succeeded at
every level \(n \le 15\), so those entries are **unconditional**; the level \(n = 23\)
(degree 24) is GRH-conditional (§6); the hyperbolic levels \(5 \le n \le 21\) are
unconditional except the degree-24 field at \(n = 21\). The numerics verify the proof; the proof does
not rest on them.

Notation as in Paper II and the Euler document: \(K = \mathbb{Q}(i)\), \(\chi_{-4}\) its
character, \(N_e(n) = n\prod_{p\mid n}(1-\chi_{-4}(p)/p)\), \(h = \tfrac12N_e(n)\) for \(n \ge 2\);
\(\Lambda_\mathfrak{c} \subseteq \mathbb{Z}[i]\) the primitive index-\(n\) sublattice of the class
\(\mathfrak{c}\); \(M(n) = \prod_\mathfrak{c}G_\mathfrak{c}\) the \(\Delta\)-mass; for \(\sigma \in \mathrm{Gal}(H_n/K)\)
with Artin class \(\mathfrak{a} = \mathfrak{c}(\sigma)\) we write \(\sigma = \sigma_\mathfrak{a}\), so that
(Paper II, Theorem 2.6(1))
$$
\sigma_\mathfrak{a}(G_\mathfrak{c}) = G_{\mathfrak{a}^{-1}\mathfrak{c}}, \qquad
\overline{G_\mathfrak{c}} = G_{\mathfrak{c}^{-1}}, \qquad
\sum_\mathfrak{c}\chi(\mathfrak{c})\log|G_\mathfrak{c}| = -12\,L'(0,\chi)\quad(\chi \ne 1).
$$
Throughout, \(f(\mathfrak{c}) := \log|G_\mathfrak{c}|\) and \(S_\chi := \sum_\mathfrak{c}\chi(\mathfrak{c})f(\mathfrak{c})\).

## 1. Places, the regulator, and the index as a ratio of covolumes

**1.1 The field.** \(H_n \supset K\) is totally complex of degree \(2h\), Galois over
\(\mathbb{Q}\) with group \(\mathrm{Pic}(\mathcal{O}_n)\rtimes\langle\iota\rangle\), \(\iota\) acting by
inversion (generalized dihedral). Its unit rank is \(r_1 + r_2 - 1 = h - 1\). Fix the
embedding \(H_n \subset \mathbb{C}\) in which the \(G_\mathfrak{c}\) are computed (the one restricting
to \(i \mapsto i\) on \(K\)). The embeddings of \(H_n\) into \(\mathbb{C}\) are the
\(\sigma_\mathfrak{a}\) (restricting to the identity on \(K\)) and the \(\iota\sigma_\mathfrak{a}\)
(restricting to conjugation on \(K\)), \(\mathfrak{a} \in \mathrm{Pic}(\mathcal{O}_n)\); a complex place is
a conjugate pair \(\{\tau, \iota\tau\}\), and each pair contains exactly one
\(\sigma_\mathfrak{a}\). Hence

> **Lemma 1.1 (places \(\leftrightarrow\) Pic).** \(\mathfrak{a} \mapsto v_\mathfrak{a} := \{\sigma_\mathfrak{a}, \iota\sigma_\mathfrak{a}\}\)
> is a bijection \(\mathrm{Pic}(\mathcal{O}_n) \to \{\text{archimedean places of } H_n\}\), and
> \(\log|x|_{v_\mathfrak{a}} = 2\log|\sigma_\mathfrak{a}(x)|\) in the normalization \(|x|_v = |\tau x|^2\)
> of a complex place. In particular, for \(\mathfrak{c} \ne 1\),
> $$
> \log\bigl|G_\mathfrak{c}/G_1\bigr|_{v_\mathfrak{a}} \;=\; 2\bigl(f(\mathfrak{a}^{-1}\mathfrak{c}) - f(\mathfrak{a}^{-1})\bigr).
> $$

**1.2 Regulators and indices.** For a number field \(F\) with \(r = r_1 + r_2 - 1\) and a
subgroup \(U \subset \mathcal{O}_F^\times\) with elements \(u_1, \dots, u_r\), put
$$
R(u_1,\dots,u_r) := \bigl|\det\bigl(\log|u_j|_{v}\bigr)_{v \in S,\ 1\le j\le r}\bigr| ,
$$
\(S\) any set of \(r\) of the \(r+1\) archimedean places (the choice is irrelevant: by the
product formula the rows over all \(r+1\) places sum to zero); this is the
regulator \(R_F\) when the \(u_j\) are fundamental units. The logarithmic embedding
\(\ell: \mathcal{O}_F^\times \to \mathbb{R}^{r+1}\), \(u \mapsto (\log|u|_v)_v\), has kernel \(\mu(F)\) and
image a lattice \(\Lambda_F\) of rank \(r\) in the trace-zero hyperplane; \(R(u_1,\dots,u_r)\) is
the covolume of the sublattice \(\ell(U)\) in the same normalization in which
\(R_F\) is the covolume of \(\Lambda_F\) (both minors are the covolume divided by
\(\sqrt{r+1}\)). If \(\ell(U)\) has rank \(r\), an index of full-rank lattices is the ratio
of covolumes, and \(\ell\) induces \(\mathcal{O}_F^\times/\mu(F)U \cong \Lambda_F/\ell(U)\):
$$
\bigl[\mathcal{O}_F^\times : \mu(F)\,U\bigr] \;=\; \frac{R(u_1,\dots,u_r)}{R_F}
\qquad\text{(Washington, Lemma 4.15).}
\tag{1.1}
$$
If moreover \(U\) is generated by \(u_1,\dots,u_r\), then \(U \cap \mu(F) = 1\) and \(U\) is free of
rank \(r\): a torsion element \(\prod u_j^{e_j}\) has \(\sum e_j\ell(u_j) = 0\), forcing \(e = 0\).

**1.3 The regulator of \(\mathcal{V}_n\) is a group determinant.** Take \(F = H_n\),
\(U = \mathcal{V}_n\) with generators \(v_\mathfrak{c} = G_\mathfrak{c}/G_1\) (\(\mathfrak{c} \ne 1\)), and omit the place
\(v_1\) of the identity embedding. By Lemma 1.1 the matrix is
$$
M \;=\; \Bigl(2\bigl(f(\mathfrak{a}^{-1}\mathfrak{c}) - f(\mathfrak{a}^{-1})\bigr)\Bigr)_{\mathfrak{a}, \mathfrak{c} \ne 1}, \qquad
R(\mathcal{V}_n) := |\det M| .
$$

## 2. The Dedekind determinant

> **Lemma 2.1 (Washington, *Introduction to Cyclotomic Fields*, Lemma 5.26).** Let \(G\) be
> a finite abelian group and \(g: G \to \mathbb{C}\) any function. Then
> $$
> \det\bigl(g(\sigma\tau^{-1}) - g(\sigma)\bigr)_{\sigma,\tau \in G\setminus\{1\}}
> \;=\; \prod_{\chi \ne 1}\;\sum_{\sigma \in G}\chi(\sigma)\,g(\sigma).
> $$

*Proof (for completeness).* Let \(A = (g(\sigma\tau^{-1}))_{\sigma,\tau \in G}\). For a character
\(\chi\), \(\sum_\tau g(\sigma\tau^{-1})\chi(\tau) = \chi(\sigma)\sum_\rho g(\rho)\bar\chi(\rho)\), so the
characters are eigenvectors and \(\det A = \prod_\chi S_{\bar\chi}(g) = \prod_\chi S_\chi(g)\),
\(S_\chi(g) := \sum_\sigma\chi(\sigma)g(\sigma)\). Subtract the column \(\tau = 1\) from every other
column (determinant unchanged): the new entries are \(g(\sigma\tau^{-1}) - g(\sigma)\) for
\(\tau \ne 1\) and \(g(\sigma)\) for \(\tau = 1\). Every column \(\tau \ne 1\) now sums to
\(\sum_\sigma g(\sigma\tau^{-1}) - \sum_\sigma g(\sigma) = 0\), so replacing the row \(\sigma = 1\) by the
sum of all rows (determinant unchanged) makes it \((S_1(g), 0, \dots, 0)\), and expanding
along it gives \(\det A = S_1(g)\cdot\det\bigl(g(\sigma\tau^{-1}) - g(\sigma)\bigr)_{\sigma,\tau\ne1}\).
For \(S_1(g) \ne 0\) this is the claim; both sides are polynomials in the values of
\(g\), so the identity holds in general. \(\square\)

> **Proposition 2.2 (the regulator of \(\mathcal{V}_n\)).** With \(f(\mathfrak{c}) = \log|G_\mathfrak{c}|\),
> $$
> \det M \;=\; 2^{\,h-1}\prod_{\chi \ne 1}S_\chi \;=\; 2^{\,h-1}\prod_{\chi\ne1}\bigl(-12\,L'(0,\chi)\bigr),
> \qquad
> R(\mathcal{V}_n) \;=\; 24^{\,h-1}\prod_{\chi\ne1}|L'(0,\chi)| .
> $$

*Proof.* Put \(g(x) := f(x^{-1})\). Since the group is abelian,
\(f(\mathfrak{a}^{-1}\mathfrak{c}) - f(\mathfrak{a}^{-1}) = g(\mathfrak{a}\mathfrak{c}^{-1}) - g(\mathfrak{a})\), so \(M = 2\,(g(\sigma\tau^{-1}) - g(\sigma))_{\sigma = \mathfrak{a},\ \tau = \mathfrak{c}}\)
is exactly the matrix of Lemma 2.1 for the function \(g\) (the relabelling
\(\sigma \mapsto \sigma^{-1}\) is absorbed into \(g\); no row or column permutation, hence no
sign, is involved). Therefore \(\det M = 2^{h-1}\prod_{\chi\ne1}S_\chi(g)\), and
\(S_\chi(g) = \sum_\mathfrak{a}\chi(\mathfrak{a})f(\mathfrak{a}^{-1}) = S_{\bar\chi}(f)\); as \(\chi \mapsto \bar\chi\)
permutes the nontrivial characters, \(\det M = 2^{h-1}\prod_{\chi\ne1}S_\chi(f)\). (Here in
fact \(g = f\), because \(|G_{\mathfrak{c}^{-1}}| = |\overline{G_\mathfrak{c}}| = |G_\mathfrak{c}|\).) The second
form is Paper II, Theorem 2.6(2). \(\square\)

Each \(S_\chi\) is real: \(f(\mathfrak{c}^{-1}) = f(\mathfrak{c})\) gives \(S_\chi = S_{\bar\chi} = \overline{S_\chi}\).

## 3. Nonvanishing: \(\mathcal{V}_n\) has full rank

The determinant of §2 is nonzero iff no \(L'(0,\chi)\) vanishes. This is where the
imprimitivity multipliers enter, and we take the occasion to record the Hecke
recursion of the Euler document as an identity of Dirichlet series, which gives
the multipliers at every \(s\).

**3.1 Primitive levels.** For \(m \mid n\) let \(\pi_m: \mathrm{Pic}(\mathcal{O}_n) \to \mathrm{Pic}(\mathcal{O}_m)\)
be the projection \([\Lambda] \mapsto [\mathcal{O}_m\Lambda]\); a character \(\chi\) *factors through
level \(m\)* if it is trivial on \(\ker\pi_m\), i.e. \(\chi = \chi_m\circ\pi_m\). Ideal-theoretically
(Cox, *Primes of the form \(x^2+ny^2\)*, Prop. 7.22 and Thm. 9.18),
\(\mathrm{Pic}(\mathcal{O}_m) \cong I_K(m)/P_{K,\mathbb{Z}}(m)\), where \(I_K(m)\) are the fractional ideals
prime to \(m\) and \(P_{K,\mathbb{Z}}(m)\) the principal ideals \((\alpha)\) with \(\alpha \equiv a \bmod m\mathcal{O}_K\),
\(a \in \mathbb{Z}\), \(\gcd(a,m) = 1\); and \(\pi_m\) is induced by \(I_K(n) \subset I_K(m)\). A character
\(\chi\) of \(\mathrm{Pic}(\mathcal{O}_n)\) is a Hecke character on \(I_K(n)\); its *conductor*
\(\mathfrak{f}_\chi \mid n\mathcal{O}_K\) is the smallest modulus \(\mathfrak{f}\) with \(\chi\) trivial on
\(P_{K,1}(\mathfrak{f})\cap I_K(n)\) (principal ideals \((\alpha)\), \(\alpha \equiv 1 \bmod^\times\mathfrak{f}\)), and
\(L_{\mathrm{prim}}(s,\chi) = \sum_{(\mathfrak{a},\mathfrak{f}_\chi)=1}\chi(\mathfrak{a})N\mathfrak{a}^{-s}\) is the Hecke \(L\)-function
of the primitive character.

> **Lemma 3.1 (primitive level = conductor).** Let \(p \mid n\). If \(\mathfrak{f}_\chi \mid \tfrac np\mathcal{O}_K\),
> then \(\chi\) factors through level \(n/p\). Consequently, if \(\chi\) factors through no
> level \(n/p\), every prime of \(K\) above every \(p \mid n\) divides \(\mathfrak{f}_\chi\), and
> \(L(s,\chi) = L_{\mathrm{prim}}(s,\chi)\). Moreover the set of levels \(m \mid n\) through which
> \(\chi\) factors is closed under \(\gcd\), so \(\chi\) has a well-defined **primitive level**
> \(m_\chi\), and \(L(s,\chi_{m_\chi}) = L_{\mathrm{prim}}(s,\chi)\).

*Proof.* Let \((\gamma) \in P_{K,\mathbb{Z}}(n/p)\cap I_K(n)\): \(\gamma \equiv c \bmod \tfrac np\mathcal{O}_K\) with
\(c \in \mathbb{Z}\), \(\gcd(c, n/p) = 1\), and \((\gamma)\) prime to \(n\). Changing \(c\) by a multiple of
\(n/p\) we may take \(\gcd(c,n) = 1\) (if \(p \nmid n/p\) the residues \(c + k\tfrac np\) cover all
classes mod \(p\); if \(p \mid n/p\) nothing is needed). Then \((\gamma) = (\gamma/c)\cdot(c)\) with
\((c) \in P_{K,\mathbb{Z}}(n)\) — \(\chi((c)) = 1\) — and \(\gamma/c \equiv 1 \bmod^\times \tfrac np\mathcal{O}_K\), so
\((\gamma/c) \in P_{K,1}(\mathfrak{f}_\chi)\) and \(\chi((\gamma/c)) = 1\). Hence \(\chi\) is trivial on
\(P_{K,\mathbb{Z}}(n/p)\cap I_K(n)\), which is \(\ker\pi_{n/p}\) under the isomorphism
\(I_K(n)/(P_{K,\mathbb{Z}}(n/p)\cap I_K(n)) \cong I_K(n/p)/P_{K,\mathbb{Z}}(n/p)\) (every class contains
ideals prime to \(n\)). For the consequence: if \(\chi\) factors through no \(n/p\), then for
each \(p \mid n\) some prime above \(p\) divides \(\mathfrak{f}_\chi\); ring class characters satisfy
\(\chi(\bar{\mathfrak{a}}) = \chi(\mathfrak{a})^{-1}\) (as \(\mathfrak{a}\bar{\mathfrak{a}} = (N\mathfrak{a})\) lies in
\(P_{K,\mathbb{Z}}\)), so \(\mathfrak{f}_\chi\) is conjugation-stable and both primes above a split \(p\)
divide it. Then the ideals prime to \(\mathfrak{f}_\chi\) are the ideals prime to \(n\), and
Paper II, Lemma 6.3 (the conductor-supported invertible \(\mathcal{O}_n\)-ideals contribute
nothing when \(\chi\) is nontrivial on every \(\ker\pi_{n/p}\)) gives
\(L(s,\chi) = \sum_{(\mathfrak{b},n)=1}\bar\chi(\mathfrak{b})N\mathfrak{b}^{-s} = L_{\mathrm{prim}}(s,\chi)\) (recall
\(L(s,\chi) = L(s,\bar\chi)\)). Finally, if \(\chi\) factors through \(m_1\) and \(m_2\), it is trivial
on \(P_{K,\mathbb{Z}}(m_1)\cap I_K(n)\) and on \(P_{K,\mathbb{Z}}(m_2)\cap I_K(n)\), hence on their product,
which contains \(P_{K,\mathbb{Z}}(\gcd(m_1,m_2))\cap I_K(n)\): for \((\gamma)\) with \(\gamma \equiv c \bmod \gcd(m_1,m_2)\),
\(c \in \mathbb{Z}\), write \(\gamma = (\gamma/c)\cdot c\) as above with \(\gamma/c \equiv 1 \bmod^\times \gcd(m_1, m_2)\),
and use \(P_{K,1}(\gcd(m_1,m_2)) = P_{K,1}(m_1)P_{K,1}(m_2)\) (the local statement
\((1 + \mathfrak{m}_1)(1+\mathfrak{m}_2) = 1 + \mathfrak{m}_1 + \mathfrak{m}_2\) at each prime, with weak approximation).
The last equality of the lemma is then the previous statement at the level \(m_\chi\),
where \(\chi_{m_\chi}\) factors through no proper divisor. \(\square\)

**3.2 The Hecke recursion at every \(s\).** Recall the fiber lemma
([schmidt-euler-system.md](schmidt-euler-system.md), Lemma 1.2): for a prime \(\ell\) and
\(\Lambda = \Lambda_\mathfrak{c}\), the \(\ell+1\) sublattices of index \(\ell\) are the primitive
representatives of the fiber \(\pi^{-1}(\mathfrak{c}) \subset \mathrm{Pic}(\mathcal{O}_{n\ell})\) together with the
\(\mathcal{O}_n\)-stable ones: \(\mathfrak{l}_n\Lambda = \lambda\Lambda_{[\mathfrak{l}]\mathfrak{c}}\) and
\(\bar{\mathfrak{l}}_n\Lambda\) (split \(\ell = \lambda\bar\lambda\)), none (inert), \(\mathfrak{p}_n\Lambda\)
(\(\ell = 2\)), \(\ell\Lambda^+\) with \(\Lambda^+\) the primitive index-\(n/\ell\) lattice of \(\pi(\mathfrak{c})\)
(\(\ell \mid n\)). The form of a class is \(x \mapsto N(x)\) on its lattice (discriminant \(-4n^2\)),
so \(\zeta_\mathfrak{c}(s) = \sum'_{x\in\Lambda_\mathfrak{c}}N(x)^{-s}\) and \(L(s,\chi) = \tfrac12\sum_\mathfrak{c}\chi(\mathfrak{c})\zeta_\mathfrak{c}(s)\).

> **Lemma 3.2 (Dirichlet-series form of the norm relations).** For a nontrivial \(\chi\) of
> \(\mathrm{Pic}(\mathcal{O}_n)\) and \(\chi^{(n\ell)} = \chi\circ\pi\),
> $$
> L(s,\chi^{(n\ell)}) \;=\;
> \begin{cases}
> \bigl(1 + \ell^{1-2s} - \ell^{-s}(\chi(\mathfrak{l}) + \bar\chi(\mathfrak{l}))\bigr)\,L(s,\chi) & \ell \text{ split},\ \ell\nmid n,\\
> \bigl(1 + \ell^{1-2s}\bigr)\,L(s,\chi) & \ell \text{ inert},\ \ell\nmid n,\\
> \bigl(1 + 2^{1-2s} - 2^{-s}\chi(\mathfrak{p})\bigr)\,L(s,\chi) & \ell = 2 \nmid n,\\
> \bigl(1 + \ell^{1-2s}\bigr)\,L(s,\chi) - \ell^{-2s}\,m_\pi\,L(s,\chi^{(n/\ell)}) & \ell \mid n,
> \end{cases}
> $$
> the last term present only if \(\chi\) factors through level \(n/\ell\), with
> \(m_\pi = |\ker(\mathrm{Pic}(\mathcal{O}_n)\to\mathrm{Pic}(\mathcal{O}_{n/\ell}))|\). At \(s = 0\) the multipliers are
> \(P_\ell(\chi) = \ell+1-\chi(\mathfrak{l})-\bar\chi(\mathfrak{l})\), \(\ell+1\), \(3-\chi(\mathfrak{p})\), and
> \((\ell+1, m_\pi)\): Theorem 2 of the Euler document, now at every \(s\).

*Proof.* Sum the theta series \(\theta_{\Lambda'}(\tau) = \sum_{x\in\Lambda'}q^{N(x)}\) over the
\(\ell+1\) sublattices: a point \(x \in \Lambda\) lies in \(\ell+1\) of them if \(x \in \ell\Lambda\) and in
exactly one otherwise, so \(\sum_{\Lambda'}\theta_{\Lambda'}(\tau) = \theta_\Lambda(\tau) + \ell\,\theta_\Lambda(\ell^2\tau)\),
i.e. \(\sum_{\Lambda'}\zeta_{\Lambda'}(s) = (1 + \ell^{1-2s})\zeta_\mathfrak{c}(s)\). The stable sublattices
contribute \(\theta_{\lambda\Lambda_{[\mathfrak{l}]\mathfrak{c}}}(\tau) = \theta_{[\mathfrak{l}]\mathfrak{c}}(\ell\tau)\), i.e.
\(\ell^{-s}\zeta_{[\mathfrak{l}]\mathfrak{c}}(s)\) (and the same for \(\bar{\mathfrak{l}}\), resp. \(\mathfrak{p}\)), and
\(\theta_{\ell\Lambda^+}(\tau) = \theta^{(n/\ell)}_{\pi(\mathfrak{c})}(\ell^2\tau)\), i.e. \(\ell^{-2s}\zeta^{(n/\ell)}_{\pi(\mathfrak{c})}(s)\).
Subtract them, multiply by \(\tfrac12\chi(\mathfrak{c})\) and sum over \(\mathfrak{c}\), using
\(\sum_\mathfrak{c}\chi(\mathfrak{c})\zeta_{[\mathfrak{l}]\mathfrak{c}} = \bar\chi(\mathfrak{l})\cdot 2L(s,\chi)\) (re-index by
\(\mathfrak{c}' = [\mathfrak{l}]\mathfrak{c}\)), \(\bar\chi(\bar{\mathfrak{l}}) = \chi(\mathfrak{l})\), and for \(\ell \mid n\)
\(\sum_\mathfrak{c}\chi(\mathfrak{c})\zeta^{(n/\ell)}_{\pi(\mathfrak{c})} = \sum_{\mathfrak{c}^+}\zeta^{(n/\ell)}_{\mathfrak{c}^+}\sum_{\pi(\mathfrak{c}) = \mathfrak{c}^+}\chi(\mathfrak{c})\),
which is \(0\) unless \(\chi\) is trivial on \(\ker\pi\) and then \(m_\pi\cdot 2L(s,\chi^{(n/\ell)})\).
\(\square\)

> **Corollary 3.3 (the multipliers are positive integers).** Let \(m = m_\chi\) be the
> primitive level of \(\chi\) and \(n = m\prod_\ell\ell^{k_\ell}\). Then
> \(L(s,\chi) = C_\chi(s)\,L_{\mathrm{prim}}(s,\chi)\) with \(C_\chi(s) = \prod_\ell c_{\ell,k_\ell}(s)\) a
> product of polynomials in \(\ell^{-s}\) — in particular entire — and
> \(C_\chi(0) = \prod_\ell a_{\ell,k_\ell}\), where
> $$
> a_{\ell,0} = 1,\qquad a_{\ell,1} = \begin{cases}P_\ell(\chi) & \ell \nmid m\\ \ell + 1 & \ell \mid m\end{cases},\qquad
> a_{\ell,k+1} = (\ell+1)\,a_{\ell,k} - m_{\ell,k}\,a_{\ell,k-1},\quad
> m_{\ell,1} = \begin{cases}\ell - \chi_{-4}(\ell) & \ell\nmid m\\ \ell & \ell \mid m\end{cases},\ \ m_{\ell,k} = \ell\ (k\ge2).
> $$
> Every \(a_{\ell,k}\) is a real algebraic integer \(\ge 2\) for \(k \ge 1\); hence
> \(C_\chi(0) > 0\), \(L'(0,\chi) = C_\chi(0)L'_{\mathrm{prim}}(0,\chi) \ne 0\), and
> \(\prod_{\chi\ne1}C_\chi(0)\) is a positive rational integer (it is Galois-invariant).

*Proof.* Adjoin the primes one at a time, starting from the primitive level, where
\(L(s,\chi_m) = L_{\mathrm{prim}}(s,\chi)\) (Lemma 3.1); at each step Lemma 3.2 applies with
\(\chi(\mathfrak{l})\) evaluated at the current level, which is the value at level \(m\) (the
classes of \(\mathfrak{l}\cap\mathcal{O}\) are compatible with the projections), and the last term is
present from the second \(\ell\)-step on with the stated \(m_\pi\) (\(|\ker| = N_e(n\ell)/N_e(n)\)).
Since \(L(0,\chi) = 0\) at every level (Paper II, Prop. 2.2), differentiating at \(s = 0\)
gives the recursion for the values \(L'(0,\cdot)\) with the coefficients
\(\ell + 1\) and \(m_{\ell,k}\). Positivity: \(a_{\ell,1} \ge \ell - 1 \ge 2\) in every case
(\(|\chi(\mathfrak{l})| = 1\); split \(\ell \ge 5\), inert \(\ell \ge 3\), \(3 - \chi(\mathfrak{p}) \in \{2,4\}\) since
\([\mathfrak{p}_n]^2 = [(2)] = 1\)); \(a_{\ell,2} - a_{\ell,1} = \ell a_{\ell,1} - m_{\ell,1} > 0\) in every case
(\(\ell(\ell-1) - (\ell-1)\), \(\ell(\ell+1) - (\ell+1)\), \(2\cdot2 - 2\), \(\ell(\ell+1) - \ell\)); and for \(k \ge 2\),
\(a_{\ell,k+1} - a_{\ell,k} = \ell(a_{\ell,k} - a_{\ell,k-1})\), so the sequence increases. Nonvanishing
of \(L'_{\mathrm{prim}}(0,\chi)\): by the functional equation \(L'_{\mathrm{prim}}(0,\chi) = (\sqrt{|d_{\mathfrak{f}}|}/2\pi)\,L_{\mathrm{prim}}(1,\chi)\)
up to a positive constant, and \(L_{\mathrm{prim}}(1,\chi) \ne 0\) for a nontrivial Hecke
character (Hecke). Galois invariance: an automorphism of \(\mathbb{Q}(\mu_h)\) sends \(\chi\) to
\(\chi^t\) and \(C_\chi(0)\), an integer polynomial in the values of \(\chi\), to \(C_{\chi^t}(0)\);
the product over all \(\chi \ne 1\) is a fixed positive algebraic integer, hence a positive
rational integer. \(\square\)

Together with Proposition 2.2: \(\det M \ne 0\), so \(\ell(\mathcal{V}_n)\) has rank \(h - 1\),
\(\mathcal{V}_n\) is free of rank \(h-1\) and \(\mathcal{V}_n\cap\mu(H_n) = 1\) (§1.2). The primitive
levels and multipliers of every character at \(n \le 23\) are listed in §6; e.g. at
\(n = 9\) the quadratic character comes from level \(3\) with \(C = 4\) (\(\ell = 3 \mid 3\)), at
\(n = 15\) the two pulled-back quadratic characters have \(C = 8\) (\(\ell = 5\) split,
\(\chi(\mathfrak{l}) = -1\)) and \(C = 4\) (\(\ell = 3\) inert), and \(\prod_\chi C_\chi(0) = 1, 1, 1, 4, 1, 1, 32, 1\)
at \(n = 3, 5, 7, 9, 11, 13, 15, 23\).

## 4. The class number formula, the torsion, and the proof of Theorem 1

**4.1 The class number formula at \(s = 0\).** By the Artin formalism for the abelian
extension \(H_n/K\),
$$
\zeta_{H_n}(s) \;=\; \prod_{\chi}L_{\mathrm{prim}}(s,\chi) \;=\; \zeta_K(s)\prod_{\chi\ne1}L_{\mathrm{prim}}(s,\chi),
$$
the product over the characters of \(\mathrm{Gal}(H_n/K) = \mathrm{Pic}(\mathcal{O}_n)\), each
\(L_{\mathrm{prim}}\) the Hecke \(L\)-function of the primitive character. Dedekind's formula
at \(s = 0\) (Neukirch, *Algebraic Number Theory*, Cor. VII.5.11) reads
\(\zeta_F(s) = -\tfrac{h_FR_F}{w_F}s^{\,r_F} + O(s^{r_F+1})\) with \(r_F = r_1 + r_2 - 1\), the
regulator in the normalization of §1.2. For \(K = \mathbb{Q}(i)\): \(r_K = 0\), \(\zeta_K(0) = -\tfrac{h_K}{w_K} = -\tfrac14\).
For \(H_n\): \(r = h - 1\). Each \(L_{\mathrm{prim}}(s,\chi)\), \(\chi \ne 1\), has a simple zero at
\(s = 0\) (Corollary 3.3), so comparing the coefficients of \(s^{h-1}\),
$$
-\frac{h_{H_n}R_{H_n}}{w_{H_n}} \;=\; -\frac14\prod_{\chi\ne1}L'_{\mathrm{prim}}(0,\chi),
\qquad\text{i.e.}\qquad
\prod_{\chi\ne1}L'_{\mathrm{prim}}(0,\chi) \;=\; \frac{4\,h_{H_n}R_{H_n}}{w_{H_n}} \;>\; 0 .
\tag{4.1}
$$
(The product is automatically positive; indeed every factor is: for \(\chi \ne \bar\chi\)
the two conjugate factors are equal and real, and for real \(\chi\) the genus
factorization \(L_{\mathrm{prim}} = L(s,\chi_{d_1})L(s,\chi_{d_2})\) has positive values at \(s = 1\).)

**4.2 The roots of unity of \(H_n\).**

> **Lemma 4.1.** \(\mu(H_n) = \mu_{w}\) with \(w = 4\cdot 2^{[4\mid n]}\cdot 3^{[3\mid n]}\): \(\mu_4\) always,
> \(\sqrt{-3} \in H_n\) iff \(3 \mid n\), \(\sqrt 2 \in H_n\) iff \(4 \mid n\), and nothing else.

*Proof.* \(H_n\) is generalized dihedral over \(\mathbb{Q}\) with abelianization
\(\mathrm{Pic}(\mathcal{O}_n)/\mathrm{Pic}(\mathcal{O}_n)^2\times\langle\iota\rangle\) of exponent \(2\); so if
\(\zeta_m \in H_n\) then \((\mathbb{Z}/m)^\times\) has exponent \(\le 2\), i.e. \(m \mid 24\). As \(i \in H_n\),
\(\mu(H_n) = \mu_{4\cdot a\cdot b}\) with \(a = 2\) iff \(\zeta_8 \in H_n\) iff \(\sqrt2 \in H_n\), and \(b = 3\) iff
\(\zeta_3 \in H_n\) iff \(\sqrt{-3} \in H_n\) iff \(\sqrt3 \in H_n\). Now an abelian extension \(L/K\)
lies in the ring class field \(H_n\) iff it is generalized dihedral over \(\mathbb{Q}\) and its
conductor divides \(n\mathcal{O}_K\) (Cox, Thm. 9.18; in the ideal-theoretic description
\(L \subset H_n\) iff the Artin symbols of the rational principal ideals \((a)\), \(a\) prime to
\(n\), are trivial on \(L\), and for a generalized dihedral \(L\) one has
\(((a),L/K) = \sigma_\mathfrak{p}\sigma_{\bar{\mathfrak{p}}} = \sigma_\mathfrak{p}\sigma_\mathfrak{p}^{-1} = 1\) for
\((a) = \mathfrak{p}\bar{\mathfrak{p}}\), and similarly \(= 1\) in the inert and ramified cases, since
\(\iota\) acts by inversion). \(K(\sqrt3) = K(\sqrt{-3})\) is biquadratic over \(\mathbb{Q}\), hence
generalized dihedral, unramified over \(K\) outside \(3\) (as \(\mathbb{Q}(\sqrt{-3})/\mathbb{Q}\) is
unramified at \(2\)) and ramified at \(3\): its conductor is \((3)\), which divides
\(n\mathcal{O}_K\) iff \(3 \mid n\). \(K(\sqrt2) = \mathbb{Q}(\zeta_8)\) has relative discriminant
\((1+i)^4 = (4)\) over \(K\) (\(2^8/2^2 = 2^4\) in norm), hence conductor \((4)\), which divides
\(n\mathcal{O}_K\) iff \(4 \mid n\). \(\square\)

**4.3 Proof of Theorem 1.** By §1.2–§3, \(\mathcal{V}_n\) is free of rank \(h-1\), \(\mathcal{V}_n\cap\mu(H_n) = 1\),
and (1.1) with Proposition 2.2 gives
$$
\bigl[\mathcal{O}_{H_n}^\times:\mu(H_n)\mathcal{V}_n\bigr] = \frac{R(\mathcal{V}_n)}{R_{H_n}}
= \frac{24^{h-1}\prod_{\chi\ne1}|L'(0,\chi)|}{R_{H_n}}
= \frac{24^{h-1}\prod_{\chi\ne1}C_\chi(0)\prod_{\chi\ne1}L'_{\mathrm{prim}}(0,\chi)}{R_{H_n}}
= \frac{4\cdot24^{h-1}}{w_{H_n}}\,h_{H_n}\prod_{\chi\ne1}C_\chi(0),
$$
using Corollary 3.3 (\(|L'(0,\chi)| = C_\chi(0)|L'_{\mathrm{prim}}(0,\chi)|\) with \(C_\chi(0) > 0\), and
\(L'_{\mathrm{prim}}(0,\chi) > 0\)) and (4.1). The torsion is Lemma 4.1, and
\([\mathcal{O}_{H_n}^\times : \mathcal{V}_n] = [\mathcal{O}_{H_n}^\times:\mu\mathcal{V}_n]\cdot[\mu\mathcal{V}_n:\mathcal{V}_n] = w_{H_n}\cdot[\mathcal{O}_{H_n}^\times:\mu\mathcal{V}_n]\)
since \(\mathcal{V}_n\cap\mu = 1\). \(\blacksquare\)

**Anchors, by hand.** At \(n = 3\): \(H_3 = \mathbb{Q}(\zeta_{12})\), \(w = 12\), \(h_{H_3} = 1\), and
\(D_3 = x^2 - 378x + 729\) has roots \(189 \pm 108\sqrt3 = 27\,\varepsilon_{12}^{\pm2}\), so
\(\mathcal{V}_3 = \langle\varepsilon_{12}^4\rangle\); the unit index of \(\mathbb{Q}(\zeta_{12})\) over
\(\mathbb{Q}(\sqrt3)\) is \(2\) (Washington, Cor. 4.13: \(12\) is not a prime power), so
\([\mathcal{O}^\times:\mu_{12}\langle\varepsilon_{12}^4\rangle] = 2\cdot4 = 8 = 4\cdot24/12\). At \(n = 5\):
\(H_5 = \mathbb{Q}(i,\sqrt5)\), \(w = 4\), \(D_5 = x^2 - 322x + 1\) with roots \(161 \pm 72\sqrt5 = \varepsilon_5^{\pm12}\),
\(\mathcal{V}_5 = \langle\varepsilon_5^{24}\rangle\), index \(24\,Q\) with \(Q\) the unit index of the
biquadratic field, and \(h_{H_5} = \tfrac12Q\,h(-4)h(5)h(-20) = Q\) (the class number formula
for imaginary bicyclic biquadratic fields); PARI: \(h_{H_5} = 1 = Q\), index \(24\).

## 5. Layers: the index in every subfield cut out by a subgroup

Let \(A \le \mathrm{Pic}(\mathcal{O}_n)\) be a subgroup, \(B = \mathrm{Pic}(\mathcal{O}_n)/A\), and \(H^A \subset H_n\)
the fixed field of \(A\), so \(\mathrm{Gal}(H^A/K) = B\), \([H^A:K] = |B|\), and \(H^A\) is again
totally complex generalized dihedral, of unit rank \(|B| - 1\). For \(b \in B\) put
$$
\theta^{(b)} \;:=\; \prod_{\mathfrak{c}\in bA}G_\mathfrak{c} \;=\; N_{H_n/H^A}(G_{\mathfrak{c}_0})\quad(\mathfrak{c}_0 \in bA),
\qquad
\mathcal{V}^A := \bigl\langle\theta^{(b)}/\theta^{(b')}\bigr\rangle = N_{H_n/H^A}(\mathcal{V}_n) \subset \mathcal{O}_{H^A}^\times ,
$$
the second equality because \(\sigma_\mathfrak{a}(G_{\mathfrak{c}_0}) = G_{\mathfrak{a}^{-1}\mathfrak{c}_0}\) runs over the
coset as \(\mathfrak{a}\) runs over \(A = \mathrm{Gal}(H_n/H^A)\).

> **Theorem 2 (layer indices).** \(\mathcal{V}^A\) is free of rank \(|B| - 1\), meets \(\mu(H^A)\)
> trivially, and
> $$
> \bigl[\mathcal{O}_{H^A}^\times : \mu(H^A)\,\mathcal{V}^A\bigr]
> \;=\; \frac{24^{\,|B|-1}\prod_{\psi\in\widehat B\setminus1}|L'(0,\psi)|}{R_{H^A}}
> \;=\; \frac{4\cdot 24^{\,|B|-1}}{w_{H^A}}\;h_{H^A}\prod_{\psi\in\widehat B\setminus1}C_\psi(0),
> $$
> the products over the characters of \(\mathrm{Pic}(\mathcal{O}_n)\) trivial on \(A\).

*Proof.* Word for word the proof of Theorem 1 with \(\mathrm{Pic}(\mathcal{O}_n)\) replaced by \(B\)
and \(f\) by \(f_B(b) := \sum_{\mathfrak{c}\in bA}f(\mathfrak{c}) = \log|\theta^{(b)}|\): the places of \(H^A\)
correspond to \(B\) (Lemma 1.1 for \(H^A\)); \(\log|\sigma_b(\theta^{(b')}/\theta^{(1)})|\) is
\(f_B(b^{-1}b') - f_B(b^{-1})\); the character sums are
\(\sum_b\psi(b)f_B(b) = \sum_\mathfrak{c}\psi(\mathfrak{c})f(\mathfrak{c}) = S_\psi = -12L'(0,\psi)\) for \(\psi \in \widehat B\);
Lemma 2.1 and Corollary 3.3 give \(R(\mathcal{V}^A) = 24^{|B|-1}\prod_{\psi\ne1}|L'(0,\psi)| \ne 0\);
and \(\zeta_{H^A} = \zeta_K\prod_{\psi\in\widehat B\setminus1}L_{\mathrm{prim}}(s,\psi)\) gives (4.1) for
\(H^A\). \(\square\)

**5.1 The mass-normalized coset units.** When \(|M(n)|^{1/|B|} \in \mathbb{Z}\) (always the case for
\(|B| \in \{2, 3\}\) at odd \(n\) with \(3 \mid h\) resp. \(2 \mid h\): every exponent of the mass law
carries the factor \(6\), and \(N_e\) of the cofactor is even), the numbers
\(\theta^{(b)}_u := \theta^{(b)}/|M(n)|^{1/|B|}\) are units of \(H^A\) by the per-class valuation
law (Paper II, Prop. 6.2 for \(|B| = 3\); the same argument for any \(|B|\)), with
\(\prod_b\theta^{(b)}_u = \pm1\), and \(\mu\langle\theta^{(b)}_u : b\rangle \supset \mu\mathcal{V}^A\) with index
exactly \(|B|\) (the log-lattice of the \(\theta_u^{(b)}\) is \(\mathbb{Z}^B/\mathbb{Z}(1,\dots,1)\) and the ratio
lattice is the image of the sum-zero vectors, of index \(|B|\)). Hence
$$
\bigl[\mathcal{O}_{H^A}^\times : \mu(H^A)\langle\theta^{(b)}_u\rangle_b\bigr] \;=\; \frac{4\cdot24^{|B|-1}}{|B|\,w_{H^A}}\,h_{H^A}\prod_{\psi\ne1}C_\psi(0).
$$

**5.2 The quadratic layer** (\(A = \ker\chi_2\), \(B = \mathbb{Z}/2\), \(H^A = K(\sqrt{d_2}) = \mathbb{Q}(i,\sqrt{d_2})\)
for the genus factorization \(-4n^2 = d_1d_2\) of Paper II, Prop. 6.1, \(d_2 > 0\)). Theorem 2
gives \([\mathcal{O}_{K(\sqrt{d_2})}^\times:\mu\mathcal{V}^A] = (96/w)\,h(K(\sqrt{d_2}))\,C_{\chi_2}(0)\). The *real*
statement: \(\theta^{(2)}_u := \theta^{(1)}/|M(n)|^{1/2}\) is fixed by complex conjugation
(\(\overline{G_\mathfrak{c}} = G_{\mathfrak{c}^{-1}}\) and \(\mathfrak{c}^{-1} \in A\) iff \(\mathfrak{c} \in A\)), so it is a unit of the
real quadratic field \(\mathbb{Q}(\sqrt{d_2})\), and
\(\log|\theta^{(2)}_u| = \tfrac12S_{\chi_2} = -6L'(0,\chi_2) = -6\,m_{\chi}\log\varepsilon_{d_2}\) with
\(m_\chi = \tfrac{2h(d_1^*)}{w(d_1^*)}h(d_2)C(0)\) (Paper II, Prop. 6.1):

> **Corollary 5.1 (real quadratic layer).** \(\theta^{(2)}_u = \pm\varepsilon_{d_2}^{-6m_\chi}\), so
> \(\bigl[\mathbb{Z}[\varepsilon_{d_2}]^\times : \langle-1,\theta^{(2)}_u\rangle\bigr] = 6\,m_\chi\)
> (\(= 2, 12, 6, 8, 6, 12\) at \(n = 3, 5, 7, 9, 11, 13\); at \(n = 15\): \(48, 24, 16\) for the three real
> characters — the one from level \(5\) (\(\mathbb{Q}(\sqrt5)\)), the primitive genus character
> (\(\mathbb{Q}(\sqrt{15})\)), and the one from level \(3\) (\(\mathbb{Q}(\sqrt3)\)).

(Here \(\varepsilon_{d_2}\) is the fundamental unit of \(\mathbb{Q}(\sqrt{d_2})\); the script checks
\(\theta^{(2)}_u\) against PARI's fundamental unit, exactly.)

**5.3 The cubic and sextic layers** (\(3 \mid h\); \(A = \ker\chi_3\), \(B = \mathbb{Z}/3\)). Paper II,
§6.2, works in the real cubic field \(L_3 = F\cap\mathbb{R}\) of the sextic \(F := H^A = K(\theta_u)\),
with \(\theta_u = \theta^{(1)}_u\): \([\mathcal{O}_{L_3}^\times:\langle-1,\theta_u\rangle] = 8h_{L_3}C_n(0)\).
Theorem 2 gives the sextic statement. Note that \(F\) is **not** a CM field: \(L_3\) has
signature \((1,1)\), so \(F\) is totally complex of unit rank \(2\) while \(\mathcal{O}_{L_3}^\times\)
has rank \(1\) — there is no finite "unit index of \(F\) over \(L_3\)"; the finite-index
statement is about the ratio group \(\mathcal{V}^A = \langle\theta_u^{(b)}/\theta_u^{(b')}\rangle\) of
rank \(2\):

> **Corollary 5.2 (sextic layer).** With \(w_F = 4\) (Lemma 4.1 applied to \(F \subset H_n\):
> \(\sqrt3, \sqrt2 \notin F\) since \([F:K] = 3\) is odd),
> $$
> \bigl[\mathcal{O}_F^\times:\mu_4\mathcal{V}^A\bigr] = 576\,h_F\,C_{\chi_3}(0)^2,\qquad
> \bigl[\mathcal{O}_F^\times:\mu_4\langle\theta_u,\theta_u'\rangle\bigr] = 192\,h_F\,C_{\chi_3}(0)^2 .
> $$

At \(n = 13\): \(h_F = 3\) (\(= h_{L_3}\)), exponent vectors of \(\theta_u, \theta_u', \theta_u''\) in
PARI's fundamental-unit basis of \(F\): \((24,0)\), \((0,24)\), \((-24,-24)\), index \(576\) for
the conjugates, \(1728 = 576\cdot3\) for the ratio group (§6). The factor
\(576 = 24^2\) is the \(|B| = 3\) instance of \(24^{|B|-1}\), exactly as \(8 = 24/3\) in the real
cubic layer is the weight bookkeeping \(12/(3/2)\) of Paper II.

**5.4 The full index versus the layers.** The index of Theorem 1 does not factor over
the layers (the lattice \(\ell(\mathcal{V}_n)\) is not the sum of its eigen-pieces), but the
*regulator* does: \(R(\mathcal{V}_n) = \prod_{\chi\ne1}|S_\chi|\cdot 2^{h-1}\) is the product over the
character orbits, and each orbit's contribution is the regulator of the
corresponding layer up to the powers of \(2\) and \(3\). The class number \(h_{H_n}\)
absorbs the class numbers of all the subfields through the Brauer–Kuroda relations;
e.g. at \(n = 13\), \(h_{H_{13}} = 3 = h_F = h_{L_3}\), and at \(n = 23\), \(h_{H_{23}} = 12\) with
\(\mathrm{Cl}(H_{23}) \cong \mathbb{Z}/6\times\mathbb{Z}/2\) (§6).

**5.5 The \(24\)-th-root saturation: index exactly \(h_{H_n}\prod C_\chi(0)\).** The exact
exponent matrices \(E_n\) of the units \(G_\mathfrak{c}/G_1\) in PARI's fundamental-unit bases
(§6) show more than the index. Let
$$
\mathcal{W}_n \;:=\; \bigl\{u \in \mathcal{O}_{H_n}^\times : u^{24} \in \mu(H_n)\mathcal{V}_n\bigr\}
$$
be the \(24\)-th-root saturation of \(\mu\mathcal{V}_n\) (the natural object, since
\(\Delta = \eta^{24}\) and the ratios \(\Delta(\Lambda)/\Delta(\Lambda')\) are \(24\)-th powers of
\(\eta\)-quotients). In terms of the Smith invariants \(d_1, \dots, d_{h-1}\) of \(E_n\),
\([\mathcal{W}_n:\mu\mathcal{V}_n] = \prod_i\gcd(d_i, 24)\) and
\(\mathcal{O}_{H_n}^\times/\mathcal{W}_n \cong \bigoplus_i\mathbb{Z}/(d_i/\gcd(d_i,24))\). At every computed level:

| \(n\) | \(w\) | \(E_n \subset\) | Smith invariants of \(E_n\) | \([\mathcal{W}_n:\mu\mathcal{V}_n]\) | \(\mathcal{O}^\times_{H_n}/\mathcal{W}_n\) | \(h_{H_n}\prod C_\chi(0)\) |
|---|---|---|---|---|---|---|
| 3 | 12 | \(8\mathbb{Z}\) | \((8)\) | \(8\) | \(1\) | \(1\) |
| 5 | 4 | \(24\mathbb{Z}\) | \((24)\) | \(24\) | \(1\) | \(1\) |
| 7 | 4 | \(24\mathbb{Z}\) | \((24,24,24)\) | \(24^3\) | \(1\) | \(1\) |
| 9 | 12 | \(8\mathbb{Z}\) | \((96,24,24,24,8)\) | \(24^5/3\) | \(\mathbb{Z}/4\) | \(4\) |
| 11 | 4 | \(24\mathbb{Z}\) | \((24,24,24,24,24)\) | \(24^5\) | \(1\) | \(1\) |
| 13 | 4 | \(24\mathbb{Z}\) | \((72,24,24,24,24)\) | \(24^5\) | \(\mathbb{Z}/3\) | \(3\) |
| 15 | 12 | \(8\mathbb{Z}\) | \((384,48,48,24,24,24,8)\) | \(24^7/3\) | \(\mathbb{Z}/16\times\mathbb{Z}/2\times\mathbb{Z}/2\) | \(64\) |
| 23 | 4 | \(24\mathbb{Z}\) | \((144,48,24,\dots,24)\) | \(24^{11}\) | \(\mathbb{Z}/6\times\mathbb{Z}/2\) | \(12\) |

So at the levels with \(w = 4\) every \(G_\mathfrak{c}/G_1\) is a \(24\)-th power in
\(\mathcal{O}_{H_n}^\times\) times a root of unity (at \(n = 13, 23\) even an exact \(24\)-th power:
all torsion components vanish), at the levels with \(w = 12\) an \(8\)-th power, with
exactly one Smith invariant not divisible by \(3\); and in every case
$$
\bigl[\mathcal{O}_{H_n}^\times : \mathcal{W}_n\bigr] \;=\; h_{H_n}\prod_{\chi\ne1}C_\chi(0)
\qquad\text{(certified at } n = 3, 5, 7, 9, 11, 13, 15, 23\text{)},
$$
with \(\mathcal{O}^\times_{H_n}/\mathcal{W}_n \cong \mathrm{Cl}(H_n)\) as abstract groups at the levels where
all characters are primitive and \(h_{H_n} > 1\) (\(\mathbb{Z}/3\) at \(n = 13\), \(\mathbb{Z}/6\times\mathbb{Z}/2\) at
\(n = 23\), matching PARI's \(\mathrm{Cl}(H_{23}) \cong [6,2]\)). Given Theorem 1, the displayed
identity is equivalent to \([\mathcal{W}_n:\mu\mathcal{V}_n] = 4\cdot24^{h-1}/w_{H_n}\), i.e. to a
statement about which \(24\)-th roots of the \(\Delta\)-quotients lie in \(H_n\) — the
Siegel-unit refinement of the \(\Delta\)-data — which we state as:

> **Conjecture 5.3 (saturated Kubert–Lang form).** For every \(n \ge 2\),
> \([\mathcal{O}_{H_n}^\times:\mathcal{W}_n] = h_{H_n}\prod_{\chi\ne1}C_\chi(0)\); at levels all of whose
> characters are primitive, \(\mathcal{O}_{H_n}^\times/\mathcal{W}_n \cong \mathrm{Cl}(H_n)\) as
> \(\mathbb{Z}[\mathrm{Gal}(H_n/K)]\)-modules.

This is the Gras-type statement for the Schmidt \(\Delta\)-units; the first half should
follow from Theorem 1 and the theory of Siegel units in ring class fields
(Kubert–Lang, Ch. 11–12; Ramachandra), the second is of the Gras/Rubin type and is
open. Note that Theorem 1 already implies the *order* of the \(p\)-part of
\(\mathcal{O}^\times/\mathcal{W}_n\) for every \(p \ge 5\): it equals that of \(h_{H_n}\prod C_\chi(0)\),
since the saturation only changes \(2\)- and \(3\)-parts.

**5.6 The hyperbolic side: the odd-unit index.** The hyperbolic units of Paper II,
Theorem 4.2, \(R_f = r_0^6\Delta(\mathfrak{b}_1)/\Delta(\mathfrak{r}^{-1}\mathfrak{b}_1)\), \(f \in \mathrm{Cl}(D)\),
\(D = 1 - n^2\), live in the ring class field \(H\) of the order of discriminant \(D\) over
\(K' = \mathbb{Q}(\sqrt D)\), with \(\sigma_\mathfrak{a}(R_f) = R_{\mathfrak{a}^{-1}f}\), \(\overline{R_f} = R_{f^{-1}}\)
and \(R_{\mathfrak{r}f} = R_f^{-1}\). Let \(\tau = \sigma_\mathfrak{r}\) (an involution of \(H/K'\)),
\(H^+ = H^\tau\) (degree \(h\) over \(\mathbb{Q}\), totally complex — it is not a real subfield),
\(E = \mathcal{O}_H^\times/\mu\), \(E^\pm = \ker(1 \mp \tau) \subset E\) (the "even" and "odd" units;
\(E^+\) is the image of \(\mathcal{O}_{H^+}^\times\)), and \(Q^- := [E : E^+E^-]\), a power of \(2\). The
character sums vanish on the even characters and are \(-24L'(0,\chi)\) on the odd ones
(\(\chi(\mathfrak{r}) = -1\); Paper II, Theorem 2.10), so \(\langle R_f\rangle \subset E^-\) has rank
\(h/2 = \operatorname{rank}E^-\).

> **Theorem 3 (hyperbolic odd index).** For every odd \(n \ge 3\),
> $$
> \bigl[E^- : \langle R_f : f \in \mathrm{Cl}(D)\rangle\bigr]
> \;=\; 24^{h/2}\cdot\frac{2^{h/2-1}}{Q^-}\cdot\frac{h_H\,w_{H^+}}{h_{H^+}\,w_H}\cdot\prod_{\chi\ \mathrm{odd}}C_\chi(0),
> $$
> where \(C_\chi(0) = L'(0,\chi)/L'_{\mathrm{prim}}(0,\chi)\) (\(=1\) whenever \(D\) is a fundamental
> discriminant, e.g. \(n = 21\)).

*Proof.* Write \(f(x) = \log|R_x|\), an odd function on \(G = \mathrm{Cl}(D)\) with respect to
\(C = \langle\mathfrak{r}\rangle\), and let \(T\) be a transversal of \(C\). The log-vectors of \(R_f\)
at the places \(v_\mathfrak{a}\) are odd under \(\mathfrak{a} \mapsto \mathfrak{r}\mathfrak{a}\); the odd
subspace of \(\mathbb{R}^G\) has dimension \(h/2\) and the coordinates indexed by \(T\) are
coordinates on it. For \(R^-(U) := |\det(2\log|\sigma_a(u_j)|)_{a\in T,\,j}|\) (a set of \(h/2\)
elements \(u_j\) of \(E^-\)) the index of full-rank sublattices of \(E^-\) is again the ratio
of the \(R^-\)'s. *(i) The twisted group determinant.* For \(a \in T\) and an odd
character \(\chi\), \(\sum_{b\in T}f(a^{-1}b)\chi(b) = \tfrac12\sum_{b\in G}f(a^{-1}b)\chi(b) = \tfrac12\chi(a)S_{\bar\chi}(f)\)
(the two elements of a coset contribute equally), and the odd characters restricted to
\(T\) are a basis of the functions on \(T\); hence \(\det(f(a^{-1}b))_{a,b\in T} = 2^{-h/2}\prod_{\chi\ \mathrm{odd}}S_\chi(f)\)
and \(R^-(\langle R_f\rangle) = \bigl|\prod_{\chi\ \mathrm{odd}}S_\chi\bigr| = 24^{h/2}\prod_{\chi\ \mathrm{odd}}|L'(0,\chi)| \ne 0\).
*(ii) The relative class number formula.* \(\zeta_H/\zeta_{H^+} = \prod_{\chi\ \mathrm{odd}}L_{\mathrm{prim}}(s,\chi)\)
(the characters of \(\mathrm{Gal}(H^+/K')\) are the even ones), and \(r_H - r_{H^+} = h/2\), so
\(\prod_{\chi\ \mathrm{odd}}L'_{\mathrm{prim}}(0,\chi) = \dfrac{h_HR_H/w_H}{h_{H^+}R_{H^+}/w_{H^+}}\).
*(iii) Regulators.* In the orthogonal coordinates \((x_{v}+x_{\tau v})/\sqrt2\), \((x_v - x_{\tau v})/\sqrt2\)
on \(\mathbb{R}^{\mathrm{places}(H)}\) the lattice \(E^+E^-\) splits, and the covolume computation gives
\(R_H/R_{H^+} = 2^{h/2-1}R^-(E^-)/Q^-\) (the factor \(2^{h/2-1}\) is the CM-field factor
\(2^{r}/Q\) of Washington, Prop. 4.16, in this non-CM situation; both \(H\) and \(H^+\)
are totally complex, so no \(\delta_v\)-mismatch occurs). *(iv)* Assemble:
\([E^-:\langle R_f\rangle] = R^-(\langle R_f\rangle)/R^-(E^-)\), insert (i)–(iii) and
\(L'(0,\chi) = C_\chi(0)L'_{\mathrm{prim}}(0,\chi)\). \(\square\)

The same statement holds in every layer \(F = H^A\), \(A \le \mathrm{Cl}(D)\) with \(\mathfrak{r} \notin A\),
for the coset units \(u_b = N_{H/F}(R_f)\) and the odd characters of \(\mathrm{Cl}(D)/A\).

**The record** (script phase H; \(\mathrm{Cl} = \mathrm{Cl}(D)\), \(\tau\)-action on \(\mu(H)\) shown as
\(\zeta \mapsto \zeta^t\); all fields `bnfcertify` = 1 except \(H\) at \(n = 21\), degree 24, GRH):

| \(n\) | \(D\) | \(\mathrm{Cl}\) | \(h_H\), \(\mathrm{Cl}(H)\) | \(w_H\) | \(h_{H^+}\), \(\mathrm{Cl}(H^+)\) | \(w_{H^+}\) | \(t\) | \(Q^-\) | \(\prod_{\mathrm{odd}}C\) | \([E^-:\langle R_f\rangle]\) |
|---|---|---|---|---|---|---|---|---|---|---|
| 5 | \(-24\) | \(\mathbb{Z}/2\) | 1 | 6 | 2, \(\mathbb{Z}/2\) | 2 | 5 | 1 | 1 | \(4 = 24/6\) |
| 7 | \(-48 = 4^2(-3)\) | \(\mathbb{Z}/2\) | 1 | 12 | 1 | 6 | 7 | 1 | 1 | \(12 = 24/2\) |
| 9 | \(-80 = 2^2(-20)\) | \(\mathbb{Z}/4\) | 1 | 4 | 1 | 4 | 1 | 2 | 1 | \(576 = 24^2\) |
| 11 | \(-120\) | \((\mathbb{Z}/2)^2\) | 2, \(\mathbb{Z}/2\) | 6 | 4, \(\mathbb{Z}/4\) | 2 | 5 | 2 | 1 | \(96 = 24^2/6\) |
| 13 | \(-168\) | \((\mathbb{Z}/2)^2\) | 1 | 6 | 2, \(\mathbb{Z}/2\) | 6 | 1 | 2 | 1 | \(288 = 24^2/2\) |
| 15 | \(-224 = 2^2(-56)\) | \(\mathbb{Z}/4\times\mathbb{Z}/2\) | 1 | 8 | 2, \(\mathbb{Z}/2\) | 2 | 7 | 8 | 16 | \(663552 = 2\cdot24^4\) |
| 17 | \(-288 = 6^2(-8)\) | \((\mathbb{Z}/2)^2\) | 1 | 24 | 1 | 8 | 17 | 2 | 4 | \(768 = \tfrac43\,24^2\) |
| 19 | \(-360 = 3^2(-40)\) | \(\mathbb{Z}/4\times\mathbb{Z}/2\) | 1 | 6 | 1 | 6 | 1 | 4 | 1 | \(663552 = 2\cdot24^4\) |
| 21 | \(-440\) | \(\mathbb{Z}/6\times\mathbb{Z}/2\) | 4, \(\mathbb{Z}/4\) | 2 | 4, \((\mathbb{Z}/2)^2\) | 2 | 1 | 32 | 1 | \(191102976 = 24^6\) |

Every row satisfies Theorem 3 exactly (the exact index from `bnfisunit` against the
right side), the regulator relation \(R_H/R_{H^+} = 2^{h/2-1}R^-(E^-)/Q^-\) (spare \(\ge 119\)),
and the relative class number formula \(\prod_{\mathrm{odd}}L'_{\mathrm{prim}}(0,\chi) = (hR/w)_H/(hR/w)_{H^+}\)
(spare \(\ge 118\)). Three readings. (i) The powers of \(2\) and \(3\) are governed by the
roots of unity of \(H\) and \(H^+\) and by \(Q^-\); \(Q^- = 2^{h/2-1}\) at eight of the nine
levels (so that the \(2\)-adic factor of Theorem 3 is \(1\)) but \(Q^- = 2^{h/2-2}\) at \(n = 19\):
\(Q^-\) is a genuine invariant, not a bookkeeping constant. (ii) At the fundamental
levels with \(w_H = w_{H^+}\) — \(n = 9, 19, 21\) — the index is \(24^{h/2}\) times
\(h_H/h_{H^+}\) times \(2^{h/2-1}/Q^-\): \(24^2\), \(2\cdot24^4\), \(24^6\). (iii) The multipliers of
the imprimitive levels are computed through exact projections \(\mathrm{Cl}(D) \to \mathrm{Cl}(f'^2d_K)\)
(ideal extension in the order of conductor \(f'\)) and the independent evaluation at
the smaller discriminant: at \(n = 15\) the two odd quartic characters come from
\(\mathrm{Cl}(-56)\) with \(C = 4 = 3 - \chi(\mathfrak{p}_2)\) each (\(2\) ramified in \(\mathbb{Q}(\sqrt{-14})\),
\(\chi(\mathfrak{p}_2) = -1\)), at \(n = 17\) one odd character comes from \(\mathrm{Cl}(-72)\) with
\(C = 4\) — the hyperbolic instances of Lemma 3.2, whose proof is the same for any
imaginary quadratic field.

**The sextic layer at \(n = 21\).** For the odd character \(\chi\) of order \(6\) with kernel
\(A\) (order \(2\)), the coset units \(u_b = N_{H/F}(R_f)\), \(F = H^A\) of degree \(12\), are the
roots of the certified palindromic sextic
\(x^6 - 87350782811055827117922566\,x^5 - 114499870290142503112954895345\,x^4 - 3230730316507362342294534850844180\,x^3 - \cdots + 1\);
\(F\): \(h_F = 8\) (\(\mathrm{Cl}(F) \cong \mathbb{Z}/8\)), \(w_F = 2\); \(F^+\) (degree 6): \(h = 4\), \(w = 2\);
\(Q^-_F = 4 = 2^{3-1}\); both `bnfcertify` = 1. The exact index is
$$
\bigl[E^-(F) : \langle u_b\rangle\bigr] \;=\; 27648 \;=\; 2\cdot24^{3} \;=\; 24^3\cdot\frac{2^2}{Q^-_F}\cdot\frac{h_F}{h_{F^+}},
$$
the first genuinely hyperbolic Robert-index datum, **unconditional**: \(24^3\) times the
relative class number \(h_F/h_{F^+} = 2\).

## 6. Machine verification

`python3 scripts/robert_index_full.py --selftest` (Euclidean levels \(n = 3, 5, 7, 9, 11, 13, 15\)
and hyperbolic levels \(n = 5, 7, \dots, 21\); add `--with-23` for the degree-\(24\) Euclidean
level; 78 s in total on one core, plus 55 s for \(n = 23\), of which 40 s are PARI's
construction of \(H_{23}\)). Environment: Python 3 with `mpmath`
(150 digits) and `sympy` (Smith normal forms), and **PARI/GP 2.15.4** (`gp`, installed
with `apt-get install pari-gp`), called through `gp -q` with generated scripts. For
each level the script asserts:

1. **(V1)** the determinant identity of Proposition 2.2 from the \(h\) numbers
   \(f(\mathfrak{c}) = \log|G_\mathfrak{c}|\): \(\det\bigl(2(f(\mathfrak{a}^{-1}\mathfrak{c}) - f(\mathfrak{a}^{-1}))\bigr)_{\mathfrak{a},\mathfrak{c}\ne1} = 2^{h-1}\prod_{\chi\ne1}S_\chi\)
   *with sign* (relative spare \(\ge 150\) digits); this checks Lemma 2.1 and the
   normalization of §1 independently of PARI. The class groups, products, inverses,
   kernels and characters are the exact HNF arithmetic of `schmidt_euler_system.py`.
2. **(V2)** \(S_\chi = -12L'(0,\chi)\) for every nontrivial \(\chi\), with \(L'(0,\chi)\) from the
   independent incomplete-gamma evaluation of Paper II, Lemma 2.3 (exact
   representation numbers of the lattices, \(218\)–\(4000\) terms; relative spare \(\ge 148\)).
3. **(V3)** the primitive level \(m_\chi\) of every character (exactly, from the kernels of
   the projections; the set of levels through which \(\chi\) factors is asserted to be
   \(\gcd\)-closed), the multiplier \(C_\chi(0)\) from the recursion of Corollary 3.3, and,
   for every imprimitive \(\chi\), the identity \(L'(0,\chi) = C_\chi(0)\,L'(0,\chi_{m_\chi})\) with
   both sides from the independent evaluation (spare \(\ge 148\)); \(\prod_\chi C_\chi(0)\)
   certified as an integer.
4. **(V4)** PARI: \(H_n = \) `polredbest(polcompositum(y^2+1, polclass(-4n^2)))`,
   `bnfinit(·, 1)` at 120 digits: \(h_{H_n}\), \(\mathrm{Cl}(H_n)\), \(w_{H_n}\) (= `bnf.tu[1]`, checked
   against Lemma 4.1 through `nfroots` for \(\sqrt2, \sqrt3\)), \(R_{H_n}\); `bnfcertify` at
   degree \(\le 16\) (returned \(1\) at every such level — the unit groups and class
   numbers at \(n \le 15\) are **unconditional**; \(n = 23\) is GRH-conditional). Then
   \(24^{h-1}\prod|L'(0,\chi)|/R_{H_n}\) is certified as an integer with the absolute-error
   criterion and \(\ge 105\) spare digits (\(\ge 40\) demanded), and asserted equal to
   \((4\cdot24^{h-1}/w_{H_n})\,h_{H_n}\prod_\chi C_\chi(0)\).
5. **(V5)** the exact index: the certified integer polynomial \(D_n(x) = \prod_\mathfrak{c}(x - G_\mathfrak{c})\)
   (integer coefficients with \(\ge 82\) spare digits; \(D_n(0) = (-1)^hM(n)\) checked
   against the mass law; the polynomials of Paper II, Thm 2.6 reproduced) is factored in
   \(H_n\) by `nfroots` (all \(h\) roots found), the units \(G_\mathfrak{c}/G_1\) are written in
   PARI's fundamental-unit basis by `bnfisunit`, and \(|\det|\) of the integer exponent
   matrix \(E_n\) is asserted equal to the index of (V4). The Smith normal form of \(E_n\),
   the torsion components, and the saturation quotient of §5.5 are recorded.
6. **(V6)** the layers: for every real character, \(\theta^{(2)}_u = \prod_{\ker}G/|M|^{1/2}\) is
   certified as a root of an integer polynomial \(x^2 - tx \pm 1\), the field
   \(\mathbb{Q}(\theta^{(2)}_u)\) is built in PARI, and `bnfisunit` gives the exact exponent
   \(\pm6m_\chi\) against PARI's (certified) fundamental unit; for the cubic character
   (\(n = 9, 11, 13, 23\)) the coset cubic and the \(\theta_u\)-cubic are certified,
   \(L_3 = \mathbb{Q}(\theta_u)\) and \(F = K(\theta_u)\) are built, \(h_{L_3}, R_{L_3}, h_F, R_F\)
   computed (`bnfcertify` = 1 for both), \(\theta_u = \pm\eta^{\pm8h_{L_3}C_n(0)}\) read off exactly
   (an independent confirmation of Paper II, Theorem 6.6 — the fundamental units there
   were obtained by root descent, here by PARI), and the exponent vectors of the three
   conjugates of \(\theta_u\) in \(F\) give the sextic indices of Corollary 5.2 exactly.

7. **(H1–H4)** the hyperbolic phase, odd \(n = 5, \dots, 21\): the units \(R_f\) from the
   lattice lemma (Paper II, Lemma 4.1), the \(R\)-polynomial re-certified against the
   record of `schmidt_units.py`; the twisted group determinant over a transversal of
   \(\langle\mathfrak{r}\rangle\) against \(\prod_{\mathrm{odd}}S_\chi\) (signed, spare \(\ge 150\)) with
   \(S_\chi = -24L'(0,\chi)\) (independent evaluation) and the vanishing of the even sums;
   the primitive conductor of every odd character through exact projections to the
   orders of conductor \(f' \mid f\) and \(C_\chi(0) = L'(0,\chi)/L'(0,\chi_{f'})\) (both
   independent evaluations), \(\prod C\) certified as an integer; PARI: \(H\) (from
   `polclass`), \(h_H\), \(w_H\), \(R_H\), the \(h\) roots \(R_f\) as units, \(\tau\) as the unique
   automorphism of \(H/K'\) inverting every \(R_f\), \(E^\pm\) as kernels on the unit lattice,
   the exact index \([E^-:\langle R_f\rangle]\), \(H^+\) generated by \(\sqrt D + \sum c^i(R_{f_i} + R_{f_i}^{-1})\),
   \(h_{H^+}, w_{H^+}, R_{H^+}\), \(Q^-\) (with the action of \(\tau\) on \(\mu(H)\)), the regulator
   relation of step (iii), the relative class number formula of step (ii), and the
   identity of Theorem 3 as an equality of integers; at \(n = 21\) also the sextic layer.
   `bnfcertify` succeeded for every field of degree \(\le 16\).

The record (\(\mathrm{Pic} = \mathrm{Pic}(\mathcal{O}_n)\); "index" \(= [\mathcal{O}_{H_n}^\times:\mu(H_n)\mathcal{V}_n]\);
"spare" = spare digits of the integer certification of \(24^{h-1}\prod|L'|/R_{H_n}\)):

| \(n\) | \(\mathrm{Pic}\) | \(h\) | \(w\) | \(h_{H_n}\) | \(\mathrm{Cl}(H_n)\) | \(\prod C_\chi(0)\) | index | \(= 2^a3^b\cdot\) | \(R_{H_n}\) | spare | status |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 3 | \(\mathbb{Z}/2\) | 2 | 12 | 1 | 1 | 1 | \(8\) | \(2^3\) | 1.3169578969 | 119 | certified |
| 5 | \(\mathbb{Z}/2\) | 2 | 4 | 1 | 1 | 1 | \(24\) | \(2^3\,3\) | 0.9624236501 | 118 | certified |
| 7 | \(\mathbb{Z}/4\) | 4 | 4 | 1 | 1 | 1 | \(13824\) | \(2^9\,3^3\) | 6.1139246305 | 116 | certified |
| 9 | \(\mathbb{Z}/6\) | 6 | 12 | 1 | 1 | 4 | \(10616832\) | \(2^{17}\,3^4\) | 98.483917427 | 113 | certified |
| 11 | \(\mathbb{Z}/6\) | 6 | 4 | 1 | 1 | 1 | \(7962624\) | \(2^{15}\,3^5\) | 446.47890419 | 113 | certified |
| 13 | \(\mathbb{Z}/6\) | 6 | 4 | 3 | \(\mathbb{Z}/3\) | 1 | \(23887872\) | \(2^{15}\,3^6\) | 393.80763206 | 112 | certified |
| 15 | \(\mathbb{Z}/4\times\mathbb{Z}/2\) | 8 | 12 | 2 | \(\mathbb{Z}/2\) | 32 | \(97844723712\) | \(2^{27}\,3^6\) | 4348.5075613 | 109 | certified |
| 23 | \(\mathbb{Z}/12\) | 12 | 4 | 12 | \(\mathbb{Z}/6\times\mathbb{Z}/2\) | 1 | \(18260173718028288\) | \(2^{35}\,3^{12}\) | \(1.0838443350\cdot10^9\) | 105 | GRH |

The defining polynomials (PARI's `polredbest`): \(y^4 - y^2 + 1\) (\(n = 3\)),
\(y^4 + 3y^2 + 1\) (\(5\)), \(y^8 - 2y^7 + 2y^6 + 4y^5 - y^4 + 4y^3 + 2y^2 - 2y + 1\) (\(7\)),
\(y^{12} - 2y^9 + 2y^6 + 2y^3 + 1\) (\(9\)), and at \(n = 11, 13, 15, 23\) the polynomials printed
by the script (degrees \(12, 12, 16, 24\)). Layer data (all `bnfcertify` = 1):

| \(n\) | quadratic: \(\mathbb{Q}(\sqrt{d_2})\), \(\theta^{(2)}_u = \pm\varepsilon^{e}\), index | cubic: \(L_3\), \(h_{L_3}\), \(\theta_u = \pm\eta^{e}\) | sextic \(F\): \(h_F\), \([\mathcal{O}_F^\times:\mu_4\mathcal{V}^A]\) |
|---|---|---|---|
| 3 | \(\mathbb{Q}(\sqrt3)\), \(e = 2\), 2 | — | — |
| 5 | \(\mathbb{Q}(\sqrt5)\), \(e = -12\), 12 | — | — |
| 7 | \(\mathbb{Q}(\sqrt7)\), \(e = 6\), 6 | — | — |
| 9 | \(\mathbb{Q}(\sqrt3)\), \(e = 8\), 8 | \(y^3 - 3y - 4\), \(1\), \(e = 8\) | \(1\), \(576\) |
| 11 | \(\mathbb{Q}(\sqrt{11})\), \(e = -6\), 6 | \(y^3 - y^2 + 4y + 2\), \(1\), \(e = 8\) | \(1\), \(576\) |
| 13 | \(\mathbb{Q}(\sqrt{13})\), \(e = -12\), 12 | \(y^3 - y^2 - 4y + 12\), \(3\), \(e = 24\) | \(3\), \(1728\) |
| 15 | \(\mathbb{Q}(\sqrt5)\) \(e = -48\), 48; \(\mathbb{Q}(\sqrt{15})\) \(e = 24\), 24; \(\mathbb{Q}(\sqrt3)\) \(e = 16\), 16 | — | — |
| 23 | \(\mathbb{Q}(\sqrt{23})\), \(e = 18\), 18 | \(y^3 - y^2 + 8y - 6\), \(2\), \(e = 16\) | \(4\), \(2304\) |

Guard rails: precision set in `main()`; every integer read from a real number
passes the absolute-error criterion with \(\ge\max(20,\mathrm{dps}/5)\) spare digits
(\(\ge 40\) for the index; observed \(\ge 82\) everywhere); exact arithmetic for the class
groups, characters, kernels, exponent matrices and Smith forms; no PSLQ. PARI's
`bnfinit` is GRH-conditional; the status column says where `bnfcertify` removed the
condition. A separate `bnfcertify` of the degree-\(24\) field \(H_{23}\) (same session,
realprecision 60) did not finish within 90 minutes of CPU time and was stopped; the
\(n = 23\) row — its class number \(12\), unit group and regulator — therefore remains
GRH-conditional, while the Dedekind-determinant side of that row (the \(\Delta\)-data,
the \(L'\)-values, the multipliers) is unconditional.

## 7. What is proved, what is certified, what is open

**Proved (every level \(n \ge 2\)).** Lemma 1.1; Proposition 2.2; Lemma 3.1 (primitive
level = conductor; the \(\gcd\)-closure); Lemma 3.2 (the Hecke recursion of the Euler
document as an identity of Dirichlet series, at every \(s\)); Corollary 3.3 (the
multipliers are positive algebraic integers, \(\prod_\chi C_\chi(0) \in \mathbb{Z}_{>0}\), and
\(L'(0,\chi) \ne 0\)); Lemma 4.1 (\(\mu(H_n)\)); **Theorem 1**; **Theorem 2** (layers) with
Corollaries 5.1–5.2; **Theorem 3** (the hyperbolic odd index, with the exact
\(2\)-adic invariant \(Q^-\)). Paper II's Conjecture 6.7 is thereby a theorem in its general
form, and its cubic-layer form (Theorem 6.6 there) is contained in Theorem 2 up to the
passage from \(F\) to its real subfield.

**Certified.** The record of §6 at \(n = 3, 5, 7, 9, 11, 13, 15\) (unconditional: PARI's
`bnfcertify`) and \(n = 23\) (GRH-conditional); the exact exponent matrices and Smith
forms; the saturation identity \([\mathcal{O}^\times_{H_n}:\mathcal{W}_n] = h_{H_n}\prod C_\chi(0)\) at all eight
levels (§5.5); the independent PARI confirmation of the eight-level cubic table of
Paper II at its four primitive levels; the hyperbolic record of §5.6 at the nine odd
levels \(5 \le n \le 21\) (unconditional except \(H\) at \(n = 21\)) with the sextic-layer index
\(27648\) at \(n = 21\) (unconditional).

**Open.** Conjecture 5.3 (the \(24\)-th-root saturation; the Galois-module statement
\(\mathcal{O}^\times/\mathcal{W}_n \cong \mathrm{Cl}(H_n)\) at primitive levels — Gras-type); the
comparison with Stark's units of \(H_n\) (Stark; Hajir–Rodriguez Villegas at prime
conductor), whose group has index \(h_{H_n}\) in \(\mathcal{O}_{H_n}^\times\) modulo torsion: Theorem 1
says \(\mathcal{V}_n\) has index \(4\cdot24^{h-1}/w\) times that, and §5.5 suggests the \(24\)-th
roots of the \(\Delta\)-quotients *are* a Stark-type system; the Kolyvagin bound
([schmidt-euler-system.md](schmidt-euler-system.md) §5), for which Theorem 1 supplies
the exact index that the bound should reproduce at every \(p \ge 5\); the law of the
\(2\)-adic invariant \(Q^-\) of Theorem 3 (\(= 2^{h/2-1}\) at eight levels, \(2^{h/2-2}\) at \(n = 19\)).

## 8. Outlook

- **The hyperbolic index beyond the first levels.** Theorem 3 reduces the hyperbolic
  Robert index to \(h_H/h_{H^+}\), the roots of unity, the multipliers and the \(2\)-adic
  invariant \(Q^- = [E:E^+E^-]\) of the quadratic extension \(H/H^{\sigma_\mathfrak{r}}\); the
  data (\(Q^- = 2^{h/2-1}\) except at \(n = 19\)) ask for the cohomological law of \(Q^-\)
  (\(\widehat H^0\)/\(H^1\) of \(\langle\sigma_\mathfrak{r}\rangle\) on the units, i.e. ambiguous
  classes of \(H/H^+\)), and for the levels \(n \ge 23\) with larger class groups.
- **Siegel units.** Prove Conjecture 5.3's first half: identify the \(24\)-th roots of the
  \(G_\mathfrak{c}/G_1\) with Siegel–Ramachandra units of the ring class field (Kubert–Lang,
  Ch. 11–12) and compute their torsion ambiguity, which is exactly the \(w/4\)-defect
  of the \(w = 12\) levels.
- **Other fields.** Theorem 1 is field-independent: for any imaginary quadratic \(K\) and
  any order \(\mathcal{O}\), the \(\Delta\)-quotients of proper \(\mathcal{O}\)-ideals have index
  \((4\cdot24^{h-1}/w)\,h_{H}\prod C_\chi(0)\cdot\)(the \(w_K/2\)-factor of \(\zeta_K(0) = -h_K/w_K\)) in the
  units of the ring class field — the constant \(4 = w_K\) of Theorem 1 is \(\zeta_K(0)^{-1}\).
  This is the index formula to test at \(\mathbb{Q}(\sqrt{-2})\) and \(\mathbb{Q}(\sqrt{-3})\)
  (item 3 of [PROGRAM.md](PROGRAM.md)).
- **Paper II.** The theorem replaces Conjecture 6.7; the verification table of §6 goes
  into §6 of the paper (done in this session's port).
