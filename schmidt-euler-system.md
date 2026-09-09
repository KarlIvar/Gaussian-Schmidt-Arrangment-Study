# The Schmidt \(\Delta\)-system as an Euler system: norm relations over the ring class tower of \(\mathbb{Q}(i)\)

This document carries out item 1 of the "beyond the Gaussian packing" program
of [outlook.md](outlook.md) §6: the Euclidean \(\Delta\)-data of
[phase-kronecker-limit.md](phase-kronecker-limit.md) Theorem 2 (Paper II,
Theorem 2.6),
$$
G_\mathfrak{c} \;=\; n^{12}\,\frac{\Delta(\Lambda_\mathfrak{c})}{\Delta(\mathbb{Z}[i])},
\qquad \mathfrak{c} \in \mathrm{Pic}(\mathcal{O}_n),\quad \mathcal{O}_n = \mathbb{Z} + n\mathbb{Z}[i],
$$
are **norm-compatible along the ring class tower** \(H_n \subset H_{n\ell}\), with
Euler factors of Heegner-point shape and the Eisenstein eigenvalue \(\ell + 1\)
in place of \(a_\ell\) (Theorem 1). Two consequences are proved: a Hecke-type
recursion for the Kronecker-limit values \(L'(0,\chi)\) of pulled-back
characters (Theorem 2), which explains every "imprimitive Euler multiplier"
recorded in [phase-kronecker-limit.md](phase-kronecker-limit.md) §4 and
[schmidt-units.md](schmidt-units.md) §5.2 — the correction \(1 + 3^{1-2s}\) at
\(n = 9\), the two unfactored real characters at \(n = 15\), and the column
\(C_n(0) = 2, 2, 2, 4\) of the Robert-index table — and the behaviour of the
cubic Robert index along the tower (Theorem 3):
$$
\bigl[\mathcal{O}_{L_3}^\times : \langle -1, \theta_u^{(n\ell)}\rangle\bigr]
\;=\; P_\ell(\chi_3)\,\bigl[\mathcal{O}_{L_3}^\times : \langle -1, \theta_u^{(n)}\rangle\bigr],
\qquad
P_\ell = \ell + 1 - \chi_3(\mathfrak{l}) - \chi_3(\mathfrak{l})^{-1}\ (\text{split}),\ \ \ell + 1\ (\text{inert}),\ \ 2\ (\ell = 2),
$$
with the Heegner recursion \((\ell+1)\,\log|\theta_u^{(n)}| - m\,\log|\theta_u^{(n/\ell)}|\)
at \(\ell \mid n\). The first step of Kolyvagin's descent for this system is
recorded (§5); the class-group bound it should produce is left open, with the
precise remaining steps listed. Everything displayed is machine-verified by
[scripts/schmidt_euler_system.py](scripts/schmidt_euler_system.py) — exact
Hermite-normal-form class-group arithmetic, \(\Delta\)-values at 80 digits, every
identity certified with \(\ge 77\) spare digits — under the guard rails of
[CLAUDE.md](CLAUDE.md); §7 states what is proved, what is certified, and what is
open.

Notation as in Paper II (`papers/2-schmidt-elliptic-units`): \(K = \mathbb{Q}(i)\),
\(\chi = \chi_{-4}\), \(N_e(n) = n\prod_{p\mid n}(1 - \chi(p)/p)\),
\(h(n) := |\mathrm{Pic}(\mathcal{O}_n)| = \tfrac12 N_e(n)\) for \(n \ge 2\); \(H_n\) the
ring class field of \(\mathcal{O}_n\); \(\Lambda_\mathfrak{c} \subseteq \mathbb{Z}[i]\) the
primitive index-\(n\) sublattice of the class \(\mathfrak{c}\), unique up to
\(\Lambda \mapsto i\Lambda\) (Paper I, Theorem 7.6); \(M(n) = \prod_\mathfrak{c} G_\mathfrak{c}\)
the \(\Delta\)-mass (Paper I, Theorem 7.18). Two facts from Paper II are used
throughout:

- **Galois action** (Theorem 2.6(1)): for \(\sigma \in \mathrm{Gal}(\bar{\mathbb{Q}}/K)\)
  with Artin class \(\mathfrak{c}(\sigma) \in \mathrm{Pic}(\mathcal{O}_n)\),
  \(\sigma(G_\mathfrak{c}) = G_{\mathfrak{c}(\sigma)^{-1}\mathfrak{c}}\); and
  \(\overline{G_\mathfrak{c}} = G_{\mathfrak{c}^{-1}}\).
- **Per-class valuations** (Theorem 4.11): \(v_\mathfrak{P}(G_\mathfrak{c}) = w_p(k)\) is
  independent of \(\mathfrak{c}\) at every place. Hence:

> **Corollary 0.1 (the Schmidt \(\Delta\)-units).** For every level \(n \ge 2\) and
> all \(\mathfrak{c}, \mathfrak{c}' \in \mathrm{Pic}(\mathcal{O}_n)\), the ratio
> \(v_{\mathfrak{c},\mathfrak{c}'} := G_\mathfrak{c}/G_{\mathfrak{c}'} = \Delta(\Lambda_\mathfrak{c})/\Delta(\Lambda_{\mathfrak{c}'})\)
> is a **unit** of \(H_n\). The group
> \(\mathcal{V}_n := \langle G_\mathfrak{c}/G_{\mathfrak{c}'}\rangle \subset \mathcal{O}_{H_n}^\times\)
> is \(\mathrm{Gal}(H_n/\mathbb{Q})\)-stable, and the cubic coset units
> \(\theta_u\) of Paper II §6 are elements of \(\mathcal{V}_n\).

(Immediate: the ideal \((G_\mathfrak{c})\) is the same for every class, and
\(\prod_\mathfrak{c} G_\mathfrak{c}/G_1 = M(n)/G_1^{\,h}\) shows \(\theta_u^{\,3} = \theta^3/M(n)
\in \mathcal{V}_n\); since \(\mathcal{V}_n\) is saturated in the unit group only up
to the index of Paper II's Conjecture 6.7, "\(\theta_u \in \mathcal{V}_n\)" is meant
up to the roots of unity \(\pm1\).) These are the classical \(\Delta\)-quotients
of proper \(\mathcal{O}_n\)-ideals; what is new below is the exact norm relation
they satisfy along the tower, and what it computes.

## 1. The tower and the fiber lemma

For a prime \(\ell\), \(\mathcal{O}_{n\ell} \subset \mathcal{O}_n\) gives \(H_n \subset H_{n\ell}\),
the projection \(\pi: \mathrm{Pic}(\mathcal{O}_{n\ell}) \to \mathrm{Pic}(\mathcal{O}_n)\),
\(\mathfrak{c}' \mapsto [\mathcal{O}_n\Lambda_{\mathfrak{c}'}]\), and
$$
\mathrm{Gal}(H_{n\ell}/H_n) \;\cong\; \ker\pi, \qquad
|\ker\pi| = \frac{N_e(n\ell)}{N_e(n)} =
\begin{cases} \ell - \chi(\ell), & \ell \nmid n,\\ \ell, & \ell \mid n\end{cases}
\qquad (n \ge 2).
$$

**Lemma 1.1 (local structure).** Let \(\Lambda\) be a primitive proper
\(\mathcal{O}_n\)-lattice of index \(n\) in \(\mathbb{Z}[i]\). For every prime \(p\),
\(\Lambda_p = \epsilon_p\,\mathcal{O}_{n,p}\) with \(\epsilon_p \in \mathbb{Z}_p[i]^\times\);
\(n\mathbb{Z}[i] \subseteq \Lambda\); and \(\Lambda \cap i\Lambda = n\mathbb{Z}[i]\).

*Proof.* The first statement is in the proof of Paper II, Theorem 4.11
(\(\Lambda_p\) is an invertible \(\mathcal{O}_{n,p}\)-module inside \(\mathbb{Z}_p[i]\) of the
right index, so a unit multiple of \(\mathcal{O}_{n,p}\)). Since \(n\mathbb{Z}[i] \subset
\mathcal{O}_n\) and \(\mathbb{Z}[i]\Lambda = \mathbb{Z}[i]\), \(n\mathbb{Z}[i] = n\mathbb{Z}[i]\Lambda
\subseteq \mathcal{O}_n\Lambda = \Lambda\). Finally \(\Lambda \cap i\Lambda\) is
\(i\)-stable, hence a \(\mathbb{Z}[i]\)-ideal containing \(n\mathbb{Z}[i]\), and locally
\(\mathcal{O}_{n,p} \cap i\mathcal{O}_{n,p} = \{a + nbi : n \mid a\} = n\mathbb{Z}_p[i]\), so it
has index \(n\) in \(\Lambda\), i.e. index \(n^2\) in \(\mathbb{Z}[i]\): it is \(n\mathbb{Z}[i]\). \(\square\)

**Lemma 1.2 (the fiber lemma).** Let \(\mathfrak{c} \in \mathrm{Pic}(\mathcal{O}_n)\),
\(\Lambda = \Lambda_\mathfrak{c}\), and \(\ell\) prime. Among the \(\ell + 1\)
sublattices \(\Lambda' \subset \Lambda\) of index \(\ell\):

1. *(\(\ell \nmid n\))* Exactly \(1 + \chi(\ell)\) are \(\mathcal{O}_n\)-stable, namely
   \(\Lambda \cap \mathfrak{l} = \mathfrak{l}_n\Lambda\) for the prime ideals
   \(\mathfrak{l} = \lambda\mathbb{Z}[i]\) of norm \(\ell\), where \(\mathfrak{l}_n := \mathfrak{l} \cap \mathcal{O}_n\)
   is the invertible \(\mathcal{O}_n\)-ideal of norm \(\ell\) with \(\mathfrak{l}_n\mathbb{Z}[i] = \mathfrak{l}\);
   they are not primitive.
2. *(\(\ell \mid n\))* Exactly one is \(\mathcal{O}_n\)-stable, namely \(\ell\,\Lambda^+\)
   with \(\Lambda^+ := \mathcal{O}_{n/\ell}\Lambda\), the primitive proper
   \(\mathcal{O}_{n/\ell}\)-lattice of index \(n/\ell\) in the class \(\pi(\mathfrak{c})\)
   (containing \(\Lambda\) with index \(\ell\)); it is not primitive.
3. The remaining \(\ell - \chi(\ell)\) (resp. \(\ell\)) sublattices are primitive proper
   \(\mathcal{O}_{n\ell}\)-lattices of index \(n\ell\), and \(\Lambda' \mapsto [\Lambda']\) is a
   **bijection** onto the fiber \(\pi^{-1}(\mathfrak{c}) \subset \mathrm{Pic}(\mathcal{O}_{n\ell})\).
4. The fiber is a single orbit of \(\mathrm{Gal}(H_{n\ell}/H_n)\), and for every
   \(\mathfrak{c}' \in \pi^{-1}(\mathfrak{c})\):
   \(\;N_{H_{n\ell}/H_n}(G_{\mathfrak{c}'}) = \prod_{\mathfrak{c}'' \in \pi^{-1}(\mathfrak{c})} G_{\mathfrak{c}''}\).

*Proof.* Index-\(\ell\) sublattices of \(\Lambda\) are the preimages of the lines of
\(\Lambda/\ell\Lambda \cong \Lambda_\ell/\ell\Lambda_\ell \cong \mathcal{O}_{n,\ell}/\ell\mathcal{O}_{n,\ell}\)
(Lemma 1.1, via \(\epsilon_\ell\)), and \(\Lambda'\) is \(\mathcal{O}_n\)-stable iff its line
is an \(\mathcal{O}_{n,\ell}\)-submodule. For \(\ell \nmid n\), \(\mathcal{O}_{n,\ell}/\ell =
\mathbb{Z}[i]/\ell\) is \(\mathbb{F}_\ell\times\mathbb{F}_\ell\), \(\mathbb{F}_{\ell^2}\), \(\mathbb{F}_2[x]/x^2\) in the
split, inert, ramified case, with \(2, 0, 1\) ideal lines; the non-ideal lines are
spanned by units \(\bar u\), \(u \in \mathbb{Z}_\ell[i]^\times\), and their preimages are
\(u(\mathbb{Z}_\ell + \ell\mathbb{Z}_\ell[i]) = u\,\mathcal{O}_{\ell,\ell}\). For \(\ell^k \| n\),
\(\mathcal{O}_{n,\ell}/\ell \cong \mathbb{F}_\ell[x]/x^2\) with \(x = \ell^ki\); the unique ideal line
\(\mathbb{F}_\ell x\) lifts to \(\epsilon(\ell\mathbb{Z}_\ell + \ell^k\mathbb{Z}_\ell[i]) = \ell\cdot\epsilon\,\mathcal{O}_{n/\ell,\ell}\),
and the \(\ell\) unit lines \(\mathbb{F}_\ell(1 + tx)\) lift to \(u\,\mathcal{O}_{n\ell,\ell}\),
\(u = 1 + t\ell^ki\) (both inclusions are one-line checks). In every case a
non-stable \(\Lambda'\) agrees with \(\Lambda\) at \(p \ne \ell\) and is a unit
multiple of \(\mathcal{O}_{n\ell,\ell}\) at \(\ell\): a locally principal, hence proper,
\(\mathcal{O}_{n\ell}\)-lattice, primitive (its \(\mathbb{Z}[i]\)-span is \(\mathbb{Z}[i]\) locally
everywhere), of index \(n\ell\). The stable ones are as stated: for \(\ell \nmid n\),
\(\mathfrak{l}_n\Lambda \subseteq \Lambda \cap \mathfrak{l}\) and both have index \(\ell\) in
\(\Lambda\) (\(\mathfrak{l}_n\) is invertible of norm \(\ell\); \(\mathbb{Z}[i]/\mathfrak{l}
\cong \mathcal{O}_n/\mathfrak{l}_n\) since \(\ell \nmid n\)), and \(\mathfrak{l}_n\Lambda \subseteq
\mathfrak{l} \ne \mathbb{Z}[i]\) is not primitive; for \(\ell \mid n\), \(\Lambda^+ =
\mathcal{O}_{n/\ell}\Lambda\) is locally \(\epsilon\,\mathcal{O}_{n/\ell,p}\), so primitive,
proper for \(\mathcal{O}_{n/\ell}\), of index \(n/\ell\), containing \(\Lambda\) with
index \(\ell\), and \(\ell\Lambda^+ \subseteq \ell\mathbb{Z}[i]\).

(3) Injectivity: if \(\Lambda'_1, \Lambda'_2 \subset \Lambda\) are non-stable in the same
class, \(\Lambda'_2 = \alpha\Lambda'_1\) with \(\alpha \in K^\times\); primitivity of
both forces \(\alpha\mathbb{Z}[i] = \mathbb{Z}[i]\), so \(\alpha \in \mu_4\), and if
\(\Lambda'_2 = i\Lambda'_1\) then \(\Lambda'_1 \subseteq \Lambda \cap i\Lambda = n\mathbb{Z}[i]\)
(Lemma 1.1), contradicting primitivity for \(n \ge 2\). Surjectivity: for
\(\mathfrak{c}' \in \pi^{-1}(\mathfrak{c})\) with primitive representative \(\Lambda'\),
\(\mathcal{O}_n\Lambda'\) is locally \(\epsilon'\mathcal{O}_{n,p}\), hence a primitive
proper \(\mathcal{O}_n\)-lattice of index \(n\) in the class \(\mathfrak{c}\), i.e.
\(\mathcal{O}_n\Lambda' = \alpha\Lambda\) with \(\alpha \in \mu_4\); replacing \(\Lambda'\)
by \(\alpha^{-1}\Lambda'\) (same class) gives \(\Lambda' \subset \Lambda\) of index
\(\ell\), non-stable because primitive. The counts agree with \(|\ker\pi|\).

(4) \(\ker\pi\) acts freely and transitively on the fiber by translation, and by the
Galois action of Paper II, Theorem 2.6(1), \(\sigma \in \mathrm{Gal}(H_{n\ell}/H_n)\)
— i.e. Artin class \(\mathfrak{c}(\sigma) \in \ker\pi\), the Artin map of \(H_n\) being
\(\pi\) composed with that of \(H_{n\ell}\) — sends \(G_{\mathfrak{c}'}\) to
\(G_{\mathfrak{c}(\sigma)^{-1}\mathfrak{c}'}\). \(\blacksquare\)

(Verified: at all 42 pairs \((n,\ell)\) of §8 the fiber sizes are as stated and the
projection is computed by exact HNF arithmetic.)

## 2. The norm relations

**Theorem 1 (norm relations).** Let \(n \ge 2\), \(\ell\) prime, \(\mathfrak{c} \in
\mathrm{Pic}(\mathcal{O}_n)\), and \(\mathfrak{c}' \in \pi^{-1}(\mathfrak{c})\). Write \(N =
N_{H_{n\ell}/H_n}\), \([\mathfrak{l}]\mathfrak{c}\) for the class of \(\mathfrak{l}_n\Lambda_\mathfrak{c}\),
\(\mathfrak{c}^+ = \pi(\mathfrak{c}) \in \mathrm{Pic}(\mathcal{O}_{n/\ell})\) when \(\ell \mid n\), and
\(G^{(1)} := 1\). Then
$$
N(G_{\mathfrak{c}'}) \;=\;
\begin{cases}
\dfrac{G_\mathfrak{c}^{\,\ell+1}}{G_{[\mathfrak{l}]\mathfrak{c}}\,G_{[\bar{\mathfrak{l}}]\mathfrak{c}}},
& \ell = \lambda\bar\lambda \text{ split},\ \ell \nmid n,\\[10pt]
\ell^{12}\,G_\mathfrak{c}^{\,\ell+1}, & \ell \text{ inert},\ \ell \nmid n,\\[6pt]
2^{6}\,\dfrac{G_\mathfrak{c}^{\,3}}{G_{[\mathfrak{p}]\mathfrak{c}}}, & \ell = 2 \nmid n,\ \ \mathfrak{p} = (1+i),\\[10pt]
(-1)^{[\ell = 2]}\;\dfrac{G_\mathfrak{c}^{\,\ell+1}}{G^{(n/\ell)}_{\mathfrak{c}^+}}, & \ell \mid n .
\end{cases}
$$

*Proof.* Step M1 of the \(\Delta\)-mass law (Paper I, Theorem 7.18; the full product
over all index-\(\ell\) sublattices of any lattice \(L\) is the constant
\(A(\ell) = (-1)^{t(\ell)}\prod_{d\mid\ell}d^{-12d}\)) gives, for prime \(\ell\),
$$
\prod_{[\Lambda:\Lambda']=\ell}\frac{\Delta(\Lambda')}{\Delta(\Lambda)} \;=\; A(\ell) = (-1)^{[\ell=2]}\,\ell^{-12\ell},
\qquad \Lambda = \Lambda_\mathfrak{c}.
$$
Split the product according to Lemma 1.2. Each non-stable \(\Lambda'\) is the primitive
representative of a fiber class \(\mathfrak{c}''\), so \(\Delta(\Lambda')/\Delta(\Lambda) =
\bigl[G_{\mathfrak{c}''}/(n\ell)^{12}\bigr]\big/\bigl[G_\mathfrak{c}/n^{12}\bigr] = G_{\mathfrak{c}''}/(\ell^{12}G_\mathfrak{c})\),
and by Lemma 1.2(4) the product over them is \(N(G_{\mathfrak{c}'})/(\ell^{12}G_\mathfrak{c})^{|\ker\pi|}\).
For a stable sublattice with \(\ell \nmid n\): \(\mathfrak{l}_n\Lambda = \lambda\Lambda_\lambda\)
with \(\Lambda_\lambda := \lambda^{-1}\mathfrak{l}_n\Lambda \subseteq \mathbb{Z}[i]\) primitive
(its \(\mathbb{Z}[i]\)-span is \(\lambda^{-1}\mathfrak{l} = \mathbb{Z}[i]\)) of index \(n\ell/\ell = n\)
in the class \([\mathfrak{l}]\mathfrak{c}\); by homogeneity
\(\Delta(\mathfrak{l}_n\Lambda)/\Delta(\Lambda) = \lambda^{-12}\,G_{[\mathfrak{l}]\mathfrak{c}}/G_\mathfrak{c}\).
For \(\ell \mid n\): \(\Delta(\ell\Lambda^+)/\Delta(\Lambda) = \ell^{-12}\Delta(\Lambda^+)/\Delta(\Lambda)
= \ell^{-12}\bigl[G^{(n/\ell)}_{\mathfrak{c}^+}/(n/\ell)^{12}\bigr]\big/\bigl[G_\mathfrak{c}/n^{12}\bigr]
= G^{(n/\ell)}_{\mathfrak{c}^+}/G_\mathfrak{c}\). Assembling:

- split, \(\ell \nmid n\) (\(|\ker\pi| = \ell - 1\)):
  \(N(G_{\mathfrak{c}'}) = \ell^{12(\ell-1)}G_\mathfrak{c}^{\ell-1}\,A(\ell)\,(\lambda\bar\lambda)^{12}\,G_\mathfrak{c}^{2}/(G_{[\mathfrak{l}]\mathfrak{c}}G_{[\bar{\mathfrak{l}}]\mathfrak{c}})\),
  and \(\ell^{12(\ell-1)}\cdot\ell^{-12\ell}\cdot\ell^{12} = 1\);
- inert (\(|\ker\pi| = \ell+1\), no stable sublattice):
  \(N(G_{\mathfrak{c}'}) = \ell^{12(\ell+1)}\ell^{-12\ell}G_\mathfrak{c}^{\ell+1} = \ell^{12}G_\mathfrak{c}^{\ell+1}\);
- \(\ell = 2 \nmid n\) (\(|\ker\pi| = 2\), one stable sublattice, \((1+i)^{12} = -64\)):
  \(N(G_{\mathfrak{c}'}) = 2^{24}G_\mathfrak{c}^{2}\cdot(-2^{-24})\cdot(-64)\,G_\mathfrak{c}/G_{[\mathfrak{p}]\mathfrak{c}} = 2^6G_\mathfrak{c}^3/G_{[\mathfrak{p}]\mathfrak{c}}\);
- \(\ell \mid n\) (\(|\ker\pi| = \ell\)):
  \(N(G_{\mathfrak{c}'}) = \ell^{12\ell}G_\mathfrak{c}^{\ell}\,A(\ell)\,G_\mathfrak{c}/G^{(n/\ell)}_{\mathfrak{c}^+}
  = (-1)^{[\ell=2]}G_\mathfrak{c}^{\ell+1}/G^{(n/\ell)}_{\mathfrak{c}^+}\). \(\blacksquare\)

Verified at every class of the 42 pairs \((n,\ell)\) listed in §8 — all four
cases, including \(\ell \mid n\) with \(n/\ell > 1\) (e.g. \(6\to18\), \(9\to27\),
\(18\to54\), \(15\to45, 75\), \(25\to125\), \(14\to98\)) — with \(\ge 77\) spare
digits at 80.

**Corollary 2.1 (Galois form; the Eisenstein–Heegner shape).** Let
\(\mathrm{Fr}_\mathfrak{l} = (\mathfrak{l}_n, H_n/K)\); by Theorem 2.6(1) of Paper II,
\(G_{[\mathfrak{l}]\mathfrak{c}} = \mathrm{Fr}_\mathfrak{l}^{-1}(G_\mathfrak{c})\), and \(\mathrm{Fr}_{\bar{\mathfrak{l}}}
= \mathrm{Fr}_\mathfrak{l}^{-1}\) (\(\mathfrak{l}_n\bar{\mathfrak{l}}_n = (\ell)\)). In additive
notation on \(H_n^\times\otimes\mathbb{Q}\) the relations read
$$
N_{H_{n\ell}/H_n}\,G_{\mathfrak{c}'} \;=\;
\begin{cases}
(\ell + 1 - \mathrm{Fr}_\mathfrak{l} - \mathrm{Fr}_\mathfrak{l}^{-1})\,G_\mathfrak{c} & (\ell \text{ split}),\\
(\ell + 1)\,G_\mathfrak{c} + 12\log\ell & (\ell \text{ inert}),\\
(3 - \mathrm{Fr}_\mathfrak{p})\,G_\mathfrak{c} + 6\log 2 & (\ell = 2),\\
(\ell + 1)\,G_\mathfrak{c} - G^{(n/\ell)}_{\mathfrak{c}^+} & (\ell \mid n),
\end{cases}
$$
and the rational constants disappear on the unit ratios of Corollary 0.1:
\(N(v_{\mathfrak{c}'_1,\mathfrak{c}'_2}) = P_\ell\cdot v_{\mathfrak{c}_1,\mathfrak{c}_2}\) with the
same operators. These are **exactly the norm relations of Heegner points**
(Gross, *Kolyvagin's work on modular elliptic curves*, Prop. 3.7:
\(a_\ell\,y_n\) at \(\ell \nmid n\) inert, \((a_\ell - \sigma_\lambda - \sigma_{\bar\lambda})y_n\)
at split \(\ell\), \(a_\ell\,y_n - y_{n/\ell}\) at \(\ell \mid n\)) with the Hecke eigenvalue
\(a_\ell\) replaced by the **Eisenstein eigenvalue** \(\ell + 1\): the Schmidt
\(\Delta\)-system is the Eisenstein member of the anticyclotomic (ring-class)
family of Euler systems over \(\mathbb{Q}(i)\). It is *not* an Euler system for
\(\mathbb{Z}_p(1)\) in Rubin's sense (those need \((1 - \mathrm{Fr}^{-1})\) along ray class
fields); the relevant machine is Kolyvagin's original descent, §5.

**Remark 2.2 (the level \(n = 1\) is the \(\Delta\)-mass law).** At \(n = 1\) the
non-stable sublattices of \(\mathbb{Z}[i]\) represent each class of \(\mathrm{Pic}(\mathcal{O}_\ell)\)
twice, so the same computation gives \(M(\ell)^2 = \ell^{12}\) (inert), \(1\) (split),
\(2^6\) (\(\ell = 2\)) — the prime levels of Paper I's mass table
(\(M(3) = 3^6, M(7) = 7^6, M(5) = 1, M(2) = 2^3\)).

## 3. The Hecke recursion for \(L'(0,\chi)\) and the imprimitive Euler factors

Let \(\chi\) be a nontrivial character of \(\mathrm{Pic}(\mathcal{O}_n)\) and
\(\chi^{(n\ell)} = \chi\circ\pi\) its pullback to level \(n\ell\); write
\(S_\chi(m) = \sum_{\mathfrak{c}\in\mathrm{Pic}(\mathcal{O}_m)}\chi(\mathfrak{c})\log|G_\mathfrak{c}|
= -12\,L'(0,\chi)\) (Paper II, Theorem 2.6(2)), and \(m_\pi := |\ker(\mathrm{Pic}(\mathcal{O}_n)
\to \mathrm{Pic}(\mathcal{O}_{n/\ell}))|\) when \(\ell \mid n\).

> **Theorem 2 (Hecke recursion).** With \(P_\ell(\chi) = \ell + 1 - \chi(\mathfrak{l}) - \chi(\mathfrak{l})^{-1}\)
> (split), \(\ell + 1\) (inert), \(3 - \chi(\mathfrak{p})\) (\(\ell = 2\)):
> $$
> L'(0,\chi^{(n\ell)}) = P_\ell(\chi)\,L'(0,\chi) \quad (\ell \nmid n), \qquad
> L'(0,\chi^{(n\ell)}) = (\ell+1)\,L'(0,\chi) - m_\pi\,L'(0,\chi^{(n/\ell)}) \quad (\ell \mid n),
> $$
> where the last term is present only if \(\chi\) factors through \(\mathrm{Pic}(\mathcal{O}_{n/\ell})\)
> (and \(\chi^{(n/\ell)}\) denotes the character it comes from).

*Proof.* Sum \(\chi(\mathfrak{c})\log|\cdot|\) of Theorem 1 over \(\mathfrak{c}\):
\(S_{\chi}(n\ell) = \sum_\mathfrak{c}\chi(\mathfrak{c})\log|N(G_{\mathfrak{c}'})|\). The rational
constants are killed by \(\sum\chi = 0\); \(\sum_\mathfrak{c}\chi(\mathfrak{c})\log|G_{[\mathfrak{l}]\mathfrak{c}}|
= \chi(\mathfrak{l})^{-1}S_\chi(n)\) by re-indexing; and for \(\ell \mid n\),
\(\sum_\mathfrak{c}\chi(\mathfrak{c})\log|G^{(n/\ell)}_{\pi(\mathfrak{c})}| = \sum_{\mathfrak{c}^+}
\log|G^{(n/\ell)}_{\mathfrak{c}^+}|\sum_{\mathfrak{c}\in\pi^{-1}(\mathfrak{c}^+)}\chi(\mathfrak{c})\), which is
\(0\) if \(\chi\) is nontrivial on \(\ker\pi\) and \(m_\pi S_{\chi^{(n/\ell)}}(n/\ell)\)
otherwise. Divide by \(-12\). \(\blacksquare\)

**Corollary 3.1 (the imprimitive Euler multipliers, proved).** In Paper II's
notation \(L(s,\chi^{(m)}) = L_{\mathrm{prim}}(s,\chi)\,C_m(s)\) (Lemma 6.3 there),
the values \(C_m(0)\) along any pullback chain are the products of the multipliers
of Theorem 2. In particular the following facts, recorded as *certified* or
*observed* in [phase-kronecker-limit.md](phase-kronecker-limit.md) §4(a) and
Paper II (Prop. 3.1, Remark 3.2, Theorem 6.6), are theorems:

| pullback | case | multiplier | consequence |
|---|---|---|---|
| \(\chi_2\): \(3 \to 9\) | \(\ell = 3 \mid 3\), \(\chi_2\) primitive at 3 | \(\ell + 1 = 4\) | \(L'(0,\chi_2^{(9)}) = 4\cdot\tfrac13\log\varepsilon_{12}\): the correction \(1 + 3^{1-2s}\), \(C(0) = 4\) |
| \(\chi_2\): \(3 \to 15\) | \(\ell = 5\) split, \(\chi_2(\mathfrak{l}_3) = -1\) | \(5 + 1 + 2 = 8\) | \(L'(0,\chi_{(2,1)}^{(15)}) = \tfrac83\log\varepsilon_{12}\) — the "unfactored" character of Remark 3.2 |
| \(\chi_2\): \(5 \to 15\) | \(\ell = 3\) inert | \(4\) | \(L'(0,\chi_{(0,1)}^{(15)}) = 8\log\varepsilon_5\) — the other one |
| \(\chi_3\): \(9 \to 18\), \(11 \to 22\), \(13 \to 26\) | \(\ell = 2\), \(\chi_3(\mathfrak{p}) = 1\) | \(2\) | \(C_n(0) = 2\) |
| \(\chi_3\): \(9 \to 27\) | \(\ell = 3 \mid 9\), \(\chi_3\) primitive at 9 | \(4\) | \(C_{27}(0) = 4\) |
| \(\chi_3\): \(9 \to 45\), \(11 \to 55\), \(13 \to 65\) | \(\ell = 5\) split, \(\chi_3(\mathfrak{l}) \ne 1\) | \(7\) | new |
| \(\chi_3\): \(9 \to 63\) | \(\ell = 7\) inert | \(8\) | new |
| \(\chi_3\): \(11 \to 33\), \(13 \to 39\), \(23 \to 69\) | \(\ell = 3\) inert | \(4\) | new |

(All entries verified directly as identities between the character sums
\(S_\chi\) at the two levels, and the first three against the proved genus closed
forms.) The "two-prime factorization shapes" left open in
[phase-kronecker-limit.md](phase-kronecker-limit.md) §8 (open 2) are thus not
genus factorizations at all: the two characters are pullbacks, and their
\(L'\)-values are Eisenstein multiples of the lower level's.

Note that for a cubic character and \(\ell \ne 3\) the pullback is automatic:
\(\ker\pi\) has order \(\ell - \chi(\ell)\) or \(\ell\), prime to \(3\) unless
\(\ell \equiv \chi(\ell) \bmod 3\) (i.e. \(\ell \equiv 1, 11 \bmod 12\)) or \(\ell = 3\) with
\(\ell^2 \mid n\); this is why the primitive cubic levels are \(9 = 3^2\) and
\(11, 13, 23 \equiv \pm1 \bmod 12\).

## 4. The Robert index along the tower

Let \(\chi_3\) be a cubic character of \(\mathrm{Pic}(\mathcal{O}_n)\) with kernel \(K_n\),
\(L_3\) the real cubic field it cuts out, \(\theta^{(m)} = \prod_{\mathfrak{c}\in K_m}G_\mathfrak{c}\)
the principal coset product at a level \(m\) of the pullback chain, and
\(\theta_u^{(m)} = \theta^{(m)}/|M(m)|^{1/3}\) the coset unit of Paper II,
Proposition 6.2 (the same field \(L_3\) at every level of the chain, since the
pulled-back character cuts out the same extension).

> **Theorem 3 (index along the tower).** With \(P_\ell = P_\ell(\chi_3)\) as in Theorem 2
> (\(P_2 = 2\), since \(\chi_3(\mathfrak{p}) = 1\)):
> $$
> \log|\theta_u^{(n\ell)}| = P_\ell\,\log|\theta_u^{(n)}| \quad (\ell \nmid n),\qquad
> \log|\theta_u^{(n\ell)}| = (\ell+1)\log|\theta_u^{(n)}| - m_\pi\log|\theta_u^{(n/\ell)}| \quad (\ell \mid n),
> $$
> the last term present iff \(\chi_3\) factors through \(\mathrm{Pic}(\mathcal{O}_{n/\ell})\).
> Consequently
> \(\bigl[\mathcal{O}_{L_3}^\times : \langle -1, \theta_u^{(n\ell)}\rangle\bigr]
> = P_\ell\,\bigl[\mathcal{O}_{L_3}^\times : \langle -1, \theta_u^{(n)}\rangle\bigr]\) for
> \(\ell \nmid n\), and the Heegner recursion for \(\ell \mid n\); Paper II's column
> \(C_n(0) = 2, 2, 2, 4\) at \(n = 18, 22, 26, 27\) is proved, and the index at every
> level of the four pullback chains below is \(8\,h_{L_3}\) times the tabulated
> multiplier.

*Proof.* \(K_{n\ell} = \pi^{-1}(K_n)\), so \(\theta^{(n\ell)} = \prod_{\mathfrak{c}\in K_n}N(G_{\mathfrak{c}'})\).
Insert Theorem 1 and use: \(\prod_{\mathfrak{c}\in K_n}G_{[\mathfrak{l}]\mathfrak{c}}\) is the coset
product \(\theta_{[\mathfrak{l}]K_n}\), equal to \(\theta^{(n)}\) if \(\chi_3(\mathfrak{l}) = 1\) and
to one of the two conjugate coset products otherwise, with
\(\theta\,\theta'\,\theta'' = M(n)\); for \(\ell \mid n\),
\(\prod_{\mathfrak{c}\in K_n}G^{(n/\ell)}_{\pi(\mathfrak{c})} = (\theta^{(n/\ell)})^{m_\pi}\) if \(\chi_3\)
factors and \(= M(n/\ell)^{m_\pi/3}\) otherwise (then \(K_n\) surjects onto
\(\mathrm{Pic}(\mathcal{O}_{n/\ell})\) with fibers of size \(m_\pi/3\)). This gives
\(\theta^{(n\ell)} = \theta^{\ell-1}\) or \(\theta^{\ell+2}/M(n)\) (split, according
to \(\chi_3(\mathfrak{l}) = 1\) or not), \(\ell^{12h/3}\theta^{\ell+1}\) (inert),
\(2^{2h}\theta^{2}\) (\(\ell = 2\); \([\mathfrak{p}]\) is \(2\)-torsion, so in \(K_n\)),
and \(\pm\theta^{\ell+1}/(\theta^{(n/\ell)})^{m_\pi}\) resp. \(\pm\theta^{\ell+1}/M(n/\ell)^{m_\pi/3}\)
(\(\ell \mid n\)). The mass law (Paper I, Theorem 7.18) gives the exponent
identities \(v_p(M(n\ell)) = (\ell - \chi(\ell))\,v_p(M(n))\) for \(p \mid n\),
\(v_\ell(M(n\ell)) = 12h, 6h, 0\) (inert, \(\ell = 2\), split) when \(\ell \nmid n\), and
\(v_p(M(n\ell)) = (\ell+1)v_p(M(n)) - m_\pi v_p(M(n/\ell))\) at every \(p\) when
\(\ell \mid n\) (at \(p = \ell\): \((\ell+1)(\ell^k - 1) - m_\pi(\ell^{k-1} - 1) = \ell^{k+1} - 1\)
with \(m_\pi = \ell\) for \(k \ge 2\), \(= \ell - \chi(\ell)\) for \(k = 1\); at \(p \ne \ell\):
\(N_e\) multiplies by \(m_\pi\) and then by \(\ell\)). Substituting \(\log|\theta| =
\log|\theta_u| + \tfrac13\log|M|\) at each level, the mass terms cancel identically
in every case and the displayed relations remain. The index statement is
Paper II, Theorem 6.6: the index is \(\log|\theta_u|/R_{L_3}\). \(\blacksquare\)

**The computed chains** (index multipliers relative to the base level; each
certified as an integer with \(\ge 79\) spare digits at 80; base indices
\(8, 8, 24, 16\) from Paper II):

| base \(n\) (\(h_{L_3}\)) | level | step | multiplier | index |
|---|---|---|---|---|
| 9 (1) | 18 | \(\ell = 2\) | 2 | 16 |
| | 36 | \(2 \mid 18\): \(3\log\theta_u^{(18)} - 2\log\theta_u^{(9)}\) | 4 | 32 |
| | 27 | \(3 \mid 9\), \(\chi_3\) primitive at 9 | 4 | 32 |
| | 81 | \(3 \mid 27\): \(4\log\theta_u^{(27)} - 3\log\theta_u^{(9)}\) | 13 | 104 |
| | 45 | \(\ell = 5\) split, \(\chi_3(\mathfrak{l}) = \omega\) | 7 | 56 |
| | 63 | \(\ell = 7\) inert | 8 | 64 |
| 11 (1) | 22 | \(\ell = 2\) | 2 | 16 |
| | 44 | \(2 \mid 22\) | 4 | 32 |
| | 33 | \(\ell = 3\) inert | 4 | 32 |
| | 55 | \(\ell = 5\) split, \(\chi_3(\mathfrak{l}) = \omega\) | 7 | 56 |
| 13 (3) | 26 | \(\ell = 2\) | 2 | 48 |
| | 52 | \(2 \mid 26\) | 4 | 96 |
| | 39 | \(\ell = 3\) inert | 4 | 96 |
| | 65 | \(\ell = 5\) split, \(\chi_3(\mathfrak{l}) = \omega^2\) | 7 | 168 |
| 23 (2) | 46 | \(\ell = 2\) | 2 | 32 |
| | 69 | \(\ell = 3\) inert | 4 | 64 |

The multipliers \(1, 4, 13\) along \(9, 27, 81\) are \(1 + 3 + 9\)-type sums, the
Eisenstein Hecke eigenvalues of \(T_{3^k}\) on the ring class tower: the
\(3\)-adic ladder of the index is \((3^{k+1} - 1)/2\).

## 5. Toward Kolyvagin: the derivative classes

Fix an odd prime \(p \ge 5\), a level \(n \ge 2\), and an inert prime \(\ell \nmid n\)
with \(\ell \equiv -1 \pmod p\). Then \(G := \mathrm{Gal}(H_{n\ell}/H_n)\) is cyclic of
order \(\ell + 1 \equiv 0 \pmod p\), and \(\ell\) splits completely in \(H_n/K\)
(the ideal \((\ell)\) is principal in \(\mathcal{O}_n\)). Let \(\sigma\) generate \(G\),
\(D_\ell = \sum_{i=1}^{\ell}i\sigma^i\) Kolyvagin's derivative, so that
\((\sigma - 1)D_\ell = \ell + 1 - N_G\).

**Lemma 5.1 (descent of the derivative classes).** Let \(v = v_{\mathfrak{c}'_1,\mathfrak{c}'_2}
\in \mathcal{V}_{n\ell}\) be a unit ratio and \(w = v_{\mathfrak{c}_1,\mathfrak{c}_2} \in \mathcal{V}_n\)
the ratio of the projected classes. Then \((\sigma - 1)D_\ell v = (v/w)^{\ell+1}\) is a
\(p\)-th power in \(H_{n\ell}^\times\), the class of \(D_\ell v\) in
\(H_{n\ell}^\times/(H_{n\ell}^\times)^p\) is \(G\)-invariant, and it descends to a
unique class \(\kappa_{n,\ell}(v) \in H_n^\times/(H_n^\times)^p\).

*Proof.* Theorem 1 (inert case) on the ratios gives \(N_G(v) = w^{\ell+1}\), so
\((\sigma - 1)D_\ell v = v^{\ell+1}/N_G(v) = (v/w)^{\ell+1}\), a \(p\)-th power as
\(p \mid \ell + 1\). For the descent: \(H_{n\ell}\) is generalized dihedral over
\(\mathbb{Q}\), so its maximal subfield abelian over \(\mathbb{Q}\) is an elementary
\(2\)-extension, which cannot contain \(\mathbb{Q}(\zeta_p)\) for \(p \ge 5\); hence
\(\mu_p(H_{n\ell}) = 1\), and inflation–restriction in Kummer theory identifies
\(H_n^\times/(H_n^\times)^p = H^1(H_n, \mu_p)\) with
\(\bigl(H^1(H_{n\ell},\mu_p)\bigr)^G = \bigl(H_{n\ell}^\times/(H_{n\ell}^\times)^p\bigr)^G\).
\(\square\)

This is the first step of Kolyvagin's argument for the Eisenstein–Heegner system.
What remains — and is **not** done here — is the analysis of \(\kappa_{n,\ell}(v)\)
at the primes of \(H_n\) above \(\ell\) (its singular part is governed by the
reduction of \(w\) modulo \(\ell\), through the finite–singular comparison of
Kolyvagin/Gross), from which one bounds the \(p\)-part of the class group of
\(H_n\) (or of the \(\chi\)-components on which \(\mathcal{V}_n\) has rank one) by the
index of \(\mathcal{V}_n\); the expected statement is the anticyclotomic analogue of
the Thaine–Rubin bound, with Küçüksakallı's class numbers of ring class fields
of prime conductor as the test data. The comparison of \(\mathcal{V}_n\) with Stark's
units of \(H_n\) (index \(h_{H_n}\)) would give the relation between the two
indices.

## 6. Comparison with the classical elliptic units

The units \(v_{\mathfrak{c},\mathfrak{c}'} = \Delta(\Lambda_\mathfrak{c})/\Delta(\Lambda_{\mathfrak{c}'})\) are
\(\Delta\)-quotients of proper \(\mathcal{O}_n\)-ideals, the raw material of every
construction of elliptic units in ring class fields (Siegel–Ramachandra–Robert
via Siegel functions; Stark; Hajir–Rodriguez Villegas; Schertz); the norm relation
of Theorem 1 is the distribution relation of \(\Delta\) (the step M1 of Paper I)
read along the tower, and the same distribution relation underlies the
norm-compatibility of the classical elliptic units (de Shalit, *Iwasawa theory of
elliptic curves with complex multiplication*, Ch. II §2; Rubin, Invent. Math. 103
(1991), for ray class fields). What this document adds is not the existence of a
norm-compatible \(\Delta\)-system but its exact form on the ring class tower of
\(\mathbb{Q}(i)\) with the Eisenstein–Heegner operators, and the two consequences
(Theorems 2 and 3) that turn certified entries of Paper II into theorems. The
exact relation between \(\mathcal{V}_n\) and Stark's unit of \(H_n\) is open (§5).

## 7. What is proved, what is certified, what is open

**Proved (all levels \(n \ge 2\), all primes \(\ell\)).** Corollary 0.1; Lemmas 1.1–1.2;
Theorem 1; Corollary 2.1; Theorem 2 and Corollary 3.1 (hence: the conductor
correction at \(n = 9\), the two \(n = 15\) identifications of Paper II Remark 3.2,
and the Euler multipliers \(C_n(0)\) at every pullback level); Theorem 3;
Lemma 5.1.

**Certified (exact HNF arithmetic; \(\Delta\)-identities with \(\ge 77\) spare digits
at 80).** The 42 norm-relation instances; the 10 recursion instances of Theorem 2;
the 16 tower levels of Theorem 3 (index multipliers as certified integers).

**Open.** The Kolyvagin bound (§5); the comparison with Stark's units and the
full index of \(\mathcal{V}_n\) in \(\mathcal{O}_{H_n}^\times\) (item 2 of the program,
[outlook.md](outlook.md) §6); the transport of Theorem 1 to other imaginary
quadratic fields (item 3 there), where the fiber lemma needs the cusps
\(\leftrightarrow\) ideal classes bookkeeping and the Frobenius factors depend on
splitting in \(K\).

## 8. Machine verification

`python3 scripts/schmidt_euler_system.py --selftest` (~2 min on one core; 80
digits). The class groups \(\mathrm{Pic}(\mathcal{O}_m)\) are built exactly: canonical
primitive index-\(m\) sublattices in Hermite normal form (each class as
\(\min(\Lambda, i\Lambda)\)), products, inverses, orders, projections
\(\Lambda \mapsto \mathcal{O}_{n}\Lambda\), ideal twists \(\Lambda \mapsto
\lambda^{-1}(\Lambda\cap\lambda\mathbb{Z}[i])\), kernels and cubic characters — no
numerics anywhere in the bookkeeping. The script asserts:

1. **(T1/T2)** at the 42 pairs \((n,\ell)\): \(n \in \{3,7,3,9,11,2,6,4\}\) with split
   \(\ell \in \{5, 13\}\); \(n \in \{3,5,7,2,4,11,13,10,5\}\) with inert \(\ell \in \{3,7\}\);
   \(n \in \{3,5,7,9,11,13,15,25\}\) with \(\ell = 2\); and the \(\ell \mid n\) pairs
   \((2,2),(3,3),(5,5),(4,2),(6,2),(6,3),(9,3),(10,2),(10,5),(12,2),(15,3),(15,5),(8,2),
   (18,3),(25,5),(14,2),(14,7)\) (top level 125): the fiber sizes of Lemma 1.2 and
   the four-case relation of Theorem 1, every class, relative error
   \(\le 10^{-77}\);
2. **(T3)** the recursion of Theorem 2 for \(\chi_2\) at \(3 \to 9\), \(3 \to 15\),
   \(5 \to 15\) — against the proved closed forms \(-4\log\varepsilon_{12}\),
   \(-16\log\varepsilon_{12}\), \(-32\log\varepsilon_{12}\), \(-24\log\varepsilon_5\),
   \(-96\log\varepsilon_5\) of the sums \(S_\chi = -12L'\) — and for \(\chi_3\) at
   \(9\to45, 63, 18\), \(11\to33, 55\), \(13\to65, 26\), with the multiplier computed from
   the class of \(\mathfrak{l}_n\) (exact) and compared to the certified ratio of the sums;
3. **(T4)** Theorem 3 on the four chains of §4 (16 levels, up to 81), each
   multiplier certified as an integer, the \(\ell \mid n\) steps against the Heegner
   recursion with the exact \(m_\pi = h(n)/h(n/\ell)\).

Guard rails: precision set in `main()`; absolute/relative-error certification with
\(\ge\max(20,\mathrm{dps}/5)\) spare digits (observed \(\ge 77\)); no PSLQ.

## 9. Outlook

- **Kolyvagin's bound** (§5): the localization analysis of \(\kappa_{n,\ell}\) and the
  resulting bound on the \(p\)-part of \(\mathrm{Cl}(H_n)\) in terms of
  \([\mathcal{O}_{H_n}^\times : \mathcal{V}_n]\); test against class numbers of ring
  class fields of \(\mathbb{Q}(i)\) at prime conductor.
- **The full index** \([\mathcal{O}_{H_n}^\times : \mathcal{V}_n]\) at \(n = 9, 11, 13\)
  (degree-12 fields; PARI's unit group), predicted by Theorem 3 and the class number
  formula to be of Kubert–Lang shape — this would settle Paper II's Conjecture 6.7
  at those levels.
- **Porting to Paper II.** Remark 3.2, the conductor correction of Prop. 3.1, the
  column \(C_n(0)\) of Theorem 6.6 and the pullback clause of Conjecture 6.7 can be
  upgraded from certified to proved by citing Theorems 2–3; the natural place is a
  short new subsection of §6.
- **Other fields**: the fiber lemma is local and transports to any imaginary
  quadratic \(K\) once the primitive-lattice dictionary is set up cusp by cusp.
