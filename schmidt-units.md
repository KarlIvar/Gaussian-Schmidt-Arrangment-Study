# Schmidt units: the unit theorem, per-class valuations, and the Robert index

This document carries out the Schmidt-unit program of [outlook.md](outlook.md)
§2.7 and §3.4(a), starting from Theorems 2 and 5 of
[phase-kronecker-limit.md](phase-kronecker-limit.md) and the first-power
descent of [first-power-descent.md](first-power-descent.md). Results:

1. **The unit theorem** (§2): the \(\mathfrak{r}\)-twisted \(\Delta\)-ratios
   \(R_f\) are algebraic **units**, for every class of every odd level —
   proved, via a new lattice-level closed form
   \(R_f = r_0^6\,\Delta(\mathfrak{b}_1)/\Delta(\mathfrak{r}^{-1}\mathfrak{b}_1)\)
   and a local rigidity argument for isogenies attached to invertible
   ambiguous ideals. The certified palindromic integer polynomials of
   phase-kronecker-limit Theorem 5 are thereby theorems, and are extended to
   all odd \(n \le 21\) and to the imprimitive strata (where \(R \equiv 1\)
   exactly when the induced twist ideal is principal).
2. **The per-class valuation law** (§3): the Euclidean \(\Delta\)-data
   \(G_\mathfrak{c} = n^{12}\Delta(\Lambda_\mathfrak{c})/\Delta(\mathbb{Z}[i])\)
   have \(P\)-adic valuation **independent of the class**, with explicit value
   \(w_p(k)\) — the per-class refinement of the \(\Delta\)-mass law
   ([euclidean-moduli-invariants.md](euclidean-moduli-invariants.md) Thm 4),
   proved without quasi-canonical-lifting inputs and verified by Newton
   polygons at every computed level. This settles the split-prime part of
   euclidean §6.1 (the \(5^{4k}\)-ladder) exactly, and decomposes the
   recorded inert ladders (\(7^3, 7^6, 7^9, 7^{13}\) at \(n=7\); the 3-adic
   ladder at \(n=9\)) into certified single-slope ingredients.
3. **First-power Schmidt units** (§4): the canonical sixth root
   \(w_f = u_f\,(\gamma_2^2\gamma_3)(\tau_{\mathfrak{r}f})/(\gamma_2^2\gamma_3)(\tau_f)\)
   satisfies \(w_f^6 = R_f\), \(w_fw_{\mathfrak{r}f} = 1\) and
   \(\overline{w_f} = w_{f^{-1}}\) **exactly at first power**; its Galois
   cocycle lies in \(\mu_6\), and the minimal exponent \(m(n)\) with
   \(\{w_f^{m}\}\) Galois-stable is certified for all odd \(n \le 35\). At
   \(n = 3, 7, 13, 15, 25\) the descent is complete (\(m = 1\)): e.g. at
   \(n = 7\) the first-power Schmidt units are
   \(w = -\varepsilon_{12}^{\pm1}\), the fundamental unit of
   \(\mathbb{Q}(\sqrt3)\).
4. **The Robert index** (§5): the quadratic-layer projections of the unit
   systems are \(-24\,m_\chi\log\varepsilon_{d_2}\) (hyperbolic) and
   \(-12\,m_\chi\log\varepsilon_{d_2}\) (Euclidean) with
   \(m_\chi = \tfrac{2h(d_1)}{w(d_1)}h(d_2)C(0)\) — class numbers appear as
   indices; and on the cubic layer (all eight Euclidean levels with a cubic
   character, \(n = 9, 11, 13, 18, 22, 23, 26, 27\)) the \(\Delta\)-coset
   unit \(\theta_u\) satisfies
   $$
   \bigl[\mathcal{O}_{L_3}^\times : \langle -1, \theta_u\rangle\bigr]
   \;=\; 8\,h_{L_3}\,C_n(0),
   $$
   \(C_n(0) = 1\) at primitive \(\chi_3\): with \(h_{L_3} = 1, 1, 3, 2\) at
   \(n = 9, 11, 13, 23\) — **the class number of the cubic field is the
   Robert index**, the exact analogue of "cyclotomic units have index
   \(h^+\)", confirmed out-of-sample at \(n = 23\). The fundamental units
   are computed and certified (\(x^3+15x^2+57x-1\), \(x^3-25x^2+201x-1\),
   \(x^3-x^2+9x-1\), \(x^3-49x^2+601x-1\)).

Everything displayed is machine-verified by
[scripts/schmidt_units.py](scripts/schmidt_units.py) under the certification
guard rails of [CLAUDE.md](CLAUDE.md); §6 states exactly what is proved, what
is certified, and what failed. Notation as in
[first-power-descent.md](first-power-descent.md): \(n \ge 3\) odd,
\(D = 1-n^2\), \(r_0 = \tfrac{n-1}2\), \(s_0 = \tfrac{n+1}2\), \(\mathcal{O}\)
the order of discriminant \(D\), \(\mathfrak{r} = [r_0, \omega_0]\)
(\(\omega_0 = \sqrt D/2\)) the invertible ambiguous twist ideal with
\(\bar{\mathfrak{r}} = \mathfrak{r}\), \(\mathfrak{r}^2 = (r_0)\),
\(N\mathfrak{r} = r_0\); per primitive class \(f\), \(\mathfrak{b}_1 = [1, m_1]\)
with \(m_1\) the CM point of the reduced form, \(\beta_1 = j(f)\),
\(\beta_2 = j(\mathfrak{r}f)\), and
\(R_f = u_f^6\,\beta_2^4(\beta_2-1728)^3/\bigl(\beta_1^4(\beta_1-1728)^3\bigr)\).

## 1. The lattice form of the twisted \(\Delta\)-ratio

**Lemma 1.1 (lattice lemma).** For every primitive class \(f\) of every odd
level \(n \ge 3\),
$$
R_f \;=\; r_0^{\,6}\;\frac{\Delta(\mathfrak{b}_1)}{\Delta(\mathfrak{r}^{-1}\mathfrak{b}_1)} .
$$
In particular the \(j\)-dressing of \(R_f\) is pure bookkeeping: the twisted
ratio **is** the Siegel-type \(\Delta\)-quotient along the invertible
ambiguous twist, normalized by its norm.

*Proof.* From the proven closed form
\(u_f = -\varepsilon\mu^{-2}h_2(\mathfrak{b}_1)/h_2(\mathfrak{b}_2)\) and the
weight algebra \(h_2^6 = c\,j^4(j-1728)^3\Delta\)
([moduli-invariants.md](moduli-invariants.md) §5.5),
$$
R_f = \varepsilon^6\mu^{-12}\,\frac{\Delta(\mathfrak{b}_1)}{\Delta(\mathfrak{b}_2)} .
$$
Write \(\mathfrak{b}_2 = \nu\,\mathfrak{r}^{-1}\mathfrak{b}_1\)
([first-power-descent.md](first-power-descent.md) §3); homogeneity gives
\(\Delta(\mathfrak{b}_2) = \nu^{-12}\Delta(\mathfrak{r}^{-1}\mathfrak{b}_1)\),
so \(R_f = \varepsilon^6(\nu^2/\mu^2)^6\,\Delta(\mathfrak{b}_1)/\Delta(\mathfrak{r}^{-1}\mathfrak{b}_1)\).
By Theorem 3.3 there (\(\omega_f = \varepsilon\nu^2/(r_0\mu^2) = 1\)),
\(\varepsilon\nu^2/\mu^2 = r_0\), whence \(\varepsilon^6(\nu^2/\mu^2)^6 = r_0^6\). \(\square\)

(Verified classwise at every odd \(n \le 21\), max relative deviation
\(6.5\cdot10^{-419}\) at 420 digits; the twist lattice
\(\mathfrak{r}^{-1}\mathfrak{b}_1\) is computed by exact HNF arithmetic in
\(K\).)

At \(n = 3\): \(r_0 = 1\), \(\mathfrak{r} = \mathcal{O}\), and \(R = 1\)
exactly — the mandatory anchor.

## 2. The unit theorem

> **Theorem 1 (unit theorem).** For every odd \(n \ge 3\) and every class
> \(f\) of discriminant \(1-n^2\) — primitive or not, with \(R_f\) defined for
> imprimitive classes by the lattice formula of Lemma 1.1 —
> \(v_\mathfrak{P}(R_f) = 0\) at every finite place \(\mathfrak{P}\) of
> \(\bar{\mathbb{Q}}\): \(R_f\) is an **algebraic unit**. Consequently the
> level polynomial \(\prod_f(x - R_f)\) (product over primitive classes) is a
> **monic integer polynomial with constant term \(\pm1\)**, palindromic by law
> 2 — the certified observations of phase-kronecker-limit Theorem 5 are
> theorems, for every \(n\).

The proof runs through five lemmas. Fix a prime \(p\) and a place
\(\mathfrak{P} \mid p\) of \(\bar{\mathbb{Q}}\), with \(v\) the valuation
normalized by \(v(p) = 1\).

**Lemma 2.1 (tangent-scalar expression).** Let \(\mathfrak{a}\) be a lattice
with CM by \(\mathcal{O}\), \(\mathfrak{c}\) an invertible
\(\mathcal{O}\)-ideal, and \(\varphi: E_{\mathfrak{a}} \to
E_{\mathfrak{c}^{-1}\mathfrak{a}}\) the isogeny that is \(z \mapsto z\)
analytically. Then
$$
v\!\left(\frac{\Delta(\mathfrak{a})}{\Delta(\mathfrak{c}^{-1}\mathfrak{a})}\right)
\;=\; -12\,v(\delta_\varphi),
$$
where \(\delta_\varphi\) is the scalar by which \(\varphi\) pulls back Néron
differentials of good-reduction models at \(\mathfrak{P}\).

*Proof.* Choose Weierstrass models \((E_i, \omega_i)\) over \(\bar{\mathbb{Q}}\)
whose period lattices are \(\lambda\mathfrak{a}\), \(\mu\mathfrak{c}^{-1}\mathfrak{a}\);
then \(\Delta(E_1,\omega_1) = \lambda^{-12}\Delta(\mathfrak{a})\),
\(\Delta(E_2,\omega_2) = \mu^{-12}\Delta(\mathfrak{c}^{-1}\mathfrak{a})\), and
\(\varphi^*\omega_2 = t\,\omega_1\) with \(t = \mu/\lambda \in \bar{\mathbb{Q}}\)
(the isogeny and the differentials are algebraic), so
\(\Delta(\mathfrak{a})/\Delta(\mathfrak{c}^{-1}\mathfrak{a}) =
t^{-12}\,\Delta(E_1,\omega_1)/\Delta(E_2,\omega_2)\). CM curves have
potentially good reduction; over a finite extension \(L/\mathbb{Q}_p\) choose
minimal models: their discriminants are units, and writing
\(\omega_i = u_i\,\omega_{i,\min}\) one gets
\(v(\Delta(E_i, \omega_i)) = -12v(u_i)\) and
\(\delta_\varphi = t\,u_1/u_2\); the \(u_i\) cancel and the display follows. \(\square\)

**Lemma 2.2 (multiplicativity; prime-to-\(p\) triviality).** Tangent scalars
multiply along composites, and \(v(\delta_\varphi) = 0\) whenever
\(\deg\varphi\) is prime to \(p\) (both \(\delta_\varphi\) and
\(\delta_{\hat\varphi}\) are integral — isogenies extend to Néron models — and
\(\delta_\varphi\delta_{\hat\varphi} = \deg\varphi\)).

**Lemma 2.3 (the \(p\)-part of the twist).** \(\mathfrak{r} =
\mathfrak{b}\mathfrak{c}\) with \(\mathfrak{b}, \mathfrak{c}\) invertible
\(\mathcal{O}\)-ideals, \(N\mathfrak{b} = p^{k}\) (\(k = v_p(r_0)\)),
\(N\mathfrak{c}\) prime to \(p\), and \(\mathfrak{b}\) inherits both
structural properties: \(\bar{\mathfrak{b}} = \mathfrak{b}\) and
\(\mathfrak{b}^2 = (p^{k})\).

*Proof.* Invertible fractional \(\mathcal{O}\)-ideals decompose along primes
(localization); let \(\mathfrak{b}\) be the ideal with
\(\mathfrak{b}_q = \mathcal{O}_q\) for \(q \ne p\) and
\(\mathfrak{b}_p = \mathfrak{r}_p\). Conjugation acts locally, and
\(\bar{\mathfrak{r}} = \mathfrak{r}\) gives \(\bar{\mathfrak{b}} = \mathfrak{b}\).
Locally \((\mathfrak{b}^2)_p = (\mathfrak{r}^2)_p = r_0\mathcal{O}_p =
p^k\mathcal{O}_p\) and \((\mathfrak{b}^2)_q = \mathcal{O}_q\), which is the
local data of \((p^k)\). \(\square\)

**Lemma 2.4 (rigidity).** The quantity \(d_p := v(\delta_{\varphi})\) for the
\(\mathfrak{b}\)-twist isogeny \(\varphi: E_{\mathfrak{a}} \to
E_{\mathfrak{b}^{-1}\mathfrak{a}}\) is the **same for every proper
\(\mathcal{O}\)-ideal \(\mathfrak{a}\)**.

*Proof.* \(v(\delta_\varphi)\) is computed on formal groups: it is the tangent
valuation of the quotient map \(G \to G/C\) of \(p\)-divisible groups over
\(\mathcal{O}_L\), where \(G = E_{\mathfrak{a}}[p^\infty]\) and \(C\) is the
schematic closure of \(\ker\varphi = \mathfrak{b}_p^{-1}\mathfrak{a}_p/\mathfrak{a}_p\)
(the quotient elliptic curve is the good model of
\(E_{\mathfrak{b}^{-1}\mathfrak{a}}\) by the Néron property, and
\(\operatorname{Lie}\) of an elliptic curve equals \(\operatorname{Lie}\) of
its \(p\)-divisible group). So \(d_p\) depends only on the isomorphism class
of the pair \((G, C)\). Given two classes, choose an invertible
\(\mathcal{O}\)-ideal \(\mathfrak{q}\) of norm prime to \(p\) connecting them;
the connecting isogeny has degree prime to \(p\), hence induces an
isomorphism of \(p\)-divisible groups, is \(\mathcal{O}_p\)-linear on Tate
modules (\(T_p E_{\mathfrak{a}} = \mathfrak{a}\otimes\mathbb{Z}_p\), free of
rank one over \(\mathcal{O}_p\) since proper ideals of quadratic orders are
invertible), and identifies the two kernels
\(\mathfrak{b}_p^{-1}\mathfrak{a}_p/\mathfrak{a}_p\). Homothety adjusts
representatives within a class. \(\square\)

**Proposition 2.5 (balance).** \(2\,d_p = k\).

*Proof.* Compose the \(\mathfrak{b}\)-twist at \(\mathfrak{a}\) with the
\(\mathfrak{b}\)-twist at \(\mathfrak{b}^{-1}\mathfrak{a}\):
since \(\mathfrak{b}^2 = (p^k)\), \(\mathfrak{b}^{-2}\mathfrak{a} =
p^{-k}\mathfrak{a}\), and the composite followed by the isomorphism
\(z \mapsto p^k z\) is the endomorphism \([p^k]\) of \(E_{\mathfrak{a}}\), of
tangent scalar \(p^{k}\). By Lemma 2.2 the tangent valuations add to \(k\);
by Lemma 2.4 both steps contribute the same \(d_p\) (the twist of a proper
ideal is proper). \(\square\)

*Proof of Theorem 1.* By Lemma 1.1 and Lemmas 2.1–2.3,
\(v(R_f) = 6v(r_0) - 12\,v(\delta_{\varphi_{\mathfrak{b}}}) = 6k - 12\,d_p
= 6k - 6k = 0\) at every \(\mathfrak{P} \mid p\) with \(p \mid r_0\), and
\(v(R_f) = 0\) trivially at \(p \nmid r_0\) (prime-to-\(p\) twist, Lemma 2.2).
For an imprimitive class of content \(g\), \(\mathfrak{b}_1\) is proper over
the order \(\mathcal{O}'\) of discriminant \(D/g^2\) and the twist is by
\(\mathfrak{r}' := \mathfrak{r}\mathcal{O}'\) — invertible (extension of an
invertible ideal), ambiguous, with \(\mathfrak{r}'^2 = r_0\mathcal{O}'\) — so
the identical argument runs over \(\mathcal{O}'\). Finally the polynomial
statement: rationality of \(\prod(x - R_f)\) is phase-kronecker-limit Theorem
5 (from the first-power descent); the roots being units makes the
coefficients rational algebraic integers, i.e. integers, with
\(\prod_f R_f = 1\) for \(n \ge 5\) (law 2 pairs \(f\) with \(\mathfrak{r}f\))
and \(R = 1\) at \(n = 3\); palindromy is law 2. \(\blacksquare\)

**Remark (what the proof did *not* need).** No Kronecker limit formula, no
Siegel units, and — despite the Deuring dictionary that motivated it — no
quasi-canonical-lifting valuation formulas: Gross's theory enters only as the
conceptual reading (the twist is a *horizontal* move, and horizontal moves
preserve local \(\Delta\)-data). The inputs are potentially good reduction,
Néron differentials, and the rigidity of \(p\)-divisible-group quotients.

**The extended certified record.** All level polynomials
\(\prod_{f\ \mathrm{prim}}(x - R_f)\), certified integer, palindromic,
constant term \(1\), and **irreducible over \(\mathbb{Q}\)** (single Galois
orbit; exact factorization), \(\ge 237\) spare digits:

| \(n\) | \(h\) | \(\prod_f(x - R_f)\) |
|---|---|---|
| 3 | 1 | \(x - 1\) |
| 5 | 2 | \(x^2 - 34x + 1\) |
| 7 | 2 | \(x^2 - 2702x + 1\) |
| 9 | 4 | \(x^4 - 339524x^3 - 95354x^2 - 339524x + 1\) |
| 11 | 4 | \(x^4 - 56529284x^3 + 1538876166x^2 - 56529284x + 1\) |
| 13 | 4 | \(x^4 - 11382984004x^3 + 885435408006x^2 - 11382984004x + 1\) |
| 15 | 8 | \(x^8 - 2628641876392x^7 - 21595933374628x^6 - 1373071731101336x^5 + 9740462908109254x^4 - (\text{sym})\) |
| 17 | 4 | \(x^4 - 673122277718404x^3 + 8553847041196806x^2 - 673122277718404x + 1\) |
| 19 | 8 | \(x^8 - 186863535844922888x^7 + 44665915402536486036508x^6 + 131633377547326082495944x^5 + 179879784238619113420870x^4 + (\text{sym})\) |
| 21 | 12 | \(x^{12} - 55348592922774901452x^{11} + 87282798056992201360611266x^{10} + \cdots\) (script prints all coefficients) |

**The imprimitive strata.** With \(R\) defined by the lattice formula, every
stratum polynomial is again integer, palindromic, with constant term
\(\pm1\) — the unit property holds on all strata — and
$$
R \equiv 1 \text{ on the content-}g\text{ stratum}
\iff \mathfrak{r}\mathcal{O}' \text{ is principal in } \mathcal{O}'.
$$
The direction \(\Leftarrow\) is proved: \(\mathfrak{r}\mathcal{O}' = (\lambda)\)
forces \(\lambda^2 = r_0\zeta\) with \(\zeta \in \mathcal{O}'^\times\)
torsion, so \(\Delta(\mathfrak{r}^{-1}\mathfrak{b}_1) =
\lambda^{12}\Delta(\mathfrak{b}_1) = r_0^6\Delta(\mathfrak{b}_1)\) and
\(R = 1\). The direction \(\Rightarrow\) is observed (all computed strata).
Data:

| \(n\) | content | stratum disc | \(\mathfrak{r}\mathcal{O}'\) | stratum polynomial |
|---|---|---|---|---|
| 7 | 2 | \(-12\) | \((\sqrt{-3}\,)\) principal | \(x - 1\) |
| 7 | 4 | \(-3\) | principal | \(x - 1\) |
| 9 | 2 | \(-20\) | \((2)\) principal | \((x-1)^2\) |
| 15 | 2 | \(-56\) | \([7, \sqrt{-14}]\) **non**principal | \(x^4 - 1988x^3 - 3194x^2 - 1988x + 1\) |
| 17 | 2 | \(-72\) | \([4, 6\sqrt{-2}]\) **non**principal | \(x^2 - 9602x + 1\) |
| 17 | 3 | \(-32\) | \((2\sqrt{-2}\,)\) principal | \((x-1)^2\) |
| 17 | 6 | \(-8\) | \((2\sqrt{-2}\,)\) principal | \(x - 1\) |
| 19 | 3 | \(-40\) | \((3)\) principal | \((x-1)^2\) |

So "what breaks on the imprimitive strata" is: **nothing** — the unit theorem
persists; what changes is that the twist can *die* (principal induced ideal),
collapsing the stratum to \(R \equiv 1\).

## 3. The per-class valuation law (Euclidean corollary)

Recall \(G_\mathfrak{c} = n^{12}\Delta(\Lambda_\mathfrak{c})/\Delta(\mathbb{Z}[i])\)
over \(\mathrm{Cl}(\mathcal{O}_n)\), \(\mathcal{O}_n = \mathbb{Z}+n\mathbb{Z}[i]\),
with \(D_n(x) = \prod_\mathfrak{c}(x - G_\mathfrak{c}) \in \mathbb{Z}[x]\) and
\(|D_n(0)| = |M(n)|\) the \(\Delta\)-mass
(phase-kronecker-limit Theorem 2; euclidean Theorem 4).

> **Theorem 2 (per-class valuations).** For every level \(n \ge 2\), every
> prime \(p\) with \(p^k \| n\), and every place \(\mathfrak{P} \mid p\)
> (normalized \(v(p) = 1\)), the valuation \(v_\mathfrak{P}(G_\mathfrak{c})\)
> is **independent of the class \(\mathfrak{c}\)**, equal to
> $$
> w_p(k) \;=\;
> \begin{cases}
> 0, & p \text{ split } (p \equiv 1 \bmod 4),\\[2pt]
> \dfrac{12\,(p^k-1)}{(p-1)\,p^{k-1}(p+1)}, & p \text{ inert } (p \equiv 3 \bmod 4),\\[6pt]
> 6\,(2^k-1)/2^k, & p = 2;
> \end{cases}
> $$
> and \(v_\mathfrak{P}(G_\mathfrak{c}) = 0\) at \(p \nmid n\).

*Proof.* Write \(\psi_\mathfrak{c}: E_{\Lambda_\mathfrak{c}} \to E_0 :=
E_{\mathbb{Z}[i]}\) for the isogeny \(z \mapsto z\) (degree \(n\)); by Lemma
2.1, \(v(G_\mathfrak{c}) = 12k - 12\,v(\delta)\) where \(\delta\) is the
tangent scalar of the \(p\)-primary part of \(\psi_\mathfrak{c}\) (Lemma 2.2).
By duality it suffices to control \(\hat\psi\)-data on the **fixed** curve
\(E_0\): the \(p\)-primary kernel of the dual is the subgroup
\(C_\mathfrak{c} = p^{-k}\Lambda_{\mathfrak{c},p}/\mathbb{Z}_p[i]
\subset E_0[p^\infty]\). Local structure: \(\Lambda_{\mathfrak{c},p}\) is an
invertible \(\mathcal{O}_{n,p}\)-module inside \(\mathbb{Z}_p[i]\) of index
\(p^k = [\mathbb{Z}_p[i] : \mathcal{O}_{n,p}]\), hence
\(\Lambda_{\mathfrak{c},p} = \epsilon\,\mathcal{O}_{n,p}\) with
\(\epsilon \in K_p^\times\); the containment and index force
\(\epsilon \in \mathbb{Z}_p[i]^\times\). Since
\(\mathbb{Z}_p[i] \subset \operatorname{End}(E_0[p^\infty])\) (the closure of
\(\mathbb{Z}[i]\) in the \(p\)-adically complete endomorphism ring),
\(\epsilon\) is an *automorphism of the \(p\)-divisible group* carrying
\(C_{\mathfrak{c}'}\) to \(C_{\mathfrak{c}}\): the pairs
\((E_0[p^\infty], C_\mathfrak{c})\) are isomorphic for all classes, so
\(v(\delta_{\hat\psi_\mathfrak{c}})\), hence \(v(G_\mathfrak{c})\), is
class-independent. The value is then pinned by the \(\Delta\)-mass law:
\(h\cdot w = v_p\bigl(\prod_\mathfrak{c} G_\mathfrak{c}\bigr) = v_p(M(n))\),
and \(v_p(M(n))/h = w_p(k)\) by euclidean Theorem 4 together with
\(h = \tfrac12N_e(p^k)N_e(n/p^k)\). \(\blacksquare\)

**Verification (exact).** The multiset of \(\mathfrak{P}\)-adic root
valuations of \(D_n\) is its Newton polygon at \(p\). For every computed
level — \(n \le 18\), \(21\), \(25\), \(27\), \(49\) — and every
\(p \mid n\), the polygon has a **single slope of the predicted value**
across the full degree \(h\):

| \(n\) | polygon slopes | \(n\) | polygon slopes |
|---|---|---|---|
| 2 | \(2{:}\,3\) | 12 | \(2{:}\,9/2;\ 3{:}\,3\) |
| 4 | \(2{:}\,9/2\) | 14 | \(2{:}\,3;\ 7{:}\,3/2\) |
| 7 | \(7{:}\,3/2\) | 16 | \(2{:}\,45/8\) |
| 8 | \(2{:}\,21/4\) | 18 | \(2{:}\,3;\ 3^2{:}\,4\) |
| 9 | \(3^2{:}\,4\) | 27 | \(3^3{:}\,13/3\) |
| 11 | \(11{:}\,1\) | 49 | \(7^2{:}\,12/7\) |

(split primes: slope 0 throughout, e.g. \(n = 5, 13, 25\)). The denominators
of the slopes (\(p^{k-1}(p+1)\) inert, \(2^k\) ramified) are exactly the
ramification degrees of Gross's quasi-canonical liftings of level \(k\) — the
per-class law *is* the quasi-canonical valuation ladder, obtained here from
rigidity plus the elementary mass recursion.

**Corollary 3.1 (the split ladder of euclidean §6.1).** Let \(p \equiv 1
\bmod 4\), \(p \nmid 6\), \(p^k \| n\), \(n \ne p^k\). Then for **every**
class,
\(v_\mathfrak{P}(u_\mathfrak{c}^2) = -4k\).
*Proof.* \(u^2 = -12\beta(\beta-1728)\,g_2(\Lambda_\mathfrak{c})/g_2(\mathbb{Z}[i])\)
(euclidean §5.6) and \((g_2\text{-quot})^3 = (\beta/1728)\,G_\mathfrak{c}/n^{12}\),
so \(v(u^2) = \tfrac43v(\beta) + v(\beta - 1728) + \tfrac13v(G) - 4k\).
Theorem 2 gives \(v(G) = 0\). At an ordinary place, \(\beta \equiv 0\) would
force the reduction's endomorphism ring to contain orders of two distinct
imaginary quadratic fields (impossible), and \(\beta \equiv 1728\) would
force equal prime-to-\(p\) conductors \(n/p^k = 1\) (Deuring), excluded.
\(\square\)
This proves the exactly arithmetic \(5^{4k}\)-denominator ladder observed at
\(n = 15\) (with the witness \(v_5(H_{-900}(0)) = v_5(H_{-900}(1728)) = 0\)
certified). At non-split \(p\) the same reduction expresses the
\(\lambda_n\)-fine-structure through \(w_p(k)\) plus the \(\beta\)-side
collision valuations — the remaining open ingredient is exactly the
conductor-degenerate Gross–Zagier data of outlook §2.8, no longer the
\(\Delta\)-part.

**Corollary 3.2 (the inert ladders of euclidean §5.5, decomposed).** At
\(n = 7\) (\(p = 7\), \(k = 1\)) and \(n = 9\) (\(p = 3\), \(k = 2\)) every
ingredient of \(v_P(u^2_\mathfrak{c}) = \tfrac43v(\beta) + v(\beta-1728) +
\tfrac13v(G) - 4k\) is certified to be a *single* Newton-polygon slope —
the \(\beta\)-collision valuations are class-independent too — giving the
constant per-class values
$$
n=7:\ v_7(u^2_\mathfrak{c}) = \tfrac43\cdot0 + \tfrac14 + \tfrac13\cdot\tfrac32 - 4
= -\tfrac{13}4,
\qquad
n=9:\ v_3(u^2_\mathfrak{c}) = \tfrac43\cdot\tfrac12 + \tfrac12 + \tfrac13\cdot4 - 8
= -\tfrac{11}2 .
$$
The recorded denominator ladders follow: at \(n = 7\) the integrality bound
\(7^{\lfloor 13k/4\rfloor} = 7^3, 7^6, 7^9, 7^{13}\) is met **exactly** by
the data of euclidean §5.5; at \(n = 9\) the bound
\(3^{\lfloor 11k/2\rfloor}\) is attained at the coset-structured
\(k = 3, 6\) (\(3^{16}, 3^{33}\)), with structured cancellation at the
other \(k\) (\(3^2, 3^6, 3^{19}, 3^{22}\)). So the \(\lambda_n\)-question is
reduced, at every computed level, to one datum: the per-class
\(\beta\)-collision slopes (certified constant here; their general law is
outlook §2.8).

## 4. First-power Schmidt units

Set \(\gamma_2 = E_4/\eta^8\), \(\gamma_3 = E_6/\eta^{12}\) (canonical
holomorphic Weber functions, \(\gamma_2^3 = j\), \(\gamma_3^2 = j - 1728\))
and define, per primitive class with reduced CM points \(\tau_f\),
$$
w_f \;:=\; u_f\cdot
\frac{(\gamma_2^2\gamma_3)(\tau_{\mathfrak{r}f})}{(\gamma_2^2\gamma_3)(\tau_f)} .
$$

**Proposition 4.1 (exact first-power laws).** For every odd \(n\) and every
primitive class:
1. \(w_f^{\,6} = R_f\) (from \((\gamma_2^2\gamma_3)^6 = j^4(j-1728)^3\));
2. \(w_f\,w_{\mathfrak{r}f} = 1\) (law 2 plus telescoping of the
   \(\gamma\)-ratio along the 2-torsion twist);
3. \(\overline{w_f} = w_{f^{-1}}\) whenever neither reduced representative is
   a boundary form (\(b = a\) or \(a = c\)): there
   \(\tau_{f^{-1}} = -\bar\tau_f\) exactly and the real \(q\)-expansions give
   the mirror law;
4. \(\sigma(w_f) = \zeta(\sigma, f)\,w_{f^{e(\sigma)}\mathfrak{c}(\sigma)}\)
   with \(\zeta(\sigma,f)^6 = 1\), for every
   \(\sigma \in \mathrm{Gal}(\bar{\mathbb{Q}}/\mathbb{Q})\) (sixth powers are
   strictly equivariant by Theorem 1 + first-power descent).

Moreover \(w_f\) has the **Siegel-unit shape** predicted by outlook §2.2:
with \(\tau_0\) an adapted uniformizer of the pair
\((\mathfrak{b}_1, \mathfrak{r}^{-1}\mathfrak{b}_1)\)
(first-power-descent Lemma 1.1), \(j' = -2\pi i\,\gamma_2^2\gamma_3\,\eta^4\)
turns \(u_f = -\tfrac1{r_0}\,j'(\tau_0)/j'(r_0\tau_0)\) into
$$
w \;=\; -\,\frac1{r_0}\Bigl(\frac{\eta(\tau_0)}{\eta(r_0\tau_0)}\Bigr)^{\!4}
$$
— a weight-0 \(\eta\)-quotient on \(\Gamma_0(r_0)\) evaluated at a Heegner
point, well-defined up to the \(\mu_6\)-valued \(\eta^4\)-multiplier of the
basis choice. The cocycle \(\zeta(\sigma,f)\) *is* this multiplier system.

**The coherence table.** Let \(m(n) \mid 6\) be minimal with
\(\{w_f^{m}\}\) Galois-stable (equivalently: \(\prod_f(x - w_f^m)\)
certified integer). For all odd \(3 \le n \le 35\):

| \(n\) | 3 | 5 | 7 | 9 | 11 | 13 | 15 | 17 | 19 | 21 | 23 | 25 | 27 | 29 | 31 | 33 | 35 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| \(m(n)\) | 1 | 3 | 1 | 3 | 6 | 1 | 1 | 3 | 2 | 3 | 6 | 1 | 6 | 6 | 6 | 6 | 6 |

At the fully coherent levels the **first-power Schmidt units** exist
unconditionally, with certified irreducible minimal polynomials:
$$
\begin{aligned}
n=3:&\quad x + 1 \qquad (w = -1)\\
n=7:&\quad x^2 + 4x + 1 \qquad (w = -\varepsilon_{12}^{\mp1} = -2\pm\sqrt3\,)\\
n=13:&\quad x^4 + 50x^3 + 123x^2 + 50x + 1\\
n=15:&\quad x^8 + 112x^7 - 630x^6 + 1568x^5 - 2109x^4 + 1568x^3 - 630x^2 + 112x + 1\\
n=25:&\quad x^8 + 13346x^7 + 24067x^6 - 4876x^5 + 11653x^4 - 4876x^3 + 24067x^2 + 13346x + 1
\end{aligned}
$$
The \(n = 7\) case is emblematic: the phase of a level-7 Schmidt circle,
\(\gamma\)-normalized, **is** the fundamental unit of \(\mathbb{Q}(\sqrt3)\)
up to sign. At the partially coherent levels the intermediate polynomials are
equally clean, e.g. \(w^3\): \(x^2+6x+1\) at \(n=5\)
(\(w^3 = -\varepsilon_8^{\pm2}\)), \(x^4+584x^3+766x^2+584x+1\) at \(n=9\),
\(x^4 + 25944604x^3 + 99499206x^2 + 25944604x + 1\) at \(n=17\); \(w^2\) at
\(n=19\): \(x^8 - 571772x^7 + 36149386x^6 - 42809072x^5 + 82773139x^4 -
(\text{sym})\).

**The failure pattern (certified negative results).**
- No congruence condition in \(n\) alone explains the table: the sign part
  \(m_2\) and cube part \(m_3\) vary independently (\(m_2 = 2\) exactly at
  \(n \in \{11, 19, 23, 27, 29, 31, 33, 35\}\), \(m_3 = 3\) exactly at
  \(n \in \{5, 9, 11, 17, 21, 23, 27, \dots, 35\}\) in the computed range),
  and every candidate congruence law formulated on the data through \(n = 27\)
  was falsified at \(29 \le n \le 35\).
- At boundary reduced forms the mirror law genuinely fails, by the explicit
  factor \((\gamma_2^2\gamma_3)(\tau+1) = -\zeta_3\,(\gamma_2^2\gamma_3)(\tau)\)
  resp. \((\gamma_2^2\gamma_3)(-1/\tau) = -(\gamma_2^2\gamma_3)(\tau)\)
  (observed as non-real \(w\) on ambiguous boundary classes at \(n = 31\)).
- Re-evaluating the \(\gamma\)'s on a Schertz 6-system (representatives with
  \(\gcd(a, 6) = 1\), \(12 \mid b\)) does **not** restore first-power
  coherence (it can even break reality, and at \(n = 9\) produces
  \(\mathbb{Z}[\sqrt5]\)-coefficients instead); nor does the
  \(\eta\)-quotient form with HNF/SNF-adapted \(\tau_0\) (the
  \(\Gamma_0(r_0)\)-multiplier moves the failure elsewhere: \(m\) drops to
  \(1\) at \(n = 15\) but rises at \(n = 7\)). The complete determination of
  \(\zeta(\sigma, f)\) is the Kubert–Lang/Schertz multiplier bookkeeping for
  this *hybrid* object (a class function \(u_f\) times a ratio of level-6
  class invariants), and remains open — the certified \(m(n)\) table is the
  constraint any such computation must reproduce.

## 5. The Robert index

### 5.1 Quadratic layer

For an odd real (genus) character \(\chi\) of \(\mathrm{Cl}(1-n^2)\) with
proved factorization \(L(s,\chi) = L(s,(\tfrac{d_1}\cdot))L(s,(\tfrac{d_2}\cdot))C(s)\)
(phase-kronecker-limit §6), Theorem 5 there gives
$$
\sum_f \chi(f)\,\log|R_f| \;=\; -24\,L'(0,\chi)
\;=\; -24\,\frac{2h(d_1)}{w(d_1)}\,h(d_2)\,C(0)\,\log\varepsilon_{d_2},
$$
and analogously \(\sum_\mathfrak{c}\chi(\mathfrak{c})\log|G_\mathfrak{c}| =
-12\,(\cdots)\log\varepsilon_{d_2}\) on the Euclidean side: the
\(\chi\)-eigenprojection of the Schmidt-unit system against the fundamental
unit of its real quadratic field is an explicit integer built from the two
**class numbers** \(h(d_1), h(d_2)\). The certified values (residuals
\(\le 10^{-249}\)):

| aspect | \(n\) | \(\chi\) | projection |
|---|---|---|---|
| hyp | 11 | \((0,1)\) | \(-48\log\varepsilon_8\) |
| hyp | 11 | \((1,0)\) | \(-16\log\varepsilon_{40}\) |
| hyp | 13 | both | \(-24\log\varepsilon_{24}\), \(-24\log\varepsilon_{21}\) |
| hyp | 15 | \((0,1)\), \((2,1)\) | \(-24\log\varepsilon_{28}\), \(-12\log\varepsilon_{56}\) |
| euc | 3, 5, 7 | \(\chi_2\) | \(-4\log\varepsilon_{12}\), \(-24\log\varepsilon_5\), \(-12\log\varepsilon_{28}\) |
| euc | 9, 11, 13 | \(\chi_2\) | \(-16\log\varepsilon_{12}\), \(-12\log\varepsilon_{44}\), \(-24\log\varepsilon_{13}\) |

### 5.2 Cubic layer: the index is \(8\,h_{L_3}\)

At a Euclidean level whose class group has 3-torsion
(\(n = 9, 11, 13, 18, 22, 23, 26, 27\) among \(n \le 31\)) let \(\chi_3\) be
a cubic character, \(L_3\) the real cubic subfield of the ring class field
\(H_n\) cut out by \(\ker\chi_3\), and
\(\theta = \prod_{\mathfrak{c}\in\ker\chi_3} G_\mathfrak{c}\) the principal
\(\Delta\)-mass coset product (the real root of the certified coset cubics of
phase-kronecker-limit §3).

**Proposition 5.1.** \(\theta_u := \theta/|M(n)|^{1/3}\) is an algebraic
**unit** in \(L_3\), and \(\log|\theta_u| = -8\,L'(0,\chi_3)\) at the real
embedding.

*Proof.* Integrality both ways: by Theorem 2 the valuation of
\(G_\mathfrak{c}\) at every place is the constant \(w_p(k)\), so
\(v(\theta) = \tfrac{h}3 w_p(k) = \tfrac13v(M(n))\) at every place over
\(p \mid n\) and \(0\) elsewhere — exactly the valuations of
\(|M(n)|^{1/3} \in \mathbb{Z}\) (the mass is a cube at every level with
3-torsion: \(3 \mid h\) forces \(3 \mid v_p(M)\) through
\(h = \tfrac12N_e(p^k)N_e(n')\); e.g. \(3^{24}, 11^6, 1, 2^{36}3^{48},
2^{36}11^{12}, 23^6, 2^{36}, 3^{78}\)). The logarithm: writing
\(a = \log|\theta|\), \(b\) for the common log-modulus of the two conjugate
coset products, the uniform Stark law
\(-12L'(0,\chi_3) = \sum\chi_3\log|G|\)
(phase-kronecker-limit Thm 2, proved) reads \(a - b = -12L'\), while
\(a + 2b = \log|M(n)|\); solving,
\(\log|\theta_u| = a - \tfrac13\log|M| = -8L'(0,\chi_3)\). \(\square\)

**Proposition 5.2.** \(L_3\) has one real and one complex place (unit rank
1), and \(h_{L_3}R_{L_3}\,C_n(0) = L'(0,\chi_3)\), where \(L(s,\chi_3)\) is
the Epstein sum at level \(n\) and \(C_n(s)\) is the finite Euler product
relating it to the primitive \(L\)-function of \(\chi_3\)
(\(C_n \equiv 1\) when \(\chi_3\) is primitive — at \(n = 9, 11, 13, 23\),
since the class groups of the proper divisor levels have no 3-torsion except
along the displayed pullbacks).

*Proof.* \(\zeta_{L_3}(s) = \zeta(s)\,L_{\mathrm{prim}}(s,\chi_3)\)
(dihedral Artin formalism for the ring class extension), and
\(L(s,\chi_3) = L_{\mathrm{prim}}(s,\chi_3)C_n(s)\). \(L(0,\chi_3) = 0\),
\(L'(0,\chi_3) \neq 0\) and \(\zeta(0) = -\tfrac12\), so \(\zeta_{L_3}\)
vanishes to first order at \(s = 0\): \(r_1 + r_2 - 1 = 1\) with a real
embedding present, forcing \((r_1, r_2) = (1,1)\). The class number formula
at \(s = 0\), \(\zeta_{L_3}'(0) = -h_LR_L/w_L\) with \(w_L = 2\), gives
\(h_LR_L = -2\zeta'_{L_3}(0) = L'_{\mathrm{prim}}(0,\chi_3)\), and
\(L'(0,\chi_3) = C_n(0)L'_{\mathrm{prim}}(0,\chi_3)\) since
\(L_{\mathrm{prim}}(0) = 0\). \(\square\)

> **Theorem 3 (Robert index, cubic layer).** At every Euclidean level with a
> cubic character (\(n = 9, 11, 13, 18, 22, 23, 26, 27\)):
> $$
> \bigl[\mathcal{O}_{L_3}^\times : \langle -1,\ \theta_u\rangle\bigr]
> \;=\; \frac{\log|\theta_u|}{R_{L_3}} \;=\; 8\,\frac{L'(0,\chi_3)}{R_{L_3}}
> \;=\; 8\,h_{L_3}\,C_n(0),
> $$
> where \(C_n(0) = 1\) when \(\chi_3\) is primitive (its conductor is not a
> proper divisor level) and \(C_n(0)\) is the integer imprimitivity Euler
> multiplier otherwise. All quantities computed and certified:
>
> | \(n\) | \(\chi_3\) | fundamental unit of \(L_3\) | \(h_{L_3}\) | \(C_n(0)\) | index |
> |---|---|---|---|---|---|
> | 9 | primitive | \(x^3+15x^2+57x-1\) | 1 | 1 | **8** |
> | 11 | primitive | \(x^3-25x^2+201x-1\) | 1 | 1 | **8** |
> | 13 | primitive | \(x^3-x^2+9x-1\) | 3 | 1 | **24** |
> | 23 | primitive | \(x^3-49x^2+601x-1\) | 2 | 1 | **16** |
> | 18 | from \(n=9\) | \(x^3+15x^2+57x-1\) | 1 | 2 | 16 |
> | 22 | from \(n=11\) | \(x^3-25x^2+201x-1\) | 1 | 2 | 16 |
> | 26 | from \(n=13\) | \(x^3-x^2+9x-1\) | 3 | 2 | 48 |
> | 27 | from \(n=9\) | \(x^3+15x^2+57x-1\) | 1 | 4 | 32 |
>
> (\(\theta_u\) minimal polynomials in the script output; e.g.
> \(x^3 - 11708931x^2 + 115597311109635x - 1\) at \(n=9\),
> \(x^3 - 28994720086003708422147x^2 + \cdots x - 1\) at \(n=23\).) The
> identity index \(= 8L'/R_{L_3}\) follows from Propositions 5.1–5.2; the
> content of the computation is the fundamental unit (hence \(R_{L_3}\) and
> the certified integer \(L'(0,\chi_3)/R_{L_3}\), \(\ge 249\) spare digits).
> Fundamentality is rigorous: the descent tests every \(k\)-th root
> (\(k \le 40\)) for certified integrality of the minimal polynomial, and
> Friedman's unconditional regulator bound \(R > 0.2052\) caps any residual
> index by \(R_{L_3}/0.2052 < 40\).

Three structural readings. (i) At \(n = 13\) and \(n = 23\) the class
numbers \(h_{L_3} = 3, 2\) of the cubic fields are forced into the index —
precisely the Kubert–Lang/Robert phenomenon ("cyclotomic units have index
\(h^+\)"), realized by the Schmidt disks of curvature \(2n\). The
\(n = 23\) row is a genuine out-of-sample confirmation: a **new** cubic
field, computed after the law \(8h_{L_3}\) was formulated on \(9, 11, 13\).
(ii) The pullback levels are a strong internal consistency check: the
*identical* fundamental units and regulators recur at \(n = 9, 18, 27\) and
at \(n = 11, 22\) and \(n = 13, 26\), while the coset units \(\theta_u\)
differ — the index grows by exactly the local Euler value \(C_n(0)\)
(observed: \(2\) per added prime, \(4\) at the square-conductor step
\(9 \to 27\) — the same local Euler phenomenon as phase-kronecker-limit §8,
opens 2). (iii) The universal factor \(8 = 12/(3/2)\) is the weight
bookkeeping of the \(\Delta\)-normalization (the coset object is a
12th-power-level datum, projected with the \((3/2)\)-coefficient of the
Stark relation).

**Conjecture (Robert index for Schmidt units).** For every Euclidean level
whose cubic character is primitive, \([\mathcal{O}_{L_3}^\times : \langle -1,
\theta_u\rangle] = 8\,h_{L_3}\) (imprimitive levels: times the local Euler
multiplier \(C_n(0)\)); and more generally, for every ring class field
\(H_n\), the full group generated by the \(\{G_\mathfrak{c}\}\)
modulo its mass-normalization sits inside the unit/\(S\)-unit group of
\(H_n\) with index a product of class numbers of the subfields cut out by
the character group, up to powers of 2 and 3 — the Kubert–Lang ch. 12–13
shape. (Certified at all eight computed cubic levels; the quadratic layer's
class-number factors are §5.1.)

**Remark (hyperbolic side).** The hyperbolic \(R_f\) satisfy
\(R_{\mathfrak{r}f} = 1/R_f\), so their projections vanish on every
character with \(\chi(\mathfrak{r}) = +1\); since cubic characters are even
on the 2-torsion class \(\mathfrak{r}\), **the hyperbolic cubic layer
degenerates identically** — the hyperbolic Robert index lives entirely on
the odd characters (quadratic layer of §5.1, and the sextic
\(\chi_2\chi_3\)-layer at levels like \(n = 21\), where the odd-part coset
products of the \(R_f\) are units by Theorem 1; their index theory is the
natural continuation).

## 6. What is proved, what is certified, what failed

**Proved (all levels).**
- Lemma 1.1 (lattice form of \(R_f\)); Theorem 1 (unit theorem, primitive
  and imprimitive classes; monic integer palindromic level polynomials with
  constant term \(\pm1\)); the collapse \(R \equiv 1\) on strata with
  principal induced twist.
- Theorem 2 (per-class valuation law) and Corollary 3.1 (split ladder,
  \(n \neq p^k\)).
- Proposition 4.1 (first-power laws for \(w_f\); \(\mu_6\)-boundedness of
  the cocycle); the \(\eta\)-quotient closed form.
- Propositions 5.1–5.2 and the index identity of Theorem 3 (given the
  computed \(R_{L_3}\)).

**Certified (exact integer arithmetic or \(\ge 237\) spare digits at
250–420 digits), not proved in general.**
- The displayed \(R\)-polynomials (\(n \le 21\)) and stratum polynomials
  match Theorem 1's guarantees and are additionally **irreducible** (single
  Galois orbit) — irreducibility for general \(n\) remains open exactly as
  in moduli-invariants §5.10 (the distinctness subgroup \(T = 1\)).
- The Newton-polygon verification of Theorem 2 at \(n \le 18, 21, 25, 27,
  49\).
- The \(m(n)\) coherence table and the coherent-level \(w\)-polynomials of
  §4; the fundamental units, \(h_{L_3}\)-values, Euler multipliers and
  indices of Theorem 3 at all eight cubic levels (fundamentality rigorous
  via Friedman's bound; \(n = 23\) out-of-sample).
- \(R \equiv 1 \Rightarrow\) principal induced twist (converse direction of
  the stratum criterion).

**Failed / open.**
- No congruence law in \(n\) for \(m(n)\) survived the data through
  \(n = 35\); the Galois cocycle of \(w_f\) awaits the Kubert–Lang
  multiplier computation (outlook §2.2 route), for which the certified
  table is the target.
- The non-split part of the \(\lambda_n\)-law still needs the
  conductor-degenerate GZ valuations of the \(\beta\)-dressing (outlook
  §2.8); the \(\Delta\)-part is now closed.
- The general Robert-index conjecture beyond the computed levels, and the
  hyperbolic sextic layer (\(n = 21\)).

## 7. Machine verification

All displayed statements: `python3 scripts/schmidt_units.py --selftest`
(phases 1–4; ~4 minutes). Per phase the script asserts:

1. the lattice lemma classwise (max deviation \(< 10^{-3\,\mathrm{dps}/5}\));
   integer/palindromy/constant-term certification of every \(R\)-polynomial
   (regression-checked against the recorded coefficients) and of every
   imprimitive stratum polynomial; the \(R \equiv 1\) collapse on the listed
   strata; irreducibility via exact factorization;
2. the single-slope Newton polygons of the certified \(D_n\) at every
   \(p \mid n\), all levels including \(n = 49\) at 460 digits, against
   \(w_p(k)\); the split-ladder witness at \(n = 15\); the inert-ladder
   decompositions of Corollary 3.2 (denominator ladders re-derived from
   certified slopes);
3. the laws \(w^6 = R\), \(w_fw_{\mathfrak{r}f} = 1\), the off-boundary
   mirror law, and the \(m(n)\) table (each entry re-certified);
4. the quadratic-layer projections against
   \(\log\varepsilon_{d_2}\) (proved KLF factorizations as input); at all
   eight cubic levels the unit-normalized coset cubics, the \(k\)-th-root
   descent with the Friedman cap, \(L'(0,\chi_3)/R_{L_3}\) as certified
   integers (independent Epstein evaluation of \(L'\)), the Stark relation
   \(\log|\theta_u| = -8L'\), the index \(= 8L'/R_{L_3}\), and the
   regression record of Theorem 3 (identical fundamental units across the
   pullback families).

Guard rails: precision set after imports; absolute-error certification with
\(\ge\max(20, \mathrm{dps}/5)\) spare digits; no PSLQ anywhere.

## 8. Outlook

- **The cocycle computation** (the one gap in §4): implement Shimura
  reciprocity for \(\gamma_2, \gamma_3\) (Gee–Stevenhagen /
  Schertz N-systems) for the hybrid object \(w_f\) and derive \(m(n)\);
  the certified table is a sharp target, and the \(\eta^4\)-multiplier form
  suggests the answer is a Dedekind-sum formula in \((r_0, s_0)\).
- **Hyperbolic sextic Robert index**: at \(n = 21\) (and any level with
  6-part in the class group) the odd-character coset products of the
  \(R_f\) are units in the sextic subfield; compute their index against the
  rank-2 unit lattice — the first genuinely hyperbolic index datum.
- **Genus-refined GZ** (outlook §2.8) now closes the *entire*
  \(\lambda_n\)-question: only \(\beta\)-collision valuations remain.
- **First-power Euclidean variant**: the same \(\gamma\)-normalization
  applied to the Euclidean phases (where the square level is already
  cocycle-free) should give first-power lemniscatic units and drop the
  factor 8 in Theorem 3 toward the sharp Robert normalization.
- **Imprimitive strata**: the principality criterion suggests the full
  imprimitive class formula (outlook §2.5) can be organized by
  \(\mathfrak{r}\mathcal{O}'\)-classes; the \(R\)-collapse is its shadow.
