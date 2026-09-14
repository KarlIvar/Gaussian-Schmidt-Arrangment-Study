# The program: Schmidt circles as an elliptic-unit and Euler-system carrier

This file is the working program of the repository. It replaces the diffuse
outlook of earlier sessions ([outlook.md](outlook.md) remains as the historical
master outlook) with one line of research and the six items that continue it.
Everything in the main tree of the repository serves this line; everything else
has been moved to [extra/](extra/README.md).

## 0. Where we stand

Two papers and one further document establish the line, each with a named
verification script and a proved / certified ledger:

- **Paper I** (`papers/1-schmidt-circles/`): the level stratification of the
  Gaussian Schmidt arrangement, the three class-number censuses, the
  involution class formula \(\hat\sigma[f] = [\mathfrak{r}_n][f]^{-1}\) with its
  canonical ambiguous twist class, the phase \(u_f = \Phi_y/\Phi_x(\beta_1,\beta_2)\)
  with its laws and first-power Galois theory, and the Euclidean phase with the
  \(\Delta\)-mass law.
- **Paper II** (`papers/2-schmidt-elliptic-units/`): the Kronecker limit formulas
  \(S(\chi) = -2L'(0,\chi) + \tfrac23\Sigma_0 + \tfrac12\Sigma_{1728}\) (Euclidean)
  and \(S(\chi) = -4L'(0,\chi) + \tfrac43\Sigma_0 + \Sigma_{1728}\) on odd characters
  (hyperbolic; \(\varepsilon\) and \(\mu\) cancel); genus characters in closed form
  with the field \(\mathbb{Q}(\sqrt n)\); the unit theorem for the twisted
  \(\Delta\)-ratios \(R_f = r_0^6\Delta(\mathfrak{b}_1)/\Delta(\mathfrak{r}^{-1}\mathfrak{b}_1)\)
  on all strata; the class-independent valuations of the Euclidean
  \(\Delta\)-data; the first-power units \(w_f\); the Robert index
  \([\mathcal{O}_{L_3}^\times : \langle -1, \theta_u\rangle] = 8h_{L_3}C_n(0)\).
- **[schmidt-euler-system.md](schmidt-euler-system.md)**: the Euclidean
  \(\Delta\)-data \(G_\mathfrak{c} = n^{12}\Delta(\Lambda_\mathfrak{c})/\Delta(\mathbb{Z}[i])\)
  are norm-compatible along the ring class tower \(H_n \subset H_{n\ell}\) with
  Heegner-shaped Euler factors and the Eisenstein eigenvalue \(\ell+1\); the
  Hecke recursion \(L'(0,\chi^{(n\ell)}) = P_\ell(\chi)L'(0,\chi)\); the Robert
  index along the tower; the first Kolyvagin step.

The line, in one sentence: *a circle of the arrangement is a Heegner point
with an Atkin–Lehner partner; its phase is the derivative of the modular
correspondence there; the absolute value of the phase is a twisted
\(\Delta\)-quotient; the character sums of \(\log|u|\) are \(L'(0,\chi)\); the
twisted quotients are units; their coset products have class-number index; and
along the conductor tower they form an Euler system of Eisenstein type.*

## 1. What carries the proofs beyond \(\mathbb{Z}[i]\)

Field-independent (they transport to any imaginary quadratic \(K\) and, in
part, beyond): the Kronecker limit formula and the independent
incomplete-gamma evaluation (Paper II, Prop. 2.2, Lemma 2.3); the rigidity of
\(p\)-divisible-group quotients and the balance argument (Paper II, §4); the
local argument of the per-class valuation law (Paper II, Thm. 4.11); the
cubic-layer index identity (Paper II, Props. 6.2–6.4); the fiber lemma and the
distribution relation behind the norm relations (schmidt-euler-system.md, §1–2).

Gaussian-specific: the \(\Delta\)-mass stratification (M2 of Paper I, Thm. 7.18)
and the "every class exactly twice" dictionary (Paper I, Lemma 7.5) use that
\(\mathbb{Z}[i]\) is a principal ideal domain with unit group \(\mu_4\); the
hyperbolic side uses the level coupling \(D = -4r_0s_0\), \(s_0 = r_0 + 1\), and the
invertible ambiguous twist \(\mathfrak{r}\) with \(\mathfrak{r}^2 = (r_0)\).

## 2. The items

### Item 1 — The Schmidt \(\Delta\)-system as an Euler system over \(\mathbb{Q}(i)\). Status: done, continuation staged.

*Question.* Do the \(\Delta\)-data of the levels \(n \mid n\ell \mid \dots\) form a
norm-compatible system, and what does the system compute?

*Done* ([schmidt-euler-system.md](schmidt-euler-system.md)): Theorem 1 (norm
relations of Heegner shape with \(a_\ell \mapsto \ell+1\), proved from Paper I's
full-mass identity and a fiber lemma; 42 pairs verified to level 125),
Theorem 2 (Hecke recursion for \(L'(0,\chi)\) under pullback: every imprimitive
Euler multiplier of Paper II — the correction \(1 + 3^{1-2s}\) at \(n = 9\), the two
\(n = 15\) \(\varepsilon\)-ratios, the column \(C_n(0) = 2,2,2,4\) — is a theorem),
Theorem 3 (the Robert index multiplies by \(P_\ell(\chi_3)\) along the tower;
multipliers \(2, 4, 7, 8, 13\) certified on four chains to level 81), Lemma 5.1
(Kolyvagin's derivative classes descend to \(H_n^\times/p\) for \(p \ge 5\)).

*Remaining* ([prompts/08-euler-system-kolyvagin.md](prompts/08-euler-system-kolyvagin.md)):
the localization of the derivative classes at the primes above \(\ell\) and the
resulting bound on the \(p\)-part of \(\mathrm{Cl}(H_n)\) in terms of the index
of the \(\Delta\)-units — the anticyclotomic Thaine–Rubin statement, with
Küçüksakallı's class numbers of ring class fields of prime conductor as test
data; the comparison with Stark's unit of \(H_n\).

*Risk.* The system is of Heegner shape, not of Rubin's \((1 - \mathrm{Fr}^{-1})\)
shape; Kolyvagin's original descent is the right machine, and the
localization step is where the anticyclotomic setting differs from the
cyclotomic one.

### Item 2 — The full Robert index. Status: **done** ([robert-index-full.md](robert-index-full.md)).

*Question.* What is \([\mathcal{O}_{H_n}^\times : \mathcal{V}_n]\) for the group
\(\mathcal{V}_n\) generated by all ratios \(G_\mathfrak{c}/G_{\mathfrak{c}'}\) (units by
Paper II, Thm. 4.11)?

*Why it should work.* Paper II settles one cubic layer; Theorem 3 of the Euler
document controls every pullback layer; the class number formula predicts a
Kubert–Lang-shaped product over the characters of \(\mathrm{Pic}(\mathcal{O}_n)\).

*First move.* \(n = 9, 11, 13\) (degree-12 fields): the unit group from PARI
(`bnfinit`), the \(\Delta\)-units from `scripts/schmidt_euler_system.py`, the
index decomposed by characters (quadratic, cubic, sextic layers); then \(n = 23\)
(degree 24) and the hyperbolic sextic layer at \(n = 21\).

*Done.* Theorem 1 of [robert-index-full.md](robert-index-full.md): for every
\(n \ge 2\), \([\mathcal{O}_{H_n}^\times:\mu(H_n)\mathcal{V}_n] = 24^{h-1}\prod_{\chi\ne1}|L'(0,\chi)|/R_{H_n}
= (4\cdot24^{h-1}/w_{H_n})\,h_{H_n}\prod_{\chi\ne1}C_\chi(0)\) — a Dedekind group determinant
against the class number formula at \(s = 0\); the layer theorem gives the index in
every subfield cut out by a subgroup (the cubic \(8h_{L_3}C_n(0)\) of Paper II and the
sextic \(576h_FC_n(0)^2\)); verified with PARI/GP at \(n = 3, \dots, 15\) (unconditional)
and \(23\) (GRH) including the exact exponent matrices; Paper II's Conjecture 6.7
replaced by the theorem (§6.3). The hyperbolic side is Theorem 3 there: the index
of the \(R_f\) in the odd units of the class field of discriminant \(1-n^2\) is
\(24^{h/2}(2^{h/2-1}/Q^-)(h_Hw_{H^+})/(h_{H^+}w_H)\prod_{\mathrm{odd}}C_\chi(0)\), exact at the
nine odd levels \(5 \le n \le 21\) (\(24^6\) at \(n = 21\); sextic layer \(2\cdot24^3\)).
*Remaining*: the \(24\)-th-root saturation conjecture (index exactly
\(h_{H_n}\prod C_\chi(0)\), certified at all eight levels; Gras-type module statement)
and the law of the \(2\)-adic invariant \(Q^-\).

### Item 3 — Other imaginary quadratic fields. Status: **done at class number one** ([other-fields.md](other-fields.md)); class number two experimental.

*Question.* Transport the whole line to Stange's arrangements \(\mathcal{S}_K\).

*Why it should work.* Everything listed as field-independent in §1; the fiber
lemma is local.

*First move.* Class number one first: \(\mathbb{Q}(\sqrt{-2})\), then
\(\mathbb{Q}(\sqrt{-3})\) (the \(\mu_6\) torsion is the cleanest test of the
phase's root-of-unity bookkeeping), then \(\mathbb{Q}(\sqrt{-7})\): redo the
classification, the level, the twist ideal and the Euclidean dictionary
\(\Lambda \subset \mathcal{O}_K\); predict from genus theory that the real field of
the Euclidean genus characters is \(\mathbb{Q}(\sqrt f)\) or
\(\mathbb{Q}(\sqrt{f|D_K|})\) and test it. Class number above one afterwards: the
cusps of the Bianchi orbifold correspond to the ideal classes, the level at one
cusp sees only \(\ker(\mathrm{Pic}(\mathcal{O}_f) \to \mathrm{Pic}(\mathcal{O}_K))\)
(Stange's theorem), the full unit system lives on the pairs \(\Lambda \subset L\)
over all cusps, and the Frobenius factors of the norm relations depend on the
splitting in \(K\).

*Deliverable.* The Kronecker limit formulas, the unit theorem and the Euler
system for each \(K\); the Robert index at the first cubic levels.

*Done* ([other-fields.md](other-fields.md), `scripts/other_fields.py`): for
\(\mathbb{Q}(\sqrt{-2}), \mathbb{Q}(\sqrt{-3}), \mathbb{Q}(\sqrt{-7}), \mathbb{Q}(\sqrt{-11})\) the whole line —
classification by \(\beta \equiv 1 \bmod \sqrt{d_K}\) with the orientation and level laws
(\(2\alpha \equiv -2s \bmod |d_K|\), half-integral levels for \(d_K \equiv 1 \bmod 4\)), the
discriminant \(4(\alpha^2-1)/d_K\), the class formula \(\hat\sigma[f] = [\mathfrak{r}_\alpha][f]^{-s}\)
with the twist ideal of norm \((2\alpha-2)/\gcd(2\alpha-2,|d_K|)\), the unit theorem, the
Euclidean dictionary \(N_e = \tfrac{w_K}{2}h\), the mass law with \(\tfrac{24}{w_Ke_p}\), the
per-class valuations, the limit formulas with the genus fields \(\mathbb{Q}(\sqrt p)\) /
\(\mathbb{Q}(\sqrt{p|d_K|})\), the norm relations with \(\pi^{12}\) and the sign \((-1)^{[\ell=2]}\),
the torsion lemma and the index \(\tfrac{w_K24^{h-1}}{w_{H_n}}h_{H_n}\prod C\) — proved, and
certified at 19 levels with `bnfcertify` (saturation \(= h_{H_n}\prod C\) throughout). At
\(\mathbb{Q}(\sqrt{-5})\), \(n = 2, 3\): the single cusp sees the kernel (Stange), the two-cusp data
\(n^{12}\Delta(\Lambda)/\Delta(\mathcal{O}_K\Lambda)\) satisfy the limit formula with a correction on the
characters pulled back from \(\mathrm{Pic}(\mathcal{O}_K)\), and the index carries \(w_K/h_K\) and a
base-character factor (experimental). *Remaining*: the cusp bookkeeping at \(h_K > 1\)
and a correction-free normalization; the non-Euclidean class-number-one fields; the
sign law of \(M(n)\); the phase \(u = \Phi_y/\Phi_x\) over \(K\).

### Item 4 — Horizontal families: Atkin–Lehner-coupled Heegner points. Status: **done** ([horizontal-families.md](horizontal-families.md)).

*Question.* The unit theorem needs only an invertible ambiguous ideal
\(\mathfrak{r}\) with \(\mathfrak{r}^2\) principal in an order of an imaginary
quadratic field. For which pairs \((N, D)\) is \(N^6\Delta(\mathfrak{b})/\Delta(\mathfrak{r}^{-1}\mathfrak{b})\)
a unit, and which arrangements realize which slices?

*Why it should work.* The proof of Paper II, Thm. 4.2 is verbatim general; the
Gaussian arrangement is the slice \(D = -4N(N+1)\).

*Progress.* The general theorem is Theorem 4 of [other-fields.md](other-fields.md)
(the proof of Paper II §4 for an arbitrary invertible ambiguous twist with principal
square), and the arrangements realize different slices of one family: the unit
polynomials agree on equal pairs \((D, [\mathfrak{r}])\) — \(D = -24\) at \(\mathbb{Q}(i)\) level 5 and
\(\mathbb{Q}(\sqrt{-2})\) level 7, \(D = -20\) at \(\mathbb{Q}(\sqrt{-3})\) level 4 and \(\mathbb{Q}(\sqrt{-7})\) level 6,
\(D = -36\) at \(\mathbb{Q}(\sqrt{-7})\) level 8 and \(\mathbb{Q}(\sqrt{-11})\) level 10.

*First move.* State and prove the general theorem; evaluate the limit formulas
and the index on a few other slices (\(D = -4N(N+k)\), the odd discriminants
\(1 - n^2\) of the even levels of \(i\mathcal{S}\)); identify the arrangements or
Kleinian-group orbits that carry them.

*Done* ([horizontal-families.md](horizontal-families.md), `scripts/horizontal_families.py`):
the slice theorem (the primitive ambiguous ideals of \(\mathcal{O}_D\) are \((r_0,0,s_0)\) and
\((r_0,r_0,\tfrac{r_0+s_0}4)\), \(2^\mu\) of them, in pairs \(\mathfrak{r}\mathfrak{s}=(\theta)\) ↔ \(\mathrm{Cl}(D)[2]\);
the unit system has rank \(h/2\)); the levels of \(\mathcal{S}_K\) at the principal cusp of every \(K\)
are \(k=t|d_K|\pm2\) resp. \(4mt\pm2\) and carry the pairs \(\{t,t|d_K|\pm4\}\) resp. \(\{t,mt\pm1\}\)
(twist norm \(t\)): a slice is realized by a class-number-one field iff its pair satisfies
\(|d_K|=(s_0\mp4)/r_0\), \(|r_0-s_0|\in\{1,2\}\) (\(\mathbb{Q}(i)\)) or \(s_0=2r_0\pm1\) (\(\mathbb{Q}(\sqrt{-2})\));
\(D=-84\) is realized with its three nontrivial classes at three fields, twelve slices with
\(|D|\le120\) at no primitive stratum; the strata of conductor \(f\) carry \((D/f^2,[\mathfrak{r}\mathcal{O}_{D/f^2}])\),
the levels of a fixed \(D'\) are the powers of the fundamental Pell solution and the stratum
twist class is \(c^j\) (odd-prime rule proved: \(p^{v_p(D')}\mid N\mathfrak{r}'\iff\eta\equiv1\bmod\mathfrak{p}\)),
so each \(K\) reaches exactly one nontrivial class per \(D'\) and \((-160,[(5,0,8)])\) none; the
even levels of \(i\mathcal{S}\) (class formula with \((n-1,n-1,\tfrac n2)\), \(x^2-7x+1\) at \(n=4\), index,
phase) close the even-level ledger item; the second-kind orbit \(\mathcal{S}_K^\perp\) is classified
(discriminants \(4-y^2|d_K|\), no twist); the genus law of the twist, Theorem 3 of the index
document and the phase \(u=\Phi_y/\Phi_x\) hold on every slice (12 resp. 16 slices certified).
*Remaining*: the \(2\)-adic stratum law, the non-principal cusps, a carrier for the
unreached slices, saturation and \(Q^-\) at \(|\mathrm{Cl}(D)[2]|=8\), a limit theory for
\(\mathcal{S}_K^\perp\) (§13 there).

### Item 5 — Other special unit systems from the same mechanism. Status: open.

(a) Stark recognition of the cubic and quartic coset objects of Paper II
(Remark 2.7, §3.4) as Stark units of the corresponding subfields.
(b) The spherical aspect ([spherical-moduli-invariants.md](spherical-moduli-invariants.md)):
the virtual partners \(m_\ell/u^2\), units twisted by the negative-Pell unit, a
system with no classical name; transport the Kronecker-limit machinery once
the spherical \(u^2\)-laws are proved.
(c) The first-power units \(w_f\) after the cocycle \(m(n)\) is settled (Paper II,
Open problem 5.4): the factor 8 in the index should drop.
(d) \(p\)-adic logarithms of the Euler system as a route to the anticyclotomic
\(p\)-adic \(L\)-function of \(\mathbb{Q}(i)\) — speculative.

### Item 6 — Beyond imaginary quadratic fields. Status: speculative.

Arithmetic Kleinian groups from quaternion algebras have no cusps, hence no
\(\Delta\) and no units, but keep Heegner points on Shimura curves: the phase
becomes a derivative of a Shimura-curve correspondence and the limit formula a
Gross–Zagier height formula. Function fields offer Drinfeld modular units,
Gekeler's \(\Delta\) and Hayes' elliptic units. Both need the level notion
rebuilt first.

## 3. Order and dependencies

1 (remaining part) and 2 reuse existing scripts and are single-session tasks;
2 depends on nothing new. 3 at class number one is done (one session, one
script); 3 at class number above one depends on the cusp bookkeeping. 4 is done (one session,
one script); its remaining parts depend on the cusp bookkeeping of 3. 5(a) depends on Paper II §3;
5(b) on the spherical proofs; 5(c) on the cocycle; 5(d) on 1. 6 is the deep end.

Recommended order: 1 (remaining), 2, 3 at class number one, 4, 5, 3 at class
number above one, 6.

## 4. Repository layout for the program

Main tree (everything here feeds the line):

| file | role in the program |
|---|---|
| [circle-classification.md](circle-classification.md), [euclidean-counting.md](euclidean-counting.md), [hyperbolic-counting.md](hyperbolic-counting.md) | the arrangement, the level, the counts \(N_e(n) = 2h(-4n^2)\) and \(3H(n^2-1)\) that index the unit systems |
| [involution.md](involution.md), [class-formula-proof.md](class-formula-proof.md) | the twist class \(\mathfrak{r}_n\) and its proof — the source of the unit theorem's twist |
| [moduli-invariants.md](moduli-invariants.md), [first-power-descent.md](first-power-descent.md) | the phase, its closed form, the Norm Lemma, \(\omega \equiv 1\) — the lattice lemma's inputs |
| [euclidean-moduli-invariants.md](euclidean-moduli-invariants.md) | the Euclidean phase and the \(\Delta\)-mass law, whose step M1 is the distribution relation of the Euler system |
| [phase-kronecker-limit.md](phase-kronecker-limit.md), [schmidt-units.md](schmidt-units.md), [schmidt-euler-system.md](schmidt-euler-system.md) | the line proper |
| [spherical-moduli-invariants.md](spherical-moduli-invariants.md) | item 5(b) |
| `papers/` | Papers I and II with their bibliographies and literature records |
| `scripts/` | the verification scripts of the above (`--selftest` each) |
| [other-fields.md](other-fields.md) | item 3: the line over \(\mathbb{Q}(\sqrt{-2}), \mathbb{Q}(\sqrt{-3}), \mathbb{Q}(\sqrt{-7}), \mathbb{Q}(\sqrt{-11})\), and \(\mathbb{Q}(\sqrt{-5})\) at \(n = 2, 3\) |
| [horizontal-families.md](horizontal-families.md) | item 4: the slices \((D,[\mathfrak{r}])\), their carriers (primitive and imprimitive strata, \(i\mathcal{S}\), class-number-two cusps, the second-kind orbit), the index and the phase on slices |
| `prompts/` | 01, 04, 05 (done or staged on the line), 06–07 (the papers), 08 (item 1's continuation), 09 (item 2, done), 10 (item 3, done at class number one), 11 (item 4, done) |
| [outlook.md](outlook.md) | the historical master outlook; its §6 summarizes this program |

Moved to `extra/` ([extra/README.md](extra/README.md)): the monoid and operation
layer, the spectral layer, the circle-language composition, the phase atlas
with the sign law and the DIT comparison, their scripts, figures and prompts.

## 5. Working rules

The house rules of [CLAUDE.md](CLAUDE.md) apply unchanged: every finite claim is
machine-verified by a named script; documents state proved / certified /
experimental; precision is set after imports; integers and rationals are
accepted only with \(\ge\max(20,\mathrm{dps}/5)\) spare digits in absolute error;
no PSLQ-derived claim without exact re-verification; new documents follow the
statement–proof style with a verification section and an outlook, and are
registered in CLAUDE.md's document map and this file.
