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

### Item 2 — The full Robert index. Status: staged ([prompts/09-full-robert-index.md](prompts/09-full-robert-index.md)).

*Question.* What is \([\mathcal{O}_{H_n}^\times : \mathcal{V}_n]\) for the group
\(\mathcal{V}_n\) generated by all ratios \(G_\mathfrak{c}/G_{\mathfrak{c}'}\) (units by
Paper II, Thm. 4.11)?

*Why it should work.* The regulator of \(\mathcal{V}_n\) is a Dedekind group
determinant over \(\mathrm{Pic}(\mathcal{O}_n)\), equal to \(24^{h-1}\prod_{\chi\ne1}|L'(0,\chi)|\)
by Paper II, Thm 2.6; the class number formula for \(H_n\) at \(s = 0\) then gives
\([\mathcal{O}_{H_n}^\times : \mu(H_n)\mathcal{V}_n] = (4\cdot24^{h-1}/w_{H_n})\,h_{H_n}\prod_{\chi\ne1}C_\chi(0)\)
for every \(n\) — the Kubert–Lang shape of Paper II's Conjecture 6.7, with the
hand-checkable anchor \(n = 3\) (index 8 in \(\mathbb{Q}(\zeta_{12})\)). Prompt 09
states the proof to be written and the PARI verification.

*First move.* \(n = 9, 11, 13\) (degree-12 fields): the unit group from PARI
(`bnfinit`), the \(\Delta\)-units from `scripts/schmidt_euler_system.py`, the
index decomposed by characters (quadratic, cubic, sextic layers); then \(n = 23\)
(degree 24) and the hyperbolic sextic layer at \(n = 21\).

*Deliverable.* A theorem settling Paper II's Conjecture 6.7 at the computed
levels, with the powers of 2 and 3 explicit.

### Item 3 — Other imaginary quadratic fields. Status: open.

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

### Item 4 — Horizontal families: Atkin–Lehner-coupled Heegner points. Status: open.

*Question.* The unit theorem needs only an invertible ambiguous ideal
\(\mathfrak{r}\) with \(\mathfrak{r}^2\) principal in an order of an imaginary
quadratic field. For which pairs \((N, D)\) is \(N^6\Delta(\mathfrak{b})/\Delta(\mathfrak{r}^{-1}\mathfrak{b})\)
a unit, and which arrangements realize which slices?

*Why it should work.* The proof of Paper II, Thm. 4.2 is verbatim general; the
Gaussian arrangement is the slice \(D = -4N(N+1)\).

*First move.* State and prove the general theorem; evaluate the limit formulas
and the index on a few other slices (\(D = -4N(N+k)\), the odd discriminants
\(1 - n^2\) of the even levels of \(i\mathcal{S}\)); identify the arrangements or
Kleinian-group orbits that carry them.

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
2 depends on nothing new. 3 at class number one depends on nothing new either,
but is a multi-session rebuild; 3 at class number above one depends on the
cusp bookkeeping. 4 depends on Paper II §4 only. 5(a) depends on Paper II §3;
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
| `prompts/` | 01, 04, 05 (done or staged on the line), 06–07 (the papers), 08 (item 1's continuation), 09 (item 2) |
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
