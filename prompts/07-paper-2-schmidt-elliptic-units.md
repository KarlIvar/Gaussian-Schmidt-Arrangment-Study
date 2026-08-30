# Paper II — *Schmidt circles and elliptic units: Kronecker limit formulas and the Robert index*

Read `CLAUDE.md` first, then **Paper I** — assumed already written at
`papers/1-schmidt-circles/schmidt-circles.tex` (prompt 06); this paper cites
it as [SchmidtI] and must not restate its proofs. Primary sources:
[phase-kronecker-limit.md](../phase-kronecker-limit.md) (all of it),
[schmidt-units.md](../schmidt-units.md) (all of it),
[euclidean-moduli-invariants.md](../euclidean-moduli-invariants.md) §5.5–5.6
(for statements imported through Paper I),
[moduli-invariants.md](../moduli-invariants.md) §5.4–5.7 (background for the
closed forms). Reusable machinery: `scripts/phase_klf.py`,
`scripts/schmidt_units.py`. Follow the certification guard rails of CLAUDE.md.

## Mission

Write the second, **narrow** paper: one theme, pursued to its quantitative
end — *the phases of the Schmidt arrangement form an elliptic-unit system,
and the system's index is a class number*. Everything not on the direct line
from the classwise \(|u|\)-law to the Robert index is excluded. The paper's
three headline theorems:

1. **Kronecker limit formulas** for the character sums of \(\log|u|\), in
   both aspects, with genus-character closed forms;
2. **the unit theorem**: the \(\mathfrak{r}\)-twisted \(\Delta\)-ratios
   \(R_f = r_0^6\Delta(\mathfrak{b}_1)/\Delta(\mathfrak{r}^{-1}\mathfrak{b}_1)\)
   are algebraic units (every class, every stratum), together with the
   per-class valuation law for the Euclidean \(\Delta\)-data;
3. **the Robert index**:
   \([\mathcal{O}_{L_3}^\times : \langle-1, \theta_u\rangle] = 8\,h_{L_3}\,C_n(0)\)
   on the cubic layer — the Schmidt-arrangement analogue of "cyclotomic units
   have index \(h^+\)", with the out-of-sample confirmation at \(n = 23\).

Deliverable: `papers/2-schmidt-elliptic-units/schmidt-elliptic-units.tex`
+ `references.bib`, same conventions as Paper I (deliver source and
syntax-validate if no TeX toolchain is available).

## Contents (be exactly this picky)

**§1 Introduction.** The arc in one page: the phase \(u_f\) of [SchmidtI]
is, classwise, a twisted \(\Delta\)-quotient; its character sums compute
\(L'(0,\chi)\); the twisted ratios are genuine units; and the group they
generate sits in the unit group of the ring class tower with index a class
number. Position against Siegel–Ramachandra–Robert and Kubert–Lang; state
what is proved vs certified up front.

**§2 The classwise laws and the master identities.** The Epstein
\(L\)-functions of a level and \(L'(0,\chi) = -\tfrac1{12}\sum\chi\log g\)
(Prop 1.1 of phase-kronecker-limit, with the independent incomplete-gamma
evaluation as the verification anchor); the classwise laws (Lemma 2.1
Euclidean, Lemma 5.1 hyperbolic — the cancellation of \(\varepsilon\) and
\(\mu\) is a highlight); the two master identities
\(S(\chi) = -2L' + \tfrac23\Sigma_0 + \tfrac12\Sigma_{1728}\) and the
odd/even hyperbolic dichotomy. The trivial-character anchor via the
\(\Delta\)-mass law (statement imported from [SchmidtI]).

**§3 Genus characters in closed form.** The Sturm-proved
\(R_\chi = 2(\mathrm{conv}*\mathrm{corr})\) factorizations and the resulting
\(L'(0,\chi_2)\) closed forms (both aspects; the conductor Euler corrections
at square conductors; the "real quadratic field is \(\mathbb{Q}(\sqrt n)\)"
phenomenon); the exact \(\Sigma\)-dressing factorizations in
\(\mathbb{Q}(\sqrt d)\) with Gross–Zagier split-prime support — present the
main tables, cite `phase_klf.py --selftest` for the full record; include the
certified PSLQ **non-fits** for non-real characters as honest negatives (the
Stark-regime discussion), briefly.

**§4 The unit theorem and the per-class valuation law.** From
schmidt-units.md, in full: the lattice lemma; Theorem 1 with its five-lemma
proof (tangent scalars, the \(p\)-part of the twist, rigidity, balance) —
this is the paper's technical core, keep the proof complete; the level
polynomials (integer, palindromic, constant term 1, irreducible, \(n \le
21\)); the imprimitive strata with the principality criterion. Then Theorem
2 (per-class \(w_p(k)\)) with the Newton-polygon verification table and the
ladder corollaries (the split \(5^{4k}\) ladder proved; the inert
\(7^3,7^6,7^9,7^{13}\) decomposition) — noting the Deuring/quasi-canonical
reading as interpretation, not input.

**§5 First-power Schmidt units.** One compact section: the definition of
\(w_f\), the exact first-power laws, the \(\eta\)-quotient closed form
\(w = -\tfrac1{r_0}(\eta(\tau_0)/\eta(r_0\tau_0))^4\), the coherent levels
with their minimal polynomials (\(w = -\varepsilon_{12}^{\pm1}\) at \(n=7\)
is the emblem), and the certified \(m(n)\) table with the cocycle stated as
the open problem (the boundary-form mechanism and the failed repairs in two
sentences, not more).

**§6 The Robert index.** Quadratic layer: the eigenprojection table
\(-24m_\chi\log\varepsilon_{d_2}\) / \(-12m_\chi\log\varepsilon_{d_2}\) with
\(m_\chi\) built from \(h(d_1)h(d_2)C(0)\). Cubic layer: Propositions 5.1–5.2
of schmidt-units.md (\(\theta_u\) is a unit; \(h_LR_LC_n(0) = L'(0,\chi_3)\)),
the certified fundamental units via root-descent + Friedman's bound, and
Theorem 3 with the full eight-level table (primitive levels \(9, 11, 13, 23\)
with \(h_{L_3} = 1,1,3,2\); pullback levels with the Euler multiplier; the
identical fundamental units across pullback families as the internal
consistency check). Close with the index conjecture, stated carefully
(Kubert–Lang shape, controlled 2,3-torsion).

**§7 Machine verification + outlook.** Appendix listing every script,
statement, precision, spare digits, runtime. Outlook limited to the four
genuine continuations: the \(w_f\)-cocycle (Kubert–Lang bookkeeping), the
genus-refined Gross–Zagier exponent law (outlook §2.8), the hyperbolic
sextic-layer index, Stark recognition for order \(>2\) characters.

**Excluded — do not include:** everything already in [SchmidtI] beyond
statement-level citations (the classification, counts, class formula,
composition, moduli coordinates, first-power descent of \(u\), Euclidean
irreducibility, \(\Delta\)-mass proof); the monoid/operation layer; the
spectral layer; all spherical material (census and phase alike); the sign
law and DIT comparison (Paper I's); exact GZ denominator valuations (open);
the \(\lambda_n\) non-split fine structure beyond §4's corollaries; any
merely-observed pattern not already flagged in the sources.

## Method and integrity

- Every theorem environment keeps the repo's proved/certified label; the
  Robert-index table is a certified computation resting on two proved
  propositions — present it exactly that way. No PSLQ-derived claim appears
  without its exact-arithmetic re-verification (there are none in the
  sources; keep it so).
- Re-run `scripts/phase_klf.py --selftest` and
  `scripts/schmidt_units.py --selftest` before finalizing; if any recorded
  constant disagrees, stop and reconcile before writing.
- Cite [SchmidtI] by its actual theorem numbers — read the .tex, do not
  guess. If Paper I ported the phase-atlas material, cite its sign law where
  the ambiguous-class signs arise in §3; otherwise omit signs entirely.
- **Literature diligence**: Siegel, Ramachandra, Robert (*Unités
  elliptiques*), Kubert–Lang ch. 11–13, Schertz (*Complex Multiplication*),
  Gross (*On canonical and quasi-canonical liftings*), Gross–Zagier (*On
  singular moduli*), Friedman (regulator bound — verify the exact constant
  0.2052 and statement before citing), Stark's papers, Meyer/Vlasenko–Zagier
  (higher Kronecker limit formulas) — the index theorem's novelty claim must
  survive a genuine search for prior ring-class elliptic-unit index
  computations over \(\mathbb{Q}(i)\). Record verdicts as for Paper I.
- House bookkeeping: CLAUDE.md paper-planning + outlook.md §4/§5 updates.

## Deliverables

- `papers/2-schmidt-elliptic-units/schmidt-elliptic-units.tex` +
  `references.bib`, compiling or syntax-validated, citing [SchmidtI].
- Machine-verification appendix; both selftests re-run and logged.
- CLAUDE.md + outlook.md updates; literature-diligence summary.
- A closing statement of what is new (the two limit formulas, the unit
  theorem, the per-class law, the index \(8h_{L_3}\)) and what remains open,
  each with its repo anchor.
