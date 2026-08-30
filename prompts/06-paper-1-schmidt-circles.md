# Paper I — *Schmidt circles of the Gaussian integers: classification, class numbers, and the phase invariant*

Read `CLAUDE.md` first. Primary sources, in reading order:
[circle-classification.md](../circle-classification.md),
[euclidean-counting.md](../euclidean-counting.md),
[hyperbolic-counting.md](../hyperbolic-counting.md),
[involution.md](../involution.md),
[class-formula-proof.md](../class-formula-proof.md),
[circle-composition.md](../circle-composition.md),
[spherical-moduli-invariants.md](../spherical-moduli-invariants.md) §1–2 only,
[moduli-invariants.md](../moduli-invariants.md),
[first-power-descent.md](../first-power-descent.md),
[euclidean-moduli-invariants.md](../euclidean-moduli-invariants.md).
Follow the certification guard rails of CLAUDE.md throughout.

## Mission

Write the first of the project's two papers: the **broad structural article** on
the Gaussian Schmidt arrangement. It contains every foundational and
modular-invariant result, and *only* results with complete proofs (per-level
computational verifications are stated as such, with the verifying script
named). Paper II (prompt 07) will build on it; do **not** include Paper II's
material (see the exclusion list).

Deliverable: a self-contained LaTeX manuscript in `papers/1-schmidt-circles/`
(`schmidt-circles.tex` + `references.bib`; `amsart` or `article` with `amsthm`),
written for a strong number-theory journal readership. If no TeX toolchain is
available in the environment, deliver the source and validate it by a careful
syntax pass; content must not depend on compilation.

## Contents (be exactly this picky)

**§1 Introduction.** The arrangement \(\mathcal{S} =
\mathrm{PSL}_2(\mathbb{Z}[i])\cdot\hat{\mathbb{R}}\); the three headline
structures (classification, class-number counts in three geometries, the phase
invariant); the novelty discussion (see "Sign law and DIT" below); a roadmap.
State plainly which results are theorems and which are per-level certified
computations.

**§2 Classification and the level.** The congruence classification of
\(\mathcal{S}\) (necessity + descent sufficiency), the single equation
\(u^2 + v(v+1) = nm\), the \(\mathrm{PGL}_2\)-variant
\(\mathcal{S} \sqcup i\mathcal{S}\); the level \(\alpha\) with its three
characterizations (imaginary part / \(\coth\) of hyperbolic radius; trace on
the Cartan embedding; inversive product with \(\hat{\mathbb{R}}\)).

**§3 Counting in three geometries.** (a) \(N_e(n) = n\prod_{p\mid n}(1 -
\chi_{-4}(p)/p)\) and \(\sum_{n\le X}N_e(n) = X^2/2G + O(X\log X)\) (Catalan's
constant); (b) level-\(n\) circles = hyperbolic circles of radius
\(\operatorname{arcoth} n\) at CM points of disc \(1-n^2\), incidence lemma,
weighted ideal-triangle count \(3H(n^2-1)\); (c) the spherical census: radius
quantization \(\cot\theta = \ell \in \mathbb{Z}\) and
\(N_{\mathrm{sph}}(\ell) = \tfrac13 r_3(\ell^2+1) = 4H(4(\ell^2+1))\) with the
free \(V_4\)-action (spherical §1–2 **only** — the proved census, not the
shape bijection, not the spherical phase). Present the triptych as one
statement: one arrangement, three geometries, three class-number counts.

**§4 The involution and the class formula.** \(\sigma(X) = \bar X^{-1}\)
identified (unitary involution, Cartan embedding, inversion in the circle);
\(\operatorname{tr} = -2\alpha\); circle \(\to\) closed geodesic of length
\(2\operatorname{arccosh} n\); the class formula
\(\hat\sigma[f] = [\mathfrak{r}_n][f]^{-1}\) on **primitive** classes with the
full proof (Lemma A's closed-form \(X\) realizing a prescribed circle is a key
reusable tool — display it). Imprimitive strata: one remark, labeled
verified/unproved, no more.

**§5 Composition in circle language.** Inverse = mirror; composition =
CRT/magnification; 2-torsion = mirror-symmetric circles; the matrix recipe
\(W_X = \tfrac i2(Z_X + nI)\).

**§6 The phase invariant.** The six-coordinate system
\((\alpha, \beta_1, \beta_2, \arg u)\) on
\(\Gamma\backslash\mathrm{SL}_2(\mathbb{C})/\Gamma\); invariance and the
functional equations (laws 1–2); \(\hat\sigma\)-pairs as Heegner pairs on
\(X_0(\tfrac{n-1}2)\times_{X(1)}X_0(\tfrac{n+1}2)\) with discriminant coupled
to the level. Then the first-power descent as the section's theorem:
\(\omega_f \equiv 1\), hence \(u_f = \Phi_y/\Phi_x(\beta_1,\beta_2)\) — the
derivative of the modular correspondence at its Heegner pair — and the full
dihedral Galois law at first power. Consequences: integer level polynomials
\(Q_n\), irreducibility (proved per computed level, \(n \le 17\); general
\(n\) reduced to distinctness \(T = 1\) — state as open), and the
Gross–Zagier **support** theorem for the denominators (valuations stay open —
say so; exact valuations and the unit theory are Paper II).

**§7 The Euclidean phase.** \(N_e(n) = 2h(-4n^2)\); \(j\)-products
\(H_{-4n^2}^2\); the root-free \(u^2 = -12\beta(\beta-1728)g_2\)-form and the
cocycle-free translation law; **monic** integer level polynomials;
irreducibility for **every** \(n\) (Theorem 5 with Lemma T and archimedean
dominance — the contrast with §6's per-level status is a highlight); the
\(\Delta\)-mass law (Theorem 4) with its elementary proof. The first-power
polynomial \(P^{(2)}_n(x^2)\): certified \(n \le 16\), open in general —
one clearly-labeled paragraph.

**Sign law and DIT (conditional).** `phase-atlas.md` with
`scripts/make_phase_atlas.py` and `scripts/dit_comparison.py` may still live
only on the branch `claude/phase-atlas-schmidt-circles-gc1nfc`
(check `git branch -r | grep phase-atlas`); that branch predates the current
main, so do **not** merge it wholesale — if the files are absent from main,
cherry-pick or port the document and its two scripts only, re-run their
selftests, and then include: (a) the proved divisor-class sign law of its
Theorem 1 in §6; (b) the certified DIT non-relation as the introduction's
novelty statement (the phase is not expressible through
Duke–Imamoğlu–Tóth cycle-integral invariants at safe PSLQ parameters); (c) a
phase-portrait figure if the atlas figures port cleanly. If porting fails,
proceed without, and keep the novelty discussion qualitative.

**Excluded — do not include:** the entire monoid/operation layer
(half-plane-monoid, atomic-census, product-cocycle, middle-kernel, the
\(3/\pi\) split, Apollonian addresses), the entire spectral layer
(spectral-geometry, spectral-outlook, trace-formula anchors), the spherical
shape bijection and all spherical phase laws, Zagier-trace slices, the
Kronecker limit formula, elliptic units, \(R_f\)/\(w_f\)/Robert index (all
Paper II), and every merely-certified spherical/imprimitive statement beyond
the flagged remarks.

## Method and integrity

- Import proofs from the repo documents, rewriting for journal prose; never
  weaken a "verified to \(n \le X\)" into a theorem or vice versa. Keep the
  proved / certified separation explicit in every theorem environment
  (certified results go in `Proposition ... (computational)` or tables with
  the verifying script cited).
- **Verification pass**: every numerical claim that survives into the paper
  must be re-runnable. The known gap: `scripts/first_power_descent.py` is
  referenced by first-power-descent.md §6 but was never committed — rewrite
  it (its content is specified there: exact \(\Phi_m\) construction,
  \(\omega_f = 1\) exact checks odd \(n \le 21\), \(u = \Phi_y/\Phi_x\)
  numerics, exact \(\Pi_n\) with irreducibility) and commit it **before**
  citing it. Re-run `--selftest` of every script the paper cites and record
  runtimes in a short "machine verification" appendix.
- **Literature diligence** (required before the introduction is finalized):
  Stange's Schmidt-arrangement papers (how much of §2–3 is folklore there —
  cite precisely and delimit what is new), Elstrodt–Grunewald–Mennicke
  (Bianchi background), Cox (ring class fields), Gross–Zagier *On singular
  moduli* (for §6's support theorem), Sarnak (reciprocal geodesics, for §4's
  geodesic remark), Duke–Imamoğlu–Tóth (for the novelty discussion). Record
  the searches and verdicts in the PR/commit message or a `NOTES.md` beside
  the paper.
- House bookkeeping: add the paper to CLAUDE.md (paper-planning section) and
  mark outlook.md §4 progress. Keep the paper's numbering of theorems stable
  once written — Paper II will cite them.

## Deliverables

- `papers/1-schmidt-circles/schmidt-circles.tex` + `references.bib`
  (+ figures if the atlas port succeeds), compiling or syntax-validated.
- `scripts/first_power_descent.py` recovered, with `--selftest`.
- A short machine-verification appendix in the paper listing script,
  statement verified, precision, and runtime.
- CLAUDE.md + outlook.md updates; a summary of the literature-diligence
  verdicts.
- State plainly at the end: what the paper claims as new, and where each
  claim is proved in the repo.
