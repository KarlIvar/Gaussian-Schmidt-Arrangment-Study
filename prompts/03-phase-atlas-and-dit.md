# Prompt: the phase atlas, and the Duke–Imamoğlu–Tóth comparison

Read `CLAUDE.md` first, then [moduli-invariants.md](../moduli-invariants.md) §4–5,
[first-power-descent.md](../first-power-descent.md), [hyperbolic-counting.md](../hyperbolic-counting.md),
[euclidean-moduli-invariants.md](../euclidean-moduli-invariants.md) §5.2, and
[outlook.md](../outlook.md) §1.1, §1.6, §2.4, §2.6. Follow the CLAUDE.md guard rails.

One session, one dataset, two deliverables: computing \(u_f\) for every class at many
levels is the expensive step, and both halves of this task consume exactly that data.

## Part A — the phase atlas

Render the Schmidt circles of the ideal triangle \(T = (0, 1, \infty)\) colored by the
phase — "the phase portrait of the Schmidt arrangement" — the natural signature figure
for the phase paper, and the cheapest way to expose patterns the formulas hide.

**Computation.** For each odd level \(3 \le n \le 41\) (skip nothing; include a couple of
larger spot levels, e.g. \(n = 101\), if cheap), enumerate the level-\(n\) circles in
\(T\) with `scripts/alpha_circles.py` machinery, map each circle to its form class, and
compute \(u_f\) per class at 60+ digits by the numerically cheap closed form
\(u_f = -\varepsilon\,\mu^{-2}\,h_2(\mathfrak{b}_1)/h_2(\mathfrak{b}_2)\)
(moduli-invariants §5.5; compute \(j' = -2\pi i E_4^2E_6/\Delta\) via theta constants /
Eisenstein series — never trust a library default). Cross-check at \(n \le 17\) against
the exact \(\Phi_y/\Phi_x\) route of first-power-descent.md — any mismatch is a bug in
the atlas, full stop.

**Rendering** (matplotlib; PNGs into `figures/`, generator into
`scripts/make_phase_atlas.py`):

1. Per-level panels of \(T\) (half-plane and the \(D_3\)-symmetric disk model, as in
   `alpha_circles.py`): fill color = hue by \(\arg u_f\), and a companion panel
   shaded by \(\log|u_f|\) (normalize per level: \(\log|u_f|/\log|u_1|\), the principal
   class being the max by far).
2. Overlays: mark ambiguous classes (where \(u_f \in \mathbb{R}\)) with their **sign**;
   draw the mirror lines \(\operatorname{Re}z \in \{0, \tfrac12, 1\}\) and \(|z|=1\), \(|z-1|=1\)
   (circle-composition.md §2 ties 2-torsion to them); connect or co-mark
   \(\mathfrak{r}_n\)-twin pairs (\(u_{\mathfrak{r}f} = 1/u_f\): antipodal hues, inverted radii —
   verify this is visible).
3. A multi-level contact sheet (all \(n \le 41\)) for pattern-hunting.
4. Optional Euclidean companion (euclidean 6.5): disks of the unit square colored the
   same way, with the **proved** R/I center-criterion (euclidean §5.2) overlaid as a
   correctness check on the rendering pipeline.

**Free rider — the sign law data (outlook §1.1 first move).** While every \(u_f\) is in
memory, extend the ambiguous-class sign table from \(n \le 17\) to all odd \(n \le 41\),
and test the stated hypotheses: a 2-adic / conductor-aware genus character explaining
the anomalous patterns at \(n = 11, 15\) (identical class groups at \(n = 11, 13\),
different sign behavior — that is the discriminating pair). Record the table and the
verdict; if a candidate character works at every level, state it as a conjecture with
the certified table as evidence.

## Part B — the Duke–Imamoğlu–Tóth comparison (outlook §2.6)

**Question.** Is \(u_f\) (or \(\log|u_f|\), or \(\arg u_f\), or the pair-sums
\(S_x = u_f + 1/u_f\)) expressible through DIT-type real-quadratic invariants? The
level-\(n\) circle data lives on the real quadratic side too: the attached closed
geodesic has trace \(-2n\), length \(2\log\varepsilon_n\), \(\varepsilon_n = n+\sqrt{n^2-1}\) —
discriminant \(4(n^2-1)\). **Either outcome is valuable**: a match imports the DIT
machinery (cycle integrals, mock modular forms, linking numbers) into the phase
theory; a certified mismatch is the novelty claim the phase paper's introduction needs.

**Do.**

1. Compute, for \(4 \le d = 4(n^2-1)\), \(n = 3, \dots, 17\), the DIT invariants of
   discriminant \(d\): cycle integrals of \(j\) (per class and the trace
   \(\mathrm{Tr}_d(j)\)), cycle integrals of \(\log|\Delta|\)/Eisenstein type (equivalently
   Zagier-reduced continued-fraction / Dedekind-sum data of the class), and the
   Rademacher symbol / linking number of the associated modular knot. Implement from
   the definitions (Zagier reduction of indefinite forms; numerical integration along
   the geodesic at 50+ digits) — do not trust unverified published tables beyond the
   anchors you can re-derive.
2. Compare against the phase data, per class and in aggregate: PSLQ (safe parameters)
   of \(\log|u_f|\) and \(S_x\) against cycle-integral bases; test \(\arg u_f\) against
   linking-number / Rademacher data; test the \(\hat\sigma\)-pair structure against DIT's
   form-class pairings. Any candidate identity must transfer to a level not used in
   the fit.
3. Probe the cusp degeneration (outlook §2.4) at the cheap end: as the phase theory
   degenerates toward the Ford stratum, check numerically whether the
   \(\arg u\)-cocycle data reduces to classical Dedekind sums — one afternoon, and it
   would explain the \(\mu_{12}\) bookkeeping if true.
4. Write the verdict either way: a certified identity, or a documented non-match
   (bases tried, precisions, margins) — the latter is the citable statement "\(u_f\) is
   not a linear shadow of known cycle-integral invariants".

## Deliverables

- `scripts/make_phase_atlas.py` (figures + selftest that re-checks the twin and
  conjugation laws pixel-independently on the computed data) and
  `scripts/dit_comparison.py` (invariant computations + comparison log with
  certification margins).
- Figures in `figures/` (per-level panels, contact sheet, Euclidean companion if done).
- One house-style document `phase-atlas.md` containing: the atlas with observations
  (each observed pattern either matched to a proved law with a cross-reference, or
  flagged as new), the extended sign table and the sign-law verdict, and the DIT
  comparison section with its verdict. Statements labeled proved / certified /
  observed.
- Update CLAUDE.md (ledger, scripts table). If a genuinely new pattern shows up in the
  atlas, add it to the open-problems ledger with a pointer to the figure.
