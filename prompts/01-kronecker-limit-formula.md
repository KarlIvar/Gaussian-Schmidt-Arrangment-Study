# Prompt: the Kronecker limit formula for the phase — \(\log|u_f|\) against \(L'(0,\chi)\)

Read `CLAUDE.md` first, then [moduli-invariants.md](../moduli-invariants.md) (§4–5),
[first-power-descent.md](../first-power-descent.md),
[euclidean-moduli-invariants.md](../euclidean-moduli-invariants.md) (§5, outlook 6.6), and
[outlook.md](../outlook.md) §2.3. Follow the certification guard rails in CLAUDE.md
(safe PSLQ, no import-time `mp.dps`, absolute-error criterion) without exception.

## Mission

Decide, first experimentally and then by proof, whether the character sums of the
log-absolute-values of the phase invariants evaluate to explicit combinations of
derivatives of \(L\)-functions at \(s = 0\) — i.e. whether the phases form a
**geometrically defined elliptic-unit / Stark-adjacent system**. Both aspects:

- **Euclidean (do this first)**: \(u_\mathfrak{c}\), \(\mathfrak{c} \in \mathrm{Cl}(\mathcal{O}_n)\),
  \(\mathcal{O}_n = \mathbb{Z}+n\mathbb{Z}[i]\), with the proved closed forms
  \(u_\mathfrak{c}^2 = -12\,\beta(\beta-1728)\,g_2(\Lambda_\mathfrak{c})/g_2(\mathbb{Z}[i])\) and
  \(u^6 = -\beta^4(\beta-1728)^3\,\Delta(\Lambda_\mathfrak{c})/\Delta(\mathbb{Z}[i])\). Here the
  conductor aspect over the fixed field \(\mathbb{Q}(i)\) puts you squarely in
  Robert / Kubert–Lang (ch. 12–13) / Schertz elliptic-unit territory — the strongest
  available literature.
- **Hyperbolic (the headline if it works)**: \(u_f\), \(f \in \mathrm{Cl}(1-n^2)\), with
  \(u_f = \Phi_y/\Phi_x(\beta_1,\beta_2)\) and
  \(u_f^6 = \varepsilon^6\mu^{-12}\,\frac{\beta_1^4(\beta_1-1728)^3}{\beta_2^4(\beta_2-1728)^3}\,
  \frac{\Delta(\tau_1)}{\Delta(\tau_2)}\), \(|\mu|^2 = \varepsilon q_1/q_2\).

## Structural facts to use (all proved in the repo)

1. Hyperbolic: \(u_{\mathfrak{r}f}u_f = 1\) forces
   \(S(\chi) := \sum_f \chi(f)\log|u_f| = 0\) **unless \(\chi(\mathfrak{r}_n) = -1\)** — only
   characters odd on the twist class carry information. \(u_{f^{-1}} = \bar u_f\) gives
   \(S(\chi) = S(\bar\chi)\).
2. Euclidean: the **trivial-character sum is already a theorem** — assemble
   \(\sum_\mathfrak{c} 2\log|u_\mathfrak{c}|\) from Theorem 2 (\(j\)-product \(= H_{-4n^2}^2\)) and
   the \(\Delta\)-mass law (Theorem 4). Use it as the mandatory sanity anchor before
   fitting any nontrivial character.
3. The \(\Delta\)-part of \(\log|u|\) is exactly what the (second) Kronecker limit
   formula computes: for a nontrivial ring class character \(\chi\),
   \(\sum_\mathfrak{c}\chi(\mathfrak{c})\log\bigl(N\mathfrak{a}_\mathfrak{c}^{6}|\Delta(\mathfrak{a}_\mathfrak{c})|\bigr)\)
   is proportional to \(L'(0,\chi)\) of the associated ring-class-field \(L\)-function.
   The \(j\)-value parts \(\sum\chi(\mathfrak{c})\log|j(\mathfrak{c})|\),
   \(\sum\chi(\mathfrak{c})\log|j(\mathfrak{c})-1728|\) are Gross–Zagier collision quantities
   (finite places: \(\log p\) over the GZ sets identified in moduli-invariants §5.7 and
   euclidean §5.5).
4. For **genus characters** \(\chi \leftrightarrow (d_1, d_2)\):
   \(L(s,\chi) = L(s,\chi_{d_1})L(s,\chi_{d_2})\), so \(L'(0,\chi)\) evaluates in closed form
   through \(h(d)\log\varepsilon_d\) (real \(d\)) and \(L(0,\chi_d) = 2h(d)/w(d)\) (imaginary
   \(d\)). The certified pair-sum fields (\(\mathbb{Q}(\sqrt5)\) at \(n = 9, 11\);
   \(\mathbb{Q}(\sqrt{14})\) at \(n = 13\); \(\mathbb{Q}(\sqrt2)\) at \(n = 17\)) tell you which
   fundamental units to expect.

## Plan

**Phase 1 — experiment (Euclidean).** Compute all \(u_\mathfrak{c}\) at 200+ digits for
\(n = 9, 11, 13\) (reuse `scripts/euclidean_moduli_invariants.py` machinery). For every
character \(\chi\) of \(\mathrm{Cl}(\mathcal{O}_n)\) form \(S(\chi)\) and PSLQ against the basis
\(\{L'(0,\chi),\ h\log\varepsilon_d\ (d \mid \text{disc data}),\ \log p\ (p \mid 2n\ \text{and}\
p \in \operatorname{supp}H_{-4n^2}(0)\cup\operatorname{supp}H_{-4n^2}(1728)),\ \log 2\}\).
Compute \(L'(0,\chi)\) via the genus factorization when available, else numerically
(mpmath, functional equation), certified to the working precision. **Fit at one level,
verify the same-shaped identity at a level not used in the fit** — a fit that does not
transfer is discarded.

**Phase 2 — proof (Euclidean).** Assemble the identity from: Kronecker's limit formula
for ring class characters of orders of \(\mathbb{Q}(i)\) (Robert/Schertz normalization),
the GZ factorizations of \(H_{-4n^2}(0)\) and \(H_{-4n^2}(1728)\), and the proved closed
form of \(u^2\). Target theorem: an exact formula for \(S(\chi)\) for every character,
making \(\{u_\mathfrak{c}\}\) an explicit elliptic-unit system up to the identified
\(S\)-integer dressing.

**Phase 3 — hyperbolic transfer.** Repeat the experiment at hyperbolic levels
\(n = 9, 11, 13, 15\) (only \(\chi(\mathfrak{r}) = -1\) characters). The expected new
ingredient is \(\log\varepsilon_{n^2-1}\) (from \(\varepsilon^6\mu^{-12}\), whose absolute
values are controlled by the Norm Lemma) alongside \(L'(0,\chi)\) of ring class
characters of disc \(1-n^2\). If a clean law appears, prove it by the same assembly,
using first-power-descent.md's algebraic closed form.

**Optional extension.** The spherical aspect (spherical-moduli-invariants.md outlook
item 7): \(\log|u^2|\) for disc \(-4(\ell^2+1)\) with the \(\varepsilon_\ell^4\) cap-swap as
archimedean correction. Only if Phases 1–2 succeed.

## Deliverables

- A house-style document `phase-kronecker-limit.md`: statements, proofs (or precise
  conjectural identities with certified numerical evidence), machine-verification
  section, outlook. State plainly which parts are proved / certified / failed.
- `scripts/phase_klf.py`, standalone from repo root, with a selftest re-verifying every
  displayed identity; certification margins printed.
- Update CLAUDE.md: status ledger entry and scripts table row.
- **Negative results are results**: if some \(S(\chi)\) matches nothing in the basis,
  record the certified non-fit (basis used, precision, margins) — that constrains the
  theory and belongs in the document.
