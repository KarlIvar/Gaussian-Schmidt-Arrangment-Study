# Prompt: genus-refined Gross–Zagier — the exponent law of the phase dressing

Read `CLAUDE.md` first, then [phase-kronecker-limit.md](../phase-kronecker-limit.md)
§4(b), §6 and §8 (the \(\Sigma\)-tables and their open exponent law),
[moduli-invariants.md](../moduli-invariants.md) §5.7 and §5.9 (denominators and the
factored leading coefficients \(d_n'\)),
[euclidean-moduli-invariants.md](../euclidean-moduli-invariants.md) §5.5 and §6.1,
and [outlook.md](../outlook.md) §1.1, §1.3, §2.8. Follow the certification guard
rails in CLAUDE.md without exception. Reusable machinery: `scripts/phase_klf.py`
(`identify_real_charsum`, `factor_ratio` — exact \(\mathbb{Q}(\sqrt d)\)
factorizations with Hensel valuations and canonical generators; class groups,
characters, GZ prime sets via `gz_denominators.gz_primes`).

## Mission

Prove the **exponent law of the Gross–Zagier dressing**: the certified tables
express every genus-character \(j\)-sum as
$$
\Sigma_0(\chi) = k\log\varepsilon_d + \sum_p e_p\,\lambda_p,\qquad
\lambda_p = \log|\pi_p/\pi_p'|,
$$
with \(e_p \in \{\pm3\}\) and split support in \(\mathrm{GZ}(D,-3)\) (and even
\(e_p\), support \(\mathrm{GZ}(D,-4)\), for \(\Sigma_{1728}\)). The multiplicities,
supports, and *which* prime of \(\mathbb{Q}(\sqrt d)\) carries which sign are
currently exact data with no law. Derive them from Gross–Zagier's singular-moduli
factorization, genus-character refined — and harvest the three corollaries that
have been open across the repo:

1. **exact denominator valuations** of the hyperbolic phase (moduli §5.7 / outlook
   1.3): the primes and exponents of the leading coefficients \(d_n'\) of \(Q_n\);
2. **the \(\lambda_n\)-valuation law** (euclidean §6.1) — note prompt 04 attacks
   this by volcano geometry; the two routes must agree, which is a cross-check, not
   a redundancy;
3. **the sign law of \(u_f\) on ambiguous classes** (outlook 1.1) — the
   \(n = 11, 15\) anomaly.

## Structural facts to use

1. The data side is **finished and exact**: for every real character of every
   computed level (Euclidean \(n \le 15\), hyperbolic \(n \le 15\)),
   `factor_ratio` returns sign, unit power \(k\), and all \((p, e_p, \pi_p)\),
   verified by exact multiplication in \(\mathbb{Q}(\sqrt d)\). The \(\pm3\)/even
   pattern is the \((4,3)\)-exponent shadow of
   \(u^6 = -\beta^4(\beta-1728)^3\cdot(\Delta\text{-quotient})\).
2. The **unit power \(k\) is not an independent unknown**: the Kronecker limit
   theorems (phase-kronecker-limit Thms 1/3/4) determine
   \(\log|A/B| = \Sigma_x(\chi)\) and the \(L'\)-closed forms determine the
   \(\varepsilon\)-content of the whole \(S(\chi)\); once the finite exponents
   \(e_p\) are proved, \(k\) follows. Only the finite places need new arithmetic.
3. Our products are \(\prod_\mathfrak{c}(j(\mathfrak{c}) - j_0)^{\chi(\mathfrak{c})}\)
   with \(j_0 \in \{0, 1728\}\) — Gross–Zagier differences against the discs
   \(-3\) and \(-4\), but with the **second disc non-fundamental**
   (\(-4n^2\), resp. \(1-n^2\) of conductor 2 at some levels) and **not coprime**
   to the first at \(2\) or \(3\). The classical GZ paper (§5–6, the
   \(\epsilon(\mathfrak{n})\)-formula) handles coprime fundamental pairs and is
   *already genus-refined* — their local multiplicities are indexed by genus
   characters. The extension to arbitrary pairs is **Lauter–Viray**, *On singular
   moduli for arbitrary discriminants* — the key external reference; specialize it
   rather than re-derive.
4. The prime-by-prime shape to prove: for \(p\) split in \(\mathbb{Q}(\sqrt d)\)
   (the genus field piece cut by \(\chi\)), the collision multiplicities of the two
   \(\chi\)-cosets at the two primes \(\mathfrak{p}, \mathfrak{p}'\) above \(p\)
   are separated exactly by the genus character value on the collision ideal — this
   is what \(e_p \ne 0\) measures; inert/ramified \(p\) cancel (proved already in
   phase-kronecker-limit: \(v_\mathfrak{p}(A/A') = 0\) there). The GZ counting
   quantities are the representation numbers of \(\tfrac{e|D| - x^2}4\)-forms that
   `gz_denominators.gz_primes` already enumerates — refine that enumeration to
   count *with* the genus character weight and compare against the exact \(e_p\).
5. For the **sign law**: signs of \(u_f\) on ambiguous classes are computable in
   **exact arithmetic, no numerics**, via the first-power closed form
   \(u_f = \Phi_y/\Phi_x(\beta_1, \beta_2)\) evaluated in
   \(F_1 = \mathbb{Q}[t]/(H_D)\) with real-embedding sign tracking
   (first-power-descent Thm 4.2 route). NOTE: `scripts/first_power_descent.py` is
   the repo's known missing script (CLAUDE.md "Known gap") — rebuilding this exact
   pipeline is part of Phase 1 and **closes that gap as a bonus deliverable**.

## Plan

**Phase 1 — data extension, all exact.**
(a) Rebuild the exact \(\Phi_{r_0}\)-pipeline (integer \(q\)-expansions, Hilbert
class polynomials, \(u(t) \in F_1\)) as `scripts/first_power_descent.py`, matching
the documented interface of first-power-descent.md §6; use it to compute the exact
sign of \(u_f\) on every ambiguous class for all odd \(n \le 41\) — extending the
outlook 1.1 table with proofs instead of numerics.
(b) Extend the \(\Sigma\)-factorization tables: Euclidean \(n \le 20\), hyperbolic
odd \(n \le 21\) (raise dps as needed; certify with the standard margins).
(c) For each split \(p\) in each table row, compute the genus-weighted GZ counting
numbers (the \(\epsilon\)-weighted representation counts of
\(\tfrac{e|D|-x^2}4\)) and fit the empirical law \(e_p = \) (weighted count
difference) — this is the conjecture to prove, stated precisely from data.

**Phase 2 — proof.** Specialize Lauter–Viray to the pairs \((-3, D)\), \((-4, D)\)
with \(D \in \{-4n^2, 1-n^2\}\); refine by the character: the \(\chi\)-weighted sum
collapses the per-class valuations onto the two primes of \(\mathbb{Q}(\sqrt d)\)
with genus-character signs. Deliver the exponent theorem: closed formula for
\(e_p\) (and hence \(k\), by fact 2). Where Lauter–Viray's hypotheses genuinely
fail (if anywhere at our levels), do the local Deuring-lift counting by hand at
that prime and say so.

**Phase 3 — the three corollaries.**
(a) Exact denominator theorem: predict the full factorization of every \(d_n'\)
(moduli §5.9 table, \(n \le 17\)) and of the Euclidean \(P^{(2)}_n(0)\)-supports;
verify against the certified tables **exactly**.
(b) The \(\lambda_n\)-law (euclidean §6.1) — compare with prompt 04's
volcano-route answer if that project has run; the agreement is a theorem-level
cross-check.
(c) The sign law: from the Phase 1 exact sign table plus the proved exponent
theorem, isolate what distinguishes \(n = 11, 15\) (the prediction to test: a
2-adic genus character of the conductor-2 part that the crude prime-disc
decomposition misses). A clean statement here completes the real-quadratic
description of \(u\) on ambiguous classes.

**Verification discipline.** Fit at Euclidean \(n = 9, 11, 13\), verify the proved
law at levels not used in the fit (\(n = 7, 15, 17\)–\(20\); hyperbolic
\(17, 19, 21\)). A law that does not transfer is discarded; a certified non-fit is
recorded with basis, precision, and margins.

## Deliverables

- House-style `gz-dressing.md`: the exponent theorem with proof, the three
  corollaries with their status (proved / certified / failed), the extended exact
  sign table, machine-verification section, outlook.
- `scripts/gz_dressing.py` (+ the rebuilt `scripts/first_power_descent.py`),
  standalone from repo root, `--selftest` re-verifying every displayed identity.
- CLAUDE.md: status ledger (close §5.7-valuations, §6.1, and 1.1 if achieved;
  remove the "Known gap" note once `first_power_descent.py` is restored); scripts
  table rows; outlook.md: mark 1.1/1.3/2.8 accordingly.
