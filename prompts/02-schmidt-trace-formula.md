# Prompt: the Schmidt relative trace formula for \((\mathrm{PSL}_2(\mathbb{Z}[i]), \mathrm{PSL}_2(\mathbb{Z}))\)

Read `CLAUDE.md` first, then [spectral-outlook.md](../spectral-outlook.md) §0–§2 (the
program and its risks, stated honestly there), [spectral-geometry.md](../spectral-geometry.md)
(the verified state of the art — especially §2 Eisenstein theory, §9 periods/distinction,
§13 counting, §14 gaps G2–G3), [hyperbolic-counting.md](../hyperbolic-counting.md), and
[involution.md](../involution.md) §2. Follow the CLAUDE.md guard rails.

## Mission

Establish the exact identity ("Schmidt trace formula") obtained by integrating the
automorphic kernel \(K(x,y) = \sum_{\gamma\in\mathrm{PSL}_2(\mathbb{Z}[i])} k(x,\gamma y)\) of the
Picard orbifold \(M = \mathrm{PSL}_2(\mathbb{Z}[i])\backslash\mathbb{H}^3\) over \(Y \times Y\),
\(Y\) the immersed modular surface (area \(\pi/3\)), and decomposing \(\gamma\) into
\(H\)-double cosets, \(H = \mathrm{PSL}_2(\mathbb{Z})\):

$$
\sum_j |P_Y(\phi_j)|^2\,h(t_j) + (\text{Eisenstein})
\;=\;
(\text{2D trace formula of } Y \text{ from } HeH)
+ \sum_{n\ge2} c\,H(n^2-1)\,F_h(\operatorname{arccosh} n)
+ (\alpha = 1\ \text{regularized}),
$$

whose geometric side is **exact integer arithmetic** — the multiplicities \(H(n^2-1)\) at
distances \(\operatorname{arccosh} n\) are proved in hyperbolic-counting.md; the bookkeeping
constant \(c\) must be derived, not guessed (orientations; the \(\mathrm{PSL}/\mathrm{PGL}\)
strata — even levels live in \(i\mathcal{S}\)).

Do the work in stages; each stage is a self-contained, committable deliverable.

## Stage 0 — numerics baseline (half a day)

Extend `scripts/hurwitz_sum.py` to \(X = 10^4\) via the class number formula (not form
enumeration). Fit the second-order term of
\(\sum_{n\le X}H(n^2-1) - \tfrac{\pi}{12G}X^2\) and test the error-exponent hypotheses
\(X^{1}\) vs \(X^{4/3}\) (spectral-outlook §1 deliverable (b), §9.1). PSLQ any stable
\(X\)-coefficient against \(\{\pi/G,\ \zeta'(2)/\zeta(2),\ \log\text{-constants}\}\) under the
safety rules. This data disciplines every later claim about the spectral side.

## Stage 1 — the Eisenstein period lemma (gap G3; the key self-contained lemma)

Compute in closed form the regularized period \(P_Y^{\mathrm{reg}}(E(\cdot, s))\) of the
Picard Eisenstein series over \(Y\), by unfolding over \(H\backslash\Gamma/\Gamma_\infty\)
with Zagier-style renormalization (the period diverges — the cusp of \(Y\) sits inside
the cusp of \(M\)). Expected shape: a ratio of completed \(\zeta\)-functions; the
binary-Hermitian-form zeta values of Elstrodt–Grunewald–Mennicke (Math. Ann. 277
(1987)) and Flórez–Karabulut–Vu are the predicted answer's home. Mandatory anchors:

- residue at \(s = 1\) of the Jacquet–Lai shape \(\operatorname{area}(Y)/(2\operatorname{vol} M)\)-type
  (\(\operatorname{vol} M = G/3\), scattering \(\varphi(s) = \tfrac{\pi}{s}\zeta_K(s)/\zeta_K(s+1)\));
- numerical verification of the closed form at several real \(s > 1\) by direct
  high-precision summation (spectral-outlook §9.4).

This lemma alone is publishable groundwork and also unlocks the point-to-plane
lattice count (Stage 4 fallback).

## Stage 2 — the compact-support identity (derivation + machine check)

For \(k\) of small support only finitely many levels contribute. Derive the identity by
the classical pre-trace argument plus the incidence bookkeeping of
hyperbolic-counting.md (this is where the constant \(c\) gets pinned). Machine-check the
**geometric side** exactly: evaluate \(\int_{Y\times Y}K\) numerically for 2–3 explicit
small-support kernels and match \(\sum_n c\,H(n^2-1)F_h(\operatorname{arccosh} n)\) term by
term using `alpha_circles.py` data; check the \(h \to\) counting limit against
`hurwitz_sum.py`. Also derive the **plane-pair spherical transform** \(F_h\) on
\(\mathbb{H}^3\) explicitly (Fock coordinates; elementary) — needed both here and in
Stage 3, and worth its own verified subsection.

## Stage 3 — the full identity (port Lekkas–Petridis)

Their 2D double-coset decomposition, orbital transforms, and error management
(arXiv:2509.12902) are the blueprint; Martin–McKee–Wambach (IJNT 7 (2011)) gives the
period-spectrum \(\leftrightarrow\) ortholength-spectrum principle. The genuinely new
technical content, in order of risk:

1. the \(\alpha = 1\) Ford stratum (tangent planes \(\leftrightarrow\) unipotent double
   cosets): a parabolic-type term to be regularized **together with** the continuous
   spectrum — this is where Zagier/Arthur-style truncation enters, and where the
   Stage 1 lemma is consumed;
2. the regularized Eisenstein \(\times\) Eisenstein inner products;
3. the \(\mathrm{SL}_2\)-packet subtleties: work classically with the four symmetry
   classes \(D, G, C, H\) of Bianchi forms rather than adelically, per
   Anandavardhanan–Prasad (the period is not always factorizable and can vanish on
   abstractly-distinguished representations — do not import \(\mathrm{GL}_2\)
   statements blindly). Predict which classes can carry \(P_Y \ne 0\) (base-change
   lifts with central character \(\chi_{-4}\), by Flicker–Zinoviev) and state the
   prediction as a testable claim.

Deliverable: the exact identity with every term defined, plus internal consistency
checks (small-support positivity; the counting limit; the constant-eigenfunction term
reproducing \(\tfrac{\pi}{12G} = \tfrac{\operatorname{area}(Y)^2}{4\pi\operatorname{vol}(M)}\) — anchor 2
of spectral-outlook §0).

## Stage 4 — consequences (stretch; do not block Stages 1–3 on it)

(a) \(\sum_{n\le X}H(n^2-1) = \tfrac{\pi}{12G}X^2 + O(X^{\theta})\) with an explicit
spectral \(\theta\) (first target \(4/3+\varepsilon\); the Eisenstein contribution is the
likely bottleneck — compare against Stage 0 data). (b) The relative Weyl law
\(\sum_{t_j\le T}|P_Y(\phi_j)|^2\). (c) Package the free corollary: the ortholength
spectrum of \(Y\) in \(M\) is \(\{2\operatorname{arccosh} n\}\) with \(H(n^2-1)\)-multiplicities —
the first exactly computed ortholength spectrum in dimension 3 (spectral-outlook §5).

**Fallback if Stage 3 stalls**: the point-to-plane lattice count for the Picard group
(spectral-outlook §2(b)) using only the Stage 1 lemma — a clean standalone paper, and
nothing from Stages 0–2 is wasted.

## Deliverables and conventions

- One house-style document per completed stage (or one growing
  `schmidt-trace-formula.md` with staged sections), each with statements, proofs,
  machine-verification section, and honest status labels (proved / certified /
  derived-modulo-X).
- Scripts: extend `hurwitz_sum.py`; add `scripts/eisenstein_period.py` (Stage 1) and
  `scripts/rtf_compact_support.py` (Stage 2) with selftests.
- Update CLAUDE.md (ledger + scripts table). Cite only literature verified per the
  spectral-geometry.md citation policy; mark anything unverified `[UNVERIFIED]`.
- Do not silently change conventions: eigenvalue normalization \(\lambda = 1+t^2\),
  \(\operatorname{vol}(M) = G/3\), \(\operatorname{area}(Y) = \pi/3\) as in spectral-geometry.md.
