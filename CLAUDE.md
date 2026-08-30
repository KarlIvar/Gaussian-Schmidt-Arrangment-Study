# CLAUDE.md — Gaussian Schmidt Arrangement Study

## What this project is

A research monorepo on the **Gaussian Schmidt arrangement**
\(\mathcal{S} = \mathrm{PSL}_2(\mathbb{Z}[i])\cdot\hat{\mathbb{R}}\) — the orbit of the
extended real line under the Bianchi group — produced by several agent sessions. The
goal is one or more papers with new results. The repo consists of self-contained
research documents (Markdown + LaTeX) in the root, verification scripts in
`scripts/`, and figures in `figures/`. **Every claim with finite content is
machine-verified** by a named script; documents state their verification status
explicitly (proved / certified / experimental).

Three literature probes (see [spectral-geometry.md](spectral-geometry.md) §14) confirm
the core results are **not in the literature**: only 3 published papers mention
"Schmidt arrangement" (Stange et al.), none touch the level stratification, the
involution class formula, the phase invariant, or spectral geometry.

## The objects, in one page

- **Classification** ([circle-classification.md](circle-classification.md)): a circle
  of radius \(1/2n\) with center \(\zeta/2n\) lies in \(\mathcal{S}\) iff
  \(\zeta \in \mathbb{Z}[i]\), \(\zeta \equiv i \pmod 2\), \(|\zeta|^2 \equiv 1 \pmod{4n}\).
  Lines: \(\operatorname{Im} z = k\). The \(\mathrm{PGL}_2\)-arrangement is
  \(\mathcal{S} \sqcup i\mathcal{S}\) (disjoint).
- **The level** \(\alpha\): for a circle in \(\mathbb{H}\), \(\alpha = y = \operatorname{Im}\zeta
  = \coth(\text{hyperbolic radius})\) — an integer \(n \ge 2\); also
  \(\alpha = -\tfrac12\operatorname{tr}(X\bar X^{-1})\) (trace on the Cartan embedding), and the
  inversive product with \(\hat{\mathbb{R}}\).
- **Level-\(n\) circles \(\leftrightarrow\) quadratic forms**: \(\omega_{q,x} \leftrightarrow
  f_\omega = (q,-x,m)\), positive definite of discriminant \(D = 1-n^2\); hyperbolic
  centers are CM points of disc \(1-n^2\); classes carry the class group
  \(\mathrm{Cl}(1-n^2)\).
- **The involution** \(\sigma(X) = \bar X^{-1}\): the unitary (second-kind) involution
  adjoint to \(M_0 = \binom{0\ i}{-i\ 0}\); \(X\bar X^{-1}\) = Cartan embedding = inversion
  in the circle \(X(\hat{\mathbb{R}})\). On classes:
  \(\hat\sigma[f] = [\mathfrak{r}_n][f]^{-1}\), where
  \(\mathfrak{r}_n = [(\tfrac{n-1}2, 0, \tfrac{n+1}2)]\) is a canonical ambiguous class —
  formally the Atkin–Lehner action on Heegner points.
- **The phase** \(u_f = \varepsilon\Theta_f\) (\(\varepsilon = n+\sqrt{n^2-1}\)): the sixth
  coordinate on double cosets \(\Gamma\backslash\mathrm{SL}_2(\mathbb{C})/\Gamma\),
  \(\Gamma = \mathrm{SL}_2(\mathbb{Z})\). Closed form:
  \(u_f = \Phi_y/\Phi_x(\beta_1,\beta_2)\) — the derivative of the modular correspondence
  \(X_0(\tfrac{n-1}2)\) at the Heegner pair. Algebraic, norm-1, dihedral Galois action at
  first power, irreducible level polynomials, Gross–Zagier-prime denominators.
- **The monoid** \(\Omega = \{X \in \mathrm{SL}_2(\mathbb{Z}[i]) : X(\mathbb{H})\subseteq\mathbb{H}\}\):
  rigid (unique factorization), units \(\mathrm{SL}_2(\mathbb{Z})\), atoms = disks of the
  Apollonian gasket in the strip; every form class of disc \(1-n^2\) gets an
  **Apollonian address** (word in atom classes).
- **The Euclidean theory** (translations \(\times\ \mathrm{SL}_2(\mathbb{Z})\) instead of
  \(\Gamma\times\Gamma\)): disks of curvature \(2n\) mod translation \(\leftrightarrow\)
  \(\mathrm{Cl}(\mathbb{Z}+n\mathbb{Z}[i])\), so \(N_e(n) = 2h(-4n^2)\); its phase
  \(u = \Theta/\Omega_{\text{lemn}}\) has **monic** integer level polynomials, irreducible
  for **every** \(n\) (proved).
- **The spherical theory** (finite \(\Gamma_{\mathrm{sph}} = \mathrm{SU}(2)\cap\mathrm{SL}_2(\mathbb{Z}[i])\)
  \(\times\ \mathrm{SL}_2(\mathbb{Z})\)): on the Riemann sphere, radii quantize to
  \(\cot\theta = \ell = q+m \in \mathbb{Z}\); census \(N_{\mathrm{sph}}(\ell) = \tfrac13 r_3(\ell^2+1)
  = 4H(4(\ell^2+1))\); shape classes sweep **all** form classes of disc \(-4(\ell^2+1)\)
  exactly once (imprimitive strata included); the phase is normalized by the
  negative-Pell unit \(\varepsilon_\ell = \ell+\sqrt{\ell^2+1}\), with cap-swap law
  \(u^2_- = \varepsilon^4 u^2_+\), level norm \(u^2\sigma(u^2) = m_\ell\) (GZ-supported), and the
  **half-orbit phenomenon**: for \(\ell \ge 2\) the phases are only half of their Galois
  root system — the other half are Pell-twisted virtual partners \(m_\ell/u^2\) realized
  by no circle, with a unit-valued (infinite-order) Galois cocycle.

Key constants: Catalan's \(G = L(2,\chi_{-4})\); \(\operatorname{vol}(\mathrm{PSL}_2(\mathbb{Z}[i])\backslash\mathbb{H}^3) = G/3\);
\(\operatorname{area}(Y_{\mathrm{mod}}) = \pi/3\); Apollonian dimension \(\delta = 1.3056867\ldots\)

## Document map (read in this order)

Foundational layer (elementary, self-contained):

| doc | headline result |
|---|---|
| [circle-classification.md](circle-classification.md) | Diophantine classification of \(\mathcal{S}\); single equation \(u^2+v(v+1)=nm\); \(\mathrm{PGL}_2\) variant |
| [euclidean-counting.md](euclidean-counting.md) | \(N_e(n) = n\prod_{p\mid n}(1-\chi_{-4}(p)/p)\); \(\sum_{n\le X}N_e(n) = X^2/2G + O(X\log X)\) |
| [hyperbolic-counting.md](hyperbolic-counting.md) | level-\(n\) circles = hyperbolic circles of radius \(\operatorname{arcoth} n\) at CM points of disc \(1-n^2\); weighted ideal-triangle count \(= 3H(n^2-1)\) |
| [involution.md](involution.md) | \(\sigma\) identified (unitary involution, Cartan embedding, Gelfand pair); trace \(=-2\alpha\); 3 twisted classes of det \(-1\) Hermitian forms; class formula stated; circle \(\to\) closed geodesic of length \(2\operatorname{arccosh} n\) |
| [class-formula-proof.md](class-formula-proof.md) | full proof of \(\hat\sigma[f]=[\mathfrak{r}_n][f]^{-1}\); Lemma A = closed-form \(X\in\mathrm{SL}_2(\mathbb{Z}[i])\) realizing any prescribed circle |
| [circle-composition.md](circle-composition.md) | Gauss composition in circle language: inverse = mirror; composition = CRT/magnification; 2-torsion = mirror-symmetric circles; matrix recipe \(W_X = \tfrac i2(Z_X+nI)\) |

Phase layer (CM machinery):

| doc | headline result |
|---|---|
| [moduli-invariants.md](moduli-invariants.md) | the 6-coordinate system \((\alpha,\beta_1,\beta_2,\arg u)\); \(\hat\sigma\)-pairs on the fiber product \(X_0(\tfrac{n-1}2)\times_{X(1)}X_0(\tfrac{n+1}2)\), discriminant coupled to level; laws \(u_{f^{-1}}=\bar u_f\), \(u_fu_{\mathfrak{r}f}=1\); Theorem A (dihedral equivariance of \(u^{12}\)); GZ-prime denominators (§5.7); level polynomials \(Q_n\) (§5.9) irreducible at computed levels (§5.10) |
| [first-power-descent.md](first-power-descent.md) | \(u_f = \Phi_y/\Phi_x(\beta_1,\beta_2)\) (derivative of the modular correspondence); Galois law at **first power**; \(\omega_f\equiv1\); irreducibility re-proved by exact arithmetic, no PSLQ |
| [euclidean-moduli-invariants.md](euclidean-moduli-invariants.md) | \(N_e(n)=2h(-4n^2)\); \(j\)-values \(=H_{-4n^2}^2\); trace slice \(t(4n^2)\); lemniscatic phase; \(\Delta\)-mass law (Thm 4, proved); \(P^{(2)}_n\) irreducible for **every** \(n\) (Thm 5, proved) |
| [phase-kronecker-limit.md](phase-kronecker-limit.md) | character sums of \(\log\lvert u\rvert\) = elliptic-unit theory: \(S(\chi) = -2L'(0,\chi) + \tfrac23\Sigma_0 + \tfrac12\Sigma_{1728}\) (Euclidean, proved); hyperbolic \(S(\chi) = 0\) for \(\chi(\mathfrak{r})=+1\), \(= -4L'(0,\chi)+\tfrac43\Sigma_0+\Sigma_{1728}\) for odd \(\chi\) (proved, \(\varepsilon\) and \(\mu\) cancel); genus-character closed forms \(c\, h\log\varepsilon_d\) with the Euclidean field \(=\mathbb{Q}(\sqrt n)\); exact GZ-supported \(S\)-unit dressing; certified \(\Delta\)-mass polynomials \(D_n\) |
| [spherical-moduli-invariants.md](spherical-moduli-invariants.md) | third geometry: \(\cot\theta = \ell\); census \(4H(4(\ell^2+1))\); shape polynomial \(\prod_{f^2\mid\ell^2+1}H_{-4(\ell^2+1)/f^2}\) (each stratum once); trace slice \(t(4\ell^2+4)\); Pell-unit phase, cap-swap \(u^2_-=\varepsilon^4u^2_+\), level norm \(m_\ell\) with GZ valuation law; **half-orbit phenomenon** (phases = half the Galois root system, golden-ratio cocycle at \(\ell=2\)) |

Monoid / operation layer:

| doc | headline result |
|---|---|
| [half-plane-monoid.md](half-plane-monoid.md) | \(\Omega\) is rigid; laminarity (inversive products are odd integers); atoms = Apollonian gasket; \(\alpha\)-spectrum = Apollonian orbit of \((-1,1,1,1)\) |
| [atomic-census.md](atomic-census.md) | cell decomposition of \(\Gamma\backslash\Omega/\Gamma\); Apollonian addresses; level pairing \(\alpha(XY)=\langle M_{X^{-1}},M_Y\rangle\); superadditivity \(\varepsilon_{\alpha(XY)}\ge\varepsilon_{\alpha(X)}\varepsilon_{\alpha(Y)}\); pure-twist \([f]\mapsto[\mathfrak{r}_nf]\) + word reversal; Ford census (\(\varphi(c)\) at \(n=2c^2+1\)); extremal depth \(\tfrac{n+1}2\); the \(3/\pi\) density split |
| [product-cocycle.md](product-cocycle.md) | Gram rigidity (product's \(\beta_1\) remembers only \(X\), \(\beta_2\) only \(Y\)); aligned stratum splices invariants; phase cocycle \(\Theta(W)\Theta(X)\Theta(Y)=KQ^2\) (anti-additive, explicit geometric coboundary) |
| [middle-kernel.md](middle-kernel.md) | \(\sum_\gamma f(X\gamma Y)\) for radial \(f\) = automorphic kernel of the modular surface at a Heegner pair; \(\alpha(X\gamma Y)=ab+\sqrt{(a^2-1)(b^2-1)}\cosh d(z_X,\gamma z_Y)\); exact rank-one pole; obstruction = cusp forms; Green's-function/GKZ connection |

Spectral layer:

| doc | headline result |
|---|---|
| [spectral-geometry.md](spectral-geometry.md) | verified literature survey of the Picard orbifold (spectra, PGT, periods, sup-norms, thin groups, CFs); **seven verified gaps G1–G7** the project can fill |
| [spectral-outlook.md](spectral-outlook.md) | the program: Schmidt relative trace formula (flagship), effective counting, Eichler–Selberg over \(\mathbb{Z}[i]\), ortholength spectrum, phase as spectral object, transfer operators; anchor \(\sum_{n\le X}H(n^2-1)\sim\frac{\pi}{12G}X^2 = \frac{\operatorname{area}(Y)^2}{4\pi\operatorname{vol}(M)}X^2\) |

Planning: [outlook.md](outlook.md) — the master outlook (small/medium/large questions, consolidation plan, recommended order). Each of atomic-census, product-cocycle, middle-kernel, euclidean-moduli-invariants also carries its own outlook section.

## Main theorems at a glance

1. **Classification** of \(\mathcal{S}\) (necessity + descent sufficiency). Proved.
2. **Euclidean count** \(N_e(n)\) multiplicative; Catalan asymptotics. Proved.
3. **Hyperbolic count** \(3H(n^2-1)\) in the ideal triangle (incidence lemma). Proved.
4. **Class formula** \(\hat\sigma[f] = [\mathfrak{r}_n][f]^{-1}\) on primitive classes. Proved
   (+ verified \(n\le41\)); imprimitive strata at the ramified prime 2 still ad hoc.
5. **Composition/inversion recipes** in circle language = Gauss composition. Proved.
6. **Unique factorization** of \(\Omega\); atoms = Apollonian gasket. Proved.
7. **Cells and addresses**; superadditivity; pure twist + word reversal; Ford census
   closed forms; depth bound \(\tfrac{n+1}2\) with exact attainment. Proved.
8. **Six-invariant coordinates**; functional equations of \(\Theta\); simultaneous
   modular equations for \(\hat\sigma\)-pairs. Proved.
9. **Phase laws**: \(u_{f^{-1}} = \bar u_f\), \(u_fu_{\mathfrak{r}f} = 1\). Proved.
10. **First-power dihedral Galois law** \(\sigma(u_f) = u_{f^{e(\sigma)}\mathfrak{c}(\sigma)}\)
    via \(u_f = \Phi_y/\Phi_x\). Proved. Level polynomial \(Q_n\) irreducible: proved for
    every **computed** level (\(n\le17\)); general \(n\) needs only distinctness (\(T=1\)).
11. **GZ denominators**: support of denominators = Gross–Zagier primes of
    \((1-n^2,-3)\cup(1-n^2,-4)\). Support proved via §5.7 mechanism; exact valuations open.
12. **Euclidean structure** \(N_e(n) = 2h(-4n^2)\); \(j\)-product \(=H_{-4n^2}^2\). Proved.
13. **\(\Delta\)-mass law** (Euclidean Thm 4): closed form for
    \(\prod_\mathfrak{c} n^{12}\Delta(\Lambda_\mathfrak{c})/\Delta(\mathbb{Z}[i])\); split primes
    contribute nothing. Proved.
14. **Euclidean irreducibility** (Thm 5): \(P^{(2)}_n\), \(P^{(6)}_n\) irreducible for
    **every** \(n\) (cocycle-free square-level Shimura translation + archimedean dominance
    of the principal class). Proved. First-power \(P^{(2)}_n(x^2)\): certified \(n\le16\), open in general.
15. **Product laws**: Gram rigidity, CM-circle law, phase cocycle \(KQ^2\). Proved.
16. **Middle-kernel collapse** onto the modular surface; rank-one pole \(\tfrac6s f_s\otimes f_s\);
    no-go for exact factorization (= existence of cusp forms). Proved.
17. **The \(3/\pi\) split**: Ford stratum has density \(3/\pi\) among all classes
    (deep/Apollonian stratum \(1-3/\pi\)). Derivation solid, one equidistribution step
    not yet rigorous; numerics to \(X=1601\).
18. **Spherical census** \(N_{\mathrm{sph}}(\ell) = \tfrac13 r_3(\ell^2+1) = 4H(4(\ell^2+1))\);
    free \(V_4\)-action; resultant identity (phase never vanishes). Proved. Shape-class
    bijection onto **all** reduced forms of disc \(-4(\ell^2+1)\): verified exactly
    \(\ell\le20\), proof (Gauss-quaternionic) outstanding.
19. **Spherical phase laws**: cap-swap \(u^2_- = \varepsilon_\ell^4 u^2_+\), mirror, level norm
    \(u^2\sigma(u^2)=m_\ell\) with valuation law \(v_p(m_\ell) = 2v_p(4\ell^2)+\tfrac43v_p(H_0)+v_p(H_{1728})\),
    and the half-orbit phenomenon with unit cocycle. Certified at computed levels
    (\(\ell\le6\)); proofs open (the hyperbolic \(\xi\)-torsion/Kronecker step provably fails here).

## Status ledger (open problems, deduplicated)

- **Kronecker limit formula for the phases: settled** (outlook §2.3 / euclidean §6.6,
  [phase-kronecker-limit.md](phase-kronecker-limit.md)): \(S(\chi) = \sum\chi\log|u|\)
  evaluates through \(L'(0,\chi)\) in both aspects (proved); genus characters in full
  closed form (Euclidean field \(=\mathbb{Q}(\sqrt n)\)); remaining open there:
  the exponent law of the \(\Sigma\)-dressings (genus-refined Gross–Zagier), the
  local Euler factor at square conductors, Stark closed forms for order \(>2\)
  characters, and the Robert-unit index of \(\{G_\mathfrak{c}\}\).
- **Hyperbolic first-power loose ends**: \(T = 1\) (distinctness of \(u_f^{12}\)) for general
  \(n\) — the only gap between per-level and all-\(n\) irreducibility. (The Euclidean
  analogue was closed by archimedean dominance; try the same at level \(1-n^2\).)
- **Exact denominator/\(\lambda_n\) valuations** (hyperbolic §5.7, Euclidean §6.1) —
  Gross–Zagier multiplicities / quasi-canonical lifting calculus.
- **Sign law of \(u_f\) on ambiguous classes** (outlook §1.1; the Euclidean R/I criterion
  is proved — hyperbolic analogue open). Data: \(n=11,15\) anomalous.
- **Imprimitive strata** of the class formula at the ramified prime 2 (verified, unproved).
- **\(3/\pi\) split**: make the incomplete-period equidistribution step rigorous.
- **Euclidean first power**: irreducibility of \(P^{(2)}_n(x^2)\) for all \(n\)
  (\(\pm\)-cocycle, weight-2 refinement of Lemma T).
- **Even levels / \(i\mathcal{S}\)**: whole \(\alpha\)-theory run only at odd \(n\);
  transport expected clean (odd discriminants).
- **\(\alpha=2\) regularization** (\(j'(\rho)=0\)).
- **Spherical opens**: the shape bijection unconditionally (Gauss-quaternionic count);
  cap-swap law and \(T\)-norm lemma proofs; the \(\ell=4\) coefficient field; the cocycle
  calculus (what replaces Kronecker's theorem when the cocycle lives in a real field);
  the geometric carrier of the virtual partners \(m_\ell/u^2\).
- The **large programs**: see outlook.md §3, spectral-outlook.md, and the per-document
  outlooks (trace formula, ortholength spectrum, Kronecker limit formula for
  \(\log|u|\), heights/unlikely intersections, other fields, transfer operators,
  bijective Kronecker–Hurwitz).

## Scripts and verification

Environment: python3 with `mpmath` and `sympy` (a venv was used in the sessions;
`matplotlib` needed only for figure scripts). Everything runs standalone from the repo
root, e.g. `python3 scripts/alpha_circles.py --selftest`.

| script | verifies | doc |
|---|---|---|
| `verify_classification.py` | orbit BFS = congruence classes; \(N_e\) closed form; \(X^2/2G\) | circle-classification, euclidean-counting |
| `alpha_circles.py` | level algorithm; \(3H(n^2-1)\), \(n\le40\) (+ spot 301, 1000); figures | hyperbolic-counting |
| `hurwitz_sum.py` | \(\operatorname{vol} = G/3\); \(\sum H(n^2-1)\sim \pi X^2/12G\) | spectral-outlook |
| `involution_experiments.py`, `involution_classmap.py` | \(\sigma\)-identities; class map for all classes, odd \(n\le41\); 3 twisted classes | involution |
| `proof_check.py` | every lemma of the class-formula proof (169 classes) | class-formula-proof |
| `composition_check.py`, `matrix_composition_check.py` | composition/inversion recipes (1913 pairs); matrix-level recipes | circle-composition |
| `moduli_invariants.py` | \(\Theta\) invariance, laws, Hilbert polys, isogenies, traces | moduli-invariants |
| `uf_integer_polynomial.py`, `uf_irreducibility.py` | certified \(Q_n,\Psi_n\) (\(n\le17\)); double-route irreducibility; GZ tagging | moduli-invariants §5.9–5.10 |
| `gz_denominators.py` | GZ-prime identification of denominators | moduli-invariants §5.7 |
| `omega.py`, `omega_verify.py` | monoid membership, factorization, atoms (41 checks) | half-plane-monoid |
| `atom_invariants.py`, `apollonian_chain_demo.py` | atom tables; chain expansion demo | half-plane-monoid |
| `atomic_census.py` | cells, addresses, census, twist/reversal, densities (4414 checks) | atomic-census |
| `product_cocycle.py` | Gram rigidity, cocycle (1205 checks) | product-cocycle |
| `middle_kernel.py` | master formula, layer integral, counting law (654 checks) | middle-kernel |
| `euclidean_moduli_invariants.py` | Euclidean structure, ring class polys, phases, \(P^{(2)},P^{(6)}\) | euclidean-moduli-invariants |
| `mass_law_and_irreducibility.py` | Theorem 4 (M1–M4), Theorem 5 ingredients, exact factorizations | euclidean-moduli-invariants |
| `spherical_moduli_invariants.py` | census, shape bijection (\(\ell\le20\) exact), resultant identity (symbolic), phase laws, level norms | spherical-moduli-invariants |
| `phase_klf.py` | Kronecker-limit character sums: master identities, Epstein \(L'(0,\chi)\), genus factorizations with conductor Euler factors, exact \(\mathbb{Q}(\sqrt d)\) coset factorizations, \(D_n\)/\(R\)-polynomials, safe-PSLQ non-fits (`--selftest`: 250–400 digits) | phase-kronecker-limit |
| `make_composition_figure.py`, `make_omega_figure.py` | figures | — |

**Known gap**: `scripts/first_power_descent.py` is referenced by
[first-power-descent.md](first-power-descent.md) §6 but was **never committed** — the
document's results are independently confirmed by `uf_integer_polynomial.py` /
`uf_irreducibility.py` agreement, but the script should be recovered or rewritten
before publication.

### Certification guard rails (learned the hard way — do not regress)

1. **PSLQ safety**: never trust a multi-term fit when
   (number of terms) × (coefficient digits) approaches the working precision. All
   surviving claims use safe parameters, and the key ones were re-derived by exact
   arithmetic.
2. **Never let a module import set `mp.dps`** — an import-time assignment once silently
   capped "high-precision" runs at ~54 digits. Set precision after imports; compute
   \(j'\) exactly via theta constants / \(-2\pi i E_4^2E_6/\Delta\).
3. **Absolute-error criterion**: accept an integer/rational only with
   \(\ge\max(20,\mathrm{dps}/5)\) spare digits in absolute error (a huge value with an
   \(O(1)\) fractional part is rejected regardless of relative precision).
4. Prefer exact integer/rational arithmetic (Hermite normal forms, `fractions`,
   resultants) over numerics wherever the statement is finite.

## Conventions

- \(\Gamma = \mathrm{SL}_2(\mathbb{Z})\) (or \(\mathrm{PSL}_2\)); the Bianchi group is
  \(\mathrm{SL}_2(\mathbb{Z}[i])\) (or \(\mathrm{PSL}_2\), stated locally). Watch the
  \(\mathrm{PSL}\) vs \(\mathrm{PGL}\) distinction — it is real (odd vs even levels,
  \(\mathcal{S}\) vs \(i\mathcal{S}\)).
- Circles are Hermitian matrices \(M = \binom{A\ B}{\bar B\ C}\), \(\det = -1\);
  \(M_{g\mathcal{C}} = (g^{-1})^\dagger M g^{-1}\); \(\hat{\mathbb{R}} \leftrightarrow M_0 = \binom{0\ i}{-i\ 0}\).
- Two normalizations of the class formula coexist: exterior orientation + reflection
  gives \([\mathfrak{r}_n][f]^{-1}\) (class-formula-proof.md); the inner-disk
  normalization of \(\Omega\) gives the pure twist \([\mathfrak{r}_nf]\)
  (atomic-census.md Thm 6). They differ by exactly one mirror. Decisive test cases:
  \(n=9\) (\(\mathbb{Z}/4\)) and \(n=11\) (\((\mathbb{Z}/2)^2\)).
- New documents should follow the house style: statement–proof Markdown with
  machine-verification section, cross-links to sibling documents, and an outlook.

## Next-step prompts

`prompts/` holds ready-to-run session prompts for the highest-value next steps:
[01-kronecker-limit-formula.md](prompts/01-kronecker-limit-formula.md) (\(\log|u|\) vs
\(L'(0,\chi)\) — the likely headline theorem), [02-schmidt-trace-formula.md](prompts/02-schmidt-trace-formula.md)
(the flagship spectral identity, staged), [03-phase-atlas-and-dit.md](prompts/03-phase-atlas-and-dit.md)
(the atlas figure + sign-law data + Duke–Imamoğlu–Tóth novelty check, one dataset).

## Paper planning

The consolidation plan lives in [outlook.md](outlook.md) §4–5 (two-paper split:
*Counting and composing Schmidt circles* and *The phase of a Schmidt circle*) and
[spectral-outlook.md](spectral-outlook.md) §10 (spectral paper ladder, starting with
the ortholength-spectrum note). The Euclidean theory
([euclidean-moduli-invariants.md](euclidean-moduli-invariants.md)) is fully proved and
publishable as a standalone companion or as part of the phase paper. Before drafting:
the literature diligence pass of outlook.md §4 (Stange, Sarnak, Duke–Imamoğlu–Tóth,
Gross–Zagier, Kubert–Lang/Schertz, Mertens/Bringmann–Kane, Vlasenko–Zagier,
Elstrodt–Grunewald–Mennicke, Parkkonen–Paulin).
