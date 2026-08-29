# The Schmidt arrangement as spectral geometry: a research outlook

Companion to [spectral-geometry.md](spectral-geometry.md) (the state of the art, with all citations; its §14 lists the seven verified gaps referenced below as **G1–G7**). This document lays out a program for studying the Gaussian Schmidt arrangement through the spectral geometry of the Picard orbifold \(M = \mathrm{PSL}_2(\mathbb{Z}[i])\backslash\mathbb{H}^3\). As far as three independent literature probes can tell, *no published work connects Schmidt arrangements to Laplace spectra, trace formulas, or Maass forms* — the combination below is unclaimed territory, with Elstrodt–Grunewald–Mennicke (binary Hermitian forms and Eisenstein series) and Parkkonen–Paulin (Hermitian-form counting via common perpendiculars) the two bodies of work any write-up must position itself against.

Format follows [outlook.md](outlook.md): the question, why it should work, the first move.

## 0. The dictionary

Everything the project has proved has a name on the spectral side:

| project object | spectral-geometric object |
|---|---|
| circle \(\omega \in \mathcal{S}\), Hermitian matrix \(M_\omega\) | totally geodesic plane \(P_\omega \subset \mathbb{H}^3\); all of them project onto the immersed modular surface \(Y = \mathrm{PSL}_2(\mathbb{Z})\backslash\mathbb{H}^2 \looparrowright M\), area \(\pi/3\) |
| \(\alpha(\omega) = n\) (inversive product with \(\hat{\mathbb{R}}\)) | \(\cosh d(P_0, P_\omega) = n\): the **distance spectrum** of the plane orbit is \(\{\operatorname{arccosh} n\}\) |
| level-\(n\) classes \(\leftrightarrow\) forms of disc \(1-n^2\); weighted count \(H(n^2-1)\) per copy of \(Y\) | multiplicity of the distance \(\operatorname{arccosh} n\) in the **ortholength-type spectrum** of \((Y, Y)\) in \(M\) |
| Cartan image \(Y_X = X\bar X^{-1}\), \(\operatorname{tr} = -2\alpha\); \(-Y_X\) loxodromic of length \(2\operatorname{arccosh} n\) | the **common perpendicular** geodesic of \(P_0\) and \(P_\omega\), doubled: a closed geodesic meeting \(Y\) orthogonally |
| \(\sigma(X) = \bar X^{-1}\), Gelfand pair, \(\hat\sigma[f] = [\mathfrak{r}_n][f]^{-1}\) | the involution behind the Jacquet–Lai relative trace formula: its geometric side is parametrized by \(g \mapsto g\bar g^{-1}\) (Hilbert 90) — *our Cartan embedding, adelically* |
| \(N_e\)-counting constant \(1/2G\) | \(1/(6\operatorname{vol} M)\), by Humbert: \(\operatorname{vol}(M) = G/3\) |
| the \(\Omega\)-monoid's atoms (Apollonian gasket) | a thin subgroup with its own (Lax–Phillips) spectral theory: \(\lambda_0 = \delta(2-\delta)\), \(\delta = 1.3056867280\ldots\) |
| infinitude of \(\mathcal{S}\) | arithmeticity of \(M\) (Bader–Fisher–Miller–Stover): the arrangement is a *visual certificate* of arithmeticity |

**Two numerical anchors** (session of 2026-08-28, [scripts/hurwitz_sum.py](scripts/hurwitz_sum.py)):

1. \(\operatorname{vol}(M) = \tfrac{8\zeta_K(2)}{4\pi^2} = \tfrac{G}{3} = 0.305321\ldots\) exactly.
2. The level census satisfies, numerically to three digits up to \(X = 1200\) with residuals of size only \(O(X)\):
$$
\sum_{n \le X} H(n^2 - 1) \;\sim\; \frac{\pi}{12\,G}\,X^2
\;=\; \boxed{\ \frac{\operatorname{area}(Y)^2}{4\pi\,\operatorname{vol}(M)}\;X^2\ }
$$
The constant has the exact volume-ratio form of an equidistribution/relative-trace main term — the \(\operatorname{area}^2/\operatorname{vol}\) is the contribution of the constant eigenfunction \(\phi_0 = 1/\sqrt{\operatorname{vol}}\) through its period \(P_Y(\phi_0) = \operatorname{area}(Y)/\sqrt{\operatorname{vol}}\), and \(1/4\pi\) is a universal plane-pair transform factor. *Derivation of \(\pi/12G\):* at level \(n\), the number of circles with half-curvature \(q\) is the number of solutions of \(x^2 \equiv 1-n^2 \pmod{4q}\); averaged over \(n\), this has mean \(N_e(q)/q\), whose mean over \(q\) is \(1/G\) ([euclidean-counting.md](euclidean-counting.md) §3); the ideal-triangle truncation profile integrates to \(\int_0^\infty \varphi(u)\,du = \pi/2\) where \(\varphi(u) = 1 - \sqrt{\max(0, 1 - u^{-2})}\). This is provable elementarily with a modest error term; **the spectral program's first quantitative goal is the error term beyond it** (§1–2).

## 1. The relative trace formula for \((\mathrm{PSL}_2(\mathbb{Z}[i]), \mathrm{PSL}_2(\mathbb{Z}))\) — the flagship

**The question.** Establish the exact identity ("Schmidt trace formula"): for a point-pair invariant \(k\) on \(\mathbb{H}^3\) with transform \(h\),
$$
\underbrace{\sum_j |P_Y(\phi_j)|^2\, h(t_j) \;+\; (\text{Eisenstein})}_{\text{spectral side}}
\;=\;
\underbrace{(\text{2D trace formula of } Y)\;+\;\sum_{n \ge 2} c\,H(n^2-1)\, F_h(\operatorname{arccosh} n) \;+\; (\alpha = 1\ \text{regularized})}_{\text{geometric side}}
$$
obtained by integrating the automorphic kernel \(K(x,y) = \sum_{\gamma\in\Gamma}k(x, \gamma y)\) over \(Y \times Y\) and decomposing \(\Gamma\) into \(H\)-double cosets: \(P_Y(\phi) = \int_Y \phi\), \(F_h\) the plane-pair orbital transform, \(c\) an explicit bookkeeping constant (orientations, the two \(\mathrm{PSL}/\mathrm{PGL}\) strata). The identity double coset \(HeH = H\) reproduces the full Selberg trace formula of the modular surface with the restricted 3D kernel; the level-\(n\) cosets contribute through the teammate-proved multiplicity \(H(n^2-1)\); the Ford stratum \(\alpha = 1\) (tangent planes) is a parabolic-type term to be regularized together with the continuous spectrum.

**Why it should work.**
- Every ingredient exists (survey §9, §13): Jacquet–Lai is this identity adelically; Lapid–Rogawski/Chaudouard give the Eisenstein-side technology; Zagier's renormalization handles the divergent \(Y\)-periods of Eisenstein series; the complete 2D template — double cosets of two hyperbolic subgroups, with error terms — is Lekkas–Petridis (2025), and the period-spectrum \(\leftrightarrow\) ortholength-spectrum principle is Martin–McKee–Wambach (2011), both for surfaces.
- **Our geometric side is exact integer arithmetic.** In every 2D instance the ortholengths are transcendental data; here the distance multiplicities are Hurwitz class numbers and the distances are \(\operatorname{arccosh}(n)\), \(n \in \mathbb{Z}\). This is the rare relative trace formula whose geometric side is completely explicit — the analogue of what makes the Eichler–Selberg formula (where \(H(4m-t^2)\) plays this role) more useful than the bare Selberg formula.
- The constant-eigenfunction term already *matches the data* (§0, anchor 2).
- The spectral side is sparse and structured: \(P_Y(\phi) \neq 0\) forces \(\phi\) **distinguished**, and by Flicker + Flicker–Zinoviev (valid for \(\mathbb{Q}(i)/\mathbb{Q}\)) the distinguished cuspidal spectrum consists of base-change lifts with central character \(\chi_{-4}\) — a density-zero subfamily. The error term in the census is governed by a thin, explicitly describable family of \(L\)-values.

**Deliverables.** (a) The exact identity; (b) \(\sum_{n\le X} H(n^2-1) = \tfrac{\pi}{12G}X^2 + O(X^{?})\) with a power-saving spectral exponent (the observed residuals \(O(X)\) suggest the truth is \(X^{1+\varepsilon}\); by analogy with Good's \(X^{2/3}\) in 2D, a first target is \(X^{4/3+\varepsilon}\), with the Eisenstein contribution the likely bottleneck); (c) the **relative Weyl law** \(\sum_{t_j \le T}|P_Y(\phi_j)|^2 \sim c' T^{?}\) — the quantitative measure of how distinguished the Picard spectrum is, the 3D analogue of Martin–McKee–Wambach's period asymptotics; (d) as a by-product, the first exact **ortholength spectrum** in dimension 3 (§5).

**First moves.**
1. *The Eisenstein period in closed form* (gap **G3**, self-contained lemma): compute \(P_Y^{\mathrm{reg}}(E(\cdot, s))\) by unfolding over \(H\backslash\Gamma/\Gamma_\infty\) with Zagier regularization. Sanity anchors: the residue at \(s = 1\) must be \(\operatorname{area}(Y)/(2\operatorname{vol} M)\)-type (Jacquet–Lai's residue mechanism), and the answer should be the ratio of \(\zeta\)-functions that EGM's Hermitian-form special values (Math. Ann. 277 (1987)) and Flórez–Karabulut–Vu predict.
2. *Derive the compact-support identity.* With \(k\) of small support, only finitely many levels contribute; the derivation is the classical pre-trace argument plus the incidence bookkeeping already proved in [hyperbolic-counting.md](hyperbolic-counting.md). Machine-check it: geometric side from `alpha_circles.py`, spectral side deferred (needs periods), but internal consistency (small-support positivity, the \(h \to\) counting limit against `hurwitz_sum.py`) is checkable now.
3. *Port Lekkas–Petridis.* Their double-coset decomposition, orbital transforms, and large-sieve error management are the closest blueprint; the new technical content is (i) the plane-pair spherical transform on \(\mathbb{H}^3\) (Fock coordinates; elementary), (ii) the cusp of \(Y\) inside the cusp of \(M\) (the \(\alpha = 1\) Ford stratum ↔ the unipotent double cosets), (iii) the \(\mathrm{SL}_2\)-packet subtleties (Anandavardhanan–Prasad): work classically with the four symmetry classes \(D, G, C, H\) of Bianchi forms rather than adelically — predict and then verify which classes can carry a nonzero \(Y\)-period.

**Risks, honestly.** The \(\mathrm{SL}_2\) period is not always factorizable and can vanish identically on abstractly-distinguished representations (A–P 2006) — the classical (symmetry-class) formulation sidesteps interpretation but the arithmetic of *which* forms appear must be settled by the local analysis at \(1+i\) and \(\infty\); the regularized Eisenstein × Eisenstein inner products are delicate (this is where Zagier-type renormalization or truncation à la Arthur enters); the main term must be assembled from *both* the residual and continuous spectra — the numerics of §0 constrain but do not prove the bookkeeping.

## 2. Effective circle counting: the plane-pair and point-plane lattice problems

**The question(s).** (a) For a nice region \(E \subset \mathbb{C}\): \(\#\{\omega \in \mathcal{S} : \operatorname{curv} \le T,\ \text{center} \in E\} = \tfrac{T^2}{8G}\operatorname{Leb}(E) + O_E(T^{2-\eta})\) with explicit \(\eta > 0\). (b) The point-to-plane problem for the Picard group: \(\#\{\gamma \in \Gamma : \cosh d(P, \gamma P_0) \le X\}\) — extend Laaksonen's cocompact result to the cofinite case. (c) The plane-pair problem — this is §1(b) again, seen as lattice-point counting on \(H\backslash G/H\).

**Why it should work.** The soft main term is known in two independent ways — Oh–Shah's lattice case (\(cT^2\) with Lebesgue equidistribution of circle centers; mixing, no error term) and the teammate's *exact* multiplicative formula in the translation-averaged regime — and the missing effective input is exactly the spectral theory surveyed: \(\lambda_1(M) \approx 44.85\) is enormous, the Eisenstein series is explicit through \(\zeta_{\mathbb{Q}(i)}\), and the modern Kuznetsov toolbox over \(\mathbb{Z}[i]\) (Bruggeman–Motohashi; Qi's symmetric-square large sieve) is state-of-the-art. The de Sitter picture makes (a) a counting problem for the three \(\Gamma\)-orbits of integral points on the Hermitian determinant-\((-1)\) quadric (the three twisted classes of [involution.md](involution.md) §3!) — a Duke–Rudnick–Sarnak/Eskin–McMullen affine symmetric variety, never effectivized in this case. Parkkonen–Paulin's exponential-mixing error terms (BAPP monograph) give a fallback exponent; the spectral route should beat it and make it explicit.

**First move.** Do (b) first — it is a clean, self-contained paper ("the hyperbolic lattice-point problem relative to a geodesic plane for the Picard group"): spectral expansion of \(\sum_\gamma k(d(P, \gamma P_0))\) against the *period-integrated* automorphic kernel; the new arithmetic input is the expansion of \(P_{P_0}\)-type integrals of Eisenstein series (the same Lemma as §1 first-move 1). Then (a) by summing over the orbit structure, and (c) via §1.

## 3. Class-number identities: Eichler–Selberg over \(\mathbb{Z}[i]\) and the modularity of the census

**The question.** Three nested targets. (i) Prove the Kronecker–Hurwitz relation \(\sum_t H(4m - t^2) = \sum_{d|m}\max(d, m/d)\) bijectively on the arrangement ([outlook.md](outlook.md) 3.3) — now with the added interpretation that both sides are geometric-side terms of §1-type formulas. (ii) Decide whether \(\mathcal{H}_{\mathrm{Sch}}(\tau) := \sum_{n\ge 2} H(n^2-1)\,q^n\) (and its refinement by symmetry classes and by \(\hat\sigma\)-orbits) is a mock/quasi-modular object — the "shifted-square slice" program of [outlook.md](outlook.md) 3.2, which the Kudla–Millson philosophy upgrades to a precise prediction: pairing the theta lift for \(\mathrm{O}(3,1)\) against the \(Y\)-cycle should produce a modular generating series whose coefficients are exactly the level counts (verified absent as a publication, **G6**-adjacent). (iii) The full Bianchi Eichler–Selberg formula: traces of Hecke operators \(T_\mathfrak{n}\) on Bianchi forms expanded in class numbers of binary Hermitian forms (**G6**), with the Schmidt census as its geometric model — the 3D sibling of what Zagier's appendix does over \(\mathbb{Q}\), and of the quartic-CM-field Hilbert version of Kuga–Seymour-Howell–Wakatsuki (2026).

**Why it should work.** (ii) is a finite experiment: Mertens/Bringmann–Kane technology handles exactly such \(H\)-sums with quadratic constraints, and Walker (ANT 2025) shows the correlation machinery is current. (iii) has all local ingredients in EGM (Math. Ann. 277) + Flórez–Karabulut–Vu, and the trace-formula side in EGM98.

**First move.** PSLQ/Sturm-style test of \(\mathcal{H}_{\mathrm{Sch}}\) against low-weight bases on \(\Gamma_0(4)\)-type groups (one afternoon with the census data). Whatever the outcome informs (i) and (iii).

## 4. Periods of Maass forms over the Schmidt planes; the sup-norm anomaly

**The questions.**
1. Pin down the distinguished spectrum at level 1: which of the four symmetry classes \(D, G, C, H\) support \(P_Y \neq 0\), and the precise base-change description (lifts with central character \(\chi_{-4}\); the local condition at the ramified prime \(1+i\) is where \(\hat\sigma\)'s twist \(\mathfrak{r}_n\) should reappear as an Atkin–Lehner-type sign — conjectural, testable).
2. Prove period bounds/asymptotics for \(P_Y(\phi_j)\) on the *non-compact* Picard orbifold (Hou's amplification results are compact-only), and the restricted-QUE statement: does \(|\phi_j|^2\big|_Y\,\) equidistribute on \(Y\)?
3. Interpret the sup-norm anomaly through the arrangement: the \(\lambda^{1/4}\)-power growth of eigenfunctions on arithmetic 3-manifolds (Rudnick–Sarnak, Milićević) exists *exactly* on manifolds with geodesic surfaces — BHM call the coincidence "extremely intriguing". The Schmidt arrangement is the organizing skeleton of those special points; mapping where Then's large eigenfunction values sit relative to the circles/planes of \(\mathcal{S}\) is an unexplored and cheap experiment that could turn "intriguing coincidence" into geometry.

**Why it should work.** For (1)–(2) the analytic tools (amplifier over \(\mathbb{Z}[i]\), Kuznetsov, the explicit Eisenstein theory) all exist; the modular surface is the simplest possible cycle. For (3) the data (13,950 eigenfunctions) already exists in the quantum-chaos literature.

**First move.** Numerical periods: recompute the first \(\sim\)20 Hecke–Maass forms on \(M\) (Hejhal-style over \(\mathbb{Z}[i]\), or obtain Then's data) and integrate over \(Y\); compare with the base-change predictions. This simultaneously seeds §1's spectral side.

## 5. The ortholength spectrum, and reciprocal geodesics in dimension 3

**The question.** Package what is already proved as: *the ortholength spectrum of the modular surface in the Picard orbifold is \(\{2\operatorname{arccosh} n : n \ge 2\}\) with multiplicities given by \(H(n^2-1)\)-data* (with the precise stabilizer weights of [hyperbolic-counting.md](hyperbolic-counting.md) and the \(\hat\sigma\)-pairing of [involution.md](involution.md) describing the unoriented count). To our knowledge this would be the **first exactly computed ortholength spectrum of a geodesic surface in any hyperbolic 3-manifold** — the objects Martin–McKee–Wambach's RTF and Parkkonen–Paulin's counting theory quantify abstractly. Then develop the 3D theory of **reciprocal geodesics** (verified gap **G5**): classify and count loxodromic classes of \(\mathrm{PSL}_2(\mathbb{Z}[i])\) conjugate to their inverses (or to their Galois conjugates — several flavors exist in 3D), in the style of Sarnak's Clay paper, where the counts were governed by genus theory; here the \(\sigma\)-twisted classes with \(\operatorname{tr} = -2n\) are counted by \(\mathrm{Cl}(1-n^2)/\)-data with the \(\hat\sigma\)-fixed classes \([f]^2 = [\mathfrak{r}_n]\) as the distinguished stratum.

**Why it should work.** The class-formula machinery is *done*; what remains is translation (twisted conjugacy vs. ordinary conjugacy bookkeeping — the \(\hat\sigma\)-orbit invariance already proved in [involution.md](involution.md) §4 is the key step) plus positioning against Sarnak/Erlandsson–Souto/Parkkonen–Paulin (strongly reversible geodesics, JIMJ 2026 — their "\(G\)-reciprocal Hermitian forms" from 2011 is the same idea one level down).

**First move.** A short standalone note: statement, proof via the Cartan bridge, comparison table with Sarnak's 2D counts (\(\Pi_{\langle\phi_R\rangle}(x)\sim \tfrac38 x\) etc. vs. our \(\sum H(n^2-1) \sim \tfrac{\pi}{12G}X^2\)). This is the cheapest publishable unit in the whole program.

## 6. The phase \(u_f\) as a spectral object (speculative)

The proved law \(\Theta_f\,\Theta_{\mathfrak{r}f} = e^{-\ell}\) says the phase is a canonical *half-multiplier* along the closed geodesic attached to the circle — cocycle data over the geodesic flow. Three directions:

1. **Twisted relative trace formula**: insert weights \(u_f^k\) (or the characters \(\chi\) of \(\mathrm{Cl}(1-n^2)\)) into §1's geometric side. The character-averaged geometric side is precisely the teammate's planned \(\sum_f \chi(f)\log|u_f|\)-type sums ([outlook.md](outlook.md) 2.3, Kronecker limit formula); the RTF would provide the *dual, spectral* evaluation — first-derivative Eisenstein periods, in the territory of Duke–Imamoğlu–Tóth ([outlook.md](outlook.md) 2.6) and of the 2024 Herrero–Imamoğlu–von Pippich–Schwagenscheidt "twisted traces on \(\mathbb{H}^3\)" (the closest published relative).
2. **Equidistribution of \(\arg u_f\)** ([outlook.md](outlook.md) 3.4) as a *sparse-family Duke problem*: the CM points of discriminant \(1-n^2\) at level \(n\) form a sparse (discriminant \(\asymp\) square) family — precisely the frontier of current equidistribution technology. Flag: hard; but even conditional/partial results would be novel.
3. **Holonomy refinement of §5**: the complex length of the doubled perpendicular is real (\(\theta = 0\)); \(u_f\) looks like the natural "next coefficient" — a torsion-like secondary invariant. Whether \(u_f\) appears in a geometric expansion (linking numbers, eta-invariants, Frahm–Spilioti-type twisted Ruelle zeta values) is wide open and worth one exploratory afternoon against the DIT dictionary.

## 7. Transfer operators: a Mayer theorem for the gasket and for the arrangement

**The question.** (a) Express the Selberg zeta function of the Apollonian/\(\Omega\)-side (the Whitehead-adjacent, infinite-volume manifold) as a Fredholm determinant of the **atom transfer operator**
$$
(\mathcal{L}_s f)(z) \;=\; \sum_{A\ \mathrm{atom}} |A'(z)|^{s}\, f(A z),
$$
using the unique-factorization theorem of [half-plane-monoid.md](half-plane-monoid.md) — rigidity of \(\Omega\) is exactly the "no overcounting" that makes symbolic codings clean. (b) The full Picard version via A. Schmidt's continued fractions: Nakada (1988) already realized the natural extension as a cross-section of the geodesic flow on \(M\); what is missing (**G4**, confirmed absent for *every* Bianchi group) is the Mayer-type determinant identity. The 2D blueprint — including the treatment of cusps by "cuspidal acceleration" — is the Pohl–Wabnitz Memoir (2026).

**Why it should work.** The alphabet is arithmetically meaningful (atoms are classified by \((\alpha, [f])\)), the coding is proven rigid, the target constant \(\delta = 1.3056867280\ldots\) is known to 128 rigorous digits (Vytnova–Wormell) for validation, and \(\lambda_0 = \delta(2-\delta)\) links the leading zero to Kontorovich–Oh's base eigenfunction. Payoffs: resonances of the gasket manifold; Lewis–Zagier-type *period functions over \(\mathbb{Z}[i]\)* (new objects); a new algorithm for \(\delta\); and the natural home for the \(\alpha\)-spectrum question below.

**Adjacent question.** The \(\alpha\)-spectrum of the atoms ([half-plane-monoid.md](half-plane-monoid.md) Q1) is a local–global problem for a *semigroup* orbit — exactly the setting where Rickards–Stange (Duke 2025) found reciprocity obstructions after HKRS. First move: compute the spectrum to \(10^6\) and test for missing quadratic families \(u\cdot k^2\) before conjecturing density-one.

## 8. Inverse spectral questions: hearing the arrangement

**The question.** The Laplace spectrum of an arithmetic hyperbolic 3-manifold determines its commensurability class (Reid; Chinburg–Hamilton–Long–Reid via the length *set*). Does the **plane-distance census** — the map \(n \mapsto\) (weighted count of level-\(n\) classes), i.e. \(H(n^2-1)\) for \(\mathbb{Q}(i)\) — likewise determine the Bianchi orbifold among Bianchi orbifolds? Concretely: compute the level function of \(\mathcal{S}_K\) for other imaginary quadratic \(K\) (Stange's framework; over \(K\) the invariant will be \(\alpha \in \tfrac{1}{\sqrt{|d_K|}}\mathbb{Z}\)-normalized) and identify the class-number family that replaces \(H(n^2-1)\) — this is also the natural first step of [outlook.md](outlook.md) 3.6. A "can one hear the Schmidt arrangement?" theorem would tie the project to the isospectrality literature (Vignéras orbifolds, Bartel–Page) from a genuinely new side: the *relative* (surface-anchored) spectrum instead of the absolute one.

## 9. Numerics roadmap (cheap, high information)

1. **Census asymptotics**: extend [scripts/hurwitz_sum.py](scripts/hurwitz_sum.py) to \(X = 10^4\) via the class-number formula (not form enumeration), fit the second-order term, and test the error exponent \(X^{1}\) vs \(X^{4/3}\) hypothesis of §1. PSLQ the \(X\)-coefficient against \(\{\pi/G,\ \log,\ \zeta'(2)/\zeta(2)\}\)-type constants (PSLQ safety rules of [outlook.md](outlook.md) §4 apply).
2. **Modularity test** of \(\sum H(n^2-1)q^n\) and its symmetry-class refinements (§3(ii)).
3. **Maass data**: implement/obtain the first \(\sim\)20 Picard Hecke–Maass forms; compute \(Y\)-periods; test the distinction predictions (§4.1); map sup-norm peaks against the arrangement (§4.3).
4. **Eisenstein period**: verify the closed form of \(P^{\mathrm{reg}}_Y(E(\cdot,s))\) numerically once derived (anchor: residue at \(s=1\)).
5. **Atom transfer operator**: finite-rank Fredholm determinant numerics; compare against \(\delta\) to many digits (§7).
6. **Phase atlas** ([outlook.md](outlook.md) 1.6) extended by spectral overlays: circles colored by \(\arg u_f\) *and* by the value of low Maass forms at the hyperbolic center — one plot, possibly two discoveries.

## 10. Recommended order and paper map

1. **(§5)** The ortholength-spectrum note — nearly free, stakes the claim, introduces the dictionary to the spectral community.
2. **(§1 first moves 1–2 + §2(b))** The Eisenstein-period lemma and the point-to-plane Picard lattice count — one coherent medium paper; unlocks everything quantitative.
3. **(§1 full + §2(a))** The Schmidt trace formula and effective circle counting — the flagship ("Spectral geometry of the Schmidt arrangement", natural paper III alongside the two planned in [outlook.md](outlook.md) §4).
4. **(§3(ii) + §9.2)** The modularity experiment — an afternoon; a positive result redirects §3 to Kudla–Millson.
5. **(§7)** The gasket Mayer theorem — the deep dynamical prize, independent of 1–4.
6. **(§4, §6, §8)** The long arc: periods/QUE, the phase spectrally, inverse problems.

Items 1–2 are low-risk and self-contained; item 3 is where the project's exact arithmetic (multiplicities \(H(n^2-1)\), the involution \(\hat\sigma\), the \(\mathfrak{r}_n\)-twist) should produce theorems that the general theory (Oh–Shah, Parkkonen–Paulin, EGM) provably cannot see — none of those frameworks knows the level stratification exists.
