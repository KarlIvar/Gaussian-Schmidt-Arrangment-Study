# Outlook: where the Gaussian Schmidt arrangement project can go next

## 0. Where we stand

Proved and machine-certified across the documents:

1. **Diophantine classification** of the circles of \(\mathcal{S} = \mathrm{PSL}_2(\mathbb{Z}[i])\hat{\mathbb{R}}\) ([circle-classification.md](circle-classification.md)); Euclidean counting \(N_e(n) = n\prod(1-\chi_{-4}(p)/p)\) with Catalan-constant asymptotics ([euclidean-counting.md](euclidean-counting.md)).
2. **Hyperbolic counting**: circles at level \(\alpha = n\) are hyperbolic circles of radius \(\operatorname{arcoth} n\) at CM points of discriminant \(1-n^2\); weighted count \(3H(n^2-1)\) in the ideal triangle ([hyperbolic-counting.md](hyperbolic-counting.md)).
3. **The involution** \(\sigma(X) = \bar X^{-1}\): class formula \(\hat\sigma[f] = [\mathfrak{r}_n][f]^{-1}\), proved via the explicit unitary basis ([class-formula-proof.md](class-formula-proof.md)); circle-level class group operations ([circle-composition.md](circle-composition.md)).
4. **The phase invariant** \(\Theta\), the sixth coordinate on double cosets; functional equations (mirror, inversion = \(e^{-\ell}\)); the normalized \(u_f = \varepsilon\Theta_f\); full generalized-dihedral Galois equivariance of \(u^{12}\) via Shimura–Siegel reciprocity, with the \(\xi\)-torsion closed by the cocycle-ideal lemma; explicit fields for the pair-sums via genus characters; denominators = Gross–Zagier primes of \((1-n^2, -3)\) and \((1-n^2, -4)\) ([moduli-invariants.md](moduli-invariants.md)); and, as a corollary of the equivariance, the phase units of a level are the roots of **one polynomial with integer coefficients** (§5.9 there: unconditionally in degree \(12h\), certified in the sharp degree \(h\) for odd \(n \le 17\)), and that polynomial is **irreducible** at every computed level — one level, one Galois orbit, the degree-\(h\) polynomial being the common minimal polynomial of the level's phase units (§5.10 there).
5. **The Euclidean moduli theory** (translations \(\times\ \mathrm{SL}_2(\mathbb{Z})\) instead of \(\mathrm{SL}_2(\mathbb{Z}) \times \mathrm{SL}_2(\mathbb{Z})\)): disks of curvature \(2n\) mod translation \(\leftrightarrow\) \(\mathrm{Cl}(\mathbb{Z} + n\mathbb{Z}[i])\) (so \(N_e(n) = 2h(-4n^2)\)); \(j\)-values = ring class polynomial squared; trace slice \(t(4n^2)\); the Euclidean sixth invariant \(\Theta = j'(w)/c^2\) and its lemniscatic phase units with monic integer level polynomials; the \(\Delta\)-mass law **proved** (Theorem 4 there: constancy of the full Hecke mass + PID stratification + generating-function identities + the conjugation sign analysis), and the level polynomials \(P^{(2)}_n, P^{(6)}_n\) **proved irreducible for every \(n\)** (Theorem 5 there: the root-free closed form \(u^2 = -12\beta(\beta-1728)\,g_2(\Lambda)/g_2(\mathbb{Z}[i])\), a cocycle-free square-level Shimura translation law, and an archimedean dominance lemma killing the coincidence subgroup — so the square-level analogues of both hyperbolic open hypotheses, first-power descent and \(T = 1\), are closed on the Euclidean side); \(H_{-4n^2}\) irreducible classically, and the first-power disk polynomial \(P^{(2)}(x^2)\) irreducible at every computed level ([euclidean-moduli-invariants.md](euclidean-moduli-invariants.md), with its own outlook §6 there).

6. **The Kronecker limit formula for the phases** ([phase-kronecker-limit.md](phase-kronecker-limit.md)):
   \(S(\chi) = \sum\chi\log|u| = -2L'(0,\chi) + \tfrac23\Sigma_0 + \tfrac12\Sigma_{1728}\)
   (Euclidean) and \(S(\chi) = -4L'(0,\chi) + \tfrac43\Sigma_0 + \Sigma_{1728}\) for
   \(\chi(\mathfrak{r}_n) = -1\), \(S \equiv 0\) for even \(\chi\) (hyperbolic — the
   Norm Lemma cancels \(\varepsilon\) and \(\mu\) in \(|u_f|\)); both **proved**.
   Genus characters in closed form \(c\,h(d)\log\varepsilon_d\) (factorizations
   Sturm-proved, conductor Euler corrections included), with the Euclidean genus
   field equal to **\(\mathbb{Q}(\sqrt n)\)**; the \(\Delta\)-data uniformized by
   certified integer polynomials \(D_n\) (constant term \(\pm M(n)\)); the
   hyperbolic \(\mathfrak{r}\)-twisted \(\Delta\)-ratios \(R_f\) certified to be
   **units** (palindromic integer level polynomials with constant term 1); the
   GZ-supported \(j\)-dressing factored exactly in \(\mathbb{Q}(\sqrt d)\);
   certified non-fits for all order-\(>2\) characters (Stark regime).

What follows are next steps, roughly by size. Each: the question, why it should work, the first move.

## 1. Small (days to weeks)

**1.1 The sign of \(u_f\).** Fresh data (session of 2026-08-28): on ambiguous classes, \(\operatorname{sign}(u_f) = -1\) identically for \(n = 3, 5, 7, 9, 13, 17\), but at \(n = 11, 15\) it equals \(-\psi(x)\), \(\psi\) the quadratic character of the \(\langle\mathfrak{r}\rangle\)-coset:

| \(n\) | 3 | 5 | 7 | 9 | 11 | 13 | 15 | 17 |
|---|---|---|---|---|---|---|---|---|
| sign pattern | − | −− | −− | −− | −++− (by coset) | −−−− | −+−+ (by coset) | −−−− |

Note \(n = 11\) and \(13\) have identical class groups \((\mathbb{Z}/2)^2\) but different sign behavior. Candidate explanations: a 2-adic genus character (the crude prime-disc decomposition mishandles conductor-2 parts), or the sign of the Gross–Zagier monomial. First move: recompute with proper conductor-aware characters; extend to \(n \le 41\). A clean sign law would complete the real-quadratic description of \(u\) on ambiguous classes.
*New tool (from 2.3):* `phase_klf.py` computes the exact sign of every genus
coset-product factorization \(A/B = \pm\varepsilon^k\prod(\pi/\pi')^{e}\) in
\(\mathbb{Q}(\sqrt d)\) — the sign of \(u_f\) on ambiguous classes should be
exactly this \(\pm\) traced through the closed form; likely resolvable together
with 2.8.

**1.2 Exact formulas for the counts via Eichler's relations.** \(H(n^2-1)\) is the \(t = \pm1\) term of the classical relation \(\sum_{t^2\le 4m} H(4m-t^2) = \ldots\) Combined with Mertens/Bringmann–Kane-type mock-modular class number identities, one should get exact recurrences for the triangle counts \(3H(n^2-1)\) — closed-form arithmetic for our census. First move: assemble the relation for \(m = (n^2-1)/4 + \) shifts and isolate the diagonal.

**1.3 Exact denominator valuations.** §5.7 identified the support and the squared-vs-first-power pattern; Gross–Zagier gives multiplicities of each \(\frac{e|D|-x^2}{4}\). Trace them through the \((4,3)\)-exponents and the 6th root to prove \(v_p(\operatorname{den} u)\) exactly. Also: nobody has looked at the **numerators** — their primes are not GZ-constrained and may encode height data.

**1.4 Even levels.** The whole \(\alpha\)-theory was run at odd \(n\) (family \(\mathcal{S}\)). Even \(n\) lives in \(i\mathcal{S}\), discriminants \(1-n^2 \equiv 3 \pmod 4\), twist ideals of odd norm \(\frac{n\mp1}{2}\)-adjusted. Everything should transport with cleaner 2-adic behavior (odd discriminant!) — a good consistency laboratory. First move: adapt `build_P` to the \(i\mathcal{S}\)-parity and rerun the pipelines.

**1.5 Regularizing \(\alpha = 2\).** The level \(n = 2\) (the inscribed circle; CM point \(\rho\)) is excluded because \(j'(\rho) = 0\). Both \(m_1, m_2\) are \(\rho\)-points, so the *ratio* of leading Laurent coefficients of any kernel defines a regularized \(\Theta(2)\). Compute it; it should be a particularly clean number (the class group is trivial and \(\varepsilon = 2+\sqrt3\)).

**1.6 A phase atlas.** Render the circles of the ideal triangle colored by \(\arg u\) and \(\log|u|\) across levels — a "phase portrait of the Schmidt arrangement". Cheap, likely to expose patterns (e.g. the sign dichotomy of 1.1 geometrically), and the natural figure for any eventual paper.

**1.7 The \(\zeta_8\)-levels.** \(n^2 - 1 = 2\square\) (\(n = 3, 17, 99, \ldots\) — NSW-adjacent) are exactly where \(B = \mathbb{Q}(\zeta_8)\) and the \(u^6\)-translation may pick up a sign. At \(n = 17\) the data shows all signs \(-\): check whether the metaplectic sign actually occurs, or the \(\mu_8\)-caveat is vacuous.

## 2. Medium (weeks to months)

**2.1 First-power descent — the isogeny-differential route (recommended).** The remaining \(\mu_{12}\)-cocycle would evaporate if \(u\) has a purely algebraic-geometric expression. Candidate: the \(\hat\sigma\)-pair carries two cyclic isogenies \(\phi_\mp: E_2 \to E_1'\)-conjugates of degrees \(\frac{n\mp1}{2}\) whose composite is multiplication by \(\omega\) (norm \(\frac{n^2-1}{4}\), and \(\omega^2 = D/4\)); their pullbacks of Néron differentials give scalars \(\lambda_\mp\) with \(\lambda_-\lambda_+ \sim \omega\). Conjecture: \(u\) is a monomial in \(\varepsilon, \lambda_\mp, \mu\) — testable by PSLQ in an afternoon. If true, \(u\) is Galois-theoretically transparent at first power (isogenies are algebraic), law 3 closes at \(u\)-level with no Siegel-unit bookkeeping, and the invariant gets a moduli-theoretic definition. If false, fall back on:

**2.2 First-power descent — the Kubert–Lang route.** Write the kernel as \(\eta^4\cdot(\text{level-6 weight-0 function})\); \(u\) becomes a Siegel–Ramachandra-type invariant; the known first-power Galois action of Siegel units (Kubert–Lang ch. 11, Schertz) resolves the sixth root. Routine, longer, certain to work.

**2.3 Kronecker limit formula for the phase: \(\log|u|\) and \(L'\)-values — SETTLED**
(see §0.6 and [phase-kronecker-limit.md](phase-kronecker-limit.md)): the clean
identity exists and is proved in both aspects, so \(u\) *is* a geometrically
defined Stark-adjacent unit system, up to the exactly identified GZ dressing.
The answer raises four sharp new questions, distributed below:
(i) *why* the twisted \(\Delta\)-ratios \(R_f\) are units — a per-class
valuation-cancellation statement across the \(\mathfrak{r}\)-twist (→ 2.7/3.4);
(ii) the exponent law of the GZ dressing (→ 2.8, new);
(iii) the conceptual local Euler factor at square/composite conductors
(phase-kronecker-limit §8, opens 2);
(iv) what geometric structure at Euclidean conductor \(n\) knows the Pell unit
of \(\mathbb{Q}(\sqrt n)\) — a real-quadratic shadow in the conductor aspect,
mirroring \(\varepsilon_{n^2-1}\) in the discriminant aspect (→ 3.4).

**2.4 Cusp degeneration: the phase at \(n \to 1\) and Dedekind sums.** As \(n \to 1\), \(\varepsilon \to 1\), the geodesic shrinks into the cusp, and the circles become Ford circles. The degenerate limit of \(u\) (and of the \(\zeta\)-cocycle, which is built from \(\eta\)-multipliers) should be classical Dedekind-sum data. Making this precise would tie the phase theory to the classical modular-transformation literature and likely *explain* the \(\mu_{12}\).

**2.5 The imprimitive strata.** The class formula's behaviour on non-invertible classes (content \(g > 1\), ramified prime 2) is verified but unproved. The tools now available (cocycle ideal, \(\mathcal{O}_B\)-extension trick) should handle it: non-invertibility of \(\mathfrak{a}_f\) obstructs only one index count, and passing to \(\mathcal{O}_B\) may again dissolve it.

**2.6 Comparison with Duke–Imamoğlu–Tóth.** DIT's cycle integrals of \(j\) along real-quadratic geodesics, and their linking numbers of modular knots, mix real-quadratic geodesics with modular data in the same way our bridge does (\(\operatorname{tr} Z = -2n\), length \(2\log\varepsilon\), CM points of \(1-n^2\)). Determine whether \(u\) (or \(\log|u|\), or \(\arg u\)) is expressible through DIT invariants — either outcome is valuable: a match imports their machinery; a mismatch means \(u\) is genuinely new.

**2.7 Kernel optimization: "Schmidt units".** The denominators come from the kernel's zeros at the elliptic points. Search weight-2 kernels (eta-quotients on \(\Gamma_0(\ell)\), Weber-function combinations) minimizing or eliminating denominators — ideally making \(\varepsilon\Theta\) a genuine algebraic *unit*. Also yields smaller class-polynomial-style certificates.
*Update after 2.3:* the target object very likely already exists — the
\(\mathfrak{r}\)-twisted \(\Delta\)-ratio
\(R_f = u_f^6\,\beta_2^4(\beta_2-1728)^3/\beta_1^4(\beta_1-1728)^3\) is
**certified to be a unit** at every computed level (palindromic integer level
polynomials with constant term 1; [phase-kronecker-limit.md](phase-kronecker-limit.md) Thm 5).
What remains: (a) prove unitness in general — a per-class statement that the
\(\mathfrak{r}\)-twist preserves every local \(\Delta\)-valuation, i.e. the
per-class refinement of the \(\Delta\)-mass law, natural territory for Gross's
quasi-canonical-lifting valuations (this would simultaneously settle euclidean
§6.1); (b) descend from \(R_f\) (a sixth-power-level object) to a first-power
kernel-level unit \(\varepsilon\Theta_g\).

**2.8 Genus-refined Gross–Zagier: the exponent law of the dressing (new, from 2.3).**
The certified tables of [phase-kronecker-limit.md](phase-kronecker-limit.md)
§4(b)/§6 express every genus-character \(j\)-sum as
\(k\log\varepsilon_d + \sum_p e_p\lambda_p\) with \(e_p \in \{\pm3\}\)
(\(\Sigma_0\)) resp. even (\(\Sigma_{1728}\)), supported exactly on
\(\mathrm{GZ}(D,-3)\) resp. \(\mathrm{GZ}(D,-4)\). Gross–Zagier's *On singular
moduli* computes genus-character-weighted prime factorizations of
\((j-j')\)-norms via their \(\epsilon(\mathfrak{n})\)-formula; matching it
against the tables should **prove** the exponents and signs. One unified
finite-place paper would settle at once: the exact denominator valuations
(1.3, moduli §5.7), the \(\lambda_n\)-law (euclidean §6.1), and plausibly the
sign anomaly of 1.1 — the exact-factorization machinery of `phase_klf.py`
already computes the sign of every \(A/B\)-decomposition, so the data side is
finished; only the GZ bookkeeping is missing.

## 3. Large (a paper or program each)

**3.1 The arithmetic length spectrum of the Bianchi orbifold.** The Cartan bridge sends level-\(n\) circle classes to closed geodesics of length \(2\log\varepsilon_n\) in \(\mathrm{PSL}_2(\mathbb{Z}[i])\backslash\mathbb{H}^3\). Make the correspondence exact (twisted/\(\sigma\)-conjugacy classes versus geodesics), and prove a multiplicity statement: the length \(2\log(n+\sqrt{n^2-1})\) occurs with multiplicity governed by \(H(n^2-1)\)-data. This is a Bianchi analogue of the classical class-number–geodesic dictionary — with the striking twist that *real* lengths carry *imaginary*-quadratic class numbers. Then: run it through the Selberg/Bianchi trace formula for spectral consequences (average of \(H(n^2-1)\); error terms against the spectrum).

**3.2 The weight-3/2 slicing program.** \(\sum_n t(n^2-1)q^n\) and \(\sum_n H(n^2-1)q^n\) are shifted-square slices of Zagier's and Cohen–Eisenstein's weight-3/2 forms. Identify these slices exactly (theta decomposition / Shimura lift; Mertens' mock-modular class number identities), obtaining closed generating-function identities for the Schmidt census — and, in the other direction, a *circle-geometric interpretation* of mock modularity.

**3.3 Kronecker–Hurwitz relations, bijectively.** Prove \(\sum_t H(4m - t^2) = \sum_{d\mid m}\max(d, m/d)\) by an explicit bijection/involution on Schmidt circles across \(\alpha\)-levels (the right-hand side counts cusp data — in the arrangement, tangencies to the horizontal lines). The arrangement is the first geometric model where all the classes in the relation coexist as concrete circles; a bijective proof would be a genuinely new take on Eichler–Selberg-type identities.

**3.4 Stark-adjacent theory and equidistribution of phases.** *2.3 has
succeeded — this program is now unlocked, and it is the most promising next
step.* Develop \(u\) systematically:
(a) **the Robert index** — determine the index of the subgroups generated by
the Euclidean \(\Delta\)-data \(\{G_\mathfrak{c}\}\) and the hyperbolic units
\(\{R_f\}\) inside the elliptic-unit / full \(S\)-unit groups of their ring
class fields (Robert; Kubert–Lang ch. 12–13): an index formula with class
numbers appearing as indices would be the *quantitative* form of "the Schmidt
arrangement is an elliptic-unit system", the analogue of the cyclotomic-unit
index \(= h^+\);
(b) **Stark recognition for the non-real characters** — the certified coset
cubics/quartics of phase-kronecker-limit §3 pin \(L'(0,\chi_3)\) etc. as logs
of explicit cubic algebraic numbers; identify them as the Stark units of the
dihedral cubic/quartic subfields (Schertz, *Complex Multiplication* ch. 6–7 —
the abelian-over-imaginary-quadratic Stark conjecture is a theorem via
elliptic units, so this is recognition, not conjecture);
(c) unit-normalization at first power (2.7) and the distribution of
\(\arg u_f\) over classes as \(n \to \infty\) (between Duke's CM-point
equidistribution and closed-geodesic equidistribution; the proved fiber-rate
\(\sqrt{\alpha^2-1}\) gives the deterministic part);
(d) the \(\mathbb{Q}(\sqrt n)\)-mystery of 2.3(iv): find the geometric carrier
of the Pell unit at Euclidean conductor \(n\) — both aspects now attach a
*real* quadratic unit to a level, and neither has a geometric explanation.

**3.5 Unlikely intersections for the coupled fiber product.** The \(\hat\sigma\)-pairs are CM points on \(X_0(\frac{n-1}2)\times_{X(1)}X_0(\frac{n+1}2)\) with **discriminant coupled to the level** — exactly the shape quantified by André–Oort/Zilber–Pink. Two directions: heights of the coupled points (Gross–Zagier/Kudla program: our denominators-as-GZ-primes is the finite-place shadow of a height formula — find the archimedean side), and whether the Schmidt family is, in a precise sense, the *complete* solution of its unlikely-intersection problem.
*Update after 2.3: the requested archimedean side is now computed* — the KLF
character sums are precisely the archimedean heights of the
\(\chi\)-components (note \(\log|u_f|\) is a combination of automorphic
Green's-function-type values \(\log|\Delta|\), \(\log|j - j(\text{elliptic})|\)
at the Heegner pair \((f, \mathfrak{r}f)\), i.e. along the graph of the
Atkin–Lehner twist). What remains is assembling KLF (archimedean) + 2.8
(finite places) into one height identity for the phase — with
[middle-kernel.md](middle-kernel.md)'s Green's-function collapse as the bridge
to the trace-formula program.

**3.6 Other imaginary quadratic fields.** Everything transports to Schmidt arrangements over \(\mathcal{O}_K\) for other \(K\) (Stange's general setting): the parity/congruence classification, the \(\alpha\)-invariant, the twist \(\mathfrak{r}\), the phase. Comparative questions: the Eisenstein case \(\mathbb{Z}[\omega]\) has \(\mu_6\)-units (richer torsion in the phase — a clean test of the \(\mu(B)\)-analysis); class number \(> 1\) fields (\(\mathcal{O}_{-5}\)) break which steps? The Catalan constant becomes \(L(2, \chi_{d_K})\); the unit \(\varepsilon\) is field-independent — meaning the *geodesic side is universal* while the arithmetic side varies. A uniform treatment would be the definitive paper on Schmidt-arrangement arithmetic.

**3.7 Tangency arithmetic and local-global.** Derive the Diophantine criterion for tangency between Schmidt circles (inversive pairing \(= \pm1\)), classify the tangency graph per level pair, and pose the Apollonian-style local-global question for curvatures in tangency components — importing Bourgain–Kontorovich technology into the Schmidt setting.

## 4. Consolidation

- **Two papers**: (I) *Counting and composing Schmidt circles* — classification, \(N_e\)/Catalan, \(3H(n^2-1)\), class formula for \(\hat\sigma\), circle-language composition (self-contained, elementary tools); (II) *The phase of a Schmidt circle* — \(\Theta\), functional equations, dihedral reciprocity, GZ-denominators (CM machinery). The documents already contain the proofs; the work is selection and polish.
- **Software**: consolidate the nine scripts into one module with a test suite; port the exact-ideal parts to Sage/PARI to replace PSLQ certification by exact arithmetic where possible. Guard rails learned this session: never trust a PSLQ fit with (terms)×(coefficient digits) near the working precision, and never let module import set `mp.dps`.
- **Literature diligence** before writing: Stange (Schmidt arrangements — how much of the classification/counting is folklore there), Sarnak (reciprocal geodesics), Duke–Imamoğlu–Tóth (cycle integrals), Gross–Zagier (singular moduli), Kubert–Lang/Schertz (Siegel units in ring class fields), Mertens/Bringmann–Kane (class number identities), Vlasenko–Zagier and Andersen (higher Kronecker limit formulas along real quadratic data — possibly the closest existing relatives of \(\log|u|\)).

## 5. Recommended order

(2.1 and 2.3 are done — first-power-descent.md and phase-kronecker-limit.md.)

1. (2.7 + 3.4a) **the Schmidt-unit program**: prove the unit property of the
   \(\mathfrak{r}\)-twisted \(\Delta\)-ratios \(R_f\) (per-class valuation
   cancellation — also closes euclidean §6.1), descend to a first-power unit,
   and compute the Robert index — turns the KLF answer into a structure
   theorem, and is the natural headline of the phase paper;
2. (2.8) **genus-refined Gross–Zagier** for the dressing exponents — one
   finite-place paper closing 1.1, 1.3 and the valuation questions
   simultaneously; the exact data is already computed;
3. (3.5 + euclidean 6.12) **the height formula for the phase** — archimedean
   side now in hand (KLF), finite side from 2.8; assembling them gives the
   Gross–Zagier-style identity along the Atkin–Lehner graph;
4. (3.1) the geodesic multiplicity statement — the deepest structural payoff
   of the whole bridge, now with the ortholength/spectral notes of
   spectral-outlook.md as staging;
5. literature pass (§4) before drafting — add Robert, Schertz
   (*Complex Multiplication*), and Stark's original papers to the list.
