# Outlook: where the Gaussian Schmidt arrangement project can go next

## 0. Where we stand

Proved and machine-certified across the six documents:

1. **Diophantine classification** of the circles of \(\mathcal{S} = \mathrm{PSL}_2(\mathbb{Z}[i])\hat{\mathbb{R}}\) ([circle-classification.md](circle-classification.md)); Euclidean counting \(N_e(n) = n\prod(1-\chi_{-4}(p)/p)\) with Catalan-constant asymptotics ([euclidean-counting.md](euclidean-counting.md)).
2. **Hyperbolic counting**: circles at level \(\alpha = n\) are hyperbolic circles of radius \(\operatorname{arcoth} n\) at CM points of discriminant \(1-n^2\); weighted count \(3H(n^2-1)\) in the ideal triangle ([hyperbolic-counting.md](hyperbolic-counting.md)).
3. **The involution** \(\sigma(X) = \bar X^{-1}\): class formula \(\hat\sigma[f] = [\mathfrak{r}_n][f]^{-1}\), proved via the explicit unitary basis ([class-formula-proof.md](class-formula-proof.md)); circle-level class group operations ([circle-composition.md](circle-composition.md)).
4. **The phase invariant** \(\Theta\), the sixth coordinate on double cosets; functional equations (mirror, inversion = \(e^{-\ell}\)); the normalized \(u_f = \varepsilon\Theta_f\); full generalized-dihedral Galois equivariance of \(u^{12}\) via Shimura–Siegel reciprocity, with the \(\xi\)-torsion closed by the cocycle-ideal lemma; explicit fields for the pair-sums via genus characters; denominators = Gross–Zagier primes of \((1-n^2, -3)\) and \((1-n^2, -4)\) ([moduli-invariants.md](moduli-invariants.md)); and, as a corollary of the equivariance, the phase units of a level are the roots of **one polynomial with integer coefficients** (§5.9 there: unconditionally in degree \(12h\), certified in the sharp degree \(h\) for odd \(n \le 17\)), and that polynomial is **irreducible** at every computed level — one level, one Galois orbit, the degree-\(h\) polynomial being the common minimal polynomial of the level's phase units (§5.10 there).
5. **The Euclidean moduli theory** (translations \(\times\ \mathrm{SL}_2(\mathbb{Z})\) instead of \(\mathrm{SL}_2(\mathbb{Z}) \times \mathrm{SL}_2(\mathbb{Z})\)): disks of curvature \(2n\) mod translation \(\leftrightarrow\) \(\mathrm{Cl}(\mathbb{Z} + n\mathbb{Z}[i])\) (so \(N_e(n) = 2h(-4n^2)\)); \(j\)-values = ring class polynomial squared; trace slice \(t(4n^2)\); the Euclidean sixth invariant \(\Theta = j'(w)/c^2\) and its lemniscatic phase units with monic integer level polynomials; the \(\Delta\)-mass law **proved** (Theorem 4 there: constancy of the full Hecke mass + PID stratification + generating-function identities + the conjugation sign analysis), and the level polynomials \(P^{(2)}_n, P^{(6)}_n\) **proved irreducible for every \(n\)** (Theorem 5 there: the root-free closed form \(u^2 = -12\beta(\beta-1728)\,g_2(\Lambda)/g_2(\mathbb{Z}[i])\), a cocycle-free square-level Shimura translation law, and an archimedean dominance lemma killing the coincidence subgroup — so the square-level analogues of both hyperbolic open hypotheses, first-power descent and \(T = 1\), are closed on the Euclidean side); \(H_{-4n^2}\) irreducible classically, and the first-power disk polynomial \(P^{(2)}(x^2)\) irreducible at every computed level ([euclidean-moduli-invariants.md](euclidean-moduli-invariants.md), with its own outlook §6 there).

What follows are next steps, roughly by size. Each: the question, why it should work, the first move.

## 1. Small (days to weeks)

**1.1 The sign of \(u_f\).** Fresh data (session of 2026-08-28): on ambiguous classes, \(\operatorname{sign}(u_f) = -1\) identically for \(n = 3, 5, 7, 9, 13, 17\), but at \(n = 11, 15\) it equals \(-\psi(x)\), \(\psi\) the quadratic character of the \(\langle\mathfrak{r}\rangle\)-coset:

| \(n\) | 3 | 5 | 7 | 9 | 11 | 13 | 15 | 17 |
|---|---|---|---|---|---|---|---|---|
| sign pattern | − | −− | −− | −− | −++− (by coset) | −−−− | −+−+ (by coset) | −−−− |

Note \(n = 11\) and \(13\) have identical class groups \((\mathbb{Z}/2)^2\) but different sign behavior. Candidate explanations: a 2-adic genus character (the crude prime-disc decomposition mishandles conductor-2 parts), or the sign of the Gross–Zagier monomial. First move: recompute with proper conductor-aware characters; extend to \(n \le 41\). A clean sign law would complete the real-quadratic description of \(u\) on ambiguous classes.

**1.2 Exact formulas for the counts via Eichler's relations.** \(H(n^2-1)\) is the \(t = \pm1\) term of the classical relation \(\sum_{t^2\le 4m} H(4m-t^2) = \ldots\) Combined with Mertens/Bringmann–Kane-type mock-modular class number identities, one should get exact recurrences for the triangle counts \(3H(n^2-1)\) — closed-form arithmetic for our census. First move: assemble the relation for \(m = (n^2-1)/4 + \) shifts and isolate the diagonal.

**1.3 Exact denominator valuations.** §5.7 identified the support and the squared-vs-first-power pattern; Gross–Zagier gives multiplicities of each \(\frac{e|D|-x^2}{4}\). Trace them through the \((4,3)\)-exponents and the 6th root to prove \(v_p(\operatorname{den} u)\) exactly. Also: nobody has looked at the **numerators** — their primes are not GZ-constrained and may encode height data.

**1.4 Even levels.** The whole \(\alpha\)-theory was run at odd \(n\) (family \(\mathcal{S}\)). Even \(n\) lives in \(i\mathcal{S}\), discriminants \(1-n^2 \equiv 3 \pmod 4\), twist ideals of odd norm \(\frac{n\mp1}{2}\)-adjusted. Everything should transport with cleaner 2-adic behavior (odd discriminant!) — a good consistency laboratory. First move: adapt `build_P` to the \(i\mathcal{S}\)-parity and rerun the pipelines.

**1.5 Regularizing \(\alpha = 2\).** The level \(n = 2\) (the inscribed circle; CM point \(\rho\)) is excluded because \(j'(\rho) = 0\). Both \(m_1, m_2\) are \(\rho\)-points, so the *ratio* of leading Laurent coefficients of any kernel defines a regularized \(\Theta(2)\). Compute it; it should be a particularly clean number (the class group is trivial and \(\varepsilon = 2+\sqrt3\)).

**1.6 A phase atlas.** Render the circles of the ideal triangle colored by \(\arg u\) and \(\log|u|\) across levels — a "phase portrait of the Schmidt arrangement". Cheap, likely to expose patterns (e.g. the sign dichotomy of 1.1 geometrically), and the natural figure for any eventual paper.

**1.7 The \(\zeta_8\)-levels.** \(n^2 - 1 = 2\square\) (\(n = 3, 17, 99, \ldots\) — NSW-adjacent) are exactly where \(B = \mathbb{Q}(\zeta_8)\) and the \(u^6\)-translation may pick up a sign. At \(n = 17\) the data shows all signs \(-\): check whether the metaplectic sign actually occurs, or the \(\mu_8\)-caveat is vacuous.

## 2. Medium (weeks to months)

**2.1 First-power descent — the isogeny-differential route (recommended).** The remaining \(\mu_{12}\)-cocycle would evaporate if \(u\) has a purely algebraic-geometric expression. Candidate: the \(\hat\sigma\)-pair carries two cyclic isogenies \(\phi_\mp: E_2 \to E_1'\)-conjugates of degrees \(\frac{n\mp1}{2}\) whose composite is multiplication by \(\omega\) (norm \(\frac{n^2-1}{4}\), and \(\omega^2 = D/4\)); their pullbacks of Néron differentials give scalars \(\lambda_\mp\) with \(\lambda_-\lambda_+ \sim \omega\). Conjecture: \(u\) is a monomial in \(\varepsilon, \lambda_\mp, \mu\) — testable by PSLQ in an afternoon. If true, \(u\) is Galois-theoretically transparent at first power (isogenies are algebraic), law 3 closes at \(u\)-level with no Siegel-unit bookkeeping, and the invariant gets a moduli-theoretic definition. If false, fall back on:

**2.2 First-power descent — the Kubert–Lang route.** Write the kernel as \(\eta^4\cdot(\text{level-6 weight-0 function})\); \(u\) becomes a Siegel–Ramachandra-type invariant; the known first-power Galois action of Siegel units (Kubert–Lang ch. 11, Schertz) resolves the sixth root. Routine, longer, certain to work.

**2.3 Kronecker limit formula for the phase: \(\log|u|\) and \(L'\)-values.** From the closed form, \(\log|u_f|\) is a linear combination of \(\log|\Delta(\mathfrak{a})|\)-terms, logs of singular-moduli differences, and \(\log\)-rationals. Character sums \(\sum_f \chi(f)\log|u_f|\) should therefore evaluate to explicit combinations of \(L'(0,\chi)\) (Kronecker limit formula) and \(\log p\) over GZ-primes (Gross–Zagier factorization). First move: compute the sums numerically for all characters at \(n = 9, 11, 13\) and PSLQ against \(\{L'(0,\chi), \log \varepsilon, \log p\}\). A clean identity here makes \(u\) a *geometrically defined Stark-adjacent unit system* — in my view the strongest single result within reach.

**2.4 Cusp degeneration: the phase at \(n \to 1\) and Dedekind sums.** As \(n \to 1\), \(\varepsilon \to 1\), the geodesic shrinks into the cusp, and the circles become Ford circles. The degenerate limit of \(u\) (and of the \(\zeta\)-cocycle, which is built from \(\eta\)-multipliers) should be classical Dedekind-sum data. Making this precise would tie the phase theory to the classical modular-transformation literature and likely *explain* the \(\mu_{12}\).

**2.5 The imprimitive strata.** The class formula's behaviour on non-invertible classes (content \(g > 1\), ramified prime 2) is verified but unproved. The tools now available (cocycle ideal, \(\mathcal{O}_B\)-extension trick) should handle it: non-invertibility of \(\mathfrak{a}_f\) obstructs only one index count, and passing to \(\mathcal{O}_B\) may again dissolve it.

**2.6 Comparison with Duke–Imamoğlu–Tóth.** DIT's cycle integrals of \(j\) along real-quadratic geodesics, and their linking numbers of modular knots, mix real-quadratic geodesics with modular data in the same way our bridge does (\(\operatorname{tr} Z = -2n\), length \(2\log\varepsilon\), CM points of \(1-n^2\)). Determine whether \(u\) (or \(\log|u|\), or \(\arg u\)) is expressible through DIT invariants — either outcome is valuable: a match imports their machinery; a mismatch means \(u\) is genuinely new.

**2.7 Kernel optimization: "Schmidt units".** The denominators come from the kernel's zeros at the elliptic points. Search weight-2 kernels (eta-quotients on \(\Gamma_0(\ell)\), Weber-function combinations) minimizing or eliminating denominators — ideally making \(\varepsilon\Theta\) a genuine algebraic *unit*. Also yields smaller class-polynomial-style certificates.

## 3. Large (a paper or program each)

**3.1 The arithmetic length spectrum of the Bianchi orbifold.** The Cartan bridge sends level-\(n\) circle classes to closed geodesics of length \(2\log\varepsilon_n\) in \(\mathrm{PSL}_2(\mathbb{Z}[i])\backslash\mathbb{H}^3\). Make the correspondence exact (twisted/\(\sigma\)-conjugacy classes versus geodesics), and prove a multiplicity statement: the length \(2\log(n+\sqrt{n^2-1})\) occurs with multiplicity governed by \(H(n^2-1)\)-data. This is a Bianchi analogue of the classical class-number–geodesic dictionary — with the striking twist that *real* lengths carry *imaginary*-quadratic class numbers. Then: run it through the Selberg/Bianchi trace formula for spectral consequences (average of \(H(n^2-1)\); error terms against the spectrum).

**3.2 The weight-3/2 slicing program.** \(\sum_n t(n^2-1)q^n\) and \(\sum_n H(n^2-1)q^n\) are shifted-square slices of Zagier's and Cohen–Eisenstein's weight-3/2 forms. Identify these slices exactly (theta decomposition / Shimura lift; Mertens' mock-modular class number identities), obtaining closed generating-function identities for the Schmidt census — and, in the other direction, a *circle-geometric interpretation* of mock modularity.

**3.3 Kronecker–Hurwitz relations, bijectively.** Prove \(\sum_t H(4m - t^2) = \sum_{d\mid m}\max(d, m/d)\) by an explicit bijection/involution on Schmidt circles across \(\alpha\)-levels (the right-hand side counts cusp data — in the arrangement, tangencies to the horizontal lines). The arrangement is the first geometric model where all the classes in the relation coexist as concrete circles; a bijective proof would be a genuinely new take on Eichler–Selberg-type identities.

**3.4 Stark-adjacent theory and equidistribution of phases.** Assuming 2.3 succeeds: develop \(u\) systematically — integrality after kernel optimization (2.7), index of the subgroup generated by \(\{u_f\}\) in the relevant unit/S-unit groups, and the distribution of \(\arg u_f\) over classes as \(n \to \infty\) (sitting between Duke's CM-point equidistribution and closed-geodesic equidistribution; the proved fiber-rate \(\sqrt{\alpha^2-1}\) gives the deterministic part).

**3.5 Unlikely intersections for the coupled fiber product.** The \(\hat\sigma\)-pairs are CM points on \(X_0(\frac{n-1}2)\times_{X(1)}X_0(\frac{n+1}2)\) with **discriminant coupled to the level** — exactly the shape quantified by André–Oort/Zilber–Pink. Two directions: heights of the coupled points (Gross–Zagier/Kudla program: our denominators-as-GZ-primes is the finite-place shadow of a height formula — find the archimedean side), and whether the Schmidt family is, in a precise sense, the *complete* solution of its unlikely-intersection problem.

**3.6 Other imaginary quadratic fields.** Everything transports to Schmidt arrangements over \(\mathcal{O}_K\) for other \(K\) (Stange's general setting): the parity/congruence classification, the \(\alpha\)-invariant, the twist \(\mathfrak{r}\), the phase. Comparative questions: the Eisenstein case \(\mathbb{Z}[\omega]\) has \(\mu_6\)-units (richer torsion in the phase — a clean test of the \(\mu(B)\)-analysis); class number \(> 1\) fields (\(\mathcal{O}_{-5}\)) break which steps? The Catalan constant becomes \(L(2, \chi_{d_K})\); the unit \(\varepsilon\) is field-independent — meaning the *geodesic side is universal* while the arithmetic side varies. A uniform treatment would be the definitive paper on Schmidt-arrangement arithmetic.

**3.7 Tangency arithmetic and local-global.** Derive the Diophantine criterion for tangency between Schmidt circles (inversive pairing \(= \pm1\)), classify the tangency graph per level pair, and pose the Apollonian-style local-global question for curvatures in tangency components — importing Bourgain–Kontorovich technology into the Schmidt setting.

## 4. Consolidation

- **Two papers**: (I) *Counting and composing Schmidt circles* — classification, \(N_e\)/Catalan, \(3H(n^2-1)\), class formula for \(\hat\sigma\), circle-language composition (self-contained, elementary tools); (II) *The phase of a Schmidt circle* — \(\Theta\), functional equations, dihedral reciprocity, GZ-denominators (CM machinery). The documents already contain the proofs; the work is selection and polish.
- **Software**: consolidate the nine scripts into one module with a test suite; port the exact-ideal parts to Sage/PARI to replace PSLQ certification by exact arithmetic where possible. Guard rails learned this session: never trust a PSLQ fit with (terms)×(coefficient digits) near the working precision, and never let module import set `mp.dps`.
- **Literature diligence** before writing: Stange (Schmidt arrangements — how much of the classification/counting is folklore there), Sarnak (reciprocal geodesics), Duke–Imamoğlu–Tóth (cycle integrals), Gross–Zagier (singular moduli), Kubert–Lang/Schertz (Siegel units in ring class fields), Mertens/Bringmann–Kane (class number identities), Vlasenko–Zagier and Andersen (higher Kronecker limit formulas along real quadratic data — possibly the closest existing relatives of \(\log|u|\)).

## 5. Recommended order

1. (2.1) the isogeny-differential test for \(u\) — one afternoon, potentially reshapes everything downstream;
2. (2.3) the \(L'\)-character-sum experiment — the likely headline theorem;
3. (1.1)–(1.3) the small certainties (sign law, Eichler formulas, valuations) to round out paper I;
4. (3.1) the geodesic multiplicity statement — the deepest structural payoff of the whole bridge;
5. literature pass (§4) before drafting.
