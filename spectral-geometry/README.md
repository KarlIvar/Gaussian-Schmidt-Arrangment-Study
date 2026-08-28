# Spectral geometry of the Gaussian Schmidt arrangement

This folder connects the objects built in the companion documents — the classification congruences, the counting functions \(N_e(n)\) and \(3H(n^2-1)\), the invariant \(\alpha\), the Cartan matrices \(X\bar X^{-1}\), the level polynomials — to the **spectral geometry of the Bianchi orbifold**
$$
M \;=\; \mathrm{PSL}_2(\mathbb{Z}[i]) \backslash \mathcal{H}^3 .
$$
The key is to distinguish one surface: the totally geodesic plane of \(\mathcal{H}^3\) **above the real line**, \(P_{\hat{\mathbb{R}}} = \{(x, 0, t)\}\). Its \(\Gamma\)-orbit has the Schmidt arrangement as boundary shadow, its stabilizer is \(\mathrm{PGL}_2(\mathbb{Z})\), and its image in \(M\) is an immersed modular orbifold \(Y\) of area \(\pi/6\). Every counting theory of the project turns out to be a metric statement about the pair \((M, Y)\).

## The documents

**[orthogeodesics.md](orthogeodesics.md)** — *the geometric/length-spectrum side.*
- The dictionary: circles = sheets of \(Y\); the hyperbolic plane of [hyperbolic-counting.md](../hyperbolic-counting.md) **is** \(Y\) (canonical identification); the two families \(\mathcal{S}\), \(i\mathcal{S}\) = two immersed surfaces \(Y, Y'\), each with laminar (non-crossing) lifts, crossing each other only at right angles.
- \(\alpha = \cosh(\text{ortho-distance})\), and **the CM point \(m_1\) is the foot of the common perpendicular** on the surface (exact feet, distances, orthogonality — proved and machine-verified).
- The Cartan matrix \(Z = X\bar X^{-1}\) is the holonomy of that perpendicular's closed double: axis \((m_1, \bar m_1)\), length \(2\operatorname{arccosh} n = 2\log\varepsilon_n\), characteristic polynomial \(\lambda^2 \mp 2n\lambda + 1\) — verified in exact arithmetic for all 160 circles of all odd levels \(n \le 15\).
- **The ortho-length spectrum theorem**: the perpendicular lengths between the surfaces are \(\{\operatorname{arccosh} n\}_{n \ge 2}\) (odd \(n\): \(Y\!\leftrightarrow\!Y\); even: \(Y\!\leftrightarrow\!Y'\)) with Hurwitz-class-number multiplicities \(H(n^2-1)\) — real quadratic lengths carrying imaginary quadratic class numbers, as anticipated in [outlook.md](../outlook.md) 3.1.
- Growth \(\sum_{n\le X} H(n^2-1) \approx 0.286\,X^2\) (sieved to \(10^6\); the constant exceeds the random-discriminant heuristic \(\pi/12\zeta(3)\) by a family singular series \(\approx 1.31\), left open), and the \(\sigma\)-geodesics form an **entropy-\(\tfrac12\)** family: \(\asymp e^{L}\) against the full \(e^{2L}/2L\).
- The strata are decorated by the project's polynomials (length polynomial, Hilbert class polynomial of the feet, Zagier trace, level polynomial \(Q_n\)); irreducibility of \(Q_n\) says Galois permutes the ortho-arcs of each length transitively.

**[eisenstein-catalan.md](eisenstein-catalan.md)** — *the analytic/automorphic side.*
- Conventions and the scattering term of \(M\): \(\varphi(s) = \frac{\pi}{s}\frac{\zeta_K(s)}{\zeta_K(s+1)}\) (derived, and verified to \(10^{-7}\) by brute-force lattice summation of the Eisenstein series); Humbert's \(\operatorname{vol}(M) = G/3\) verified by direct integration; \(\operatorname{Res}_{s=1}\varphi = 1/(2\operatorname{vol} M)\).
- The counting Dirichlet series rewritten as \(\sum N_e(n) n^{-s} = \zeta(s-1)\zeta(s)/\zeta_K(s)\): a zeta ratio of scattering (Asai) type.
- **The Catalan volume identity** (closing Remark 1 of [euclidean-counting.md](../euclidean-counting.md)):
  $$\sum_{n \le X} N_e(n) \sim \frac{X^2}{\pi}\cdot\frac{\operatorname{area}(Y)}{\operatorname{vol}(M)} = \frac{X^2}{2G},$$
  proved through an incidence kernel \(\Theta_s\) summing a distance-weight to the plane orbit: its cusp constant term is an exact series over the circle count (with the elementary "kinematic lemma" \(\int_0^\infty J_s(x - \frac1{4x})dx = C_s\) reconciling the two evaluations), its global mass is \(\operatorname{area}/\operatorname{vol}\), and expanding horospheres equidistribute.
- **The Mellin/period identity**: \(\int_0^\infty \operatorname{Hem}_s(t)\,t^{w-2}dt = \pi K_s(w)\,\zeta(w)\zeta(w+1)/\zeta_K(w+1)\) with explicit Gamma factors — the regularized pairing of the surface kernel with the Eisenstein series, whose residue at \(w = 1\) matches the Maass–Selberg residue exactly (\(\pi C_s/2G\) both ways). Cusp forms sit far away (\(\lambda_1 \approx 44.85\), no exceptional eigenvalues), so the error theory is governed by \(\zeta_K\).
- The thin companion: the **atoms** of [half-plane-monoid.md](../half-plane-monoid.md) count with exponent \(\delta = 1.30568\ldots\) (exact Descartes enumeration fits \(1.3033\)) — base eigenvalue \(\lambda_0 = \delta(2-\delta)\) instead of a volume: the lattice/thin dichotomy inside one arrangement.

## Verification

Two new scripts in [../scripts](../scripts), following the project's machine-verification convention:

| script | checks | status |
|---|---|---|
| [spectral_verify.py](../scripts/spectral_verify.py) | volume, \(\zeta_K\), scattering, Dirichlet series, kinematic lemma, horosphere limit, Mellin identity, residues | 23/23 pass |
| [orthospectrum_verify.py](../scripts/orthospectrum_verify.py) | exact holonomy/axis (160 circles), orthogeodesic geometry, \(H(n^2-1)\) counts, Hurwitz growth + entropy, gasket exponent | 8/8 pass |

Run: `python3 scripts/spectral_verify.py` (≈ 1–2 min; `--fast` skips the lattice sum) and `python3 scripts/orthospectrum_verify.py` (seconds). Both need `numpy` and `mpmath`.

## Literature notes (for the diligence pass of [outlook.md](../outlook.md) §4)

The dictionary "Schmidt arrangement = boundary of a plane orbit" is implicit in K. Stange's papers (*The Apollonian structure of Bianchi groups*, arXiv:1505.03121; *Visualising the arithmetic of imaginary quadratic fields*), which should be cited as the origin of the geodesic-plane viewpoint. Spectral background: Elstrodt–Grunewald–Mennicke, *Groups Acting on Hyperbolic Space* (Eisenstein series, fundamental domain, volumes); Humbert's volume formula; Then, *Arithmetic quantum chaos of Maass waveforms* (math-ph/0305048) for the Picard-group eigenvalue data quoted (\(r_1 = 6.6221193\ldots\)). Counting/equidistribution: Duke–Rudnick–Sarnak, Eskin–McMullen (orbit counting via mixing), Oh–Shah and Kontorovich–Oh (circle counting, thin groups), McMullen (the dimension \(\delta = 1.30568\ldots\)). Periods: Flicker (Asai/distinction), Zagier (Rankin–Selberg regularization), Kudla–Millson (cycle theta lifts). Reciprocal geodesics: Sarnak. Prime geodesic theorem for \(\mathrm{PSL}_2(\mathbb{Z}[i])\): arXiv:1903.05111. These attributions are recalled from the literature and were spot-checked online during this session, but the full diligence pass before any write-up for publication still applies.
