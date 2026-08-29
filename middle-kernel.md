# The middle kernel: \(\Gamma\)-bi-invariant sums \(S_f(X,Y) = \sum_\gamma f(X\gamma Y)\) collapse onto the modular surface

Companion to [atomic-census.md](atomic-census.md) and [product-cocycle.md](product-cocycle.md). Question addressed (posed 2026-08-29): for \(\Gamma\)-bi-invariant \(f : \mathrm{SL}_2(\mathbb{C}) \to V\), find \(f\) — ideally with \(\dim V < 6\) — for which
$$
S_f(X, Y) \;=\; \sum_{\gamma \in \Gamma} f(X\gamma Y)
$$
has a simple form in the \((\alpha, \beta_1, \beta_2, u)\)-coordinates of the two double cosets; best of all, \(S_f(X,Y)\) a function of \(f(X), f(Y)\) alone.

**Answer in one paragraph.** Take \(f\) *radial*: \(f = k\circ\alpha\) (\(\dim V = 1\); these are exactly the \(\mathrm{SL}_2(\mathbb{R})\)-bi-invariant functions). Then \(S_f\) is nothing other than **the automorphic kernel of the modular surface evaluated at the pair of Heegner points carrying \(\beta_2(X)\) and \(\beta_1(Y)\)** — the whole middle sum collapses onto \(\mathrm{PSL}_2(\mathbb{Z})\backslash\mathbb{H}^2\). Consequently: (i) \(S_f\) depends on precisely four of the twelve invariants, namely \((\alpha_X, \beta_2(X); \alpha_Y, \beta_1(Y))\) — the phases \(u\) and the outer moduli \(\beta_1(X), \beta_2(Y)\) drop out identically; (ii) for the canonical kernel \(f_s = \varepsilon_\alpha^{-s}(\alpha^2-1)^{-1/2}\) the sum has the exact rank-one leading behaviour \(S_{f_s} \sim \frac{6}{s}\,f_s(X)\,f_s(Y)\) at the pole \(s = 0\), with class-independent residue; (iii) exact factorization at fixed \(s\) is *impossible* — the obstruction is exactly the cuspidal-plus-Eisenstein spectrum of the modular surface, and the class-dependence of \(S_f\) is a Maass-form Weyl sum at Heegner points (Duke's theorem territory). Everything is machine-verified: [scripts/middle_kernel.py](scripts/middle_kernel.py), **654 checks passing**.

## 1. The radial class and the master formula

\(\Gamma\)-bi-invariance alone leaves the 6-dimensional space \(\Gamma\backslash G/\Gamma\) of values, so a generic \(f\) gives \(S_f\) depending on everything. The natural sub-class is the *maximally symmetric* one: \(f = k\circ\alpha\), which is bi-invariant under all of \(H = \mathrm{SL}_2(\mathbb{R})\) (since \(\alpha(g) = \langle M_0, M_g\rangle\) and \(H\) stabilises \(M_0\)) — the analogue of the spherical Hecke subalgebra. For these:

> **Theorem 1 (master formula).** Let \(X, Y \in \Omega\) with levels \(a = \alpha(X) \ge 2\), \(b = \alpha(Y) \ge 2\), and set
> $$
> z_X := m_1(\sigma(X)) \quad (\text{the CM point of the class } [\mathfrak{r}_a f_X],\ \text{i.e. the point carrying } \beta_2(X)),
> \qquad
> z_Y := m_1(Y) \quad (\text{carrying } \beta_1(Y)).
> $$
> Then for every \(\gamma \in \Gamma\):
> $$
> \boxed{\;\alpha(X\gamma Y) \;=\; ab \;+\; \sqrt{(a^2-1)(b^2-1)}\;\cosh d_{\mathbb{H}}\bigl(z_X,\ \gamma\, z_Y\bigr)\;}
> $$
> where \(d_{\mathbb{H}}\) is hyperbolic distance in \(\mathbb{H}^2\) and \(\gamma\) acts by Möbius transformations. Hence for radial \(f = k \circ \alpha\):
> $$
> S_f(X, Y) \;=\; 2\sum_{\gamma \in \mathrm{PSL}_2(\mathbb{Z})} K_{a,b}\bigl(d_{\mathbb{H}}(z_X, \gamma z_Y)\bigr),
> \qquad
> K_{a,b}(d) := k\bigl(ab + \sqrt{(a^2-1)(b^2-1)}\cosh d\bigr)
> $$
> — **the classical automorphic point-pair kernel of the modular surface, evaluated at the Heegner pair \((z_X, z_Y)\).**

*Proof.* By the product-level formula ([atomic-census.md](atomic-census.md) Prop. 4), \(\alpha(X\gamma Y) = \langle M_{X^{-1}}, \gamma M_Y\rangle\). By [involution.md](involution.md) §6, \(M_{X^{-1}} = -\overline{M_{\sigma X}}\); writing \(M_{\sigma X} = (q_1, x_1 + i a, m_1)\) and \(\gamma M_Y = (q_2, x_2 + ib, m_2)\) (both honest \(\Omega\)-disks), the pairing of \(-\bar M\) against \(M'\) evaluates to
\(\langle -\bar M_{\sigma X}, \gamma M_Y\rangle = 2ab - \langle M_{\sigma X}, \gamma M_Y\rangle\),
and Lemma 2 of [product-cocycle.md](product-cocycle.md) gives \(\langle M_{\sigma X}, \gamma M_Y\rangle = ab - \cosh d(z_X, \gamma z_Y)\sqrt{(a^2-1)(b^2-1)}\). (The factor \(2\) in \(S_f\) is \(-I \in \Gamma\).) \(\square\)

Superadditivity (\(\alpha(X\gamma Y) \ge ab + \sqrt{(a^2-1)(b^2-1)}\)) is now simply \(\cosh d \ge 1\), with the aligned/equality case \(d = 0\), i.e. \(\gamma z_Y = z_X\). Verified to 45 digits on 120 random products — including the pleasing rigidity that \(ab + \sqrt{(a^2-1)(b^2-1)}\cosh d(z_X, \gamma z_Y)\) is always an **odd integer**: the distance spectrum between the two Heegner orbits is quantized by the arrangement.

> **Corollary 2 (what \(S_f\) depends on).** For radial \(f\), whenever the sum converges absolutely (\(k(n) \ll n^{-1-\epsilon}\) suffices, by the hyperbolic circle problem; levels \(\ge 2\) so the relevant stabilisers are finite):
> $$
> S_f(X, Y) \;=\; \mathcal{K}_f\bigl(\alpha_X,\ \beta_2(X)\ ;\ \alpha_Y,\ \beta_1(Y)\bigr)
> $$
> — a function of four of the twelve coordinates. The phases \(u_X, u_Y\) and the outer moduli \(\beta_1(X), \beta_2(Y)\) drop out **identically**, and only the "inner" pair \((\beta_2(X), \beta_1(Y))\) survives, exactly matching the marginal-rigidity picture of [product-cocycle.md](product-cocycle.md). The twist symmetry \(S_f(X, Y) = S_f(\sigma Y, \sigma X)\) holds on the nose.

## 2. Spectral form, the \(\varepsilon\)-kernel, and the exact rank-one limit

The spectral expansion of automorphic kernels on \(Y_{\mathrm{mod}} = \mathrm{PSL}_2(\mathbb{Z})\backslash\mathbb{H}^2\) now applies verbatim:
$$
\tfrac12\,S_f(X,Y) \;=\;
\underbrace{\frac{1}{\mathrm{vol}(Y_{\mathrm{mod}})}\int_{\mathbb{H}^2} K_{a,b}(d(z, w))\,d\mu(w)}_{\text{rank-one main term}}
\;+\; \sum_j \hat h_{a,b}(t_j)\,\phi_j(z_X)\overline{\phi_j(z_Y)}
\;+\; \frac{1}{4\pi}\int_{\mathbb{R}} \hat h_{a,b}(t)\,E(z_X, \tfrac12{+}it)\overline{E(z_Y, \tfrac12{+}it)}\,dt,
$$
with \(\hat h_{a,b}\) the Selberg transform of \(K_{a,b}\) and \(\phi_j\) the level-one Maass forms. Three consequences.

**(a) The canonical kernel factors in the main term.** Take
$$
f_s(g) \;:=\; \frac{\varepsilon_{\alpha(g)}^{-s}}{\sqrt{\alpha(g)^2 - 1}}, \qquad \varepsilon_n = n + \sqrt{n^2-1} = e^{\operatorname{arccosh} n}.
$$
The layer integral evaluates exactly (substituting \(n = ab + P\cosh\rho\), \(P = \sqrt{(a^2-1)(b^2-1)}\), and using \(\varepsilon_{ab+P} = \varepsilon_a\varepsilon_b\) — the multiplicativity of \(\varepsilon\) at the aligned point):
$$
\int_{\mathbb{H}^2} K_{a,b}\,d\mu \;=\; \frac{2\pi}{P}\int_{ab+P}^{\infty} k_s(n)\,dn \;=\; \frac{2\pi}{s}\,\frac{\varepsilon_a^{-s}\varepsilon_b^{-s}}{P},
$$
so with \(\mathrm{vol}(Y_{\mathrm{mod}}) = \pi/3\):
$$
\boxed{\;\tfrac12\,S_{f_s}(X, Y) \;=\; \frac{6}{s}\, f_s(X)\, f_s(Y) \;+\; \bigl(\text{cuspidal} + \text{Eisenstein at } (z_X, z_Y)\bigr).\;}
$$
Moreover the counting main term is class-independent, so the pole is exact:
$$
\lim_{s \to 0^+} s\cdot\tfrac12 S_{f_s}(X,Y) \;=\; \frac{6}{\sqrt{(\alpha_X^2-1)(\alpha_Y^2-1)}} \;=\; 6\,f_0(X)\,f_0(Y)
$$
— **at the leading pole, \(S_f\) is exactly a rank-one function of \((f_0(X), f_0(Y))\), for a \(\dim V = 1\) kernel.** This is the precise sense in which the requested factorization holds.

**(b) No-go at fixed \(s\).** Exact factorization \(S_f(X,Y) = F(f(X))G(f(Y))\) for all pairs would force the automorphic kernel to be rank one, i.e. \(\hat h_{a,b}(t_j) = 0\) for all Maass forms and \(\hat h_{a,b} \equiv 0\) on the critical line; for an admissible kernel the Selberg transform is analytic in a strip, so it would vanish identically and \(k = 0\). **The obstruction to the ideal answer is precisely the existence of cusp forms on the modular surface.** Numerically the obstruction is very real: at levels \((9,9)\), changing the class of \(X\) from principal to \((3,2,7)\) changes \(S_{f_2}\) by \(18\%\) and \(S_{f_3}\) by \(38\%\).

**(c) The failure is Duke-equidistribution arithmetic.** The class-dependence of \(S_f\) sits entirely in \(\phi_j(z_X)\overline{\phi_j(z_Y)}\) and the Eisenstein values — Maass forms at Heegner points. Averaging over the class group on either side (levels fixed) turns these into Weyl sums over the full Heegner orbit of discriminant \(1 - \alpha^2\), which decay with power savings by Duke's theorem: **class-averaged, the operation is asymptotically rank one**,
$$
\frac{1}{h_X h_Y}\sum_{[f_X], [f_Y]} S_{f_s} \;\longrightarrow\; \frac{12}{s}\,f_s\otimes f_s \quad (\alpha_X\ \text{or}\ \alpha_Y \to \infty).
$$

**Numerics** (levels \((9,9)\), complete orbit enumeration to pairing \(4000\)): counting law \(\#\{\gamma : \alpha(X\gamma Y) \le T\}\) vs \(6T/P\): ratios \(0.92, 0.96, 0.97\) at \(T = 1500, 3000, 4000\) (the expected \(O(T^{-1/3})\)-flavoured convergence); rank-one ratios \(S/\text{main} = 0.58\)–\(0.69\) at \(s = 2\) (spectral terms of the same order at fixed \(s\), as theory predicts — the factorization is exact only at the pole and after class-averaging).

## 3. Distinguished kernel choices

- **Counting kernel** \(k = \mathbf{1}_{[1,T]}\): \(S_f = \#\{\gamma : \alpha(X\gamma Y) \le T\}\), the hyperbolic circle problem at Heegner pairs: \(= \frac{12\,T}{\sqrt{(a^2-1)(b^2-1)}} + O_{a,b}(T^{2/3})\) with Selberg's classical error exponent.
- **Resolvent/Green kernels**: choosing \(K_{a,b}(d) = G_s(\cosh d)\) (the free resolvent) makes \(S_f(X,Y)\) the **automorphic Green's function \(G_s(z_X, z_Y)\) at a pair of CM points of two different discriminants** \(1-a^2\), \(1-b^2\). At integral \(s = j \ge 2\) these are the *higher Green's functions* of Gross–Kohnen–Zagier, whose values at CM pairs are logarithms of algebraic numbers (conjectured by GZ, proved by Yingkun Li, arXiv:1810.13214 / Invent. Math.). So there exists a \(\Gamma\)-bi-invariant \(f\) for which \(S_f(X, Y)\) is, up to explicit periods, \(\log|\text{algebraic}|\) — a "simple form" of a deep kind, and exactly the kind of object the phase program ([outlook.md](outlook.md) 2.3, \(\log|u_f|\) vs \(L'\)-values) is built from. The two programs plausibly meet here.
- **Legendre kernels** \(k = Q_{s-1}\circ(\cdot)\): make \(\hat h\) fully explicit (standard transforms), convenient for the meromorphic continuation of the "monoid zeta functions" below.

## 4. What this buys for the atoms of \(\Omega\)

1. **The multiplication table of atoms has modular-surface asymptotics.** For atom classes \([A], [B]\) at levels \(a, b\), the cell sizes of [atomic-census.md](atomic-census.md) obey
   \(\#\{\text{double cosets in } \Gamma A\Gamma B\Gamma \text{ of level} \le T\} = \frac{c_{A,B}\,T}{\sqrt{(a^2-1)(b^2-1)}} + O(T^{2/3})\),
   with \(c_{A,B}\) explicit (torsion/fiber bookkeeping) — e.g. the depth-2 constants observed in the census are the \(6T/P\) law in disguise.
2. **The atomic transfer operator is an automorphic object.** The entries of the operator \(\mathsf{T}\) of [atomic-census.md](atomic-census.md) §8.1 are the multiplicities \(r(n) = \#\{\gamma : \cosh d(z_A, \gamma z_B) = \frac{n - ab}{P}\}\) — hyperbolic lattice counts between Heegner points. Their generating Dirichlet series \(\sum_n r(n)\varepsilon_n^{-s}\) (the natural "monoid zeta" of a cell) inherits meromorphic continuation and a spectral pole structure from the modular surface: **any Mayer-type Fredholm determinant for \(\Omega\) (gap G4) must contain the Maass spectrum of \(\mathrm{PSL}_2(\mathbb{Z})\)** — a sharp, previously invisible constraint on that program.
3. **The \(\alpha\)-spectrum question in a new metric.** Which levels \(n\) occur in a cell is now: which quantized distances \(\frac{n-ab}{P}\) occur between two Heegner orbits — a question about representation of integers by the pairing form of two CM points, where genus theory and Duke-type equidistribution both apply. This reframes the local–global problem Q1 of [half-plane-monoid.md](half-plane-monoid.md) (the atoms' \(\alpha\)-spectrum) inside a classical counting framework.
4. **The aligned stratum is the Heegner-coincidence stratum**: \(\varepsilon\)-multiplicative products exist iff the two Heegner orbits meet (\(\gamma z_Y = z_X\)), i.e. iff \(\beta_2(X) = \beta_1(Y)\) — e.g. the \((A_7, A_7)\)-cell attains its floor \(n = 97\) because both points are the elliptic point \(\rho\).

## 5. Verification

[scripts/middle_kernel.py](scripts/middle_kernel.py) (exact arithmetic via [scripts/omega.py](scripts/omega.py); mpmath at 60 digits; run inside the mpmath venv): the master formula and odd-integrality on 120 random \(\Omega\)-products (45 digits); \(\varepsilon_{ab+P} = \varepsilon_a\varepsilon_b\); the layer integral against quadrature (30 digits); complete orbit enumeration at Heegner pairs (coprime-pair ellipse method) with the counting law and integrality of all pairing values; the rank-one comparison and the class-dependence witness. **654 checks, all passing.**

## 6. Outlook

1. **Weight-\(m\) kernels: making \(S_f\) see the phase.** The radial class is blind to \(u\) — by design. The next class up: \(f\) bi-invariant but equivariant along the fiber circle of [moduli-invariants.md](moduli-invariants.md) §1 (the rotation in which \(\arg\Theta\) moves linearly). Fourier-decomposing in that angle should produce *weight-\(2m\) automorphic kernels* at the Heegner pair, with \(S_f\) picking up factors \(e^{im\arg u_X}\), \(e^{-im\arg u_Y}\) and values of weight-\(2m\) cusp forms at CM points — i.e. a kernel-theoretic home for the phase invariant, and a route from \(u_f\) to holomorphic-form periods.
2. **The RTF connection.** Integrating \(S_f(X, Y)\) over the Schmidt classes at fixed levels (with the \(H(n^2-1)\)-multiplicities) is exactly a geometric-side assembly of the relative trace formula of [spectral-outlook.md](spectral-outlook.md) §1; the present note is its "fixed double coset" fiber. The class-average statement (c) of §2 is the first quantitative instance.
3. **Green's values program.** Compute \(S_f\) for the higher-Green kernel at small Schmidt pairs and PSLQ against \(\log\) of algebraic numbers (Li's theorem gives the target field); compare with the Gross–Zagier-prime structure of the \(u\)-denominators ([moduli-invariants.md](moduli-invariants.md) §5.7). A match would fuse the phase arithmetic and the kernel arithmetic into one height-theoretic statement.
4. **Effective class-dependence.** The 18–38% class-swings observed at \((9,9)\) are Maass values at Heegner points; making outlook item (3) of [atomic-census.md](atomic-census.md) §8 quantitative (the deep-stratum density \(1 - 3/\pi\)) plausibly needs exactly these Weyl sums — the two error analyses should be run together.
