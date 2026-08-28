# The distinguished surface and its orthospectrum: Schmidt circles as geodesic planes in the Bianchi orbifold

Every companion document of this project lives on the boundary \(\hat{\mathbb{C}}\) of hyperbolic \(3\)-space. This one moves the whole picture *inside* \(\mathcal{H}^3\) and into the quotient
$$
M \;=\; \mathrm{PSL}_2(\mathbb{Z}[i]) \backslash \mathcal{H}^3,
$$
the **Bianchi orbifold** of the Gaussian integers. The organizing object is a single distinguished surface: the totally geodesic plane of \(\mathcal{H}^3\) lying **above the real line**. Its \(\Gamma\)-orbit has the Schmidt arrangement as boundary shadow, and under this dictionary the invariants of the companion documents become metric data of \(M\):

| boundary object (companion docs) | interior object (this doc) |
|---|---|
| circle \(\omega \in \mathcal{S}\) | totally geodesic plane (hemisphere) \(P_\omega \subset \mathcal{H}^3\); a sheet of the immersed surface \(Y \subset M\) |
| the invariant \(\alpha(\omega) = n\) ([hyperbolic-counting.md](../hyperbolic-counting.md)) | \(\cosh\) of the distance between \(P_{\hat{\mathbb{R}}}\) and \(P_\omega\) |
| the hyperbolic center \(m_1\) (a CM point of discriminant \(1-n^2\)) | the **foot** of the common perpendicular of the two planes, as a point of the surface \(Y\) |
| the Cartan matrix \(Z = X\bar X^{-1}\) ([involution.md](../involution.md)) | the holonomy of the closed geodesic doubling that perpendicular: a loxodromic with axis \((m_1, \bar m_1)\) and length \(2\operatorname{arccosh} n = 2\log\varepsilon_n\) |
| the triangle count \(3H(n^2-1)\) ([hyperbolic-counting.md](../hyperbolic-counting.md)) | the multiplicity of \(\operatorname{arccosh} n\) in the **ortho-length spectrum** of \(Y\) |
| level polynomials \(Q_n\), Hilbert class polynomials \(H_{1-n^2}\), traces \(t(n^2-1)\) ([moduli-invariants.md](../moduli-invariants.md)) | spectral decorations of the ortho-length strata |
| the atoms of \(\Omega\) ([half-plane-monoid.md](../half-plane-monoid.md)) | the boundary planes of an infinite-covolume (Apollonian) sub-quotient; counting exponent \(\delta = 1.30568\ldots\) |

All finite statements below are machine-verified — most of them in **exact arithmetic** — by [scripts/orthospectrum_verify.py](../scripts/orthospectrum_verify.py) (8/8 checks pass). The companion document [eisenstein-catalan.md](eisenstein-catalan.md) treats the analytic side (Eisenstein series, scattering, and the volume-theoretic origin of the Catalan constant).

## 1. The plane above the real line

Work in the upper half-space model \(\mathcal{H}^3 = \{(z, t): z \in \mathbb{C},\ t > 0\}\), \(ds^2 = (|dz|^2 + dt^2)/t^2\), on which \(\Gamma := \mathrm{PSL}_2(\mathbb{Z}[i])\) acts by orientation-preserving isometries extending the Möbius action on \(\hat{\mathbb{C}}\). Totally geodesic planes of \(\mathcal{H}^3\) are the vertical half-planes and hemispheres; each is determined by its boundary circle or line in \(\hat{\mathbb{C}}\). Write \(P_\mathcal{C}\) for the plane over a circle \(\mathcal{C}\). The **distinguished plane** is
$$
P_{\hat{\mathbb{R}}} \;=\; \{(x, 0, t) : x \in \mathbb{R},\ t > 0\},
$$
the vertical plane above the real line. Since the Schmidt arrangement is by definition \(\mathcal{S} = \Gamma\cdot\hat{\mathbb{R}}\):

> The map \(\gamma \mapsto \gamma(P_{\hat{\mathbb{R}}}) = P_{\gamma(\hat{\mathbb{R}})}\) is a bijection from \(\Gamma/\Gamma_H\) onto the planes over the circles and lines of \(\mathcal{S}\): **the Schmidt arrangement is the boundary shadow of one \(\Gamma\)-orbit of totally geodesic planes.**

**The stabilizer.** \(\Gamma_H := \operatorname{Stab}_\Gamma(P_{\hat{\mathbb{R}}}) \cong \mathrm{PGL}_2(\mathbb{Z})\). *Proof.* A Möbius transformation preserving \(\hat{\mathbb{R}}\) is represented by a real matrix \(g\), so an element of \(\Gamma_H\) is a class \(\lambda g\) with \(g\) real and \(\lambda^2 \det g = 1\). If \(\det g = 1\) then \(\lambda = \pm 1\) and \(\lambda g \in \mathrm{SL}_2(\mathbb{Z}[i]) \cap \mathrm{GL}_2(\mathbb{R}) = \mathrm{SL}_2(\mathbb{Z})\). If \(\det g = -1\) then \(\lambda = \pm i\), and \(ig\) has entries in \(i\mathbb{Z} \subset \mathbb{Z}[i]\): these classes do occur (e.g. \(z \mapsto -z\), the class of \(\operatorname{diag}(i, -i)\)). So \(\Gamma_H = \mathrm{PSL}_2(\mathbb{Z}) \sqcup i\!\cdot\!\mathrm{GL}_2(\mathbb{Z})^{\det = -1} \cong \mathrm{PGL}_2(\mathbb{Z})\). \(\square\)

Note what the two components do to \(\mathcal{H}^3 \cong P_{\hat{\mathbb{R}}} \times \mathbb{R}\) (normal coordinates): \(\mathrm{PSL}_2(\mathbb{Z})\) preserves each side of the plane and acts on it in the usual way; the second component acts on the plane by orientation-*reversing* isometries and **swaps the two sides** (e.g. \(z \mapsto -z\) exchanges the upper and lower half-planes of \(\mathbb{C}\)) — it must, since the total action on \(\mathcal{H}^3\) is orientation-preserving.

**The immersed surface.** The image of \(P_{\hat{\mathbb{R}}}\) in \(M\) is an immersed totally geodesic \(2\)-orbifold
$$
Y \;=\; \Gamma_H \backslash P_{\hat{\mathbb{R}}} \;\cong\; \mathrm{PGL}_2(\mathbb{Z})\backslash \mathbb{H}^2, \qquad \operatorname{area}(Y) = \frac{\pi}{6}
$$
(half of \(\operatorname{area}(\mathrm{PSL}_2(\mathbb{Z})\backslash\mathbb{H}^2) = \pi/3\)): the modular orbifold folded along its mirror symmetry — a \((2, 3, \infty)\) triangle with mirror boundary. It is non-compact: it runs up the cusp of \(M\).

**The canonical identification.** The teammates' "hyperbolic plane" — the upper half-plane \(\{\operatorname{Im} w > 0\} \subset \mathbb{C}\) of [hyperbolic-counting.md](../hyperbolic-counting.md), where the hyperbolic centers, the ideal triangle and the CM points live — sits on the *boundary* of \(\mathcal{H}^3\). It is canonically the same surface as \(P_{\hat{\mathbb{R}}}\):
$$
\iota:\ \{\operatorname{Im} w > 0\} \longrightarrow P_{\hat{\mathbb{R}}}, \qquad
\iota(w) = (\operatorname{Re} w,\ 0,\ \operatorname{Im} w),
$$
an isometry which is \(\mathrm{PSL}_2(\mathbb{R})\)-equivariant (a real Möbius matrix acts on quaternions \(x + tj\) with real \(x\) exactly as on the complex number \(x + it\)). Geometrically, \(\iota(w)\) is the point where the geodesic of \(\mathcal{H}^3\) with ideal endpoints \(w, \bar w\) crosses \(P_{\hat{\mathbb{R}}}\). **All hyperbolic-plane statements of the companion documents are therefore statements about the distinguished surface \(Y\) itself.**

**The second orbit, and orthogonal crossings.** The planes over the rotated family \(i\mathcal{S}\) (the \(x\)-odd–\(y\)-even circles and the vertical lines, cf. [circle-classification.md](../circle-classification.md)) form a *second* \(\Gamma\)-orbit, \(\Gamma \cdot P_{i\hat{\mathbb{R}}}\), with conjugate stabilizer; its image is a second immersed copy \(Y' \subset M\) of the same orbifold. The parity computation of [half-plane-monoid.md](../half-plane-monoid.md) (Prop. 5) now reads geometrically:

- within each family the inversive product is an **odd integer**, so two sheets of \(Y\) (or of \(Y'\)) are equal, tangent, or disjoint — the lifts of \(Y\) form a laminar family, never crossing transversally;
- between the families the inversive product is an **even integer**, so a sheet of \(Y\) and a sheet of \(Y'\) either avoid each other (\(|\langle\cdot,\cdot\rangle| \geq 2\)) or cross **orthogonally** (\(\langle\cdot,\cdot\rangle = 0\); inversive product zero means perpendicular circles). Example: the planes over \(\operatorname{Im} z = 0\) and \(\operatorname{Re} z = 0\) meet at right angles along the \(t\)-axis, a cusp-to-cusp geodesic of \(M\).

So \(M\) carries exactly two \(\Gamma\)-classes of "Schmidt planes"; each immersed surface is self-tangent only, and the two cross each other only at right angles, along geodesics with both endpoints in the cusp.

## 2. \(\alpha\) is the ortho-distance; the CM point is the foot

Fix a level-\(n\) circle \(\omega\) (curvature \(2q\), center \(\frac{x + ni}{2q}\), \(x^2 + n^2 - 1 = 4qm\), \(N := n^2 - 1\)), and recall its hyperbolic center \(m_1 = \frac{x + i\sqrt N}{2q}\), the CM point of discriminant \(1 - n^2\) attached to the form \((q, -x, m)\). [half-plane-monoid.md](../half-plane-monoid.md) §7 already identified \(\alpha = \langle M_0, M_\omega\rangle\) with \(\cosh\) of the distance in \(\mathcal{H}^3\) between the planes over \(\hat{\mathbb{R}}\) and over \(\omega\). Here is the perpendicular itself.

> **Proposition (the orthogeodesic).** Let \(\gamma^*\) be the geodesic of \(\mathcal{H}^3\) with ideal endpoints \(m_1\) and \(\bar m_1\) — the semicircle \(\{(\tfrac{x}{2q} + iy,\ t) : y^2 + t^2 = \tfrac{N}{4q^2}\}\). Then \(\gamma^*\) is the common perpendicular of \(P_{\hat{\mathbb{R}}}\) and \(P_\omega\). Its feet are
> $$
> F_1 = \Bigl(\frac{x}{2q},\ 0,\ \frac{\sqrt N}{2q}\Bigr) = \iota(m_1) \in P_{\hat{\mathbb{R}}},
> \qquad
> F_2 = \Bigl(\frac{x}{2q},\ \frac{N}{2qn},\ \frac{\sqrt N}{2qn}\Bigr) \in P_\omega,
> $$
> and the ortho-distance is
> $$
> \boxed{\ \cosh d(F_1, F_2) \;=\; n \;=\; \alpha(\omega).\ }
> $$
> In particular **the CM point \(m_1\) is precisely the foot of the orthogeodesic on the distinguished surface**, via the canonical identification \(\iota\) of §1.

*Proof.* Orthogonality by symmetry, twice. (i) The reflection \(z \mapsto \bar z\) of \(\hat{\mathbb{C}}\) extends to the reflection of \(\mathcal{H}^3\) through \(P_{\hat{\mathbb{R}}}\); it swaps the endpoints \(m_1 \leftrightarrow \bar m_1\), hence preserves \(\gamma^*\), and \(\gamma^*\) meets the mirror (at \(F_1\), where \(y = 0\)); an invariant geodesic through the mirror and not contained in it is orthogonal to it. (ii) The anti-Möbius inversion \(\tau_\omega\) in \(\omega\) satisfies \(\tau_\omega(z) = Z(\bar z)\) with \(Z = X\bar X^{-1}\) ([involution.md](../involution.md) §2), and \(Z\) fixes both \(m_1\) and \(\bar m_1\) (§3 below), so \(\tau_\omega(m_1) = Z(\bar m_1) = \bar m_1\) and \(\tau_\omega(\bar m_1) = m_1\): the inversion in \(\omega\) also preserves \(\gamma^*\), swapping its endpoints. And \(\gamma^*\) genuinely crosses the hemisphere: \(m_1\) lies inside the disk of \(\omega\) (distance \(\frac{n - \sqrt N}{2q} < \frac{1}{2q}\) from the center, as \(n - \sqrt{n^2-1} = \varepsilon_n^{-1} < 1\)) while \(\bar m_1\) lies outside. A geodesic invariant under the reflection in a sphere and crossing it is orthogonal to it.

Feet and distance by direct computation. Intersecting \(y^2 + t^2 = \frac{N}{4q^2}\) with the hemisphere \((y - \tfrac{n}{2q})^2 + t^2 = \tfrac{1}{4q^2}\) (both at \(\operatorname{Re} z = \tfrac{x}{2q}\)) gives \(y = \tfrac{N}{2qn}\), \(t = \tfrac{\sqrt N}{2qn}\), which is \(F_2\); \(F_1\) is the \(y = 0\) point. With \(\cosh d = 1 + \frac{|\Delta z|^2 + \Delta t^2}{2 t_1 t_2}\):
$$
|\Delta z|^2 + \Delta t^2
= \frac{N^2}{4q^2n^2} + \frac{N(n-1)^2}{4q^2n^2}
= \frac{N\,(n^2 - 1 + n^2 - 2n + 1)}{4q^2n^2}
= \frac{N(n-1)}{2q^2 n},
\qquad
2t_1t_2 = \frac{N}{2q^2n},
$$
so \(\cosh d = 1 + (n - 1) = n\). (One can also check the orthogonality at \(F_2\) in coordinates: the tangent of \(\gamma^*\) there is proportional to \((0, t, -y) \propto (0, \sqrt N, -N)\), the Euclidean radial direction of the hemisphere is \((0, y - \tfrac{n}{2q}, t) \propto (0, -1, \sqrt N)\), and these are parallel; since the model is conformal, Euclidean and hyperbolic angles agree.) \(\square\)

Machine verification: exact (40-digit) checks of the feet, the distance, the orthogonality, and \(|\langle M_0, M_\omega\rangle| = n\), for 42 circles across \(n \in \{3,4,5,7,9,12,15\}\) — including the even levels, where \(P_\omega\) is a sheet of \(Y'\).

**The tangency stratum.** \(\alpha = 1\) (the Ford circles) is the degenerate case: \(P_\omega\) is tangent to \(P_{\hat{\mathbb{R}}}\) at a rational boundary point, the perpendicular has length \(\operatorname{arccosh} 1 = 0\), and the "arc" escapes into the cusp. The ortho-length spectrum below starts at \(n = 2\).

## 3. The Cartan matrix is the holonomy of the orthogeodesic

[involution.md](../involution.md) proved that \(Z = X \bar X^{-1}\) (for \(X \in \mathrm{SL}_2(\mathbb{Z}[i])\) with \(X(\hat{\mathbb{R}}) = \omega\)) is, as an isometry of \(\mathcal{H}^3\), the composition of the inversions in \(P_\omega\) and in \(P_{\hat{\mathbb{R}}}\), with \(\operatorname{tr} Z = \pm 2n\). The product of the reflections in two planes at distance \(d\) is the hyperbolic translation of length \(2d\) along their common perpendicular; combining with §2:

> **Proposition (exact holonomy).** Let \(\varepsilon_n = n + \sqrt{n^2 - 1}\) and \(m_1\) as above. Then, exactly in the biquadratic field \(\mathbb{Q}(i, \sqrt{N})\),
> $$
> Z \binom{m_1}{1} = \lambda \binom{m_1}{1}, \qquad \lambda = \pm\varepsilon_n^{\pm 1},
> $$
> so the axis of the loxodromic \(Z \in \Gamma\) is exactly \(\gamma^* = (m_1, \bar m_1)\), the orthogeodesic of §2, with translation length
> $$
> \ell_n = 2\log \varepsilon_n = 2\operatorname{arccosh} n
> $$
> and characteristic ("length") polynomial \(\lambda^2 \mp 2n\lambda + 1\). The invariant trace field of the geodesic is the **real** quadratic field \(\mathbb{Q}(\sqrt{n^2-1})\).

This is Lemma A of [moduli-invariants.md](../moduli-invariants.md) §5 read geometrically; [scripts/orthospectrum_verify.py](../scripts/orthospectrum_verify.py) re-verifies it independently, in exact arithmetic, for **all 160 circles of all odd levels \(n \le 15\)**: the matrix \(X\) is produced by the descent of [circle-classification.md](../circle-classification.md), the eigenvector identity is checked in \(\mathbb{Q}(i)[\sqrt N]\), and \(\lambda\) always comes out \(\pm(n \pm \sqrt N)\).

In \(M\), the picture is: the orthogeodesic arc from \(Y\) to \(Y\) (odd \(n\); to \(Y'\) for even \(n\)) of length \(\operatorname{arccosh} n\), and the **closed geodesic** obtained by doubling it across the two mirrors — the axis of \(Z\) — of length \(2\operatorname{arccosh} n\), which crosses the immersed surface orthogonally at both feet. We call these the **\(\sigma\)-geodesics** of \(M\), since they are exactly the closed geodesics produced by the Cartan/Galois involution \(\sigma(X) = \bar X^{-1}\) of [involution.md](../involution.md); they are the Bianchi analogue of Sarnak's reciprocal geodesics on the modular surface (geodesics carrying a distinguished involution), with the roles twisted: *real* lengths \(2\log\varepsilon_n\), *imaginary* quadratic class data.

Two bookkeeping facts from the companion documents transfer directly:

- distinct circles at level \(n\) give the same closed geodesic exactly along the fibers of the involution \(\hat\sigma[f] = [\mathfrak{r}_n][f]^{-1}\) ([involution.md](../involution.md) §4, Corollary), so closed \(\sigma\)-geodesics of length \(2\operatorname{arccosh} n\) are indexed by \(\hat\sigma\)-orbits of classes;
- the arcs themselves (based at the surface) are indexed by the classes with no folding — see the theorem below.

## 4. The ortho-length spectrum: multiplicities are Hurwitz class numbers

Cut \(M\) along the immersed surfaces and count the perpendiculars. An **ortho-arc at level \(n\)** is a \(\Gamma_H\)-orbit of common perpendiculars from \(P_{\hat{\mathbb{R}}}\) to a plane of the arrangement at inversive distance \(n\); by §2–3 these correspond to level-\(n\) circles modulo \(\Gamma_H\). Since \(z \mapsto -z \in \Gamma_H\) identifies the upper-half-plane circles with the lower ones, ortho-arcs correspond to upper level-\(n\) circles modulo \(\mathrm{PSL}_2(\mathbb{Z})\) — exactly the objects counted in [hyperbolic-counting.md](../hyperbolic-counting.md). Its theorem ("weighted count in the ideal triangle \(= 3H(n^2-1)\)", where the triangle is three copies of the modular fundamental domain) becomes:

> **Theorem (arithmetic ortho-length spectrum).** In the Bianchi orbifold \(M\), the ortho-length spectrum of the pair of immersed modular surfaces \((Y, Y')\) is
> $$
> \{\operatorname{arccosh} n \ :\ n = 2, 3, 4, \ldots\},
> $$
> with \(Y\!\leftrightarrow\!Y\) perpendiculars at odd \(n\) and \(Y\!\leftrightarrow\!Y'\) perpendiculars at even \(n\), and the level-\(n\) ortho-arcs based at the distinguished sheet, counted with orbifold weights \(1/|\mathrm{Stab}|\), have total mass
> $$
> \boxed{\;H(n^2 - 1)\;}
> $$
> — the Hurwitz class number of discriminant \(1 - n^2\). The feet of the level-\(n\) arcs on the surface are precisely the CM points of discriminant \(1 - n^2\), and the closed \(\sigma\)-geodesics of length \(2\operatorname{arccosh} n\) correspond to \(\hat\sigma\)-orbits \(\{[f], [\mathfrak{r}_n][f]^{-1}\}\) of form classes.

The multiplicity table is the \(H(n^2-1)\) row of [hyperbolic-counting.md](../hyperbolic-counting.md): \(\tfrac13, 1, 2, 2, 2, \tfrac{10}3, 5, 6, 3, 4, 10, \ldots\) for \(n = 2, 3, \ldots\) (the fractional weights at \(n = 2, 7, 26, \ldots\) are the arcs whose foot is the orbifold point of order 3 — the Pell levels \(n^2 - 3f^2 = 1\)). Machine check: the weighted count equals \(3H(n^2-1)\) per triangle, i.e. \(H(n^2-1)\) per fundamental domain, for all \(2 \le n \le 20\).

This is the precise form of outlook item 3.1 ("the arithmetic length spectrum of the Bianchi orbifold"): the striking feature promised there — **real quadratic lengths carrying imaginary quadratic class numbers** — is now a theorem, with the mechanism identified: the lengths are forced by the *distance to a fixed totally geodesic surface* (a real-quadratic quantity, \(\varepsilon_n = e^{\ell_n/2}\)), while the multiplicity counts the *feet on that surface* (CM points, an imaginary-quadratic count).

## 5. Growth, entropy, and equidistribution of the orthospectrum

**Quadratic growth of multiplicities.** Summing the theorem over levels (bulk Hurwitz sieve up to \(10^6\), cross-checked against the exact routine of [scripts/alpha_circles.py](../scripts/alpha_circles.py)):
$$
S_H(X) := \sum_{2 \le n \le X} H(n^2 - 1) \;\sim\; C\,X^2,
\qquad
\frac{S_H(X)}{X^2} = 0.2823,\ 0.2846,\ 0.2859,\ 0.2858 \quad (X = 125, 250, 500, 1000).
$$
For comparison, if the discriminants \(n^2 - 1\) behaved like random admissible discriminants of the same size, the average of \(H\) would give \(\pi/(12\zeta(3)) = 0.21778\ldots\); the measured constant exceeds this by a stable factor \(\approx 1.312\) — the **singular series of the family \(n^2 - 1 = (n-1)(n+1)\)**, whose exact evaluation we leave open (§8.1).

**Half the entropy.** The closed \(\sigma\)-geodesics of length \(\le L\) number
$$
\Pi_\sigma(L) \;=\; \sum_{2\operatorname{arccosh} n \le L} (\text{mult}) \;\asymp\; S_H(\cosh(L/2)) \;\asymp\; e^{L},
$$
and the sieve data confirms the exponential rate directly: \(d\log \Pi_\sigma/dL = 1.007,\ 1.003\) across \(L = 10 \to 12 \to 14\). Against the prime geodesic theorem for \(\mathrm{PSL}_2(\mathbb{Z}[i])\) — all closed geodesics number \(\sim e^{2L}/(2L)\), topological entropy \(2\) for \(\mathcal{H}^3\) (see e.g. *The Prime Geodesic Theorem for \(\mathrm{PSL}_2(\mathbb{Z}[i])\) and Spectral Exponential Sums*, arXiv:1903.05111) — the \(\sigma\)-family is an exponentially rich but exponentially thin family of **exactly half the entropy**:
$$
\Pi_\sigma(L) \asymp e^{(h_{\mathrm{top}}/2)\,L}, \qquad h_{\mathrm{top}} = 2.
$$
This is the precise 3-dimensional analogue of the density of reciprocal geodesics on the modular surface (Sarnak: \(\asymp e^{L/2}\) against \(e^L\) — again the square root), as befits geodesics forced to hit a fixed totally geodesic object: one crossing condition costs half the entropy.

**Equidistribution of the feet.** The level-\(n\) feet are the CM points of discriminant \(1 - n^2\) on \(Y\). By Duke's theorem (equidistribution of Heegner points; for the non-fundamental discriminants occurring here one needs the extensions via subconvexity — literature pass before citing in print, cf. [outlook.md](../outlook.md) §4), the feet equidistribute on \(Y\) with respect to hyperbolic area as \(n \to \infty\): **the orthospectrum of the Bianchi orbifold rains down uniformly on the distinguished surface.** Quantitatively, level \(n\) drops \(\asymp n\) feet of arcs of length \(\operatorname{arccosh} n \approx \log 2n\) — the surface is struck ever more densely by ever longer perpendiculars.

## 6. Spectral decorations of the strata

Each stratum \(\{\text{ortho-length } \operatorname{arccosh} n\}\) of the orthospectrum carries, from the companion documents, a stack of canonical invariants — all indexed by the arcs of that stratum:

1. **The length polynomial** \(\lambda^2 - 2n\lambda + 1\): the characteristic polynomial of the holonomy; roots \(\varepsilon_n^{\pm 1}\); one per stratum.
2. **The Hilbert class polynomial** \(H_{1-n^2}(x) = \prod_{\text{feet}} (x - j(\text{foot}))\): the \(j\)-coordinates of the (primitive-class) feet, an integer polynomial per stratum ([moduli-invariants.md](../moduli-invariants.md) §2).
3. **The Zagier trace** \(\mathrm{Tr}_n = \sum_{\text{feet}} \frac{j(\text{foot}) - 744}{w} = t(n^2 - 1)\): the total \(j\)-observable of the stratum; assembling over strata slices Zagier's weight-\(3/2\) form along \(d = n^2 - 1\) ([moduli-invariants.md](../moduli-invariants.md) §3).
4. **The phase units and the level polynomial**: each arc (class \(f\)) carries the unit \(u_f = \varepsilon_n \Theta_f\); the stratum carries the single integer polynomial \(Q_n\) of degree \(h(1-n^2)\) whose roots they are, irreducible at every computed level ([moduli-invariants.md](../moduli-invariants.md) §5.9–5.10). Via Theorem A there (dihedral Galois equivariance), irreducibility says: **the absolute Galois group permutes the ortho-arcs of a given length transitively** — a Linnik-flavored algebraic transitivity underneath the analytic equidistribution of §5.

The involution \(\hat\sigma\) acts on each stratum (twisted inversion \([\mathfrak{r}_n][f]^{-1}\)), pairing the arcs that double to the same closed geodesic; the phase laws \(u_{\hat\sigma f} = \bar u_f^{-1}\) and \(u_{\mathfrak{r} f} u_f = 1\) are then statements about how the decoration transforms along the geodesic's two crossings of the surface.

**The frame.** All of this is the *geometric side* of a relative trace formula. \((\mathrm{SL}_2(\mathbb{C}), \mathrm{SL}_2(\mathbb{R}))\) is a symmetric — indeed Gelfand — pair (as [involution.md](../involution.md) §5 identified), \(Y \subset M\) is the associated cycle, and the double cosets \(\Gamma_H\backslash\Gamma/\Gamma_H\), which are exactly the circle classes at the various levels, index the geometric (orbital-integral) side; the multiplicities \(H(n^2-1)\) are its arithmetic content. The *spectral side* is periods over \(Y\) of the automorphic spectrum of \(M\) — Eisenstein series and Bianchi–Maass forms — and is taken up in [eisenstein-catalan.md](eisenstein-catalan.md), where the first payoff (the Catalan constant of the Euclidean count as a volume ratio) is proved and machine-verified.

## 7. The infinite-covolume stratum: atoms

The atoms of the half-plane monoid \(\Omega\) ([half-plane-monoid.md](../half-plane-monoid.md), Thm. 9) are the planes over the strip Apollonian gasket — the boundary of an *infinite-covolume* quotient sitting inside the surface system. Counting them per period by exact Descartes recursion on the Hermitian triples (\(D_j' = 2(D_k + D_l + D_m) - D_j\), quadruples deduplicated modulo \(z \mapsto z + 1\)):
$$
N(K) = \#\{\text{atoms of curvature} \le K\}\text{ per period}: \quad
N(10^2, 10^3, 10^4, 3\cdot 10^4) = 49,\ 951,\ 19299,\ 80993,
$$
with fitted exponent \(\delta_{\mathrm{fit}} = 1.3033\) — against the Hausdorff dimension of the Apollonian gasket \(\delta = 1.30568672\ldots\) (McMullen's eigenvalue algorithm; high-precision value from the recent literature). This answers the growth part of Question 2 of [half-plane-monoid.md](../half-plane-monoid.md): **the atom count is governed not by the volume of \(M\) but by the base eigenvalue**
$$
\lambda_0 = \delta(2 - \delta) \approx 0.9066
$$
of the Laplacian on the infinite-covolume manifold of the Apollonian group (Patterson–Sullivan; the counting asymptotics \(N(K) \sim c\,K^\delta\) is Kontorovich–Oh). The contrast with the full arrangement is the lattice/thin-group dichotomy in one picture: the full circle count grows like \(K^2\) with a **volume-ratio** constant ([eisenstein-catalan.md](eisenstein-catalan.md) §3), the atom count like \(K^\delta\) with a **Patterson–Sullivan** constant, and both are counts of the *same* geometric objects, restricted to different sub-orbits.

## 8. Machine verification

`python3 scripts/orthospectrum_verify.py` — 8 checks, all passing:

1. **Exact holonomy** (§3): for all 160 circles of all odd levels \(n \le 15\): \(X\) from the descent, \(\operatorname{tr}(X\bar X^{-1}) = \pm 2n\) exactly, and the eigenvector identity \(Z(m_1, 1)^{\mathsf T} = \pm(n \pm \sqrt N)(m_1, 1)^{\mathsf T}\) verified in exact arithmetic in \(\mathbb{Q}(i, \sqrt N)\); translation length \(2\log\varepsilon_n = 2\operatorname{arccosh} n\).
2. **Orthogeodesic geometry** (§2): feet on the respective planes, \(\cosh d(F_1, F_2) = n\), orthogonality at both feet, \(|\langle M_0, M_\omega\rangle| = n\); 42 circles including even levels; 40 digits.
3. **Counts** (§4): weighted triangle count \(= 3H(n^2-1)\), \(2 \le n \le 20\).
4. **Growth and entropy** (§5): bulk Hurwitz sieve to \(10^6\) (cross-checked against `alpha_circles.hurwitz` on 25 random \(N\) and all \(N = n^2-1\), \(n \le 30\)); \(S_H(X)/X^2\) stable \(\approx 0.286\); entropy slopes \(\to 1\).
5. **Atoms** (§7): exact Descartes/strip-gasket enumeration (every generated disk verified to be a Schmidt disk in \(\mathbb{H}\)); \(N(2) = 1\), \(N(8) = 3\); fitted \(\delta = 1.3033 \approx 1.30568\).

## 9. Questions

1. **The singular series of \(n^2 - 1\).** Evaluate \(C = \lim S_H(X)/X^2\) in closed form (the measured \(C \approx 0.2858\), an enhancement \(\approx 1.312\) over the random-discriminant heuristic \(\pi/(12\zeta(3))\)). The natural route: average \(L(1, \chi_{1-n^2})\) over the family, or extract the diagonal from Eichler/Kronecker–Hurwitz relations (cf. [outlook.md](../outlook.md) 1.2) — which would simultaneously give the *exact* second-order terms of the ortho-length counting function.
2. **Primitivity.** Which \(\sigma\)-geodesics are primitive in \(\Gamma\) (not proper powers), and what is the exact relation between the arc count \(H(n^2-1)\) and the multiplicity of \(2\operatorname{arccosh} n\) in the honest (primitive) length spectrum of \(M\)? The \(\hat\sigma\)-fold is understood ([involution.md](../involution.md)); powers and coincidences of lengths across different \(n\) (\(\varepsilon_n\) powers, e.g. \(\varepsilon_{2n^2-1} = \varepsilon_n^2\)) are the remaining bookkeeping.
3. **The crossing geodesics.** The locus \(Y \cap Y'\) (right-angle crossings, §1) consists of cusp-to-cusp geodesics; classify them arithmetically (they are indexed by pairs of orthogonal circles \(\langle M, M'\rangle = 0\), one from each family) and fit them into the trace-formula frame as the "singular" orbital terms.
4. **Orthospectrum identities.** For manifolds with totally geodesic *boundary*, the Basmajian and Bridgeman–Kane identities express boundary area, respectively volume, as sums over the orthospectrum. Here \(Y\) is immersed, not boundary; is there an identity summing an explicit function of \(\{\operatorname{arccosh} n\}\) with multiplicities \(H(n^2-1)\) against \(\operatorname{area}(Y)\) or \(\operatorname{vol}(M)\)? The Fubini computations of [eisenstein-catalan.md](eisenstein-catalan.md) §3 are a smoothed version of exactly this and suggest the identity exists in a regularized form.
5. **Kronecker–Hurwitz, geometrically** (sharpening [outlook.md](../outlook.md) 3.3): the class-number relations now assert linear identities between multiplicities of *different ortho-lengths*; a bijective proof would be an orthospectrum recursion — plausibly a cut-and-paste of arcs across the tangency stratum \(\alpha = 1\).
