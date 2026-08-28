# Eisenstein series, the scattering matrix, and the Catalan volume: the spectral home of the Euclidean count

[euclidean-counting.md](../euclidean-counting.md) ended on a promissory note (Remark 1): Catalan's constant \(G = L(2, \chi_{-4})\) appeared there "purely arithmetically", while "the same constant is famously the hyperbolic-volume quantum of the Bianchi orbifold \(\ldots\) the agreement is of course not a coincidence and is worth revisiting". This document revisits it and closes it. The mechanism is the distinguished surface of [orthogeodesics.md](orthogeodesics.md): the counting function of the Schmidt arrangement is the shadow, along the cusp, of the orbit of the plane above the real line, and its leading constant is forced by equidistribution to be a **volume ratio**:

$$
\boxed{\;\sum_{n \le X} N_e(n) \;\sim\; \frac{X^2}{\pi}\cdot\frac{\operatorname{area}(Y)}{\operatorname{vol}(M)}
\;=\; \frac{X^2}{\pi}\cdot\frac{\pi/6}{G/3} \;=\; \frac{X^2}{2G}\;}
$$

— the elementary theorem of [euclidean-counting.md](../euclidean-counting.md), with its Catalan constant now identified as \(\operatorname{vol}(M) = G/3\) (Humbert). Catalan's constant enters the project in three *a priori* different ways, and the theorems below identify all three:

| appearance | where | mechanism |
|---|---|---|
| \(L(2, \chi_{-4})\) from local densities of \(x^2 + y^2 \equiv 1 \bmod 4n\) | [euclidean-counting.md](../euclidean-counting.md) | Euler product of \(N_e\) |
| \(\operatorname{vol}(\mathrm{PSL}_2(\mathbb{Z}[i])\backslash\mathcal{H}^3) = G/3\) | Humbert's formula | \(\zeta_K(2) = \tfrac{\pi^2}{6}G\) |
| \(\operatorname{Res}_{s=1}\varphi(s) = \tfrac{3}{2G} = \tfrac{1}{2\operatorname{vol}(M)}\), \(\varphi\) = scattering term | §1 below | Maass–Selberg |

Everything quantitative below is machine-verified by [scripts/spectral_verify.py](../scripts/spectral_verify.py) (23/23 checks pass), including a brute-force lattice-sum verification of the scattering formula to \(10^{-7}\).

## 1. The spectrum of the Bianchi orbifold: conventions and the scattering term

On \(\mathcal{H}^3 = \{(z, t)\}\) with \(ds^2 = (|dz|^2 + dt^2)/t^2\), \(dV = t^{-3}\,dx\,dy\,dt\), the Laplacian \(\Delta = t^2(\partial_x^2 + \partial_y^2 + \partial_t^2) - t\,\partial_t\) has \(-\Delta\, t^{1+s} = (1 - s^2)\,t^{1+s}\): spectral parameter \(s\), eigenvalue \(\lambda = 1 - s^2\), continuous spectrum \([1, \infty)\) (\(s \in i\mathbb{R}\)).

\(\Gamma = \mathrm{PSL}_2(\mathbb{Z}[i])\) has one cusp (class number one), with stabilizer \(\Gamma_\infty = \{z \mapsto z + \mu\} \sqcup \{z \mapsto -z + \mu\}\), \(\mu \in \mathbb{Z}[i]\): the cusp cross-section is the torus \(\mathbb{C}/\mathbb{Z}[i]\) folded by \(z \mapsto -z\), of area \(\tfrac12\). The Eisenstein series is
$$
E(P, s) \;=\; \sum_{\gamma \in \Gamma_\infty\backslash\Gamma} t(\gamma P)^{1+s},
\qquad
t(\gamma P) = \frac{t}{|cz + d|^2 + |c|^2 t^2},
$$
absolutely convergent for \(\operatorname{Re} s > 1\). Cosets \(\Gamma_\infty\backslash\Gamma\) correspond to coprime bottom rows \((c, d)\) modulo the units \(\mu_4\), and removing coprimality gives the summation-friendly form used by the verification script:
$$
E(P, s) \;=\; \frac{1}{4\,\zeta_K(1+s)} \sum_{(c,d) \in \mathbb{Z}[i]^2 \smallsetminus 0} \Bigl(\frac{t}{|cz+d|^2 + |c|^2t^2}\Bigr)^{1+s},
\qquad
\zeta_K(s) = \zeta(s)\,L(s, \chi_{-4})
$$
the Dedekind zeta function of \(K = \mathbb{Q}(i)\).

**The constant term.** Write \(a_0(t) = \int_{[0,1]^2} E((z,t), s)\,dz\). The \(c = 0\) cosets (\((0,1), (0,i)\) mod \(\pm\), i.e. one coset of \(\Gamma_\infty\)) give \(t^{1+s}\). For \(c \neq 0\), summing \(d\) over each coprime residue \(d_0 \bmod c\) and unfolding the \(z\)-integral over \(\mathbb{C}\):
$$
\int_{\mathbb{C}} \Bigl(\frac{t}{|c|^2(|w|^2 + t^2)}\Bigr)^{1+s} dw
= \frac{\pi}{s}\, t^{1-s}\, |c|^{-2(1+s)},
\qquad
\sum_{\text{ideals } \mathfrak{c} \neq 0} \frac{\varphi_K(\mathfrak{c})}{N\mathfrak{c}^{\,1+s}} = \frac{\zeta_K(s)}{\zeta_K(s+1)} .
$$

> **Scattering term.**
> $$
> a_0(t) = t^{1+s} + \varphi(s)\, t^{1-s}, \qquad
> \boxed{\;\varphi(s) = \frac{\pi}{s}\,\frac{\zeta_K(s)}{\zeta_K(s+1)}\;}
> $$
> (the non-constant Fourier modes are the usual exponentially decaying \(K\)-Bessel terms with Gaussian divisor-sum coefficients; Elstrodt–Grunewald–Mennicke, ch. 8). Machine check: brute-force lattice summation of \(E\) at \((s, t) = (2, 2)\) and \((1.5, 2.5)\), Fourier-averaged over a \(4\times4\) grid, reproduces \(\varphi(s)\) with relative error \(< 10^{-6}\).

**The Catalan chain.** Three classical evaluations line up:
$$
\zeta_K(2) = \zeta(2) L(2, \chi_{-4}) = \frac{\pi^2}{6}\,G;
\qquad
\operatorname{vol}(M) = \frac{|d_K|^{3/2}\zeta_K(2)}{4\pi^2} = \frac{8 \cdot \frac{\pi^2}{6}G}{4\pi^2} = \frac{G}{3}
$$
(Humbert; verified both symbolically and by direct numerical integration of the Elstrodt–Grunewald–Mennicke fundamental domain \(\{|{\operatorname{Re} z}| \le \tfrac12,\ 0 \le \operatorname{Im} z \le \tfrac12,\ |z|^2 + t^2 \ge 1\}\), which matches \(G/3 = 0.30532186\ldots\) to 14 digits); and, with \(\operatorname{Res}_{s=1}\zeta_K = \tfrac{\pi}{4}\) (class number formula),
$$
\operatorname{Res}_{s=1} \varphi(s) \;=\; \frac{\pi \cdot \frac\pi4}{\zeta_K(2)} \;=\; \frac{3}{2G}
\;=\; \frac{\text{area of the cusp cross-section}}{\operatorname{vol}(M)} \;=\; \frac{1/2}{G/3},
$$
the Maass–Selberg relation: \(E(P, s)\) has its pole at \(s = 1\) with residue the **constant function** \(3/(2G)\). This is the spectral face of the Catalan quantum, and it is the residue that will match the counting in §4.

## 2. The counting Dirichlet series is a zeta ratio of scattering type

[euclidean-counting.md](../euclidean-counting.md) computed \(\sum N_e(n)\,n^{-s} = \zeta(s-1)/L(s, \chi_{-4})\). Multiplying numerator and denominator by \(\zeta(s)\):

> **Lemma 1.**
> $$
> \sum_{n \ge 1} \frac{N_e(n)}{n^{s}} \;=\; \frac{\zeta(s-1)\,\zeta(s)}{\zeta_K(s)} .
> $$

The Dedekind zeta function of \(\mathbb{Z}[i]\) — the arithmetic engine of the scattering matrix — sits in the **denominator** of the counting series; the numerator is rational-integer zeta data. This is precisely the fingerprint of a \(\mathbb{Q}\)-cycle inside a \(\mathbb{Q}(i)\)-manifold (an Asai-type ratio; see §4), and the pole at \(s = 2\), of residue \(\zeta(2)/\zeta_K(2) = 1/G\), is the counting constant. (Sieve check at \(s = 3, 4\) to \(10^{-5}\), plus the mean-value check \((2G/X^2)\sum_{n \le 10^6} N_e(n) = 1.000001\).)

## 3. The plane-orbit kernel: two evaluations of one integral, and the closure of Remark 1

Attach to the orbit of the distinguished plane the **incidence kernel**
$$
\Theta_s(P) \;=\; \sum_{\mathfrak{p}\, \in\, \Gamma\cdot P_{\hat{\mathbb{R}}}} \bigl(\cosh \operatorname{dist}(P, \mathfrak{p})\bigr)^{-(2+s)}
\qquad (s > 0),
$$
a \(\Gamma\)-invariant function on \(M\) summing a decaying weight of the distance to every sheet of the immersed surface \(Y\). Two exact computations, one along the cusp and one over the fundamental domain, evaluate its mass in unrelated-looking ways; equating them is the theorem. Throughout, \(C_s := \int_0^\infty \cosh^{-s}u\,du = \tfrac{\sqrt\pi}{2}\Gamma(\tfrac s2)/\Gamma(\tfrac{s+1}2)\) and \(J_s(v) := \int_v^\infty (1 + u^2)^{-(2+s)/2}\,du\).

**(a) The cusp expansion is the circle count.** The signed distance from \(P = (z, t)\) to the plane over a circle (center \(z_0\), radius \(r\)) is \(\sinh u = \frac{|z - z_0|^2 + t^2 - r^2}{2rt}\); to the plane over the line \(\operatorname{Im} z = k\) it is \(\sinh u = \frac{\operatorname{Im} z - k}{t}\). Integrating over one period \(z \in [0,1]^2\) and using the classification (per unit area, \(N_e(n)\) circles of curvature \(2n\); one line per integer height):

> **Proposition A (exact constant term).** For every \(t > 0\),
> $$
> \operatorname{CT}(t) := \int_{[0,1]^2} \Theta_s((z,t))\,dz
> \;=\; \underbrace{\frac{\sqrt{\pi}\,\Gamma(\tfrac{1+s}{2})}{\Gamma(\tfrac{2+s}{2})}\; t}_{\text{lines}}
> \;+\; \underbrace{\pi t \sum_{n \ge 1} \frac{N_e(n)}{n}\, J_s\!\Bigl(nt - \frac{1}{4nt}\Bigr)}_{\text{hemispheres}} .
> $$
> *Proof.* Positivity permits term-by-term integration. For one hemisphere, substituting \(w = |z - z_0|^2\) and then \(v = \frac{w + t^2 - r^2}{2rt}\) gives \(\int_{\mathbb{C}}(\cosh u)^{-(2+s)}dz = 2\pi r t\, J_s(v_0)\) with \(v_0 = \frac{t^2 - r^2}{2rt}\), and \(r = \frac{1}{2n}\) makes \(v_0 = nt - \frac{1}{4nt}\). For the lines, \(\int_{\mathbb{R}} (1 + (y/t)^2)^{-(2+s)/2} dy = t\, B(\tfrac12, \tfrac{1+s}{2})\). \(\square\)

**(b) The global mass is an area-to-volume ratio.** Unfold over \(\Gamma_H = \operatorname{Stab}(P_{\hat{\mathbb{R}}})\) and use Fermi coordinates off the plane (\(dV = \cosh^2 u \, dA\, du\)), folding by the side-swapping component of \(\Gamma_H\):
$$
\int_{\Gamma\backslash\mathcal{H}^3} \Theta_s \, dV
= \int_{\Gamma_H\backslash\mathcal{H}^3} (\cosh u)^{-(2+s)}\,dV
= \operatorname{area}\bigl(\mathrm{PSL}_2(\mathbb{Z})\backslash\mathbb{H}^2\bigr)\int_0^\infty \cosh^{-s}u\,du
= \frac{\pi}{3}\,C_s,
$$
so the average of \(\Theta_s\) over \(M\) is
$$
\langle \Theta_s \rangle = \frac{(\pi/3)\,C_s}{\operatorname{vol}(M)} = \frac{\pi C_s}{G}
\;=\; 2\,C_s\,\frac{\operatorname{area}(Y)}{\operatorname{vol}(M)}
\qquad\Bigl(\operatorname{area}(Y) = \frac{\pi}{6},\ \operatorname{vol}(M) = \frac{G}{3}\Bigr).
$$

**(c) The hinge.** The horospheres \(t = \mathrm{const}\) project to expanding closed horospheres in \(M\) as \(t \to 0\), and \(\operatorname{CT}(t)\) is exactly the horosphere average of \(\Theta_s\). Equidistribution of expanding closed horospheres in finite-volume hyperbolic manifolds (mixing; Eskin–McMullen) predicts \(\operatorname{CT}(t) \to \langle\Theta_s\rangle\). On the other side, Proposition A turns this limit into a statement about the mean of \(N_e(n)/n\). The two are reconciled by a clean Fubini identity:

> **Lemma 2 (kinematic lemma).** For every \(s > 0\),
> $$
> \int_0^\infty J_s\!\Bigl(x - \frac{1}{4x}\Bigr) dx \;=\; C_s .
> $$
> *Proof.* Substitute \(x = \tfrac12 e^\theta\), so \(x - \tfrac{1}{4x} = \sinh\theta\), and \(J_s(\sinh\theta) = \int_\theta^\infty \cosh^{-(1+s)}\varphi\, d\varphi\) (via \(u = \sinh\varphi\)). Then
> $$
> \int_{\mathbb{R}} J_s(\sinh\theta)\,\frac{e^\theta}{2}\,d\theta
> = \frac12 \int_{\mathbb{R}} \cosh^{-(1+s)}\varphi \Bigl(\int_{-\infty}^{\varphi} e^\theta d\theta\Bigr) d\varphi
> = \frac12 \int_{\mathbb{R}} \cosh^{-(1+s)}\varphi\;(\cosh\varphi + \sinh\varphi)\, d\varphi = C_s. \; \square
> $$

> **Proposition B (unconditional horosphere limit).** \(\displaystyle \lim_{t \to 0^+} \operatorname{CT}(t) = \frac{\pi\,C_s}{G}\).
> *Proof.* The line term is \(O(t)\). For the hemisphere term, [euclidean-counting.md](../euclidean-counting.md) gives \(\sum_{n \le X} N_e(n)/n = X/G + O(\log^2 X)\) by partial summation from its mean-count theorem; summing \(\frac{N_e(n)}{n}\, h(nt)\) against the \(C^1\), integrably-decaying kernel \(h(x) = J_s(x - \frac{1}{4x})\) by partial summation gives \(\pi t \cdot \frac1G \cdot \frac1t \int_0^\infty h + o(1) \to \frac{\pi}{G}\,C_s\) by Lemma 2. \(\square\)

Numerics (both \(s = 2\), target \(\pi/G = 3.429815\ldots\), and \(s = 1\), target \(\pi^2/2G = 5.387541\ldots\)): \(\operatorname{CT}(t)\) at \(t = 0.2, 0.1, 0.05, 0.02\) marches to the target with visibly linear-in-\(t\) error. Notice \(\pi/G\) is also the mean total circumference per unit area per level found in [euclidean-counting.md](../euclidean-counting.md) — the \(t \to 0\) horosphere sees each plane with weight proportional to its boundary length, which is that statement in kernel form.

Proposition B is proved *from* the elementary count; running the logic in the other direction, the horosphere equidistribution theorem plus Lemma 2 forces the mean of \(N_e(n)/n\) to be \(\langle\Theta_s\rangle \cdot G/\pi \cdot \ldots\) — i.e. forces the counting constant. Either way the constant is identified:

> **Theorem (the Catalan volume identity; Remark 1 of [euclidean-counting.md](../euclidean-counting.md) closed).** The Euclidean counting constant of the Gaussian Schmidt arrangement is a volume ratio:
> $$
> \sum_{n \le X} N_e(n) \;=\; \frac{X^2}{2G} + O(X\log X), \qquad
> \frac{1}{2G} \;=\; \frac{1}{\pi}\cdot \frac{\operatorname{area}(Y)}{\operatorname{vol}(M)} .
> $$
> Equivalently: the number of Schmidt circles of curvature \(\le \kappa\) per unit area is \(\sim \frac{\kappa^2}{4\pi}\cdot\frac{\operatorname{area}(Y)}{\operatorname{vol}(M)}\), and in the Duke–Rudnick–Sarnak/Eskin–McMullen normalization (invariant measure \(dz_0\,dr/r^3\) on the space of circles) the orbit-counting constant is \(\frac{\operatorname{area}(Y)}{2\pi\operatorname{vol}(M)} = \frac{1}{4G}\). Catalan's constant appears **because** \(\operatorname{vol}(M) = G/3\): the local solution densities of \(x^2 + y^2 \equiv 1 \bmod 4n\) and Humbert's volume integral are two computations of the same Tamagawa-type quantity, matched through the surface \(Y\).

(The \(\mathrm{PGL}_2\)-variant doubles consistently: the arrangement \(\mathcal{S} \sqcup i\mathcal{S}\) doubles the plane-orbit — both immersed surfaces \(Y, Y'\) of [orthogeodesics.md](orthogeodesics.md) now contribute — and \(N_e^{\mathrm{PGL}} = 2N_e\).)

## 4. The Mellin transform along the cusp: a period identity of Asai type

The whole \(t\)-dependence of §3(a) compresses into Gamma factors. Let
$$
K_s(w) := \int_0^\infty J_s\!\Bigl(x - \frac{1}{4x}\Bigr)\,x^{w-1}\,dx
\;=\; \frac{2^{\,s-w}\,\Gamma\bigl(\tfrac{1+s+w}{2}\bigr)\Gamma\bigl(\tfrac{1+s-w}{2}\bigr)}{w\,\Gamma(1+s)}
\qquad (0 < \operatorname{Re} w < 1 + s),
$$
by the same substitution as Lemma 2 followed by parts and Euler's Beta integral; \(K_s(1) = C_s\) is the duplication formula (Lemma 2 again).

> **Theorem (Mellin/period identity).** For \(1 < \operatorname{Re} w < 1 + s\), with \(\operatorname{Hem}_s(t)\) the hemisphere part of Proposition A,
> $$
> \int_0^\infty \operatorname{Hem}_s(t)\; t^{w-2}\, dt
> \;=\; \pi\, K_s(w)\; \frac{\zeta(w)\,\zeta(w+1)}{\zeta_K(w+1)} .
> $$
> *Proof.* Tonelli termwise; in the \(n\)-th term substitute \(t = x/n\):
> \(\pi \sum_n \frac{N_e(n)}{n}\, n^{-w} \int_0^\infty J_s(x - \tfrac{1}{4x})\, x^{w-1} dx = \pi K_s(w) \sum_n N_e(n)\, n^{-(w+1)}\), and apply Lemma 1 at \(w + 1\). \(\square\)
>
> Machine check: direct numerical \(t\)-integration of the exact series at \((s, w) = (2, 2.5)\) and \((2, 1.5)\) against the closed form (relative error \(3\cdot10^{-7}\), \(6\cdot10^{-5}\)).

**Interpretation.** By Rankin–Selberg unfolding in Zagier's regularized sense, the Mellin transform of a \(\Gamma\)-invariant function's constant term *is* its pairing with the Eisenstein series: \(\langle \Theta_s, E(\cdot, w)\rangle^{\mathrm{reg}} = \tfrac12\int_0^\infty \operatorname{CT}(t)\, t^{w-2}dt\) (the \(\tfrac12\) is the half-square cusp cross-section; the pure-power line term of \(\operatorname{CT}\) contributes \(0\) after regularization). So the theorem computes the regularized pairing of the **surface kernel with the Eisenstein series**:
$$
\bigl\langle \Theta_s,\, E(\cdot, w) \bigr\rangle^{\mathrm{reg}}
= \frac{\pi}{2}\, K_s(w)\,\frac{\zeta(w)\,\zeta(w+1)}{\zeta_K(w+1)} .
$$
Unfolding the left side instead over \(\Gamma_H\) exhibits it as a transversally smoothed **period of \(E(\cdot, w)\) over the surface \(Y\)** — and the right side has exactly the shape predicted by the theory of such periods: for automorphic forms on \(\mathrm{GL}_2\) over \(K = \mathbb{Q}(i)\), the period over a \(\mathrm{GL}_2/\mathbb{Q}\)-cycle is governed by the **Asai \(L\)-function** (Flicker's period criterion: distinguished \(\iff\) Asai pole), and for the Eisenstein spectrum the Asai \(L\)-function degenerates into products of \(\zeta\) and \(L(\cdot, \chi_{-4})\) — the ratio \(\zeta(w)\zeta(w+1)/\zeta_K(w+1) = \zeta(w)/L(w+1, \chi_{-4})\) above. (A second, purely computational route to the same arithmetic: restrict the Fourier expansion of \(E\) to the plane \(z \in \mathbb{R}\) and integrate over \(x \bmod 1\) — only the purely imaginary frequencies \(\mu \in i\mathbb{Z}\) survive, and the resulting Dirichlet series of Gaussian divisor sums along \(\mathbb{Z}\) reassembles into the same ratio. We leave the exact unsmoothed period — expected \(= \Gamma\text{-factors} \times \zeta(w)\zeta(w+1)/\zeta_K(w+1)\) on the nose — as the natural sequel; see §7.)

**Pole bookkeeping (all verified numerically).** In the band, the only pole of the right side is at \(w = 1\), from \(\zeta(w)\):
$$
\operatorname{Res}_{w=1} \bigl\langle \Theta_s, E(\cdot, w)\bigr\rangle^{\mathrm{reg}}
= \frac{\pi}{2} K_s(1) \frac{\zeta(2)}{\zeta_K(2)} = \frac{\pi C_s}{2G}
\qquad\text{vs.}\qquad
\bigl\langle \Theta_s,\ \operatorname{Res}_{w=1}E \bigr\rangle = \frac{3}{2G}\cdot\frac{\pi}{3}C_s = \frac{\pi C_s}{2G} .
$$
The Maass–Selberg residue of §1 and the counting residue of §3 agree exactly — the two faces of the Catalan constant are the same number reached through the same pole. Continuing left, at \(w = 0\) the factors \(1/w\) (from \(K_s\)) and \(\zeta(w+1)\) produce a double pole: the small-\(t\) expansion of \(\operatorname{Hem}_s\) therefore carries a \(t\log(1/t)\) term after the constant — the hemisphere family "feels" the line stratum logarithmically. Further left, the poles at the zeros of \(\zeta_K(w+1)\) govern the error term in Proposition B: the equidistribution rate along the cusp is controlled by the zero-free region of \(\zeta_K\) — i.e. by the **continuous spectrum** of \(M\), since:

**Cusp forms are (essentially) invisible here.** The discrete spectrum of \(M\) pairs with \(\Theta_s\) through periods of Bianchi–Maass forms over \(Y\) — the base-change/distinction question of Flicker's criterion. For \(\mathrm{PSL}_2(\mathbb{Z}[i])\) the cuspidal spectrum is far away in any case: the first cusp form has spectral parameter \(r = 6.6221193\ldots\), \(\lambda = 1 + r^2 \approx 44.85\) (numerics of Steil and of Then, *Arithmetic quantum chaos of Maass waveforms*, math-ph/0305048, consistent with the lower bound \(\lambda_1 > 2\pi^2/3 \approx 6.58\)) — in particular **no exceptional eigenvalues**: nothing between the constant and the continuous threshold \(\lambda = 1\). This is why the teammates' elementary error term \(O(X\log X)\) is not merely consistent with the spectral picture but ahead of the crude spectral bound: the arithmetic (Möbius) structure of \(N_e\) beats generic equidistribution rates, exactly as in classical lattice-point refinements.

## 5. The thin companion: atoms count with the resonance, not the volume

Restricting the same count from the full arrangement to the **atoms** of the half-plane monoid — the strip Apollonian gasket of [half-plane-monoid.md](../half-plane-monoid.md) — replaces the lattice \(\Gamma\) by the thin Apollonian group, of infinite covolume. The counting law changes shape entirely:
$$
N_{\mathrm{atoms}}(K) \;\sim\; c\,K^{\delta}, \qquad \delta = 1.30568672\ldots
$$
(Kontorovich–Oh; \(\delta\) = Hausdorff dimension of the gasket, McMullen), with the constant now built from the Patterson–Sullivan measure of the base eigenfunction: the bottom of the spectrum of the infinite-covolume quotient is \(\lambda_0 = \delta(2 - \delta) \approx 0.9066 < 1\), an actual \(L^2\)-eigenvalue below the threshold — the "resonance" that replaces \(\operatorname{vol}(M)^{-1}\). Our exact Descartes enumeration ([orthogeodesics.md](orthogeodesics.md) §7) fits \(\delta_{\mathrm{fit}} = 1.3033\) over \(K \le 3\cdot 10^4\). Lattice versus thin, in one arrangement: exponent \(2\) with the volume-ratio constant of §3 for all circles; exponent \(\delta\) with a spectral-gap constant for the maximal (atom) circles.

## 6. Machine verification

`python3 scripts/spectral_verify.py` — 23 checks, all passing (add `--fast` to skip the slower lattice sum):

1. Humbert volume by direct integration of the fundamental domain, \(= G/3\) to 14 digits; the formula \(8\zeta_K(2)/4\pi^2 = G/3\).
2. \(\zeta_K(3)\) lattice sum vs. \(\zeta(3)L(3,\chi_{-4})\).
3. The scattering term \(\varphi(s) = (\pi/s)\zeta_K(s)/\zeta_K(s+1)\) from brute-force summation of \(E\) (relative error \(\sim 10^{-7}\) at two \((s,t)\)); \(\operatorname{Res}_{s=1}\varphi = 3/(2G) = 1/(2\operatorname{vol})\); \(\operatorname{Res}_{s=1}\zeta_K = \pi/4\).
4. Lemma 1 against the sieved \(N_e\) (to \(10^6\)); the mean-count normalization.
5. Lemma 2 and the \(K_s(w)\) evaluation at five \((s,w)\); \(K_s(1) = C_s\).
6. Proposition A/B: \(\operatorname{CT}(t) \to \pi C_s/G\) along \(t = 0.2, 0.1, 0.05, 0.02\), for the two elementary kernels \(s = 2, 1\).
7. The Mellin/period identity at \((s,w) = (2, 2.5), (2, 1.5)\) by direct \(t\)-integration.
8. The residue bookkeeping \(\tfrac12 \pi K_s(1)\zeta(2)/\zeta_K(2) = \operatorname{Res} E \cdot \int_M \Theta_s = \pi C_s/2G\), and \(1/(2G) = \tfrac1\pi \operatorname{area}(Y)/\operatorname{vol}(M)\).

## 7. Questions and next steps

1. **The exact period.** Remove the transverse smoothing: compute \(\int^{\mathrm{reg}}_Y E(P, w)\,dA\) exactly (Fourier route sketched in §4: frequencies \(\mu \in i\mathbb{Z}\), Gaussian divisor sums along \(\mathbb{Z}\)) and identify the precise \(\Gamma\)-factors in front of \(\zeta(w)\zeta(w+1)/\zeta_K(w+1)\); the \(K_s(w)\)-family computed here determines them up to the inversion of an explicit \(u\)-transform. This would be the clean Hecke-type formula for the Bianchi–Eisenstein period over the modular surface.
2. **Functional equation.** \(K_s(w)\) is even in \(w\) up to the \(1/w\) factor, and \(E(\cdot, w)\) has its functional equation through \(\varphi(w)\); combine them into the expected functional equation of the two-variable kernel pairing, and identify what \(w \mapsto -w\) does to the geometric side (a horosphere \(t \mapsto\) dual-cusp inversion).
3. **Error terms.** Make Proposition B effective: the explicit-formula version of the Mellin identity converts the zero-free region of \(\zeta_K\) into a rate in \(\operatorname{CT}(t) - \pi C_s/G\), to compare with the elementary \(O(X \log X)\); conversely, the elementary count gives (very modest) zero-density-flavored information about \(\zeta_K\) through the identity. Amusing either way.
4. **Cusp-form periods and base change.** Evaluate \(\langle \Theta_s, u_j \rangle\) for the first Bianchi–Maass forms numerically (eigenvalue data exists: \(r = 6.6221\ldots\)); Flicker's criterion predicts vanishing except for the (base-change) part of the spectrum — a directly testable instance of the period criterion on the first genuinely 3-dimensional arithmetic manifold.
5. **The relative trace formula.** §3–4 are the spectral side of the \((\mathrm{PSL}_2(\mathbb{C}), \mathrm{PGL}_2(\mathbb{R}))\)-relative trace formula whose geometric side is the orthospectrum theorem of [orthogeodesics.md](orthogeodesics.md) §4: matching the two term-by-term (orbital integrals at level \(n\) against \(H(n^2-1)\)) is outlook 3.1 of [outlook.md](../outlook.md) in its natural final form, and would give the exact second-order terms of both counting theories at once.
6. **Kudla–Millson.** The pairing \(\langle\Theta_s, \cdot\rangle\) is a crude (scalar) shadow of the Kudla–Millson theta lift attached to the cycle \(Y\); upgrading it should recover the weight-\(3/2\) slicing program of [outlook.md](../outlook.md) 3.2 — with \(\sum_n H(n^2-1) q^n\) and \(\sum_n t(n^2-1)q^n\) as the cycle-integral coefficients along the \(\alpha\)-strata.
