# The phase atlas: the Schmidt arrangement colored by its phase, the sign law, and the Duke–Imamoğlu–Tóth comparison

This document renders the level-\(n\) Schmidt circles of the ideal triangle
\(T = (0, 1, \infty)\) **colored by the phase invariant** \(u_f\) of
[moduli-invariants.md](moduli-invariants.md) §4 — the phase portrait of the Gaussian
Schmidt arrangement — and harvests the two results the one expensive dataset (every
\(u_f\), every odd level \(3 \le n \le 41\), plus \(n = 101\)) makes cheap:

1. **The sign law of \(u_f\) on ambiguous classes** (outlook.md §1.1). The
   genus-character hypothesis is **refuted** — at \(n = 31\) and \(41\) the sign is not
   multiplicative on \(\mathrm{Cl}[2]\), so *no* character of any kind produces it —
   and replaced by a **proved archimedean law**: on every divisor-type ambiguous class
   the sign is an explicit gcd inequality (Theorem 1), read off from which side of the
   \(E_6\)-zero at \(\tau = i\) the \(\mathfrak{r}_n\)-twisted CM point falls.
2. **The DIT comparison** (outlook.md §2.6). The real-quadratic cycle-integral
   invariants of discriminant \(4(n^2-1)\) — cycle integrals of \(j\) and of the
   Kronecker-limit kernel, and Rademacher symbols / linking numbers, all implemented
   from the definitions — admit **no linear relation** with \(\log|u_f|\), \(\arg u_f\)
   or the pair-sums at safe PSLQ parameters: a certified mismatch (§4), which is the
   novelty statement the phase paper's introduction needs. The cheap end of the cusp
   degeneration (outlook.md §2.4) is also negative: \(\arg u_f/\pi\) is irrational at
   every non-ambiguous class computed — no Dedekind-sum quantization away from the
   ambiguous locus.

Status labels: **proved** (mathematical proof, machine-verified), **certified**
(numerical statement with stated precision and safety margins), **observed**
(pattern in the data, no proof claimed). Scripts:
[scripts/make_phase_atlas.py](scripts/make_phase_atlas.py),
[scripts/dit_comparison.py](scripts/dit_comparison.py).

## 1. The dataset, and why it is trustworthy

For each level the circles of \(T\) come from the congruence sweep of
[hyperbolic-counting.md](hyperbolic-counting.md) §4; a circle \((q, x, m)\) is colored
through its class \(f = \overline{(q, -x, m)}\) of discriminant \(D = 1 - n^2\). Every
primitive class receives its phase unit \(u_f = \varepsilon\,\Theta_f\),
\(\varepsilon = n + \sqrt{n^2-1}\), computed by **two independent routes** at two
precisions (100 and 140 digits, \(j' = -2\pi i\,E_4^2E_6/\Delta\) via theta constants,
`mp.dps` set after imports — the guard rails of CLAUDE.md):

- **Route A** (canonical matrix): \(u = \varepsilon\, j'(m_1)X'(m_2)/\overline{j'(\bar m_2)}\)
  on the Lemma-A matrix of [class-formula-proof.md](class-formula-proof.md) — the
  closed form \(-\varepsilon\mu^{-2}h_2(\mathfrak{b}_1)/h_2(\mathfrak{b}_2)\) of
  moduli-invariants.md §5.5.
- **Route B** (derivative of the modular correspondence,
  [first-power-descent.md](first-power-descent.md) Prop. 1.3 + Thm. 3.3):
  \(u = -r_0\,h_2(\mathfrak{b}_1)/h_2(\mathfrak{r}^{-1}\mathfrak{b}_1)\),
  \(r_0 = \tfrac{n-1}2\), the twisted lattice produced by exact integer HNF
  arithmetic and both kernel values evaluated after exact
  \(\mathrm{SL}_2(\mathbb{Z})\)-reduction. Analytically this **is**
  \(\Phi_y/\Phi_x(\beta_1, \beta_2)\), evaluated through the uniformization of
  \(X_0(r_0)\) instead of through the integer coefficients of \(\Phi_{r_0}\).

**Certification** (`python3 scripts/make_phase_atlas.py --selftest`, all odd
\(3 \le n \le 41\) and \(n = 101\)):

- routes A and B agree to \(\ge 138\) digits at every primitive class of every level;
- the proved laws are re-checked on the data: \(u_{f^{-1}} = \bar u_f\) and
  \(u_{\mathfrak{r}f}u_f = 1\) to \(\ge 138\) digits, reality on ambiguous classes,
  \(|u_1|\) maximal on the principal class;
- **exactness anchors**: at every odd \(n \le 17\), every computed \(u_f\) is a root of
  the *published integer level polynomial* \(Q_n\) of moduli-invariants.md §5.9
  (residual \(\le 10^{-140}\)) — and those polynomials were produced for \(n \le 13\)
  by the exact rational arithmetic of first-power-descent.md Thm. 4.2; at \(n = 5\)
  the phase is additionally re-derived as \(\Phi_y/\Phi_x(\beta_1,\beta_2)\) from the
  **classical integer modular polynomial** \(\Phi_2\) (itself re-certified against
  \(\Phi_2(j(2\tau), j(\tau)) = 0\) at random points to 140 digits), agreeing to
  \(10^{-140}\);
- the weighted circle count re-verifies \(3H(n^2-1)\) at every atlas level.

Any mismatch anywhere in this chain would have been a bug in the atlas; there is none.

**Imprimitive classes.** \(\Theta\) is a well-defined class invariant on the
imprimitive strata too (the two-sided invariance never used primitivity), *except* on
the elliptic-core classes \(g\cdot(1,1,1)\), \(g\cdot(1,0,1)\), where it is \(0/0\)
(\(j'(\rho) = j'(i) = 0\) — the \(\alpha = 2\) phenomenon of outlook.md §1.5). The
atlas computes imprimitive phases with a two-precision stability check and draws the
elliptic-core classes neutrally. See §2.3 for what the imprimitive phases turn out to
do — one of the two genuinely new patterns the atlas exposed.

## 2. The atlas

Figures (in [figures/](figures/)): per-level panels `phase-atlas-n{3..41}.png`
(disk model and half-plane detail, hue \(= \arg u_f\), companion shade
\(= \log|u_f|/\log|u_1|\)), the contact sheets `phase-atlas-contact.png` (all levels
and \(n = 101\), hue) and `phase-atlas-contact-logabs.png` (shade), and the Euclidean
companion `phase-atlas-euclidean.png`. Overlays: the mirror geodesics
\(\operatorname{Re}z \in \{0, \tfrac12, 1\}\), \(|z| = 1\), \(|z - 1| = 1\); the sign
of \(u_f\) on ambiguous classes; thin segments joining \(\mathfrak{r}_n\)-twin
representatives. The selftest re-checks every law *on the data that is drawn*, so the
figures cannot silently drift from the mathematics.

### 2.1 Observed patterns matched to proved laws

| observation (visible in the figures) | law |
|---|---|
| every circle of an ambiguous class is pure red or pure cyan (hue only at \(\arg u \in \{\pi, 0\}\)); all *edge* circles are of this kind | \(u_{f^{-1}} = \bar u_f\) (law 1, **proved**): ambiguous \(\Rightarrow u \in \mathbb{R}\); edge circles have \(x \in \{0, 2q, 2m\}\), giving 2-torsion forms — the mirror-line geometry of [circle-composition.md](circle-composition.md) §2 |
| mirror-image circles (about \(\operatorname{Re}z = \tfrac12\)) carry complex-conjugate hues | mirror = class inversion + law 1 (**proved**) |
| \(\mathfrak{r}_n\)-twin pairs: **mirrored** hues (\(\arg \mapsto -\arg\)) and **inverted** shades (\(\log|u| \mapsto -\log|u|\)) | \(u_{\mathfrak{r}f}u_f = 1\) (law 2, **proved**). Note: the twins are *hue-mirrored*, not hue-antipodal — antipodal hues belong to the Euclidean disk pairs \(D, -D\) (\(u(-D) = -u(D)\)), clearly visible in the companion figure |
| the first non-red/cyan hues appear at \(n = 9\), and interior hue diversity grows with \(n\) | first class group with 4-torsion (\(\mathbb{Z}/4\)); complex \(u\) first at \(n = 9\) (moduli-invariants.md §4, **certified record**) |
| levels whose panel is *entirely* red/cyan: \(n = 3, 5, 7, 13, 17, 29\) (all circles real-phased; at \(n = 7\) plus the lone gray elliptic-core circle) | exactly the levels \(\le 41\) with 2-torsion class group \(\mathrm{Cl} = (\mathbb{Z}/2)^k\) — all classes ambiguous (\(n = 29\): \((\mathbb{Z}/2)^3\), eight real phases, striking in the sheet); the imprimitive phases there are real too (§2.3) |
| in the shade panels the two extreme shades sit at the cusps of \(T\) (the principal class, and its \(\mathfrak{r}_n\)-twin with the opposite extreme); the interior classes are pale (\(\lvert u\rvert \approx 1\)), so the shade organizes radially, extreme at the cusps and washing out toward the centroid | \(|u_1|\) maximal (certified at every level); growth \(\log|u_1| = \pi\sqrt{n^2-1}\,(1 + o(1))\): observed \(126.52\) vs \(\pi\sqrt{1680} = 128.78\) at \(n = 41\), \(316.85\) vs \(317.30\) at \(n = 101\) — the leading term is the \(q^{-1}\)-pole of \(F = E_4^2E_6/\Delta\) at the principal CM point of height \(\sqrt N/2\); more generally \(\log\lvert u_f\rvert\) tracks the height of the CM point (largest curvature drop across the twin) (**observed**, leading term derivable) |
| the gray (excluded) circle at the centroid of \(T\) at \(n = 7\) | the elliptic-core class \(4\cdot(1,1,1)\) — precisely the centroid circle of the Pell levels \(n^2 - 3f^2 = 1\) of [hyperbolic-counting.md](hyperbolic-counting.md) §5, where the phase degenerates to \(0/0\) (outlook.md §1.5) |

### 2.2 The sign geography (new, and the seed of §3)

In the hue panels the ambiguous classes split \(T\) into a **red majority and cyan
minority**, and the cyan classes sit in a structured way: never the principal class,
never the \(\mathfrak{r}_n\)-class, and always in twin-symmetric pairs. The sign table
of §3 quantifies this; the b = 0 part is now a theorem (Theorem 1): cyan = the
divisor classes whose \(\mathfrak{r}_n\)-twisted CM point falls *below* the unit
circle. **Observed** corollary visible in the sheets: levels \(n \equiv \pm 3
\pmod 8\) with exactly four ambiguous classes always show the coset pattern
\(-++-\) or \(----\), and the two cyan classes, when present, are the two
"mixed-divisor" classes.

### 2.3 The imprimitive strata carry their own phase theory (new)

The atlas computes the imprimitive phases as a by-product, and they are not noise
(**observed**, stability-checked to 100+ digits, `--selftest` prints the values):

- content-\(g\) classes with tiny core class number give **exact roots of unity**:
  \(u_{2(1,0,3)} = 1\) at \(n = 7\); \(u_{2(1,0,5)} = u_{2(2,2,3)} = -1\) at
  \(n = 9\); \(u_{3(1,0,8)} = u_{3(3,2,3)} = 1\) and \(u_{6(1,0,2)} = 1\) at
  \(n = 17\); \(u_{3(1,0,10)} = u_{3(2,0,5)} = -1\) at \(n = 19\) (all equalities
  to \(10^{-139}\), conjecturally exact);
- richer cores reproduce the *full law structure inside the stratum*: at \(n = 15\)
  (core disc \(-56\)) the four content-2 classes come as a twin-like pair
  \(u \approx 254720.64\), \(u' = 1/u\) **and** a conjugate pair on the unit circle;
  at \(n = 17\) (core disc \(-72\)): \(u \approx -1229104.6\) with partner \(1/u\).

So the level-\(n\) phase restricted to the content-\(g\) stratum appears to be a
phase theory of the non-maximal order \(\mathcal{O}_{D/g^2}\) — with its own twist
class and its own \(\pm1\) anchors — even though \(D/g^2\) is *not* of the form
\(1 - m^2\). None of this is covered by the proved laws (whose derivations use
invertibility of \(\mathfrak{r}_n\) and primitivity); it is a new open stratum,
recorded in the ledger.

### 2.4 The Euclidean companion

`phase-atlas-euclidean.png` colors the curvature-\(2n\) disks of the unit square by
\(\arg u\) for \(n = 7, 12\) (the lemniscatic phase of
[euclidean-moduli-invariants.md](euclidean-moduli-invariants.md) §5). The proved R/I
center criterion of §5.2 there is overlaid letter-by-letter and re-verified against
the computed phases at \(n \le 13\) (deviations \(\le 10^{-60}\)) — an end-to-end
correctness check on the rendering pipeline, run through an *independent* phase
theory. The disk pairs \(D, -D\) show the antipodal-hue law \(u(-D) = -u(D)\); the
R-pair at half-integer height shows red/cyan adjacency (\(\pm\) real values), the
I-disks the purple/chartreuse pair (\(\pm\) imaginary).

## 3. The sign of \(u_f\) on ambiguous classes: character refuted, archimedean law proved

Setup: on an ambiguous (2-torsion) class, \(u_f \in \mathbb{R}^\times\) (law 1);
\(\operatorname{sign}(u_f)\) is constant on \(\langle\mathfrak{r}_n\rangle\)-cosets
(law 2). Outlook.md §1.1 recorded the patterns for \(n \le 17\) and asked whether a
conductor-aware genus character \(\psi\) explains them via
\(\operatorname{sign}(u_f) = -\psi(f)\). The full table
(`python3 scripts/make_phase_atlas.py --signs`; classes in reduced-form order,
certified from 140-digit phases):

| \(n\) | ambiguous classes and signs | pattern |
|---|---|---|
| 3 | \((1,0,2)^-\) | \(-\) |
| 5 | \((1,0,6)^-, (2,0,3)^-\) | \(--\) |
| 7 | \((1,0,12)^-, (3,0,4)^-\) | \(--\) |
| 9 | \((1,0,20)^-, (4,0,5)^-\) | \(--\) |
| 11 | \((1,0,30)^-, (2,0,15)^+, (3,0,10)^+, (5,0,6)^-\) | \(-++-\) |
| 13 | \((1,0,42)^-, (2,0,21)^-, (3,0,14)^-, (6,0,7)^-\) | \(----\) |
| 15 | \((1,0,56)^-, (4,4,15)^+, (7,0,8)^-, (8,8,9)^+\) | \(-+-+\) |
| 17 | \((1,0,72)^-, (4,4,19)^-, (8,0,9)^-, (8,8,11)^-\) | \(----\) |
| 19 | \((1,0,90)^-, (2,0,45)^+, (5,0,18)^+, (9,0,10)^-\) | \(-++-\) |
| 21 | \((1,0,110)^-, (2,0,55)^-, (5,0,22)^-, (10,0,11)^-\) | \(----\) |
| 23 | \((1,0,132)^-, (3,0,44)^+, (4,0,33)^+, (11,0,12)^-\) | \(-++-\) |
| 25 | \((1,0,156)^-, (3,0,52)^-, (4,0,39)^-, (12,0,13)^-\) | \(----\) |
| 27 | \((1,0,182)^-, (2,0,91)^+, (7,0,26)^+, (13,0,14)^-\) | \(-++-\) |
| 29 | \((1,0,210)^-, (2,0,105)^-, (3,0,70)^+, (5,0,42)^+, (6,0,35)^+, (7,0,30)^-, (10,0,21)^+, (14,0,15)^-\) | \(--+++-+-\) |
| 31 | \((1,0,240)^-, (3,0,80)^-, (4,4,61)^+, (5,0,48)^-, (12,12,23)^-, (15,0,16)^-, (16,16,19)^+, (17,14,17)^-\) | \(--+---+-\) |
| 33 | \((1,0,272)^-, (4,4,69)^-, (16,0,17)^-, (16,16,21)^-\) | \(----\) |
| 35 | \((1,0,306)^-, (2,0,153)^+, (9,0,34)^+, (17,0,18)^-\) | \(-++-\) |
| 37 | \((1,0,342)^-, (2,0,171)^-, (9,0,38)^-, (18,0,19)^-\) | \(----\) |
| 39 | \((1,0,380)^-, (4,0,95)^+, (5,0,76)^+, (19,0,20)^-\) | \(-++-\) |
| 41 | \((1,0,420)^-, (3,0,140)^+, (4,0,105)^-, (5,0,84)^-, (7,0,60)^+, (12,0,35)^-, (15,0,28)^-, (20,0,21)^-\) | \(-+--+---\) |
| 101 | \((1,0,2550)^-, (2,0,1275)^-, (3,0,850)^+, (6,0,425)^+, (17,0,150)^+, (25,0,102)^-, (34,0,75)^+, (50,0,51)^-\) | \(--+++-+-\) |

### 3.1 The character hypothesis is refuted (certified)

**Claim (certified).** *No* map of the form \(\operatorname{sign}(u_f) = -\psi(f)\)
with \(\psi\) a homomorphism \(\mathrm{Cl}[2] \to \{\pm1\}\) — in particular no genus
character, conductor-aware or otherwise — reproduces the table.

At \(n = 41\), write \(\psi := -\operatorname{sign}(u)\) (so a character
explanation needs \(\psi\) multiplicative). Exact Gauss composition gives
\([(3,0,140)]\cdot[(4,0,105)] = [(12,0,35)]\) in \(\mathrm{Cl}(-1680)\), but
\(\psi([3]) = -1\), \(\psi([4]) = +1\), \(\psi([12]) = +1\): the product is
\(-1 \ne +1\). Likewise \([3]\cdot[5] = [15]\) there, and at \(n = 31\) the triple
\([(3,0,80)]\cdot[(4,4,61)] = [(12,12,23)]\).
The script checks every product of ambiguous classes at every level:
\(-\operatorname{sign}\) *is* multiplicative at all levels \(\le 29\), at 33–39, and
at 101 — which is why the hypothesis survived the \(n \le 17\) data — and fails
exactly at \(n = 31, 41\). The splitting-character matches that do exist at
individual levels (e.g. \(\chi_5\) at \(n = 11\), \(\chi_{-7}\) at \(n = 29\),
\(\chi_{17}\) at \(n = 101\)) are accidents of small 2-torsion. The
\(n = 11\) vs \(13\) discriminating pair of outlook.md §1.1 is subsumed: the sign is
not a character at all.

### 3.2 The archimedean law (proved for the divisor classes)

What the sign actually is: an inequality, not a character. Write
\(r_0 = \tfrac{n-1}2\), \(s_0 = \tfrac{n+1}2\), \(N = n^2-1 = 4r_0s_0\). The
ambiguous classes with \(b = 0\) are the **divisor classes**
\(f = (a, 0, c)\), \(ac = r_0s_0\), \(\gcd(a, c) = 1\): each prime power of
\(r_0s_0\) goes wholly to \(a\) or to \(c\), so
\[
a = a_r a_s, \quad c = c_r c_s, \qquad r_0 = a_r c_r,\ s_0 = a_s c_s,
\qquad a_r = \gcd(a, r_0),\ a_s = \gcd(a, s_0), \dots
\]

> **Theorem 1 (sign law on divisor classes).** For every odd \(n \ge 3\) and every
> ambiguous class \(f = (a, 0, c)\) of discriminant \(1 - n^2\),
> \[
> \operatorname{sign}(u_f) \;=\; -\operatorname{sign}\bigl(a_r c_s - a_s c_r\bigr),
> \]
> and \(a_rc_s = a_sc_r\) never happens. Equivalently: \(u_f > 0\) iff the twisted
> divisor \(a_rc_s\) (the \(r_0\)-part of \(a\) times the \(s_0\)-part of \(c\)) is
> \(< \sqrt{N}/2\).

*Proof.* By [first-power-descent.md](first-power-descent.md) Thm. 3.3 and Lemma 1.1,
\(u_f = -r_0\,h_2(\mathfrak{b}_1)/h_2(\mathfrak{r}^{-1}\mathfrak{b}_1)\) with
\(\mathfrak{b}_1 = [1, m_1]\), \(m_1 = i\sqrt N/(2a)\).

*(i) The twisted lattice, exactly.* In the basis \((1, \omega_0)\),
\(\omega_0 = i\sqrt N/2\), the module \((2a)\,\mathfrak{r}\mathfrak{b}_1\) is
generated by \((2ar_0, 0), (0, 2r_0), (0, 2a), (-2r_0s_0, 0)\). The
\(\omega_0\)-column gcd is \(2\gcd(r_0, a) = 2a_r\) (contributed by vectors with
vanishing first coordinate, so no mixing), and the first-coordinate gcd of the
remaining generators is \(\gcd(2ar_0, 2r_0s_0) = 2r_0\gcd(a, s_0) = 2r_0a_s\). Hence
\(\mathfrak{r}^{-1}\mathfrak{b}_1 = \tfrac1{r_0}\mathfrak{r}\mathfrak{b}_1 =
\lambda\,[1, \tau_2]\) with \(\lambda = \tfrac{a_s}{a} \in \mathbb{Q}_{>0}\) and
\[
\tau_2 = i y_2, \qquad
y_2 = \frac{2a_r}{2r_0a_s}\cdot\frac{\sqrt N}{2}
    = \frac{a_r\sqrt{r_0s_0}}{r_0\,a_s}
    = \sqrt{\frac{a_rc_s}{a_sc_r}}.
\]
(The script asserts this closed form against the machine HNF at every class.)

*(ii) Positivity of the kernel above \(i\).* \(F = E_4^2E_6/\Delta\) is real on the
imaginary axis; \(E_4(iy) > 0\) and \(\Delta(iy) > 0\) termwise, and
\(E_6(iy) = 1 - 504\sum\sigma_5(m)e^{-2\pi my}\) is strictly increasing in \(y\) with
\(E_6(i) = 0\) (weight-6 \(S\)-invariance), so \(F(iy) > 0\) for \(y > 1\).

*(iii) Assembling.* \(m_1 = iy_1\) with \(y_1 = \sqrt N/(2a) > 1\) (reduced,
\(a < c\) strictly since \(\gcd(a,c) = 1\) excludes \(a = c\)), so
\(h_2(\mathfrak{b}_1) = C\,F(iy_1)\) with \(F(iy_1) > 0\). By homogeneity
\(h_2(\mathfrak{r}^{-1}\mathfrak{b}_1) = \lambda^{-2}C\,F(iy_2)\), and by the
weight-2 functional equation \(F(iy_2) = -y_2^{-2}F(i/y_2)\) when \(y_2 < 1\). Thus
\(\operatorname{sign}(u_f) = -\operatorname{sign} F(iy_2) = -1\) if \(y_2 > 1\) and
\(+1\) if \(y_2 < 1\), i.e. \(-\operatorname{sign}(a_rc_s - a_sc_r)\) by (i).
Finally \(a_rc_s = a_sc_r\) would force (pairwise coprimality of the four factors)
\(a_r = c_r = 1\) and \(a_s = c_s = \sqrt{s_0}\) with \(r_0 = 1\), i.e. \(n = 3\),
\(s_0 = 2 = a_s^2\) — impossible. \(\blacksquare\)

Machine verification: the predicted sign matches the computed one at **every**
divisor-type ambiguous class of every odd \(n \le 41\) and of \(n = 101\) — 81
classes — and the closed-form HNF is asserted exactly at each
(`--signs`). The theorem explains the refutation of §3.1 structurally: an inequality
between divisors is not multiplicative. It also recovers the old observations: the
principal class has \((a_r, a_s) = (1,1)\), \(y_2 = \sqrt{s_0/r_0} > 1\), sign
\(-\) always; the \(\mathfrak{r}_n\)-class \((r_0, 0, s_0)\) has
\(y_2 = \sqrt{r_0s_0} > 1\), sign \(-\); and a \(+\) requires genuinely mixed
divisors — impossible
when \(r_0\) or \(s_0\) is \(1\) or when every ambiguous \(a\) is a pure
\(r_0\)- or \(s_0\)-part, which is exactly what happens at
\(n = 3, 5, 7, 9, 13, 17, 21, 25, 33, 37\).

**The geometric reading.** \(\tau_2\) is the CM point of the
\(\mathfrak{r}_n\)-twisted lattice, and the sign records **on which side of the
\(E_6\)-zero locus (the \(\mathrm{SL}_2(\mathbb{Z})\)-orbit of \(i\)) the twisted
point falls** — equivalently, the orientation sign
\(\operatorname{sign}((c\tau_2 + d)^2)\) of its reduction cocycle. This is precisely
"the sign data lives in the geometry of \(\Phi_{r_0}\) at real points"
(first-power-descent.md §5), now made exact; and it connects to the
denominator law: the same two elliptic fixed points whose Deuring collisions
produce the Gross–Zagier denominators (moduli-invariants.md §5.7) govern the
archimedean sign.

### 3.3 What remains open: the 2-adic ambiguous classes

At levels \(n \equiv \pm1 \pmod 8\) (and only there) some ambiguous classes have
\(b \ne 0\): types \((a, a, c)\) (CM point on \(\operatorname{Re} = -\tfrac12\)) and
\((a, b, a)\) (on \(|\tau| = 1\)) — the classes over the ramified prime 2. There the
per-factor decomposition still holds (the script prints
\(E_6\)-region and cocycle signs per class; on the \(\operatorname{Re} = \tfrac12\)
line \(\Delta < 0\) and the region sign flips), but a uniform closed form in the
shape of Theorem 1 is not yet extracted: at \(n = 15\): \((4,4,15)^+, (8,8,9)^+\);
at \(n = 17\): \((4,4,19)^-, (8,8,11)^-\); at \(n = 31\): \((4,4,61)^+,
(12,12,23)^-, (16,16,19)^+, (17,14,17)^-\); at \(n = 33\): \((4,4,69)^-,
(16,16,21)^-\). **Open** (ledger): the 2-adic sign law — same mechanism, needs the
HNF closed form on the ramified stratum.

## 4. The Duke–Imamoğlu–Tóth comparison

**The bridge.** The level-\(n\) data lives on the real-quadratic side through the
Cartan geodesic: trace \(-2n\), length \(2\log\varepsilon_n\), discriminant
\(d = 4(n^2-1) = -4D\). The alignment is exact: the minimal solution of
\(t^2 - du^2 = 4\) is \((t, u) = (2n, 1)\), so **the fundamental automorph of every
form of discriminant \(d\) is \(\varepsilon_n\) itself** — the unit that normalizes
the phase. Note the discriminant pairing this puts on the table: the phase couples
\(D = 1 - n^2 < 0\) with \(d = -4D > 0\), exactly the index shape \((D, d)\) of the
two-sign Fourier coefficients in DIT's Katok–Sarnak framework.

**What was computed** (`python3 scripts/dit_comparison.py`, from the definitions —
no tables trusted): for \(n = 3, \dots, 17\) odd, all
\(\mathrm{SL}_2(\mathbb{Z})\)-classes of forms of discriminant \(d\) (Gauss-reduced
cycles; imprimitive classes included and marked), and per class \(A\):

- the **cycle integral** \(C_A(f) = \int_{\Gamma_Q\backslash S_Q} f(z)\,
  \tfrac{\sqrt d\,dz}{Q(z,1)} = \int_1^{\varepsilon_A^2} f(\sigma(iy))\tfrac{dy}y\)
  for \(f = j - 744\) and for the Kronecker-limit kernel
  \(f = \log(\operatorname{Im}(z)^6|\Delta(z)|)\) (write \(E_A\)), where
  \(\varepsilon_A\) is the **primitive** stabilizer eigenvalue of the class — for the
  imprimitive classes with smaller core discriminant the level geodesic *wraps* the
  primitive one (at \(n = 7\) the content-4 class rides the core-disc-12 geodesic,
  \(\varepsilon_{12}^2 = \varepsilon_7\); at \(n = 9\) the content-8 class rides the
  **golden geodesic** of disc 5, \(\varepsilon_9 = \bigl(\tfrac{3+\sqrt5}2\bigr)^3\)),
  and the integral runs over one primitive period, as in DIT;
- the **Rademacher symbol** \(\Psi(\gamma_A)\) of the automorph (= the linking number
  of the modular knot with the trefoil), computed **twice**: exact Dedekind sums
  (\(\Phi = \tfrac{a+d}c - 12\operatorname{sign}(c)\,s(d, |c|)\),
  \(\Psi = \Phi - 3\operatorname{sign}(c(a+d))\)) and the minus-continued-fraction
  cycle (\(\Psi = \sum(b_i - 3)\)); the two agree at every class of every level.

**Validation of the DIT side**: \(C_A(1) = 2\log\varepsilon_A\) to \(10^{-70}\);
every value re-derived under random \(\mathrm{SL}_2(\mathbb{Z})\)-transforms of the
form and under an independent second parametrization (\(d\theta/\sin\theta\) along
the semicircle), agreeing to \(10^{-140}\)–\(10^{-200}\) (the working precision is
elevated internally by \(\pi\sqrt d/\ln 10\) digits, the growth of the \(j\)-kernel
at the top of the tallest geodesic, so that 150 digits stay certified);
Dedekind-sum reciprocity checked against the brute-force definition; \(\Psi\)
conjugation-invariance checked on random conjugates; hand-proved anchors
\(\Psi = 0\) (golden geodesic, \(d = 5\)) and \(\Psi([[4,-1],[1,0]]) = 1\)
(\(d = 12\)) reproduced by both routes. Sample \(\Psi\)-vectors:
\((11, 1, 4, -1, 1, -1, -4, -11)\) at \(d = 192\),
\((31, 3, 7, 9, 2, 1, 1, -1, -1, 14, -2, -9, 0, -7, 4, -3, -4, -14, -31)\) at
\(d = 1152\) — always antisymmetric under class inversion.

Sample values (full tables in the script log): at \(n = 3\) (\(d = 32\), 3 classes):
\(\mathrm{Tr}_d(j - 744) = -318.54299833\ldots\), \(\sum_A E_A = -84.12550259\ldots\),
\(\Psi = (3, 0, -3)\); at \(n = 5\) (\(d = 96\), 6 classes):
\(\mathrm{Tr}_d(j-744) = -809.26549902\ldots\), \(\sum_A E_A = -244.95086317\ldots\),
\(\Psi = (7, 1, -1, 2, -2, -7)\).

### 4.1 The comparison battery and its verdict

PSLQ with the certification guard rails (a candidate is accepted only when
\((\text{terms})\times(\text{coefficient digits})\) is far below the working
precision **and** the residual sits at the noise floor; every run logs its margin).
Working precision 150 digits; phase data at 170. Bases per level, with the
conjugate-pair degeneracy \(C_{A^{-1}} = \overline{C_A}\), \(E_{A^{-1}} = E_A\)
removed:

| test | basis | coeff. bound | outcome |
|---|---|---|---|
| \(\log\lvert u_f \rvert\) (one per twin pair) | \([1, \log\varepsilon, \mathrm{Tr}_d(j{-}744), \sum E_A, \pi^2]\) | \(10^8\) | **no relation**, every class, every level |
| \(\log\lvert u_f \rvert\) | \([1, \log\varepsilon]\) | \(10^{10}\) | **no relation** (so \(u_f\) is not \(\pm\,\mathrm{rational}\times\varepsilon^k\) in disguise) |
| \(\log\lvert u_f \rvert\) | per-pair \(\operatorname{Re}C_A(j{-}744)\) and \(E_A\) | \(10^6\) | **no relation** (run wherever the deduplicated basis stays \(\le 12\) terms, i.e. \(n \le 15\)) |
| \(\arg u_f\) (non-ambiguous) | \([\pi, \log\varepsilon]\) | \(10^{10}\) | **no relation** |
| \(\arg u_f\) | per-pair \(\operatorname{Im}C_A(j{-}744)\), \(\pi\) | \(10^6\) | **no relation** |
| pair-sums \(S_x\) (real) | aggregate basis | \(10^6\) | **no relation** (the sole hit is \(n = 3\): \(S + 2 = 0\), i.e. the proved anchor \(u = -1\)) |

> **Verdict (certified).** At \(3 \le n \le 17\), neither \(\log|u_f|\) nor
> \(\arg u_f\) nor the pair-sums \(S_x\) is a \(\mathbb{Q}\)-linear combination of
> the discriminant-\(4(n^2-1)\) cycle integrals of \(j\), the Kronecker-limit
> integrals \(E_A\), \(\log\varepsilon_n\), \(\pi\), \(\pi^2\) and \(1\), with
> coefficients of height up to the stated bounds, at 150-digit precision. **The
> Schmidt phase is not a linear shadow of the known cycle-integral invariants of its
> own real-quadratic discriminant.**

A by-product worth recording (**observed**, spare \(\ge 140\) digits): the
\(E_A\)-basis itself carries integer relations *across strata* — e.g. at \(d = 192\):
\(E_{(-12,12,1)} + 6E_{(-8,8,4)} - 5E_{(-6,12,2)} + E_{(-4,12,3)} = 0\), mixing the
primitive stratum with the wrapped core-disc-12 and core-disc-48 classes. This is
genuine mathematics, not noise: the \(E_A\) are logarithms of algebraic
(Kronecker-limit) invariants, which can be multiplicatively dependent. The battery
quotients such internal relations out before testing the phase target, and logs each
one.

This is the expected outcome, and the structural reason deserves record: \(u_f\) is
**algebraic** (proved, first-power-descent.md), while cycle integrals of \(j\) are
expected to be transcendental and are not known (or expected) to be logarithms of
algebraic numbers; a nontrivial linear identity would have been a transcendence
miracle. The value of the certified mismatch is the licence it gives the phase
paper: \(u_f\) is genuinely *new* real-quadratic-adjacent data at
\(d = 4(n^2-1)\), not repackaged DIT material. What the two theories *do* share —
exactly \(\varepsilon_n\), the coupled pair \((D, -4D)\), and (Theorem 1) the
\(E_6\)-geometry at the elliptic point — marks the structured place to look next:
DIT-type **mixed** objects indexed by a pair of discriminants of opposite sign
(Katok–Sarnak two-sign coefficients, biharmonic/local polynomials), not the
one-discriminant invariants tested here.

### 4.2 The Rademacher/linking side

\(\Psi\) on the trace-\(2n\) classes comes out antisymmetric under class inversion
(\(\Psi(A^{-1}) = -\Psi(A)\); e.g. \((3, 0, -3)\) at \(n = 3\),
\((7, 1, -1, 2, -2, -7)\) at \(n = 5\)), so \(\sum_A \Psi_A = 0\) identically — while
\(\operatorname{sign}(u_f)\) is defined on the *imaginary*-quadratic 2-torsion, is
not antisymmetric, and (§3.1) is not even a character. A per-class matching
\(\operatorname{sign}(u_f) \leftrightarrow (-1)^{\Psi}\)-type is already
dimensionally impossible in general (\(h(1-n^2) \ne h(4(n^2-1))\) as multisets of
classes); no aggregate correlation survives the battery either. **Observed**: no
Rademacher shadow of the phase.

### 4.3 The cusp-degeneration probe (outlook.md §2.4, the cheap end)

If the \(\arg u\)-cocycle degenerated to classical Dedekind-sum data at the compact
levels, \(\arg u_f/\pi\) would show small-denominator rationals. Probe
(`--cusp`, 110 digits, all odd \(n \le 41\), 86 non-ambiguous classes): **no**
\(\arg u_f/\pi\) is rational with denominator \(\le 10^4\); the closest approach in
the CF quality measure \(\operatorname{err}\cdot q^2\) is \(\approx 5\cdot10^{-3}\)
(\(n = 39\), \(\arg u/\pi \approx \pm7/86\), err \(7\cdot10^{-7}\) — far above the
\(10^{-110}\) noise floor, i.e. genuinely irrational at this precision). The phase
angle quantizes on the ambiguous locus (\(\{0, \pi\}\), proved) and **nowhere else**;
whatever the \(n \to 1\) Ford degeneration of the \(\mu_{12}\)-bookkeeping is, it is
not visible as Dedekind rationality at compact levels. The probe does not touch the
genuine degeneration regime (a family \(n \to 1\) does not exist through odd integer
levels); the honest formulation of §2.4 remains the \(\eta\)-multiplier limit, still
open.

## 5. Machine verification

| check | command | result |
|---|---|---|
| routes A = B (138+ digits), laws 1–2, ambiguity/reality, principal max, published \(Q_n\) roots (\(n \le 17\)), \(\Phi_2\) anchor, Euclidean R/I, \(3H(n^2-1)\) | `python3 scripts/make_phase_atlas.py --selftest` | all pass, every odd \(3 \le n \le 41\) and \(n = 101\) |
| figures | `python3 scripts/make_phase_atlas.py --figures` | 20 per-level panels + 2 contact sheets + Euclidean companion into `figures/` |
| sign table, character refutation, Theorem-1 rule, HNF closed form, \(E_6\)/cocycle decomposition | `python3 scripts/make_phase_atlas.py --signs` | table of §3; multiplicativity fails exactly at \(n = 31, 41\); Theorem-1 rule matches all 81 divisor classes |
| DIT selftest (Dedekind reciprocity vs brute force, \(\Psi\) anchors \(d = 5, 12\), conjugation invariance, \(C(1) = 2\log\varepsilon\), two parametrizations, \(\mathrm{SL}_2\)-invariance) | `python3 scripts/dit_comparison.py --selftest` | all pass |
| the comparison battery of §4.1 | `python3 scripts/dit_comparison.py` | margins logged per test; no non-trivial relation |
| cusp probe | `python3 scripts/dit_comparison.py --cusp` | no rational \(\arg u/\pi\), 86 classes |

Precisions: phases at 100 and 140 digits (two-precision drift \(\le 10^{-95}\));
DIT integrals certified at 150 digits (internally elevated by
\(\pi\sqrt d/\ln 10\) digits against the kernel's growth), validations to
\(10^{-140}\)–\(10^{-200}\); PSLQ at 150 digits with the (terms) × (digits) safety
rule of CLAUDE.md. Across the full battery (\(n = 3, \dots, 17\)): 52 certified
non-relations, zero unsafe or rejected candidates, and the only target-involving
hit is the proved \(n = 3\) anchor \(S + 2 = 0\).

## 6. Outlook

1. **The 2-adic sign law** (§3.3): redo the Theorem-1 HNF computation on the
   \((a,a,c)\) and \((a,b,a)\) ambiguous families; expect a closed form in the
   \(E_6\)-region/orientation data of the \(\operatorname{Re} = \tfrac12\) and
   \(|\tau| = 1\) walls. Together with Theorem 1 this would close outlook.md §1.1
   completely.
2. **The imprimitive phase strata** (§2.3): identify the stratum theory — twist
   class, \(\pm1\) anchors, unit-circle classes — presumably the phase theory of the
   suborder \(\mathcal{O}_{D/g^2}\) with its own \(\varepsilon\); prove the observed
   exact values \(u = \pm1\). This is fresh structure the atlas exposed.
3. **Mixed two-discriminant DIT objects**: the certified mismatch (§4.1) rules out
   the one-discriminant linear shadow; the coupled pair \((1-n^2,\ 4(n^2-1))\) with
   shared unit \(\varepsilon_n\) is exactly the index set of Katok–Sarnak-type mixed
   coefficients — the right comparison object for a second pass.
4. **\(E_6\)-geometry unification**: Theorem 1 (archimedean sign from the
   \(E_6\)-zero), the GZ-denominator law (finite places from Deuring collisions with
   \(j = 0, 1728\), moduli-invariants.md §5.7), and the \(\tau_2 \to i\) collision of
   the principal class as \(n \to \infty\) all point at one height-type statement on
   the coupled fiber product — a concrete entry into outlook.md §3.5.
5. **Atlas continuation**: even levels (\(i\mathcal{S}\), outlook.md §1.4) and the
   \(\alpha = 2\)/elliptic-core regularization (outlook.md §1.5) would fill the gray
   circles of the figures.
