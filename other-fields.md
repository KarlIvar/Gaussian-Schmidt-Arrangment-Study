# The Schmidt line over other imaginary quadratic fields

This document carries out item 3 of [PROGRAM.md](PROGRAM.md): the line of
Papers I and II — classification, level, involution and class formula, the
twisted \(\Delta\)-units, the Euclidean \(\Delta\)-data with their mass law and
per-class valuations, the Kronecker limit formulas, the Euler system over the
ring class tower, and the full Robert index — is transported from
\(\mathbb{Q}(i)\) to Stange's arrangements
$$
\mathcal{S}_K \;=\; \mathrm{PSL}_2(\mathcal{O}_K)\cdot\hat{\mathbb{R}},
\qquad K = \mathbb{Q}(\sqrt{d_K}),\quad \mathcal{O}_K = \mathbb{Z}[\omega],
$$
for the class-number-one fields \(K = \mathbb{Q}(\sqrt{-2})\), \(\mathbb{Q}(\sqrt{-3})\),
\(\mathbb{Q}(\sqrt{-7})\) and \(\mathbb{Q}(\sqrt{-11})\) (the Gaussian field is re-run as the
anchor), and, at the first two levels, to the class-number-two field
\(\mathbb{Q}(\sqrt{-5})\). Everything finite is verified by
[scripts/other_fields.py](scripts/other_fields.py) (`--selftest`; exact
\(\mathcal{O}_K\)-lattice arithmetic, certified \(\Delta\)-values, PARI/GP for the
class fields) under the guard rails of [CLAUDE.md](CLAUDE.md); §10 states what is
proved for every \(K\), what is proved for the computed \(K\), and what is
certified only.

**The results in one paragraph.** The circles of \(\mathcal{S}_K\) are the Hermitian
matrices \(\binom{\sqrt{|d_K|}q\ \ i\beta}{-i\bar\beta\ \ \sqrt{|d_K|}m}\) with
\(\beta \in \mathcal{O}_K\), \(\beta \equiv 1 \pmod{\sqrt{d_K}\mathcal{O}_K}\) and
\(N(\beta) = 1 + |d_K|qm\) (necessity for every \(K\), sufficiency by descent for the
five Euclidean fields). For \(K \ne \mathbb{Q}(i)\) every circle carries a **unique
orientation**, the level \(\alpha = -\operatorname{sgn}(q)\operatorname{Re}\beta\) lies in
\(\tfrac12\mathbb{Z}\) with \(2\alpha \equiv -2\operatorname{sgn}(q) \pmod{|d_K|}\), and the levels
alternate in orientation. Level-\(\alpha\) circles are the forms of discriminant
\(D_K(\alpha) = 4(\alpha^2-1)/d_K\), counted by \(3H(|D_K(\alpha)|)\). The involution
\(\sigma(X) = \bar X^{-1}\) acts by \(\hat\sigma[f] = [\mathfrak{r}_\alpha][f]^{-s}\), \(s\) the
orientation sign, with the ambiguous twist ideal \(\mathfrak{r}_\alpha\) of norm
\(r_0 = (2\alpha-2)/\gcd(2\alpha-2,|d_K|)\) (the "\(\alpha-1\) side" of
\(N(\theta_D) = r_0s_0\)); the proof is Paper I's with the Gram form
\(g_s(u) = \alpha N(u) + s\operatorname{Re}(u^2)\). The twisted ratios
\(R_f = r_0^6\Delta(\mathfrak{b}_f)/\Delta(\mathfrak{r}_\alpha^{-1}\mathfrak{b}_f)\) are units
(Paper II's proof verbatim), with \(\sum\chi\log|R_f| = -24L'(0,\chi)\) on odd
characters. On the Euclidean side, disks of curvature \(\sqrt{|d_K|}n\) modulo
translation are the primitive index-\(n\) sublattices of \(\mathcal{O}_K\),
\(N_e(n) = \tfrac{w_K}{2}h(\mathcal{O}_n)\); the \(\Delta\)-data
\(G_\mathfrak{c} = n^{12}\Delta(\Lambda_\mathfrak{c})/\Delta(\mathcal{O}_K)\) have integer level
polynomials, the mass law
\(v_p(M(n)) = \tfrac{24}{w_Ke_p}\tfrac{p^k-1}{p-1}N_e(n/p^k)\) at the non-split
primes, and class-independent valuations
\(w_p(k) = 12(p^k-1)/(e_p(p-1)N_e(p^k))\); the Kronecker limit formula
\(\sum\chi\log|G_\mathfrak{c}| = -12L'(0,\chi)\) holds, the genus character at a prime
level \(p\) has real field \(\mathbb{Q}(\sqrt p)\) (\(p \equiv 1 \bmod 4\)) or
\(\mathbb{Q}(\sqrt{p|d_K|})\) (\(p \equiv 3 \bmod 4\)); the norm relations of the
Euler system hold with the Frobenius factors of \(K\), the constant \(\pi^{12}\) at
the ramified prime and a global sign \((-1)^{[\ell=2]}\); and the full Robert index is
$$
\bigl[\mathcal{O}_{H_n}^\times:\mu(H_n)\mathcal{V}_n\bigr]
= \frac{w_K\cdot24^{h-1}}{w_{H_n}}\,h_{H_n}\prod_{\chi\ne1}C_\chi(0),
$$
with \(\mu(H_n)\) from a conductor lemma — certified with PARI (`bnfcertify`) at
every computed level. At \(\mathbb{Q}(\sqrt{-5})\) the single cusp sees only
\(\ker(\mathrm{Pic}(\mathcal{O}_n)\to\mathrm{Pic}(\mathcal{O}_K))\); the two-cusp data
\(n^{12}\Delta(\Lambda)/\Delta(\mathcal{O}_K\Lambda)\) satisfy the limit formula with a
correction on the characters pulled back from \(\mathrm{Pic}(\mathcal{O}_K)\), and the
index acquires the factor \(w_K/h_K\) together with a base-character correction
(experimental).

**Notation.** \(d_K\) the discriminant, \(\sqrt{d_K} := i\sqrt{|d_K|}\),
\(\omega = \sqrt{d_K}/2\) (\(d_K \equiv 0 \bmod 4\)) or \((1+\sqrt{d_K})/2\)
(\(d_K \equiv 1 \bmod 4\)); \(w_K = |\mathcal{O}_K^\times|\), \(h_K\); \(\chi_K\) the
Kronecker symbol of \(d_K\); for \(x \in \mathcal{O}_K\), \(\operatorname{im}_K(x) := 2\operatorname{Im}(x)/\sqrt{|d_K|} \in \mathbb{Z}\)
is the \(\omega\)-coordinate. Orders \(\mathcal{O}_n = \mathbb{Z}+n\mathcal{O}_K\) of discriminant \(n^2d_K\),
\(N_e(n) = n\prod_{p\mid n}(1-\chi_K(p)/p)\), \(h = h(\mathcal{O}_n) = |\mathrm{Pic}(\mathcal{O}_n)|\), \(H_n\) the ring
class field. Circles are Hermitian matrices \(M = \binom{A\ B}{\bar B\ C}\) of determinant
\(-1\), \(M_{g\mathcal{C}} = (g^{-1})^\dagger Mg^{-1}\), \(\hat{\mathbb{R}} \leftrightarrow M_0 = \binom{0\ i}{-i\ 0}\),
the interior \(\{M[z] < 0\}\) for \(A > 0\) (Paper I §2.1). Forms \((a,b,c)\) of discriminant
\(D\) correspond to the ideals \(\mathfrak{a}_f = \mathbb{Z}a + \mathbb{Z}\tfrac{b+\sqrt D}{2}\) (Paper I's
convention: a positively oriented basis, \(\operatorname{Im}(\bar\alpha\beta) > 0\), has norm form
\(f\)); \(\theta_D\) denotes \(\sqrt D\) (odd \(d_K\)) resp. \(\sqrt D/2\) (even \(d_K\)). Status
conventions as in the papers: **proved**, **certified** (finite statement verified by
the script), **experimental**.

## 1. Classification, level, orientation

For \(X = \binom{a\ b}{c\ d} \in \mathrm{SL}_2(\mathcal{O}_K)\), Paper I's two-line computation gives
$$
M_{X(\hat{\mathbb{R}})} = \begin{pmatrix} 2\operatorname{Im}(c\bar d) & i(a\bar d - b\bar c)\\ \overline{i(a\bar d-b\bar c)} & 2\operatorname{Im}(a\bar b)\end{pmatrix}.
$$
Since \(\operatorname{Im}(\mathcal{O}_K) = \tfrac{\sqrt{|d_K|}}{2}\mathbb{Z}\), the curvature and co-curvature lie in
\(\sqrt{|d_K|}\mathbb{Z}\) (Stange, IMRN 2018) and we write
$$
M = \begin{pmatrix}\sqrt{|d_K|}\,q & i\beta\\ -i\bar\beta & \sqrt{|d_K|}\,m\end{pmatrix},
\qquad q = \operatorname{im}_K(c\bar d),\quad \beta = a\bar d - b\bar c,\quad m = \operatorname{im}_K(a\bar b),
$$
so that \(\det M = -1\) reads \(N(\beta) = 1 + |d_K|qm\). The circle has radius \(1/(\sqrt{|d_K|}|q|)\)
and centre \(-i\beta/(\sqrt{|d_K|}q)\); translation by \(\lambda \in \mathcal{O}_K\) acts by
\(\beta \mapsto \beta + \sqrt{d_K}\,q\lambda\), \(m \mapsto m + \operatorname{im}_K(\bar\lambda\beta) + qN(\lambda)\), and
\(z \mapsto -1/z\) by \((q,\beta,m)\mapsto(m,\bar\beta,q)\).

> **Theorem 1 (classification, level, orientation).**
> 1. *(Necessity, every \(K\).)* Every circle of \(\mathcal{S}_K\) is of the form above with
>    \(\beta \equiv 1 \pmod{\sqrt{d_K}\mathcal{O}_K}\). The lines of \(\mathcal{S}_K\) are
>    \(\operatorname{Im}(u^2z) \in \tfrac{\sqrt{|d_K|}}{2}\mathbb{Z}\), \(u \in \mathcal{O}_K^\times\): one direction
>    for \(w_K \in \{2,4\}\), three directions for \(\mathbb{Q}(\sqrt{-3})\).
> 2. *(Sufficiency, \(d_K \in \{-3,-4,-7,-8,-11\}\).)* Conversely every such Hermitian matrix
>    (with \(q \ne 0\), or a line as in 1) lies in the orbit of \(M_0\): descent by translation
>    and inversion terminates.
> 3. *(Orientation.)* \(-M_0\) lies in the orbit of \(M_0\) iff \(K = \mathbb{Q}(i)\). Hence for
>    \(K \ne \mathbb{Q}(i)\) every circle of \(\mathcal{S}_K\) carries a unique orientation, and the sign of
>    \(q\) is an invariant of the circle.
> 4. *(Level.)* For a circle in \(\mathbb{H}\), \(\alpha := \operatorname{Im}(z_0)/r = -\operatorname{sgn}(q)\operatorname{Re}\beta\)
>    (\(= \coth\) of the hyperbolic radius, the inversive product with \(\hat{\mathbb{R}}\), and
>    \(-\tfrac12\operatorname{tr}(X\bar X^{-1})\) up to orientation, exactly as in Paper I, Prop. 2.5).
>    Since \(\operatorname{Re}\beta \in 1 + \tfrac{|d_K|}{2}\mathbb{Z}\), the level is a half-integer with
>    $$2\alpha \;\equiv\; -2\operatorname{sgn}(q) \pmod{|d_K|},$$
>    so for \(K \ne \mathbb{Q}(i)\) the level determines the orientation sign \(s = s(\alpha)\), and
>    for \(\mathbb{Q}(i)\) both orientations occur at every odd level. The levels are all
>    \(\alpha \ge 2\) with \(2\alpha \equiv \pm2 \bmod |d_K|\) (integers for \(d_K \equiv 0 \bmod 4\)).

*Proof.* (1) \(\bar x \equiv x \pmod{\sqrt{d_K}\mathcal{O}_K}\) for all \(x \in \mathcal{O}_K\) (as \(\omega - \bar\omega = \sqrt{d_K}\)),
so \(\beta = a\bar d - b\bar c \equiv ad - bc = 1\). Lines: \(q = 0\) means \(c/d \in \mathbb{Q}\), so
\(c = pu\), \(d = q'u\) with \(p, q'\) coprime integers and \(u\) a unit (\(h_K = 1\) is not needed:
\(\gcd(c,d) = 1\) as ideals and \(c\bar d\) real); \(ad - bc = 1\) gives \(u^{-1} = aq' - bp\) and
\(B = i\bar u u^{-1} = i\bar u^2\), so the line is \(2\operatorname{Im}(u^2z) + C = 0\) with \(C \in \sqrt{|d_K|}\mathbb{Z}\);
conversely \(\binom{\bar u\ b}{0\ u}\) realizes every such line. (2) Translate so that
\(N(\beta)\) is minimal: the lattice \(\sqrt{d_K}q\mathcal{O}_K\) has covering radius
\(|q|\sqrt{|d_K|}\rho_K\), \(\rho_K^2 = \tfrac{1+m}{4}\) for \(d_K = -4m\) and \(\tfrac{(1+m)^2}{16m}\) for
\(d_K = -m\), i.e. \(\rho_K^2 = \tfrac12, \tfrac34, \tfrac13, \tfrac47, \tfrac9{11} < 1\) for the five fields, so
\(N(\beta) \le q^2|d_K|\rho_K^2 < q^2|d_K| + 1\) and \(|m| = (N(\beta)-1)/(|d_K||q|) < |q|\). Inverting
strictly decreases \(|q|\); the walk ends at \(q = 0\), where \(\beta\) is a unit \(\equiv 1\), hence
\(\beta = \bar u^2\) for a unit \(u\) (the units \(\equiv 1 \bmod \sqrt{d_K}\) are \(\pm1\) for \(\mathbb{Q}(i)\),
\(\{1,\zeta_3,\zeta_3^2\}\) for \(\mathbb{Q}(\sqrt{-3})\), \(\{1\}\) otherwise, and these are exactly the
\(\bar u^2\)), a line of \(\mathcal{S}_K\) by (1). Reversing the walk exhibits \(M\) in the orbit. (3) The
orientation-reversing stabilizers of \(\hat{\mathbb{R}}\) in \(\mathrm{SL}_2(\mathbb{C})\) are the \(ih\),
\(h \in \mathrm{GL}_2(\mathbb{R})\), \(\det h = -1\); \(ih \in M_2(\mathcal{O}_K)\) forces the entries of \(h\) into
\(\mathcal{O}_K \cap i\mathbb{R}\cdot(-i) = \sqrt{|d_K|/4}\,\mathbb{Z}\) resp. \(\sqrt{|d_K|}\mathbb{Z}\), and \(\det h = -1\) then
forces \(|d_K|/4 = 1\). (4) The formula \(\alpha = \operatorname{Im}(z_0)/r\) is Paper I's; with the centre
\(-i\beta/(\sqrt{|d_K|}q)\) it reads \(-\operatorname{Re}(\beta)\operatorname{sgn}(q)\). From \(\beta = 1 + \sqrt{d_K}\gamma\),
\(\operatorname{Re}\beta = 1 - \sqrt{|d_K|}\operatorname{Im}\gamma \in 1 + \tfrac{|d_K|}{2}\mathbb{Z}\). \(\square\)

**Certified** (script, phase A): for each of the five fields the breadth-first orbit of
\(M_0\) (translations by \(\pm1, \pm\omega\), inversion; about 900 circles) has exactly the
residues \(\beta \bmod q\sqrt{d_K}\mathcal{O}_K\) allowed by the congruences for every \(|q| \le 5\)
(their number is \(N_e(|q|)\)), and its level set is the predicted one. The level sets:

| \(K\) | levels \(2\alpha\) (orientation sign) |
|---|---|
| \(\mathbb{Q}(i)\) | \(6, 10, 14, \dots\) (\(\pm\), both) — the odd levels of Paper I |
| \(\mathbb{Q}(\sqrt{-2})\) | \(6\,(+), 10\,(-), 14\,(+), 18\,(-), \dots\): all odd \(\alpha \ge 3\), \(s = +\) iff \(\alpha \equiv 3 \bmod 4\) |
| \(\mathbb{Q}(\sqrt{-3})\) | \(4\,(+), 5\,(-), 7\,(+), 8\,(-), 10\,(+), 11\,(-), \dots\): \(2\alpha \not\equiv 0 \bmod 3\), \(s = +\) iff \(2\alpha \equiv 1 \bmod 3\) |
| \(\mathbb{Q}(\sqrt{-7})\) | \(5\,(+), 9\,(-), 12\,(+), 16\,(-), 19\,(+), 23\,(-), \dots\): \(2\alpha \equiv \mp2 \bmod 7\) |
| \(\mathbb{Q}(\sqrt{-11})\) | \(9\,(+), 13\,(-), 20\,(+), 24\,(-), 31\,(+), 35\,(-), \dots\): \(2\alpha \equiv \mp2 \bmod 11\) |

Half-integral levels occur exactly for \(d_K \equiv 1 \bmod 4\); the level \(\alpha = 2\) (the
inscribed circle of the ideal triangle, \(D = -4\)) belongs to \(\mathbb{Q}(\sqrt{-3})\), and
\(\alpha = 5/2\) (\(D = -3\), the centroid) to \(\mathbb{Q}(\sqrt{-7})\).

## 2. The hyperbolic dictionary and the class formula

Write \(k = 2\alpha \in \mathbb{Z}\) and \(D = D_K(\alpha) := 4(\alpha^2-1)/d_K = -(k-2)(k+2)/|d_K|\); by
Theorem 1(4) this is an integer, and \(D \equiv k^2 \equiv b^2 \pmod 4\) below.

> **Theorem 2 (level-\(\alpha\) circles are the forms of discriminant \(D_K(\alpha)\)).** Fix
> a level \(\alpha\) with orientation sign \(s\). Writing \(\beta = -s(\alpha + b\sqrt{d_K}/2)\), the
> level-\(\alpha\) circles in \(\mathbb{H}\) are exactly
> $$\omega_{a,b,c}:\quad q = sa,\ m = sc,\ \ \text{radius } \tfrac{1}{\sqrt{|d_K|}a},\ \text{centre } \tfrac{-b + i\alpha}{\sqrt{|d_K|}a},\qquad b^2 - 4ac = D,\ a \ge 1,$$
> and \(\omega_{a,b,c} \mapsto (a,b,c)\) is a bijection onto the positive definite forms of
> discriminant \(D\) (imprimitive included), equivariant for \(\mathrm{PSL}_2(\mathbb{Z})\); the
> hyperbolic centre is the CM point \(\tfrac{-b+\sqrt D}{2a}\) of the field \(\mathbb{Q}(\sqrt D)\),
> and the hyperbolic radius is \(\operatorname{arcoth}\alpha\). The weighted number of level-\(\alpha\)
> circles with hyperbolic centre in the ideal triangle is \(3H(|D|)\).

*Proof.* Identical to Paper I, Prop. 2.6 and Thm 3.2: with \(x = \operatorname{im}_K\beta = -sb\) the norm
condition \(N(\beta) = \alpha^2 + |d_K|x^2/4 = 1 + |d_K|ac\) is \(b^2 - 4ac = D\); the congruence
\(\beta \equiv 1\) is the level condition of Theorem 1(4) together with \(b \equiv k \bmod 2\), which
is forced by \(D \equiv k^2 \bmod 4\); the hyperbolic centre is \(x_0 + i\sqrt{y_0^2 - r^2} = \tfrac{-b + i\sqrt{|D|}}{2a}\);
and the incidence lemma of Paper I gives \(3H\). \(\square\)

**The twist ideal.** Put \(\rho = \alpha + s\), \(\rho' = \alpha - s\) (so \(\{\rho,\rho'\} = \{\alpha\pm1\}\)) and
$$
t_0 := \frac{2\rho}{|d_K|} = \frac{k + 2s}{|d_K|} \in \mathbb{Z}, \qquad
r_0 := \frac{k-2}{\gcd(k-2,|d_K|)},\quad s_0 := \frac{k+2}{\gcd(k+2,|d_K|)}\ \ (d_K \text{ odd}),
$$
$$
r_0 := \frac{\alpha-1}{\gcd(\alpha-1,|d_K|/2)},\quad s_0 := \frac{\alpha+1}{\gcd(\alpha+1,|d_K|/2)}\ \ (d_K \text{ even}),
$$
so that \(r_0s_0 = N(\theta_D)\) (\(= |D|\) for odd \(d_K\), \(|D|/4\) for even \(d_K\)), and \(t_0 = s_0\)
when \(s = +1\), \(t_0 = r_0\) when \(s = -1\) (the \(t_0\)-ideal is the one carrying the Gram form). Let
$$
\mathfrak{r}_\alpha := \bigl(r_0,\ \tfrac{r_0+\sqrt D}{2}\bigr)\ \text{(form } (r_0, r_0, \tfrac{r_0+s_0}{4})\text{)},\qquad
\mathfrak{r}_\alpha := \bigl(r_0,\ \tfrac{\sqrt D}{2}\bigr)\ \text{(form } (r_0,0,s_0)\text{)}
$$
for odd resp. even \(d_K\), and \(\mathfrak{s}_\alpha\) the same with \(r_0 \leftrightarrow s_0\). Both are
invertible ambiguous ideals of the order \(\mathcal{O}_D\) with \(\mathfrak{r}_\alpha^2 = (r_0)\),
\(\mathfrak{s}_\alpha^2 = (s_0)\), \(\mathfrak{r}_\alpha\mathfrak{s}_\alpha = (\theta_D)\), and \([\mathfrak{r}_\alpha] = [\mathfrak{s}_\alpha]\)
(their forms are primitive: for odd \(d_K\), \(r_0 + s_0 \equiv 0 \bmod 4\) and \((r_0+s_0)/4\) is
prime to \(r_0\) — a direct check on the \(2\)-part using \(\gcd(k-2,k+2) \mid 4\); certified at
every listed level). For \(\mathbb{Q}(i)\) this is Paper I's \(\mathfrak{r}_n = (\tfrac{n-1}2, \tfrac{\sqrt D}{2})\).

> **Theorem 3 (class formula).** Let \(\alpha\) be a level of \(\mathcal{S}_K\) with orientation
> sign \(s\), and \(f\) a primitive form of discriminant \(D = D_K(\alpha)\). Then the circle of
> \(\sigma(X) = \bar X^{-1}\) has the negative orientation and lies in the lower half-plane when
> \(s = +1\) (level \(-\alpha\)) and in \(\mathbb{H}\) when \(s = -1\) (level \(\alpha\)); reflecting it into
> \(\mathbb{H}\) in the first case, its class is
> $$\hat\sigma[f] \;=\; [\mathfrak{r}_\alpha]\,[f]^{-s} .$$
> For \(\mathbb{Q}(i)\) with the positive orientation this is Paper I's \([\mathfrak{r}_n][f]^{-1}\); with the
> negative orientation it is the pure twist \([\mathfrak{r}_nf]\) of the inner-disk normalization.

The proof is that of Paper I §4.3 with three general lemmas; we display the changes.
Fix \(f = (a,b,c)\), \(\beta = -s(\alpha + b\sqrt{d_K}/2)\), and write elements of \(\mathcal{O}_K\) as
\(u = x + y\tfrac{\sqrt{d_K}}{2}\) (\(x \in \tfrac12\mathbb{Z}\) or \(\mathbb{Z}\), \(y \in \mathbb{Z}\)).

**Lemma A′ (the explicit unimodular basis).** Let
\(\mathcal{K}_f := \{u \in \mathcal{O}_K : \bar u - \beta u \in a\sqrt{d_K}\,\mathcal{O}_K\}\). Then
\([\mathcal{O}_K : \mathcal{K}_f] = a\), \(a\mathcal{O}_K \subseteq \mathcal{K}_f\), and for a basis \((u_1,u_2)\) of
\(\mathcal{K}_f\) with \(\operatorname{im}_K(u_1\bar u_2) = sa\) the numbers
\(v_k := s(\bar u_k - \beta u_k)/(a\sqrt{d_K}) \in \mathcal{O}_K\) give
\(P = \binom{u_1\ v_1}{u_2\ v_2} \in \mathrm{SL}_2(\mathcal{O}_K)\) with \(P^\dagger M_0P = M_{\omega_{a,b,c}}\),
i.e. \(X := P^{-1}\) realizes the circle of \(f\).

*Proof.* The Gram entries are \(h_0(p_1,p_1) = 2\operatorname{Im}(u_1\bar u_2) = \sqrt{|d_K|}sa\),
\(h_0(p_1,p_2) = i(\bar u_1v_2 - \bar u_2v_1)\) and \(\det P = u_1v_2 - u_2v_1\); solving the two
linear equations \(u_1v_2 - u_2v_1 = 1\), \(\bar u_1v_2 - \bar u_2v_1 = \beta\) gives exactly the
displayed \(v_k\) (the determinant of the system is \(2i\operatorname{Im}(u_1\bar u_2) = sa\sqrt{d_K}\)),
and \(h_0(p_2,p_2)\) is then forced by \(\det(P^\dagger M_0P) = -1\). For the index: multiplying
by \(s\), \(s(\bar u - \beta u) = X' + Y'\tfrac{\sqrt{d_K}}{2}\) with
$$X' = \rho x - \tfrac{|d_K|}{4}by,\qquad Y' = \rho'y + bx,$$
and \(\sqrt{d_K}\mathcal{O}_K = \{-\tfrac{|d_K|}{2}y' + 2x'\tfrac{\sqrt{d_K}}{2}\}\); so \(u \in \mathcal{K}_f\) iff
\(X' \in \tfrac{a|d_K|}{2}\mathbb{Z}\) and \(Y' \in 2a\mathbb{Z}\) (even \(d_K\)), resp. \(Y' \in a\mathbb{Z}\) with
\(Y'/a \equiv 2X'/(a|d_K|) \bmod 2\) (odd \(d_K\)) — Paper I's two congruences. Since
\(\bar v \equiv v \equiv \beta v \bmod \sqrt{d_K}\), \(a\mathcal{O}_K \subseteq \mathcal{K}_f\). The map
\(\gamma' := (\bar u - \beta u)/\sqrt{d_K}\) is \(\mathbb{Z}\)-linear of determinant \((N(\beta)-1)/|d_K| = ac\);
in the coordinates \((x, y)\) its matrix is \(s\binom{b/2\ \ \rho'/2}{-t_0\ \ b/2}\) (even \(d_K\); the odd
case is the same with the doubled coordinates), whose content is \(1\) because an odd prime
dividing \(t_0\) and \(\rho'/2\) would divide \(\rho\) and \(\rho'\), hence \(\rho - \rho' = 2s\), while
\(2 \mid t_0, \rho'/2\) is excluded by \(\rho - \rho' = 2s\) at \(|d_K| \in \{4, 8\}\). So the Smith
invariants of \(\gamma'\) are \((1, ac)\), and \(\mathcal{K}_f/a\mathcal{O}_K = \ker(\gamma' \bmod a)\) has order
\(\gcd(1,a)\gcd(ac,a) = a\). \(\square\)

**Lemma B′ (the Gram form of \(\sigma\)).** With \(g_s(u) := \alpha N(u) + s\operatorname{Re}(u^2)\) and its
polarization \(g_s(u,u')\), the circle of \(\sigma(X)\) is \(-\bar N\), \(N = X^\dagger M_0X\), with
$$N_{11} = \frac{2g_s(u_2)}{a\sqrt{|d_K|}},\qquad N_{22} = \frac{2g_s(u_1)}{a\sqrt{|d_K|}},\qquad
\operatorname{Re}(\beta') = \operatorname{Re}(\beta),$$
so the image has curvature \(-2g_s(u_2)/(a\sqrt{|d_K|}) < 0\) and level \(-s\alpha\); after the
reflection \(z \mapsto \bar z\) when \(s = +1\) (which replaces \(\beta' \mapsto -\bar\beta'\) and
preserves the reading of the form), the image form is
$$f' \;=\; \frac{2}{a|d_K|}\Bigl(g_s(u_2),\ -2g_s(u_1,u_2),\ g_s(u_1)\Bigr) \;=\; \frac{2}{a|d_K|}\,g_s\big|_{\mathcal{K}_f}
\ \text{ in the basis } (u_2, -u_1).$$

*Proof.* \(X = \binom{v_2\ -v_1}{-u_2\ u_1}\), so \(N_{11} = 2\operatorname{Im}(u_2\bar v_2)\) and \(N_{22} = 2\operatorname{Im}(u_1\bar v_1)\)
(Paper I, Lemma B), and \(u_k\bar v_k = is(u_k^2 - \bar\beta|u_k|^2)/(a\sqrt{|d_K|})\) has imaginary part
\(s(\operatorname{Re}(u_k^2) - \operatorname{Re}(\beta)|u_k|^2)/(a\sqrt{|d_K|}) = g_s(u_k)/(a\sqrt{|d_K|})\), using \(\operatorname{Re}\beta = -s\alpha\).
The level identity \(\operatorname{Re}(a\bar d - c\bar b) = \operatorname{Re}(a\bar d - b\bar c)\) is immediate; the middle
coefficient is the polarization by the same computation. In the Gaussian case
\(g_+ = (n+1)s^2 + (n-1)t^2\) is Paper I's \(g\). \(\square\)

**Lemma C′ (\(\mathcal{K}_f\) is an ideal and factors).** Let \(\iota: \mathcal{O}_K \to K' = \mathbb{Q}(\sqrt D)\),
\(\iota(x + y\tfrac{\sqrt{d_K}}{2}) := t_0x + \tfrac{y}{2}\sqrt D\). Then
\(N_{K'}(\iota(u)) = t_0\cdot\tfrac{2}{|d_K|}g_s(u)\); \(\mathfrak{t} := \iota(\mathcal{O}_K)\) is the ideal
\(\mathfrak{s}_\alpha\) (\(s = +1\)) resp. \(\mathfrak{r}_\alpha\) (\(s = -1\)) of norm \(t_0\); \(\iota\) intertwines
multiplication by \(\sqrt D\) with \(\phi(x + y\tfrac{\sqrt{d_K}}{2}) = -\rho'y + 2t_0x\tfrac{\sqrt{d_K}}{2}\); and
for primitive \(f\),
$$\iota(\mathcal{K}_f) \;=\; \mathfrak{t}\,\mathfrak{a}_f,\qquad \mathfrak{a}_f = \mathbb{Z}a + \mathbb{Z}\tfrac{b+\sqrt D}{2}.$$

*Proof.* \(g_s(u) = \rho x^2 + \rho'|d_K|y^2/4 = \rho\,(x^2 + \tfrac{\rho'^2}{|D|}y^2)\) since
\(|D| = 4\rho\rho'/|d_K|\), and \(N(t_0x + \tfrac y2\sqrt D) = t_0^2x^2 + \tfrac{|D|}{4}y^2 = \tfrac{t_0^2}{\rho}g_s(u)\)
using \(\tfrac{|D|}{4} = \tfrac{t_0\rho'}{2}\); this is the norm identity. \(\mathfrak{t}\) is
\(\mathbb{Z}t_0 + \mathbb{Z}\tfrac{t_0+\sqrt D}{2}\) (odd \(d_K\)) or \(\mathbb{Z}t_0 + \mathbb{Z}\tfrac{\sqrt D}{2}\) (even \(d_K\)), the
ideal of norm \(t_0\) of the form \((t_0,t_0,\cdot)\) resp. \((t_0,0,\cdot)\). The intertwining is the
computation \(\sqrt D(t_0x + \tfrac y2\sqrt D) = -\tfrac{|D|}{2}y + t_0x\sqrt D\) with \(\tfrac{|D|}{2t_0} = \rho'\).
*Stability.* For \(u \in \mathcal{K}_f\) the coordinates \((X',Y')\) of Lemma A′ transform under \(\phi\) as
$$(X', Y') \;\longmapsto\; \bigl(-\rho Y',\ \tfrac{4\rho'}{|d_K|}X'\bigr),$$
(using \(\tfrac{|d_K|}{2}t_0 = \rho\)), which preserves the divisibility conditions
(\(2a\rho = a|d_K|t_0\), \(\tfrac{4\rho'}{|d_K|}\cdot\tfrac{a|d_K|}{2} = 2a\rho'\); for odd \(d_K\) the parity
condition is preserved because \(2\rho' \equiv t_0 \equiv k \bmod 2\)); hence \(\phi(\mathcal{K}_f) \subseteq \mathcal{K}_f\)
and \(\iota(\mathcal{K}_f)\) is \(\mathcal{O}_D\)-stable (the extra generator \(\tfrac{D+\sqrt D}{2}\) for odd \(D\) is
a parity check of the same kind). *Containment.* \(\mathfrak{t}\mathfrak{a}_f = a\mathfrak{t} + \tfrac{b+\sqrt D}{2}\mathfrak{t}
= \iota(a\mathcal{O}_K) + \iota(\psi(\mathcal{O}_K))\) with \(\psi = \tfrac{b+\phi}{2}\); \(a\mathcal{O}_K \subseteq \mathcal{K}_f\) by
Lemma A′, and for \(u = x + y\tfrac{\sqrt{d_K}}{2}\) the coordinates of \(\psi(u)\) are
\(X' = -\tfrac{|d_K|ac}{2}y\), \(Y' = 2acx\) (using \(4\rho\rho' + |d_K|b^2 = 4|d_K|ac\) and
\(2t_0\rho' + b^2 = 4ac\)), which satisfy the conditions. *Equality.* Both sides have index \(a\) in
\(\mathfrak{t}\): \([\mathfrak{t}:\iota(\mathcal{K}_f)] = [\mathcal{O}_K:\mathcal{K}_f] = a\) and \([\mathfrak{t}:\mathfrak{t}\mathfrak{a}_f] = N(\mathfrak{a}_f) = a\)
because \(\mathfrak{a}_f\) is proper (\(f\) primitive). \(\square\)

*Proof of Theorem 3.* By Lemmas B′ and C′, \(f'\) is the norm form of the proper ideal
\(\iota(\mathcal{K}_f) = \mathfrak{t}\mathfrak{a}_f\) of norm \(t_0a\) in the basis \((\iota(u_2), \iota(-u_1))\), whose
orientation is \(\operatorname{Im}(\overline{\iota(u_2)}\iota(-u_1)) = -\tfrac{t_0\sqrt{|D|}}{2}\operatorname{im}_K(u_1\bar u_2) = -\tfrac{t_0\sqrt{|D|}}{2}sa\):
negative for \(s = +1\), positive for \(s = -1\). Hence \([f'] = [\mathfrak{t}\mathfrak{a}_f]^{-s} = [\mathfrak{t}]^{-s}[f]^{-s}
= [\mathfrak{r}_\alpha][f]^{-s}\), as \([\mathfrak{t}] = [\mathfrak{r}_\alpha] = [\mathfrak{s}_\alpha]\) is \(2\)-torsion. \(\blacksquare\)

**Certified** (phase B): at the first ten levels of each field (all classes, both
orientations for \(\mathbb{Q}(i)\)), the circle of \(\sigma(X)\) for the descent matrix \(X\) has the
stated orientation and level, its class is \([\mathfrak{r}_\alpha][f]^{-s}\) at every primitive class,
and Lemmas A′ (index, integrality of \(v_k\), \(\det P = 1\), \(P^{-1}\) realizes the circle), B′
(the image form equals \(\tfrac{2}{a|d_K|}g_s\) in the basis \((u_2,-u_1)\)) and C′
(\(\iota(\mathcal{K}_f) = \mathfrak{t}\mathfrak{a}_f\) as lattices in \(K'\), and the norm identity) hold at every
class — \(30, 26, 26, 30\) primitive classes at \(d_K = -3,-7,-8,-11\) and \(98\) at \(\mathbb{Q}(i)\). The
twist classes at the first levels (\(D\), \(h\), the reduced form of \(\mathfrak{r}_\alpha\), \(r_0\), \(s_0\)):

| \(K\) | \(2\alpha\) | \(s\) | \(D\) | \(h\) | \([\mathfrak{r}_\alpha]\) | \(r_0\) | \(s_0\) |
|---|---|---|---|---|---|---|---|
| \(\mathbb{Q}(\sqrt{-2})\) | 6, 10, 14, 18, 22, 26, 30 | \(+,-,+,-,+,-,+\) | \(-4,-12,-24,-40,-60,-84,-112\) | 1,1,2,2,2,4,2 | \(1, 1, (2,0,3), (2,0,5), (3,0,5), (3,0,7), (4,0,7)\) | 1,1,3,2,5,3,7 | 1,3,2,5,3,7,4 |
| \(\mathbb{Q}(\sqrt{-3})\) | 4, 5, 7, 8, 10, 11, 13, 14, 16, 17 | \(+,-,+,-,\dots\) | \(-4,-7,-15,-20,-32,-39,-55,-64,-84,-95\) | 1,1,2,2,2,4,4,2,4,8 | \(1,1,(2,1,2),(2,2,3),(3,2,3),(3,3,4),(4,3,4),(4,4,5),(5,4,5),(5,5,6)\) | 2,1,5,2,8,3,11,4,14,5 | 2,7,3,10,4,13,5,16,6,19 |
| \(\mathbb{Q}(\sqrt{-7})\) | 5, 9, 12, 16, 19, 23, 26, 30 | \(+,-,\dots\) | \(-3,-11,-20,-36,-51,-75,-96,-128\) | 1,1,2,2,2,2,4,4 | \(1,1,(2,2,3),(2,2,5),(3,3,5),(3,3,7),(4,4,7),(4,4,9)\) | 3,1,10,2,17,3,24,4 | 1,11,2,18,3,25,4,32 |
| \(\mathbb{Q}(\sqrt{-11})\) | 9, 13, 20, 24, 31, 35 | \(+,-,\dots\) | \(-7,-15,-36,-52,-87,-111\) | 1,2,2,2,6,8 | \(1,1,(2,2,5),(2,2,7),(3,3,8),(3,3,10)\) | 7,1,18,2,29,3 | 1,15,2,26,3,37 |

(At \(2\alpha = 13\) over \(\mathbb{Q}(\sqrt{-11})\), \(D = -15\), the twist is trivial although \(h = 2\):
\(\mathfrak{r} = (1,\cdot)\) is principal and \(\hat\sigma[f] = [f]\) — the first level where the
involution fixes every class.)

## 3. The unit theorem

> **Theorem 4 (unit theorem over \(K\)).** For every level \(\alpha\) and every form \(f\) of
> discriminant \(D = D_K(\alpha)\) (primitive or not),
> $$R_f \;:=\; r_0^{\,6}\,\frac{\Delta(\mathfrak{b}_f)}{\Delta(\mathfrak{r}_\alpha^{-1}\mathfrak{b}_f)},\qquad \mathfrak{b}_f = \mathfrak{a}_f,$$
> is an algebraic unit; \(R_{\mathfrak{r}f} = R_f^{-1}\), \(R_{f^{-1}} = \overline{R_f}\), the level polynomial
> \(\prod_f(x - R_f)\) over the primitive classes is a palindromic integer polynomial with
> constant term \(1\), and \(R_f\) is unchanged if \(\mathfrak{r}_\alpha\) is replaced by \(\mathfrak{s}_\alpha\). For
> every nontrivial character \(\chi\) of \(\mathrm{Cl}(D)\),
> $$\sum_f\chi(f)\log|R_f| \;=\; \bigl(1 - \chi(\mathfrak{r}_\alpha)\bigr)\bigl(-12L'(0,\chi)\bigr)
> \;=\; \begin{cases}-24L'(0,\chi), & \chi(\mathfrak{r}_\alpha) = -1,\\ 0, & \chi(\mathfrak{r}_\alpha) = +1.\end{cases}$$

*Proof.* The unit statement is Paper II, Theorem 4.2, whose proof uses only that the twist
is an invertible ambiguous ideal with \(\mathfrak{r}^2 = (r_0)\): the \(p\)-part lemma, the rigidity of
\(p\)-divisible-group quotients and the balance \(2d_p = k\) are verbatim, at every prime
including the ramified primes of \(\mathbb{Q}(\sqrt D)\) and \(p = 2, 3\), because they are statements
about the invertible \(\mathcal{O}_D\)-ideal \(\mathfrak{r}\) alone; on an imprimitive stratum the twist is
\(\mathfrak{r}\mathcal{O}'\). The laws: \(\mathfrak{r}^{-1}\mathfrak{b}_{\mathfrak{r}f} = \mathfrak{r}^{-2}\mathfrak{b}_f = r_0^{-1}\mathfrak{b}_f\)
up to the homothety fixing the class, so \(R_{\mathfrak{r}f} = r_0^6\Delta(\mathfrak{r}^{-1}\mathfrak{b}_f)/\Delta(r_0^{-1}\mathfrak{b}_f)
= R_f^{-1}\); conjugation is \(\bar{\mathfrak{r}} = \mathfrak{r}\). Independence of the choice: \(\mathfrak{s}^{-1} = \mathfrak{r}/\theta_D\), so
\(\Delta(\mathfrak{s}^{-1}\mathfrak{b}) = \theta_D^{12}\Delta(\mathfrak{r}\mathfrak{b}) = \theta_D^{12}r_0^{-12}\Delta(\mathfrak{r}^{-1}\mathfrak{b})\) and
\(s_0^6\Delta(\mathfrak{b})/\Delta(\mathfrak{s}^{-1}\mathfrak{b}) = s_0^6r_0^{12}(r_0s_0)^{-6}\Delta(\mathfrak{b})/\Delta(\mathfrak{r}^{-1}\mathfrak{b}) = R_f\)
(\(\theta_D^{12} = N(\theta_D)^6\)). The character sum: \(\log|R_f| = \log g(f) - \log g(\mathfrak{r}f)\)
because the covolumes \(r_0^{-1}\) and the factor \(r_0^6\) cancel, and Paper II, Prop. 2.2
(valid for every discriminant) gives \(\sum\chi\log g = -12L'(0,\chi)\). \(\square\)

**Certified** (phase C, 200 digits): at the levels \(2\alpha \le 18\) of each field the level
polynomials are integer, palindromic, with constant term \(1\) (spare digits \(\ge 190\)); the
imprimitive strata are units as well; the two laws and the character-sum identity hold
with the independent evaluation of \(L'(0,\chi)\) (Paper II, Lemma 2.3). The first
polynomials:

| \(K\) | \(2\alpha\) | \(D\) | \(\prod_f(x - R_f)\) |
|---|---|---|---|
| \(\mathbb{Q}(\sqrt{-2})\) | 6, 10 | \(-4, -12\) | \(x - 1\) |
| | 14 | \(-24\) | \(x^2 - 34x + 1\) |
| | 18 | \(-40\) | \(x^2 - 322x + 1\) |
| | 22 | \(-60\) | \(x^2 - 15127x + 1\) |
| \(\mathbb{Q}(\sqrt{-3})\) | 4 | \(-4\) | \(x + 1\) (\(R = -1\): \(\mathfrak{r} = (2, 1+i) = (1+i)\) is principal with \((1+i)^{12} = -2^6 = -r_0^6\)) |
| | 7 | \(-15\) | \(x^2 - 7x + 1\) |
| | 8 | \(-20\) | \(x^2 + 18x + 1\) |
| | 10 | \(-32\) | \(x^2 + 198x + 1\) |
| | 11 | \(-39\) | \(x^4 - 679x^3 - 831x^2 - 679x + 1\) |
| | 13 | \(-55\) | \(x^4 - 9959x^3 - 13359x^2 - 9959x + 1\) |
| | 14 | \(-64\) | \(x^2 + 39202x + 1\) |
| | 16 | \(-84\) | \(x^4 + 686308x^3 + 7365318x^2 + 686308x + 1\) |
| \(\mathbb{Q}(\sqrt{-7})\) | 12 | \(-20\) | \(x^2 + 18x + 1\) |
| | 16 | \(-36\) | \(x^2 + 194x + 1\) |
| | 19 | \(-51\) | \(x^2 - 4354x + 1\) |
| \(\mathbb{Q}(\sqrt{-11})\) | 13 | \(-15\) | \((x-1)^2\) (trivial twist) |
| | 20 | \(-36\) | \(x^2 + 194x + 1\) |
| | 24 | \(-52\) | \(x^2 + 1298x + 1\) |

**Horizontal coincidences (item 4 of the program).** \(R_f\) depends only on the pair
\((D, [\mathfrak{r}_\alpha])\). The slice \(D = -24\), \([\mathfrak{r}] = [(2,0,3)]\) is the level \(n = 5\) of
\(\mathbb{Q}(i)\) and the level \(\alpha = 7\) of \(\mathbb{Q}(\sqrt{-2})\), and the two unit polynomials coincide
(\(x^2 - 34x + 1\), Paper II's \(n = 5\)); \(D = -20\), \([(2,2,3)]\) is \(\alpha = 4\) of \(\mathbb{Q}(\sqrt{-3})\) and
\(\alpha = 6\) of \(\mathbb{Q}(\sqrt{-7})\) (both \(x^2 + 18x + 1\)); \(D = -36\), \([(2,2,5)]\) is \(\alpha = 8\) of
\(\mathbb{Q}(\sqrt{-7})\) and \(\alpha = 10\) of \(\mathbb{Q}(\sqrt{-11})\) (both \(x^2 + 194x + 1\)). The arrangements
realize different slices of one horizontal family; which pairs \((D, \mathfrak{r})\) are realized by
which \(K\) is exactly the question of item 4.

## 4. The Euclidean dictionary, the \(\Delta\)-data, the mass law

> **Theorem 5 (Euclidean dictionary; \(h_K = 1\)).**
> 1. Disks \(X(\mathbb{H})\) of curvature \(\sqrt{|d_K|}\,n\) (\(n \ge 1\), the orientation with \(q > 0\))
>    modulo translation by \(\mathcal{O}_K\) are in bijection, via the bottom row \((c,d) \mapsto \Lambda = \mathbb{Z}c + \mathbb{Z}d\),
>    with the primitive index-\(n\) sublattices \(\Lambda \subseteq \mathcal{O}_K\) (\(\mathcal{O}_K\Lambda = \mathcal{O}_K\)); every
>    such \(\Lambda\) is a proper \(\mathcal{O}_n\)-lattice (conductor \(=\) curvature), and
>    \(\Lambda \mapsto [\Lambda] \in \mathrm{Pic}(\mathcal{O}_n)\) is surjective with fibers \(\{u\Lambda : u \in \mathcal{O}_K^\times/\pm1\}\)
>    of size \(w_K/2\) for \(n \ge 2\). Hence
>    $$N_e(n) \;=\; \tfrac{w_K}{2}\,h(\mathcal{O}_n) \qquad(n \ge 2),$$
>    \(= 2h\) at \(\mathbb{Q}(i)\), \(3h\) at \(\mathbb{Q}(\sqrt{-3})\), \(h\) at \(\mathbb{Q}(\sqrt{-2}), \mathbb{Q}(\sqrt{-7}), \mathbb{Q}(\sqrt{-11})\).
>    (For \(K \ne \mathbb{Q}(i)\) the unoriented circles of curvature \(\sqrt{|d_K|}n\) number \(2N_e(n) = w_Kh\)
>    per translation class, each with its unique orientation; the two orientation signs of a lattice
>    are exchanged by \(z \mapsto -z\), which normalizes \(\mathrm{PSL}_2(\mathcal{O}_K)\) and preserves \(\mathcal{S}_K\).)
> 2. The \(\Delta\)-data \(G_\mathfrak{c} := n^{12}\Delta(\Lambda_\mathfrak{c})/\Delta(\mathcal{O}_K)\) are well defined on
>    \(\mathrm{Pic}(\mathcal{O}_n)\) (\(u^{12} = 1\) for every unit), are algebraic integers of \(H_n\) with
>    \(\sigma_\mathfrak{a}(G_\mathfrak{c}) = G_{\mathfrak{a}^{-1}\mathfrak{c}}\), \(\overline{G_\mathfrak{c}} = G_{\mathfrak{c}^{-1}}\), and
>    \(D_n(x) := \prod_\mathfrak{c}(x - G_\mathfrak{c}) \in \mathbb{Z}[x]\) with \(D_n(0) = (-1)^hM(n)\), \(M(n) := \prod_\mathfrak{c}G_\mathfrak{c}\).
> 3. *(Mass law.)* \(|M(n)| = \prod_{p^k\parallel n,\ p\text{ not split}}p^{\,\frac{24}{w_Ke_p}\cdot\frac{p^k-1}{p-1}\,N_e(n/p^k)}\),
>    \(e_p = 2\) at the ramified prime, \(1\) at inert primes; split primes contribute nothing.
> 4. *(Per-class valuations.)* For every place \(\mathfrak{P} \mid p\), \(p^k \parallel n\),
>    \(v_\mathfrak{P}(G_\mathfrak{c}) = w_p(k) := \dfrac{12\,(p^k-1)}{e_p\,(p-1)\,N_e(p^k)}\) independently of \(\mathfrak{c}\)
>    (\(= 0\) at split \(p\)): \(\dfrac{12(p^k-1)}{(p-1)p^{k-1}(p+1)}\) inert, \(\dfrac{6(p^k-1)}{(p-1)p^k}\) ramified.
>    Hence \(\mathcal{V}_n := \langle G_\mathfrak{c}/G_{\mathfrak{c}'}\rangle \subset \mathcal{O}_{H_n}^\times\).
> 5. *(Sign; certified.)* \(\operatorname{sgn}M(n) = \operatorname{sgn}\bigl(\Delta(\mathcal{O}_K)\bigr)^{h}\cdot\epsilon_K(n)\), where
>    \(\Delta(\mathcal{O}_K) < 0\) iff \(d_K \equiv 1 \bmod 4\) (\(\operatorname{Re}\omega = \tfrac12\)) and \(\epsilon_K(n) = -1\) iff
>    \(n\) is a power of the ramified prime of \(K\) (with \(n \ge 4\) when \(K = \mathbb{Q}(i)\)).

*Proof.* (1) Paper I, Lemmas 7.4–7.6 verbatim: a basis \((c,d)\) with \(\gcd(c,d) = 1\) completes
to \(\mathrm{SL}_2(\mathcal{O}_K)\) and the completions form one translation class; the index is
\(|\operatorname{im}_K(c\bar d)| = q\); the conductor computation is local at each \(p \mid n\) (primitivity
excludes exactly the ideal lines of \(\mathcal{O}_K/p\), which re-derives \(N_e\)); for a class pick
\(\mathfrak{a}\) prime to \(n\) and a generator of \(\mathfrak{a}\mathcal{O}_K\) (\(h_K = 1\)), and \(u\Lambda = \Lambda\) iff
\(u \in \mathcal{O}_n^\times = \{\pm1\}\) for \(n \ge 2\) (\(\zeta_6, \zeta_3 \notin \mathbb{Z}+n\mathcal{O}_K\)). (2) Paper II,
Theorem 2.6: the transport lemma is Shimura reciprocity for the isobaric ratio
\(n^{12}\Delta(\Lambda)/\Delta(\mathcal{O}_K)\), and the ambiguity \(\Lambda \mapsto u\Lambda\) is invisible. (3) Steps
M1 and M3 of Paper I, Thm 7.18 are field-independent (the full mass \(A(n) = (-1)^{t(n)}\prod_{d\mid n}d^{-12d}\);
\(v_p(\gamma(p^j)^{12}) = 6j(j+1), 6j[2\mid j], 6j\) for split, inert, ramified \(p\), \(e_p\) being the
ramification index; \(r * N_e = \sigma\) with \(r = 1 * \chi_K\)). In M2 the stratification
\(\Lambda = \delta\Lambda''\) uses \(h_K = 1\) and \(\gamma(m)^{12}\) is well defined since \(u^{12} = 1\); the product
\(Q(n)\) over the \(N_e(n)\) primitive lattices is \(\bigl(\prod_\mathfrak{c}\Delta(\Lambda_\mathfrak{c})/\Delta(\mathcal{O}_K)\bigr)^{w_K/2}\),
so \(v_p(M(n)) = 12hk + \tfrac{2}{w_K}v_p(Q(n))\) with \(v_p(Q)\) given by the same recursion as at
\(\mathbb{Q}(i)\) (where \(v_p(M) = 12hk + \tfrac12v_p(Q)\)); substituting \(h = \tfrac{2}{w_K}N_e(n)\) the terms
\(12hk\) cancel and the exponent \(\tfrac{24}{w_Ke_p}\tfrac{p^k-1}{p-1}N_e(n/p^k)\) remains. (4) Paper II,
Theorem 4.11 with \(\mathcal{O}_{K,p}\) in place of \(\mathbb{Z}_p[i]\): \(\Lambda_{\mathfrak{c},p} = \epsilon\mathcal{O}_{n,p}\) with
\(\epsilon \in \mathcal{O}_{K,p}^\times\) an automorphism of \(E_0[p^\infty]\), so the pairs
\((E_0[p^\infty], C_\mathfrak{c})\) are isomorphic; the value is \(v_p(M(n))/h\). (5) is certified at
\(n \le 12\) in every field; the mechanism is Paper I's M4 with the extra sign of the reference
lattice for \(d_K \equiv 1 \bmod 4\). \(\square\)

**Certified** (phase D, 150 digits, \(2 \le n \le 9\) per field, and \(n \le 12\) in exploration):
\(D_n \in \mathbb{Z}[x]\), \(|D_n(0)|\) equal to the mass law, the sign law, a single Newton slope
\(w_p(k)\) across the full degree at every non-split \(p \mid n\) and slope \(0\) at split \(p\),
\(\overline{G_\mathfrak{c}} = G_{\mathfrak{c}^{-1}}\). First polynomials and masses:

| \(K\) | \(n\) | \(h\) | \(D_n(x)\) | \(\lvert M(n)\rvert\), sign | slopes |
|---|---|---|---|---|---|
| \(\mathbb{Q}(\sqrt{-2})\) | 2 | 2 | \(x^2 + 112x - 64\) | \(2^6\), \(-\) | \(2: 3\) |
| | 3 | 2 | \(x^2 - 98x + 1\) | \(1\) | \(3: 0\) |
| | 4 | 4 | \(x^4 + 2816x^3 - 3963904x^2 + 5881593856x - 2^{18}\) | \(2^{18}\), \(-\) | \(2^2: 9/2\) |
| | 5, 7, 8 | 6, 8, 8 | | \(5^{12}, 7^{12}, 2^{42}\) (\(-\) at 8) | \(5: 2\); \(7: 3/2\); \(2^3: 21/4\) |
| \(\mathbb{Q}(\sqrt{-3})\) | 2 | 1 | \(x + 16\) | \(2^4\), \(-\) | \(2: 4\) |
| | 3 | 1 | \(x - 9\) | \(3^2\) | \(3: 2\) |
| | 4 | 2 | \(x^2 + 3328x + 4096\) | \(2^{12}\) | \(2^2: 6\) |
| | 6, 8, 9 | 3, 4, 3 | | \(2^{12}3^6\) (\(-\)), \(2^{28}\), \(3^8\) | \(2:4, 3:2\); \(2^3: 7\); \(3^2: 8/3\) |
| \(\mathbb{Q}(\sqrt{-7})\) | 2 | 1 | \(x + 1\) (\(G = -1\)) | \(1\), \(-\) | \(2: 0\) |
| | 3 | 4 | \(x^4 - 756x^3 - 223074x^2 - 16671501x + 3^{12}\) | \(3^{12}\) | \(3: 3\) |
| | 4 | 2 | \(x^2 + 4048x + 1\) | \(1\) | \(2: 0\) |
| | 7, 9 | 7, 12 | | \(7^6\), \(3^{48}\) | \(7: 6/7\); \(3^2: 4\) |
| \(\mathbb{Q}(\sqrt{-11})\) | 2 | 3 | \(x^3 + 48x^2 + 33536x + 4096\) | \(2^{12}\), \(-\) | \(2: 4\) |
| | 3 | 2 | \(x^2 - 2114x + 1\) | \(1\) | \(3: 0\) |
| | 4 | 6 | \(x^6 + 9984x^5 + \cdots + 2^{36}\) | \(2^{36}\) | \(2^2: 6\) |
| | 8, 11 | 12, 11 | | \(2^{84}\), \(11^6\) | \(2^3: 7\); \(11: 6/11\) |

(\(\mathbb{Q}(\sqrt{-2})\), \(n = 2\): \(H_2 = \mathbb{Q}(\zeta_8)\), the roots are \(G = 8(1+\sqrt2)^{-3}\) and
\(-8(1+\sqrt2)^{3}\), so \(\mathcal{V}_2 = \langle(1+\sqrt2)^6\rangle\) up to torsion and the index of §7 is \(6\), the
fundamental unit of \(\mathbb{Q}(\zeta_8)\) being \(1+\sqrt2\).)

## 5. The Kronecker limit formula and the genus characters

> **Theorem 6.** For \(n \ge 2\) and every nontrivial character \(\chi\) of \(\mathrm{Pic}(\mathcal{O}_n)\),
> $$\sum_{\mathfrak{c}}\chi(\mathfrak{c})\log|G_\mathfrak{c}| \;=\; -12\,L'(0,\chi),$$
> \(L(s,\chi) = \tfrac12\sum_\mathfrak{c}\chi(\mathfrak{c})\zeta_{Q_\mathfrak{c}}(s)\) the Epstein \(L\)-function of discriminant
> \(n^2d_K\) (\(w = 2\)). At a prime level \(n = p \nmid d_K\), \(\mathrm{Pic}(\mathcal{O}_p)\) has exactly one
> quadratic character (\(h_K = 1\)), it is primitive, and with \(\{d_1, d_2\} = \{p^*, d_Kp^*\}\),
> \(d_2 > 0\), \(p^* = (-1)^{(p-1)/2}p\),
> $$L'(0,\chi_2) \;=\; \frac{2h(d_1)}{w(d_1)}\,h(d_2)\log\varepsilon_{d_2},\qquad
> \text{real field } \mathbb{Q}(\sqrt{d_2}) = \begin{cases}\mathbb{Q}(\sqrt p), & p \equiv 1 \bmod 4,\\ \mathbb{Q}(\sqrt{p|d_K|}), & p \equiv 3 \bmod 4.\end{cases}$$
> For \(\mathbb{Q}(i)\) both cases give \(\mathbb{Q}(\sqrt p)\) — Paper II's "\(\mathbb{Q}(\sqrt n)\) phenomenon" is the
> case \(4p \sim p\) of the general rule; for \(\mathbb{Q}(\sqrt{-2})\) it is \(\mathbb{Q}(\sqrt{2p})\) at \(p \equiv 3 \bmod 4\).

*Proof.* The first statement is Paper II, Prop. 2.2 (valid for every discriminant) with the
classwise identity \(\log|G_\mathfrak{c}| = \log g_\mathfrak{c} + 6\log n - \log|\Delta_q(\omega)|\cdot\)const. At a prime
level the genus theory of the order \(\mathcal{O}_p\) (Cox, Thm 3.15 for orders; \(\mu = 2\) assigned
characters when \(|d_K|\) is a prime power) has one nontrivial genus character, the one attached
to the factorization \(p^2d_K = d_1d_2\) into the two fundamental discriminants \(p^*\) and
\(d_Kp^*\) (\(d_2 = 9\)-type factorizations give the trivial character); it is primitive by
[robert-index-full.md](robert-index-full.md) Lemma 3.1, so \(L(s,\chi_2) = L(s,\chi_{d_1})L(s,\chi_{d_2})\)
and Dirichlet's formulas \(L(0,\chi_{d_1}) = 2h(d_1)/w(d_1)\), \(L'(0,\chi_{d_2}) = h(d_2)\log\varepsilon_{d_2}\)
give the value. \(\square\)

**Certified** (phase E): the limit formula at every nontrivial character of every level
\(n \le 9\) of each field (independent incomplete-gamma evaluation, spare \(\ge 148\) at 150
digits); the genus closed form at the prime levels with PARI's class numbers and
fundamental units (`qfbclassno`, `quadunit`): \(\mathbb{Q}(\sqrt{-2})\): \(p = 3, 5, 7\) →
\(\mathbb{Q}(\sqrt{24}), \mathbb{Q}(\sqrt5), \mathbb{Q}(\sqrt{56})\) with \((d_1,d_2) = (-3,24), (-40,5), (-7,56)\);
\(\mathbb{Q}(\sqrt{-3})\): \(p = 5, 7\) → \(\mathbb{Q}(\sqrt5), \mathbb{Q}(\sqrt{21})\); \(\mathbb{Q}(\sqrt{-7})\): \(p = 3, 5\) →
\(\mathbb{Q}(\sqrt{21}), \mathbb{Q}(\sqrt5)\) (with \(h(-35) = 2\)); \(\mathbb{Q}(\sqrt{-11})\): \(p = 3, 5, 7\) →
\(\mathbb{Q}(\sqrt{33}), \mathbb{Q}(\sqrt5), \mathbb{Q}(\sqrt{77})\) (with \(h(-55) = 4\)).

## 6. The Euler system

> **Theorem 7 (norm relations over \(K\)).** Let \(n \ge 2\), \(\ell\) prime, \(\pi: \mathrm{Pic}(\mathcal{O}_{n\ell}) \to \mathrm{Pic}(\mathcal{O}_n)\).
> 1. *(Fiber lemma.)* Of the \(\ell+1\) index-\(\ell\) sublattices of \(\Lambda_\mathfrak{c}\), exactly
>    \(1 + \chi_K(\ell)\) are \(\mathcal{O}_n\)-stable when \(\ell \nmid n\) (\(\mathfrak{l}_n\Lambda\) for the primes \(\mathfrak{l} \mid \ell\))
>    and exactly one when \(\ell \mid n\) (\(\ell\Lambda^+\)); the others are the primitive representatives
>    of the fiber \(\pi^{-1}(\mathfrak{c})\), bijectively, and \(|\ker\pi| = \ell - \chi_K(\ell)\) resp. \(\ell\).
> 2. With \(N = N_{H_{n\ell}/H_n}\), \([\mathfrak{l}]\mathfrak{c}\) the class of \(\mathfrak{l}_n\Lambda_\mathfrak{c}\), \(\mathfrak{c}^+ = \pi(\mathfrak{c})\),
>    \(\mathfrak{p} = (\pi)\) the prime above the ramified \(\ell\) and \(G^{(1)} := 1\):
>    $$N(G_{\mathfrak{c}'}) = (-1)^{[\ell=2]}\cdot\begin{cases}
>    G_\mathfrak{c}^{\,\ell+1}/(G_{[\mathfrak{l}]\mathfrak{c}}G_{[\bar{\mathfrak{l}}]\mathfrak{c}}), & \ell \text{ split},\ \ell\nmid n,\\
>    \ell^{12}\,G_\mathfrak{c}^{\,\ell+1}, & \ell \text{ inert},\ \ell\nmid n,\\
>    \pi^{12}\,G_\mathfrak{c}^{\,\ell+1}/G_{[\mathfrak{p}]\mathfrak{c}}, & \ell \text{ ramified},\ \ell\nmid n,\\
>    G_\mathfrak{c}^{\,\ell+1}/G^{(n/\ell)}_{\mathfrak{c}^+}, & \ell \mid n,\end{cases}$$
>    where \(\pi^{12} \in \mathbb{Z}\) is independent of the generator: \(-64\) at \(\mathbb{Q}(i)\) and \(+64\) at
>    \(\mathbb{Q}(\sqrt{-2})\) (\(\ell = 2\)), \(3^6, 7^6, 11^6\) at \(\mathbb{Q}(\sqrt{-3}), \mathbb{Q}(\sqrt{-7}), \mathbb{Q}(\sqrt{-11})\).
> 3. *(Hecke recursion.)* For \(\chi \ne 1\) on \(\mathrm{Pic}(\mathcal{O}_n)\) and \(\ell \nmid n\),
>    \(L(s,\chi^{(n\ell)}) = \bigl(1 + \ell^{1-2s} - \ell^{-s}(\chi(\mathfrak{l}) + \bar\chi(\mathfrak{l}))\bigr)L(s,\chi)\) (split),
>    \((1 + \ell^{1-2s})L(s,\chi)\) (inert), \((1 + \ell^{1-2s} - \ell^{-s}\chi(\mathfrak{p}))L(s,\chi)\) (ramified), and
>    the \(\ell \mid n\) relation of [robert-index-full.md](robert-index-full.md) Lemma 3.2; at \(s = 0\)
>    the multipliers are \(P_\ell(\chi) = \ell + 1 - \chi(\mathfrak{l}) - \bar\chi(\mathfrak{l})\), \(\ell+1\), \(\ell + 1 - \chi(\mathfrak{p})\).

*Proof.* (1) is [schmidt-euler-system.md](schmidt-euler-system.md) Lemma 1.2, which is local:
the lines of \(\mathcal{O}_{n,\ell}/\ell\) and the unit lines; the only global input is the injectivity
step, where \(\Lambda'_2 = u\Lambda'_1\) with \(u \in \mathcal{O}_K^\times \setminus\{\pm1\}\) forces
\(\Lambda'_1 \subseteq \Lambda\cap u^{-1}\Lambda = n\mathcal{O}_K\) — for \(u = \zeta_6, \zeta_3\) the local computation
\(\mathcal{O}_{n,p}\cap u\mathcal{O}_{n,p} = n\mathcal{O}_{K,p}\) is the same as for \(u = i\) (the \(\omega\)-coordinate of \(u\)
is \(\pm1\)). (2) is the proof of Theorem 1 there with \(A(\ell) = (-1)^{[\ell=2]}\ell^{-12\ell}\): the sign of
\(A(2)\) is global and appears in every case (at \(\mathbb{Q}(i)\) only the ramified and the \(2 \mid n\) cases
exist), and at a ramified \(\ell\) the stable sublattice contributes
\(\Delta(\mathfrak{p}_n\Lambda)/\Delta(\Lambda) = \pi^{-12}G_{[\mathfrak{p}]\mathfrak{c}}/G_\mathfrak{c}\). (3) is the theta-series
identity of Lemma 3.2 there, verbatim. \(\square\)

**Certified** (phase F, 80 digits): the fiber lemma (stable sublattices counted and the
non-stable ones matched to the fiber) and the norm relation at every class of \(\ge 20\)
pairs \((n,\ell)\) per field covering the four cases (\(26, 27, 25, 21\) pairs at \(d_K = -8, -3, -7, -11\),
top level \(60\); the ramified prime tested at \(n = 3, 5, 7, 9, 11\) over \(\mathbb{Q}(\sqrt{-2})\), the inert
\(\ell = 2\) with its sign over \(\mathbb{Q}(\sqrt{-3})\) and \(\mathbb{Q}(\sqrt{-11})\), the split \(\ell = 2\) over
\(\mathbb{Q}(\sqrt{-7})\)), with \(\ge 77\) spare digits; the Hecke recursion at the first level with \(h \ge 2\)
for two or three primes, the multipliers certified as integers: \(P_3 = 6\), \(P_5 = 6\) at \(2 \to 6, 10\)
over \(\mathbb{Q}(\sqrt{-2})\) (\(3\) split with \(\chi(\mathfrak{l}) = -1\), \(5\) inert); \(P_3 = 5 = 3 + 1 - \chi(\mathfrak{p})\) at
\(4 \to 12\) over \(\mathbb{Q}(\sqrt{-3})\) (the ramified multiplier with \(\chi(\mathfrak{p}) = -1\)); \(P_2 \in \{3, 5\}\) at
\(3 \to 6\) over \(\mathbb{Q}(\sqrt{-7})\); \(P_3 = 5\), \(P_5 = 7\) at \(2 \to 6, 10\) over \(\mathbb{Q}(\sqrt{-11})\).

## 7. The full Robert index

> **Lemma 7.1 (roots of unity of \(H_n\)).** \(\mu(H_n) = \mu_w\) with \(w = 2\cdot a\cdot b\), where
> \(b = 3\) iff \(\zeta_3 \in H_n\) iff \(K = \mathbb{Q}(\sqrt{-3})\) or \(3 \mid n\); and \(a \in \{1, 2, 4\}\) with
> \(\zeta_4 \in H_n\) iff \(K = \mathbb{Q}(i)\), or \(K = \mathbb{Q}(\sqrt{-2})\) and \(2 \mid n\), or \(4 \mid n\);
> \(\zeta_8 \in H_n\) iff \(K = \mathbb{Q}(i)\) and \(4 \mid n\), or \(K = \mathbb{Q}(\sqrt{-2})\) and \(2 \mid n\), or \(8 \mid n\).
> So \(w_{H_n} = 2\cdot4^{[2\mid n]}3^{[3\mid n]}\) at \(\mathbb{Q}(\sqrt{-2})\), \(6\cdot2^{[4\mid n]}2^{[8\mid n]}\) at
> \(\mathbb{Q}(\sqrt{-3})\), \(2\cdot2^{[4\mid n]}2^{[8\mid n]}3^{[3\mid n]}\) at \(\mathbb{Q}(\sqrt{-7}), \mathbb{Q}(\sqrt{-11})\).

*Proof.* As [robert-index-full.md](robert-index-full.md) Lemma 4.1: \(\zeta_m \in H_n\) forces
\(m \mid 24\); \(K(\zeta_m) \subseteq H_n\) iff the conductor of \(K(\zeta_m)/K\) divides \(n\mathcal{O}_K\)
(generalized dihedral is automatic). The conductors, from the discriminants of the
biquadratic fields \(\mathbb{Q}(\sqrt{d_K}, \sqrt{-3})\), \(\mathbb{Q}(\sqrt{d_K}, i)\), \(\mathbb{Q}(\sqrt{d_K},\sqrt2)\)
(\(d_{L} = d_1d_2d_3\)) divided by \(d_K^2\): \(K(\sqrt{-3})/K\) has conductor \((3)\) for \(K \ne \mathbb{Q}(\sqrt{-3})\);
\(K(i)/K\) has conductor \((4)\) except at \(\mathbb{Q}(\sqrt{-2})\), where \(K(i) = \mathbb{Q}(\zeta_8)\) has relative
discriminant of norm \(2^8/2^6 = 4\), i.e. conductor \(\mathfrak{p}_2^2 = (2)\), and \(\zeta_8 \in K(i)\);
\(K(\sqrt2)/K\) has conductor \((8)\) for \(d_K\) odd. \(\square\)

> **Theorem 8 (the full Robert index over \(K\), \(h_K = 1\)).** For every \(n \ge 2\),
> \(\mathcal{V}_n\) is free of rank \(h - 1\), meets \(\mu(H_n)\) trivially, and
> $$\bigl[\mathcal{O}_{H_n}^\times : \mu(H_n)\mathcal{V}_n\bigr] \;=\; \frac{24^{h-1}\prod_{\chi\ne1}|L'(0,\chi)|}{R_{H_n}}
> \;=\; \frac{w_K\cdot24^{h-1}}{w_{H_n}}\;h_{H_n}\prod_{\chi\ne1}C_\chi(0),$$
> with the multipliers \(C_\chi(0)\) of Theorem 7(3) (positive integers in product). The layer
> theorem holds for every subgroup \(A \le \mathrm{Pic}(\mathcal{O}_n)\) with the same constant \(w_K/w_{H^A}\);
> on a cubic layer with \(|M(n)|\) a cube, \([\mathcal{O}_{L_3}^\times:\langle-1,\theta_u\rangle] = 8h_{L_3}C_{\chi_3}(0)\)
> in the real cubic field \(L_3\) and \([\mathcal{O}_F^\times:\mu_F\mathcal{V}^A] = 576\,h_FC_{\chi_3}(0)^2\) in the
> sextic \(F = K(\theta_u)\) (\(w_F = w_K\)), independently of \(K\).

*Proof.* [robert-index-full.md](robert-index-full.md) Theorem 1 verbatim: the Dedekind
determinant (Lemma 2.1 there), the nonvanishing through the multipliers (Cor. 3.3, with
the ramified multiplier \(\ell + 1 - \chi(\mathfrak{p}) \in \{\ell, \ell+2\}\)), and the class number formula
\(\zeta_{H_n} = \zeta_K\prod_{\chi\ne1}L_{\mathrm{prim}}(s,\chi)\) at \(s = 0\) with
\(\zeta_K(0) = -h_K/w_K = -1/w_K\): the constant \(4\) of the Gaussian theorem is \(w_K\). For the cubic
layer \(\log|\theta_u| = -8L'(0,\chi_3)\) and \(\zeta_{L_3} = \zeta\cdot L_{\mathrm{prim}}(s,\chi_3)\) give
\(L'_{\mathrm{prim}}(0,\chi_3) = h_{L_3}R_{L_3}\) with no root-of-unity constant; \(\mu(F) = \mu(K)\) since
\([F:K]\) is odd. \(\square\)

**Certified with PARI/GP** (phase G; `polclass(n^2d_K)`, `polcompositum` with \(y^2 - d_K\),
`bnfinit` at 120 digits, `bnfcertify` at degree \(\le 16\), `nfroots` of \(D_n\), `bnfisunit`
exponent matrices, Smith forms): at every listed level the determinant identity (V1), the
limit formula (V2), the primitive levels and multipliers (V3), \(w_{H_n}\) against Lemma 7.1
(through `nfroots` of \(x^2+1, x^2+3, x^2-2\)), the index as a certified integer (spare
\(\ge 105\)) equal to the formula (V4), the exact index \(|\det E_n|\) (V5), the \(24\)-th-root
saturation \([\mathcal{O}^\times:\mathcal{W}_n] = h_{H_n}\prod C_\chi(0)\) (Conjecture 5.3 of the index document,
now seen away from \(\mathbb{Q}(i)\)), the quadratic layers (\(\theta^{(2)}_u = \pm\varepsilon_{d_2}^{e}\), \(|e| = 6L'(0,\chi)/\log\varepsilon\))
and one cubic layer per field (V6). The record is in §8 and §11.

## 8. The four fields

**\(\mathbb{Q}(\sqrt{-2})\)** (\(d_K = -8\), \(w_K = 2\), \(\mathcal{O}_K = \mathbb{Z}[\sqrt{-2}]\), \(2\) ramified, \(3\) split, \(5, 7\) inert).
The cleanest transport: all odd levels \(\alpha \ge 3\), orientation by \(\alpha \bmod 4\); twist
\((r_0, 0, s_0)\) with \(\{r_0,s_0\} = \{\tfrac{\alpha-1}{2},\tfrac{\alpha+1}{4}\}\) or \(\{\tfrac{\alpha-1}{4},\tfrac{\alpha+1}{2}\}\);
\(N_e(n) = h(-8n^2)\); mass exponents \(\tfrac{12}{e_p}\tfrac{p^k-1}{p-1}N_e(n/p^k)\), sign \(-\) exactly at
\(n = 2^k\); genus fields \(\mathbb{Q}(\sqrt6), \mathbb{Q}(\sqrt5), \mathbb{Q}(\sqrt{14})\) at \(p = 3, 5, 7\); \(\pi^{12} = 64\) at
\(\ell = 2\) with the global sign, so \(N(G_{\mathfrak{c}'}) = -64\,G_\mathfrak{c}^3/G_{[\mathfrak{p}]\mathfrak{c}}\) — the negative of
the Gaussian relation; \(\mu(H_n) = \mu_{8}\) at even \(n\) (\(H_2 = \mathbb{Q}(\zeta_8)\), index \(6 = 2\cdot24/8\)).
The cubic level \(n = 5\) (\(\mathrm{Pic} = \mathbb{Z}/6\)): index \(24^5 = 7962624\), \(h_{H_5} = 1\), \(w = 2\);
\(L_3 = \mathbb{Q}[y]/(y^3 - y^2 + 2y + 2)\), \(h_{L_3} = 1\), \(\theta_u = \pm\eta^8\); sextic index \(576\).

**\(\mathbb{Q}(\sqrt{-3})\)** (\(d_K = -3\), \(w_K = 6\), \(3\) ramified, \(2, 5, 11\) inert, \(7, 13\) split). The
sharpest test of the root-of-unity bookkeeping: levels \(2\alpha \not\equiv 0 \bmod 3\) including the
half-integral ones, lines in three directions, \(N_e(n) = 3h(-3n^2)\), mass exponents
\(\tfrac{4}{e_p}\tfrac{p^k-1}{p-1}N_e(n/p^k)\) (\(M(2) = -16\), \(M(3) = 9\): \(G_1 = 3^{12}\Delta(3\omega)/\Delta(\omega) = 9\)),
\(\Delta(\mathcal{O}_K) < 0\) so the sign law carries \((-1)^h\); the ramified prime gives the per-class
valuations \(w_3(k) = 3(3^k-1)/3^k\); \(\mu_6 \subset \mu(H_n)\) always, \(\mu_{12}\) at \(4 \mid n\), \(\mu_{24}\)
at \(8 \mid n\); the \(w_K = 6\) in the index formula: at \(n = 4\) (\(\mathrm{Pic} = \mathbb{Z}/2\), \(H_4 = \mathbb{Q}(\zeta_{12})\))
the index is \(6\cdot24/12 = 12\). The first cubic level is \(n = 6\) (\(\mathrm{Pic} = \mathbb{Z}/3\), \(|M(6)| = 2^{12}3^6\)
a cube): \(L_3 = \mathbb{Q}(\sqrt[3]2)\) (PARI: \(y^3 - 2\)), \(h_{L_3} = 1\), \(\theta_u = \theta/144 = \pm\eta^8\),
\(H_6 = K(\sqrt[3]2)\) with \(h = 1\), \(w = 6\), full index \(6\cdot24^2/6 = 576\), sextic index \(576\).

**\(\mathbb{Q}(\sqrt{-7})\)** (\(d_K = -7\), \(w_K = 2\), \(7\) ramified, \(2, 11\) split, \(3, 5\) inert). Half-integral
levels \(5/2, 9/2, 6, 8, 19/2, \dots\); the level \(\alpha = 5/2\) has \(D = -3\) (the centroid);
\(N_e(n) = h(-7n^2)\), \(M(2) = -1\) (\(G_1 = -1\) at \(n = 2\): \(2\) splits, and \(\Delta(\mathcal{O}_K) < 0\)); the
ramified prime contributes \(w_7(1) = 6/7\) (\(M(7) = 7^6\), \(h = 7\)); genus fields
\(\mathbb{Q}(\sqrt{21})\) at \(p = 3\) (with \(h(-3) = 1\), \(w = 6\)) and \(\mathbb{Q}(\sqrt5)\) at \(p = 5\) (\(h(-35) = 2\));
the cubic level \(n = 5\) (\(\mathrm{Pic} = \mathbb{Z}/6\)).

**\(\mathbb{Q}(\sqrt{-11})\)** (\(d_K = -11\), \(w_K = 2\), \(11\) ramified, \(2, 7\) inert, \(3, 5\) split; the
out-of-sample check). Levels \(2\alpha = 9, 13, 20, 24, 31, 35, \dots\); the trivial twist at
\(D = -15\) (\(2\alpha = 13\)); \(N_e(n) = h(-11n^2)\); \(M(2) = -2^{12}\) with \(h = 3\) — the first cubic
level is \(n = 2\): \(H_2\) of degree \(6\), \(L_3 = \mathbb{Q}[y]/(y^3 - y^2 - y - 1)\) (the tribonacci field),
\(h_{L_3} = 1\), \(\theta_u = \theta/16 = \pm\eta^8\), full index \(576 = 2\cdot24^2/2\); genus fields
\(\mathbb{Q}(\sqrt{33}), \mathbb{Q}(\sqrt5), \mathbb{Q}(\sqrt{77})\) at \(p = 3, 5, 7\).

**The index record** (all levels `bnfcertify = 1`, i.e. unconditional; "sat" is
\([\mathcal{O}^\times_{H_n}:\mathcal{W}_n]\), the index modulo \(24\)-th roots, equal to \(h_{H_n}\prod C\) in every row):

| \(K\) | \(n\) | \(\mathrm{Pic}(\mathcal{O}_n)\) | \(h\) | \(w_{H_n}\) | \(h_{H_n}\) | \(\prod C_\chi(0)\) | index | \(= \tfrac{w_K24^{h-1}}{w}h\prod C\) | sat | \(\mathcal{O}^\times/\mu\mathcal{V}_n\) |
|---|---|---|---|---|---|---|---|---|---|---|
| \(\mathbb{Q}(\sqrt{-2})\) | 2 | \(\mathbb{Z}/2\) | 2 | 8 | 1 | 1 | 6 | \(2\cdot24/8\) | 1 | \(\mathbb{Z}/6\) |
| | 3 | \(\mathbb{Z}/2\) | 2 | 6 | 1 | 1 | 8 | \(2\cdot24/6\) | 1 | \(\mathbb{Z}/8\) |
| | 4 | \(\mathbb{Z}/4\) | 4 | 8 | 1 | 3 (\(\chi_2\) from level 2, \(\ell = 2 \mid m\)) | 10368 | \(2\cdot24^3\cdot3/8\) | 3 | \(\mathbb{Z}/6\times\mathbb{Z}/24\times\mathbb{Z}/72\) |
| | 5 | \(\mathbb{Z}/6\) | 6 | 2 | 1 | 1 | 7962624 | \(24^5\) | 1 | \((\mathbb{Z}/24)^5\) |
| | 6 | \((\mathbb{Z}/2)^2\) | 4 | 24 | 1 | 24 (\(4\) from level 3: \(\ell = 2\) ramified, \(\chi(\mathfrak{p}) = -1\); \(6\) from level 2: \(\ell = 3\) split, \(\chi(\mathfrak{l}) = -1\)) | 27648 | \(2\cdot24^3\cdot24/24\) | 24 | \(\mathbb{Z}/2\times\mathbb{Z}/24\times\mathbb{Z}/576\) |
| \(\mathbb{Q}(\sqrt{-3})\) | 4 | \(\mathbb{Z}/2\) | 2 | 12 | 1 | 1 | 12 | \(6\cdot24/12\) | 1 | \(\mathbb{Z}/12\) |
| | 5 | \(\mathbb{Z}/2\) | 2 | 6 | 1 | 1 | 24 | \(6\cdot24/6\) | 1 | \(\mathbb{Z}/24\) |
| | 6 | \(\mathbb{Z}/3\) | 3 | 6 | 1 | 1 | 576 | \(6\cdot24^2/6\) | 1 | \((\mathbb{Z}/24)^2\) |
| | 7 | \(\mathbb{Z}/2\) | 2 | 6 | 1 | 1 | 24 | \(6\cdot24/6\) | 1 | \(\mathbb{Z}/24\) |
| | 8 | \((\mathbb{Z}/2)^2\) | 4 | 24 | 1 | 3 (from level 4, \(\ell = 2 \mid m\)) | 10368 | \(6\cdot24^3\cdot3/24\) | 3 | \(\mathbb{Z}/6\times\mathbb{Z}/24\times\mathbb{Z}/72\) |
| \(\mathbb{Q}(\sqrt{-7})\) | 3 | \(\mathbb{Z}/4\) | 4 | 6 | 1 | 1 | 4608 | \(2\cdot24^3/6\) | 1 | \(\mathbb{Z}/8\times(\mathbb{Z}/24)^2\) |
| | 4 | \(\mathbb{Z}/2\) | 2 | 4 | 1 | 1 | 12 | \(2\cdot24/4\) | 1 | \(\mathbb{Z}/12\) |
| | 5 | \(\mathbb{Z}/6\) | 6 | 2 | 1 | 1 | 7962624 | \(24^5\) | 1 | \((\mathbb{Z}/24)^5\) |
| | 6 | \(\mathbb{Z}/4\) | 4 | 6 | 1 | 45 (\(3, 5, 3\) from level 3: \(\ell = 2\) split, \(\chi(\mathfrak{l}) = \pm i, -1\)) | 207360 | \(2\cdot24^3\cdot45/6\) | 45 | \(\mathbb{Z}/8\times\mathbb{Z}/72\times\mathbb{Z}/360\) |
| \(\mathbb{Q}(\sqrt{-11})\) | 2 | \(\mathbb{Z}/3\) | 3 | 2 | 1 | 1 | 576 | \(2\cdot24^2/2\) | 1 | \((\mathbb{Z}/24)^2\) |
| | 3 | \(\mathbb{Z}/2\) | 2 | 6 | 1 | 1 | 8 | \(2\cdot24/6\) | 1 | \(\mathbb{Z}/8\) |
| | 4 | \(\mathbb{Z}/6\) | 6 | 4 | 1 | 9 (\(3\cdot3\): both cubic characters from level 2, \(\ell = 2 \mid m\)) | 35831808 | \(2\cdot24^5\cdot9/4\) | 9 | \(\mathbb{Z}/12\times(\mathbb{Z}/24)^2\times(\mathbb{Z}/72)^2\) |
| | 5 | \(\mathbb{Z}/4\) | 4 | 2 | 2 (\(\mathrm{Cl} = \mathbb{Z}/2\)) | 1 | 27648 | \(2\cdot24^3\cdot2/2\) | 2 | \((\mathbb{Z}/24)^2\times\mathbb{Z}/48\) |
| \(\mathbb{Q}(i)\) (anchor) | 5 | \(\mathbb{Z}/2\) | 2 | 4 | 1 | 1 | 24 | Paper II | 1 | \(\mathbb{Z}/24\) |

The level \(n = 5\) of \(\mathbb{Q}(\sqrt{-11})\) is the first with \(h_{H_n} > 1\) away from \(\mathbb{Q}(i)\): the
saturated quotient \(\mathcal{O}^\times/\mathcal{W}_n \cong \mathbb{Z}/2 \cong \mathrm{Cl}(H_5)\), as Conjecture 5.3 of the index
document predicts. Cubic layers: \(\mathbb{Q}(\sqrt{-7})\), \(n = 5\): \(L_3 = \mathbb{Q}[y]/(y^3 - 5y - 5)\), \(h_{L_3} = 1\),
\(\theta_u = \pm\eta^8\), sextic index \(576\); \(\mathbb{Q}(\sqrt{-11})\), \(n = 4\) (pullback of the level-2 cubic
character, \(C = 3\)): the same tribonacci \(L_3\), \(\theta_u = \pm\eta^{24} = \pm\eta^{8\cdot1\cdot3}\), sextic index
\(5184 = 576\cdot9\) — the tower theorem of the Euler document with multiplier \(\ell + 1 = 3\) at \(2 \to 4\).
Quadratic layers \(\theta^{(2)}_u = \pm\varepsilon_{d_2}^{e}\), \(|e| = 6L'(0,\chi)/\log\varepsilon_{d_2}\) exactly at every real
character (e.g. \(|e| = 3\) in \(\mathbb{Q}(\sqrt2)\) at \(\mathbb{Q}(\sqrt{-2})\), \(n = 2\); \(12\) in \(\mathbb{Q}(\sqrt5)\) at every
level-5 field; \(10\) in \(\mathbb{Q}(\sqrt{21})\) at \(\mathbb{Q}(\sqrt{-7})\), \(n = 6\)).

## 9. Class number two: \(\mathbb{Q}(\sqrt{-5})\) (experimental)

\(K = \mathbb{Q}(\sqrt{-5})\), \(d_K = -20\), \(\mathcal{O}_K = \mathbb{Z}[\sqrt{-5}]\), \(h_K = 2\), \(w_K = 2\), the non-principal
class represented by \(L = \mathfrak{p}_2 = (2, 1+\sqrt{-5})\), \(L^2 = (2)\). The Bianchi orbifold has two
cusps; the single orbit \(\mathcal{S}_K\) of \(\hat{\mathbb{R}}\) is the cusp of \(\mathcal{O}_K\).

**Stange's statement, verified.** The primitive index-\(n\) sublattices of \(\mathcal{O}_K\) — the
disks of curvature \(\sqrt{20}\,n\) of \(\mathcal{S}_K\) modulo translation — number \(N_e(n) = h(\mathcal{O}_n)/h_K\)
and represent exactly the kernel of \(\mathrm{Pic}(\mathcal{O}_n) \to \mathrm{Pic}(\mathcal{O}_K)\), each class once
(\(w_K/2 = 1\)), a subgroup of index \(2\): at \(n = 2\) (\(\mathrm{Pic} = \mathbb{Z}/4\)) and \(n = 3\)
(\((\mathbb{Z}/2)^2\)) the two kernel classes are closed under products (certified).

**The two-cusp dictionary.** The proper \(\mathcal{O}_n\)-lattices \(\Lambda\) with \(\mathcal{O}_K\Lambda \in \{\mathcal{O}_K, L\}\)
and \([\mathcal{O}_K\Lambda:\Lambda] = n\) represent all of \(\mathrm{Pic}(\mathcal{O}_n)\) (each class once), and
$$G_\mathfrak{c} := n^{12}\,\frac{\Delta(\Lambda)}{\Delta(\mathcal{O}_K\Lambda)}$$
is homogeneous of degree \(0\), so well defined on classes, with \(\sigma_\mathfrak{a}(G_\mathfrak{c}) = G_{\mathfrak{a}^{-1}\mathfrak{c}}\)
(the reference lattice moves with the class: \(\mathcal{O}_K(\mathfrak{a}^{-1}\Lambda) = \mathfrak{a}^{-1}\mathcal{O}_K\Lambda\)) and
\(\overline{G_\mathfrak{c}} = G_{\mathfrak{c}^{-1}}\); \(D_n \in \mathbb{Z}[x]\): \(D_2 = x^4 - 1056x^3 - 47744x^2 - 1263616x + 4096\),
\(D_3 = x^4 - 10116x^3 + 7304966x^2 - 3005316x + 1\) (certified).

**The corrected limit formula.** Since \(\log|G_\mathfrak{c}| = 6\log n + \log g(\Lambda_\mathfrak{c}) - \log g(\mathcal{O}_K\Lambda_\mathfrak{c})\),
$$\sum_\mathfrak{c}\chi(\mathfrak{c})\log|G_\mathfrak{c}| = -12L'(0,\chi) + 12\,|\ker\pi|\,L'(0,\chi_K)\cdot[\chi = \chi_K\circ\pi],$$
where \(\chi_K\) is the character of \(\mathrm{Pic}(\mathcal{O}_K)\) and \(L'(0,\chi_K)\) its Epstein value at level \(1\)
(\(w = 2\)); and \(L'(0,\chi_K\circ\pi) = C\cdot L'(0,\chi_K)\) with \(C = 4\) at \(n = 2\) (\(2\) ramified,
\(\chi_K(\mathfrak{p}_2) = -1\)) and \(C = 6\) at \(n = 3\) (\(3 = \mathfrak{p}\bar{\mathfrak{p}}\) with non-principal factors). Certified
at both levels (spare \(\ge 149\)). Consequently the character sum of the pulled-back
character is \(-12L'(0,\chi)(1 - |\ker\pi|/C)\), nonzero at \(n = 2, 3\) (\(C \ne |\ker\pi| = 2\)).

**The index.** With the Dedekind determinant and \(\zeta_K(0) = -h_K/w_K = -1\):
$$\bigl[\mathcal{O}_{H_n}^\times:\mu(H_n)\mathcal{V}_n\bigr] = \frac{24^{h-1}w_K}{h_K\,w_{H_n}}\,h_{H_n}
\prod_{\chi\ne1}C_\chi(0)\prod_{\chi = \chi_K\circ\pi}\Bigl|1 - \frac{|\ker\pi|}{C_\chi(0)}\Bigr| ,$$
certified with PARI (`bnfcertify = 1`): \(n = 2\): \(H_2\) of degree \(8\), \(h = 1\), \(w = 4\), index
\(6912 = 24^3\cdot\tfrac14\cdot4\cdot\tfrac12\); \(n = 3\): degree \(8\), \(h = 1\), \(w = 12\), index
\(4608 = 24^3\cdot\tfrac1{12}\cdot6\cdot\tfrac23\). The single-cusp data alone (the kernel classes) have
rank \(h/h_K - 1 = 1\) only; the full system has rank \(h - 1 = 3\). Everything in this section is
**experimental** in the sense of the program: the cusp bookkeeping (that the second cusp's
circles are the lattices \(\Lambda \subset L\)) is not proved here, and the normalization
\(\Delta(\mathcal{O}_K\Lambda)\) — forced by homogeneity — is what produces the base-character correction.

## 10. Uniform statements and their status

| statement | proved for every \(K\) | proved for the computed \(K\) | certified only |
|---|---|---|---|
| Thm 1: congruence description, lines, orientation, level set | necessity, orientation, level, \(2\alpha \equiv -2s\) | sufficiency (descent: the five Euclidean fields) | — |
| Thm 2: level-\(\alpha\) circles = forms of disc \(4(\alpha^2-1)/d_K\); census \(3H\) | yes (given Thm 1) | | |
| Thm 3: \(\hat\sigma[f] = [\mathfrak{r}_\alpha][f]^{-s}\), Lemmas A′–C′ | yes (all \(K\), all levels) | | imprimitive classes |
| Thm 4: unit theorem for \(R_f\), laws, KLF on odd characters | yes | | polynomials at \(2\alpha \le 18\) |
| Thm 5: \(N_e = \tfrac{w_K}{2}h\), \(D_n \in \mathbb{Z}[x]\), mass exponents, per-class law | yes for \(h_K = 1\) | | the sign law |
| Thm 6: KLF; genus field \(\mathbb{Q}(\sqrt p)\) / \(\mathbb{Q}(\sqrt{p|d_K|})\) at prime levels | yes | | |
| Thm 7: fiber lemma, norm relations, Hecke recursion | yes | | 26 + 27 + 25 + 21 pairs |
| Lemma 7.1: \(\mu(H_n)\) | | the five fields | |
| Thm 8: full index \(\tfrac{w_K24^{h-1}}{w_{H_n}}h_{H_n}\prod C\); layers | yes for \(h_K = 1\) | | saturation \(= h\prod C\) |
| §9: two-cusp dictionary, corrected KLF, index with \(w_K/h_K\) | | | experimental (\(n = 2, 3\)) |

**Where the expected statements of the prompt failed, and the corrected statements.**
(i) The level is not an integer for \(d_K \equiv 1 \bmod 4\): \(2\alpha \in \mathbb{Z}\) with
\(2\alpha \equiv \pm2 \bmod |d_K|\), and the sign fixes the orientation. (ii) The twist is not
\((\tfrac{\alpha-1}{2}, 0, \tfrac{\alpha+1}{2})\)-shaped in general: for odd \(d_K\) it is the ambiguous form
\((r_0, r_0, \cdot)\) with \(r_0 = (2\alpha-2)/\gcd(2\alpha-2,|d_K|)\), of norm \(2(\alpha-1)\) or \(2(\alpha-1)/|d_K|\),
and the class formula has the exponent \(-s\), so at the negatively oriented levels it is the pure
twist. (iii) The norm relation carries the global sign \((-1)^{[\ell=2]}\) also at inert and split
\(2\). (iv) \(\mu(H_n)\) at \(\mathbb{Q}(\sqrt{-2})\) is \(\mu_8\) already at even \(n\) (not \(\mu_4\)), because
\(K(i) = \mathbb{Q}(\zeta_8)\). (v) At \(\mathbb{Q}(\sqrt{-5})\) the limit formula holds uncorrected only off
the characters of the base class group.

## 11. Machine verification

`python3 scripts/other_fields.py --selftest` (needs `mpmath`, `sympy`, PARI/GP 2.15.4;
timings at the end of this section). Phases, with the guard rails of CLAUDE.md (precision
set in `main()`; integers accepted only with \(\ge\max(20,\mathrm{dps}/5)\) spare digits in
absolute error; exact HNF arithmetic with respect to \((1,\omega)\) for lattices, class groups,
projections, ideal twists and characters; PARI results labelled by `bnfcertify`):

- **A** (60 digits): orbit BFS versus congruences for \(|q| \le 5\); the level set and its signs.
- **B**: descent to an explicit \(X\) for every circle of the first ten levels; the class
  formula at every primitive class; Lemmas A′, B′, C′ at every class (\(\iota(\mathcal{K}_f) = \mathfrak{t}\mathfrak{a}_f\)
  as exact lattices in \(\mathbb{Q}(\sqrt D)\)).
- **C** (200 digits): the units \(R_f\) by the lattice formula with the exact twist lattice
  \(\mathfrak{r}^{-1}\mathfrak{b}\); integer palindromic level polynomials with constant term \(1\), the
  imprimitive strata, the twist and conjugation laws, \(\sum\chi\log|R| = -24L'\) against the
  independent evaluation.
- **D/E** (150 digits): \(D_n \in \mathbb{Z}[x]\) at \(2 \le n \le 9\); \(|D_n(0)|\) = the mass law; the sign law;
  Newton polygons (one slope \(w_p(k)\)); the limit formula at every character; the genus
  closed forms at prime levels with PARI's `qfbclassno`/`quadunit`.
- **F** (80 digits): the fiber lemma and the four-case norm relation at \(\ge 20\) pairs per
  field (\(26, 27, 25, 21\)); the Hecke recursion at two or three primes.
- **G** (150 digits; PARI 120 digits): the index at \(n = 2,3,4,5,6\) (\(\mathbb{Q}(\sqrt{-2})\)),
  \(4,5,6,7,8\) (\(\mathbb{Q}(\sqrt{-3})\)), \(3,4,5,6\) (\(\mathbb{Q}(\sqrt{-7})\)), \(2,3,4,5\) (\(\mathbb{Q}(\sqrt{-11})\)) and the
  Gaussian anchor \(n = 5\) (index \(24\), Paper II); V1–V6 as in [robert-index-full.md](robert-index-full.md);
  regression records `INDEX_RECORD`.
- **H** (150 digits): \(\mathbb{Q}(\sqrt{-5})\) at \(n = 2, 3\).

Runtimes of the selftest (one core, 212 s in total; PARI 2.15.4): geometry \(< 1\) s per field;
units \(4, 22, 4, 2\) s; Euclidean \(18, 10, 17, 22\) s; Euler \(22, 5, 23, 19\) s; index \(12, 8, 11, 9\) s
(and \(1\) s for the Gaussian anchor); \(\mathbb{Q}(\sqrt{-5})\) \(4\) s (fields in the order
\(-8, -3, -7, -11\)). Every `bnfcertify` returned \(1\) (all class fields have degree \(\le 12\)), so
the class numbers, unit groups and regulators of §7–§9 are **unconditional**. Observed spare
digits: \(\ge 190\) (units), \(\ge 102\) (\(D_n\)), \(\ge 148\) (limit formulas), \(\ge 77\) (norm relations),
\(\ge 105\) (indices). The regression record `INDEX_RECORD` in the script fixes \(h, w_{H_n}, \prod C\),
the index and the Smith invariants at all \(19\) index levels.

## 12. Literature diligence

Searches ran through the web-search relay; arXiv, journal and aggregator full texts
were not reachable from this environment (blocked by the egress proxy), so the
statements below are verified from the abstracts and the bibliographic records that the
relay returned, as in the papers' `NOTES.md`.

1. **K. E. Stange**, *Visualising the arithmetic of imaginary quadratic fields*, IMRN 2018,
   no. 12, 3908–3938 (arXiv:1410.0417). Abstract, verified: the orbit of \(\hat{\mathbb{R}}\) under
   \(\mathrm{PSL}_2(\mathcal{O}_K)\) is the Schmidt arrangement \(\mathcal{S}_K\); the curvatures of its circles
   are integer multiples of \(\sqrt{-\Delta}\); the curvatures of tangent circles are
   described by the norm form; the circles are in bijection with certain ideal classes in
   orders of \(K\), the conductor being a certain multiple of the curvature, which counts
   circles with class numbers; \(\mathcal{S}_K\) is connected iff \(\mathcal{O}_K\) is Euclidean. The
   bijection "circles of curvature \(f\sqrt{-\Delta}\) up to translation and rotation by
   \(180^\circ\) \(\leftrightarrow\) \(\ker(\mathrm{Pic}(\mathcal{O}_f)\to\mathrm{Pic}(\mathcal{O}_K))\)" is the form recorded
   in Paper I's `NOTES.md`; we could not read the theorem numbers. **Credited** for:
   curvature integrality (Theorem 1's curvature statement), the class-number census
   (Theorem 5(1) is this bijection made explicit through primitive sublattices, with the
   factor \(w_K/2\) and the orientation bookkeeping), and the single-cusp statement of §9
   (verified at \(\mathbb{Q}(\sqrt{-5})\), \(n = 2, 3\)).
2. **K. E. Stange**, *The Apollonian structure of Bianchi groups*, Trans. AMS 370 (2018)
   6169–6219 (arXiv:1505.03121). Abstract, verified: \(K\)-Bianchi circles, the Schmidt
   arrangement as a geometric realization of the arithmetic of \(K\), the decomposition
   (with orientations) into primitive integral \(K\)-Apollonian packings, the \(K\)-Apollonian
   groups, a local-to-global conjecture on curvatures; for \(K \ne \mathbb{Q}(\sqrt{-3})\) circles meet
   only tangentially. **Credited** for the definition of \(\mathcal{S}_K\) and the oriented viewpoint.
3. **J. Rickards, K. E. Stange**, *Eisenstein circle packings and the Eisenpint Schmidt
   arrangement*, arXiv:2605.16053 (May 2026). Abstract, verified: at \(\mathbb{Q}(\sqrt{-3})\) circles
   meet at angles \(\pi/3, 2\pi/3\); a modified arrangement (the Eisenpint arrangement) and
   Eisenstein circle packings, strong approximation, congruence obstructions, density-one
   local-global, quadratic reciprocity obstructions. **Adjacent**: it concerns packings and
   curvatures at \(\mathbb{Q}(\sqrt{-3})\), not the level stratification, the involution, the units or
   the index; the three-direction line structure of §1 is consistent with its setting.
4. **D. Cox**, *Primes of the form \(x^2+ny^2\)*, 2nd ed. (2013), §7 (orders, \(\mathrm{Pic}(\mathcal{O}_n)\)),
   §9 (ring class fields), Thm 3.15 (genus theory); **Neukirch** Cor. VII.5.11 (class number
   formula at \(s = 0\)); **Washington** Lemmas 4.15, 5.26; **de Shalit** Ch. II §2 and **Rubin**
   (1991) for the classical norm-compatibility of elliptic units; **Robert**, **Kubert–Lang**,
   **Schertz** for the index-equals-class-number shape — as positioned in Paper II's notes.
   None of these treats the Schmidt-geometric dictionary at a general \(K\), the explicit
   twist ideal \(\mathfrak{r}_\alpha\), the constants \(w_K/w_{H_n}\), or the index formulas.

**Verdict.** New here: the level, its half-integrality and orientation law; the explicit
congruence classification with its descent proof for the five Euclidean fields (in
Stange's framework); the discriminant \(4(\alpha^2-1)/d_K\); the class formula with the twist
ideal \(\mathfrak{r}_\alpha\) and the orientation exponent; the unit theorem and the horizontal
coincidences; the mass law with the constant \(24/(w_Ke_p)\); the genus-field law; the norm
relations with \(\pi^{12}\) and the global sign; the torsion lemma and the index formula with
\(w_K\); the two-cusp correction. Not new: Stange's curvature integrality and class-number
bijection; the classical elliptic-unit theory of ring class fields.

## 13. Outlook

- **Non-Euclidean class number one** (\(d_K = -19, -43, -67, -163\)): necessity, the level, the
  dictionary and everything on the Euclidean side hold; sufficiency of the congruence
  classification (Theorem 1(2)) needs the class number of the indefinite unimodular
  Hermitian lattice \((\mathcal{O}_K^2, h_0)\) in place of the descent. The unit and index theorems
  do not depend on it.
- **Class number \(> 1\)**: prove the cusp bookkeeping of §9 (the circles of the cusp of
  \(L\) are the lattices \(\Lambda \subset L\) with \(\mathcal{O}_K\Lambda = L\), realized by
  \(\mathrm{SL}(L\oplus\mathcal{O}_K)\)); decide whether a normalization of the reference lattice
  removes the base-character correction (a Galois-equivariant choice of
  \(\Delta(L)^{1/h_K}\)-type data), and compute \(\mathbb{Q}(\sqrt{-6})\), \(\mathbb{Q}(\sqrt{-15})\).
- **Item 4 (horizontal families)**: the coincidences of §3 show that the pair \((D,[\mathfrak{r}])\)
  is the invariant; enumerate which pairs each \(K\) realizes (\(D = -(k-2)(k+2)/|d_K|\) with the
  twist of norm \(r_0\)) and prove the unit theorem for an arbitrary invertible ambiguous
  twist with principal square directly — the proof of Theorem 4 already does.
- **The phase** \(u_f = \Phi_y/\Phi_x\) was not transported: the closed form, the Norm Lemma
  and \(\omega_f = 1\) over \(K\) would identify \(R_f\) with the twisted \(|u_f|^6\) as in Paper II
  Lemma 4.1; the correspondence is \(X_0(t_0)\times_{X(1)}X_0(\cdot)\) with the two norms
  \(r_0, s_0\) in place of \(\tfrac{n\mp1}{2}\).
- **The sign law of \(M(n)\)** (Theorem 5(5)): prove it by Paper I's M4 with the sign of
  \(\Delta(\mathcal{O}_K)\) and the conjugation-stable representatives \((1-\zeta)\Lambda\) at \(\mathbb{Q}(\sqrt{-3})\).
- **Saturation**: the Gras-type Conjecture 5.3 of the index document holds at every level
  here; the \(\mathbb{Q}(\sqrt{-3})\) levels with \(w_{H_n} = 12, 24\) are the place to study the
  \(3\)-part of the \(24\)-th-root defect.
