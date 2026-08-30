# The Kronecker limit formula for the phase: character sums of \(\log|u|\) and \(L'(0,\chi)\)

This document answers the question of [outlook.md](outlook.md) §2.3 and
[euclidean-moduli-invariants.md](euclidean-moduli-invariants.md) §6.6: **do the
character sums of the log-absolute-values of the phase invariants evaluate to
explicit combinations of derivatives of \(L\)-functions at \(s = 0\)?**

The answer is **yes, as a theorem, in both aspects**. For every character
\(\chi\) of the relevant class group,
$$
\text{Euclidean:}\quad S(\chi) \;=\; -2\,L'(0,\chi) \;+\; \tfrac23\,\Sigma_0(\chi) \;+\; \tfrac12\,\Sigma_{1728}(\chi),
\qquad
\text{hyperbolic:}\quad S(\chi) \;=\;
\begin{cases}
0, & \chi(\mathfrak{r}_n) = +1,\\[2pt]
-4\,L'(0,\chi) + \tfrac43\,\Sigma_0(\chi) + \Sigma_{1728}(\chi), & \chi(\mathfrak{r}_n) = -1,
\end{cases}
$$
where \(S(\chi) = \sum \chi\,\log|u|\), \(\Sigma_x(\chi) = \sum \chi\,\log|j - x|\)
over the classes of the level, and \(L(s,\chi)\) is the form-class (Epstein)
\(L\)-function of the discriminant of the level. The \(L'\)-part is the
Kronecker limit formula; the \(j\)-parts are Gross–Zagier collision quantities,
identified **exactly** below (real quadratic \(S\)-numbers for genus characters,
with certified unit powers and split-prime data). Everything displayed is
machine-verified by [scripts/phase_klf.py](scripts/phase_klf.py) at 250 digits
(400 at Euclidean \(n = 15\)) under the certification guard rails of
[CLAUDE.md](CLAUDE.md); the sections state explicitly what is proved, what is
certified, and what failed to fit (§8).

Two structural surprises. First, on the hyperbolic side the expected new
ingredient \(\log\varepsilon_{n^2-1}\) (from \(\varepsilon^6\mu^{-12}\)) **cancels
identically** — the Norm Lemma kills both \(\varepsilon\) and \(\mu\) in
\(|u_f|\), leaving a scale-invariant two-CM-point law (Lemma 5.1). Second, on
the Euclidean side the real quadratic field carrying the genus-character sum at
conductor \(n\) is **\(\mathbb{Q}(\sqrt n)\)** itself: the phases of the
curvature-\(2n\) disks know the fundamental unit of \(\mathbb{Q}(\sqrt n)\) (§4).

Throughout: \(K = \mathbb{Q}(i)\); Euclidean level \(n\) means
\(\mathcal{O}_n = \mathbb{Z} + n\mathbb{Z}[i]\), disc \(D = -4n^2\), phases
\(u_\mathfrak{c} = \Theta/\Omega\) with \(|u|\) well-defined on
\(\mathrm{Cl}(\mathcal{O}_n)\) ([euclidean-moduli-invariants.md](euclidean-moduli-invariants.md) §5);
hyperbolic level \(n\) (odd) means disc \(D = 1-n^2\), phases
\(u_f = \varepsilon\Theta_f\) on primitive classes
([moduli-invariants.md](moduli-invariants.md) §4–5), twist class
\(\mathfrak{r}_n = [(\tfrac{n-1}2, 0, \tfrac{n+1}2)]\). In both cases the unit
count of the order is \(w = 2\). We write \(\Delta_q = \eta^{24}\),
\(\tau_\mathfrak{c}\) for the CM point of the reduced form
\(Q_\mathfrak{c} = (A, B, C)\) of a class, \(y_\mathfrak{c} = \operatorname{Im}\tau_\mathfrak{c}
= \sqrt{|D|}/2A\), and
$$
g_\mathfrak{c} \;:=\; y_\mathfrak{c}^{\,6}\,\bigl|\Delta_q(\tau_\mathfrak{c})\bigr|
$$
— the scale-invariant Delta-datum of the class (\(g = \mathrm{covol}^6|\Delta|\)
of any lattice in the class, up to one global constant).

## 1. The \(L\)-function of a level and the limit formula

**Definition.** For a discriminant \(D < 0\) with class group \(\mathrm{Cl}(D)\)
(primitive forms) and a character \(\chi\), set
$$
L(s,\chi) \;:=\; \frac1w \sum_{\mathfrak{c}} \chi(\mathfrak{c})\, \zeta_{Q_\mathfrak{c}}(s),
\qquad
\zeta_Q(s) = \sideset{}{'}\sum_{(x,y)\in\mathbb{Z}^2} Q(x,y)^{-s} .
$$
For \(m \geq 1\) let \(r_\mathfrak{c}(m)\) be the number of representations of
\(m\) by \(Q_\mathfrak{c}\) and \(R_\chi(m) = \sum_\mathfrak{c}\chi(\mathfrak{c})r_\mathfrak{c}(m)\).
Since \((x,y)\mapsto(x,-y)\) bijects representations by \(Q\) and \(Q^{-1}\),
\(R_\chi = R_{\bar\chi}\) and \(L(s,\chi) = L(s,\bar\chi)\) is real.

**Proposition 1.1 (functional equation and the limit value).** Write
\(\gamma(s) = \bigl(\tfrac{\sqrt{|D|}}{2\pi}\bigr)^{s}\Gamma(s)\). Then
\(\gamma(s)\,\zeta_{Q}(s) = \gamma(1-s)\,\zeta_{Q}(1-s)\) for every class, so for
**nontrivial** \(\chi\) the function \(L(s,\chi)\) is entire with
\(L(0,\chi) = 0\), and
$$
L'(0,\chi) \;=\; \frac{\sqrt{|D|}}{2\pi}\,L(1,\chi)
\;=\; -\,\frac{1}{12}\sum_{\mathfrak{c}}\chi(\mathfrak{c})\,\log g_\mathfrak{c}.
$$

*Proof.* With \(\tau_Q = \tfrac{-B + i\sqrt{|D|}}{2A}\) one has
\(A|x + y\tau_Q|^2 = Q(x,-y)\), whence
\(\zeta_Q(s) = \bigl(\tfrac{2}{\sqrt{|D|}}\bigr)^{s} E(\tau_Q, s)\) for the
nonholomorphic Eisenstein series \(E(\tau,s) = \sum' y^s/|m\tau+n|^{2s}\)
(using \(A\,y_{\tau} = \sqrt{|D|}/2\)). The completed series
\(E^*(\tau,s) = \pi^{-s}\Gamma(s)E(\tau,s)\) satisfies \(E^* (\tau,s) = E^*(\tau,1-s)\),
and \(\gamma(s)\zeta_Q(s) = E^*(\tau_Q,s)\) — the functional equation, class by
class. For nontrivial \(\chi\) the poles cancel in the character sum, so
\(\Gamma(s)L(s,\chi)\) is finite at \(s = 0\); since \(\Gamma(s)\sim 1/s\),
\(L(0,\chi) = 0\) and \(L'(0,\chi) = \lim_{s\to0}\gamma(s)L(s,\chi)
= \gamma(1)L(1,\chi)\). The first Kronecker limit formula
\(E(\tau, s) = \tfrac{\pi}{s-1} + 2\pi\bigl(\gamma_E - \log 2 - \log(\sqrt{y}\,|\eta(\tau)|^2)\bigr) + O(s-1)\)
gives, after the \(\chi\)-sum kills all constants and with
\(\sqrt{y}|\eta|^2 = (y^6|\Delta_q|)^{1/12} = g^{1/12}\),
\(L(1,\chi) = -\tfrac1w\cdot\tfrac{2}{\sqrt{|D|}}\cdot\tfrac{\pi}{6}\sum\chi\log g\);
multiply by \(\gamma(1) = \sqrt{|D|}/2\pi\) and use \(w = 2\). \(\square\)

**Independent evaluation (used in every check below).** Unfolding
\(\gamma(s)\zeta_Q(s)\) against the theta series and applying
\(\theta_Q(1/t) = t\,\theta_{Q^{-1}}(t)\) yields, for nontrivial \(\chi\), the
rapidly convergent
$$
L'(0,\chi) \;=\; \frac1w \sum_{m\geq1} R_\chi(m)\left[\frac{e^{-\alpha m}}{\alpha m} + E_1(\alpha m)\right],
\qquad \alpha = \frac{2\pi}{\sqrt{|D|}},
$$
(\(E_1\) the exponential integral). This computes \(L'(0,\chi)\) from the exact
integer representation numbers alone — independent of every modular quantity —
and is the anchor against which all limit-formula identities are verified
(residuals \(\le 10^{-247}\) at 250 digits).

## 2. Euclidean: the master identity

**Lemma 2.1 (classwise law).** For every class \(\mathfrak{c}\) of every level
\(n\ge2\), with \(\beta = j(\mathfrak{c})\),
$$
|u_\mathfrak{c}|^6 \;=\; \bigl|\beta^4(\beta - 1728)^3\bigr|\;
\frac{g_\mathfrak{c}}{n^6\,|\Delta_q(i)|}.
$$

*Proof.* \(u^6 = -\beta^4(\beta-1728)^3\,\Delta(\Lambda_\mathfrak{c})/\Delta(\mathbb{Z}[i])\)
(proved, [euclidean-moduli-invariants.md](euclidean-moduli-invariants.md) §5.1).
Writing \(\Lambda = c\,[1, w]\), \(|c|^2 y_w = n = \mathrm{covol}(\Lambda)\) and
\(|\Delta(\Lambda)| = |c|^{-12}|\Delta_q(w)|(2\pi)^{12}\), so
\(|\Delta(\Lambda)/\Delta(\mathbb{Z}[i])| = y_w^6|\Delta_q(w)|/(n^6|\Delta_q(i)|)\),
and \(y_w^6|\Delta_q(w)| = g_\mathfrak{c}\) by \(\mathrm{SL}_2(\mathbb{Z})\)-invariance. \(\square\)

(Verified classwise, constant included, to \(10^{-247}\) at every computed level.)

> **Theorem 1 (Euclidean master identity).** For every nontrivial character
> \(\chi\) of \(\mathrm{Cl}(\mathcal{O}_n)\), with
> \(S(\chi) = \sum_\mathfrak{c}\chi(\mathfrak{c})\log|u_\mathfrak{c}|\) and
> \(\Sigma_x(\chi) = \sum_\mathfrak{c}\chi(\mathfrak{c})\log|j(\mathfrak{c}) - x|\):
> $$
> \boxed{\;S(\chi) \;=\; -2\,L'(0,\chi) \;+\; \tfrac{2}{3}\,\Sigma_0(\chi) \;+\; \tfrac{1}{2}\,\Sigma_{1728}(\chi).\;}
> $$
> For the trivial character (the mandatory anchor; a corollary of Theorem 2 and
> the \(\Delta\)-mass law of
> [euclidean-moduli-invariants.md](euclidean-moduli-invariants.md) Thm 4):
> $$
> 6\sum_\mathfrak{c}\log|u_\mathfrak{c}| \;=\; 4\log|H_{-4n^2}(0)| + 3\log|H_{-4n^2}(1728)| + \log|M(n)| - 12\,h\log n .
> $$

*Proof.* Take \(\log\) of Lemma 2.1, sum against \(\chi\) (the
\(\mathfrak{c}\)-independent constant \(-6\log n - \log|\Delta_q(i)|\) dies for
\(\chi \neq 1\)), and substitute
\(\sum\chi\log g = -12L'(0,\chi)\) from Proposition 1.1. For the trivial
character, sum the same identity over all classes:
\(\sum\log|\beta| = \log|H(0)|\), \(\sum\log|\beta-1728| = \log|H(1728)|\), and
\(\sum_\mathfrak{c}(\log g_\mathfrak{c} - \log|\Delta_q(i)| - 6\log n)
= \log|M(n)| - 12h\log n\) by definition of the mass
\(M(n) = \prod_\mathfrak{c} n^{12}\Delta(\Lambda_\mathfrak{c})/\Delta(\mathbb{Z}[i])\),
whose closed form is Theorem 4 there. \(\blacksquare\)

Both statements verified at \(n = 3, 5, 7, 9, 11, 13, 15\), every character,
with residuals \(\le 10^{-247}\) (250 digits; the identity itself and the
independent \(L'\)-evaluation of §1 are separate computations).

## 3. Euclidean: the elliptic-unit uniformization

The limit formula becomes an algebraic statement through one certified integer
polynomial per level.

> **Theorem 2 (the \(\Delta\)-mass polynomial).** Set
> \(G_\mathfrak{c} := n^{12}\,\Delta(\Lambda_\mathfrak{c})/\Delta(\mathbb{Z}[i])\)
> (an algebraic integer, nonzero, of absolute value \(n^6 g_\mathfrak{c}/|\Delta_q(i)|\)).
> Then:
> 1. \(\sigma\bigl(G_\mathfrak{c}\bigr) = G_{\mathfrak{c}(\sigma)^{-1}\mathfrak{c}}\)
>    for \(\sigma \in \mathrm{Gal}(\bar{\mathbb{Q}}/K)\), and
>    \(\overline{G_\mathfrak{c}} = G_{\mathfrak{c}^{-1}}\): the multiset
>    \(\{G_\mathfrak{c}\}\) is Galois-stable, so
>    \(D_n(x) := \prod_\mathfrak{c}(x - G_\mathfrak{c}) \in \mathbb{Z}[x]\),
>    with \(|D_n(0)| = |M(n)|\) given by the \(\Delta\)-mass law.
> 2. For **every** nontrivial \(\chi\):
>    \(\;-12\,L'(0,\chi) = \sum_\mathfrak{c}\chi(\mathfrak{c})\log|G_\mathfrak{c}|\).
>
> So all the \(L'(0,\chi)\) of a level are character combinations of the logs
> of the roots of one integer polynomial — the phases' \(\Delta\)-data
> \(\{G_\mathfrak{c}\}\) form an **elliptic-unit system for the ring class
> tower of \(\mathbb{Q}(i)\)** in the classical (Siegel–Robert) sense, realized
> geometrically by the curvature-\(2n\) disks.

*Proof.* (1) \(G_\mathfrak{c} = F(\Lambda_\mathfrak{c}; \mathbb{Z}[i])\) with
\(F(\Lambda; L) = n^{12}\Delta(\Lambda)/\Delta(L)\), isobaric of **joint weight
zero**; Lemma T of
[euclidean-moduli-invariants.md](euclidean-moduli-invariants.md) §5.6 gives
\(\sigma(F(\Lambda;\mathbb{Z}[i])) = F(s^{-1}\Lambda; s^{-1}\mathbb{Z}[i])\), and
rescaling by a generator of \((s^{-1}\mathbb{Z}[i])^{-1}\) (exactly as in the
translation law there, but now with **no ambiguity at all**, the weight being
zero rather than four) lands on the primitive representative of the class
\(\mathfrak{c}(\sigma)^{-1}\mathfrak{c}\). Mirror: \(\Delta(\bar\Lambda) = \overline{\Delta(\Lambda)}\)
and \(\bar\Lambda\) represents \(\mathfrak{c}^{-1}\). Rationality and (by
integrality of \(n^{12}\Delta\)-quotients, Lang, *Elliptic Functions* ch. 12)
integrality of the coefficients follow; \(D_n(0) = \pm\prod G_\mathfrak{c} = \pm M(n)\).
(2) \(\log|G_\mathfrak{c}| = \log g_\mathfrak{c} + 6\log n - \log|\Delta_q(i)|\),
so the character sum reproduces Proposition 1.1. \(\blacksquare\)

Certified polynomials (integer coefficients, \(\ge 190\) spare digits at 250;
constant terms match \(\pm M(n)\) exactly):

| \(n\) | \(D_n(x)\) |
|---|---|
| 3 | \(x^2 - 378x + 729\) |
| 5 | \(x^2 - 322x + 1\) |
| 7 | \(x^4 + 58604x^3 + 5618502582x^2 - 191605917748x + 117649\) |
| 9 | \(x^6 + 777114x^5 + 469570343535x^4 + 88033943990172204x^3 + 4975970841324612804303x^2 - 6464018128844293554150x + 3^{24}\) |
| 11 | \(x^6 - 2940366x^5 + 2977596073815x^4 + 116499046538184860x^3 + 38813935158542440699935x^2 - 1046289507684469650030x + 11^{6}\) |
| 13 | \(x^6 + 7210938x^5 + 12700046391567x^4 - 1529217345755339156x^3 + 61048319821249936271631x^2 - 22813565655771828294x + 1\) |
| 15 | \(x^8 - 243432x^7 - 49497426055620x^6 - \cdots + 3^{24}\) |

(Full coefficient lists are printed by the script; at \(n = 9, 11, 13\) the
cubic **coset products** of the \(G\)'s along the order-3 character —
\(x^3 - 76822296291x^2 + 4976085199686658256835x - 3^{24}\) at \(n=9\),
\(x^3 + 277946532501x^2 + 38813938298955894507411x - 11^6\) at \(n=11\),
\(x^3 + 446643445245x^2 + 61048319249786206560771x - 1\) at \(n=13\) — are
also individually certified, giving the Stark-type closed form
\(L'(0,\chi_3) = -\tfrac1{12}\bigl[\tfrac32\log|G_0^{(3)}| - \tfrac12\log|C_\Delta(0)|\bigr]\)
for the cubic characters, verified to \(10^{-247}\).)

## 4. Euclidean: genus characters in closed form

\(\mathrm{Cl}(\mathcal{O}_n)\) is cyclic of order 6 at \(n = 9, 11, 13\) (order
2 at \(n = 3, 5\); \(\mathbb{Z}/4\) at \(n = 7\); \(\mathbb{Z}/4\times\mathbb{Z}/2\)
at \(n = 15\)), so there is one nontrivial real character \(\chi_2\) per level
(three at \(n=15\)). For it, everything in Theorem 1 evaluates in closed form.

**(a) The \(L'\)-part.** The exact representation-number identity
$$
R_{\chi_2}(m) \;=\; 2\,\bigl(\mathrm{conv}_{d_1,d_2} * \,\mathrm{corr}\bigr)(m),
\qquad
\mathrm{conv}_{d_1,d_2}(m) = \sum_{e \mid m}\Bigl(\tfrac{d_1}{e}\Bigr)\Bigl(\tfrac{d_2}{m/e}\Bigr)
$$
holds — verified as an identity of integers for all \(m \le 300\) at every
level — for the decompositions and finite conductor Euler corrections

| \(n\) | \(d_1\cdot d_2 = -4n^2\) | correction | \(L'(0,\chi_2)\) |
|---|---|---|---|
| 3 | \((-3)\cdot 12\) | — | \(\tfrac13\log\varepsilon_{12}\) |
| 5 | \((-20)\cdot 5\) | — | \(2\log\varepsilon_{5}\) |
| 7 | \((-7)\cdot 28\) | — | \(\log\varepsilon_{28}\) |
| 9 | \((-27)\cdot 12\) | \(1 + 3^{1-2s}\) | \(\tfrac43\log\varepsilon_{12}\) |
| 11 | \((-11)\cdot 44\) | — | \(\log\varepsilon_{44}\) |
| 13 | \((-52)\cdot 13\) | — | \(2\log\varepsilon_{13}\) |
| 15 | \((-15)\cdot 60\) | — | \(4\log\varepsilon_{60}\) |

with \(\varepsilon_{12} = 2+\sqrt3\), \(\varepsilon_5 = \tfrac{1+\sqrt5}2\),
\(\varepsilon_{28} = 8+3\sqrt7\), \(\varepsilon_{44} = 10+3\sqrt{11}\),
\(\varepsilon_{13} = \tfrac{3+\sqrt{13}}2\), \(\varepsilon_{60} = 4+\sqrt{15}\).
(At \(n = 15\), \(\mathrm{Cl}(\mathcal{O}_{15}) = \mathbb{Z}/4\times\mathbb{Z}/2\)
has three real characters; the one factoring as displayed is again the
\(\mathbb{Q}(\sqrt{15})\)-character. The other two admit **no** finite-Euler
factorization over the divisor pairs of \(-900\) — an honest negative result of
the discovery procedure — but their values are cleanly identified (certified
ratio, 398 spare digits): \(L'(0,\chi) = 8\log\varepsilon_5\) and
\(\tfrac83\log\varepsilon_{12}\), the units of the two proper subfield levels
\(5\) and \(3\); their factorizations presumably need Euler data at both
conductor primes simultaneously.) Consequently
\(L(s,\chi_2) = L(s,(\tfrac{d_1}{\cdot}))L(s,(\tfrac{d_2}{\cdot}))\cdot C(s)\)
as Dirichlet series (Kronecker symbols, possibly imprimitive), and since the
positive part always has \(L(0) = 0\),
$$
L'(0,\chi_2) \;=\; \frac{2h(d_1)}{w(d_1)}\; h(d_2)\,\log\varepsilon_{d_2}\; \cdot\, C(0).
$$
**The real quadratic field is \(\mathbb{Q}(\sqrt n)\) in every case** —
\(d_2\) is the discriminant of \(\mathbb{Q}(\sqrt n)\) — and the closed forms
match the independent \(L'\)-computation of §1 to \(10^{-249}\) or better.

*Proof.* The coefficient identity is verified exactly for all \(m\le300\) at
every level. This is in fact a **proof**, by the Sturm bound: the left side is
the \(q\)-expansion of the genus-character theta combination
\(\tfrac12\sum_\mathfrak{c}\chi(\mathfrak{c})\theta_{Q_\mathfrak{c}}\), a
holomorphic weight-1 modular form of level \(|D| = 4n^2\) and character
\((\tfrac{D}{\cdot})\) (Hecke–Schoeneberg, valid for arbitrary
discriminants); the right side is the weight-1 Eisenstein series
\(E_1(\chi_{d_1}, \chi_{d_2})\) plus its \(p^a\)-dilations weighted by the
correction coefficients, of level dividing \(|D|\) as well (e.g. at \(n = 9\):
\(E + 3E(9\tau)\) with \(E\) of level 36). A weight-1 form with odd character
that vanishes beyond the Sturm bound \(\mu(\Gamma_0(|D|))/12\) — here at most
\(180\) (at \(n=15\)), far below the certified 300 — vanishes identically, and
no constant discrepancy can survive at weight 1. \(\square\)

(The conductor correction at \(n = 9\) — square conductor, and only there —
is thereby also proved as an identity; what remains open is only its
*conceptual* derivation from the local proper-\(\mathcal{O}_9\)-ideal
structure at 3, and the general-\(n\) law. The composite square-free conductor
\(n = 15\) needs **no** correction for its \(\mathbb{Q}(\sqrt{15})\)-character,
exactly as this reading predicts — verified at 400 digits.)

**(b) The \(j\)-dressing.** For the real character, split the classes by
\(\chi_2 = \pm1\) and let \(A, B\) be the products of \(j(\mathfrak{c})\)
(resp. of \(j(\mathfrak{c}) - 1728\)) over the two cosets. Then \(A + B\) and
\(AB\) are **certified integers**, \(A, B = \tfrac{s \pm t\sqrt{d_2}}2\) are
conjugate quadratic integers in the same field \(\mathbb{Q}(\sqrt n\,)\), and
the norm-one ratio \(A/B\) factors **exactly** (verified in
\(\mathbb{Q}(\sqrt{d_2})\)-arithmetic, valuations computed 2-adically/Hensel):
$$
\Sigma_0(\chi_2) = \log\Bigl|\frac{A}{B}\Bigr|
= k\,\log\varepsilon_{d_2} + \sum_p e_p\,\lambda_p,
\qquad
\lambda_p := \log\Bigl|\frac{\pi_p}{\pi_p'}\Bigr|,
$$
over **split** primes \(p\) of \(\mathbb{Q}(\sqrt{d_2})\), \(\pi_p\) an explicit
generator of \(\mathfrak{p}_p^{m_p}\) (\(m_p\) the class-order of the prime;
\(e_p \in \tfrac1{m_p}\mathbb{Z}\)), canonically normalized so that
\(0 \le \lambda_p < 2\log\varepsilon_{d_2}\) — the generators are printed by
the script and each identity is re-multiplied out exactly. The split-prime
support lies in \(\mathrm{GZ}(-4n^2,-3)\cup\{p \mid 2n\}\) for
\(\Sigma_0\) and \(\mathrm{GZ}(-4n^2,-4)\cup\{p\mid 2n\}\) for \(\Sigma_{1728}\)
— the same Gross–Zagier sets that govern the denominators in
[moduli-invariants.md](moduli-invariants.md) §5.7. The certified factorizations
(coefficients in the canonical normalization):

| \(n\) | \(\Sigma_0(\chi_2)\) | \(\Sigma_{1728}(\chi_2)\) |
|---|---|---|
| 3 | \(16\log\varepsilon_{12} - 3\lambda_{11} - 3\lambda_{23}\) | \(4\log\varepsilon_{12} + 2\lambda_{11}\) |
| 5 | \(30\log\varepsilon_{5} - 3\lambda_{11} + 3\lambda_{59} + 3\lambda_{71}\) | \(28\log\varepsilon_{5} + 2\lambda_{11} + 2\lambda_{19}\) |
| 7 | \(24\log\varepsilon_{28} - 3\lambda_{47} - 3\lambda_{83} - 3\lambda_{131}\) | \(18\log\varepsilon_{28} + 2\lambda_{3} - 2\lambda_{19}\) |
| 9 | \(16\log\varepsilon_{12} + 3\lambda_{11} - 3\lambda_{23} + 3\lambda_{47} + 3\lambda_{179} - 3\lambda_{227} - 3\lambda_{239}\) | \(24\log\varepsilon_{12} - 2\lambda_{11} - 4\lambda_{23}\) |
| 11 | \(12\log\varepsilon_{44} - 3\lambda_{107} + 3\lambda_{167} + 3\lambda_{263} + 3\lambda_{347} - 3\lambda_{359}\) | \(4\log\varepsilon_{44} + 6\lambda_{19} + 2\lambda_{43}\) |
| 13 | \(18\log\varepsilon_{13} - 3\lambda_{23} + 3\lambda_{107} + 3\lambda_{251} - 3\lambda_{311} + 3\lambda_{443} - 3\lambda_{491} + 3\lambda_{503}\) | \(28\log\varepsilon_{13} - 2\lambda_{3} - 4\lambda_{23} - 2\lambda_{43}\) |
| 15 \((2,0)\) | \(42\log\varepsilon_{60} - 3\lambda_{11} + 3\lambda_{59} + 3\lambda_{71} + 3\lambda_{191} + 3\lambda_{419} + 3\lambda_{479} + 3\lambda_{659}\) | \(52\log\varepsilon_{60} + 2\lambda_{11} + 2\lambda_{43} + 2\lambda_{59}\) |
| 15 \((0,1)\) | \(132\log\varepsilon_{5} + 3\lambda_{11} + 3\lambda_{59} + 3\lambda_{71} + 3\lambda_{191} - 3\lambda_{419} - 3\lambda_{479} - 3\lambda_{659}\) | \(124\log\varepsilon_{5} - 2\lambda_{11} - 4\lambda_{19} + 4\lambda_{31} - 2\lambda_{59}\) |
| 15 \((2,1)\) | \(62\log\varepsilon_{12} - 9\lambda_{11} - 6\lambda_{47} - 3\lambda_{59} - 3\lambda_{71} - 3\lambda_{191} + 3\lambda_{419} + 3\lambda_{479} - 3\lambda_{659}\) | \(16\log\varepsilon_{12} + 6\lambda_{11} + 4\lambda_{23} + 4\lambda_{47} + 2\lambda_{59}\) |

(the \(\pm3\)-multiples in \(\Sigma_0\) and even coefficients in
\(\Sigma_{1728}\) are the \((4,3)\)-exponent shadow of
\(u^6 = -\beta^4(\beta-1728)^3\cdot\Delta\)-quotient, i.e. of §5.7 there; note
that at \(n = 15\) all three real characters get exact \(\Sigma\)-data even
though two of them lack the finite genus factorization for \(L'\)).

> **Theorem 3 (closed form for the genus character).** For the real character
> \(\chi_2\) of \(\mathrm{Cl}(\mathcal{O}_n)\), \(n \le 13\) (and the
> \(\mathbb{Q}(\sqrt{15})\)-character at \(n = 15\); the other two real
> characters there carry the same \(\Sigma\)-closed forms with the certified
> \(\varepsilon\)-ratio standing in for the \(L'\)-factorization),
> $$
> S(\chi_2) \;=\; -2\,\frac{2h(d_1)}{w(d_1)}h(d_2)C(0)\,\log\varepsilon_{d_2}
> \;+\; \tfrac23\Sigma_0(\chi_2) + \tfrac12\Sigma_{1728}(\chi_2)
> $$
> with the tabulated exact \(\Sigma\)-factorizations: an explicit element of
> \(\mathbb{Q}\,\log\varepsilon_{\mathbb{Q}(\sqrt n)} + \sum_p \mathbb{Q}\,\lambda_p\)
> over the Gross–Zagier split primes. Verified to \(10^{-247}\); every
> ingredient exact. The phases are thus an elliptic-unit system **up to the
> identified \(\mathrm{GZ}\)-supported \(S\)-integer dressing**, precisely as
> conjectured in outlook §2.3.

**(c) Negative results (certified non-fits).** For the order-3 and order-6
characters at \(n = 9, 11, 13\) and the order-4 characters at \(n = 7, 15\),
safe PSLQ (250 digits, coefficient height \(10^6\), tolerance \(10^{-200}\))
finds **no** relation between \(S(\chi)\), \(L'(0,\chi)\), \(\{\log p\}\) over
\(p \mid 2n\,H(0)H(1728)\), and \(\log\varepsilon_{d_2}\): the dressing of the
non-real characters genuinely lives in the cubic/quartic subfields of the ring
class field (its exact description is the certified coset-cubic data of §3),
not in the rank-one basis of the naive conjecture. This is the expected
Stark-regime behaviour, and the master identity of Theorem 1 is exactly the
part that survives.

## 5. Hyperbolic: cancellation of \(\varepsilon\) and the odd-character law

**Lemma 5.1 (scale-invariant closed form).** For every primitive class \(f\) of
disc \(1 - n^2\), with \(\beta_1 = j(f)\), \(\beta_2 = j(\mathfrak{r}_nf)\) and
\(g(f) = y_f^6|\Delta_q(\tau_f)|\):
$$
|u_f|^6 \;=\;
\frac{\bigl|\beta_1^4(\beta_1-1728)^3\bigr|}{\bigl|\beta_2^4(\beta_2-1728)^3\bigr|}
\cdot \frac{g(f)}{g(\mathfrak{r}_nf)} .
$$
In particular **neither \(\varepsilon = n+\sqrt{n^2-1}\) nor \(\mu\) survives**
in \(|u_f|\).

*Proof.* From the proved closed form \(u_f = -\varepsilon\mu^{-2}h_2(\mathfrak{b}_1)/h_2(\mathfrak{b}_2)\)
and the Norm Lemma \(|\mu|^2 = \varepsilon q_1/q_2\)
([moduli-invariants.md](moduli-invariants.md) §5.5):
\(|u_f| = (q_2/q_1)\,|h_2(\mathfrak{b}_1)/h_2(\mathfrak{b}_2)|\). The weight algebra
\(h_2^6 = c\,j^4(j-1728)^3\Delta\cdot(\text{const})\) and
\(|\Delta(\mathfrak{b}_i)| = g(\mathfrak{b}_i)/\mathrm{covol}(\mathfrak{b}_i)^6\) with
\(\mathrm{covol}(\mathfrak{b}_1)/\mathrm{covol}(\mathfrak{b}_2) = q_2/q_1\)
(both centers have \(\operatorname{Im} = \sqrt{n^2-1}/2q_i\)) give
\(|h_2(\mathfrak{b}_1)/h_2(\mathfrak{b}_2)|^6 = (q_1/q_2)^6\,[\beta\text{-dressing}]\cdot g_1/g_2\);
assemble. Finally \([\mathfrak{b}_2] = [\mathfrak{r}_n][\mathfrak{b}_1]\)
([class-formula-proof.md](class-formula-proof.md)) identifies \(\beta_2, g_2\)
as the class functions at \(\mathfrak{r}_nf\). \(\square\)

(Verified classwise at \(n = 9, 11, 13, 15\), max deviation \(10^{-247}\).)

> **Theorem 4 (hyperbolic master identity).** Let \(\chi\) be a character of
> \(\mathrm{Cl}(1-n^2)\) (primitive classes), \(S(\chi) = \sum_f\chi(f)\log|u_f|\),
> \(\Sigma_x(\chi) = \sum_f\chi(f)\log|\beta_1(f) - x|\). Then
> $$
> S(\chi) = 0 \quad\text{if } \chi(\mathfrak{r}_n) = +1, \qquad
> \boxed{\;S(\chi) \;=\; -4\,L'(0,\chi) \;+\; \tfrac43\,\Sigma_0(\chi) \;+\; \Sigma_{1728}(\chi)\;}
> \quad\text{if } \chi(\mathfrak{r}_n) = -1 .
> $$

*Proof.* Sum \(6\log|u_f|\) from Lemma 5.1 against \(\chi\); re-indexing the
\(\mathfrak{r}_nf\)-terms by \(f \mapsto \mathfrak{r}_nf\) multiplies their sum
by \(\chi(\mathfrak{r}_n)^{-1} = \chi(\mathfrak{r}_n)\) (the twist class is
2-torsion), so
\(6S(\chi) = (1 - \chi(\mathfrak{r}_n))\bigl[4\Sigma_0 + 3\Sigma_{1728} + \Sigma_\gamma\bigr]\)
with \(\Sigma_\gamma(\chi) = \sum\chi(f)\log g(f) = -12L'(0,\chi)\) by
Proposition 1.1 (disc \(1-n^2\), \(w = 2\)). \(\blacksquare\)

The vanishing statement is of course law 2 (\(u_{\mathfrak{r}f}u_f = 1\)) seen
through the limit formula; both branches verified at \(n = 9, 11, 13, 15\), all
characters, residuals \(\le 10^{-247}\) (even-character sums vanish to
\(10^{-240}\) or better).

> **Theorem 5 (the twisted \(\Delta\)-ratio polynomial).** Set
> \(R_f := u_f^{6}\,\dfrac{\beta_2^4(\beta_2-1728)^3}{\beta_1^4(\beta_1-1728)^3}\)
> — a class function with \(|R_f| = g(f)/g(\mathfrak{r}_nf)\) (Lemma 5.1). By
> the first-power descent
> ([first-power-descent.md](first-power-descent.md) Thm 3.4) the multiset
> \(\{R_f\}\) is Galois-stable, so \(\prod_f(x - R_f) \in \mathbb{Q}[x]\); and
> for every odd character,
> $$
> -24\,L'(0,\chi) \;=\; \sum_f \chi(f)\,\log|R_f| .
> $$
> Certified (all levels, \(\ge 230\) spare digits): the coefficients are
> **integers**, the polynomials are **palindromic with constant term 1** —
> so the \(R_f\) are algebraic **units**, inverted in pairs by the
> \(\mathfrak{r}_n\)-twist:
> $$
> \begin{aligned}
> n=9:&\quad x^4 - 339524x^3 - 95354x^2 - 339524x + 1\\
> n=11:&\quad x^4 - 56529284x^3 + 1538876166x^2 - 56529284x + 1\\
> n=13:&\quad x^4 - 11382984004x^3 + 885435408006x^2 - 11382984004x + 1\\
> n=15:&\quad x^8 - 2628641876392x^7 - 21595933374628x^6 - 1373071731101336x^5\\
> &\qquad + 9740462908109254x^4 - (\text{sym})
> \end{aligned}
> $$
> with the \(L'\)-identity verified to \(10^{-247}\). This is the hyperbolic
> counterpart of Theorem 2, and sharper: the \(L'(0,\chi)\) of the odd
> characters are \(\chi\)-combinations of logs of genuine **units** — the
> \(\mathfrak{r}_n\)-twisted \(\Delta\)-ratios that the phase itself carries —
> with no \(S\)-integer dressing at all in the \(\Delta\)-part.

*Proof of the displayed identity.* \(\sum_f\chi(f)[\log g(f) - \log g(\mathfrak{r}_nf)]
= (1-\chi(\mathfrak{r}))\Sigma_\gamma(\chi) = 2\Sigma_\gamma(\chi) = -24L'(0,\chi)\)
for odd \(\chi\). Galois stability of \(\{R_f\}\): \(u_f\), \(\beta_1\),
\(\beta_2\) are strictly dihedrally equivariant (first-power descent), so
\(\sigma(R_f) = R_{f^{e(\sigma)}\mathfrak{c}(\sigma)}\). \(\square\)

## 6. Hyperbolic: genus characters in closed form

At \(n = 11, 13, 15\) some odd characters are real (genus) characters; the same
machinery as §4 then evaluates everything. The discovered exact
\(R_\chi = 2\,\mathrm{conv}\) decompositions (certified for all
\(m \le 300\) and **proved** by the same Sturm argument as §4(a) — the
bounds here are \(\le 32\); no conductor corrections arise for any odd
character, only for the even ones noted below):

| \(n\) | odd \(\chi\) | \(d_1 \cdot d_2 = 1-n^2\) | \(L'(0,\chi)\) | field |
|---|---|---|---|---|
| 11 | \(\chi_{(0,1)}\) | \((-15)\cdot 8\) | \(2\log(1{+}\sqrt2)\) | \(\mathbb{Q}(\sqrt2)\) |
| 11 | \(\chi_{(1,0)}\) | \((-3)\cdot 40\) | \(\tfrac23\log(3{+}\sqrt{10})\) | \(\mathbb{Q}(\sqrt{10})\) |
| 13 | \(\chi_{(0,1)}\) | \((-7)\cdot 24\) | \(\log(5{+}2\sqrt6)\) | \(\mathbb{Q}(\sqrt6)\) |
| 13 | \(\chi_{(1,0)}\) | \((-8)\cdot 21\) | \(\log\tfrac{5+\sqrt{21}}2\) | \(\mathbb{Q}(\sqrt{21})\) |
| 15 | \(\chi_{(0,1)}\) | \((-8)\cdot 28\) | \(\log(8{+}3\sqrt7)\) | \(\mathbb{Q}(\sqrt7)\) |
| 15 | \(\chi_{(2,1)}\) | \((-4)\cdot 56\) | \(\tfrac12\log(15{+}4\sqrt{14})\) | \(\mathbb{Q}(\sqrt{14})\) |

The **even** real characters pair with the complementary fields — at \(n = 11\)
the even character factors (certified, same procedure) as
\(L(s,(\tfrac{-24}\cdot))\,L(s,(\tfrac5\cdot))\): the field is \(\mathbb{Q}(\sqrt5)\),
exactly the invariant field of the pair-sums in
[moduli-invariants.md](moduli-invariants.md) §5.8; at \(n = 13\) as
\(L(s,(\tfrac{-3}\cdot))\,L(s,(\tfrac{56}\cdot))\), giving \(\mathbb{Q}(\sqrt{14})\)
— again the §5.8 field. So the parity \(\chi(\mathfrak{r}_n) = \pm1\) splits
the genus fields of \(1-n^2\) into the **pair-sum fields** (even; carry the
algebra of \(u + 1/u\)) and the **\(\log|u|\)-fields** (odd; carry the
character sums) — complementary halves of the same genus theory. At \(n = 9\),
disc \(-80\) (conductor 2), the single real character is even; its Epstein
\(L\) factors as
\(L(s,(\tfrac{-16}\cdot))L(s,(\tfrac5\cdot))\,(1 + 2^{-s} + 2\cdot4^{-s})\)
(certified), and at \(n = 15\) the even \(\chi_{(2,0)}\) as
\(L(s,(\tfrac{-7}\cdot))L(s,(\tfrac{32}\cdot))\,(1 - 2^{-s} + 2\cdot4^{-s})\) —
two more instances of the conductor Euler phenomenon.

The \(j\)-dressing of the odd real characters factors exactly as in §4(b), now
with coset products taken over \(\ker\chi\) — certified quadratic-integer
conjugates in the matching field. Where a split prime is non-principal the
generator is taken for \(\mathfrak{p}^{m}\) (\(m\) the class order) and the
exponent is \(e/m\); at \(n = 11\) in \(\mathbb{Q}(\sqrt{10})\) (class number
2), \(\pi_3 = 1 + \tfrac12\sqrt{40}\) generates \(\mathfrak{p}_3^2\):

| \(n\) | \(\chi\) | \(\Sigma_0(\chi)\) | \(\Sigma_{1728}(\chi)\) |
|---|---|---|---|
| 11 | \((-15)\cdot8\) | \(30\log\varepsilon_{8} + 3\lambda_{41} + 3\lambda_{89}\) | \(36\log\varepsilon_{8} + 2\lambda_{71}\) |
| 11 | \((-3)\cdot40\) | \(14\log\varepsilon_{40} + 3\lambda_{3} - 3\lambda_{41} - 3\lambda_{89}\) | \(16\log\varepsilon_{40} - 2\lambda_{71}\) |
| 13 | \((-7)\cdot24\) | \(18\log\varepsilon_{24} - 3\lambda_{5} + 3\lambda_{101}\) | \(12\log\varepsilon_{24} + 4\lambda_{19} + 2\lambda_{47} + 2\lambda_{167}\) |
| 13 | \((-8)\cdot21\) | \(24\log\varepsilon_{21} - 3\lambda_{5} - 3\lambda_{101}\) | \(24\log\varepsilon_{21} - 2\lambda_{47} - 2\lambda_{167}\) |
| 15 | \((-8)\cdot28\) | \(30\log\varepsilon_{28} - 6\lambda_{29} + 3\lambda_{47} - 6\lambda_{53} + 3\lambda_{167}\) | \(20\log\varepsilon_{28} + 2\lambda_{31} - 2\lambda_{47} - 2\lambda_{103} - 2\lambda_{199} - 2\lambda_{223}\) |
| 15 | \((-4)\cdot56\) | \(6\log\varepsilon_{56} - 6\lambda_{11} + 3\lambda_{47} - 3\lambda_{167}\) | \(2\log\varepsilon_{56} + 4\lambda_{11} + 2\lambda_{31} + 4\lambda_{43} - 2\lambda_{47} + 2\lambda_{103} - 2\lambda_{199} - 2\lambda_{223}\) |

(\(\lambda_p\) as in §4 with the canonical generators printed by the script;
split-prime supports again land in the \(\mathrm{GZ}(1-n^2, -3)\) resp.
\(\mathrm{GZ}(1-n^2,-4)\) sets together with \(p \mid 2(n^2-1)\).) Combining
with the table of \(L'\)-values, **\(S(\chi)\) is in closed form for every odd
real character of every computed level**, verified to \(10^{-247}\).

**The non-real odd characters.** At \(n = 9\) (\(\mathrm{Cl} = \mathbb{Z}/4\))
the two odd characters are quartic, and their sums collapse by law 2 to
\(S(\chi_4) = 2\log|u_1|\) — the master identity is then precisely the
classwise Lemma 5.1 at the principal class, and
\(-12L'(0,\chi_4) = \log|R\)-data\(|\) is a Stark-type evaluation in the quartic
ring class field: certified through Theorem 5's polynomial, while safe PSLQ
confirms (certified non-fit, as in §4(c)) that no rank-one
\(\{\log\varepsilon, \log p\}\)-expression exists. Same at \(n = 15\) for its
quartic odd pair.

## 7. The certified numerical record

Character-sum and \(L'\)-values (20 digits shown; computed at 250 digits, 400
at Euclidean \(n=15\); \(S(\chi) = S(\bar\chi)\), one value per conjugate
pair). Euclidean:

| \(n\) | \(\mathrm{ord}\,\chi\) | \(S(\chi)\) | \(L'(0,\chi)\) |
|---|---|---|---|
| 3 | 2 | 10.119185642834999673 | 0.43898596564160556954 |
| 5 | 2 | 16.401110452980199012 | 0.96242365011920689500 |
| 7 | 2 | 50.898733868194595161 | 2.7686593833135738327 |
| 7 | 4 | 22.684295755688513572 | 1.4860221248769270925 |
| 9 | 2 | 18.137584684821439049 | 1.7559438625664222782 |
| 9 | 3 | 70.144229691414912554 | 4.0476417264498781634 |
| 9 | 6 | 34.382429251891407159 | 2.1364594185869239177 |
| 11 | 2 | 50.213109260828194576 | 2.9932228461263808979 |
| 11 | 3 | 86.033492426467556694 | 5.3026856597888262543 |
| 11 | 6 | 27.769444924657409111 | 2.3032180902028192162 |
| 13 | 2 | 23.802362926332249068 | 2.3895264345742186082 |
| 13 | 3 | 101.24102291036467905 | 6.5582440777920829100 |
| 13 | 6 | 50.399596052674761580 | 3.3904644593763457318 |
| 15 | 2 \((0,1)\) | 62.742974077314449734 | 3.8496946004768275800 |
| 15 | 2 \((2,0)\) | 123.17609219442540179 | 8.2537482755822421869 |
| 15 | 2 \((2,1)\) | 48.769335467586396360 | 3.5118877251328445563 |
| 15 | 4 \((1,0)\) | 132.06845592069919001 | 8.2852528669601990706 |
| 15 | 4 \((1,1)\) | 39.877919196363264725 | 3.4800432906097017186 |

Hyperbolic (odd characters; the even ones have \(S = 0\) to \(10^{-240}\) or
better):

| \(n\) | \(\chi\) | \(S(\chi)\) | \(L'(0,\chi)\) |
|---|---|---|---|
| 9 | order 4 | 45.337283759443440194 | 1.0612750619050356520 |
| 11 | \((0,1)\) | 71.156367508225883636 | 1.7627471740390860505 |
| 11 | \((1,0)\) | 46.591472482040984242 | 1.2122976394880445490 |
| 13 | \((0,1)\) | 86.596498648832811158 | 2.2924316695611776878 |
| 13 | \((1,0)\) | 57.828173780263506754 | 1.5667992369724110787 |
| 15 | \((0,1)\) | 96.655012640029353845 | 2.7686593833135738327 |
| 15 | \((2,1)\) | 42.411391969797829839 | 1.7000422070566697504 |
| 15 | order 4 | 101.30861654091162593 | 2.5318972768037898006 |

(Two spot identities visible to the eye: Euclidean \(n=7\), order 2, and
hyperbolic \(n=15\), \(\chi_{(0,1)}\), share
\(L' = \log(8+3\sqrt7) = 2.76865938\ldots\) — the same Kronecker-symbol pair
\((-7)\cdot28\) resp. \((-8)\cdot28\)-data arriving at the same
\(\varepsilon_{28}\); and the Euclidean \(n=15\) values
\(8.2537\ldots = 4\log\varepsilon_{60}\), \(3.8497\ldots = 8\log\varepsilon_5\),
\(3.5119\ldots = \tfrac83\log\varepsilon_{12}\).)

## 8. What is proved, what is certified, what failed

**Proved** (for all levels of each aspect):
- Proposition 1.1, Lemma 2.1, Theorem 1, Theorem 2, Lemma 5.1, Theorem 4,
  Theorem 5 (the polynomial rationality granting the first-power descent of
  [first-power-descent.md](first-power-descent.md), which is proved).
- The genus factorizations of the §4/§6 tables, conductor Euler corrections
  included, at every displayed level — by the Sturm-bound argument of §4(a)
  (verified coefficients \(m \le 300\) against Sturm bounds \(\le 180\));
  hence all displayed \(L'(0,\chi_2)\) closed forms are theorems.

**Certified (exact integer/ideal arithmetic or \(\ge170\) spare digits), not
proved in general:**
- the \(D_n\)/\(R\)-polynomial coefficient values (each certified;
  \(D_n \in \mathbb{Z}[x]\) and \(\prod(x - R_f) \in \mathbb{Q}[x]\) are
  proved; the observed unit property of the \(R_f\) — palindromy, constant
  term 1, integrality — is certified per level, unproved in general);
- the two \(\varepsilon\)-ratio identifications at the composite conductor
  \(n = 15\) (\(L' = 8\log\varepsilon_5\), \(\tfrac83\log\varepsilon_{12}\);
  398 spare digits, but no factorization behind them yet);
- each coset-product factorization of §4(b)/§6 **is** individually proved
  (verified by exact arithmetic in its quadratic field); what is missing is a
  general law for the exponents \(k, e_p\) — a genus-refined Gross–Zagier
  statement.

**Certified negative results:** the \(S(\chi)\) of non-real characters match
nothing in the basis \(\{L'(0,\chi), \log p\ (p \mid 2nH(0)H(1728)),
\log\varepsilon\}\) (safe PSLQ, height \(10^6\), 250–400 digits) — the
dressing is genuinely of higher degree, as the certified cubics quantify. A
fit of the naive basis was also attempted and refuted for the hyperbolic
quartic characters at \(n = 9, 15\); and two of the three real characters at
Euclidean \(n = 15\) admit **no** finite-Euler genus factorization over the
divisor pairs of \(-900\).

**Open:**
1. The exponent law in the \(\Sigma\)-tables (genus-refined GZ; the
   \(3 \mid \Sigma_0\)-exponents and even \(\Sigma_{1728}\)-exponents are the
   \((4,3)\)-shadow, but the signs/unit powers await the height-theoretic
   derivation — outlook item for §1.3/§6.1 of the parent documents).
2. The general local Euler factor at conductor primes: a conceptual
   (proper-local-ideal) derivation of the corrections at square conductors,
   and the two-prime factorization shapes behind the certified
   \(\varepsilon\)-ratios at composite conductors (\(n = 15\)).
3. Stark-type closed forms for the non-real characters (order 3/4/6): the
   certified cubic/quartic coset data pins the objects; recognizing them as
   Robert/Stark units of the cubic and quartic subfields is the natural
   continuation (Schertz, *Complex Multiplication*, ch. 6–7).
4. The spherical aspect (disc \(-4(\ell^2+1)\), cap-swap \(\varepsilon_\ell^4\)
   as archimedean correction) — the machinery transports verbatim once the
   spherical \(u^2\)-laws are proved.
5. Even hyperbolic levels / \(i\mathcal{S}\), and the index of the subgroup
   generated by the \(\{G_\mathfrak{c}\}\) in the full elliptic-unit group of
   the ring class field (Robert index questions).

## 9. Machine verification

All statements above: `python3 scripts/phase_klf.py --selftest` (runs both
aspects at the default 250 digits over Euclidean \(n = 9, 11, 13, 7, 5, 3\) and
hyperbolic \(n = 9, 11, 13, 15\); Euclidean \(n = 15\) is a separate 400-digit
run, `python3 scripts/phase_klf.py --dps 400 euclid 15`). The script asserts,
per level:

1. the classwise laws (Lemma 2.1 / Lemma 5.1) including constants;
2. \(\Sigma_\gamma(\chi) = -12L'(0,\chi)\) against the independent
   incomplete-gamma evaluation of §1 (this is the KLF check proper);
3. the master identities (Theorems 1 and 4), every character;
4. the vanishing \(S(\chi) = 0\) for even hyperbolic characters;
5. the trivial-character anchor against the certified \(H_{-4n^2}\) and the
   \(\Delta\)-mass closed form;
6. the \(D_n\) certification with \(D_n(0) = \pm M(n)\), and the uniform Stark
   law \(-12L' = \sum\chi\log|G|\) for all characters at once;
7. the genus factorizations \(R_\chi = 2(\mathrm{conv}*\mathrm{corr})\) as exact
   integer identities (\(m \le 300\)) and the resulting \(L'\) closed forms;
8. the coset-product integer certifications and the **exact**
   \(\mathbb{Q}(\sqrt d)\)-verification of every \(A/B\)-factorization (unit
   power and split-prime exponents re-multiplied symbolically, valuations via
   Hensel lifts — no PSLQ anywhere in the exact chain);
9. the safe-PSLQ fits/non-fits of §4(c)/§6 (information-theoretic safety
   margin enforced: a relation is accepted only when its information content is
   \(\le \mathrm{dps}/4\) digits and the residual is at rounding level).

Guard rails: precision is set after imports; every integer/rational
certification uses the absolute-error criterion with \(\ge\max(20,\mathrm{dps}/5)\)
spare digits; all multi-term fits are safe-parameter; every PSLQ discovery that
is *used* is re-verified in exact arithmetic. Residual summary of the recorded
runs: identities \(\le 2\cdot10^{-248}\) at 250 digits and
\(\le 3\cdot10^{-398}\) at 400 (asserted below \(10^{-150}\) resp.
\(10^{-240}\)); even-character sums vanish to \(10^{-249}\) or exactly;
certifications carry 170–240+ spare digits (\(\ge 100\) at the 400-digit
level).
