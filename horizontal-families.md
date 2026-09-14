# Horizontal families: Atkin–Lehner-coupled Heegner points and the slices of the Schmidt unit systems

This document carries out item 4 of [PROGRAM.md](PROGRAM.md). The unit theorem of
Paper II (Theorem 4.2; [other-fields.md](other-fields.md) Theorem 4 for an arbitrary
twist) needs only an invertible ambiguous ideal \(\mathfrak{r}\) of an order
\(\mathcal{O}_D\) of an imaginary quadratic field: for every such pair the twisted
\(\Delta\)-ratios
$$
R_\mathfrak{b} \;=\; N(\mathfrak{r})^6\,\frac{\Delta(\mathfrak{b})}{\Delta(\mathfrak{r}^{-1}\mathfrak{b})},
\qquad \mathfrak{b} \in \mathrm{Pic}(\mathcal{O}_D),
$$
are units of the ring class field \(H_D\). We call the pair \((D, [\mathfrak{r}])\) a
**slice** and \(\{R_\mathfrak{b}\}\) its **unit system**. The Gaussian arrangement
\(\mathcal{S}\) realizes the slices \(D = -4r_0(r_0+1)\), \([\mathfrak{r}] = [(r_0, 0, r_0+1)]\),
one per odd level \(n = 2r_0+1\) (Papers I–II); its \(\mathrm{PGL}_2\)-companion
\(i\mathcal{S}\) realizes \(D = 1-n^2\) at the even levels; Stange's arrangements
\(\mathcal{S}_K\) realize \(D = -(k^2-4)/|d_K|\) at the level \(2\alpha = k\)
([other-fields.md](other-fields.md)). Here the family of all slices is described, the
arrangements (primitive strata, imprimitive strata, the even levels of
\(i\mathcal{S}\), the principal cusps of class-number-two fields, and the second-kind
orbit of the imaginary axis) are matched to the slices they carry, and the limit
formulas, the odd index and the phase are transported to arbitrary slices. In the
language of the program's one-sentence line: the pair
\((\mathfrak{b}, \mathfrak{r}^{-1}\mathfrak{b})\) is a Heegner point on
\(X_0(r_0)\times_{X(1)}X_0(s_0)\) (\(\mathfrak{r}\mathfrak{s} = (\theta)\), \(r_0s_0 = N(\theta)\)),
swapped by the Atkin–Lehner involutions, and an arrangement is a way of coupling the two
levels \(r_0, s_0\) (\(s_0 = r_0+1\) at \(\mathbb{Q}(i)\); in general
\(g_+s_0 - g_-r_0 = 4\) with \(g_-g_+ \in \{|d_K|, 4|d_K|\}\)). Everything finite is
verified by [scripts/horizontal_families.py](scripts/horizontal_families.py)
(`--selftest`, §11) under the guard rails of [CLAUDE.md](CLAUDE.md); §10 separates
proved, certified and experimental, and records where the expected statements of the
session prompt failed.

**The results in one paragraph.** *Slices* (§1): the primitive ambiguous invertible
ideals of \(\mathcal{O}_D\) are the ideals of the forms \((r_0, 0, s_0)\) and
\((r_0, r_0, \tfrac{r_0+s_0}{4})\); they number \(2^\mu\) (Gauss's \(\mu\)), come in
pairs \(\{\mathfrak{r}, \mathfrak{s}\}\) with \(\mathfrak{r}\mathfrak{s} = (\theta)\), \(\theta\)
purely imaginary of norm \(r_0s_0\), the pairs are in bijection with \(\mathrm{Cl}(D)[2]\),
and \([\mathfrak{r}] = 1\) iff \(\min(r_0, s_0) = 1\); the unit system depends only on the
slice and has rank \(h/2\) (\([\mathfrak{r}] \ne 1\)) or \(0\). *Realization* (§2): the
levels of \(\mathcal{S}_K\) in \(\mathbb{H}\) are \(k = t|d_K| \pm 2\) (odd \(d_K\)) and
\(k = 4mt \pm 2\) (\(d_K = -4m\)), \(t \ge 1\), and the \(t\)-th level pair carries the
slices with pairs \(\{t, t|d_K| \pm 4\}\) resp. \(\{t, mt \pm 1\}\): **the twist of the
level \(t|d_K|\pm2\) has norm \(t\)**; \(i\mathcal{S}\) at \(n = 2t\) carries
\(\{2t-1, 2t+1\}\). Hence \(\mathcal{S}_K\) has a level of discriminant \(D\) iff
\(|d_K||D|+4\) is a square, a slice with pair \(\{r_0, s_0\}\) is realized by the odd
\(d_K\) with \(|d_K| = (s_0 \mp 4)/r_0\) or \((r_0 \mp 4)/s_0\), by \(\mathbb{Q}(i)\) iff
\(|r_0 - s_0| \in \{1, 2\}\) (\(\mathcal{S}\), \(i\mathcal{S}\)) and by
\(\mathbb{Q}(\sqrt{-2})\) iff \(s_0 = 2r_0 \pm 1\) or \(r_0 = 2s_0 \pm 1\) — with the
type (\(D = -r_0s_0\) or \(-4r_0s_0\)) matching. Equal slices give equal unit systems
(\(D = -15\) at four arrangements, \(-20, -24, -36\) at two), different classes of the
same \(D\) give different systems with different odd characters (\(D = -84\) at three
fields, \(-160, -180, -195, -224, -255, -260, -276, -288\) at two). The same theorems hold
at the principal cusp of every \(K\): the class-number-two fields realize new slices
(\((-168, [(3,0,14)])\) at \(\mathbb{Q}(\sqrt{-5})\), \((-88, [(2,0,11)])\) at
\(\mathbb{Q}(\sqrt{-6})\), \((-168, [(2,0,21)])\) at \(\mathbb{Q}(\sqrt{-10})\)). *Strata* (§3):
the stratum of conductor \(f\) at level \(k\) carries the slice
\((D/f^2, [\mathfrak{r}_\alpha\mathcal{O}_{D/f^2}])\); \(D'\) occurs at
\(\mathcal{S}_K\) iff \(k^2 - |d_K||D'|f^2 = 4\), infinitely often unless \(K = \mathbb{Q}(\sqrt{D'})\),
and along the powers \(\eta^j\) of the fundamental solution the stratum twist class is
\(c^j\): trivial for even \(j\), and for odd \(j\) the class \(c\) of the pair whose
\(r_0\)-side contains \(p^{v_p(D')}\) exactly for the odd primes \(p\) with
\(\eta \equiv 1 \bmod \mathfrak{p}\) (proved; the \(2\)-part certified). So every
\(\mathcal{S}_K\) realizes, on its strata of discriminant \(D'\), exactly one class of
\(\mathrm{Cl}(D')[2]\) besides the trivial one — e.g. \((-56, [(2,0,7)])\) at nine fields,
but \((-160, [(5,0,8)])\) at none. *The even levels of \(i\mathcal{S}\)* (§5): the class
formula \(\hat\sigma[f] = [(n-1, n-1, \tfrac n2)][f]^{-s}\) is proved by Paper I's Lemma A
on the coset \(\det X = i\) with the odd-\(D\) Gram form \(nN(u) - s\operatorname{Im}(u^2)\);
the unit polynomial at \(n = 4\) is \(x^2 - 7x + 1\), as predicted by the coincidence with
\(\mathbb{Q}(\sqrt{-3})\); the ledger item "even levels / \(i\mathcal{S}\)" is closed. *The
second-kind orbit* \(\mathcal{S}_K^\perp = \mathrm{PSL}_2(\mathcal{O}_K)\cdot i\hat{\mathbb{R}}\)
(§6) is classified for the five Euclidean fields: its level-\(y\) circles are all forms
of discriminant \(4 - y^2|d_K|\), a family disjoint from that of \(\mathcal{S}_K\), but it
has no level-preserving involution and hence carries discriminants, not slices.
*Genus* (§7): the genus vector of \(\mathfrak{r}_\alpha\) is \((g_+/p)\) at the odd primes
of \(r_0\) and \((-g_-/p)\) at those of \(s_0\). *Index* (§8): Theorem 3 of
[robert-index-full.md](robert-index-full.md) holds on every slice, exact and unconditional
at twelve slices. *Phase* (§9): \(u_\mathfrak{b} = \Phi_y/\Phi_x(\beta_1,\beta_2)\),
\(\Phi = \Phi_{r_0}\), is defined on every slice with \(D < -4\) (the branch through the
Atkin–Lehner pair is unique for the smaller-norm twist), satisfies the first-power
Galois law and \(u^6 = R\cdot j\)-dressing, and its exact slice polynomials are
irreducible at all sixteen computed slices.

**Notation.** \(D < 0\) a discriminant, \(\mathcal{O}_D\) the order, \(\mathrm{Cl}(D)\)
its class group of order \(h\), \(H_D\) the ring class field over \(K' = \mathbb{Q}(\sqrt D)\).
A form \((a, b, c)\) has the ideal \(\mathfrak{a} = [a, \tfrac{-b+\sqrt D}{2}]\) (the convention of
`scripts/other_fields.py`; composition of forms is composition of ideals). A
**primitive ambiguous ideal** is a lattice \([a, \tfrac{-b+\sqrt D}{2}]\) (not divisible by a
rational integer \(> 1\)) that is an invertible \(\mathcal{O}_D\)-ideal equal to its conjugate.
\(\theta_D := \sqrt D\) (\(D\) odd) or \(\sqrt D/2\) (\(D\) even). For an imaginary quadratic
field \(K\): \(d_K\), \(\mathcal{O}_K = \mathbb{Z}[\omega]\), \(\sqrt{d_K} = i\sqrt{|d_K|}\),
\(\mathcal{S}_K = \mathrm{PSL}_2(\mathcal{O}_K)\cdot\hat{\mathbb{R}}\), levels \(2\alpha = k\),
\(D_K(\alpha) = -(k^2-4)/|d_K|\), the twist ideal \(\mathfrak{r}_\alpha\) of norm
\(r_0 = (k-2)/g_-\), partner \(\mathfrak{s}_\alpha\) of norm \(s_0 = (k+2)/g_+\),
\(g_\mp = \gcd(k \mp 2, |d_K|)\), all as in [other-fields.md](other-fields.md) §2.
Status conventions as in the papers: **proved**, **certified** (finite statement
verified by the script), **experimental**.

## 1. Slices

> **Theorem 1 (slices).** Let \(D < -4\) be a discriminant.
> 1. The primitive ambiguous invertible ideals of \(\mathcal{O}_D\) are exactly the
>    ideals of the primitive forms \((r_0, 0, s_0)\) with \(4r_0s_0 = |D|\) and
>    \((r_0, r_0, \tfrac{r_0+s_0}{4})\) with \(r_0s_0 = |D|\), \(\gcd(r_0, s_0) = 1\) in both
>    cases (write \(\mathfrak{r} = \mathfrak{r}(r_0, s_0)\) for the ideal, of norm \(r_0\)).
> 2. There are \(2^\mu\) of them, \(\mu\) the number of assigned characters of \(D\)
>    (the number of odd primes dividing \(D\), plus \(0, 1, 1, 2\) for
>    \(D/4 \equiv 3, \{1, 2\}, 4, 0 \pmod{4, 4, 8, 8}\) when \(D \equiv 0 \bmod 4\)).
> 3. They come in pairs: \(\mathfrak{r}(r_0, s_0)\) and \(\mathfrak{s} = \mathfrak{r}(s_0, r_0)\)
>    satisfy \(\mathfrak{r}\mathfrak{s} = (\theta)\) with \(\theta = \sqrt D/2\) resp. \(\sqrt D\)
>    purely imaginary, \(N(\theta) = r_0s_0\); and \(\mathfrak{r}^2 = (r_0)\), \(\mathfrak{s}^2 = (s_0)\).
> 4. Two primitive ambiguous ideals in the same class are the two members of one
>    pair; the map pair \(\mapsto\) class is a bijection onto \(\mathrm{Cl}(D)[2]\),
>    which has order \(2^{\mu-1}\).
> 5. \([\mathfrak{r}] = 1\) iff \(\min(r_0, s_0) = 1\).
> 6. The unit system \(\{R_\mathfrak{b}\}\) of \(\mathfrak{r}\) is unchanged under
>    \(\mathfrak{r} \mapsto \mathfrak{s}\) and under rational scaling, so it depends only on the
>    slice \((D, [\mathfrak{r}])\); the \(R_\mathfrak{b}\) are Galois conjugate
>    (\(\sigma_\mathfrak{a}R_\mathfrak{b} = R_{\mathfrak{a}^{-1}\mathfrak{b}}\)), and the group they
>    generate has rank \(h/2\) if \([\mathfrak{r}] \ne 1\) and rank \(0\) (all
>    \(R_\mathfrak{b} = 1\)) if \([\mathfrak{r}] = 1\).

*Proof.* (1) A primitive ideal has a basis \([a, \tfrac{-b+\sqrt D}{2}]\) with
\(b^2 - 4ac = D\); it is invertible iff \(\gcd(a, b, c) = 1\), and its conjugate is
\([a, \tfrac{b+\sqrt D}{2}]\), which is the same lattice iff \(b \equiv -b \bmod 2a\), i.e.
\(a \mid b\), i.e. \(b \in \{0, a\} \bmod 2a\). With \(b = 0\): \(4ac = |D|\); with
\(b = a\): \(a(4c - a) = |D|\); primitivity of the form is \(\gcd(a, c) = 1\), which for
\(b = 0\) is \(\gcd(r_0, s_0) = 1\) and for \(b = a\), \(s_0 := 4c - a\), is
\(\gcd(r_0, s_0) = 1\) (a prime dividing \(a\) and \(c\) divides \(4c - a\), and conversely
a prime dividing \(a\) and \(4c-a\) divides \(4c\), hence \(c\) if odd; \(2 \mid a\) and
\(2 \mid s_0\) with \(s_0 = 4c - a\) forces \(2 \mid c\) when \(4 \mid a\), while
\(a \equiv 2 \bmod 4\) gives \(s_0 \equiv 2 \bmod 4\) and then \(\gcd(a, c) = 1\) iff \(c\)
odd). (2) *Odd \(D = -m\)*: only the second type; \(a\) runs over the divisors of \(m\),
\(a + m/a \equiv 0 \bmod 4\) is automatic (\(a\cdot m/a \equiv 3 \bmod 4\)), and
\(\gcd(a, c) = 1\) iff \(\gcd(a, m/a) = 1\): \(2^{\omega(m)}\) unitary divisors, \(\mu = \omega(m)\).
*Even \(D = -4n\)*: the first type gives the \(2^{\omega(n)}\) unitary factorizations \(n = ac\).
For the second type \(a\) is even, \(a = 2a'\), \(a'(2c - a') = n\) with \(c\) odd and
\(\gcd(a', c) = 1\); writing \(n = a'e\), \(c = (a'+e)/2\) must be an odd integer:
for \(n \equiv 1 \bmod 4\) every unitary factorization works (\(a' \equiv e \bmod 4\)),
for \(n \equiv 3 \bmod 4\) none (\(a' + e \equiv 0 \bmod 4\)), for \(n \equiv 2 \bmod 4\) none
(\(a', e\) would both have to be even), for \(n \equiv 4 \bmod 8\) none (\(a' = 2a''\),
\(e = 2e''\), \(a''e'' = n/4\) odd, \(a''+e''\) even), and for \(8 \mid n\) exactly the
\(2^{\omega(n/4)} = 2^{\omega(n)}\) unitary factorizations of \(n/4\). Adding the two types
gives \(2^{\omega(n)}\cdot(1, 2, 1, 1, 2)\) in the five cases, i.e. \(2^\mu\) with the
\(\mu\) of Cox, Thm 3.15. (3) \([r_0, \tfrac{\sqrt D}{2}]\cdot[s_0, \tfrac{\sqrt D}{2}]
= [r_0s_0, \tfrac{\sqrt D}{2}]\) (using \(\gcd(r_0, s_0) = 1\) and \(D/4 = -r_0s_0\)) \(= \tfrac{\sqrt D}{2}\mathcal{O}_D\);
for the second type the four products of generators divided by \(\sqrt D\) are
\(-\sqrt D\), \(\tfrac{r_0 - \sqrt D}{2}\), \(\tfrac{s_0 - \sqrt D}{2}\) and \(\tfrac{r_0+s_0}{4}\) (using
\(r_0s_0 = |D|\), \(r_0, s_0\) odd), all in \(\mathcal{O}_D\), so the product is contained in
\(\sqrt D\mathcal{O}_D\) and has the same norm \(r_0s_0 = |D|\).
\(\mathfrak{r}^2 = (r_0)\) is the same computation with \(s_0 = r_0\) formally, i.e.
\([r_0^2, r_0\tfrac{\pm r_0 + \sqrt D}{2}, (\tfrac{r_0+\sqrt D}{2})^2] = r_0\mathcal{O}_D\).
(4) If \(\mathfrak{r}' = \lambda\mathfrak{r}\) with both ambiguous then
\(\lambda\mathfrak{r} = \bar\lambda\mathfrak{r}\), so \(\lambda/\bar\lambda \in \mathcal{O}_D^\times = \{\pm1\}\)
(\(D < -4\)): \(\lambda\) is rational, and primitivity gives \(\mathfrak{r}' = \mathfrak{r}\); or
\(\lambda\) is purely imaginary, \(\lambda = q\theta\) with \(q \in \mathbb{Q}\), and
\(\mathfrak{r}' = q\theta\mathfrak{r} = q\,r_0\,\mathfrak{s}\) is a rational multiple of \(\mathfrak{s}\), hence
\(= \mathfrak{s}\). So pair \(\mapsto\) class is injective; its image lies in \(\mathrm{Cl}(D)[2]\)
by (3); both sets have \(2^{\mu-1}\) elements (Gauss: the number of genera is
\(2^{\mu-1}\), and \(|\mathrm{Cl}[2]| = |\mathrm{Cl}/\mathrm{Cl}^2|\)). (5) If \(r_0 = 1\) then
\(\mathfrak{r} = \mathcal{O}_D\); if \([\mathfrak{r}] = 1\), \(\mathfrak{r} = (\lambda)\) with \((\bar\lambda) = (\lambda)\),
so \(\lambda\) is rational or purely imaginary as in (4), giving \(\mathfrak{r} = \mathcal{O}_D\) or
\(\mathfrak{s} = \mathcal{O}_D\). (6) Independence of the partner and of scaling is Theorem 4 of
[other-fields.md](other-fields.md); the Galois law is Paper II's (Theorem 2.9 / the
Shimura reciprocity of \(\Delta\)-quotients). For the rank: the places of \(H_D\) are
permuted simply transitively by \(\mathrm{Cl}(D)\) ([robert-index-full.md](robert-index-full.md)
Lemma 1.1), so the rank of \(\langle R_\mathfrak{b}\rangle\) is the rank of the group matrix
\((\log|R_{\mathfrak{a}^{-1}\mathfrak{b}}|)_{\mathfrak{a},\mathfrak{b}}\), which is the number of
characters \(\chi\) with \(S_\chi := \sum_\mathfrak{b}\chi(\mathfrak{b})\log|R_\mathfrak{b}| \ne 0\);
by Theorem 4 of other-fields.md, \(S_\chi = (1 - \chi(\mathfrak{r}))(-12L'(0,\chi))\) and
\(L'(0,\chi) \ne 0\) for \(\chi \ne 1\) ([robert-index-full.md](robert-index-full.md)
Cor. 3.3), so the count is the number of \(\chi\) with \(\chi(\mathfrak{r}) = -1\): a coset of
the index-two subgroup \(\{\chi(\mathfrak{r}) = 1\}\) when \([\mathfrak{r}] \ne 1\), empty otherwise
(and then \(\mathfrak{r} = \mathcal{O}_D\) or \(\mathfrak{s} = \mathcal{O}_D\) by (5), so every
\(R_\mathfrak{b} = 1\), directly or by the partner independence). \(\blacksquare\)

**Certified** (phase S): at all \(198\) discriminants \(-400 \le D < -4\) the
enumeration of the two types gives \(2^\mu\) primitive forms whose ideals are ambiguous
(exact HNF arithmetic in \(\mathbb{Q}(\sqrt D)\)) with \(\mathfrak{r}^2 = (r_0)\),
\(\mathfrak{r}\mathfrak{s} = (\theta)\); the \(398\) pairs (\(200\) with nontrivial class) hit
\(\mathrm{Cl}(D)[2]\) bijectively, \([\mathfrak{r}] = 1\) iff \(\min = 1\), and the ambiguous forms
of the third kind \((a, b, a)\) are never ambiguous ideals. Examples:

| \(D\) | \(\mathrm{Cl}(D)\) | pairs \(\{r_0, s_0\}\) \(\to\) class |
|---|---|---|
| \(-15\) | \(\mathbb{Z}/2\) | \(\{1,15\}\to1\), \(\{3,5\}\to(2,1,2)\) |
| \(-20\) | \(\mathbb{Z}/2\) | \(\{1,5\}\to1\), \(\{2,10\}\to(2,2,3)\) |
| \(-56\) | \(\mathbb{Z}/4\) | \(\{1,14\}\to1\), \(\{2,7\}\to(2,0,7)\) |
| \(-84\) | \((\mathbb{Z}/2)^2\) | \(\{1,21\}\to1\), \(\{2,42\}\to(2,2,11)\), \(\{3,7\}\to(3,0,7)\), \(\{6,14\}\to(5,4,5)\) |
| \(-160\) | \((\mathbb{Z}/2)^2\) | \(\{1,40\}\to1\), \(\{4,40\}\to(4,4,11)\), \(\{5,8\}\to(5,0,8)\), \(\{8,20\}\to(7,6,7)\) |
| \(-168\) | \((\mathbb{Z}/2)^2\) | \(\{1,42\}\to1\), \(\{2,21\}\to(2,0,21)\), \(\{3,14\}\to(3,0,14)\), \(\{6,7\}\to(6,0,7)\) |

The two types coexist at even \(D\) (\(-20\): \(\{1,5\}\) of the first, \(\{2,10\}\) of
the second type), and the pair does not determine \(D\): \(\{3, 5\}\) is a pair of
\(D = -15\) (second type) and of \(D = -60\) (first type). The type is part of the slice
(it is read off from \(r_0s_0 \in \{|D|, |D|/4\}\)).

**Remark (what is classical).** (1)–(5) are Gauss's theory of ambiguous forms
(*Disquisitiones*, Art. 257–262: the forms \((a, b, c)\) with \(a \mid b\), their number,
and the fact that a class of order \(\le 2\) contains exactly two of them up to the
normalization \(b \in \{0, a\}\)) in ideal language; the count (2) is the standard one
behind Cox's Theorem 3.15. The unit property (6) at the maximal order is the classical
statement that \(\Delta(\mathfrak{a})/\Delta(\mathfrak{b})\) lies in the Hilbert class field and
generates \((\mathfrak{b}\mathfrak{a}^{-1})^{12}\) (Lang, *Elliptic Functions*, Ch. 12; see §12): with
\(\mathfrak{a} = \mathfrak{b}\), \(\mathfrak{b}' = \mathfrak{r}^{-1}\mathfrak{b}\) it generates \(\mathfrak{r}^{-12} = (r_0)^{-6}\), so
\(R_\mathfrak{b}\) is a unit. For non-maximal orders the same conclusion is Paper II's
five-lemma proof (other-fields.md Theorem 4), which uses only the invertibility of
\(\mathfrak{r}\); what is new in this section is the bookkeeping of the pairs and of the type,
and the rank statement. Nothing here uses an arrangement.

## 2. Realization at primitive strata

For a level \(2\alpha = k\) of \(\mathcal{S}_K\) write \(g_\mp = \gcd(k \mp 2, |d_K|)\),
\(r_0 = (k-2)/g_-\), \(s_0 = (k+2)/g_+\) as in other-fields.md; then \(g_+s_0 - g_-r_0 = 4\)
and \(g_-g_+ \in \{|d_K|, 4|d_K|\}\), the two cases giving the two types of Theorem 1(1).

> **Theorem 2 (the levels of \(\mathcal{S}_K\) and their slices).** Let \(K\) be any
> imaginary quadratic field (no restriction on \(h_K\); at \(h_K > 1\) the statements
> concern the orbit of \(\hat{\mathbb{R}}\), the principal cusp).
> 1. *(Odd \(d_K = -p\).)* The levels of \(\mathcal{S}_K\) in \(\mathbb{H}\) are exactly
>    \(k = tp \pm 2\), \(t \ge 1\) (\(k \ge 3\)); the level \(tp + 2\) has orientation sign
>    \(s = -1\), the level \(tp - 2\) has \(s = +1\); the level \(tp \pm 2\) carries the
>    discriminant \(D = -t(tp \pm 4)\) with the twist \(\mathfrak{r}_\alpha\) of the pair
>    \(\{t, tp \pm 4\}\) (second type, form \((t, t, \tfrac{t(p+1)\pm4}{4})\)): **the twist of
>    the \(t\)-th level pair has norm \(t\)**.
> 2. *(Even \(d_K = -4m\), \(m \ge 2\).)* The levels are \(k = 4mt \pm 2\), \(t \ge 1\), with
>    \(s = \mp1\), discriminant \(D = -4t(mt \pm 1)\) and the twist of the pair
>    \(\{t, mt \pm 1\}\) (first type, form \((t, 0, mt \pm 1)\)).
> 3. *(\(\mathbb{Q}(i)\).)* \(\mathcal{S}\): \(k = 2n\), \(n = 2t+1\) odd, both orientations, pair
>    \(\{t, t+1\}\) (first type, \(D = -4t(t+1) = 1 - n^2\)); \(i\mathcal{S}\): \(k = 2n\),
>    \(n = 2t\) even, both orientations, pair \(\{2t-1, 2t+1\}\) (second type,
>    \(D = 1 - n^2\), form \((n-1, n-1, \tfrac n2)\)) — §5.
> 4. In every case the level-\(k\) circles in \(\mathbb{H}\) are all the forms of discriminant
>    \(D\) (census \(3H(|D|)\)) and the involution acts by
>    \(\hat\sigma[f] = [\mathfrak{r}_\alpha][f]^{-s}\); the twist is trivial exactly at \(t = 1\).
> 5. *(Criterion.)* \(\mathcal{S}_K\) has a level of discriminant \(D\) at a primitive stratum
>    iff \(|d_K||D| + 4 = k^2\) is a perfect square (then \(k \equiv \pm2 \bmod |d_K|\) is
>    automatic for the nine class-number-one fields; at \(\mathbb{Q}(i)\), \(k \equiv 0 \bmod 4\)
>    gives \(i\mathcal{S}\)). A slice with pair \(\{r_0, s_0\}\) is realized at a primitive
>    stratum of the arrangement of a class-number-one field iff: second type and
>    \(|d_K| = (s_0 \mp 4)/r_0\) or \((r_0 \mp 4)/s_0 \in \{3, 7, 11, 19, 43, 67, 163\}\); or
>    second type and \(|r_0 - s_0| = 2\) (\(i\mathcal{S}\)); or first type and
>    \(s_0 = r_0 \pm 1\) (\(\mathcal{S}\)) or \(s_0 = 2r_0 \pm 1\) or \(r_0 = 2s_0 \pm 1\)
>    (\(\mathbb{Q}(\sqrt{-2})\)).

*Proof.* The levels are the \(k \ge 3\) with \(2\alpha \equiv -2s \bmod |d_K|\)
(other-fields.md, Theorem 1(4), whose necessity part holds for every \(K\)); for odd
\(p = |d_K|\) this is \(k = tp \pm 2\), for \(|d_K| = 4m\) it is \(k = 4mt \pm 2\) (for \(m = 1\)
the two parametrizations coincide and one takes \(k = 4t+2\)). Every form of discriminant
\(D_K(\alpha)\) is realized at level \(k\) by the explicit \(P \in \mathrm{SL}_2(\mathcal{O}_K)\) of
Lemma A′ of other-fields.md, which uses neither \(h_K = 1\) nor the Euclidean property
(it only needs the lattice \(\mathcal{K}_f\), whose index is computed by a Smith form); the
class formula and its Lemmas B′, C′ are likewise general. The twist data: at
\(k = tp + 2\), \(g_- = \gcd(tp, p) = p\), \(g_+ = \gcd(tp+4, p) = 1\), so \(r_0 = t\),
\(s_0 = tp + 4\), \(g_-g_+ = p\) (second type), \(D = -r_0s_0\); at \(k = tp - 2\),
\(g_- = 1\), \(g_+ = p\), \(r_0 = tp - 4\), \(s_0 = t\); \(t(p+1) \pm 4 \equiv 0 \bmod 4\) since
\(p \equiv 3 \bmod 4\). At \(k = 4mt + 2\): \(g_- = 4m\), \(g_+ = 4\gcd(mt+1, m) = 4\), so
\(r_0 = t\), \(s_0 = mt+1\), \(g_-g_+ = 16m = 4|d_K|\) (first type); at \(k = 4mt - 2\)
symmetrically. \(\mathbb{Q}(i)\): Paper I and §5. (4) is Theorems 2–3 of other-fields.md
and Theorem 4 below; trivial twist iff \(\min(r_0, s_0) = 1\) iff \(t = 1\) (Theorem 1(5)).
(5) restates (1)–(3). \(\blacksquare\)

**Certified** (phases R and X′): (R1) the criterion at every \(D \ge -300\) and every
class-number-one field; (R5) the \(t\)-parametrization at every level \(k \le 40|d_K|\)
and the converse criterion at all \(287\) pairs with \(|D| \le 300\); the class
formula with Lemmas A′–C′ through the explicit \(P\) (no descent) at the first four
levels of \(\mathbb{Q}(\sqrt{-19}), \mathbb{Q}(\sqrt{-43}), \mathbb{Q}(\sqrt{-67}), \mathbb{Q}(\sqrt{-163})\)
(X4; e.g. \(k = 324\) at \(\mathbb{Q}(\sqrt{-163})\): \(D = -644\), \(h = 16\), twist \((2,2,81)\))
and at \(14\) Euclidean levels against the descent route of `other_fields.py` (X3). The
first slices of each field (\(k\): \(D\), twist \(\sim\) class):

| \(K\) | first levels |
|---|---|
| \(\mathbb{Q}(\sqrt{-3})\) | 4: \(-4\) (trivial); 5: \(-7\) (trivial); 7: \(-15\) \((5,5,2)\sim(2,1,2)\); 8: \(-20\) \((2,2,3)\); 10: \(-32\) \((8,8,3)\sim(3,2,3)\); 11: \(-39\) \((3,3,4)\); 13: \(-55\) \((11,11,4)\sim(4,3,4)\); 14: \(-64\) \((4,4,5)\); 16: \(-84\) \((14,14,5)\sim(5,4,5)\); 17: \(-95\) \((5,5,6)\) |
| \(\mathbb{Q}(i)\) | 6: \(-8\) (trivial); 8: \(-15\) \((3,3,2)\sim(2,1,2)\) [\(i\mathcal{S}\)]; 10: \(-24\) \((2,0,3)\); 12: \(-35\) \((5,5,3)\sim(3,1,3)\); 14: \(-48\) \((3,0,4)\); 16: \(-63\) \((7,7,4)\sim(4,1,4)\); 18: \(-80\) \((4,0,5)\); 20: \(-99\) \((9,9,5)\sim(5,1,5)\); 22: \(-120\) \((5,0,6)\); 24: \(-143\) \((11,11,6)\sim(6,1,6)\) |
| \(\mathbb{Q}(\sqrt{-7})\) | 9: \(-11\) (trivial); 12: \(-20\) \((10,10,3)\sim(2,2,3)\); 16: \(-36\) \((2,2,5)\); 19: \(-51\) \((17,17,5)\sim(3,3,5)\); 23: \(-75\) \((3,3,7)\); 26: \(-96\) \((24,24,7)\sim(4,4,7)\); 30: \(-128\) \((4,4,9)\); 33: \(-155\) \((31,31,9)\sim(5,5,9)\); 37: \(-195\) \((5,5,11)\); 40: \(-228\) \((38,38,11)\sim(6,6,11)\) |
| \(\mathbb{Q}(\sqrt{-2})\) | 6: \(-4\); 10: \(-12\) (trivial); 14: \(-24\) \((3,0,2)\sim(2,0,3)\); 18: \(-40\) \((2,0,5)\); 22: \(-60\) \((5,0,3)\sim(3,0,5)\); 26: \(-84\) \((3,0,7)\); 30: \(-112\) \((7,0,4)\sim(4,0,7)\); 34: \(-144\) \((4,0,9)\); 38: \(-180\) \((9,0,5)\sim(5,0,9)\); 42: \(-220\) \((5,0,11)\) |
| \(\mathbb{Q}(\sqrt{-11})\) | 9: \(-7\) (trivial); 13: \(-15\) (trivial); 20: \(-36\) \((18,18,5)\sim(2,2,5)\); 24: \(-52\) \((2,2,7)\); 31: \(-87\) \((29,29,8)\sim(3,3,8)\); 35: \(-111\) \((3,3,10)\); 42: \(-160\) \((40,40,11)\sim(4,4,11)\); 46: \(-192\) \((4,4,13)\); 53: \(-255\) \((51,51,14)\sim(5,5,14)\); 57: \(-295\) \((5,5,16)\) |
| \(\mathbb{Q}(\sqrt{-19})\) | 17: \(-15\) (trivial); 21: \(-23\) (trivial); 36: \(-68\) \((34,34,9)\sim(2,2,9)\); 40: \(-84\) \((2,2,11)\); 55: \(-159\) \((53,53,14)\sim(3,3,14)\); 59: \(-183\) \((3,3,16)\); 74: \(-288\) \((72,72,19)\sim(4,4,19)\); 78: \(-320\) \((4,4,21)\) |
| \(\mathbb{Q}(\sqrt{-43})\) | 41: \(-39\) (trivial); 45: \(-47\) (trivial); 84: \(-164\) \((2,2,21)\); 88: \(-180\) \((2,2,23)\); 127: \(-375\) \((3,3,32)\); 131: \(-399\) \((3,3,34)\) |
| \(\mathbb{Q}(\sqrt{-67})\) | 65: \(-63\) (trivial); 69: \(-71\) (trivial); 132: \(-260\) \((2,2,33)\); 136: \(-276\) \((2,2,35)\); 199: \(-591\); 203: \(-615\) |
| \(\mathbb{Q}(\sqrt{-163})\) | 161: \(-159\) (trivial); 165: \(-167\) (trivial); 324: \(-644\) \((2,2,81)\); 328: \(-660\) \((2,2,83)\) |

(the \(t = 2\) level pair of every field carries the norm-\(2\) twist \((2, 2, \cdot)\)
or \((2, 0, \cdot)\) — the "\(2\)-twist family" \(D = -2(2|d_K| \pm 4) = -4(|d_K| \pm 2)\)
resp. \(-8(2m \pm 1)\)).

**Coincidences and their classes** (R3; \(|D| \le 300\): \(61\) discriminants are realized
at primitive strata, \(18\) at two or three fields):

| \(D\) | \(\mathrm{Cl}(D)\) | realizations (field, \(k\): class) |
|---|---|---|
| \(-15\) | \(\mathbb{Z}/2\) | \(\mathbb{Q}(\sqrt{-3})\) 7: \((2,1,2)\); \(\mathbb{Q}(i)\) 8 (\(i\mathcal{S}\), \(n = 4\)): \((2,1,2)\); \(\mathbb{Q}(\sqrt{-11})\) 13, \(\mathbb{Q}(\sqrt{-19})\) 17: trivial |
| \(-20\) | \(\mathbb{Z}/2\) | \(\mathbb{Q}(\sqrt{-3})\) 8, \(\mathbb{Q}(\sqrt{-7})\) 12: \((2,2,3)\) |
| \(-24\) | \(\mathbb{Z}/2\) | \(\mathbb{Q}(i)\) 10, \(\mathbb{Q}(\sqrt{-2})\) 14: \((2,0,3)\) |
| \(-36\) | \(\mathbb{Z}/2\) | \(\mathbb{Q}(\sqrt{-7})\) 16, \(\mathbb{Q}(\sqrt{-11})\) 20: \((2,2,5)\) |
| \(-39, -63, -159\) | \(\mathbb{Z}/4, \mathbb{Z}/4, \mathbb{Z}/10\) | one field carries the nontrivial class (\(\mathbb{Q}(\sqrt{-3})\) 11, \(\mathbb{Q}(i)\) 16, \(\mathbb{Q}(\sqrt{-19})\) 55) and \(\mathbb{Q}(\sqrt{-43})\) 41, \(\mathbb{Q}(\sqrt{-67})\) 65, \(\mathbb{Q}(\sqrt{-163})\) 161 the trivial one |
| \(-84\) | \((\mathbb{Z}/2)^2\) | \(\mathbb{Q}(\sqrt{-3})\) 16: \((5,4,5)\); \(\mathbb{Q}(\sqrt{-2})\) 26: \((3,0,7)\); \(\mathbb{Q}(\sqrt{-19})\) 40: \((2,2,11)\) — the three nontrivial classes, one per field |
| \(-160\) | \((\mathbb{Z}/2)^2\) | \(\mathbb{Q}(\sqrt{-3})\) 22: \((7,6,7)\); \(\mathbb{Q}(\sqrt{-11})\) 42: \((4,4,11)\) |
| \(-180\) | \((\mathbb{Z}/2)^2\) | \(\mathbb{Q}(\sqrt{-2})\) 38: \((5,0,9)\); \(\mathbb{Q}(\sqrt{-43})\) 88: \((2,2,23)\) |
| \(-195\) | \((\mathbb{Z}/2)^2\) | \(\mathbb{Q}(i)\) 28 (\(i\mathcal{S}\), \(n = 14\)): \((7,1,7)\); \(\mathbb{Q}(\sqrt{-7})\) 37: \((5,5,11)\) |
| \(-224\) | \(\mathbb{Z}/4\times\mathbb{Z}/2\) | \(\mathbb{Q}(\sqrt{-3})\) 26: \((8,8,9)\); \(\mathbb{Q}(i)\) 30 (\(n = 15\)): \((7,0,8)\) |
| \(-255\) | \(\mathbb{Z}/6\times\mathbb{Z}/2\) | \(\mathbb{Q}(i)\) 32 (\(i\mathcal{S}\)): \((8,1,8)\); \(\mathbb{Q}(\sqrt{-11})\) 53: \((5,5,14)\) |
| \(-260\) | \(\mathbb{Z}/4\times\mathbb{Z}/2\) | \(\mathbb{Q}(\sqrt{-3})\) 28: \((9,8,9)\); \(\mathbb{Q}(\sqrt{-67})\) 132: \((2,2,33)\) |
| \(-276\) | \(\mathbb{Z}/4\times\mathbb{Z}/2\) | \(\mathbb{Q}(\sqrt{-7})\) 44: \((6,6,13)\); \(\mathbb{Q}(\sqrt{-67})\) 136: \((2,2,35)\) |
| \(-288\) | \((\mathbb{Z}/2)^2\) | \(\mathbb{Q}(i)\) 34 (\(n = 17\)): \((8,0,9)\); \(\mathbb{Q}(\sqrt{-19})\) 74: \((4,4,19)\) |

With \(|\mathrm{Cl}(D)[2]| = 2\) all nontrivial realizations of a \(D\) agree (forced); with
\(|\mathrm{Cl}(D)[2]| \ge 4\) different fields hit different classes, as Theorem 2(5)
predicts: at \(D = -84\) the pairs \(\{6, 14\}\), \(\{3, 7\}\), \(\{2, 42\}\) give
\(|d_K| = (14 + 4)/6 = 3\), \(m = (7 - 1)/3 = 2\), \(|d_K| = (42 - 4)/2 = 19\) — and no other
class-number-one field.

**Unrealized slices** (R4). Among the slices with nontrivial class and \(|D| \le 120\),
twelve are realized by no primitive stratum of the nine fields:
\(-56\ (2,0,7)\), \(-72\ (2,0,9)\), \(-88\ (2,0,11)\), \(-91\ (5,3,5)\), \(-96\ (3,0,8)\) and
\((5,2,5)\), \(-100\ (2,2,13)\), \(-104\ (2,0,13)\), \(-115\ (5,5,7)\), \(-116\ (2,2,15)\),
\(-120\ (2,0,15)\) and \((3,0,10)\). The eight discriminants
\(-56, -72, -88, -91, -100, -104, -115, -116\) have no realization at all (the list of the
session prompt); \(D = -96\) and \(-120\) are realized at one class only
(\(\mathbb{Q}(\sqrt{-7})\) \(k = 26\): \((4,4,7)\); \(\mathbb{Q}(i)\) \(n = 11\): \((5,0,6)\)). Their
carriers are the imprimitive strata (§3) and the principal cusps of the class-number-two
fields (below).

**The principal cusp at class number two** (X5, certified). Theorem 2 holds verbatim at
the principal cusp of every \(K\) (the orbit of \(\hat{\mathbb{R}}\)): at
\(d_K = -20, -24, -40, -52\) the class formula with Lemmas A′–C′ holds at every class of
the first four levels \(k = 4mt \pm 2\), and the \(t = 2\) levels carry new slices:

| \(K\) | \(t = 1\): \(k\), \(D\) (trivial twist) | \(t = 2\): \(k\), \(D\), twist |
|---|---|---|
| \(\mathbb{Q}(\sqrt{-5})\) | 18: \(-16\); 22: \(-24\) | 38: \(-72\), \((2,0,9)\); 42: \(-88\), \((2,0,11)\) |
| \(\mathbb{Q}(\sqrt{-6})\) | 22: \(-20\); 26: \(-28\) | 46: \(-88\), \((2,0,11)\); 50: \(-104\), \((2,0,13)\) |
| \(\mathbb{Q}(\sqrt{-10})\) | 38: \(-36\); 42: \(-44\) | 78: \(-152\), \((2,0,19)\); 82: \(-168\), \((2,0,21)\) |
| \(\mathbb{Q}(\sqrt{-13})\) | 50: \(-48\); 54: \(-56\) | 102: \(-200\), \((2,0,25)\); 106: \(-216\), \((2,0,27)\) |

and \(k = 58\) at \(\mathbb{Q}(\sqrt{-5})\) (\(t = 3\), \(s = +1\)) carries \(D = -168\) with the
twist \((14, 0, 3) \sim (3, 0, 14)\): with \(\mathbb{Q}(i)\) at \(n = 13\) (\((6,0,7)\)) and
\(\mathbb{Q}(\sqrt{-10})\) at \(k = 82\) (\((2,0,21)\)) all three nontrivial classes of
\(\mathrm{Cl}(-168) \cong (\mathbb{Z}/2)^2\) are realized. The residues
\(k^2 \equiv 4 \bmod |d_K|\) with \(k \not\equiv \pm2\) (\(k \equiv \pm8 \bmod 20\),
\(\pm10 \bmod 24\), \(\pm18 \bmod 40\), \(\pm24 \bmod 52\)) are **not** levels of the
principal cusp (Theorem 1(4) of other-fields.md forbids them); the levels
\(k = 58\) at \(\mathbb{Q}(\sqrt{-6})\) (\(D = -140\)) and \(k = 98\) at \(\mathbb{Q}(\sqrt{-10})\)
(\(D = -240\)) expected in the session prompt are of this kind, and remain candidates for
the non-principal cusp (§13).

## 3. Imprimitive strata: the Pell equation and the stratum twist class

At a level \(k\) with \(D = f^2D'\) the stratum of conductor \(f\) consists of the forms
\(f\cdot f'\), \(f'\) primitive of discriminant \(D'\); their ideals are the proper
\(\mathcal{O}_{D'}\)-ideals \(f\cdot\mathfrak{a}_{f'}\), and Paper II's imprimitive clause
(other-fields.md Theorem 4) twists them by \(\mathfrak{r}_\alpha\mathcal{O}_{D'}\).

> **Theorem 3 (strata and the Pell equation).** Let \(K\) be an imaginary quadratic field
> and \(D'\) a discriminant with \(|d_K||D'|\) not a square (i.e. \(K \ne \mathbb{Q}(\sqrt{D'})\)).
> 1. \(\mathcal{S}_K\) has a stratum of discriminant \(D'\) at the level \(k\) with conductor
>    \(f\) iff \(k^2 - |d_K||D'|f^2 = 4\); the solutions \((k_j, f_j)\), \(j \ge 1\), are the powers
>    \(\eta^j = (k_j + f_j\sqrt{m})/2\), \(m = |d_K||D'|\), of the fundamental norm-one unit
>    \(\eta > 1\) of the order of discriminant \(m\) (\(m \equiv 0, 1 \bmod 4\)) or \(4m\): infinitely
>    many levels.
> 2. The stratum carries the slice \((D', [\mathfrak{r}'_j])\) with
>    \(\mathfrak{r}_\alpha\mathcal{O}_{D'} = c_j\mathfrak{r}'_j\), \(c_j \in \mathbb{Z}_{>0}\) the content and
>    \(\mathfrak{r}'_j\) a primitive ambiguous invertible ideal of \(\mathcal{O}_{D'}\) of norm
>    \(r_0(k_j)/c_j^2\); the values of \(R\) on the stratum are the values of the slice
>    \((D', [\mathfrak{r}'_j])\), class by class.
> 3. *(The law of the stratum twist.)* For an odd prime \(p \mid D'\) and every \(j\):
>    \(p^{v_p(D')} \,\|\, N(\mathfrak{r}'_j)\) if \(p \mid k_j - 2\) (equivalently
>    \(\eta^j \equiv 1 \bmod \mathfrak{p}\), \(\mathfrak{p}\) the prime above \(p\) of \(\mathbb{Q}(\sqrt m)\)),
>    and \(p \nmid N(\mathfrak{r}'_j)\) if \(p \mid k_j + 2\). The side of each odd prime is the
>    same for all odd \(j\) and is the "\(k-2\) side" for all even \(j\). Consequently, for
>    odd \(D'\): \([\mathfrak{r}'_j] = c^j\) with \(c = [\mathfrak{r}'_1]\) the class of the pair whose
>    \(r_0\)-side is \(\prod_{p \mid k_1 - 2}p^{v_p(D')}\) — a homomorphism from the solution
>    group \(\langle\eta\rangle\) to \(\mathrm{Cl}(D')[2]\), trivial on even \(j\). For even \(D'\) the
>    same statement holds with the \(2\)-part of \(N(\mathfrak{r}'_j)\) determined as in
>    (P4) below (certified).
> 4. Hence the strata of \(\mathcal{S}_K\) of discriminant \(D'\) realize exactly two slices,
>    the trivial one and \((D', c(K, D'))\), \(c(K, D') = [\mathfrak{r}'_1]\) — and only one if
>    \(c(K, D') = 1\).

*Proof.* (1) \(D = -(k^2-4)/|d_K| = f^2D'\) is the displayed equation; \((k + f\sqrt m)/2\)
has norm \(1\); for \(m \equiv 0, 1 \bmod 4\) it lies in the order of discriminant \(m\)
(\(k \equiv f \bmod 2\) is forced by \(k^2 \equiv mf^2 \bmod 4\)), for \(m \equiv 2, 3 \bmod 4\)
both \(k, f\) are even and \((k + f\sqrt m)/2 \in \mathbb{Z}[\sqrt m]\); the norm-one units
\(> 1\) of a real quadratic order form the cyclic group \(\langle\eta\rangle\). (2) is
other-fields.md Theorem 4 (imprimitive clause): the twist on a proper \(\mathcal{O}_{D'}\)-ideal
is \(\mathfrak{r}\mathcal{O}_{D'}\), invertible ambiguous with square \(r_0\mathcal{O}_{D'}\) and norm
\(r_0\) (norms are local and \(\mathfrak{r}\) is locally principal), and \(R\) is invariant under
rational scaling of the twist (Theorem 1(6)); the classwise identity is the homogeneity
\(\Delta(f\mathfrak{a}) = f^{-12}\Delta(\mathfrak{a})\) on both terms. (3) Let \(p\) be odd,
\(p \mid D'\). Since \(\gcd(k-2, k+2) \mid 4\), \(p\) divides exactly one of \(k \pm 2\); mod
\(\mathfrak{p}\) one has \(\sqrt m \equiv 0\), so \(\eta \equiv k_1/2\) and \(\eta^j \equiv (k_1/2)^j\)
with \(k_1/2 \equiv \pm1\) (as \(k_1^2 \equiv 4 \bmod p\)); then \(k_j = \operatorname{Tr}(\eta^j)
\equiv 2(\pm1)^j\), which proves the parity statement. If \(p \mid k_j + 2\) then
\(p \nmid r_0\) and \(p \nmid N(\mathfrak{r}'_j)\). If \(p \mid k_j - 2\): \(v_p(k_j - 2) = v_p(|d_K|) + v_p(D') + 2v_p(f_j)\)
and \(g_-\) absorbs exactly \(p^{v_p(d_K)}\) (\(d_K\) is squarefree at odd \(p\)), so
\(v_p(r_0) = v_p(D') + 2\varphi\), \(\varphi := v_p(f_j)\). The content \(c_j\) satisfies
\(v_p(c_j) \ge \varphi\): \(\mathfrak{r} = \mathbb{Z}r_0 + \mathbb{Z}\tfrac{b + f_j\sqrt{D'}}{2}\) with
\(b \in \{0, r_0\}\) lies in \(p^\varphi\mathcal{O}_{D'}\). And \(v_p(c_j) \le \varphi\): an element
\(\tfrac{x + y\sqrt{D'}}{2}\) of \(p^{\varphi+1}\mathcal{O}_{D'}\) has \(p^{\varphi+1} \mid y\) (\(p\) odd),
whereas \(\tfrac{b + f_j\sqrt{D'}}{2} \in \mathfrak{r}\mathcal{O}_{D'}\) has \(y = f_j\) with
\(v_p(f_j) = \varphi\). So \(v_p(N\mathfrak{r}'_j) = v_p(r_0) - 2\varphi = v_p(D')\). (4) follows
from (3): the class depends only on the parity of \(j\), and at even \(j\) the pair is
\(\{N(\theta'), 1\}\) for odd \(D'\) (all odd primes on the \(r_0\)-side); the certified
\(2\)-part completes the even case. \(\blacksquare\)

**Certified** (phase P; the twist extension \(\mathfrak{r}\mathcal{O}_{D'}\), its content, the
primitive part and its class computed by exact HNF arithmetic; fundamental units from
PARI's `quadunit`): (P1) the ten instances predicted while the prompt was written are
confirmed — \(\mathbb{Q}(\sqrt{-3})\) \(k = 26\) (\(D = -224 = 2^2\cdot(-56)\)): \((8,8,9)\) extends
to \(2\cdot(2,0,7)\); \(\mathbb{Q}(\sqrt{-2})\) \(k = 22\) (\(D = -60\)): \((5,0,3)\) extends to
\((5,5,2) \sim (2,1,2)\); \(\mathbb{Q}(i)\) \(n = 9\): \((4,0,5)\) extends to \(2\cdot(1,0,5)\)
(trivial, \(R \equiv 1\)); \(n = 17\), \(f = 2\): \((2,0,9)\); \(\mathbb{Q}(\sqrt{-19})\) \(k = 74\), \(f = 2\):
trivial; \(\mathbb{Q}(\sqrt{-7})\) \(k = 54\), \(f = 2\): \((2,0,13)\); \(\mathbb{Q}(i)\) \(n = 51\), \(f = 5\):
trivial; \(\mathbb{Q}(\sqrt{-3})\) \(k = 52\), \(f = 3\): \((2,2,13)\); \(\mathbb{Q}(\sqrt{-2})\) \(k = 198\),
\(f = 7\): trivial; \(\mathbb{Q}(i)\) \(n = 15\), \(f = 2\): \((2,0,7)\). (P2–P4) the law at
\(14 \times 9\) pairs \((D', K)\), \(D' \in \{-15, -20, -24, -36, -56, -72, -84, -88, -91, -100, -104, -115, -116, -160\}\),
four solutions each (\(k_j\) up to \(133\) digits): the class is \(c^j\) at every pair; the
odd part of \(N(\mathfrak{r}'_j)\) is \(\prod_{p \text{ odd},\ v_p(k-2) > v_p(k+2)}p^{v_p(D')}\)
always, and on the pairs of the first type \((r'_0, 0, s'_0)\) the same rule holds at
\(p = 2\) with the exponent \(v_2(|D'|/4)\) (for the second type at even \(D'\) the \(2\)-part
splits between the two sides, e.g. \(\{2, 58\}\) at \(D' = -116\), \(\{8, 20\}\) at \(-160\)).
(P5) The odd-\(j\) class \(c(K, D')\):

| \(D'\) | \(\mathrm{Cl}[2]\) nontrivial classes | \(c(K, D')\) for \(K = \mathbb{Q}(\sqrt{-3}), \mathbb{Q}(i), \mathbb{Q}(\sqrt{-7}), \mathbb{Q}(\sqrt{-2}), \mathbb{Q}(\sqrt{-11}), \mathbb{Q}(\sqrt{-19}), \mathbb{Q}(\sqrt{-43}), \mathbb{Q}(\sqrt{-67}), \mathbb{Q}(\sqrt{-163})\) |
|---|---|---|
| \(-56\) | \((2,0,7)\) | nontrivial at all nine (\(f_1 = 2, 2, 10, 12, 1716, 42, 28, 560, 74740\)) |
| \(-72\) | \((2,0,9)\) | \(A, A, A, 1, A, 1, 1, A, 1\) (\(A = (2,0,9)\); \(\mathbb{Q}(i)\): \(n = 17\), \(f = 2\)) |
| \(-84\) | \((5,4,5), (3,0,7), (2,2,11)\) | \((5,4,5), (3,0,7), (3,0,7), (3,0,7), (5,4,5), (2,2,11), (3,0,7), (3,0,7), (3,0,7)\) |
| \(-88\) | \((2,0,11)\) | nontrivial except at \(\mathbb{Q}(\sqrt{-7})\) |
| \(-91\) | \((5,3,5)\) | trivial at \(\mathbb{Q}(\sqrt{-3}), \mathbb{Q}(\sqrt{-43})\), nontrivial at the other seven (\(\mathbb{Q}(\sqrt{-2})\): \(k = 54\), \(f = 2\)) |
| \(-100\) | \((2,2,13)\) | nontrivial at \(\mathbb{Q}(\sqrt{-3}), \mathbb{Q}(\sqrt{-7}), \mathbb{Q}(\sqrt{-43}), \mathbb{Q}(\sqrt{-67}), \mathbb{Q}(\sqrt{-163})\), trivial at the rest |
| \(-104\) | \((2,0,13)\) | trivial at \(\mathbb{Q}(\sqrt{-3}), \mathbb{Q}(i), \mathbb{Q}(\sqrt{-43})\), nontrivial at the rest |
| \(-115\) | \((5,5,7)\) | trivial at \(\mathbb{Q}(\sqrt{-11}), \mathbb{Q}(\sqrt{-19})\), nontrivial at the rest |
| \(-116\) | \((2,2,15)\) | nontrivial at \(\mathbb{Q}(\sqrt{-3}), \mathbb{Q}(\sqrt{-11}), \mathbb{Q}(\sqrt{-19}), \mathbb{Q}(\sqrt{-43}), \mathbb{Q}(\sqrt{-163})\), trivial at the rest |
| \(-160\) | \((7,6,7), (4,4,11), (5,0,8)\) | \((7,6,7), 1, (7,6,7), 1, (4,4,11), (4,4,11), (7,6,7), (7,6,7), (7,6,7)\) — **\((5,0,8)\) at none** |

So every nontrivial slice with \(|D'| \le 120\) is carried by an imprimitive stratum of
some class-number-one field (and \((-56, [(2,0,7)])\) by all nine), but not every slice by
every field: \((D', [\mathfrak{r}'])\) is reached by the strata of \(\mathcal{S}_K\) iff
\([\mathfrak{r}'] = c(K, D')\), and the class \((5, 0, 8)\) of \(D' = -160\) is reached by no
stratum of any of the nine fields (nor by their primitive strata: its pair \(\{5, 8\}\)
satisfies none of the congruences of Theorem 2(5)). The unit system of that slice is
nevertheless computed in §8–§9 — the slice theory is larger than the arrangements.
**Exact identification** (U2, \(200\) digits): on the strata of the ten instances the
values \(R\) coincide class by class with the values of the slice \((D', [\mathfrak{r}'])\)
computed directly (spare digits \(\ge 200\)), and the stratum polynomials agree with the
primitive-stratum polynomials of the other realizations: \(x^2 - 7x + 1\) at
\(D' = -15\) (\(\mathbb{Q}(\sqrt{-3})\) \(k = 7\), \(i\mathcal{S}\) \(n = 4\), the \(f = 2\) stratum of
\(\mathbb{Q}(\sqrt{-2})\) \(k = 22\)); \(x^4 - 1988x^3 - 3194x^2 - 1988x + 1\) at \(D' = -56\) (the
\(f = 2\) strata of \(\mathbb{Q}(\sqrt{-3})\) \(k = 26\) and \(\mathbb{Q}(i)\) \(n = 15\), matching Paper II's
imprimitive record); \(x^2 - 9602x + 1\) at \((-72, (2,0,9))\); \(x^2 + 103682x + 1\) at
\((-100, (2,2,13))\); \((x-1)^2\) on the trivial stratum of \(n = 9\).

## 4. Units on slices

**Theorem 4 (units on a slice; = other-fields.md Theorem 4).** For every slice
\((D, [\mathfrak{r}])\): the \(R_\mathfrak{b}\) are units of \(H_D\), \(R_{\mathfrak{r}\mathfrak{b}} = R_\mathfrak{b}^{-1}\),
\(R_{\bar{\mathfrak{b}}} = \overline{R_\mathfrak{b}}\), the slice polynomial \(\prod_\mathfrak{b}(x - R_\mathfrak{b})\) is
a palindromic integer polynomial with constant term \(1\), and
\(\sum_\mathfrak{b}\chi(\mathfrak{b})\log|R_\mathfrak{b}| = -24L'(0,\chi)\) on the odd characters
(\(\chi(\mathfrak{r}) = -1\)), \(0\) on the even ones. (Proved: the proof never used an
arrangement. Certified here, phase U, \(200\) digits: at \(54\) slices — the first
levels of the nine fields, even levels of \(i\mathcal{S}\) included — with the twist and
conjugation laws, the independence of the partner \(\mathfrak{s}\), and the limit formula
against the independent incomplete-gamma evaluation at every character, spare
\(\ge 169\).)

**Equal slices, equal systems; different classes, different systems** (U1). The
polynomials of the \(54\) slices agree whenever the slices agree — \(D = -15\), class
\((2,1,2)\): \(x^2 - 7x + 1\) at \(\mathbb{Q}(\sqrt{-3})\) \(k = 7\), \(i\mathcal{S}\) \(n = 4\) (and the
\(f = 2\) stratum of \(\mathbb{Q}(\sqrt{-2})\) \(k = 22\)); \(-20\): \(x^2 + 18x + 1\) at
\(\mathbb{Q}(\sqrt{-3})\) \(k = 8\), \(\mathbb{Q}(\sqrt{-7})\) \(k = 12\); \(-24\): \(x^2 - 34x + 1\) at \(\mathbb{Q}(i)\)
\(n = 5\), \(\mathbb{Q}(\sqrt{-2})\) \(k = 14\); \(-36\): \(x^2 + 194x + 1\) at \(\mathbb{Q}(\sqrt{-7})\) \(k = 16\),
\(\mathbb{Q}(\sqrt{-11})\) \(k = 20\) — and differ, with different sets of odd characters,
whenever the classes differ:

| \(D\) | slice | realization | \(\prod(x - R)\) |
|---|---|---|---|
| \(-84\) | \((5,4,5)\) | \(\mathbb{Q}(\sqrt{-3})\) \(k = 16\) | \(x^4 + 686308x^3 + 7365318x^2 + 686308x + 1\) |
| \(-84\) | \((3,0,7)\) | \(\mathbb{Q}(\sqrt{-2})\) \(\alpha = 13\) | \(x^4 - 297220x^3 + 7312902x^2 - 297220x + 1\) |
| \(-84\) | \((2,2,11)\) | \(\mathbb{Q}(\sqrt{-19})\) \(\alpha = 20\) | \(x^4 + 27940x^3 + 76614x^2 + 27940x + 1\) |
| \(-160\) | \((7,6,7)\) | \(\mathbb{Q}(\sqrt{-3})\) \(k = 22\) | \(x^4 + 5677418956x^3 + 13748395686x^2 + \cdots\) |
| \(-160\) | \((4,4,11)\) | \(\mathbb{Q}(\sqrt{-11})\) \(k = 42\) | \(x^4 + 2146623116x^3 + 4535235366x^2 + \cdots\) |
| \(-160\) | \((5,0,8)\) | none | \(x^4 - 4064541764x^3 + 12286753926x^2 - \cdots\) |
| \(-168\) | \((6,0,7)\) | \(\mathbb{Q}(i)\) \(n = 13\) | Paper II: \(x^4 - 11382984004x^3 + 885435408006x^2 - \cdots\) |
| \(-168\) | \((3,0,14)\) | \(\mathbb{Q}(\sqrt{-5})\) \(k = 58\) | \(x^4 - 844926404x^3 + 885289852806x^2 - \cdots\) |
| \(-180\) | \((5,0,9)\), \((2,2,23)\) | \(\mathbb{Q}(\sqrt{-2})\) 38, \(\mathbb{Q}(\sqrt{-43})\) 88 | \(x^4 - 28049681284x^3 + \cdots\), \(x^4 + 22199076x^3 + 48146246x^2 + \cdots\) |
| \(-195\) | \((7,1,7)\), \((5,5,11)\) | \(i\mathcal{S}\) 14, \(\mathbb{Q}(\sqrt{-7})\) 37 | \(x^4 - 174683640964x^3 + 2849307736326x^2 - \cdots\), \(x^4 - 112143790724x^3 + 2842988279046x^2 - \cdots\) |

The three \(D = -84\) systems have odd characters \(\{\chi_{(0,1)}, \chi_{(1,0)}\}\),
\(\{\chi_{(0,1)}, \chi_{(1,1)}\}\), \(\{\chi_{(1,0)}, \chi_{(1,1)}\}\) of
\(\mathrm{Cl}(-84) \cong (\mathbb{Z}/2)^2\): each nontrivial class \([\mathfrak{r}]\) is odd for exactly the
two characters not killing it, so the three unit systems are the three "halves" of the
same character theory. The phases \(\Theta = u/\varepsilon\) of these realizations
differ as well: the Cartan unit \(\varepsilon_K(\alpha) = \alpha + \sqrt{\alpha^2 - 1}\) lies in
\(\mathbb{Q}(\sqrt{|d_K||D|})\) — \(\mathbb{Q}(\sqrt7)\), \(\mathbb{Q}(\sqrt{42})\), \(\mathbb{Q}(\sqrt{399})\) at the three
realizations of \(D = -84\) — while \(u\) itself is slice-intrinsic (§9).

**Genus coincidences across discriminants** (observed, explained by §7). The
quadratic slices \(D = -35\) (\(i\mathcal{S}\), \(n = 6\)) and \(D = -40\) (\(\mathbb{Q}(\sqrt{-2})\), \(k = 18\))
have the same polynomial \(x^2 - 322x + 1\), i.e. \(R = \varepsilon_5^{\pm12}\) in both
ring class fields, which both contain \(\mathbb{Q}(\sqrt5)\); and \(x^2 - 7x + 1\) (\(D = -15\)),
\(x^2 + 18x + 1\) (\(D = -20\)), \(x^2 - 34x + 1\) (\(D = -24\)) are \(\varepsilon_5^{\pm4}\),
\(-\varepsilon_5^{\pm6}\), \(\varepsilon_8^{\pm4}\): by the closed form of §7,
\(\log|R_1| = -12L'(0,\chi) = -12\tfrac{2h(d_1)}{w(d_1)}h(d_2)\log\varepsilon_{d_2}\cdot C(0)\) at
an \(h = 2\) slice with \(D = d_1d_2\), so slices with the same \((w(d_1), h(d_1), d_2)\) have
the same unit up to sign (\(D = -35 = (-7)\cdot5\) and \(-40 = (-8)\cdot5\): \(w = 2\), \(h = 1\)).

## 5. The even levels of \(i\mathcal{S}\)

\(i\mathcal{S} = \mathrm{PGL}_2(\mathbb{Z}[i])\cdot\hat{\mathbb{R}} \setminus \mathcal{S}\) is the orbit of the
imaginary axis, \(i\mathcal{S} = \mathrm{PSL}_2(\mathbb{Z}[i])\cdot i\hat{\mathbb{R}}\), and also the set of
circles \(X(\hat{\mathbb{R}})\), \(X \in \mathrm{GL}_2(\mathbb{Z}[i])\), \(\det X = i\). Its circles are the
Hermitian matrices \(\binom{2a\ \ i\beta}{-i\bar\beta\ \ 2c}\) with \(\beta \equiv i \bmod 2\) — centre
\((x + ni)/(2a)\) with \(x\) odd, \(n\) even, \(x^2 + n^2 \equiv 1 \bmod 4a\)
([circle-classification.md](circle-classification.md), Remark) — so its levels are the even
\(n \ge 2\), with both orientations, and the level-\(n\) circles are the forms
\((a, b, c)\) of discriminant \(1 - n^2\) (odd), \(b = -x\), imprimitive included, with the
census \(3H(n^2 - 1)\) ([hyperbolic-counting.md](hyperbolic-counting.md)). The involution
\(\sigma(X) = \bar X^{-1}\) preserves the coset \(\det X = i\) (\(\det\bar X^{-1} = 1/\overline{\det X} = i\))
and the level (\(-\tfrac12\operatorname{tr}(X\bar X^{-1})\) is conjugation invariant).

> **Theorem 5 (class formula of \(i\mathcal{S}\)).** Let \(n \ge 4\) be even, \(D = 1 - n^2\),
> \(f = (a, b, c)\) a primitive form of discriminant \(D\) with circle of orientation sign
> \(s\) (\(\beta = -s(n + bi)\), \(q = sa\)). The circle of \(\sigma(X)\) has negative orientation,
> level \(n\), lies in the lower half-plane iff \(s = +1\), and after reflection its class is
> $$\hat\sigma[f] \;=\; [\mathfrak{r}_n]\,[f]^{-s},\qquad
> \mathfrak{r}_n = \bigl(n-1,\ n-1,\ \tfrac n2\bigr) \sim \bigl(\tfrac n2,\ 1,\ \tfrac n2\bigr),$$
> the ambiguous class of the pair \(\{n-1, n+1\}\), \(\mathfrak{r}_n\mathfrak{s}_n = (\sqrt D)\).

The proof is Paper I's (Lemmas A–C of [class-formula-proof.md](class-formula-proof.md),
in the general form of other-fields.md) on the coset \(\det X = i\), with the odd-\(D\)
Gram form. Fix \(f\), \(\beta\), and put \(\delta' := -i\).

**Lemma A′′.** \(\mathcal{K}_f := \{u \in \mathbb{Z}[i] : \delta'\bar u - \beta u \in 2a\mathbb{Z}[i]\}\) has index
\(a\) in \(\mathbb{Z}[i]\) and contains \(a\mathbb{Z}[i]\); for a basis \((u_1, u_2)\) with
\(\operatorname{Im}(u_1\bar u_2) = sa\), \(v_k := s(\delta'\bar u_k - \beta u_k)/(2ia) \in \mathbb{Z}[i]\),
\(P = \binom{u_1\ v_1}{u_2\ v_2}\) has \(\det P = \delta'\) and \(P^\dagger M_0P = M_f\); \(X := P^{-1}\)
(\(\det X = i\)) realizes the circle.

*Proof.* Solving \(u_1v_2 - u_2v_1 = \delta'\), \(\bar u_1v_2 - \bar u_2v_1 = \beta\) (the two
Gram conditions \(\det P = \delta'\), \(h_0(p_1, p_2) = i\beta\)) gives the displayed \(v_k\)
(the determinant of the system is \(2i\operatorname{Im}(u_1\bar u_2) = 2isa\)), and
\(h_0(p_2,p_2)\) is then forced by \(\det(P^\dagger M_0P) = -|\det P|^2 = -1\). The map
\(u \mapsto \delta'\bar u - \beta u\) has real determinant \(N(\beta) - |\delta'|^2 = 4ac\) and
is \(\equiv 0 \bmod 2\) (\(\delta' \equiv \beta \equiv i\), \(\bar u \equiv u \bmod 2\)), so
\(\psi := \tfrac12(\delta'\bar u - \beta u)\) is an integral map of determinant \(ac\) with
\(\mathcal{K}_f = \psi^{-1}(a\mathbb{Z}[i])\); in the coordinates \(u = x + iy\) its matrix is
\(\binom{sn/2\ \ -(sb+1)/2}{(sb-1)/2\ \ sn/2}\), of content \(1\) (a prime dividing \(n/2\) and both
\((b\pm1)/2\) would divide \(1\)), so the Smith invariants are \((1, ac)\) and the index of
the preimage of \(a\mathbb{Z}[i]\) is \(\gcd(1, a)\gcd(ac, a) = a\). \(\square\)

**Lemma B′′.** With \(g_s(u) := nN(u) - s\operatorname{Im}(u^2)\) (\(= n(x^2+y^2) - 2sxy\), the form
\(2\cdot(\tfrac n2, -s, \tfrac n2)\) of discriminant \(D\) in \((x, y)\)), \(N = X^\dagger M_0X\) has
\(N_{11} = g_s(u_2)/a\), \(N_{22} = g_s(u_1)/a\), \(\operatorname{Re}\beta' = \operatorname{Re}\beta\); the
circle of \(\sigma(X)\) is \(-\bar N\), of negative orientation, below \(\hat{\mathbb{R}}\) iff
\(s = +1\) (exactly as in Theorem 3 of other-fields.md), and after the reflection when
\(s = +1\) the image form is
\(\tfrac{1}{2a}\bigl(g_s(u_2), -2g_s(u_1,u_2), g_s(u_1)\bigr) = \tfrac{1}{2a}g_s|_{\mathcal{K}_f}\) in the
basis \((u_2, -u_1)\).

*Proof.* \(X = \delta'^{-1}\operatorname{adj}P\) and \(|\delta'| = 1\), so \(N = (\operatorname{adj}P)^\dagger M_0(\operatorname{adj}P)\)
as in Lemma B′, and \(2\operatorname{Im}(u_k\bar v_k) = \tfrac{s}{a}\operatorname{Re}(\bar\delta'u_k^2 - \bar\beta|u_k|^2)
= \tfrac1a(nN(u_k) + s\operatorname{Re}(iu_k^2)) = g_s(u_k)/a\), using \(\operatorname{Re}\beta = -sn\). \(\square\)

**Lemma C′′.** \(\iota(x + iy) := x\tfrac n2 + y\tfrac{-s + \sqrt D}{2}\) satisfies
\(N(\iota(u)) = \tfrac n4g_s(u)\), \(\mathfrak{t} := \iota(\mathbb{Z}[i]) = [\tfrac n2, \tfrac{-s+\sqrt D}{2}]\) is
the ideal of the form \((\tfrac n2, s, \tfrac n2)\), \(\iota\) intertwines multiplication by
\(\sqrt D\) with the integral map \(x + iy \mapsto (sx - ny) + i(nx - sy)\), and for
primitive \(f\), \(\iota(\mathcal{K}_f) = \mathfrak{t}\mathfrak{a}_f\), \(\mathfrak{a}_f = [a, \tfrac{b+\sqrt D}{2}]\).

*Proof.* The norm identity is \((\tfrac n2x - \tfrac s2y)^2 - \tfrac D4y^2 = \tfrac n4(nx^2 - 2sxy + ny^2)\)
(\(1 - D = n^2\)); the intertwining is the computation of §2 of other-fields.md with
\(\varphi^2 = D\); \(\mathcal{K}_f\) is stable under this map (the two divisibility conditions are
preserved — checked as in Lemma C′), so \(\iota(\mathcal{K}_f)\) is an \(\mathcal{O}_D\)-module;
\(a\mathfrak{t} = \iota(a\mathbb{Z}[i]) \subseteq \iota(\mathcal{K}_f)\) and \(\tfrac{b+\sqrt D}{2}\mathfrak{t} \subseteq \iota(\mathcal{K}_f)\)
by the same coordinate computation; both sides have index \(a\) in \(\mathfrak{t}\). \(\square\)

*Proof of Theorem 5.* By Lemmas B′′–C′′ the image form is the norm form of the proper
ideal \(\iota(\mathcal{K}_f) = \mathfrak{t}\mathfrak{a}_f\) in the basis \((\iota(u_2), \iota(-u_1))\), of
orientation \(-s\), so \([f'] = [\mathfrak{t}\mathfrak{a}_f]^{-s} = [\mathfrak{t}][f]^{-s}\); and
\([\mathfrak{t}] = [(\tfrac n2, s, \tfrac n2)] = [(\tfrac n2, 1, \tfrac n2)] = [(n-1, n-1, \tfrac n2)]\):
the reduction of \((n-1, n-1, \tfrac n2)\) is \((\tfrac n2, 1, \tfrac n2)\) (swap, then
\(b \mapsto b + n\)), a \(2\)-torsion class. \(\blacksquare\)

**Certified** (phase E): at every primitive class and both orientations of
\(n = 4, 6, \dots, 16\) (\(72\) pairs (class, orientation)): the explicit \(P\) of Lemma A′′
(index \(a\), \(\det P = -i\), \(P^\dagger M_0P = M_f\), \(X = P^{-1}\) of determinant \(i\) realizing
the circle), the Gram identities of Lemma B′′, the identities of Lemma C′′
(\(N\circ\iota = \tfrac n4g_s\), \(\iota(\mathcal{K}_f) = \mathfrak{t}\mathfrak{a}_f\) as exact lattices, \([\mathfrak{t}] = [\mathfrak{r}_n]\)),
the class of \(\sigma(X)(\hat{\mathbb{R}})\) equal to \([\mathfrak{r}_n][f]^{-s}\), reflection exactly when
\(s = +1\); and the BFS orbit of \(i\hat{\mathbb{R}}\) equals the congruence set of the Remark of
circle-classification.md for every curvature \(\le 6\). The slices of \(i\mathcal{S}\) with their
unit polynomials (phase U) and odd index (phase I, §8):

| \(n\) | \(D\) | \(\mathrm{Cl}(D)\) | \(\mathfrak{r}_n\) | \(\prod_f(x - R_f)\) | \([E^-:\langle R_f\rangle]\) |
|---|---|---|---|---|---|
| 4 | \(-15\) | \(\mathbb{Z}/2\) | \((2,1,2)\) | \(x^2 - 7x + 1\) | \(4 = 24/6\) |
| 6 | \(-35\) | \(\mathbb{Z}/2\) | \((3,1,3)\) | \(x^2 - 322x + 1\) | \(12 = 24/2\) |
| 8 | \(-63 = 3^2(-7)\) | \(\mathbb{Z}/4\) | \((4,1,4)\) | \(x^4 - 31279x^3 + 45681x^2 - 31279x + 1\) | \(576 = 24^2\) |
| 10 | \(-99 = 3^2(-11)\) | \(\mathbb{Z}/2\) | \((5,1,5)\) | \(x^2 - 4468994x + 1\) | \(8 = 24/3\) |
| 12 | \(-143\) | \(\mathbb{Z}/10\) | \((6,1,6)\) | \(x^{10} - 813095735x^9 + 41713983345x^8 - \cdots\) | — |
| 14 | \(-195\) | \((\mathbb{Z}/2)^2\) | \((7,1,7)\) | \(x^4 - 174683640964x^3 + 2849307736326x^2 - \cdots\) | \(576 = 24^2\) |
| 16 | \(-255\) | \(\mathbb{Z}/6\times\mathbb{Z}/2\) | \((8,1,8)\) | (class map only) | — |

The prediction \(x^2 - 7x + 1\) at \(n = 4\) (from the coincidence with
\(\mathbb{Q}(\sqrt{-3})\), \(k = 7\)) is confirmed; \(n = 6\) coincides with the \(\mathbb{Q}(\sqrt{-2})\)
slice \(D = -40\) as a unit (§4), not as a slice. This closes the ledger item "even
levels / \(i\mathcal{S}\)": the whole \(\alpha\)-theory of Papers I–II (classification, census,
class formula, unit theorem, limit formula, odd index, phase) runs at the even levels
with the twist \((n-1, n-1, \tfrac n2)\) of norm \(n - 1\) in place of \((\tfrac{n-1}2, 0, \tfrac{n+1}2)\).

## 6. The second-kind orbit \(\mathcal{S}_K^\perp\)

For any \(K\), \(\sqrt{d_K}\hat{\mathbb{R}} = i\hat{\mathbb{R}}\), and the **second-kind orbit** is
\(\mathcal{S}_K^\perp := \mathrm{PSL}_2(\mathcal{O}_K)\cdot i\hat{\mathbb{R}}\) (\(= i\mathcal{S}\) at \(\mathbb{Q}(i)\)). The
stabilizer of \(i\hat{\mathbb{R}}\) in \(\mathrm{SL}_2(\mathcal{O}_K)\) is the group of matrices
\(\binom{a\ \ \beta}{\gamma\ \ d}\) with \(a, d \in \mathbb{Z}\), \(\beta, \gamma \in \mathcal{O}_K \cap i\mathbb{R}\), conjugate
to \(\Gamma_0(|d_K|)\) (odd \(d_K\)) or \(\Gamma_0(m)\) (\(d_K = -4m\)). The circle
\(X(i\hat{\mathbb{R}})\), \(X = \binom{a\ b}{c\ d}\), has the Hermitian matrix
\(\binom{A\ \ B}{\bar B\ \ C}\), \(A = -\operatorname{Tr}(c\bar d)\), \(B = a\bar d + b\bar c\), \(C = -\operatorname{Tr}(a\bar b)\),
\(N(B) = 1 + AC\): the orbit of the second unimodular Hermitian form
\(h_1 = \operatorname{Tr}(x\bar y)\) (\(M_1 = \binom{0\ 1}{1\ 0}\)) instead of \(h_0 = 2\operatorname{Im}(x\bar y)\) (\(M_0\));
in the normalization of \(\mathcal{S}_K\) it is \(\beta = iB\), \(\sqrt{|d_K|}q = -A\). At
\(\mathbb{Q}(i)\) the two forms are isometric (\(\operatorname{diag}(i,1)\)), which is why
\(i\mathcal{S}\) is a rotated \(\mathcal{S}\); for \(K \ne \mathbb{Q}(i)\) they are not.

> **Theorem 6 (\(\mathcal{S}_K^\perp\)).** 1. *(Necessity, every \(K\).)* Every circle of
> \(\mathcal{S}_K^\perp\) has \(A, C \in \operatorname{Tr}(\mathcal{O}_K)\) (\(= \mathbb{Z}\) for odd \(d_K\), \(2\mathbb{Z}\) for
> even), \(B \in \mathcal{O}_K\), \(N(B) = 1 + AC\), and at \(\mathbb{Q}(i)\) moreover \(B \equiv 1 \bmod 2\).
> Its lines are \(\operatorname{Re}(\bar u z) = C/2\), \(u \in \mathcal{O}_K^\times\).
> 2. *(Sufficiency, the five Euclidean fields.)* Conversely every such matrix is in the
> orbit: the descent (translate \(B \mapsto B - \lambda A\), invert \((A,B,C)\mapsto(C,-\bar B,A)\))
> terminates at a line or at the unit circle \(\pm(1, 0, -1)\), which lies in the orbit for
> odd \(d_K\) (\(X_0 = \binom{-1\ -\bar\omega}{1\ -\omega}\)).
> 3. *(Levels and forms.)* Both orientations occur at every level. The level of a circle in
> \(\mathbb{H}\) is \(\alpha = y\sqrt{|d_K|}/2\), \(y = -\operatorname{sgn}(A)\operatorname{im}_K(B) \in \mathbb{Z}_{>0}\)
> (irrational for \(K \ne \mathbb{Q}(i)\)), its hyperbolic centre is the CM point of the form
> \((A, 2\operatorname{Re}B, C)\) of discriminant \(4 - y^2|d_K|\) (odd \(d_K\)), resp. of
> \(\tfrac12(A, 2\operatorname{Re}B, C)\) of discriminant \(1 - my^2\) with \(y\) even (\(d_K = -4m\)), and
> every form of that discriminant occurs, with both orientations: the level-\(y\)
> circles in the ideal triangle number \(3H(y^2|d_K| - 4)\) resp. \(3H(my^2 - 1)\).
> 4. *(No twist.)* The adjoint involution of \(h_1\), \(\tau(X) = M_1X^\dagger M_1 = E\bar X^{-1}E\)
> (\(E = \operatorname{diag}(1,-1)\)), is well defined on circles of \(\mathcal{S}_K^\perp\) and preserves
> the level with respect to the imaginary axis, but not the level \(y\) with respect to
> \(\hat{\mathbb{R}}\); for \(K \ne \mathbb{Q}(i)\) there is no level-preserving involution of the type
> \(X \mapsto J\bar X^{-1}J'\) with \(J, J' \in \mathrm{GL}_2(\mathcal{O}_K)\), and the circles of
> \(\mathcal{S}_K^\perp\) are Heegner points without an Atkin–Lehner partner:
> \(\mathcal{S}_K^\perp\) realizes discriminants, not slices.

*Proof.* (1) The formulas for \(A, B, C\) are the entries of \((X^{-1})^\dagger M_1X^{-1}\);
\(N(B) - AC = |ad - bc|^2 = 1\); at \(\mathbb{Q}(i)\), \(B \equiv ad + bc = 1 + 2bc \equiv 1 \bmod 2\)
since \(\bar x \equiv x \bmod 2\). Lines: \(A = 0\) forces \(N(B) = 1\). (2) After the
translation \(N(B) \le A^2\rho_K^2\) with \(\rho_K^2 < 1\), so \(|C| = (N(B)-1)/|A| < |A|\)
unless \(N(B) = 0\), \(|A| = 1\): the unit circle. Terminal lines: \(B = u^2\) is realized
by \(\binom{u\ b}{0\ \bar u}\), \(B = -u^2\) by \(\binom{b\ -u}{\bar u\ 0}\) (every unit is \(\pm\) a
square), with \(C = \mp\operatorname{Tr}(u\bar b)\) arbitrary in \(\operatorname{Tr}(\mathcal{O}_K)\). (3) The
inversion maps \((0, 1, 0)\) to \((0, -1, 0)\), so \(-M_1\) is in the orbit and every circle
occurs with both orientations. The centre is \(-B/A\), the radius \(1/|A|\), so
\(\alpha = -\operatorname{sgn}(A)\operatorname{Im}B\) and \(\operatorname{Im}B \in \tfrac{\sqrt{|d_K|}}{2}\mathbb{Z}\); the
hyperbolic centre \((-2\operatorname{Re}B + i\sqrt{y^2|d_K| - 4})/(2A)\) is the root of
\(Az^2 + 2\operatorname{Re}(B)z + C\) (discriminant \(4(\operatorname{Re}B)^2 - 4(N(B)-1) = 4 - 4(\operatorname{Im}B)^2\)).
For even \(d_K\), \(A, C, 2\operatorname{Re}B\) are even and \(N(B) \equiv 1 \bmod 4\) forces \(y\) even.
Given a form of the right discriminant, \(B\) is reconstructed with
\(2\operatorname{Re}B = b\), \(\operatorname{im}_K B = \mp y\) (the parity condition \(b \equiv y \bmod 2\) holds
because \(D \equiv y^2 \bmod 4\)), and (2) realizes the circle. (4) For \(\gamma\) in the
stabilizer, \(\bar\gamma = E\gamma E\), so \(\tau(X\gamma) = \tau(X)\gamma^{-1}\) and the circle
\(\tau(X)(i\hat{\mathbb{R}})\) is well defined; it is \(-\overline{X^{-1}(i\hat{\mathbb{R}})}\), the reflection in
\(i\hat{\mathbb{R}}\) of the circle of \(X^{-1}\) — the exact analogue of
\(\sigma(X)(\hat{\mathbb{R}}) = \overline{X^{-1}(\hat{\mathbb{R}})}\) in \(\mathcal{S}_K\). The level with respect to
\(\hat{\mathbb{R}}\) is \(\operatorname{tr}(XE\bar X^{-1})/2i\), the level with respect to \(i\hat{\mathbb{R}}\) is
\(-\tfrac12\operatorname{tr}(XE\bar X^{-1}E)\); the second is preserved by \(\tau\) and by
\(X \mapsto X^{-1}\), the first is not (\(\operatorname{tr}(EX\bar X^{-1}) \ne \operatorname{tr}(E\bar X^{-1}X)\) in
general). A level-preserving \(J\bar X^{-1}J'\) would need \(\bar J^{-1}J = -E\) and
\(J'E\bar J'^{-1} = -I\) up to scalars, i.e. \(J = \operatorname{diag}(\sqrt{d_K}, 1)\)-type matrices,
which lie in \(\mathrm{GL}_2(\mathcal{O}_K)\) only at \(\mathbb{Q}(i)\). In the frame of the imaginary axis
(rotate by \(\zeta_8^{-1}\)) the arrangement is \(\Gamma'_K\cdot\hat{\mathbb{R}}\) with
\(\Gamma'_K = D_8^{-1}\mathrm{SL}_2(\mathcal{O}_K)D_8\), its levels \(-\operatorname{sgn}(A)\operatorname{Re}B \in \tfrac12\mathbb{Z}\)
are rational, but its hyperbolic centres \((-2\operatorname{Im}B + i\sqrt{k'^2 - 4})/(2A)\),
\(k' = 2\operatorname{Re}B\), are quadratic over \(\mathbb{Q}(\sqrt{|d_K|})\) and not CM points over \(\mathbb{Q}\)
for \(K \ne \mathbb{Q}(i)\). \(\blacksquare\)

**Certified** (phase E, X1–X2): for \(K = \mathbb{Q}(i), \mathbb{Q}(\sqrt{-3}), \mathbb{Q}(\sqrt{-7}), \mathbb{Q}(\sqrt{-2}),
\mathbb{Q}(\sqrt{-11})\) the BFS orbit of \(i\hat{\mathbb{R}}\) (\(658\)–\(3464\) circles) equals the
congruence set of (1) for every curvature \(|A| \le 5\) or \(6\); every form of the
discriminants of (3), imprimitive included and with both orientations, is descended to
an explicit \(X\) at the first levels (\(y \le 6\), resp. \(y \le 8\) even); and the
\(\hat{\mathbb{R}}\)-levels of \(\tau(X)(i\hat{\mathbb{R}})\) scatter (at \(\mathbb{Q}(\sqrt{-3})\), \(y = 4\):
\(\{0, 2, 4, 8\}\); at \(\mathbb{Q}(\sqrt{-2})\), \(y = 6\): \(\{0, 8, 16, \infty\}\)).

| \(K\) | levels \(y\) and discriminants \(4 - y^2|d_K|\) (odd \(d_K\)) / \(1 - my^2\) (even) |
|---|---|
| \(\mathbb{Q}(\sqrt{-3})\) | 2: \(-8\); 3: \(-23\); 4: \(-44\); 5: \(-71\); 6: \(-104\); 7: \(-143\); 8: \(-188\) |
| \(\mathbb{Q}(\sqrt{-7})\) | 1: \(-3\); 2: \(-24\); 3: \(-59\); 4: \(-108\); 5: \(-171\); 6: \(-248\) |
| \(\mathbb{Q}(\sqrt{-11})\) | 1: \(-7\); 2: \(-40\); 3: \(-95\); 4: \(-172\); 5: \(-271\) |
| \(\mathbb{Q}(\sqrt{-2})\) | 2: \(-7\); 4: \(-31\); 6: \(-71\); 8: \(-127\) |
| \(\mathbb{Q}(\sqrt{-19})\) | 1: \(-15\); 2: \(-72\); 3: \(-167\) |

The family \(|D| + 4 = y^2|d_K|\) is disjoint from the family \(|d_K||D| + 4 = k^2\) of
\(\mathcal{S}_K\) (at \(\mathbb{Q}(\sqrt{-3})\): \(-8, -23, -44, \dots\) against \(-4, -7, -15, -20, \dots\)),
and it reaches discriminants of the unrealized list — \(-72\) at \(\mathbb{Q}(\sqrt{-19})\),
\(-104\) at \(\mathbb{Q}(\sqrt{-3})\) — but without a twist: the answer to the question of item
8(b) of the session prompt is that the second-kind orbit is a genuine second arrangement,
classified above, whose circles carry the Heegner points of a new family of discriminants
and no Atkin–Lehner coupling.

## 7. Genus structure of a slice

For a slice \((D, [\mathfrak{r}])\) the genus characters of \(\mathrm{Cl}(D)\) are the
quadratic characters; those with \(\chi(\mathfrak{r}) = -1\) contribute the closed forms of
Paper II §3 to the unit system, the others nothing. The genus of \([\mathfrak{r}]\) is read off
from the assigned characters (Cox, Thm 3.15: \((\cdot/p)\) at the odd primes \(p \mid D\), and
\(\delta, \varepsilon, \delta\varepsilon\) at \(2\) according to \(D/4 \bmod 8\)) evaluated at a
number represented by the twist form and prime to \(D\).

> **Theorem 7 (genus of the twist).** At the level \(k\) of \(\mathcal{S}_K\) (any \(K\)),
> with \(g_\mp = \gcd(k \mp 2, |d_K|)\): for an odd prime \(p \mid D\),
> $$\chi_p(\mathfrak{r}_\alpha) = \Bigl(\frac{g_+}{p}\Bigr)\ \text{ if } p \mid r_0,\qquad
> \chi_p(\mathfrak{r}_\alpha) = \Bigl(\frac{-g_-}{p}\Bigr)\ \text{ if } p \mid s_0 .$$
> For odd \(d_K\) this reads: \((|d_K|/p) = (d_K/p)\cdot(-1/p)\)-type values —
> \(\chi_p(\mathfrak{r}) = (|d_K|/p)\) at \(p \mid r_0\) and \((-1/p)\) at \(p \mid s_0\) when
> \(k \equiv -2 \bmod |d_K|\), and \(\chi_p(\mathfrak{r}) = 1\) at \(p \mid r_0\),
> \((-|d_K|/p) = (d_K/p)\) at \(p \mid s_0\) when \(k \equiv 2 \bmod |d_K|\). At \(\mathbb{Q}(i)\)
> (odd \(n\), \(g_\pm = 4\)): \(\chi_p(\mathfrak{r}_n) = 1\) at \(p \mid r_0\), \((-1/p)\) at \(p \mid s_0\)
> — Paper II's Proposition 3.1 — and at the even levels (\(g_\pm = 2\)): \((2/p)\) at
> \(p \mid n - 1\), \((-2/p)\) at \(p \mid n + 1\). The \(2\)-adic characters of \(\mathfrak{r}_\alpha\)
> depend only on \(k \bmod 32\) (even \(d_K\)) resp. \(k \bmod 16|d_K|\) (odd \(d_K\)).
> Consequently \([\mathfrak{r}_\alpha]\) lies in the principal genus (no genus character is odd,
> and the closed forms contribute nothing) iff all these symbols are \(+1\).

*Proof.* The twist form represents \(s_0\) (at \((0, 1)\)) and \(r_0\) (at \((1, 0)\)); for
\(p \mid r_0\) use \(s_0\): \(s_0g_+ = k + 2 \equiv 4 \bmod p\) (as \(p \mid k - 2\)), so
\((s_0/p) = (4g_+^{-1}/p) = (g_+/p)\); for \(p \mid s_0\) use \(r_0\): \(r_0g_- = k - 2 \equiv -4\),
so \((r_0/p) = (-g_-/p)\). For the \((r_0, r_0, c)\) type the represented values are \(c \equiv s_0/4\)
at \(p \mid r_0\) and \(r_0\) at \(p \mid s_0\), giving the same symbols. \(g_\pm \in \{1, |d_K|\}\)
for odd \(d_K\) and \(\{4, 16\}\) resp. \(\{2, 4\}\)-type values at \(\mathbb{Q}(i)\) give the
displayed cases; the \(2\)-adic characters are functions of the represented values mod
\(8\), which depend on \(k\) modulo a fixed power of \(2\) times \(|d_K|\). \(\blacksquare\)

**Certified** (phase G): at the \(80\) levels of the nine fields with \(|D| \le 300\): the
assigned characters are class functions, multiplicative, with kernel \(\mathrm{Cl}(D)^2\)
and image of order \(2^{\mu-1}\); the genus vector of \(\mathfrak{r}_\alpha\) obeys the odd-prime law
at every level; the \(2\)-adic characters are functions of \(k \bmod 32\) resp.
\(k \bmod 16|d_K|\) (e.g. at \(\mathbb{Q}(i)\): \(\delta(\mathfrak{r}) = +\) for \(k \equiv 2, 18 \bmod 32\),
\(-\) for \(14, 30\); \(\varepsilon(\mathfrak{r}) = +\) for \(2, 30\), \(-\) for \(10, 22\); \(\delta\varepsilon = +\)
at \(6\), \(-\) at \(26\)). **The odd genus characters and their real fields**: for each odd
quadratic \(\chi = \chi_{d_1,d_2}\) (\(D = d_1d_2\), \(d_2 > 0\)) the closed form
\(L'(0,\chi) = \tfrac{2h(d_1^*)}{w(d_1^*)}h(d_2^*)\log\varepsilon_{d_2^*}\cdot C(0)\) was
verified with PARI's `qfbclassno`/`quadunit` against the independent evaluation at the
slices \(D = -84\) (three classes), \(-120\), \(-168\) (two classes), \(-160\) (two classes),
with \(C(0) = 1\) for the primitive characters and \(C(0) = 4\) for
\(\chi_{-32,5}\) of \(D = -160\) (conductor \(2\) over \(-40\)); the real fields
\(\mathbb{Q}(\sqrt{d_2})\) of the odd characters:

| slice | odd genus characters \((d_1, d_2)\) |
|---|---|
| \(-84\), \((5,4,5)\) | \((-7, 12)\), \((-3, 28)\) — \(\mathbb{Q}(\sqrt3), \mathbb{Q}(\sqrt7)\) |
| \(-84\), \((3,0,7)\) | \((-7, 12)\), \((-4, 21)\) — \(\mathbb{Q}(\sqrt3), \mathbb{Q}(\sqrt{21})\) |
| \(-84\), \((2,2,11)\) | \((-3, 28)\), \((-4, 21)\) — \(\mathbb{Q}(\sqrt7), \mathbb{Q}(\sqrt{21})\) |
| \(-120\), \((5,0,6)\) (\(n = 11\)) | \((-15, 8)\), \((-3, 40)\) — \(\mathbb{Q}(\sqrt2), \mathbb{Q}(\sqrt{10})\) (Paper II) |
| \(-168\), \((6,0,7)\) (\(n = 13\)) | \((-7, 24)\), \((-8, 21)\) — \(\mathbb{Q}(\sqrt6), \mathbb{Q}(\sqrt{21})\) (Paper II) |
| \(-168\), \((3,0,14)\) (\(\mathbb{Q}(\sqrt{-5})\)) | \((-7, 24)\), \((-3, 56)\) — \(\mathbb{Q}(\sqrt6), \mathbb{Q}(\sqrt{14})\) |
| \(-160\), \((7,6,7)\) | \((-32, 5)\) [\(C(0) = 4\)], \((-4, 40)\) — \(\mathbb{Q}(\sqrt5), \mathbb{Q}(\sqrt{10})\) |
| \(-160\), \((4,4,11)\) | \((-4, 40)\), \((-20, 8)\) — \(\mathbb{Q}(\sqrt{10}), \mathbb{Q}(\sqrt2)\) |

so the three fields \(\mathbb{Q}(\sqrt3), \mathbb{Q}(\sqrt7), \mathbb{Q}(\sqrt{21})\) of the genus theory of
\(-84\) are distributed two per slice, each field on exactly two of the three slices —
the "complementary halves" of Paper II §3.3 in the horizontal direction. The real
fields occurring at the first levels of each \(K\) (all odd genus characters,
\(|D| \le 300\)): \(\mathbb{Q}(\sqrt d)\) for \(d = 5, 8, 12, 13, 17, 21, 24, 28, 29, 33, 40, 53, 56, 57, 60, 65, 69, 76, 85, 88, 92\)
(e.g. \(\mathbb{Q}(\sqrt{13}), \mathbb{Q}(\sqrt5)\) at \(-195\); \(\mathbb{Q}(\sqrt{85}), \mathbb{Q}(\sqrt5)\) at \(-255\);
\(\mathbb{Q}(\sqrt{69})\) at \(-207\) and \(-276\)).

## 8. The odd index on a slice

> **Theorem 8 (odd index on a slice; = [robert-index-full.md](robert-index-full.md)
> Theorem 3 with an arbitrary twist).** Let \((D, [\mathfrak{r}])\) be a slice with
> \([\mathfrak{r}] \ne 1\), \(H = H_D\), \(\tau = \sigma_\mathfrak{r}\), \(H^+ = H^\tau\), \(E^\pm\) the even and odd
> units modulo torsion, \(Q^- = [E : E^+E^-]\). Then
> $$\bigl[E^- : \langle R_\mathfrak{b} : \mathfrak{b} \in \mathrm{Cl}(D)\rangle\bigr]
> \;=\; 24^{h/2}\cdot\frac{2^{h/2-1}}{Q^-}\cdot\frac{h_Hw_{H^+}}{h_{H^+}w_H}\cdot\prod_{\chi(\mathfrak{r}) = -1}C_\chi(0).$$

*Proof.* Verbatim the proof of Theorem 3 of robert-index-full.md: the twisted group
determinant over a transversal of \(\langle\mathfrak{r}\rangle\), the relative class number
formula for \(H/H^+\), and the regulator relation \(R_H/R_{H^+} = 2^{h/2-1}R^-(E^-)/Q^-\);
nothing there uses \(s_0 = r_0 + 1\). \(\blacksquare\)

**Certified with PARI/GP** (phase I; `polclass`, `bnfinit` at \(120\) digits, `bnfcertify`
(all degrees \(\le 8\): unconditional), `nfroots` of the slice polynomial, `bnfisunit`,
`nfgaloisconj` for \(\tau\), Smith forms; the multipliers through exact projections
\(\mathrm{Cl}(D) \to \mathrm{Cl}(f'^2d_{K'})\)):

| slice | realization | \(\mathrm{Cl}(D)\) | \(h_H\), \(w_H\) | \(h_{H^+}\), \(w_{H^+}\) | \(Q^-\) | \(\prod_{\rm odd}C\) | \([E^-:\langle R\rangle]\) | Smith form |
|---|---|---|---|---|---|---|---|---|
| \(-15\), \((2,1,2)\) | \(i\mathcal{S}\) 4, \(\mathbb{Q}(\sqrt{-3})\) 7 | \(\mathbb{Z}/2\) | 1, 6 | 2, 2 | 1 | 1 | \(4 = 24/6\) | \([4]\) |
| \(-35\), \((3,1,3)\) | \(i\mathcal{S}\) 6 | \(\mathbb{Z}/2\) | 1, 2 | 2, 2 | 1 | 1 | \(12 = 24/2\) | \([12]\) |
| \(-63\), \((4,1,4)\) | \(i\mathcal{S}\) 8 | \(\mathbb{Z}/4\) | 1, 6 | 1, 6 | 2 | 1 | \(576 = 24^2\) | \([24, 24]\) |
| \(-99\), \((5,1,5)\) | \(i\mathcal{S}\) 10 | \(\mathbb{Z}/2\) | 1, 6 | 1, 2 | 1 | 1 | \(8 = 24/3\) | \([8]\) |
| \(-195\), \((7,1,7)\) | \(i\mathcal{S}\) 14 | \((\mathbb{Z}/2)^2\) | 4, 6 | 4, 6 | 2 | 1 | \(576 = 24^2\) | \([48, 12]\) |
| \(-84\), \((5,4,5)\) | \(\mathbb{Q}(\sqrt{-3})\) 16 | \((\mathbb{Z}/2)^2\) | 1, 12 | 2, 4 | 2 | 1 | \(96 = 24^2/6\) | \([24, 4]\) |
| \(-84\), \((3,0,7)\) | \(\mathbb{Q}(\sqrt{-2})\) 26 | \((\mathbb{Z}/2)^2\) | 1, 12 | 2, 6 | 2 | 1 | \(144 = 24^2/4\) | \([24, 6]\) |
| \(-84\), \((2,2,11)\) | \(\mathbb{Q}(\sqrt{-19})\) 40 | \((\mathbb{Z}/2)^2\) | 1, 12 | 2, 2 | 2 | 1 | \(48 = 24^2/12\) | \([24, 2]\) |
| \(-160\), \((7,6,7)\) | \(\mathbb{Q}(\sqrt{-3})\) 22 | \((\mathbb{Z}/2)^2\) | 1, 8 | 2, 2 | 2 | 4 | \(288 = 24^2/2\) | \([48, 6]\) |
| \(-160\), \((4,4,11)\) | \(\mathbb{Q}(\sqrt{-11})\) 42 | \((\mathbb{Z}/2)^2\) | 1, 8 | 1, 2 | 2 | 1 | \(144 = 24^2/4\) | \([24, 6]\) |
| \(-160\), \((5,0,8)\) | none | \((\mathbb{Z}/2)^2\) | 1, 8 | 2, 4 | 2 | 4 | \(576 = 24^2\) | \([48, 12]\) |
| \(-168\), \((3,0,14)\) | \(\mathbb{Q}(\sqrt{-5})\) 58 | \((\mathbb{Z}/2)^2\) | 1, 6 | 2, 2 | 2 | 1 | \(96 = 24^2/6\) | \([24, 4]\) |

Every row satisfies Theorem 8 exactly (index from `bnfisunit` against the right side;
the regulator relation and the relative class number formula to \(\ge 119\) spare digits).
Readings: the three slices of \(D = -84\) share \(H\) (\(h_H = 1\), \(w_H = 12\)) but have
different \(H^+ = H^{\sigma_\mathfrak{r}}\) — \(w_{H^+} = 4, 6, 2\) (\(H^+ \ni i\), \(\ni\zeta_3\), neither)
— so the index takes the three values \(96, 144, 48 = 24^2\cdot w_{H^+}/w_H\cdot\tfrac{h_H}{h_{H^+}}\);
the multiplier \(C(0) = 4\) of the imprimitive character of \(-160\) (conductor \(2\) over
\(-40\)) enters at the two slices where that character is odd; \(Q^- = 2^{h/2-1}\) at all
twelve slices (the \(2\)-adic factor of Theorem 8 is \(1\) throughout — the exception
\(Q^- = 2^{h/2-2}\) of Paper II's \(n = 19\) has no analogue here). The \(24\)-saturated index
(the index divided by \(\prod\gcd(s_i, 24)\) over the Smith invariants \(s_i\)) is \(1\) except
at the three slices with a nontrivial \(h_H/h_{H^+}\) or multiplier (\(2\) at
\(-195\ (7,1,7)\), \(-160\ (7,6,7)\), \(-160\ (5,0,8)\)) — consistent with the Gras-type
reading of robert-index-full.md Conjecture 5.3 on the odd part. The larger \(2\)-torsion
\(|\mathrm{Cl}(D)[2]| = 8\) (\(D = -480\), three fields) is prepared in the script
(`index --with-480`, degree-\(16\) fields) but was not run in the selftest.

## 9. The phase on a slice

Let \((D, [\mathfrak{r}])\) be a slice, \(D < -4\), \(\mathfrak{r}\) **the twist of smaller norm** of the pair
(\(r_0 < s_0\)), and for \(\mathfrak{b} \in \mathrm{Cl}(D)\) put
\((\beta_1, \beta_2) := (j(\mathfrak{b}), j(\mathfrak{r}^{-1}\mathfrak{b}))\), a point of the fiber product
\(X_0(r_0)\times_{X(1)}X_0(s_0)\) (\(\mathfrak{r}^{-1}\mathfrak{b} = r_0^{-1}\mathfrak{r}\mathfrak{b}\) is cyclically
\(r_0\)-isogenous to \(\mathfrak{b}\), and \(s_0\)-isogenous through \(\mathfrak{s}\)).

> **Theorem 9 (the phase on a slice).** 1. *(Uniqueness of the branch.)* The
> number of cyclic \(r_0\)-isogenies from \(E_\mathfrak{b}\) to a curve with \(j = \beta_2\), modulo
> automorphisms, is the number of proper representations of \(r_0\) by the twist form
> modulo automorphs; for the smaller-norm twist and \(D < -4\) it is \(1\) (the second
> type: \(4r_0 = r_0(2x+y)^2 + s_0y^2\) forces \(y = 0\) when \(s_0 > r_0\) unless
> \(s_0 \le 4\), i.e. \(D = -3\); the first type: \(r_0 = r_0x^2 + s_0y^2\) forces \(y = 0\)).
> Hence \((\beta_1, \beta_2)\) is a smooth point of \(\Phi_{r_0} = 0\) with
> \(\Phi_x\Phi_y \ne 0\), and
> $$u_\mathfrak{b} \;:=\; -r_0\frac{h_2(\mathfrak{b})}{h_2(\mathfrak{r}^{-1}\mathfrak{b})} \;=\; \frac{\Phi_y}{\Phi_x}(\beta_1, \beta_2),
> \qquad \Phi = \Phi_{r_0},$$
> is well defined on every slice.
> 2. \(u_\mathfrak{b} \in H_D\), and \(\sigma(u_\mathfrak{b}) = u_{\mathfrak{b}^{e(\sigma)}\mathfrak{c}(\sigma)}\) for every
> \(\sigma \in \mathrm{Gal}(\bar{\mathbb{Q}}/\mathbb{Q})\) (first-power dihedral law).
> 3. *(\(j\)-dressing.)* \(u_\mathfrak{b}^6 = R_\mathfrak{b}\cdot\dfrac{\beta_1^4(\beta_1 - 1728)^3}{\beta_2^4(\beta_2 - 1728)^3}\).
> 4. *(Partner.)* The partner \(\mathfrak{s}\) gives \(-u_\mathfrak{b}\): \(-s_0h_2(\mathfrak{b})/h_2(\mathfrak{s}^{-1}\mathfrak{b}) = -u_\mathfrak{b} = \Phi^{(s_0)}_y/\Phi^{(s_0)}_x(\beta_1,\beta_2)\).
> 5. The slice polynomial \(\Pi(x) = \prod_\mathfrak{b}(x - u_\mathfrak{b}) \in \mathbb{Q}[x]\) is a power of an
> irreducible polynomial, irreducible iff the \(u_\mathfrak{b}\) are pairwise distinct.
> On the Gaussian slices \(u_f = \varepsilon\Theta_f\) with \(\omega_f = 1\): the circle
> geometry produces \(u\) up to the Cartan unit \(\varepsilon_K(\alpha)\) of the realization.

*Proof.* (1) is Lemma 1.2 of [first-power-descent.md](first-power-descent.md) with its
hypothesis \(s_0 > r_0\) replaced by the representation count stated (the kernels of the
isogenies in question are \(E[\mathfrak{c}']\) for invertible \(\mathfrak{c}'\) of norm \(r_0\) with
\([\mathfrak{c}'] = [\mathfrak{r}]\), i.e. proper representations of \(r_0\) by the form of \([\mathfrak{r}]\)); the
formula is Proposition 1.3 there. (2) is Proposition 2.1 there (Shimura reciprocity for
\(j\), \(\Phi \in \mathbb{Z}[x,y]\), \(\bar{\mathfrak{r}} = \mathfrak{r}\)); it needs no \(\varepsilon\) and no \(\omega_f\).
(3) is the weight algebra \(h_2^6 = c\,j^4(j-1728)^3\Delta\) applied to \(u^6 = r_0^6h_2(\mathfrak{b})^6/h_2(\mathfrak{r}^{-1}\mathfrak{b})^6\)
(Paper II, Lemma 4.1, in the slice-intrinsic form). (4): \(\mathfrak{s}^{-1}\mathfrak{b} = (r_0/\theta)\mathfrak{r}^{-1}\mathfrak{b}\)
and \(h_2\) has weight \(2\), so \(h_2(\mathfrak{s}^{-1}\mathfrak{b}) = (\theta/r_0)^2h_2(\mathfrak{r}^{-1}\mathfrak{b}) = -(s_0/r_0)h_2(\mathfrak{r}^{-1}\mathfrak{b})\).
(5) is Theorem 4.1 of first-power-descent.md (single Galois orbit). \(\blacksquare\)

**Certified** (phase F, \(150\) digits; exact \(\Phi_m\) from `first_power_descent.py` for
\(m \le 12\) and from PARI's `polmodular` for prime \(m \le 40\); \(H_D\) by certified
rounding; the exact \(\Pi\) through \(\mathbb{Q}[t]/(H_D)\)): at the sixteen slices
\(D = -15, -35, -63, -99, -195\) (\(i\mathcal{S}\)), \(-24, -20, -36\), \(-84\) (three), \(-160\) (three),
\(-168\ (3,0,14)\), \(-56\ (2,0,7)\): \(\Phi_{r_0}(\beta_1,\beta_2) = 0\) (spare \(\ge 149\)),
the \(j\)-dressing identity between the two independently computed sides (spare
\(\ge 146\)), the partner law where \(\Phi_{s_0}\) is available (\(s_0 \le 12\) or prime),
and the exact \(\Pi\), irreducible and squarefree, matching the numerical roots (spare
\(\ge 115\)). Where several primitive invertible ideals of norm \(r_0\) exist the
pairing gcd \(\gcd(H_D(y), \Phi_{r_0}(t, y))\) over \(\mathbb{Q}[t]/(H_D)\) has degree \(> 1\) — at
\(D = -160\), \((7,6,7)\), \(r_0 = 8\), the ideals \((8,8,7)\) and \((8,0,5)\) — and the
Atkin–Lehner branch is isolated exactly on the fiber product, by the gcd with
\(\Phi_{s_0}\) or by dividing out the branch of the other ambiguous ideal through its own
partner (\(\Phi_5\) here); at all other slices \(\mathfrak{r}\) is the only primitive ideal of
norm \(r_0\) and the gcd has degree \(1\). First slice polynomials (primitive integer form):

| slice | \(\Pi\) |
|---|---|
| \(-15\), \((3,3,2)\), \(r_0 = 3\) | \(11x^2 + 4783x + 11\) |
| \(-35\), \((5,5,3)\), \(r_0 = 5\) | \(589x^2 + 221979023x + 589\) |
| \(-99\), \((9,9,5)\), \(r_0 = 9\) | \(95948x^2 + 25717824888795379x + 95948\) |
| \(-20\), \((2,2,3)\), \(r_0 = 2\) | \(209x^2 - 538338x + 209\) |
| \(-36\), \((2,2,5)\), \(r_0 = 2\) | \(5819x^2 - 144395894x + 5819\) |
| \(-24\), \((2,0,3)\) | \(6647x^2 + 30594194x + 6647\) (Paper I, \(n = 5\)) |
| \(-56\), \((2,0,7)\) | \(85646669273x^4 + 21816145776449948x^3 + 43628432475562262x^2 + \cdots\) |

(The \(j\)-dressing identity is what makes the phase computable from the \(\Delta\)-quotient
route, and conversely; on a slice without arrangement, e.g. \((-160, [(5,0,8)])\), \(u\) is
defined by (1) alone.)

## 10. Status: proved, certified, experimental; what failed

| statement | status |
|---|---|
| Thm 1 (slices): types, count \(2^\mu\), pairs, \(\mathrm{Cl}[2]\), \([\mathfrak{r}] = 1\) iff \(\min = 1\), rank \(h/2\) | proved (Gauss + Theorem 4 of other-fields.md); certified \(-400 \le D < -4\) |
| Thm 2 (levels and slices of \(\mathcal{S}_K\), \(t\)-parametrization, criterion) | proved for every \(K\) (principal cusp); certified at 9 fields, \(|D| \le 300\), and at \(d_K = -20, -24, -40, -52\) |
| Thm 3 (strata; Pell; class \(= c^j\); odd-prime rule) | proved; the \(2\)-part of the stratum twist certified (\(14 \times 9\) pairs, \(4\) solutions each) |
| Thm 4 (units on slices, laws, KLF) | proved (other-fields.md Thm 4); certified at 54 slices + 7 strata, 200 digits |
| Thm 5 (\(i\mathcal{S}\) class formula, Lemmas A′′–C′′) | proved; certified \(n \le 16\), both orientations |
| Thm 6 (\(\mathcal{S}_K^\perp\)): necessity; sufficiency; levels/forms; no twist | necessity proved (all \(K\)); sufficiency proved for the five Euclidean fields; the no-twist statement proved for involutions of the type \(J\bar X^{-1}J'\), \(J, J' \in \mathrm{GL}_2(\mathcal{O}_K)\) |
| Thm 7 (genus of the twist) | proved; certified at 80 levels; closed forms at 8 slices with PARI |
| Thm 8 (odd index on a slice) | proved (Theorem 3 of robert-index-full.md); exact and unconditional at 12 slices |
| Thm 9 (phase on a slice) | proved; certified at 16 slices, exact \(\Pi\) irreducible |
| the \(2\)-adic law of the stratum twist; the saturation of the odd index; \(Q^- = 2^{h/2-1}\) | certified only |
| the non-principal cusps at \(h_K > 1\) (\(k \equiv \pm10 \bmod 24\) etc.) | open (not realized here) |

**Where the expected statements of the session prompt failed, and the corrected
statements.** (i) The list of "slices with a nontrivial class realized by no primitive
stratum" is twelve slices at ten discriminants, not eight discriminants: \(D = -96\) and
\(-120\) are realized at one class and not at the others. (ii) \(k = 58\) at
\(\mathbb{Q}(\sqrt{-6})\) and \(k = 98\) at \(\mathbb{Q}(\sqrt{-10})\) are not levels of the principal cusp
(\(k \not\equiv \pm2 \bmod |d_K|\)); the principal cusps realize the slices of Theorem 2 at
\(k = 4mt \pm 2\), e.g. \((-88, (2,0,11))\) and \((-104, (2,0,13))\) at \(\mathbb{Q}(\sqrt{-6})\). (iii) The
stratum twist class is a homomorphism \(j \mapsto c^j\) as expected, but "the genus of the
fundamental solution" is made precise as the sign pattern \(\eta \equiv \pm1 \bmod \mathfrak{p}\)
at the odd primes of \(D'\); it is not true that every slice is realized by every
\(\mathcal{S}_K\): each \((K, D')\) realizes exactly one nontrivial class on its strata, and
\((-160, (5,0,8))\) is realized by no stratum of any class-number-one field. (iv) The
second-kind orbit has irrational levels and no twist; it carries discriminants
\(4 - y^2|d_K|\), not slices. (v) The uniqueness of the branch holds on every slice
(\(D < -4\), smaller-norm twist), but the exact \(\Pi\) route needs the fiber product
when \(\mathfrak{r}\) is not the only primitive ideal of norm \(r_0\). (vi) All predictions about
\(i\mathcal{S}\) (twist \((n-1,n-1,\tfrac n2)\), \(x^2 - 7x + 1\) at \(n = 4\)) and the ten stratum
instances of the prompt were confirmed.

## 11. Machine verification

`python3 scripts/horizontal_families.py --selftest` (needs `mpmath`, `sympy`, PARI/GP
2.15.4; \(270\) s on one core; every phase also runs standalone:
`slices`, `realize`, `pell`, `units`, `even`, `other`, `genus`, `index [--with-480]`, `phase`).
Precision is set in `main()`; integers and rationals are accepted only with
\(\ge\max(20,\mathrm{dps}/5)\) spare digits in absolute error; forms, ideals, class groups,
characters, stratum twists and Pell solutions are exact (HNF/Smith arithmetic); PARI
results are labelled by `bnfcertify` (all \(1\) here); no PSLQ.

- **S** (\(0.2\) s): Theorem 1 at \(-400 \le D < -4\).
- **R** (\(< 1\) s): the criterion, the first slices of the nine fields, the coincidence
  table, the unrealized slices, the \(t\)-parametrization and its converse at 287 pairs.
- **P** (\(1\) s; PARI `quadunit`): the ten instances, the law at \(14 \times 9\) pairs with
  \(4\) solutions each, the realization-by-strata table.
- **U** (\(232\) s, \(200\) digits): units at 54 slices (integer palindromic polynomials,
  laws, partner independence, KLF at every character), the regression record
  `UNIT_RECORD`, the equal-slice / different-class comparisons, the seven strata.
- **E/X** (\(< 1\) s): \(i\mathcal{S}\) at \(n \le 16\) with Lemmas A′′–C′′; \(\mathcal{S}_K^\perp\) at
  five fields (BFS, descent, forms census, \(\tau\)-levels).
- **X′** (\(< 1\) s): the class formula with Lemmas A′–C′ through the explicit \(P\) at the
  non-Euclidean fields and at the principal cusps of \(d_K = -20, -24, -40, -52\).
- **G** (\(5\) s, \(120\) digits; PARI): genus vectors at 80 levels; the closed forms.
- **I** (\(27\) s, \(150\) digits; PARI \(120\) digits): the odd index at 12 slices, regression
  record `INDEX_RECORD_EXPECTED`.
- **F** (\(5\) s, \(150\) digits): the phase at 16 slices, exact \(\Pi\).

Observed spare digits: \(\ge 169\) (units), \(\ge 197\) (KLF), \(\ge 119\) (index identities),
\(\ge 146\) (\(j\)-dressing), \(\ge 115\) (exact vs numerical \(\Pi\)).

## 12. Literature diligence

Searches ran through the web-search relay; every full-text fetch attempted in this
session (Springer, arXiv, Numdam/EuDML mirrors, the Leiden, MIT, BU and Chicago
course pages) was blocked by the egress proxy, so — as in the papers' `NOTES.md` — the
statements below are verified only as far as the relay's bibliographic records go, and
the precise chapter references are given from memory and marked as such.

1. **C. F. Gauss**, *Disquisitiones Arithmeticae* (1801), Art. 257–262 (formae ancipites:
   the forms \((a, b, c)\) with \(a \mid b\); the number of ambiguous classes) and Art. 305
   (the \(2\)-torsion). Theorem 1(1)–(5) is this theory in ideal language, with Cox's
   normalization; **not new**. **D. A. Cox**, *Primes of the form \(x^2 + ny^2\)*, 2nd ed.,
   Thm 3.15 (the assigned characters, \(2^{\mu-1}\) genera) and §7 (orders, proper ideals):
   used for \(\mu\) and for the genus theory of §7 (relay: bibliographic record only).
2. **S. Lang**, *Elliptic Functions*, 2nd ed., GTM 112, Ch. 12 (the "\(\Delta\)-quotients"
   chapter: for fractional ideals \(\mathfrak{a}, \mathfrak{b}\) of the maximal order,
   \(\Delta(\mathfrak{a})/\Delta(\mathfrak{b})\) lies in the Hilbert class field and generates the ideal
   \((\mathfrak{b}\mathfrak{a}^{-1})^{12}\); from memory — the chapter could not be read). This is the
   classical statement behind Theorem 1(6) at the maximal order (\(D\) fundamental): the
   unit property of \(R_\mathfrak{b}\) is a corollary there. For non-maximal orders Paper II's
   five lemmas (other-fields.md Theorem 4) are an independent proof; whether Schertz's
   *Complex Multiplication* (New Math. Monographs 15, CUP 2010, Ch. 6–7: elliptic units,
   ring class fields) states the ideal factorization for proper ideals of orders could
   not be checked (relay: table of contents only, as in Paper II's notes). **Credited**:
   the unit property; **new**: the slice bookkeeping (pairs, types), the realization
   theory, the stratum law, the index and the phase on slices.
3. **D. Kubert, S. Lang**, *Modular Units*, Grundlehren 244 (1981): \(\Delta(N\tau)/\Delta(\tau)\)
   is a modular unit (an \(S\)-unit at CM points, \(S\) the primes above \(N\)) — the
   \(12\)-th power of the Siegel-function products; the relay returned only the bibliographic
   record and Rohrlich's survey title. Our \(r_0^6\Delta(\mathfrak{b})/\Delta(\mathfrak{r}^{-1}\mathfrak{b})\) is the
   specialization of \(N^6\Delta(N\tau)/\Delta(\tau)\)-type quotients along the ambiguous
   ideal, whose \(S\)-part cancels because \(\mathfrak{r}^2 = (r_0)\) (Paper II's balance lemma).
4. **K. E. Stange**, IMRN 2018 and Trans. AMS 2018 (verified from abstracts in
   other-fields.md §12): the arrangements \(\mathcal{S}_K\) and their curvature integrality. The
   second-kind orbit \(\mathcal{S}_K^\perp\) (the orbit of the imaginary axis, the form \(h_1\))
   does not appear there; the relay found no paper on it.
5. **Ramachandra, Robert, Gillard–Robert, Schertz** (as positioned in Paper II's notes):
   elliptic units and their index; the relay's records (the Hasse–Ramachandra and
   Schertz-conjecture papers, Jung–Koo–Shin's ring class invariants) confirm that
   \(\Delta\)-quotients at proper ideals of orders generate ring class fields (Schertz's
   argument), which is the field statement \(R_\mathfrak{b} \in H_D\), not the slice theory.

**Verdict.** Not new: Gauss's ambiguous forms; the unit property of \(\Delta\)-quotients along
ambiguous ideals (maximal order: Lang; orders: Paper II / Schertz-type arguments). New
here: the slice bookkeeping with the type, the \(t\)-parametrization of the levels of
\(\mathcal{S}_K\) and the realization criterion, the coincidence and non-realization tables,
the stratum law along the Pell solutions, the even levels of \(i\mathcal{S}\) with the odd-\(D\)
Gram form, the classification of the second-kind orbit and its no-twist property, the
genus law of the twist, and the index and phase theorems on slices.

## 13. Outlook

- **The non-principal cusps.** At \(h_K > 1\) the residues \(k^2 \equiv 4 \bmod |d_K|\),
  \(k \not\equiv \pm2\), carry the discriminants \(-(k^2-4)/|d_K|\) (\(-3, -7, -39, -51, -115, \dots\)
  at \(\mathbb{Q}(\sqrt{-5})\); \(-140\) at \(k = 58\) of \(\mathbb{Q}(\sqrt{-6})\), \(-240\) at \(k = 98\) of
  \(\mathbb{Q}(\sqrt{-10})\)); the conjecture of the prompt is that they are the levels of the
  orbit of \(\hat{\mathbb{R}}\) under \(\mathrm{SL}(L \oplus \mathcal{O}_K)\), \(L\) a non-principal class. The
  cusp bookkeeping of other-fields.md §13 is the prerequisite; the slice theory says
  what to expect (the pair \(\{t, mt \pm 1\}\)-type formulas with \(L\) twisted in).
- **The \(2\)-adic stratum law.** Prove the \(2\)-part of Theorem 3(3): the two types at
  even \(D'\) and the splitting of the \(2\)-part (\(\{2, 58\}\) at \(-116\)); with it the
  homomorphism \(j \mapsto c^j\) is a theorem for even \(D'\) as well.
- **Slices without arrangement.** \((-160, [(5,0,8)])\) has a unit system, an index
  (\(576\)) and a phase (irreducible quartic) but no geometric carrier among the
  arrangements of this document; the quaternionic carriers of item 6 of the program, or
  the non-principal cusps, are the candidates.
- **The second-kind orbit.** \(\mathcal{S}_K^\perp\) has its own hyperbolic dictionary
  (discriminants \(4 - y^2|d_K|\), irrational levels \(y\sqrt{|d_K|}/2\)); its Euclidean side
  (curvature-\(A\) disks modulo translation, \(A \in \operatorname{Tr}(\mathcal{O}_K)\)) and its
  relation to the isometry classes of the unimodular Hermitian forms \(h_0, h_1\) over
  \(\mathcal{O}_K\) are untouched; a Kronecker-limit theory for the \(\Delta\)-data of the
  \(\mathcal{S}_K^\perp\)-circles would be a Heegner-point theory without Atkin–Lehner partner.
- **The saturation and \(Q^-\) on slices.** \(Q^- = 2^{h/2-1}\) at all twelve slices and the
  \(24\)-saturated index \(\in \{1, 2\}\): the \(D = -480\) slices (\(|\mathrm{Cl}[2]| = 8\), three
  fields, degree-\(16\) fields) are the next test (`index --with-480`).
- **Paper III.** The slice theory (§1–§3, §7–§9) with the even levels (§5) is the natural
  third paper of the program: *Slices of the Schmidt unit systems*.
