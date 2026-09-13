# Prompt 10 — Item 3: the line over other imaginary quadratic fields

Item 3 of [PROGRAM.md](../PROGRAM.md): transport the line — classification, level,
involution, the Euclidean \(\Delta\)-data, the Kronecker limit formulas, the unit
theorem, the Euler system and the Robert index — from \(\mathbb{Q}(i)\) to Stange's
arrangements \(\mathcal{S}_K = \mathrm{PSL}_2(\mathcal{O}_K)\cdot\hat{\mathbb{R}}\) for other imaginary
quadratic fields \(K\), class number one first, then class number \(> 1\).

Read `CLAUDE.md` and [PROGRAM.md](../PROGRAM.md) §1 (what is field-independent) and
item 3; then Paper I (`papers/1-schmidt-circles/`) §2–§4 and §7 (classification,
level, census, involution, the Euclidean dictionary and the \(\Delta\)-mass law),
Paper II (`papers/2-schmidt-elliptic-units/`) §2, §4 and §6.3,
[schmidt-euler-system.md](../schmidt-euler-system.md) §1–§2 (the fiber lemma and the
norm relations) and [robert-index-full.md](../robert-index-full.md) (the index
theorem, its §8 "other fields"). Literature to read *and verify before use*:
K. E. Stange, *Visualising the arithmetic of imaginary quadratic fields* (IMRN 2018)
and *The Apollonian structure of Bianchi groups* (TAMS 2018) — the definition of
\(\mathcal{S}_K\), the curvature and centre description of its circles, and the
statement about which ideal classes a single cusp sees; credit them for everything
they contain (the Gaussian classification of Paper I §2 is already credited to
Stange). Record the diligence in the new document's notes, as the papers' `NOTES.md`
do.

Machinery that transports with small changes: `scripts/phase_klf.py` (`group_data`,
`characters`, `rep_numbers`, `epstein_Lprime0` are generic in the discriminant),
`scripts/robert_index_full.py` (the Dedekind determinant, the multiplier recursion,
the PARI glue — replace `y^2+1` by the polynomial of \(K\) and \(-4n^2\) by \(n^2d_K\)
in `polclass`, the torsion lemma by its \(K\)-version). Machinery that must be
rewritten: the lattice arithmetic of `scripts/schmidt_euler_system.py` (Hermite
normal forms with respect to \((1, \omega)\), \(\mathcal{O}_K = \mathbb{Z}[\omega]\), instead of
\((1, i)\); the projections \(\Lambda \mapsto \mathcal{O}_m\Lambda\); the ideal twists), the
\(\Delta\)-values (\(G_\mathfrak{c} = n^{12}\Delta(\Lambda_\mathfrak{c})/\Delta(\mathcal{O}_K)\), reference lattice
\(\mathcal{O}_K\)), and the hyperbolic `R_lattice` of `scripts/schmidt_units.py`
(the twist ideal changes). Follow the certification guard rails of CLAUDE.md.

## Established (do not redo)

Field-independent, to be *cited* and reused (PROGRAM.md §1): the Kronecker limit
formula and the independent incomplete-gamma evaluation (Paper II, Prop. 2.2,
Lemma 2.3 — valid for every discriminant); the rigidity and balance arguments of
the unit theorem (Paper II §4); the local argument of the per-class valuation law
(Paper II, Thm 4.11); the fiber lemma (schmidt-euler-system.md, Lemma 1.2 — local at
each prime); the Dirichlet-series form of the Hecke recursion
(robert-index-full.md, Lemma 3.2 — a theta-series identity over the index-\(\ell\)
sublattices); the index theorem (robert-index-full.md, Theorem 1) whose proof uses
only the Galois action \(\sigma_\mathfrak{a}(G_\mathfrak{c}) = G_{\mathfrak{a}^{-1}\mathfrak{c}}\), the
conjugation law, the class-independent valuations, \(S_\chi = -12L'(0,\chi)\) and
the class number formula, with \(\zeta_K(0) = -h_K/w_K\) in place of \(-\tfrac14\).

Gaussian-specific, to be *redone*: the \(\Delta\)-mass stratification (Paper I, Thm 7.18,
step M2) and the "every class exactly twice" dictionary (Paper I, Lemma 7.5) use that
\(\mathbb{Z}[i]\) is a PID with unit group \(\mu_4\); the hyperbolic coupling \(D = -4r_0s_0\),
\(s_0 = r_0 + 1\), and the twist ideal \(\mathfrak{r}\) with \(\mathfrak{r}^2 = (r_0)\) use the
Gaussian level; the torsion lemma \(\mu(H_n) = \mu_{4\cdot2^{[4|n]}3^{[3|n]}}\) uses \(i \in K\).

## Expected statements (derive, then prove or certify; hedge where they fail)

Notation: \(K = \mathbb{Q}(\sqrt{d_K})\), \(\mathcal{O}_K = \mathbb{Z}[\omega]\), \(w_K = |\mathcal{O}_K^\times|\), \(h_K\);
orders \(\mathcal{O}_n = \mathbb{Z} + n\mathcal{O}_K\) of discriminant \(n^2d_K\); ring class fields \(H_n\);
\(h = |\mathrm{Pic}(\mathcal{O}_n)|\).

1. **Level.** For \(X \in \mathrm{SL}_2(\mathcal{O}_K)\) the Cartan trace \(\alpha(X) = -\tfrac12\mathrm{tr}(X\bar X^{-1})\)
   is real and lies in \(\tfrac12\mathrm{tr}(\mathcal{O}_K) \subseteq \tfrac12\mathbb{Z}\); it equals
   \(\coth\) of the hyperbolic radius of the circle \(X(\hat{\mathbb{R}})\) when that circle
   lies in the upper half-plane. From \(M_{X(\hat{\mathbb{R}})} = (X^{-1})^\dagger M_0X^{-1}\)
   with \(X = \binom{a\ b}{c\ d}\): curvature \(2|\operatorname{Im}(c\bar d)|\), centre
   \(-i(a\bar d - b\bar c)/(2\operatorname{Im}(c\bar d))\), so curvatures lie in
   \(\sqrt{|d_K|}\,\mathbb{Z}\) (Stange) and \(\alpha = \mp\operatorname{Re}(a\bar d - b\bar c)\), an
   integer when \(d_K \equiv 0 \bmod 4\) and possibly a half-integer when
   \(d_K \equiv 1 \bmod 4\). Determine the set of levels that occur, for each \(K\).
2. **Hyperbolic dictionary.** The hyperbolic centre of a level-\(\alpha\) circle is
   \(x_0 + i\,r\sqrt{\alpha^2-1}\), a CM point of the field \(\mathbb{Q}(\sqrt{d_K(\alpha^2-1)})\)
   — *not* of \(K\). Expected: level-\(\alpha\) circles modulo \(\mathrm{SL}_2(\mathbb{Z})\times\mathrm{SL}_2(\mathbb{Z})\)
   correspond to the primitive (and imprimitive) form classes of a discriminant
   \(D_K(\alpha) \in d_K(\alpha^2-1)\,\mathbb{Q}^{\times2}\), with the census a Hurwitz class
   number as in Paper I §3. Find \(D_K(\alpha)\) exactly (it should be
   \(-4r_0s_0\)-like, with \(r_0, s_0\) built from \(\alpha \mp 1\) and \(d_K\)), and the
   ambiguous twist class \(\mathfrak{r}_\alpha\) of the involution
   \(\sigma(X) = \bar X^{-1}\) (class formula \(\hat\sigma[f] = [\mathfrak{r}_\alpha][f]^{-1}\);
   \(\mathfrak{r}_\alpha^2\) principal). The unit theorem for
   \(R_f = N(\mathfrak{r})^6\Delta(\mathfrak{b})/\Delta(\mathfrak{r}^{-1}\mathfrak{b})\) then follows from
   item 4 of PROGRAM.md (an invertible ambiguous twist with principal square is all
   the proof of Paper II §4 needs); check the ramified primes of the new CM field and
   \(p = 2, 3\) by hand.
3. **Euclidean dictionary.** Disks of a fixed curvature modulo translations by
   \(\mathcal{O}_K\) should correspond to the primitive proper \(\mathcal{O}_n\)-sublattices of
   \(\mathcal{O}_K\) of index \(n\), each class of \(\mathrm{Pic}(\mathcal{O}_n)\) counted \(w_K/2\) times
   (\(\Lambda\) and \(u\Lambda\), \(u \in \mathcal{O}_K^\times/\pm1\)): expected census
   \(N_e(n) = \tfrac{w_K}2\,|\mathrm{Pic}(\mathcal{O}_n)|\) — \(2h\) at \(\mathbb{Q}(i)\) (Paper I), \(3h\) at
   \(\mathbb{Q}(\sqrt{-3})\), \(h\) at \(\mathbb{Q}(\sqrt{-2}), \mathbb{Q}(\sqrt{-7})\). The Euclidean phase is
   normalized by the period of the CM curve with CM by \(\mathcal{O}_K\) (lemniscatic at
   \(\mathbb{Q}(i)\), equianharmonic at \(\mathbb{Q}(\sqrt{-3})\)); the \(\Delta\)-data
   \(G_\mathfrak{c} = n^{12}\Delta(\Lambda_\mathfrak{c})/\Delta(\mathcal{O}_K)\) are algebraic integers of \(H_n\)
   with the Galois law of Paper II Thm 2.6 (the transport lemma is Shimura
   reciprocity, field-independent), and \(D_n(x) = \prod(x - G_\mathfrak{c}) \in \mathbb{Z}[x]\).
   The \(\Delta\)-mass \(M(n) = \prod G_\mathfrak{c}\): redo step M2 of Paper I Thm 7.18 with the
   splitting of the primes \(p \mid n\) in \(K\) (split primes contribute nothing;
   inert and ramified ones do — get the exponents).
4. **Kronecker limit formulas.** \(\sum_\mathfrak{c}\chi(\mathfrak{c})\log|G_\mathfrak{c}| = -12L'(0,\chi)\)
   with the Epstein \(L\)-function of discriminant \(n^2d_K\) (\(w = 2\) for \(n \ge 2\);
   avoid \(n = 1\) at \(\mathbb{Q}(i), \mathbb{Q}(\sqrt{-3})\) where the order has \(w = 4, 6\)). Genus
   characters: for \(n^2d_K = d_1d_2\) the real field of a Euclidean genus character
   is \(\mathbb{Q}(\sqrt{d_2})\) with \(d_2 > 0\) built from the primes dividing \(n\) and the
   \(2\)-part of \(d_K\); predict \(\mathbb{Q}(\sqrt n)\) or \(\mathbb{Q}(\sqrt{n|d_K|})\) or
   \(\mathbb{Q}(\sqrt{n|d_K|/4})\) from genus theory of the order and test at the first
   levels (Paper II Prop. 3.1's method, Sturm bound).
5. **Per-class valuations and units.** \(v_\mathfrak{P}(G_\mathfrak{c})\) independent of \(\mathfrak{c}\)
   (the local argument of Paper II Thm 4.11 with \(\mathcal{O}_{K,p}\) in place of \(\mathbb{Z}_p[i]\));
   Newton polygons of the certified \(D_n\) at every \(p \mid n\) as the check. Hence
   \(\mathcal{V}_n = \langle G_\mathfrak{c}/G_{\mathfrak{c}'}\rangle \subset \mathcal{O}_{H_n}^\times\).
6. **Euler system.** Fiber lemma over \(\mathcal{O}_K\): of the \(\ell+1\) index-\(\ell\)
   sublattices of a primitive proper \(\mathcal{O}_n\)-lattice, \(1 + \chi_K(\ell)\) are
   \(\mathcal{O}_n\)-stable for \(\ell \nmid n\) (\(\chi_K\) the Kronecker symbol of \(d_K\)) and one
   for \(\ell \mid n\); norm relations with \(A(\ell)\) and the factor \(\pi^{12}\) at a
   ramified \(\ell = \mathfrak{p}^2\) (at \(\mathbb{Q}(i)\): \((1+i)^{12} = -64\)); Hecke recursion at
   every \(s\) (robert-index-full.md Lemma 3.2 verbatim). Verify at \(\ge 20\) pairs
   \((n, \ell)\) per field.
7. **The full Robert index.** For \(h_K = 1\):
   $$
   \bigl[\mathcal{O}_{H_n}^\times : \mu(H_n)\mathcal{V}_n\bigr] \;=\; \frac{24^{h-1}\prod_{\chi\ne1}|L'(0,\chi)|}{R_{H_n}}
   \;=\; \frac{w_K\cdot24^{h-1}}{w_{H_n}}\,h_{H_n}\prod_{\chi\ne1}C_\chi(0)
   $$
   (the \(4\) of \(\mathbb{Q}(i)\) is \(w_K\); at \(\mathbb{Q}(\sqrt{-3})\), \(w_K = 6\) and \(\mu_6 \subset \mu(H_n)\)
   always, so the powers of \(3\) rearrange). Prove \(\mu(H_n)\) by the criterion
   "\(K(\sqrt d) \subset H_n\) iff the conductor of \(K(\sqrt d)/K\) divides \(n\mathcal{O}_K\)"
   (robert-index-full.md Lemma 4.1). Verify with PARI at the first levels with
   \(h \ge 2\) (`polclass(n^2*d_K)`, `polcompositum` with the polynomial of \(K\),
   `bnfinit` at 120 digits, `bnfcertify` where it finishes, `nfroots` of \(D_n\),
   `bnfisunit` exponent matrices, Smith forms); test the saturation observation
   (index modulo \(24\)-th roots \(= h_{H_n}\prod C\)) and the layer theorem at a cubic
   level of each field. This is the cleanest test of Conjecture 5.3 of
   robert-index-full.md away from \(\mathbb{Q}(i)\).
8. **Class number \(> 1\)** (\(K = \mathbb{Q}(\sqrt{-5})\), \(\mathbb{Q}(\sqrt{-6})\), \(\mathbb{Q}(\sqrt{-15})\)):
   the Bianchi orbifold has \(h_K\) cusps, one per ideal class \([L]\); the arrangement at
   the cusp of \(L\) is built from \(\mathrm{SL}(L \oplus \mathcal{O}_K)\)-type matrices. Verify
   Stange's statement that a single cusp sees only \(\ker(\mathrm{Pic}(\mathcal{O}_n) \to \mathrm{Pic}(\mathcal{O}_K))\);
   set up the Euclidean dictionary on pairs \(\Lambda \subset L\) (\(\Lambda\) a proper
   \(\mathcal{O}_n\)-lattice of index \(n\) in the ideal \(L\)) over all cusps, so that the
   \(\Delta\)-data \(G_{\mathfrak{c}} = n^{12}\Delta(\Lambda)/\Delta(L)\) are indexed by all of
   \(\mathrm{Pic}(\mathcal{O}_n)\); then items 4–7 with \(\zeta_K(0) = -h_K/w_K\): the index acquires
   the factor \(w_K/h_K\). Expected: the single-cusp data have rank \(h/h_K - 1\) only;
   the full system is needed for full rank.

## Mission

1. Class number one: \(K = \mathbb{Q}(\sqrt{-2})\) (the cleanest: \(w_K = 2\), \(d_K \equiv 0 \bmod 4\)),
   then \(\mathbb{Q}(\sqrt{-3})\) (\(\mu_6\): the sharpest test of every root-of-unity
   bookkeeping — phase, \(\Delta\)-mass, \(\mu(H_n)\)), then \(\mathbb{Q}(\sqrt{-7})\)
   (\(d_K \equiv 1 \bmod 4\), half-integral levels?), with \(\mathbb{Q}(\sqrt{-11})\) as an
   out-of-sample check. Items 1–7 above for each, in one document
   `other-fields.md` (house style: statement–proof, one section per field plus a
   uniform-statement section, verification section, outlook), verified by
   `scripts/other_fields.py --selftest` (generic \(\mathcal{O}_K\)-lattice HNF arithmetic;
   the \(\Delta\)-data; KLF against the independent evaluation; norm relations; the
   index with PARI; regression records; runtimes recorded).
2. Class number two: \(K = \mathbb{Q}(\sqrt{-5})\) at the first two levels — the two-cusp
   dictionary and the index with the factor \(w_K/h_K = 1\); label everything
   experimental until the cusp bookkeeping is proved.
3. State the uniform theorems (level, dictionary, KLF, unit theorem, norm relations,
   index) with the field-dependent constants displayed (\(w_K\), \(h_K\), \(d_K\)'s
   \(2\)-part, the splitting of the conductor primes), and say for each which parts
   are proved for all \(K\), which are proved for the computed \(K\), and which are
   certified only.
4. Update CLAUDE.md (document map, a theorem 23 entry, ledger, scripts table),
   PROGRAM.md (item 3 status, and item 4 if the twist ideal of expected statement 2
   settles it), outlook.md; commit with clear messages and push.

## Guard rails

As in CLAUDE.md: precision set after imports; integers and rationals accepted only
with \(\ge\max(20,\mathrm{dps}/5)\) spare digits in absolute error; exact arithmetic
for lattices, class groups, characters and projections; no PSLQ-derived claim
without exact re-verification; PARI outputs labelled GRH-conditional unless
`bnfcertify` returned 1; proved / certified / experimental kept distinct. Do not
claim Stange's classification or the classical elliptic-unit theory of ring class
fields (Robert, Kubert–Lang, Schertz) as new; what is new is the Schmidt-geometric
dictionary at each \(K\), the explicit constants, and the index formulas. Where an
expected statement fails (levels, census, twist class, genus field), record the
failure and the corrected statement — the failures at \(\mathbb{Q}(\sqrt{-3})\) and at
\(d_K \equiv 1 \bmod 4\) are the information.
