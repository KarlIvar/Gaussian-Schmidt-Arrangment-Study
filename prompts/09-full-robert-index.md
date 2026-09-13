# Prompt 09 — The full Robert index of the Schmidt Δ-units

**Status: done** — results in [robert-index-full.md](../robert-index-full.md)
(Theorem 1: the full index; Theorem 2: the layer indices; Theorem 3: the hyperbolic
odd index), verified by `scripts/robert_index_full.py --selftest`, and ported to
Paper II §6.3 (Theorems 6.7, 6.12, 6.13, Conjecture 6.16). The prompt as it was run:

Item 2 of [PROGRAM.md](../PROGRAM.md). Read `CLAUDE.md`, then Paper II
(`papers/2-schmidt-elliptic-units/`, §2 and §6 in full),
[schmidt-euler-system.md](../schmidt-euler-system.md) (all of it), and
[schmidt-units.md](../schmidt-units.md) §5. Machinery:
`scripts/schmidt_euler_system.py` (exact class groups of \(\mathcal{O}_n\),
projections, characters, the \(\Delta\)-data \(G_\mathfrak{c}\) and coset units),
`scripts/phase_klf.py` (the independent incomplete-gamma evaluation of the Epstein
\(L'(0,\chi)\), 250 digits), `scripts/schmidt_units.py` (fundamental units of the
cubic layer, root descent, Friedman cap). Follow the certification guard rails of
CLAUDE.md throughout.

## Established (do not redo)

- \(G_\mathfrak{c} = n^{12}\Delta(\Lambda_\mathfrak{c})/\Delta(\mathbb{Z}[i])\), \(\mathfrak{c} \in \mathrm{Pic}(\mathcal{O}_n)\),
  are algebraic integers of \(H_n\) with \(\sigma(G_\mathfrak{c}) = G_{\mathfrak{c}(\sigma)^{-1}\mathfrak{c}}\),
  \(\overline{G_\mathfrak{c}} = G_{\mathfrak{c}^{-1}}\), and
  \(\sum_\mathfrak{c}\chi(\mathfrak{c})\log|G_\mathfrak{c}| = -12\,L'(0,\chi)\) for every nontrivial
  \(\chi\), where \(L(s,\chi)\) is the Epstein \(L\)-function of the level (Paper II, Thm 2.6).
- \(v_\mathfrak{P}(G_\mathfrak{c})\) is independent of \(\mathfrak{c}\) (Paper II, Thm 4.11), so
  \(\mathcal{V}_n := \langle G_\mathfrak{c}/G_{\mathfrak{c}'}\rangle = \langle G_\mathfrak{c}/G_1 : \mathfrak{c} \ne 1\rangle
  \subset \mathcal{O}_{H_n}^\times\) is a Galois-stable group of units
  (schmidt-euler-system.md, Cor. 0.1).
- \(L(s,\chi) = L_{\mathrm{prim}}(s,\chi)\,C_\chi(s)\) with \(C_\chi \equiv 1\) when \(\chi\) is
  nontrivial on every kernel \(\mathrm{Pic}(\mathcal{O}_n) \to \mathrm{Pic}(\mathcal{O}_{n/p})\)
  (Paper II, Lemma 6.3), and \(C_\chi(0)\) along pullbacks is the Eisenstein multiplier
  \(P_\ell(\chi)\) (schmidt-euler-system.md, Thm 2).
- Cubic layer: \([\mathcal{O}_{L_3}^\times : \langle -1, \theta_u\rangle] = 8h_{L_3}C_n(0)\)
  (Paper II, Thm 6.6; multipliers along the tower proved in schmidt-euler-system.md, Thm 3).
- Quadratic layer: \(\sum\chi_2\log|G| = -12\,m_\chi\log\varepsilon_{d_2}\),
  \(m_\chi = \tfrac{2h(d_1^*)}{w(d_1^*)}h(d_2)C(0)\) (Paper II, Prop. 6.1).

## The expected theorem (prove it)

\(H_n\) is totally complex of degree \(2h\), \(h = h(-4n^2)\), with unit rank \(h - 1\);
its \(h\) complex places correspond to \(\mathrm{Gal}(H_n/K) = \mathrm{Pic}(\mathcal{O}_n)\)
(one representative per coset of the complex conjugation \(\iota\)). Since
\(\log|\sigma(G_\mathfrak{c}/G_1)| = \log|G_{\mathfrak{c}(\sigma)^{-1}\mathfrak{c}}| - \log|G_{\mathfrak{c}(\sigma)^{-1}}|\),
the regulator of \(\mathcal{V}_n\) (basis \(G_\mathfrak{c}/G_1\), \(\mathfrak{c} \ne 1\); normalization
\(2\log|\cdot|\) at complex places, one place omitted) is a Dedekind group
determinant (Washington, *Introduction to Cyclotomic Fields*, Lemma 5.26:
\(\det(f(\sigma\tau^{-1}) - f(\sigma))_{\sigma,\tau\ne1} = \prod_{\chi\ne1}\sum_\sigma\chi(\sigma)f(\sigma)\)),
hence
$$
R(\mathcal{V}_n) \;=\; \prod_{\chi \ne 1}\Bigl|\sum_\mathfrak{c}\chi(\mathfrak{c})\,2\log|G_\mathfrak{c}|\Bigr|
\;=\; 24^{\,h-1}\prod_{\chi\ne1}|L'(0,\chi)| .
$$
The class number formula at \(s = 0\) for \(H_n\), with
\(\zeta_{H_n}(s) = \zeta_K(s)\prod_{\chi\ne1}L_{\mathrm{prim}}(s,\chi)\) and \(\zeta_K(0) = -\tfrac14\),
gives \(h_{H_n}R_{H_n} = \tfrac{w_{H_n}}{4}\prod_{\chi\ne1}L'_{\mathrm{prim}}(0,\chi)\). Therefore,
since an index of full-rank lattices is the ratio of covolumes,
$$
\boxed{\;\bigl[\mathcal{O}_{H_n}^\times : \mu(H_n)\,\mathcal{V}_n\bigr]
\;=\; \frac{24^{\,h-1}\prod_{\chi\ne1}|L'(0,\chi)|}{R_{H_n}}
\;=\; \frac{4\cdot 24^{\,h-1}}{w_{H_n}}\; h_{H_n}\prod_{\chi\ne1}C_\chi(0)\;}
$$
— the Robert index of the Schmidt \(\Delta\)-units is the class number of the ring
class field times an explicit power of 2 and 3 and the imprimitivity multipliers:
the Kubert–Lang shape of Paper II, Conjecture 6.7, for every \(n\). Anchors that
can be checked by hand before any computation: at \(n = 3\), \(H_3 = \mathbb{Q}(\zeta_{12})\),
\(w = 12\), \(h = 1\), \(\mathcal{V}_3 = \langle\pm\varepsilon_{12}^{4}\rangle\) (from
\(\sum\chi_2\log|G| = -4\log\varepsilon_{12}\) and \(G_1G_t = 3^6\)), and the unit
index of \(\mathbb{Q}(\zeta_{12})\) over \(\mathbb{Q}(\sqrt3)\) is 2, so the index is
\(4\cdot2 = 8 = 4\cdot24/12\); at \(n = 5\), \(H_5 = \mathbb{Q}(i,\sqrt5)\), \(w = 4\),
\(\mathcal{V}_5 = \langle\pm\varepsilon_5^{24}\rangle\), index \(24\,h_{H_5}\) (with
\(h_{H_5} = Q\), the unit index of the CM field, consistent either way).

Do the proof carefully: (i) the correspondence places \(\leftrightarrow\) Pic and the
regulator normalization; (ii) the determinant lemma with the relabelling
\(\sigma \mapsto \sigma^{-1}\) (signs only); (iii) \(L'(0,\chi)\) real and nonzero
(\(L(s,\chi) = L(s,\bar\chi)\), \(L(1,\chi) \ne 0\)), so \(\mathcal{V}_n\) has full rank;
(iv) the class number formula with the exact \(\zeta_K(0)\); (v) the torsion:
the index is of \(\mu(H_n)\mathcal{V}_n\); record \(\mu(H_n)\) (\(\mu_4\) always;
\(\mu_{12}\) iff \(\sqrt{-3} \in H_n\), i.e. iff \(K(\sqrt3)\) is a genus subfield,
e.g. \(n = 3, 9, 15\); \(\mu_8\) iff \(\sqrt2 \in H_n\)).

## Mission

1. **Prove the theorem** above as a statement about all \(n \ge 2\), in a new
   document `robert-index-full.md` in the house style (statement–proof, verification
   section, outlook). Derive the layer indices from it or alongside it: the real
   quadratic layer \([\mathcal{O}_{\mathbb{Q}(\sqrt{d_2})}^\times : \langle -1, \theta_u^{(2)}\rangle] = 6m_\chi\)
   (from \(\log|\theta_u^{(2)}| = \tfrac12\sum\chi_2\log|G|\)), the cubic layer of Paper II,
   and the CM sextic subfield (unit index \(Q \in \{1,2\}\) over \(L_3\)).
2. **Verify with PARI/GP.** Install with `apt-get install -y pari-gp` (apt works in
   this environment; record the gp version in the verification section) and call
   `gp -q` from a script `scripts/robert_index_full.py --selftest`. For each level
   \(n \in \{3, 5, 7, 9, 11, 13, 15\}\) (degrees \(4, 4, 8, 12, 12, 12, 16\)) and, if feasible,
   \(n = 23\) (degree 24): build \(H_n\) as `polredbest(polcompositum(x^2+1, polclass(-4*n^2))[1])`
   (or via `bnrclassfield` of the ring class group, whichever gives the smaller
   polynomial), run `bnfinit(pol, 1)` at `\p 120` for \(h_{H_n}\), \(R_{H_n}\) and
   `bnf.tu` (this gives \(w\)), and `bnfcertify` where it finishes — label every
   quantity GRH-conditional if it does not. Then certify, with the absolute-error
   criterion and \(\ge 40\) spare digits, that \(24^{h-1}\prod_{\chi\ne1}|L'(0,\chi)|/R_{H_n}\)
   is an integer (the \(L'\) from `phase_klf.py`'s independent evaluation, or from the
   character sums of the \(G_\mathfrak{c}\) of `schmidt_euler_system.py`, which is the same
   number by Thm 2.6), and that it equals \((4\cdot24^{h-1}/w)\,h_{H_n}\prod_\chi C_\chi(0)\)
   with \(\prod C_\chi(0)\) computed from Lemma 6.3 / Thm 2 of the Euler document
   (e.g. \(4\) at \(n = 9\), \(1\) at \(n = 11, 13\)). Also verify \(R(\mathcal{V}_n)\) directly
   from the \(h\) numbers \(\log|G_\mathfrak{c}|\) as a determinant, against the character-sum
   product — this checks the determinant lemma and the normalization independently.
3. **Express the units** (optional but valuable): with `bnfisunit`, write
   \(G_\mathfrak{c}/G_1\) in PARI's fundamental-unit basis at \(n = 9, 11, 13\) (the elements
   must be given exactly; use the certified integer polynomials \(D_n\) of Paper II
   Thm 2.6 to identify \(G_\mathfrak{c}\) as roots inside \(H_n\), e.g. via `nfroots`), and
   read off the \(p\)-parts of the index and which fundamental units are missed.
4. **The hyperbolic sextic layer at \(n = 21\)** (stretch): discriminant \(-440\),
   \(\mathrm{Pic} \cong \mathbb{Z}/2\times\mathbb{Z}/6\); for the odd characters of order 6 the
   coset products of the units \(R_f\) (Paper II, Thm 4.2) over \(\ker\chi\) lie in the
   subfield of the ring class field cut out by \(\ker\chi\); identify that field with
   PARI, compute its unit group, and the index of the coset units — the first
   genuinely hyperbolic index datum (Paper II, Remark 6.8).
5. **Port to Paper II**: replace Conjecture 6.7 by the theorem, add the verification
   table (levels, \(w\), \(h_{H_n}\), \(R_{H_n}\), index, GRH status), and re-run both
   selftests of Paper II before touching the text.

## Guard rails

Precision after imports; absolute-error certification with spare digits for every
integer read off from a real number; no PSLQ; exact arithmetic for the class
groups and characters (reuse `ClassGroup` of `schmidt_euler_system.py`); state
plainly which PARI outputs are GRH-conditional and whether `bnfcertify`
succeeded; keep the proved / certified labels distinct; do not claim the theorem
from the numerics — the numerics verify the proof. Update CLAUDE.md (document
map, theorem 22, ledger, scripts table), PROGRAM.md (item 2 status) and
outlook.md; commit with clear messages and push.

## What the session found beyond the prompt (for the record)

- The sextic \(F = K(\theta_u)\) is totally complex but **not CM** (\(L_3\) has a complex
  place), so "unit index \(Q\) over \(L_3\)" is not finite; the correct statement is the
  layer theorem (Theorem 2), \([\mathcal{O}_F^\times:\mu_4\mathcal{V}^A] = 576\,h_F\,C_n(0)^2\).
- Two lemmas were needed that the prompt did not foresee: primitive level = conductor
  (closing a gap in the use of Paper II's Lemma 6.3) and the Hecke recursion as an
  identity of Dirichlet series (positivity and integrality of the multipliers).
- The exact exponent matrices show that modulo \(24\)-th roots the index is exactly
  \(h_{H_n}\prod C_\chi(0)\), with \(\mathcal{O}^\times/\mathcal{W}_n \cong \mathrm{Cl}(H_n)\) at \(n = 13, 23\)
  (Conjecture 5.3 there, Gras type).
- The hyperbolic index has a general theorem (Theorem 3) with a genuine \(2\)-adic
  invariant \(Q^- = [E:E^+E^-]\); \(24^6\) at \(n = 21\), and \(27648\) in its sextic layer.
