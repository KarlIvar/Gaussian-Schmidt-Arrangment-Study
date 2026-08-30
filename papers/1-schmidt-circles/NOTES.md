# Literature-diligence record for Paper I

Date of pass: 2026-08-30 (session of the paper draft). Searches were run
through a web-search relay (arXiv full text was not directly reachable from
this environment; statements below distinguish what was verified from
abstracts/summaries vs. what could be read verbatim). This complements the
earlier probes recorded in [spectral-geometry.md](../../spectral-geometry.md)
§14 (arXiv metadata, citation graphs of Stange's papers, full-text probes),
which found exactly three published papers mentioning "Schmidt arrangement"
and none touching the level stratification, the involution class formula, or
the phase.

## The Schmidt-arrangement literature (complete, to our knowledge)

1. **A. L. Schmidt**, *Diophantine approximation of complex numbers*, Acta
   Math. 134 (1975) 1-85. The circles appear as the Farey-type subdivision
   underlying Schmidt's complex continued fractions; no classification, no
   counting.
2. **K. E. Stange**, *The Apollonian structure of Bianchi groups*, Trans.
   AMS 370 (2018) 6169-6219 (arXiv:1505.03121). Defines Schmidt arrangements
   S_K for imaginary quadratic K, proves tangency structure (for K != Q(sqrt-3)
   circles meet only tangentially), decomposes S_K into K-Apollonian packings.
   Its Theorem 1.3 (cited as such by Rickards-Stange 2026) identifies the
   scaled Gaussian Schmidt arrangement with the Apollonian super-packing of
   Graham-Lagarias-Mallows-Wilks-Yan, containing one copy of every primitive
   integral Apollonian packing up to symmetry.
3. **K. E. Stange**, *Visualising the arithmetic of imaginary quadratic
   fields*, IMRN 2018, no. 12, 3908-3938 (arXiv:1410.0417). Proves: the
   curvatures of K-Bianchi circles are integer multiples of sqrt(-Delta);
   tangent curvatures are governed by the norm form; **circles of curvature
   f sqrt(-Delta), up to translation and rotation by 180 degrees, are in
   bijection with ker(Pic(O_f) -> Pic(O_K))**; S_K is connected iff O_K is
   Euclidean.
4. **J. Rickards, K. E. Stange**, *Eisenstein circle packings and the
   Eisenpint Schmidt arrangement*, arXiv:2605.16053 (May 2026). Eisenstein
   case: congruence obstructions, strong approximation, density-one
   local-global, reciprocity obstructions. Recalls the Gaussian congruence
   description of S_{Q(i)} (curvature 2p', center 2t' + (2s'+1)i with a
   divisibility condition of exactly our Theorem-2.1 shape).

Adjacent: **D. E. Martin**, *A geometric study of circle packings and ideal
class groups* (arXiv:2202.10530), circle-packing avatars of class groups in
Schmidt-arrangement-like settings; **GLMWY**, *Apollonian circle packings:
geometry and group theory II. Super-Apollonian group and integral packings*,
Discrete Comput. Geom. 35 (2006) 1-36 (strongly integral super-packings:
integral curvatures AND curvature x centers; eight geometric strongly
integral super-packings; contains every integral Apollonian packing);
expository: Stange's short exposition *Visualizing imaginary quadratic
fields*; Fuchs' surveys; the AMS Visual Insight post (2015).

## Verdicts, section by section

- **Section 2 (classification).** NOT claimed new. The congruence
  classification of S_{Q(i)} is equivalent to GLMWY's description of the
  strongly integral super-packing combined with Stange's identification
  (Trans. AMS Thm 1.3), and a congruence form is stated in Rickards-Stange
  (2026). Curvature 2Z and the Hermitian-matrix framework are Stange's. The
  paper presents the classification with a short self-contained proof
  (descent), the single-equation form u^2 + v(v+1) = nm, and the explicit
  PSL_2-vs-PGL_2 parity split S vs iS, and says exactly this in the text.
- **Section 3(a) (Euclidean count).** The bijection behind N_e(n) = 2h(-4n^2)
  is Stange's IMRN theorem specialized to K = Q(i) (where the Pic-kernel is
  the whole ring class group); the paper credits this prominently. What we
  add: the explicit primitive-sublattice dictionary with conductor =
  curvature (used heavily later), the multiplicative closed form
  N_e(n) = n prod(1 - chi(p)/p) with the mean-value asymptotic
  sum N_e ~ X^2/(2G) (Catalan's constant) -- elementary, and we found no
  prior statement; presented as a proposition with proof.
- **Section 3(b) (hyperbolic level count 3H(n^2-1)).** NEW. No paper
  stratifies a Schmidt arrangement by the inversive level; H(n^2-1) is not
  the answer to any indexed counting problem in the literature (searches on
  "Hurwitz class number n^2-1", hyperbolic-circle counts at CM points,
  arccoth radii: nothing).
- **Section 3(c) (spherical census 4H(4(l^2+1))).** NEW as a census of a
  Schmidt arrangement; reduces to Gauss's r_3 = 12 H(4n) (cite Gauss,
  Grosswald) by a parity argument.
- **Section 4 (involution and class formula).** NEW. No Cartan-embedding
  involution on circle classes in the literature; the class formula
  sigma-hat[f] = [r_n][f]^{-1} appears nowhere. The closed-geodesic remark
  connects to Sarnak's reciprocal geodesics (Clay Proc. 7 (2007)) --
  cited as the 2D analogue of trace-(-2n) classes; the 3D twisted version is
  new (spectral-geometry.md gap 5).
- **Section 5 (composition in circle language).** The underlying algebra is
  Dirichlet/Gauss composition (classical; cite Cox); the circle-geometric
  recipes (inverse = mirror, composition = CRT/magnification, 2-torsion =
  mirror symmetry, the matrix read-off W_X) are new as statements about the
  arrangement.
- **Section 6 (phase).** NEW. The six-coordinate system, the laws, the
  first-power descent u = Phi_y/Phi_x at the Heegner pair, the level
  polynomials and the GZ support theorem have no antecedent. The
  Gross-Zagier support input is *On singular moduli* (Crelle 355 (1985)
  191-220). The certified non-relation with Duke-Imamoglu-Toth cycle
  integrals (Ann. of Math. 173 (2011) 947-981; Katok-Sarnak, Israel J. Math.
  84 (1993) for the two-sign framework) is the paper's novelty statement,
  with the PSLQ battery of scripts/dit_comparison.py re-run for this draft.
- **Section 7 (Euclidean phase).** The j-value theorem (H_{-4n^2}^2), trace
  slice t(4n^2) (cite Zagier, *Traces of singular moduli*, 2002), the
  lemniscatic phase, monic level polynomials, all-n irreducibility
  (Theorem 5) and the Delta-mass law are NEW; the input bijection is
  Stange's as in 3(a).

## Standard references fixed for the bibliography

Cox (Primes of the form x^2+ny^2, 2nd ed., Wiley 2013) for ring class
fields, reduced forms, composition, Cl <-> Pic dictionaries, Thm 11.36
(reciprocity for j of an order); Elstrodt-Grunewald-Mennicke (Groups Acting
on Hyperbolic Space, Springer 1998) for Bianchi groups and binary Hermitian
forms; Shimura (Introduction..., Princeton 1971) for the main theorem of CM
and weight-2 algebraicity; Lang (Elliptic Functions, 2nd ed., GTM 112) for
Delta-quotients and isogenies; Zagier (CRAS 281 (1975)) and Cohen (Math.
Ann. 217 (1975)) for Hurwitz class numbers and the weight-3/2 series;
Grosswald (Sums of Squares, 1985) for r_3 = 12H; Duke-Imamoglu-Toth (Ann.
Math. 173 (2011)); Katok-Sarnak (1993); Sarnak (Reciprocal geodesics, Clay
Math. Proc. 7 (2007) 217-237); Kubert-Lang (Modular Units, 1981) --
mentioned for the Siegel-unit route that the first-power descent avoids;
GLMWY (DCG 35 (2006) 1-36); Rickards-Stange (arXiv:2605.16053); A. L.
Schmidt (Acta Math. 134 (1975)); mpmath and SymPy for the verification
appendix.

## What Paper I claims as new (final list, mirrored in the introduction)

1. The level stratification of S and the three-geometry census triptych:
   the hyperbolic 3H(n^2-1) (Thm in section 3) and the spherical
   4H(4(l^2+1)) (section 3); the Euclidean member restates Stange + the
   closed form and Catalan asymptotics.
2. The identification of sigma(X) = conj(X)^{-1} (unitary involution,
   Cartan embedding = inversion, trace = -2 alpha) and the class formula
   sigma-hat[f] = [r_n][f]^{-1} with its proof (Lemma A's closed-form
   section included).
3. Composition and inversion of Schmidt circles in circle language.
4. The six-invariant moduli system with the phase, its functional
   equations, the Heegner-pair structure on X_0((n-1)/2) x_{X(1)}
   X_0((n+1)/2) with discriminant coupled to level.
5. The first-power descent u_f = Phi_y/Phi_x(beta_1, beta_2) and the
   dihedral Galois law at first power; integer level polynomials Q_n,
   irreducible at every computed level; the GZ support theorem for
   denominators; the divisor-class sign law.
6. The certified DIT non-relation (novelty statement).
7. The Euclidean phase theory: H^2 j-products, lemniscatic normalization,
   monic integer level polynomials irreducible for every n (Theorem 5),
   the Delta-mass law (Theorem 4), the reality/center criterion.

Where each claim is proved in the repo: see the Machine verification
appendix of the paper and the source documents cited beside each theorem.
