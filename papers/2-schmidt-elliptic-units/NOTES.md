# Literature-diligence record for Paper II

Date of pass: 2026-09-04 (session of the Paper II draft). Searches were run
through a web-search relay; full texts at arXiv, Springer, Numdam, EuDML
and zbMATH were not reachable from this environment, so the statements
below distinguish what was verified from abstracts, reviews and citing
papers (bibliographic data, statements of main theorems as quoted by
others) from what could be read verbatim (nothing in full). This
complements the Paper I record in
[papers/1-schmidt-circles/NOTES.md](../1-schmidt-circles/NOTES.md) and the
probes of [spectral-geometry.md](../../spectral-geometry.md) §14.

## The classical elliptic-unit literature, as positioned in the paper

1. **C. L. Siegel**, *Lectures on advanced analytic number theory*, notes
   by S. Raghavan, Tata Institute, Bombay 1961 (lectures 1959/60). Kronecker
   limit formulas and their use in class field theory of imaginary
   quadratic fields. Bibliographic data verified.
2. **K. Ramachandra**, *Some applications of Kronecker's limit formulas*,
   Ann. of Math. (2) 80 (1964) 104–148. Construction of units (Ramachandra
   units) generating ray class fields via Siegel functions. Verified.
3. **G. Robert**, *Unités elliptiques*, Bull. Soc. Math. France, Mém. 36
   (1973). Defines the elliptic units of an abelian extension of an
   imaginary quadratic field and proves the index in the full unit group is
   the class number up to an explicit factor (from reviews/citations).
   Verified.
4. **R. Gillard, G. Robert**, *Groupes d'unités elliptiques*, Bull. SMF 107
   (1979) 305–317: two groups of elliptic units per abelian extension, index
   related to the class number by a Leopoldt-type formula (from the
   abstract). Verified.
5. **D. Kubert, S. Lang**, *Modular units*, Grundlehren 244 (1981), ch.
   11–13: Siegel/Klein forms, elliptic units, index formulas. Cited for the
   multiplier bookkeeping (open problem of §5) and the shape of the index
   conjecture (§6).
6. **H. Oukhaba**, *Index formulas for ramified elliptic units*, Compositio
   137 (2003) 1–22: Sinnott-type index formulas. Verified.
7. **H. M. Stark**, *L-functions at s = 1. IV. First derivatives at s = 0*,
   Adv. Math. 35 (1980) 197–235: the conjecture over an imaginary quadratic
   base is a theorem through elliptic units; for ring class fields of prime
   conductor the Stark unit (with the roots of unity) has index equal to
   the class number of the ring class field (as quoted by 8 and 9).
   Verified via 8/9.
8. **F. Hajir, F. Rodriguez Villegas**, *Explicit elliptic units, I*, Duke
   Math. J. 90 (1997) 495–521. Explicit elliptic units in ring class fields
   of prime conductor and their relation to Stark units. Verified.
9. **Ö. Küçüksakallı**, *Class numbers of ring class fields of prime
   conductor*, Acta Arith. 153 (2012) 251–269: uses Stark's elliptic units
   (index = class number of the ring class field) to compute class numbers
   of ring class fields over imaginary quadratic fields. Verified from the
   abstract. **This is the closest prior work**: an index-equals-class-number
   statement for elliptic units in ring class fields. It concerns Stark's
   units in the full ring class field of prime conductor; it does not
   concern the Δ-coset units of the Schmidt data, cubic subfields, or Q(i)
   specifically.
10. **R. Schertz**, *Complex multiplication*, New Math. Monographs 15, CUP
    2010, ch. 6–7: elliptic units and class numbers, Stark units of ring
    class fields, N-systems. Verified (table of contents).
11. **B. H. Gross**, *On canonical and quasi-canonical liftings*, Invent.
    Math. 84 (1986) 321–326. Verified. Used only as the conceptual reading
    of the per-class valuation law (denominators p^{k-1}(p+1), 2^k).
12. **B. H. Gross, D. B. Zagier**, *On singular moduli*, Crelle 355 (1985)
    191–220. Verified (Paper I). Reference for the collision sets GZ(D,-3),
    GZ(D,-4) and for the open exponent law (outlook §2.8).
13. **E. Friedman**, *Analytic formulas for the regulator of a number
    field*, Invent. Math. 98 (1989) 599–622. **Constant verified**: the
    regulator of every number field of degree ≥ 2 satisfies R > 0.2052,
    the minimum being attained at a unique sextic field of discriminant
    -10051 (as quoted in several recent regulator papers; the paper also
    determines all fields with R/w < 1/8). The script and the paper use
    R > 0.2052 as an unconditional cap on the residual index of the
    root descent (largest ratio R_L/0.2052 < 31.2 at n = 23, below the
    tested k ≤ 40).
14. **C. Meyer**, *Die Berechnung der Klassenzahl abelscher Körper über
    quadratischen Zahlkörpern*, Akademie-Verlag 1957: class numbers of
    abelian extensions of quadratic fields via Kronecker limit formulas.
    Verified.
15. **M. Vlasenko, D. Zagier**, *Higher Kronecker "limit" formulas for real
    quadratic fields*, Crelle 679 (2013) 23–64. Verified. Closest relatives
    of log|u_f| on the real-quadratic side; no overlap (Paper I's DIT
    non-relation covers the cycle-integral side).
16. **J. Sturm**, LNM 1240 (1987) 275–280 (Sturm bound); **B. Schoeneberg**,
    Math. Ann. 116 (1939) 511–523 (theta series of binary forms as weight-1
    forms); **Diamond–Shurman** §4.8 (weight-1 Eisenstein series with two
    characters, vanishing constant term). Verified; used in the proof of the
    genus factorizations.

Searches for "Schmidt arrangement" together with "elliptic units",
"Kronecker limit", "Robert index" return nothing relevant (only generic
Kronecker-limit-formula literature). Searches for elliptic-unit index
computations in ring class fields of Q(i), or for cubic subfields of ring
class fields with index 8h, return nothing beyond items 7–9 (general
imaginary quadratic base, Stark units, prime conductor).

## Verdicts, section by section

- **§2 (master identities).** The Kronecker limit formula itself is
  classical (1, 14). NEW: that the geometric phase of the Schmidt
  arrangement is classwise a twisted Δ-quotient whose character sums are
  exactly -2L'(0,χ) (Euclidean) / -4L'(0,χ) on odd characters (hyperbolic)
  plus GZ dressings; the cancellation of ε = n + sqrt(n²-1) and μ; the
  Δ-mass polynomial D_n as an elliptic-unit system for the ring class tower
  of Q(i).
- **§3 (genus characters).** The factorization of an Epstein L-function of
  a genus character into two Dirichlet L-functions is classical
  (Dirichlet/Kronecker genus theory); the conductor Euler corrections at
  square conductors and the identification "the real field is Q(sqrt n)"
  in the Euclidean aspect, and the exact GZ-supported Σ-factorizations,
  are new as statements about the arrangement; the Sturm-bound proof is
  standard technique.
- **§4 (unit theorem, per-class law).** NEW. Elliptic units of ring class
  fields are classical (3–10), but the statement that the Δ-quotient along
  the invertible ambiguous twist r of a *non-maximal* order is a unit on
  every stratum, proved by rigidity of p-divisible-group quotients without
  quasi-canonical liftings, and the class-independent valuation w_p(k) of
  the Euclidean Δ-data, have no antecedent we could find. The Deuring/Gross
  reading (11) is cited as interpretation only.
- **§5 (first-power units).** NEW as an object; the multiplier problem is
  the classical Kubert–Lang/Schertz bookkeeping (5, 10), stated as open.
- **§6 (Robert index).** The shape "index = class number × explicit
  factor" is Robert–Gillard–Robert–Kubert–Lang–Stark (3–5, 7); the
  prime-conductor ring-class statement "index = h" for Stark's units is
  known (7–9). NEW: the Δ-coset units θ_u of the Schmidt data, the explicit
  index 8 h_{L_3} C_n(0) in the unit group of the real cubic subfield over
  Q(i), certified at eight levels with the out-of-sample n = 23, and the
  pullback consistency (identical fundamental units, Euler multipliers 2,
  2, 2, 4). The relation between θ_u and the Stark unit of H_n is stated
  as an open comparison, not claimed.

## Corrections and additions made while writing

- schmidt-units.md §5.2 listed the Δ-mass at n = 18 as 2^18·3^48; the
  Δ-mass law (Paper I, Thm 7.18; euclidean-moduli-invariants.md Thm 4)
  gives 2^36·3^48 (p = 2: 3·N_e(9) = 36). Corrected in schmidt-units.md and
  used in the paper; nothing else depends on it (|M(18)|^{1/3} = 2^12·3^16
  is an integer either way).
- The genus closed form L'(0,χ_2) = (2h(d_1)/w(d_1)) h(d_2) log ε_{d_2}
  C(0) of phase-kronecker-limit.md §4(a) must be read with h, w of the
  *fundamental* discriminant d_1^* underlying the Kronecker symbol (d_1/·):
  at n = 9, d_1 = -27, d_1^* = -3, w = 6, giving (1/3)·4·log ε_12 =
  (4/3) log ε_12 as tabulated. The paper states this explicitly (Prop. 3.1).
- The paper adds a short lemma (Lemma 6.3) proving that the Epstein
  L-function of a ring class character equals the primitive Hecke
  L-function whenever the character is nontrivial on every kernel
  Cl(O_n) → Cl(O_{n/p}), i.e. C_n ≡ 1 at the primitive cubic levels; the
  source document asserted this without proof. At the pullback levels
  C_n(0) is read off from the certified integer L'(0,χ_3)/R_{L_3} = h C_n(0).

## Standard references fixed for the bibliography

Cox (ring class fields), Lang (Elliptic functions, ch. 12 and 20),
Neukirch (class number formula at s = 0), Diamond–Shurman (weight-1
Eisenstein series), Sturm, Schoeneberg, Siegel, Ramachandra, Robert,
Gillard–Robert, Kubert–Lang, Oukhaba, Stark, Hajir–Rodriguez Villegas,
Küçüksakallı, Schertz, Gross, Gross–Zagier, Lauter–Viray, Friedman, Meyer,
Vlasenko–Zagier, Schmidt, Stange (TAMS, IMRN), mpmath, SymPy, and Paper I
as [SchmidtI].

## What Paper II claims as new (final list, mirrored in the introduction)

1. The two Kronecker limit formulas for the phase (Theorems 2.5, 2.9) with
   the trivial-character anchor, the Δ-mass polynomial and the uniform
   Stark law (Theorem 2.6), and the twisted-ratio polynomial (Theorem 2.10).
2. The genus closed forms (Propositions 3.1, 3.6; Theorem 3.4), the field
   Q(sqrt n), the exact Σ-factorizations with GZ split-prime support, and
   the certified non-fits for non-real characters.
3. The unit theorem on all strata (Theorem 4.2) with the lattice lemma and
   the principality criterion; the per-class valuation law (Theorem 4.11)
   with the split ladder (Corollary 4.13) and the inert-ladder
   decompositions.
4. The first-power Schmidt units w_f, their exact laws, the η-quotient form,
   the coherence table and w = -(2 ± sqrt 3) at n = 7.
5. The Robert index 8 h_{L_3} C_n(0) on the cubic layer (Theorem 6.6), with
   Propositions 6.2, 6.4 proved and the fundamental units rigorous through
   Friedman's bound; the quadratic-layer projections (Proposition 6.1).
