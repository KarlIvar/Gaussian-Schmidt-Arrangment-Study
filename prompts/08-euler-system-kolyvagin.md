# Prompt 08 — The Schmidt Δ-system as an Euler system: Kolyvagin's bound and the full Robert index

Read `CLAUDE.md`, then [schmidt-euler-system.md](../schmidt-euler-system.md)
(all of it), [schmidt-units.md](../schmidt-units.md) §3 and §5, and Paper II
(`papers/2-schmidt-elliptic-units/`) §6. Machinery:
`scripts/schmidt_euler_system.py` (exact class groups of \(\mathcal{O}_n\),
projections, ideal twists, coset units along the tower),
`scripts/schmidt_units.py` (fundamental units, indices).

## Established (do not redo)

Theorem 1 of schmidt-euler-system.md: the \(\Delta\)-data
\(G_\mathfrak{c} = n^{12}\Delta(\Lambda_\mathfrak{c})/\Delta(\mathbb{Z}[i])\) satisfy
Heegner-type norm relations with the Eisenstein eigenvalue \(\ell+1\) along
\(H_n \subset H_{n\ell}\); Theorem 2: the Hecke recursion for \(L'(0,\chi)\) under
pullback (all imprimitive Euler multipliers proved); Theorem 3: the cubic
Robert index multiplies by \(P_\ell(\chi_3)\) along the tower (chains from
\(9, 11, 13, 23\) computed to level 81); Lemma 5.1: the Kolyvagin derivative
classes \(\kappa_{n,\ell}(v) \in H_n^\times/p\) exist for inert
\(\ell \equiv -1 \bmod p\), \(p \ge 5\).

## Mission

1. **Kolyvagin's bound.** Carry out the localization analysis of
   \(\kappa_{n,\ell}(v)\) at the primes of \(H_n\) above \(\ell\) (finite–singular
   comparison; the singular part is governed by \(w \bmod \ell\)) and prove the
   anticyclotomic Thaine–Rubin-type bound: the \(p\)-part of a suitable
   \(\chi\)-component of \(\mathrm{Cl}(H_n)\) is bounded by the corresponding
   component of \([\mathcal{O}_{H_n}^\times : \mathcal{V}_n]\). Follow Gross's
   *Kolyvagin's work on modular elliptic curves* and Rubin's *Euler Systems*
   Ch. 3–4 for the structure; note the system is of Heegner shape, not of
   the \((1 - \mathrm{Fr}^{-1})\) shape. Test the bound numerically against
   class numbers of ring class fields of \(\mathbb{Q}(i)\) of prime conductor
   (Küçüksakallı, Acta Arith. 153 (2012)); PARI/GP (`bnfinit`) may be
   installed with apt for the class-number side — record any such use in the
   verification section.
2. **The full index** is now its own handoff,
   [09-full-robert-index.md](09-full-robert-index.md) (the analytic index
   formula \([\mathcal{O}_{H_n}^\times : \mu\mathcal{V}_n] = (4\cdot24^{h-1}/w)\,h_{H_n}\prod C_\chi(0)\)
   to prove, and its PARI verification); do it first if both are in one session,
   since the Kolyvagin bound is stated in terms of that index.
3. **Stark comparison.** Express Stark's unit of \(H_n\) (prime conductor
   \(n = 11, 13, 23\); Hajir–Rodriguez Villegas' explicit form) in terms of
   \(\mathcal{V}_n\), or show it is not in \(\mathcal{V}_n\); relate the two
   indices \(h_{H_n}\) and \(8h_{L_3}\).
4. Port Theorems 2–3 into Paper II §6 (a short subsection) and upgrade
   Remark 3.2, Prop. 3.1's correction, Theorem 6.6's \(C_n(0)\) column and
   Conjecture 6.7's pullback clause from certified to proved.

## Guard rails

As in CLAUDE.md: precision after imports, absolute-error certification with
spare digits, exact arithmetic wherever finite, no PSLQ-derived claim
without exact re-verification; every theorem environment keeps its
proved/certified label. Update CLAUDE.md (document map, theorems, ledger,
scripts) and outlook.md §6.
