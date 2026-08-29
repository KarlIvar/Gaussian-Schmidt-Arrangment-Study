# Euclidean moduli invariants of Schmidt disks: ring class fields of \(\mathbb{Q}(i)\), a sixth invariant, and lemniscatic phase units

This is the Euclidean companion of [moduli-invariants.md](moduli-invariants.md). There we studied functions of \(X \in \mathrm{SL}_2(\mathbb{C})\) invariant under \(X \mapsto \gamma X \gamma'\) for \(\gamma, \gamma' \in \mathrm{SL}_2(\mathbb{Z})\) — the *hyperbolic* moduli problem, whose values on the Bianchi group turned out to be singular moduli of discriminant \(1-n^2\) and the phase units \(u_f\). Here we replace the left factor by the **Euclidean symmetries of the plane inside the Bianchi group**: the translations
$$
N \;=\; \Bigl\{ T_\alpha = \begin{pmatrix} 1 & \alpha \\ 0 & 1 \end{pmatrix} : \alpha \in \mathbb{Z}[i] \Bigr\},
$$
and study functions of \(X\) invariant under
$$
X \;\longmapsto\; T_\alpha\, X\, \gamma', \qquad \alpha \in \mathbb{Z}[i],\ \gamma' \in \mathrm{SL}_2(\mathbb{Z}).
$$
Two invariants are visible at once: the Cartan-type matrix \(X \mapsto X\bar X^{-1}\) (equivalently the circle \(X(\hat{\mathbb{R}})\), which right multiplication fixes exactly and left translation merely translates), and \(X \mapsto j(-X^{-1}(\infty))\), defined whenever \(X(\mathbb{H})\) is a bounded disk. Together they have \(3 + 2 = 5\) real parameters against \(\dim_{\mathbb{R}}\mathrm{SL}_2(\mathbb{C}) = 6\). This document: (1) finds the missing sixth invariant — again a phase, with an even cleaner mechanism than the hyperbolic \(\Theta\); (2) computes the value set of \(j(-X^{-1}(\infty))\) over the Schmidt disks of curvature \(2n\) in the unit square — the answer is the **square of the ring class polynomial of the order \(\mathbb{Z} + n\mathbb{Z}[i]\)**; (3) does the same for the sixth invariant — after normalizing by a **lemniscatic period**, its values at a level are the roots of one monic integer polynomial, which is moreover **irreducible** at every computed level, as are all the level polynomials of this study (§5.6), while the total \(\Delta\)-quotient of a level obeys an exact product law, proved in §5.5 (Theorem 4); (4) closes with a research outlook. All numerics: [scripts/euclidean_moduli_invariants.py](scripts/euclidean_moduli_invariants.py) and [scripts/mass_law_and_irreducibility.py](scripts/mass_law_and_irreducibility.py) (mpmath, 120–840 digits, absolute-error certification; the proof ingredients of Theorem 4 and all factorizations in exact arithmetic).

Throughout, \(K = \mathbb{Q}(i)\), \(\mathcal{O}_K = \mathbb{Z}[i]\), and \(\mathcal{O}_n = \mathbb{Z} + n\mathbb{Z}[i]\) is the order of conductor \(n\) and discriminant \(-4n^2\), with class group \(\mathrm{Cl}(\mathcal{O}_n)\), class number \(h(\mathcal{O}_n)\), and ring class polynomial \(H_{-4n^2}(x) = \prod_{\mathfrak{c}}(x - j(\mathfrak{c})) \in \mathbb{Z}[x]\).

## 1. The Euclidean moduli problem and the five known invariants

**Normalization.** Let \(X = \begin{pmatrix} a & b \\ c & d \end{pmatrix} \in \mathrm{SL}_2(\mathbb{C})\) with \(c \neq 0\), and suppose \(X(\mathbb{H})\) is a **bounded disk** \(D\). Writing the boundary circle's Hermitian matrix as in [circle-classification.md](circle-classification.md), boundedness is equivalent to
$$
\kappa \;:=\; 2\operatorname{Im}(\bar c d) \;>\; 0
\qquad\Longleftrightarrow\qquad
w \;:=\; \frac{d}{c} \;=\; -X^{-1}(\infty) \;\in\; \mathbb{H},
$$
and then \(\kappa\) is the curvature of \(\partial D\). The pole \(X^{-1}(\infty) = -d/c = -w\) lies in the lower half-plane; \(w\) is its reflection through \(0\). For \(X \in \mathrm{SL}_2(\mathbb{Z}[i])\) this is the positive-curvature normalization of a **Schmidt disk**: curvature \(\kappa = 2n\), center \(\zeta/2n\) with \(\zeta = i(a\bar d - b\bar c) \equiv i \pmod 2\).

**The first invariant.** \(Z := X \bar X^{-1}\) satisfies \(Z\bar Z = 1\), and \(v \mapsto Z\bar v\) is the inversion in the circle \(\omega = X(\hat{\mathbb{R}}) = \partial D\) (as in moduli-invariants.md §1); so \(Z\) *is* the circle, with orientation. Under the two actions:
$$
Z(X\gamma') = Z(X) \quad (\gamma' \text{ real — exact invariance}), \qquad
Z(T_\alpha X) = T_\alpha\, Z(X)\, T_{-\bar\alpha}
$$
— the twisted conjugation by which the inversion in \(\omega\) becomes the inversion in \(\omega + \alpha\). So the honest bi-invariant content of \(Z\) is **the circle modulo \(\mathbb{Z}[i]\)-translation**: the curvature \(\kappa\) and the center mod \(\mathbb{Z}[i]\) — three real parameters. (For Schmidt disks: \(\kappa = 2n\) and \(\zeta \bmod 2n\mathbb{Z}[i]\), the classification data of [circle-classification.md](circle-classification.md).)

**The second invariant.** Left translation by *any* \(\alpha \in \mathbb{C}\) leaves the bottom row of \(X\) untouched, hence fixes \(w\) exactly. Right multiplication by \(\gamma' = \begin{pmatrix} p & q \\ r & s\end{pmatrix} \in \mathrm{SL}_2(\mathbb{Z})\) sends \((c, d) \mapsto (cp + dr,\, cq + ds)\), i.e.
$$
w \;\mapsto\; \frac{sw + q}{rw + p} = \gamma^\vee(w), \quad \gamma^\vee = \begin{pmatrix} s & q \\ r & p \end{pmatrix} \in \mathrm{SL}_2(\mathbb{Z}),
\qquad
c \;\mapsto\; c\,(p + rw).
$$
So \(w\) moves by an \(\mathrm{SL}_2(\mathbb{Z})\)-Möbius map and
$$
\beta(X) \;:=\; j(w) \;=\; j\bigl(-X^{-1}(\infty)\bigr)
$$
is a bi-invariant — two more real parameters, defined exactly when \(X(\mathbb{H})\) is bounded. Total: five.

**What the five see, and what they miss.** The double coset \(N X\,\mathrm{SL}_2(\mathbb{Z})\) of a generic \(X\) is completely described by two independent pieces of data:

- the **oriented lattice** \(\Lambda = \mathbb{Z}c + \mathbb{Z}d \subset \mathbb{C}\) (the bottom row modulo \(\pm\mathrm{SL}_2(\mathbb{Z})\)-basis change; orientation \(\operatorname{Im}(\bar c d) > 0\)) — 4 real parameters;
- the **center of the disk modulo \(\mathbb{Z}[i]\)** (the top row is determined by the bottom row up to \((a, b) \mapsto (a + tc, b + td)\), \(t \in \mathbb{C}\), which translates the disk by \(t\); the left action quotients \(t\) by \(\mathbb{Z}[i]\)) — 2 real parameters.

Of the lattice \(\Lambda\), the five invariants see the covolume (\(= \kappa/2\), via \(\kappa = 2|c|^2\operatorname{Im} w\)) and the homothety class (\(= j(w)\), as \(\Lambda = c\,(\mathbb{Z} + \mathbb{Z}w)\) up to ordering) — but **not the phase of the homothety** \(c\). That one missing real dimension is the sixth invariant.

## 2. The sixth invariant: the phase of the homothety

> **Definition.** For \(X\) as above (bounded-disk normalization), set
> $$
> \Theta(X) \;:=\; \frac{j'(w)}{c^2} \;=\; -\operatorname*{Res}_{z = X^{-1}(\infty)} X(z)\; \cdot\; j'\bigl(-X^{-1}(\infty)\bigr).
> $$
> (The residue of the Möbius map at its pole is \(\mathrm{Res}_{-d/c}\,\frac{az+b}{cz+d} = -1/c^2\); this is the coordinate-free form.) More generally \(j'\) may be replaced by any meromorphic weight-2 kernel \(g\) for \(\mathrm{SL}_2(\mathbb{Z})\) with real \(q\)-coefficients, as in moduli-invariants.md §5.

**Proposition (invariance).** \(\Theta(T_\alpha X \gamma') = \Theta(X)\) for every \(\alpha \in \mathbb{C}\) (not just \(\mathbb{Z}[i]\)) and every \(\gamma' \in \mathrm{SL}_2(\mathbb{Z})\).

*Proof.* Left translations do not touch \((c,d)\). On the right, by the transformation above and \(j'(\gamma^\vee w) = (rw + p)^2 j'(w)\) (weight 2), numerator and denominator of \(j'(w)/c^2\) pick up the same factor \((p + rw)^2\). \(\square\)

**Lattice reading.** Let \(\hat h_2\) be the weight-2 homogeneous extension of \(j'\) to oriented lattices: \(\hat h_2(\mathbb{Z}c + \mathbb{Z}d) := j'(d/c)/c^2\) (basis-independence *is* the right-invariance above; \(\hat h_2(\lambda\Lambda) = \lambda^{-2}\hat h_2(\Lambda)\)). Then
$$
\Theta(X) = \hat h_2(\Lambda_X), \qquad \Lambda_X = \mathbb{Z}c + \mathbb{Z}d :
$$
**the sixth invariant is the weight-2 value at the actual lattice**, not just at its homothety class — precisely the phase datum the five invariants discard. This is the exact structural analogue of the hyperbolic \(\Theta = j'(m_1)\,X'(m_2)\,\overline{j'(\bar m_2)}^{-1}\): there, two kernels (one per \(\mathrm{SL}_2(\mathbb{Z})\)-factor) glued by a derivative; here, one kernel glued by a residue — the left factor is a translation group, whose derivative is \(1\) and which needs no kernel.

**Proposition (fiber; \(\arg\Theta\) is the missing coordinate).** Fix the five invariants. The identity component of the fiber through \(X\) is \(X\,\{h_t\}\), where \(h_t \subset \mathrm{SL}_2(\mathbb{R})\) is the rotation group fixing the pole \(X^{-1}(\infty) = -w\) (equivalently the pair \(\{-w, -\bar w\}\)). Along it
$$
|\Theta(Xh_t)| \ \text{is constant}, \qquad \arg \Theta(Xh_t) = \arg\Theta(X) + 2t ,
$$
and \(|\Theta| = \operatorname{Im}(w)\,|j'(w)|\,/\,n\) is a function of the five invariants. Hence \((\kappa,\ \text{center mod } \mathbb{Z}[i],\ \beta,\ \arg\Theta)\) — six real dimensions — specify the double coset up to finite ambiguity.

*Proof.* A real \(h\) preserves \(\hat{\mathbb{R}}\), hence the circle exactly, and \(w(Xh) = -h^{-1}(-w)\); so the fiber condition is \(h(-w) = -w\), the elliptic subgroup rotating about \(-\bar w \in \mathbb{H}\). Its generator \(K\) has left eigenvector \((1, w)\) — for \(K_p\) fixing \(p \in \mathbb{H}\) the row eigenvectors are \((1, -\bar p)\) and \((1, -p)\), and \(-\overline{(-\bar w)} = w\) — with eigenvalue of modulus 1, so the bottom row \((c,d) = c(1,w)\) transforms as \(c \mapsto c e^{-it}\): \(w\) fixed, \(|c|\) fixed, \(\Theta = j'(w)/c^2 \mapsto e^{2it}\Theta\). The formula for \(|\Theta|\) follows from \(\kappa = 2n = 2|c|^2\operatorname{Im}(w)\), and \(y\,|j'|\) is \(\mathrm{SL}_2(\mathbb{Z})\)-invariant, hence a function of the \(j(w)\)-class. \(\square\)

(Compare the hyperbolic fiber: real rotations about the axis \(\{m_2, \bar m_2\}\), \(|\Theta|\) constant, \(\arg\Theta\) linear at rate \(\sqrt{\alpha^2 - 1}\). The Euclidean rate is the constant \(2\) — the weight.)

**The Borel refinement.** The full stabilizer of \(\infty\) in \(\mathrm{PSL}_2(\mathbb{Z}[i])\) is generated by \(N\) and \(z \mapsto -z\) (the class of \(\operatorname{diag}(i, -i)\)), which sends \(D \mapsto -D\) and \(\Theta \mapsto -\Theta\). So \(\Theta\) is invariant under the translations exactly and anti-invariant under the extra unit; \(\Theta^2\) is invariant under the whole Borel. This sign is not a defect — it will match the two-to-one structure of §3 perfectly.

**Remark (the involution \(\sigma\) changes sides).** In the hyperbolic study \(\sigma(X) = \bar X^{-1}\) acted on the double cosets. Here \(\sigma(T_\alpha X \gamma') = \bar\gamma'^{-1}\,\sigma(X)\,T_{-\bar\alpha}\): \(\sigma\) exchanges the Euclidean moduli problem with its transpose (left \(\mathrm{SL}_2(\mathbb{Z})\), right translations). One checks that \(\sigma X\) is again in bounded-disk normalization, with lattice \(\overline{\mathbb{Z}a + \mathbb{Z}c}\) — the *columns* replacing the rows. The two-sided theory pairing a disk with its \(\sigma\)-disk lives on the product of the two transposed problems; see the outlook (§6).

## 3. Schmidt disks of curvature \(2n\): the exact structure

Now let \(X \in \mathrm{SL}_2(\mathbb{Z}[i])\), \(D = X(\mathbb{H})\) a Schmidt disk of curvature \(2n\). The unit square \([0,1) + [0,1)i\) is a fundamental domain for the translations, so "disks in the unit square" = "disks mod \(\mathbb{Z}[i]\)", and [euclidean-counting.md](euclidean-counting.md) counts them: \(N_e(n) = n\prod_{p \mid n \text{ odd}}(1 - \chi_{-4}(p)/p)\).

**Lemma 1 (disk \(\leftrightarrow\) lattice).** The map \(D \mapsto \Lambda_D = \mathbb{Z}c + \mathbb{Z}d\) is a bijection
$$
\{\text{Schmidt disks of curvature } 2n\} / \mathbb{Z}[i]\text{-translation}
\;\xrightarrow{\ \sim\ }\;
\{\Lambda \subseteq \mathbb{Z}[i] \text{ of index } n \text{ with } \Lambda\,\mathbb{Z}[i] = \mathbb{Z}[i]\},
$$
onto the **primitive index-\(n\) sublattices** (primitive: the \(\mathbb{Z}[i]\)-ideal generated by \(\Lambda\) is everything).

*Proof.* \(D\) determines \(X\) up to \(X \mapsto T_\lambda X h\), \(h \in \pm\mathrm{SL}_2(\mathbb{Z})\) (the Möbius stabilizer of \(\mathbb{H}\) inside \(\mathrm{PSL}_2(\mathbb{Z}[i])\) is \(\mathrm{PSL}_2(\mathbb{Z})\)), so the oriented lattice \(\Lambda\) is well-defined; it has index \([\mathbb{Z}[i] : \Lambda] = |\operatorname{Im}(\bar c d)| = n\) and is primitive because \(ad - bc = 1 \in (c, d)\). Conversely an oriented basis \((c,d)\) of a primitive \(\Lambda\) has \(\gcd_{\mathbb{Z}[i]}(c,d) = 1\), so it completes to \(X \in \mathrm{SL}_2(\mathbb{Z}[i])\); the completions form exactly one translation class, since replacing the top row by \((a + tc, b + td)\) shifts \(\zeta = i(a\bar d - b\bar c)\) by \(2nt\), i.e. the center by \(t\). \(\square\)

**Lemma 2 (conductor \(=\) curvature).** Every primitive index-\(n\) sublattice \(\Lambda \subseteq \mathbb{Z}[i]\) is a proper \(\mathcal{O}_n\)-lattice: its multiplier ring is exactly \(\mathbb{Z} + n\mathbb{Z}[i]\).

*Proof.* Primitivity forces the quotient \(\mathbb{Z}[i]/\Lambda\) to be cyclic (\(\Lambda \subseteq p\mathbb{Z}[i]\) would violate \(\Lambda\mathbb{Z}[i] = \mathbb{Z}[i]\)), so \(\Lambda = \ker\phi\) for a surjection \(\phi: \mathbb{Z}[i] \to \mathbb{Z}/n\), \(\phi(x + yi) = sx + ty\), \(\gcd(s, t, n) = 1\). A direct computation (eliminate the scalar by which \(x + yi\) acts on the cyclic quotient) gives
$$
x + yi \in \mathcal{O}(\Lambda) \iff n \mid y\,(s^2 + t^2),
$$
so the conductor is \(n/\gcd(n, s^2+t^2)\). For \(p \mid n\): \(p \mid s^2 + t^2\) iff the line \(\ker(\phi \bmod p)\) is stable under multiplication by \(i\), i.e. is a *prime ideal line* of \(\mathbb{Z}[i]/p\) (at inert \(p\) there is none and \(p \mid s^2+t^2\) would force \(p \mid \gcd(s,t)\)). Primitivity excludes exactly the ideal lines, so \(\gcd(n, s^2+t^2) = 1\). \(\square\)

(The same count — \(p+1\) lines minus \(0/1/2\) ideal lines at inert/ramified/split \(p\) — re-derives \(N_e(n)\) of euclidean-counting.md from scratch.)

**Lemma 3 (every class, exactly twice).** For \(n \ge 2\), the map \(\Lambda \mapsto [\Lambda] \in \mathrm{Cl}(\mathcal{O}_n)\) from primitive index-\(n\) sublattices is surjective with all fibers \(\{\Lambda, i\Lambda\}\) of size two; the corresponding disk pairs are \(\{D, -D\}\). For \(n = 1\) it is the bijection \(\{\text{Ford disk}\} \to \mathrm{Cl}(\mathbb{Z}[i])\).

*Proof.* Given a proper \(\mathcal{O}_n\)-ideal class, pick an integral representative \(\mathfrak{a}\) prime to \(n\), and let \(\alpha\) generate \(\mathfrak{a}\mathcal{O}_K\) (\(\mathbb{Z}[i]\) is a PID). Then \(\Lambda = \alpha^{-1}\mathfrak{a} \subseteq \alpha^{-1}\mathfrak{a}\mathcal{O}_K = \mathcal{O}_K\) is primitive, proper for \(\mathcal{O}_n\), of index
$$
[\mathcal{O}_K : \alpha^{-1}\mathfrak{a}] = \frac{[\mathcal{O}_K : \mathfrak{a}]}{N(\mathfrak{a}\mathcal{O}_K)} = \frac{n\,N_{\mathcal{O}_n}(\mathfrak{a})}{N(\mathfrak{a}\mathcal{O}_K)} = n,
$$
the last step because extension of ideals prime to the conductor preserves norms; the construction is independent of all choices up to the unit \(\mu_4 = \{\pm 1, \pm i\}\), and \(-\Lambda = \Lambda\), \(i\Lambda \neq \Lambda\) (else \(i \in \mathcal{O}(\Lambda) = \mathcal{O}_n\)). Conversely any primitive \(\Lambda\) arises this way from its own class. Finally \(i\Lambda\) is the lattice of the disk \(-D\): \(\operatorname{diag}(i,-i)\,X\) has bottom row \((-ic, -id)\) and Möbius action \(z \mapsto -z\). Homothety by \(i \in K^\times\) does not change the class, so the fiber over each class is exactly \(\{\Lambda, i\Lambda\} \leftrightarrow \{D, -D\}\). \(\square\)

> **Theorem 1 (structure).** For \(n \ge 2\), the assignment \(D \mapsto [\Lambda_D]\) induces a bijection
> $$
> \{\text{Schmidt disks of curvature } 2n \text{ in } [0,1)+[0,1)i\}\,/\,(z \mapsto -z)
> \;\xrightarrow{\ \sim\ }\; \mathrm{Cl}(\mathbb{Z} + n\mathbb{Z}[i]),
> $$
> and in particular
> $$
> \boxed{\;N_e(n) \;=\; 2\,h(-4n^2) \quad (n \ge 2), \qquad N_e(1) = h(-4) = 1.\;}
> $$
> Under it, complex conjugation of disks (\(D \mapsto \bar D\)) is inversion in the class group, and the tangent-to-\(\hat{\mathbb{R}}\) disks (\(\zeta = \pm i\), the "Ford row") are the principal class, with \(w \sim ni\).

The count identity reconciles [euclidean-counting.md](euclidean-counting.md) with the class number formula for orders: \(h(\mathcal{O}_n) = \tfrac{n}{2}\prod_{p \mid n}\bigl(1 - \chi_{-4}(p)/p\bigr)\) for \(n > 1\). *The Euclidean census of the Gaussian Schmidt arrangement is literally a class number census of the non-maximal orders of \(\mathbb{Q}(i)\)* — the counting document's multiplicative function, seen arithmetically.

**Machine verification** (`python3 scripts/euclidean_moduli_invariants.py`, experiment A): for all \(n \le 24\), exact set equality of \(\zeta\)-classes with the classification congruence classes; every class of disc \(-4n^2\) hit exactly twice via the exact class map \(\Lambda = \mathbb{Z}(b+ai) + \mathbb{Z}d \mapsto\) reduced form of \((d^2, -2bd, a^2+b^2)\); conductor exactly \(n\) for every primitive lattice; \(i\)-pairing as stated.

## 4. The values of \(\beta = j(-X^{-1}(\infty))\): ring class polynomials, squared

The point \(w = d/c\) is a Gaussian rational in \(\mathbb{H}\), hence a CM point of \(K = \mathbb{Q}(i)\); by Lemma 2 its order is exactly \(\mathcal{O}_n\), and by Theorem 1 the disks of a level sweep out each class twice. Since \(j(w) = j([\Lambda_D])\):

> **Theorem 2 (singular moduli of the Euclidean family).** For \(n \ge 2\),
> $$
> \prod_{\substack{D \text{ Schmidt disk of curvature } 2n \\ \text{center} \in [0,1)+[0,1)i}} \bigl(x - j(-X_D^{-1}(\infty))\bigr) \;=\; H_{-4n^2}(x)^{2} \;\in\; \mathbb{Z}[x],
> $$
> the square of the ring class polynomial of \(\mathbb{Z} + n\mathbb{Z}[i]\) (degree \(2h(-4n^2) = N_e(n)\)); for \(n = 1\) it is \(H_{-4}(x) = x - 1728\). In particular the values are **algebraic integers forming full Galois orbits, and the answer to "are they the roots of a rational polynomial?" is yes — an integer polynomial, each root twice.**

The certified table for small \(n\) (experiment C; every coefficient integer-certified with hundreds of spare digits, and the disk-side product checked against \(H^2\) exactly):

| \(n\) | disc | \(h\) | \(H_{-4n^2}(x)\) |
|---|---|---|---|
| 1 | \(-4\) | 1 | \(x - 1728\) |
| 2 | \(-16\) | 1 | \(x - 287496\) |
| 3 | \(-36\) | 2 | \(x^2 - 153542016\,x - 1790957481984\) |
| 4 | \(-64\) | 2 | \(x^2 - 82226316240\,x - 7367066619912\) |
| 5 | \(-100\) | 2 | \(x^2 - 44031499226496\,x - 292143758886942437376\) |
| 6 | \(-144\) | 4 | \(x^4 - 23578503968570400x^3 + 269499185406087942528x^2 + 490453856866850787293184x + 571751321233328637579104256\) |
| 7 | \(-196\) | 4 | \(x^4 - 12626092121367165696x^3 - 44864481851299856707307347968x^2 + 250850701957837760512539510177792x - 2108010653658430719613224868701536256\) |

The anchors are classical: the Ford disks see \(j(i) = 1728\); the curvature-4 disks see \(j(2i) = 66^3\); the curvature-6 disks see \(j(3i)\) and \(j(\tfrac{1+3i}{2})\), the conjugate pair generating \(\mathbb{Q}(\sqrt3)\).

**The duality with the hyperbolic study.** At hyperbolic level \(n\), the \(\beta\)-values were singular moduli of discriminant \(1 - n^2\): *fixed (maximal) order type, field varying with the level*. At Euclidean curvature \(2n\) they are singular moduli of discriminant \(-4n^2\): **fixed field \(\mathbb{Q}(i)\), conductor equal to the level**. The Schmidt arrangement thus sees, along its two natural gradings, the two orthogonal directions of the theory of complex multiplication — discriminant aspect and conductor aspect.

**Corollary (traces: the even-square slice).** Let \(\mathrm{Tr}^E_n := \sum_{\mathfrak{c} \in \mathrm{Cl}(\mathcal{O}_n)}(j(\mathfrak{c}) - 744)\) for \(n \ge 2\), and \(\mathrm{Tr}^E_1 := \tfrac{1728 - 744}{2} = 492\). Then
$$
\sum_{g \mid n} \mathrm{Tr}^E_{n/g} \;=\; t(4n^2),
$$
Zagier's trace of singular moduli. Verified for all \(n \le 13\) (experiment C), e.g. \(\mathrm{Tr}^E_2 = 286752\) and \(t(16) = 286752 + 492 = 287244\). So the Euclidean census reads the coefficients of Zagier's weight-\(3/2\) form along the **even squares** \(d = 4n^2\) — the complementary diagonal to the hyperbolic slice \(d = n^2 - 1\) of moduli-invariants.md §3. (The two slicings meet the two natural quadratic progressions \(d = \square\) and \(d = \square - 1\); their interaction with the Hurwitz–Kronecker relations is taken up in the outlook.)

| \(n\) | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|
| \(\mathrm{Tr}^E_n\) | 286752 | 153540528 | 82226314752 | 44031499225008 | 23578503968567424 | 12626092121367162720 |

## 5. The phase on the Bianchi group: lemniscatic normalization and one integer polynomial per level

### 5.1 The period, and algebraicity

On integral \(X\) the five invariants produced ring class integers. The sixth is transcendental — but by exactly one period. Define
$$
\Omega \;:=\; \frac{\Gamma(\tfrac14)^4}{8\pi^2} \;=\; \frac{\varpi^2}{\pi} \;=\; \Bigl((2\pi)^6\,\eta(i)^{24}\Bigr)^{1/6} \;=\; 2.18843961522647\ldots,
\qquad
u(D) \;:=\; \frac{\Theta(X_D)}{\Omega},
$$
where \(\varpi = 2\int_0^1 dt/\sqrt{1-t^4} = \Gamma(\tfrac14)^2/\sqrt{8\pi}\) is the **lemniscate constant** — the period of the CM curve \(\mathbb{C}/\mathbb{Z}[i]\), whose homothety class is the unique class at \(n = 1\). (All three expressions for \(\Omega\) are verified to 100+ digits; experiment B.) From \(j'^{\,6} = -(2\pi)^6 j^4 (j - 1728)^3\,\Delta\) (with \(\Delta = \eta^{24}\)):
$$
u(D)^6 \;=\; -\,\beta^4\,(\beta - 1728)^3\; \frac{\Delta(\Lambda_D)}{\Delta(\mathbb{Z}[i])},
\qquad \beta = j(w_D),
$$
where \(\Delta(\Lambda)\) is the weight-12 lattice discriminant. The \(\Delta\)-quotient of commensurable CM lattices is algebraic and \(n^{12}\Delta(\Lambda)/\Delta(\mathbb{Z}[i])\) is an algebraic integer (classical; Lang, *Elliptic Functions* ch. 12), so **\(u(D)\) is an algebraic number and \(n^2 u(D)\) an algebraic integer**. The identity is verified numerically at every level (experiment B).

**Anchor (\(n = 2\)): the phase is an integer times \(i\).** For the curvature-4 disks, \(\beta = 66^3\), \(\beta - 1728 = 2^3 3^6 7^2\), and \(\Delta(\Lambda)/\Delta(\mathbb{Z}[i]) = \eta(2i)^{24}/\eta(i)^{24} = 2^{-9}\) (Weber), giving \(u^2 = -\beta(\beta - 1728)\,\gamma_2(2i)/8 = -2^4 3^{10} 7^2 11^4\) exactly:
$$
u \;=\; \pm\,2^2 3^5\, 7 \cdot 11^2\, i \;=\; \pm\,823284\,i .
$$
This is the Euclidean analogue of the hyperbolic anchor "\(u = -1\) exactly at \(n = 3\)" — with the characteristic difference that the Euclidean phase is not a unit but an \(S\)-integer (see §5.5).

### 5.2 Functional equations and the reality law

**Proposition (laws).** For every Schmidt disk \(D\) of curvature \(2n \ge 4\):
1. *(units)* \(u(-D) = -u(D)\); the two disks of a ring class carry opposite phases, so \(u^2\) is well-defined on \(\mathrm{Cl}(\mathcal{O}_n)\): write \(v_\mathfrak{c}^2\).
2. *(mirror)* \(u(\bar D) = \overline{u(D)}\), and \([\Lambda_{\bar D}] = [\Lambda_D]^{-1}\): complex conjugation of disks is inversion in the class group. (In lattice terms \(u(\bar\Lambda) = -\overline{u(\Lambda)}\), the sign because \(\bar D \leftrightarrow i\bar\Lambda\).)
3. *(reality quantization)* On an ambiguous class (\(\mathfrak{c}^2 = 1\)), \(u \in \mathbb{R} \cup i\mathbb{R}\); precisely, with \(\zeta = x + yi\) the curvature-center of either disk of the class:
$$
u \in \mathbb{R} \iff n \mid y \ (\text{center at half-integer height}), \qquad
u \in i\mathbb{R} \iff n \mid x \ (\text{center on the lines } \operatorname{Re} = 0, \tfrac12),
$$
and for \(n \ge 2\) exactly one of the two holds on each ambiguous class. Since \(y\) is odd, **even levels have no real classes**; at every odd \(n \ge 3\) the class of \(\zeta = (n+1) + ni\) is real — at all computed odd levels it is the norm-2 class, over the ramified prime.

*Proof.* (1) is the Borel sign of §2. (2): \(\bar D = Y(\mathbb{H})\) for \(Y = \bar X \operatorname{diag}(i,-i)\), whose invariant point is \(w_Y = -\bar w\); with the real-kernel identity \(j'(-\bar\tau) = -\overline{j'(\tau)}\) the factors combine to \(\Theta(Y) = \overline{\Theta(X)}\); and \(\bar\Lambda\) represents the conjugate, i.e. inverse, class. (3): ambiguity means \(\bar D \equiv D\) or \(\bar D \equiv -D\) mod translation, i.e. \(\bar\zeta \equiv \pm\zeta \pmod{2n}\), i.e. \(n \mid y\) resp. \(n \mid x\); then (2) and (1) give \(\bar u = u\) resp. \(\bar u = -u\). Both congruences together would force \(n^2 \mid x^2 + y^2 \equiv 1 \pmod{4n}\), impossible for \(n \ge 2\). Solvability of \(x^2 \equiv 1 - n^2 \pmod{4n}\) at \(y = n\): \(x = n + 1\) works for every odd \(n\). \(\square\)

The observed reality patterns over \(\mathrm{Cl}(\mathcal{O}_n)\) (R real, I imaginary, C complex; experiment D — the criterion of (3) is asserted class-by-class):

| \(n\) | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pattern | I | IR | II | IR | IICC | IRCC | IICC | IRCCCC | IICC | IRCCCC | IICCCICI | IRCCCC |

(classes in reduced-form order; principal class always I; the norm-2 class R at odd levels — so at \(n = 12\) all four ambiguous classes are I and the two conjugate pairs are C. This table is the Euclidean analogue of the hyperbolic sign table of outlook.md §1.1 — with the dichotomy now *proved* by the center criterion.)

### 5.3 One level, one monic integer polynomial (unconditional at sixth powers)

> **Theorem 3.** Fix \(n \ge 2\) and let \(v_\mathfrak{c}\) (\(\mathfrak{c} \in \mathrm{Cl}(\mathcal{O}_n)\)) be the phases, each taken at either disk of its class. Then
> $$
> P^{(6)}_n(y) \;:=\; \prod_{\mathfrak{c} \in \mathrm{Cl}(\mathcal{O}_n)} \bigl(y - (n^2 v_\mathfrak{c})^6\bigr) \;\in\; \mathbb{Z}[y],
> $$
> **monic** of degree \(h(-4n^2)\), and the product over all \(2h\) disks of the level is \(P^{(6)}_n(y)^2\). In particular the sixth powers of the normalized phases of a level are the roots of one monic integer polynomial.

*Proof.* (i) *Rationality.* For \(\tau \in \mathbb{H}\) and a cyclic index-\(n\) sublattice \(\Lambda_{a,b,d} = \mathbb{Z}(a\tau + b) + \mathbb{Z}d \subset \mathbb{Z}\tau + \mathbb{Z}\) (\(ad = n\), \(\gcd(a,b,d) = 1\)), the quantity \(F_{a,b,d}(\tau) := (n^2)^6\hat h_2(\Lambda_{a,b,d})^6 / \bigl((2\pi)^6\Delta(\tau)\bigr) = -n^{12} j_w^4 (j_w - 1728)^3 \Delta(w)/(d^{12}\Delta(\tau))\), \(w = (a\tau+b)/d\), is one member of a Hecke orbit: \(\mathrm{SL}_2(\mathbb{Z})\) permutes the \(F_{a,b,d}\), so every elementary symmetric function of them (and jointly of them with the \(j_w\)) is a meromorphic modular function, holomorphic off the cusp, with rational \(q\)-expansion (the \(\zeta_d\)-coefficients cancel in the symmetrization over \(b\), as for the classical modular polynomials) — hence a polynomial in \(j(\tau)\) with rational coefficients. Evaluating at \(\tau = i\) (\(j = 1728 \in \mathbb{Q}\)) makes all symmetric functions of the cyclic family rational, jointly with the \(\beta\)'s. (ii) *The primitive cut.* By Lemma 2 the primitive sublattices are exactly the cyclic ones whose \(\beta\) has discriminant \(-4n^2\); the discriminant of \(\beta\)'s order is a Galois invariant, so \(\mathrm{Gal}(\bar{\mathbb{Q}}/\mathbb{Q})\), which permutes the joint pairs \((\beta_\Lambda, F_\Lambda)\), preserves the primitive sub-multiset; its symmetric functions are therefore rational. Each class appears twice with the same sixth power (\((\pm v)^6\)), giving \(P^{(6)}_n \in \mathbb{Q}[y]\) and the square statement. (iii) *Integrality.* \(\beta\) and \(\beta - 1728\) are algebraic integers and \(n^{12}\Delta(\Lambda)/\Delta(\mathbb{Z}[i])\) is an algebraic integer; a monic rational polynomial with algebraic-integer roots lies in \(\mathbb{Z}[y]\). \(\square\)

Contrast with the hyperbolic Theorem C (moduli-invariants.md §5.9): same shape — *one level, one integer polynomial* — but the Euclidean polynomial is **monic**: the hyperbolic denominators (Gross–Zagier primes of \((1-n^2, -3), (1-n^2, -4)\)) came from kernel values at *both* CM points in a ratio; here there is one kernel value and its elliptic-point factors \(\beta^4(\beta-1728)^3\) sit in the *numerator*. All denominators of the Euclidean phase live over \(n\) (§5.5), and \(n^2\) clears them.

### 5.4 The certified sharpening: squares, minimal normalizer, irreducibility, Galois groups

Experimentally the descent from sixth powers to **squares** holds at every computed level, with a smaller normalizer than \(n^2\). Let \(\lambda_n\) be the least positive integer with \(\prod_\mathfrak{c}\bigl(y - (\lambda_n v_\mathfrak{c})^2\bigr) \in \mathbb{Z}[y]\); then define
$$
P^{(2)}_n(y) \;:=\; \prod_{\mathfrak{c} \in \mathrm{Cl}(\mathcal{O}_n)} \bigl(y - (\lambda_n v_\mathfrak{c})^2\bigr).
$$
Certified at 140–600 digits with the absolute-error criterion (experiment D; every coefficient carries 100+ spare digits; \(P^{(6)}\) is then derived *exactly* from \(P^{(2)}\) by the cube resultant and re-checked numerically; the levels \(14 \le n \le 16\), at up to 840 digits, are covered by the irreducibility suite of §5.6):

$$
\begin{aligned}
P^{(2)}_2(y) &= y + 677796544656 = y + 2^4 3^{10} 7^2 11^4 \qquad (\lambda_2 = 1)\\
P^{(2)}_3(y) &= y^2 + 194359980365709312\,y - 61348044280308881216765952 \qquad (\lambda_3 = 1)\\
P^{(2)}_4(y) &= y^2 + 891727674920361943695108\,y + 669500278835555301992309603844 \qquad (\lambda_4 = 4)\\
P^{(2)}_5(y) &= y^2 + 9988440693113556636537936936960\,y \\&\qquad - 566463402951991282866091504109681092004916756480 \qquad (\lambda_5 = 25)\\
P^{(2)}_6(y) &= y^4 + 4582717575900659305402000363330368\,y^3 \\&\qquad + 351294264522418974183125673664523963340288\,y^2 \\&\qquad - 806641252061379699326345644173807440652498026496\,y \\&\qquad + 458309874496608718457444408199953336549693979660189696 \qquad (\lambda_6 = 1)\\
P^{(2)}_7(y) &= y^4 + 3155154778580067195294294125640706443706368\,y^3\\&\qquad - 197111694939076777568905182387467955698550112233175852728818073600\,y^2\\&\qquad + 4575429520322458394762189283839982543295080609171830501600451900661064269824\,y\\&\qquad - 239522477510167975017898749934852425525535101345499925886422989252893280414570855268352 \quad (\lambda_7 = 49)
\end{aligned}
$$

| \(n\) | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| \(h(-4n^2)\) | 1 | 2 | 2 | 2 | 4 | 4 | 4 | 6 | 4 | 6 | 8 | 6 |
| \(\lambda_n\) | 1 | 1 | 4 | 25 | 1 | 49 | \(2^5\) | \(3^3\) | 25 | \(11^2\) | 4 | \(13^2\) |
| \(\mathrm{Cl}(\mathcal{O}_n)\) | 1 | \(\mathbb{Z}/2\) | \(\mathbb{Z}/2\) | \(\mathbb{Z}/2\) | \(\mathbb{Z}/4\) | \(\mathbb{Z}/4\) | \(\mathbb{Z}/4\) | \(\mathbb{Z}/6\) | \(\mathbb{Z}/4\) | \(\mathbb{Z}/6\) | \(\mathbb{Z}/4{\times}\mathbb{Z}/2\) | \(\mathbb{Z}/6\) |
| \(\mathrm{Gal}(P^{(2)}_n)\) | 1 | \(\mathbb{Z}/2\) | \(\mathbb{Z}/2\) | \(\mathbb{Z}/2\) | \(D_4\) | \(D_4\) | \(D_4\) | \(D_6\) | \(D_4\) | \(D_6\) | order 16, even | \(D_6\) |

Three structural facts, holding at every computed level:

- **Irreducibility.** \(P^{(2)}_n\) is irreducible over \(\mathbb{Q}\) of degree exactly \(h\): one level, one Galois orbit — so \(P^{(2)}_n\) is the *minimal* polynomial of every \((\lambda_n v_\mathfrak{c})^2\), and \([\mathbb{Q}(v_\mathfrak{c}^2) : \mathbb{Q}] = h(-4n^2)\). (The analogue of moduli-invariants.md §5.10, now visible already at square level.)
- **The Galois group is the generalized dihedral group of \(\mathrm{Cl}(\mathcal{O}_n)\)** — \(\mathrm{Cl} \rtimes \{\pm1\}\) in its action on \(\mathrm{Cl}\): \(D_4\) at the \(\mathbb{Z}/4\)-levels, \(D_6\) at the \(\mathbb{Z}/6\)-levels (sympy `galois_group`), and at \(n = 12\) (\(\mathrm{Cl} = \mathbb{Z}/4\times\mathbb{Z}/2\), degree 8) the discriminant of \(P^{(2)}_{12}\) is a perfect square, exactly as predicted: translations and inversion of \(\mathbb{Z}/4\times\mathbb{Z}/2\) are all even permutations, while every group with a cyclic factor \(\mathbb{Z}/4\) or \(\mathbb{Z}/6\) contains an odd translation cycle — and indeed those discriminants are non-squares. This is the ring-class-field Galois structure (translation by \(\mathrm{Cl}\), inversion by the mirror law) acting on the phases, i.e. the Euclidean counterpart of hyperbolic Theorem A.
- **First powers.** The disk-level multiset \(\{\lambda_n u(D)\}\) is \(\{\pm\lambda_n v_\mathfrak{c}\}\), the root set of the even polynomial \(P^{(2)}_n(x^2) \in \mathbb{Z}[x]\) of degree \(2h = N_e(n)\). So the answer to "are the sixth-invariant values roots of a rational polynomial?" is: *after dividing by the single transcendental \(\Omega = \varpi^2/\pi\) and clearing \(\lambda_n \mid n^2\) — yes, of a monic integer polynomial.*

The descent \(u^6 \to u^2\) has a visible mechanism when \(3 \nmid n\): \(u^2 = -\beta\,(\beta-1728)\,\gamma_2(w)\,R_4(\Lambda)\) with \(\gamma_2 = j^{1/3}\) the Weber function (in the ring class field precisely when \(3 \nmid \operatorname{disc}\)) and \(R_4(\Lambda) = \eta(w)^8/(c^4\,\eta(i)^8)\) the weight-4 \(\eta^8\)-quotient (\(R_4^3 = \Delta\)-quotient). Remarkably the certified integrality holds **also at \(3 \mid n\)** (\(n = 3, 6, 9, 12\)), where \(\gamma_2(w)\) and \(R_4(\Lambda)\) individually leave the ring class field but their product does not — the Euclidean analogue of the hyperbolic first-power descent, and open in the same way (§6).

### 5.5 Prime structure: \(S_n\)-integrality, collision numerators, and the \(\Delta\)-mass

**Denominators.** The unnormalized symmetric functions \(e_k(\{v^2\})\) are rationals whose denominators are supported **only at primes dividing \(n\)** (verified at every level — e.g. at \(n = 7\) the denominators of \(e_1, \dots, e_4\) are \(7^3, 7^6, 7^9, 7^{13}\)); this is the \(S_n\)-integrality forced by the \(\Delta\)-quotient. The minimal normalizer \(\lambda_n\) (table above) shows a fine structure begging for a valuation law: trivial at \(n = 2, 3, 6\), equal to \(p^2\) at \(n = 5, 7, 10, 11, 13\), but \(2^2\) at \(n \in \{4, 12\}\), \(2^5\) at \(n = 8\), \(3^3\) at \(n = 9\) — the analogue of the exact-valuation problem of moduli-invariants.md §5.7, now purely local at \(n\).

**Numerators.** The support of \(P^{(2)}_n(0)\) always lies in \(\{p \mid 2n\} \cup \operatorname{supp} H_{-4n^2}(0) \cup \operatorname{supp} H_{-4n^2}(1728)\) (experiment D'), with exponents transferred from \(u^6 = -\beta^4(\beta - 1728)^3\cdot(\Delta\text{-quotient})\): away from \(2n\), a prime enters \(P^{(2)}_n(0)\) with \(\tfrac43\) times its \(H(0)\)-exponent plus its \(H(1728)\)-exponent (e.g. at \(n = 13\): \(503^3 \| H(0) \Rightarrow 503^4 \| P^{(2)}(0)\); \(47^4 \| H(1728) \Rightarrow 47^4 \| P^{(2)}(0)\)). The two value-norms are collision quantities in the Gross–Zagier sense: \(H(0) = \pm\prod_\mathfrak{c} j(\mathfrak{c})\) measures collisions of the conductor-\(n\) points with \(j = 0\) (discriminant pair \((-4n^2, -3)\), coprime — genuine GZ), while \(H(1728) = \pm\prod_\mathfrak{c}(j(\mathfrak{c}) - 1728)\) measures collisions with \(j(i)\) — a **same-field, conductor-degenerate** Gross–Zagier situation (\(-4 \mid -4n^2\)), i.e. arithmetic intersection with the \(n = 1\) stratum of the same family.

**The \(\Delta\)-mass.** The total \(\Delta\)-quotient of a level is startlingly clean. Define \(M(n) := \prod_{\mathfrak{c}} n^{12}\,\Delta(\Lambda_\mathfrak{c})/\Delta(\mathbb{Z}[i])\). Computed directly (independently of the polynomials) for \(2 \le n \le 16\) and \(n = 18, 21, 25, 27, 49\):

| \(n\) | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 25 | 27 | 49 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| \(M(n)\) | \(2^3\) | \(3^6\) | \(-2^9\) | \(1\) | \(2^{12}3^{12}\) | \(7^6\) | \(-2^{21}\) | \(3^{24}\) | \(2^{12}\) | \(11^6\) | \(2^{36}3^{24}\) | \(1\) | \(2^{24}7^{12}\) | \(3^{24}\) | \(-2^{45}\) | \(1\) | \(3^{78}\) | \(7^{48}\) |

> **Theorem 4 (the \(\Delta\)-mass law).** For every \(n \ge 2\),
> $$
> M(n) \;=\; \varepsilon(n) \prod_{\substack{p^k \| n \\ p \text{ not split}}} p^{\;\frac{6}{e_p}\cdot\frac{p^k - 1}{p - 1}\; N_e(n/p^k)},
> \qquad
> \varepsilon(n) = \begin{cases} -1 & n = 2^k,\ k \ge 2,\\ +1 & \text{otherwise,}\end{cases}
> $$
> with \(e_2 = 2\) (ramified), \(e_p = 1\) (inert); **split primes contribute nothing at all** — e.g. \(M(5) = M(13) = M(25) = 1\) exactly.

*Proof.* Four steps; each is machine-verified separately in [scripts/mass_law_and_irreducibility.py](scripts/mass_law_and_irreducibility.py).

**(M1) The full mass is a constant.** For \(\tau \in \mathbb{H}\) let \(A(n)(\tau) := \prod_{[L_\tau : \Lambda] = n} \Delta(\Lambda)/\Delta(L_\tau)\), the product over **all** \(\sigma(n)\) index-\(n\) sublattices of \(L_\tau = \mathbb{Z}\tau + \mathbb{Z}\). The family of sublattices is lattice-intrinsic and \(\Delta\)-ratios are homothety-invariant, so \(A(n)\) is \(\mathrm{SL}_2(\mathbb{Z})\)-invariant; in the Hermite parametrization \(\Lambda_{a,b,d} = \mathbb{Z}(a\tau+b) + \mathbb{Z}d\) (\(ad = n\), \(0 \le b < d\)) each factor \(d^{-12}\Delta_q(\tfrac{a\tau+b}{d})/\Delta_q(\tau)\) is holomorphic and nonvanishing on \(\mathbb{H}\) and meromorphic in \(q\) at the cusp. So \(A(n)\) is a rational function of \(j\) with no zeros or poles in \(\mathbb{C}\) — a constant, equal to its \(q \to 0\) limit. The factor \((a,b,d)\) contributes leading term \(d^{-12}\zeta_d^b\, q^{a/d - 1}\); the \(q\)-exponents cancel (\(\sum_{ad = n} d\,(\tfrac{a}{d} - 1) = \sum (a - d) = 0\)) and \(\zeta_d^{0 + 1 + \cdots + (d-1)} = e^{\pi i (d-1)} = (-1)^{d-1}\), so
$$
A(n) \;=\; (-1)^{t(n)} \prod_{d \mid n} d^{-12d}, \qquad t(n) = \#\{d \mid n : d \text{ even}\}.
$$

**(M2) Stratification over \(\mathbb{Z}[i]\), and the recursion.** Now put \(\tau = i\). For \(\Lambda \subseteq \mathbb{Z}[i]\) of index \(n\), the ideal \(\mathfrak{d} := \Lambda\,\mathbb{Z}[i]\) is principal, \(\mathfrak{d} = (\delta)\), and \(\Lambda = \mathfrak{d}\,\Lambda''\) with \(\Lambda'' = \mathfrak{d}^{-1}\Lambda\) **primitive** of index \(n/N\mathfrak{d}\); the assignment \(\Lambda \mapsto (\mathfrak{d}, \Lambda'')\) is a bijection, and \(\Delta(\delta\Lambda'') = \delta^{-12}\Delta(\Lambda'')\). Writing \(Q(k) := \prod_{\Lambda \text{ prim., index } k} \Delta(\Lambda)/\Delta(\mathbb{Z}[i])\), \(r(m) := \#\{\text{ideals of norm } m\}\), and \(\gamma(m) :=\) a generator of \(\prod_{N\mathfrak{d} = m}\mathfrak{d}\) — so \(\gamma(m)^{12}\) is well-defined, and rational because the ideal product is conjugation-stable (\(\bar\gamma = u\gamma\), \(u \in \mu_4\)) — the constant of M1 factors as
$$
A(n) \;=\; \prod_{m \mid n} \gamma(m)^{-12\,N_e(n/m)}\; Q(n/m)^{\,r(m)}.
$$
The \(m = 1\) stratum carries \(Q(n)\) to the first power, so this recursion determines \(Q(n) \in \mathbb{Q}\) inductively from \(Q(1) = 1\) — *the rationality of the mass needs no CM theory* — and \(Q(n) = \bigl[\prod_\mathfrak{c}\Delta(\Lambda_\mathfrak{c})/\Delta(\mathbb{Z}[i])\bigr]^2 \ge 0\) because \(\Delta(i\Lambda) = \Delta(\Lambda)\) pairs the two lattices of each class. Hence \(|M(n)| = n^{12h}\sqrt{Q(n)}\) is pinned by the recursion: it suffices to check that the claimed value satisfies it.

**(M3) The exponents.** Both sides are units away from \(n\). Fix \(p \mid n\), write \(n = p^k n'\) (\(p \nmid n'\)). Every ingredient factors along \(p\): \(v_p\bigl(\gamma(p^j m')^{12}\bigr) = r(m')\,g_j\) with
$$
g_j := v_p\bigl(\gamma(p^j)^{12}\bigr) = \begin{cases} 6j(j+1) & p \text{ split}\\ 6j\,[\,2 \mid j\,] & p \text{ inert}\\ 6j & p = 2, \end{cases}
$$
\(v_p(A(n)) = -12\bigl(\sum_{j \le k} j p^j\bigr)\sigma(n')\), and \(r, N_e\) are multiplicative with \(r * N_e = \sigma\) (Dirichlet convolution — the counting shadow of the stratification, since \(r = 1 * \chi_{-4}\) and \(N_e = \mathrm{Id} * \mu\chi_{-4}\)). Substituting the ansatz \(v_p\bigl(Q(p^{j} m)\bigr) = v_p\bigl(Q(p^{j})\bigr) N_e(m)\) (\(p \nmid m\)) — which the claimed closed form satisfies, since \(h(p^j m) = \tfrac12 N_e(p^j)N_e(m)\) — the recursion factors as (\(p\)-power recursion) \(\times\) (\(r * N_e = \sigma\) over \(n'\)), and everything reduces to \(n = p^k\). There, in generating functions \(\sum_j (\cdot)\,T^j\), the recursion is the single identity
$$
\alpha + G\,V \;=\; R\,S, \qquad
V = \frac{1 - \chi T}{1 - pT},\quad R = \frac{1}{(1-T)(1-\chi T)},\quad \alpha = \frac{-12\,pT}{(1-T)(1-pT)^2},
$$
\(\chi = \chi_{-4}(p)\), \(G = \sum g_jT^j\), and \(S = \sum v_p(Q(p^j))T^j\) the claimed series. All three cases are one-line rational-function checks; e.g. split (\(\chi = 1\), \(G = \tfrac{12T}{(1-T)^3}\), \(S = \tfrac{-12(p-1)T}{(1-pT)^2}\)):
$$
\alpha + GV = \frac{12T}{(1-T)(1-pT)}\Bigl[\frac{-p}{1-pT} + \frac{1}{1-T}\Bigr] = \frac{-12(p-1)T}{(1-T)^2(1-pT)^2} = RS;
$$
inert (\(\chi = -1\), \(G = \tfrac{12T^2}{(1-T^2)^2}\), \(S = \tfrac{12T}{(1-T)(1-pT)} - \tfrac{12(1+p)T}{(1-pT)^2}\)) and ramified (\(p = 2\), \(\chi = 0\), \(G = \tfrac{6T}{(1-T)^2}\), \(S = \tfrac{6T}{(1-T)(1-2T)} - \tfrac{24T}{(1-2T)^2}\)) reduce to the numerator identities \(-p(1-T^2) + T(1-pT) = T - p\) and \(-4(1-T) + (1-2T) = 2T - 3\). This proves \(|M(n)|\) equals the claimed product.

**(M4) The sign.** In \(M(n) = \prod_\mathfrak{c} R_\mathfrak{c}\), conjugation pairs the factor of \(\mathfrak{c}\) with that of \(\mathfrak{c}^{-1}\) (\(\Delta(\bar\Lambda) = \overline{\Delta(\Lambda)}\)), so \(\operatorname{sign} M(n) = \prod_{\text{ambiguous}} \operatorname{sign} R_\mathfrak{c}\), each such factor being real. On an **I-class** every lattice is conjugation-stable (\(\bar\Lambda = \Lambda\), from \(\bar D \equiv -D\)), hence of Hermite shape \(\langle d, b + ai\rangle\) with \(b \in \{0, d/2\}\); then \(\Delta(\Lambda) = d^{-12}(2\pi)^{12}\Delta_q(\tfrac{b}{d} + \tfrac{a}{d}i)\) is real with the sign of \(q_w = e^{2\pi i b/d}e^{-2\pi a/d}\): **positive for \(b = 0\), negative for \(b = d/2\)**. On an **R-class** (\(\bar\Lambda = i\Lambda\); odd \(n\)), the lattice \((1-i)\Lambda\) is conjugation-stable, and it cannot be rectangular (\(b = 0\)): \((1-i)\Lambda = \mathbb{Z}r + \mathbb{Z}si\) would force \(r, s\) even and \(\Lambda \subseteq (1+i)\mathbb{Z}[i]\), contradicting primitivity; so \(\Delta((1-i)\Lambda) < 0\) and \(\Delta(\Lambda) = (1-i)^{12}\Delta((1-i)\Lambda) = -64\,\Delta((1-i)\Lambda) > 0\). Since the sign is constant on a class, the two I-flavors separate classes; completing \(\langle d, ai\rangle\) gives \(\zeta = (1 - 2au)i\) (\(ua + vd = 1\)), so \(b = 0\) lattices lie in the \(x = 0\) classes, and a count (\(\#\{b = 0 \text{ primitive lattices}\} = 2^{\omega(n)}\) = twice the number of \(x = 0\) classes) forces \(b = d/2 \iff x \equiv n\). Finally the \(x = n\) classes are counted by the odd solutions of \(y^2 \equiv 1 - n^2 \pmod{4n}\) mod \(2n\), up to sign: zero unless \(4 \mid n\) (a mod-8 obstruction), and for \(n = 2^k m\) (\(k \ge 2\), \(m\) odd) equal to \(2^{\omega(m)}\) by CRT (four square roots of \(1\) mod \(2^{k+2}\), \(2^{\omega(m)}\) mod \(m\)). So \(\operatorname{sign} M(n) = (-1)^{2^{\omega(m)}}\) for \(4 \mid n\) and \(+1\) otherwise: \(-1\) exactly when \(m = 1\), i.e. \(n = 2^k\), \(k \ge 2\). \(\blacksquare\)

**Machine verification** (`python3 scripts/mass_law_and_irreducibility.py`): M1 at three random \(\tau\) per level to 30+ digits; the recursion of M2 solved in exact rational arithmetic reproduces the closed form for **all \(n \le 60\)** (with the \(\gamma^{12}\)-sign bookkeeping matching \((-1)^{t(n)}\)); the three generating-function identities of M3 verified symbolically; the sign analysis of M4 — reality and sign of every ambiguous class, Hermite flavors, the \(x = n\) count \(2^{\omega(m)}\) — checked class by class for all \(n \le 20\). Together with the direct evaluation table above (through \(n = 49\)), the law is proved and doubly verified.

**Remark (Deuring reading).** The proof is elementary, but the shape of the exponents keeps its conceptual reading: split \(p\) = ordinary reduction of the lemniscatic curve (the class-group product cancels along the horizontal isogeny volcano), non-split \(p\) = supersingular reduction, where the geometric series \(\tfrac{p^k-1}{p-1}\) is precisely the growth of quasi-canonical-lifting valuations (Gross). Matching the proof's global bookkeeping to the per-class local valuations remains open — see §6.1.

### 5.6 Irreducibility: one level, one Galois orbit — at first power

Every polynomial this study introduces is now settled (exact integer factorization over \(\mathbb{Q}\), all levels \(2 \le n \le 16\); `python3 scripts/mass_law_and_irreducibility.py irred`):

| polynomial | degree | irreducible? |
|---|---|---|
| \(H_{-4n^2}(x)\) | \(h\) | **yes, for every \(n \ge 1\)** — a theorem of ring class field theory: \(j(\mathcal{O}_n)\) generates the ring class field, of degree \(h(-4n^2)\) over \(K\), and being real it has the same degree over \(\mathbb{Q}\), so \(H_{-4n^2}\) is its minimal polynomial over \(\mathbb{Q}\) (and over \(K\)) — Cox, *Primes of the form \(x^2 + ny^2\)*, §9; checked exactly for \(n \le 16\) |
| \(H_{-4n^2}(x)^2\) (the disk \(j\)-product of Theorem 2) | \(2h\) | **no, by construction** — its irreducible factorization is the two copies of \(H_{-4n^2}\) |
| \(P^{(2)}_n(y)\) | \(h\) | **yes at every computed level** \(2 \le n \le 16\), with exact squarefreeness certificates \(\gcd\bigl(P^{(2)}, P^{(2)\prime}\bigr) = 1\) |
| \(P^{(6)}_n(y)\) | \(h\) | **yes at every computed level** \(2 \le n \le 16\) (and squarefree): cubing the square-phases loses no degree, \(\mathbb{Q}\bigl((n^2v_\mathfrak{c})^6\bigr) = \mathbb{Q}\bigl((\lambda_n v_\mathfrak{c})^2\bigr)\) of degree \(h\) |
| \(P^{(2)}_n(x^2)\) (the first-power disk polynomial) | \(2h\) | **yes at every computed level** \(2 \le n \le 16\) |

The last line is the strongest statement: the \(2h = N_e(n)\) normalized phases \(\pm\lambda_n v_\mathfrak{c}\) of a level — one per disk in the unit square — form a **single Galois orbit**, \(P^{(2)}_n(x^2)\) is their common minimal polynomial, and
$$
[\mathbb{Q}(\lambda_n v_\mathfrak{c}) : \mathbb{Q}] \;=\; 2h(-4n^2) \;=\; N_e(n)
\qquad\text{for every disk of every computed level:}
$$
passing from \(v^2\) to \(v\) always doubles the degree (no phase is rational over its square's field), and the \(\pm\)-pair ambiguity of §5.2 is exactly the quadratic step. This sharpens the hyperbolic picture (moduli-invariants.md §5.10), where the first-power orbit statement needed the \(\mathfrak{r}\)-pair bookkeeping; here the disk polynomial itself is irreducible.

The computation extends the \(\lambda_n\)-table: \(\lambda_{14} = 7^2\), \(\lambda_{15} = 5^2\), \(\lambda_{16} = 2^7\) — and at the split level \(n = 15\) the unnormalized denominators of \(e_1, \dots, e_8\) form the exactly arithmetic ladder \(5^4, 5^8, \dots, 5^{32}\) (compare the irregular inert and 2-adic ladders of §5.5), a clean hint for the valuation problem of §6.1.

**Conditional theorem (the analogue of moduli-invariants.md Theorem D).** Grant the square-level translation law of §6.7 (\(\sigma\) acts on \(\{(\lambda_nv_\mathfrak{c})^2\}\) through \(\mathfrak{c} \mapsto \mathfrak{c}^{e(\sigma)}\mathfrak{c}(\sigma)\)). Since translations act transitively on \(\mathrm{Cl}(\mathcal{O}_n)\), all roots of \(P^{(2)}_n\) are then conjugate, the coincidence set \(\{\mathfrak{t} : (\lambda v_\mathfrak{t})^2 = (\lambda v_1)^2\}\) is a subgroup \(T\), and \(P^{(2)}_n = m^{|T|}\) for an irreducible \(m\): **\(P^{(2)}_n\) is irreducible iff it is squarefree**, a finite exact computation — certified above at every level. So irreducibility for all \(n\) rests only on the translation law plus per-level squarefreeness, exactly as in the hyperbolic §5.10.

### 5.7 What changed relative to the hyperbolic phase

| | hyperbolic (moduli-invariants.md) | Euclidean (this document) |
|---|---|---|
| left/right groups | \(\mathrm{SL}_2(\mathbb{Z}) \times \mathrm{SL}_2(\mathbb{Z})\) | translations \(\mathbb{Z}[i]\) \(\times\) \(\mathrm{SL}_2(\mathbb{Z})\) |
| CM data of a level | disc \(1 - n^2\), varying field | fixed field \(\mathbb{Q}(i)\), conductor \(n\) |
| count of a level | \(3H(n^2-1)\) (Hurwitz-weighted) | \(N_e(n) = 2h(-4n^2)\) (exact) |
| \(j\)-polynomial | Hilbert class polynomial \(H_{1-n^2}\) | ring class polynomial squared \(H_{-4n^2}^2\) |
| trace slice | \(t(n^2 - 1)\) | \(t(4n^2)\) |
| sixth invariant | two kernels + derivative; fiber rate \(\sqrt{n^2-1}\) | one kernel + residue; fiber rate \(2\) |
| normalization | \(\varepsilon = n + \sqrt{n^2-1}\) (real quadratic unit) | \(\Omega = \varpi^2/\pi\) (lemniscatic period) |
| level polynomial | integer, non-monic (GZ denominators) | **monic** integer, denominators over \(n\) only |
| twist law | \(u_f u_{\mathfrak{r}f} = 1\) (norm-1 pairs) | none needed: absolute values via \(\Delta\)-mass |
| sign/reality law | experimental sign table (outlook 1.1) | **proved** center-criterion (\(n \mid x\) vs \(n \mid y\)) |

## 6. Research outlook

### Small (days to weeks)

**6.1 Exact valuations of \(\lambda_n\): the per-class refinement of Theorem 4.** Theorem 4 pins the *total* \(p\)-valuation of a level's \(\Delta\)-quotients; the open local problem is its distribution over the classes — the per-coefficient denominators of \(e_k(\{v^2\})\) (data: \(7^3, 7^6, 7^9, 7^{13}\) at \(n = 7\); \(3^2, 3^6, 3^{16}, 3^{19}, 3^{22}, 3^{33}\) at \(n = 9\); the exactly arithmetic \(5^{4k}\) at the split level \(n = 15\)), hence the minimal normalizers \(\lambda_n\). This is a Newton-polygon-over-classes statement, and the natural tool is exactly the quasi-canonical valuation calculus of Gross that the Remark after Theorem 4 points at: compute \(v_\mathfrak{p}(\Delta(\Lambda_\mathfrak{c})/\Delta(\mathbb{Z}[i]))\) class by class for \(n = 7, 9, 15\) and match against the filtration level of \(\Lambda_\mathfrak{c} \otimes \mathbb{Z}_p\); Theorem 4 then serves as the mass check on any candidate law.

**6.2 The square-descent at \(3 \mid n\).** \(P^{(2)}\)-integrality is certified at \(n = 3, 6, 9, 12\) although \(\gamma_2(w)\) leaves the ring class field there. Identify the invariant combination \(\gamma_2(w) R_4(\Lambda)\) in Weber/Schertz terms (the \(\mathfrak{f}\)-function calculus for conductor \(3 \mid n\)) and prove the descent — the exact analogue of the hyperbolic first-power problem, one power lower. Also settle whether first powers (not just squares) satisfy an integer polynomial after resolving the \(\pm\)-pairing by orientation data.

**6.3 The real-class census.** By §5.2 the real classes are those meeting the half-height row \(y = n\); their number is the number of solutions of \(x^2 \equiv 1 - n^2 \pmod{4n}\) up to sign — a genus-theory count. Data shows exactly one real class (the norm-2 class) at every odd \(n \le 13\); decide whether that persists (it should fail once \(\mathrm{Cl}(\mathcal{O}_n)[2]\) is large enough) and derive the full R/I distribution on \(\mathrm{Cl}[2]\) from genus characters — the proved counterpart of the hyperbolic sign question (outlook.md 1.1).

**6.4 The \(i\mathcal{S}\) companion and \(\mathrm{PGL}_2\).** Everything here used \(\zeta \equiv i \pmod 2\). The rotated family \(i\mathcal{S}\) (\(x\) odd, \(y\) even) doubles the census (euclidean-counting.md Remark 2); transporting the phase theory to it (and to the full \(\mathrm{PGL}_2\)-arrangement) should exchange the roles of the two reality lines and give the "missing" real classes at even \(n\). A clean consistency laboratory, like even levels were for the hyperbolic side.

**6.5 A lemniscatic phase atlas.** Plot the disks of the unit square colored by \(\arg u\) and \(\log|u|\) across curvatures — the Euclidean companion of outlook.md 1.6, likely exposing the R/I geometry of §5.2 and the class-group translation structure visually.

### Medium (weeks to months)

**6.6 Kronecker limit formula: \(\log|u|\) against \(L'\)-values of ring class characters of \(\mathbb{Q}(i)\).** \(\log|u_\mathfrak{c}|\) is, by the closed form, a linear combination of \(\log|\Delta(\Lambda_\mathfrak{c})/\Delta(\mathbb{Z}[i])|\) and logs of ring class integers. Character sums \(\sum_\mathfrak{c} \chi(\mathfrak{c}) \log|u_\mathfrak{c}|\) should evaluate to \(L'(0, \chi)\)-combinations by the Kronecker limit formula for orders — making the Euclidean phases a **lemniscatic elliptic-unit system in the ring class towers of \(\mathbb{Q}(i)\)** (Siegel–Robert units of the lemniscatic curve, with the \(\beta^4(\beta-1728)^3\)-dressing). This is the same experiment as outlook.md 2.3 but in the conductor aspect over a fixed field, where the elliptic-unit literature (Robert; Kubert–Lang ch. 12–13; Schertz) is strongest. PSLQ against \(\{L'(0,\chi), \log p\ (p \mid 2n)\}\) at \(n = 9, 11, 13\) is an afternoon.

**6.7 The translation law, made explicit.** Theorem 3's mechanism predicts \(\sigma_\mathfrak{c}(v^6_{\mathfrak{c}'}) = v^6_{\mathfrak{c}^{-1}\mathfrak{c}'}\) (Artin-normalized) and inversion by the mirror law — the generalized dihedral picture the computed Galois groups already display. Prove it via Shimura reciprocity for lattice functions (the adelic action on \(\hat h_2^6/\Delta\)-ratios), then push to \(v^2\) with the Weber multiplier bookkeeping of 6.2. This upgrades §5.4 from "observed groups" to a theorem, and is strictly easier than the hyperbolic §5.5–5.8 chain (no \(\varepsilon\), no \(\mu\)-cocycle, class group abelian over a fixed field).

**6.8 The transpose pairing and a two-sided theory.** The involution \(\sigma(X) = \bar X^{-1}\) maps our moduli problem to its transpose (§2), pairing each disk with a \(\sigma\)-disk whose lattice is built from the *columns*. On the product of the two problems \(\sigma\) is an honest involution; the composite invariants (e.g. \(\Theta(X)\overline{\Theta(\sigma X)}\), which unwinds to \(-j'(d/c)\,j'(a/c)/c^4\)) couple two conductors at once — the Euclidean shadow of the hyperbolic \(\hat\sigma\)-pairing and its simultaneous modular equations. Work out the class formula linking \([\Lambda_{\text{rows}}]\) and \([\Lambda_{\text{cols}}]\); the determinant relation \(ad - bc = 1\) should make it a conductor-mixing analogue of \(\hat\sigma[f] = [\mathfrak{r}_n][f]^{-1}\).

**6.9 The center as a torsion point: division values.** The left invariant (center mod \(\mathbb{Z}[i]\)) of a curvature-\(2n\) disk is a \(2n\)-torsion point \(\zeta/2n\) of the lemniscatic curve \(\mathbb{C}/\mathbb{Z}[i]\). The classification congruence \(x^2 + y^2 \equiv 1 \pmod{4n}\) selects a specific Galois-stable set of torsion points, and \(\wp(\zeta/2n;\, \mathbb{Z}[i])\) are division values living in ray class fields of \(\mathbb{Q}(i)\) — Abel's and Eisenstein's lemniscate division theory. Combine with §4: the full invariant of a disk is (torsion point, ring class point, phase) — a point on a mixed Shimura-type object (universal elliptic curve over \(X_0\)-data for \(\mathbb{Q}(i)\)). Identify the image: which (torsion, class) pairs occur is exactly the classification congruence — a reciprocity statement waiting to be named.

**6.10 Hecke comparison: the primitive part of \(\Phi_n(x, 1728)\).** The cyclic index-\(n\) family gives \(\Phi_n(x, 1728)\) (modular polynomial); Theorem 2 identifies its conductor-\(n\) part as \(H_{-4n^2}(x)^{2/\ldots}\)-structure and the lower strata as smaller conductors. Write out \(\Phi_n(x, 1728) = \prod_{f} H_{-4f^2}(x)^{e(f,n)}\) with explicit exponents (local line counts of Lemma 2), and similarly the \(\Delta\)-modular equation at \(\tau = i\) against the mass law — turning §5.5 into statements about classical modular polynomials at the elliptic point.

### Large (a paper or program each)

**6.11 Equidistribution of the phase in the conductor aspect.** As \(n \to \infty\), the normalized angles \(\arg v_\mathfrak{c}\) over \(\mathrm{Cl}(\mathcal{O}_n)\) sit on top of the conductor-aspect equidistribution of the CM points \(w_\mathfrak{c}\) (known by subconvexity/duke-type results for orders). The phase adds the homothety angle — a torus coordinate above the modular curve. Formulate and prove equidistribution of \((w_\mathfrak{c}, \arg v_\mathfrak{c})\) on the unit tangent-like bundle; the deterministic fiber rate 2 (§2) is the vertical part. This is the Euclidean twin of outlook.md 3.4, in the aspect where the harmonic analysis is cleanest.

**6.12 Heights: the mass as an arithmetic intersection number.** \(\log M(n)\) is a Faltings-height difference between the conductor-\(n\) locus and the lemniscatic point, and Theorem 4 — proved here by elementary modular means — has exactly the shape of an intersection with the supersingular locus (split primes: zero; non-split: geometric series). Re-derive it from the Gross–Keating/quasi-canonical side (thereby solving 6.1's per-class refinement on the way), and in the other direction interpret \(H_{-4n^2}(1728)\) (the same-field collision) via degenerate Gross–Zagier — the arithmetic self-intersection of the family with its own \(n = 1\) stratum. Together with 6.6 this would give the Euclidean phase a complete archimedean-plus-finite height story, with Theorem 4 as its proven finite-place skeleton.

**6.13 Other fields, and the two aspects at once.** Over \(\mathcal{O}_K\) for general imaginary quadratic \(K\) (Stange's arrangements), the Euclidean study transports verbatim: disks of curvature \(2n\)-analogues \(\leftrightarrow\) \(\mathrm{Cl}(\mathbb{Z} + n\mathcal{O}_K)\), period \(\Omega_K\) from the Chowla–Selberg formula, \(\mu_6\)-phases in the Eisenstein case (sharper torsion tests). Combined with the hyperbolic study one obtains, for each \(K\), *both* CM aspects — discriminant (hyperbolic grading) and conductor (Euclidean grading) — realized by one circle packing. The uniform statement ("a Schmidt arrangement is a geometric model of the full CM theory of its field: orders in one direction, discriminants in the other") is the definitive form of the project, and the two-aspect interaction (6.8, 6.10) its novel content.

**6.14 The weight-\(3/2\) two-slice.** The hyperbolic census slices Zagier's form along \(d = n^2 - 1\), the Euclidean along \(d = 4n^2\). Shifted squares and even squares are precisely the two quadratic progressions occurring in the Eichler–Selberg/Hurwitz–Kronecker relations \(\sum_t H(4m - t^2)\); a bijective proof of those relations *on the arrangement* (outlook.md 3.3) would now have both slice families available as concrete disks and circles, with the \(t = \pm 2n'\)-terms as Euclidean strata. A single geometric identity tying the two censuses would be a genuinely new take on class number relations.

## 7. Files

- [scripts/euclidean_moduli_invariants.py](scripts/euclidean_moduli_invariants.py) — all experiments. `python3 scripts/euclidean_moduli_invariants.py` runs A (exact structure, \(n \le 24\)), B (invariance, fiber, period, closed form), C (ring class polynomials, \(H^2\), traces, \(n \le 13\)); `... phase` runs D (laws, reality, \(P^{(2)}, P^{(6)}\), \(\lambda_n\), factorizations, Galois groups), D' (prime tagging), D'' (the \(\Delta\)-mass, including \(n = 25, 27, 49\)); `... all` runs everything. Requires mpmath, sympy.
- [scripts/mass_law_and_irreducibility.py](scripts/mass_law_and_irreducibility.py) — proof verification for Theorem 4 (`... mass`: M1 at random \(\tau\); the exact recursion vs. the closed form, \(n \le 60\); the symbolic generating-function identities; the class-by-class sign analysis, \(n \le 20\)) and the irreducibility suite of §5.6 (`... irred`: \(H\), \(P^{(2)}\), \(P^{(6)}\), \(P^{(2)}(x^2)\) factored exactly with squarefreeness certificates, \(2 \le n \le 16\)).
- Certification policy: integers/rationals are accepted only with \(\ge \max(20, \mathrm{dps}/5)\) spare digits in the **absolute-error** sense (a huge value with \(O(1)\) fractional part is rejected no matter how many digits it has — the guard rail learned in moduli-invariants.md §5.6, one trap deeper).
- Precisions: structure exact (integer arithmetic); analytic experiments at 120–840 digits per level (table `DPS` in the script); Theorem 4's recursion and all factorizations in exact arithmetic.
