# extra/ — material not on the program's line

Everything in this folder is finished, machine-verified work of the project
that is **not needed for the program** of [../PROGRAM.md](../PROGRAM.md) (the
elliptic-unit / Euler-system line from the classification to the Robert index
and beyond). It is kept intact, with its scripts and figures, and its
cross-links have been rewritten for the new location.

| document | layer | why it is here |
|---|---|---|
| [circle-composition.md](circle-composition.md) | foundational | Gauss composition in circle language (Paper I §5); the program computes class groups from lattices instead |
| [half-plane-monoid.md](half-plane-monoid.md), [atomic-census.md](atomic-census.md), [product-cocycle.md](product-cocycle.md), [middle-kernel.md](middle-kernel.md) | monoid / operation layer | the monoid \(\Omega\), Apollonian addresses, the product cocycle and the middle kernel — excluded from both papers by design |
| [spectral-geometry.md](spectral-geometry.md), [spectral-outlook.md](spectral-outlook.md) | spectral layer | the Picard-orbifold survey (its §14 holds the early literature probes cited by the papers' NOTES) and the trace-formula program |
| [phase-atlas.md](phase-atlas.md) | phase layer | phase portraits, the divisor-class sign law and the Duke–Imamoğlu–Tóth non-relation — Paper I §6 material, not on the unit line |
| [prompts/02-schmidt-trace-formula.md](prompts/02-schmidt-trace-formula.md), [prompts/03-phase-atlas-and-dit.md](prompts/03-phase-atlas-and-dit.md) | prompts | the spectral prompt (staged) and the atlas prompt (done) |

Scripts (`extra/scripts/`): `omega.py`, `omega_verify.py`, `atom_invariants.py`,
`apollonian_chain_demo.py`, `atomic_census.py`, `product_cocycle.py`,
`middle_kernel.py`, `hurwitz_sum.py`, `composition_check.py`,
`matrix_composition_check.py`, `make_composition_figure.py`,
`make_omega_figure.py`, `make_phase_atlas.py`, `dit_comparison.py`. Run them
from the repository root, e.g. `python3 extra/scripts/atomic_census.py`; the
ones that import the main verification scripts (`composition_check`,
`matrix_composition_check`, `make_phase_atlas`, `dit_comparison`) put
`../scripts` on their path themselves.

Figures (`extra/figures/`): the phase-atlas pages and `omega-atoms.png`. The
three figures used by Paper I (`alpha15-disk.png`, `composition-n9.png`,
`phase-atlas-n9.png`) stay in the main `figures/` folder, which is what the
paper's `\graphicspath` points to; `make_phase_atlas.py --figures` writes to
`extra/figures/` by default.

Nothing here is deprecated: the literature verdicts of Paper I's NOTES cite
`spectral-geometry.md` §14, and Paper I's sections 5–6 rest on
`circle-composition.md` and `phase-atlas.md`.
