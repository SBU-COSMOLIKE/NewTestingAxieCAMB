# Strategy: axionHMcode boost factor → AxiECAMB via Cobaya (Cocoa)

Created 2026-07-01 from the prompt draft `prompt_draft_hmcode.tex`. Read in order.
All file:line references verified against the working tree on 2026-07-01.

- [01-objective-and-deliverables.md](01-objective-and-deliverables.md) — what we are building and the constraints
- [02-key-discovery-no-circular-dependency.md](02-key-discovery-no-circular-dependency.md) — the PR#480 mechanism already solves the CAMB↔axionHMcode cycle
- [03-proposed-architecture.md](03-proposed-architecture.md) — the Theory-class design, data flow, and staging plan
- [04-verified-facts.md](04-verified-facts.md) — inventory of code facts with file:line citations
- [05-axionhmcode-api.md](05-axionhmcode-api.md) — axionHMcode call sequence, dictionaries, conventions
- [06-risks-and-open-questions.md](06-risks-and-open-questions.md) — decisions that need the PI, and known traps
- [07-validation-plan.md](07-validation-plan.md) — regime-complete validation matrix
- [08-implementation-phases.md](08-implementation-phases.md) — phased plan with checkpoints (no code written yet)
- [09-mass-prior-yaml-recipes.md](09-mass-prior-yaml-recipes.md) — PI-decided mass window (1e-25..1e-23 eV) and README-ready yaml recipes: fixed mass vs sampled log-mass
- [10-collaborator-guidance-2026-07-01.md](10-collaborator-guidance-2026-07-01.md) — collaborator design guidance: dome default + exposed nuisance params, smooth-component physics for DE-like/large-fax, z above the KG→EFA switch (R10), drag-and-drop axionHMcode constraint

Related prior work: the AxiECAMB port to CAMB 1.6.7 (see auto-memory `axiecamb-port-project`,
this repo's `README.rst`, and the `hmcode` branch with the halofit-level axion-aware HMcode).
Key papers live in `../../papers/` (rayne/papers/, outside this git repo — large PDFs, not
committed; see 01 for the arXiv ids to re-download).

Location note: this folder lives inside the New_AxiECAMB repo (moved from rayne/strategy/
on 2026-07-01 so it can be committed); working-tree paths like `cocoa/...`, `axionHMcode/...`
in these notes are relative to the parent `rayne/` directory, one level up from the repo root.
