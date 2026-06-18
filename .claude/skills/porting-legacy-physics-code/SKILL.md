---
name: porting-legacy-physics-code
description: Port physics modifications from a legacy scientific-code fork (e.g. an old CAMB/Boltzmann-code fork) onto a modern upstream codebase. Use when migrating physics modules between code versions, modernizing Fortran/C research code, validating such a port, or wiring the result into a sampler framework (Cobaya/CosmoMC/MontePython). Covers verbatim-numerics discipline, regime-complete ratio validation, bit-identical upstream limits, and inline change-marking.
---

# Porting legacy physics code onto a modern codebase

Distilled from a real port: AxiECAMB (ultralight-axion physics, built on CAMB Nov-2013)
onto modern CAMB 1.6.7 (class-based Fortran + Python wrapper), validated to 0.01–0.1%
against the original and wired into pristine Cobaya. Examples below reference that port;
the method is general.

## Core principles

1. **Triangulate, never read the fork raw.** Obtain the *pristine base* the fork was
   built on. The physics is `diff(pristine_base, fork)` — usually 10–30× smaller than
   the fork itself. Measure real size with whitespace-insensitive diffs before planning.
2. **Classify every change before porting any change.** Tag each diff hunk:
   `[PHYSICS]` (port verbatim), `[PLUMBING]` (re-derive for the new architecture),
   `[ACCURACY]` (port if physics-motivated), `[OBSOLETE]` (superseded upstream),
   `[COSMETIC]` (ignore). Most of a legacy fork is not physics.
3. **Find the chokepoints in the modern code.** Modern architectures centralize what old
   ones scattered. One background-density hook (`grho_no_de` → `dtauda`) made ~90% of the
   fork's recombination/reionization/distance edits obsolete — they flowed through
   automatically. Count chokepoints before counting changes.
4. **Numerics verbatim, plumbing restructured.** Copy Butcher tableaux, tolerances,
   magic constants, iteration orders, even crude spline boundary conditions exactly —
   reproducing the original is the only ground truth you have. Restructure *only* state:
   globals → state objects, sentinels → logicals at the API surface (keep internal
   semantics), `stop` → recoverable errors.
5. **Gate everything; prove the null limit.** Every modification must be behind a
   component-active check so that with the new physics off, output is **byte-identical**
   to upstream. This is your strongest regression test and it is binary.
6. **Mark every hunk.** One greppable tag (e.g. `!AxiECAMB`) on every change site, and a
   per-hunk audit script (`diff` vs pristine, assert each hunk contains the tag). The
   diff-vs-upstream then doubles as the change documentation.

## Workflow

### Phase 0 — Baseline the toolchain
Build and run the *unmodified* modern code first. Record build quirks (serial-make
races, SDK paths, compiler env) before they can be confused with your bugs. Save
reference outputs for later byte-comparisons.

### Phase 1 — Inventory the fork
Diff every file against the pristine base. For each modified/new file produce a written
report: every non-cosmetic change with file:line, **verbatim equations** (never
paraphrase math), units, and a classification tag. Check the build system for what is
*actually compiled* — forks ship dead files (see pitfalls). Parallel analysis agents
work well here; make the reports the durable artifact (they become the implementation
spec and, later, developer documentation).

### Phase 2 — Map the modern architecture
For each class of change, find its landing zone in the modern code: parameter structs
and their language-binding mirrors, background hooks, perturbation-equation extension
points, switch/approximation machinery, build system, Python/ini surfaces. Identify
existing precedents to copy (modern codes usually contain an analog of what you need —
e.g. quintessence solvers, mid-evolution variable switches for tight coupling).

### Phase 3 — Design document
A change-by-change mapping: old location → new location → decision. Include an explicit
**"NOT ported (with reasons)"** list — it is as valuable as the port itself. Decide the
component architecture here (own component vs. hijacking an existing slot; see
patterns.md "Semantic honesty of slots").

### Phase 4 — Implement in compile checkpoints
Order work so the code compiles and passes the null-limit test at named milestones
(e.g. background-only first, then perturbations, then sources/integrators, then
wrappers). Build after every checkpoint; run the byte-identical null test every time.

### Phase 5 — Validate by regimes and by ratios
- **Ratios, not absolutes**: compare (modified/base) observable ratios between old and
  new code, so version-to-version baseline differences cancel. Define a quantitative
  target (we used ≤0.1% on C_ℓ ratios, ≤0.01% on P(k) ratios).
- **Regime-complete matrix**: one test per physical branch and per code path — every
  switch position (early/late/never), every parameterization branch, extreme prior
  corners (including pathological ones like negative effective Λ).
- **Frameworks end-to-end**: run the upstream test suite; for sampler integration,
  execute a real evaluate + a micro-MCMC, not just an import.

### Phase 6 — Document and review
Ship the Phase-1 reports and Phase-3 design as developer docs. Write a critical-review
list: every constant you could not re-derive, every parameter-space discontinuity, every
accuracy knob that does not scale — ranked for human (author-level) attention. You are
the person who has read this code most closely; say what you didn't like.

## Files

- [patterns.md](patterns.md) — the load-bearing patterns, each with the concrete case
  that earned it.
- [pitfalls.md](pitfalls.md) — the catalog of traps, each with the bite it gave or
  nearly gave.
