# Patterns

Each pattern: the problem, the move, and the concrete case from the AxiECAMB → CAMB 1.6.7
port that earned it a place here.

## 1. Triangulation diff

**Problem:** A research fork is thousands of lines of re-indented, comment-laden code;
reading it raw conflates physics with noise.
**Move:** Acquire the exact pristine base the fork branched from. Generate per-file
unified diffs; measure with `diff -w -B | grep -c '^[<>]'` to size *real* changes.
**Case:** `cmbmain.f90` showed 3124 changed lines raw; 543 non-whitespace; ~163 real
after de-commenting — and exactly one new mechanism (a boundary term). Plans made from
the raw number would have been 20× oversized.

## 2. Classification tags as decisions

**Problem:** Without forced classification, every change defaults to "port it," and you
inherit ten years of workarounds for problems the modern code no longer has.
**Move:** Tag every hunk `[PHYSICS]/[PLUMBING]/[ACCURACY]/[OBSOLETE]/[COSMETIC]` during
analysis, before implementation. Each tag is a decision with a default action.
**Case:** The fork's entire recombination file (2929 diff lines) reduced to *one*
`[PHYSICS]` item (an axion dHdz term); everything else was `[OBSOLETE]` because the
modern code routes H(z) generically. An elaborate xe-smoothing block looked like physics
and was a disabled no-op (`[OBSOLETE]`, dead upstream comment confirmed it).

## 3. Chokepoint accounting

**Problem:** Old architectures forced the fork to patch H(z), Ω_m, etc. in a dozen
places; porting each patch 1:1 multiplies risk.
**Move:** Before porting, list where the modern code *centralizes* each quantity. Port
to the chokepoint once; then *verify consumers inherit it* (don't assume — read each).
**Case:** Adding the axion density to `grho_no_de` made distances, thermal history,
reionization, tensor background, *and* the CosmoRec/HyRec wrappers axion-aware with zero
additional edits (the wrappers pass CAMB's tabulated H(z), not parameters). The same
audit found the one place that did NOT inherit: recfast's analytic dHdz.

## 4. Verbatim numerics, restructured state

**Problem:** "Improving" numerics while porting destroys your only ground truth.
**Move:** Transcribe integrator tableaux, tolerances, bracket logic, iteration order,
spline boundary conditions, even crude finite-difference BCs, exactly — and preserve
behavior quirks (e.g. a fixed-point variable that's stateful across calls). Restructure
only ownership: module globals → component/state members.
**Case:** The 16-stage RK8 tableau, the 1e-6 shooting tolerance, the wEFA fixed point
(tol 1e-7, max 30), and the deliberately-crude equality-spline BCs were copied digit for
digit. Result: a_osc matched to 5 digits, σ8 to 4, on first successful run.

## 5. Ratio validation

**Problem:** Old-base and new-base codes differ at the percent level for *baseline*
physics (recombination versions, ν treatment); absolute comparisons can't separate
port errors from version differences.
**Move:** Validate `O_modified/O_base` ratios computed *within* each code, old vs new.
Baseline differences cancel; the new physics is isolated. Set numeric targets.
**Case:** TT C_ℓ suppression ratios agreed to ≤0.01–0.1% across ℓ=2–2600 in all regimes
while absolute baselines differed by ~0.2%. The residual ratio difference appeared
equally in a regime where the new physics was negligible — diagnosing it as baseline
non-cancellation, not port error.

## 6. Regime-complete test matrix

**Problem:** Branchy physics (switches, parameterization modes) hides untested paths.
**Move:** Enumerate every physical regime and every input branch; one matched old/new
test per cell. Include hostile prior corners.
**Case:** Five cells: switch before recombination, after recombination (activates a
boundary term that is exp(−κ)-suppressed otherwise — invisible in the "easy" cell!),
no-switch (DE-like), fraction-parameterization branch, and 100%-replacement. A
negative-Λ corner (frozen field supplying 90% of DE) was added later as a stress test.

## 7. Gated changes + byte-identical null limit

**Problem:** Diffuse small edits silently shift baseline results for users who never
enable the new physics.
**Move:** Gate every edit on component activity. CI test: with the component off,
outputs must be byte-identical to unmodified upstream (`diff` the output files, not
"close enough").
**Case:** Proven twice — for the main port (ΛCDM limit) and again for the HMcode
adaptation (ΛCDM + HMcode limit). The byte test caught nothing both times *because* the
discipline was applied from the start; its value is that it is binary and free.

## 8. Inline markers + per-hunk audit

**Problem:** Reviewers (and future you) cannot find what changed, and prose docs drift.
**Move:** One tag (`!AxiECAMB`) on every change site. Then a script that diffs each file
against pristine upstream, splits hunks, and asserts every hunk contains the tag.
**Case:** The audit found 9 unmarked hunks out of 74 on first run (use-statements,
one-line density terms, block tails). Cost to fix: minutes. Result: `grep -n AxiECAMB`
*is* the change map.

## 9. Reports as spec, then as documentation

**Problem:** Analysis effort evaporates if it lives in a conversation or a head.
**Move:** Write per-subsystem analysis reports with verbatim equations and line refs to
durable files. Implement *from the reports*, not from re-reading the fork. Afterwards,
convert the reports into the repo's developer documentation.
**Case:** Seven analysis reports (~230 KB) drove the implementation, then became the
README's "Port documentation" section nearly unchanged. The "hardest item" (a 130-line
variable-projection at a mid-evolution switch) was ported from the report's quoted code
without reopening the original file.

## 10. Semantic honesty of component slots

**Problem:** It is tempting to hijack an existing extension slot (e.g. the dark-energy
class) for new physics because the hooks are free.
**Move:** Ask what every *consumer* of the slot believes it means (derived parameters,
w(a) reported to nonlinear modules, density-row outputs). If the new component is
sometimes-matter, a dark-energy slot lies to all of them. Prefer a first-class component
even at the cost of explicit wiring.
**Case:** The axion as a separate `CAMBparams%Axion` kept Ω_de, halofit's
`Effective_w_wa`, and background-density outputs honest, and preserved the user's
freedom to combine axions with fluid/PPF dark energy — at the cost of ~10 explicit,
auditable hook sites in the equations module.

## 11. Three meanings of Ω_m (audit all consumers)

**Problem:** "Matter density" is overloaded: expansion source, P(k)/σ normalization,
and halo-collapse density are *different quantities* once a component half-clusters.
**Move:** Grep every consumer of the density variable in the target module; assign each
consumer to a meaning; decide each meaning explicitly per regime. Watch for *closure
remainders* — expressions like `(1 − Ω_m − Ω_v)` that silently assume the budget is
exhausted.
**Case:** HMcode's `Hubble2` assigns its closure remainder a curvature-like (1+z)²
scaling — correct for every model upstream supports, silently wrong the moment a new
component enters the budget. The fix audit found seven distinct consumers (expansion,
acceleration, Ω_m(z), Ω_cold(z), f_ν, σ(R)↔M mass mapping, Dolag ΛCDM reference), each
needing its own regime decision.

## 12. Additive exact component, not wholesale replacement

**Problem:** When feeding exact background tables into a fitting-calibrated module
(halo models, emulators), replacing the whole internal background changes the module's
behavior for *standard* cosmologies and breaks its calibration conventions.
**Move:** Add the new component as a separate exact term; leave the calibrated
parametric pieces untouched; remove only the new component's share from closure
remainders. Gate on activity.
**Case:** HMcode's internal expansion deliberately omits radiation ("consistent with
simulations") and normalizes growth as g(a)=a in matter domination. Wholesale
CAMB-exact H(a) tables would have broken both. The additive Ω_ax(a), w_ax(a) term kept
upstream behavior bit-identical for non-axion runs and exact for axion runs.

## 13. Error semantics follow the caller

**Problem:** Old code's error handling encodes its *original calling context*; the
modern context may invert the right behavior.
**Move:** For every error path you port, ask: who calls this now, and what do they do
with the failure? Per-*sample* rejection (MCMC driver discards the point) is not
per-*trial* rejection (a root-finder probing extreme intermediate values *expects*
some trials to fail and needs a graceful verdict, not an abort).
**Case:** The fork's collapse flag was a CosmoMC sample-rejection signal. Modern
θ★→H0 shooting (brentq inside CAMB) evaluates bracket endpoints at extreme H0 where
deliberately under-dense bisection candidates have no expanding solution. Promoting the
flag to a hard error — which felt like an improvement — killed every chain point. The
correct port: failed trial ⇒ "under-dense, push the bracket", hard error only for an
invalid *converged* solution.

## 14. Adapt your package to the framework, not the framework

**Problem:** Patching the sampler framework (Cobaya/CosmoMC) creates a fork users must
maintain forever.
**Move:** Find the framework's discovery hooks and implement them in *your* package.
Test against the pristine framework from PyPI, not your local copy.
**Case:** Cobaya's CAMB wrapper discovers parameters via the camb package's own
`get_valid_numerical_params` and sets them via `set_params`. Adding the axion setter to
those two chains (ordered *before* `set_cosmology`, so θ→H0 shooting sees the axion)
made unmodified Cobaya work. Zero framework patches; verified by running evaluate and a
micro-MCMC under pristine Cobaya.

## 15. The critical-review deliverable

**Problem:** The porter ends up knowing the code better than anyone but the authors —
and that knowledge evaporates if the port ends at "tests pass."
**Move:** Deliver a ranked list: (a) constants you could not re-derive from the paper,
(b) discontinuities in parameter space (gates, hard regime boundaries) that samplers
and emulators will differentiate through, (c) accuracy knobs that don't scale with the
global accuracy setting, (d) failure modes that are silent. Flag the top items for the
original authors explicitly.
**Case:** Out of the whole port, exactly two equations resisted re-derivation (an 11/10
opacity compensation and a 5/4 sound-speed coefficient with visible calibration
history). Naming precisely those two — rather than vague unease — produced a
substantive author discussion that adjudicated each suggestion on physics grounds.
