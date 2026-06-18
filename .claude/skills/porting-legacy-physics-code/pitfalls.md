# Pitfalls

Each entry: the trap, how it bites, and the concrete encounter in the AxiECAMB →
CAMB 1.6.7 port (bitten, or disarmed just in time).

## Units and conventions

### P1. Mid-routine base changes
A table abscissa held ln(a) during integration and was converted **in place** to
log10(a) before output (`loga_table = dlog10(dexp(loga_table))`), with `log(10)`
Jacobians appearing only in the spline boundary derivatives. Every downstream lookup
assumes log10. Flagged by analysis as "the single most likely source of porting bugs" —
correctly. **Defense:** write a units ledger per array (name, units, base, valid range)
before porting; keep the original convention even if ugly, and comment it at the
declaration.

### P2. Misleading names encode units
A function named `littleh` returned **conformal aH in units of 100 km/s/Mpc** — neither
little-h nor H. Three different Hubble normalizations (H/H0, aH/100km/s/Mpc, 1/Mpc)
coexisted, with conversion factors like `hnot/littlehfunc` scattered at call sites.
**Defense:** never trust names; derive the unit of every quantity from one verified
call chain, then write it down.

### P3. Framework conventions differ by factors of 100
The chain sampled 100·θ★ while the modern API expects θ★ itself; neither the pristine
nor the forked sampler wrapper divides by 100. Caught only because a test point
produced "theta looks wrong" instead of silently shooting H0 to a wrong root.
**Defense:** for every "same-named" parameter crossing a framework boundary, test one
point with a known answer before wiring the chain.

## State, sentinels, and control flow

### P4. Magic sentinel values as control flow
`a_osc = 15` meant "no oscillation"; `omaxh2_guess = 42` meant "not computed";
`> 1.0` meant "recompute" (and 42 > 1 makes the sentinels overlap by design). A
later consumer capped a_osc to 1 — *two files away*. Porting any one site without the
others breaks the protocol invisibly. **Defense:** map the full lifecycle of every
sentinel across files before touching it; expose logicals at the new API surface but
keep internal semantics until validation passes.

### P5. Stateful scratch inside parameter structs
A fixed-point coefficient (wEFA_c) was initialized per solve, mutated inside an inner
routine, and **fed back** through the parameter struct into subsequent iterations —
bit-level reproducibility depends on call order. **Defense:** reproduce the call order
exactly first, validate, and only then consider refactoring (we kept it).

### P6. Garbage-but-unused fields that downstream code later uses
In the no-switch regime the final matching ran at the a_osc=15 *sentinel*, splining
1.2 decades beyond the table — producing garbage values that were "never used"… except
one fed a perturbation rescale factor, which only worked because the garbage happened
to stay positive. **Defense:** for every "unused in this regime" field, grep all
consumers; replace regime-invalid values with explicit safe ones.

## Dead code and build truth

### P7. The plausible stale twin
The fork shipped `cmbmainOMP.f90` — looking like "the OpenMP variant" of the main loop.
It was an older snapshot missing the newest physics (the boundary-term machinery) and
was **not in the Makefile**. Using it as a porting reference would have silently dropped
a published feature. **Defense:** the Makefile is the only source of truth for what
exists; diff twin files against each other before trusting either.

### P8. Files that never compiled
A modified `writefits.f90` contained syntax errors (orphan `end if`, undeclared
symbols) — it had never been compiled since the edits (its target was disabled).
Analysis time spent on it: wasted unless you check compilability first.

### P9. Wrappers shipped but superseded
The fork carried a verbatim DVERK clone added "solely to thread scalars" — modern code
passes state objects, making the entire mechanism obsolete. ~90% of that file's diff
was this plus re-indentation. **Defense:** classify before reading line-by-line.

## Numerics

### P10. Discontinuities in the background propagate everywhere
A derivative kink in dtauda at the switch degrades every Romberg integral crossing it.
The fork patched **three** sites ad hoc (time, sound horizon at z★, reionization) and
missed others (r_drag — a BAO observable!, optical depths, damping scale). **Defense:**
fix kinks at one helper (split-at-kink integrator) and route *every* background
integral through it; grep the integrator name to find them all.

### P11. Accuracy knobs that don't scale
The fork's ULA-specific accuracy settings (table resolution = 100·dfac, fixed-step RK8
with the embedded error estimate computed and *discarded*, hand-tuned time-step
windows) do **not** respond to the global AccuracyBoost — so standard convergence tests
exercise everything except the new physics. **Defense:** document which knobs are
frozen; if possible add one boost that scales the new machinery's resolution/tolerances
together.

### P12. Hard regime boundaries inside the prior
A binary DM-like/DE-like classification at m/H0 = 10 switches the *definitions* of
σ8, halofit's Ω_m, and the total transfer function discontinuously — inside the chain's
prior range. Current likelihoods happened to be continuous through it; any P(k)-based
likelihood would not be, and **emulators trained on the output will differentiate
through the step**. **Defense:** inventory every gate/threshold; verify which
observables jump; document the boundary as part of the model definition.

### P13. Trial failures are not sample failures
See patterns.md #13. Promoting a per-sample rejection flag to a hard error — a
seemingly safe "improvement" — broke θ★→H0 shooting at its bracket endpoints, where
under-dense *trial* candidates legitimately have no expanding solution (the Friedmann
RHS, the would-be (aH)², goes negative: the trial universe recollapses before reaching
that scale factor; physical H² is never negative). The original's NaN-circulation
through the bisection was sloppier but *less wrong* than the strict port.

## Language-binding mirrors

### P14. Positional ctypes mirroring fails silently
The Python mirror of the Fortran parameter struct is matched **by position, not name**;
one missing or misplaced field shifts every later field with no error. Appending new
fields at the end (both sides) minimizes risk; the framework's struct-size handshake
does not catch ordering errors. **Defense:** append-only placement, then run the
upstream Python test suite — it is the only layout guard.

### P15. Non-idempotent input transformations break round-trips
Reading `Hinf` as log10(GeV) and *storing* the converted ratio made write-ini→read-ini
non-idempotent: the suite's round-trip test failed with a doubly-converted value.
**Defense:** store raw inputs; convert at the point of use. Any `x = f(x_input)` at
read time is a round-trip bug waiting for whichever field has a non-involutive f.

## Process

### P16. Log spam is a correctness issue in samplers
A once-per-evaluation warning becomes millions of lines in an MCMC. Per-process
`logical, save` once-flags for any diagnostic on the likelihood path.

### P17. Build-system races and environment look like code bugs
Parallel make raced on Fortran module files (twice, two layers); conda gfortran on
macOS fails to link without `SDKROOT`. Both were diagnosed *before* any code edits —
because the pristine baseline was built first (patterns: Phase 0). Record quirks in
the repo docs.

### P18. Validating with the easy case only
The first axion test (switch before recombination) agreed to 0.01% — *with the
boundary-term machinery entirely missing*, because exp(−κ)≈0 suppressed it there.
Only the post-recombination-switch regime exercised that code. **Defense:** for every
mechanism you port, identify the regime where it is *largest* and test there
(patterns #6).

### P19. The closure remainder assumption
Any expression of the form `rest = 1 − Σ(known components)` inherits a scaling
assumption about `rest`. Upstream HMcode's remainder is curvature-like (1+z)² — valid
for its entire supported model space, silently wrong the moment your port adds a
component to the budget. Your new component breaks *other people's* closure
assumptions in modules you never edited. **Defense:** after adding any component, grep
the codebase for `1 -` / `1._dl -` budget expressions and audit each.
