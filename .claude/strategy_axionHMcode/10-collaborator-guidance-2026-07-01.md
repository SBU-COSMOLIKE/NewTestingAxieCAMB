---
name: axionhmcode-collaborator-guidance
description: Collaborator design guidance (2026-07-01) on the open questions — version default, validity-domain physics, switch-redshift risk, ratio convention, compatibility constraint — with epistemic status and follow-up leads
metadata:
  type: project
---

# Collaborator design guidance (2026-07-01)

Design guidance gathered from group physics discussion on 2026-07-01, addressing the open
questions in file 06. Recorded here in summary form; further, more specific input from
the collaboration is expected and is tracked under "Pending".

## Epistemic status

This guidance reflects physics judgment, not hands-on experience with the axionHMcode
implementation, and the collaboration flagged exactly that limitation when giving it.
Consequence: every decision sourced from this file (dome default, exposed nuisance
parameters, the smooth-component reasoning, risk R10) is working guidance, not
code-verified fact. It must be re-verified against the actual code and papers during
Phases 0-1, and any disagreement between this guidance and the code/papers is itself a
finding to surface to the group — not a conflict to resolve silently in either direction.

## Version default and nuisance parameters (question 1)

dome is the most recent recalibration and is adopted as the default `version`. Exposing
the calibration parameters is considered useful: the Dentler nuisance parameters
(alpha_1, alpha_2, gamma_1, gamma_2) become Theory-block parameters that can be sampled
(defaults = the calibration values; the example MCMC shows the Gaughan-et-al prior
blocks, commented or active per run).

## Validity-domain physics (question 2)

- Conditions (a) DE-like axion and (b) large fax are not fundamentally out of domain:
  a proper smooth-vs-clustered decomposition of the axion would cover them, though
  warnings remain appropriate. For (a), the axion has horizon-scale smoothness, so two
  routes exist in principle: treat the DE-like axion analogously to w0wa smooth dark
  energy, or as a smooth component of mixed dark matter. Either route requires pinning
  down what "matter" means in every convention layer (CAMB Transfer_tot, axionHMcode
  Eq. 9 weights, the boost denominator) and how that maps to the Weyl potential sourcing
  the lensing. That is a physics extension, not plumbing — DEFERRED pending further
  collaborator input; the current build keeps the hard error for (a) as a pragmatic
  guard, explicitly documented as not-fundamental (see [[axionhmcode-architecture]]
  "Regime gating").
- (b) stays warn-under-`strict: False`, consistent with this guidance.
- (c),(d) warn — but with one real breakdown identified: redshifts above the KG→EFA
  switch (z > z_osc, i.e. a < a_osc, before the axion oscillates and becomes DM-like).
  There the axion transfer function is not a DM density contrast in the halo-model sense
  at all. Logged as risk R10 in [[axionhmcode-open-questions]]: the Theory class must
  check the top of its z grid against z_osc = 1/a_osc - 1 (available as
  `results.Params.Axion.a_osc`, camb/axion.py:36). For the target window
  (m ~ 1e-25..1e-23 eV) oscillation begins deep in the radiation era (z_osc >> 1e4), far
  above any lensing-source redshift, so this assertion should never fire in production —
  but it must exist, because a light-mass misconfiguration would silently feed garbage
  otherwise.

## Ratio convention (question 6) and compatibility constraint (question 4)

- The Eq.9-consistent ratio proposal (numerator AND denominator in axionHMcode's own
  decomposition) is assessed as plausible — adopted as the working convention, pending
  further collaborator input.
- Question 4 gains a design constraint: DRAG-AND-DROP COMPATIBILITY with future
  axionHMcode updates. The Theory class must call axionHMcode only through its public
  entry points (func_power_spec_dic-equivalent construction, HMCode_param_dic,
  func_axion_param_dic, func_full_halo_model_ax) and never modify axionHMcode files;
  if a group fork is created it is pin-only (no code changes). The Eq.9-internal ratio
  definition supports this: all conventions live on our side of the interface.

## Follow-up leads on the switch problem (R10 area)

Directions raised in the same discussion, recorded as investigation leads (not decisions):

- Phase-1 check: how the public AxiCAMB (github.com/adammoss/AxiCAMB, the code base of
  arXiv:2605.12054) handles redshifts near/above the EFA switch in its axionHMcode
  nonlinear interface — MCMCs were shipped with it, so either the problem was hit and
  solved there, or its z range never reaches the switch.
- Possible resolution: output instantaneous density transfer functions in the KG phase,
  accepting a definitional discontinuity at the switch. Phase-1 sub-task: establish what
  New_AxiECAMB's Transfer_axion actually returns for z > z_osc (field phase) — the port
  docs document careful source-discontinuity handling at the switch (README.rst:2843
  ultra-fine time window; :4114/:4259 metric-delta source pieces; :4890 dtauda kink
  split), but the transfer-output definition in the KG phase is unverified.
- A "blend" across the switch (an idea from the original port-era discussions; not found
  in the committed port docs — no written spec exists) could make weak-lensing observables
  behave continuously across it. Would need a written definition before implementation.

None of this blocks the target window (z_osc >> any lensing redshift for
m ~ 1e-25..1e-23 eV); it matters if the pipeline is ever pushed to lighter masses.

## Pending

- Further collaborator input (expected: the smooth-component conventions for
  DE-like/large-fax axions and their Weyl-potential mapping; possibly ratio-convention
  details for question 6).
- The switch-problem leads above (Phase 1).
- Dome flag configuration question (found 2026-07-02 on full re-read of 2409.11469;
  PI decision: keep released-code defaults, documented in
  axionhmcode_boost/README.md "The dome flag configuration" section). Dome et al.
  Sec. 4.4 state the calibration ran with the halo bloating eta (their Eq. 57) and
  the perturbative two-halo damping (their Eq. 55) active, so the Table 5 alpha fits
  were calibrated jointly with them; the released code's README instead recommends
  only `alpha = True` and `concentration_param = True` for the dome version and
  leaves eta/two-halo damping off. Our defaults follow the released code and
  reproduce Gaughan et al. Figs. 1-3 (V3). Question for collaborators/upstream:
  should dome runs set `eta_given=True` and `two_halo_damping=True` to match the
  Sec. 4.4 configuration? If yes, no code change:
  `model_flags: {eta_given: True, two_halo_damping: True}` in the yaml. Expected
  size of the difference: a few percent in the quasi-linear regime, inside the
  model's 10-20% accuracy band.
