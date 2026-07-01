---
name: axionhmcode-collaborator-guidance
description: Wayne Hu's 2026-07-01 email comments on the open questions — verbatim record plus resulting decisions; further specifics expected from Dan
metadata:
  type: project
---

# Collaborator guidance (Wayne Hu, email relayed by PI 2026-07-01)

Verbatim record of Wayne's comments on the open questions in file 06 (email of
2026-07-01; further specifics expected from Dan). Interleaved formatting from the email
is flattened here; each block is Wayne replying to the numbered question.

## On question 1 (basic vs dome, nuisance parameters)

> Dome I think is the most recent recalibration. Exposing params is useful.

DECISION: `version: dome` is the default; the Dentler nuisance parameters
(alpha_1, alpha_2, gamma_1, gamma_2) are exposed as parameters of the Theory block so
they CAN be sampled (defaults = the calibration values; the example MCMC shows the
Gaughan-et-al prior blocks, commented or active per run).

## On question 2 (out-of-validity policy)

> in principle there doesn't need to be restrictions for (a), (b) if the axions can be
> treated as a smooth vs clustered component though warnings may still be good. For (a)
> the axion has a horizon scale smoothness so in principle could be treated in the same
> way as w0wa but could also be a smooth component of mixed dark matter so long as we
> know what the conventions for "matter" are and how they are related to the Weyl
> potential for lensing. (c),(d) warn but there may be a problem with redshifts above
> the switch.

Interpretation and consequences:

- (a) DE-like and (b) large fax are not fundamentally out of domain — a proper
  smooth-vs-clustered decomposition would cover them. For (a) two routes exist in
  principle: treat the DE-like axion like w0wa smooth dark energy, or as a smooth
  component of mixed dark matter. Either requires pinning down what "matter" means in
  every convention layer (CAMB Transfer_tot, axionHMcode Eq. 9 weights, the boost
  denominator) and how that maps to the Weyl potential sourcing the lensing. That is a
  physics extension, not plumbing — DEFERRED until Dan's specifics; the current build
  keeps the hard error for (a) as a pragmatic guard, explicitly documented as
  not-fundamental (see [[axionhmcode-architecture]] "Regime gating").
- (b) stays warn-under-`strict: False` (unchanged, now with Wayne's blessing).
- (c),(d): warn — BUT Wayne flags a real breakdown: redshifts above the KG→EFA switch
  (z > z_osc, i.e. a < a_osc, before the axion oscillates and becomes DM-like). There the
  axion transfer function is not a DM density contrast in the halo-model sense at all.
  Logged as new risk R10 in [[axionhmcode-open-questions]]: the Theory class must check
  the top of its z grid against z_osc = 1/a_osc - 1 (available as
  `results.Params.Axion.a_osc`, camb/axion.py:36). For the target window
  (m ~ 1e-25..1e-23 eV) oscillation begins deep in the radiation era (z_osc >> 1e4),
  far above any lensing-source redshift, so this is an assertion that should never fire
  in production — but it must exist, because a light-mass misconfiguration would
  silently feed garbage otherwise.

## On question 6 (Eq.9-consistent ratio) and question 4 (staging)

> plausible, possibly the only way to fix that if we want a drag and drop compatibility
> for axionHMCode updates in 4.

- Question 6: the Eq.9-consistent ratio proposal (numerator AND denominator in
  axionHMcode's own decomposition) is deemed plausible — adopted as the working
  convention, pending Dan's specifics.
- Question 4 gains a design constraint: DRAG-AND-DROP COMPATIBILITY with future
  axionHMcode updates. The Theory class must call axionHMcode only through its public
  entry points (func_power_spec_dic-equivalent construction, HMCode_param_dic,
  func_axion_param_dic, func_full_halo_model_ax) and never modify axionHMcode files;
  if a group fork is created it is pin-only (no code changes). The Eq.9-internal ratio
  definition supports this: all conventions live on our side of the interface.

## Pending

- Dan to provide more specifics (expected: the smooth-component conventions for
  DE-like/large-fax axions and their Weyl-potential mapping; possibly ratio-convention
  details for question 6).
