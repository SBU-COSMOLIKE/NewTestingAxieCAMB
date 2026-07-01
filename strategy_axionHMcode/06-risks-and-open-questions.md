---
name: axionhmcode-open-questions
description: PI decisions needed (version default, out-of-validity policy, mass priors, staging, performance budget, neutrino ratio convention) and known risks
metadata:
  type: project
---

# Risks and open questions

## Questions for the PI (blocking design decisions)

1. **basic vs dome default, and nuisance parameters.** ANSWERED as working default
   (collaborator guidance 2026-07-01, [[axionhmcode-collaborator-guidance]] — note the
   epistemic caveat there): dome is the default (most recent recalibration); the
   Dentler nuisance parameters (alpha_1, alpha_2, gamma_1, gamma_2) are exposed as
   Theory-block parameters so they can be sampled (defaults = calibration values; example
   MCMC shows the Gaughan-et-al prior blocks). Verify against code/papers in Phases 0-1.
2. **Out-of-validity policy.** ANSWERED (PI, 2026-07-01): a single boolean `strict` flag
   written in the theory-block yaml section. `strict: True` = hard error on any
   out-of-calibration condition; `strict: False` (default) = warn + extrapolate for
   (b) fax beyond 0.5/0.3, (c) z outside 1-8 for dome, (d) m far from 1e-24.5 eV for dome.
   Condition (a), DE-like axion (is_de_like True), is a hard error regardless of the flag
   in the CURRENT build — per collaborator guidance (2026-07-01,
   [[axionhmcode-collaborator-guidance]]) this restriction is pragmatic, not fundamental:
   a smooth-vs-clustered decomposition (DE-like axion treated like w0wa, or as a smooth
   mixed-DM component) could lift (a) and (b) once the "matter" conventions and their
   Weyl-potential mapping are pinned down — a physics extension deferred pending further
   collaborator input. Full semantics: [[axionhmcode-architecture]] "Regime gating".
3. **Mass prior for the new MCMC yaml.** ANSWERED (PI, 2026-07-01): target window
   m ~ 1e-25..1e-23 eV. Document BOTH modes in a README: the EVALUATE example carries the
   log-mass prior block (logmx uniform in [-25,-23]) with the mass pinned in the evaluate
   sampler's override (single cosmology); the MCMC example fixes the mass, with yaml
   comments showing how to switch to the uniform log-mass prior. Full README-ready recipes:
   [[axionhmcode-mass-prior-recipes]] (file 09).
4. **Staging decisions.** INTERIM ANSWER (PI, 2026-07-01): symlinks from the current
   positions into external_modules/code (created and verified — see
   [[axionhmcode-architecture]] "Staging into Cocoa"); no Cocoa scripts touched until the
   pipeline works. STILL OPEN for the permanent mechanism: (a) clone script with pinned
   commit for New_AxiECAMB; (b) axionHMcode vendored vs cloned, upstream vs group fork
   (fork recommended — pinning + we may need to touch numba/global-state details).
   Part (c) is ANSWERED (PI directive, 2026-07-01): the Theory block lives as its own
   folder inside the New_AxiECAMB repo, committable alongside strategy_axionHMcode/
   (see [[axionhmcode-architecture]] "Staging into Cocoa"). NEW CONSTRAINT (collaborator
   guidance, 2026-07-01): whatever the permanent mechanism, preserve DRAG-AND-DROP
   compatibility with future axionHMcode updates — call it only through public entry
   points, never modify its files; any group fork is pin-only.
5. **Performance budget.** ANSWERED (PI, 2026-07-01): accuracy is the worry; performance is
   not a constraint. Consequences: (a) z-grid policy = dense — evaluate axionHMcode at every
   redshift CAMB actually uses (the full nonlinear-lensing-augmented PK_redshifts grid);
   (b) any grid thinning or interpolation shortcut must be justified by an explicit
   convergence test (validation V6), never by cost; (c) if the pipeline later proves too
   costly for production MCMC, the PI will use THIS code as the training-data generator for
   ML emulators (emulmps-style, dropping into the same get_non_linear_ratio interface) —
   an emulator is a later phase, not part of this build.
6. **Ratio convention with massive neutrinos.** WORKING CONVENTION ADOPTED (assessed
   plausible in collaborator guidance, 2026-07-01, [[axionhmcode-collaborator-guidance]];
   final details pending further collaborator input): define the boost with numerator AND
   denominator in axionHMcode's own Eq. 9
   decomposition (P_NL^Eq9 / P_L^Eq9), so the convention cancels at linear order, and CAMB
   multiplies its own P_L^tot(incl. nu). This also preserves drag-and-drop compatibility
   with future axionHMcode updates (all conventions live on our side of the interface). The
   residual is how nu suppression enters the halo model — same limitation Gaughan et al.
   accepted. Confirm.

## Known risks / traps (mitigations in [[axionhmcode-architecture]])

- R1 — linear-P(k) extraction: get_linear_matter_power_spectrum on transfers-only results
  triggers calc_power_spectra → error stop (ratio unset). Use get_matter_transfer_data only.
- R2 — transfer normalization mismatch between CAMB's MatterTransferData and axionCAMB's
  file convention would silently rescale the boost. Killed by validation V1 (numeric check),
  not by reading docs.
- R3 — which transfer variables include the axion in each regime (Transfer_tot kluge:
  axion in tot iff DM-like; Transfer_nonu semantics for the axion unverified). Must be
  established in Phase 1 with a direct dump per regime.
- R4 — z-grid content under NonLinear_both: PK_redshifts is augmented internally by CAMB for
  nonlinear lensing; actual node count/range must be logged at runtime (drives cost).
  High-z clamping reuses the last z slice — pad the grid so the boundary slice is ≈ 1.
- R5 — power-law primordial assumption: axionHMcode's primordial_PS is pure power law
  (no running). Fine for the planned chains; document as a limitation.
- R6 — axionHMcode internal growth G_a and Omega_w_0 = 1 - Omega_m_0 assume flat LCDM
  background without the axion's early-DE phase. This is part of the calibrated model —
  keep verbatim, note in the README report (same "three Omega_m's" caveat family as the
  halofit→HMcode analysis in the prompt's Previous Attempt Notes).
- R7 — GetNonLinRatios_All hard-stops: any requirement that requests nonlinear velocity
  spectra (var pairs with velocities under NonLinear_pk) will kill the run. Not needed for
  TT/TE/EE/phiphi; document.
- R8 — global/module state and MPI: axionHMcode is functional (dict in, array out) — no
  module-level caches spotted so far; confirm before declaring MPI-safe. numba JIT cost is
  paid once per process.
- R9 — RESOLVED (was wrong; PI correction 2026-07-01): Cocoa DOES pin the cobaya commit —
  `export COBAYA_GIT_COMMIT="899f30a49f85de610dac321e91a1af50018e56aa"` at
  `cocoa/Cocoa/set_installation_options.sh:220` (consumed by setup_cobaya.sh:122-124;
  unset in flags_impl_unset_keys.sh:73). That is exactly the commit verified to contain the
  full PR#480 machinery, so the mechanism is stable under reinstalls. Residual note: if the
  group ever bumps COBAYA_GIT_COMMIT, re-verify the use_non_linear_ratio code paths
  ([[axionhmcode-verified-facts]]) as part of the upgrade.
- R10 — redshifts above the KG→EFA switch (collaborator guidance, 2026-07-01,
  [[axionhmcode-collaborator-guidance]]): for z > z_osc (a < a_osc) the axion has not yet
  begun oscillating and is not DM-like — its transfer function is not a halo-model density
  contrast, so feeding those redshifts to axionHMcode is meaningless, not merely
  uncalibrated. The Theory class must check the top of its z grid against
  z_osc = 1/a_osc - 1 (`results.Params.Axion.a_osc`, camb/axion.py:36) and hard-error if
  the grid reaches it (regardless of `strict`). For the target window
  (m ~ 1e-25..1e-23 eV) z_osc >> 1e4, far above any lensing-source redshift, so the
  assertion should never fire in production — it guards against light-mass
  misconfigurations silently producing garbage.
