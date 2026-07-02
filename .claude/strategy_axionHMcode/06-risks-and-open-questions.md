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
4. **Staging decisions.** FULLY ANSWERED. (a)+(b) permanent mechanism (PI, 2026-07-02):
   Cocoa scripts `installation_scripts/setup_axie_camb.sh` (clones the AxiECAMB repo,
   pinned via AXIE_CAMB_GIT_COMMIT, into `external_modules/code/axiecamb` — CANONICAL
   NAME, lowercase — applying Makefile/compiler patches from
   `cocoa_installation_libraries/axiecamb_changes/`; also clones UPSTREAM
   SophieMLV/axionHMcode, pinned via AXION_HMCODE_GIT_COMMIT, into
   `external_modules/code/axionHMcode`) and `compile_axie_camb.sh` (setup.py build with
   the RECOMBINATION_FILES variants), gated by INSTALL_AXIE_CAMB_V2 in
   set_installation_options.sh. axionHMcode needs NO compile step: pure Python, no
   packaging, deps already in the Cocoa layers; the only compilation is numba JIT once
   per process (~7 s; upstream sets no cache=True, and drag-and-drop forbids adding it).
   Part (c) (PI directive, 2026-07-01): the Theory block lives as its own folder inside
   the New_AxiECAMB repo. CONSTRAINT (collaborator guidance, 2026-07-01): drag-and-drop
   compatibility with axionHMcode updates — public entry points only, never modify its
   files. Interim mac symlinks: the `AxiECAMB` symlink collides case-insensitively with
   an `axiecamb` clone on macOS — remove it before running the setup script there.
5. **Performance budget.** ANSWERED (PI, 2026-07-01): accuracy is the worry; performance is
   not a constraint. Consequences: (a) z-grid policy = dense — evaluate axionHMcode at every
   redshift CAMB actually uses (the full nonlinear-lensing-augmented
   `results.transfer_redshifts` grid — see R4);
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
- R4 — RESOLVED by review pass 2 (2026-07-01, [[axionhmcode-verified-facts]]): the
  nonlinear-lensing grid is built by GetComputedPKRedshifts (fortran/results.f90:1168) —
  nint(50*NL_Boost) redshifts linear in [0, 10] ([0, 15] if NL_Boost >= 2.5), merged with
  the user PK_redshifts into `results.transfer_redshifts` (the ONLY Python-visible union;
  Params.Transfer.PK_redshifts does NOT include it — the trivial-test pattern would
  silently miss the lensing grid). Theory class reads np.array(results.transfer_redshifts),
  sorted ascending. High-z clamping reuses the last z slice — pad so the boundary
  slice is ≈ 1. Phase-0 runtime log remains as confirmation.
- R5 — power-law primordial assumption: axionHMcode's primordial_PS is pure power law
  (no running). Fine for the planned chains; document as a limitation.
- R6 — axionHMcode internal growth G_a and Omega_w_0 = 1 - Omega_m_0 assume flat LCDM
  background without the axion's early-DE phase. This is part of the calibrated model —
  keep verbatim, note in the README report (same "three Omega_m's" caveat family as the
  halofit→HMcode analysis in the prompt's Previous Attempt Notes).
- R7 — DOWNGRADED by review pass 2: GetNonLinRatios_All (which ExternalNonLinearRatio
  error-stops on) is reachable only through the 21cm power path (fortran/results.f90:4153;
  guarded, single-redshift). TT/TE/EE/phiphi and Pk_grid never call it. Document as
  "21cm outputs unsupported with the external ratio".
- R8 — global/module state and MPI: axionHMcode is functional (dict in, array out) — no
  module-level caches spotted so far; confirm before declaring MPI-safe. numba part:
  available from the conda layer (0.60), but Cocoa's .local pip layer overlays
  numpy==1.26.3 on top of it (PI correction 2026-07-01 — two-layer environment,
  [[axionhmcode-verified-facts]]); compatible on paper, smoke-test inside the active
  (cocoa)(.local) env in Phase 0. New pip deps, if ever needed, go into
  setup_pip_core_packages.sh (--prefix ${ROOTDIR}/.local). JIT cost paid once per process.
- R11 — cobaya treats EVERY `get_*` method on a Theory subclass as a providable product
  (theory.py:173 → tools.py:937-948). The boost class must expose exactly one such method
  (get_non_linear_ratio); all helpers get underscore-prefixed names, or cobaya's
  dependency resolver will see phantom products.
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
