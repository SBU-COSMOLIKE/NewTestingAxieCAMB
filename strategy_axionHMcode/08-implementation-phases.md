---
name: axionhmcode-implementation-phases
description: Phased implementation plan with checkpoints; Phase 0 benchmarks gate the design decisions
metadata:
  type: project
---

# Implementation phases (no code written as of 2026-07-01)

Each phase ends at a checkpoint that must pass before the next phase starts (porting-skill
discipline: compile checkpoints + null tests every time).

## Phase 0 — baselines and benchmarks (no new code)

- Two-track environment policy (PI, 2026-07-01): for TESTING (Phases 0-2) a plain
  `pip install "cobaya>=3.6.2"` environment suffices — no Cocoa needed. 3.6.2 is the
  version at Cocoa's pinned commit (899f30a4 = the 3.6.2 version bump, 2026-03-27), and
  the last camb.py change before it is PR#480 itself (commit 975a9413), so pip cobaya
  >= 3.6.2 has the identical use_non_linear_ratio mechanism (the port was already
  validated this way — EXAMPLE_EVALUATE1.yaml header, cobayapristine env). Cocoa's
  (cocoa)(.local) environment becomes mandatory at Phase 3 (staging, Cocoa likelihoods,
  end-to-end V5). Before Phase 3: run the environment smoke tests inside (cocoa)(.local)
  — python 3.10, `import numba` + trivial @njit against the overlaid numpy 1.26.3
  ([[axionhmcode-verified-facts]] "two-layered"). Write all code 3.10-compatible from the
  start regardless of the testing env.
- Build New_AxiECAMB main branch (serial make; miniforge gfortran 14.3; SDKROOT quirks —
  [[axiecamb-port-project]]).
- Run the existing python test suite + EXAMPLE_EVALUATE1.yaml to confirm the starting state.
- Run cobaya's test_trivial_non_linear_ratio against stock CAMB, then V0a against
  New_AxiECAMB.
- Benchmark axionHMcode: single-z, single-cosmology wall time (basic and dome; with numba
  warm), on the notebook's own example inputs. Informational only — question 5 is answered
  (accuracy first, dense z-grid; performance never drives design). The number sets MCMC
  wall-time expectations and tells us when the future emulator phase becomes attractive.
- Log what `results.transfer_redshifts` actually contains under NonLinear_both with lensed
  Cls requested — confirming the GetComputedPKRedshifts analysis (R4, resolved on paper:
  ~50*NL_Boost nodes linear in [0,10]; [[axionhmcode-verified-facts]] review pass 2).
- Checkpoint: V0a passes; benchmark numbers recorded here.

## Phase 1 — standalone prototype (pyCAMB + axionHMcode, no Cobaya)

- One script: AxiECAMB run → get_matter_transfer_data → power_spec_dic → axionHMcode →
  ratio grid → set on ExternalNonLinearRatio → calc_power_spectra → lensed spectra.
  This is the prompt's "if we were calling pyCAMB directly" snippet made real (note the
  actual API is params.NonLinearModel = ExternalNonLinearRatio(); .set_ratio(...);
  there is no results.set_nonlin_ratio).
- Establish transfer-variable semantics per regime (V1, risks R2/R3), including what
  Transfer_axion returns for z > z_osc (KG phase) — R10 lead in
  [[axionhmcode-collaborator-guidance]].
- Check how the public AxiCAMB (github.com/adammoss/AxiCAMB) handles redshifts near/above
  the EFA switch in its axionHMcode interface (R10 lead).
- Reproduce Gaughan figures (V3) from this script.
- Checkpoint: V1 + V3 pass; conventions frozen and written into
  [[axionhmcode-architecture]].

## Phase 2 — the Theory class

- Write AxionHMcodeBoost(Theory) per [[axionhmcode-architecture]]: emulmps style (2-space,
  concise), Python 3.10 (the Cocoa environment python — no 3.11+ features),
  get_requirements = CAMB_transfers, get_non_linear_ratio(results) wrapping the
  Phase-1 pipeline, class options (version, grids, the `strict` validity flag,
  nuisance values).
- PI directive (2026-07-01): deliver it as its own folder inside the New_AxiECAMB repo
  (committable, sibling of strategy_axionHMcode/) — module + yaml defaults/examples in one
  place, referenced from yamls via python_path through the external_modules symlink.
- Unit-level test mirroring TrivialNonLinearRatio but with the real class on a fixed
  cosmology; compare against the Phase-1 script outputs (must be identical — same code).
- Checkpoint: cobaya get_model + loglikes runs on a laptop-scale likelihood set.

## Phase 3 — Cocoa staging + yamls

- Stage New_AxiECAMB at external_modules/code/AxiECAMB and axionHMcode per the PI's
  answer to open question 4; write the EVALUATE and MCMC yamls (mass window
  1e-25..1e-23 eV, use_non_linear_ratio top-level, theory block for the boost class):
  EVALUATE = log-mass prior + evaluate-override pin; MCMC = fixed mass + commented
  log-mass alternative; README from the recipes in [[axionhmcode-mass-prior-recipes]].
- If (and only if) yaml python_path loading fails inside Cocoa: fall back to the single
  appended-patch route in setup_cobaya.sh (03, "Fallback").
- Checkpoint: V5 evaluate passes inside Cocoa environment (start_cocoa.sh sourced).

## Phase 4 — validation battery

- Full V0-V5 matrix from [[axionhmcode-validation-plan]]; record all numbers in a
  validation report file next to these strategy notes.
- Checkpoint: acceptance targets met or deviations explained and signed off by PI.

## Phase 5 — documentation and review

- Append the project report to New_AxiECAMB/README.rst (what was built, conventions chosen,
  validation numbers, limitations R5/R6, out-of-validity policy).
- Critical-review list for the PI: every convention we could not verify against a reference,
  every extrapolation region, accuracy knobs that do not scale.
- Update these strategy files to their as-built state; update auto-memory.
