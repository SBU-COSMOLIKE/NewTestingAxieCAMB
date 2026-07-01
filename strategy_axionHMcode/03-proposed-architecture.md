---
name: axionhmcode-architecture
description: Design of the AxionHMcodeBoost Theory class, data flow, ratio definition, and Cocoa staging
metadata:
  type: project
---

# Proposed architecture

## One new Theory class, zero Cobaya patches

A single Python file defining (working name) `class AxionHMcodeBoost(Theory)`, referenced
from the yaml with an explicit class path — no file needs to live inside the cobaya tree:

```yaml
theory:
  camb:
    path: ./external_modules/code/AxiECAMB
    use_non_linear_ratio: True        # top-level camb option, NOT inside extra_args
    extra_args:
      nonlinear: NonLinear_both
      halofit_version: mead2020       # value irrelevant once use_non_linear_ratio is on;
                                      # ExternalNonLinearRatio is forced (camb.py:331-334)
      ...
  axionhmcode_boost.AxionHMcodeBoost:
    python_path: ./external_modules/code/AxiECAMB   # or wherever the file lives
    version: dome            # 'basic' | 'dome'
    strict: False            # PI directive 2026-07-01 (see Regime gating): False =
                             # warn + extrapolate outside calibration; True = hard error
    ...
```

Note the prompt's yaml example places `use_non_linear_ratio: True` inside `extra_args` —
that is the wrong level; it is a class attribute of the camb theory block
(cobaya `theories/camb/camb.yaml:26`; the upstream test passes it top-level).

## Class skeleton (contracts only — no code yet)

- Language target: Python 3.10 — the Cocoa environment python (cocoapy310*.yml; PI
  directive 2026-07-01). Write 3.10-compatible code only; no 3.11+ features.
- `initialize()`: append the axionHMcode checkout to `sys.path`, import its modules once,
  read class options (version basic/dome, mass-grid size, the `strict` validity flag,
  optional fixed nuisance values alpha_1, alpha_2, gamma_1, gamma_2).
- `get_requirements()`: `{"CAMB_transfers": None}`. Nothing else — every cosmological number
  is read from the results object at call time, so the class needs no sampled params of its
  own (unless we later expose the Dentler nuisance parameters as sampled — see
  [[axionhmcode-open-questions]]).
- `get_non_linear_ratio(self, results)` — the whole physics pipeline:
  1. Cosmology from `results.Params`: h = H0/100, ombh2, omch2, As/ns/pivot from
     `InitPower` (power-law assumed — document), axion m_ax/omaxh2 and regime flag
     `is_de_like` from `results.Params.Axion` (`New_AxiECAMB/camb/axion.py:33`;
     read from results.Params, the state copy — [[axiecamb-port-project]]).
  2. `tdata = results.get_matter_transfer_data()`; slice `delta_cdm`, `delta_baryon`,
     `delta_axion` (Transfer_axion = 14, `camb/model.py:56`), `delta_tot` at each z.
     CAMB's matter transfer convention (T/k^2, k in h/Mpc) matches what axionHMcode's
     `load_transfer_from_file` reads from axionCAMB output files — verify numerically
     once (validation V1) before trusting.
  3. Build `power_spec_dic` per z by reusing axionHMcode's own `transfer_to_PS` and the
     cold-combination formula from `func_power_spec_dic`
     (`axionHMcode/axionCAMB_and_lin_PS/lin_power_spectrum.py:15,51`) — reuse their code,
     do not re-derive (numerics-verbatim discipline).
  4. Loop over the z grid (axionHMcode computes ONE redshift per call: `cosmo_dic['z']`
     is scalar): build cosmo_dic (including their internal growth `G_a` exactly as
     `load_cosmology.py:205` does), `HMCode_param_dic`, `func_axion_param_dic`, then
     `func_full_halo_model_ax(...)` → P_NL(k) at that z.
  5. ratio_z = sqrt(P_NL / P_L) with numerator and denominator in the SAME decomposition
     (axionHMcode Eq. 9 weights) so the neutrino convention cancels — see
     [[axionhmcode-open-questions]] item 6.
  6. Return `{"k_h": k, "z": zgrid, "ratio": ratio}` (ratio shape (nz, nk), z ascending).

## Grid policy

- k: the transfer k grid (or a log subset); CAMB clamps outside the grid
  (ExternalNonLinearRatio.f90:81-85), so cover [k_min_transfer, kmax_yaml≈10 h/Mpc].
- z: ACCURACY-FIRST POLICY (PI, 2026-07-01 — performance is explicitly not a constraint):
  evaluate axionHMcode at every redshift CAMB actually uses, i.e. the full
  `results.Params.Transfer.PK_redshifts[:PK_num_redshifts]` grid (the trivial-test pattern,
  test_cosmo_multi_theory.py:287-292), which under NonLinear_both is augmented internally
  for the lensing sources — log its actual content at runtime in Phase 0 to confirm
  coverage. No coarse-grid/interpolation shortcut unless a convergence test (V6) proves it
  changes lensed C_ell and C_L^phiphi below the numerical targets; cost alone never
  justifies thinning. If per-step cost ever becomes prohibitive for production MCMC, the
  escape hatch is an ML emulator trained on THIS pipeline's output (same
  get_non_linear_ratio interface), as a later phase.
  High-z boundary must decay toward 1 (or be padded with 1) because clamping reuses the last
  z slice at all higher redshifts.

## Regime gating (PI decision 2026-07-01: the `strict` yaml flag)

One boolean class option, `strict`, written in the theory-block yaml section (see the
snippet above). Semantics over the out-of-validity conditions
(a) DE-like axion, (b) fax beyond basic/dome calibration, (c) z outside 1-8 for dome,
(d) m far from 1e-24.5 eV for dome:

- `strict: True` — any of (a)-(d) is a hard error (stop_at_error semantics): the run halts
  with a message naming the violated condition and its bound.
- `strict: False` (default) — (b)-(d) log a warning and extrapolate, matching how
  Gaughan et al. deliberately ran extrapolations; the warning states the condition,
  the sampled value, and the calibration bound.
- Condition (a) (`results.Params.Axion.is_de_like == True`, m/H0 < 10) is a hard error
  REGARDLESS of the flag in the current build. Per collaborator guidance (2026-07-01,
  [[axionhmcode-collaborator-guidance]]) this is pragmatic, not fundamental: a
  smooth-vs-clustered decomposition (DE-like axion as w0wa-like smooth DE, or as a smooth
  mixed-DM component) could lift (a) and (b) once the "matter" conventions and their
  Weyl-potential mapping are defined — deferred pending further collaborator input. Until
  then, use the halofit-original path (existing EXAMPLE yamls) for DE-like masses.
- Additional unconditional check (R10, from collaborator guidance): hard-error if the top of the z grid
  reaches z_osc = 1/a_osc - 1 (`results.Params.Axion.a_osc`) — above the KG→EFA switch the
  axion is not yet DM-like and its transfer function is not a halo-model density contrast.
  Never fires in the target mass window (z_osc >> 1e4); guards light-mass misconfigurations.

## Staging into Cocoa

PI decision 2026-07-01: for now, symlinks — get the pipeline working BEFORE touching any
Cocoa script. Created (relative links, verified resolving):

    cocoa/Cocoa/external_modules/code/AxiECAMB    -> ../../../../New_AxiECAMB
    cocoa/Cocoa/external_modules/code/axionHMcode -> ../../../../axionHMcode

Notes on the interim setup:
- The link is named `AxiECAMB` (matching the yaml `path: ./external_modules/code/AxiECAMB`)
  but points at `New_AxiECAMB` (the CAMB-1.6.7 port) — NOT at `rayne/AxiECAMB/` (the
  original Nov13 fork). Do not confuse the two.
- Re-running Cocoa setup scripts with OVERWRITE flags could in principle disturb
  external_modules content; the links are cheap to recreate (two ln -s commands above).
- Permanent staging (Cocoa-style clone scripts with pinned commits; group fork of
  axionHMcode recommended) is deferred until the pipeline works end-to-end —
  [[axionhmcode-open-questions]] item 4 stays open for the permanent mechanism.
- `axionHMcode` has no setup.py/pyproject (plain package folder with `__init__.py`) — the
  Theory class inserts its parent dir into sys.path. It uses numba — env dependency to
  check in Phase 0.
- PI DIRECTIVE (2026-07-01): the Theory block, once built, is added as its OWN FOLDER inside
  the New_AxiECAMB repo (committable, like strategy_axionHMcode/) — e.g.
  `New_AxiECAMB/axionhmcode_theory/` (final name TBD) holding the Theory-class module plus
  its yaml defaults/examples. The yaml then references it via
  `python_path: ./external_modules/code/AxiECAMB/<folder>` (works through the symlink).
  Do NOT scatter it as loose files in the repo root.

## Fallback (only if the zero-patch route fails)

If yaml `python_path` class loading proves incompatible with some Cocoa wrapper detail, the
fallback is ONE new-file "patch": add the theory module under
`cobaya_changes/cobaya/theories/<name>/` and one `cppatch`/`cppatchfolder` block APPENDED at
the end of `setup_cobaya.sh` (after all existing blocks, per the ordering constraint). The
existing patches are never touched either way.
