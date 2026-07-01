---
name: axionhmcode-verified-facts
description: Code facts verified 2026-07-01 with file:line citations (cobaya PR#480 present, CAMB-side classes, transfer variables, traps)
metadata:
  type: project
---

# Verified facts (2026-07-01)

Everything below was read from the working tree, not assumed. Re-verify line numbers after
any git pull.

## Cobaya (cocoa/Cocoa/cobaya — HEAD 899f30a4, PINNED by Cocoa)

The checkout is not a floating master: Cocoa pins it via
`export COBAYA_GIT_COMMIT="899f30a49f85de610dac321e91a1af50018e56aa"`
(`cocoa/Cocoa/set_installation_options.sh:220`, consumed in setup_cobaya.sh:122-124).
All facts below are therefore stable under Cocoa reinstalls.

- `use_non_linear_ratio: bool` class attribute: `cobaya/theories/camb/camb.py:263`; yaml
  default False at `cobaya/theories/camb/camb.yaml:26`. It is a TOP-LEVEL camb option, not
  an extra_args entry.
- initialize() forces `extra_args["non_linear_model"] = camb.nonlinear.ExternalNonLinearRatio`
  when the flag is on: camb.py:331-334.
- must_provide adds `non_linear_ratio` requirement iff flag and needs_perts: camb.py:640-641.
- calculate() order (camb.py:676-731): get_CAMB_transfers → set InitPower on results (:714)
  → `provider.get_non_linear_ratio(results)` (:717) → `results.Params.NonLinearModel.set_ratio`
  (:718-722) → `results.power_spectra_from_transfer()` (:731). Ratio provider is called with
  the live results object — this kills the circular dependency.
- Reference provider+test: `TrivialNonLinearRatio` (requires only "CAMB_transfers") and
  `test_trivial_non_linear_ratio`: `tests/test_cosmo_multi_theory.py:279-338`. Ratio dict
  format: `{"k_h": 1d, "z": 1d ascending, "ratio": 2d (nz, nk)}`; ratio = sqrt(P_NL/P_L)
  (test asserts pk_nonlin = ratio^2 * pk_lin at :306).
- `CambTransfers` helper (camb.py:1162): with non_linear_sources uses
  `get_transfer_functions(camb_params, only_time_sources=True)` (:1219-1221); else full
  transfers or background only (:1223-1227).
- Cocoa patch inventory: `cocoa/cocoa_installation_libraries/cobaya_changes/` — mostly
  full-file copies (`cppatch`) plus three real patch files (planck_clik.patch applied with -R,
  model.patch, InstallableLikelihood.patch), applied in `setup_cobaya.sh`; then
  `pip install --editable`. Cocoa's replacement camb.yaml keeps `use_non_linear_ratio` →
  the mechanism survives Cocoa's patching. New patches would be appended after existing
  blocks in setup_cobaya.sh.

## New_AxiECAMB (CAMB 1.6.7 + ULA port, [[axiecamb-port-project]])

- Fortran: `fortran/ExternalNonLinearRatio.f90` — `TExternalNonLinearRatio` extends
  TNonLinearModel; TInterpGrid2D over (k/h, z); clamps k and z to grid bounds (:81-85);
  initializes nonlin_ratio to 1 (:78); `error stop` if ratio not set (:69-70);
  GetNonLinRatios_All (velocity corrections) unsupported → `NonLinear_pk` var pairs with
  velocities would hard-stop (:96).
- Python: `camb/nonlinear.py:322` `ExternalNonLinearRatio` with `set_ratio(k_h, z, ratio)`
  (:355-372; validates ratio.shape == (len(z), len(k_h))) and `clear_ratio`.
- Transfer variables: `Transfer_axion = 14` == `delta_axion` (`camb/model.py:56,100`;
  `fortran/results.f90:3282,3287`). `Transfer_nonu = 8` exists; whether the axion is inside
  `Transfer_tot` follows the Hlozek clustering kluge (axion in tot iff DM-like) — re-verify
  which variables include the axion in each regime before building spectra (Phase 1 check).
- Axion state on results: `results.Params.Axion` exposes `is_de_like` (m/H0 < 10), `a_osc`,
  cycle-averaged quantities (`camb/axion.py:33-52`); set_params signature
  m_ax/omaxh2/omdah2/axfrac/dfac (:78).
- TRAP: `get_linear_matter_power_spectrum` calls `calc_power_spectra` when the results are
  transfers-only (`camb/results.py:835-836`) → would invoke the unset ExternalNonLinearRatio
  → Fortran error stop. Use `get_matter_transfer_data()` (results.py:768) instead.
- `only_time_sources=True` semantics (`camb/camb.py:58-72`): skips CMB l,k transfer functions
  and nonlinear scaling only; matter transfers remain available; results reusable with
  different initial power AND consistent nonlinear lensed spectra — exactly our use case.
- `EXAMPLE_EVALUATE1.yaml` / `EXAMPLE_MCMC1.yaml` exist and run pristine Cobaya with no
  wrapper patch (axion params discovered via camb.get_valid_numerical_params); they sample
  logmx in [-34,-31] (DE-like) with halofit_version: original.

## Cocoa staging status

- `cocoa/Cocoa/external_modules/code/` contains CAMB, FAST-PT, emulators, etc. — but NO
  AxiECAMB/New_AxiECAMB and NO axionHMcode yet (verified 2026-07-01). Staging is part of
  the work.

## axionHMcode

See [[axionhmcode-api]]. Cloned from SophieMLV/axionHMcode; includes the Dome et al. update
(version 'basic' | 'dome'), numba (@njit in halo_bias.py, variance.py), Dentler nuisance
parameters alpha_1/alpha_2/gamma_1/gamma_2 via cosmo_dic keys, no packaging files.
