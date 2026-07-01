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

- UPDATE (2026-07-01, later same day): interim relative symlinks now in place —
  `external_modules/code/AxiECAMB -> ../../../../New_AxiECAMB` and
  `external_modules/code/axionHMcode -> ../../../../axionHMcode`
  ([[axionhmcode-architecture]] "Staging into Cocoa"). Permanent mechanism still open.

## axionHMcode

See [[axionhmcode-api]]. Cloned from SophieMLV/axionHMcode; includes the Dome et al. update
(version 'basic' | 'dome'), numba (@njit in halo_bias.py, variance.py), Dentler nuisance
parameters alpha_1/alpha_2/gamma_1/gamma_2 via cosmo_dic keys, no packaging files.

## Review pass 2 (2026-07-01) — additional verified facts

- Provider dispatch: `Theory.get_can_provide_methods` (cobaya/theory.py:173) uses
  `get_class_methods(..., start="get_", not_base=Theory)` (tools.py:937-948): EVERY method
  named `get_*` defined on the Theory subclass (with self as first arg) is registered as a
  providable quantity. Implementation rule: name all internal helpers `_something`, never
  `get_something`, or cobaya will treat them as products.
- Nonlinear-lensing z grid: `GetComputedPKRedshifts` (fortran/results.f90:1168) — with
  NonLinear_lens/both + DoLensing, NLL grid = nint(10*5*NL_Boost) redshifts LINEAR in
  [0, maxRedshift], maxRedshift = 10 (15 when NL_Boost >= 2.5),
  NL_Boost = AccuracyBoost*NonlinSourceBoost; also raises Transfer kmax to
  max(kmax, 5*NL_Boost) (internal units — confirm in Phase 0). Master array = union with
  PK_redshifts, exposed to Python ONLY as `results.transfer_redshifts`
  (camb/results.py:226; order not guaranteed ascending — sort before use).
  Params.Transfer.PK_redshifts does NOT include the NLL grid. At AccuracyBoost=1.5
  (group default) expect ~75 z nodes.
- `GetNonLinRatios_All` (which ExternalNonLinearRatio error-stops on) is reachable only via
  the 21cm power path (fortran/results.f90:4153, guarded by NonLinear /= None and /= Lens
  inside the 21cm PK function, single-redshift only). TT/TE/EE/phiphi and Pk_grid never hit
  it — risk R7 downgraded to "21cm outputs unsupported".
- Python environment is TWO-LAYERED (PI correction 2026-07-01; cocoa/README.md "(cocoa)(.local)
  is a feature, not a bug"): conda env cocoapy310 (yml files; includes numba 0.60) PLUS a
  Cocoa-private pip prefix at `${ROOTDIR}/.local` installed by
  `installation_scripts/setup_pip_core_packages.sh` and activated by `start_cocoa.sh`.
  The .local layer SHADOWS conda packages — notably `numpy==1.26.3` (or 1.23.5 under
  COCOA_FORCE_NUMPY_1_23), plus mpi4py 4.0.3, setuptools, emcee, sacc, jax 0.4.18, etc.
  numba 0.60 (conda layer) supports numpy 1.22-2.0, so the overlay is compatible on paper,
  but the Phase-0 smoke test (`import numba` + a trivial @njit call) must run inside the
  ACTIVE (cocoa)(.local) environment, not the bare conda env. If axionHMcode ever needs a
  new pip package or a re-pin, the sanctioned mechanism is adding it to the arrays in
  setup_pip_core_packages.sh (pip --prefix ${ROOTDIR}/.local + sentinel-hash caching) —
  never an ad-hoc global pip install. Reloading (.local) requires re-sourcing
  start_cocoa.sh after any set_installation_options.sh edit.
- InitPower attribute names on the ctypes class (camb/initialpower.py:128-136): `ns`,
  `pivot_scalar` (Mpc^-1), `As` — read off results.Params.InitPower after the cobaya
  wrapper sets them (camb.py:714, before the ratio callback at :717).
- axionHMcode notebook ground truth (example_file.ipynb): mass grid
  `M_arr = np.logspace(cosmo_dic['M_min'], cosmo_dic['M_max'], 100)` — M_min/M_max are
  LOG10 EXPONENTS (7, 18 in input_file.txt); the linear defaults in load_cosmology.py
  (1e8/1e17) are inconsistent with this usage — an upstream quirk; our cosmo_dic must use
  the exponent convention. `func_full_halo_model_ax` returns a tuple; element [0] is the
  total nonlinear P(k). One power_spec_dic instance is passed everywhere (no separate
  sigma dict in the notebook flow).
- Port doc gap (candidate Phase-5 fix): the MatterTransferData class docstring in
  camb/results.py lists transfer indices only through 13 — Transfer_axion = 14 is missing
  from the docs (the code itself is correct).
