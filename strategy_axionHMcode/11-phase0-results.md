---
name: axionhmcode-phase0-results
description: Phase 0 measured results (2026-07-01) — V0a machine-precision pass, R4 grid confirmed at runtime, transfer-convention V1 pass, axionHMcode timings, dependency findings
metadata:
  type: project
---

# Phase 0 results (measured 2026-07-01)

Test environment: `rayne/cobaya_test_env` (venv, python 3.12.11 from miniforge) —
cobaya 3.6.2, numpy 1.26.4, scipy 1.13.1, numba 0.60.0, astropy 6.1.7. Versions chosen to
mirror the Cocoa layers (numpy/numba identical pins). Local camb = New_AxiECAMB (checked
out on branch `hmcode`, camblib.so prebuilt; ExternalNonLinearRatio verified present on
`main` too). Scripts: scratchpad phase0_v0a.py, phase0_axionhmcode_bench.py (to be
migrated into the Theory-block folder's tests in Phase 2).

## V0a — trivial-ratio null test: PASS at machine precision

cobaya get_model + loglikes with `use_non_linear_ratio: True`, camb path = New_AxiECAMB,
TrivialNonLinearRatio ratio = 2: Pk_grid(nonlinear) = 4 x Pk_grid(linear) with max
relative deviation 1.1e-15. The full chain (wrapper -> provider callback -> Fortran
TExternalNonLinearRatio -> power_spectra_from_transfer) works against the port.

## R4 — lensing z-grid: CONFIRMED at runtime

Lensed-Cl run (NonLinear_both, lens_potential_accuracy 1, AccuracyBoost 1):
`results.transfer_redshifts` = 50 nodes, linear, [0.0, 9.8] step 0.2 — exactly the
GetComputedPKRedshifts prediction. `Params.Transfer.PK_redshifts` = [0.] ONLY (the
trivial-test pattern would have clamped a z=0 boost slice across the whole lensing
kernel). Matter transfer data available on the transfers-only results in the lensing
path: shape (14, 217, 50) — 14 includes Transfer_axion.

## V1 (z=0 slice) — transfer convention: PASS

P_total assembled from get_matter_transfer_data + axionHMcode transfer_to_PS
(T in CAMB convention, k h/Mpc, As/ns/k_piv power law) vs
results.get_linear_matter_power_spectrum: median |dev| 2.1e-8, max 5.7e-8
(float32 transfer-data precision). Convention identity holds; R2 closed at z=0
(multi-z sweep still owed in Phase 1).

## Regime facts (m_ax = 1e-25 eV, fax = 0.1, input_file.txt cosmology)

is_de_like = False; a_osc = 3.83e-5 -> z_osc = 2.6e4 >> z_grid_max ~ 10. R10 assertion
can never fire in the target window, as designed.

## axionHMcode single-z timings (warm, numba JIT paid; M grid 100 pts, nk 222)

- dome:  first call 7.9 s (JIT), warm 2.7 s
- basic: first call 1.3 s (JIT), warm 1.3 s
- Production estimate at nz = 50 (AccuracyBoost 1): ~65 s/eval (basic) to ~135 s/eval
  (dome), single-threaded. nz = 75 at AccuracyBoost 1.5: scale by 1.5.
- The z loop is embarrassingly parallel — parallelizing it changes wall time only, not
  numbers; plan a `processes` option in the Theory class (accuracy-first compliant).
- Boost sanity at z=0 (that cosmology): B(k=0.1) = 1.06-1.07, B(1) = 7.9 (basic) /
  11.3 (dome), B(5) = 55 (basic) / 77 (dome). Large-k boost magnitudes plausible for
  z=0 halo model; low-k limit check (B -> 1 as k -> 1e-3) owed in Phase 1.

## Dependency findings

- axionHMcode needs ASTROPY (halo_model/axion_density_profile.py imports
  astropy.constants) — present in all Cocoa env files; added to the test venv (6.1.7).
- Third-party surface (full): numpy, scipy, numba, astropy. Nothing else.
- FUTURE RISK: axionHMcode has a dead `from scipy import interpolate, misc` import
  (HMcode_params.py:5; misc never used, only referenced in a comment). scipy removed
  scipy.misc in 1.14 — Cocoa's scipy 1.12 and venv 1.13.1 are safe, but any future scipy
  bump breaks the import. Drag-and-drop fix if it ever bites: pre-seed a stub
  sys.modules["scipy.misc"] in OUR wrapper before importing axionHMcode (no upstream edit).
- numpy 2.x breaks the CAMB-1.6.7-era python layer (YHe ctypes assignment,
  model.py:722) — another reason the numpy 1.26 pin matters everywhere this port runs.
