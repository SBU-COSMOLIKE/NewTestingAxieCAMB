---
name: axionhmcode-validation-plan
description: Regime-complete validation matrix for the axionHMcode boost pipeline (plumbing nulls, LCDM limits, Gaughan+26 reproduction, cross-checks, end-to-end)
metadata:
  type: project
---

# Validation plan

Discipline from the porting skill: ratios not absolutes; regime-complete matrix; framework
end-to-end; a null test that is binary. Related: [[axionhmcode-architecture]].

## V0 — plumbing nulls (before any axion physics)

- V0a: run cobaya's own `test_trivial_non_linear_ratio` pattern against New_AxiECAMB
  (LCDM params, ratio ≡ 2): Pk_grid(nonlinear=True) must equal 4 x linear. Proves the whole
  use_non_linear_ratio chain works with OUR CAMB inside the Cocoa cobaya.
- V0b (strong plumbing test): compute the HMcode-2020 ratio externally (standalone pyCAMB run
  of the same LCDM cosmology, sqrt(P_NL/P_L) on a grid), feed it through the external-ratio
  path, and compare lensed TT/EE/TE and C_L^phiphi against the SAME cosmology run with
  internal `halofit_version: mead2020`. Agreement to interpolation accuracy (<~0.1% on C_l)
  isolates grid/interpolation errors with zero new physics.
- V0c: ratio ≡ 1 with NonLinear_both vs NonLinear_none lensed spectra — quantifies the pure
  code-path difference (expect ≈ linear-lensing agreement).

## V1 — transfer-convention check (kills R2/R3)

Standalone pyCAMB (no cobaya): for one DM-like and one DE-like axion cosmology, build
power_spec_dic from `get_matter_transfer_data()` and compare power_total against
`results.get_linear_matter_power_spectrum()` from a full calc_power_spectra run
(legal there). Must match to <0.1% where power is not exponentially suppressed. Also dump
delta_tot vs (weighted cdm+baryon+axion) per regime to pin the Transfer_tot kluge and decide
the correct variables for power_cold/power_total.

## V2 — LCDM / small-fax limit

fax → 1e-4 (m DM-like): axionHMcode P_NL vs CAMB HMcode-2020 P_NL (Dome update claims
HMcode-2020 agreement in LCDM — README items 5-6). Quantify; this is the fax=0 limit-recovery
test in the camb-dev sense. Note Gaughan et al. report the basic/dome fax→0 limit is NOT
exactly their LCDM reference (Delta chi2 ~ -0.5..-0.8) — do not expect perfection; record
the offset.

## V3 — reproduce Gaughan, Green & Moss (2605.12054)

The anchor validation. fax = 0.3, m in {1e-23, 1e-24, 1e-25} eV, z in {0, 2}, Planck-2018
fiducial (their Sec. III): reproduce
- Fig. 1: P_NL^axion/P_NL^LCDM for axionHMcode basic and dome (and naive-HMcode curves from
  our V0b machinery for context);
- Fig. 2: lensed TT/EE percent differences vs LCDM (basic ~ -1% at l~3000 TT; dome
  oscillatory up to ~3%);
- Fig. 3: C_L^phiphi differences (5-15% at L >~ 500 between prescriptions).
Their code base is AxiCAMB (standard-EFA lineage) vs our AxiECAMB (PH EFA), so agreement
target is qualitative-to-~percent on the ratios, not exact — differences should shrink at
masses where the EFAs agree (m >~ 1e-24 eV; cf. 2501.13662).

## V4 — cross-check vs the hmcode branch

The New_AxiECAMB `hmcode` branch (axion-aware internal HMcode, items 1+2 of the earlier
proposal, [[axiecamb-port-project]]) gives an independent Fortran-side prediction. In the
DM-like regime compare its lensed spectra vs the axionHMcode-boost pipeline: expected
disagreement at the few-to-10% level in P (they are different halo models — internal HMcode
has no Jeans-scale halo physics), but background/growth-driven differences should be small.
Any order-of-magnitude surprise flags a plumbing bug.

## V5 — end-to-end in Cocoa

- Evaluate yaml (EXAMPLE_EVALUATE1 + boost theory block, DM-like mass): finite loglike,
  timing log (per-component), stable across repeated evaluations (determinism).
- Micro-MCMC (few hundred steps, MPI): no crashes at prior corners; per-step cost within
  the budget from open question 5; derived params sane.
- Prior-corner sweep: fax near prior edge, m at both ends of the window, extreme omaxh2 —
  the pipeline must fail loudly or extrapolate per the configured policy, never silently.

## V6 — numerical convergence of the boost grid (accuracy-first policy)

Performance is not a constraint (PI, 2026-07-01), so the grids must be demonstrated
converged, not merely affordable. At a representative DM-like point (fax=0.3, m=1e-24 eV)
and one extreme corner (fax near prior edge):
- double nk; double nz (if any thinning below the full PK_redshifts grid is ever proposed,
  this test decides it);
- refine the halo-mass integration grid M and widen [M_min, M_max] one decade each way;
- tighten axionHMcode's internal root-finding/integration tolerances where exposed.
Converged means lensed TT/EE/TE shift < 0.1% and C_L^phiphi shifts < 0.5% — an order of
magnitude below the model-to-model differences the science case rests on (Gaughan et al.:
~1-3% in C_ell, 5-15% in C_L^phiphi).

## Acceptance targets

- V0a/V0b binary pass (interpolation-level agreement).
- V1 <0.1% linear-power reconstruction.
- V3 reproduces the qualitative ordering and magnitudes of all three Gaughan figures.
- V5 runs an MCMC without babysitting.
- V6 convergence targets met with the production grid settings.
