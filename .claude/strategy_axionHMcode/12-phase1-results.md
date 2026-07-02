---
name: axionhmcode-phase1-results
description: Phase 1 measured results (2026-07-01) — R3 semantics settled, boost denominator derived and verified, V0c machine-precision null, first physics spectra, kmax insensitivity
metadata:
  type: project
---

# Phase 1 results (measured 2026-07-01)

Scripts: scratchpad phase1_semantics.py, phase1_lensedcls.py (migrated into the Theory
folder in Phase 2). Cosmology: input_file.txt-like (fax_dm = 0.1, m = 1e-25 eV, h = 0.674).

## R3 — transfer-variable semantics: SETTLED (numerically, both regimes)

- DM-like (m = 1e-25 eV): Transfer_tot = density-weighted cdm+baryon+AXION+massive-nu
  (median match 1.8e-6; every other combination fails at >= 3e-4).
  Transfer_nonu = cdm+baryon ONLY — the axion is NOT in nonu (match 3.2e-8).
  So: P_cold for axionHMcode = (omega_d T_cdm + omega_b T_b)/omega_db (== T_nonu);
  CAMB's total includes axion and nu.
- DE-like (m = 1e-33 eV): Transfer_tot = cdm+baryon+nu — axion REMOVED from tot
  (match 2.0e-6; cb+ax fails at 8e-2). Hlozek kluge confirmed as ported. KG-phase
  T_ax/T_cdm ~ -1e-4 at k = 0.01 falling to -6e-11 at k = 9 (horizon-scale smooth,
  finite, sign-flipped) — answers the R10-lead question about what Transfer_axion
  returns in the field phase. a_osc = 1.0 exactly for DE-like.

## Boost definition: derived closed form, verified

The Eq.9 model-linear-limit denominator is a perfect square, so the sqrt CAMB wants is
exact:

    sqrt(P_L_eq9) = (O_db/O_m) sqrt(P_cold) + (O_ax/O_m) [fc sqrt(P_cold)
                    + (1-fc) sqrt(P_ax)]           (fc = frac_cluster, scalar per z)

with B = P_NL_eq9 / P_L_eq9. Verified: B(k <= 0.01) = 1 within 6e-4 (dome) / 2e-5 (basic)
at z = 0 and 2 (residual = the damped one-halo (k/k*)^4 tail — physical, keep).
frac_cluster: 0.47 (z=0) -> 0.13-0.15 (z=2). Sample values, this cosmology:
z=0 dome B(0.5)=2.83, B(1)=11.1; z=2 dome B(0.5)=1.81, B(1)=3.17; basic systematically
lower. Reasoning recorded here because the denominator choice is load-bearing:
dividing by CAMB's P_L_tot instead would put a spurious few-1e-3 low-k tilt into the
lensed spectra (neutrino + composition mismatch); with the Eq.9-limit denominator the
low-k boost is 1 by construction and CAMB's own linear spectrum is untouched there.

## kmax truncation: NOT a concern

Boost recomputed from inputs truncated at k <= 10 h/Mpc vs k <= 100: relative change
<= 1.6e-4 for all k <= 8 (dome, z=0). The sigma(M)/HMF tails do not propagate to the
CMB-relevant boost. Production yamls can keep kmax ~ 10; V6 re-checks per point.

## V0c — code-path null: PASS at machine precision

ratio == 1 through (get_transfer_functions(only_time_sources=True) -> set_ratio ->
power_spectra_from_transfer) vs direct NonLinear_none get_results: lensed TT/EE agree to
~5e-13, C_L^phiphi to ~1e-15. Zero code-path systematic.

## First physics spectra (dome, fax_dm = 0.1, m = 1e-25 eV; boost vs linear lensing)

- C_L^phiphi: +1.9% (L=100), +19% (L=500), +45% (L=1000), +98% (L=2000).
- Lensed TT/EE: oscillatory, |dTT| < 0.4% to l = 2000, +1.1% at l = 2400 — the standard
  nonlinear lensing-smoothing signature. Magnitudes consistent with the known size of
  nonlinear corrections to CMB lensing.
- Boost grid on the production lensing grid (nz = 50, nk = 220): 75 s (1.5 s/z, dome,
  warm) single-threaded.

## Misc

- The Fortran warning "mismatch in integrated times (CAMB: CalcScalarSources)"
  (cmbmain.f90:1071) is upstream CAMB's own end-of-integration tolerance notice
  (|tau - tau_transfer| > 5e-5 Mpc at the last step for some k) — fractional effect
  ~1e-8; benign; present during the port's 0.1% validation. Not axion physics.
- V3 (full Gaughan figure reproduction, multiple masses + LCDM references) deferred to
  Phase 4 as planned.
