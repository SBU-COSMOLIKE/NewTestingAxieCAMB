---
name: axionhmcode-api
description: axionHMcode call sequence, input dictionaries, unit conventions, versions, nuisance parameters
metadata:
  type: project
---

# axionHMcode API (as cloned in axionHMcode/)

Physics: Vogt et al. 2209.13445 (basic), Dome et al. 2409.11469 (dome calibration + speedups).
Mixed-dark-matter halo model: P = (Ocb/Om)^2 Pcb + 2 OcbOax/Om^2 Pcb,ax + (Oax/Om)^2 Pax,
with the axion split into halo-bound (clustered) and linear (unclustered) parts below/above
the Jeans scale (their Eq. 9; README "How Does the Code Work?").

## Call sequence (from README + example notebook + signatures)

1. `cosmo_dic` — plain dict (bypass `load_cosmology_input`; build it directly). Required keys
   (load_cosmology.py): omega_b_0, omega_d_0 (CDM ONLY — axion excluded), omega_ax_0, m_ax,
   h, z (SCALAR — one redshift per call), ns, As, k_piv (1/Mpc), M_min, M_max (LOG10
   exponents of Msun/h — RESOLVED by the notebook: `np.logspace(M_min, M_max, 100)` with
   7/18 from input_file.txt; the linear defaults 1e8/1e17 in load_cosmology.py are
   inconsistent with that usage — upstream quirk, do not copy them),
   transfer_kmax (default 1e3), version ('basic'|'dome'), plus derived omega_db_0/Omega_*_0,
   Omega_w_0 = 1 - Omega_m_0, and G_a = func_D_z_unnorm_int(z, Omega_m_0, Omega_w_0)
   (load_cosmology.py:204-205) — their internal scale-independent growth; keep verbatim, it
   is part of the calibrated model even though it is LCDM-approximate for axion cosmologies.
   Optional nuisance keys: alpha_1, alpha_2, gamma_1, gamma_2 (Dentler et al. 2111.01199
   Eq. 36; priors used by Gaughan et al.: [0.6,2.0], [1.43,2.54], [5.0,45.0], [-0.37,-0.23]).
2. `power_spec_dic` — from transfer functions via `func_power_spec_dic`
   (axionCAMB_and_lin_PS/lin_power_spectrum.py:51): keys k, power_total, power_CDM,
   power_baryon, power_cold (T_cold = (Ob*T_b + Od*T_cdm)/Odb), power_axion.
   `transfer_to_PS` (lin_power_spectrum.py:15): P(k) = T^2 (k h)^4 P_prim(k) 2 pi^2 / k^3,
   T in the axionCAMB/CAMB convention (transfer divided by k^2, Mpc^2), k in h/Mpc,
   P in (Mpc/h)^3. CAMB's `get_matter_transfer_data()` provides the same convention —
   verify numerically once (V1). The example notebook passes ONE power_spec_dic instance
   everywhere (the separate sigma-dict mentioned in docstrings is not exercised) —
   mirror the notebook.
3. `hmcode_dic = HMCode_param_dic(cosmo_dic, k, PS_cold)` (halo_model/HMcode_params.py:107).
4. `axion_dic = func_axion_param_dic(M, cosmo_dic, power_spec_dic, hmcode_dic,
   concentration_param=...)` (axion_functions/axion_params.py:10) — cut mass M_cut, M_ax(M),
   clustered fraction, central density.
5. `P_NL = func_full_halo_model_ax(M, power_spec_dic, cosmo_dic, hmcode_dic, axion_dic,
   alpha=..., eta_given=..., one_halo_damping=True, two_halo_damping=..., 
   concentration_param=..., full_2h=...)` (halo_model/PS_nonlin_axion.py:13).
   M is the halo-mass integration grid (Msun/h; notebook uses 100 logspaced points; the
   inline comment notes the HMcode reference used 129 and prefers 1025 — a V6 convergence
   knob). Returns a tuple — element [0] is the total nonlinear P(k).

## Flag conventions

- basic version: alpha=False (i.e. alpha=1 smoothing), HMcode-2020 params off; calibrated
  fax<0.5, wide mass range 1e-33..1e-21 eV.
- dome version: alpha=True and concentration_param=True recommended (README); calibrated
  0.01<fax<0.3, 1<z<8, m near 1e-24.5 eV; full_2h=False was used in calibration.
- README warns HMcode-2020 parameter switches are NOT calibrated for MDM — keep the
  per-version defaults from the example notebook; do not invent flag combinations.

## Performance notes

- Performance is NOT a design constraint (PI, 2026-07-01) — accuracy is. One halo-model
  evaluation per redshift (z scalar in cosmo_dic); the z-grid policy is dense (every z CAMB
  uses, [[axionhmcode-architecture]]). The Phase-0 single-z benchmark is informational only
  (sets MCMC wall-time expectations), not a design gate.
- numba @njit in halo_bias.py and variance.py (Dome update) — first call pays JIT
  compilation. RESOLVED: numba 0.60 is pinned in all Cocoa env files (cocoapy310*.yml);
  no new dependency needed.
- Gaughan et al. ran full MCMCs with this code inside AxiCAMB (accuracy=2.5, m/H*=50),
  so per-step cost is known to be tractable; ACT DR6 used an emulator instead. If cost ever
  becomes prohibitive, the PI's plan is to train ML emulators USING this pipeline as the
  training-data generator (emulmps-style, swapping the internals of get_non_linear_ratio,
  same Theory interface) — a later phase, not part of this build.

## Gotchas

- `load_cosmology_input` uses eval() on the input file — another reason to build cosmo_dic
  directly in the Theory class.
- Massive neutrinos are absent from axionHMcode's matter budget (omega_m = b + CDM + ax);
  see [[axionhmcode-open-questions]] on the ratio convention with mnu = 0.06 fixed.
- Units: all k in h/Mpc, all P in (Mpc/h)^3, masses in Msun/h; As dimensionless at k_piv
  in 1/Mpc (their primordial_PS multiplies k by h — lin_power_spectrum.py:12).
