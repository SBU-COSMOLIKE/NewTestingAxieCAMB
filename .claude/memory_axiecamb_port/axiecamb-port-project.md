---
name: axiecamb-port-project
description: Ongoing project — port AxiECAMB (ultralight-axion CAMB Nov13 fork) to modern CAMB 1.6.7 in New_AxiECAMB/
metadata: 
  node_type: memory
  type: project
  originSessionId: ce7b9b4b-a767-493b-bf09-b247a44c375f
---

Started 2026-06-09; **port completed same day**. AxiECAMB (arXiv:2412.15192; axionCAMB/CAMB-Nov13 base) ported to modern CAMB 1.6.7 in `New_AxiECAMB/`, Fortran + Python. Validated vs original: axion/LCDM suppression ratios agree to ≤0.1% (TT, ℓ=2–2600) and ≤0.01% (P(k) where not exponentially dead) across DM-like pre/post-recombination-switch and DE-like regimes; LCDM limit bit-identical to CAMB 1.6.7; python test suite passes. Key files: fortran/AxionBackground.f90 (new TAxionModel component in CAMBparams%Axion), edits in results/equations/cmbmain/recfast/halofit/camb.f90, camb/axion.py, README_AxiECAMB.md, inifiles/params_axion.ini. Derived axion outputs read from `results.Params.Axion` (state copy), not the input pars. Isocurvature force-disabled (v1.0 parity).

Key facts:
- Public repo: `New_AxiECAMB/` is the local clone of `git@github.com:SBU-COSMOLIKE/NewTestingAxieCAMB.git` (branches `main` + `hmcode`). Its readme was Cocoa-restructured 2026-07-02: `README.md` is the user-facing readme (overview/run/appendices) and `PORT_DEVELOPER_GUIDE.rst` holds the generated developer dossier (Phase-1 reports/Phase-3 design turned into docs; ~335 KB); `.port_analysis/` holds the source diffs/reports/rst.
- Working-tree layout under `rayne/`: `AxiECAMB/` = original Nov13 fork; `OLDCAMB/` = pristine CAMB Nov13 baseline for diffing; `CAMB/` = pristine modern 1.6.7 baseline; `New_AxiECAMB/` = the port; `cobaya/` = Cocoa Cobaya CAMB-wrapper fork.
- `OLDCAMB/` = pristine CAMB Nov13 baseline for diffing; real (non-whitespace) AxiECAMB changes: equations_ppf.f90 ~1420 lines, modules.f90 ~720, cmbmain.f90 ~543; new files axion_background.F90 (1429 ln), recfast_axion.F90, inidriver_axion.F90, halofit_ppf.f90.
- Physics: KG equation for ULA background+perturbations until m/H=10, then effective fluid (time-averaged oscillations); m/H0<10 ⇒ DE-like, KG to today. Isocurvature params exist but disabled in v1.0. Params: m_ax (eV), omaxh2 / (omdah2, axfrac, use_axfrac).
- Analysis diffs & agent reports in `.port_analysis/` (diffs/ and reports/).
- Toolchain: gfortran 14.3 at /Users/vivianmiranda/miniforge/bin/gfortran; CAMB makefiles race with -j4 (forutils + main); build serially. Python 3.9.6 at /usr/bin/python3.

Cobaya integration (2026-06-10): camb.set_params chain + get_valid_numerical_params now route m_ax/omaxh2/omdah2/axfrac/dfac (Axion.set_params called before set_cosmology so theta→H0 shooting includes the axion). Chain yaml at rayne/EXAMPLE_MCMC_axions_newaxiecamb.yaml (replaces the emulator blocks of arXiv:2510.14957's EXAMPLE_MCMC7): sample thetastar100 (CAMB wants θ*, not 100θ*; lambda /100), logA→As, omegaaxh2→omaxh2, logmx→m_ax lambdas with drop, use_renames for omegabh2/omegach2, extra_args theta_H0_range [40,130] + halofit_version original + nnu 3.046. User rule: NO .f90 edits without telling them + reason first. Cobaya fork wrapper at rayne/cobaya/theories/camb (Cocoa); tested locally with cobayapristine env (~1.5s/eval).

Branch `hmcode` (2026-06-10, requested by Wayne Hu, items 1+2 of my HMcode proposal): halofit.f90 only. Exact axion background (Ω_ax(a), w_ax(a) tables from GrhoAx/FieldValsAta) added as separate component to HMcode's internal Hubble2/AH (axion share removed from the (1+z)² closure remainder), hence exact growth ODE; DM-like axions count as fully-clustering cold matter in Omega_m_hm/Omega_cold_hm/f_nu/EH99-Tcb/cosmic_density/Dolag-LCDM-reference (no Jeans suppression — input P(k) carries scale dependence); DE-like → expansion only. All gated on HasAxion(): LCDM+mead byte-identical to upstream (verified). New mead vs classic-original axion suppression agrees 1-8%; consistency fix worth ≤2.5% at f_ax=10%; negative-Λ extreme (ωax=0.39 DE-like) stable. recfast→CosmoRec needs NO axion work (cosmorec.f90:160 passes Hz array from dtauda; hyrec exports dtauda callback).

Related: [[full-file-access-granted]].
