---
name: axionhmcode-objective
description: Main objective — lensed CMB + lensing spectra from AxiECAMB with axionHMcode nonlinear boost, via Cobaya inside Cocoa
metadata:
  type: project
---

# Objective and deliverables

## Main objective (from prompt_draft_hmcode.tex)

Compute lensed C_ell^TT, C_ell^EE, C_ell^TE and the lensing potential C_L^phiphi with
AxiECAMB (the ported ultralight-axion CAMB in `New_AxiECAMB/`, see [[axiecamb-port-project]]),
where the nonlinear matter power spectrum comes from axionHMcode (`axionHMcode/`) instead of
CAMB's internal halofit/HMcode. The coupling is through the boost factor
B(k,z) = P_NL(k,z)/P_L(k,z); CAMB's `ExternalNonLinearRatio` ingests x = sqrt(B) via
`set_ratio(k_h, z, ratio)`. Everything must run through the Cobaya interface (MCMC-able),
inside the Cocoa framework (`cocoa/`).

Physics motivation and validation target: Gaughan, Green & Moss 2026 (arXiv:2605.12054)
integrated axionHMcode into AxiCAMB natively and showed the nonlinear prescription
(naive HMcode vs axionHMcode basic vs DOME) shifts ULA constraints by more than the
statistical preference itself for m ~ 1e-25..1e-23 eV. We are building the same physics
pipeline, but for AxiECAMB (Passaglia-Hu EFA) and through Cobaya Theory blocks.

## Deliverables

1. A Cobaya `Theory` class (Python) that wraps axionHMcode and provides `non_linear_ratio`
   to the stock Cobaya CAMB wrapper (`use_non_linear_ratio: True`). Template for mechanics:
   `TrivialNonLinearRatio` in cobaya's `tests/test_cosmo_multi_theory.py:279`. Template for
   group style: `class emulmps(Theory)` from CosmoLike/emulators_code (copy in scratchpad;
   re-fetch from https://github.com/CosmoLike/emulators_code/blob/main/emulmps/emulmps.py).
2. YAML interfaces: an EVALUATE and an MCMC yaml modeled on `New_AxiECAMB/EXAMPLE_EVALUATE1.yaml`
   and `EXAMPLE_MCMC1.yaml`, targeting m ~ 1e-25..1e-23 eV (PI, 2026-07-01; the existing
   examples sample logmx in [-34,-31] — DE-like — where axionHMcode does not apply).
   EVALUATE example: log-mass prior block + mass pinned in the evaluate override.
   MCMC example: fixed mass, with comments showing the uniform log-mass alternative.
   Plus README instructions covering both modes — recipes ready in
   [[axionhmcode-mass-prior-recipes]] (file 09).
3. Staging into Cocoa: New_AxiECAMB at `external_modules/code/AxiECAMB` (the path used in the
   prompt's yaml example) and axionHMcode importable by the Theory class.
4. Validation battery (see [[axionhmcode-validation-plan]]).
5. Append a project report to `New_AxiECAMB/README.rst` (do not rewrite existing content).

## Constraints (user-imposed)

- Any Cobaya modification must be a patch applied inside
  `cocoa/Cocoa/installation_scripts/setup_cobaya.sh`, AFTER the existing patches in
  `cocoa/cocoa_installation_libraries/cobaya_changes/`; existing patches must not be edited.
  (Finding: zero new Cobaya patches are expected to be needed — see
  [[axionhmcode-no-circular-dependency]].)
- Solution with the least patches preferred; a pure Theory-block solution is "ideal"
  per the prompt — and it is achievable.
- Group rules: no `.f90` edits without telling the user first with a reason
  ([[axiecamb-port-project]]); Theory-class Python in group style (2-space indent, concise).

## Paper map (rayne/papers/ — sibling of this repo, not committed)

- arxiv_2412.15192.pdf — AxiECAMB physics (anchor)
- arxiv_2605.12054.pdf — Gaughan/Green/Moss: nonlinear modelling role; primary validation target
- arxiv_2409.11469.pdf — Dome et al.: DOME calibration of axionHMcode
- arxiv_2209.13445.pdf — Vogt et al.: the original axionHMcode paper (added; the prompt cited
  only the DOME paper for axionHMcode)
- arxiv_2501.13662.pdf — Moss et al.: PH-EFA AxiCAMB
- arxiv_2201.10238.pdf — Passaglia & Hu EFA (theory base of AxiECAMB)
- arxiv_astro-ph_0003365.pdf — Hu, Barkana & Gruzinov fuzzy DM
- arxiv_2005.05290.pdf — Cobaya paper
- CAMB_notes.pdf — CAMB notes
