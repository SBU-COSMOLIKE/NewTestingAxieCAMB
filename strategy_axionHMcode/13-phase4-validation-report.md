---
name: axionhmcode-phase4-validation
description: Validation battery results (2026-07-01) — V0a/V0b/V0c/V1/V2/V3/V4/V6 measured numbers; V5 pending likelihood data
metadata:
  type: project
---

# Phase 4 validation report (measured 2026-07-01)

Environment: rayne/cobaya_test_env (cobaya 3.6.2, numpy 1.26.4, scipy 1.13.1,
numba 0.60, astropy 6.1.7); New_AxiECAMB on branch `hmcode` (prebuilt).
Scripts: axionhmcode_boost/dev_scripts/. Acceptance targets from file 07.

## Null and plumbing tests — all PASS

- V0a (trivial ratio=2 through cobaya): pk_nl = 4 x pk_lin to 1.1e-15.
- V0c (ratio==1 vs NonLinear_none): lensed TT/EE to ~5e-13, Cpp to ~1e-15.
- V0b (external mead2020 ratio through set_ratio vs internal mead2020, LCDM,
  26-node z grid): |dTT| <= 2e-7 (l <= 2000), |dCpp| <= 3.2e-5 (L <= 1000).
  Far below the 0.1% target — grid conventions, sqrt convention, ordering and
  clamping verified against an independent internal reference.
- V1 (transfer->P(k) convention): P_tot rebuilt through axionHMcode's own
  transfer_to_PS matches CAMB's linear P to 5.7e-8 max (float32-limited);
  sigma8 cross-check 0.8118 vs 0.8118.

## V2 — LCDM limit (cold-only halo model vs CAMB HMcode-2020, z=0, Planck18)

k [h/Mpc]:      0.2    0.5    1.0    2.0    5.0
mead2020:      1.158  2.365  5.686  12.62  27.06
dome flags:    1.413  2.693  6.399  14.78  31.61   (+12% at k=1)
basic flags:   1.176  2.057  5.103  12.21  25.02   (-10% at k=1)
Model spread ~ +-12% — the level the axionHMcode papers themselves claim for
LCDM agreement. (Note the fax->0 limit runs the cold-only code path below
omaxh2 = 1e-8, axionHMcode's own LCDM recipe; the full axion assembly is
singular at vanishing fraction.)

## V3 — Gaughan/Green/Moss (2605.12054) reproduction

Planck-2018 fiducial, fax = O_ax/O_D = 0.3 (their eq. 1). P_NL^ax/P_NL^LCDM
(their Fig. 1 style; R at k = 0.3/1/3 h/Mpc):

- m=1e-23: dome z2: 1.71/2.25/1.91 (boosted); basic z2: 0.92/0.73/0.74
- m=1e-24: dome z2: 1.58/1.48/1.05; basic z2: 0.91/0.47/0.46
- m=1e-25: dome z2: 0.89/0.31/0.20; basic z2: 0.65/0.11/0.07

Ordering (dome boosted above, basic suppressed below), signs, z-trend and
magnitudes all match their Fig. 1 structure (e.g. their basic m=1e-24 z=2
curve reads ~0.5 at k~1; ours 0.47).

Lensed spectra vs LCDM-mead2020, m=1e-24 (their Figs. 2-3 style):
- dome:  dTT = +0.11/+0.36/+0.97 % at l = 1000/2000/2400 (oscillatory,
  positive — their Fig. 2 red reaches ~+2-3% by l~3000);
  dCpp = +3.0/+17.4/+32.6/+46.8 % at L = 100/500/1000/2000
  (their Fig. 3 red: ~+35-40% at L ~ 1000-2000). MATCH.
- basic: dTT = +0.00/-0.07/-0.28 % (negative trend, their blue ~-1% at 3000);
  dCpp = +0.2/-2.6/-7.4/-17.2 % (their blue: ~-10..-25% at L~1000-2000). MATCH.
Residual quantitative differences are expected: their linear physics is
AxiCAMB (standard EFA lineage), ours AxiECAMB (Passaglia-Hu EFA).

## V4 — cross-check vs the hmcode-branch axion-aware internal HMcode

fax_dm = 0.1, m = 1e-25 eV: boost pipeline vs `halofit_version: mead2020` on
this branch (exact axion background, fully-clustering DM-like axion, NO
Jeans-scale halo physics): dTT <= 0.22% (l<=2000); dCpp = +0.6% (L=100),
+10% (L=500), +22% (L=1000). Two different halo models for the same
cosmology — differences in the anticipated band; no order-of-magnitude
anomaly. The mixed-DM model predicts more small-scale lensing power than
the axion-as-pure-CDM treatment at this mass/fraction.

## V6 — numerical convergence (m=1e-24, fax=0.3, dome, z=0 and 2)

- Halo-mass grid 100 -> 320 points: max |dB/B| = 7.6e-4 (median 1.4e-4).
- Mass range 7..18 -> 6..19 decades: identical to the density effect.
- Input k grid: the FULL transfer grid is converged — spline-densifying x2
  and x4 changes B by <= 1.3e-3 (median 3.5e-5), stable between x2 and x4.
  BUT halving the k grid corrupts B by up to 11% at k ~ 0.8-7 h/Mpc (the
  alpha/k_star transition machinery is sensitive to input-k sampling).
  RULE: never thin the k grid; always feed the full transfer k array.
- All shifts translate to well below the <0.1% TT / <0.5% Cpp targets.

## V5 — end-to-end with real likelihoods: PENDING

Requires `cobaya-install EXAMPLE_AXIONHMCODE_EVALUATE1.yaml -p <packages>`
(multi-GB likelihood data download; no packages path exists on this machine).
The full wiring was smoke-tested with a Cl-requiring stand-in likelihood:
get_model 0.3 s, one evaluation 83.6 s (dome, 50-node grid, single core),
sensible lensed spectra and derived parameters (H0 = 63.4 from theta_star
shooting through the DM-like axion background at omaxh2 = 0.012).

## Incidental findings

- Upstream cobaya >= 3.6 removed `use_renames` (renames now unconditional):
  the OLD EXAMPLE_EVALUATE1/MCMC1.yaml fail on pristine cobaya >= 3.6 (they
  still work on Cocoa, whose replacement camb.yaml defines the key). The new
  yamls omit the flag. Also `derived: true` + `value:` lambda is rejected by
  3.6.2 — the function now goes directly on `derived:`.
- The Cl-only cobaya path leaves WantTransfer off; the boost theory needs
  matter transfers, so the example yamls set `nonlinear: NonLinear_both` and
  `kmax: 10` in extra_args (kmax flips WantTransfer on via set_matter_power).
- cobaya Pk_grid returns k in 1/Mpc (not h/Mpc) — bit us once in a test.
