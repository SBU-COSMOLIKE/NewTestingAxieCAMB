---
name: axionhmcode-performance-anatomy
description: Profiled explanation of why axionHMcode costs seconds per redshift while Fortran HMcode costs milliseconds — nested root-finds recomputing untabulated growth and sigma(M) integrals; upstream speedup estimates
metadata:
  type: project
---

# Performance anatomy: why axionHMcode is ~1000x slower than HMcode

Measured 2026-07-02 (cProfile, one dome redshift, z=0, nk=212, M-grid 100;
script: axionhmcode_boost/dev_scripts/profile_axionhmcode.py). Context: the ACT DR6
axion analysis flagged the same slowness and built the axionEMU emulator around it.

## Where the time goes (stage level, per redshift)

| stage | dome | basic |
|---|---|---|
| HMCode_param_dic (cold HMcode parameters) | 0.000 s | 0.000 s |
| func_axion_param_dic (axion cut mass, soliton central density, frac_cluster) | 2.35 s | 1.22 s |
| func_full_halo_model_ax (the actual halo-model sums) | 0.43 s | 0.33 s |

The cold-HMcode parameter stage is instantaneous (already numba-vectorized by the Dome
update). ~85% of the runtime is the AXION parameter stage, ~15% the honest halo-model
integrals.

## The mechanism (function level, one dome redshift = 2.7M function calls)

| function | calls | tottime |
|---|---|---|
| func_D_z_norm (growth factor, re-integrated per call incl. its constant D(0) normalization) | 109,334 | 0.86 s |
| func_nu = delta_c/sigma(M) (full 212-point variance integral per call, no table) | 99,055 | 0.77 s |
| func_conc_param / formation-z machinery (cold_density_profile.py) | 99,052 | 1.88 s cum |
| scipy brentq root-finds | 6,451 | 1.56 s cum |
| scipy quad/dblquad (growth integrals, G_a-type) | 626 | 1.17 s cum |
| func_E_z | 246,666 | 0.04 s |

The call cascade (verified by reading the sources):

    func_cut_mass_axion_halo:            brentq over log10 M in [7, 17]
      -> func_jeans_virial_ratio(M_trial)
        -> func_halo_jeans_kscale(M_trial)
          -> func_conc_param(M_trial)     [off-grid M -> no reuse possible]
            -> func_formation_z_HMcode2020: ANOTHER brentq over z_f in [z, 100]
              -> each iteration: func_D_z_norm(z_f)  [re-integrates growth,
                 re-computes the constant D(0) normalization every call]
              -> and func_nu(f*M)          [re-integrates sigma(M) over the
                 whole k grid every call]

and the same pattern again inside the soliton central-density solver
(axion_density_profile: MaxofMc / func_central_density_param). Root-finders nested
inside root-finders, with numerical integrals recomputed from scratch at the bottom of
every innermost iteration, and zero memoization anywhere: D(0) — a run constant — is
recomputed 109k times per redshift; sigma(M) — a smooth one-argument function — is
re-integrated 99k times instead of being splined once on ~100 masses.

## Why Fortran HMcode does the same class of physics in ~milliseconds

Not primarily the language. HMcode-2020 (in CAMB) tabulates sigma(M) once per redshift
on a small mass grid and interpolates; obtains the formation redshift by INVERSE
INTERPOLATION of the tabulated growth (no root-finder); uses closed-form fitted
delta_c/Delta_v; and its c(M) is then a one-line formula over those tables. One pass
over small tables per redshift. axionHMcode's extra axion physics (which halos can hold
axions; how dense the soliton core is) is genuinely implicit — it needs SOME iteration —
but the current implementation solves each implicit condition by brute force at every
trial point of every enclosing solver. The ~1000x gap factorizes roughly as: ~1 order of
magnitude interpreter/call overhead x ~2 orders of magnitude redundant recomputation.

## Timing accounting in cobaya (the "camb: 84 s" confusion)

Cobaya bills wall time to the component whose calculate() is executing. The boost
theory's own calculate() is a no-op (measured 1.6e-6 s); all axionHMcode work happens in
get_non_linear_ratio(results), which the CAMB wrapper invokes from inside ITS OWN
calculate() (camb.py:717). So in `timing: true` output, ~50 z x ~1.6 s of axionHMcode
appears under "camb", while "camb.transfers" (~0.86 s) is the actual Boltzmann solve and
"axionhmcode_boost" shows microseconds. CAMB itself is not slow.

## Speedup implementation (fork_axionHMcode, PI sign-off 2026-07-02)

IMPLEMENTED in the group fork (github.com/SBU-COSMOLIKE/axionHMcode, folder
rayne/fork_axionHMcode; all edits fenced with `# VM-SPEEDUP`, grep lists every site):

1. `cosmology/fast_tables.py` (new): growth D(z) table (4001 nodes, node values from
   the ORIGINAL upstream quadrature), integrated-growth G(z) table (cumulative
   trapezoid replacing scipy.dblquad), sigma(M) log-log table (321 nodes from upstream
   func_sigma_M); per-evaluation caches ride inside cosmo_dic ('_vm_' keys).
2. `func_z_formation`: monotone inverse interpolation of the growth table replaces the
   brentq cascade; clamp branches replicated exactly (decided, never interpolated).
3. Exact memoization of `func_conc_param` (scalar) and `NFW_profile` (per mass +
   r-grid fingerprint) — the central-density root solver re-evaluated identical
   quantities ~30x per mass.
4. `func_D_z_unnorm_int`: G-table instead of dblquad (0.15 s -> ~0 per call).

MEASURED (dev_scripts/fork_validate.py, fork vs unmodified upstream a85ba26):
- Numerics: max |dB/B| = 1.6e-5 across dome/basic x z=0/2 x three cosmologies incl.
  the LCDM limit (gate was 1e-4); deviations IDENTICAL before/after the memoization
  round (proving those are exact); growth/G primitives <= 1.6e-6 vs original
  quadratures; boost pytest suite 4/4 against the fork.
- Speed per redshift: dome z=0 3.5 -> 1.07 s (3.3x), basic 2.0 -> 0.62 s (3.2x),
  z=2 ~1.6x, LCDM path 13-22x. Per likelihood evaluation (50-node grid): dome
  ~135 s -> ~52 s, before `processes: N` division.
- REMAINING structure (profiled): ~0.7 s/z = the soliton central-density solver
  (optimize.root per halo mass; ~28 profile evaluations each even with cached NFW —
  simpson-call overhead dominates); ~0.4 s/z = the (M,k) profile Fourier transform in
  func_full_halo_model_ax. Untouched deliberately: replacing the root solver changes
  its guess-based no-solution rejection semantics (behavior, not just precision), and
  reimplementing scipy's non-uniform Simpson weights exactly is version-fragile.
  Next-round candidates if more speed is needed; else the emulator route (ACT's
  choice) remains the production path for large campaigns.

## Round-2 direction (PI, 2026-07-02 evening)

The version-fragility objection above is scoped correctly but narrowly: it applies to
CLONING scipy's simpson weights on a fixed externally-given grid (the even-N last
interval policy changed in scipy 1.11 and the escape keyword was removed in 1.14; on
our even-length geomspace grids the old and new defaults differ by ~1e-6-1e-7
relative). It does NOT apply to Gauss-Legendre quadrature, whose nodes and weights
are canonical mathematics — the PI's preferred pattern, used throughout cosmolike's
cosmo2D.c (gsl_integration_glfixed_table + create_cosmo_nodes: evaluate all expensive
functions once at the GL nodes, hot loops become weighted sums, OpenMP + SIMD outside;
node count = Ntable.high_def_integration is the accuracy knob). Python translation:
numpy-vectorized integrand evaluation at GL nodes plays SIMD, `processes:` plays
OpenMP. A GL rewrite of the soliton r-integrals (in ln r; the integrands are smooth)
would cut integrand evaluations ~1000-2000 -> ~30-100 nodes, not just the 52x
per-call bookkeeping overhead — but it changes WHERE upstream samples the integrand,
so it moves those integrals from bit-exact to converged-with-a-knob and needs the
fork_validate gate re-run plus a 1-vs-2 convergence check.

That knob now exists: `accuracy_boost` yaml option in the boost theory (PI request,
2026-07-02; CAMB-style single multiplier). Scales m_grid_points and, via
fast_tables.set_accuracy_boost() (forwarded by _compute_row; upstream checkouts
ignore it), the fork's growth/G/sigma table node counts. Measured boost 1 vs 2
(dome z=0, gaughan/inputfile/lcdm): max |dB/B| = 0.7-1.1e-3, dominated by the
halo-mass grid — consistent with V6's M-grid sensitivity; fork tables contribute
<~2e-6. Cost scales ~linearly.

## Round-2 OUTCOME (2026-07-02/03, PI go-ahead "and yes")

The GL rewrite of func_ax_halo_mass was implemented, measured, and REJECTED by
the gate — replaced by a design that keeps the speed and passes:

1. GL in ln r (128 nodes) failed fork_validate: gaughan_z0_dome 4.3e-3,
   inputfile_z2_dome 6.8e-3 (gate 1e-4). Node-count scan (128/256/512/1024)
   showed non-monotone jumps up to 1.9e-2 at k ~ 8-16 h/Mpc -> not quadrature
   error but discrete flips: func_dens_profile_ax composes soliton-vs-NFW by
   sign changes of the sampled difference, so the composed profile depends on
   node positions; grazing crossovers flip with the grid. Behavior change, not
   precision. Grid-independent crossover detection would fail the gate the same
   way (upstream itself grid-snaps).
2. ADOPTED design (fast_tables.geom_simpson_grid): keep upstream's exact
   np.geomspace(1e-15, r_vir, 2000) nodes (bit-identical composition, snapping,
   rejection) and integrate the identical samples with a precomputed canonical
   weight vector: unevenly-spaced composite Simpson parabola pairs + trapezoid
   on the leftover final interval. Textbook math, no scipy call in the hot
   loop, no version-policy dependence; deviation from scipy's Cartwright tail
   confined to the outermost interval, far below the gate.
3. The 2000-point radial grid is PINNED — excluded from accuracy_boost: with
   the grid scaled, boost 1-vs-2 showed up to 6.7e-2 (gaughan_z0 dome) from
   crossover re-snapping + rejection flips. That sensitivity is a property of
   the released calibrated model (Dome calibration ran with this grid), so the
   grid is model definition, not a fork convergence parameter. With the grid
   pinned, boost 1-vs-2 is back to 0.7-1.1e-3 (M-grid dominated). The ~7e-2
   high-k grid sensitivity of the model itself is worth flagging to
   collaborators/upstream alongside the file-10 flag-configuration question.
4. Gate PASSED: worst max|dB/B| = 1.59e-5 across all 12 cases (same level as
   round 1); primitives unchanged (<=1.6e-6); pytest 4/4.
5. Timings (warm, single core): dome z=0 3.5 -> 0.73 s (4.9x), basic z=0
   2.0 -> 0.49 s (4.0x), z=2 ~2.1-2.4x, LCDM 15-22x. Dome evaluation
   ~135 (upstream) -> ~36 s at 50 z-nodes, before process division.
6. Parallelism = OMP_NUM_THREADS, always (PI decision 2026-07-03: the yaml
   must not be able to overrule OpenMP — the `processes` option was removed
   entirely after one day; workers = OMP_NUM_THREADS, fallback 1, no knob).
   Note: Cocoa's axicambv2 EXAMPLE_EVALUATE2/MCMC2 yamls still carry
   `processes: 1` lines — they must be deleted when the PI repins axiecamb
   (report-only tree). Full record: fork README appendix A.7-A.11.

Remaining bottleneck structure (fork, dome z=0 ~0.73 s wall; re-profiled
2026-07-03 post-round-2): the central-density solver stage 0.46 s cum (of
which 0.26 s = ~26 hybr objective evaluations per mass x 100 masses, and
0.19 s = the guess construction: 21-point adaptive quad of the unit-amplitude
soliton + 1000-pt simpson of the NFW part) and the (M,k) profile Fourier
transform 0.38 s cum (0.17 s = building the sin(kr)/(kr) kernel on 212x2000
arrays per mass, ~0.1 s composed-profile arrays, ~0.1 s simpson reduction).
CORRECTION to the earlier note: the FT grid is NOT branch-free — its
integrand is the same composed profile (crossover detection on the samples),
so that grid is behavior-tied and must stay pinned like A.7; safe FT
optimizations are same-samples reductions only. The solver stage is
semantics-locked under the match-upstream gate (hybr failure behavior +
|guess-rho_c|>100 rejection are part of the effective calibrated model);
gate-safe crumbs: exact memoization of func_delta_char/func_r_vir per mass
would cut the guess-quad overhead (~0.19 -> ~0.05 s) without changing any
value.

## Round-3 OUTCOME (2026-07-03, PI: "yes, don't forget to pad")

FFTLog j0 transform implemented for func_dens_profile_ax_kspace
(fast_tables.fftlog_j0_grid/fftlog_j0_eval), on the pinned upstream samples
(composition bit-identical); the cfftlog/cosmo2D.c algorithm with Bessel
order fixed to zero (no per-k kernel loop; their per-ell phase-2 has no
analog here). Design settled by prototype iteration (dev_scripts/
fork_fftlog_edge_scan.py, fork_fftlog_interp_scan.py):
- bare FFTLog missed by 5.4e-3: the truncation edge at r_vir (where
  F = rho r^3 peaks) integrates O(dlnr) differently in a Fourier basis than
  in simpson. Fix: fill the padding with the constant F(r_vir) (continuous,
  bias-killed by e^{-36q} before the wrap) and subtract its tail
  analytically, F(rvir)[sin(a)/a - Ci(a)], a = k r_vir, via sici.
- remaining 3.3e-5 at high k was linear interpolation of the oscillatory
  I(ln k); 4-point Lagrange (cubic, closed-form on the uniform ln k grid)
  -> max |dI|/I(0) = 5e-6 (1.7e-5 on a second grid), insensitive to bias
  q (0.9), window (0.25), n_fft (4096).
- padding does double duty (PI emphasis): power-of-two FFT length AND
  extends the reciprocal k grid ~36 e-folds below 1/r_vir, so one
  transform covers k down to 1e-4.
- FFT count: 2/mass = 200/z = 1e4/eval; ~0.11 ms/mass eval + 0.36 ms/mass
  grid setup (Mellin loggammas, cached per grid in cosmo_dic). Plan reuse:
  pocketfft's per-length twiddle cache (all transforms share one length) =
  the static-FFTW-plan trick of cosmo2D.c; not additionally threaded
  (process-level z-loop fork owns the OMP_NUM_THREADS core budget).
Gate PASSED: worst max|dB/B| = 1.59e-5 (unchanged — still table-carried);
1-vs-2 accuracy_boost unchanged (0.7-1.1e-3); pytest 4/4. Timings: dome
z=0 2.9-3.5 -> 0.34 s (~8x), basic ~7x, z=2 ~4x, LCDM 16-24x; dome eval
~135 -> ~17 s single-threaded. Remaining cost: the central-density solver
structure only (semantics-locked, see above) + the gate-safe
delta_char/r_vir memo crumb.

FFT-backend benchmark (2026-07-03, PI asked about pyfftw/FFTW plan reuse;
script dev_scripts/fork_fftw_vs_pocketfft_bench.py): 4096-pt rfft+irfft
pair on this mac — pyfftw MEASURE-plan reuse 15.1 us (plan 516 ms once),
ESTIMATE 16.0 us (plan 1.3 ms), scipy pocketfft 29.5 us, numpy pocketfft
43.2 us; batched (100,4096) forward 614/870/1136 us. FFTW is ~2x faster
per transform, but the whole FFT share is ~3 ms/z (~1% of 0.34 s/z), so
end-to-end gain ~0.1% — dependency not justified (macOS build was DIY: no
py312/arm64 wheel; built 0.15.0 against miniforge's libfftw3 with cython
+ a filtered lib dir to dodge Apple clang's missing -fopenmp; left
installed in cobaya_test_env for reference only). ADOPTED instead: swap
np.fft -> scipy.fft in fftlog_j0_eval (same pocketfft library, bit-
identical results, ~1.5x on the FFT share for free); gate re-run
byte-identical, pytest 4/4. If FFTW is ever wanted: 15-line optional
backend shim in fast_tables (try-import pyfftw, plan cache keyed by
shape, wisdom export across workers), zero hard dependency.
