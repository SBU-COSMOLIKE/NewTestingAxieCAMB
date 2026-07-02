# AxiECAMB ⊕ axionHMcode via Cobaya — developer guide

This is the complete technical record of how the axionHMcode mixed-dark-matter
nonlinear boost was merged with AxiECAMB through a Cobaya `Theory` block
(built and validated 2026-07-01/02). It plays the same role for the boost that
[PORT_DEVELOPER_GUIDE.rst](PORT_DEVELOPER_GUIDE.rst) plays for the CAMB-1.6.7 port:
everything a developer needs to audit, maintain, or extend the integration —
the mechanism, the conventions and why they were chosen, the implementation
walkthrough, the traps, and every measured validation number.

The raw working documents behind this guide (design notes, decision log, phase-by-phase
measured results) live in [`.claude/strategy_axionHMcode/`](.claude/strategy_axionHMcode/00-INDEX.md).
The implementation is [`axionhmcode_boost/`](axionhmcode_boost/README.md) (user-level
options are documented in its README). File:line references were verified against:
this repository at commit `91771d5`; Cobaya 3.6.2 (= Cocoa's pinned commit `899f30a4`);
the axionHMcode checkout at upstream commit `a85ba26`.

# Table of contents
1. [The problem and the objective](#problem)
2. [The central mechanism: why there is no circular dependency](#mechanism)
3. [End-to-end data flow](#dataflow)
4. [The boost definition (the load-bearing convention)](#convention)
5. [Transfer-variable semantics (verified, both regimes)](#transfers)
6. [The redshift grid](#zgrid)
7. [Implementation walkthrough: AxionHMcodeBoost](#implementation)
8. [Traps catalog](#traps)
9. [Validation battery — measured results](#validation)
10. [Performance](#performance)
11. [Validity domain and known limitations](#validity)
12. [Cocoa integration](#cocoa)
13. [Decision log](#decisions)
14. [File inventory](#files)

# 1. The problem and the objective <a name="problem"></a>

Goal: compute lensed C_l^TT/TE/EE and the lensing potential C_L^phiphi for
ultralight-axion cosmologies where the nonlinear matter power comes from
**axionHMcode** (Vogt et al. [arXiv:2209.13445](https://arxiv.org/abs/2209.13445);
Dome et al. recalibration [arXiv:2409.11469](https://arxiv.org/abs/2409.11469))
instead of halofit/HMcode — through the Cobaya interface, so MCMC runs are possible,
with the fewest possible modifications to Cobaya (ideally none). Science context:
[arXiv:2605.12054](https://arxiv.org/abs/2605.12054) showed the choice of nonlinear
prescription in the m_ax ~ 1e-25..1e-23 eV window shifts ULA constraints by more than
the statistical preference itself.

The apparent blocker: axionHMcode is a halo model, not an emulator. It needs the
linear P(k,z) from CAMB to produce P_NL(k,z), and CAMB then needs
sqrt(P_NL/P_L) back before it can compute lensed spectra. Naively that is a
CAMB → axionHMcode → CAMB sandwich — a cyclic dependency that Cobaya's
requirement resolver would reject.

# 2. The central mechanism: why there is no circular dependency <a name="mechanism"></a>

Two facts, both invisible from the yaml, dissolve the cycle.

**Fact 1 — the `camb` block is two graph nodes, not one.** At initialization the
Cobaya CAMB wrapper spawns a helper component named `camb.transfers`
(`CAMB.get_helper_theories`, cobaya `theories/camb/camb.py:1049-1061`; helper class
`CambTransfers` at `camb.py:1162`). The helper runs the Boltzmann solve once per point
and provides the quantity `CAMB_transfers` (the transfers-level `CAMBdata`). The main
`camb` node never runs the Boltzmann code: its `calculate()` takes the cached transfers
and finishes with `power_spectra_from_transfer()` — folding the primordial spectrum into
the precomputed transfer functions and applying the nonlinear scaling. The split is
physically exact because the transfer functions T(k,z) and the CMB time sources do not
depend on A_s, n_s, or the nonlinear model. It also preserves Cobaya's fast/slow
caching: the helper result is reused when only primordial or nonlinear parameters
change (docstring, `camb.py:1164-1166`).

**Fact 2 — the linear P(k) never enters the dependency graph.** With
`use_non_linear_ratio: True` (a top-level `camb` option, `camb.yaml:26` — *not* an
`extra_args` entry), the wrapper:

- forces the nonlinear model to `ExternalNonLinearRatio` (`camb.py:331-334`);
- declares a `non_linear_ratio` requirement (`camb.py:640-641`);
- and, inside its own `calculate()` (`camb.py:676`), after setting InitPower on the
  results object (`:714`) and **before** `power_spectra_from_transfer()` (`:731`), executes

      non_linear_ratio = self.provider.get_non_linear_ratio(results)   # camb.py:717
      results.Params.NonLinearModel.set_ratio(k_h, z, ratio)           # camb.py:718-722

The provider's `get_non_linear_ratio(self, results)` receives the live `CAMBdata`
**as a function argument** — plain Python method dispatch at runtime, invisible to the
dependency resolver. The boost theory therefore declares only `CAMB_transfers` as a
requirement, and the graph is a linear chain:

    camb.transfers  →  AxionHMcodeBoost  →  camb (final spectra)  →  likelihoods

Cobaya matches the provider automatically: any `get_X` method on a `Theory` subclass is
registered as providing `X` (`cobaya/theory.py:173` → `tools.py:937-948`).

Proof the pattern is sanctioned: `test_trivial_non_linear_ratio` in Cobaya's own suite
(`tests/test_cosmo_multi_theory.py:279-338`) builds exactly this sandwich. The mechanism
is in every Cobaya >= 3.6.2: the version-bump commit to 3.6.2 (`899f30a4`, which Cocoa
pins) is the first release whose last `camb.py` change is PR #480 itself (`975a9413`).
Consequently **zero Cobaya patches** were needed — plain `pip install "cobaya>=3.6.2"`
carries the identical machinery, and Cocoa's replacement `camb.yaml` retains the key.
The Fortran side needed no changes either: `ExternalNonLinearRatio`
(`fortran/ExternalNonLinearRatio.f90`, python wrapper `camb/nonlinear.py:322`)
was already part of the port.

# 3. End-to-end data flow <a name="dataflow"></a>

One likelihood evaluation:

1. `camb.transfers` runs the Boltzmann solve. For lensed Cls the wrapper uses
   `get_transfer_functions(camb_params, only_time_sources=True)` (`camb.py:1219-1221`);
   `only_time_sources` skips only the CMB l,k transfer functions and the nonlinear
   scaling — **matter transfers remain available** (`camb/camb.py:58-72`), which is what
   makes step 3 possible.
2. The main `camb` node sets the sampled primordial parameters on the results object.
3. `AxionHMcodeBoost.get_non_linear_ratio(results)` runs:
   reads the cosmology off `results.Params` (h, ombh2, omch2; axion m_ax/omaxh2 and
   regime flags from `results.Params.Axion`); reads the matter transfer functions via
   `results.get_matter_transfer_data()`; builds the linear component spectra with
   axionHMcode's own `transfer_to_PS`; loops the halo model over every redshift in
   `results.transfer_redshifts`; returns `{"k_h", "z", "ratio"}` with
   ratio = sqrt(P_NL/P_L), shape (nz, nk), z ascending.
4. The wrapper calls `set_ratio` — the Fortran `TExternalNonLinearRatio` stores the grid
   in a bilinear `TInterpGrid2D`, initializes `nonlin_ratio = 1`, and clamps requested
   (k, z) outside the grid to the boundary (`ExternalNonLinearRatio.f90:78-87`).
5. `power_spectra_from_transfer()` produces P_NL(k,z), lensed TT/TE/EE and C_L^phiphi
   with the boost applied consistently to the lensing sources.

# 4. The boost definition (the load-bearing convention) <a name="convention"></a>

axionHMcode assembles its nonlinear total in the Eq. 9 decomposition
(`halo_model/PS_nonlin_axion.py`, final assembly at the end of
`func_full_halo_model_ax`):

    P_NL = (O_db/O_m)^2 P_cold_NL + 2 (O_db O_ax / O_m^2) P_cross + (O_ax/O_m)^2 P_ax_tot

The **denominator** of the boost is the same assembly's own linear limit (one-halo terms
damped away at low k, two-halo terms → their linear inputs), which collapses to a
perfect square:

    sqrt(P_L_eq9) = (O_db/O_m) sqrt(P_cold) + (O_ax/O_m) [ fc sqrt(P_cold) + (1-fc) sqrt(P_ax) ]

with `fc = frac_cluster` (scalar per redshift, from `func_axion_param_dic`). So the
returned quantity is exactly

    sqrt(B) = sqrt(P_NL_eq9) / sqrt(P_L_eq9)

Why this denominator and not CAMB's own linear total: (i) B → 1 at low k **by
construction** (measured: within 6e-4 for dome, 2e-5 for basic, at z = 0 and 2 — the
residual is the damped one-halo (k/k_star)^4 tail, which is physical); (ii) CAMB's
P_L^tot includes massive neutrinos and the exact Boltzmann component composition, which
axionHMcode's budget does not — dividing by it would imprint a spurious few-1e-3 tilt
on the lensed spectra at low k. With the Eq. 9-internal ratio, all convention choices
live on our side of the interface and CAMB's own linear spectrum is untouched where the
model is linear. (Assessed plausible in collaborator guidance, 2026-07-01; adopted as
the working convention.)

**LCDM fallback.** The full axion assembly is singular at vanishing axion fraction
(the `frac_cluster` normalization integrates over an empty mass range). Below
`omaxh2 = 1e-8` the class switches to axionHMcode's own LCDM recipe — the cold-only
halo model `func_non_lin_PS_matter` (`halo_model/PS_nonlin_cold.py:15`, which explicitly
handles f_ax < 0.01 in its dome-alpha branch) with a tiny placeholder fraction, and
sqrt(B) = sqrt(P_cold_NL / P_cold). This is what makes the fax → 0 limit (validation V2)
well defined.

# 5. Transfer-variable semantics (verified, both regimes) <a name="transfers"></a>

Established numerically (not from documentation) by density-weighted recombination of
the component transfer functions:

| variable | DM-like axion (m/H0 >= 10) | DE-like axion (m/H0 < 10) |
|---|---|---|
| `Transfer_tot` | cdm + baryon + **axion** + massive nu (match 1.8e-6) | cdm + baryon + nu — axion **removed** (match 2.0e-6) |
| `Transfer_nonu` | cdm + baryon only — axion **not** included (3.2e-8) | cdm + baryon (3.3e-8) |
| `Transfer_axion` (=14, `camb/model.py:56`) | axion density contrast | KG-phase field output: finite, smooth, T_ax/T_cdm ~ −1e-4 at k = 0.01 falling to −6e-11 at k = 9 |

So `P_cold` for axionHMcode is the (omega_c T_cdm + omega_b T_b)/omega_cb combination
(identical to `Transfer_nonu`), and the axion enters CAMB's total exactly per the
Hlozek clustering convention the port carries.

**Transfer → P(k) identity.** axionHMcode natively consumes axionCAMB-style transfer
functions and converts with `transfer_to_PS`
(`axionCAMB_and_lin_PS/lin_power_spectrum.py:15`):
P(k) = T² (k h)⁴ P_prim(k) 2π²/k³, T in the CAMB convention (divided by k², Mpc² units),
k in h/Mpc, P in (Mpc/h)³. Feeding it CAMB's `get_matter_transfer_data()` arrays
reproduces `results.get_linear_matter_power_spectrum` to max 5.7e-8 (median 2.1e-8;
float32 transfer-storage precision), and the top-hat sigma8 from the reconstructed
P_tot matches `results.get_sigma8_0()` to four digits (0.8118 vs 0.8118). The
convention identity is exact; the same primordial parameters (As, ns, `pivot_scalar`)
are read off `results.Params.InitPower` (`camb/initialpower.py:128-136`), which the
wrapper sets *before* the ratio callback.

# 6. The redshift grid <a name="zgrid"></a>

The boost must cover every redshift CAMB will query. Under `NonLinear_lens/both` with
lensing on, `GetComputedPKRedshifts` (`fortran/results.f90:1168`) builds an internal
nonlinear-lensing grid of `nint(50 × AccuracyBoost × NonlinSourceBoost)` nodes,
**linear** in z on [0, 10] ([0, 15] when the boost product >= 2.5), merged with the
user's PK redshifts. Measured at AccuracyBoost 1: 50 nodes, [0.0, 9.8], step 0.2.

The union is exposed to Python **only** as `results.transfer_redshifts`
(`camb/results.py:226`; sort ascending before use — Fortran stores it descending).
`Params.Transfer.PK_redshifts` does **not** contain it: in the measured lensed-Cl run it
held just `[0.]`. The upstream `TrivialNonLinearRatio` test reads `PK_redshifts` and
gets away with it only because that test requests no lensed Cls; copying that pattern
here would have clamped a single z = 0 boost slice across the entire lensing kernel —
silently. The Fortran clamps (k, z) outside the supplied grid to the boundary, so the
top of the grid is where the boost has naturally decayed toward 1.

# 7. Implementation walkthrough: AxionHMcodeBoost <a name="implementation"></a>

One module, [`axionhmcode_boost/axionhmcode_boost.py`](axionhmcode_boost/axionhmcode_boost.py)
(Python 3.10, the Cocoa environment python). Structure:

- **`initialize()`** — resolves `axionhmcode_path`, pre-seeds a `scipy.misc` stub if
  scipy >= 1.14 already removed it (axionHMcode carries a dead
  `from scipy import interpolate, misc` at `halo_model/HMcode_params.py:5`; the shim
  keeps upstream files untouched), inserts the checkout into `sys.path`, imports the
  axionHMcode modules once (so fork workers inherit the compiled numba state), resolves
  the per-version call flags, and logs the dome z-calibration note once.
- **`get_requirements()`** — `{"CAMB_transfers": None}`, plus the four Dentler nuisance
  names when `sample_nuisance: True` (they then flow in as ordinary sampled parameters
  read via `provider.get_param`).
- **`get_non_linear_ratio(results)`** — the pipeline of section 3, plus the guards of
  section 11. This is deliberately the **only** public `get_*` method on the class:
  cobaya registers every `get_*` method as a providable product, so helpers are
  underscore-prefixed.
- **`_compute_row(payload)`** — one redshift of the grid, module-level so
  `multiprocessing` fork workers can run it (`processes: N` divides wall time with
  bit-identical numerics; keep 1 under MPI unless the node layout is understood).
  Per redshift it builds `cosmo_dic` exactly as axionHMcode's `load_cosmology_input`
  does — including their internal scale-independent growth
  `G_a = func_D_z_unnorm_int(z, Omega_m, Omega_w)` (`load_cosmology.py:204-205`),
  kept verbatim because it is part of the calibrated model — then calls, in order,
  `HMCode_param_dic` → `func_axion_param_dic` → `func_full_halo_model_ax`
  (or `func_non_lin_PS_matter` in LCDM mode) and forms sqrt(B) with the Eq. 9-limit
  denominator. **Only public axionHMcode entry points are used and no upstream file is
  modified** (the drag-and-drop constraint: axionHMcode updates must drop in).

Version flag mapping (from the upstream README and example notebook; overridable via
`model_flags` for V6-style experiments):

| flag | `dome` (default) | `basic` |
|---|---|---|
| `alpha` (1h/2h smoothing) | True | False |
| `concentration_param` | True | False |
| `full_2h` | False (calibration choice) | True |
| `one_halo_damping` | True | True |
| `two_halo_damping` / `eta_given` | False | False |

Axion parameters are taken from `results.Params.Axion` (the state copy):
`omaxh2` directly, or `omdah2 × axfrac` when `use_axfrac` is set; regime and switch
information from `is_de_like`, `a_osc`, `m_ovH0` (`camb/axion.py:22-52`).

Tests: [`axionhmcode_boost/tests/test_boost.py`](axionhmcode_boost/tests/test_boost.py)
(pytest; axion-run sanity, LCDM limit, DE-like hard error, strict gating — fast
Pk-grid-only configuration). The measurement scripts behind every number in this guide
are in [`axionhmcode_boost/dev_scripts/`](axionhmcode_boost/dev_scripts/).

# 8. Traps catalog <a name="traps"></a>

Each of these bit (or nearly bit) during development; all are guarded or documented.

1. **Never call the power-spectrum getters inside the provider.**
   `results.get_linear_matter_power_spectrum(..., have_power_spectra=False)` on the
   transfers-only results triggers `calc_power_spectra` (`camb/results.py:835-836`),
   which applies the nonlinear model — `ExternalNonLinearRatio` with no ratio set —
   and the Fortran hard-stops (`ExternalNonLinearRatio.f90:69-70`). Use
   `get_matter_transfer_data()` (`results.py:768`) exclusively.
2. **Never declare `Pk_grid` in the provider's requirements.** `Pk_grid` is provided by
   the final `camb` node, which itself requires `non_linear_ratio` — that would create
   the real cycle the mechanism exists to avoid.
3. **The z grid is `results.transfer_redshifts`, not `PK_redshifts`** (section 6).
4. **`WantTransfer` stays off on a Cl-only path.** A likelihood that requests only Cls
   leaves matter transfers uncomputed and `get_matter_transfer_data()` raises. The
   example yamls set `nonlinear: NonLinear_both` and `kmax: 10` in `extra_args`
   (`kmax` flips `WantTransfer` on via `set_matter_power`).
5. **`use_non_linear_ratio` is a top-level camb option**, not an `extra_args` entry
   (inside `extra_args` it reaches CAMBparams and fails).
6. **Never thin the input k grid.** The halo model's alpha/k_star machinery is
   sensitive to input-k sampling: halving the transfer k grid corrupts the boost by up
   to 11% at k ~ 0.8-7 h/Mpc, while the full grid is converged to <= 1.3e-3 against
   ×2/×4 spline densification (validation V6). Always feed the full transfer k array.
7. **cobaya >= 3.6 API changes**: `use_renames` was removed (renames are unconditional;
   the pre-boost `EXAMPLE_EVALUATE1/MCMC1.yaml` need that key deleted on pristine
   cobaya >= 3.6 — they still work on Cocoa, whose replacement camb.yaml defines the
   key); derived-parameter lambdas go directly on `derived:` (the old
   `derived: true` + `value:` combination is rejected).
8. **`Pk_grid` returns k in 1/Mpc**, not h/Mpc (`set_ratio` wants h/Mpc — the class
   passes the transfer k/h array, which is already correct).
9. **numpy >= 2 breaks this port's python layer** (ctypes scalar assignment,
   `camb/model.py:722`); **scipy >= 1.14 breaks axionHMcode's import** (dead
   `scipy.misc`; shimmed by the class). Pin numpy < 2, scipy < 1.14 — the Cocoa layers
   already do (numpy 1.26.x, scipy 1.12).
10. **axionHMcode quirks**: `M_min`/`M_max` in `cosmo_dic` are **log10 exponents**
    (the example notebook does `np.logspace(M_min, M_max, 100)`; the linear defaults in
    `load_cosmology.py` are an upstream inconsistency — do not copy them);
    `func_full_halo_model_ax` returns a tuple, element `[0]` is the total;
    `cosmo_dic["z"]` is a scalar — one halo-model evaluation per redshift;
    astropy is a real dependency (`halo_model/axion_density_profile.py`).
11. **21cm outputs are unsupported** with the external ratio: `GetNonLinRatios_All`
    error-stops, and its only call site is the 21cm power path
    (`fortran/results.f90:4153`). TT/TE/EE/phiphi and Pk_grid never reach it.

# 9. Validation battery — measured results <a name="validation"></a>

Environment: venv with cobaya 3.6.2, numpy 1.26.4, scipy 1.13.1, numba 0.60,
astropy 6.1.7 (mirroring the Cocoa pins). Full write-ups:
[`13-phase4-validation-report.md`](.claude/strategy_axionHMcode/13-phase4-validation-report.md)
(and files 11-12 for phases 0-1). Reproduction scripts: `axionhmcode_boost/dev_scripts/`.

| test | result |
|---|---|
| V0a — trivial ratio = 2 through cobaya | pk_NL = 4 × pk_L to **1.1e-15** |
| V0c — ratio ≡ 1 vs `NonLinear_none` | lensed TT/EE to ~5e-13, C_L^phiphi to ~1e-15 |
| V0b — external mead2020 ratio vs internal mead2020 (LCDM) | dTT <= 2e-7 (l <= 2000), dCpp <= **3.2e-5** (L <= 1000) |
| V1 — transfer→P(k) identity | max 5.7e-8; sigma8 0.8118 = 0.8118 |
| V2 — LCDM limit vs CAMB HMcode-2020, z=0 | R(k=1): mead 5.69, basic 5.10 (−10%), dome 6.40 (+12%) — the papers' own claimed agreement level |
| V3 — Gaughan/Green/Moss reproduction | see below |
| V4 — vs this branch's axion-aware internal HMcode (fax_dm = 0.1, m = 1e-25) | dTT <= 0.22%; dCpp +0.6/+10/+22% at L = 100/500/1000 — two different halo models, expected band |
| V6 — numerical convergence (m = 1e-24, fax = 0.3, dome) | M-grid 100→320: max dB/B 7.6e-4; k ×2/×4 densification: <= 1.3e-3; k halving: up to **11%** (hence trap 6) |
| V5 — full-likelihood evaluate/MCMC | **pending** likelihood data (`cobaya-install`); wiring smoke-tested end-to-end at 83.6 s/eval |

**V3 in detail** (Planck-2018 fiducial, fax = Ω_ax/Ω_D = 0.3, the setup of
[arXiv:2605.12054](https://arxiv.org/abs/2605.12054)). Boost-ratio curves
P_NL^axion/P_NL^LCDM at z = {0, 2} for m = 1e-23/-24/-25 reproduce their Fig. 1
structure: dome boosted above unity, basic suppressed below, suppression deepening
toward lower mass and higher z — e.g. basic m = 1e-24, z = 2 gives R(k=1) = 0.47 vs
their ~0.5. Lensed-spectra differences vs LCDM (m = 1e-24): dome dTT oscillatory,
+0.97% at l = 2400 (their Fig. 2 red reaches ~+2-3% by l ~ 3000); dome dCpp
+3.0/+17.4/+32.6/+46.8% at L = 100/500/1000/2000 vs their Fig. 3 read-offs
~+4/+13/+26/+38; basic dCpp −2.6/−7.4/−17.2% at L = 500/1000/2000 vs their
~−4/−10/−19. Overlay figures (our curves vs approximate read-offs from their published
figures): `axionhmcode_boost/dev_scripts/gaughan_comparison_fig1.png` and
`gaughan_comparison_fig23.png` (regenerate with `gaughan_comparison_plot.py`). Residual
dome offsets at the highest L and at m = 1e-23 (1.5 decades above the dome calibration
pivot) are consistent with genuine inter-code differences: their linear inputs are
AxiCAMB (standard-EFA lineage), ours AxiECAMB (Passaglia-Hu EFA), and the two differ
most in exactly this mass range.

# 10. Performance <a name="performance"></a>

Measured warm (numba JIT paid; first call adds ~7 s per process): **1.3 s/redshift
(basic), 2.7 s/redshift (dome)**, nk ~ 220. A lensed-Cl evaluation needs the full
lensing grid — 50 nodes at AccuracyBoost 1 (75 at 1.5) — giving ~65-135 s per
likelihood evaluation single-threaded; the smoke-tested end-to-end figure was 83.6 s
(dome, 50 nodes). `processes: N` forks the z loop (identical numerics, wall time / N).
Accuracy was prioritized by design decision (see the decision log): no grid thinning
is permitted on cost grounds; if production MCMC cost becomes prohibitive, the intended
path is training an ML emulator on this pipeline's output, dropping into the same
`get_non_linear_ratio` interface.

# 11. Validity domain and known limitations <a name="validity"></a>

Hard errors, regardless of the `strict` option:

- **DE-like axions** (m/H0 < 10, `results.Params.Axion.is_de_like`): the mixed-DM halo
  model is undefined (the axion is not in the matter budget; no Jeans-scale
  clustered/unclustered split). Use `halofit_version: original` for those masses.
  Collaborator guidance (2026-07-01) notes this restriction is pragmatic, not
  fundamental — a smooth-vs-clustered decomposition (DE-like axion as w0wa-like smooth
  DE, or as a smooth mixed-DM component) could lift it once the "matter" conventions
  and their Weyl-potential mapping are pinned down; deferred pending further input.
- **z grids reaching 20% of z_osc** (the KG→EFA switch): above the switch the axion
  transfer function is a field-phase quantity, not a halo-model density contrast.
  Never triggered in the target window (z_osc >= 2.6e4 for m >= 1e-25 eV vs grid
  max ~10); the guard exists so a light-mass misconfiguration fails loudly.

Gated by `strict` (False = warn once + extrapolate, the 2605.12054 practice;
True = hard error): axion fraction beyond the version's calibration (dome:
f_ax = Ω_ax/Ω_m <= 0.3; basic: f_ax,dm <= 0.5) and dome masses outside
[1e-25.5, 1e-23.5] eV (calibration pivot 1e-24.5). The lensing z grid inherently
exceeds dome's 1 < z < 8 calibration — logged once at initialization, not gated.

Modeling limitations (documented, accepted): power-law primordial spectra only
(axionHMcode's internal P_prim); massive neutrinos outside axionHMcode's budget
(handled by the Eq. 9-internal ratio convention, section 4); axionHMcode's internal
growth `G_a` and Ω_w = 1 − Ω_m assume a flat ΛCDM background without the axion's
early-DE phase (part of the calibrated model, kept verbatim). Target mass window:
m_ax ~ 1e-25..1e-23 eV.

# 12. Cocoa integration <a name="cocoa"></a>

- `installation_scripts/setup_axie_camb.sh` clones this repository (pinned by
  `AXIE_CAMB_GIT_COMMIT`) into `external_modules/code/axiecamb` (canonical name,
  lowercase) with compiler patches from `cocoa_installation_libraries/axiecamb_changes/`,
  and clones upstream axionHMcode (pinned by `AXION_HMCODE_GIT_COMMIT`) into
  `external_modules/code/axionHMcode`. `compile_axie_camb.sh` builds AxiECAMB
  (`setup.py build` with the RECOMBINATION_FILES variants). Both are gated by
  `INSTALL_AXIE_CAMB_V2` in `set_installation_options.sh`.
- **axionHMcode needs no compile step**: pure Python, no packaging, dependencies
  already present in Cocoa's layered environment (conda `cocoapy310` provides
  numba 0.60 and astropy; the `.local` pip prefix overlays numpy 1.26.x). The only
  compilation is the per-process numba JIT; upstream sets no `cache=True`, and the
  drag-and-drop constraint forbids adding it.
- Cobaya itself is pinned by Cocoa at `COBAYA_GIT_COMMIT = 899f30a4` — exactly the
  3.6.2 version bump verified to contain the full mechanism. If that pin is ever moved,
  re-verify the `use_non_linear_ratio` code paths cited in section 2.

# 13. Decision log <a name="decisions"></a>

Recorded fully in [`.claude/strategy_axionHMcode/`](.claude/strategy_axionHMcode/00-INDEX.md)
(files 06 and 10); the operative outcomes:

- **`version: dome` is the default** (most recent recalibration), with the Dentler
  nuisance parameters (alpha_1, alpha_2, gamma_1, gamma_2; arXiv:2111.01199) exposed —
  fixed values via class options or sampled via `sample_nuisance: True`.
- **Out-of-calibration policy = the `strict` yaml flag** (section 11); DE-like and
  switch-redshift conditions hard-error unconditionally.
- **Mass modes**: evaluate example carries the sampled-log-mass block (prior
  logmx ∈ [−25, −23]) with the mass pinned in the `evaluate` override; the MCMC example
  fixes the mass (one chain per mass, the practice of arXiv:2301.08361 and 2605.12054)
  with the sampled-mass alternative in comments.
- **Accuracy first, performance never drives design**: the boost is evaluated on the
  full CAMB grid; any thinning must pass the V6 convergence test; the cost escape hatch
  is a future emulator trained on this pipeline (same interface).
- **Boost convention** = Eq. 9-internal ratio (section 4), assessed plausible in
  collaborator guidance; it is also what makes drag-and-drop axionHMcode updates
  possible (all conventions on our side of the interface).
- **Drag-and-drop constraint**: axionHMcode is called through public entry points only,
  never modified; the pinned upstream checkout can be bumped independently.

# 14. File inventory <a name="files"></a>

Added by this project (all in this repository):

| path | content |
|---|---|
| `axionhmcode_boost/axionhmcode_boost.py` | the `AxionHMcodeBoost` Theory class |
| `axionhmcode_boost/README.md` | user-level options, yaml recipes, cost table |
| `axionhmcode_boost/tests/test_boost.py` | pytest suite (4 tests) |
| `axionhmcode_boost/dev_scripts/` | the 10 measurement/validation scripts + 2 comparison figures behind sections 5-6 and 9 |
| `EXAMPLE_AXIONHMCODE_EVALUATE1.yaml` | single-point evaluation (mode A: log-mass prior + override) |
| `EXAMPLE_AXIONHMCODE_MCMC1.yaml` | MCMC (mode B: fixed mass; sampled-mass blocks in comments) |
| `.claude/strategy_axionHMcode/00-13` | the raw working documents: strategy, verified facts, decision log, phase results, validation report |
| this file | the consolidated developer guide |

Nothing was modified in: Cobaya (zero patches), the AxiECAMB Fortran, or axionHMcode.
