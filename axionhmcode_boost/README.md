# axionhmcode_boost

Cobaya `Theory` class feeding the axionHMcode mixed-dark-matter nonlinear boost
B(k,z) = P_NL/P_L to AxiECAMB through cobaya's `use_non_linear_ratio` mechanism
(cobaya >= 3.6.2). The Boltzmann code runs once per likelihood evaluation; the
boost is computed from its linear transfer functions at every redshift CAMB
uses (including the internal ~50-node nonlinear-lensing grid over 0 <= z <~ 10)
and applied consistently to P(k), lensed TT/TE/EE, and C_L^phiphi.

Physics: Vogt et al. (arXiv:2209.13445, `version: basic`) and the Dome et al.
recalibration (arXiv:2409.11469, `version: dome`, default). Validation targets
and measured results live in `../.claude/strategy_axionHMcode/` (files 07, 11, 12).

## Requirements

- The AxiECAMB port as the CAMB module (this repository; needs the
  `Transfer_axion` variable and `ExternalNonLinearRatio`).
- cobaya >= 3.6.2 (first release with `use_non_linear_ratio`, PR #480).
- An unmodified axionHMcode checkout (github.com/SophieMLV/axionHMcode).
  It is called only through public entry points, so upstream updates drop in.
- numpy < 2 (the CAMB-1.6.7 python layer is not numpy-2 compatible), scipy
  < 1.14 (axionHMcode carries a dead `scipy.misc` import; this class shims it
  if scipy already removed it), numba, astropy.

## Minimal yaml

```yaml
theory:
  camb:
    path: ./external_modules/code/axiecamb
    use_non_linear_ratio: True      # top-level camb option, NOT extra_args
    extra_args:
      num_massive_neutrinos: 1
      nnu: 3.046
      nonlinear: NonLinear_both
      lens_potential_accuracy: 4
  axionhmcode_boost.AxionHMcodeBoost:
    python_path: ./external_modules/code/axiecamb/axionhmcode_boost
    axionhmcode_path: ./external_modules/code/axionHMcode
    version: dome                   # 'dome' (default) | 'basic'
    strict: False                   # True: hard-error outside calibration
```

## Options

| option | default | meaning |
|---|---|---|
| `axionhmcode_path` | (required) | path to the axionHMcode checkout |
| `version` | `dome` | calibration: `dome` (Dome et al.) or `basic` (Vogt et al.) |
| `strict` | `False` | out-of-calibration fax/mass: `False` = warn once + extrapolate (Gaughan et al. practice), `True` = hard error |
| `sample_nuisance` | `False` | expose `alpha_1 alpha_2 gamma_1 gamma_2` (Dentler et al. 2111.01199) as sampled input parameters |
| `alpha_1 ... gamma_2` | `None` | fixed nuisance values when not sampling (None = axionHMcode internal defaults) |
| `m_grid_points` | 100 | halo-mass integration grid size |
| `m_min_exponent`, `m_max_exponent` | 7, 18 | log10 of the mass grid bounds (Msun/h) |
| `model_flags` | `None` | expert overrides of the per-version axionHMcode call flags |
| `accuracy_boost` | 1 | CAMB-style multiplier on the internal discretizations: scales `m_grid_points` and, when `axionhmcode_path` is the SBU-COSMOLIKE fork, its growth/G/sigma(M) table node counts. Rerun the same yaml at 1 and 2 and compare boost grids to check convergence (measured 1-vs-2 shift: max dB/B ~ 1e-3, dominated by the halo-mass grid) |

Parallelism (deliberately not a yaml option): the redshift loop forks one
worker per `OMP_NUM_THREADS` core (identical numerics, wall time only;
falls back to 1 when unset). The environment's core budget is always
respected — the same contract as every OpenMP code in Cocoa, and correct
under MPI, where each rank forks within its own allocation. There is no
yaml override.

Nuisance sampling: with `sample_nuisance: True`, add the four parameters to the
yaml `params` block (Gaughan et al. flat priors: alpha_1 [0.6, 2.0], alpha_2
[1.43, 2.54], gamma_1 [5.0, 45.0], gamma_2 [-0.37, -0.23]).

## The dome flag configuration: released code vs. calibration paper

For `version: dome` this theory calls axionHMcode with `alpha=True,
concentration_param=True, one_halo_damping=True, full_2h=False,
eta_given=False, two_halo_damping=False`. These defaults follow the released
code: the upstream README recommends exactly `alpha = True` and
`concentration_param = True` for the dome version, notes that the calibration
"was performed with `full_2h = False`", and warns that the HMcode-2020
parameter switches are not calibrated to mixed axion/cold-dark-matter
cosmologies. This configuration also reproduces Figs. 1-3 of Gaughan, Green &
Moss ([arXiv:2605.12054](https://arxiv.org/abs/2605.12054)), the published
constraint pipeline built on this code (validation test V3 in the
[developer guide](../BOOST_DEVELOPER_GUIDE.md)).

There is, however, a documented tension with the calibration paper itself.
Dome et al. ([arXiv:2409.11469](https://arxiv.org/abs/2409.11469), Sec. 4.4)
state that the recalibrated model activates by default the halo bloating
parameter eta (their Eq. 57) and the perturbative two-halo damping (their
Eq. 55), alongside the one-halo damping and the alpha smoothing — meaning the
fitted alpha values of their Table 5 were calibrated with eta and the two-halo
damping active. Running those fitted alpha values with the two switches off,
as the released code recommends, is therefore a small internal inconsistency
in the upstream recommendation. The expected size of the difference is a few
percent in the quasi-linear regime, inside the model's quoted 10-20% accuracy
band.

> [!NOTE]
> Our defaults deliberately follow the released code — that keeps this
> pipeline consistent with how the ecosystem actually runs axionHMcode and
> with the published Gaughan et al. results. The question has been put to the
> collaborators/upstream. If it is resolved in favor of the paper
> configuration, no code change is needed; flip the two switches from the
> yaml:
>
> ```yaml
> axionhmcode_boost.AxionHMcodeBoost:
>   # ...
>   model_flags: {eta_given: True, two_halo_damping: True}
> ```

## Axion mass: the two run modes

`logmx` = log10(m_ax/eV); target window for the boost pipeline is
m_ax ~ 1e-25 .. 1e-23 eV (arXiv:2605.12054). See `EXAMPLE_AXIONHMCODE_*.yaml`
in the repository root for complete files.

Mode A — evaluate (single cosmology): keep the sampled-log-mass parameter block
and pin the mass in the evaluate sampler override:

```yaml
params:
  logmx:
    prior: {min: -25.0, max: -23.0}
    ref: {dist: norm, loc: -24.0, scale: 0.1}
    proposal: 0.1
    latex: \log{m_{\rm ax}}
    drop: true
  m_ax:
    value: 'lambda logmx: 10**logmx'
    latex: m_\mathrm{ax}
sampler:
  evaluate:
    override:
      logmx: -24.0
```

Mode B — MCMC (the committed example): fix the mass, one chain per mass (the
Gaughan/Rogers practice; the mass posterior is poorly constrained). To sample
the mass instead, replace the fixed `m_ax` with the two blocks from Mode A.

```yaml
params:
  m_ax:
    value: 1.0e-24
    latex: m_\mathrm{ax}
```

## Validity domain and hard limits

- DE-like axions (m/H0 < 10, `results.Params.Axion.is_de_like`) hard-error
  regardless of `strict`: the mixed-dark-matter halo model is undefined there.
  Use the halofit path (existing EXAMPLE yamls) for those masses.
- A z grid reaching 20% of z_osc (the KG->EFA switch) hard-errors regardless
  of `strict` — never triggered in the target window (z_osc >~ 1e4).
- dome calibration: fax = O_ax/O_m in [0.01, 0.3], m near 1e-24.5 eV,
  1 < z < 8. fax/mass violations follow `strict`; the z range is inherently
  exceeded by the lensing grid and is logged once as a documented limitation.
- Power-law primordial spectra only (axionHMcode's internal P_prim assumption).
- Massive neutrinos are outside axionHMcode's matter budget. The boost is
  defined inside the model's own Eq. 9 decomposition (numerator and
  denominator), so the convention cancels at linear order and the boost -> 1
  at low k by construction; CAMB's own linear total (which includes nu) is
  untouched there.

## Cost

Measured single core, warm numba, per redshift:

| checkout | basic | dome |
|---|---|---|
| upstream SophieMLV | ~1.3 s | ~2.7 s (3.5 s at z=0) |
| SBU-COSMOLIKE fork (VM-SPEEDUP rounds 1-3, 2026-07-02/03) | ~0.2 s | ~0.34 s |

The fork changes are validated against upstream to max |dB/B| = 1.6e-5
(see the fork README appendix for the mathematics and the protocol).
Lensed-Cl runs need the full internal lensing grid (50 nodes at
AccuracyBoost 1, 75 at 1.5): ~17 s (fork) / ~135 s (upstream) per dome
likelihood evaluation single-threaded. The z-loop fork parallelism divides
the wall time without changing any number, with one worker per
`OMP_NUM_THREADS` core (always; no yaml override). If this is still too
slow for production MCMC, the intended path is training an ML emulator on
this pipeline's output (same interface).

## Tests and development scripts

- `tests/test_boost.py` — pytest suite (axion run sanity, LCDM limit, DE-like
  hard error, strict gating). Set `AXIONHMCODE_PATH` if axionHMcode is not at
  `../axionHMcode` relative to the repository parent.
- `dev_scripts/` — the Phase 0/1 measurement scripts behind strategy files
  11-12 (trivial-ratio null test, timing benchmark, regime semantics, and the
  standalone lensed-Cl prototype).
