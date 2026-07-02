# AxiECAMB (modern-CAMB port)

Ultralight-axion (ULA) physics on modern CAMB 1.6.7, with a Fortran driver, a Python wrapper, a Cobaya interface, and an axionHMcode nonlinear-boost add-on.

# Table of contents
1. [Overview](#overview)
2. [The physics in brief](#the-physics-in-brief)
3. [Installation and compilation](#installation-and-compilation)
4. [Running AxiECAMB](#running-axiecamb)
    1. [Fortran driver and .ini files](#fortran-driver-and-ini-files)
    2. [Python wrapper](#python-wrapper)
    3. [Cobaya interface (axion MCMC)](#cobaya-interface-axion-mcmc)
    4. [axionHMcode nonlinear boost (Cobaya theory block)](#axionhmcode-nonlinear-boost-cobaya-theory-block)
5. [Warnings and known differences](#warnings-and-known-differences)

### Appendices (port internals)

6. [Appendix A: What was ported and where](#appendix-a-what-was-ported-and-where)
7. [Appendix B: Validation against the original AxiECAMB](#appendix-b-validation-against-the-original-axiecamb)
8. [Appendix C: the axionHMcode boost — implementation details](#appendix-c-the-axionhmcode-boost--implementation-details)
9. [Appendix D: About the underlying CAMB](#appendix-d-about-the-underlying-camb)
10. [Port developer guide (separate file)](PORT_DEVELOPER_GUIDE.rst)

# Overview <a name="overview"></a>

**AxiECAMB** is the ultralight-axion (ULA) effective method of [arXiv:2412.15192](https://arxiv.org/abs/2412.15192) (Liu, Hu et al.), ported from its original CAMB-Nov13 base onto **modern CAMB 1.6.7** (Fortran + Python wrapper). It computes the CMB, lensing, matter power spectra and background evolution for cosmologies with a ULA component, whose behaviour is set by the axion mass `m_ax`: dark-matter-like for `m/H0 >= 10` and dark-energy-like below that.

This repository provides three ways to run the code and one add-on:

- a **Fortran/.ini** driver and a **Python** wrapper (both below);
- a **Cobaya** interface for MCMC that works with *unmodified* Cobaya;
- an **axionHMcode nonlinear-boost** Cobaya theory block (`axionhmcode_boost/`) that feeds a mixed-dark-matter halo-model boost B(k,z) = P_NL/P_L into the lensed CMB and lensing-potential spectra.

New users should read [The physics in brief](#the-physics-in-brief) and [Running AxiECAMB](#running-axiecamb); the appendices document the port internals (what was changed, validation), and the full [developer guide](PORT_DEVELOPER_GUIDE.rst) lives in a separate file.

> [!NOTE]
> When using the axion module, please cite [arXiv:2412.15192](https://arxiv.org/abs/2412.15192) and Passaglia & Hu 2022 ([arXiv:2201.10238](https://arxiv.org/abs/2201.10238)), on which the effective method builds. The original AxiECAMB heavily modified [axionCAMB](https://github.com/dgrin1/axionCAMB) (Hlozek et al., arXiv:1410.2896). For the axionHMcode boost, also cite Vogt et al. ([arXiv:2209.13445](https://arxiv.org/abs/2209.13445)) and Dome et al. ([arXiv:2409.11469](https://arxiv.org/abs/2409.11469)).

# The physics in brief <a name="the-physics-in-brief"></a>

The axion has a quadratic potential and is evolved in synchronous gauge:

- **Background**: the exact Klein-Gordon (KG) equation is solved (16-stage 8th-order fixed-step Runge-Kutta in ln a) from deep radiation domination until m = dfac*H, with the initial field value found by shooting to match the requested relic abundance. At the switch the field is projected onto WKB cos/sin amplitudes and matched onto an **effective fluid** (EFA) whose density follows a^-3 with an exp[3 int w dln a] residual, w(a) = wEFA_c (H/m)^2; the matching coefficients (\<H\>, wEFA_c) are iterated to self-consistency. dfac (default 10) is retuned internally: the oscillation phase at the switch is targeted to 2 beta = 7.08 pi for light DM-like axions, and the switch is pushed out of the recombination window z in (800, 1300).
- **Perturbations**: exact KG before the switch (with a per-k conditioning rescale of delta-phi), the (1+w)-weighted GDM effective fluid after, with sound speed cs^2 = (sqrt(1+kappa)-1)^2/kappa + (5/4)(H/am)^2, kappa = k^2/(a^2 m^2). At the switch the KG variables are projected onto (delta_ax, u_ax) preserving velocity and shear continuity; the residual metric jump is absorbed into eta (sub-horizon) or carried as a delta-function boundary term in the temperature line-of-sight integral (super-horizon).
- **DM vs DE**: for m/H0 >= 10 the axion is dark-matter-like — it counts in the matter transfer functions, sigma_8, the equality redshift, CosmoMC theta, and halofit Omega_m. For m/H0 < 10 (m <~ 1.4e-32 eV) it is dark-energy-like: KG is solved to a = 1, there is no fluid switch, and the matter transfer excludes the axion.

# Installation and compilation <a name="installation-and-compilation"></a>

You need `gfortran` (and a working C compiler). Build `forutils` first if your tree does not ship it prebuilt.

**Fortran executable** (run in `fortran/`):

    make camb

**Python library** — build `camblib.so` into `camb/` (run in `fortran/`):

    make python PYCAMB_OUTPUT_DIR=../camb/

(equivalently, `python setup.py make` from the repository root).

> [!NOTE]
> Inside [Cocoa](https://github.com/CosmoLike/cocoa), installation and compilation are automated by `installation_scripts/setup_axie_camb.sh` (clones this repo, pinned by commit, into `external_modules/code/axiecamb` and applies the compiler patches) and `compile_axie_camb.sh` (the `setup.py build`). The axionHMcode checkout is cloned by the same setup script and needs no compilation (pure Python). Both are gated by the `INSTALL_AXIE_CAMB_V2` flag in `set_installation_options.sh`.

# Running AxiECAMB <a name="running-axiecamb"></a>

Four entry points, in increasing level of automation: the Fortran driver, the Python wrapper, the Cobaya interface for axion MCMC runs, and the axionHMcode nonlinear-boost theory block.

## Fortran driver and .ini files <a name="fortran-driver-and-ini-files"></a>

Build and run (in `fortran/`; build forutils first if needed):

    make camb
    ./camb ../inifiles/params_axion.ini

New ini keys (see `inifiles/params_axion.ini`):

| key | meaning |
|-----|---------|
| `m_ax` | ULA mass in eV (negative input = log10(m_ax/eV)) |
| `use_axfrac` | `T`: use (`omdah2`, `axfrac`); `F`: use `omaxh2` (+ usual `omch2`) |
| `omaxh2` | Omega_ax h^2 (when `use_axfrac = F`) |
| `omdah2` | total dark-matter Omega h^2 (when `use_axfrac = T`) |
| `axfrac` | axion fraction of DM (m/H0 >= 10) or of DE (m/H0 < 10) |
| `axion_dfac` | switch threshold m = dfac*H (default 10; retuned internally) |
| `axion_isocurvature`, `Hinf` | accepted but isocurvature is force-disabled (v1.0 parity) |

With `use_axfrac = T`, `omch2` may be omitted (it is derived). Constant-w dark energy (fluid or PPF) can be combined with the axion; quintessence dark-energy models cannot (the axion background solver treats DE as Lambda, as in the original).

## Python wrapper <a name="python-wrapper"></a>

```python
import camb
pars = camb.CAMBparams()
pars.set_cosmology(H0=67.32, ombh2=0.02238, omch2=0.108, mnu=0.06, tau=0.054)
pars.InitPower.set_params(As=2.1e-9, ns=0.966)
pars.set_axion(m_ax=1e-27, omaxh2=0.012)                  # or omdah2=..., axfrac=...
results = camb.get_results(pars)
Ax = results.Params.Axion        # derived quantities live on the *result* state copy
print(Ax.a_osc, Ax.dfac_used, Ax.tau_osc, Ax.m_ovH0)
```

The axion density perturbation is available as the `delta_axion` matter transfer column (`camb.model.Transfer_axion`). `delta_tot` (and hence sigma_8 and the default matter power) includes the axion when it is DM-like, excludes it when DE-like.

To build the Python library run `make python` in `fortran/` (or use `python setup.py make`), which places `camblib.so` in `camb/`.

## Cobaya interface (axion MCMC) <a name="cobaya-interface-axion-mcmc"></a>

This code works with **unmodified (pristine) Cobaya** — no patching of Cobaya's CAMB theory wrapper is needed. The axion parameters (`m_ax`, `omaxh2`, `omdah2`, `axfrac`, `dfac`) are discovered and set through this package's own `camb.get_valid_numerical_params` and `camb.set_params` hooks, which the stock wrapper already uses. Two ready-to-run examples ship in the repository root:

- `EXAMPLE_EVALUATE1.yaml` — single-point posterior evaluation;
- `EXAMPLE_MCMC1.yaml` — the corresponding MCMC (Planck lite + lowl TT/EE + DESI DR2 BAO + DES-Y5 SN + ACT DR6 lensing), sampling (logA, ns, 100theta_*, omegabh2, omegach2, tau, omegaaxh2, logmx).

Setup and run (from the repository root, with `cobaya` installed):

    cd fortran && make camb && make python PYCAMB_OUTPUT_DIR=../camb/ && cd ..
    pip install act_dr6_lenslike
    cobaya-install EXAMPLE_MCMC1.yaml -p /path/to/cobaya_packages
    cobaya-run EXAMPLE_MCMC1.yaml -p /path/to/cobaya_packages

Conventions used in the examples (see comments inside the yaml files): the chain samples `thetastar100` (= 100 theta_*) and feeds `thetastar` to CAMB via a value-lambda (CAMB's input is theta_* itself); `logA -> As`, `omegaaxh2 -> omaxh2` and `logmx -> m_ax` are standard value-lambda mappings; `theta_H0_range: [40, 130]` brackets the theta -> H0 solution over the whole prior box; `halofit_version: original` is the AxiECAMB-validated non-linear treatment; the `omegam` derived parameter excludes the axion (the DE-like convention).

## axionHMcode nonlinear boost (Cobaya theory block) <a name="axionhmcode-nonlinear-boost-cobaya-theory-block"></a>

The `axionhmcode_boost/` folder ships a Cobaya `Theory` class, `AxionHMcodeBoost`, that supplies the mixed-dark-matter nonlinear boost B(k,z) = P_NL/P_L from [axionHMcode](https://github.com/SophieMLV/axionHMcode) to AxiECAMB, so lensed TT/TE/EE and the lensing potential C_L^phiphi use a Jeans-scale-aware nonlinear prescription instead of halofit/HMcode. It runs through Cobaya's `use_non_linear_ratio` mechanism (cobaya >= 3.6.2) with **no Cobaya patches and no Fortran changes**.

Two ready-to-run examples ship in the repository root:

- `EXAMPLE_AXIONHMCODE_EVALUATE1.yaml` — single-point evaluation (log-mass prior, mass pinned in the `evaluate` sampler);
- `EXAMPLE_AXIONHMCODE_MCMC1.yaml` — the MCMC (fixed mass; comments show how to switch to a sampled log-mass).

Both target the mass window `m_ax` ~ 1e-25..1e-23 eV, where the axion Jeans scale sits in the quasi-linear regime probed by CMB lensing ([arXiv:2605.12054](https://arxiv.org/abs/2605.12054)). Setup and run:

    git clone https://github.com/SophieMLV/axionHMcode ../axionHMcode
    cobaya-install EXAMPLE_AXIONHMCODE_EVALUATE1.yaml -p /path/to/packages
    cobaya-run     EXAMPLE_AXIONHMCODE_EVALUATE1.yaml -p /path/to/packages

> [!TIP]
> Options (`version: dome|basic`, the `strict` validity flag, nuisance-parameter sampling, fork parallelism, ...) are documented in [`axionhmcode_boost/README.md`](axionhmcode_boost/README.md). The design rationale, conventions and measured validation live in `.claude/strategy_axionHMcode/` (start at `00-INDEX.md`); the implementation details are summarized in [Appendix C](#appendix-c-the-axionhmcode-boost--implementation-details).

# Warnings and known differences <a name="warnings-and-known-differences"></a>

> [!WARNING]
> Keep the following in mind when running:

- **Isocurvature is disabled** (as in AxiECAMB v1.0; the original mode-6 vector targets variables that are not evolved). Inputs are accepted and ignored with a warning.
- The **growth-rate (Transfer_f) column** of the original is disabled there and was not ported; modern CAMB's own growth outputs are available.
- The **non-linear mode** is inherited from axionCAMB and not extensively tested. Two supported paths: `halofit_version = 1` (original; the AxiECAMB-validated treatment — Takahashi was found unstable for axion models), and HMcode (`mead2020` etc.), which has been made axion-consistent: the exact axion background (KG/EFA density and equation of state) enters HMcode's internal expansion and growth as a separate component, and DM-like axions (m/H0 >= 10) count as fully clustering cold matter in the halo-model quantities (Omega_m, the sigma(R)-mass mapping, the EH99 cold ratio, f_nu and the Dolag reference) — i.e. no Jeans-scale halo suppression is modelled (the linear input P(k) carries the scale-dependence). DE-like axions affect HMcode only through the expansion. Non-axion HMcode results are bit-identical to upstream (all changes gated).
- For z > 0 transfer outputs, mind whether the requested z is before or after the switch: the axion density contrast is defined differently in the two regimes.
- P(k) at wavenumbers where the spectrum is suppressed by >~ 6 orders of magnitude differs from the original (which zeroed rather than extrapolated the dead tail); set `transfer_kmax` high enough for any application sensitive to that region.
- Default accuracy (`accuracy_boost = 1`) is what was validated in arXiv:2412.15192; higher boosts apply to the non-ULA accuracy settings only.
- CosmoMC theta counts the axion in omega_dm unconditionally (original behaviour), which is only meaningful for DM-like axions.

---

# Appendix A: What was ported and where <a name="appendix-a-what-was-ported-and-where"></a>

All modifications in the Fortran sources are marked with inline `!AxiECAMB` comments (`grep -n AxiECAMB fortran/*.f90` lists every change site).

- `fortran/AxionBackground.f90` (new): the KG background solver `w_evolve`, EFA matching `auxiIC`, phase targeting and recombination-skip dfac retuning (moved here from the old `inidriver_axion.F90` so the Python interface gets them too), as the component class `TAxionModel` stored in `CAMBparams%Axion`.
- `results.f90`: density budget/closure with the axion, the solver invocation, tau_osc, background integrals split at the dtauda kink at a_osc (applied uniformly: times, distances, sound horizons, optical depths — the original only split some), fine time-step window around tau_osc, thermo values cached at tau_osc, `Transfer_axion` column, z_eq and CosmoMC theta definitions.
- `equations.f90`: the two axion perturbation equations (KG <-> EFA), the mid-evolution switch in the `next_switch` chain with the WKB projection (`AxionSwitchKGtoEFA`), adiabatic delta-phi initial conditions, axion terms in dgrho/dgq/grho/gpres, low-k lmaxnr boost ("WH smoother"). Tensors need no axion terms: the modern tensor background comes from the dtauda-based thermo table (this also fixes an original-code issue where tensors extrapolated the field table past a_osc).
- `cmbmain.f90`: the switch boundary term in the temperature LOS integral (`deltaBCSrc` machinery, flat and curved cases), axion-aware integration start time.
- `recfast.f90`: dHdz in the tightly-coupled T_mat term includes the axion (numerical derivative of the exact H(z), stepped away from the a_osc kink).
- `halofit.f90`: axion counted in Omega_m (DM-like) or in the smooth DE (DE-like), with a warning that the non-linear mode is inherited from axionCAMB and not well tested.
- Python: `camb/axion.py` (`AxionModel`), `CAMBparams.set_axion`, transfer-name lists.

# Appendix B: Validation against the original AxiECAMB <a name="appendix-b-validation-against-the-original-axiecamb"></a>

With matched cosmologies, the axion/LCDM suppression ratios agree between the original AxiECAMB (Nov13 base) and this port to:

| case | TT C_l ratio (l=2-2600) | P(k) ratio |
|------|-------------------------|------------|
| m=1e-27 eV, 10% of DM (switch z~1341) | <= 0.01% | <= 0.005% |
| m=1e-27 eV, 100% of DM (`use_axfrac`) | <= 0.10% | <= 0.005% (where suppression < 10^3) |
| m=1e-30 eV, 10% of DM (switch z~24, boundary term active) | <= 0.08% | <= 0.01% |
| m=1.4e-33 eV (DE-like, no switch) | <= 0.09% | <= 0.01% |

Residuals at the 0.03-0.1% level are dominated by Nov13 <-> CAMB-1.6.7 baseline physics differences that do not perfectly cancel in the ratios (the absolute LCDM baselines differ by ~0.2%). The pure-LCDM limit of this code is bit-identical to unmodified CAMB 1.6.7. The standard CAMB Python test suite passes.

# Appendix C: the axionHMcode boost — implementation details <a name="appendix-c-the-axionhmcode-boost--implementation-details"></a>

This documents the second project phase built on the port: feeding the **axionHMcode** mixed-dark-matter nonlinear boost B(k,z) = P_NL/P_L (Vogt et al., [arXiv:2209.13445](https://arxiv.org/abs/2209.13445); Dome et al. recalibration, [arXiv:2409.11469](https://arxiv.org/abs/2409.11469)) into AxiECAMB through Cobaya, so lensed TT/TE/EE and the lensing potential C_L^phiphi are computed with a Jeans-scale-aware nonlinear prescription (science context: [arXiv:2605.12054](https://arxiv.org/abs/2605.12054)). Full design notes, decisions and measured validation numbers live in `.claude/strategy_axionHMcode/` (start at `00-INDEX.md`); the implementation is in [`axionhmcode_boost/`](axionhmcode_boost/README.md), with example inputs `EXAMPLE_AXIONHMCODE_EVALUATE1.yaml` and `EXAMPLE_AXIONHMCODE_MCMC1.yaml`.

### Architecture

Cobaya >= 3.6.2 splits the `camb` block into two graph nodes (`camb.transfers` runs the Boltzmann solve once; `camb` assembles spectra) and, with `use_non_linear_ratio: True`, calls `provider.get_non_linear_ratio(results)` from inside its own `calculate()`, passing the transfers-level CAMBdata as an argument. The `AxionHMcodeBoost` Theory class (`axionhmcode_boost/`) therefore requires only `CAMB_transfers` — no dependency cycle exists — and extracts the linear transfer functions itself via `get_matter_transfer_data()` (never the power-spectrum getters, which would invoke the not-yet-configured nonlinear model). The boost is evaluated at every redshift CAMB uses, read from `results.transfer_redshifts` (under nonlinear lensing this includes the internal grid of nint(50 x AccuracyBoost x NonlinSourceBoost) nodes, linear in z on [0, 10]; `Params.Transfer.PK_redshifts` does NOT contain it). No Cobaya patches are needed; Fortran side unchanged (`ExternalNonLinearRatio` was already part of the port).

### Boost convention

Numerator and denominator both live in axionHMcode's Eq. 9 decomposition. The denominator is the model's own linear limit, a perfect square:

    sqrt(P_L_eq9) = (O_db/O_m) sqrt(P_cold)
                    + (O_ax/O_m) [fc sqrt(P_cold) + (1-fc) sqrt(P_ax)]

so B -> 1 at low k by construction (verified to <= 6e-4) and CAMB's own linear total (which includes massive neutrinos, outside axionHMcode's budget) is untouched there. P_cold comes from the cdm+baryon transfer combination (== Transfer_nonu; the axion is in Transfer_tot iff DM-like — verified numerically in both regimes). Below omaxh2 = 1e-8 the class switches to axionHMcode's own LCDM recipe (cold-only halo model; the full axion assembly is singular at vanishing fraction).

### Validity domain

DE-like axions (m/H0 < 10) and z grids reaching the KG->EFA switch hard-error regardless of the `strict` yaml flag (the mixed-DM halo model is undefined there; use `halofit_version: original` for DE-like masses). Out-of-calibration axion fractions/masses warn-and-extrapolate (`strict: False`, default, the 2605.12054 practice) or hard-error (`strict: True`). The lensing z grid inherently exceeds dome's 1 < z < 8 calibration — logged once, not gated. Power-law primordial spectra only. Target mass window: m_ax ~ 1e-25..1e-23 eV.

### Validation summary

(full numbers in `.claude/strategy_axionHMcode/13-*.md`)

- Trivial-ratio and ratio==1 nulls: machine precision (1e-15/1e-13).
- External mead2020 ratio through the plumbing reproduces internal mead2020 lensed spectra to <= 3.2e-5 (Cpp, L <= 1000).
- Transfer->P(k) convention identity: 5.7e-8; sigma8 cross-check exact to 4 digits.
- LCDM limit vs CAMB HMcode-2020: -10%/+12% at k = 1 h/Mpc (basic/dome), the agreement level the axionHMcode papers themselves claim.
- Gaughan/Green/Moss reproduction (fax = 0.3, m = 1e-23/-24/-25): boost-ratio ordering, signs and magnitudes match their Fig. 1 (e.g. basic m=1e-24 z=2 R(k=1) = 0.47 vs their ~0.5); lensed-spectra differences match Figs. 2-3 (dome dCpp(L=1000) = +33% vs their ~+35-40%).
- Cross-check vs this branch's axion-aware internal HMcode: <= 0.2% in TT, +10%/+22% in Cpp at L = 500/1000 — two different halo models, as expected.
- Numerical convergence: mass-grid and k-density converged below 1.3e-3 in B (well under the 0.1% TT / 0.5% Cpp targets). Never thin the input k grid: halving it corrupts B by up to 11% at k ~ 1 (alpha/k_star sensitivity).
- Pending: the full-likelihood evaluate/MCMC (needs cobaya-install data); wiring smoke-tested at 83.6 s/eval (dome, 50-node grid, single core; `processes: N` divides wall time with identical numerics).

### Environment requirements

numpy < 2 (this port's python layer is not numpy-2 compatible), scipy < 1.14 (axionHMcode has a dead `scipy.misc` import; the Theory class shims it when absent), numba, astropy, cobaya >= 3.6.2. Note cobaya >= 3.6 removed `use_renames` (renames are unconditional now): the older `EXAMPLE_EVALUATE1/MCMC1.yaml` need that key removed to run on pristine cobaya >= 3.6 (they still work on Cocoa, whose camb.yaml defines the key).

# Appendix D: About the underlying CAMB <a name="appendix-d-about-the-underlying-camb"></a>

This code is built on CAMB 1.6.7 (Antony Lewis and Anthony Challinor, https://camb.info/): a cosmology code for calculating cosmological observables, including CMB, lensing, source count and 21cm angular power spectra, matter power spectra, transfer functions and background evolution; Python package with numerical code in modern Fortran. See the [CAMB documentation](https://camb.readthedocs.io/en/latest/) and the upstream repository at https://github.com/cmbant/CAMB. You will need gfortran installed to compile.

---

The complete port developer guide (change-by-change mapping, the modern-CAMB architecture map, and the exhaustive original-code analyses) is in **[PORT_DEVELOPER_GUIDE.rst](PORT_DEVELOPER_GUIDE.rst)**.
