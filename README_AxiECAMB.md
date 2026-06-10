# AxiECAMB (modern-CAMB port)

This is **AxiECAMB** — the ultralight-axion (ULA) effective method of
[arXiv:2412.15192](https://arxiv.org/abs/2412.15192) (Liu, Hu et al.) — ported from its
original CAMB-Nov13 base onto **modern CAMB 1.6.7**, including the Python wrapper.

When using the axion module, please cite [arXiv:2412.15192](https://arxiv.org/abs/2412.15192)
(and Passaglia & Hu 2022, [arXiv:2201.10238](https://arxiv.org/abs/2201.10238), on which the
method builds). The original AxiECAMB heavily modified
[axionCAMB](https://github.com/dgrin1/axionCAMB) (Hlozek et al., arXiv:1410.2896).

## Method

The axion has a quadratic potential and is evolved in synchronous gauge:

- **Background**: the exact Klein-Gordon (KG) equation is solved (16-stage 8th-order
  fixed-step Runge-Kutta in ln a) from deep radiation domination until m = dfac·H, with
  the initial field value found by shooting to match the requested relic abundance.
  At the switch the field is projected onto WKB cos/sin amplitudes and matched onto an
  **effective fluid** (EFA) whose density follows a⁻³ with an exp[3∫w dln a] residual,
  w(a) = wEFA_c (H/m)²; the matching coefficients (⟨H⟩, wEFA_c) are iterated to
  self-consistency. dfac (default 10) is retuned internally: the oscillation phase at
  the switch is targeted to 2β = 7.08π for light DM-like axions, and the switch is
  pushed out of the recombination window z ∈ (800, 1300).
- **Perturbations**: exact KG before the switch (with a per-k conditioning rescale of
  δφ), the (1+w)-weighted GDM effective fluid after, with sound speed
  cs² = (√(1+κ)−1)²/κ + (5/4)(H/am)², κ = k²/(a²m²). At the switch the KG variables are
  projected onto (δ_ax, u_ax) preserving velocity and shear continuity; the residual
  metric jump is absorbed into η (sub-horizon) or carried as a delta-function boundary
  term in the temperature line-of-sight integral (super-horizon).
- **DM vs DE**: for m/H0 ≥ 10 the axion is dark-matter-like — it counts in the matter
  transfer functions, σ₈, the equality redshift, CosmoMC θ, and halofit Ω_m. For
  m/H0 < 10 (m ≲ 1.4e-32 eV) it is dark-energy-like: KG is solved to a = 1, there is no
  fluid switch, and the matter transfer excludes the axion.

## Usage (Fortran / .ini)

```
make camb                     # in fortran/ (build forutils first if needed)
./camb ../inifiles/params_axion.ini
```

New ini keys (see `inifiles/params_axion.ini`):

| key | meaning |
|---|---|
| `m_ax` | ULA mass in eV (negative input = log10(m_ax/eV)) |
| `use_axfrac` | T: use (`omdah2`, `axfrac`); F: use `omaxh2` (+ usual `omch2`) |
| `omaxh2` | Ω_ax h² (when `use_axfrac = F`) |
| `omdah2` | total dark-matter Ω h² (when `use_axfrac = T`) |
| `axfrac` | axion fraction of DM (m/H0 ≥ 10) or of DE (m/H0 < 10) |
| `axion_dfac` | switch threshold m = dfac·H (default 10; retuned internally) |
| `axion_isocurvature`, `Hinf` | accepted but isocurvature is force-disabled (v1.0 parity) |

With `use_axfrac = T`, `omch2` may be omitted (it is derived). Constant-w dark energy
(fluid or PPF) can be combined with the axion; quintessence dark-energy models cannot
(the axion background solver treats DE as Λ, as in the original).

## Usage (Python)

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

The axion density perturbation is available as the `delta_axion` matter transfer
column (`camb.model.Transfer_axion`). `delta_tot` (and hence σ₈ and the default
matter power) includes the axion when it is DM-like, excludes it when DE-like.

## What was ported and where

- `fortran/AxionBackground.f90` (new): the KG background solver `w_evolve`, EFA
  matching `auxiIC`, phase targeting and recombination-skip dfac retuning (moved here
  from the old `inidriver_axion.F90` so the Python interface gets them too), as the
  component class `TAxionModel` stored in `CAMBparams%Axion`.
- `results.f90`: density budget/closure with the axion, the solver invocation,
  τ_osc, background integrals split at the dtauda kink at a_osc (applied uniformly:
  times, distances, sound horizons, optical depths — the original only split some),
  fine time-step window around τ_osc, thermo values cached at τ_osc, `Transfer_axion`
  column, z_eq and CosmoMC θ definitions.
- `equations.f90`: the two axion perturbation equations (KG ↔ EFA), the mid-evolution
  switch in the `next_switch` chain with the WKB projection (`AxionSwitchKGtoEFA`),
  adiabatic δφ initial conditions, axion terms in dgrho/dgq/grho/gpres, low-k lmaxnr
  boost ("WH smoother"). Tensors need no axion terms: the modern tensor background
  comes from the dtauda-based thermo table (this also fixes an original-code issue
  where tensors extrapolated the field table past a_osc).
- `cmbmain.f90`: the switch boundary term in the temperature LOS integral
  (`deltaBCSrc` machinery, flat and curved cases), axion-aware integration start time.
- `recfast.f90`: dHdz in the tightly-coupled T_mat term includes the axion (numerical
  derivative of the exact H(z), stepped away from the a_osc kink).
- `halofit.f90`: axion counted in Ω_m (DM-like) or in the smooth DE (DE-like), with a
  warning that the non-linear mode is inherited from axionCAMB and not well tested.
- Python: `camb/axion.py` (`AxionModel`), `CAMBparams.set_axion`, transfer-name lists.

## Validation against the original AxiECAMB

With matched cosmologies, the axion/ΛCDM suppression ratios agree between the original
AxiECAMB (Nov13 base) and this port to:

| case | TT C_ℓ ratio (ℓ=2–2600) | P(k) ratio |
|---|---|---|
| m=1e-27 eV, 10% of DM (switch z≈1341) | ≤ 0.01% | ≤ 0.005% |
| m=1e-27 eV, 100% of DM (`use_axfrac`) | ≤ 0.10% | ≤ 0.005% (where suppression < 10³) |
| m=1e-30 eV, 10% of DM (switch z≈24, boundary term active) | ≤ 0.08% | ≤ 0.01% |
| m=1.4e-33 eV (DE-like, no switch) | ≤ 0.09% | ≤ 0.01% |

Residuals at the 0.03–0.1% level are dominated by Nov13 ↔ CAMB-1.6.7 baseline physics
differences that do not perfectly cancel in the ratios (the absolute ΛCDM baselines
differ by ~0.2%). The pure-ΛCDM limit of this code is bit-identical to unmodified
CAMB 1.6.7. The standard CAMB Python test suite passes.

## Warnings / known differences (carried over or documented)

- **Isocurvature is disabled** (as in AxiECAMB v1.0; the original mode-6 vector
  targets variables that are not evolved). Inputs are accepted and ignored with a
  warning.
- The **growth-rate (Transfer_f) column** of the original is disabled there and was
  not ported; modern CAMB's own growth outputs are available.
- The **non-linear mode** is inherited from axionCAMB and not extensively tested;
  `halofit_version = 1` (original) is suggested — Takahashi was found unstable for
  axion models.
- For z > 0 transfer outputs, mind whether the requested z is before or after the
  switch: the axion density contrast is defined differently in the two regimes.
- P(k) at wavenumbers where the spectrum is suppressed by ≳6 orders of magnitude
  differs from the original (which zeroed rather than extrapolated the dead tail);
  set `transfer_kmax` high enough for any application sensitive to that region.
- Default accuracy (`accuracy_boost = 1`) is what was validated in arXiv:2412.15192;
  higher boosts apply to the non-ULA accuracy settings only.
- CosmoMC θ counts the axion in ω_dm unconditionally (original behaviour), which is
  only meaningful for DM-like axions.

Intermediate analysis artifacts from the porting work (full diff inventories of the
original vs CAMB Nov13, design document) are in `../.port_analysis/`.
