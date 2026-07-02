---
name: axionhmcode-mass-prior-recipes
description: PI-decided mass window (1e-25..1e-23 eV) and README-ready yaml recipes for fixed-mass and sampled-log-mass runs
metadata:
  type: project
---

# Mass prior recipes (README-ready instructions)

PI decision (2026-07-01): the boost pipeline targets m_ax ~ 1e-25 .. 1e-23 eV — the window
where the axion Jeans scale sits in the quasi-linear regime probed by CMB lensing
(arXiv:2605.12054) and the axion is deep DM-like for AxiECAMB (m/H0 >> 10).
Two run modes must both be documented in the final example files; the recipes below are
the content for that README (drop-in, up to file-name changes).

Convention used throughout (same as EXAMPLE_EVALUATE1.yaml): `logmx` = log10(m_ax / eV) is
the sampled/user-facing parameter; `m_ax` (eV) is what CAMB receives, via the value lambda.

## Mode A — EVALUATE yaml (single cosmology)

Directive: keep the parameter block in sampled-log-mass form (prior on `logmx` over the
target window), and pin the mass in the `evaluate` sampler's `override` section — we only
run one cosmology in an evaluate, so the prior is never actually explored, but the block
is then copy-paste ready for a sampling run.

```yaml
params:
  # ... other cosmological parameters as in EXAMPLE_EVALUATE1.yaml ...
  logmx:
    prior:
      min: -25.0
      max: -23.0
    ref:
      dist: norm
      loc: -24.0
      scale: 0.1
    proposal: 0.1
    latex: \log{m_{\rm ax}}
    drop: true
  m_ax:
    value: 'lambda logmx: 10**logmx'
    latex: m_\mathrm{ax}

sampler:
  evaluate:
    override:
      # single-point evaluation: the mass is fixed HERE, not in the prior block
      logmx: -24.0
      # ... overrides for the remaining sampled params ...
```

## Mode B — MCMC yaml (the committed example: fixed mass)

Directive: the example MCMC fixes the mass (one chain per mass, the Gaughan/Rogers
practice — the mass posterior is poorly constrained and multimodal-prone, so production
science runs scan a small set of fixed masses). The yaml carries comments showing exactly
how to switch to a uniform log-mass prior.

```yaml
params:
  # ... other cosmological parameters ...
  m_ax:
    value: 1.0e-24        # fixed mass for this chain; run one chain per target mass
    latex: m_\mathrm{ax}
  # --- To sample the mass instead (uniform prior in log10 m_ax/eV), replace the
  # --- fixed m_ax block above with the two blocks below:
  # logmx:
  #   prior:
  #     min: -25.0
  #     max: -23.0
  #   ref:
  #     dist: norm
  #     loc: -24.0
  #     scale: 0.1
  #   proposal: 0.1
  #   latex: \log{m_{\rm ax}}
  #   drop: true
  # m_ax:
  #   value: 'lambda logmx: 10**logmx'
  #   latex: m_\mathrm{ax}
```

## Notes that belong in the README next to the recipes

- Regime sanity: the whole window is DM-like (`results.Params.Axion.is_de_like` False);
  the old EXAMPLE yamls' window logmx in [-34,-31] is DE-like and must NOT be reused with
  the boost theory — a DE-like point hard-errors regardless of the `strict` flag
  ([[axionhmcode-architecture]] "Regime gating").
- dome-version calibration is centred at m = 1e-24.5 eV; at the window edges (especially
  1e-23 eV, 1.5 decades up) dome results are extrapolations — governed by the theory-block
  `strict` flag (False = warn + extrapolate, True = hard error). basic version covers the
  full window.
- dome is the default version (collaborator guidance 2026-07-01 — most recent
  recalibration). Sampled
  log-mass runs cross the dome calibration boundary continuously: under `strict: False`
  they warn + extrapolate (treat such chains as sensitivity tests, Gaughan et al.
  Sec. IV B); basic remains available for wide-mass scans via `version: basic`.
- Axion density parameterization: the port natively supports either `omaxh2` directly or
  (`omdah2`, `axfrac`) with use_axfrac ([[axiecamb-port-project]]). The Gaughan-style
  science parameterization {Omega_D h^2 total, f_ax} maps to (omdah2, axfrac) — recommended
  for production since it keeps the CMB-friendly total-dark-matter degeneracy direction;
  the fixed-mass example can use either.
- The DE-like-window example yamls (EXAMPLE_EVALUATE1/EXAMPLE_MCMC1) remain valid for the
  halofit-original path; the new boost examples are separate files, not edits to them.
