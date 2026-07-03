---
name: axionhmcode-readability-pass
description: Plan for the code-simplification and formal-documentation pass over the boost theory and the fork's VM-SPEEDUP additions (PI directive 2026-07-03), including the target inventory and the re-validation protocol
metadata:
  type: project
---

# Readability pass: simplify the Python, formalize the documentation

PI directive (2026-07-03, conditional on the round-4 V5 postprocess landing
green, which it did — see file 14): simplify Python that is hard for a human
whose main language is not Python, in the spirit of the earlier PyTorch
teaching project, and provide formal, good documentation. The working rules
live in the auto-memory note `readable-python-code-rules` and are summarized
here so the pass is reproducible.

## Rules

1. Anonymous lambdas doing real work become named functions with docstrings.
2. Comprehensions and chained one-liners become plain loops and named
   intermediate variables (with physical names; units in comments) unless
   the plain form has a real cost.
3. Vectorized numpy idioms (slice arithmetic, fancy indexing, FFT
   conjugation conventions, id()-based cache fingerprints) are either
   replaced by the plain form when that is cost-free, or kept with a
   derivation comment a non-Python reader can follow line by line.
4. Every function gets a formal docstring: purpose, arguments with units
   and array shapes, return value, and the equation/paper it implements.
5. Numerics must not change: after each file, re-run the fork_validate gate
   (byte-identical expectation on the legacy path) and the boost pytest
   suite; default-solver behavior re-checked with dev_scripts/
   fork_aggressive_map.py.
6. VM-SPEEDUP fences stay greppable; boost code keeps the group 2-space
   style, fork code keeps upstream's 4-space style.

## Target inventory (ordered)

1. `New_AxiECAMB/axionhmcode_boost/axionhmcode_boost.py` — the file the PI
   reads most. Known offenders: the `t2p` lambda in `_compute_row`; the
   payload-dict plumbing (documented but dense); dict comprehensions in
   `get_requirements`/`_nuisance_values`; the derived-parameter lambdas are
   yaml-side and exempt. Add a module-level "reading guide" section.
2. `fork_axionHMcode/cosmology/fast_tables.py` — the weight-vector slice
   arithmetic in `geom_simpson_grid` (w[0:m:2] pair assembly), the FFTLog
   conjugation identity and the windowed-slice index gymnastics in
   `fftlog_j0_eval`, the 4-point Lagrange weights, cache fingerprints
   (id(PS) semantics deserve a plain-language paragraph).
3. `fork_axionHMcode/halo_model/axion_density_profile.py` (fenced blocks) —
   `_ax_halo_mass_smooth` cell-correction algebra; the bracket-expansion
   loop in the default solver.
4. `fork_axionHMcode/halo_model/cold_density_profile.py` (fenced blocks) —
   memo-key construction; the formation-z branch logic.
5. `dev_scripts/` — bring the measurement scripts up to the same standard
   (they are read as records).
6. Documentation deliverable: expanded docstrings throughout + a short
   CODE_TOUR.md in axionhmcode_boost/ walking a non-Python reader through
   one full boost evaluation, file by file, function by function.

## Status

- [ ] File 1 (boost theory)
- [ ] File 2 (fast_tables)
- [ ] File 3 (axion_density_profile fences)
- [ ] File 4 (cold_density_profile fences)
- [ ] File 5 (dev_scripts)
- [ ] CODE_TOUR.md
