---
name: axionhmcode-no-circular-dependency
description: The CAMB↔axionHMcode circular dependency is already solved by Cobaya PR#480's get_non_linear_ratio(results) mechanism
metadata:
  type: project
---

# Key discovery: the circular dependency does not exist

The prompt's central open question — "CAMB (output P_L) → axionHMcode (output sqrt(P_NL/P_L))
→ CAMB again (output lensed spectra), which I don't know how to implement via Theory blocks" —
is already solved by the design of Cobaya PR#480, and every needed piece is present in the
working tree. No new dependency-graph tricks and (most likely) no Cobaya patches are needed.

## Why there is no cycle (the "sandwich" question, answered in full)

The natural mental model — an axionHMcode Theory block sandwiched between two CAMB theory
blocks, "CAMB → boost → CAMB again" — would indeed be rejected: Cobaya resolves
`get_requirements()`/`provides` declarations into a directed acyclic graph, and a theory that
requires P_L from camb while camb requires the boost from it is a true cycle. The sandwich
works anyway for two reasons, neither visible from the yaml.

### Reason 1 — the `camb:` block is two graph nodes, not one

When Cobaya initializes the CAMB wrapper, `CAMB.get_helper_theories()` (camb.py:1049-1061)
spawns a separate component named `camb.transfers` (class `CambTransfers`, camb.py:1162).
That helper runs the actual Boltzmann solve — the expensive stage — and provides
`CAMB_transfers` (the CAMBdata object holding time sources and matter transfer functions).
The main `camb` node never runs the Boltzmann code: its `calculate()` takes the cached
transfers object and finishes with `power_spectra_from_transfer()` — the cheap stage (fold
the primordial spectrum into precomputed transfers, apply the nonlinear rescaling, produce
lensed C_ell and P(k)). "CAMB is called again" never happens; CAMB is internally a two-stage
code and Cobaya exposes each stage as its own node.

The physics that makes the split exact: transfer functions T(k,z) do not depend on As, ns,
or the nonlinear model; P_L(k,z) = P_prim(k) x T^2(k,z), assembled after the fact. This is
also why caching survives: `camb.transfers` caches its result when only primordial or
nonlinear parameters change (docstring, camb.py:1164-1166), preserving Cobaya's fast/slow
parameter split with the boost in place.

Precedent: this split is not new in PR#480 — it is the same machinery behind
`external_primordial_pk`, which the group already used for GSR (`GSRPrimordialPk` is a
Theory sandwiched between `camb.transfers` and `camb` in exactly this way; see the camb-dev
skill notes). PR#480 added a second quantity, `non_linear_ratio`, flowing through the same gap.

### Reason 2 — P_L never enters the dependency graph at all

The ratio provider declares ONLY `CAMB_transfers` as a requirement (like
`TrivialNonLinearRatio` at `tests/test_cosmo_multi_theory.py:280-281`). It must NOT request
`Pk_grid` — that quantity is provided by the final `camb` node, so requesting it would
recreate the cycle for real (see trap T2 below). The dependency chain is linear:

    camb.transfers  →  AxionHMcodeBoost  →  camb (final spectra)  →  likelihoods

P_L reaches the provider by ordinary argument passing, invisible to the dependency resolver:
inside `CAMB.calculate()` (camb.py:676), after primordial power is set on the results object
(line 714) and BEFORE `power_spectra_from_transfer()` (line 731), the wrapper calls

    non_linear_ratio = self.provider.get_non_linear_ratio(results)   # camb.py:717
    results.Params.NonLinearModel.set_ratio(k_h, z, ratio)           # camb.py:718-722

i.e. the provider's `get_non_linear_ratio(self, results)` method receives the live CAMBdata
object as a call argument — a plain Python method dispatch at runtime, not a graph edge.
The provider extracts the linear spectrum from `results` itself, runs axionHMcode, and
returns `{"k_h": ..., "z": ..., "ratio": ...}` with ratio = sqrt(P_NL/P_L), shape
(len(z), len(k_h)). Cobaya auto-detects that a Theory provides `non_linear_ratio` from the
existence of the `get_non_linear_ratio` method.

Proof the pattern is allowed: `test_trivial_non_linear_ratio`
(test_cosmo_multi_theory.py:329) constructs precisely this sandwich and passes in the
upstream test suite. The Boltzmann solve executes exactly once per MCMC step.

## Everything is already installed

- The Cocoa-cloned Cobaya (`cocoa/Cocoa/cobaya`, recent master, HEAD 899f30a4 "version bump")
  contains the full PR#480 machinery: `use_non_linear_ratio` class attr (camb.py:263),
  auto-selection of `ExternalNonLinearRatio` as the nonlinear model (camb.py:331-334),
  the must_provide hook (camb.py:640-641), the calculate-time call (camb.py:716-722),
  and the reference test `test_trivial_non_linear_ratio` (test_cosmo_multi_theory.py:329).
- Cocoa's full-file replacement `cobaya_changes/cobaya/theories/camb/camb.yaml` retains
  `use_non_linear_ratio: False` (line 26) — the Cocoa patching pipeline does not break it.
- New_AxiECAMB already carries the CAMB-side implementation, both Fortran
  (`New_AxiECAMB/fortran/ExternalNonLinearRatio.f90`, `TExternalNonLinearRatio`) and Python
  (`New_AxiECAMB/camb/nonlinear.py:322`, `ExternalNonLinearRatio.set_ratio`).

## The two traps

T1 — inside `get_non_linear_ratio(results)` the results object is transfers-only. Calling
`results.get_linear_matter_power_spectrum(..., have_power_spectra=False)` triggers
`self.calc_power_spectra(params)` (`New_AxiECAMB/camb/results.py:835-836`), which applies the
nonlinear model — `ExternalNonLinearRatio` with no ratio set yet — and the Fortran hard-stops:
`error stop 'ExternalNonLinearRatio: ratio not set'` (ExternalNonLinearRatio.f90:69-70).

The safe extraction is `results.get_matter_transfer_data()` (results.py:768), which only reads
the already-computed matter transfer functions. This is doubly convenient because axionHMcode's
native input IS transfer functions in the axionCAMB convention (see [[axionhmcode-api]]).
Matter transfers ARE available in the lensing case: with `non_linear_sources=True` the helper
calls `get_transfer_functions(camb_params, only_time_sources=True)` (camb.py:1219-1221), and
`only_time_sources` skips only the CMB l,k transfer functions and the nonlinear scaling — not
the matter transfers (`New_AxiECAMB/camb/camb.py:58-72`).

T2 — the boost theory must never declare `Pk_grid` (or any power-spectrum product of the
final `camb` node) in `get_requirements()`. Those quantities are provided by `camb`, which
itself requires `non_linear_ratio` from the boost theory — declaring them turns the imagined
sandwich cycle into a real one and Cobaya will refuse to build the model. The only legal
requirement is `CAMB_transfers`; everything else comes off the `results` argument at call
time (Reason 2 above).

Related: [[axionhmcode-architecture]], [[axionhmcode-verified-facts]].
