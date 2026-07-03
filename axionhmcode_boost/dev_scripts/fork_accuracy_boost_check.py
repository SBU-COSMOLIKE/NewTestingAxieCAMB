"""Convergence check for the new accuracy_boost flag: compute the boost row
at accuracy_boost = 1 (validated defaults) and 2 (doubled M-grid + doubled
fork table nodes) on the saved validation transfers, fork checkout.
"""
import sys
import time
import numpy as np

SCRATCH = ("/private/tmp/claude-501/-Users-vivianmiranda-data-research-WayneHu"
           "-rayne/4e32f7cb-e470-4d65-8ad6-bf48eb7553b7/scratchpad")
CAMB_PATH = "/Users/vivianmiranda/data/research/WayneHu/rayne/New_AxiECAMB"
FORK = "/Users/vivianmiranda/data/research/WayneHu/rayne/fork_axionHMcode"
sys.path.insert(0, CAMB_PATH)
sys.path.insert(0, f"{CAMB_PATH}/axionhmcode_boost")
sys.path.insert(0, FORK)

from axionhmcode_boost import _compute_row, _VERSION_FLAGS

d = np.load(f"{SCRATCH}/fork_val_transfers.npz")
NUIS = {n: None for n in ("alpha_1", "alpha_2", "gamma_1", "gamma_2")}

print(f"{'case':22s} {'max|dB/B| (2 vs 1)':>19s} {'t(1)':>7s} {'t(2)':>7s}")
for base in ["gaughan_z0", "inputfile_z0", "lcdm_z0"]:
  meta = d[f"case_{base}"]
  rows, times = {}, {}
  for boost in (1.0, 2.0):
    n_m = int(round(100 * boost))
    payload = {
      "z": float(meta[0]), "h": float(meta[5]), "omega_b": float(meta[1]),
      "omega_d": float(meta[2]), "omega_ax": float(meta[3]),
      "m_ax": float(meta[4]), "lcdm_mode": bool(meta[8]),
      "ns": float(meta[6]), "As": float(meta[7]), "k_piv": 0.05,
      "version": "dome", "flags": _VERSION_FLAGS["dome"],
      "M_arr": np.logspace(7, 18, n_m), "nuisance": NUIS,
      "m_min_exponent": 7, "m_max_exponent": 18,
      "accuracy_boost": boost,
      "k": d[f"k_{base}"], "T_cdm": d[f"Tc_{base}"], "T_b": d[f"Tb_{base}"],
      "T_ax": d[f"Ta_{base}"], "T_tot": d[f"Tt_{base}"]}
    row = _compute_row(payload)          # first call pays numba JIT
    t0 = time.time()
    row = _compute_row(payload)
    times[boost] = time.time() - t0
    rows[boost] = row
  B1, B2 = rows[1.0]**2, rows[2.0]**2
  dev = float(np.max(np.abs(B2 / B1 - 1)))
  print(f"{base:22s} {dev:19.2e} {times[1.0]:7.2f} {times[2.0]:7.2f}")
