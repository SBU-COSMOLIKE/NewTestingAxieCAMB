"""Diagnose the GL gate failure: recompute the failing dome rows with
increasing GL node counts and compare each to the saved upstream rows.
Convergence toward upstream = quadrature error; a persistent jump = a
flipped no-solution rejection in the central-density solver.
"""
import sys
import numpy as np

SCRATCH = ("/private/tmp/claude-501/-Users-vivianmiranda-data-research-WayneHu"
           "-rayne/4e32f7cb-e470-4d65-8ad6-bf48eb7553b7/scratchpad")
CAMB_PATH = "/Users/vivianmiranda/data/research/WayneHu/rayne/New_AxiECAMB"
FORK = "/Users/vivianmiranda/data/research/WayneHu/rayne/fork_axionHMcode"
sys.path.insert(0, CAMB_PATH)
sys.path.insert(0, f"{CAMB_PATH}/axionhmcode_boost")
sys.path.insert(0, FORK)

from axionhmcode_boost import _compute_row, _VERSION_FLAGS
from cosmology import fast_tables

d = np.load(f"{SCRATCH}/fork_val_transfers.npz")
up = np.load(f"{SCRATCH}/fork_val_upstream.npz")
NUIS = {n: None for n in ("alpha_1", "alpha_2", "gamma_1", "gamma_2")}

for base in ["gaughan_z0", "inputfile_z2", "gaughan_z2"]:
  meta = d[f"case_{base}"]
  k = d[f"k_{base}"]
  B_up = up[f"{base}_dome"]**2
  print(f"--- {base} (dome) ---")
  for n_gl in (128, 256, 512, 1024):
    fast_tables._GL_R_NODES = n_gl
    payload = {
      "z": float(meta[0]), "h": float(meta[5]), "omega_b": float(meta[1]),
      "omega_d": float(meta[2]), "omega_ax": float(meta[3]),
      "m_ax": float(meta[4]), "lcdm_mode": bool(meta[8]),
      "ns": float(meta[6]), "As": float(meta[7]), "k_piv": 0.05,
      "version": "dome", "flags": _VERSION_FLAGS["dome"],
      "M_arr": np.logspace(7, 18, 100), "nuisance": NUIS,
      "m_min_exponent": 7, "m_max_exponent": 18,
      "k": k, "T_cdm": d[f"Tc_{base}"], "T_b": d[f"Tb_{base}"],
      "T_ax": d[f"Ta_{base}"], "T_tot": d[f"Tt_{base}"]}
    B = _compute_row(payload)**2
    dev = np.abs(B / B_up - 1)
    i = int(np.argmax(dev))
    print(f"  n_gl={n_gl:5d}: max|dB/B| = {dev.max():.2e} "
          f"at k = {k[i]:.3f} h/Mpc")
