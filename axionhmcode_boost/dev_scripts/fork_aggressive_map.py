"""Aggressive-mode difference map: fork with aggressive_optimization=True
vs unmodified upstream (saved rows), on the validation cases. Also reports
the per-evaluation solver diagnostics and timings."""
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
up = np.load(f"{SCRATCH}/fork_val_upstream.npz")
NUIS = {n: None for n in ("alpha_1", "alpha_2", "gamma_1", "gamma_2")}

print(f"{'case':26s} {'max|dB/B| vs upstream':>22s} {'at k':>8s} {'t':>7s}")
for base in ["gaughan_z0", "gaughan_z2", "inputfile_z0", "inputfile_z2",
             "lcdm_z0", "lcdm_z2"]:
  meta = d[f"case_{base}"]
  k = d[f"k_{base}"]
  for version in ("dome", "basic"):
    payload = {
      "z": float(meta[0]), "h": float(meta[5]), "omega_b": float(meta[1]),
      "omega_d": float(meta[2]), "omega_ax": float(meta[3]),
      "m_ax": float(meta[4]), "lcdm_mode": bool(meta[8]),
      "ns": float(meta[6]), "As": float(meta[7]), "k_piv": 0.05,
      "version": version, "flags": _VERSION_FLAGS[version],
      "M_arr": np.logspace(7, 18, 100), "nuisance": NUIS,
      "m_min_exponent": 7, "m_max_exponent": 18,
      "accuracy_boost": 1.0, "legacy_root_finder": False,
      "k": k, "T_cdm": d[f"Tc_{base}"], "T_b": d[f"Tb_{base}"],
      "T_ax": d[f"Ta_{base}"], "T_tot": d[f"Tt_{base}"]}
    row = _compute_row(payload)          # warm
    t0 = time.time()
    row = _compute_row(payload)
    dt = time.time() - t0
    B_ag, B_up = row**2, up[f"{base}_{version}"]**2
    dev = np.abs(B_ag / B_up - 1)
    i = int(np.argmax(dev))
    print(f"{base+'_'+version:26s} {dev.max():22.2e} {k[i]:8.3f} {dt:7.3f}")
