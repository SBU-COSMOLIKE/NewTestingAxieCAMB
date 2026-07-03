"""Simulate consecutive MCMC steps to see whether the current pattern
(fresh pool per call, parent never warmed) re-pays the numba JIT every
step, and whether pre-warming the parent once fixes it.

Fresh process: the parent has NOT compiled numba (mirrors the real theory,
whose initialize() only imports).
"""
import sys
import time
import numpy as np
import multiprocessing as mp

SCRATCH = ("/private/tmp/claude-501/-Users-vivianmiranda-data-research-WayneHu"
           "-rayne/4e32f7cb-e470-4d65-8ad6-bf48eb7553b7/scratchpad")
CAMB_PATH = "/Users/vivianmiranda/data/research/WayneHu/rayne/New_AxiECAMB"
FORK = "/Users/vivianmiranda/data/research/WayneHu/rayne/fork_axionHMcode"
sys.path.insert(0, CAMB_PATH)
sys.path.insert(0, f"{CAMB_PATH}/axionhmcode_boost")
sys.path.insert(0, FORK)

from axionhmcode_boost import _compute_row, _VERSION_FLAGS

d = np.load(f"{SCRATCH}/fork_val_transfers.npz")
meta = d["case_gaughan_z0"]


def make_payload():
  return {
    "z": float(meta[0]), "h": float(meta[5]), "omega_b": float(meta[1]),
    "omega_d": float(meta[2]), "omega_ax": float(meta[3]),
    "m_ax": float(meta[4]), "lcdm_mode": bool(meta[8]),
    "ns": float(meta[6]), "As": float(meta[7]), "k_piv": 0.05,
    "version": "dome", "flags": _VERSION_FLAGS["dome"],
    "M_arr": np.logspace(7, 18, 100),
    "nuisance": {n: None for n in ("alpha_1", "alpha_2", "gamma_1", "gamma_2")},
    "m_min_exponent": 7, "m_max_exponent": 18, "accuracy_boost": 1.0,
    "legacy_root_finder": False,
    "k": d["k_gaughan_z0"], "T_cdm": d["Tc_gaughan_z0"],
    "T_b": d["Tb_gaughan_z0"], "T_ax": d["Ta_gaughan_z0"],
    "T_tot": d["Tt_gaughan_z0"]}


NZ, NW = 50, 8
ctx = mp.get_context("fork")


def one_step():
  payloads = [make_payload() for _ in range(NZ)]
  t0 = time.time()
  with ctx.Pool(NW) as pool:
    pool.map(_compute_row, payloads)
  return time.time() - t0


print("CURRENT pattern (parent never warmed):")
for step in range(1, 4):
  print(f"  step {step}: {one_step():.2f} s")

print("\nAfter warming the parent once (single serial _compute_row):")
t0 = time.time(); _compute_row(make_payload()); tw = time.time() - t0
print(f"  parent warm-up call: {tw:.2f} s (one-time)")
for step in range(1, 4):
  print(f"  step {step}: {one_step():.2f} s")
