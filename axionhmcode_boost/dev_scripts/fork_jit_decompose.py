"""Decompose the single-evaluate boost cost: cold numba JIT vs warm work,
and whether pre-warming the PARENT process (so forked workers inherit the
compiled numba state) removes the per-pool recompile.
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
base = "gaughan_z0"
meta = d[f"case_{base}"]


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
    "k": d[f"k_{base}"], "T_cdm": d[f"Tc_{base}"], "T_b": d[f"Tb_{base}"],
    "T_ax": d[f"Ta_{base}"], "T_tot": d[f"Tt_{base}"]}


NZ, NW = 50, 8

# (1) cold vs warm single call -> the numba JIT floor
p = make_payload()
t0 = time.time(); _compute_row(p); t_cold = time.time() - t0
t0 = time.time(); _compute_row(p); t_warm = time.time() - t0
print(f"single _compute_row: cold {t_cold:.2f} s, warm {t_warm:.3f} s "
      f"-> JIT floor ~ {t_cold - t_warm:.2f} s")
print(f"serial {NZ} nodes (warm): {NZ * t_warm:.1f} s")

payloads = [make_payload() for _ in range(NZ)]

# (2) pool with a COLD parent (parent never compiled -> each worker JITs)
#     NOTE: this parent is already warm from step (1); to emulate a truly
#     cold parent we time a fresh subprocess below instead. Here we measure
#     the WARM-parent pool: workers fork an already-compiled parent.
ctx = mp.get_context("fork")
t0 = time.time()
with ctx.Pool(NW) as pool:
  rows = pool.map(_compute_row, payloads)
t_pool_warm_parent = time.time() - t0
print(f"{NW}-worker pool, warm parent (forks inherit compiled numba): "
      f"{t_pool_warm_parent:.2f} s")
