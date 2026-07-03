"""Child: compute boost rows with a given axionHMcode checkout.
argv: <checkout_path> <transfers.npz> <out.npz>
"""
import sys
import time
import numpy as np

CHECKOUT, TRANSFERS, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
CAMB_PATH = "/Users/vivianmiranda/data/research/WayneHu/rayne/New_AxiECAMB"
sys.path.insert(0, CAMB_PATH)
sys.path.insert(0, f"{CAMB_PATH}/axionhmcode_boost")
sys.path.insert(0, CHECKOUT)

from axionhmcode_boost import _compute_row, _VERSION_FLAGS

d = np.load(TRANSFERS)
NUIS = {n: None for n in ("alpha_1", "alpha_2", "gamma_1", "gamma_2")}
results, times = {}, {}

for case in [c for c in d.files if c.startswith("case_")]:
  meta = d[case]  # [z, omega_b, omega_d, omega_ax, m_ax, h, ns, As, lcdm]
  base = case.replace("case_", "")
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
      "k": k, "T_cdm": d[f"Tc_{base}"], "T_b": d[f"Tb_{base}"],
      "T_ax": d[f"Ta_{base}"], "T_tot": d[f"Tt_{base}"]}
    t0 = time.time()
    row = _compute_row(payload)
    dt = time.time() - t0
    # warm rerun for a fair timing (JIT paid on first)
    t0 = time.time()
    row = _compute_row(payload)
    times[f"{base}_{version}"] = time.time() - t0
    results[f"{base}_{version}"] = row

np.savez(OUT, **results,
         **{f"time_{key}": val for key, val in times.items()})

# primitive checks (fork only has fast paths; recompute originals inline)
if OUT.endswith("_fork.npz"):
  from scipy import integrate
  from cosmology.basic_cosmology import func_E_z
  from cosmology.overdensities import func_D_z_unnorm
  from cosmology import fast_tables
  Om, Ow = 0.30964, 0.69036
  print("primitive checks (fork vs original quadratures):")
  for z in [0.0, 1.0, 5.0, 9.8]:
    f = lambda y, x: func_E_z(x, Om, Ow) / (1 + x) * (1 + y) \
        / func_E_z(y, Om, Ow)**3
    G_ref = 2.5 * Om * integrate.dblquad(f, z, 10000, lambda x: x, 10000)[0]
    G_fast = fast_tables.G_integral_fast(z, Om, Ow)
    Dn_ref = func_D_z_unnorm(z, Om, Ow) / func_D_z_unnorm(0.0, Om, Ow)
    Dn_fast = fast_tables.D_norm_fast(z, Om, Ow)
    print(f"  z={z:4.1f}: dG/G = {G_fast/G_ref-1:+.2e}, "
          f"dDn/Dn = {Dn_fast/Dn_ref-1:+.2e}")
