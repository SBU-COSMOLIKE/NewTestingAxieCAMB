"""Per-halo diagnosis of aggressive vs strict solver outcomes on the
gaughan_z0 dome case: which masses flip soliton status, the diagnostics
counters, and where the negative one-halo cross term comes from."""
import sys
import numpy as np

SCRATCH = ("/private/tmp/claude-501/-Users-vivianmiranda-data-research-WayneHu"
           "-rayne/4e32f7cb-e470-4d65-8ad6-bf48eb7553b7/scratchpad")
CAMB_PATH = "/Users/vivianmiranda/data/research/WayneHu/rayne/New_AxiECAMB"
FORK = "/Users/vivianmiranda/data/research/WayneHu/rayne/fork_axionHMcode"
sys.path.insert(0, CAMB_PATH)
sys.path.insert(0, f"{CAMB_PATH}/axionhmcode_boost")
sys.path.insert(0, FORK)

from axionCAMB_and_lin_PS import lin_power_spectrum
from cosmology.overdensities import func_D_z_unnorm_int
from cosmology import fast_tables
from halo_model import HMcode_params
from axion_functions import axion_params
from axionhmcode_boost import _VERSION_FLAGS

d = np.load(f"{SCRATCH}/fork_val_transfers.npz")
base = "gaughan_z0"
meta = d[f"case_{base}"]
z, h = float(meta[0]), float(meta[5])
cd = {
  "M_min": 7, "M_max": 18, "transfer_kmax": float(d[f"k_{base}"].max()) * h,
  "version": "dome",
  "omega_b_0": float(meta[1]), "omega_d_0": float(meta[2]),
  "omega_ax_0": float(meta[3]),
  "omega_db_0": float(meta[2]) + float(meta[1]),
  "omega_m_0": float(meta[1]) + float(meta[2]) + float(meta[3]),
  "m_ax": float(meta[4]), "h": h, "z": z,
  "ns": float(meta[6]), "As": float(meta[7]), "k_piv": 0.05,
}
for key in ("b", "d", "ax", "db", "m"):
  cd[f"Omega_{key}_0"] = cd[f"omega_{key}_0"] / h**2
cd["Omega_w_0"] = 1 - cd["Omega_m_0"]
cd["G_a"] = func_D_z_unnorm_int(z, cd["Omega_m_0"], cd["Omega_w_0"])

k = d[f"k_{base}"]
t2p = lambda T: lin_power_spectrum.transfer_to_PS(k, T, cd)
cold_T = (cd["Omega_b_0"] * d[f"Tb_{base}"]
          + cd["Omega_d_0"] * d[f"Tc_{base}"]) / cd["Omega_db_0"]
pd = {"k": k, "power_total": t2p(d[f"Tt_{base}"]),
      "power_CDM": t2p(d[f"Tc_{base}"]), "power_baryon": t2p(d[f"Tb_{base}"]),
      "power_cold": t2p(cold_T), "power_axion": t2p(d[f"Ta_{base}"])}
M_arr = np.logspace(7, 18, 100)
hmc = HMcode_params.HMCode_param_dic(cd, k, pd["power_cold"])

results = {}
for mode in (False, True):
  fast_tables.set_legacy_root_finder(not mode)
  cdm = dict(cd)   # fresh cosmo_dic so '_vm_' caches and diag are per mode
  axd = axion_params.func_axion_param_dic(M_arr, cdm, pd, hmc,
                                          concentration_param=True)
  dens = np.asarray(axd["central_dens"] if "central_dens" in axd
                    else axd.get("dens_param", np.nan))
  results[mode] = (axd, dens, cdm.get("_vm_agg_diag"))
  print(f"mode aggressive={mode}: keys {sorted(axd.keys())[:8]}...")

axd_s, dens_s, _ = results[False]
axd_a, dens_a, diag = results[True]
print("\ndiagnostics (aggressive):", diag)
print(f"M_cut = {axd_s['M_cut']:.3e} vs aggressive {axd_a['M_cut']:.3e}")
print(f"frac_cluster strict {axd_s['frac_cluster']:.6f} "
      f"aggressive {axd_a['frac_cluster']:.6f}")

if dens_s.shape == dens_a.shape:
  zero_s, zero_a = dens_s == 0, dens_a == 0
  both = ~zero_s & ~zero_a
  print(f"\nhalos with soliton: strict {np.sum(~zero_s)}, "
        f"aggressive {np.sum(~zero_a)}")
  print(f"flips strict-0 -> aggressive-solved: {np.sum(zero_s & ~zero_a)}")
  print(f"flips strict-solved -> aggressive-0: {np.sum(~zero_s & zero_a)}")
  if np.any(both):
    rel = np.abs(dens_a[both] / dens_s[both] - 1)
    print(f"among both-solved: max |d rho_c/rho_c| = {rel.max():.2e}")
  flip_idx = np.where(zero_s != zero_a)[0]
  for i in flip_idx[:12]:
    print(f"  M = {M_arr[i]:.3e}  strict rho_c = {dens_s[i]:.3e}  "
          f"aggressive rho_c = {dens_a[i]:.3e}")
