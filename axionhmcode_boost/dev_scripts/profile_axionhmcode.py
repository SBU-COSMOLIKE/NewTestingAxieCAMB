"""Deep performance study: where do axionHMcode's seconds actually go?

Stage-level timings (HMCode_param_dic / func_axion_param_dic /
func_full_halo_model_ax), then cProfile of one full dome redshift, plus a
comparison against CAMB's internal Fortran HMcode-2020 cost for ALL redshifts.
"""
import sys
import time
import cProfile
import pstats
import io

import numpy as np

AXHM = "/Users/vivianmiranda/data/research/WayneHu/rayne/axionHMcode"
CAMB_PATH = "/Users/vivianmiranda/data/research/WayneHu/rayne/New_AxiECAMB"
BOOST = f"{CAMB_PATH}/axionhmcode_boost"
for p in (AXHM, CAMB_PATH, BOOST):
  sys.path.insert(0, p)

import camb
from camb import model as cm
from axionCAMB_and_lin_PS import lin_power_spectrum
from cosmology.overdensities import func_D_z_unnorm_int
from halo_model import HMcode_params, PS_nonlin_axion
from axion_functions import axion_params
from axionhmcode_boost import _compute_row, _VERSION_FLAGS

fid = dict(H0=67.32, ombh2=0.022383, As=2.101e-9, ns=0.96605,
           num_massive_neutrinos=1, mnu=0.06, nnu=3.046)
OMD, FAX, M_AX = 0.12011, 0.3, 1e-24
omaxh2 = OMD * FAX
omch2 = OMD - omaxh2
NUIS = {n: None for n in ("alpha_1", "alpha_2", "gamma_1", "gamma_2")}

# ---- transfers once --------------------------------------------------------
pars = camb.set_params(omch2=omch2, omaxh2=omaxh2, m_ax=M_AX,
                       redshifts=[0.0], kmax=10.0, WantTransfer=True, **fid)
pars.NonLinear = cm.NonLinear_none
res = camb.get_results(pars)
td = np.asarray(res.get_matter_transfer_data().transfer_data, dtype=np.float64)
k = td[0, :, 0]
print(f"nk = {len(k)}")


def payload(version):
  return {"z": 0.0, "h": fid["H0"] / 100, "omega_b": fid["ombh2"],
          "omega_d": omch2, "omega_ax": omaxh2, "m_ax": M_AX,
          "lcdm_mode": False, "ns": fid["ns"], "As": fid["As"],
          "k_piv": 0.05, "version": version,
          "flags": _VERSION_FLAGS[version],
          "M_arr": np.logspace(7, 18, 100), "nuisance": NUIS,
          "m_min_exponent": 7, "m_max_exponent": 18,
          "k": k, "T_cdm": td[1, :, 0], "T_b": td[2, :, 0],
          "T_ax": td[13, :, 0], "T_tot": td[6, :, 0]}


# warm the numba JIT
_compute_row(payload("dome"))

# ---- stage-level timing -----------------------------------------------------
def build_inputs(version):
  cd = {"M_min": 7, "M_max": 18, "transfer_kmax": float(k.max()) * 0.6732,
        "version": version, "omega_b_0": fid["ombh2"], "omega_d_0": omch2,
        "omega_ax_0": omaxh2, "omega_db_0": omch2 + fid["ombh2"],
        "omega_m_0": fid["ombh2"] + omch2 + omaxh2, "m_ax": M_AX,
        "h": fid["H0"] / 100, "z": 0.0, "ns": fid["ns"], "As": fid["As"],
        "k_piv": 0.05}
  for key in ("b", "d", "ax", "db", "m"):
    cd[f"Omega_{key}_0"] = cd[f"omega_{key}_0"] / cd["h"]**2
  cd["Omega_w_0"] = 1 - cd["Omega_m_0"]
  cd["G_a"] = func_D_z_unnorm_int(0.0, cd["Omega_m_0"], cd["Omega_w_0"])
  t2p = lambda T: lin_power_spectrum.transfer_to_PS(k, T, cd)
  cold_T = (cd["Omega_b_0"] * td[2, :, 0] + cd["Omega_d_0"] * td[1, :, 0]) \
      / cd["Omega_db_0"]
  pd = {"k": k, "power_total": t2p(td[6, :, 0]), "power_CDM": t2p(td[1, :, 0]),
        "power_baryon": t2p(td[2, :, 0]), "power_cold": t2p(cold_T),
        "power_axion": t2p(td[13, :, 0])}
  return cd, pd


for version in ("dome", "basic"):
  cd, pd = build_inputs(version)
  flags = _VERSION_FLAGS[version]
  M_arr = np.logspace(7, 18, 100)
  t0 = time.time()
  hmc = HMcode_params.HMCode_param_dic(cd, k, pd["power_cold"])
  t1 = time.time()
  axd = axion_params.func_axion_param_dic(
    M_arr, cd, pd, hmc, concentration_param=flags["concentration_param"])
  t2 = time.time()
  out = PS_nonlin_axion.func_full_halo_model_ax(
    M_arr, pd, cd, hmc, axd, **flags)
  t3 = time.time()
  print(f"{version:6s}: HMCode_param_dic = {t1-t0:6.3f} s | "
        f"axion_param_dic = {t2-t1:6.3f} s | full_halo_model = {t3-t2:6.3f} s"
        f" | total = {t3-t0:6.3f} s")

# ---- cProfile of one full dome row ------------------------------------------
print("\n" + "=" * 78)
print("cProfile, one dome redshift (top 22 by tottime)")
print("=" * 78)
prof = cProfile.Profile()
prof.enable()
_compute_row(payload("dome"))
prof.disable()
s = io.StringIO()
ps = pstats.Stats(prof, stream=s).sort_stats("tottime")
ps.print_stats(22)
lines = s.getvalue().split("\n")
for ln in lines[:34]:
  print(ln[:150])

# call counts of key primitives
print("=" * 78)
print("call counts of key primitives (one dome redshift)")
print("=" * 78)
s2 = io.StringIO()
ps2 = pstats.Stats(prof, stream=s2).sort_stats("ncalls")
ps2.print_stats("quad|brentq|sici|D_z_unnorm|dens_profile|halo_mass_function"
                "|variance|func_sigma|conc|simpson")
for ln in s2.getvalue().split("\n"):
  if any(w in ln for w in ("quad", "brentq", "sici", "D_z_unnorm",
                           "dens_profile", "halo_mass", "variance", "sigma",
                           "conc", "simpson")) and "/" in ln:
    print(ln[:150])
