"""Phase 0: axionHMcode single-z benchmark, fed by New_AxiECAMB transfers.

Cosmology mirrors axionHMcode's input_file.txt: omega_b = 0.02237,
total dark = 0.12 with ax_fraction = 0.1 (omega_ax = 0.012, omega_cdm = 0.108),
m_ax = 1e-25 eV, h = 0.674. One redshift (z = 0), dome and basic versions.
"""
import sys
import time
import numpy as np

AXHM = "/Users/vivianmiranda/data/research/WayneHu/rayne/axionHMcode"
CAMB_PATH = "/Users/vivianmiranda/data/research/WayneHu/rayne/New_AxiECAMB"
sys.path.insert(0, AXHM)
sys.path.insert(0, CAMB_PATH)

import camb
from camb import model as camb_model
from axionCAMB_and_lin_PS import lin_power_spectrum
from cosmology.overdensities import func_D_z_unnorm_int
from halo_model import HMcode_params
from axion_functions import axion_params
from halo_model import PS_nonlin_axion

h = 0.674
omega_b, total_dark, ax_frac = 0.02237, 0.12, 0.1
omega_ax = total_dark * ax_frac
omega_d = total_dark - omega_ax
m_ax, As, ns, k_piv = 1e-25, 2.1e-9, 0.9655, 0.05
z_eval = 0.0

t0 = time.time()
pars = camb.set_params(
  H0=100 * h, ombh2=omega_b, omch2=omega_d, omaxh2=omega_ax, m_ax=m_ax,
  As=As, ns=ns, num_massive_neutrinos=1, mnu=0.06, nnu=3.046,
  redshifts=[z_eval], kmax=100.0, WantTransfer=True,
  NonLinear=camb_model.NonLinear_none)
results = camb.get_results(pars)
td = results.get_matter_transfer_data()
print(f"AxiECAMB run: {time.time()-t0:.2f} s; "
      f"transfer_data shape {np.asarray(td.transfer_data).shape}")
print(f"axion regime: is_de_like = {results.Params.Axion.is_de_like}, "
      f"a_osc = {results.Params.Axion.a_osc:.3e} "
      f"(z_osc = {1/results.Params.Axion.a_osc - 1:.3e})")

tdata = np.asarray(td.transfer_data, dtype=np.float64)
k = tdata[camb_model.Transfer_kh - 1, :, 0]
T_cdm = tdata[camb_model.Transfer_cdm - 1, :, 0]
T_b = tdata[camb_model.Transfer_b - 1, :, 0]
T_ax = tdata[camb_model.Transfer_axion - 1, :, 0]
T_tot = tdata[camb_model.Transfer_tot - 1, :, 0]

cosmo_dic = {
  "M_min": 7, "M_max": 18, "transfer_kmax": 100,
  "omega_b_0": omega_b, "omega_d_0": omega_d, "omega_ax_0": omega_ax,
  "omega_db_0": omega_d + omega_b, "omega_m_0": omega_b + total_dark,
  "m_ax": m_ax, "h": h, "z": z_eval, "ns": ns, "As": As, "k_piv": k_piv,
}
for key in ["b", "d", "ax", "db", "m"]:
  cosmo_dic[f"Omega_{key}_0"] = cosmo_dic[f"omega_{key}_0"] / h**2
cosmo_dic["Omega_w_0"] = 1 - cosmo_dic["Omega_m_0"]
cosmo_dic["G_a"] = func_D_z_unnorm_int(
  z_eval, cosmo_dic["Omega_m_0"], cosmo_dic["Omega_w_0"])

power_spec_dic = {
  "k": k,
  "power_total": lin_power_spectrum.transfer_to_PS(k, T_tot, cosmo_dic),
  "power_CDM": lin_power_spectrum.transfer_to_PS(k, T_cdm, cosmo_dic),
  "power_baryon": lin_power_spectrum.transfer_to_PS(k, T_b, cosmo_dic),
  "power_cold": lin_power_spectrum.transfer_to_PS(
    k, (cosmo_dic["Omega_b_0"] * T_b + cosmo_dic["Omega_d_0"] * T_cdm)
    / cosmo_dic["Omega_db_0"], cosmo_dic),
  "power_axion": lin_power_spectrum.transfer_to_PS(k, T_ax, cosmo_dic),
}

kh_camb, z_camb, pk_camb = results.get_linear_matter_power_spectrum(
  hubble_units=True, k_hunit=True)
pk_interp = np.interp(k, kh_camb, pk_camb[0])
dev = np.abs(power_spec_dic["power_total"] / pk_interp - 1)
print(f"V1 check (P_tot from transfers vs CAMB linear P): "
      f"median dev {np.median(dev):.2e}, max {dev.max():.2e}")

M_arr = np.logspace(cosmo_dic["M_min"], cosmo_dic["M_max"], 100)

for version, flags in [
    ("dome", dict(alpha=True, eta_given=False, one_halo_damping=True,
                  two_halo_damping=False, concentration_param=True,
                  full_2h=False)),
    ("basic", dict(alpha=False, eta_given=False, one_halo_damping=True,
                   two_halo_damping=False, concentration_param=False,
                   full_2h=True))]:
  cosmo_dic["version"] = version
  times = []
  for rep in range(3):
    t0 = time.time()
    hmcode_dic = HMcode_params.HMCode_param_dic(
      cosmo_dic, power_spec_dic["k"], power_spec_dic["power_cold"])
    axion_dic = axion_params.func_axion_param_dic(
      M_arr, cosmo_dic, power_spec_dic, hmcode_dic,
      concentration_param=flags["concentration_param"])
    PS_nl = PS_nonlin_axion.func_full_halo_model_ax(
      M_arr, power_spec_dic, cosmo_dic, hmcode_dic, axion_dic, **flags)
    times.append(time.time() - t0)
  P_NL = np.asarray(PS_nl[0])
  boost = P_NL / power_spec_dic["power_total"]
  i01 = np.argmin(np.abs(k - 0.1))
  i1 = np.argmin(np.abs(k - 1.0))
  i5 = np.argmin(np.abs(k - 5.0))
  print(f"{version:6s}: t_first = {times[0]:.2f} s (JIT), "
        f"t_warm = {min(times[1:]):.2f} s | "
        f"B(k=0.1) = {boost[i01]:.3f}, B(1) = {boost[i1]:.3f}, "
        f"B(5) = {boost[i5]:.3f}")

print("nz for production lensing grid ~50-75 -> per-eval cost estimate:")
print(f"  ~{50*min(times[1:]):.0f}-{75*min(times[1:]):.0f} s/likelihood-eval "
      f"at warm single-z cost {min(times[1:]):.2f} s")
