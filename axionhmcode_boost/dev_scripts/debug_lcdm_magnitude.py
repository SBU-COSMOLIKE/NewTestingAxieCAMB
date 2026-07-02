"""Debug: axionHMcode cold-model LCDM boost vs CAMB HMcode-2020 (mead2020)."""
import sys
import numpy as np

AXHM = "/Users/vivianmiranda/data/research/WayneHu/rayne/axionHMcode"
CAMB_PATH = "/Users/vivianmiranda/data/research/WayneHu/rayne/New_AxiECAMB"
sys.path.insert(0, AXHM)
sys.path.insert(0, CAMB_PATH)

import camb
from camb import model as cm
from axionCAMB_and_lin_PS import lin_power_spectrum
from cosmology.overdensities import func_D_z_unnorm_int
from halo_model import HMcode_params, PS_nonlin_cold

h = 0.6732
pars = camb.set_params(
  H0=100 * h, ombh2=0.022383, omch2=0.12011, As=2.101e-9, ns=0.96605,
  num_massive_neutrinos=1, mnu=0.06, nnu=3.046,
  redshifts=[0.0], kmax=100.0, WantTransfer=True,
  NonLinear=cm.NonLinear_pk, halofit_version="mead2020")
res = camb.get_results(pars)
kh_l, _, pk_l = res.get_linear_matter_power_spectrum(hubble_units=True,
                                                     k_hunit=True)
kh_n, _, pk_n = res.get_nonlinear_matter_power_spectrum(hubble_units=True,
                                                        k_hunit=True)
R_camb = pk_n[0] / pk_l[0]
print(f"sigma8 (CAMB) = {res.get_sigma8_0():.4f}")

td = np.asarray(res.get_matter_transfer_data().transfer_data,
                dtype=np.float64)
k = td[cm.Transfer_kh - 1, :, 0]
T_c, T_b = td[cm.Transfer_cdm - 1, :, 0], td[cm.Transfer_b - 1, :, 0]
T_tot = td[cm.Transfer_tot - 1, :, 0]

cd = {"M_min": 7, "M_max": 18, "transfer_kmax": 100, "version": "dome",
      "omega_b_0": 0.022383, "omega_d_0": 0.12011, "omega_ax_0": 1e-20 * h**2,
      "m_ax": 1e-25, "h": h, "z": 0.0, "ns": 0.96605, "As": 2.101e-9,
      "k_piv": 0.05}
cd["omega_db_0"] = cd["omega_b_0"] + cd["omega_d_0"]
cd["omega_m_0"] = cd["omega_db_0"] + cd["omega_ax_0"]
for key in ["b", "d", "ax", "db", "m"]:
  cd[f"Omega_{key}_0"] = cd[f"omega_{key}_0"] / h**2
cd["Omega_w_0"] = 1 - cd["Omega_m_0"]
cd["G_a"] = func_D_z_unnorm_int(0.0, cd["Omega_m_0"], cd["Omega_w_0"])

t2p = lambda T: lin_power_spectrum.transfer_to_PS(k, T, cd)
cold_T = (cd["Omega_b_0"] * T_b + cd["Omega_d_0"] * T_c) / cd["Omega_db_0"]
P_cold, P_tot = t2p(cold_T), t2p(T_tot)

# sigma8 from our linear power (top-hat, R = 8 Mpc/h)
x = k * 8.0
W = 3 * (np.sin(x) - x * np.cos(x)) / x**3
sig8 = np.sqrt(np.trapz(P_tot * W**2 * k**2, k) / (2 * np.pi**2))
print(f"sigma8 (our P_tot integral) = {sig8:.4f}")

M_arr = np.logspace(7, 18, 100)
hmc = HMcode_params.HMCode_param_dic(cd, k, P_cold)
print("hmcode_dic:", {key: (float(val[0]) if hasattr(val, '__len__')
                            else round(float(val), 4))
                      for key, val in hmc.items()})

combos = {
  "dome-like  (a=T,1h=T,2h=F,cp=T,f2=F)": dict(
    alpha=True, one_halo_damping=True, two_halo_damping=False,
    concentration_param=True, full_2h=False),
  "basic-like (a=F,1h=T,2h=F,cp=F,f2=T)": dict(
    alpha=False, one_halo_damping=True, two_halo_damping=False,
    concentration_param=False, full_2h=True),
  "hm2020ish  (a=T,1h=T,2h=T,cp=T,f2=F)": dict(
    alpha=True, one_halo_damping=True, two_halo_damping=True,
    concentration_param=True, full_2h=False),
}
iq = {kq: np.argmin(np.abs(k - kq)) for kq in [0.2, 0.5, 1.0, 2.0, 5.0]}
print("\nk:          " + "  ".join(f"{kq:>7.1f}" for kq in iq))
print("R mead2020: " + "  ".join(
  f"{np.interp(k[i], kh_n, R_camb):7.3f}" for kq, i in iq.items()))
for label, fl in combos.items():
  out = PS_nonlin_cold.func_non_lin_PS_matter(
    M_arr, k, P_cold, cd, hmc, cd["Omega_db_0"], eta_given=False,
    ax_one_halo=False, axion_dic=None, **fl)
  R = np.asarray(out[0]) / P_cold
  print(f"{label}: " + "  ".join(f"{R[i]:7.3f}" for kq, i in iq.items()))
