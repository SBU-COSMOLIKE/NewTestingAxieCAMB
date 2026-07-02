"""Phase 1 part 3: end-to-end standalone lensed spectra with the axionHMcode
boost through ExternalNonLinearRatio (the prompt's pyCAMB snippet, made real).

Run A: NonLinear_none, full get_results          -> linear-lensing reference
Run C: transfers path + ratio == 1               -> code-path null (V0c)
Run B: transfers path + real axionHMcode boost   -> the production physics
"""
import sys
import time
import numpy as np

AXHM = "/Users/vivianmiranda/data/research/WayneHu/rayne/axionHMcode"
CAMB_PATH = "/Users/vivianmiranda/data/research/WayneHu/rayne/New_AxiECAMB"
sys.path.insert(0, AXHM)
sys.path.insert(0, CAMB_PATH)

import camb
from camb import model as cm
from camb.nonlinear import ExternalNonLinearRatio
from axionCAMB_and_lin_PS import lin_power_spectrum
from cosmology.overdensities import func_D_z_unnorm_int
from halo_model import HMcode_params
from axion_functions import axion_params
from halo_model import PS_nonlin_axion

h, omega_b, total_dark, ax_frac = 0.674, 0.02237, 0.12, 0.1
omega_ax = total_dark * ax_frac
omega_d = total_dark - omega_ax
As, ns, k_piv = 2.1e-9, 0.9655, 0.05
m_ax = 1e-25
LMAX = 2500

VERSION = "dome"
FLAGS = dict(alpha=True, eta_given=False, one_halo_damping=True,
             two_halo_damping=False, concentration_param=True, full_2h=False)


def make_pars(nonlinear):
  pars = camb.set_params(
    H0=100 * h, ombh2=omega_b, omch2=omega_d, omaxh2=omega_ax, m_ax=m_ax,
    As=As, ns=ns, num_massive_neutrinos=1, mnu=0.06, nnu=3.046,
    redshifts=[0.0], kmax=10.0, WantTransfer=True,
    lmax=LMAX, lens_potential_accuracy=1)
  pars.NonLinear = nonlinear
  return pars


def spectra(res):
  cl = res.get_lensed_scalar_cls(lmax=LMAX, CMB_unit="muK")
  pp = res.get_lens_potential_cls(lmax=2000)
  return cl, pp


def compute_boost_grid(res):
  """axionHMcode sqrt(P_NL/P_L_eq9) on the full transfer_redshifts grid."""
  zgrid = np.array(res.transfer_redshifts)
  order = np.argsort(zgrid)
  td = np.asarray(res.get_matter_transfer_data().transfer_data,
                  dtype=np.float64)
  k = td[cm.Transfer_kh - 1, :, 0]
  nz, nk = len(zgrid), len(k)
  sqrtB = np.ones((nz, nk))
  t0 = time.time()
  for row, iz in enumerate(order):
    z = float(zgrid[iz])
    cd = {
      "M_min": 7, "M_max": 18, "transfer_kmax": 100, "version": VERSION,
      "omega_b_0": omega_b, "omega_d_0": omega_d, "omega_ax_0": omega_ax,
      "omega_db_0": omega_d + omega_b, "omega_m_0": omega_b + total_dark,
      "m_ax": m_ax, "h": h, "z": z, "ns": ns, "As": As, "k_piv": k_piv,
    }
    for key in ["b", "d", "ax", "db", "m"]:
      cd[f"Omega_{key}_0"] = cd[f"omega_{key}_0"] / h**2
    cd["Omega_w_0"] = 1 - cd["Omega_m_0"]
    cd["G_a"] = func_D_z_unnorm_int(z, cd["Omega_m_0"], cd["Omega_w_0"])
    t2p = lambda T: lin_power_spectrum.transfer_to_PS(k, T, cd)
    T_c = td[cm.Transfer_cdm - 1, :, iz]
    T_b = td[cm.Transfer_b - 1, :, iz]
    T_a = td[cm.Transfer_axion - 1, :, iz]
    cold_T = (cd["Omega_b_0"] * T_b + cd["Omega_d_0"] * T_c) / cd["Omega_db_0"]
    pd = {"k": k, "power_total": t2p(td[cm.Transfer_tot - 1, :, iz]),
          "power_CDM": t2p(T_c), "power_baryon": t2p(T_b),
          "power_cold": t2p(cold_T), "power_axion": t2p(T_a)}
    M_arr = np.logspace(cd["M_min"], cd["M_max"], 100)
    hmc = HMcode_params.HMCode_param_dic(cd, k, pd["power_cold"])
    axd = axion_params.func_axion_param_dic(
      M_arr, cd, pd, hmc, concentration_param=FLAGS["concentration_param"])
    out = PS_nonlin_axion.func_full_halo_model_ax(
      M_arr, pd, cd, hmc, axd, **FLAGS)
    fc = axd["frac_cluster"]
    wdb = cd["Omega_db_0"] / cd["Omega_m_0"]
    wax = cd["Omega_ax_0"] / cd["Omega_m_0"]
    sqc, sqa = np.sqrt(pd["power_cold"]), np.sqrt(pd["power_axion"])
    sqrtB[row] = np.sqrt(np.asarray(out[0])) / (
      wdb * sqc + wax * (fc * sqc + (1 - fc) * sqa))
  dt = time.time() - t0
  z_sorted = zgrid[order]
  print(f"  boost grid: nz = {nz}, nk = {nk}, z in "
        f"[{z_sorted[0]:.2f}, {z_sorted[-1]:.2f}], took {dt:.1f} s "
        f"({dt/nz:.2f} s/z)")
  print(f"  sqrtB(z=0) at k~1: {sqrtB[0][np.argmin(np.abs(k-1)):][0]:.3f}; "
        f"sqrtB(z_max) range [{sqrtB[-1].min():.4f}, {sqrtB[-1].max():.4f}]")
  return k, z_sorted, sqrtB


print("RUN A: NonLinear_none reference (linear lensing)")
res_a = camb.get_results(make_pars(cm.NonLinear_none))
cl_a, pp_a = spectra(res_a)

print("RUN C: transfers path, ratio == 1 (code-path null, V0c)")
pars_c = make_pars(cm.NonLinear_both)
pars_c.NonLinearModel = ExternalNonLinearRatio()
res_c = camb.get_transfer_functions(pars_c, only_time_sources=True)
kk = np.logspace(-4, 2, 50)
zz = np.linspace(0, 12, 25)
res_c.Params.NonLinearModel.set_ratio(kk, zz, np.ones((len(zz), len(kk))))
res_c.power_spectra_from_transfer()
cl_c, pp_c = spectra(res_c)

print("RUN B: transfers path + axionHMcode dome boost")
pars_b = make_pars(cm.NonLinear_both)
pars_b.NonLinearModel = ExternalNonLinearRatio()
res_b = camb.get_transfer_functions(pars_b, only_time_sources=True)
k_h, z_grid, sqrtB = compute_boost_grid(res_b)
res_b.Params.NonLinearModel.set_ratio(k_h, z_grid, sqrtB)
res_b.power_spectra_from_transfer()
cl_b, pp_b = spectra(res_b)

print()
print("V0c (C vs A), lensed spectra rel. differences:")
for ell in [500, 1000, 2000]:
  print(f"  l={ell}: dTT = {cl_c[ell,0]/cl_a[ell,0]-1:+.2e}, "
        f"dEE = {cl_c[ell,1]/cl_a[ell,1]-1:+.2e}")
for L in [100, 500, 1000]:
  print(f"  L={L}: dClpp = {pp_c[L,0]/pp_a[L,0]-1:+.2e}")

print()
print("PHYSICS (B vs A): effect of the axionHMcode boost on lensed spectra:")
for ell in [500, 1000, 1500, 2000, 2400]:
  print(f"  l={ell}: dTT = {cl_b[ell,0]/cl_a[ell,0]-1:+.3%}, "
        f"dEE = {cl_b[ell,1]/cl_a[ell,1]-1:+.3%}")
for L in [100, 500, 1000, 2000]:
  print(f"  L={L}: dClpp = {pp_b[L,0]/pp_a[L,0]-1:+.3%}")
