"""Phase 4b (V3): Gaughan/Green/Moss-style comparison.
Planck 2018 fiducial, fax = O_ax/O_D = 0.3 (their eq. 1), masses
1e-23/-24/-25 eV. Part 1: P_NL^axion / P_NL^LCDM ratio curves at z = 0, 2
(their Fig. 1 lower panels). Part 2: lensed TT/EE and Cpp differences vs
LCDM for m = 1e-24 (their Figs. 2-3), dome and basic.
"""
import sys
import numpy as np

AXHM = "/Users/vivianmiranda/data/research/WayneHu/rayne/axionHMcode"
CAMB_PATH = "/Users/vivianmiranda/data/research/WayneHu/rayne/New_AxiECAMB"
BOOST = f"{CAMB_PATH}/axionhmcode_boost"
for p in (AXHM, CAMB_PATH, BOOST):
  sys.path.insert(0, p)

import camb
from camb import model as cm
from camb.nonlinear import ExternalNonLinearRatio
from axionhmcode_boost import _compute_row, _VERSION_FLAGS

LMAX = 2500
OMD = 0.12011
FAX = 0.3
fid = dict(H0=67.32, ombh2=0.022383, As=2.101e-9, ns=0.96605,
           num_massive_neutrinos=1, mnu=0.06, nnu=3.046,
           lmax=LMAX, lens_potential_accuracy=1)
NUIS = {n: None for n in ("alpha_1", "alpha_2", "gamma_1", "gamma_2")}


def payload(z, k, td, iz, m_ax, omch2, omaxh2, version):
  return {"z": z, "h": fid["H0"] / 100, "omega_b": fid["ombh2"],
          "omega_d": omch2, "omega_ax": omaxh2, "m_ax": m_ax,
          "lcdm_mode": False, "ns": fid["ns"], "As": fid["As"],
          "k_piv": 0.05, "version": version,
          "flags": _VERSION_FLAGS[version],
          "M_arr": np.logspace(7, 18, 100), "nuisance": NUIS,
          "m_min_exponent": 7, "m_max_exponent": 18,
          "k": k, "T_cdm": td[1, :, iz], "T_b": td[2, :, iz],
          "T_ax": td[13, :, iz], "T_tot": td[6, :, iz]}


omaxh2 = OMD * FAX
omch2_ax = OMD - omaxh2

print("=" * 72)
print("V3 part 1: P_NL^axion / P_NL^LCDM at z = 0, 2 (Gaughan Fig. 1 style)")
print(f"  fax = O_ax/O_D = {FAX}, omaxh2 = {omaxh2:.4f}")
print("=" * 72)
pars_l = camb.set_params(omch2=OMD, redshifts=[2.0, 0.0], kmax=20.0,
                         WantTransfer=True, **fid)
pars_l.NonLinear = cm.NonLinear_pk
pars_l.NonLinearModel = camb.nonlinear.Halofit(halofit_version="mead2020")
res_l = camb.get_results(pars_l)
kh_l, z_l, pk_lcdm = res_l.get_nonlinear_matter_power_spectrum(
  hubble_units=True, k_hunit=True)

for m_ax in [1e-23, 1e-24, 1e-25]:
  pars_a = camb.set_params(omch2=omch2_ax, omaxh2=omaxh2, m_ax=m_ax,
                           redshifts=[2.0, 0.0], kmax=20.0,
                           WantTransfer=True, **fid)
  pars_a.NonLinear = cm.NonLinear_none
  res_a = camb.get_results(pars_a)
  td = np.asarray(res_a.get_matter_transfer_data().transfer_data,
                  dtype=np.float64)
  k = td[0, :, 0]
  zs = np.array(res_a.Params.Transfer.PK_redshifts[
    :res_a.Params.Transfer.PK_num_redshifts])
  for version in ["dome", "basic"]:
    line = f"m={m_ax:0.0e} {version:6s}"
    for z_t in [0.0, 2.0]:
      iz = int(np.argmin(np.abs(zs - z_t)))
      pl = payload(z_t, k, td, iz, m_ax, omch2_ax, omaxh2, version)
      sqrtB = _compute_row(pl)
      # P_NL^axion in the Eq.9 decomposition = (sqrtB * sqrt(P_L_eq9))^2;
      # rebuild P_L_eq9 exactly as _compute_row does
      from axionCAMB_and_lin_PS import lin_power_spectrum
      from cosmology.overdensities import func_D_z_unnorm_int
      cd = {"omega_b_0": fid["ombh2"], "omega_d_0": omch2_ax,
            "omega_ax_0": omaxh2, "h": fid["H0"] / 100, "z": z_t,
            "ns": fid["ns"], "As": fid["As"], "k_piv": 0.05}
      cd["omega_db_0"] = cd["omega_b_0"] + cd["omega_d_0"]
      cd["omega_m_0"] = cd["omega_db_0"] + cd["omega_ax_0"]
      for key in ["b", "d", "ax", "db", "m"]:
        cd[f"Omega_{key}_0"] = cd[f"omega_{key}_0"] / cd["h"]**2
      # frac_cluster is inside _compute_row; instead of re-deriving, use the
      # CAMB linear total as the anchor: P_NL_ax ~= B * P_L_tot is NOT exact
      # (Eq.9 vs CAMB composition differ at the sub-percent level at high k),
      # but for the Fig.1-style ratio it is the same convention CAMB applies.
      iz_pk = int(np.argmin(np.abs(z_l - z_t)))
      kh_a, z_a, pk_a_lin = res_a.get_linear_matter_power_spectrum(
        hubble_units=True, k_hunit=True)
      p_nl_ax = sqrtB**2 * pk_a_lin[int(np.argmin(np.abs(z_a - z_t)))]
      p_nl_l = np.interp(k, kh_l, pk_lcdm[iz_pk])
      r = p_nl_ax / p_nl_l
      vals = " ".join(
        f"R({kq:g},z={z_t:.0f})={np.interp(kq, k, r):.3f}"
        for kq in [0.3, 1.0, 3.0])
      line += "  " + vals
    print(line)

print()
print("=" * 72)
print("V3 part 2: lensed TT/EE and Cpp vs LCDM, m = 1e-24 (Figs. 2-3 style)")
print("=" * 72)
pars_ref = camb.set_params(omch2=OMD, redshifts=[0.0], kmax=10.0,
                           WantTransfer=True, **fid)
pars_ref.NonLinear = cm.NonLinear_both
pars_ref.NonLinearModel = camb.nonlinear.Halofit(halofit_version="mead2020")
res_ref = camb.get_results(pars_ref)
cl_ref = res_ref.get_lensed_scalar_cls(lmax=LMAX, CMB_unit="muK")
pp_ref = res_ref.get_lens_potential_cls(lmax=2000)

for version in ["dome", "basic"]:
  pars_b = camb.set_params(omch2=omch2_ax, omaxh2=omaxh2, m_ax=1e-24,
                           redshifts=[0.0], kmax=10.0, WantTransfer=True,
                           **fid)
  pars_b.NonLinear = cm.NonLinear_both
  pars_b.NonLinearModel = ExternalNonLinearRatio()
  res_b = camb.get_transfer_functions(pars_b, only_time_sources=True)
  zg = np.array(res_b.transfer_redshifts)
  order = np.argsort(zg)
  td = np.asarray(res_b.get_matter_transfer_data().transfer_data,
                  dtype=np.float64)
  k = td[0, :, 0]
  rows = [_compute_row(payload(float(zg[iz]), k, td, iz, 1e-24, omch2_ax,
                               omaxh2, version)) for iz in order]
  res_b.Params.NonLinearModel.set_ratio(k, zg[order], np.vstack(rows))
  res_b.power_spectra_from_transfer()
  cl_b = res_b.get_lensed_scalar_cls(lmax=LMAX, CMB_unit="muK")
  pp_b = res_b.get_lens_potential_cls(lmax=2000)
  tt = " ".join(f"dTT({l})={cl_b[l,0]/cl_ref[l,0]-1:+.3%}"
                for l in [1000, 2000, 2400])
  ee = " ".join(f"dEE({l})={cl_b[l,1]/cl_ref[l,1]-1:+.3%}"
                for l in [1000, 2000])
  pp = " ".join(f"dCpp({L})={pp_b[L,0]/pp_ref[L,0]-1:+.3%}"
                for L in [100, 500, 1000, 2000])
  print(f"{version:6s}: {tt}\n        {ee}\n        {pp}")
