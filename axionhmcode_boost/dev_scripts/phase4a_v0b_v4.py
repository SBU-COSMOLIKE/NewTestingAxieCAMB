"""Phase 4a: V0b plumbing null (external mead2020 ratio must reproduce
internal mead2020 lensed spectra) and V4 cross-check (boost pipeline vs the
hmcode-branch axion-aware internal HMcode).
"""
import sys
import numpy as np

AXHM = "/Users/vivianmiranda/data/research/WayneHu/rayne/axionHMcode"
CAMB_PATH = "/Users/vivianmiranda/data/research/WayneHu/rayne/New_AxiECAMB"
BOOST = f"{CAMB_PATH}/axionhmcode_boost"
sys.path.insert(0, AXHM)
sys.path.insert(0, CAMB_PATH)
sys.path.insert(0, BOOST)

import camb
from camb import model as cm
from camb.nonlinear import ExternalNonLinearRatio
from axionhmcode_boost import _compute_row, _VERSION_FLAGS

LMAX = 2500
fid = dict(H0=67.32, ombh2=0.022383, As=2.101e-9, ns=0.96605,
           num_massive_neutrinos=1, mnu=0.06, nnu=3.046,
           lmax=LMAX, lens_potential_accuracy=1)


def lensed(res):
  cl = res.get_lensed_scalar_cls(lmax=LMAX, CMB_unit="muK")
  pp = res.get_lens_potential_cls(lmax=2000)
  return cl, pp


def report(tag, cl_a, pp_a, cl_b, pp_b, ells=(500, 1000, 2000),
           Ls=(100, 500, 1000)):
  parts = [f"dTT(l={l}) = {cl_b[l,0]/cl_a[l,0]-1:+.2e}" for l in ells]
  parts += [f"dCpp(L={L}) = {pp_b[L,0]/pp_a[L,0]-1:+.2e}" for L in Ls]
  print(f"{tag}: " + ", ".join(parts))


print("=" * 72)
print("V0b: external mead2020 ratio through set_ratio vs internal mead2020")
print("=" * 72)
zs = list(np.linspace(0, 10, 26))
pars_ref = camb.set_params(omch2=0.12011, redshifts=sorted(zs, reverse=True),
                           kmax=10.0, WantTransfer=True, **fid)
pars_ref.NonLinear = cm.NonLinear_both
pars_ref.NonLinearModel = camb.nonlinear.Halofit(halofit_version="mead2020")
res_ref = camb.get_results(pars_ref)
cl_ref, pp_ref = lensed(res_ref)
kh, z_pk, pk_l = res_ref.get_linear_matter_power_spectrum(hubble_units=True,
                                                          k_hunit=True)
_, _, pk_n = res_ref.get_nonlinear_matter_power_spectrum(hubble_units=True,
                                                         k_hunit=True)
sqrt_ratio = np.sqrt(pk_n / pk_l)

pars_ext = camb.set_params(omch2=0.12011, redshifts=sorted(zs, reverse=True),
                           kmax=10.0, WantTransfer=True, **fid)
pars_ext.NonLinear = cm.NonLinear_both
pars_ext.NonLinearModel = ExternalNonLinearRatio()
res_ext = camb.get_transfer_functions(pars_ext, only_time_sources=True)
res_ext.Params.NonLinearModel.set_ratio(kh, z_pk, sqrt_ratio)
res_ext.power_spectra_from_transfer()
cl_ext, pp_ext = lensed(res_ext)
report("V0b (ext vs int mead2020)", cl_ref, pp_ref, cl_ext, pp_ext)

print()
print("=" * 72)
print("V4: boost pipeline vs hmcode-branch axion-aware internal HMcode")
print("    (axion fax_dm = 0.1, m = 1e-25 eV; different halo models -> ")
print("     expect few-percent-level Cpp differences, not agreement)")
print("=" * 72)
ax = dict(omch2=0.108, omaxh2=0.012, m_ax=1e-25)

pars_int = camb.set_params(redshifts=[0.0], kmax=10.0, WantTransfer=True,
                           **fid, **ax)
pars_int.NonLinear = cm.NonLinear_both
pars_int.NonLinearModel = camb.nonlinear.Halofit(halofit_version="mead2020")
res_int = camb.get_results(pars_int)
cl_int, pp_int = lensed(res_int)

pars_b = camb.set_params(redshifts=[0.0], kmax=10.0, WantTransfer=True,
                         **fid, **ax)
pars_b.NonLinear = cm.NonLinear_both
pars_b.NonLinearModel = ExternalNonLinearRatio()
res_b = camb.get_transfer_functions(pars_b, only_time_sources=True)

zgrid = np.array(res_b.transfer_redshifts)
order = np.argsort(zgrid)
td = np.asarray(res_b.get_matter_transfer_data().transfer_data,
                dtype=np.float64)
k = td[0, :, 0]
h = fid["H0"] / 100
rows = []
for iz in order:
  rows.append(_compute_row({
    "z": float(zgrid[iz]), "h": h, "omega_b": fid["ombh2"],
    "omega_d": ax["omch2"], "omega_ax": ax["omaxh2"], "m_ax": ax["m_ax"],
    "lcdm_mode": False, "ns": fid["ns"], "As": fid["As"], "k_piv": 0.05,
    "version": "dome", "flags": _VERSION_FLAGS["dome"],
    "M_arr": np.logspace(7, 18, 100),
    "nuisance": {n: None for n in
                 ("alpha_1", "alpha_2", "gamma_1", "gamma_2")},
    "m_min_exponent": 7, "m_max_exponent": 18,
    "k": k, "T_cdm": td[1, :, iz], "T_b": td[2, :, iz],
    "T_ax": td[13, :, iz], "T_tot": td[6, :, iz]}))
res_b.Params.NonLinearModel.set_ratio(k, zgrid[order], np.vstack(rows))
res_b.power_spectra_from_transfer()
cl_b, pp_b = lensed(res_b)
report("V4 (boost vs internal axion-HMcode)", cl_int, pp_int, cl_b, pp_b)
