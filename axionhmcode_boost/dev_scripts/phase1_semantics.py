"""Phase 1 part 1-2: transfer-variable semantics per regime (R3), KG-phase
axion transfer (R10 lead), Eq.9 boost construction with low-k -> 1 check,
and kmax-truncation sensitivity.
"""
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
from halo_model import HMcode_params
from axion_functions import axion_params
from halo_model import PS_nonlin_axion

h, omega_b, total_dark, ax_frac = 0.674, 0.02237, 0.12, 0.1
omega_ax = total_dark * ax_frac
omega_d = total_dark - omega_ax
As, ns, k_piv = 2.1e-9, 0.9655, 0.05
omega_nu = 0.06 * (3.046 / 3) ** 0.75 / 94.0708


def run_axiecamb(m_ax, redshifts, kmax=100.0):
  pars = camb.set_params(
    H0=100 * h, ombh2=omega_b, omch2=omega_d, omaxh2=omega_ax, m_ax=m_ax,
    As=As, ns=ns, num_massive_neutrinos=1, mnu=0.06, nnu=3.046,
    redshifts=sorted(redshifts, reverse=True), kmax=kmax, WantTransfer=True,
    NonLinear=cm.NonLinear_none)
  return camb.get_results(pars)


def transfers_at(results, z):
  td = np.asarray(results.get_matter_transfer_data().transfer_data,
                  dtype=np.float64)
  npk = results.Params.Transfer.PK_num_redshifts
  zs = np.array(results.Params.Transfer.PK_redshifts[:npk])
  iz = np.argmin(np.abs(zs - z))
  assert abs(zs[iz] - z) < 1e-8, (zs, z)
  out = {name: td[getattr(cm, f"Transfer_{name}") - 1, :, iz]
         for name in ["kh", "cdm", "b", "axion", "tot", "nonu", "nu"]}
  return out


print("=" * 72)
print("PART 1 (R3): what do Transfer_tot / Transfer_nonu contain, per regime?")
print("=" * 72)
for m_ax, label in [(1e-25, "DM-like"), (1e-33, "DE-like")]:
  res = run_axiecamb(m_ax, [0.0])
  ax = res.Params.Axion
  tr = transfers_at(res, 0.0)
  # candidate compositions (density-weighted transfer functions)
  cb = (omega_d * tr["cdm"] + omega_b * tr["b"]) / (omega_d + omega_b)
  cba = (omega_d * tr["cdm"] + omega_b * tr["b"] + omega_ax * tr["axion"]) / (
    omega_d + omega_b + omega_ax)
  cban = (omega_d * tr["cdm"] + omega_b * tr["b"] + omega_ax * tr["axion"]
          + omega_nu * tr["nu"]) / (omega_d + omega_b + omega_ax + omega_nu)
  cbn = (omega_d * tr["cdm"] + omega_b * tr["b"] + omega_nu * tr["nu"]) / (
    omega_d + omega_b + omega_nu)

  def match(name, cand, ref):
    dev = np.median(np.abs(cand / ref - 1))
    return f"{name}: {dev:.2e}"

  print(f"\nm_ax = {m_ax:g} eV ({label}): is_de_like = {ax.is_de_like}, "
        f"a_osc = {ax.a_osc:.3e}")
  print("  Transfer_tot  vs " + " | ".join([
    match("cb", cb, tr["tot"]), match("cb+ax", cba, tr["tot"]),
    match("cb+ax+nu", cban, tr["tot"]), match("cb+nu", cbn, tr["tot"])]))
  print("  Transfer_nonu vs " + " | ".join([
    match("cb", cb, tr["nonu"]), match("cb+ax", cba, tr["nonu"])]))
  ka = tr["kh"]
  for kq in [0.01, 0.1, 1.0, 10.0]:
    i = np.argmin(np.abs(ka - kq))
    print(f"  k={ka[i]:7.3f} h/Mpc: T_ax/T_cdm = {tr['axion'][i]/tr['cdm'][i]:+.4e}")

print()
print("=" * 72)
print("PART 2: Eq.9 boost with model-linear-limit denominator; low-k -> 1")
print("=" * 72)


def make_cosmo_dic(z, version):
  d = {
    "M_min": 7, "M_max": 18, "transfer_kmax": 100, "version": version,
    "omega_b_0": omega_b, "omega_d_0": omega_d, "omega_ax_0": omega_ax,
    "omega_db_0": omega_d + omega_b, "omega_m_0": omega_b + total_dark,
    "m_ax": 1e-25, "h": h, "z": z, "ns": ns, "As": As, "k_piv": k_piv,
  }
  for key in ["b", "d", "ax", "db", "m"]:
    d[f"Omega_{key}_0"] = d[f"omega_{key}_0"] / h**2
  d["Omega_w_0"] = 1 - d["Omega_m_0"]
  d["G_a"] = func_D_z_unnorm_int(z, d["Omega_m_0"], d["Omega_w_0"])
  return d


def make_power_dic(tr, cosmo_dic):
  k = tr["kh"]
  t2p = lambda T: lin_power_spectrum.transfer_to_PS(k, T, cosmo_dic)
  cold_T = (cosmo_dic["Omega_b_0"] * tr["b"] + cosmo_dic["Omega_d_0"]
            * tr["cdm"]) / cosmo_dic["Omega_db_0"]
  return {"k": k, "power_total": t2p(tr["tot"]), "power_CDM": t2p(tr["cdm"]),
          "power_baryon": t2p(tr["b"]), "power_cold": t2p(cold_T),
          "power_axion": t2p(tr["axion"])}


VERSION_FLAGS = {
  "dome": dict(alpha=True, eta_given=False, one_halo_damping=True,
               two_halo_damping=False, concentration_param=True,
               full_2h=False),
  "basic": dict(alpha=False, eta_given=False, one_halo_damping=True,
                two_halo_damping=False, concentration_param=False,
                full_2h=True),
}


def boost_at_z(tr, z, version):
  cd = make_cosmo_dic(z, version)
  pd = make_power_dic(tr, cd)
  flags = VERSION_FLAGS[version]
  M_arr = np.logspace(cd["M_min"], cd["M_max"], 100)
  hmc = HMcode_params.HMCode_param_dic(cd, pd["k"], pd["power_cold"])
  axd = axion_params.func_axion_param_dic(
    M_arr, cd, pd, hmc, concentration_param=flags["concentration_param"])
  out = PS_nonlin_axion.func_full_halo_model_ax(
    M_arr, pd, cd, hmc, axd, **flags)
  P_NL = np.asarray(out[0])
  fc = axd["frac_cluster"]
  wdb = cd["Omega_db_0"] / cd["Omega_m_0"]
  wax = cd["Omega_ax_0"] / cd["Omega_m_0"]
  sqc, sqa = np.sqrt(pd["power_cold"]), np.sqrt(pd["power_axion"])
  sqrt_PL_eq9 = wdb * sqc + wax * (fc * sqc + (1 - fc) * sqa)
  sqrt_B = np.sqrt(P_NL) / sqrt_PL_eq9
  return pd["k"], sqrt_B, fc


res = run_axiecamb(1e-25, [0.0, 2.0])
for z in [0.0, 2.0]:
  tr = transfers_at(res, z)
  for version in ["dome", "basic"]:
    k, sB, fc = boost_at_z(tr, z, version)
    B = sB**2
    idx = {kq: np.argmin(np.abs(k - kq)) for kq in [1e-3, 1e-2, 0.1, 0.5, 1, 5]}
    print(f"z={z:.0f} {version:6s} frac_cluster={float(fc):.4f} | "
          + " ".join(f"B({kq:g})={B[i]:.4f}" for kq, i in idx.items()))

print()
print("=" * 72)
print("PART 2b: kmax-truncation sensitivity (k <= 10 vs k <= 100 h/Mpc inputs)")
print("=" * 72)
tr0 = transfers_at(res, 0.0)
k_full, sB_full, _ = boost_at_z(tr0, 0.0, "dome")
tr_cut = {key: (val[tr0["kh"] <= 10.0] if key != "kh" else val[tr0["kh"] <= 10.0])
          for key, val in tr0.items()}
k_cut, sB_cut, _ = boost_at_z(tr_cut, 0.0, "dome")
B_full_on_cut = np.interp(k_cut, k_full, sB_full**2)
dev = np.abs(sB_cut**2 / B_full_on_cut - 1)
for kq in [0.1, 0.5, 1.0, 3.0, 8.0]:
  i = np.argmin(np.abs(k_cut - kq))
  print(f"k={k_cut[i]:6.2f}: boost rel. change from k-truncation = {dev[i]:.3e}")
