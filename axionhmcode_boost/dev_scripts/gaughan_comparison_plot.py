"""Comparison figures: our AxiECAMB+axionHMcode pipeline vs approximate
read-offs from Gaughan/Green/Moss 2605.12054 Figs. 1-3.
Planck 2018 fiducial, fax = O_ax/O_D = 0.3.
"""
import sys
import numpy as np

AXHM = "/Users/vivianmiranda/data/research/WayneHu/rayne/axionHMcode"
CAMB_PATH = "/Users/vivianmiranda/data/research/WayneHu/rayne/New_AxiECAMB"
BOOST = f"{CAMB_PATH}/axionhmcode_boost"
for p in (AXHM, CAMB_PATH, BOOST):
  sys.path.insert(0, p)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import camb
from camb import model as cm
from camb.nonlinear import ExternalNonLinearRatio
from axionhmcode_boost import _compute_row, _VERSION_FLAGS

LMAX = 2500
OMD, FAX = 0.12011, 0.3
omaxh2 = OMD * FAX
omch2_ax = OMD - omaxh2
fid = dict(H0=67.32, ombh2=0.022383, As=2.101e-9, ns=0.96605,
           num_massive_neutrinos=1, mnu=0.06, nnu=3.046)
NUIS = {n: None for n in ("alpha_1", "alpha_2", "gamma_1", "gamma_2")}
COL = {"dome": "tab:red", "basic": "tab:blue"}


def payload(z, k, td, iz, m_ax, version):
  return {"z": z, "h": fid["H0"] / 100, "omega_b": fid["ombh2"],
          "omega_d": omch2_ax, "omega_ax": omaxh2, "m_ax": m_ax,
          "lcdm_mode": False, "ns": fid["ns"], "As": fid["As"],
          "k_piv": 0.05, "version": version,
          "flags": _VERSION_FLAGS[version],
          "M_arr": np.logspace(7, 18, 100), "nuisance": NUIS,
          "m_min_exponent": 7, "m_max_exponent": 18,
          "k": k, "T_cdm": td[1, :, iz], "T_b": td[2, :, iz],
          "T_ax": td[13, :, iz], "T_tot": td[6, :, iz]}


# ---------- Fig A data: P_NL ratio curves ----------
pars_l = camb.set_params(omch2=OMD, redshifts=[2.0, 0.0], kmax=20.0,
                         WantTransfer=True, lmax=LMAX,
                         lens_potential_accuracy=1, **fid)
pars_l.NonLinear = cm.NonLinear_pk
pars_l.NonLinearModel = camb.nonlinear.Halofit(halofit_version="mead2020")
res_l = camb.get_results(pars_l)
kh_l, z_l, pk_lcdm = res_l.get_nonlinear_matter_power_spectrum(
  hubble_units=True, k_hunit=True)

curves = {}
for m_ax in [1e-23, 1e-24, 1e-25]:
  pars_a = camb.set_params(omch2=omch2_ax, omaxh2=omaxh2, m_ax=m_ax,
                           redshifts=[2.0, 0.0], kmax=20.0,
                           WantTransfer=True, lmax=LMAX,
                           lens_potential_accuracy=1, **fid)
  pars_a.NonLinear = cm.NonLinear_none
  res_a = camb.get_results(pars_a)
  td = np.asarray(res_a.get_matter_transfer_data().transfer_data,
                  dtype=np.float64)
  k = td[0, :, 0]
  zs = np.array(res_a.Params.Transfer.PK_redshifts[
    :res_a.Params.Transfer.PK_num_redshifts])
  kh_a, z_a, pk_a_lin = res_a.get_linear_matter_power_spectrum(
    hubble_units=True, k_hunit=True)
  for version in ["dome", "basic"]:
    for z_t in [0.0, 2.0]:
      iz = int(np.argmin(np.abs(zs - z_t)))
      sqrtB = _compute_row(payload(z_t, k, td, iz, m_ax, version))
      p_nl_ax = sqrtB**2 * pk_a_lin[int(np.argmin(np.abs(z_a - z_t)))]
      p_nl_l = np.interp(k, kh_l, pk_lcdm[int(np.argmin(np.abs(z_l - z_t)))])
      curves[(m_ax, version, z_t)] = (k, p_nl_ax / p_nl_l)
  print(f"fig A: m={m_ax:g} done")

# ---------- Fig B data: lensed TT and Cpp differences, m = 1e-24 ----------
pars_ref = camb.set_params(omch2=OMD, redshifts=[0.0], kmax=10.0,
                           WantTransfer=True, lmax=LMAX,
                           lens_potential_accuracy=1, **fid)
pars_ref.NonLinear = cm.NonLinear_both
pars_ref.NonLinearModel = camb.nonlinear.Halofit(halofit_version="mead2020")
res_ref = camb.get_results(pars_ref)
cl_ref = res_ref.get_lensed_scalar_cls(lmax=LMAX, CMB_unit="muK")
pp_ref = res_ref.get_lens_potential_cls(lmax=2000)

cl_diff, pp_diff = {}, {}
for version in ["dome", "basic"]:
  pars_b = camb.set_params(omch2=omch2_ax, omaxh2=omaxh2, m_ax=1e-24,
                           redshifts=[0.0], kmax=10.0, WantTransfer=True,
                           lmax=LMAX, lens_potential_accuracy=1, **fid)
  pars_b.NonLinear = cm.NonLinear_both
  pars_b.NonLinearModel = ExternalNonLinearRatio()
  res_b = camb.get_transfer_functions(pars_b, only_time_sources=True)
  zg = np.array(res_b.transfer_redshifts)
  order = np.argsort(zg)
  td = np.asarray(res_b.get_matter_transfer_data().transfer_data,
                  dtype=np.float64)
  k = td[0, :, 0]
  rows = [_compute_row(payload(float(zg[iz]), k, td, iz, 1e-24, version))
          for iz in order]
  res_b.Params.NonLinearModel.set_ratio(k, zg[order], np.vstack(rows))
  res_b.power_spectra_from_transfer()
  cl_b = res_b.get_lensed_scalar_cls(lmax=LMAX, CMB_unit="muK")
  pp_b = res_b.get_lens_potential_cls(lmax=2000)
  ell = np.arange(2, LMAX)
  LL = np.arange(2, 2000)
  cl_diff[version] = (ell, 100 * (cl_b[2:LMAX, 0] / cl_ref[2:LMAX, 0] - 1))
  pp_diff[version] = (LL, 100 * (pp_b[2:2000, 0] / pp_ref[2:2000, 0] - 1))
  print(f"fig B: {version} done")

# ---------- approximate read-offs from Gaughan et al. figures ----------
fig1_anchors = {
  (1e-23, "dome", 2.0): [(0.3, 1.40), (10, 1.40)],
  (1e-23, "basic", 2.0): [(0.5, 0.88), (2, 0.63), (10, 0.78)],
  (1e-24, "dome", 2.0): [(0.4, 1.35), (0.8, 1.57), (2, 1.05), (7, 0.62)],
  (1e-24, "basic", 2.0): [(0.4, 0.78), (1.5, 0.47), (7, 0.44)],
  (1e-25, "dome", 2.0): [(0.25, 0.75), (0.4, 0.55)],
  (1e-25, "basic", 2.0): [(0.25, 0.62), (0.4, 0.45)],
  (1e-24, "dome", 0.0): [(1, 1.20), (3, 1.15)],
  (1e-24, "basic", 0.0): [(1, 0.88), (3, 0.85)],
  (1e-23, "dome", 0.0): [(2, 1.15)],
  (1e-23, "basic", 0.0): [(2, 0.93)],
  (1e-25, "basic", 0.0): [(1, 0.60)],
  (1e-25, "dome", 0.0): [(1, 0.85)],
}
fig3_anchors = {
  "dome": [(200, 4), (500, 13), (1000, 26), (2000, 38)],
  "basic": [(500, -4), (1000, -10), (2000, -19)],
}
fig2_anchors = {"dome": [(2000, 0.9), (2400, 1.5)],
                "basic": [(2000, -0.2), (2400, -0.4)]}

# ---------- Figure A ----------
fig, axes = plt.subplots(2, 3, figsize=(13, 6.6), sharex=True, sharey="row")
for col, m_ax in enumerate([1e-23, 1e-24, 1e-25]):
  for row, z_t in enumerate([0.0, 2.0]):
    ax = axes[row, col]
    for version in ["dome", "basic"]:
      k, r = curves[(m_ax, version, z_t)]
      mask = (k > 0.02) & (k < 20)
      ax.semilogx(k[mask], r[mask], color=COL[version], lw=1.8,
                  label=f"ours {version}" if (row, col) == (0, 0) else None)
      for kq, rq in fig1_anchors.get((m_ax, version, z_t), []):
        ax.plot(kq, rq, "o", mfc="none", mec=COL[version], ms=9, mew=1.8,
                label=("Gaughan+26 Fig.1 read-off"
                       if (row, col) == (0, 0) and version == "dome"
                       and (kq, rq) == fig1_anchors[(1e-23, "dome", 0.0)][0]
                       else None))
    ax.axhline(1, color="gray", lw=0.6, ls=":")
    ax.set_title(f"$m = 10^{{{int(np.log10(m_ax))}}}$ eV,  z = {z_t:.0f}",
                 fontsize=10)
    if row == 1:
      ax.set_xlabel(r"$k$  [$h\,{\rm Mpc}^{-1}$]")
    if col == 0:
      ax.set_ylabel(r"$P^{\rm axion}_{\rm NL}/P^{\Lambda{\rm CDM}}_{\rm NL}$")
    ax.set_ylim(0, 2.4)
h1 = [plt.Line2D([], [], color=COL["dome"], lw=1.8),
      plt.Line2D([], [], color=COL["basic"], lw=1.8),
      plt.Line2D([], [], marker="o", mfc="none", mec="k", ls="", ms=9)]
fig.legend(h1, ["ours: dome", "ours: basic",
                "Gaughan+26 Fig. 1 (approx. read-off, $\\pm\\sim$0.1)"],
           loc="center", bbox_to_anchor=(0.5, 0.955), ncol=3, frameon=False)
fig.suptitle(r"AxiECAMB + axionHMcode boost vs Gaughan/Green/Moss 2605.12054"
             r"  ($f_{\rm ax}=\Omega_{\rm ax}/\Omega_{\rm D}=0.3$)", y=1.005,
             fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.92))
out_a = "/private/tmp/claude-501/-Users-vivianmiranda-data-research-WayneHu-rayne/4e32f7cb-e470-4d65-8ad6-bf48eb7553b7/scratchpad/gaughan_comparison_fig1.png"
fig.savefig(out_a, dpi=150, bbox_inches="tight")

# ---------- Figure B ----------
fig, (axt, axp) = plt.subplots(1, 2, figsize=(11.5, 4.4))
for version in ["dome", "basic"]:
  ell, dtt = cl_diff[version]
  axt.plot(ell, dtt, color=COL[version], lw=1.4, label=f"ours {version}")
  for lq, vq in fig2_anchors[version]:
    axt.errorbar(lq, vq, yerr=0.4, fmt="o", mfc="none", mec=COL[version],
                 ecolor=COL[version], ms=8, mew=1.6, capsize=3)
axt.axhline(0, color="gray", lw=0.6, ls=":")
axt.set_xlabel(r"$\ell$")
axt.set_ylabel(r"$C_\ell^{TT,\rm axion}/C_\ell^{TT,\Lambda{\rm CDM}} - 1$ [%]")
axt.set_title(r"lensed TT,  $m=10^{-24}$ eV (their Fig. 2 centre)")
axt.legend(frameon=False, fontsize=9)

for version in ["dome", "basic"]:
  LL, dpp = pp_diff[version]
  axp.semilogx(LL, dpp, color=COL[version], lw=1.6, label=f"ours {version}")
  pts = np.array(fig3_anchors[version])
  axp.errorbar(pts[:, 0], pts[:, 1], yerr=3, fmt="o", mfc="none",
               mec=COL[version], ecolor=COL[version], ms=8, mew=1.6,
               capsize=3,
               label=(f"Gaughan+26 Fig. 3 read-off" if version == "dome"
                      else None))
axp.axhline(0, color="gray", lw=0.6, ls=":")
axp.set_xlabel(r"$L$")
axp.set_ylabel(r"$C_L^{\phi\phi,\rm axion}/C_L^{\phi\phi,\Lambda{\rm CDM}} - 1$ [%]")
axp.set_title(r"lensing potential,  $m=10^{-24}$ eV (their Fig. 3)")
axp.legend(frameon=False, fontsize=9)
fig.suptitle("markers: approximate read-offs from Gaughan/Green/Moss "
             "2605.12054 figures (error bars = reading uncertainty)",
             fontsize=10, y=1.02)
fig.tight_layout()
out_b = "/private/tmp/claude-501/-Users-vivianmiranda-data-research-WayneHu-rayne/4e32f7cb-e470-4d65-8ad6-bf48eb7553b7/scratchpad/gaughan_comparison_fig23.png"
fig.savefig(out_b, dpi=150, bbox_inches="tight")
print("saved", out_a, out_b)
