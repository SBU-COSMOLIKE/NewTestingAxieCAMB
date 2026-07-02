"""Phase 4c (V6): boost-grid numerical convergence at the representative
point m = 1e-24 eV, fax = O_ax/O_D = 0.3 (dome), z = 0 and 2.
Variations: halo-mass grid density (100 -> 320), mass-range widening
(7..18 -> 6..19), and input-k coarsening (every 2nd point).
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
from axionhmcode_boost import _compute_row, _VERSION_FLAGS

fid = dict(H0=67.32, ombh2=0.022383, As=2.101e-9, ns=0.96605,
           num_massive_neutrinos=1, mnu=0.06, nnu=3.046)
OMD, FAX, M_AX = 0.12011, 0.3, 1e-24
omaxh2 = OMD * FAX
omch2 = OMD - omaxh2
NUIS = {n: None for n in ("alpha_1", "alpha_2", "gamma_1", "gamma_2")}

pars = camb.set_params(omch2=omch2, omaxh2=omaxh2, m_ax=M_AX,
                       redshifts=[2.0, 0.0], kmax=20.0, WantTransfer=True,
                       **fid)
pars.NonLinear = cm.NonLinear_none
res = camb.get_results(pars)
td = np.asarray(res.get_matter_transfer_data().transfer_data,
                dtype=np.float64)
k = td[0, :, 0]
zs = np.array(res.Params.Transfer.PK_redshifts[
  :res.Params.Transfer.PK_num_redshifts])


def boost(z_t, m_pts=100, m_lo=7, m_hi=18, k_slice=None):
  iz = int(np.argmin(np.abs(zs - z_t)))
  sl = slice(None) if k_slice is None else k_slice
  return k[sl], _compute_row({
    "z": z_t, "h": fid["H0"] / 100, "omega_b": fid["ombh2"],
    "omega_d": omch2, "omega_ax": omaxh2, "m_ax": M_AX,
    "lcdm_mode": False, "ns": fid["ns"], "As": fid["As"], "k_piv": 0.05,
    "version": "dome", "flags": _VERSION_FLAGS["dome"],
    "M_arr": np.logspace(m_lo, m_hi, m_pts), "nuisance": NUIS,
    "m_min_exponent": m_lo, "m_max_exponent": m_hi,
    "k": k[sl], "T_cdm": td[1, sl, iz], "T_b": td[2, sl, iz],
    "T_ax": td[13, sl, iz], "T_tot": td[6, sl, iz]})


for z_t in [0.0, 2.0]:
  k_base, b_base = boost(z_t)
  B_base = b_base**2
  mask = k_base <= 10.0
  print(f"z = {z_t:.0f}:")
  for label, kwargs in [
      ("M grid 100 -> 320   ", dict(m_pts=320)),
      ("M range 7..18 -> 6..19 (320)", dict(m_pts=320, m_lo=6, m_hi=19)),
      ("k grid every 2nd point", dict(k_slice=slice(None, None, 2)))]:
    k_v, b_v = boost(z_t, **kwargs)
    B_v = np.interp(k_base, k_v, b_v**2)
    dev = np.abs(B_v / B_base - 1)[mask]
    print(f"  {label}: max |dB/B| (k<=10) = {dev.max():.3e}, "
          f"median = {np.median(dev):.3e}")
