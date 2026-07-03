"""Parent: build transfers for the validation cases, run upstream and fork
children, compare boost grids and timings.
"""
import subprocess
import sys
import numpy as np

SCRATCH = ("/private/tmp/claude-501/-Users-vivianmiranda-data-research-WayneHu"
           "-rayne/4e32f7cb-e470-4d65-8ad6-bf48eb7553b7/scratchpad")
CAMB_PATH = "/Users/vivianmiranda/data/research/WayneHu/rayne/New_AxiECAMB"
UPSTREAM = "/Users/vivianmiranda/data/research/WayneHu/rayne/axionHMcode"
FORK = "/Users/vivianmiranda/data/research/WayneHu/rayne/fork_axionHMcode"
PY = ("/Users/vivianmiranda/data/research/WayneHu/rayne/cobaya_test_env/bin/"
      "python")
sys.path.insert(0, CAMB_PATH)

import camb
from camb import model as cm

fid = dict(H0=67.32, ombh2=0.022383, As=2.101e-9, ns=0.96605,
           num_massive_neutrinos=1, mnu=0.06, nnu=3.046)
h = fid["H0"] / 100

cases = {
  # name: (omch2, omaxh2, m_ax, lcdm_mode)
  "gaughan": (0.12011 * 0.7, 0.12011 * 0.3, 1e-24, 0),
  "inputfile": (0.108, 0.012, 1e-25, 0),
  "lcdm": (0.12011, 0.0, 1e-25, 1),
}

blobs = {}
for name, (omch2, omaxh2, m_ax, lcdm) in cases.items():
  kwargs = dict(omch2=omch2, redshifts=[2.0, 0.0], kmax=10.0,
                WantTransfer=True, **fid)
  if not lcdm:
    kwargs.update(omaxh2=omaxh2, m_ax=m_ax)
  pars = camb.set_params(**kwargs)
  pars.NonLinear = cm.NonLinear_none
  res = camb.get_results(pars)
  td = np.asarray(res.get_matter_transfer_data().transfer_data,
                  dtype=np.float64)
  zs = np.array(res.Params.Transfer.PK_redshifts[
    :res.Params.Transfer.PK_num_redshifts])
  for z_t in (0.0, 2.0):
    iz = int(np.argmin(np.abs(zs - z_t)))
    base = f"{name}_z{z_t:.0f}"
    blobs[f"case_{base}"] = np.array(
      [z_t, fid["ombh2"], omch2, omaxh2 if not lcdm else 0.0, m_ax, h,
       fid["ns"], fid["As"], lcdm])
    blobs[f"k_{base}"] = td[0, :, iz]
    blobs[f"Tc_{base}"] = td[1, :, iz]
    blobs[f"Tb_{base}"] = td[2, :, iz]
    blobs[f"Ta_{base}"] = td[13, :, iz]
    blobs[f"Tt_{base}"] = td[6, :, iz]

np.savez(f"{SCRATCH}/fork_val_transfers.npz", **blobs)
print("transfers saved; running children...")

for tag, checkout in [("upstream", UPSTREAM), ("fork", FORK)]:
  proc = subprocess.run(
    [PY, f"{SCRATCH}/fork_validate_child.py", checkout,
     f"{SCRATCH}/fork_val_transfers.npz", f"{SCRATCH}/fork_val_{tag}.npz"],
    capture_output=True, text=True, timeout=1800)
  if proc.returncode != 0:
    print(f"{tag} FAILED:\n", proc.stderr[-2500:])
    sys.exit(1)
  tail = [ln for ln in proc.stdout.split("\n")
          if ln.strip() and "mismatch" not in ln and "Warning" not in ln]
  for ln in tail:
    print(f"[{tag}] {ln}")

up = np.load(f"{SCRATCH}/fork_val_upstream.npz")
fk = np.load(f"{SCRATCH}/fork_val_fork.npz")
print()
print(f"{'case':28s} {'max|dB/B|':>11s} {'t_up':>7s} {'t_fork':>8s} "
      f"{'speedup':>8s}")
worst = 0.0
for key in sorted(up.files):
  if key.startswith("time_"):
    continue
  B_up, B_fk = up[key]**2, fk[key]**2
  dev = float(np.max(np.abs(B_fk / B_up - 1)))
  worst = max(worst, dev)
  t_u, t_f = float(up[f"time_{key}"]), float(fk[f"time_{key}"])
  print(f"{key:28s} {dev:11.2e} {t_u:7.2f} {t_f:8.3f} {t_u/t_f:7.1f}x")
print(f"\nWORST max|dB/B| across all cases: {worst:.2e}")
