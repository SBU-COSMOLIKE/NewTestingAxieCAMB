"""Phase 3 smoke test: EXAMPLE_AXIONHMCODE_EVALUATE1.yaml wiring, with the
data likelihoods replaced by a lightweight Cl-requiring stand-in (the real
likelihood run needs cobaya-install'ed data; that is validation V5).
Run from the New_AxiECAMB repo root.
"""
import time
import numpy as np
from cobaya.yaml import yaml_load_file
from cobaya.model import get_model

info = yaml_load_file("EXAMPLE_AXIONHMCODE_EVALUATE1.yaml")

store = {}


def cl_probe(_self=None):
  cl = _self.provider.get_Cl(ell_factor=True)
  store["tt"] = cl["tt"]
  store["pp"] = cl["pp"]
  return 0.0


info["likelihood"] = {
  "cl_probe": {"external": cl_probe,
               "requires": {"Cl": {"tt": 2500, "ee": 2500, "pp": 2000}}}}
override = info["sampler"]["evaluate"]["override"]
info.pop("sampler")
info.pop("output")

t0 = time.time()
model = get_model(info)
print(f"get_model OK ({time.time()-t0:.1f} s) — yaml wiring, python_path "
      "class load, param chain all resolved")

t0 = time.time()
logpost = model.logposterior(override)
print(f"evaluation OK ({time.time()-t0:.1f} s): logpost = {logpost.logpost:.3f}")
print(f"lensed TT[l=1000] = {store['tt'][1000]:.2f} muK^2, "
      f"[L(L+1)]^2 Cpp/2pi [L=500] = {store['pp'][500]:.3e}")
derived = dict(zip(model.parameterization.derived_params(), logpost.derived))
print("derived:", {key: round(float(val), 4) for key, val in derived.items()})
