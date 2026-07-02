"""Tests for AxionHMcodeBoost: run with pytest, or directly as a script.

Fast configuration: Pk_grid-only likelihood at z = [0, 2] (2-redshift boost
grid, a few seconds per model). The full lensed-Cl path is exercised by the
validation battery (dev_scripts/, strategy files 11-12).
"""
import os
import sys

import numpy as np
import pytest
from cobaya.likelihood import Likelihood
from cobaya.log import LoggedError
from cobaya.model import get_model

_HERE = os.path.dirname(os.path.abspath(__file__))
CAMB_PATH = os.path.dirname(os.path.dirname(_HERE))  # the New_AxiECAMB repo
AXHM_PATH = os.environ.get(
  "AXIONHMCODE_PATH",
  os.path.join(os.path.dirname(CAMB_PATH), "axionHMcode"))
sys.path.insert(0, os.path.dirname(_HERE))
from axionhmcode_boost import AxionHMcodeBoost  # noqa: E402

base_params = {
  "ombh2": 0.02237, "omch2": 0.108, "H0": 67.4, "tau": 0.0544,
  "ns": 0.9655, "As": 2.1e-9, "mnu": 0.06, "nnu": 3.046,
}
captured = {}


class PkRatioLike(Likelihood):
  def get_requirements(self):
    return {"Pk_grid": {"z": [0.0, 2.0], "k_max": 5.0,
                        "nonlinear": [False, True]}}

  def logp(self, **params_values):
    k, z, pk_lin = self.provider.get_Pk_grid(nonlinear=False)
    _, _, pk_nl = self.provider.get_Pk_grid(nonlinear=True)
    captured.update(k=k, z=z, ratio=pk_nl / pk_lin)
    return 0


def make_info(theory_opts=None, params_extra=None):
  opts = {"external": AxionHMcodeBoost, "axionhmcode_path": AXHM_PATH,
          "version": "dome"}
  opts.update(theory_opts or {})
  params = dict(base_params)
  params.update(params_extra or {})
  return {
    "likelihood": {"pkratio": PkRatioLike},
    "theory": {"camb": {"path": CAMB_PATH, "use_non_linear_ratio": True,
                        "extra_args": {"num_massive_neutrinos": 1}},
               "axboost": opts},
    "params": params, "stop_at_error": True,
  }


def run(info):
  captured.clear()
  model = get_model(info)
  model.loglikes({})
  return captured


def test_axion_boost_runs():
  # note cobaya's Pk_grid k units are 1/Mpc, so k = 1 here is ~1.5 h/Mpc,
  # where the z=0 dome boost is ~10 (mead2020 gives ~9 at that scale)
  out = run(make_info(params_extra={"m_ax": 1e-25, "omaxh2": 0.012},
                      theory_opts={"version": "dome"}))
  k, z, ratio = out["k"], out["z"], out["ratio"]
  iz0, iz2 = np.argmin(np.abs(z - 0)), np.argmin(np.abs(z - 2))
  ilo, i1 = np.argmin(np.abs(k - 1e-3)), np.argmin(np.abs(k - 1.0))
  assert abs(ratio[iz0, ilo] - 1) < 1e-2, "boost must -> 1 at low k"
  assert 4 < ratio[iz0, i1] < 30, "z=0 boost at k=1/Mpc out of sane range"
  assert ratio[iz2, i1] < ratio[iz0, i1], "boost must grow toward z=0"
  assert np.all(np.isfinite(ratio))


def test_lcdm_limit_runs():
  out = run(make_info())  # no axion parameters at all
  k, z, ratio = out["k"], out["z"], out["ratio"]
  iz0 = np.argmin(np.abs(z - 0))
  ilo, i1 = np.argmin(np.abs(k - 1e-3)), np.argmin(np.abs(k - 1.0))
  assert abs(ratio[iz0, ilo] - 1) < 1e-2
  # cold-only halo model vs mead2020 at k = 1/Mpc: ~10.7 vs ~8.9 (+12%
  # model spread, consistent with the papers' LCDM agreement claims)
  assert 4 < ratio[iz0, i1] < 25, "LCDM-limit boost at k=1/Mpc out of range"


def test_de_like_hard_errors():
  info = make_info(params_extra={"m_ax": 1e-33, "omaxh2": 0.012})
  with pytest.raises(LoggedError, match="DE-like"):
    run(info)


def test_strict_gates_large_fax():
  extra = {"m_ax": 1e-24, "omaxh2": 0.06}
  info = make_info(params_extra=extra, theory_opts={"strict": True})
  info["params"]["omch2"] = 0.06
  with pytest.raises(LoggedError, match="strict"):
    run(info)
  info = make_info(params_extra=extra, theory_opts={"strict": False})
  info["params"]["omch2"] = 0.06
  out = run(info)  # warns, but must complete
  assert np.all(np.isfinite(out["ratio"]))


if __name__ == "__main__":
  for fn in [test_axion_boost_runs, test_lcdm_limit_runs,
             test_de_like_hard_errors, test_strict_gates_large_fax]:
    fn()
    print(f"{fn.__name__}: OK")
