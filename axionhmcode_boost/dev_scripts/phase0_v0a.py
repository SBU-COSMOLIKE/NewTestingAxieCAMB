"""Phase 0 / V0a: trivial external-ratio null test against New_AxiECAMB via cobaya.

Part A: ratio = 2 everywhere -> Pk_grid(nonlinear) must equal 4 x Pk_grid(linear).
        (Mirrors cobaya's test_trivial_non_linear_ratio, pointed at the local port.)
Part B: lensed-Cl run with NonLinear_both and ratio = 1; from inside
        get_non_linear_ratio, log results.transfer_redshifts vs
        Params.Transfer.PK_redshifts (R4 confirmation) and confirm
        get_matter_transfer_data() works on the transfers-only results (R1).
"""
import numpy as np
from cobaya.theory import Theory
from cobaya.likelihood import Likelihood
from cobaya.model import get_model

CAMB_PATH = "/Users/vivianmiranda/data/research/WayneHu/rayne/New_AxiECAMB"
RATIO_AMP = 2.0

base_params = {
  "ombh2": 0.022274, "omch2": 0.11913, "cosmomc_theta": 0.01040867,
  "tau": 0.0639, "ns": 0.9667, "As": 2.105e-9,
}


class TrivialNonLinearRatio(Theory):
  def get_requirements(self):
    return "CAMB_transfers"

  def get_non_linear_ratio(self, results):
    k_h = np.logspace(-4, 2, 200)
    z = np.array(results.Params.Transfer.PK_redshifts[
      :results.Params.Transfer.PK_num_redshifts])
    z = np.sort(z)
    ratio = RATIO_AMP * np.ones((len(z), len(k_h)))
    return {"k_h": k_h, "z": z, "ratio": ratio}


class NonLinearRatioLike(Likelihood):
  def get_requirements(self):
    return {"Pk_grid": {"z": [0, 0.5, 1.0], "k_max": 10,
                        "nonlinear": [False, True]}}

  def logp(self, **params_values):
    k_lin, z_lin, pk_lin = self.provider.get_Pk_grid(nonlinear=False)
    k_nl, z_nl, pk_nl = self.provider.get_Pk_grid(nonlinear=True)
    np.testing.assert_allclose(k_nl, k_lin)
    np.testing.assert_allclose(z_nl, z_lin)
    dev = np.max(np.abs(pk_nl / (RATIO_AMP**2 * pk_lin) - 1))
    np.testing.assert_allclose(pk_nl, RATIO_AMP**2 * pk_lin, rtol=1e-4, atol=0)
    print(f"PART A OK: pk_nonlinear == {RATIO_AMP**2:.0f} x pk_linear "
          f"(max rel dev {dev:.2e})")
    return 0


grid_log = {}


class GridLoggerRatio(Theory):
  def get_requirements(self):
    return "CAMB_transfers"

  def get_non_linear_ratio(self, results):
    tz = np.sort(np.array(results.transfer_redshifts))
    npk = results.Params.Transfer.PK_num_redshifts
    pkz = np.sort(np.array(results.Params.Transfer.PK_redshifts[:npk]))
    td = results.get_matter_transfer_data()
    grid_log["transfer_redshifts"] = tz
    grid_log["PK_redshifts"] = pkz
    grid_log["transfer_data_shape"] = np.asarray(td.transfer_data).shape
    k_h = np.logspace(-4, 2, 100)
    return {"k_h": k_h, "z": tz, "ratio": np.ones((len(tz), len(k_h)))}


class ClLike(Likelihood):
  def get_requirements(self):
    return {"Cl": {"tt": 2000, "ee": 2000, "pp": 2000}}

  def logp(self, **params_values):
    cl = self.provider.get_Cl(ell_factor=True)
    print(f"PART B: got lensed Cls, TT[l=1000] = {cl['tt'][1000]:.4f} muK^2")
    return 0


print("=" * 70)
print("PART A: trivial ratio = 2, Pk_grid check")
print("=" * 70)
info_a = {
  "likelihood": {"like": NonLinearRatioLike},
  "theory": {
    "camb": {"path": CAMB_PATH, "use_non_linear_ratio": True},
    "my_nonlin": TrivialNonLinearRatio,
  },
  "params": base_params,
  "stop_at_error": True,
}
model_a = get_model(info_a)
model_a.loglikes({})

print("=" * 70)
print("PART B: lensed Cls with NonLinear_both, ratio = 1, grid logging")
print("=" * 70)
info_b = {
  "likelihood": {"like": ClLike},
  "theory": {
    "camb": {"path": CAMB_PATH, "use_non_linear_ratio": True,
             "extra_args": {"nonlinear": "NonLinear_both",
                            "lens_potential_accuracy": 1}},
    "my_nonlin": GridLoggerRatio,
  },
  "params": base_params,
  "stop_at_error": True,
}
model_b = get_model(info_b)
model_b.loglikes({})

tz = grid_log["transfer_redshifts"]
pkz = grid_log["PK_redshifts"]
print(f"transfer_redshifts: n = {len(tz)}, range [{tz.min():.3f}, {tz.max():.3f}]")
print(f"  first 6: {np.round(tz[:6], 3)}  last 3: {np.round(tz[-3:], 3)}")
print(f"PK_redshifts (Params.Transfer): n = {len(pkz)}, values = {np.round(pkz, 3)}")
print(f"matter transfer_data shape (nvar, nk, nz): {grid_log['transfer_data_shape']}")
print("R4 CONFIRMED" if len(tz) > len(pkz) else "R4 NOT CONFIRMED — investigate")
