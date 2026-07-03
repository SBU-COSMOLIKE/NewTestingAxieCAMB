"""Cobaya Theory class providing the axionHMcode nonlinear boost to AxiECAMB.

Provides `non_linear_ratio` = sqrt(P_NL/P_L) to the stock cobaya CAMB wrapper
(theory block option `use_non_linear_ratio: True`, cobaya >= 3.6.2). The CAMB
wrapper calls get_non_linear_ratio(results) from inside its own calculate(),
handing over the transfers-level CAMBdata; this class extracts the linear
transfer functions, runs axionHMcode (Vogt et al. 2209.13445; Dome et al.
2409.11469) at every redshift CAMB uses, and returns the boost grid.

Boost convention: numerator and denominator both live in axionHMcode's Eq. 9
decomposition. The denominator is the model's own linear limit, which is a
perfect square, so the returned ratio -> 1 at low k by construction and CAMB's
linear spectrum is untouched there:

  sqrt(P_L_eq9) = (O_db/O_m) sqrt(P_cold)
                  + (O_ax/O_m) [fc sqrt(P_cold) + (1 - fc) sqrt(P_ax)]

Validity policy (`strict` option): DE-like axions (m/H0 < 10) and a z grid
reaching the KG->EFA switch always hard-error (the mixed-dark-matter halo
model is undefined there). Axion fractions or masses outside the version's
calibration warn-and-extrapolate when strict is False (the Gaughan et al.
2605.12054 practice) and hard-error when strict is True. The lensing z grid
inevitably spans redshifts outside dome's 1 < z < 8 calibration; that is a
documented limitation, logged once, not gated by strict.

Requires the AxiECAMB port (Transfer_axion variable) as the camb module and an
unmodified axionHMcode checkout (called only through public entry points, so
axionHMcode updates drop in).
"""

import os
import sys
import types

import numpy as np
from cobaya.log import LoggedError
from cobaya.theory import Theory

# matter transfer indices of the AxiECAMB port (camb/model.py); the shape of
# transfer_data is guarded at runtime before these are used
_T_KH, _T_CDM, _T_B, _T_TOT, _T_AXION = 1, 2, 3, 7, 14

# per-version axionHMcode call flags, following the example notebook and the
# upstream README (dome: alpha + concentration calibrated, full_2h=False used
# in the calibration; basic: plain sum, full two-halo term)
_VERSION_FLAGS = {
  "dome": dict(alpha=True, eta_given=False, one_halo_damping=True,
               two_halo_damping=False, concentration_param=True,
               full_2h=False),
  "basic": dict(alpha=False, eta_given=False, one_halo_damping=True,
                two_halo_damping=False, concentration_param=False,
                full_2h=True),
}

_NUISANCE = ("alpha_1", "alpha_2", "gamma_1", "gamma_2")

# dome calibration domain (Dome et al. 2409.11469): fax = O_ax/O_m in
# [0.01, 0.3], m near 1e-24.5 eV; basic (Vogt et al.): fax_dm < 0.5
_DOME_FAX_M_MAX = 0.3
_DOME_LOG10M_RANGE = (-25.5, -23.5)
_BASIC_FAX_DM_MAX = 0.5

# z grids reaching this fraction of z_osc are refused (R10): above the
# KG->EFA switch the axion transfer is not a halo-model density contrast
_ZOSC_MARGIN = 0.2

# omega_ax h^2 below this is treated as LCDM: axionHMcode's own LCDM recipe
# (tiny axion fraction) is used so the boost reduces to the cold halo model
_OMAXH2_LCDM = 1e-8


def _compute_row(payload):
  """One redshift of the boost grid; module-level so fork workers can run it.

  All axionHMcode calls go through its public entry points only.
  """
  from axionCAMB_and_lin_PS import lin_power_spectrum
  from cosmology.overdensities import func_D_z_unnorm_int
  from halo_model import HMcode_params
  from axion_functions import axion_params
  from halo_model import PS_nonlin_axion
  from halo_model import PS_nonlin_cold

  # accuracy_boost hook: only the SBU-COSMOLIKE fork has fast_tables;
  # unmodified upstream has no tables to scale, so the hook is a no-op there
  try:
    from cosmology import fast_tables
  except ImportError:
    pass
  else:
    fast_tables.set_accuracy_boost(payload.get("accuracy_boost", 1.0))

  z = payload["z"]
  h = payload["h"]
  cd = {
    "M_min": payload["m_min_exponent"], "M_max": payload["m_max_exponent"],
    "transfer_kmax": float(payload["k"].max()) * h,
    "version": payload["version"],
    "omega_b_0": payload["omega_b"], "omega_d_0": payload["omega_d"],
    "omega_ax_0": payload["omega_ax"],
    "omega_db_0": payload["omega_d"] + payload["omega_b"],
    "omega_m_0": payload["omega_b"] + payload["omega_d"]
                 + payload["omega_ax"],
    "m_ax": payload["m_ax"], "h": h, "z": z,
    "ns": payload["ns"], "As": payload["As"], "k_piv": payload["k_piv"],
  }
  for key in ("b", "d", "ax", "db", "m"):
    cd[f"Omega_{key}_0"] = cd[f"omega_{key}_0"] / h**2
  cd["Omega_w_0"] = 1 - cd["Omega_m_0"]
  cd["G_a"] = func_D_z_unnorm_int(z, cd["Omega_m_0"], cd["Omega_w_0"])
  for name, value in payload["nuisance"].items():
    if value is not None:
      cd[name] = value

  k = payload["k"]
  t2p = lambda T: lin_power_spectrum.transfer_to_PS(k, T, cd)
  cold_T = (cd["Omega_b_0"] * payload["T_b"]
            + cd["Omega_d_0"] * payload["T_cdm"]) / cd["Omega_db_0"]
  pd = {"k": k, "power_total": t2p(payload["T_tot"]),
        "power_CDM": t2p(payload["T_cdm"]),
        "power_baryon": t2p(payload["T_b"]),
        "power_cold": t2p(cold_T),
        "power_axion": t2p(payload["T_ax"])}

  flags = payload["flags"]
  M_arr = payload["M_arr"]
  hmc = HMcode_params.HMCode_param_dic(cd, k, pd["power_cold"])

  if payload["lcdm_mode"]:
    # axionHMcode's own LCDM recipe: cold-only halo model (the full axion
    # assembly is singular at vanishing axion fraction)
    out = PS_nonlin_cold.func_non_lin_PS_matter(
      M_arr, k, pd["power_cold"], cd, hmc, cd["Omega_db_0"],
      alpha=flags["alpha"], eta_given=flags["eta_given"], ax_one_halo=False,
      one_halo_damping=flags["one_halo_damping"],
      two_halo_damping=flags["two_halo_damping"],
      concentration_param=flags["concentration_param"],
      full_2h=flags["full_2h"], axion_dic=None)
    return np.sqrt(np.asarray(out[0]) / pd["power_cold"])

  axd = axion_params.func_axion_param_dic(
    M_arr, cd, pd, hmc, concentration_param=flags["concentration_param"])
  out = PS_nonlin_axion.func_full_halo_model_ax(
    M_arr, pd, cd, hmc, axd, **flags)

  fc = float(axd["frac_cluster"])
  wdb = cd["Omega_db_0"] / cd["Omega_m_0"]
  wax = cd["Omega_ax_0"] / cd["Omega_m_0"]
  sqc = np.sqrt(pd["power_cold"])
  sqa = np.sqrt(pd["power_axion"])
  sqrt_pl_eq9 = wdb * sqc + wax * (fc * sqc + (1 - fc) * sqa)
  return np.sqrt(np.asarray(out[0])) / sqrt_pl_eq9


class AxionHMcodeBoost(Theory):
  """axionHMcode boost provider for `camb: {use_non_linear_ratio: True}`."""

  # path to the unmodified axionHMcode checkout (required)
  axionhmcode_path: str = ""
  # 'dome' (default; most recent recalibration) or 'basic'
  version: str = "dome"
  # True: hard-error outside the calibration domain; False: warn + extrapolate
  strict: bool = False
  # expose alpha_1/alpha_2/gamma_1/gamma_2 as sampled input parameters
  sample_nuisance: bool = False
  # fixed nuisance values used when sample_nuisance is False (None = omit,
  # axionHMcode then uses its internal calibrated defaults)
  alpha_1: float | None = None
  alpha_2: float | None = None
  gamma_1: float | None = None
  gamma_2: float | None = None
  # halo-mass integration grid: logspaced 10**m_min .. 10**m_max Msun/h
  m_grid_points: int = 100
  m_min_exponent: float = 7.0
  m_max_exponent: float = 18.0
  # CAMB-style single accuracy multiplier: scales the halo-mass grid point
  # count and, when the checkout is the SBU-COSMOLIKE fork, its internal
  # growth/G/sigma(M) table node counts. Rerun the same yaml with 1 and 2
  # and compare the boost grids to check convergence; 1 = the validated
  # defaults, so results are unchanged unless the flag is set.
  accuracy_boost: float = 1.0
  # override individual axionHMcode call flags (expert use; V6 territory)
  model_flags: dict | None = None
  # Parallelism note (not a yaml option by design): the redshift loop forks
  # one worker per OMP_NUM_THREADS core (identical numerics, wall time only;
  # falls back to 1 when unset). The environment's core budget is always
  # respected — the same contract as every OpenMP code in Cocoa, and correct
  # per MPI rank (each rank forks within its own allocation). There is
  # deliberately no yaml override (PI decision, 2026-07-03).

  def initialize(self):
    if not self.axionhmcode_path:
      raise LoggedError(
        self.log, "Set axionhmcode_path to the axionHMcode checkout "
        "(e.g. ./external_modules/code/axionHMcode).")
    path = os.path.abspath(os.path.expanduser(self.axionhmcode_path))
    if not os.path.isfile(os.path.join(path, "halo_model",
                                       "PS_nonlin_axion.py")):
      raise LoggedError(
        self.log, "axionhmcode_path=%s does not look like an axionHMcode "
        "checkout (halo_model/PS_nonlin_axion.py not found).", path)
    # axionHMcode carries a dead `from scipy import misc` import; scipy
    # removed scipy.misc in 1.14 -- pre-seed a stub so the import survives
    # without touching upstream files
    try:
      import scipy.misc  # noqa: F401
    except Exception:
      import scipy
      stub = types.ModuleType("scipy.misc")
      sys.modules["scipy.misc"] = stub
      scipy.misc = stub
    if path not in sys.path:
      sys.path.insert(0, path)
    # import once here so failures surface at initialization, and so fork
    # workers inherit the compiled numba state
    from halo_model import PS_nonlin_axion  # noqa: F401
    from axion_functions import axion_params  # noqa: F401

    self.version = str(self.version).lower()
    if self.version not in _VERSION_FLAGS:
      raise LoggedError(self.log, "version must be one of %s, got %r",
                        sorted(_VERSION_FLAGS), self.version)
    self._flags = dict(_VERSION_FLAGS[self.version])
    self._flags.update(self.model_flags or {})
    self.accuracy_boost = float(self.accuracy_boost)
    if not self.accuracy_boost > 0:
      raise LoggedError(self.log, "accuracy_boost must be > 0, got %r",
                        self.accuracy_boost)
    n_m = max(20, int(round(self.m_grid_points * self.accuracy_boost)))
    self._M_arr = np.logspace(self.m_min_exponent, self.m_max_exponent, n_m)
    if self.accuracy_boost != 1.0:
      self.log.info(
        "accuracy_boost=%g: halo-mass grid %d points (base %d); fork table "
        "node counts scaled by the same factor when the checkout is the "
        "SBU-COSMOLIKE fork.", self.accuracy_boost, n_m, self.m_grid_points)
    try:
      env_threads = int(os.environ.get("OMP_NUM_THREADS", "") or 1)
    except ValueError:
      env_threads = 1
    self._n_processes = max(1, env_threads)
    self.log.info(
      "z-loop fork parallelism: %d worker(s) from OMP_NUM_THREADS%s.",
      self._n_processes,
      "" if "OMP_NUM_THREADS" in os.environ else " (unset, defaulting to 1)")
    self._warned = set()
    if self.version == "dome":
      self.log.info(
        "dome calibration covers 1 < z < 8; the CAMB lensing grid also needs "
        "z outside that range, where the boost is an extrapolation of the "
        "calibrated fitting formulae (not gated by `strict`).")

  def get_requirements(self):
    reqs = {"CAMB_transfers": None}
    if self.sample_nuisance:
      reqs.update({name: None for name in _NUISANCE})
    return reqs

  def get_non_linear_ratio(self, results):
    p = results.Params
    ax = p.Axion
    h = p.H0 / 100.0
    omega_b, omega_c = p.ombh2, p.omch2

    if ax.active and ax.is_de_like:
      raise LoggedError(
        self.log, "DE-like axion (m/H0 = %.3g < 10): the mixed-dark-matter "
        "halo model is undefined here regardless of `strict`. Use the "
        "halofit path (halofit_version: original, use_non_linear_ratio "
        "off) for DE-like masses.", ax.m_ovH0)

    if ax.active and ax.use_axfrac:
      omega_ax = ax.omdah2 * ax.axfrac
    elif ax.active:
      omega_ax = ax.omaxh2
    else:
      omega_ax = 0.0
    lcdm_mode = omega_ax < _OMAXH2_LCDM
    if lcdm_mode:
      # LCDM limit: axionHMcode's own recipe (tiny fraction, dummy mass,
      # cold-only halo model in _compute_row)
      omega_ax, m_ax = 1e-20 * h**2, 1e-25
    else:
      m_ax = ax.m_ax
      self._check_validity(omega_ax, omega_c, omega_b, m_ax)

    ip = p.InitPower
    for attr in ("As", "ns", "pivot_scalar"):
      if not hasattr(ip, attr):
        raise LoggedError(
          self.log, "InitPower has no %r: only power-law primordial spectra "
          "are supported (axionHMcode's internal P_prim is a power law).",
          attr)

    z_raw = np.asarray(results.transfer_redshifts, dtype=np.float64)
    if z_raw.size == 0:
      raise LoggedError(self.log, "results has no transfer redshifts.")
    if ax.active and 0 < ax.a_osc < 1:
      z_osc = 1.0 / ax.a_osc - 1.0
      if z_raw.max() > _ZOSC_MARGIN * z_osc:
        raise LoggedError(
          self.log, "z grid reaches z_max = %.3g > %.2g * z_osc (z_osc = "
          "%.3g): above the KG->EFA switch the axion transfer function is "
          "not a halo-model density contrast (hard error regardless of "
          "`strict`).", z_raw.max(), _ZOSC_MARGIN, z_osc)

    td = np.asarray(results.get_matter_transfer_data().transfer_data,
                    dtype=np.float64)
    if td.shape[0] < _T_AXION:
      raise LoggedError(
        self.log, "transfer_data has %d variables; the AxiECAMB port "
        "(Transfer_axion = %d) is required.", td.shape[0], _T_AXION)
    if td.shape[2] != z_raw.size:
      raise LoggedError(
        self.log, "transfer_data nz = %d does not match transfer_redshifts "
        "(%d).", td.shape[2], z_raw.size)

    order = np.argsort(z_raw)
    k = td[_T_KH - 1, :, 0]
    nuisance = self._nuisance_values()
    payloads = []
    for iz in order:
      payloads.append({
        "z": float(z_raw[iz]), "h": h, "omega_b": omega_b,
        "omega_d": omega_c, "omega_ax": omega_ax, "m_ax": m_ax,
        "lcdm_mode": lcdm_mode,
        "ns": ip.ns, "As": ip.As, "k_piv": ip.pivot_scalar,
        "version": self.version, "flags": self._flags,
        "M_arr": self._M_arr, "nuisance": nuisance,
        "m_min_exponent": self.m_min_exponent,
        "m_max_exponent": self.m_max_exponent,
        "accuracy_boost": self.accuracy_boost,
        "k": k, "T_cdm": td[_T_CDM - 1, :, iz], "T_b": td[_T_B - 1, :, iz],
        "T_ax": td[_T_AXION - 1, :, iz], "T_tot": td[_T_TOT - 1, :, iz],
      })

    n_proc = min(self._n_processes, len(payloads))
    if n_proc > 1:
      import multiprocessing
      ctx = multiprocessing.get_context("fork")
      with ctx.Pool(n_proc) as pool:
        rows = pool.map(_compute_row, payloads)
    else:
      rows = [_compute_row(payload) for payload in payloads]

    ratio = np.vstack(rows)
    if not np.all(np.isfinite(ratio)):
      raise LoggedError(self.log, "non-finite boost values in the grid.")
    return {"k_h": k, "z": z_raw[order], "ratio": ratio}

  def _nuisance_values(self):
    if self.sample_nuisance:
      return {name: float(self.provider.get_param(name))
              for name in _NUISANCE}
    return {name: getattr(self, name) for name in _NUISANCE}

  def _check_validity(self, omega_ax, omega_c, omega_b, m_ax):
    if self.version == "dome":
      fax_m = omega_ax / (omega_ax + omega_c + omega_b)
      if fax_m > _DOME_FAX_M_MAX:
        self._out_of_domain(
          "dome_fax", f"fax = O_ax/O_m = {fax_m:.3f} > {_DOME_FAX_M_MAX} "
          "(dome calibration bound)")
      log10m = np.log10(m_ax)
      if not _DOME_LOG10M_RANGE[0] <= log10m <= _DOME_LOG10M_RANGE[1]:
        self._out_of_domain(
          "dome_mass", f"log10(m_ax/eV) = {log10m:.2f} outside "
          f"[{_DOME_LOG10M_RANGE[0]}, {_DOME_LOG10M_RANGE[1]}] (dome "
          "calibration is centred on -24.5)")
    else:
      fax_dm = omega_ax / (omega_ax + omega_c)
      if fax_dm > _BASIC_FAX_DM_MAX:
        self._out_of_domain(
          "basic_fax", f"fax_dm = O_ax/O_dm = {fax_dm:.3f} > "
          f"{_BASIC_FAX_DM_MAX} (basic halo-model bound)")

  def _out_of_domain(self, key, message):
    if self.strict:
      raise LoggedError(
        self.log, "outside axionHMcode calibration and strict=True: %s",
        message)
    if key not in self._warned:
      self._warned.add(key)
      self.log.warning(
        "outside axionHMcode calibration, extrapolating (strict=False): "
        "%s. Further occurrences of this condition are not repeated.",
        message)
