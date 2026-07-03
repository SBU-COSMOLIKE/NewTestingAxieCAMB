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

Reading guide (the order in which cobaya exercises this file):

1. `AxionHMcodeBoost.initialize` — runs once at startup: checks the
   axionHMcode checkout, imports it, resolves the yaml options, and decides
   the worker count from OMP_NUM_THREADS.
2. `AxionHMcodeBoost.get_requirements` — tells cobaya that this theory needs
   the CAMB transfer functions (and nothing else; see the developer guide
   for why requesting the power spectrum here would create a real
   dependency cycle).
3. `AxionHMcodeBoost.get_non_linear_ratio` — called by the CAMB wrapper on
   every likelihood evaluation with the freshly computed transfers. It
   validates the axion regime, builds one work package ("payload") per
   redshift, runs `_compute_row` on each (in parallel worker processes when
   OMP_NUM_THREADS > 1), and returns the boost grid.
4. `_compute_row` — the physics: one redshift of the boost. Module-level
   (outside the class) because forked worker processes must be able to call
   it without carrying the cobaya machinery.
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
  """Compute one redshift of the boost grid: sqrt(P_NL/P_L) on the k grid.

  Arguments:
    payload — a plain dictionary with everything this redshift needs
      (cosmological densities in omega = Omega h^2 form, the axion mass in
      eV, the primordial amplitude/tilt/pivot, the transfer-function columns
      on the k grid in h/Mpc, the halo-mass grid in Msun/h, and the yaml
      options). A plain dictionary — no cobaya objects — so that forked
      worker processes can receive it cheaply.

  Returns:
    A 1-D numpy array over k: sqrt(B) with B = P_NL / P_L_eq9, the square
    root because cobaya's non_linear_ratio convention is the ratio of
    amplitudes, not of powers.

  All axionHMcode calls go through its public entry points only (the
  drag-and-drop constraint: upstream updates must drop in unchanged).
  """
  from axionCAMB_and_lin_PS import lin_power_spectrum
  from cosmology.overdensities import func_D_z_unnorm_int
  from halo_model import HMcode_params
  from axion_functions import axion_params
  from halo_model import PS_nonlin_axion
  from halo_model import PS_nonlin_cold

  # fork hooks: only the SBU-COSMOLIKE fork has fast_tables; unmodified
  # upstream has neither tables to scale nor a solver-mode switch, so both
  # hooks are no-ops there (initialize() has already hard-errored if the
  # default solver was requested without the fork). Payloads without the
  # key — the standalone validation scripts — default to the legacy solver
  # so that fork_validate keeps certifying the bit-faithful path.
  try:
    from cosmology import fast_tables
  except ImportError:
    pass
  else:
    fast_tables.set_accuracy_boost(payload.get("accuracy_boost", 1.0))
    fast_tables.set_legacy_root_finder(
      payload.get("legacy_root_finder", True))

  # --- step 1: the cosmology dictionary, in axionHMcode's own conventions.
  # Lower-case omega_X_0 are physical densities (omega = Omega h^2);
  # upper-case Omega_X_0 are density fractions. Suffixes: b = baryons,
  # d = cold dark matter only, ax = axion, db = cdm + baryons ("cold"),
  # m = total matter, w = dark energy (flatness: Omega_w = 1 - Omega_m).
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
  cd["Omega_b_0"] = cd["omega_b_0"] / h**2
  cd["Omega_d_0"] = cd["omega_d_0"] / h**2
  cd["Omega_ax_0"] = cd["omega_ax_0"] / h**2
  cd["Omega_db_0"] = cd["omega_db_0"] / h**2
  cd["Omega_m_0"] = cd["omega_m_0"] / h**2
  cd["Omega_w_0"] = 1 - cd["Omega_m_0"]
  # G_a is the HMcode-2020 integrated growth (eq. A5), a single number per
  # (cosmology, redshift) that axionHMcode expects precomputed in cosmo_dic
  cd["G_a"] = func_D_z_unnorm_int(z, cd["Omega_m_0"], cd["Omega_w_0"])
  # optional Dentler et al. nuisance parameters: axionHMcode activates them
  # simply by their presence as keys in the cosmology dictionary
  for name, value in payload["nuisance"].items():
    if value is not None:
      cd[name] = value

  # --- step 2: linear power spectra from the transfer-function columns.
  k = payload["k"]

  def power_from_transfer(transfer):
    """Linear P(k) in (Mpc/h)^3 from one transfer-function column, using
    axionHMcode's own primordial power-law conventions (upstream function
    transfer_to_PS; k in h/Mpc, transfers in the CAMB output convention)."""
    return lin_power_spectrum.transfer_to_PS(k, transfer, cd)

  # "cold" = the density-weighted combination of cdm and baryons, the field
  # the halo model builds its halos from (Dome et al., footnote to eq. 37)
  cold_transfer = (cd["Omega_b_0"] * payload["T_b"]
                   + cd["Omega_d_0"] * payload["T_cdm"]) / cd["Omega_db_0"]
  pd = {"k": k,
        "power_total": power_from_transfer(payload["T_tot"]),
        "power_CDM": power_from_transfer(payload["T_cdm"]),
        "power_baryon": power_from_transfer(payload["T_b"]),
        "power_cold": power_from_transfer(cold_transfer),
        "power_axion": power_from_transfer(payload["T_ax"])}

  # --- step 3: the halo model. HMCode_param_dic derives the cold
  # HMcode-2020 parameters (sigma(M), formation redshifts, damping scales);
  # func_axion_param_dic derives the axion quantities (cut mass, soliton
  # central densities, clustered fraction); func_full_halo_model_ax
  # assembles the nonlinear total power spectrum.
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
    # out is a tuple; element [0] is the total nonlinear P(k)
    return np.sqrt(np.asarray(out[0]) / pd["power_cold"])

  axd = axion_params.func_axion_param_dic(
    M_arr, cd, pd, hmc, concentration_param=flags["concentration_param"])
  out = PS_nonlin_axion.func_full_halo_model_ax(
    M_arr, pd, cd, hmc, axd, **flags)

  # --- step 4: the boost denominator — the model's own linear limit
  # (the "Eq. 9 convention" of the developer guide). It is a perfect
  # square by construction, assembled from the same ingredients the halo
  # model uses, so B -> 1 at low k exactly:
  #   sqrt(P_L) = w_db sqrt(P_cold) + w_ax [fc sqrt(P_cold)
  #                                         + (1 - fc) sqrt(P_axion)]
  # with w_db, w_ax the density weights and fc the clustered axion fraction.
  fc = float(axd["frac_cluster"])
  weight_db = cd["Omega_db_0"] / cd["Omega_m_0"]
  weight_ax = cd["Omega_ax_0"] / cd["Omega_m_0"]
  sqrt_power_cold = np.sqrt(pd["power_cold"])
  sqrt_power_axion = np.sqrt(pd["power_axion"])
  sqrt_linear_eq9 = (weight_db * sqrt_power_cold
                     + weight_ax * (fc * sqrt_power_cold
                                    + (1 - fc) * sqrt_power_axion))
  return np.sqrt(np.asarray(out[0])) / sqrt_linear_eq9


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
  # Solver mode (PI decision 2026-07-03: the re-engineered solver is the
  # default). False (default): bracketed brentq with residual
  # classification, closest-achievable acceptance for unreachable targets,
  # per-evaluation diagnostics, and a continuously interpolated soliton/NFW
  # crossover cell — requires the SBU-COSMOLIKE fork; validated by
  # dome-version agreement with upstream (<= 7e-5 in B), attributed
  # basic-version differences, and posterior-point delta-chi2 (fork README
  # appendix A.12). True: the released code's solver verbatim —
  # optimize.root(hybr) with the |guess - rho_c| > 100 acceptance net —
  # bit-faithful to upstream (max |dB/B| = 1.6e-5, the fork_validate gate);
  # kept for upstream comparison, and the only mode a plain upstream
  # checkout supports. Keep the flag visible in every yaml and constant
  # within a chain.
  legacy_root_finder: bool = False
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
    self.legacy_root_finder = bool(self.legacy_root_finder)
    if not self.legacy_root_finder:
      try:
        from cosmology import fast_tables  # noqa: F401
      except ImportError:
        raise LoggedError(
          self.log, "the default solver requires the SBU-COSMOLIKE "
          "axionHMcode fork (cosmology/fast_tables.py not found in %s). "
          "Point axionhmcode_path at the fork, or set "
          "legacy_root_finder: true to run the plain upstream checkout "
          "with its verbatim solver.", path)
      self.log.info(
        "solver: default (bracketed brentq + residual classification, "
        "interpolated crossover cell); set legacy_root_finder: true for "
        "the released code's verbatim behavior.")
    else:
      self.log.info(
        "solver: legacy_root_finder (released axionHMcode verbatim; "
        "bit-faithful to upstream).")

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
    # becomes True after the first evaluation compiles numba in this parent
    # process (see get_non_linear_ratio); keeps the ~6 s JIT a one-time cost
    self._warmed = False
    if self.version == "dome":
      self.log.info(
        "dome calibration covers 1 < z < 8; the CAMB lensing grid also needs "
        "z outside that range, where the boost is an extrapolation of the "
        "calibrated fitting formulae (not gated by `strict`).")

  def get_requirements(self):
    """Declare what this theory needs from the rest of the cobaya graph.

    Only CAMB_transfers (the Boltzmann solve) — never the power spectrum,
    which would create a genuine dependency cycle — plus, when nuisance
    sampling is on, the four Dentler parameters as sampled inputs.
    """
    requirements = {"CAMB_transfers": None}
    if self.sample_nuisance:
      for name in _NUISANCE:
        requirements[name] = None
    return requirements

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
        "legacy_root_finder": self.legacy_root_finder,
        "k": k, "T_cdm": td[_T_CDM - 1, :, iz], "T_b": td[_T_B - 1, :, iz],
        "T_ax": td[_T_AXION - 1, :, iz], "T_tot": td[_T_TOT - 1, :, iz],
      })

    n_proc = min(self._n_processes, len(payloads))
    if n_proc > 1:
      import multiprocessing
      ctx = multiprocessing.get_context("fork")
      # numba compiles axionHMcode's just-in-time functions the first time
      # they run in a process, at a one-time cost of ~6 s. Worker processes
      # are created by forking this (parent) process, and fork copies the
      # parent's memory -- including numba's compiled machine code -- so if
      # the parent has already compiled, the workers inherit it for free.
      # On the first evaluation we therefore compute one redshift here in
      # the parent (which compiles numba) and keep its result, so the pool
      # only does the remaining redshifts: the first evaluation costs the
      # same as before, but every later evaluation forks an already-warm
      # parent and skips the recompile. Without this, each evaluation forks
      # fresh un-compiled workers and re-pays the ~6 s compile every step
      # (measured: ~10 s/step cold vs ~5 s/step warmed).
      if not self._warmed:
        first_row = _compute_row(payloads[0])
        self._warmed = True
        with ctx.Pool(n_proc) as pool:
          rest_rows = pool.map(_compute_row, payloads[1:])
        rows = [first_row] + rest_rows
      else:
        with ctx.Pool(n_proc) as pool:
          rows = pool.map(_compute_row, payloads)
    else:
      rows = [_compute_row(payload) for payload in payloads]

    ratio = np.vstack(rows)
    if not np.all(np.isfinite(ratio)):
      raise LoggedError(self.log, "non-finite boost values in the grid.")
    return {"k_h": k, "z": z_raw[order], "ratio": ratio}

  def _nuisance_values(self):
    """Current values of the four Dentler nuisance parameters, as a
    name -> value dictionary (value None means "omit", in which case
    axionHMcode falls back to its internal calibrated defaults).

    Sampled mode reads them from the cobaya provider (they change every
    likelihood evaluation); fixed mode reads the yaml options.
    """
    values = {}
    if self.sample_nuisance:
      for name in _NUISANCE:
        values[name] = float(self.provider.get_param(name))
    else:
      for name in _NUISANCE:
        values[name] = getattr(self, name)
    return values

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
