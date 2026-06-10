===========================
AxiECAMB (modern-CAMB port)
===========================

This is **AxiECAMB** — the ultralight-axion (ULA) effective method of
`arXiv:2412.15192 <https://arxiv.org/abs/2412.15192>`_ (Liu, Hu et al.) — ported from its
original CAMB-Nov13 base onto **modern CAMB 1.6.7**, including the Python wrapper.

When using the axion module, please cite `arXiv:2412.15192 <https://arxiv.org/abs/2412.15192>`_
(and Passaglia & Hu 2022, `arXiv:2201.10238 <https://arxiv.org/abs/2201.10238>`_, on which the
method builds). The original AxiECAMB heavily modified
`axionCAMB <https://github.com/dgrin1/axionCAMB>`_ (Hlozek et al., arXiv:1410.2896).

Method
======

The axion has a quadratic potential and is evolved in synchronous gauge:

- **Background**: the exact Klein-Gordon (KG) equation is solved (16-stage 8th-order
  fixed-step Runge-Kutta in ln a) from deep radiation domination until m = dfac*H, with
  the initial field value found by shooting to match the requested relic abundance.
  At the switch the field is projected onto WKB cos/sin amplitudes and matched onto an
  **effective fluid** (EFA) whose density follows a^-3 with an exp[3 int w dln a]
  residual, w(a) = wEFA_c (H/m)^2; the matching coefficients (<H>, wEFA_c) are iterated
  to self-consistency. dfac (default 10) is retuned internally: the oscillation phase
  at the switch is targeted to 2 beta = 7.08 pi for light DM-like axions, and the
  switch is pushed out of the recombination window z in (800, 1300).
- **Perturbations**: exact KG before the switch (with a per-k conditioning rescale of
  delta-phi), the (1+w)-weighted GDM effective fluid after, with sound speed
  cs^2 = (sqrt(1+kappa)-1)^2/kappa + (5/4)(H/am)^2, kappa = k^2/(a^2 m^2). At the
  switch the KG variables are projected onto (delta_ax, u_ax) preserving velocity and
  shear continuity; the residual metric jump is absorbed into eta (sub-horizon) or
  carried as a delta-function boundary term in the temperature line-of-sight integral
  (super-horizon).
- **DM vs DE**: for m/H0 >= 10 the axion is dark-matter-like — it counts in the matter
  transfer functions, sigma_8, the equality redshift, CosmoMC theta, and halofit
  Omega_m. For m/H0 < 10 (m <~ 1.4e-32 eV) it is dark-energy-like: KG is solved to
  a = 1, there is no fluid switch, and the matter transfer excludes the axion.

Usage (Fortran / .ini)
======================

Build and run (in ``fortran/``; build forutils first if needed)::

    make camb
    ./camb ../inifiles/params_axion.ini

New ini keys (see ``inifiles/params_axion.ini``):

==============================  ==================================================================
key                             meaning
==============================  ==================================================================
``m_ax``                        ULA mass in eV (negative input = log10(m_ax/eV))
``use_axfrac``                  T: use (``omdah2``, ``axfrac``); F: use ``omaxh2`` (+ usual ``omch2``)
``omaxh2``                      Omega_ax h^2 (when ``use_axfrac = F``)
``omdah2``                      total dark-matter Omega h^2 (when ``use_axfrac = T``)
``axfrac``                      axion fraction of DM (m/H0 >= 10) or of DE (m/H0 < 10)
``axion_dfac``                  switch threshold m = dfac*H (default 10; retuned internally)
``axion_isocurvature``, Hinf    accepted but isocurvature is force-disabled (v1.0 parity)
==============================  ==================================================================

With ``use_axfrac = T``, ``omch2`` may be omitted (it is derived). Constant-w dark
energy (fluid or PPF) can be combined with the axion; quintessence dark-energy models
cannot (the axion background solver treats DE as Lambda, as in the original).

Usage (Python)
==============

.. code-block:: python

    import camb
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=67.32, ombh2=0.02238, omch2=0.108, mnu=0.06, tau=0.054)
    pars.InitPower.set_params(As=2.1e-9, ns=0.966)
    pars.set_axion(m_ax=1e-27, omaxh2=0.012)                  # or omdah2=..., axfrac=...
    results = camb.get_results(pars)
    Ax = results.Params.Axion        # derived quantities live on the *result* state copy
    print(Ax.a_osc, Ax.dfac_used, Ax.tau_osc, Ax.m_ovH0)

The axion density perturbation is available as the ``delta_axion`` matter transfer
column (``camb.model.Transfer_axion``). ``delta_tot`` (and hence sigma_8 and the
default matter power) includes the axion when it is DM-like, excludes it when DE-like.

To build the Python library run ``make python`` in ``fortran/`` (or use
``python setup.py make``), which places ``camblib.so`` in ``camb/``.

Usage (Cobaya)
==============

This code works with **unmodified (pristine) Cobaya** — no patching of Cobaya's CAMB
theory wrapper is needed. The axion parameters (``m_ax``, ``omaxh2``, ``omdah2``,
``axfrac``, ``dfac``) are discovered and set through this package's own
``camb.get_valid_numerical_params`` and ``camb.set_params`` hooks, which the stock
wrapper already uses. Two ready-to-run examples ship in the repository root:

- ``EXAMPLE_EVALUATE1.yaml`` — single-point posterior evaluation;
- ``EXAMPLE_MCMC1.yaml`` — the corresponding MCMC (Planck lite + lowl TT/EE + DESI
  DR2 BAO + DES-Y5 SN + ACT DR6 lensing), sampling
  (logA, ns, 100theta_*, omegabh2, omegach2, tau, omegaaxh2, logmx).

Setup and run (from the repository root, with ``cobaya`` installed)::

    cd fortran && make camb && make python PYCAMB_OUTPUT_DIR=../camb/ && cd ..
    pip install act_dr6_lenslike
    cobaya-install EXAMPLE_MCMC1.yaml -p /path/to/cobaya_packages
    cobaya-run EXAMPLE_MCMC1.yaml -p /path/to/cobaya_packages

Conventions used in the examples (see comments inside the yaml files): the chain
samples ``thetastar100`` (= 100 theta_*) and feeds ``thetastar`` to CAMB via a
value-lambda (CAMB's input is theta_* itself); ``logA -> As``,
``omegaaxh2 -> omaxh2`` and ``logmx -> m_ax`` are standard value-lambda mappings;
``theta_H0_range: [40, 130]`` brackets the theta -> H0 solution over the whole prior
box; ``halofit_version: original`` is the AxiECAMB-validated non-linear treatment;
the ``omegam`` derived parameter excludes the axion (the DE-like convention).

What was ported and where
=========================

All modifications in the Fortran sources are marked with inline ``!AxiECAMB`` comments
(``grep -n AxiECAMB fortran/*.f90`` lists every change site).

- ``fortran/AxionBackground.f90`` (new): the KG background solver ``w_evolve``, EFA
  matching ``auxiIC``, phase targeting and recombination-skip dfac retuning (moved here
  from the old ``inidriver_axion.F90`` so the Python interface gets them too), as the
  component class ``TAxionModel`` stored in ``CAMBparams%Axion``.
- ``results.f90``: density budget/closure with the axion, the solver invocation,
  tau_osc, background integrals split at the dtauda kink at a_osc (applied uniformly:
  times, distances, sound horizons, optical depths — the original only split some),
  fine time-step window around tau_osc, thermo values cached at tau_osc,
  ``Transfer_axion`` column, z_eq and CosmoMC theta definitions.
- ``equations.f90``: the two axion perturbation equations (KG <-> EFA), the
  mid-evolution switch in the ``next_switch`` chain with the WKB projection
  (``AxionSwitchKGtoEFA``), adiabatic delta-phi initial conditions, axion terms in
  dgrho/dgq/grho/gpres, low-k lmaxnr boost ("WH smoother"). Tensors need no axion
  terms: the modern tensor background comes from the dtauda-based thermo table (this
  also fixes an original-code issue where tensors extrapolated the field table past
  a_osc).
- ``cmbmain.f90``: the switch boundary term in the temperature LOS integral
  (``deltaBCSrc`` machinery, flat and curved cases), axion-aware integration start
  time.
- ``recfast.f90``: dHdz in the tightly-coupled T_mat term includes the axion
  (numerical derivative of the exact H(z), stepped away from the a_osc kink).
- ``halofit.f90``: axion counted in Omega_m (DM-like) or in the smooth DE (DE-like),
  with a warning that the non-linear mode is inherited from axionCAMB and not well
  tested.
- Python: ``camb/axion.py`` (``AxionModel``), ``CAMBparams.set_axion``, transfer-name
  lists.

Validation against the original AxiECAMB
========================================

With matched cosmologies, the axion/LCDM suppression ratios agree between the original
AxiECAMB (Nov13 base) and this port to:

==========================================================  =========================  ====================================
case                                                        TT C_l ratio (l=2-2600)    P(k) ratio
==========================================================  =========================  ====================================
m=1e-27 eV, 10% of DM (switch z~1341)                       <= 0.01%                   <= 0.005%
m=1e-27 eV, 100% of DM (``use_axfrac``)                     <= 0.10%                   <= 0.005% (where suppression < 10^3)
m=1e-30 eV, 10% of DM (switch z~24, boundary term active)   <= 0.08%                   <= 0.01%
m=1.4e-33 eV (DE-like, no switch)                           <= 0.09%                   <= 0.01%
==========================================================  =========================  ====================================

Residuals at the 0.03-0.1% level are dominated by Nov13 <-> CAMB-1.6.7 baseline physics
differences that do not perfectly cancel in the ratios (the absolute LCDM baselines
differ by ~0.2%). The pure-LCDM limit of this code is bit-identical to unmodified
CAMB 1.6.7. The standard CAMB Python test suite passes.

Warnings / known differences (carried over or documented)
==========================================================

- **Isocurvature is disabled** (as in AxiECAMB v1.0; the original mode-6 vector
  targets variables that are not evolved). Inputs are accepted and ignored with a
  warning.
- The **growth-rate (Transfer_f) column** of the original is disabled there and was
  not ported; modern CAMB's own growth outputs are available.
- The **non-linear mode** is inherited from axionCAMB and not extensively tested;
  ``halofit_version = 1`` (original) is suggested — Takahashi was found unstable for
  axion models.
- For z > 0 transfer outputs, mind whether the requested z is before or after the
  switch: the axion density contrast is defined differently in the two regimes.
- P(k) at wavenumbers where the spectrum is suppressed by >~ 6 orders of magnitude
  differs from the original (which zeroed rather than extrapolated the dead tail);
  set ``transfer_kmax`` high enough for any application sensitive to that region.
- Default accuracy (``accuracy_boost = 1``) is what was validated in arXiv:2412.15192;
  higher boosts apply to the non-ULA accuracy settings only.
- CosmoMC theta counts the axion in omega_dm unconditionally (original behaviour),
  which is only meaningful for DM-like axions.

About the underlying CAMB
=========================

This code is built on CAMB 1.6.7 (Antony Lewis and Anthony Challinor,
https://camb.info/): a cosmology code for calculating cosmological observables,
including CMB, lensing, source count and 21cm angular power spectra, matter power
spectra, transfer functions and background evolution; Python package with numerical
code in modern Fortran. See the
`CAMB documentation <https://camb.readthedocs.io/en/latest/>`_ and the upstream
repository at https://github.com/cmbant/CAMB. You will need gfortran installed to
compile.


Port documentation (developer guide)
====================================

The sections below are the complete working documents produced while making this
port (June 2026). They are included verbatim so that anyone can learn how the port
was done, audit any individual change, or extend the code. Reading order:

1. *Port design* — the master change-by-change mapping from the original AxiECAMB
   to this code base, including the architecture decision and the list of things
   deliberately not ported.
2. *Modern CAMB architecture map* — where every kind of change lands in CAMB 1.6.7
   (parameter plumbing, background, perturbation hooks, python mirroring, build).
3. *Original-code analyses* — exhaustive inventories of every non-cosmetic change
   the original AxiECAMB made relative to its CAMB-Nov13 base, with verbatim
   equations and line references, one document per subsystem.

Conventions used throughout these documents:

- File/line references of the form ``AxiECAMB/equations_ppf.f90:3392`` or
  ``OLDCAMB/recfast.f90:828`` refer to the *original* AxiECAMB release and its
  pristine CAMB-Nov13 base, which live outside this repository (they were sibling
  directories during the port). References into ``fortran/model.f90`` etc. refer to
  pristine CAMB 1.6.7, i.e. (up to the port's own marked insertions) the files in
  this repository. Every modification made in this repository is marked with an
  inline ``!AxiECAMB`` comment, so ``grep -n AxiECAMB fortran/*.f90`` enumerates the
  realized version of everything described below.
- Classification tags used in the analyses:
  **[PHYSICS]** axion physics that had to be ported exactly;
  **[PLUMBING]** wiring re-derived for the class-based architecture;
  **[ACCURACY]** accuracy/sampling tweaks (ported when physics-motivated);
  **[OBSOLETE]** only meaningful in the Nov13 structure or superseded in 1.6.7;
  **[COSMETIC]** whitespace/comments (ignored).
- "RL"/"DG"/"RH"/"DM" in quoted comments are the original AxiECAMB/axionCAMB
  authors' initials; "Nov13" means the November 2013 CAMB release that
  axionCAMB/AxiECAMB were built on.


Port design: change-by-change mapping
-------------------------------------

Date: 2026-06-09. Sources: ``.port_analysis/reports/*.md`` (7 analysis reports), AxiECAMB originals.

0. Architecture decision
~~~~~~~~~~~~~~~~~~~~~~~~

The axion is implemented as a **new standalone component class** ``TAxionModel`` (module
``AxionBackground``, file ``fortran/AxionBackground.f90``), held in a new allocatable slot
``CAMBparams%Axion``, appended at the END of CAMBparams (and at the end of the python
``_fields_`` mirror). It is NOT placed in the DarkEnergy slot, because:

- the axion is dark *matter* for m/H0 ≥ 10 — Omega_de, halofit w(a), background density
  rows, and theta calculations stay honest;
- the user keeps full DE freedom (fluid / PPF with constant w + cs2_lam, like AxiECAMB's
  equations_ppf heritage). Quintessence/EarlyQuintessence DE + axion is rejected in
  Validate (circular background solves).

The component owns: all user inputs, the KG background tables, switch state, and the
EFA matching coefficients. Internal solver state lives in the object (thread/state safe;
no module globals). The driver-level dfac orchestration (phase targeting at 2β=7.08π +
recombination skip) moves INSIDE ``TAxionModel`` so the python path gets it too.

Activity: ``CP%Axion`` allocated always (default inactive); ``Axion%active`` true when
configured with m_ax>0 and nonzero abundance. All hooks guard on
``allocated(CP%Axion) .and. CP%Axion%active`` via helper ``State%HasAxion()``.

1. New file: fortran/AxionBackground.f90 (port of axion_background.F90)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Faithful port (identical numerics) of ``w_evolve``, ``derivs_bg``, ``auxiIC``, ``lh``,
``next_step``, ``get_phase_info`` + the driver loops, restructured as type-bound procedures:

.. code-block:: fortran

    type, extends(TCambComponent) :: TAxionModel
      ! --- inputs (python-mirrored) ---
      logical  :: active = .false.
      real(dl) :: m_ax = 0          ! eV (negative input interpreted as log10 by readers)
      real(dl) :: omaxh2 = 0        ! used if use_axfrac=F
      logical  :: use_axfrac = .false.
      real(dl) :: omdah2 = 0        ! total DM (m/H0>=10) ; with axfrac
      real(dl) :: axfrac = 0
      logical  :: axion_isocurvature = .false.  ! force-disabled (v1.0 parity)
      real(dl) :: Hinf = 13.7       ! log10 GeV input; stored as Hinf/Mpl after read
      real(dl) :: dfac = 10         ! switch threshold m = dfac*H (retuned internally)
      ! --- derived/state (python-mirrored, read-only for users) ---
      real(dl) :: m_ovH0, H0_eV, a_osc, tau_osc, aeq, aeq_LCDM, phiinit
      real(dl) :: ah_osc, ahosc_ETA, A_coeff, A_coeff_alt
      real(dl) :: tvarphi_c, tvarphi_s, tvarphi_cp, tvarphi_sp
      real(dl) :: wEFA_c, rhorefp_hsq, Prefp, wcorr_coeff
      real(dl) :: dfac_skip, a_skip, a_skipst
      real(dl) :: opac_tauosc, expmmu_tauosc      ! filled by Thermo_Init
      real(dl) :: amp_i, r_val, alpha_ax          ! isocurvature bookkeeping (inactive)
      real(dl) :: omegaax                          ! Omega_ax (final)
      logical  :: is_de_like                       ! m_ovH0 < 10
      logical  :: has_switch                       ! a_osc <= 1 and not de_like
      integer  :: ntable
      real(dl), allocatable :: loga_table(:)       ! log10(a) (kept; documented)
      real(dl), allocatable :: phinorm_table(:), phidotnorm_table(:)
      real(dl), allocatable :: phinorm_table_ddlga(:), phidotnorm_table_ddlga(:)
      real(dl), allocatable :: rhoaxh2ovrhom_logtable(:), rhoaxh2ovrhom_logtable_buff(:)
      real(dl) :: grhocrit_ovh2 = 0  ! grhocrit/h^2 cached for unit conversion
      real(dl) :: a_table_min = 0
    contains
      Solve (master: w_evolve + phase targeting + recomb-skip; called from SetParams)
      GrhoAx(a) -> 8 pi G rho_ax a^4  [Mpc^-2]  (table for a<=a_osc; analytic EFA beyond;
                  0 for a < a_table_min*?? -> follow AxiECAMB: clamp to table start value)
      FieldValsAta(a, v1, v2)  (spline lookup of phinorm/phidotnorm at log10 a)
      RhoaxH2AtA(a) -> Omega_ax(a) h^2  (raw units helper)
      ReadParams(Ini), Validate, SelfPointer, PythonClass
    end type

Solve() inputs come as plain arguments from results.f90 (no CAMBdata dependency, so the
module compiles before model.f90): omegah2_regm(b+c_eff), omegah2_rad parts, omk, hnot,
omegah2_lambda, nu data (eigenstates, masses from State%nu_masses, lhsq factors from
grhormass), TCMB. Massive-ν background via ``ThermalNuBackground%rho`` (module
MassiveNu, available pre-model). The Friedmann assembled internally MUST equal
CAMBdata's H(a): radiation = (grhog+grhornomass)/grhocrit·h², massive ν via
grhormass·rhonu, Λ = grhov/grhocrit·h². The DE is treated as Λ inside the solver
(AxiECAMB did the same even with PPF w; document).

Numerics kept verbatim: RK8 fixed step Butcher tableau, v2 attractor IC (1/5 factor),
ntable = nint(dfac*100)+1, shooting bisection tol 1e-6 (nphi=150 cap), auxiIC fixed
point tol 1e-7/30 iters with stateful wEFA_c (init 9/8 per Solve), natural-spline a_osc
root, adaptive a_final = 1.1·a_osc, log10 base change, analytic spline BCs, aeq +
aeq_LCDM, sentinels replaced by has_switch/is_de_like logicals at the API surface but
behavior identical. Phase targeting + recomb skip (a_skip=1/801, a_skipst=1/1301,
2β target 7.08π, tol 0.02π, dfac<23 gates) ported from inidriver_axion.F90:508-636.
Failure modes -> call GlobalError(error_unsupported_params) instead of stop.

2. model.f90 / model.py
~~~~~~~~~~~~~~~~~~~~~~~

- Append ``class(TAxionModel), allocatable :: Axion`` after CustomSources (last member).
- Validate: axion + (TEarlyQuintessence/TQuintessence DE) -> error; axion +
  Do21cm/CustomSources -> warning (untested); axion_isocurvature -> warn + disable
  (AxiECAMB v1.0 parity).
- python: new ``camb/axion.py`` with ``AxionModel(F2003Class)``,
  ``_fortran_class_module_="AxionBackground"``, ``_fortran_class_name_="TAxionModel"``;
  fields mirrored positionally (inputs + derived + AllocatableArrayDouble tables).
  model.py: import, append ``("Axion", AllocatableObject(AxionModel))`` at end of
  ``_fields_``; allocate default in CAMBparams.__init__ like other slots; convenience
  ``set_axion(...)`` on CAMBparams.

3. results.f90
~~~~~~~~~~~~~~

- SetParams (CAMBdata_SetParams):

  a) after h2 known and BEFORE Omega_de closure: if axion active resolve budget:

     - H0_eV, m_ovH0, is_de_like = m_ovH0 < 10
     - DM-like + use_axfrac: omaxh2 = axfrac*omdah2 ; CP%omch2 = (1-axfrac)*omdah2
     - DM-like + .not.use_axfrac: omaxh2 as given (omch2 as given)
     - DE-like + use_axfrac: CP%omch2 = omdah2 ;
       omegaax = axfrac*(1 - (ombh2+omch2+omnuh2)/h2 - omk - (grhog+grhornomass+Σgrhormass)/grhocrit)
     - DE-like + .not.use_axfrac: omegaax = omaxh2/h2
     - Omega_de = 1 - ... - omegaax  (closure includes axion)

  b) z_eq: add grhoax0=grhocrit*omegaax when DM-like (AxiECAMB modules.f90:2917 parity).
  c) after massive-nu init, BEFORE DarkEnergy%Init: call Axion%Solve(...) passing
     background pieces; on badflag -> GlobalError.
  d) after tau0: if has_switch: Axion%tau_osc = DeltaTime(0,a_osc) (split-aware).

- grho_no_de: ``+ this%grhoax(a)`` term via Axion%GrhoAx (a^4 units) — single chokepoint
  feeding dtauda everywhere (recfast/thermo/distances automatic).
- DeltaTime (CAMBdata_DeltaTime): if a_osc strictly inside (a1,a2): integrate
  [a1, a_osc*(1-tol_kink)] + [a_osc, a2] (kink in dtauda at the switch). This single
  generic fix replaces AxiECAMB's three scattered splits (DeltaTime, GYr, rs) and also
  covers reionization tau integrals that use DeltaTime; the z-space optical depth
  integral (reion_doptdepth_dz) gets the same split if a_osc in range.
- CosmomcTheta: omdmh2 += omaxh2 when DM-like.
- SetTimeSteps: refinement window around tau_osc (port AxiECAMB modules.f90:3044-3073:
  dtauosc = tau_osc/int(6000*dfac), ±6.5 steps fine window + coarser taurend→tau_osc).
- Thermo_Init end: cache opacity and exp(-kappa) at tau_osc into Axion%opac_tauosc /
  expmmu_tauosc (modern interpolation replaces Nov13 ThermoSplineOut).
- Transfer module: Transfer_axion = 14, Transfer_max bumped, name tag 'axion'.
- MatterPowerData_k extrapolation beyond kmax: when axion active, clamp (no power-law
  extrapolation of exponentially suppressed P(k)) [AxiECAMB ACCURACY change].
- a_verydom guard / RedshiftAtTimeArr ``om``: include axion as matter when DM-like and
  a > a_osc? — AxiECAMB did not change these (verify equations report); skip unless
  report says otherwise.

4. equations.f90 (details per the equations_ppf report below)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- EvolutionVars: ``axion_ix``, ``AxionIsFluid`` flag, cached tau_osc.
- SetupScalarArrayIndices: +2 equations when axion active (both phases use 2).
- GaugeInterface_EvolveScal: new switch at tau_osc -> set AxionIsFluid, re-setup
  indices, copy variables, apply the KG→EFA transformation (verbatim matching formulas
  from equations_ppf report), ind=1 dverk restart.
- derivs/output: axion background density+pressure into grho/gpres; dgrho/dgq
  contributions in KG phase (from v1_bg,v2_bg lookups + delta-phi vars) and EFA phase
  (delta_ax, u_ax with cs2(k,m,a)); metric_delta capture at tau_osc for the cmbmain
  boundary term; Transfer_axion column fill; Transfer_tot/nonu include axion when
  DM-like; initial(): adiabatic axion ICs (mode 6 NOT wired — v1.0 parity).
- derivst (tensors): axion in background grho (density+pressure) only.
- dtauda: automatic via grho_no_de.

5. cmbmain.f90
~~~~~~~~~~~~~~

- Port the switch boundary-term machinery (deltaBCSrc): capture EV%metric_delta(2) per
  k after CalcScalarSources; spline over k in InitSourceInterpolation; interpolate in
  InterpolateSources; add to temperature transfer sums (flat: sums(1) term with
  expmmu_tauosc, opac_tauosc, 11/10 factor and dotJl_osc; curved: DoRangeInt straddle).
  Gate: has_switch .and. tau_osc < tau0.
- GetTauStart: taustart = min(..., 0.3*tau_osc, tauosc-estimate, taueq) + warning if
  before background table start.
- use_cl_spline_template: force .false. when axion active (AxiECAMB removed templated
  C_l interpolation; LCDM template biases axion spectra).

6. recfast.f90
~~~~~~~~~~~~~~

- dHdz used in the tightly-coupled Tmat term: when axion active replace analytic dHdz
  with FD derivative of Hz(z) from dtauda, with the FD step directed AWAY from a_osc
  (kink). Everything else: no changes (modern recfast already uses dtauda for Hz).

7. halofit.f90
~~~~~~~~~~~~~~

- omm0: include omegaax when DM-like (matter); when DE-like add to omega_v inputs.
- Print warning when NonLinear active with axions ("inherited from axionCAMB, not
  extensively tested"; recommend halofit_version=1/original — Takahashi unstable).

8. camb.f90 (ini interface)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

- CAMB_ReadParams: allocate TAxionModel, call ``P%Axion%ReadParams(Ini)``: keys m_ax
  (negative -> 10**m_ax), use_axfrac, omaxh2, omdah2, axfrac, axion_isocurvature
  (warn-disable), alpha_ax (ignored; computed), Hinf (log10 GeV -> /Mpl with
  Mpl=2.435e18 GeV). Activity: m_ax > 0 and (use_axfrac ? axfrac*omdah2 : omaxh2) > 0
  (DE-like: axfrac>0 suffices).
- inifiles/params_axion.ini sample mirroring AxiECAMB params.ini fiducial.

9. Build & python
~~~~~~~~~~~~~~~~~

- Makefile_main: add ``AxionBackground`` to SOURCEFILES between massive_neutrinos/model
  slots (module used by model.f90 → must precede model).
- camb/axion.py + model.py changes; transfer_names + Transfer_max in model.py.
- README_AxiECAMB.md at repo root: physics, usage (ini + python), parity notes,
  warnings carried over (isocurvature disabled, growth-rate disabled, z>0 transfer
  before switch warning, nonlinear untested, accuracy_boost caveat).

10. Explicitly NOT ported (with reasons)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- cmbmainOMP.f90 (dead code, stale), writefits.f90 changes (broken/dead),
  camb.f90 call_again hack (Nov13 global-state), eV->elecV rename, ALens_Fiducial
  (already in 1.6.7), xe boxcar smoothing in inithermo (disabled no-op upstream),
  recfast recdverk + arg threading (modern structure), recfast tol 1e-3 loosening
  (undocumented accuracy regression), helium_fullreion_redshiftstart 5->7 (modern has
  user-settable 5.5 default; note in README), Transfer_f growth-rate column (disabled
  in AxiECAMB v1.0), isocurvature mode 6 wiring (disabled in v1.0; params accepted and
  warned), GetOmegak deletion (modern has no GetOmegak), neutrino bookkeeping moves
  (modern SetParams handles), P allocatable stack hack, NLL_num_redshifts (modern
  nonlinear z-grid differs; revisit only if lensing accuracy tests fail).

11. Validation plan
~~~~~~~~~~~~~~~~~~~

1. Build New_AxiECAMB (serial make; SDKROOT for conda gfortran).
2. Build AxiECAMB (Nov13, gfortran Makefile) — reference.
3. Runs at m_ax ∈ {1e-27 (RD switch), 1e-28, 1e-25 (phase-tuned), 1e-30 (MD switch),
   1e-33 (DE-like)}, axfrac ∈ {1.0, 0.1}; LCDM limit axfrac→0.
4. Compare: (a) background: a_osc, aeq, dfac after tuning, rho_ax(a) tables; (b)
   ratios (C_l axion)/(C_l LCDM) and (P(k) axion)/(P(k) LCDM) between old and new
   code (isolates axion physics from Nov13↔1.6.7 baseline differences); (c) python
   wrapper smoke test: set_axion + get_results == ini run.


Modern CAMB 1.6.7 architecture map (port landing zones)
--------------------------------------------------------

Source tree: ``/Users/vivianmiranda/data/research/WayneHu/rayne/CAMB/`` (fortran in ``fortran/``, python in ``camb/``).
All line numbers refer to that tree. Goal: tell the port designer exactly where every AxiECAMB change lands.
Tag convention used below: [PLUMBING] = architectural fact that constrains the port; [PHYSICS-LANDING] = the place
axion physics will plug in; [ACCURACY-LANDING] = where accuracy knobs live; [PRECEDENT] = an existing pattern that
already implements something analogous (copy that pattern).

1. CAMBparams (fortran/model.f90) and its python mirror (camb/model.py)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1.1 Fortran side
^^^^^^^^^^^^^^^^^

- ``fortran/model.f90:108-200`` — ``type, extends (TCAMBParameters) :: CAMBparams``. It contains (in this exact order):
  flags ``WantCls`` ... ``DoLensing`` (lines 109-120), ``NonLinear`` (121), ``type(TransferParams) :: Transfer`` (122),
  ``want_zstar/want_zdrag`` (124-125), l/eta limits (127-131), then the physical densities:

  .. code-block:: fortran

      real(dl)  :: ombh2 = 0._dl !baryon density Omega_b h^2          ! model.f90:135
      real(dl)  :: omch2 = 0._dl !cold dark matter density Omega_c h^2 ! :136
      real(dl)  :: omk = 0._dl                                         ! :137
      real(dl)  :: omnuh2 = 0._dl                                      ! :138
      real(dl)  :: H0 = 67._dl                                         ! :139

  then TCMB/Yhe/neutrino settings (140-148), then the **allocatable component classes**:

  .. code-block:: fortran

      class(TInitialPower), allocatable :: InitPower        ! model.f90:150
      class(TRecombinationModel), allocatable :: Recomb     ! :151
      class(TReionizationModel), allocatable :: Reion       ! :152
      class(TDarkEnergyModel), allocatable :: DarkEnergy    ! :153
      class(TNonLinearModel), allocatable :: NonLinearModel ! :154
      type(AccuracyParams)     :: Accuracy                  ! :155
      type(SourceTermParams)   :: SourceTerms              ! :156
      real(dl), allocatable :: z_outputs(:)                 ! :158
      integer   :: Scalar_initial_condition = 1 !adiabatic  ! :160
      real(dl), allocatable  :: InitialConditionVector(:)   ! :162

  followed by ``OutputNormalization``, ``Alens``, ``MassiveNuMethod``, ``DoLateRadTruncation``, ``Evolve_baryon_cs``,
  ``Evolve_delta_xe``, ``Evolve_delta_Ts``, ``Do21cm``, ``transfer_21cm_cl``, ``Log_lvalues``, ``use_cl_spline_template``,
  ``min_l_logl_sampling``, ``SourceWindows``, ``CustomSources`` (164-190).

- Required boilerplate per python-visible class: ``PythonClass`` (model.f90:204-207), ``SelfPointer`` (209-218),
  ``Replace`` (220-235). ``CAMBparams`` already has these; any **new** F2003 class (e.g. a TAxionModel component)
  needs its own ``SelfPointer``/``PythonClass`` (see classes.f90:45-55 comments: "All python-accessible inherited
  classes must define SelfPointer, and use @fortran_class decorator in python").

- ``CAMBparams_Validate`` model.f90:341-409 — sanity range checks (``ombh2<0.0005 .or. omch2<0 ...`` at :382).
  [PLUMBING] AxiECAMB's parameter validation (omaxh2/axfrac consistency) belongs here or in the axion class's
  ``Validate``.

1.2 How to add a parameter — the rules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

[PLUMBING] The instructions are literally in the code, ``camb/model.py:271-274``:

    "To add a new parameter, add it to the CAMBparams type in model.f90, then edit the ``_fields_`` list in the
    CAMBparams class in model.py to add the new parameter **in the corresponding location of the member list**."

ctypes layout rules (from ``camb/baseconfig.py`` ``CAMBStructureMeta.__new__``, baseconfig.py:459-591):

1. **Ordering is positional and must exactly match the Fortran declaration order** — the python ``_fields_`` tuple
   is converted 1:1 into a ctypes ``Structure._fields_`` (baseconfig.py:538 ``namespace["_fields_"] = ctypes_fields``).
   There is no name matching against Fortran; any mismatch silently shifts every later field.
2. **Booleans**: declare ``c_bool`` in python; the metaclass stores them as ``c_int`` internally
   (``BoolField``, baseconfig.py:424-434; conversion at :484-486 ``ctypes_fields.append(("_" + field_name, c_int))``).
   Fortran side must be ``logical`` (default kind, 4 bytes).
3. **Integers**: ``c_int`` ↔ Fortran default ``integer``. Optional named-value mapping via
   ``("NonLinear", c_int, {"names": NonLinear_names})`` → ``NamedIntField`` (baseconfig.py:393-421, model.py:295).
4. **Reals**: ``c_double`` ↔ ``real(dl)``.
5. **Fixed-size arrays**: ``c_double * max_nu`` with ``{"size": "nu_mass_eigenstates"}`` → ``SizedArrayField``
   (baseconfig.py:437-456); the size key must be the name of another field in the same structure (:500-505).
6. **Allocatable arrays**: Fortran ``real(dl), allocatable :: x(:)`` ↔ python ``AllocatableArrayDouble``
   (baseconfig.py:313-321); ``integer, allocatable`` ↔ ``AllocatableArrayInt`` (:302-310). These are opaque
   descriptors whose byte size is queried from Fortran at import time (``_get_fortran_sizes``, baseconfig.py:138-154)
   — so gfortran/ifort descriptor differences are handled automatically.
7. **Allocatable class members**: Fortran ``class(TDarkEnergyModel), allocatable :: DarkEnergy`` ↔
   ``("DarkEnergy", AllocatableObject(DarkEnergyModel))`` (model.py:341; machinery baseconfig.py:173-204).
   Assignment goes through ``FortranManagedField.__set__`` → ``_set_allocatable`` (Fortran ``handles`` module in
   camb_python.f90:204-213), which sources a copy of the python-held Fortran object into the param's allocatable.
8. **Nested plain structures** (TransferParams, AccuracyParams, SourceTermParams, CustomSources) are
   ``CAMB_Structure`` subclasses embedded by value (model.py:296, 343, 344, 391).
9. Fields whose python name starts with ``__`` become private/hidden (e.g. ``("__is_cosmological_constant", c_bool)``
   dark_energy.py:11), but **still occupy space and must be present**.
10. Methods callable from python are declared in ``_methods_`` (model.py:399-407) and must exist in Fortran as
    ``module procedure`` named ``<ClassName>_<MethodName>`` in module ``_fortran_class_module_`` (= ``"model"`` for
    CAMBparams, model.py:397). Lookup pattern: ``__<module>_MOD_<classname><funcname>`` (baseconfig.py:96).
11. New python classes bind with the ``@fortran_class`` decorator (baseconfig.py:820-840), which resolves the
    Fortran ``SelfPointer`` and registers the class in ``F2003Class._class_names[cls.__name__]`` (baseconfig.py:840)
    — this registry is what ``make_class_named`` / ``set_classes(dark_energy_model="...")`` use.

The full python ``_fields_`` mirror of CAMBparams is ``camb/model.py:284-392`` (order matches model.f90 exactly:
WantCls ... CustomSources). Scalar_initial_condition is mirrored as a named int:
``("scalar_initial_condition", c_int, {"names": ["initial_vector", "initial_adiabatic", "initial_iso_CDM", "initial_iso_baryon", "initial_iso_neutrino", "initial_iso_neutrino_vel"]})``
(model.py:346-359; note ``initial_vector`` = 0 so list starts at 0).
[PHYSICS-LANDING] A new mode 6 (``initial_iso_axion``) means appending a name to this list **and** raising
``initial_nummodes`` in equations.f90 (see §4.4).

[PLUMBING] **Where axion params can live**: two options. (a) Add ``m_ax, omaxh2, axfrac, use_axfrac, ...`` directly
to CAMBparams in model.f90 + model.py ``_fields_`` (simple, AxiECAMB-like). (b) Put them in a new
``TCambComponent``/``TDarkEnergyModel``-derived class allocated inside CAMBparams (the modern, preferred pattern —
TEarlyQuintessence/TAxionEffectiveFluid precedent, §3). Isocurvature params ``alpha_ax``, ``Hinf`` interact with the
initial power amplitude and initial conditions, so they likely belong in CAMBparams or InitPower.

2. Background expansion: where rho_tot/H(a) is computed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

2.1 The dtauda/grho chain (single chokepoint)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``fortran/equations.f90:4-23`` — top-level function, **outside any module**, used everywhere via
``procedure(obj_function), private :: dtauda`` interface declarations (results.f90:293, equations.f90:193,
recfast.f90:348):

.. code-block:: fortran

    function dtauda(this,a)
    use results
    use DarkEnergyInterface
    implicit none
    class(CAMBdata) :: this
    real(dl), intent(in) :: a
    real(dl) :: dtauda, grhoa2, grhov_t

    call this%CP%DarkEnergy%BackgroundDensityAndPressure(this%grhov, a, grhov_t)

    !  8*pi*G*rho*a**4.
    grhoa2 = this%grho_no_de(a) +  grhov_t * a**2
    if (grhoa2 <= 0) then
        call GlobalError('Universe stops expanding before today (recollapse not supported)', error_unsupported_params)
        dtauda = 0
    else
        dtauda = sqrt(3 / grhoa2)
    end if
    end function dtauda

``fortran/results.f90:1224-1241`` — the non-DE part:

.. code-block:: fortran

    function grho_no_de(this, a) result(grhoa2)
    !  Return 8*pi*G*rho_no_de*a**4 where rho_no_de includes everything except dark energy.
    class(CAMBdata) :: this
    ...
    grhoa2 = this%grhok * a**2 + (this%grhoc + this%grhob) * a + this%grhog + this%grhornomass
    if (this%CP%Num_Nu_massive /= 0) then
        do nu_i = 1, this%CP%nu_mass_eigenstates
            call ThermalNuBack%rho(a * this%nu_masses(nu_i), rhonu)
            grhoa2 = grhoa2 + rhonu * this%grhormass(nu_i)
        end do
    end if
    end function grho_no_de

[PHYSICS-LANDING] **If the axion is implemented as the DarkEnergy component** (or a sibling component combined
into ``BackgroundDensityAndPressure``), H(a) everywhere is automatically correct: every distance, age, thermal
history, recombination and perturbation routine flows through ``dtauda`` → ``grho_no_de + DarkEnergy``. If instead
the axion is a *separate* new component (DE stays Λ), ``grho_no_de`` (results.f90:1231) and/or ``dtauda`` itself are
the two functions to modify — but then you must also touch every place listed in §2.3.

Units: ``grho*`` quantities are 8πG·ρ·a^n/c² in Mpc⁻² with a=1 today; ``grhoa2`` is 8πGρa⁴; ``dtauda = sqrt(3/grhoa2)``.

2.2 Constant densities set in CAMBdata_SetParams (results.f90:313-...)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``fortran/results.f90:443-467``:

.. code-block:: fortran

    this%grhocrit = 3*this%CP%h0**2/c**2*1000**2 !3*h0^2/c^2 (=8*pi*G*rho_crit/c^2)
    this%grhog = kappa/c**2*4*sigma_boltz/c**3*this%CP%tcmb**4*Mpc**2
    this%grhor = 7._dl/8*(4._dl/11)**(4._dl/3)*this%grhog
    this%grhornomass=this%grhor*nu_massless_degeneracy
    this%grhormass(nu_i)=this%grhor*this%CP%Nu_mass_degeneracies(nu_i)   ! :456
    h2 = (this%CP%H0/100)**2
    this%grhoc=this%grhocrit*this%CP%omch2/h2
    this%grhob=this%grhocrit*this%CP%ombh2/h2
    this%grhok=this%grhocrit*this%CP%omk
    this%Omega_de = 1 -(this%CP%omch2 + this%CP%ombh2 + this%CP%omnuh2)/h2 - this%CP%omk  &
        - (this%grhornomass + this%grhog)/this%grhocrit                  ! :462-463
    this%grhov=this%grhocrit*this%Omega_de                               ! :464
    this%adotrad = sqrt((this%grhog+this%grhornomass+sum(this%grhormass(...)))/3)  ! :467

[PHYSICS-LANDING] **Closure/budget**: ``Omega_de`` at results.f90:462 is where "1 − everything else" is computed.
An axion with ``omaxh2`` must be subtracted here (or folded into the DE component so that grhov covers Λ+axion as
AxiECAMB's omdah2 split does). ``DarkEnergy%Init(this)`` is called at results.f90:497 *after* the grho* are set and
massive ν initialized — that is the precise hook where a KG background solver can run (TEarlyQuintessence does
exactly this, §3.3). There is **no GetOmegak function in 1.6.7**; curvature handling is ``this%CP%omk``,
``this%grhok`` (:461) and flat/closed flags at results.f90:422-433.

Other budget-sensitive places (all in results.f90):

- ``z_eq`` (defined "assuming all neutrinos massless") :476 — ``(grhob+grhoc)/(grhog+grhornomass+Σgrhormass) − 1``.
  [PHYSICS] If the axion contributes to "matter", AxiECAMB redefined z_eq; the modern analog is here and in the
  derived parameter ``ThermoDerivedParams(derived_zEQ)`` results.f90:2186-2199 (``a_eq``, ``k_EQ = 1/(a_eq*dtauda)``).
- ``CosmomcTheta`` :917-937 — uses ``omdmh2 = (this%CP%omch2+this%CP%omnuh2)`` in the Hu–Sugiyama zstar fit (:924-929).
  If omaxh2 acts as matter it must be added to ``omdmh2`` for theta-based H0 setting to be meaningful.
- ``RedshiftAtTimeArr`` :876-884 and ``initial`` (equations.f90:1830) use
  ``om = (grhob+grhoc)/sqrt(3*(grhog+grhornomass+Σgrhormass))`` — early-matter-domination expansions; an axion that
  is matter-like *before* recombination would enter ``om`` (AxiECAMB touched the analogous Nov13 code).
- ``Thermo_Init`` matter_verydom: results.f90:1814 ``a_verydom = AccuracyBoost*5*(grhog+grhornomass)/(grhoc+grhob)``.
- Feedback printout :570-579 prints ``Om_m (inc Om_u) = (ombh2+omch2+omnuh2)/h2``.
- ``GetBackgroundDensities`` results.f90:940-972 returns ``densities(8,n)`` rows
  ``[tot, K, cdm, baryon, photon, neutrino(massless), nu(massive), de]``; python mirror is
  ``model.py:125 density_names = ["tot","K","cdm","baryon","photon","neutrino","nu","de"]`` and
  ``results.py get_background_densities``. [PLUMBING] Adding an axion row means changing the ``densities(8,n)``
  shape, the row assignments (:962-969), and the python density_names list together.

2.3 Background tables / hooks for tabulated densities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``TCAMBdata`` is an empty placeholder in classes.f90:77-79; the real ``CAMBdata`` is results.f90:164-281.
  It stores: constant ``grhocrit,grhog,grhor,grhob,grhoc,grhov,grhornomass,grhok`` (:177), ``Omega_de`` (:180),
  curvature vars (:181-183), ``grhormass(max_nu)``/``nu_masses(max_nu)`` (:190-192), ``ThermoData`` (TThermoData,
  :220), ``BackgroundOutputs`` (TBackgroundOutputs ``H(:), DA(:), rs_by_D_v(:)`` results.f90:39-41), ``MT``
  (MatterTransferData :225), ``CLdata`` (:230).
- ``TThermoData`` (results.f90:63-97) holds the only generic tabulated background: log-τ-spaced arrays
  ``ScaleFactor, adot`` etc., built in ``Thermo_Init`` (results.f90:1681 ff) by direct calls to ``dtauda``:
  ``dt(i) = dtauda(State,scale_factors(i))`` (:1844), ``adot = 1/dtauda(State,a)`` (:1869), splined into
  ``ScaleFactorAtTime`` (TCubicSpline, :1846-1851).
- [PLUMBING] **There is no generic hook for an extra tabulated density component.** The supported pattern is:
  the component class itself owns its interpolation tables (TEarlyQuintessence keeps ``sampled_a, phi_a,
  phidot_a, ddphi_a, ddphidot_a, fde, ddfde`` inside the DE object — DarkEnergyQuintessence.f90:37-41, 68) and
  exposes density via ``BackgroundDensityAndPressure(grhov, a, grhov_t, w)``. The axion background table
  (KG solution before switch + w=0 fluid after) should live inside an axion class the same way.
- Derived parameters: ``nthermo_derived = 13`` constants at results.f90:43-46
  (``derived_age...derived_theta_rs_EQ``), filled at results.f90:2172-2199; python names
  ``model.py:68-82 derived_names = ["age","zstar",...,"thetarseq"]``. Adding a derived param (e.g. AxiECAMB's
  switch redshift) requires bumping ``nthermo_derived`` in both languages plus ``model.py:41 nthermo_derived = 13``.

2.4 DarkEnergy background coupling interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``fortran/DarkEnergyInterface.f90:83-103``:

.. code-block:: fortran

    subroutine BackgroundDensityAndPressure(this, grhov, a, grhov_t, w)
    !Get grhov_t = 8*pi*rho_de*a**2 and (optionally) equation of state at scale factor a
    ...
    if (this%is_cosmological_constant) then
        grhov_t = grhov * a * a
        if (present(w)) w = -1_dl
    else
        if (a > 1e-10) then
            grhov_t = grhov * this%grho_de(a) / (a * a)
        else
            grhov_t = 0._dl
        end if
        if (present(w)) w = this%w_de(a)
    end if

with overridables ``w_de(a)`` (default −1, :51-58) and ``grho_de(a)`` ("relative density 8 pi G a^4 rho_de /grhov",
default 0, :60-67). ``TDarkEnergyEqnOfState`` (:26-46) implements w0–wa and tabulated w(a):
``grho_de = a**(1-3*w_lam-3*wa) * exp(-3*wa*(1-a))`` (:223-224) and spline ``logdensity`` from
``log(rho) = -3 ∫ dlna (1+w)`` (:174-179). ``Effective_w_wa`` (:105-112) feeds halofit.

3. The component/class system and the quintessence precedent
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

3.1 classes.f90
^^^^^^^^^^^^^^^^

- ``TPythonInterfacedClass`` classes.f90:45-55 (SelfPointer/PythonClass/Replace contract).
- ``TCambComponent`` classes.f90:66-70 — adds ``ReadParams(Ini)`` and ``Validate(OK)``; every pluggable model
  (TNonLinearModel :81-87, TInitialPower :89-95, TRecombinationModel :97-108, TReionizationModel :110-116,
  TDarkEnergyModel) extends it.
- ``TClassDverk`` interface classes.f90:119-128 — class-aware dverk used by quintessence background integration
  (``procedure(TClassDverk) :: dverk`` DarkEnergyQuintessence.f90:82). There is no ``TCAMB_Calculation`` type in
  1.6.7; the orchestration object is ``CAMBdata`` (results.f90:164) driven by ``cmbmain.f90`` and ``camb.f90``
  (``CAMB_GetResults`` camb.f90:50-151).

3.2 TDarkEnergyModel perturbation hooks (DarkEnergyInterface.f90)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: fortran

    type, extends(TCambComponent) :: TDarkEnergyModel      ! DarkEnergyInterface.f90:9-24
        logical :: is_cosmological_constant = .true.
        integer :: num_perturb_equations = 0
    contains
        procedure :: Init
        procedure :: BackgroundDensityAndPressure
        procedure :: PerturbedStressEnergy !Get density perturbation and heat flux for sources
        procedure :: diff_rhopi_Add_Term
        procedure :: PerturbationInitial
        procedure :: PerturbationEvolve
        ...

Signatures (quoted):

- ``PerturbedStressEnergy(this, dgrhoe, dgqe, a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)`` (:115-127) — returns the component's ``dgrhoe = a²κδρ``, ``dgqe = a²κ(ρ+p)v`` to be added into the metric sums.
- ``PerturbationEvolve(this, ayprime, w, w_ix, a, adotoa, k, z, y)`` (:144-149) — writes ``ayprime(w_ix...)``.
- ``PerturbationInitial(this, y, a, tau, k)`` (:151-160) — initial values for the component's perturbation block
  ("can usually just set to zero").
- ``diff_rhopi_Add_Term(...) result(ppiedot)`` (:130-142) — anisotropic-stress time-derivative addition for source
  terms (used by PPF; default 0).

3.3 TEarlyQuintessence — the closest analog to the axion KG background solver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``fortran/DarkEnergyQuintessence.f90``:

- Base ``TQuintessence`` (:33-53): fields ``astart = 1e-7_dl``, ``integrate_tol = 1e-6_dl``, table arrays
  ``sampled_a, phi_a, phidot_a`` (+private ``ddphi_a, ddphidot_a``, grid bookkeeping ``npoints_linear, npoints_log,
  dloga, da, log_astart, max_a_log``), ``class(CAMBdata), pointer, private :: State``.
- ``TQuintessence_Init`` (:113-131): grabs the State pointer (``select type(State); class is (CAMBdata)`` :121-124),
  sets ``is_cosmological_constant = .false.``, ``num_perturb_equations = 2``.
- Background ODE in a (and log a wrapper :159-172), ``EvolveBackground`` :174-195 — variables ``phi=y(1)``,
  ``a² phi' = y(2)``:

  .. code-block:: fortran

      grhode=a2*(0.5d0*phidot**2 + a2*this%Vofphi(phi,0))
      tot = this%state%grho_no_de(a) + grhode
      adot=sqrt(tot/3.0d0)
      yprime(1)=phidot/adot !d phi /d a
      yprime(2)= -a2**2*this%Vofphi(phi,1)/adot

  Note it reuses ``State%grho_no_de(a)`` so the KG solve is self-consistent with all other components — exactly
  what the AxiECAMB background solver needs.
- **Grid setup & oscillation sampling** (``TEarlyQuintessence_Init`` :301-551):

  - ``npoints = 5000`` baseline log-a steps, ``min_steps_per_osc = 10`` (:66-67).
  - ``this%dloga = (-this%log_astart)/(this%npoints-1)`` (:383); log spacing from astart, switching to linear
    where the step matches: ``this%max_a_log = 1.d0/this%npoints/(exp(this%dloga)-1)`` (:386).
  - integrates with ``dverk(this,NumEqs,EvolveBackgroundLog,afrom,y,aend,this%integrate_tol,ind,c,NumEqs,w)``
    (:463); detects oscillations by sign changes of y(2):

    .. code-block:: fortran

        elseif (y(2)*lastsign < 0) then
            !derivative has changed sign. Use to probe any oscillation scale:
            da_osc = min(da_osc, exp(aend) - last_a)        ! :470-474

    and aborts log spacing when
    ``if (sampled_a(ix)*(exp(this%dloga)-1)*this%min_steps_per_osc > da_osc) exit`` (:483-486).
  - Then linear-in-a stage with
    ``this%da = min(this%max_a_log*(exp(this%dloga)-1), da_osc/this%min_steps_per_osc, (1-this%max_a_log)/(this%npoints-this%npoints_log))`` (:492-495).
  - Splines: ``call spline(this%sampled_a,this%phi_a,tot_points,splZero,splZero,this%ddphi_a)`` and same for
    ``phidot_a``, ``fde`` (:527-529).
  - **fde** (early-DE fraction) tabulated during the solve (:477-479):

    .. code-block:: fortran

        fde(ix) = 1/((this%state%grho_no_de(sampled_a(ix)) +  this%frac_lambda0*this%State%grhov*a2**2) &
            /(a2*(0.5d0* phidot_a(ix)**2 + a2*this%Vofphi(y(1),0))) + 1)

  - Peak finding by spline-derivative root (``fde_peak`` :567-598) → ``zc = 1/a_c - 1``, ``fde_zc = fdeAta(a_c)`` (:540-546).
- **Shooting/matching** (the analog of AxiECAMB's omaxh2 shooting): if ``use_zc`` (:327-381), minimizes
  ``match_fde_zc = (log(this%fde_zc)-log(fde_zc))**2 + (log(zc)-log(this%zc))**2`` (:631) over (log f, log m) with
  ``Minimize%NEWUOA(this, match_fde_zc, 2, 5, log_params, 0.8_dl, 1e-4_dl, this%DebugLevel, 500)`` (:360-361),
  using ``calc_zc_fde`` (:638-709) which forward-integrates until the fde peak. Failure sets
  ``global_error_flag = error_darkenergy`` (:364, 376). PowellMinimize.f90 supplies NEWUOA/BOBYQA; brentq is also
  available (:338-345). [PRECEDENT] AxiECAMB's "shoot for omaxh2 given m_ax" maps directly onto this structure.
- **Interpolation accessor** ``ValsAta(this,a,aphi,aphidot)`` (:206-235): index lookup is *analytic* (no search):
  log region ``ix = int((log(a)-log_astart)/dloga)+1``, linear region ``ix = npoints_log + int((a-max_a_log)/da)``;
  then in-place cubic-spline formula (:228-233). For ``a < astart``: ``aphi = phi_a(1), aphidot = 0`` (:217-220).
- **Potential** ``TEarlyQuintessence_VofPhi`` (:276-298): units comment "input variable phi is sqrt(8*Pi*G)*psi ...
  return result is in 1/Mpc^2 units":

  .. code-block:: fortran

      real(dl), parameter :: units = MPC_in_sec**2 /Tpl**2  !convert to units of 1/Mpc^2  ! :284
      theta = phi/this%f
      if (deriv==0) then
          V = units*this%m**2*this%f**2*(1 - cos(theta))**this%n + this%frac_lambda0*this%State%grhov
      else if (deriv ==1) then
          V = units*this%m**2*this%f*this%n*(1 - cos(theta))**(this%n-1)*sin(theta)
      else if (deriv ==2) then
          V = units*this%m**2*this%n*(1 - costheta)**(this%n-1)*(this%n*(1+costheta) -1)

  with ``Tpl = sqrt(kappa*hbar/c**5)`` (:29) and ``m`` in reduced-Planck units. [PHYSICS] AxiECAMB's m_ax in eV must
  be converted consistently (m[Mpc⁻¹] = m_eV/ħ·Mpc/c etc.) — this file shows the house unit conventions.
- **Quintessence perturbations** (:237-272) — field perturbations evolved in the DE block:

  .. code-block:: fortran

      clxq=ay(w_ix); vq=ay(w_ix+1)
      dgrhoe= phidot*vq +clxq*a**2*this%Vofphi(phi,1)
      dgqe= k*phidot*clxq                                   ! PerturbedStressEnergy :248-252
      ayprime(w_ix)= vq
      ayprime(w_ix+1) = - 2*adotoa*vq - k*z*phidot - k**2*clxq - a**2*clxq*this%Vofphi(phi,2)  ! :269-270

  So here ``y(w_ix)=δφ (kappa-normalized), y(w_ix+1)=δφ'``. [PRECEDENT] AxiECAMB's exact-KG perturbation phase
  maps onto this; the post-switch effective-fluid phase maps onto TDarkEnergyFluid/TAxionEffectiveFluid (§3.4).

3.4 TAxionEffectiveFluid — precedent for the post-switch effective fluid
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``fortran/DarkEnergyFluid.f90:24-41`` declares ``TAxionEffectiveFluid`` (w_n, fde_zc, zc, theta_i; cached
``a_c, pow, om, omL, acpow, freq, n``); ``num_perturb_equations = 2`` (:200). Key physics worth quoting because it is
nearly the same scale-dependent sound-speed structure AxiECAMB uses after the switch
(``TAxionEffectiveFluid_PerturbationEvolve``, :249-276):

.. code-block:: fortran

    if (this%w_n < 0.9999) then
        fac = 2*a**(2-6*this%w_n)*this%freq**2
        cs2 = (fac*(this%n-1) + k**2)/(fac*(this%n+1) + k**2)
    else
        cs2 = 1
    end if
    ...
    Hv3_over_k =  3*adotoa* y(w_ix + 1) / k
    deriv  = (acpow**2*(this%om+this%omL)+this%om*acpow-apow**2*this%omL)*this%pow &
        /((apow+acpow)*(this%omL*(apow+acpow)+this%om*(1+acpow)))   ! dw/dloga/(1+w)
    ayprime(w_ix) = -3 * adotoa * (cs2 - w) *  (y(w_ix) + Hv3_over_k) &
        -   k * y(w_ix + 1) - (1 + w) * k * z  - adotoa*deriv* Hv3_over_k
    ayprime(w_ix + 1) = -adotoa * (1 - 3 * cs2 - deriv) * y(w_ix + 1) + k * cs2 * y(w_ix)

Here ``y(w_ix)=δ``, ``y(w_ix+1)=(1+w)v``; ``PerturbedStressEnergy`` returns ``dgrhoe = ay(w_ix)*grhov_t``,
``dgqe = ay(w_ix+1)*grhov_t`` (:288-289). The standard fluid version (TDarkEnergyFluid_PerturbationEvolve,
:117-146) uses ``cs2_lam`` and ``dgqe = ay(w_ix+1) * grhov_t * (1 + w)`` (:112).
Background: ``grho_de(a) = (omL*(a^pow+acpow)+om*(1+acpow))*a^4/((a^pow+acpow)*(omL+om))`` (:234-247),
``w_de`` (:221-232). The Init shows how to use ``State%grho_no_de(a_c)`` and ``State%Omega_de`` to convert fde_zc to
Omegas (:196-199).

4. equations.f90 (GaugeInterface): EvolutionVars, indices, sums, ICs, switches
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

4.1 Structure and index assignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Module ``GaugeInterface`` equations.f90:29; gauge tag ``Eqns_name = 'cdm_gauge'`` (:42).
- Fixed indices: ``integer, parameter :: basic_num_eqns = 4`` and
  ``ix_etak=1, ix_clxc=2, ix_clxb=3, ix_vb=4`` (:46-47).
- ``type EvolutionVars`` (:66-151): per-k state holding **dynamic index assignments**: ``w_ix`` ("Index of two
  quintessence equations", :71), ``Tg_ix, reion_line_ix, xe_ix, Ts_ix`` (:72-76), ``r_ix`` (massless ν), ``g_ix``
  (photons) (:78-79), ``nvar/nvart/nvarv`` (:85), lmax bookkeeping (:88-92), ``polind`` (:94), massive-ν
  ``nu_ix(max_nu), nu_pert_ix, nq(max_nu)`` (:97-98), approximation flags
  ``MassiveNuApprox/MassiveNuApproxTime`` (:105-106), ``high_ktau_neutrino_approx`` (:109),
  ``TightCoupling, TensTightCoupling, TightSwitchoffTime`` (:115-116), ``ScalEqsToPropagate`` (:119),
  ``no_nu_multpoles, no_phot_multpoles`` (:129), ``saha, evolve_TM, evolve_baryon_cs`` (:139-141), output pointers
  ``OutputTransfer/OutputSources/CustomSources`` (:146-148).
- ``SetupScalarArrayIndices(EV, max_num_eqns)`` (:522-636) assigns indices in order: photons (g_ix=5, then
  polarization via polind), massless ν (r_ix), **then dark energy**:

  .. code-block:: fortran

      !Dark energy
      if (.not. CP%DarkEnergy%is_cosmological_constant) then
          EV%w_ix = neq + 1
          neq = neq + CP%DarkEnergy%num_perturb_equations
          maxeq = maxeq + CP%DarkEnergy%num_perturb_equations
      else
          EV%w_ix = 0
      end if                                                  ! equations.f90:555-562

  then xe/Tg/Ts blocks, then massive ν blocks (:592-629), finally ``EV%ScalEqsToPropagate = neq`` (:631).
  [PLUMBING] An axion perturbation block, if not implemented as the DarkEnergy component, would be added here
  with its own ``ax_ix`` exactly parallel to ``w_ix``. If it *is* the DE component, you just set
  ``num_perturb_equations`` (can differ before/after switch only via the copy mechanism — see §4.5 risk note).

4.2 Where dgrho/dgq sums happen (subroutine derivs, :2140 ff)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: fortran

    grhob_t=State%grhob/a
    grhoc_t=State%grhoc/a
    grhor_t=State%grhornomass/a2
    grhog_t=State%grhog/a2                                       ! :2202-2205
    if (EV%is_cosmological_constant) then
        grhov_t = State%grhov * a2 ; w_dark_energy_t = -1_dl
    else
        call State%CP%DarkEnergy%BackgroundDensityAndPressure(State%grhov, a, grhov_t, w_dark_energy_t)
    end if                                                       ! :2207-2212
    dgrho_matter=grhob_t*clxb+grhoc_t*clxc                       ! :2216
    dgq=grhob_t*vb                                               ! :2218
    if (State%CP%Num_Nu_Massive > 0) call MassiveNuVars(EV,ay,a,grhonu_t,gpres_nu,dgrho_matter,dgq, wnu_arr) ! :2223-2225
    grho_matter=grhonu_t+grhob_t+grhoc_t                         ! :2227
    grho = grho_matter+grhor_t+grhog_t+grhov_t                   ! :2228
    gpres_noDE = gpres_nu + (grhor_t + grhog_t)/3                ! :2229
    adotoa=sqrt(grho/3)  (flat)                                  ! :2232
    ...
    dgrho=dgrho + grhog_t*clxg+grhor_t*clxr                      ! :2275
    dgq=dgq + grhog_t*qg+grhor_t*qr                              ! :2278
    if (.not. EV%is_cosmological_constant) then
        call State%CP%DarkEnergy%PerturbedStressEnergy(dgrho_de, dgq_de, &
            a, dgq, dgrho, grho, grhov_t, w_dark_energy_t, gpres_noDE, etak, &
            adotoa, k, EV%Kf(1), ay, ayprime, EV%w_ix)
        dgrho = dgrho + dgrho_de
        dgq = dgq + dgq_de
    end if                                                       ! :2284-2290
    z=(0.5_dl*dgrho/k + etak)/adotoa                             ! :2294
    ayprime(ix_etak)=0.5_dl*dgq          (flat)                  ! :2298
    if (.not. EV%is_cosmological_constant) &
        call State%CP%DarkEnergy%PerturbationEvolve(ayprime, w_dark_energy_t, EV%w_ix, a, adotoa, k, z, ay)  ! :2304-2306

Output-time sums (sources/transfer): ``dgpi = grhor_t*pir + grhog_t*pig`` (:2680), DE anisotropic-stress hook in
``diff_rhopi`` via ``DarkEnergy%diff_rhopi_Add_Term(...)`` (:2689-2692), Weyl potential
``phi = -((dgrho +3*dgq*adotoa/k)/EV%Kf(1) + dgpi)/(2*k2)`` (:2693).
[PHYSICS-LANDING] AxiECAMB's axion contributions to δρ, (ρ+p)v (and zero anisotropic stress) all flow through the
single ``PerturbedStressEnergy`` call if the axion is the DE component — no scattered edits needed (contrast Nov13).

Equivalent hook sites also exist in ``output`` (around :2675-2731, used for transfer/source output) and in
``MassiveNuVarsOut`` (:1011-1071, ``dgrho = a^2 kappa \delta\rho`` comment at :1018). Tensor (``derivst``) and vector
(``derivsv``) equations do not call DarkEnergy perturbation hooks (DE has no tensor sources); the tensor background
uses ``grhov_t`` via the same BackgroundDensityAndPressure call pattern.

4.3 Initial conditions (subroutine initial, :1760-1976)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Local mode-vector indices: ``i_clxg=1,i_clxr=2,i_clxc=3,i_clxb=4,i_qg=5,i_qr=6,i_vb=7,i_pir=8,i_eta=9,
  i_aj3r=10,i_clxde=11,i_vde=12``; ``i_max=i_vde``; ``real(dl) initv(6,1:i_max), initvec(1:i_max)`` (:1771-1774).
  Note **initv is already dimensioned 6** but only rows 1–5 are filled; row 6 is free.
- Supported modes: ``initial_adiabatic=1, initial_iso_CDM=2, initial_iso_baryon=3, initial_iso_neutrino=4,
  initial_iso_neutrino_vel=5, initial_vector=0``; ``initial_nummodes = initial_iso_neutrino_vel`` (:62-64).
  Guard: ``if (CP%Scalar_initial_condition > initial_nummodes) call MpiStop('Invalid initial condition...')`` (:1839).
- Superhorizon expansion variables: ``x=k*tau``, ``omtau``, ``Rv=grhonu/(grhonu+grhog)``, ``Rc=omch2/(omch2+ombh2)``
  (:1824-1837); e.g. adiabatic ``initv(1,i_clxg)=-chi*EV%Kf(1)/3*x2*(1-omtau/5)`` (:1850). There is even a stub
  comment ``!quintessence isocurvature mode`` at :1908.
- Mode selection: ``InitVec = initv(CP%Scalar_initial_condition,:)``, adiabatic sign flip (:1916-1920); mixing via
  ``CP%InitialConditionVector(i)`` for ``initial_vector`` (:1911-1915).
- Mapping into y: ``y(ix_etak)= -InitVec(i_eta)*k/2``, ``y(ix_clxc)``, ``y(ix_clxb)``, ``y(ix_vb)``, photons
  ``y(EV%g_ix)=InitVec(i_clxg)``, ``y(EV%g_ix+1)=InitVec(i_qg)`` (:1922-1934); **dark energy block initial values**:

  .. code-block:: fortran

      if (CP%DarkEnergy%num_perturb_equations > 0) then
          call CP%DarkEnergy%PerturbationInitial(InitVec(i_clxde:i_clxde + CP%DarkEnergy%num_perturb_equations - 1), &
              a, tau,  k)
          y(EV%w_ix:EV%w_ix + CP%DarkEnergy%num_perturb_equations - 1) = &
              InitVec(i_clxde:i_clxde + CP%DarkEnergy%num_perturb_equations - 1)
      end if                                                   ! :1938-1943

  neutrinos ``y(EV%r_ix..)`` (:1950-1956), massive ν copies + ``MassiveNuApproxTime`` (:1960-1974).

[PHYSICS-LANDING] **Adding axion isocurvature (mode 6)**: (a) bump ``initial_nummodes`` to a new
``initial_iso_axion=6`` (:62-64); (b) fill ``initv(6,:)`` (dimension already 6, but increase if also keeping rows
1–5; safer: ``initv(initial_nummodes,1:i_max)``); (c) the axion field perturbation IC goes through
``DarkEnergy%PerturbationInitial`` (which receives ``a, tau, k`` — enough for the AxiECAMB δφ ∝ alpha_ax·Hinf form,
but it does *not* receive the mode number, so the axion class needs to know ``Scalar_initial_condition``, e.g. via
its State pointer ``State%CP%Scalar_initial_condition``); (d) python mirror: append ``"initial_iso_axion"`` to the
names list in model.py:350-358; (e) ini file: ``initial_condition`` is read at camb.f90:578.
[RISK] ``use_cl_spline_template`` is auto-disabled for non-adiabatic ICs (camb.f90:584).

4.4 The switch mechanism (precedent for the KG→fluid mid-evolution switch)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``recursive subroutine GaugeInterface_EvolveScal(EV,tau,y,tauend,tol1,ind,c,w)`` equations.f90:245-433 is the
canonical pattern. It computes candidate switch times, takes the min, integrates to it, rebuilds the index
layout, copies variables, optionally transforms them, then recurses:

.. code-block:: fortran

    next_switch = min(tau_switch_ktau, tau_switch_nu_massless,EV%TightSwitchoffTime, tau_switch_nu_massive, &
        tau_switch_no_nu_multpoles, tau_switch_no_phot_multpoles, tau_switch_nu_nonrel, noSwitch, &
        tau_switch_saha, tau_switch_evolve_TM)                  ! :302-304

    if (next_switch < tauend) then
        if (next_switch > tau+smallTime) then
            call GaugeInterface_ScalEv(EV, y, tau,next_switch,tol1,ind,c,w)
            if (global_error_flag/=0) return
        end if
        EVout=EV
        if (next_switch == EV%TightSwitchoffTime) then
            EVout%TightCoupling=.false.
            EVout%TightSwitchoffTime = noSwitch
            call SetupScalarArrayIndices(EVout)
            call CopyScalarVariableArray(y,yout, EV, EVout)
            EV=EVout
            y=yout
            ind=1                                              ! restart dverk
            !Set up variables with their tight coupling values
            y(EV%g_ix+2) = EV%pig
            ...
        else if (next_switch == tau_switch_nu_massive) then
            ...
            call SwitchToMassiveNuApprox(EV, a, y, nu_i)        ! :378-386
        ...
        end if
        call GaugeInterface_EvolveScal(EV,tau,y,tauend,tol1,ind,c,w)   ! recurse  :427
        return
    end if
    call GaugeInterface_ScalEv(EV,y,tau,tauend,tol1,ind,c,w)

Key ingredients for an axion switch:

- ``noSwitch = State%tau0+1`` and ``smallTime = min(tau, 1/EV%k_buf)/100`` (:257-258).
- Switch times are precomputed conformal times. Precedent for converting an a-threshold (like m/H crossing) to τ:
  ``GaugeInterface_Init`` (:472-519) computes ``nu_tau_nonrelativistic(nu_i) = DeltaTimeMaxed(0._dl, a_nonrel)`` with
  ``a_nonrel = 2.5d0/nu_mass*CP%Accuracy%AccuracyBoost`` (:511-512) using ``State%DeltaTime`` integrals — the axion
  switch time τ_osc (from a_osc where m = switch_factor·H) would be computed the same way, once, k-independent
  (or in the axion Init).
- ``CopyScalarVariableArray(y,yout,EV,EVout)`` (:638-726) maps every block old→new layout; the DE block copy is
  generic over ``num_perturb_equations`` (:649-652). Variable *transformations* at the switch (e.g. AxiECAMB's
  δφ,δφ' → δ_ax,u_ax projection) are done **after** the copy, exactly as tight-coupling seeds
  ``y(EV%g_ix+3) = (3._dl/7._dl)*y(EV%g_ix+2)*(EV%k_buf/opacity)*...`` (:330-343) or Saha seeds
  ``y(EV%xe_ix) = (1-xe)/(2-xe)*(-y(ix_clxb) + (3./2+ CB1/(CP%TCMB/a))*Delta_TM)`` (:404-415).
- ``ind=1`` resets the dverk integrator after a discontinuous change (:322, 389, 398, 406, 418).
- There is no ``EvolveStatus`` enum in 1.6.7 (that was the Nov13 mechanism); the modern equivalent is exactly this
  ``next_switch``/``EVout`` re-layout pattern plus boolean state flags in EvolutionVars.

[DESIGN OPTION] Because ``num_perturb_equations`` is global to the DarkEnergy class (and both phases use 2
equations in AxiECAMB), the *simplest* port keeps 2 DE equations throughout and performs the variable
re-interpretation inside the axion class at τ>τ_switch: ``PerturbationEvolve``/``PerturbedStressEnergy`` branch on
``a >= a_switch`` and re-interpret ``y(w_ix), y(w_ix+1)`` as (δ, (1+w)v) instead of (δφ, δφ'). The one-time variable
transformation at the switch still needs a hook; options: (i) add an axion case to GaugeInterface_EvolveScal's
switch list (clean, follows precedent, requires touching equations.f90); (ii) do the conversion inside the class
on first call past the switch (hacky: dverk history invalid — risk of step rejection; precedent says use (i) and
reset ``ind=1``).

4.5 Where EvolutionVars-driven evolution is launched
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``cmbmain.f90:728 call GetNumEqns(EV)``, ``:824 call GaugeInterface_Init``, ``:980 call initial(EV,y, taustart)``,
``:1146 call initial(EV,y, tau)`` (transfer-only pass). ``dverk`` driver ``GaugeInterface_ScalEv`` equations.f90:205-218
with the famous index-error message (:211-217).

5. Thermo/recombination: how H(z) reaches recfast
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``fortran/recfast.f90:248`` — ``RecombinationData`` holds ``class(CAMBdata), pointer :: State``; ``TRecfast``
  (recfast.f90:251 ff) extends ``TRecombinationModel`` and owns a ``Calc`` (RecombinationData) instance.
- ``TRecfast_init(this, State, WantTSpin)`` recfast.f90:558-628: ``select type(State); class is (CAMBdata):
  Calc%State => State`` (:584-586); pulls ``Calc%Tnow = State%CP%tcmb`` (:602), ``OmegaT=(omch2+ombh2)/H**2`` (:611),
  ``Calc%z_eq = State%z_eq`` (:628).
- **H(z) is obtained generically through dtauda**, recfast.f90:882:

  .. code-block:: fortran

      Hz = ainv**2/dtauda(Recomb%State,1/ainv)/MPC_in_sec

  (also :755 for τ_21). ``dtauda`` is interfaced at recfast.f90:348 (``procedure(obj_function), private :: dtauda``).

[PLUMBING] **Consequence**: any background component included in ``dtauda`` (i.e. in grho_no_de or the DarkEnergy
component) automatically propagates to recombination, thermal history (``Thermo_Init`` uses ``dtauda`` at
results.f90:1844, 1860, 1869), reionization optical depth (``reion_doptdepth_dz`` results.f90:1220), distances,
and derived parameters. AxiECAMB's manual edits to recfast's H(z) are **obsolete** in 1.6.7 — no recfast changes
needed at all (except the ``OmegaT``-style diagnostic variables at :611, which are only used in fudge formulas and
the He Saha — verify whether AxiECAMB altered those; CosmoMC-era OmegaT enters ``Trad``-independent rate terms only).

- The variants cosmorec.f90/hyrec.f90 follow the same TRecombinationModel interface (Makefile_main:7 chooses).

6. Reionization, halofit, transfer columns, initial power, ini driver
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

6.1 Reionization
^^^^^^^^^^^^^^^^^

``fortran/reionization.f90``: ``TBaseTauWithHeReionization`` (:33-60, tanh He second reionization, binary search for
zre from tau via ``SetParamsForZre``), ``TTanhReionization`` (:62-73; the tanh in ``(1+z)**Rionization_zexp``, comment
:12-14), ``TExpReionization`` (:75-87). Selected via ini ``reionization_model``/python ``reionization_model`` class
name. No axion impact expected beyond H(z) (automatic).

6.2 Halofit / nonlinear
^^^^^^^^^^^^^^^^^^^^^^^^

``fortran/halofit.f90``: version constants
``halofit_original=1, halofit_bird=2, halofit_peacock=3, halofit_takahashi=4, halofit_casarini=7, halofit_mead2016=5, halofit_halomodel=6, halofit_mead2015=8, halofit_mead2020=9, halofit_mead2020_feedback=10``,
``halofit_default=halofit_mead2020`` (:58-63). ``type, extends(TNonLinearModel) :: THalofit`` with
``halofit_version`` field (:67-68); ini key read at :258 ``this%halofit_version = Ini%Read_Int('halofit_version', halofit_default)``.
Model selected by name in camb.f90:507-519 (``'HALOFIT'`` → ``allocate(THalofit::P%NonLinearModel)``).
Python: camb/nonlinear.py (``Halofit`` class, ``halofit_version`` named field). [NOTE] AxiECAMB's Nov13 halofit edits
(if any) must be re-evaluated: 1.6.7 halofit gets Omega_m via ``CP%DarkEnergy%Effective_w_wa`` and transfer-based
sigma8; axion-as-matter affects it only through the input P(k) and omegam definitions inside halofit.f90
(grep ``omm0`` there during the port).

6.3 Transfer output columns
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fortran enumeration ``module Transfer`` results.f90:3046-3068:

  .. code-block:: fortran

      integer, parameter :: Transfer_kh =1, Transfer_cdm=2,Transfer_b=3,Transfer_g=4, &
          Transfer_r=5, Transfer_nu = 6,  &
          Transfer_tot=7, Transfer_nonu=8, Transfer_tot_de=9,  &
          Transfer_Weyl = 10, &
          Transfer_Newt_vel_cdm=11, Transfer_Newt_vel_baryon=12, &
          Transfer_vel_baryon_cdm = 13
      integer, parameter :: Transfer_max = Transfer_vel_baryon_cdm
      character(LEN=name_tag_len) :: Transfer_name_tags(Transfer_max-1) = &
          ['CDM     ', 'baryon  ', ...]

- Columns are filled in equations.f90 ``output`` (:2695-2709): ``EV%OutputTransfer(Transfer_cdm) = clxc``, ...,
  ``EV%OutputTransfer(Transfer_tot) = dgrho_matter/grho_matter``,
  ``EV%OutputTransfer(Transfer_tot_de) = dgrho/grho_matter``, ``EV%OutputTransfer(Transfer_Weyl) = k2*phi``, etc.
  Post-scaling: ``Arr(Transfer_kh+1:Transfer_max) = Arr(Transfer_kh+1:Transfer_max)/EV%k2_buf`` (:2136).
- ``transfer_power_var = Transfer_tot`` default for P(k)/σ8 (results.f90:3073-3075).
- [PHYSICS-LANDING] Adding an axion column = (1) new constant ``Transfer_axion=14`` + bump ``Transfer_max`` + extend
  ``Transfer_name_tags`` (results.f90:3051-3068); (2) set ``EV%OutputTransfer(Transfer_axion) = clxax`` in
  equations.f90 output (~:2709); (3) decide whether δ_ax enters ``Transfer_tot``/``Transfer_nonu`` numerators
  (AxiECAMB includes axions in total matter; that means editing :2702-2704 sums and likely ``grho_matter``);
  (4) python mirror: ``camb/model.py:42-55`` (``Transfer_max = ...``) **and** ``transfer_names`` list model.py:84-98
  (the string list order defines ``get_matter_transfer_data`` indexing, results.py:135-137 and var1/var2 lookup
  results.py:802-804); (5) ``MatterTransferData.TransferData(entry,k,z)`` shape flows through automatically
  (classes.f90:8-18, camb_python.f90 ``CAMBdata_MatterTransferData`` :321-339).
- 21cm aliasing trap: ``Transfer_monopole=4, Transfer_vnewt=5, Transfer_Tmat=6`` overlay columns 4-6
  (results.f90:3060) — a new column must come after 13, never reuse 4-6.

6.4 Initial power
^^^^^^^^^^^^^^^^^^

``fortran/InitialPower.f90:43-65`` ``TInitialPowerLaw`` (ns, nrun, nrunrun, nt, ntrun, r, pivot_scalar, pivot_tensor,
As, At); ``TSplinedInitialPower`` (:67-83) for tabulated P(k). ``ScalarPower(k)`` is the only thing the Boltzmann
code calls (via ``CAMBparams_PrimordialPower`` model.f90:429-452 and cmbmain). Python: camb/initialpower.py.
[PHYSICS-LANDING] AxiECAMB's isocurvature amplitude (alpha_ax·(Hinf/...)² style) would either be a new
``TInitialPower`` subclass or handled by scaling InitVec — note modern CAMB computes each IC mode's C_l with the
*same* ScalarPower; mode-dependent amplitude needs either InitialConditionVector weighting (only for
initial_vector mixed mode) or a custom power class checking ``Scalar_initial_condition``.

6.5 Ini driver / standalone fortran
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``fortran/inidriver.f90`` is now 18 lines — just ``call CAMB_CommandLineRun(InputFile)`` (inidriver.f90:15).
- All ini reading is in ``fortran/camb.f90``: ``CAMB_ReadParams(P, Ini, ErrMsg)`` camb.f90:238 ff. Dark energy model
  selection by string (camb.f90:422-436):

  .. code-block:: fortran

      DarkEneryModel = UpperCase(Ini%Read_String_Default('dark_energy_model', 'fluid'))
      if (allocated(P%DarkEnergy)) deallocate(P%DarkEnergy)
      if (DarkEneryModel == 'FLUID') then
          allocate (TDarkEnergyFluid::P%DarkEnergy)
      else if (DarkEneryModel == 'PPF') then
          allocate (TDarkEnergyPPF::P%DarkEnergy)
      else if (DarkEneryModel == 'AXIONEFFECTIVEFLUID') then
          allocate (TAxionEffectiveFluid::P%DarkEnergy)
      else if (DarkEneryModel == 'EARLYQUINTESSENCE') then
          allocate (TEarlyQuintessence::P%DarkEnergy)
      else ...
      call P%DarkEnergy%ReadParams(Ini)

  [PLUMBING] A new ``TAxiEDE``/axion class adds one ``else if`` branch here and implements
  ``ReadParams(Ini)`` (pattern: TAxionEffectiveFluid_ReadParams DarkEnergyFluid.f90:150-164 reads
  ``Ini%Read_Double('AxionEffectiveFluid_w_n')`` etc.). New top-level CAMBparams ini keys go directly in
  CAMB_ReadParams (e.g. near ``P%omnuh2 = Ini%Read_Double('omnuh2')`` camb.f90:443).
- ``initial_condition`` ini key: camb.f90:578 ``P%Scalar_initial_condition = Ini%Read_Int('initial_condition', initial_adiabatic)``;
  ``initial_vector`` string :581.
- Ini machinery (``TIniFile``, ``Ini%Read_*``) comes from **forutils** (IniObjects), not the repo-local ``inifiles/``
  directory (that only holds sample .ini files like ``inifiles/params.ini``).
- Python can also write/read ini: ``camb/_ini.py`` + ``CAMBparams.write_ini`` model.py:422-429, and
  ``camb.read_ini`` uses ``CAMB_ReadParamFile`` (camb.f90:232).

7. Build system
~~~~~~~~~~~~~~~~

7.1 Fortran
^^^^^^^^^^^^

- ``fortran/Makefile`` (146 lines): compiler detect (ifort vs gfortran, Makefile:32-127), flags
  (gfortran: ``-MMD -cpp -ffree-line-length-none -fmax-errors=4 -fopenmp``, ``FFLAGS=-O3``, Makefile:94-96),
  forutils path autodetect ``../forutils`` (Makefile:9-23), then ``include ./Makefile_main`` (:146).
- ``fortran/Makefile_main``:

  - **Source list** (Makefile_main:31-35):

    .. code-block:: make

        SOURCEFILES = constants config classes MathUtils subroutines DarkAge21cm \
                DarkEnergyInterface SourceWindows massive_neutrinos model results bessels \
                $(RECOMBINATION_FILES) $(DARKENERGY_FILES) equations \
                $(REIONIZATION_FILES) $(POWERSPECTRUM_FILES) $(NONLINEAR_FILES) \
                lensing $(BISPECTRUM) cmbmain camb camb_python

    with ``DARKENERGY_FILES ?= DarkEnergyFluid DarkEnergyPPF PowellMinimize DarkEnergyQuintessence``
    (Makefile_main:9). [PLUMBING] **Adding a new source file** (e.g. ``AxionEDE.f90``): add its stem to
    ``DARKENERGY_FILES`` (or SOURCEFILES) *in dependency-safe position* — order matters only for the first build
    since ``-MMD``/``-gen-dep`` dependency files (``-include *.d``, :171) handle rebuilds; the new module must ``use``
    only modules earlier in the compile graph (DarkEnergyInterface, results, etc., all of which precede the
    DARKENERGY_FILES slot which is expanded after ``model results`` and before ``equations``).
  - Targets: ``camb`` (static lib + executable, :94-108), ``python`` →
    ``make -C $(DLL_DIR) camblib.so F90FLAGS="$(SF90FLAGS)"`` (:100-101), ``camblib.so`` link rule (:117-120) which
    ends with ``cp $(DLL_DIR)/camblib.so $(PYCAMB_OUTPUT_DIR)`` (= ``../camb``, :25).
  - forutils is built first (``libforutils``/``libforutils_so`` :140-146) and linked via
    ``LIBLINK = -L"$(FORUTILS_DIR)" -lforutils`` (:55). Note comment :116: "cannot link the .a library, or the
    linker will strip out some things we might need" — the .so links the .o files directly.

7.2 Python build & library loading
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``setup.py:95-198 make_library()``: on unix runs
  ``make python PYCAMB_OUTPUT_DIR=../camb/ CLUSTER_SAFE=...`` (setup.py:184-188) from ``fortran/``, then verifies
  ``camb/camblib.so`` exists (:191-192); on Windows manual gfortran compile of ``SOURCEFILES`` parsed out of
  Makefile_main (:111-173). ``python setup.py make`` (MakeLibrary command, :201-212) rebuilds in place.
- ``camb/baseconfig.py:41-65``: ``DLLNAME = 'camblib.so'`` (``cambdll.dll`` on Windows), loaded from the package dir:
  ``CAMBL = osp.join(BASEDIR, DLLNAME)``; loads via ``ctypes.LibraryLoader(IfortGfortranLoader)`` (:65); fails with
  "Library file ... does not exist" otherwise. So **for the python wrapper to pick up modifications, you only
  need** ``camblib.so`` **regenerated into** ``camb/`` (``make python`` does the copy automatically). Mangling assumption:
  ``__<module>_MOD_<proc>`` (gfortran) or ``<module>_mp_<proc>_`` (ifort) — handled by IfortGfortranLoader
  (baseconfig.py:68 checks ``handles_mp_set_cls_template_``).
- The struct-size handshake at import (``_get_fortran_sizes`` baseconfig.py:138-151) plus per-class
  ``f_<method>`` imports mean a **field-order mismatch between model.f90 and model.py fails silently** — there is
  no checksum on CAMBparams layout. The test suite (``camb/tests/camb_test.py``) is the guard; run it after layout
  changes.

8. Python API surface for new axion params
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``camb.set_params(**params)`` camb/camb.py:114-203 dispatches kwargs by introspection over this fixed setter
  chain (camb.py:183-191, order matters — comment :181 "must call DarkEnergy.set_params before set_cosmology if
  setting theta rather than H0"):

  .. code-block:: python

      do_set(cp.set_accuracy)
      do_set(cp.set_classes)
      do_set(cp.DarkEnergy.set_params)
      do_set(cp.Reion.set_extra_params)
      do_set(cp.set_cosmology)
      do_set(cp.set_matter_power)
      do_set(cp.set_for_lmax)
      do_set(cp.InitPower.set_params)
      do_set(cp.NonLinearModel.set_params)

  Unknown kwargs raise ``CAMBUnknownArgumentError`` (:196-205). [PLUMBING] If axion params are fields of the
  axion DarkEnergy class with a ``set_params`` method, ``set_params(dark_energy_model='AxionEDE', m_ax=..., ...)``
  works with zero changes to camb.py (kwargs reach ``cp.DarkEnergy.set_params`` after ``set_classes`` swaps the
  class). If they are CAMBparams fields, users set them directly or you extend ``set_cosmology``
  (model.py:614-759 — signature and body quoted in §1; you would add ``m_ax=None, omaxh2=None,...`` kwargs and the
  budget logic there).
- Model selection by name: ``CAMBparams.set_classes(dark_energy_model=...)`` model.py:799-825 →
  ``self.DarkEnergy = self.make_class_named(dark_energy_model, DarkEnergyModel)``; ``make_class_named``
  (baseconfig.py:805-817) resolves through ``F2003Class._class_names`` which is populated by the ``@fortran_class``
  decorator (baseconfig.py:840) plus the alias dict at dark_energy.py:239:
  ``F2003Class._class_names.update({"fluid": DarkEnergyFluid, "ppf": DarkEnergyPPF})`` — add an
  ``"axion": AxionEDE``-style alias the same way.
- **EarlyQuintessence python class** (the template to copy), camb/dark_energy.py:

  - Base ``DarkEnergyModel`` :6-14 mirrors only ``(__is_cosmological_constant, c_bool), (__num_perturb_equations, c_int)`` — every DE subclass inherits these two hidden fields **first**.
  - ``Quintessence`` :158-189 mirrors TQuintessence fields **in declaration order including privates**:
    ``DebugLevel, astart, integrate_tol, sampled_a, phi_a, phidot_a (AllocatableArrayDouble), __npoints_linear,
    __npoints_log, __dloga, __da, __log_astart, __max_a_log, __ddphi_a, __ddphidot_a, __state (f_pointer)``.
    Note the ``class(CAMBdata), pointer :: State`` is mirrored as raw ``("__state", f_pointer)`` — pointers must be
    padded this way. ``__getstate__`` raises ("Cannot save class with splines") :188-189.
  - ``EarlyQuintessence`` :192-235: ``_fields_`` = n, f, m, theta_i, frac_lambda0, use_zc, zc, fde_zc, npoints,
    min_steps_per_osc, fde (AllocatableArrayDouble), ``__ddfde``; ``_fortran_class_name_ = "TEarlyQuintessence"``,
    ``_fortran_class_module_ = "Quintessence"`` (inherited from Quintessence :186); plus a plain-python
    ``set_params(self, n, f=0.05, m=5e-54, theta_i=0.0, use_zc=True, zc=None, fde_zc=None)``.
  - ``AxionEffectiveFluid`` :131-154 — same pattern, 4 fields, ``_fortran_class_module_ = "DarkEnergyFluid"``.
- ``CAMBdata`` python class: camb/results.py:169 ff, ``_fortran_class_module_ = "results"`` (:181), ``_methods_``
  (:236 ff) binding ``CAMBdata_*`` functions exported in camb_python.f90 (e.g. ``CAMBdata_GetBackgroundDensities``
  via ``get_background_densities``, ``CAMBdata_MatterTransferData`` camb_python.f90:321-339).

9. Risk register / surprises for the port designer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **[RISK] CAMBparams ctypes mirroring is unchecked.** Any new Fortran field must be inserted in model.py
   ``_fields_`` at the same position; allocatable/class members need the right wrapper type. Silent corruption
   otherwise. (baseconfig.py:459-538.)
2. **[RISK]** ``is_cosmological_constant`` **short-circuits.** dtauda (equations.f90:12 always calls the DE class),
   but derivs uses ``EV%is_cosmological_constant`` (set from the DE class at initial(), equations.f90:1780) to skip
   ``PerturbedStressEnergy``/``PerturbationEvolve`` (:2284, 2304) — an axion-in-DE class must set
   ``is_cosmological_constant = .false.`` in Init even when w≈−1 today.
3. **[RISK] BackgroundDensityAndPressure guard** ``if (a > 1e-10)`` **returns grhov_t=0 below a=1e-10**
   (DarkEnergyInterface.f90:95-99); quintessence overrides with its own astart=1e-7 floor
   (DarkEnergyQuintessence.f90:144-155 returns grhov_t=0, w=−1 for a<astart). The axion's early-time
   (pre-table) limit must be handled explicitly (AxiECAMB starts KG deep in radiation domination).
4. **[SURPRISE/GOOD] No edits needed in recfast/reionization/thermo for the background** — all H(a) flows through
   dtauda (§5). Most of AxiECAMB's scattered Nov13 background edits collapse into one class.
5. **[RISK] "matter" definitions are scattered**: ``grho_matter``/``dgrho_matter`` (equations.f90:2216, 2227),
   ``Transfer_tot``/``Transfer_nonu`` (equations.f90:2702-2704), z_eq (results.f90:476), CosmomcTheta omdmh2
   (results.f90:924), ``om`` superhorizon expansions (equations.f90:1830, results.f90:878), halofit's Omega_m.
   AxiECAMB treats the axion as matter post-switch; each of these needs an explicit decision.
6. **[RISK] Tensor/vector evolution uses grhov_t via BackgroundDensityAndPressure but no DE perturbations** —
   matches AxiECAMB (axion has no tensor sources), nothing to do, but verify AxiECAMB didn't add tensor terms.
7. **[PRECEDENT GAP] No existing mechanism changes** ``num_perturb_equations`` **mid-run.** The KG→fluid switch must
   either keep equation count fixed (recommended; both phases are 2-variable) with an in-class variable
   re-interpretation plus a one-time transformation hooked into GaugeInterface_EvolveScal's switch chain
   (equations.f90:302-433), or add a new EV flag + SetupScalarArrayIndices branch. Either way ``ind=1`` reset is
   mandatory after the transform.
8. **[ACCURACY]** AccuracyParams (model.f90:35-94) replaces Nov13's hardcoded boosts; AxiECAMB's hand-tuned
   sampling tweaks (e.g. extra timesteps around the switch) should map onto ``TimeStepBoost``,
   ``BackgroundTimeStepBoost``, ``IntTolBoost`` or onto class-local parameters like ``npoints``/``min_steps_per_osc``
   (TEarlyQuintessence pattern) rather than global constants.
9. **[PLUMBING] PowellMinimize.f90 (NEWUOA/BOBYQA) and brentq are already in the build** (Makefile_main:9) for
   shooting; no need to port AxiECAMB's bespoke root-finder.
10. **[NOTE] initv is dimensioned** ``initv(6,1:i_max)`` (equations.f90:1774) — one free isocurvature row already
    exists, suggesting upstream anticipated a 6th mode; still update ``initial_nummodes`` and the python names list.

10. Quick landing-zone table (AxiECAMB change → modern location)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: auto

   * - AxiECAMB concern
     - Modern landing zone
   * - m_ax, omaxh2/axfrac params
     - New ``TAxionEDE`` class fields (DarkEnergy slot) or CAMBparams (model.f90:135 ff + model.py:284 ff)
   * - KG background solve + tables
     - New class extending TQuintessence pattern: Init hook results.f90:497, solver pattern DarkEnergyQuintessence.f90:301-551
   * - omaxh2 shooting
     - NEWUOA/brentq pattern DarkEnergyQuintessence.f90:327-381, 600-636
   * - H(a) with axion
     - automatic via dtauda equations.f90:4-23 once in DE component; else grho_no_de results.f90:1231
   * - Budget closure (omdah2 split)
     - results.f90:462-464 Omega_de calculation
   * - KG perturbations (pre-switch)
     - PerturbedStressEnergy/PerturbationEvolve override (quintessence form DarkEnergyQuintessence.f90:237-272)
   * - Effective fluid + cs²(k,a) (post-switch)
     - TAxionEffectiveFluid form DarkEnergyFluid.f90:249-291
   * - Switch mechanism
     - GaugeInterface_EvolveScal next_switch chain equations.f90:245-433; τ from a via DeltaTimeMaxed pattern :460-519
   * - Axion isocurvature mode 6
     - equations.f90:62-64 + initv row 6 (:1845-1909) + PerturbationInitial + model.py:350-358 + camb.f90:578
   * - Hinf/alpha_ax amplitude
     - InitialPower subclass or InitVec scaling (see §6.4)
   * - Transfer axion column
     - results.f90:3051-3068 + equations.f90:2695-2709 + model.py:42-98
   * - Derived params (z_switch etc.)
     - nthermo_derived results.f90:43-46 + :2172-2199 + model.py:41,68-82
   * - ini keys
     - camb.f90:238-600 (CAMB_ReadParams) + class ReadParams
   * - python API
     - dark_energy.py new class + ``_class_names`` alias :239; set_params chain camb.py:183-191
   * - build
     - Makefile_main:9 DARKENERGY_FILES + ``make python`` → camb/camblib.so (setup.py:184-192)


Original-code analysis: the axion background solver (axion_background.F90)
---------------------------------------------------------------------------

**File:** ``/Users/vivianmiranda/data/research/WayneHu/rayne/AxiECAMB/axion_background.F90`` (1429 lines)

**Status:** Entirely NEW file (no OLDCAMB counterpart). Implements the ultralight-axion
background solver of arXiv:1410.2896 (Hlozek et al.) with the effective-fluid-approximation
(EFA) upgrades of arXiv:2412.15192 (AxiECAMB, "RL" = Rayne Liu edits).

1. Module and public interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Module name: ``axion_background`` (line 31). There are **no** ``private``/``public``
**statements**, so all six contained procedures are public. Only ``w_evolve`` and
``get_phase_info`` are called from outside the module (see §5); the other four are
internal helpers but technically exported.

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Procedure
     - Lines
     - Signature
     - Purpose
   * - ``w_evolve(Params, badflag)``
     - 35–1172
     - ``type(CAMBparams) :: Params`` (inout), ``integer badflag`` (out; set to 1 on
       collapse/NaN)
     - Master driver: solves KG background, shoots for ``phi_init`` to match
       ``omegah2_ax``, finds ``a_osc``, performs EFA matching, fills global
       interpolation tables and ~20 ``Params%`` fields.
   * - ``derivs_bg(a, v, dvt_dloga, omegah2_regm, omegah2_rad, omegah2_lambda, omk, hsq, maxion_twiddle, badflag, lhsqcont_massless, lhsqcont_massive, Nu_mass_eigenstates, Nu_masses)``
     - 1180–1203
     - ``v(1:2)`` in, ``dvt_dloga(1:2)`` out, rest real(dl)/int in
     - RHS of the KG system in d/d(ln a).
   * - ``get_phase_info(Params, y_beta, beta_coeff, movHETA, beta2x)``
     - 1207–1218
     - all real(dl) out except Params in
     - Analytic WKB oscillation-phase estimate at the switch; used by
       ``inidriver_axion`` to tune ``dfac``.
   * - ``auxiIC(Params, omegah2_regm, omegah2_rad, omegah2_lambda, omk, hnot, maxion_twiddle, a, v, badflag, lhsqcont_massless, lhsqcont_massive, Nu_mass_eigenstates, Numasses, littlehauxi, lhETA, A_coeff, tvarphi_c, tvarphi_cp, tvarphi_s, tvarphi_sp, rhorefp, Prefp)``
     - 1221–1327
     - ``v(1:2)`` = (v1,v2) at the switch; outputs: instantaneous conformal Hubble
       ``littlehauxi``, time-averaged ``lhETA``, ``A_coeff``, the four EFA field
       projections ``tvarphi_*``, EFA density ``rhorefp`` and pressure ``Prefp``
     - EFA matching at ``a_osc``: iteratively determines ⟨H⟩ and ``wEFA_c``, projects
       (φ, φ′) onto cos/sin WKB basis, returns cycle-averaged ρ and P.
       **Mutates** ``Params%wEFA_c``.
   * - ``lh(omegah2_regm, omegah2_rad, omegah2_lambda, omk, hsq, maxion_twiddle, a, v, littlehfunc, badflag, lhsqcont_massless, lhsqcont_massive, Nu_mass_eigenstates, nu_masses [, rho_f])``
     - 1334–1382
     - ``littlehfunc`` out; ``rho_f`` ``optional`` — if present, replaces field energy
       with given EFA density
     - Computes **conformal** ``aH/(100 km/s/Mpc)`` (despite "littleh" name; conversion
       at line 1369). Sets ``badflag=1`` if H²≤0 or NaN.
   * - ``next_step(a, v, kvec, kfinal, avec, omegah2_regm, omegah2_rad, omegah2_lambda, omk, hsq, maxion_twiddle, badflag, dloga, nstep, cmat, lhsqcont_massless, lhsqcont_massive, Nu_mass_eigenstates, Nu_masses)``
     - 1387–1424
     - ``kvec(1:2,1:nstep)``, ``kfinal(1:2)`` out; ``nstep=16``
     - One 16-stage RK evaluation: fills the stage derivatives ``kvec`` (each
       pre-multiplied by ``dloga``).

Modules used: ``ModelParams``, ``constants``, ``Precision``, ``MassiveNu`` (w_evolve,
auxiIC); helpers use ``constants``, ``Precision``, ``MassiveNu`` only. ``spline`` and
``spline_out`` come from ``subroutines.f90`` (``subroutines.f90:288`` and ``:12``).

2. Variable definitions, units, and normalization conventions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

From header comments (lines 12–25) and code:

- KG equation in conformal time: ``\ddot{φ} + 2H φ̇ + m² a² φ = 0`` (H = conformal
  Hubble; dots = d/dτ).
- **v1** (``v_vec(1,:)``, stored as ``phinorm_table``): ``phi = sqrt(3/(4πG)) v1`` ⇒
  v1 = φ·sqrt(4πG/3) = φ/(√6 M_pl) with M_pl = (8πG)^(-1/2). Hence
  ``Params%phiinit = vtwiddle_init*sqrt(6)`` (line 1151) is φ_init in reduced-Planck
  units.
- **v2** (``v_vec(2,:)``, stored as ``phidotnorm_table``): dφ/dτ = u2 with
  ``u2 = H0·u2~``, ``u2~ = sqrt(3/(4πG)) v2`` ⇒ v2 = (dφ/dτ)·sqrt(4πG/3)/H0.
  (equations_ppf.f90:3383 confirms: "dv1 = sqrt(4πG/3) dphi, dv2 = sqrt(4πG/3)/H0
  dphidot".)
- **maxion_twiddle** = ``Params%m_ovH0`` = m/H0 (dimensionless; driver sets
  ``P%m_ovH0 = P%ma/P%H0_eV``, inidriver_axion.F90:279, with ``ma`` in eV).
- **hnot** = H0/100 = ``Params%H0/100.d0``; ``hsq = hnot**2``.
- **omegah2_X** ≡ Ω_X h² for dm, b, lambda, ax, nu (lines 252–263).
- ``lh`` **output** (``littlehfunc``): conformal aH in units of 100 km/s/Mpc. So
  ``lh/(a·hnot) = H/H0`` and ``H/m = lh/(a·hnot·maxion_twiddle)``.
- **Axion energy density** (line 1024):

  .. code-block:: fortran

      rhoaxh2_ov_rhom(i)=(v_vec(2,i)/a_arr(i))**2.0d0+(maxion_twiddle*v_vec(1,i))**2.0d0

  This is the *instantaneous* Ω_ax(a)h² ≡ ρ_ax(a)·h²/ρ_crit,0 (NOT scaled to today).
  Comment lines 1032–33: "This is the axion energy density \*h²/(3H0²/8πG) … in camb's
  definitions grhoa2_axion = grhom*grhox_table_internal(i)". I.e. to get CAMB Nov13
  ``grho`` units (Mpc⁻², factor 8πG ρ a⁴-style), multiply by
  ``grhom = 3*(hsq*1.d10)/c**2`` (line 294 — this routine itself (re)sets the
  ModelParams global ``grhom``!).
- **rhorefp / Prefp**: EFA cycle-averaged density/pressure at the switch, in the same
  Ω(a)h² units; ``Params%rhorefp_ovh2 = rhorefp/hsq`` (line 1091) is Ω_ax(a_osc)
  (no h²); ``Params%Prefp = Prefp`` is left in Ω(a)h² units (line 1092, "I don't want
  to rescale with hsq anymore"). Downstream code (equations_ppf.f90:1097) divides
  ``CP%Prefp`` by ``CP%rhorefp_ovh2*(CP%H0**2/1.0d4)`` — consistent.
- **radiation**: lines 295–308 recompute ``grhog``, ``grhor`` (globals!) and build

  .. code-block:: fortran

      rhocrit=(8.0d0*const_pi*G*1.d3/(3.0d0*((1.d7/(MPC_in_sec*c*1.d2))**(2.0d0))))**(-1.0d0)
      Params%omegah2_rad=((Params%TCMB**4.0d0)/(rhocrit))/(c**2.0d0)
      Params%omegah2_rad=Params%omegah2_rad*a_rad*1.d1/(1.d4)
      lhsqcont_massless=(Params%Nu_massless_degeneracy*grhor*(c**2.0d0)/((1.d5**2.0d0)))/3.0d0
      Params%omegah2_rad=Params%omegah2_rad+lhsqcont_massless

  so ``Params%omegah2_rad`` = Ω_γ h² + Ω_ν,massless h² (massless-ν part uses
  ``Nu_massless_degeneracy`` — DG fix; original used ``Num_Nu_massless``).
- **massive neutrinos** (lines 310–325):
  ``nu_constant=(7/120)π⁴/(ζ₃·1.5)·omnuh2·(grhom/grhor)/hsq``; sets **global**
  ``Nu_masses(k)=nu_constant*Params%Nu_mass_fractions(k)/Params%Nu_mass_degeneracies(k)``
  (MassiveNu module variable) and
  ``lhsqcont_massive(k)=Params%Nu_mass_degeneracies(k)*(grhor*c**2/(1.d5**2))/3.0d0``
  (Ω h² per unit ``Nu_rho``). In ``lh``, massive-ν energy enters as
  ``sum(lhsqcont_massive*mass_correctors)/a**4`` with ``mass_correctors(i)`` from
  ``call Nu_rho(a*nu_masses(i),rhonu)``.

**[PLUMBING/OBSOLETE]** The recomputation of ``grhom/grhog/grhor``, ``rhocrit``,
``Nu_masses`` duplicates CAMB's own ``ModelParams``/``MassiveNu`` initialization
("non-trivial aspects to sharing neutrino data structures… we recompute some things
twice", lines 62–65). In modern CAMB 1.6.7 these all live in
``CAMBdata``/``TNuPerturbations`` and must NOT be re-set; read them from the state
object instead. The physics requirement is only: the background H(a) used inside this
solver must be *identical* to the one CAMB uses elsewhere.

3. Background algorithm, step by step (``w_evolve``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

3.0 Re-entrancy [PLUMBING]
^^^^^^^^^^^^^^^^^^^^^^^^^^

Lines 237–243: deallocate the 7 global tables (``loga_table``, ``phinorm_table``,
``phidotnorm_table``, ``phinorm_table_ddlga``, ``phidotnorm_table_ddlga``,
``rhoaxh2ovrhom_logtable``, ``rhoaxh2ovrhom_logtable_buff``) if allocated —
``w_evolve`` is called up to dozens of times per cosmology by the driver (§6).

3.1 Switch constant and EFA w-coefficient [PHYSICS]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Line 333: ``dfac = Params%dfac`` — the switch is at **m = dfac·H** (conformal-time
  H/a, i.e. physical H). Driver default ``P%dfac = 10._dl`` (inidriver_axion.F90:154,
  "making movH internal" — no longer read from ini). dfac is then *retuned* by the
  driver (§6).
- Line 335: ``Params%wEFA_c = 9._dl/8._dl`` — initial value of the coefficient in the
  EFA equation of state w_ax = wEFA_c·(H/m)²; later iterated to self-consistency in
  ``auxiIC``. (9/8 is the analytic RD value.)

3.2 Integrator setup [ACCURACY]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Lines 345–525: hard-coded Butcher tableau (``avec(1:16)``, ``cmat(16,16)``,
``svec(1:15)``) of **Fehlberg's classical 8th-order RK** (16 stages; NASA TR, page-75
formulas; URL in comment line 349). Used as a **fixed-step** method on a uniform ln a
grid — ``kfinal`` (the embedded lower-order estimate) is computed but never used for
step control. Comment lines 701–703: "a reasonably accurate integrator is required to
accurately obtain the adiabatic sound speed at earlier times… this integrator + 5000
grid points was necessary to avoid exciting a non-physically large low-l ISW effect."
(Note: ``ntable`` is no longer 5000; see §6 — it is ``nint(dfac*100)+1``, i.e. ~1001 by
default, because the table now ends at ~1.1·a_osc rather than a=1.)

- ``svec(15)`` third digit was fixed by DG (comment line 525) — value
  ``0.3072649547580d-1``.
- Update formula (lines 711–713):
  ``v_vec(:,i)=v_vec(:,i-1)+(svec(1)*kvec(:,1)+svec(9)*kvec(:,9)+svec(10)*kvec(:,10)+svec(11)*kvec(:,11)+svec(12)*kvec(:,12)+svec(13)*kvec(:,13)+svec(14)*kvec(:,14)+svec(15)*kvec(:,15))``
  — only stages 1, 9–15 carry weight.
- ``next_step`` builds the stages: ``vfeed(cp) = dot_product(cmat(m,1:m), kvec(cp,1:m))``
  (line 1407, "RL 011025 reverse kvec order" fix), evaluation points
  ``a*dexp(dloga*avec(m))``.

3.3 Grid in ln a [PHYSICS/ACCURACY]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Lines 534–576:

- ``a_init = min(a_rel, a_lambda, a_m, as_matt, as_rad, as_scalar)*1.d-8`` where

  .. code-block:: fortran

      as_matt=(omegah2_regm/(maxion_twiddle**2.0d0))**(1.0d0/3.0d0)
      as_rad=(Params%omegah2_rad/(maxion_twiddle**2.0d0))**(1.0d0/4.0d0)
      a_m=(Params%omegah2_rad/(omegah2_regm))
      a_rel=10.0d0
      a_lambda=(Params%omegah2_rad/omegah2_lambda)**(0.25d0)
      as_scalar=(omegah2_ax/(maxion_twiddle**2.0d0))**(1.0d0/3.0d0)

  (start 8 decades before the axion matters relative to anything).
- ``a_final = 1.0d0`` initially; uniform
  ``dloga=(log_a_final-log_a_init)/(dble(ntable-1))``; ``loga_table(i)`` filled with
  **natural log** of a (converted to log10 later, §3.8!),
  ``a_arr(i)=dexp(loga_table(i))``.
- ``omk=1.0d0-(omegah2_m+Params%omegah2_rad+omegah2_lambda+omnuh2)/hsq`` (line 536).

3.4 KG initial conditions [PHYSICS]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For each trial ``vtwiddle_init`` (lines 669–676):

.. code-block:: fortran

    v_vec(1,1)=vtwiddle_init
    v_vec(2,1)=0.0d0
    call lh(... a_arr(1), v_vec(1:2,1), littlehfunc(1), ...)
    v_vec(2,1)= - vtwiddle_init * (a_arr(1)**2.0d0) * (maxion_twiddle**2.0d0) * hnot/(5.0d0 * littlehfunc(1))

i.e. field starts at rest on the hill, then v2 is set to the **early-time attractor**
value ``v2 = −v1 a² m̃² H0/(5 aH)`` (the 1/5 factor is the RD slow-roll solution
φ′ = −m²a²φ·τ/5; works because lh = aH and in RD τ = 1/(aH)).

3.5 KG evolution equations [PHYSICS]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``derivs_bg`` (lines 1200–1202):

.. code-block:: fortran

    dvt_da(1)=v(2)*dsqrt(hsq)/(a*(lhr))
    dvt_da(2)=-2.0d0*v(2)/(a)-(maxion_twiddle**2.0d0)*a*v(1)*dsqrt(hsq)/(lhr)
    dvt_dloga(1:2)=a*dvt_da(1:2)

(``dsqrt(hsq)=hnot`` converts ``lhr`` [100 km/s/Mpc units] to aH/H0; "solution proposed
by astralsight5" comment marks a units bugfix vs original axionCAMB.)

Friedmann (``lh``, lines 1357–1369):

.. code-block:: fortran

    littlehfunc=(omegah2_regm/(a**3.0d0)+omegah2_rad/(a**4.0d0))+&
         &sum(lhsqcont_massive*mass_correctors,Nu_mass_eigenstates)/(a**4.0d0)
    littlehfunc=littlehfunc+omegah2_lambda
    ! field (or, if present(rho_f), EFA density):
    littlehfunc=littlehfunc+(maxion_twiddle*v(1))**2.0d0+((v(2)/a)**2.0d0)
    littlehfunc=littlehfunc*(a**2.0d0)+omk*hsq      ! -> conformal (aH)^2
    littlehfunc=dsqrt(littlehfunc)

``badflag=1`` if ``littlehfunc<=0`` (collapsing universe) or NaN (checked via
``.not. littlehfunc == littlehfunc``, line 1375 — ``isnan`` intrinsic deliberately
avoided, RL 010925).

3.6 Switch detection (a_osc) [PHYSICS]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- During integration: ``diagnostic(i)=dfac*littlehfunc(i)/(a_arr(i)*hnot)`` =
  dfac·H(a)/H0 (lines 688/731/1018; "RL fixed by adding hnot"). First-crossing grid
  guess (lines 734–741): when ``maxion_twiddle < diagnostic(i-1)`` and
  ``maxion_twiddle >= diagnostic(i)`` → ``aosc_guess(j)=a_arr(i)``. Sentinel
  ``15.0d0`` = "never oscillates / not yet found".
- Refinement (lines 750–776): ``f_arr = dlog(maxion_twiddle/diagnostic)`` =
  ln(m/(dfac·H)); **natural cubic spline** (``d1 = d2 = 1.0d50`` triggers natural BCs)
  of ln a vs f_arr inverted at 0:

  .. code-block:: fortran

      call spline(f_arr(1:ntable),(loga_table(1:ntable)),ntable,d1,d2,abuff(1:ntable))
      call spline_out(f_arr(1:ntable),loga_table(1:ntable),abuff(1:ntable),ntable,0.0d0,laosc)

  ``aosc = exp(laosc)``. **Exact criterion: m = dfac·H(a_osc), dfac = Params%dfac**
  (driver-controlled; 10 default, up to ~23+ after phase tuning, or ``dfac_skip`` for
  the recombination skip).
- If m̃ ≥ 10 but no crossing by a_final: ``aosc_guess(j) = 1.0_dl - 1.e-3_dl``
  (lines 798–802).
- Field values at switch (lines 804–830): spline v1 and v2 vs ln a with **analytic
  endpoint derivatives** (exact KG derivatives, not finite differences):

  .. code-block:: fortran

      d1 = v_vec(2,1)*hnot/littlehfunc(1)                       ! dv1/dlna
      d2 = v_vec(2,ntable)*hnot/littlehfunc(ntable)
      ...
      d1 = -(2.0d0*v_vec(2,1) + ((maxion_twiddle*dexp(loga_table(1)))**2.0d0)*v_vec(1,1)*hnot/littlehfunc(1))   ! dv2/dlna

  evaluate ``phiosc``, ``phidosc`` at ``laosc``.

3.7 EFA matching at the switch (``auxiIC``) [PHYSICS — core new physics vs axionCAMB]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Called per bisection trial (line 840) and once more with the converged solution
(line 1068). Inputs (φ, φ′) at a_osc; outputs everything the perturbation module needs
at the switch.

Algorithm (lines 1242–1326):

1. ``littlehauxi`` = instantaneous conformal aH at a_osc (full Friedmann with field
   energy).
2. Initialize ``lhETA = littlehauxi``,
   ``rhorefp = (v(2)/a)**2 + (maxion_twiddle*v(1))**2`` (instantaneous),
   ``tol_EFA = 1.e-7_dl``, ``maxiter = 30``,
   ``Hovm_ins = littlehauxi/(a*hnot*maxion_twiddle)``.
3. Massive-ν w: ``call Nu_background(a*numasses(i),rhonu,pnu)``; ``w_nu(i)=pnu/rhonu``.
4. Iterate (``iter_EFA = 1, maxiter``):

   .. code-block:: fortran

       Hovm_ETA = lhETA/(a*hnot*maxion_twiddle)
       w_ax = (Hovm_ETA**2.0d0)*Params%wEFA_c
       dHsqdmt_term = omegah2_regm*(-3.0d0)/a+&                          ! matter
            &omegah2_rad*(-4.0d0)/(a2)+&                                  ! radiation
            &sum(lhsqcont_massive*mass_correctors*(-3.0d0*(1.0d0 + w_nu)),Nu_mass_eigenstates)/(a2)+& ! massive nu
            &rhorefp*(-3.0d0*(w_ax + 1.0d0))*(a2)+&                       ! axion EFA
            &omk*(hnot**2.0d0)*(-2.0d0)                                   ! curvature
       dHsqdmt_term = dHsqdmt_term/(lhETA**2.0d0)                         ! dimensionless dH^2/d(mt)
       A_coeff = (-Hovm_ETA/2.0d0)*(3.0d0 - dHsqdmt_term)
       A_denom = A_coeff**2.0d0 + 3.0d0*A_coeff*Hovm_ins + 4.0d0
       tvarphi_c = v(1)
       tvarphi_cp = -3.0d0*Hovm_ins*(2.0d0*v(1) + (A_coeff + 3.0d0*Hovm_ins)*v(2)/(a*maxion_twiddle))/A_denom
       tvarphi_s = v(2)/(a*maxion_twiddle) - tvarphi_cp
       tvarphi_sp = 3.0d0*Hovm_ins*(A_coeff*v(1) - 2.0d0*v(2)/(a*maxion_twiddle))/A_denom
       rhorefp = (maxion_twiddle**2.0d0)*(tvarphi_c**2.0d0 + tvarphi_s**2.0d0 + &
            &(tvarphi_cp**2.0d0 + tvarphi_sp**2.0d0)/2.0d0 - tvarphi_c*tvarphi_sp + tvarphi_s*tvarphi_cp)
       Prefp = (maxion_twiddle**2.0d0)*(tvarphi_cp**2.0d0/2.0d0 + tvarphi_sp**2.0d0/2.0d0 -&
            & tvarphi_c*tvarphi_sp + tvarphi_s*tvarphi_cp)
       wEFA_c_upd = (Prefp/rhorefp)/((lhETA/(maxion_twiddle*hnot*a))**2._dl)
       call lh(..., lhETA_upd, ..., rhorefp)     ! <H> from Friedmann with EFA density

   Convergence test: ``abs(wEFA_c_upd/Params%wEFA_c - 1.0_dl) .lt. tol_EFA`` → exit;
   else ``Params%wEFA_c = wEFA_c_upd; lhETA = lhETA_upd`` and repeat. Warning printed
   if ``iter_EFA .eq. maxiter``.

   Note dimensional bookkeeping comments (lines 1273–1282): the a² factors per species
   and the final ``/lhETA**2`` implement dH²/d(mt) = (dH²/dlna)·(H/m) in
   (H0·m)-normalized units.
5. So: ``littlehauxi`` = instantaneous (aH)_osc → ``Params%ah_osc``; ``lhETA`` =
   ⟨aH⟩_EFA → ``Params%ahosc_ETA``; rhorefp/Prefp are cycle-averaged.

**Side effect:** ``Params%wEFA_c`` keeps its converged value between calls (it is the
*input* of the next call), so the whole shooting loop is itself a fixed-point iteration
on wEFA_c. Order of operations matters for exact reproducibility.

3.8 Density evolution after the switch [PHYSICS]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There is **no table beyond a_osc** (comment line 1095: "never use the lookup table to
define density after aosc"). Post-switch density is analytic everywhere in the code:

.. code-block:: fortran

    wcorr_coeff = Params%ahosc_ETA*aosc_guess(j)/(maxion_twiddle*hnot)            ! line 853
    omaxh2_wcorr = rhorefp*(aosc_guess(j)**3.0d0)*dexp((wcorr_coeff**2.0d0)*3.0d0*&
         &Params%wEFA_c*(1.0d0 - 1.0d0/(aosc_guess(j)**4.0d0))/4.0d0)             ! lines 857–858 (a=1)

General form (used in equations_ppf.f90:289–292, 1791–1796 etc.):
``ρ_ax(a)h² = rhorefp·(a_osc/a)³·exp(3·wEFA_c·wcorr_coeff²·(1/a⁴ − 1/a_osc⁴)/4)``
This is ρ ∝ a⁻³·exp(3∫w dln a) with w(a) = wEFA_c·(H/m)² and H ∝ a⁻² (RD scaling),
with wcorr_coeff = (⟨H⟩/m)·a_osc² so H/m = wcorr_coeff/a².

3.9 m/H0 < 10 special case (dark-energy-like) [PHYSICS]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Lines 755–761: ``if (maxion_twiddle .lt. 10._dl)`` — **no switch is performed at
all**; the KG field is evolved to a=1 and

.. code-block:: fortran

    omaxh2_guess(j) = (v_vec(2,ntable)/a_arr(ntable))**2.0d0+(maxion_twiddle*v_vec(1,ntable))**2.0d0

is the present-day density used for shooting. ``aosc_guess(j)`` stays at sentinel
15.0, so ``Params%a_osc = 15`` on convergence (line 898 — *not capped*);
``modules.f90:417`` later caps ``a_osc=1``, and the driver/``equations_ppf`` use
``CP%a_osc .le. 1`` to decide whether a switch exists. In this regime the exact field
tables (phinorm/phidotnorm) cover all of 0 < a ≤ 1 and the perturbations stay in KG
form throughout.
(Edge case: for 3 ≲ m̃ < 10 a grid crossing of m = dfac·H may have set
``aosc_guess(j)`` to a rough grid value before the m̃<10 branch was taken; that value
would be stored unrefined in ``Params%a_osc``. Given dfac ≥ 10 and H(a=1)=H0 this can
only marginally trigger near m̃→10.)

3.10 Recombination-skip support [PHYSICS/ACCURACY]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Lines 863–886: if ``a_skipst ≤ aosc < a_skip`` (driver sets ``P%a_skip = 1/801``,
``P%a_skipst = 1/1301``, inidriver_axion.F90:520–521), compute the conformal Hubble at
``a_skip`` with the EFA-extrapolated axion density passed as ``rho_f``:

.. code-block:: fortran

    call lh(..., Params%a_skip, v_vec(1:2,1), lh_skip, ..., &
         &rhorefp*((aosc_guess(j)/Params%a_skip)**3.0d0)*dexp((wcorr_coeff**2.0d0)*3.0d0*&
         &Params%wEFA_c*(1.0d0/(Params%a_skip**4.0d0) - 1.0d0/(aosc_guess(j)**4.0d0))/4.0d0))
    Params%dfac_skip = min(littlehauxi/Params%ahosc_ETA, 1._dl)*lh_skip
    Params%dfac_skip = (maxion_twiddle*Params%a_skip*hnot/Params%dfac_skip)

i.e. the dfac value that would relocate the switch to z = 800 (after recombination);
the driver re-runs ``w_evolve`` with ``P%dfac = P%dfac_skip`` (§6) to avoid placing
the KG→EFA switch inside the recombination window.

3.11 Output-table construction [PHYSICS + numerics]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After the shooting loop (lines 1014–1112):

- ``phinorm_table(i) = v_vec(1,i)``, ``phidotnorm_table(i) = v_vec(2,i)``,
  ``rhoaxh2_ov_rhom(i)`` as in §2.
- **Base change** line 1041: ``loga_table=dlog10(dexp(loga_table))`` — from here on
  ``loga_table`` is **log10(a)** (all downstream ``spline_out`` calls use
  ``dlog10(a)``).
- Spline second-derivative tables with analytic boundary derivatives (note the
  ``dlog(10)`` Jacobians):

  .. code-block:: fortran

      d1 = dlog(10._dl)*phidotnorm_table(1)*hnot/littlehfunc(1)
      ...
      call spline(loga_table(1:ntable),phinorm_table,ntable,d1,d2,phinorm_table_ddlga)
      d1 = -dlog(10._dl) * (2._dl*phidotnorm_table(1) + &
           &((maxion_twiddle*(10._dl**(loga_table(1))))**2._dl)*phinorm_table(1)*hnot/littlehfunc(1))
      ...
      call spline(loga_table(1:ntable),phidotnorm_table,ntable,d1,d2,phidotnorm_table_ddlga)

- Re-evaluate ``v1_ref``, ``v2_ref`` at ``dlog10(Params%a_osc)`` from the splines
  (lines 1063–64) and make the **final** ``auxiIC`` **call** (1068–71); then store:

  .. code-block:: fortran

      Params%ah_osc = littlehauxi
      Params%A_coeff = A_coeff
      Params%A_coeff_alt = A_coeff + 2.0d0*littlehauxi/(Params%a_osc*hnot*maxion_twiddle)
      Params%tvarphi_c/s/cp/sp = ...
      Params%rhorefp_ovh2 = rhorefp/hsq
      Params%Prefp = Prefp

  (``Params%ahosc_ETA`` was set inside ``auxiIC`` via the ``lhETA`` dummy argument.)
- Density log-table (lines 1105–1112):

  .. code-block:: fortran

      rhoaxh2ovrhom_logtable=dlog10(rhoaxh2_ov_rhom)
      d1 = -6.0_dl*(phidotnorm_table(1)/(10._dl**(loga_table(1))))**2._dl/rhoaxh2_ov_rhom(1)
      d2 = -6.0_dl*(phidotnorm_table(ntable)/(10._dl**(loga_table(ntable))))**2._dl/rhoaxh2_ov_rhom(ntable)
      call spline(loga_table,rhoaxh2ovrhom_logtable,d1,d2,rhoaxh2ovrhom_logtable_buff)

  (dlog10ρ/dlog10a = −6·KE/ρ, exact for KG.) An older block rescaling the table to w=0
  after a_osc is commented out (1098–1104) — the table is *purely* the KG solution;
  EFA extrapolation is done analytically by consumers.

3.12 Matter–radiation equality [PHYSICS]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Lines 1122–1145:

.. code-block:: fortran

    eq_arr(i)=((omegah2_regm/(a_arr(i)**3.0d0)+rhoaxh2_ov_rhom(i))/(Params%omegah2_rad/(a_arr(i)**4.0d0)))
    eq_arr=dlog(eq_arr)
    d1=(loga_table(2)-loga_table(1))/(eq_arr(2)-eq_arr(1))   ! finite-difference BCs
    ...
    call spline(eq_arr,...,loga_table,...); call spline_out(...,0.0d0,Params%aeq)
    Params%aeq=10._dl**(Params%aeq)        ! RL fixed 012524 (loga_table is log10 here)

(axions counted as matter via their *exact* KG density; massive ν counted as
radiation). Fallback if the spline breaks when ``a_osc < regzeq``:

.. code-block:: fortran

    regzeq=(Params%omegah2_rad+sum(lhsqcont_massive))/(omegah2_b+omegah2_dm+omegah2_ax)
    if (Params%a_osc.lt.regzeq) Params%aeq=regzeq

Also sets the **module-global** ``aeq_LCDM`` (ModelParams; lines 1136–1139), the
pure-LCDM analytic a_eq assuming axions scale as matter:

.. code-block:: fortran

    aeq_LCDM = (((Params%TCMB**4.0d0)/(rhocrit))/(c**2.0d0)*a_rad*1.d1/(1.d4))*&
         &(1._dl + (Params%nu_massless_degeneracy + &
         &sum(Params%Nu_mass_degeneracies(1:Params%Nu_mass_eigenstates)))*(7._dl/8._dl)*&
         &((4._dl/11._dl)**(4._dl/3._dl)))/(omegah2_b+omegah2_dm+omegah2_ax)

used for the photon-oscillation phase targeting (``get_phase_info``, driver).

3.13 Isocurvature normalization [PHYSICS]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Lines 1151–1165:

.. code-block:: fortran

    Params%phiinit=vtwiddle_init*sqrt(6.0d0)
    if (Params%axion_isocurvature) then
       Params%amp_i = Params%Hinf**2/(pi**2*Params%phiinit**2)
       Params%r_val  = 2*(Params%Hinf**2/(pi**2.*Params%InitPower%ScalarPowerAmp(1)))
       Params%alpha_ax = Params%amp_i/Params%InitPower%ScalarPowerAmp(1)
    end if
    Params%omegar=Params%omegah2_rad/hsq

``Hinf`` is read by the driver as log10(GeV) and converted to H_inf/M_planck
(``P%Hinf = (10**P%Hinf)/mplanck``, inidriver_axion.F90:316–317); ``phiinit`` is
φ/M_pl(reduced) — so ``amp_i`` = (H_inf/π φ_init)² is the axion isocurvature power
amplitude, ``r_val`` the tensor-to-scalar ratio, ``alpha_ax`` the iso fraction. These
feed the ``initial_condition = 6`` mode and CosmoMC outputs.

4. Shooting / root-finding (what is solved for)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Unknown: the initial field value** ``vtwiddle_init`` **(= v1 at a_init)**, such that
the present-day axion density ``omaxh2`` (EFA-extrapolated, or exact field value for
m̃<10) equals the input ``omegah2_ax = Params%omegaax·h²``. H0 is **not** shot for; it
is fixed input. (The original axionCAMB swept ``nphi=150`` trial values + spline
inversion — that code is commented out; RL replaced it with bracketing+bisection.)

Structure (lines 652–1006):

- 3-element working arrays ``v1_initguess(1:3)``, ``omaxh2_guess(1:3)``,
  ``aosc_guess(1:3)``.
- **Analytic first guess** (lines 598–613), using m/H\ :sub:`*` = 3 asymptotics:

  .. code-block:: fortran

      if (maxion_twiddle .lt. 3.0d0) then                 ! never oscillates (DE-like)
         v1_initguess(2) = dsqrt(omegah2_ax)/maxion_twiddle
      else if ((maxion_twiddle**2.0d0)/9.0d0 .lt. &
           &(omegah2_m**4.0d0)/(Params%omegah2_rad**3.0d0) + (Params%omegah2_rad**3.0d0)/(omegah2_m**4.0d0)) then  ! osc in MD
         v1_initguess(2) = dsqrt(omegah2_ax)/(3.0d0*dsqrt(omegah2_m)/hnot)
      else                                                 ! osc in RD
         v1_initguess(2) = dsqrt(omegah2_ax)/(((9.0d0*Params%omegah2_rad/hsq)**0.375d0)*(maxion_twiddle**0.25d0))
      end if

  bracket = [guess/2, guess·2]; middle = arithmetic mean.
- Sentinels: ``omaxh2_guess = 42.0d0`` = "not computed"; recompute condition
  ``omaxh2_guess(j) .gt. 1.0d0`` (line 660); ``aosc_guess = 15.0d0`` = "no oscillation
  found".
- Loop ``do while (iter_c < nphi)`` — ``nphi = 150`` **(modules.f90:75) reused as max
  bisection iterations**. Each pass integrates the full KG history for whichever of
  the 3 candidates needs recomputation.
- **Convergence:** ``abs(omaxh2_guess(j)/omegah2_ax - 1.0d0) .lt. 1.0d-6`` (line 893)
  → accept (``vtwiddle_init``, ``Params%a_osc = aosc_guess(j)``), set ``iter_c=-1``
  and exit.
- **Bracket expansion** (lines 925–933): if not straddling, double
  ``v1_initguess(3)`` or halve ``v1_initguess(1)`` (monotonicity assumed); once
  straddling, ``bisec_bracketed=.true.`` and standard bisection (replace the end with
  the same sign as midpoint; midpoint = mean, line 971). Bracket-lost-after-bisection
  ⇒ hard ``stop`` (line 963–965).
- **Failure path** (lines 909–923): at ``iter_c .eq. nphi-1`` print warning, accept
  middle guess, cap ``Params%a_osc`` to 1.
- **Adaptive a_final** (lines 977–993): once ``omaxh2_guess(3)/omegah2_ax < 1.1`` and
  a_osc·1.1 < 1, shrink the grid: ``a_final = aosc_guess(3)*1.1d0``, recompute
  ``dloga``, ``loga_table``, ``a_arr``. So in the oscillating case the final tables
  **end at ≈1.1·a_osc, not at a=1** — all consumers must use the analytic EFA density
  beyond a_osc.

5. Outputs, data structures and couplings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

5.1 Global tables (declared in ``modules.f90:249–255``, module ModelParams) [PLUMBING+PHYSICS]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: fortran

    integer :: ntable = 5000   ! modules.f90:249 — but ALWAYS overwritten by driver: ntable = nint(P%dfac*100) + 1
    integer, parameter :: nphi = 150          ! modules.f90:75
    real(dl) :: aeq_LCDM                      ! modules.f90:250
    real(dl), allocatable :: loga_table(:), phinorm_table(:), phidotnorm_table(:), phinorm_table_ddlga(:)
    real(dl), allocatable :: phidotnorm_table_ddlga(:), rhoaxh2ovrhom_logtable(:), rhoaxh2ovrhom_logtable_buff(:)

All length ``ntable``; abscissa ``loga_table`` = **log10(a)** spanning
[log10(a_init), log10(a_final≈1.1·a_osc or 1)]; contents:

- ``phinorm_table`` = v1 (φ in √(3/4πG) units), ``phidotnorm_table`` = v2 (φ′ in
  √(3/4πG)·H0 units), with ``_ddlga`` natural-cubic-spline 2nd-derivative buffers
  (analytic first-derivative BCs).
- ``rhoaxh2ovrhom_logtable`` = log10(Ω_ax(a)h²) with ``_buff`` spline buffer.

5.2 CAMBparams fields written by ``w_evolve``/``auxiIC``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``omegah2_rad`` (299–308), ``wEFA_c`` (335 init; mutated in auxiIC 1313),
``dfac_skip`` (877–879), ``a_osc`` (898/917), ``ahosc_ETA`` (dummy arg of auxiIC,
lines 843/1071), ``ah_osc`` (1075), ``A_coeff`` (1077), ``A_coeff_alt`` (1079),
``tvarphi_c, tvarphi_s, tvarphi_cp, tvarphi_sp`` (1081–84), ``rhorefp_ovh2`` (1091),
``Prefp`` (1092), ``aeq`` (1132/1142), ``phiinit`` (1151), ``amp_i, r_val, alpha_ax``
(1160–62), ``omegar`` (1168).
Fields read: ``H0, omegab, omegac, omegan, omegav, omegaax, m_ovH0, TCMB, dfac,
a_skip, a_skipst, Nu_mass_eigenstates, Nu_mass_fractions, Nu_mass_degeneracies,
Nu_massless_degeneracy, nu_massless_degeneracy, axion_isocurvature, Hinf,
InitPower%ScalarPowerAmp(1)``.
Declarations: modules.f90:123–127, 152, 182, 195.

5.3 Globals set as side effects [PLUMBING — must be re-derived in 1.6.7]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``grhom, grhog, grhor`` (ModelParams, modules.f90:213) at lines 294–296 —
  recomputed *before* CAMBParams_Set runs; in modern CAMB use ``State%grhom`` etc.,
  do not re-set.
- ``Nu_masses(k)`` (MassiveNu, modules.f90:220) at line 315 — requires
  ``init_massive_nu`` already called (driver does so at inidriver_axion.F90:519
  before ``w_evolve``).
- ``aeq_LCDM`` (lines 1136–1139).
- ``loga_table`` & friends (§5.1).

5.4 External call sites of public symbols
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``use axion_background``: inidriver_axion.F90:23.
- ``call w_evolve(P, badflag)``: inidriver_axion.F90:523, 536, 552, 572, 587, 611,
  627.
- ``call get_phase_info(P, y_phase, beta_coeff, movHETA_new, twobeta_new)``:
  inidriver_axion.F90:537, 553, 573, 588, 612.
- ``derivs_bg``, ``auxiIC``, ``lh``, ``next_step``: no external callers (internal
  only).
- Table consumers (representative; all via
  ``spline_out(loga_table, …, dlog10(a), …)``):

  - equations_ppf.f90:269–278 (``a_min = 10._dl**(loga_table(1))``; density lookup
    ``rhoaxh2ovrhom_logtable`` → ``grhoaxh2_ov_grhom``, with analytic EFA formula for
    a > a_osc at 289–292);
  - equations_ppf.f90:1023–24, 1770–71, 1942–43, 2540, 2987–88, 3139–40, 3785–86,
    3987–88 (``phinorm_table``/``phidotnorm_table`` → ``v1_bg``, ``v2_bg`` for KG
    perturbation sources);
  - cmbmain.f90:824–827 (taustart sanity check against ``10.0**(loga_table(1))``);
    cmbmainOMP.f90 analogous;
  - inidriver_axion.F90:726–732 (final deallocation of all 7 tables);
  - aeq_LCDM: inidriver_axion.F90:527, 533.
- Switch-parameter consumers (for context): equations_ppf.f90 reads
  ``CP%a_osc, CP%rhorefp_ovh2, CP%Prefp, CP%wEFA_c, CP%ah_osc, CP%ahosc_ETA,
  CP%A_coeff, CP%tvarphi_*`` extensively (e.g. 289–292, 1034–1104, 1790–96);
  modules.f90:414–418 (``grhoax=grhom*CP%omegaax``, ``a_osc`` capping), 603–606 &
  636–639 (``DeltaTime``/``DeltaPhysicalTimeGyr`` rombint split at ``CP%a_osc``), 440
  (Reionization_Init gets ``CP%a_osc``), 2152, 2899–2901 (sound horizon split);
  cmbmain.f90:811–812 (``taueq`` from ``CP%aeq``).

5.5 ``get_phase_info`` content [PHYSICS]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: fortran

    y_beta = Params%a_osc/aeq_LCDM
    movHETA = Params%dfac*Params%ah_osc/Params%ahosc_ETA          ! = m/<H> at switch
    beta_coeff = (4._dl*(y_beta**2 - y_beta - 2.0_dl + 2.0_dl*sqrt(1.0_dl + y_beta)))/(3._dl*(y_beta**2))
    beta2x = movHETA*beta_coeff - const_pi*3._dl*(1.0_dl + y_beta)/(4.0_dl + 3.0_dl*y_beta)

``beta2x`` ("2β") is the analytic WKB phase of the field oscillation at the switch in
a matter+radiation universe; used by the driver to pick ``dfac`` so the switch lands
at a fixed phase.

6. Driver-level orchestration around ``w_evolve`` (inidriver_axion.F90:508–636) [PLUMBING/PHYSICS]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is outside the file but inseparable from its algorithm; in CAMB 1.6.7 it must
move into the equivalent of ``CAMBdata%SetParams``/background init:

1. ``call init_massive_nu(P%omegan /=0)``; ``P%a_skip = 1/801``;
   ``P%a_skipst = 1/1301``; ``P%dfac_skip = 0``; ``P%dfac = 10`` and
   ``ntable = nint(P%dfac*100)+1`` (lines 154–155) — **table resolution is 100 points
   per dfac unit**.
2. First ``w_evolve`` (line 523).
3. **Phase targeting** (526–620): if ``P%dfac < 23 .and. P%m_ovH0 ≥ 10 .and. P%ma <
   1.e-25 .and. P%a_osc*(P%omegaax/(P%omegac+P%omegaax))/aeq_LCDM > 0.03 .and.
   P%a_osc < P%a_skipst``: target ``twobeta_tgt = 7.08_dl*const_pi``; first guess
   ``P%dfac = twobeta_tgt + 0.75π − twobeta_tgt²/(4(twobeta_tgt + 2 m̃
   aeq_LCDM^1.5/√(2(Ω_c+Ω_b+Ω_n+Ω_ax))))``; then shooting+bisection on
   ``hosc = P%ah_osc/P%a_osc`` (rescaling ``P%dfac *= (P%ah_osc/P%a_osc)/hosc_new``,
   re-running ``w_evolve`` + ``get_phase_info`` each step; tolerance
   ``beta_tol = 2.e-2*π``; ≤500 iters each phase).
4. **Recombination skip** (622–631): up to 500 times, while
   ``a_skipst ≤ P%a_osc < P%a_skip(1−1e-2)``:
   ``P%dfac = P%dfac_skip; ntable = nint(P%dfac*100)+1; call w_evolve`` — pushes the
   switch past z=800.

7. Numerical subtleties & gotchas (porting checklist)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **[PHYSICS]** ``lh`` **returns conformal aH** in 100 km/s/Mpc units, not H — every
   ``hnot/littlehfunc`` factor in the file relies on this. Name is misleading.
2. **[PHYSICS] Fixed-step RK8 in ln a** — no error control (embedded estimate
   discarded). Accuracy is governed solely by ``ntable = nint(dfac*100)+1``. The
   5000-point comment (line 702) predates the adaptive-a_final change; the effective
   resolution requirement is "≈100 steps per e-fold-equivalent of dfac" with the table
   truncated at 1.1·a_osc. In modern CAMB consider replacing with dverk/adaptive
   integrator, but verify low-ℓ ISW (the stated reason for the 8th-order method).
3. **[ACCURACY] Tolerances:** density shooting ``1e-6`` relative (line 893); EFA fixed
   point ``tol_EFA=1e-7``, ``maxiter=30``; aosc root via natural spline
   (``d1=d2=1e50``).
4. **[PHYSICS] log-base switch:** ``loga_table`` is ln(a) during integration, log10(a)
   from line 1041 onward and in all global tables. The ``dlog(10)`` factors in the
   spline boundary derivatives (1047–56) and ``Params%aeq=10**…`` exist for this
   reason. Single most likely source of porting bugs.
5. **[PHYSICS]** ``Params%wEFA_c`` **is stateful**: initialized 9/8 each ``w_evolve``
   call (line 335), then iterated inside ``auxiIC`` and *fed back* into both the
   shooting density correction and the next auxiIC call. Reproducing AxiECAMB numbers
   requires the same call order.
6. **[PHYSICS] Sentinels:** ``aosc_guess=15.0`` (no oscillation),
   ``omaxh2_guess=42.0`` (not computed), ``omaxh2_guess>1.0`` triggers recompute,
   ``a_osc=15`` may be stored in Params (capped to 1 in modules.f90:417). Replace with
   explicit logicals in the port.
7. **[PHYSICS] m̃ thresholds:** ``< 3`` → DE-like initial-guess formula; ``< 10`` →
   never switch, evolve KG to a=1 (perturbations also stay KG); ``≥ 10`` → EFA switch
   at m = dfac·H. Driver phase-tuning only for ``m_ovH0 ≥ 10`` and ``ma < 1e-25`` eV.
8. **[PHYSICS] v2 attractor IC** with the 1/5 factor (line 676) — replaces the old
   "start at rest" IC; required so the adiabatic sound speed is correct at early
   times.
9. **[PLUMBING] Side-effect setting of**
   ``grhom/grhog/grhor/Nu_masses/omegah2_rad/omegar`` duplicates standard CAMB init —
   in 1.6.7 read them from ``CAMBdata`` (they are computed in ``CAMBdata%SetParams``
   before background classes initialize). Ensure the H(a) inside the axion solver
   matches ``CAMBdata`` exactly (incl. massive-ν ``Nu_rho`` interpolation).
10. **[ACCURACY]** ``eq_arr`` **spline BCs** (lines 1127–28) are crude finite
    differences (unlike everything else, which uses analytic BCs) — fine, but don't
    "fix" silently if matching outputs.
11. **[OBSOLETE/COSMETIC] Dead variables:** ``littlehfunc_buff`` (allocated, never
    used), ``nstop``, ``H_ev``, ``d1v1/d2v1/d1v2/d2v2``, ``movH_test``, ``y_phase``
    (in w_evolve), ``dlogvtwiddle_init``, ``vtwiddle_initmax/min``,
    ``vtwiddle_initlist``, ``eq_arr_buff`` is used; large commented-out blocks (old
    nphi-sweep shooting, NaN checks, DG neutrino-fix history at 268–291).
    ``lhsqcont_massless`` is passed to ``lh``/``derivs_bg``/``next_step`` but never
    used in their bodies (it is already folded into ``omegah2_rad``).
12. **[PHYSICS] Tables end at 1.1·a_osc** in the oscillating case — any consumer in
    the new architecture must guard ``a ≤ a_osc`` and use the analytic EFA density
    (§3.8) beyond. equations_ppf.f90:269 uses ``a_min = 10**loga_table(1)`` similarly
    for the early end.
13. **[PLUMBING] badflag**: only ever set to 1 (never reset here); flags
    collapsing/NaN histories so CosmoMC-style callers can reject the sample.
14. **[ACCURACY]** ``ntable`` **lives in ModelParams as a non-parameter integer** so
    the driver can resize it between ``w_evolve`` calls; in 1.6.7 this should become a
    member of the axion component class, sized per dfac.

8. Classification summary
~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Item
     - Class
   * - KG system, Friedmann with field energy, v2 attractor IC, switch criterion
       m = dfac·H, aosc spline root, auxiIC EFA matching (tvarphi projections,
       A_coeff, ⟨H⟩/wEFA_c fixed point), post-switch ρ(a) with w-correction
       exponential, aeq/aeq_LCDM, phiinit/isocurvature amplitudes, m̃<10 DE-like
       branch, dfac_skip computation, get_phase_info formulas
     - **[PHYSICS]** — port verbatim
   * - Bisection shooting for vtwiddle_init (tolerance 1e-6, bracket ×2/÷2, nphi cap),
       adaptive a_final = 1.1·aosc, ntable = nint(dfac·100)+1, fixed-step RK8 +
       Butcher tableau, natural-spline BCs, tol_EFA/maxiter
     - **[ACCURACY]** — physics-motivated; keep values, but integrator/iteration
       structure may be modernized if validated against low-ℓ ISW
   * - Global tables in ModelParams, deallocate-on-entry, Params field wiring, driver
       loops (phase targeting, recomb skip), badflag plumbing
     - **[PLUMBING]** — re-derive for class-based CAMB (put tables + switch state in
       an axion component / CAMBdata)
   * - Recomputation of grhom/grhog/grhor/rhocrit/Nu_masses/omegah2_rad, commented-out
       nphi-sweep & DG-fix archaeology, dead variables
     - **[OBSOLETE]**
   * - Massive comment blocks, "WORK NEEDED ON COMMENTS"
     - **[COSMETIC]**


Original-code analysis: parameters and background (modules.f90)
----------------------------------------------------------------

Source diff: ``/Users/vivianmiranda/data/research/WayneHu/rayne/.port_analysis/diffs/modules.f90.diff`` (4535 lines).
File analyzed: ``/Users/vivianmiranda/data/research/WayneHu/rayne/AxiECAMB/modules.f90`` (3327 lines) vs ``/Users/vivianmiranda/data/research/WayneHu/rayne/OLDCAMB/modules.f90`` (2807 lines).

**Global note:** The entire file was re-indented (4-space → 2/3-space, emacs style). The diff therefore touches nearly every line; the overwhelming majority of hunks are [COSMETIC] whitespace. Everything substantive is itemized below with AxiECAMB line numbers. Modules in this file: ``ModelParams``, ``lvalues``, ``ModelData``, ``MassiveNu``, ``Transfer``, ``ThermoData``. ``lvalues``, ``ModelData`` (except commented timing/debug code), and ``MassiveNu`` have **no functional changes**.

Modern-CAMB mapping context: in CAMB 1.6.7, ``ModelParams`` → ``model.f90 (CAMBparams)`` + ``results.f90 (CAMBdata)``, ``Transfer``/``MatterTransferData`` → ``results.f90``, ``ThermoData`` → ``results.f90 (TThermoData)``, constants like ``Transfer_*`` → ``results.f90``.

1. Module ``ModelParams`` — constants and new globals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1.1 ``max_transfer_redshifts`` 150 → 600 — [ACCURACY]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``AxiECAMB/modules.f90:61``

.. code-block:: fortran

    integer, parameter :: max_transfer_redshifts = 600!RL 04/30/25: default 500. 600 is to accomodate accuracy_boost = 6 with do_nonlinear and should not be needed for normal use ! COSMOSIS - alter number of transfer redshifts

Needed because ``Transfer_SetForNonlinearLensing`` was changed to use 100×AccuracyBoost redshifts (see §5.4). Modern CAMB 1.6.7 has ``max_transfer_redshifts = 256`` (model.f90); if the 100× NLL change is ported, this must grow accordingly (e.g. ≥ 100*AccuracyBoost+PK redshifts).

1.2 New parameter ``nphi`` — [PHYSICS] (support constant for axion background shooting)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``AxiECAMB/modules.f90:71-75``

.. code-block:: fortran

    !Size of axion integration array (number of time slices for homogeneous scalar field evolution)
    !!integer, parameter :: ntable = 300 ! RH added this so it is callable everywhere !RL: default is 5000
    !Number of scalar initial conditions to try to build a cubic spline
    !and thus determine correct initial condition for scalar field evolution
    integer, parameter:: nphi = 150 !RL: default 150

``nphi`` = number of trial initial field values used by the ``axion_background.F90`` shooting/spline algorithm to find ``phiinit`` that produces the requested ``omegaax``. Dimensionless count.

1.3 New module-level axion globals — [PHYSICS]/[PLUMBING]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``AxiECAMB/modules.f90:223-228``

.. code-block:: fortran

    !Axions  log10 density today, a_osc =a when m=dfac*H (see background module
    !axion background.f90, drefp_hsq=axion density (not log) in same units
    !when a=a_osc -- allows simple a^-3 scaling to be applied throughout code
    real(dl) grhoax, a_osc, tau_osc, drefp_hsq !RL added tau at oscillation

- ``grhoax`` — axion analogue of ``grhoc``: ``grhoax = grhom*CP%omegaax`` (κa²ρ units, Mpc⁻²). Set in ``CAMBParams_Set`` (§2.5).
- ``a_osc`` — module copy of ``CP%a_osc`` (scale factor of KG→EFA switch where ``m = dfac*H``), clamped to ≤ 1 in ``CAMBParams_Set``.
- ``tau_osc`` — module-level conformal time at switch (Mpc). (The authoritative value is ``CP%tau_osc``, set in ``init_background`` in equations_ppf.f90.)
- ``drefp_hsq`` — reference axion density (linear, not log) at ``a = a_osc`` in the background-table units; after the switch the code scales the axion density as ``drefp_hsq*(a_osc/a)^3``. Set by the background module.

1.4 New axion background tables + ``ntable``, ``aeq_LCDM`` — [PHYSICS]/[PLUMBING]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``AxiECAMB/modules.f90:248-255``

.. code-block:: fortran

    real(dl) :: AccuracyBoost = 1._dl !1. is the default value. ...
    integer :: ntable = 5000 !RL 111123. ntable should be properly set in inidriver_axion.
    real(dl) :: aeq_LCDM !RL 031924, added for photon oscillation skipping
    real(dl), allocatable :: loga_table(:), phinorm_table(:), phidotnorm_table(:), phinorm_table_ddlga(:)
    real(dl), allocatable :: phidotnorm_table_ddlga(:), rhoaxh2ovrhom_logtable(:), rhoaxh2ovrhom_logtable_buff(:) !RL 112823
    public loga_table, phinorm_table, phidotnorm_table, phinorm_table_ddlga, phidotnorm_table_ddlga
    public rhoaxh2ovrhom_logtable, rhoaxh2ovrhom_logtable_buff !RL 112823 - replacing CP tables with public tables
    public aeq_LCDM !RL 031924

- ``ntable`` (default 5000) — runtime size of the axion background time grid; allocated/filled by ``axion_background.F90``, value set in ``inidriver_axion.F90``. (Originally these tables were members of ``CP``; RL moved them to module-level allocatables — comment "replacing CP tables with public tables".)
- ``loga_table(:)`` — log10(a) grid for the background solution.
- ``phinorm_table(:)``, ``phidotnorm_table(:)`` — normalized homogeneous scalar field φ and dφ/dt on that grid; ``*_ddlga`` arrays are spline second derivatives w.r.t. log10(a).
- ``rhoaxh2ovrhom_logtable(:)`` + ``_buff(:)`` — log10 of axion density (in h² units relative to the code's grhom convention) table + spline buffer; this is what ``dtauda``/background interpolates pre-switch.
- ``aeq_LCDM`` — matter-radiation equality scale factor of the equivalent ΛCDM model, used for the "photon oscillation skipping" (recombination-skip) machinery.

**Port note:** in modern CAMB these must NOT be globals; they belong in the ``CAMBdata``/``TCosmoComponent``-style state object (re-derive as [PLUMBING]).

1.5 ``AccuracyBoost`` default — [COSMETIC]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Still ``1._dl`` (``modules.f90:248``); only a comment added. ``lSampleBoost=1._dl``, ``lAccuracyBoost=1.``, ``HighAccuracyDefault=.false.`` all unchanged.

1.6 Derived-parameter index list 10 → 13 — [OBSOLETE]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``AxiECAMB/modules.f90:292-295``

.. code-block:: fortran

    integer, parameter :: derived_age=1, derived_zstar=2, derived_rstar=3, derived_thetastar=4, derived_DAstar = 5, &
          derived_zdrag=6, derived_rdrag=7,derived_kD=8,derived_thetaD=9, derived_zEQ =10, derived_keq =11, &
          derived_thetaEQ=12, derived_theta_rs_EQ = 13
      integer, parameter :: nthermo_derived = 13

Old code had 10 (no DAstar, keq, theta_rs_EQ). This is a back-port of later mainline CAMB; CAMB 1.6.7 already has exactly this set in ``results.f90``. **Nothing to port** except the axion modification of z_EQ itself (§6.6).

2. ``CAMBparams`` type — new/changed fields
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All in ``AxiECAMB/modules.f90:123-195``. Verbatim declarations:

.. code-block:: fortran

    real(dl)  :: omegab, omegac, omegav, omegan, omegaax, ma, m_ovH0, dfac          ! :123
    real(dl)  :: a_osc, tau_osc, opac_tauosc, expmmu_tauosc, alpha_ax, r_val, omegah2_rad, amp_i, axfrac, omegada, Hinf !RL added tau at oscillation - the correct value will be calculated at init_background in equations_ppf  ! :124
    real(dl) :: ah_osc, ahosc_ETA, A_coeff, tvarphi_c, tvarphi_cp, tvarphi_s, tvarphi_sp, wEFA_c !RL added background EFA parameters at the switch  ! :125
    real(dl) :: A_coeff_alt !RL 043024 testing with changed A coeff                  ! :126
    real(dl) :: a_skip, dfac_skip, a_skipst !RL 012524 added for skipping recombination ! :127
    real(dl)  :: H0_in_Mpc_inv, H0_eV !RL for the ease of computing KG in pert       ! :130
    real(dl)  :: ratio                                                               ! :131
    real(dl)  :: Nu_massless_degeneracy !RL pasted DG's addition                     ! :144
    logical   :: use_axfrac, axion_isocurvature                                      ! :152
    real(dl) omegak, omegar, grhor          ! (was: real(dl) omegak)                 ! :182
    real(dl) :: phiinit,aeq, ainit, lens_amp,rhorefp_ovh2, Prefp !RL added Prefp     ! :195

Field-by-field:

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Field
     - Meaning / units
     - Set where
     - Class
   * - ``omegaax``
     - Ω_axion today (density parameter)
     - inidriver from ``m_ax``, ``omaxh2``/``axfrac``
     - [PHYSICS]
   * - ``ma``
     - axion mass m_ax in **eV** (printed as ``m_ax/eV``, modules.f90:491)
     - input
     - [PHYSICS]
   * - ``m_ovH0``
     - m_ax/H0, dimensionless (mass in units of H0)
     - inidriver
     - [PHYSICS]
   * - ``dfac``
     - switch threshold: switch when ``m = dfac*H`` (paper default 10)
     - input/inidriver
     - [PHYSICS]
   * - ``a_osc``
     - scale factor at the KG→EFA switch
     - axion_background
     - [PHYSICS]
   * - ``tau_osc``
     - conformal time (Mpc) at switch; computed in ``init_background`` (equations_ppf.f90:237-249) via ``DeltaTime(0, a_osc, in_tol=1.0d-8)``; forced to ``tau0+1`` if no switch before today, capped at ``tau0``
     - equations_ppf
     - [PHYSICS]
   * - ``opac_tauosc``, ``expmmu_tauosc``
     - opacity ȧκ̇ and exp(−κ) interpolated at τ_osc, cached by ``inithermo`` via ``ThermoSplineOut`` (§6.4); used by perturbation module at switch
     - ThermoData
     - [PHYSICS]/[PLUMBING]
   * - ``alpha_ax``
     - axion isocurvature amplitude parameter α_ax
     - input
     - [PHYSICS]
   * - ``r_val``, ``ratio``
     - tensor-to-scalar ratio carried in CP "for use in axion isocurvature i.c.'s" (comment :137)
     - input
     - [PHYSICS]
   * - ``omegah2_rad``
     - ω_radiation = Ω_r h² (photons + massless ν), assembled in inidriver (inidriver_axion.F90:269)
     - inidriver
     - [PLUMBING]
   * - ``amp_i``
     - primordial scalar amplitude copy used in iso normalization
     - inidriver
     - [PLUMBING]
   * - ``axfrac``
     - Ω_ax/Ω_dark fraction (used with ``use_axfrac``)
     - input
     - [PHYSICS]
   * - ``omegada``
     - total dark (ax+cdm) density parameter, from ``omdah2``
     - input
     - [PHYSICS]
   * - ``Hinf``
     - inflationary Hubble scale (log10(H_inf/Mpl) convention handled in inidriver) for isocurvature amplitude
     - input
     - [PHYSICS]
   * - ``ah_osc``
     - (aH) at switch
     - axion_background
     - [PHYSICS]
   * - ``ahosc_ETA``
     - variant of aH at switch (ETA-time convention)
     - axion_background
     - [PHYSICS]
   * - ``A_coeff`` (+\ ``A_coeff_alt``)
     - matching coefficient of EFA at the switch (WKB joining coefficient)
     - axion_background
     - [PHYSICS] (``A_coeff_alt`` is a test variant — check usage in equations before porting)
   * - ``tvarphi_c, tvarphi_cp, tvarphi_s, tvarphi_sp``
     - cosine/sine components of the dimensionless field (tvarphi ≡ √κ φ/√6 convention) and their derivatives at the switch; used for matching exact KG solution to EFA fluid variables
     - axion_background
     - [PHYSICS]
   * - ``wEFA_c``
     - leading correction coefficient to the EFA equation of state w(a) after the switch
     - axion_background
     - [PHYSICS]
   * - ``a_skip, dfac_skip, a_skipst``
     - parameters for "recombination skipping"/photon-oscillation skipping (RL 012524) — accelerates integration of rapid oscillations
     - inidriver/equations
     - [PHYSICS] (verify usage in cmbmain/equations report)
   * - ``H0_in_Mpc_inv``
     - H0 in Mpc⁻¹ (= H0[km/s/Mpc]/c·1000)
     - inidriver
     - [PLUMBING]
   * - ``H0_eV``
     - H0 expressed in eV
     - inidriver
     - [PLUMBING]
   * - ``Nu_massless_degeneracy``
     - N_eff for massless ν; **promoted from a local variable in CAMBParams_Set to a CP field**, now set in ``inidriver_axion.f90:215,231``
     - inidriver
     - [PLUMBING] — see §3.1 risk
   * - ``use_axfrac``
     - if true, interpret inputs as (omdah2, axfrac) instead of omaxh2
     - input
     - [PHYSICS]
   * - ``axion_isocurvature``
     - enables axion isocurvature mode (initial_condition 6)
     - input
     - [PHYSICS]
   * - ``omegar`` (new on :182)
     - Ω_radiation; passed to ``Recombination_Init`` and used in equations_ppf adiabatic IC corrections (equations_ppf.f90:2670-2671)
     - inidriver
     - [PLUMBING]
   * - ``grhor`` (new on :182)
     - per-ν-species κa²ρ; duplicated into CP: ``P%grhor = (7.0d0/8.0d0)*((4.0d0/11.0d0)**(4.0d0/3.0d0))*grhog`` (inidriver_axion.F90:241) so the driver can compute ω_rad before ``CAMBParams_Set``
     - inidriver
     - [PLUMBING] (redundant with module global; re-derive cleanly in port)
   * - ``phiinit``
     - initial scalar field value found by shooting (units of the background module's normalized field)
     - axion_background
     - [PHYSICS]
   * - ``aeq``
     - matter-radiation equality scale factor (incl. axions per background solve)
     - axion_background
     - [PHYSICS]
   * - ``ainit``
     - initial scale factor of background field integration
     - axion_background
     - [PHYSICS]
   * - ``lens_amp``
     - lensing amplitude scaling (A_lens-like; check inidriver/lensing usage)
     - input
     - [PLUMBING]
   * - ``rhorefp_ovh2``
     - axion density at reference point (a_osc) over h² — companion of ``drefp_hsq``
     - axion_background
     - [PHYSICS]
   * - ``Prefp``
     - axion pressure at reference point (a_osc)
     - axion_background
     - [PHYSICS]

Commented-out (ignore, [COSMETIC]): ``RHCl_temp*`` arrays (:134-136), ``grhoax_table``/``cs2_table`` CP members (:191-194).

3. ``CAMBParams_Set`` changes (``AxiECAMB/modules.f90:308-510``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

3.1 Massive-neutrino bookkeeping block REMOVED — [PLUMBING] (high-risk)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The entire OLDCAMB block (OLDCAMB modules.f90 ~283-310) was deleted from ``CAMBParams_Set``:

- consistency check ``Num_Nu_Massive == sum(Nu_mass_numbers)``,
- conversion of massive→massless when ``Omegan==0``,
- ``nu_massless_degeneracy`` computation incl. ``share_delta_neff`` logic (``neff_i = fractional_number/(actual_massless + CP%Num_Nu_massive)`` etc.),
- ``Nu_mass_fractions`` sum check.

This logic now lives in ``inidriver_axion.f90`` (~lines 215-231) which sets ``P%Nu_massless_degeneracy`` before calling Set. Consequence inside Set:

.. code-block:: fortran

    grhornomass=grhor*CP%Nu_massless_degeneracy !RL fixed 020625      ! modules.f90:398

and the N_eff feedback print uses ``CP%Nu_massless_degeneracy`` (modules.f90:499).

**Port note:** CAMB 1.6.7 already does this correctly inside ``CAMBparams`` / ``TCAMBdata%SetParams``; do **not** port the relocation — only make sure the axion driver path doesn't bypass the standard neutrino setup. Anyone calling AxiECAMB's ``CAMBParams_Set`` without the driver gets garbage ``grhornomass`` (uninitialized field) — this is a wiring fragility, not physics.

3.2 ``GetOmegak()`` call removed — [PLUMBING]/[PHYSICS-adjacent]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``AxiECAMB/modules.f90:357-359``

.. code-block:: fortran

    !!write(*, *) 'CP%omegak before GetOmegak()', CP%omegak !RL 032724: CP%omegak is already assigned before GetOmegak(). Tested that in all cases the difference between these two results are very small (mostly 1e-16 level). Since GetOmegak() is physically wrong, I removed it and readjusted H0 instead, although in practice the effect is very small.
    !!CP%omegak = GetOmegak()

``CP%omegak`` is now assigned by the driver, computed **including radiation self-consistently**, and H0 is readjusted there to close the budget (Ω_b+Ω_c+Ω_ax+Ω_ν+Ω_r+Ω_Λ+Ω_k = 1). The omega-closure logic itself is in ``inidriver_axion.F90`` (separate report). In modern CAMB, ``omegak`` is derived in ``CAMBparams%SetCurv``-equivalent code from ``omk`` input; the axion port must add Ω_ax to the budget at that point.

3.3 grhog/grhor expressions — [COSMETIC]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Rewritten with explicit ``_dl`` exponents, numerically identical:

.. code-block:: fortran

    grhog = (kappa/c**2._dl)*4._dl*sigma_boltz/(c**3._dl)*(CP%tcmb**4._dl)*(Mpc**2._dl)   ! :385
    grhor = 7._dl/8._dl*(4._dl/11._dl)**(4._dl/3._dl)*grhog                                ! :387

3.4 Axion density + a_osc clamp — [PHYSICS]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``AxiECAMB/modules.f90:412-419``

.. code-block:: fortran

    ! Axions
    grhoax=grhom*CP%omegaax
    a_osc=CP%a_osc
    ! DM: Set a_osc=1 in case scalar field crashes for bad cosmology
    if(a_osc>1) then
       a_osc=1.
    end if

(Module-level ``a_osc`` is the clamped copy; ``CP%a_osc`` retains the raw value — note other code mixes both; the timestep logic uses ``CP%a_osc .le. 1._dl`` as the "is there a switch" test.)

3.5 ``init_massive_nu`` call moved out — [PLUMBING]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``AxiECAMB/modules.f90:435``

.. code-block:: fortran

    !call init_massive_nu(CP%omegan /=0) !RL commnented out 07/10/23
    call init_background

``init_massive_nu`` is instead called from ``inidriver_axion.F90:519`` *before* ``CAMBParams_Set``, because the axion background solver (called by the driver before Set) needs ``nu_masses``. In modern CAMB the equivalent ordering constraint is: massive-ν tables must be initialized before the axion background shooting; re-derive in the new architecture.

3.6 ``CP%tau0`` assignment removed from Set — [PLUMBING]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

OLDCAMB had ``CP%tau0=TimeOfz(0._dl)`` right after ``init_background``. In AxiECAMB, ``tau0`` is computed *inside* ``init_background`` (equations_ppf.f90:237: ``CP%tau0=TimeOfz(0._dl) !RL 061924 moved tau0 here so to initialte tau_osc``), because ``tau_osc`` must be derived there too. ``last_tau0`` caching logic unchanged (modules.f90:439,443).

3.7 ``Reionization_Init`` signature gains ``CP%a_osc`` — [PLUMBING]/[PHYSICS]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``AxiECAMB/modules.f90:440``

.. code-block:: fortran

    if (WantReion) call Reionization_Init(CP%Reion,CP%ReionHist, CP%YHe, akthom, CP%tau0, CP%a_osc, FeedbackLevel)

(old: no ``CP%a_osc`` argument). The reionization module needs a_osc to split its internal optical-depth/time integrations at the switch discontinuity (see reionization.f90 report).

3.8 Feedback prints — [COSMETIC] except noted
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``AxiECAMB/modules.f90:469-505``: added

.. code-block:: fortran

    write(*,'("H0                   = ",f10.6)') CP%H0
    write(*,'("Om_ax h^2            = ",f9.6)') CP%omegaax*(CP%H0/100)**2
    write(*,'("m_ax/eV              = ",e9.2)')  CP%ma !!/(100/.3e5)

a_eq print suppressed (commented, :485). The ν mass conversion uses constant ``elecV`` instead of ``eV`` (:501) — the constant was renamed in AxiECAMB ``constants.f90`` (name clash avoidance); [PLUMBING] detail for the constants report.

3.9 Misc declarations — [COSMETIC]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Added in Set: ``real clock_start, clock_stop``, ``integer :: i_check``, ``real(dl) dtauda / external dtauda`` (unused remnants), removed ``external GetOmegak``. Commented-out thetaMC file dump (OxFish) — ignore.

4. Background/distance functions in ``ModelParams``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

4.1 ``DeltaTime`` split at a_osc — [PHYSICS]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``AxiECAMB/modules.f90:587-612``

.. code-block:: fortran

    !RL modified for ULA switch
    if (a1 .lt. CP%a_osc .and. a2 .ge. CP%a_osc) then
       DeltaTime=rombint(dtauda,a1,CP%a_osc*(1._dl-max(atol/100.0_dl,1.d-15)),atol) + rombint(dtauda, CP%a_osc, a2, atol)
    else
       DeltaTime=rombint(dtauda,a1,a2,atol)
    end if

Rationale: ``dtauda(a)`` is (mildly) discontinuous at the KG→EFA switch; Romberg integration across the kink loses accuracy/convergence. The first integral ends at ``a_osc*(1-max(atol/100,1e-15))`` (just below the switch), the second starts exactly at ``a_osc``. **Must be ported** into modern CAMB's ``TCAMBdata%DeltaTime``/``rombint`` usage (or whatever background integrator handles dtauda) whenever an axion switch exists.

4.2 ``DeltaPhysicalTimeGyr`` split at a_osc — [PHYSICS]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``AxiECAMB/modules.f90:623-644``

.. code-block:: fortran

    !RL modified for ULA switch 043025
    if (a1 .lt. CP%a_osc .and. a2 .ge. CP%a_osc) then
       DeltaPhysicalTimeGyr = (rombint(dtda,a1,CP%a_osc*(1._dl-max(atol/100.0_dl,1.d-15)),atol) + &
            & rombint(dtda,CP%a_osc,a2,atol))*Mpc/c/Gyr
    else
       DeltaPhysicalTimeGyr = rombint(dtda,a1,a2,atol)*Mpc/c/Gyr
    end if

4.3 ``ComovingRadialDistance`` tolerance hardened — [ACCURACY]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``AxiECAMB/modules.f90:666``

.. code-block:: fortran

    ComovingRadialDistance = DeltaTime(1/(1+z),1._dl, 1.d-7) !RL reverted 062624 since thetastar may have a problem there !RL modified the tolerance here to 1.e-6 or else the code crashes for some axion masses

(old: default tolerance ``tol/1000/exp(AccuracyBoost-1)`` ≈ 1e-7 at boost 1, but explicit 1.d-7 is now forced regardless of AccuracyBoost). Physics-motivated robustness fix near the switch; port the explicit tolerance.

4.4 ``CosmomcTheta``: axions counted in ω_dm — [PHYSICS]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``AxiECAMB/modules.f90:741``

.. code-block:: fortran

    omdmh2 = (CP%omegac+CP%omegan+CP%omegaax)*(CP%h0/100.0d0)**2! RH axion

(Hu & Sugiyama z* fitting formula input; comment notes θ_MC is just a stepping parameter.) Port to modern ``CosmomcTheta()`` in results.f90.

5. Module ``Transfer`` changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

5.1 New transfer columns — [PHYSICS]/[PLUMBING]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``AxiECAMB/modules.f90:1737-1743``

.. code-block:: fortran

    integer, parameter :: Transfer_kh =1, Transfer_cdm=2,Transfer_b=3,Transfer_g=4, &
         Transfer_r=5, Transfer_nu = 6,  & !massless and massive neutrino
         Transfer_axion=7,  Transfer_f=8, & ! DM: Axions and growth rate
         Transfer_tot=9
    integer, parameter :: Transfer_max = Transfer_tot

OLDCAMB: ``Transfer_tot=7``, ``Transfer_max=7``. AxiECAMB inserts ``Transfer_axion=7`` (axion density transfer δ_ax/k²-normalized like others) and ``Transfer_f=8`` (growth rate f; filled in cmbmain — see cmbmain report) and shifts ``Transfer_tot`` to 9. ``transfer_power_var`` default still ``Transfer_tot`` (so σ₈/P(k) include axions in "total matter").

**Port note:** modern CAMB has 13 columns (``Transfer_kh..Transfer_tot=7, Transfer_nonu=8, Transfer_tot_de=9, Transfer_Weyl=10, Transfer_Newt_vel_*=11,12, Transfer_vel_baryon_cdm=13``). The axion column must be **appended** (e.g. 14) rather than inserted at 7; all places that fill ``Transfer_axion``/``Transfer_f`` and the definition of ``Transfer_tot``/``Transfer_nonu`` must be revisited: in modern CAMB, decide whether δ_ax enters ``Transfer_tot`` (it must, post-switch, for σ₈ consistency with AxiECAMB). ``Transfer_f`` duplicates modern CAMB's growth outputs (``Transfer_Newt_vel_*`` / ``MatterTransferData`` growth handling) — likely [OBSOLETE] if modern facilities suffice, but check what cmbmain stores there.

5.2 ``MatterPowerData_k``: remove −30 floor; zero high-k extrapolation — [ACCURACY] (physics-motivated)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``AxiECAMB/modules.f90:1926-1972``. Low-k branch and spline branch:

.. code-block:: fortran

    !    outpower = exp(max(-30._dl,outpower))
    outpower = exp(outpower) ! RH change for axions !RL 122024 moved location so that the second category is safer

High-k extrapolation branch (logk > log_kh(num_k)):

.. code-block:: fortran

    !Do dodgy linear extrapolation on assumption accuracy of result won't matter
    !RL 122024 this can be too dodgy for light DM masses that our results may break. Since the result won't matter we directly set it to zero
    outpower = 0._dl

Rationale: ULA P(k) is exponentially suppressed below the Jeans scale; the ``exp(max(-30,·))`` floor and log-linear extrapolation corrupt the suppressed tail. Port both changes into ``results.f90: MatterPowerData_k`` (modern code still has ``exp(max(-30._dl,outpower))`` and the dodgy extrapolation).

5.3 ``Transfer_GetMatterPowerD``: same floor removal — [ACCURACY]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``AxiECAMB/modules.f90:2078-2079``

.. code-block:: fortran

    !    outpower = exp(max(-30.d0,outpower))
    outpower = exp(outpower) ! RH change for axions

(Note: the ``outpower(il)=-30.`` assignment for k below the computed range at :2049 is retained, so those points become ``exp(-30)`` ≈ 9.4e-14 — unchanged behavior there.)

5.4 ``Transfer_SetForNonlinearLensing``: 10× more NL lensing redshifts — [ACCURACY]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``AxiECAMB/modules.f90:2219-2220``

.. code-block:: fortran

    !    P%NLL_num_redshifts =  nint(10*AccuracyBoost) ! RH original
    P%NLL_num_redshifts =  nint(100*AccuracyBoost) ! RH change for axions

Physics-motivated (axion P(k,z) shape evolves; halofit ratio interpolation in z needs finer grid), but the factor 100 is a hardcoded choice. Requires ``max_transfer_redshifts`` ≥ ~100*boost (§1.1). Port decision: keep, but consider making it conditional on ``omegaax > 0``.

5.5 ``Transfer_output_Sig8``: warning for z > z_osc — [PLUMBING] (output only)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``AxiECAMB/modules.f90:2146-2179``

.. code-block:: fortran

    real(dl) z_osc !RL 121924 for comparison
    z_osc = 1._dl/CP%a_osc - 1._dl
    if (z_osc .lt. 0._dl) then
       z_osc = 0._dl
    end if
    ...
    if (real(CP%Transfer%redshifts(j)) .gt. z_osc) then
       write(*, '(A, F9.5, A, F9.5, A)') 'Note: z = ', real(CP%Transfer%redshifts(j)), &
            &' is before the axion switch point z = ', z_osc, &
            &', hence T(k) is defined differently from that after the switch point.'
    end if

σ₈ algorithm itself unchanged (still integrates ``transfer_power_var = Transfer_tot``, which now includes axions because Transfer_tot is the 9th column filled in cmbmain). The physics of what enters Transfer_tot lives in cmbmain/equations.

5.6 ``Transfer_SaveToFiles``: 7 → 9 columns — [PLUMBING]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``AxiECAMB/modules.f90:2253``

.. code-block:: fortran

    write(fileio_unit,'(9E14.6)') MTrans%TransferData(Transfer_kh:Transfer_max,ik,i) ! axion

Output transfer file columns become: k/h, CDM, baryon, photon, massless ν, massive ν, **axion**, **growth rate f**, total. Modern CAMB writes a header + 13 columns via ``Transfer_SaveToFiles``; port = append axion column(s) with proper header label.

5.7 Timing instrumentation — [COSMETIC]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``cpu_time`` calls in ``Transfer_GetMatterPowerData``, ``MatterPowerdata_getsplines`` (prints commented out).

6. Module ``ThermoData`` changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

6.1 ``Recombination_Init`` call replaced (recfast_axion) — [PHYSICS]/[PLUMBING]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``AxiECAMB/modules.f90:2498-2503``

.. code-block:: fortran

    !axion version of recfast call
    call Recombination_Init(CP%Recomb, CP%omegac,CP%omegab,CP%Omegan,&
         CP%Omegav,CP%h0,CP%tcmb,CP%yhe,CP%omegaax,CP%omegar,CP%aeq, CP%a_osc)
    ! Data structures not passed through params because camb data structures significantly different

Old signature ended ``...,CP%yhe,CP%Num_Nu_massless + CP%Num_Nu_massive)``. New extra args: ``omegaax, omegar, aeq, a_osc``; the neutrino-count argument is **dropped** (recfast_axion builds H(z) from omegar directly). See recfast_axion.F90 report for what it does with them. In modern CAMB, ``TRecfast%Init`` receives ``State`` — the port should pass axion quantities through the state object, making this whole signature change [PLUMBING] with [PHYSICS] content inside recfast_axion.

6.2 Smoothed-xe (boxcar) machinery — **currently disabled** — [OBSOLETE]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``AxiECAMB/modules.f90:2546-2664``. A pre-pass loop integrates the Friedmann equation once to build

.. code-block:: fortran

    recdotmu_ro(i) = Recombination_xe(a)*akthom/a2
    recdotmu_ro(i) = recdotmu_ro(i)*(tau**tau_n)

then a boxcar filter of half-width ``half_w`` produces ``recxe_sm(i) = sum_sm*wt/(akthom/(scaleFactor(i)**2._dl))``, after which the original inithermo loop is repeated using ``recxe_sm`` instead of ``Recombination_xe(a)``:

.. code-block:: fortran

    xe(i) = Reionization_xe(a, tau, recxe_sm(ncount)/(tau**tau_n))   ! reionized era
    dotmu(i) = (recxe_sm(i)/(tau**tau_n) - xe(i))*akthom/a2          ! AccurateReionization branch
    xe(i) = recxe_sm(i)/(tau**tau_n)                                  ! pre-reionization
    ...
    xe(1) = recxe_sm(1)/(tauminn**tau_n)        ! modules.f90:2662
    dotmu(1) = xe(1)*akthom/a02                  ! modules.f90:2664

**BUT** the controls are hard-set to identity:

.. code-block:: fortran

    tau_n = 0 !RL 010224 ... !121524 - set tau_n = 0 to turn smoothing off since we don't need it     ! :2549
    half_w = 0 !RL 010224 !121524 - since we have recombination skip, we turn smoothing off ...        ! :2554
    wt = 1._dl/(2*half_w + 1)

With ``half_w=0, tau_n=0``, ``recxe_sm(i) ≡ Recombination_xe(a_i)`` and the result is bitwise-equivalent to stock CAMB (the original delta-function-in-ẋe problem was solved instead by the "recombination skip" in the perturbation module). **Recommendation: do not port the boxcar; treat as dead code.** Only residual effects worth noting:

- ``scaleFactor(1) = a0`` is now explicitly initialized (``modules.f90:2540``); in OLDCAMB index 1 was left unset (latent bug). Modern CAMB already initializes its scale_factor table correctly.
- the Friedmann integration is performed twice (wasted work, harmless).

6.3 Local ``vis`` renamed ``visi`` — [COSMETIC]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``modules.f90:2469`` (``real(dl) dtbdla,vfi,cf1,maxvis, visi``) — avoids shadowing the public ``vis(:)`` array so ``ThermoSplineOut`` can be validated. No behavior change.

6.4 New: ``ThermoSplineOut`` + caching of opacity/exp(−κ) at τ_osc — [PHYSICS]/[PLUMBING]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Call site after the DoThermoSpline loop, ``AxiECAMB/modules.f90:2889-2890``:

.. code-block:: fortran

    !RL 09062023
    call ThermoSplineOut(CP%tau_osc, CP%opac_tauosc, CP%expmmu_tauosc)

New subroutine ``modules.f90:3155-3224`` (active lines only):

.. code-block:: fortran

    subroutine ThermoSplineOut(tau, opac_tau, expmmu_tau)
      ...
      d=log(tau/tauminn)/dlntau+1._dl
      i=int(d)
      d=d-i
      if (i < nthermo) then
         opac_tau=dotmu(i)+d*(ddotmu(i)+d*(3._dl*(dotmu(i+1)-dotmu(i)) &
              -2._dl*ddotmu(i)-ddotmu(i+1)+d*(ddotmu(i)+ddotmu(i+1) &
              +2._dl*(dotmu(i)-dotmu(i+1)))))
         expmmu_tau=emmu(i)+d*(demmu(i)+d*(3._dl*(emmu(i+1)-emmu(i)) &
              -2._dl*demmu(i)-demmu(i+1)+d*(demmu(i)+demmu(i+1) &
              +2._dl*(emmu(i)-emmu(i+1)))))
      else
         opac_tau=dotmu(nthermo)
         expmmu_tau=emmu(nthermo)
      end if
    end subroutine ThermoSplineOut

Purpose: the perturbation module needs κ̇(τ_osc) and e^{−κ}(τ_osc) for the source/visibility handling at the switch step (used in equations_ppf/cmbmain — see those reports). **Port:** modern CAMB's ``TThermoData%Values(tau, ...)`` already provides exactly this interpolation; just call it once after ``inithermo`` and store in the axion state. No new subroutine needed → re-derive as [PLUMBING].

6.5 r_s(z*) integral split at a_osc — [PHYSICS]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``AxiECAMB/modules.f90:2896-2903``

.. code-block:: fortran

    !RL062624 This rs should be split - and 1d-6 is the atol which should be changed throughout these two lines if it needs to be changed
    if (1d-8 .le. CP%a_osc .and. 1/(z_star+1) .gt. CP%a_osc) then
       rs=rombint(dsound_da_exact,1d-8,CP%a_osc,1d-6) + &
            & rombint(dsound_da_exact, CP%a_osc*(1._dl+max(1d-6/100.0_dl,1.d-15)), 1/(z_star+1), 1d-6)
    else
       rs =rombint(dsound_da_exact,1d-8,1/(z_star+1),1d-6)
    end if

**Inconsistency to flag:** the analogous ``rdrag`` integral (``modules.f90:2913``: ``rs =rombint(dsound_da_exact,1d-8,1._dl/(z_drag+1._dl),1d-6)``) and the ``theta_rs_EQ`` integral (:2922) are NOT split at a_osc. If z_drag < z_osc the rdrag value will integrate across the kink. Decide in the port whether to split all ``dsound_da_exact``/``ddamping_da`` rombints uniformly (recommended: yes, via a helper).

6.6 z_EQ includes axions; new derived params — [PHYSICS] (the grhoax part) / [OBSOLETE] (the rest)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``AxiECAMB/modules.f90:2905-2922``

.. code-block:: fortran

    ThermoDerivedParams( derived_DAstar ) = DA/1000
    ...
    z_eq = (grhob+grhoc+grhoax)/(grhog+grhornomass+sum(grhormass(1:CP%Nu_mass_eigenstates))) -1._dl
    ThermoDerivedParams( derived_zEQ ) = z_eq
    a_eq = 1._dl/(1._dl+z_eq)
    ThermoDerivedParams( derived_kEQ ) = 1._dl/(a_eq*dtauda(a_eq))
    ThermoDerivedParams( derived_thetaEQ ) = 100._dl*timeOfz( ThermoDerivedParams( derived_zEQ ))/DA
    ThermoDerivedParams( derived_theta_rs_EQ ) = 100._dl*rombint(dsound_da_exact,1d-8,a_eq,1d-6)/DA

DAstar/kEQ/theta_rs_EQ are back-ports already present in CAMB 1.6.7 (``results.f90`` uses ``grhob+grhoc`` for zeq there). **The change to port is adding grhoax to the matter side of z_eq** — modern CAMB computes ``z_eq`` via ``1/aeq-1`` with ``aeq`` from ``grho_no_de``-style expressions; the axion must count as matter there (with the caveat the authors note: only well-defined when the axion behaves as matter, i.e. m≫H_eq). The feedback printing of z_EQ/theta_EQ is deliberately suppressed (``modules.f90:2946-2958``, DG explanation: axions can act as matter or DE depending on mass, output ambiguous) — [COSMETIC]/output policy, but worth carrying the comment.

Also the failure branch of recombination-end search now prints diagnostics (``modules.f90:2828-2831``):

.. code-block:: fortran

    call GlobalError('inithermo: failed to find end of recombination',error_reionization)
    print*, 'omaxh2, omch2, mass, H0', CP%omegaax*(CP%h0/100.)**2, CP%omegac*(CP%h0/100.)**2, CP%ma, CP%h0

[COSMETIC].
Minor format tweaks: ``r_s(zstar)`` print format ``f7.2`` → ``f7.3`` (:2940). [COSMETIC]

6.7 ``SetTimeSteps``: τ-grid refinement around τ_osc — [PHYSICS]/[ACCURACY]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``AxiECAMB/modules.f90:2967-3092``; the **active** new code (stripped of the many commented experiments):

.. code-block:: fortran

    subroutine SetTimeSteps
      real(dl) dtau0, dtauosc !RL 121323: adding dtauosc for timestep refinements around tauosc
      ...
      call Ranges_Add_delta(TimeSteps, taurst, taurend, dtaurec) !default
      ...
      call Ranges_Add_delta(TimeSteps,taurend, CP%tau0, dtau0) !default
      ...reion Ranges_Add unchanged...
    !RL 061924 Everything should happen only if we have a switch - and use a_osc from the background to definitively say whether we have a switch or not, to avoid the small but finite mis-correspondence between tau_osc and a_osc
    if (CP%a_osc .le. 1._dl) then
        if (CP%tau_osc .gt. taurst) then
             dtauosc = CP%tau_osc/int(6000._dl*CP%dfac)
             call Ranges_Add_delta(TimeSteps, max(CP%tau_osc-6.5_dl*dtauosc, taurst), min(CP%tau_osc+6.5_dl*dtauosc, CP%tau0), dtauosc)!RL
        end if
        if (CP%tau_osc .gt. taurend) then
           dtauosc = dtauosc*1000._dl
           call Ranges_Add_delta(TimeSteps, max(6.5_dl*dtauosc, taurend), min(CP%tau_osc+6.5_dl*dtauosc, CP%tau0), dtauosc)
        end if
    end if
      call Ranges_GetArray(TimeSteps)
      ...

Two refinements when a switch happens before today (``CP%a_osc <= 1``):

1. **Ultra-fine window** (only if τ_osc is after recombination start): ±6.5 steps of ``Δτ = τ_osc/int(6000·dfac)`` bracketing τ_osc (≈13 extra points, spacing ~τ_osc/60000 for dfac=10). Resolves the source discontinuity at the switch in the C_ℓ time integration.
2. **Coarse refinement** (only if τ_osc after recombination end): spacing ×1000 (i.e. ``τ_osc/int(6·dfac)``-ish) from ``max(6.5Δτ, taurend)`` up to just past τ_osc.

Subtlety: in branch 2, ``dtauosc`` relies on branch 1 having executed; safe because ``τ_osc > taurend ⇒ τ_osc > taurst``. Port into ``TThermoData%SetTimeSteps`` (results.f90) using ``TRanges%Add_delta``; same logic, state-based.
The numbers 6000, 6.5, ×1000 are hardcoded tuning ([ACCURACY], but required for correct C_ℓ; keep as-is initially and document).

6.8 ``DoThermoSpline`` — unchanged except debug comments. [COSMETIC]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

7. Items explicitly checked and NOT present / NOT changed in modules.f90
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **Initial-condition mode count / initial_iso_axion / initial vector names**: not in this file — they live in ``equations_ppf.f90`` (GaugeInterface). modules.f90 only carries the CP fields feeding them (``axion_isocurvature, alpha_ax, Hinf, r_val/ratio, amp_i``). ``CP%InitialConditionVector(1:10)`` unchanged.
- ``lSampleBoost``, ``lAccuracyBoost``, ``HighAccuracyDefault``, ``tol``, ``lmin``, ``l0max``, ``OmegaKFlat``: unchanged.
- ``initlval`` (ℓ sampling): unchanged (pure re-indent).
- ``Init_ClTransfer``, Limber machinery, ``Cl_scalar/tensor/vector`` allocation, ``output_cl_files``, ``output_lens_pot_files``, ``output_veccl_files``, ``NormalizeClsAtL``: unchanged (only commented-out full-precision debug writes, RH timing). [COSMETIC]
- ``MassiveNu`` module (Nu_init/Nu_background/Nu_rho/Nu_drho): unchanged. Only the *call site* of ``init_massive_nu`` moved (§3.5).
- ``thermo()`` interpolator, ``Thermo_OpacityToTime``, ``ddamping_da``, ``doptdepth_dz``/``optdepth``/``dragoptdepth``/``find_z``: unchanged. Note ``doptdepth_dz`` and ``ddamping_da`` use ``Recombination_xe`` and ``dtauda`` directly — with an axion switch their rombint integrands also have the dtauda kink; AxiECAMB did **not** split those (only ComovingRadialDistance/DeltaTime/DeltaPhysicalTimeGyr/rs(z*)). Same treatment question as §6.5.
- ``Transfer_Get_sigma8`` algorithm: unchanged (it inherits axion effects through ``transfer_power_var``).
- ``dsound_da_exact`` / ``dsound_da``: formulas unchanged (baryon-photon only, correct).

8. Summary of port classifications
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**[PHYSICS] (must port):**

- CP fields: omegaax, ma, m_ovH0, dfac, a_osc, tau_osc, ah_osc, ahosc_ETA, A_coeff, tvarphi_{c,cp,s,sp}, wEFA_c, phiinit, aeq, ainit, rhorefp_ovh2, Prefp, alpha_ax, Hinf, axfrac/omegada/use_axfrac, axion_isocurvature, a_skip/dfac_skip/a_skipst (§2)
- Background tables loga_table/phinorm/phidot/rhoaxh2ovrhom (+ntable, nphi, aeq_LCDM, drefp_hsq, grhoax) — as state members (§1.2-1.4)
- grhoax=grhom*omegaax + a_osc≤1 clamp (§3.4)
- DeltaTime / DeltaPhysicalTimeGyr / rs(z*) rombint splits at a_osc (§4.1, 4.2, 6.5)
- CosmomcTheta ω_dm includes ω_ax (§4.4)
- z_eq includes grhoax (§6.6)
- Recombination_Init axion arguments (content in recfast_axion) (§6.1)
- Reionization_Init gains a_osc (§3.7)
- ThermoSplineOut caching of opac/expmmu at τ_osc (re-implement via TThermoData%Values) (§6.4)
- SetTimeSteps τ_osc refinement windows (§6.7)
- Transfer_axion (and decide on Transfer_f) columns (§5.1, 5.6)

**[ACCURACY] (physics-motivated; port with care):**

- exp(max(-30,...)) floor removal + zero high-k extrapolation in MatterPowerData_k / Transfer_GetMatterPowerD (§5.2, 5.3)
- NLL_num_redshifts 10→100 ×boost + max_transfer_redshifts 600 (§5.4, §1.1)
- ComovingRadialDistance explicit 1e-7 tolerance (§4.3)

**[PLUMBING] (re-derive, don't copy):**

- Neutrino bookkeeping & init_massive_nu relocation to driver; CP%Nu_massless_degeneracy (§3.1, 3.5) — modern CAMB already handles; just preserve ordering vs axion background solve
- GetOmegak removal / omegak+H0 closure in driver (§3.2) — fold Ω_ax into modern budget logic
- tau0 computed in init_background (§3.6)
- CP%omegar/CP%grhor/omegah2_rad/H0_in_Mpc_inv/H0_eV/amp_i/ratio/lens_amp duplicated derived quantities (§2)
- 9-column transfer file output (§5.6), σ₈ z>z_osc warning (§5.5), feedback prints (§3.8)

**[OBSOLETE]:**

- 13 derived parameters (DAstar/keq/theta_rs_EQ) — already in CAMB 1.6.7 (§1.6, 6.6)
- Boxcar xe smoothing machinery — disabled (half_w=0, tau_n=0), superseded by recombination skip (§6.2)
- scaleFactor(1) init fix — already correct in modern CAMB (§6.2)
- vis→visi rename, timing code, commented experiments, full re-indent — [COSMETIC]

**Risks / surprises:**

1. ``CAMBParams_Set`` no longer self-consistent: it *requires* the driver to pre-set ``Nu_massless_degeneracy``, ``omegak``, ``omegar``, ``grhor``, and to have run the axion background solver — the call-order contract must be re-established deliberately in the class-based port.
2. ``Transfer_axion=7`` insertion renumbered ``Transfer_tot`` (7→9); any literal index assumptions elsewhere (halofit, sigma8.f90, COSMOSIS interfaces) must be audited; in modern CAMB append instead.
3. rdrag/optdepth/damping integrals are NOT split at a_osc while rs(z*) is — latent inconsistency to resolve (recommend uniform splitting helper).
4. The boxcar smoothing looks like real physics but is a disabled no-op — do not port.
5. Module-level mutable globals (``a_osc``, ``tau_osc``, ``grhoax``, tables, ``ntable``) are thread/state hazards in modern CAMB; everything must move into the state/params objects.


Original-code analysis: perturbation equations (equations_ppf.f90)
--------------------------------------------------------------------

Files:

- AxiECAMB: ``/Users/vivianmiranda/data/research/WayneHu/rayne/AxiECAMB/equations_ppf.f90`` (4175 lines)
- Base:     ``/Users/vivianmiranda/data/research/WayneHu/rayne/OLDCAMB/equations_ppf.f90`` (2920 lines)
- Diff:     ``/Users/vivianmiranda/data/research/WayneHu/rayne/.port_analysis/diffs/equations_ppf.f90.diff``

All line numbers below refer to the **AxiECAMB** file unless prefixed ``OLD:``. Most of the 5414-line
diff is re-indentation (4-space → 2-space, ``end if`` alignment) — pure [COSMETIC] noise; everything
substantive is listed here. The file also carries a very large amount of commented-out debugging
(``!RL ...``, ``!write(...)``) which is [COSMETIC] and is not itemized except where the comments document
physics choices.

0. Conventions, units and external state used by this file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The axion machinery leans on quantities computed in ``axion_background.F90`` and stored in
``ModelParams`` (``modules.f90``, see modules report). Used here:

- ``loga_table(ntable)``, ``phinorm_table``, ``phidotnorm_table`` (+ ``_ddlga`` spline buffers),
  ``rhoaxh2ovrhom_logtable`` (+ ``_buff``): background tables vs ``log10(a)``. The tables end at/near
  ``a_osc`` (the background grid endpoint is moved to ``a_osc`` during the shooting in
  ``axion_background.F90``).
- Field normalization (comment at 3383): ``v1_bg = sqrt(4*pi*G/3)*phibar`` (dimensionless) and
  ``v2_bg = sqrt(4*pi*G/3)*phibar_dot/H0``; same normalization for the perturbations ``dv1``, ``dv2``
  (with an extra numerical rescale ``EV%renorm_c`` on ``dv1``, see §5).
  ``grhoax_kg = v2_bg**2/a2 + (CP%m_ovH0*v1_bg)**2`` equals ``rho_ax/rho_crit,0 * h^2``-like fraction
  (it is divided by ``CP%H0**2/1e4 = h^2`` to get the density fraction).
- ``CP%m_ovH0`` = m/H0 (dimensionless); ``CP%H0_in_Mpc_inv`` = H0 in Mpc^-1; ``CP%H0_eV`` = H0 in eV;
  ``CP%ma`` = axion mass in eV. ``m`` in Mpc^-1 is always ``CP%m_ovH0*CP%H0_in_Mpc_inv``.
- Switch-epoch background quantities (set in ``axion_background.F90``): ``CP%a_osc``, ``CP%tau_osc``,
  ``CP%ah_osc`` and ``CP%ahosc_ETA`` (= aH at the switch, in units of 100 km/s/Mpc; ``ahosc_ETA`` is the
  "ETA" variant adopted 082924), ``CP%rhorefp_ovh2`` (axion density fraction at the reference point
  a_osc), ``CP%Prefp`` (axion pressure there), ``CP%wEFA_c`` (coefficient of the residual-w correction),
  ``CP%A_coeff``, ``CP%tvarphi_c, tvarphi_s, tvarphi_cp, tvarphi_sp`` (background cos/sin WKB
  decomposition of the field at the switch, used for the perturbation matching), ``CP%aeq``,
  ``CP%opac_tauosc``/``CP%expmmu_tauosc`` (thermo at the switch; only used in commented variants here).
- ``EV%metric_delta(2)`` is exported to cmbmain.f90:1220 (``deltaBCSrc(EV%q_ix,1,:)=EV%metric_delta``)
  where it becomes a delta-function boundary term in the line-of-sight source integral (see cmbmain
  report).

1. Background: ``dtauda``, ``grhoax_frac``, ``init_background``, ``GetOmegak``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1.1 ``GetOmegak`` deleted (lines 211–227, commented out) — [PLUMBING]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The base-code function ``GetOmegak() = 1 - (omegab+omegac+omegav+omegan)`` is entirely commented out.
RL comment (line 218): *"this is physically incorrect and we should modify H0 instead. The curvature
should always be determined at the input level."* The DG version (visible in the comments) had added
``CP%omegaax`` and a radiation/massless-neutrino term. In modern CAMB curvature closure is handled in
``CAMBparams``; nothing to port literally, but the *intent* (axion density participates in the budget,
curvature fixed at input) must be honored in the new ``SetParams`` logic.

1.2 ``init_background`` (lines 229–252) — [PHYSICS]/[PLUMBING]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: fortran

    229 subroutine init_background
    230   use LambdaGeneral
    231   use ModelParams !RL
    ...
    235   is_cosmological_constant = .not. use_tabulated_w .and. w_lam==-1_dl .and. wa_ppf==0._dl
    236
    237   CP%tau0=TimeOfz(0._dl) !RL 061924 moved tau0 here so to initialte tau_osc
    238
    239   if (CP%a_osc .gt. 1._dl) then
    240      CP%tau_osc=CP%tau0 + 1._dl !To make sure no switch happens before the present day.
    241   else
    242      CP%tau_osc=DeltaTime(0._dl, CP%a_osc, in_tol = 1.0d-8)
    243   end if
    244   !RL 062624 - add scenarios  potential problems that the background switches but the perturbation does not
    245   if (CP%a_osc .le. 1._dl .and. CP%tau_osc .gt. CP%tau0) then
    246      write(*, *) 'a_osc <= 1 and tau_osc > tau0. ...'
    249      CP%tau_osc = min(CP%tau_osc, CP%tau0)
    250   end if
    252 end  subroutine init_background

- Computes ``tau0`` *early* (needed because ``tau_osc`` is used in the perturbation switch and in
  ``EV%renorm_c``), then ``tau_osc = DeltaTime(0, a_osc, tol=1e-8)`` (tight tolerance — [ACCURACY],
  physics-motivated: the conformal time of the switch must be consistent between background and
  perturbations to ~1e-8 or the matching produces glitches).
- Guard: if ``a_osc > 1``, ``tau_osc = tau0 + 1`` (switch never happens); if rounding makes
  ``tau_osc > tau0`` while ``a_osc <= 1``, clamp to ``tau0``. [PHYSICS] (edge-case policy must be ported).
- Note ``DeltaTime`` in AxiECAMB gained an optional ``in_tol`` argument (modules.f90 change).

1.3 NEW function ``grhoax_frac(a_in)`` (lines 257–296) — [PHYSICS]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns ``rho_ax(a)/rho_crit,0`` (so that ``8πG a^4 ρ_ax = grhoax_frac(a)*grhom*a^4``):

.. code-block:: fortran

    268   a = a_in
    269   a_min = 10._dl**(loga_table(1))
    271   a2 = a**2.0d0
    273   if (a .lt. CP%a_osc) then
    274       if (a .lt. a_min) then
    275          grhoaxh2_ov_grhom = rhoaxh2ovrhom_logtable(1)
    276       else
    277          !Note that in the background the spline table is log10(rho)
    278          call spline_out(loga_table,rhoaxh2ovrhom_logtable,rhoaxh2ovrhom_logtable_buff,ntable,dlog10(a),grhoaxh2_ov_grhom)
    279       end if
    281      grhoax_frac = (10._dl**grhoaxh2_ov_grhom)/(CP%H0**2.0d0/1.0d4)
    283   else
    ...
    290      wcorr_coeff = CP%ahosc_ETA*CP%a_osc/((CP%ma/CP%H0_eV)*(CP%H0/100.0d0)) !RL082924
    292      grhoax_frac=(CP%rhorefp_ovh2)*((CP%a_osc/a)**3.0d0)*dexp((wcorr_coeff**2.0d0)*3.0d0*CP%wEFA_c*(1.0d0/(a2**2.0d0) &
    293           &- 1.0d0/(CP%a_osc**4.0d0))/4.0d0)
    294   endif

Physics:

- **Before** ``a_osc``: cubic-spline lookup of ``log10(ρ_ax h²/ρ_crit,0)`` vs ``log10 a`` from the exact
  background KG solution; for ``a < a_min`` the table is frozen at its first entry (the field is frozen,
  w = −1, constant ρ). The ``/(H0²/1e4) = /h²`` converts the table (which stores Ω h² units) to a
  density fraction.
- **After** ``a_osc``: analytic dilution ``ρ_ref (a_osc/a)³`` **times a residual-w correction**
  ``exp[3·wEFA_c·wcorr²·(a⁻⁴ − a_osc⁻⁴)/4]``. This is exactly
  ``exp[3∫ w dln a]`` for ``w(a) = wEFA_c · (wcorr_coeff)²/a⁴``, i.e. the time-averaged
  ``⟨w⟩ ∝ (H/m)²`` correction of arXiv:2412.15192 evaluated with ``ℋ_osc`` and a radiation-era a⁻⁴
  scaling. ``wcorr_coeff = ℋ(a_osc)·a_osc/(m/H0 · h)`` (dimensionless). Note this branch uses
  ``CP%ma/CP%H0_eV`` while every other site uses ``CP%m_ovH0`` — the same number computed two ways.

1.4 ``dtauda`` (lines 299–335) — [PHYSICS]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Single physics line added (the rest identical to OLD):

.. code-block:: fortran

    331   grhoa2 = grhoa2 + grhoax_frac(a)*grhom*(a2**2._dl)
    332   dtauda=sqrt(3._dl/grhoa2)

The axion enters as ``8πG ρ_ax a⁴ = grhoax_frac(a)·grhom·a⁴``. Declarations add
``grhoax_frac`` as ``external``. (``v1_test, v2_test, grhotest, i`` declared but unused — [COSMETIC].)

2. State vector & index bookkeeping
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

2.1 New initial-condition flag (lines 369–375) — [PHYSICS]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: fortran

    369   integer, parameter :: initial_adiabatic=1, initial_iso_CDM=2, &
    370        initial_iso_baryon=3,  initial_iso_neutrino=4, initial_iso_neutrino_vel=5, initial_vector = 0, &
    374        initial_iso_axion=6
    375   integer, parameter :: initial_nummodes =  initial_iso_axion !DM: added axion isocurvature

2.2 ``EvolutionVars`` additions (lines 379–450) — [PLUMBING] (+1 [OBSOLETE])
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: fortran

    390      integer a_ix !Index of the two axion fluid equations
    391      integer a_kg_ix !RL: Index of the two axion field equations - and will change to fluid equations after the switch
    442      real(dl) dgrhoec_ppf, dgqec_ppf, vTc_ppf
    447      logical oscillation_started !RL: whether the axion oscillation has started yet
    448      logical output_done !RL testing (...)
    449      real(dl) metric_delta(2) !RL 090323 adding boundary condition delta function
    450      real(dl) renorm_c !RL 050724

- ``a_ix``/``a_ix+1``: **legacy axionCAMB GDM fluid pair (clxax, v_ax) — now dead** (derivatives forced
  to 0, see §4.5); kept only because the isocurvature IC vector still writes them. [OBSOLETE]
- ``a_kg_ix``/``a_kg_ix+1``: the live pair. Before the switch: (δφ/renorm_c, δφ̇-normalized). After:
  (δ_ax, u_ax = (1+w)v_ax). [PLUMBING — in CAMB 1.6.7 this becomes the dark-matter/axion component's
  ``w_ix``-style index from ``TDarkEnergyModel``/custom component.]
- ``oscillation_started``: per-k regime flag; ``metric_delta(2)``: source-jump record (→ cmbmain);
  ``renorm_c``: per-k δφ rescale. [PLUMBING/PHYSICS]
- ``dgrhoec_ppf, dgqec_ppf, vTc_ppf``: PPF rest-frame variables, **computed in derivs (3337–3339) but
  never read anywhere** in AxiECAMB. Vestigial (from equations_cross). [OBSOLETE]
- Also ``max_l_evolve`` reduced 1024→512 (line 367) — memory/OpenMP-stack tweak, with the error
  message at 1318 extended "…may also need to increase stack size". [ACCURACY]/[PLUMBING]
  (hardcoded setting, not physics).

2.3 ``SetupScalarArrayIndices`` (lines 804–907) — [PLUMBING]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After massive neutrinos (i.e. axion equations are always the **last 4** slots):

.. code-block:: fortran

    891     !Axions - adding on to the DE parameters ...
    892     EV%a_ix=neq+1
    893     EV%a_kg_ix = neq + 3 !RL adding the new KG equation solvers
    895     neq = neq + 4 !RL
    897     maxeq = maxeq + 4 !RL

Always 4 equations added (2 dead GDM + 2 live KG/EFA) in **both** regimes — the equation *count*
does not change at the switch, only the meaning of ``y(a_kg_ix:a_kg_ix+1)``. (``neq_dummynu`` declared
for debugging — cosmetic.)

2.4 ``GetNumEqns`` (lines 1218–1396)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- New per-k initialization (lines 1270–1271) — [PLUMBING]:

  .. code-block:: fortran

      1270    EV%oscillation_started = .false. !RL
      1271    EV%metric_delta(:) = 0._dl !RL 090323

- **lmaxnr change at low k** (replaces OLD ``EV%lmaxnr=max(3,nint(min(7,nint(sqrt(scal)*150*EV%q))*lAccuracyBoost))``):

  .. code-block:: fortran

      1349          EV%lmaxnr=max(3,nint(min(8,nint(sqrt(scal)* 450 * EV%q))*lAccuracyBoost)) !RL modifying WH smoother lmaxnr
      1351          if (EV%lmaxnr < EV%lmaxnu) then
      1352                ! Nov 2020 change following Pavel Motloch report (RL added from newest CAMB)
      1353             EV%lmaxnr = EV%lmaxnu
      1355          endif

  [ACCURACY] — physics-motivated: more neutrino multipoles at low k (450 vs 150 coefficient, cap 8
  vs 7) to keep the large-scale ISW smooth across the mid-evolution switch ("WH smoother lmaxnr").
  The Motloch clamp already exists in CAMB ≥1.3; the 450-coefficient boost is an AxiECAMB tuning that
  should be re-validated (it is k-resolution/large-scale-Cl motivated by the axion switch).

3. The mid-evolution switch at τ_osc
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

3.1 Scheduling in ``GaugeInterface_EvolveScal`` (lines 512–718) — [PHYSICS]/[PLUMBING]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

New switch time, checked first:

.. code-block:: fortran

    523     real(dl) tau_switch_oscillation !RL adding axion oscillation switch time
    ...
    530     !smallTime =  min(tau, 1/EV%k_buf)/100 !The original smallTime - ...
    531     smallTime = 0._dl !RL killing smallTime
    ...
    542     !RL: axion oscillation switches - first it's noSwitch
    543     tau_switch_oscillation = noSwitch
    549     if (.not. EV%oscillation_started) then
    550        tau_switch_oscillation = CP%tau_osc
    551     end if
    ...
    580     next_switch = min(tau_switch_oscillation, tau_switch_ktau, tau_switch_nu_massless,EV%TightSwitchoffTime, &
    581          tau_switch_nu_massive, tau_switch_no_nu_multpoles, tau_switch_no_phot_multpoles, &
    582          tau_switch_nu_nonrel, noSwitch)

and the switch action (the evolver is integrated **exactly to τ_osc**, variables remapped, and dverk
restarted with ``ind = 1``):

.. code-block:: fortran

    584     if (next_switch < tauend) then
    585        if (next_switch > tau+smallTime) then
    586           call GaugeInterface_ScalEv(EV, y, tau,next_switch,tol1,ind,c,w)
    ...
    590        EVout=EV
    592        if (next_switch == tau_switch_oscillation) then
    593           EVout%oscillation_started = .true.
    594           call SetupScalarArrayIndices(EVout)
    597           call CopyScalarVariableArray(y,yout, EV, EVout)
    598           EV = EVout
    599           y = yout
    601           ind = 1
    603        else if (next_switch == EV%TightSwitchoffTime) then
    ...                                    ! (all other switches unchanged from OLD)

- The evolver does **not** stop/restart at the ODE-driver level; the standard recursive
  ``GaugeInterface_EvolveScal`` switch mechanism is reused; ``ind=1`` makes ``dverk`` re-select its step
  size from scratch after the discontinuous variable change. [PLUMBING — in CAMB 1.6.7 this maps to
  ``EvolveScal``'s ``next_switch`` chain in ``gauge_inv``/``equations.f90``.]
- ``smallTime = 0`` (line 531) replaces ``min(tau,1/k)/100``: guarantees the switch is taken even when
  the integrator is already sitting at τ_osc (no "skip switch if too close" logic). [ACCURACY]
  (deliberate — port it, otherwise per-k jitter in the switch time).
- Unused additions: ``use ModelData``, ``sources(CTransScal%NumSources)``, ``yprimetest`` declarations
  (debug leftovers) — [COSMETIC].

3.2 Variable remapping in ``CopyScalarVariableArray`` (lines 909–1149) — [PHYSICS] (core of the port)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

New declarations (914–924) include ``yprime, yprimeout`` work arrays and the matching temporaries.
At the top (1000-area): ``a=y(1)``, ``a2=a*a``, ``k=EV%k_buf``, ``k2=EV%k2_buf`` now extracted (needed for
the matching). Three ``if (.not. EV%oscillation_started .and. EVout%oscillation_started) then / end if``
**empty** blocks were inserted in the photon/neutrino/nu_pert copy sections (969, 989, 1010–1011) —
dead code, [COSMETIC] (one of them, 1010–1011, is mis-indented inside the ``has_nu_relativistic``
block — harmless but confusing).

The axion copy (verbatim, 1017–1146):

.. code-block:: fortran

    1017        ! Axions
    1018    yout(EVout%a_ix)=y(EV%a_ix)
    1019    yout(EVout%a_ix+1)=y(EV%a_ix+1)
    1020    ! RL adding KG
    1021    if (.not. EV%oscillation_started .and. EVout%oscillation_started) then !switch ...
    1023       call spline_out(loga_table,phinorm_table,phinorm_table_ddlga,ntable,dlog10(a),v1_bg)
    1024       call spline_out(loga_table,phidotnorm_table,phidotnorm_table_ddlga,ntable,dlog10(a),v2_bg)
    1028       drhoax_kg = (v2_bg*y(EV%a_kg_ix+1)/a2 + (CP%m_ovH0**2.0d0)*v1_bg*y(EV%a_kg_ix)*EV%renorm_c)*2.0d0 !RL 050324
    1029       grhoax_kg = (v2_bg)**2.0d0/a2+(CP%m_ovH0*v1_bg)**2.0d0
    1030       !For EFA, we will need to know hLdot, i.e. 2*k*z here, so some extra variables will be computed again
    1031       call derivs(EV,EV%ScalEqsToPropagate,CP%tau_osc,y,yprime)
    1035       wcorr_coeff = CP%ahosc_ETA*CP%a_osc/(CP%m_ovH0*(CP%H0/100.0d0)) !RL082924
    1037       w_ax = ((wcorr_coeff/a2)**2.0d0)*CP%wEFA_c
    1039       !Compute the U coefficient for constructing the pert EFA variables
    1040       tU = ((-CP%tvarphi_c + CP%tvarphi_sp)*2.0d0*yprime(3) &
    1041            &- 2.0d0*k2*y(EV%a_kg_ix+1)/(a2*(CP%m_ovH0**2.0d0)*CP%H0_in_Mpc_inv))/(a*CP%m_ovH0*CP%H0_in_Mpc_inv)&
    1042            &+ 6.0d0*CP%ah_osc*y(EV%a_kg_ix)*EV%renorm_c/(a*CP%m_ovH0*CP%H0/100.0d0)
    1044       !V coefficient
    1045       tV = (-(CP%tvarphi_s + CP%tvarphi_cp)*2.0d0*yprime(3) &
    1046            &+ 2.0d0*k2*y(EV%a_kg_ix)*EV%renorm_c/(a*CP%m_ovH0*CP%H0_in_Mpc_inv))/(a*CP%m_ovH0*CP%H0_in_Mpc_inv) &
    1047            &+ 6.0d0*CP%ah_osc*y(EV%a_kg_ix+1)/(a2*(CP%m_ovH0**2.0d0)*CP%H0/100.0d0)
    1049       !W coefficient
    1050       tW = CP%A_coeff**2.0d0 + 3.0d0*CP%A_coeff*CP%ah_osc/(a*CP%m_ovH0*CP%H0/100.0d0) &
    1051            &+ 2.0d0*k2/((a*CP%m_ovH0*CP%H0_in_Mpc_inv)**2.0d0) + 4.0d0
    1055       !The pert boundary conditions
    1056        tdvarphi_c = y(EV%a_kg_ix)*EV%renorm_c
    1057        tdvarphi_cp = (-2.0d0*tU - (CP%A_coeff + 3.0d0*CP%ah_osc/(a*CP%m_ovH0*CP%H0/100.0d0))*tV)/(2.0d0*tW)
    1058        tdvarphi_s = y(EV%a_kg_ix + 1)/(a*CP%m_ovH0) - tdvarphi_cp
    1059        tdvarphi_sp = (CP%A_coeff*tU - (2.0d0 + k2/((a*CP%m_ovH0*CP%H0_in_Mpc_inv)**2.0d0))*tV)/(2.0d0*tW)
    1065       !Normalized deltarho_ef
    1066        tdrho_ef = (CP%m_ovH0**2.0d0)*(CP%tvarphi_s*tdvarphi_cp - CP%tvarphi_c*tdvarphi_sp &
    1067             &+ CP%tvarphi_cp*tdvarphi_cp + CP%tvarphi_sp*tdvarphi_sp &
    1068             &+ tdvarphi_s*(2.0d0*CP%tvarphi_s + CP%tvarphi_cp) &
    1069             &+ tdvarphi_c*(2.0d0*CP%tvarphi_c - CP%tvarphi_sp))
    1071       tdP_ef_test = tdrho_ef - 2._dl*(CP%m_ovH0**2.0d0)*(tdvarphi_s*CP%tvarphi_s + tdvarphi_c*CP%tvarphi_c)
    1073       kamnorm_test = k2/(a2*((CP%m_ovH0*CP%H0_in_Mpc_inv)**2.0_dl))
    1074      if (kamnorm_test .lt. 1.e-14_dl) then
    1076          csquared_ax_test = kamnorm_test/4.0_dl + 5.0_dl*((1/(a*dtauda(a)))**2.0_dl)/(4.0_dl*(k2/kamnorm_test))
    1078    else
    1079       csquared_ax_test = (sqrt(1.0_dl + kamnorm_test) - 1.0_dl)**2.0_dl/(kamnorm_test) &
    1080            &+ 5.0_dl*((1/(a*dtauda(a)))**2.0_dl)/(4.0_dl*(k2/kamnorm_test))
    1082      end if
    1084       !Normalized u_ax_ef = (1+w)thetaax_ef/k. ...
    1086      u_ax_ef = k*CP%m_ovH0*(tdvarphi_c*(CP%tvarphi_s + CP%tvarphi_cp) &
    1087           &+ tdvarphi_s*(-CP%tvarphi_c + CP%tvarphi_sp))/(a*CP%H0_in_Mpc_inv)
    1090       u_ax_ef = u_ax_ef/(CP%rhorefp_ovh2*(CP%H0**2.0d0/1.0d4))
    1093       weight = (k/(CP%ahosc_ETA*CP%H0_in_Mpc_inv/(CP%H0/100.0d0)))**2._dl/&
    1094            &(3._dl + (k/(CP%ahosc_ETA*CP%H0_in_Mpc_inv/(CP%H0/100.0d0)))**2._dl) !RL 082024
    1097       u_ax_efa = u_ax_ef*(1._dl + w_ax)/(1._dl + (CP%Prefp/(CP%rhorefp_ovh2*(CP%H0**2.0d0/1.0d4))))
    1100       !Now the LHS of the two EoMs are assigned clxax_kg and u_ax_kg
    1102       yout(EVout%a_kg_ix) = tdrho_ef/(CP%rhorefp_ovh2*(CP%H0**2.0d0/1.0d4)) !deltaax_ef
    1103       !RL 103023: deltaax_efa needs to change in order to keep v_ef = v_efa and sigma_ef = sigma_efa
    1104       yout(EVout%a_kg_ix) = yout(EVout%a_kg_ix) + (3._dl*CP%ah_osc*CP%H0_in_Mpc_inv/(CP%H0/100.0d0))*(u_ax_ef-u_ax_efa)/k
    1106       yout(EVout%a_kg_ix+1) = u_ax_efa  !uax_ef, taken to be continuous, i.e. u_ax_efa
    1108       !Now we have finished switching to EFA for EVout, call derivs again to get variables on the EFA side
    1110       EVout%output_done = .false.
    1111       call derivs(EVout,EVout%ScalEqsToPropagate,CP%tau_osc,yout,yprimeout)
    1113       !We have y and yout variables. Construct the corresponding source boundary values
    1114       if (CP%flat) then
    1118          EVout%metric_delta(1)= ((-yprime(3)/k+3._dl*yprime(2)/k2)*(2._dl*(yprime(1)/a))/k &
    1119               &- y(2)/k) - ((-yprimeout(3)/k+3._dl*yprimeout(2)/k2)*(2._dl*(yprimeout(1)/a))/k - yout(2)/k)
    1121          EVout%metric_delta(2)= (-yprime(3)/k+3._dl*yprime(2)/k2)/k - &
    1122               &(-yprimeout(3)/k+3._dl*yprimeout(2)/k2)/k
    1124       else
    1126          EVout%metric_delta(1)= ((-yprime(3)/k+3._dl*(yprime(2)-CP%curv*(-yprime(3)/k))/k2)*&
    1127               &(2._dl*(yprime(1)/a))/(k*EV%Kf(1)) - y(2)/(k*EV%Kf(1))) -&
    1128               & ((-yprimeout(3)/k+3._dl*(yprimeout(2) - CP%curv*(-yprimeout(3)/k))/k2)*&
    1129               &(2._dl*(yprimeout(1)/a))/(k*EVout%Kf(1)) - yout(2)/(k*EV%Kf(1)))
    1131          EVout%metric_delta(2)= (-yprime(3)/k+3._dl*(yprime(2)-CP%curv*(-yprime(3)/k))/k2)/(k*EV%Kf(1)) - &
    1132               &(-yprimeout(3)/k+3._dl*(yprimeout(2) - CP%curv*(-yprimeout(3)/k))/k2)/(k*EV%Kf(1))
    1134       end if
    1136       yout(2) = (EVout%metric_delta(2)*yprimeout(1)/a)*(k*EV%Kf(1))*weight + y(2) !yout(2) is etaTEFAnew*k
    1138        EVout%metric_delta(1) = EVout%metric_delta(1) - weight*(yprimeout(1)/a)*EVout%metric_delta(2)
    1140        EVout%metric_delta(2) = EVout%metric_delta(2)*(1._dl -weight)
    1143    else
    1144       yout(EVout%a_kg_ix) = y(EV%a_kg_ix)
    1145       yout(EVout%a_kg_ix+1) = y(EV%a_kg_ix+1)
    1146    end if

Reading guide (all [PHYSICS], must be ported exactly):

- The KG solution at τ_osc is decomposed onto the WKB basis
  φ = a^{-3/2}[φ_c cos(m t) + φ_s sin(m t)] using the *background* coefficients
  ``CP%tvarphi_{c,s,cp,sp}`` and ``CP%A_coeff`` fixed by ``axion_background.F90``. ``tU,tV,tW`` and the
  resulting ``tdvarphi_{c,cp,s,sp}`` solve the 2×2 linear system for the perturbed cos/sin amplitudes
  given (δφ, δφ̇) and the metric source ḣ/2 = −yprime(3) (note ``yprime(3) = clxcdot = −k z``, so
  ``2*yprime(3)`` = −2kz = −ḣ_L; the ``CP%ah_osc`` terms are the 3ℋ/(am) pieces).
- ``tdrho_ef``, ``tdP_ef_test`` are the cycle-averaged effective-fluid δρ and δP built from those
  amplitudes; ``u_ax_ef`` the heat flux (1+w)θ/k. They are normalized by the reference density
  ``CP%rhorefp_ovh2*h²``.
- ``u_ax_efa`` rescales the flux by ``(1+w_EFA)/(1+P_ref/ρ_ref)`` so the *velocity* v (not u) is
  continuous given the slightly different w on the two sides.
- ``yout(a_kg_ix)`` gets an extra ``3ℋ(u_ef−u_efa)/k`` so that the **rest-frame density/σ are
  continuous** ("keep v_ef = v_efa and sigma_ef = sigma_efa").
- ``derivs`` is called on both sides at exactly τ_osc; ``(-yprime(3)/k + 3*yprime(2)/k2)`` is exactly
  σ (shear) reconstructed from clxcdot and etakdot (flat case; curvature version in the else),
  ``yprime(1)/a = ℋ``. So
  ``metric_delta(1) = [2ℋσ/k − η k/k]_KG − [...]_EFA`` and ``metric_delta(2) = [σ/k]_KG − [σ/k]_EFA``:
  the jumps in the two source-term combinations that multiply j_l and j_l′ in the line-of-sight
  integral.
- ``weight = (k/ℋ_osc)²/(3+(k/ℋ_osc)²)`` splits the jump: for sub-horizon k (weight→1) the η variable
  itself is shifted (``yout(2) = ... + y(2)``, line 1136) and the recorded deltas are reduced
  accordingly (1138, 1140); for super-horizon k (weight→0) η stays continuous and the full jump is
  exported via ``EV%metric_delta`` → ``deltaBCSrc`` in cmbmain (boundary terms from integrating the
  delta-function in the sources by parts).
- ``csquared_ax_test``/``tdP_ef_test``/``kamnorm_test`` at 1071–1082 are *diagnostics only* (used only in
  commented writes) — do not port. [COSMETIC]
- For all later (non-oscillation) switches the live pair is copied verbatim (1144–1145).

**Porting note**: this block is the single hardest piece. In CAMB 1.6.7 there is no generic
mid-evolution variable-meaning change; you must add an ``oscillation_started``-like flag to the
EvolutionVars analog, register τ_osc in the ``next_switch`` chain, and re-derive the η-shift/source
delta plumbing against the 1.6.7 source functions (which are assembled differently —
``output`` now computes sources via ``TSourceTermParams``-style code).

4. Perturbation equations in ``derivs`` (lines 3073–3702)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

4.1 New declarations (3106–3113) — [PLUMBING]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: fortran

    3106     ! Axions
    3107     real(dl) :: v1_bg, v2_bg, dv1, dv2, drhoax_kg, grhoax_kg, clxax_kg, clxax_kg_dot, u_ax_kg, u_ax_kg_dot
    3108     real(dl) grhoax_t, clxax, gpres_ax
    3109     real(dl) clxaxdot,v_ax,v_axdot
    3110     real(dl) w_ax, w_ax_p1, wcorr_coeff, csquared_ax,cad2
    3111     real(dl) dorp
    3112     real(dl) gr,kamnorm

4.2 Axion state extraction / background at current a (3130–3174) — [PHYSICS]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: fortran

    3130     ! Axion variables
    3131     clxax=ay(EV%a_ix)
    3132     v_ax=ay(EV%a_ix+1)
    3135     if (.not. EV%oscillation_started) then
    3136        dv1 = ay(EV%a_kg_ix) !RL 050324: not renormalizing here ... renormalize the delta rho etc.
    3137        dv2 = ay(EV%a_kg_ix+1)
    3139        call spline_out(loga_table,phinorm_table,phinorm_table_ddlga,ntable,dlog10(a),v1_bg)
    3140        call spline_out(loga_table,phidotnorm_table,phidotnorm_table_ddlga,ntable,dlog10(a),v2_bg)
    3142        w_ax_p1 = (2.0_dl*(v2_bg**2.0_dl)/a2)/((v2_bg**2.0_dl)/a2 + (CP%m_ovH0*v1_bg)**2.0_dl)
    3143        w_ax = w_ax_p1 - 1.0_dl
    3145        drhoax_kg = (v2_bg*dv2/a2 + (CP%m_ovH0**2.0d0)*v1_bg*dv1*EV%renorm_c)*2.0d0 !RL 050324
    3147        grhoax_kg = (v2_bg)**2.0d0/a2+(CP%m_ovH0*v1_bg)**2.0d0
    3148        if (v1_bg .eq. v1_bg .and. v2_bg .eq. v2_bg) then !RL inherited DG flag   [NaN guard]
    3149           dorp = grhom*grhoax_kg/(CP%H0**2.0d0/1.0d4)
    3150        else
    3151           dorp=0.0d0
    3152        end if
    3153        clxax_kg = drhoax_kg/grhoax_kg
    3154        u_ax_kg = w_ax_p1*k*dv1*EV%renorm_c/(CP%H0_in_Mpc_inv*v2_bg) !RL 050324
    3155     else !past tauosc
    3157        v1_bg = 0._dl ; v2_bg = 0._dl ; drhoax_kg = 0._dl ; grhoax_kg = 0._dl     [3157–3160]
    3163        wcorr_coeff = CP%ahosc_ETA*CP%a_osc/(CP%m_ovH0*(CP%H0/100.0d0)) !RL082924
    3164        w_ax = ((wcorr_coeff/a2)**2.0d0)*CP%wEFA_c
    3165        w_ax_p1 = 1.0_dl + w_ax
    3166        dorp=grhom*CP%rhorefp_ovh2*((CP%a_osc/a)**3.0d0)*dexp((wcorr_coeff**2.0d0)*3.0d0*&
    3167             &CP%wEFA_c*(1.0d0/(a2**2.0d0) - 1.0d0/(CP%a_osc**4.0d0))/4.0d0) !RL 110923
    3172        clxax_kg = ay(EV%a_kg_ix)
    3173        u_ax_kg = ay(EV%a_kg_ix+1)
    3174     end if

Physics content:

- Before switch: exact field expressions. δρ_ax = 2(φ̄̇ δφ̇/a² + m²φ̄ δφ) in the v-normalization
  (the explicit ``/a²`` because v2 = φ̇·a-coordinate convention of axion_background's ``v_vec``);
  δ_ax = δρ/ρ; heat flux u = (1+w)v with
  ``(ρ+P)v = 2 k φ̄̇ δφ /a²·...`` reduced to ``u = (1+w)·k·δφ/(φ̄̇)·(1/H0_in_Mpc_inv)`` (line 3154).
  ``w+1 = 2(φ̇²/a²)/(φ̇²/a²+m²φ²)`` (kinetic/total).
- After switch: EFA with w_EFA(a) = wEFA_c·(ℋ_osc a_osc/(m/H0 h))²/a⁴ (matches the background
  correction in ``grhoax_frac`` exactly — the same dorp expression, [PHYSICS] consistency requirement),
  and (δ_ax, u_ax) read directly from the state vector.

4.3 Axion in background/metric sums of derivs (3212–3239, 3374) — [PHYSICS]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: fortran

    3212     grhoax_t=dorp*a2
    3214     gpres_ax=w_ax*grhoax_t
    3227     gpres=gpres_ax
    3228     grho=grhob_t+grhoc_t+grhor_t+grhog_t+grhov_t+grhoax_t
    3236     dgrho=grhob_t*clxb+grhoc_t*clxc+grhoax_t*clxax_kg !RL replaced axion pert with kg
    3239     dgq=grhob_t*vb+grhoax_t*u_ax_kg !RL replaced axion pert with kg
    ...
    3374     gpres=(grhog_t+grhor_t)/3._dl+grhov_t*w_eff+gpres_ax   ! RL 081324: full gpres needed for cad2

- The axion contributes to ``grho``, ``gpres``, ``dgrho`` (δρ a²·8πG), ``dgq`` ((ρ+P)v a²·8πG) in both
  regimes with the *same* variable names. There is **no axion dgpi** (no anisotropic stress) and no
  contribution to ``pidot_sum/dgpi_diff`` — correct for both a scalar field at this order and the EFA.
- ``gpres`` at 3227 is temporarily only ``gpres_ax`` (radiation/DE pressure added later at 3374, and
  inside the tight-coupling block the OLD line ``gpres=gpres + grhov_t*w_eff`` is **removed**, with
  3374 placed *before* the TightCoupling ``adotdota=(adotoa*adotoa-gpres)/2`` at 3489 — net effect:
  ``adotdota`` now uses the **full** pressure including radiation even in the tight-coupling slip
  (OLD used ``(grhor+grhog)/3 + grhov*w_eff`` as well, via the two-step sum — equivalent), plus the
  axion pressure. Behavior change vs OLD: axion pressure now in ``adotdota`` — correct physics.
- These feed the standard constraint equations (unchanged in form, 3342–3352):

  .. code-block:: fortran

      3344     z=(0.5_dl*dgrho/k + etak)/adotoa
      3347        sigma=(z+1.5_dl*dgq/k2)            ! flat
      3348        ayprime(2)=0.5_dl*dgq              ! eta*k equation
      3350        sigma=(z+1.5_dl*dgq/k2)/EV%Kf(1)   ! nonflat
      3351        ayprime(2)=0.5_dl*dgq + CP%curv*z

  so the axion enters σ, z, η̇ purely through dgrho/dgq. CDM-frame (synchronous) conventions
  untouched.

4.4 The evolution equations (3384–3430) — [PHYSICS] (verbatim)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Exact KG (before τ_osc), synchronous gauge, conformal time in Mpc:**

.. code-block:: fortran

    3384     if (.not. EV%oscillation_started) then
    3391        !The untouched KG---------------
    3392        ayprime(EV%a_kg_ix) = dv2 * CP%H0_in_Mpc_inv/EV%renorm_c
    3393        ayprime(EV%a_kg_ix+1) = -2 * adotoa * dv2 - k2*dv1*EV%renorm_c/(CP%H0_in_Mpc_inv) &
    3394             &- a2*(CP%m_ovH0**2.0_dl)*dv1*EV%renorm_c*CP%H0_in_Mpc_inv - k*z*v2_bg !RL 050324

i.e. with y₁ = δφ/renorm_c, y₂ = δφ̇/H0 (a-coordinate):
δφ′ = δφ̇; δφ̇′ = −2ℋ δφ̇ − (k² + a²m²) δφ − (ḣ_L/2) φ̄̇, with ``k*z = ḣ_L/2`` (CAMB sign conventions)
and all H0 factors converting the v-normalizations to Mpc⁻¹ time units.

**Effective fluid (after τ_osc):** adiabatic sound speed, EFA sound speed, then the GDM equations:

.. code-block:: fortran

    3395     else !Now past the oscillation phase, use EFA
    3397        cad2 = w_ax*((1._dl+ gpres/grho)/w_ax_p1 +1._dl) !RL 081324 ...
    3398        kamnorm = k2/(a2*((CP%m_ovH0*CP%H0_in_Mpc_inv)**2.0_dl))
    3399      if (kamnorm .lt. 1.e-14_dl) then
    3402         csquared_ax = kamnorm/4.0_dl + 5.0_dl*(adotoa**2.0_dl)/(4.0_dl*(k2/kamnorm))
    3404       else
    3405          csquared_ax = (sqrt(1.0_dl + kamnorm) - 1.0_dl)**2.0_dl/(kamnorm) + 5.0_dl*(adotoa**2.0_dl)/(4.0_dl*(k2/kamnorm))
    3408      end if
    3424      clxax_kg_dot = -k*(u_ax_kg + z*w_ax_p1)-3.0_dl*(csquared_ax-w_ax)*adotoa*clxax_kg-&
    3425           &9.0_dl*(csquared_ax-cad2)*(adotoa**2.0_dl)*u_ax_kg/k
    3426      u_ax_kg_dot=-adotoa*u_ax_kg+3.0_dl*csquared_ax*adotoa*u_ax_kg+k*csquared_ax*clxax_kg+&
    3427           &3.0_dl*(w_ax-cad2)*adotoa*u_ax_kg
    3428        ayprime(EV%a_kg_ix) = clxax_kg_dot
    3429        ayprime(EV%a_kg_ix+1) = u_ax_kg_dot
    3430     end if

- **Sound speed** (rest-frame): with κ ≡ k²/(a²m²) (``kamnorm``, m in Mpc⁻¹),
  ``cs²(k,m,a) = (√(1+κ) − 1)²/κ + (5/4)·(ℋ/(a m))²`` — the first term is the exact relativistic
  field cs² (→ κ/4 for κ≪1, →1 for κ≫1); the second ``5ℋ²/(4a²m²)`` term is the AxiECAMB H/m
  correction (note ``k2/kamnorm = a²m²``). For κ < 1e-14 the first term is Taylor-expanded to κ/4 to
  avoid catastrophic cancellation. Commented alternatives (3/2 coefficient, 1.1 coefficient) record
  the calibration history — port the 5/4.
- **cad²** = w·[(1+gpres/grho)/(1+w) + 1] (line 3397) — an approximation to ẇ/(3ℋ(1+w))-corrected
  adiabatic sound speed built from the *total* background pressure/density (uses ``gpres`` assembled
  at 3374 including the axion); the comment notes the later PPF modification of gpres doesn't matter.
- **GDM equations** are the standard synchronous-gauge (1+w)-weighted velocity form:
  δ′ = −k(u + (1+w)z) − 3ℋ(cs²−w)δ − 9ℋ²(cs²−cad²)u/k,
  u′ = −ℋ(1−3cs²)u + k cs² δ + 3ℋ(w−cad²)u, with u=(1+w)v.

4.5 Legacy GDM pair zeroed (3469–3473) — [OBSOLETE]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: fortran

    3469     !axionCAMB's GDM variables
    3470     clxaxdot = 0. !RL: already useless, but didn't get rid of due to isocurvature
    3471     v_axdot = 0.
    3472     ayprime(EV%a_ix)=clxaxdot
    3473     ayprime(EV%a_ix+1)=v_axdot

``clxax``/``v_ax`` are frozen at their initial values and used nowhere in the metric. The huge commented
block at 3433–3467 is the original axionCAMB GDM system kept for reference. Do **not** port; but see
§5.3 — the axion isocurvature ICs still target these dead variables.

4.6 PPF block additions (3337–3339) — [OBSOLETE]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: fortran

    3337        EV%dgrhoec_ppf=dgrhoe+3._dl*(1+w_eff)*grhov_t*adotoa*vT/k
    3338        EV%dgqec_ppf=dgqe+(1+w_eff)*grhov_t*sigma
    3339        EV%vTc_ppf=vT+sigma

Computed every call, never consumed. (Rest-frame/comoving PPF variables; likely imported from
equations_cross for debugging.) The PPF Γ equation, ``dgqe/dgrhoe``, ``Fa``, ``S_Gamma`` etc. are
**bit-identical to OLD** — no axion modification of PPF dynamics. Note ``vT = dgq/(grhoT+gpres)``
with ``grhoT = grho − grhov_t`` now implicitly *includes the axion* in the "matter" total
(grho and dgq contain axion terms) — that is the correct PPF prescription (DE vs everything else).

4.7 RSA / no-multipole interactions — no formula change, axion enters implicitly
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The RSA reconstructions (3296–3300, 3312–3324; also ``output`` 2055–2061, 2076–2085) are unchanged in
form:

.. code-block:: fortran

    3298        z=(0.5_dl*dgrho/k + etak)/adotoa
    3299        dz= -adotoa*z - 0.5_dl*dgrho/k
    3300        clxr=-4*dz/k ;  qr=-4._dl/3*z ; pir=0

but since ``dgrho``, ``adotoa`` now include the axion, the RSA photon/neutrino monopoles automatically
see the axion potential. [PHYSICS — comes for free if dgrho/grho are right.]

5. Initial conditions (``initial``, lines 2430–2774)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

5.1 IC vector enlargement (2443–2448) — [PLUMBING]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: fortran

    2443     integer, parameter :: i_clxg=1,i_clxr=2,i_clxc=3, i_clxb=4, &
    2444          i_qg=5,i_qr=6,i_vb=7,i_pir=8, i_eta=9, i_aj3r=10,i_clxq=11,i_vq=12, &
    2445          i_clxax=13,i_v_ax=14, i_dphi_ax=15, i_dphidot_ax=16 !RL adding KG variables
    2447     integer, parameter :: i_max = i_dphidot_ax
    2448     real(dl) initv(6,1:i_max), InitVec(1:i_max)

(``initv`` first dimension was already 6 in OLD.) A batch of declared-but-unused helper arrays for an
aborted a_osc/aeq calculation (``lnamin, lnamax, a_arr, fax_arr, ...``, 2452–2462) — [COSMETIC].

5.2 Adiabatic mode: axion entries (2547–2584) — [PHYSICS]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Background/setup additions:

.. code-block:: fortran

    2522     om = (grhob+grhoc)/sqrt(3*(grhog+grhonu))  ! DM: for a_i<<a_osc, axions do not contribute in background expansion
    2527     Rc=(CP%omegac)/(CP%omegac+CP%omegab) ! DM: no axions in here for same reason
    2536     adotoa = 1.0_dl/(a*dtauda(a))  !RL added ... updated to using dtauda
    2540     call spline_out(loga_table,phinorm_table,phinorm_table_ddlga,ntable,dlog10(a),v1_bg)

GDM entries zero; the field ICs are set from the analytic early-time (frozen-field) solution:

.. code-block:: fortran

    2557     initv(1,i_clxax)=0
    2558     initv(1,i_v_ax)=0
    2567     dgrho = (grhonu + grhog)*initv(1,i_clxg)/a2 + (grhob + grhoc)*initv(1,i_clxb)/a !RL added, adiabatic initial pert conditions
    2569     z = (0.5_dl*dgrho/k -initv(1,i_eta)*k/2)/adotoa
    2577     initv(1,i_dphi_ax) = ((CP%m_ovH0*CP%H0_in_Mpc_inv)**2.0d0) * &
    2578          &(chi*EV%Kf(1)*(1-omtau/5)*k*x/2) * a2 * v1_bg/(210.0d0 * (adotoa**3.0d0)) !RL - clxcdotmethod
    2579     initv(1,i_dphidot_ax) = (CP%m_ovH0**2.0d0) * CP%H0_in_Mpc_inv &
    2580          &* (chi*EV%Kf(1)*(1-omtau/5)*k*x/2) * a2 * v1_bg/(35.0d0 * (adotoa**2.0d0)) !RL- clxcdotmethod

- The factor ``(chi*EV%Kf(1)*(1-omtau/5)*k*x/2)`` is the analytic ``clxcdot = −k z`` of the adiabatic
  mode (RL verified ≈ k·z to 1e-8; the ``z`` computed at 2569 is only used in commented checks —
  the "kzmethod" variants at 2582–2583 are commented out).
- Physics: deep in RD with m ≪ H the field is frozen; the leading forced response to ḣ is
  δφ = m² (kz) a² φ̄ /(210 ℋ³), δφ̇ = m² (kz) a² φ̄ /(35 ℋ²) (in Mpc/v-normalized units as written).
  These must be re-derived/checked against the new CAMB IC normalization when porting (CAMB 1.6.7
  uses the same ``initv`` chi=−1 convention, so the expressions carry over).
- The adiabatic sign flip is unchanged: ``InitVec = −InitVec`` for adiabatic (2713).

5.3 NEW axion isocurvature mode 6 (2633–2698) — [PHYSICS] (but currently incomplete — see risk)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: fortran

    2635        !axion isocurvature mode
    2636        ! DM: See pdf notes on axion theory
    2637        ! Derived by Dan Grin using a Matrix ODE formalism
    2638        ! All assume tau_i<<tau_osc
    2639        initv(6,:)=0.0d0
    2640        Ra=grhoax/(grhog+grhonu)/a_osc**3.0d0
    2641        omr=(grhog+grhonu)/(grhob+grhoc)
    2642        frac=grhoax/(grhoax+grhob+grhoc) ! Fraction of matter in axions
    2650        if(a_osc.le.CP%aeq) then
    2651           FF=1.0d0 ! DM: F(f) sets a_0 from a_osc FF=1 if aosc<aeq
    2652           Konstant=1.0d0-frac
    2653        else
    2654           FF=((1.0d0-frac)+frac*(CP%aeq/a_osc)**3.0d0)**(-1.0d0)
    2655           Konstant=(1.0d0-frac)/((1.0d0-frac)+frac*((CP%aeq/a_osc)**3.0d0))
    2656        end if
    2660        initv(6,i_clxax)=1.0d0-(x*x)/(10.0d0)&
    2661             -Konstant*(x*x)*omtau/(180.0d0)+(1.0d0/600.0d0)*(x**4.0d0)&
    2662             +(53.0d0*(Konstant**2.0d0)*(x**2.0d0)*(omtau**2.0d0))/(140.0d0*5.0d0*16.0d0)
    2663        ! Normalise to \delta_a=1.
    2664        initv(6,i_v_ax)=-x/5.0d0+(x*x*x)/30.0d0+Konstant*x*omtau/(30.0d0)&
    2665             -x*(Konstant**2.0d0)*(omtau**2.0d0)/(84.d0)&
    2666             +(11.0d0*(Konstant**3.0d0)*x*(omtau**(3.0d0)))/(42.0d0*64.0d0)&
    2667             -Konstant*(x**3.0d0)*omtau/(120.0d0)-(x**5.0d0)/(600.0d0)&
    2668             -(11.0d0*(Konstant**4.0d0)*x*(omtau**4.0d0))/(63.0d0*128.0d0)&
    2669             +(929.0d0*(Konstant**2.0d0))*(x**3.0d0)*(omtau**2.0d0)/(3780.0d0*16.0d0*5.0d0)&
    2670             +((CP%m_ovH0**2.0d0)*x*(omtau**4.0d0)/(225.0d0*5.0d0*64.0d0*CP%omegar&
    2671             *(((CP%omegab+CP%omegac+CP%omegaax)/(4.0d0*CP%omegar))**4.0d0)))*(FF**4.0d0)
    2679        AA=Ra*(omr**4.0d0)*(FF**4.0d0)
    2682        initv(6,i_clxg)= -AA/3.*omtau**4.
    2683        initv(6,i_clxr)=initv(6,i_clxg)
    2684        initv(6,i_clxb)=0.75_dl*initv(6,i_clxg)
    2685        initv(6,i_clxc)=initv(6,i_clxb)
    2686        initv(6,i_eta)=-(0.50*initv(6,i_clxg))
    2687        !DM: this eta is -2eta_s, see notes and below setting y(2)
    2689        initv(6,i_qg)=initv(6,i_clxg)*x/15.0d0
    2690        initv(6,i_qr)=initv(6,i_clxg)*x/15.0d0
    2691        initv(6,i_vb)=0.75_dl*initv(6,i_qg)
    2692        initv(6,i_pir)= (-9.0d0*Rb)*omtau*initv(6,i_clxg)*(1.0d0-Konstant)/(5.0d0*(75.0d0+4.0d0*Rv))
    2694        initv(6,i_aj3r)= initv(6,i_pir)*x/(7.0d0)
    2698     end if

Notes:

- ``grhoax``, ``a_osc`` here are the **module-level** (modules.f90) variables, not ``CP%`` ones (RL
  comment at 2640 flags this confusion). ``x=kτ``, ``omtau`` as usual. Normalized to δ_ax = 1.
- **RISK / incompleteness**: mode 6 populates the *dead* GDM pair (``i_clxax``, ``i_v_ax`` → ``y(EV%a_ix)``)
  and leaves the live KG pair (``i_dphi_ax``, ``i_dphidot_ax``) at 0. Since derivs uses only ``clxax_kg``
  (from the KG pair) in dgrho/dgq, an ``initial_condition=6`` run in AxiECAMB evolves **zero axion
  perturbation** plus the small compensating photon/neutrino entries. RL's comment at 2659
  acknowledges this: *"we might need to incorporate into the field variables too in the future. But
  stay tuned for now."* For the port: the mode-6 IC must be re-derived in terms of (δφ, δφ̇)
  (e.g. δφ = φ̄·δ_a-type relation for a frozen field) — quote the old series only as a reference for
  the compensated radiation part. The isocurvature amplitude parameters (``alpha_ax``, ``Hinf``) are
  handled in modules/cmbmain, not here.
- The IC vector is **not** sign-flipped for mode 6 (only adiabatic flips, 2713 unchanged).

5.4 State-vector initialization incl. ``renorm_c`` (2728–2737) — [PHYSICS]/[ACCURACY]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: fortran

    2730     y(EV%a_ix)=InitVec(i_clxax)
    2731     y(EV%a_ix+1)=InitVec(i_v_ax)
    2734     EV%renorm_c = sqrt(CP%rhorefp_ovh2)*CP%H0*((k*CP%tau_osc)**2._dl)/CP%m_ovH0
    2735     y(EV%a_kg_ix) = InitVec(i_dphi_ax)/EV%renorm_c !RL 050324
    2736     y(EV%a_kg_ix+1) = InitVec(i_dphidot_ax) !RL 050324

``renorm_c`` is a per-k constant rescaling of the stored δφ
(= √(ρ_ref/ρcrit h²)·H0·(kτ_osc)²/(m/H0)) chosen so the stored variable stays O(1) over the
enormous dynamic range of δφ for large m (k τ_osc large). Only ``y(a_kg_ix)`` is rescaled; ``dv2`` is
not. Every physical use multiplies back by ``EV%renorm_c`` (derivs 3145/3154/3392–3394, output
1949/1959, outtransf 2992, CopyScalarVariableArray 1028/1042/1046/1056). [ACCURACY — numerically
essential for high masses; port the idea, the exact constant is a free choice.]

6. Sources & outputs
~~~~~~~~~~~~~~~~~~~~~

6.1 ``output`` (lines 1843–2292)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Signature change — [PLUMBING]:

.. code-block:: fortran

    1843   subroutine output(EV,y,j,tau,sources,dgpi_out)
    1851     real(dl), optional, intent(out) :: dgpi_out !RL 091023 for switch boundary condition of source integration
    ...
    2285     if (present(dgpi_out)) then !RL 091023
    2286        dgpi_out = dgpi
    2289     end if

(``dgpi_out`` was meant for the commented-out variant of ``metric_delta(1)`` that includes dgpi terms;
cmbmain may call output with it — check cmbmain report. Currently the value is exported but the
consumer usage is in commented code → [OBSOLETE]-leaning, but cheap to keep.)

Axion reconstruction (1932–1989) — identical logic to derivs §4.2 (KG-side from ``y(a_kg_ix)·renorm_c``
splines; EFA-side direct), quoted at lines 1945–1959 and 1969–1988 (see §4.2 for the formulas;
``v_ax_test`` is diagnostic). Then — [PHYSICS]:

.. code-block:: fortran

    1993     grhoax_t=dorp*a2
    1996     clxax=y(EV%a_ix)            ! legacy, unused below
    1997     v_ax=y(EV%a_ix+1)
    2005     gpres_ax=w_ax*grhoax_t
    2011     grho=grhob_t+grhoc_t+grhor_t+grhog_t+grhoax_t ! RH removed the lambda term for later
    2012     gpres=(grhog_t+grhor_t)/3+gpres_ax ! RH removed the lambda pressure term for later
    2016     dgrho=grhob_t*clxb+grhoc_t*clxc+grhoax_t*clxax_kg !RL replaced with axion pert from kg
    2024     dgq=grhob_t*vb+u_ax_kg*grhoax_t ! only baryons and axions for now
    2026     if (is_cosmological_constant) then
    2027        w_eff = -1_dl
    2028        grhov_t=grhov*a2
    2029        grho = grho+grhov_t ! add in the DE
    2031     else
    2033        w_eff=w_de(a) ; grhov_t=grho_de(a)/a2 ; grho = grho+grhov_t
    2036        dgrho=dgrho+EV%dgrho_e_ppf
    2037        dgq=dgq+EV%dgq_e_ppf
    2038     end if
    2041     gpres=gpres+grhov_t*w_eff ! RH add in the DE now

The DE terms were just reordered (moved after the axion) — same totals as OLD plus axion. The
remainder of the source assembly — z, σ, ppiedot, polter, ISW, ``sources(1)``, ``sources(2)``,
lensing ``sources(3) = -2*phi*f_K(...)`` with ``phi = -(dgrho +3*dgq*adotoa/k)/(k2*EV%Kf(1)*2)`` —
is **formula-identical to OLD** (2102–2228); the axion affects them only through grho/gpres/dgrho/
dgq/adotoa/σ/z. ISW therefore automatically picks up the axion (early-ISW around the switch is
where ``metric_delta``/``deltaBCSrc`` matter).

``opac(j)``/``dopac(j)`` were renamed through locals (2071–2072 ``opacity_use = opac(j)``,
``dopacity_use = dopac(j)``) with no value change — [COSMETIC]. The ``csquared_ax`` recomputation at
2263–2271 (with a 1e-13 Taylor threshold instead of derivs' 1e-14) feeds only commented writes —
[COSMETIC]. ``qgdot = -4*dz/3`` added in the no_phot_multpoles branch (2085): present in OLD too?
**No — this is new** (OLD set ``qgdot`` only implicitly). In OLD the ``no_phot_multpoles`` branch of
``output`` did not set ``qgdot`` at all before use in sources(1) → actually OLD set it via
``qgdot =yprime(EV%g_ix+1)`` only in the else; AxiECAMB adds the RSA-consistent
``qgdot = -4*dz/3`` — [ACCURACY] bugfix-like (modern CAMB has the same term; no port action needed).

6.2 NEW ``GrowthRate`` (lines 1705–1839) — [OBSOLETE] (currently disabled)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

DM-era helper computing f = dln δ/dln a for cdm+baryon+axion (``clxtot`` from
``dgrhomat/grhomat`` with ``clxax_kg``; KG- or EFA-side construction as in §4.2; explicitly *excludes*
massive neutrinos; uses ``grhov_t=grhov*a**(-1-3*w_lam)`` i.e. constant-w DE only). Its only call site
(outtransf 3063) is commented out and ``growth = 0._dl`` is hardwired (3066). Port only if the
``Transfer_f`` column is wanted; otherwise drop.

6.3 ``outtransf`` (lines 2924–3070) — [PHYSICS]/[PLUMBING]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Signature: ``subroutine outtransf(EV, y, tau,Arr)`` (new ``tau`` argument; cmbmain call sites updated).
Axion δ reconstruction (2980–3017) identical pattern to §4.2 (KG: ``clxax_kg = drhoax_kg/grhoax_kg``
with renorm_c; EFA: ``clxax_kg = y(EV%a_kg_ix)``), then:

.. code-block:: fortran

    3019     grhoax_t=dorp*(a**2.0d0)
    3024     Arr(Transfer_axion) = clxax_kg/k2 !RL replaced the pert with that from kg
    ...
    3051     if (CP%m_ovH0 .ge. 10._dl) then !RL 070324
    3052        dgrho = dgrho+(clxc*grhoc + clxb*grhob)/a+grhoax_t*clxax_kg !RL replaced axion pert with kg
    3053        grho =  grho+(grhoc+grhob)/a+grhoax_t
    3054     else
    3056        dgrho = dgrho+(clxc*grhoc + clxb*grhob)/a
    3057        grho =  grho+(grhoc+grhob)/a
    3058     end if
    3066     growth = 0._dl
    3067     Arr(Transfer_f) = growth
    3068     Arr(Transfer_tot) = dgrho/grho/k2

- ``Transfer_axion`` (=7) and ``Transfer_f`` (=8) are new columns (modules.f90 ``Transfer_max`` reworked).
- **Clustering kluge** (Hlozek et al. 2015 bias treatment): the axion is included in
  ``Transfer_tot`` (δρ_m/ρ_m) only when oscillating-DM-like, criterion changed from the DM16
  ``CP%ma ≥ 1e-25 eV`` (commented at 3049–3050) to ``CP%m_ovH0 ≥ 10`` (m ≥ 10 H0). [PHYSICS] — port the
  criterion (it also must stay consistent with the matching kluge in halofit, see that report).

6.4 ``outputt`` (tensors, 2296–2374) — **no non-cosmetic change.**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The tensor source formulas are bit-identical to OLD (re-indent only). The axion affects tensors only
through the background in ``derivst`` (§7). ``outputv`` likewise unchanged.

7. Tensor & vector derivatives
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

7.1 ``derivst`` (lines 3928–4169) — [PHYSICS] + **a porting trap**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Only change: axion added to the background expansion (3946–3949 declarations; 3968–4009):

.. code-block:: fortran

    3968     if (.not. EV%oscillation_started) then !RL
    3987        call spline_out(loga_table,phinorm_table,phinorm_table_ddlga,ntable,dlog10(a),v1_bg)
    3988        call spline_out(loga_table,phidotnorm_table,phidotnorm_table_ddlga,ntable,dlog10(a),v2_bg)
    3989        if (v1_bg .eq. v1_bg .and. v2_bg .eq. v2_bg) then
    3991           grhoax_kg = (v2_bg)**2.0_dl/a2+(CP%m_ovH0*v1_bg)**2.0_dl
    3992           dorp = grhom*grhoax_kg/(CP%H0**2.0d0/1.0d4)
    3993        else
    3994           dorp=0.0d0
    3995        endif
    3996     else
    3999        wcorr_coeff = CP%ahosc_ETA*CP%a_osc/(CP%m_ovH0*(CP%H0/100.0d0)) !RL082924
    4001        dorp=grhom*CP%rhorefp_ovh2*((CP%a_osc/a)**3.0d0)*dexp((wcorr_coeff**2.0d0)*3.0d0*&
    4002             &CP%wEFA_c*(1.0d0/(a2**2.0d0) - 1.0d0/(CP%a_osc**4.0d0))/4.0d0)
    4004     endif
    4006     grhoax_t=dorp*(a**2.0d0)
    4009     grho=grhob_t+grhoc_t+grhor_t+grhog_t+grhov_t+grhoax_t

**TRAP**: for tensor runs ``EV%oscillation_started`` is initialized ``.false.`` in ``GetNumEqns`` (1270)
and **never switched** (only ``GaugeInterface_EvolveScal`` flips it). So ``derivst`` *always* takes the
KG-table branch, and for a > a_osc (table end) ``spline_out`` blindly extrapolates the last cubic
segment of an oscillating field (subroutines.f90:12–36 has no range clamp; the NaN guard only
catches NaN, not bad extrapolation). The tensor background after a_osc is therefore wrong/unreliable
in AxiECAMB. In the port, use the single unified ``rho_ax(a)`` background function (as ``dtauda`` does
via ``grhoax_frac``) for tensors — that automatically fixes this. Flag for validation against
AxiECAMB tensor spectra (differences expected!).

7.2 ``derivsv`` (vectors, 3706–3923) — [OBSOLETE]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Dead code: line 3729 ``stop 'ppf not implemented for vectors'`` executes before anything else (as in
OLD). The added axion block (3771–3811) mirrors derivst but contains two latent bugs (harmless since
unreachable): ``gpres_ax=w_ax_p1*grhoax_t`` (3810 — uses 1+w instead of w) and ``w_eff`` used
uninitialized at 3813. Do not port.

8. Tight coupling, RSA, approximation switches, accuracy knobs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Tight coupling**: no change to the slip/pig formulas. Only difference: ``gpres`` entering
   ``adotdota`` (3489) is now the full pressure incl. axion, assembled at 3374 instead of inside the
   TightCoupling branch (OLD line ``gpres=gpres + grhov_t*w_eff`` removed). [PHYSICS — small, correct.]
2. **RSA**: formulas unchanged; axion enters via dgrho (see §4.7). The new ``qgdot=-4*dz/3`` in
   output's no_phot_multpoles branch (2085) — [ACCURACY], already present in modern CAMB.
3. **Switch scheduling**: ``tau_switch_oscillation`` checked first; ``smallTime=0`` (531) [ACCURACY];
   ``DoLateRadTruncation`` left on (comment at 351–353 notes inidriver controls it).
4. ``max_l_evolve`` **1024→512** (367) [ACCURACY/PLUMBING — hardcoded; modern CAMB doesn't need it].
5. **lmaxnr low-k boost** ``min(8, sqrt(scal)*450*q)`` + Motloch clamp (1349–1355) [ACCURACY,
   physics-motivated for switch-era ISW smoothness; revisit in 1.6.7 which has its own defaults].
6. ``DeltaTimeMaxed`` (749–757): ``t= DeltaTime(a1,a2)`` — the ``tol`` argument is **dropped**
   (OLD passed it through). Affects massive-ν switch-time precomputation only. Looks accidental;
   [ACCURACY, minor — do not replicate].
7. ``GetNumEqns`` **lmaxnu comment tweaks** (1300–1310): net values unchanged from OLD defaults
   ([COSMETIC] besides the lmaxnr item above).
8. **dverk restarts**: every switch (incl. oscillation) sets ``ind=1``; no other integrator changes.
   The dverk source itself is unchanged (it lives elsewhere).
9. ``MassiveNuVarsOut`` (1443–1506): re-enabled ``dpnu`` bookkeeping
   (``dpnu=y(EV%nu_ix(nu_i)+1)`` at 1474; ``dpnu=dpnu/rhonu`` at 1482) — computed but not consumed
   (the ``dgp`` accumulation line was removed again). [COSMETIC/OBSOLETE]

9. PPF / dark-energy summary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``LambdaGeneral`` module: **no functional change** (re-indent only). ``w_de``, ``grho_de``,
  ``setddwa``, ``interpolrde``, ``setcgammappf``, ``cubicsplint`` identical to OLD. [COSMETIC]
- PPF perturbation (Γ) equation, ``dgrho_e_ppf/dgq_e_ppf``, ``ppiedot`` in output: identical to OLD.
- Only additions: the never-read ``EV%dgrhoec_ppf/dgqec_ppf/vTc_ppf`` (§4.6) [OBSOLETE], and the
  implicit fact that the PPF "matter" totals now include the axion (correct).
- ``is_cosmological_constant`` logic unchanged; the axion never routes through the DE machinery
  (it is its own component).

10. Porting checklist & risk register
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Item
     - Class
     - Action
   * - ``grhoax_frac`` (table < a_osc, analytic ⨯ exp-w-correction ≥ a_osc) + ``dtauda`` term
     - PHYSICS
     - Port as the component's ``grho_ax(a)``; keep ``m_ovH0`` vs ``ma/H0_eV`` consistent (one source of truth).
   * - ``tau_osc`` setup + clamps in ``init_background``; DeltaTime tol 1e-8
     - PHYSICS/ACCURACY
     - Port; in 1.6.7 do it in the component's ``Init`` after background.
   * - 4-slot state vector (2 dead GDM + δφ,δφ̇ → δ,u)
     - PLUMBING
     - Re-derive: only 2 live equations needed; drop the GDM pair.
   * - ``renorm_c`` δφ rescaling
     - ACCURACY
     - Port concept (any O(1)-keeping rescale works).
   * - KG equations (3392–3394) with ``k*z`` metric source
     - PHYSICS
     - Port verbatim (synchronous ḣ/2 = kz).
   * - EFA equations (3424–3427), cs² with +5/4(ℋ/am)² term (3398–3408), cad² (3397), w_EFA(a)=wEFA_c(wcorr/a²)²
     - PHYSICS
     - Port verbatim incl. 1e-14 Taylor branch.
   * - Switch matching in CopyScalarVariableArray (tU/tV/tW, tdvarphi_*, tdrho_ef, u_ax_ef→u_ax_efa, δ-shift for σ continuity, weight, η shift, metric_delta export)
     - PHYSICS
     - Port exactly; needs ``tvarphi_*``, ``A_coeff``, ``ah_osc``, ``ahosc_ETA``, ``Prefp``, ``rhorefp_ovh2``, ``wEFA_c`` from background module; re-derive the source-jump bookkeeping against 1.6.7's source assembly (deltaBCSrc lives in cmbmain).
   * - dgrho/dgq axion terms in derivs & output; no dgpi
     - PHYSICS
     - Port (both regimes).
   * - Adiabatic δφ,δφ̇ ICs (2577–2580)
     - PHYSICS
     - Port; verify against 1.6.7 IC normalization (chi convention identical).
   * - Iso mode 6 series (2639–2698)
     - PHYSICS (incomplete)
     - **Targets dead variables; no KG iso ICs exist.** Re-derive δφ-based iso ICs; old series useful only for the compensated radiation entries.
   * - ``Transfer_axion=δ_ax/k²``, ``Transfer_f``, clustering criterion ``m_ovH0 ≥ 10`` in Transfer_tot
     - PHYSICS/PLUMBING
     - Port column + criterion; keep consistent with halofit kluge.
   * - GrowthRate subroutine
     - OBSOLETE
     - Disabled (growth=0); port only if Transfer_f wanted.
   * - Tensor background via always-KG branch in derivst
     - PHYSICS+BUG
     - In port use unified ρ_ax(a); AxiECAMB extrapolates the field table past a_osc for tensors (validate; expect differences).
   * - derivsv axion block
     - OBSOLETE
     - Unreachable; contains bugs; drop.
   * - ``smallTime=0``, oscillation switch first, ``ind=1`` restart
     - ACCURACY/PLUMBING
     - Port the guaranteed-switch behavior.
   * - lmaxnr 450-coefficient low-k boost + Motloch clamp
     - ACCURACY
     - Motloch already in 1.6.7; re-test the 450 boost need.
   * - ``max_l_evolve=512``, DeltaTimeMaxed tol drop, opacity_use renames, dgrhoec_ppf trio, empty if-blocks, debug writes
     - OBSOLETE/COSMETIC
     - Ignore.
   * - ``GetOmegak`` removal (curvature fixed at input)
     - PLUMBING
     - Handled by 1.6.7 params logic; ensure Ω_ax in the budget.
   * - ``dgpi_out`` optional arg of output
     - PLUMBING
     - Consumer logic is commented out; only needed if reviving the dgpi-including metric_delta variant.


Original-code analysis: main loop and miscellaneous files (cmbmain.f90 and others)
-----------------------------------------------------------------------------------

Port inventory: driver/main-loop and small files (cmbmain.f90, cmbmainOMP.f90,
camb.f90, subroutines.f90, lensing.f90, bessels.f90, constants.f90, writefits.f90,
utils.F90, Makefiles).

Comparison base: ``OLDCAMB/`` (pristine CAMB Nov13) -> ``AxiECAMB/``.
All line references are to the AxiECAMB files unless prefixed ``OLDCAMB``.

Method note: the raw ``cmbmain.f90.diff`` is 3964 lines, but a whitespace-insensitive,
comment-insensitive diff (``diff -w -B`` + comment filter) reduces the file to exactly
**163 non-comment changed lines**. Every one of them is accounted for below. The same
filtering was applied to all other files in this report's scope.

0. Build system: what is actually compiled
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

AxiECAMB/Makefile (vs OLDCAMB/Makefile)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Compiler switched from ifort to gfortran:
  ``F90C = gfortran``, ``FFLAGS = -O3 -fopenmp -ffpe-summary=none`` (ifort/F90CRLINK/FISHER/MKL blocks all commented out).
  **[OBSOLETE]** — modern CAMB has its own build system; nothing to port. The only
  semi-meaningful flag is ``-ffpe-summary=none`` (suppress FP-exception summaries at exit),
  which hints the axion code generates benign FP exceptions (underflow etc.); harmless.
- FITS/HEALPIX paths commented out; ``EXTCAMBFILES`` block commented out (so it is empty). **[OBSOLETE]**

AxiECAMB/Makefile_main (vs OLDCAMB/Makefile_main)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Module selection (this defines which physics files are live):

  .. code-block:: make

      EQUATIONS     ?= equations_ppf        (was: equations)
      RECOMBINATION ?= recfast_axion        (was: recfast)
      NONLINEAR     ?= halofit_ppf          (was: halofit)
      DRIVER        ?= inidriver_axion.F90  (was: inidriver.F90)

  **[PLUMBING]** — in modern CAMB the equivalents are equations.f90 / recfast.f90 /
  halofit.f90 selected at build; the port must base its changes on those modern files
  (the \*_ppf/_axion variants are inventoried by other reports).
- New object + explicit rule:

  .. code-block:: make

      CAMBOBJ = constants.o utils.o $(EXTCAMBFILES) subroutines.o inifile.o $(POWERSPECTRUM).o $(RECOMBINATION).o \
            $(REIONIZATION).o modules.o bessels.o $(EQUATIONS).o $(NONLINEAR).o lensing.o $(BISPECTRUM).o cmbmain.o camb.o axion_background.o

      #RL added dependency fix of axion_background on modules
      axion_background.o: axion_background.F90 modules.o
                          $(F90C) $(F90FLAGS) -c axion_background.F90

  **[PLUMBING]** — the new compilation unit ``axion_background.F90`` (KG background solver)
  must be added to the modern CAMB build (fortran/Makefile) with a dependency on the
  module that holds the parameter type.
- ``F90CRLINK ?= -lstdc++`` commented out; ``clean`` also removes the ``camb`` binary. **[COSMETIC]**

KEY FACT: ``cmbmain.o`` is built from ``cmbmain.f90`` via the implicit ``%.o: %.f90`` rule.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``cmbmainOMP.f90`` is referenced **nowhere** (not in CAMBOBJ, no rule, not ``include``\ d, not
``use``\ d by any source file). **cmbmainOMP.f90 is dead code** — see section 2.

1. cmbmain.f90 — the compiled main loop
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1.1 Module-level declarations (top of module CAMBmain)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- cmbmain.f90:80 (inside ``type IntegrationVars``):

  .. code-block:: fortran

      real(dl), dimension(:,:), pointer :: Source_q, ddSource_q, metricdeltas_q !RL 090323

  (was ``:: Source_q, ddSource_q``). ``IV%metricdeltas_q(2,SourceNum)``: per-k boundary-term
  amplitudes at the KG->fluid switch, interpolated onto the integration k-grid.
  **[PHYSICS]** (data structure for the switch boundary condition, see 1.8/1.9).

- cmbmain.f90:90:

  .. code-block:: fortran

      real(dl), dimension(:,:,:), allocatable :: Src, ddSrc, deltaBCSrc, dd_deltaBCSrc
      !Sources and second derivs !RL 090323 added delta functions at the boundary to address
      !discontinuity in sources - and note that these second derivatives are for interpolations
      !across k space, not interpolations across time

  ``deltaBCSrc(Evolve_q%npoints, SourceNum, 2)`` holds, per evolved k, the coefficients of
  the delta-function (index 3rd-dim=1) and delta-prime (3rd-dim=2) pieces of the scalar
  source at tau_osc; ``dd_deltaBCSrc`` are its spline second derivatives **in k**. **[PHYSICS]**

- cmbmain.f90:102 ``qmin0=0.1_dl`` — value unchanged, comment added ("RL changed to 1000.0 for glitch testing, default is 0.1"). **[COSMETIC]**
- cmbmain.f90:117 ``fixq = 0._dl`` — value unchanged, debug values appended in comment. **[COSMETIC]**
- cmbmain.f90:129-130 ``real clock_start, clock_stop``, ``real clock_totstart, clock_totstop ! RH timing`` — declarations for (now fully commented-out) timing probes. **[COSMETIC]**

1.2 ``subroutine cmbmain`` body (lines 127-400)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All changes here are mechanical: one-line ``if (...) call X`` statements expanded to
``if (...) then / call X / end if`` blocks so that (commented-out) cpu_time probes could be
inserted, plus dozens of commented ``! RH timing`` lines and a commented note
``!!! RH axions need to be called again here !!!`` near the ClTransferToCl call (line 374;
nothing is actually called again — the executable code is identical to OLDCAMB).
Affected statements (semantics unchanged): InitializePowers guard (183-191),
SetkValuesForSources (197-204), InitTransfer (206-213), DoSourcek loop (239-252,
``!$OMP PARALLEL DO ... PRIVATE(EV, q_ix)`` retained), Transfer_Get_sigma8 guard (297-304),
GetLimberTransfers (325-332), ClTransferToCl guard (371-379).
**[COSMETIC]** — port none of it.

OMP status: every parallel region present in OLDCAMB is still present and active in
AxiECAMB cmbmain.f90 (DoSourcek loop line 239; SourceToTransfers loop line 343
``!$OMP PARAllEl DO DEFAUlT(SHARED),SHARED(TimeSteps), SCHEDUlE(STATIC,4)``; Limber Cls
line 474; TransferOut line 1316; MakeNonlinearSources line 1384; InitSourceInterpolation
line 1430 + the NEW one at 1447; InterpolateCls line 2844). No threading changes to port.

1.3 ``function GetTauStart(q)`` (lines 770-830) — axion start time **[PHYSICS]**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

New declarations:

.. code-block:: fortran

    real(dl) tauosc,om,rhomass,grhonu,taurad,tauax_init,taueq
    double precision arad,abar,taubar

New block after the massive-neutrino constraint (quoted verbatim, debug comments elided):

.. code-block:: fortran

    !Start when axions are not oscillating
    !update start time to use a_init from w_evolve routine 9/19 dgrin, comment out other checks (superfluous)
    if(CP%Omegaax>0) then
       rhomass =  sum(grhormass(1:CP%Nu_mass_eigenstates))
       grhonu=rhomass+grhornomass
       om = (grhob+grhoc)/sqrt(3.0d0*(grhog+grhonu))

       tauosc= (a_osc*0.01d0/adotrad)/om
       !DM:
       arad=((a_osc**3.0d0)*((grhog+grhonu)*0.01d0)/(grhoax))**(1.0d0/4.0d0)
       abar=((a_osc**3.0d0)*((grhob)*0.01d0)/(grhoax))**(1.0d0/3.0d0)
       taueq= (-1.0d0+dsqrt(1.0d0+CP%aeq*0.01d0/adotrad))/om/2.0d0

       !old working start time for other modes
       taustart=min(taustart,tauosc,taueq, 0.3_dl*CP%tau_osc) !RL added CP%tau_osc as a constraint (the 0.001 prefactors in the times are subject to test and change)
    end if

    GetTauStart=taustart
    !RL putting a flag to check the perturbation taustart is after the background taustart
    if (taustart .le. DeltaTime(0._dl, 10.0**(loga_table(1)), 1.0d-10)) then
       write(*, *) 'WARNING: perturbation taustart is before background taustart, &
            &the initial condition may not be valid. perttaustart, bgtaustart:', &
            &taustart, DeltaTime(0._dl, 10.0**(loga_table(1)), 1.0d-10)
    end if

Semantics: per-k integration must start (a) well before the axion field starts
oscillating — ``tauosc`` is the conformal time when a = a_osc/100 estimated from the
matter-radiation analytic solution with ``om = (grhob+grhoc)/sqrt(3*(grhog+grhonu))``
[Mpc^-1] and ``adotrad`` (radiation-era da/dtau) — (b) well before matter-radiation
equality (``taueq``, a = CP%aeq/100 via the exact RD+MD solution
tau = (sqrt(1+a/a_eq-ish)-1)/(om/2) form), and (c) no later than ``0.3*CP%tau_osc``.
``arad``, ``abar``, ``taurad``, ``taubar``, ``tauax_init`` are computed/declared but UNUSED
(dead leftovers from D. Grin; do not port). Units: all taus in Mpc; grho* are the
8*pi*G*rho_i*a^4 (radiation) / a^3 (matter) constants of CAMB Nov13.
Cross-file inputs: ``a_osc``, ``grhoax``, ``loga_table`` (module-level, modules.f90 /
axion_background), ``CP%aeq``, ``CP%tau_osc``, ``CP%Omegaax``, and ``DeltaTime(a1,a2,tol)``
with explicit tolerance 1.0d-10.
The mistakenly-named "WARNING flag" check compares the perturbation start to the first
entry of the background log-a table; this is a sanity warning only. **[PHYSICS]**
(port the min() constraint and the warning; drop dead variables).
NOTE the inconsistency to resolve at port time: comment says 0.001 prefactor but the
code uses ``*0.01d0`` for a_osc/aeq — keep code behavior (``0.01``).

1.4 ``GetSourceMem`` / ``FreeSourceMem`` (862-896) **[PHYSICS allocation]**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: fortran

    if (CP%tau_osc .lt. CP%tau0) then
       allocate(deltaBCSrc(Evolve_q%npoints,SourceNum,2)) !RL 090323, let the 3rd dimension be 2 for now (the terms that are multiplied by J_l instead of J_ldot before and after the switch)
       allocate(dd_deltaBCSrc(Evolve_q%npoints,SourceNum,2))
       deltaBCSrc = 0
       dd_deltaBCSrc = 0
    end if

and in FreeSourceMem (893):

.. code-block:: fortran

    if (allocated(deltaBCSrc))deallocate(deltaBCSrc, dd_deltaBCSrc) !RL 090323

Gate: boundary machinery only exists when the switch happens before today
(``CP%tau_osc < CP%tau0``); for very light axions that never oscillate it is fully bypassed.

1.5 ``CalcScalarSources`` (1100-1224)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- The OLDCAMB ``fixq`` debug block (writes ``evolve_q005.txt`` and stops) is fully commented
  out (lines 1130-1150). **[COSMETIC]**
- 1155: ``tol1=tol/exp(AccuracyBoost-1) !* 2e-6_dl !RL`` — value unchanged. **[COSMETIC]**
- 1207: transfer-function output call gained a ``tau`` argument:

  .. code-block:: fortran

      call outtransf(EV,y,tau, MT%TransferData(:,EV%q_ix,itf))

  (OLDCAMB: ``call outtransf(EV,y, MT%TransferData(:,EV%q_ix,itf))``).
  **[PLUMBING -> PHYSICS interface]**: ``outtransf`` lives in equations_ppf.f90 and needs the
  evaluation time to reconstruct the axion fluid/field density contrast for the new
  transfer columns (see equations report). The modern-CAMB equivalent is
  ``outtransf(EV, y, tau, Arr)`` in equations.f90 — modern CAMB **already passes tau**
  (signature ``subroutine outtransf(EV, y, tau, Arr)`` exists in CAMB 1.6.7), so this is
  free in the port.
- 1219-1221, after the time-step loop — capture of the switch boundary amplitudes for
  this k:

  .. code-block:: fortran

      if (allocated(deltaBCSrc)) then
         deltaBCSrc(EV%q_ix,1,:)=EV%metric_delta !RL 090323 adding the delta function due to source discontinuity
      end if

  **[PHYSICS]**. ``EV%metric_delta(2)`` is set inside equations_ppf.f90 at the KG->fluid
  switch (equations_ppf.f90:449 declaration; filled at ~1118-1140 of equations_ppf.f90 in
  the switch routine; see the equations report). Only source index 1 (temperature) is
  filled; indices 2..SourceNum of deltaBCSrc stay 0 (E-pol and lensing get no boundary
  term).

1.6 ``GetTransfer`` (1340-1361)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- 1347: ``atol=tol/exp(AccuracyBoost-1) !default`` — unchanged value. **[COSMETIC]**
- 1358: same signature change as 1.5: ``call outtransf(EV,y,tau,MT%TransferData(:,EV%q_ix,i))``. **[PLUMBING]**

1.7 ``InitSourceInterpolation`` (1426-1455) — spline deltaBCSrc over k **[PHYSICS]**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After the existing Src spline loop (unchanged), new parallel loop:

.. code-block:: fortran

    !RL 090323 added interpolation across ks for deltaBCSrc as well. Not varying over time
    !Note that spline only takes in 1-d arrays
    ! This needs to be parallelized separately since they don't have the same dimensions
    if (allocated(deltaBCSrc)) then
       !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(STATIC), PRIVATE(i,j), SHARED(Evolve_q)
       do i = 1, 2
          do j = 1, SourceNum
             call spline(Evolve_q%points, deltaBCSrc(1,j,i), Evolve_q%npoints, spl_large, spl_large, dd_deltaBCSrc(1,j,i))
          end do
       end do
       !$OMP END PARAllEl DO
    end if

(natural-ish spline with ``spl_large`` end conditions, same convention as the Src spline).

1.8 ``SourceToTransfers`` (590-627) + ``IntegrationVars_Init`` (1655-1663) — IV lifecycle **[PLUMBING]**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: fortran

    ! Initialize allocation status !RL 122124
    IV%Source_q => null()
    IV%ddSource_q => null()
    IV%metricdeltas_q => null()

    allocate(IV%Source_q(TimeSteps%npoints,SourceNum))
    if (CP%tau_osc .lt. CP%tau0) then
       allocate(IV%metricdeltas_q(2,SourceNum)) !RL 090323
    end if !RL 122124
    ...
    deallocate(IV%Source_q)
    if (associated(IV%metricdeltas_q)) then
       deallocate(IV%metricdeltas_q) !RL 090323 - "associated" is for pointers
    end if

and in IntegrationVars_Init (1661):

.. code-block:: fortran

    if (associated(IV%metricdeltas_q))IV%metricdeltas_q(:,:) = 0 !RL090323

Intent: each integration k gets a 2xSourceNum buffer of switch-boundary amplitudes.
In modern CAMB the analogue of IntegrationVars lives in cmbmain (still), so the
mechanics carry over directly; nullify-before-allocate was a bug fix (RL 122124) against
undefined pointer states with gfortran. **Thread-safety**: IV is private per k iteration;
deltaBCSrc/dd_deltaBCSrc are read-only during integration — safe under the existing
``!$OMP PARALLEL DO`` at line 343.

1.9 ``InterpolateSources`` (1535-1652) — k-interpolation of boundary amplitudes **[PHYSICS]**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Inserted after the standard Src k-interpolation (which already computed
klo/khi/ho/a0/b0/ho2o6/a03/b03 for IV%q on the Evolve_q grid):

.. code-block:: fortran

    ! Ensure SourceNum, klo, and khi are within valid bounds
    if (CP%WantScalars) then
       if (associated(IV%metricdeltas_q) .and. allocated(deltaBCSrc)) then
              do i=1,2
                 IV%metricdeltas_q(i,1:SourceNum) = a0*deltaBCSrc(klo,1:SourceNum,i)+ &
                      b0*deltaBCSrc(khi,1:SourceNum,i) + (a03*dd_deltaBCSrc(klo,1:SourceNum,i) &
                      +b03*dd_deltaBCSrc(khi,1:SourceNum,i))*ho2o6
              end do
       end if
    end if

i.e. the same cubic-spline-in-k evaluation used for the sources themselves.
(Only effective for scalars; tensors/vectors see zeros from IntegrationVars_Init.)

1.10 ``DoFlatIntegration`` (1722-1932) — switch boundary term, flat case **[PHYSICS — core change]**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

New declarations (1733-1739):

.. code-block:: fortran

    integer bes_ix, besix_osc, n, bes_index(IV%SourceSteps)
    real(dl) xf_osc, fac_osc, aa_osc, Jl_osc !RL090623 delt_osc,
    real(dl) dotJl_osc !RL091223
    integer i_output !RL092023
    logical source_output_done !RL092023
    real(dl) l_float, xf_n !RL 092123
    source_output_done = .False. !RL092023

(``i_output``, ``source_output_done``, ``l_float``, ``xf_n`` are debug leftovers, unused in live
code — do not port.)

Trivial precision tweak (1767): ``fac(j)=fac(j)**2._dl*aa(j)/6._dl`` (was ``fac(j)**2*aa(j)/6``).
**[COSMETIC]** (identical numerics).

The core addition, inside the l-loop (``do j=1,max_bessels_l_index``), AFTER the existing
time integral (both the SourceNum==2 fast path and the SourceNum==3 path with Limber)
and BEFORE accumulation into Delta_p_l_k (verbatim, lines 1863-1927 with debug comments
elided):

.. code-block:: fortran

    !RL 090623 adding the boundary conditions at the switch if applicable
    if (associated(IV%metricdeltas_q) .and. allocated(deltaBCSrc)) then
       !I only need to precompute the interpolation parameters once for this q
       xf_osc=abs(IV%q*(CP%tau0-CP%tau_osc))
       besix_osc=Ranges_indexOf(BessRanges,xf_osc)
       fac_osc=BessRanges%points(besix_osc+1)-BessRanges%points(besix_osc)
       aa_osc=(BessRanges%points(besix_osc+1)-xf_osc)/fac_osc
       !Calculate analytical form of the derivative of the J_l cubic spline expression before reusing fac_osc
       dotJl_osc = -IV%q*((ajl(besix_osc+1,j) - ajl(besix_osc,j))/fac_osc + &
            &fac_osc*((1._dl - 3._dl*(aa_osc**2._dl))*ajlpr(besix_osc,j) + &
            &(2._dl - 6._dl*aa_osc + 3._dl*(aa_osc**2._dl))*ajlpr(besix_osc+1,j))/6._dl)
       fac_osc=fac_osc**2._dl*aa_osc/6._dl
       !Now obtain the corresponding spline interpolated J_l
       !Note that in this loop j runs over (1,max_bessels_l_index) and doesn't represent TimeSteps anymore
       Jl_osc=aa_osc*ajl(besix_osc,j)+(1-aa_osc)*(ajl(besix_osc+1,j) - ((aa_osc+1) &
            *ajlpr(besix_osc,j)+(2-aa_osc)*ajlpr(besix_osc+1,j))*fac_osc) !cubic spline

       sums(1) = sums(1) + CP%expmmu_tauosc*Jl_osc*(IV%metricdeltas_q(1, 1) &
            &- CP%opac_tauosc*IV%metricdeltas_q(2, 1)*11._dl/10._dl) &
            &+ CP%expmmu_tauosc*dotJl_osc*IV%metricdeltas_q(2, 1) !!RL 090424
    end if

    ThisCT%Delta_p_l_k(1:SourceNum,j,IV%q_ix) = ThisCT%Delta_p_l_k(1:SourceNum,j,IV%q_ix) + sums(1:SourceNum)

Physics: the scalar temperature source has delta(tau-tau_osc) and delta'(tau-tau_osc)
pieces from the discontinuity in metric derivatives at the KG->effective-fluid switch.
Integrating S(tau) j_l(k(tau0-tau)) dtau across the switch produces a J_l(x_osc) term
with amplitude ``metricdeltas_q(1,1)`` and a dJ_l/dtau term with amplitude
``metricdeltas_q(2,1)``; the dJ_l/dtau is implemented analytically as the derivative of
the cubic-spline interpolant (chain rule dx/dtau = -q gives the leading ``-IV%q*`` factor).
The ``- CP%opac_tauosc*IV%metricdeltas_q(2,1)*11._dl/10._dl`` piece is a compensating term
(opacity correction at the switch; the 11/10 factor is hard-coded — quote and port as-is,
flagged RL 090424). Visibility factors: ``CP%expmmu_tauosc`` = e^{-kappa}(tau_osc),
``CP%opac_tauosc`` = kappa_dot(tau_osc), both precomputed via
``call ThermoSplineOut(CP%tau_osc, CP%opac_tauosc, CP%expmmu_tauosc)`` (modules.f90:2890 —
see modules report). Boundary term is added once per (l, k); it is NOT multiplied by the
time-step measure ``TimeSteps%dpoints`` (it is a point contribution, comment in source says
exactly this).
Note: the boundary term is added regardless of the ``DoInt``/Limber branch taken for the
regular integral, and only to sums(1) (temperature).

1.11 ``DoRangeInt`` (2115-2334) — switch boundary term, curved case **[PHYSICS]**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

New declaration (2139): ``real(dl) chi_osc, sh_osc, ujl_osc, dotujl_osc !RL 041824``

At the very end, after ``out = (out - sources*ujl/2)*delchi*CP%r`` (verbatim, 2315-2329):

.. code-block:: fortran

    !RL testing boundary condition 043024------
    !Note that this is boundary condition and should not be multiplied by the measure, hence it should be after the measure is multiplied
    if (CP%tau_osc .lt. CP%tau0 .and. (TimeSteps%points(nstart)-CP%tau_osc)*(TimeSteps%points(nend)-CP%tau_osc) .le. 0._dl) then
       chi_osc=(CP%tau0-CP%tau_osc)/CP%r
       call USpherBesselWithDeriv(CP%closed,chi_osc,l,nu,ujl_osc,dotujl_osc)
       sh_osc = rofChi(chi_osc)
       ujl_osc = ujl_osc/sh_osc
       dotujl_osc = -(dotujl_osc - ujl_osc*cosfunc(chi_osc))/(sh_osc*CP%r)
       out(1) = out(1) + CP%expmmu_tauosc*ujl_osc*(IV%metricdeltas_q(1, 1) &
            &- CP%opac_tauosc*IV%metricdeltas_q(2, 1)*11._dl/10._dl) &
            &+ CP%expmmu_tauosc*dotujl_osc*IV%metricdeltas_q(2, 1) !RL090424
    end if

Same physics as 1.10 with hyperspherical Bessels: ``USpherBesselWithDeriv`` returns the
raw ujl and its chi-derivative; division by ``rofChi`` and the
``-(dotujl - ujl*cosfunc(chi))/(sh*CP%r)`` conversion produce the normalized ujl and its
conformal-TIME derivative. Trigger condition: the [nstart,nend] sub-range being
integrated straddles tau_osc.
**PORT RISK**: if tau_osc coincides exactly with a TimeSteps range boundary, the ``.le. 0``
product test could fire in two adjacent ranges (double-count). Improbable but worth a
strict-inequality guard or a "added once per (l,k)" flag in the port.
**Note**: ``DoRangeIntTensor`` is untouched (no tensor boundary term).

1.12 ``InterpolateCls`` (2835-2882) — template interpolation removed **[ACCURACY/PHYSICS-motivated, must port]**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: fortran

    !! RH NB WE HAVE REMOVED THE INTERPOLATION
    do in=1,CP%InitPower%nn
       if (CP%WantScalars) then
          do i = C_Temp, C_last
             call InterpolateClArr(CTransS%ls,iCl_scalar(1,i,in),Cl_scalar(lmin, in, i), &
                  CTransS%ls%l0)

OLDCAMB called ``InterpolateClArrTemplated(CTransS%ls, ..., CTransS%ls%l0, i)``, which
blends the sparse-l interpolation with the LCDM ``highL_CL_template`` shape. With axion
physics altering the high-l damping tail, the LCDM template biases the interpolation, so
AxiECAMB falls back to plain cubic interpolation over the computed l-samples.
**Modern-CAMB mapping**: CAMB 1.6.7 ``results.f90`` ``CalcCls``/``InterpolateCls`` uses
``InterpolateClArrTemplated`` when ``CP%InitPower`` is non-standard? — concretely it calls
``this%CLData%InterpolateCls`` -> ``InterpolateClArrTemplated`` whenever
``use_spline_template`` is true. The port equivalent is to set
``CP%Accuracy / use_spline_template = .false.`` (modern CAMB exposes
``lSampleBoost``/``use_cl_spline_template``) or replicate the plain-interp call. If left
templated, lensed TT/EE residuals at high l will silently regress for large axion
fractions. Vector/tensor paths unchanged.

1.13 Things in scope that did NOT change (verified, whitespace-insensitive)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``InitTransfer`` (629-768): transfer-function k-grid construction byte-identical
  (q_switch_lowk1=0.7/taurst, dlog_lowk1=2*boost, q_switch_lowk=8/taurst, dlog_lowk=8*boost
  (\*2.5 HighAcc), q_switch_osc=min(kmax,30/taurst), d_osc=200*boost (\*1.8 HighAcc),
  q_switch_highk=min(kmax,60 or 90/taurst), dlog_osc=17*boost, dlog_highk=3*boost,
  amin=5e-5). Only commented "RL hacking for test" qmin lines added. **[COSMETIC]**
- ``SetkValuesForSources`` (975-1035): source k sampling identical;
  ``SourceAccuracyBoost = AccuracyBoost!*10._dl!*CP%dfac/10._dl !RL 102123`` — multiplier
  commented out, i.e. NOT active. **[COSMETIC]** (but evidence they experimented with
  k-sampling boosts tied to the axion ``CP%dfac`` switch parameter).
- ``SetkValuesForInt`` (1458-…): ``IntSampleBoost=AccuracyBoost!*10._dl !RL 020624`` — same story. **[COSMETIC]**
- ``GetLimberTransfers``, ``CalcLimberScalCls``, ``CalcScalCls``, ``CalcScalCls2``, ``CalcTensCls``,
  ``CalcVecCls``, ``ClTransferToCl``, ``CalcTensorSources``, ``CalcVectorSources``, ``TransferOut``,
  ``MakeNonlinearSources`` (note: nonlinear path), ``InitVars``, ``DoSourcek``,
  ``DoSourceIntegration``, ``IntegrateSourcesBessels``, ``DoRangeIntTensor``, ``GetInitPowerArrayVec``:
  zero non-comment changes.
- **Time sampling**: ``SetTimeSteps`` is NOT in cmbmain.f90 in Nov13 — it lives in
  modules.f90 (OLDCAMB modules.f90:2608). The extra time steps around tau_osc
  (``dtauosc = CP%tau_osc/int(...*CP%dfac)`` and ``Ranges_Add_delta(TimeSteps, ...)`` around
  the switch, AxiECAMB modules.f90 ~2980-3060) belong to the modules.f90 report —
  cross-reference, do not double-port.
- ALens scaling of C_phi (CalcScalCls lines 2635-2640) — already in OLDCAMB, unchanged.

2. cmbmainOMP.f90 vs cmbmain.f90 — **[OBSOLETE — do not port, do not use as reference]**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Facts established:

- Not compiled: absent from CAMBOBJ, no rule, never referenced.
- It is an EARLIER snapshot of the AxiECAMB cmbmain, not a "more parallel" variant:

  - Contains the axion GetTauStart block but the OLDER form:
    ``taustart=min(taustart,tauosc,taueq)`` (no ``0.3_dl*CP%tau_osc`` cap, no
    background-taustart warning).
  - Contains the ``call outtransf(EV,y,tau, ...)`` 4-argument form (both call sites).
  - **Lacks the entire switch boundary-condition machinery**: zero occurrences of
    ``metricdeltas_q``/``deltaBCSrc``/``tau_osc`` (except one comment) — i.e. pre-RL-090323.
  - Has the RH cpu_time/print timing probes ACTIVE (uncommented) — it would spam stdout.
  - ``SourceAccuracyBoost = AccuracyBoost`` without the commented multipliers.
  - OMP directives are the same set as OLDCAMB/cmbmain.f90 (no extra parallelism).

- The 4065-line diff between the two is therefore: live-vs-commented timing prints,
  the missing boundary-condition code, the older GetTauStart, and mass re-indentation.

Conclusion: cmbmainOMP.f90 is a stale development snapshot. Modern CAMB has its own
threading; nothing here is needed. If the port wants a reference, use ONLY cmbmain.f90.

3. camb.f90
~~~~~~~~~~~~

- ``CAMB_GetResults`` (100-251):

  1) ``type(CAMBparams), allocatable :: P !RL 111323`` + ``allocate(P)`` at entry +
     ``if (allocated(P)) deallocate(P)`` at exit. Reason: AxiECAMB's CAMBparams was at one
     point bloated with ``real(dl), dimension(100000,1,6)`` scratch arrays (now commented
     out in modules.f90), so stack copies of ``type(CAMBparams) P`` overflowed; heap
     allocation was the fix. **[OBSOLETE]** — modern CAMB's CAMBparams is a class
     (always heap); and the port must NOT put large tables inside CAMBparams anyway.
  2) **Behavioral change**: the final reset of the module flag was disabled:

     .. code-block:: fortran

         ! RH for axions 20 May at midnight    call_again = .false.

     (OLDCAMB camb.f90:215 ``call_again = .false.``; AxiECAMB camb.f90:250 commented out.
     ``call_again`` is still set ``.false.`` at the start of CAMB_GetResults, line 117, and
     ``.true.`` after each cmbmain pass.) Effect: any ``CAMBParams_Set`` executed AFTER
     CAMB_GetResults but before the next GetResults runs in "call_again" mode, skipping
     re-initialization — presumably to keep the (expensive, state-dependent) axion
     background/thermo results alive for post-processing in the driver. **[PLUMBING]** —
     re-derive intent in modern CAMB: the equivalent is making sure the axion background
     solve is not redone/clobbered between CAMBdata%CalcTransfers and later power
     computations; modern CAMB's explicit CAMBdata state largely removes the need.
     **RISK**: in Nov13 semantics this leaks call_again=.true. to subsequent independent
     runs in the same process (e.g. parameter loops) — verify the port does not need it.
  3) Everything else: commented-out ``cpu_time`` probes and ``!!!write(*, *) 'RL, if N'``
     markers; a re-indent of the ``lens_Cls`` block; no executable change. **[COSMETIC]**

- ``CAMB_TransfersToPowers``: one commented debug print. **[COSMETIC]**
- ``CAMB_ValidateParams``: trailing-whitespace only. **[COSMETIC]**
- **Transfer output columns**: NO change in camb.f90 — the extra axion transfer columns
  are implemented in modules.f90 (``Transfer_*``, ``Transfer_SaveToFiles``) and
  equations_ppf.f90 (``outtransf``); see those reports.

4. subroutines.f90
~~~~~~~~~~~~~~~~~~~

- **dverk: UNCHANGED.** The ODE integrator (subroutines.f90:407 onward) is byte-identical
  modulo whitespace — no error-handling, step-count, or tolerance changes. (AxiECAMB's
  KG/fluid switching is done in equations_ppf's wrappers around dverk, not in dverk.)
- New subroutine ``spline_out`` (AxiECAMB subroutines.f90:8-39), verbatim:

  .. code-block:: fortran

      subroutine spline_out(xarr,yarr,yarr_buff,n,x,y)
        use precision
        integer n,llo_out,lhi_out,midp
        real(dl) xarr(n),yarr(n),yarr_buff(n)
        real(dl) x,y,a0_out,b0_out,ho_out
        ! Written by Dan Grin for the Axion project
        llo_out=1
        lhi_out=n
        do while((lhi_out-llo_out).gt.1)
           midp=(llo_out+lhi_out)/2
           if (xarr(midp).gt.x) then
              lhi_out=midp
           else
              llo_out=midp
           endif
        enddo

        ho_out=xarr(lhi_out)-xarr(llo_out)
        a0_out=(xarr(lhi_out)-x)/ho_out
        b0_out=(x-xarr(llo_out))/ho_out
        y=a0_out*yarr(llo_out)+b0_out*yarr(lhi_out)+((a0_out**3.0d0-a0_out)*&
             yarr_buff(llo_out)+(b0_out**3.0d0-b0_out)*&
             yarr_buff(lhi_out))*ho_out*ho_out /6.d0
      end subroutine spline_out

  Standard NR cubic-spline *evaluation* (binary search + cubic eval; ``yarr_buff`` = y'').
  Called from axion_background.F90, equations_ppf.f90, recfast_axion.F90.
  **[PLUMBING]** — needed by the ported axion code, but in modern CAMB replace with the
  ``Interpolation`` module (``TCubicSpline%Value``) or ``classes``/``splines`` utilities instead
  of adding a free-floating routine.
- ``rombint`` failure handling hardened (subroutines.f90:203-208):

  .. code-block:: fortran

      if (i.gt.MAXITER.and.abs(error).gt.tol)  then
        write(*,*) 'Warning: Rombint failed to converge RH; '
        write (*,*)'integral, error, tol:', rombint,error, tol
        stop
      end if

  OLDCAMB only warned; AxiECAMB **stops the program**. **[ACCURACY/PLUMBING]** — decide
  policy in the port: modern CAMB should raise ``global_error_flag`` /
  ``call GlobalError(...)`` rather than ``stop``. The motivation was presumably that a
  non-converged rombint over the rapidly-oscillating axion background silently corrupts
  results. Port the *intent* (hard failure -> error) not the literal ``stop``.
- ``spline_deriv`` etc.: whitespace only. **[COSMETIC]**

5. lensing.f90 — **[OBSOLETE: already upstream]**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Changes (verbatim):

.. code-block:: fortran

    real(dl) :: ALens_Fiducial = 0._dl
    !Change from zero to set lensing smoothing by scaling amplitude of fiducial template !RL incorporating 013125
    ...
    public ... , ALens_Fiducial !RL 013125 incorporating ALens_Fiducial
    ...
    !RL incorporating ALens_Fiducial 013125
    if (ALens_Fiducial > 0) then
       do l=2, lmax
          sc = (2*l+1)/(4*pi) * 2*pi/(l*(l+1))
          Cphil3(l) =  sc * highL_CL_template(l, C_Phi) * ALens_Fiducial
       end do
    end if

This is a back-port FROM newer CAMB: CAMB 1.6.7 already contains the identical feature
(CAMB/fortran/lensing.f90:51, 68, 266-269). **Nothing to port.** Note: AxiECAMB never
reads ALens_Fiducial from the ini (inidriver only reads ``Alens``), so it is dormant.
Remaining lensing.f90 diffs are a commented debug print and re-indentation. **[COSMETIC]**

6. bessels.f90 — **[COSMETIC]**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Only change: two commented-out debug lines in GenerateBessels
(``!call bjl(1,5.2_dl,bjl_test)`` + a commented write, plus the ``!real(dl) bjl_test``
declaration comment). No functional change. Nothing to port.

7. constants.f90
~~~~~~~~~~~~~~~~~

- ``eV`` renamed to ``elecV`` (same value 1.60217646e-19):

  .. code-block:: fortran

      real(dl), parameter :: elecV = 1.60217646e-19_dl

  Users: modules.f90:502 (neutrino conv factor — pre-existing code, just renamed),
  inidriver_axion.F90:143 (``P%H0_eV = h_P*P%H0_in_Mpc_inv/(Mpc_in_sec*2._dl*const_pi*elecV)``),
  axion_background.F90. **[PLUMBING — do NOT port the rename.]** Modern CAMB keeps ``eV``
  (value updated to the exact SI 1.602176634e-19). The rename was only to avoid a name
  clash with axion-code local variables; in the port, use modern ``eV`` and rename any
  clashing locals instead. Beware the value difference (1.60217646e-19 vs
  1.602176634e-19) when doing regression comparisons at the 1e-8 level.
- New constant:

  .. code-block:: fortran

      real(dl), parameter :: const_rhocrit=(8.0d0*const_pi*G*1.d3/(3.0d0*((1.d7/(MPC_in_sec*c*1.d2))**(2.0d0))))**(-1.0d0)

  = critical density 3H^2/(8 pi G) for H = 100 km/s/Mpc in CGS-flavored units (the 1.d3,
  1.d2, 1.d7 factors convert SI G and H to cgs; result ~1.878e-29 g/cm^3 h^-2 scale).
  **Currently used ONLY inside a commented-out debug line (equations_ppf.f90:2228).**
  Check the axion_background/equations reports before porting; as of this file it is
  dead weight. **[OBSOLETE unless another report finds a live use]**
- New constant:

  .. code-block:: fortran

      real(dl), parameter :: mplanck = 2.435e18 ! the reduced planck mass in GeV

  Live use: inidriver_axion.F90:317 ``P%Hinf = (10**P%Hinf)/mplanck`` (converts input
  log10(H_inflation/GeV) to H_inf in reduced-Planck-mass units for the axion
  isocurvature amplitude). **[PHYSICS]** — port (note: default-real literal ``2.435e18``,
  not ``_dl``; make it ``2.435e18_dl`` in the port). Units: GeV.
- Everything else whitespace. **[COSMETIC]**

8. writefits.f90 — **[OBSOLETE — broken, never compiled]**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The AxiECAMB version is syntactically INVALID: an extra ``end if`` was added at line ~35
with no matching ``if``:

.. code-block:: text

      if (CP%OutputNormalization >=2) then
       fac=1
      else
       fac=OutputDenominator*CP%tcmb**2
      end if
      end if          <-- orphan

plus re-introduction of pre-Nov13 COBE-normalization code referencing a constant that no
longer exists anywhere in AxiECAMB (``outCOBE``):

.. code-block:: fortran

      if (CP%OutputNormalization == outCOBE) then
         unitstr='Kelvin-squared'
      else
         unitstr='unknown'
      end if
     ...
     COBEnorm = CP%outputNormalization==outCOBE

writefits.f90 is only compiled by the ``camb_fits`` target (Makefile_main:84-85), which is
unusable (FITSDIR/HEALPIX commented out in Makefile). Someone pasted chunks of an older
CAMB writefits and never compiled it. **Ignore entirely for the port** (modern CAMB has
no writefits).

9. utils.F90 — **[COSMETIC]**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

2396-line diff reduces to exactly 4 non-comment changed lines, all naming of END
statements:

- ``END INTERFACE`` -> ``END INTERFACE CONCAT``
- ``end subroutine`` -> ``end subroutine MpiQuietWait``

Everything else is whitespace/comment reflow. Nothing to port. (Modern CAMB replaced
utils.F90 with MiscUtils/StringUtils/etc. from forutils anyway.)

10. Cross-file contract: what the cmbmain port needs from other reports
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The cmbmain changes compile only if these exist (all defined OUTSIDE this report's files):

- ``CP%tau_osc``, ``CP%a_osc``, ``CP%aeq``, ``CP%Omegaax``, ``CP%dfac``, ``CP%opac_tauosc``,
  ``CP%expmmu_tauosc`` — CAMBparams members (modules.f90:124; opac/expmmu filled by
  ``call ThermoSplineOut(CP%tau_osc, CP%opac_tauosc, CP%expmmu_tauosc)`` at modules.f90:2890).
- module-level ``a_osc``, ``grhoax``, ``loga_table`` (axion background tables; modules.f90 /
  axion_background.F90).
- ``EV%metric_delta(2)`` — filled in equations_ppf.f90 at the switch (decl
  equations_ppf.f90:449; computation ~1118-1140 — see equations report; it already folds
  in gauge/curvature factors, including the curved-space ``EV%Kf(1)`` variants).
- ``outtransf(EV, y, tau, Arr)`` — 4-arg signature (modern CAMB already matches).
- ``DeltaTime(a1, a2, tol)`` with explicit tolerance (modern CAMB: ``DeltaTime(state, a1, a2, tol)``).
- ``ThermoSplineOut`` (modules report).
- Extra TimeSteps refinement around tau_osc: in ``SetTimeSteps`` (modules.f90, other report).
- ``spline_out`` (this report, section 4) used by axion_background/recfast_axion/equations_ppf.

11. Surprises / risks for the port designer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **cmbmainOMP.f90 is a trap**: it looks like "the OMP version" but is an obsolete
   pre-boundary-condition snapshot with debug prints live. The Makefile compiles
   cmbmain.f90. Never consult cmbmainOMP for physics.
2. **The switch boundary term is the only genuinely new physics in cmbmain** (sections
   1.4-1.11): a J_l and dJ_l/dtau point contribution at tau_osc added to the temperature
   transfer integral, in both flat (spline-interpolated j_l) and curved
   (USpherBesselWithDeriv) branches, with the hard-coded ``11/10`` opacity-compensation
   factor. In modern CAMB the corresponding insertion points are
   ``TimeSourcesToCl``/``DoFlatIntegration`` and ``DoRangeInt`` in fortran/cmbmain.f90 (still
   structurally similar), with ``BessRanges/ajl/ajlpr`` -> ``SpherBessels`` module state.
3. **Template C_l interpolation disabled** (1.12): plain ``InterpolateClArr`` instead of
   ``InterpolateClArrTemplated`` for scalars. Modern CAMB equivalent: disable the
   high-L template spline (results.f90 ``use_spline_template`` / set
   ``CTrans%ls`` interpolation accordingly) — otherwise high-l C_l for large axion
   fractions inherit LCDM template shape errors.
4. **GetTauStart cap 0.3*tau_osc** ensures every k starts in the KG regime — if the
   port keeps modern CAMB's ``GetTauStart`` structure, add the axion min() clause AND the
   background-table start check; otherwise initial conditions are applied in the wrong
   regime for high-k modes with early-oscillating (heavy) axions.
5. **rombint now stops on non-convergence** — convert to modern error handling, but do
   not drop it silently: it may be masking real convergence issues in the axion
   integrals.
6. **call_again = .false. removal in CAMB_GetResults** (3.2) is a subtle global-state
   hack ("RH for axions 20 May at midnight"); understand what the modern driver flow
   needs before replicating — most likely nothing, since CAMBdata holds state explicitly.
7. **DoRangeInt double-count edge case** (1.11): ``.le. 0`` straddle test can in principle
   fire in two adjacent ranges if tau_osc lands exactly on a range boundary.
8. **eV constant value drift**: AxiECAMB elecV=1.60217646e-19 vs modern
   eV=1.602176634e-19 — irrelevant physically, relevant for bit-level regression tests of
   m_ax conversions.
9. Many "accuracy knobs" the authors experimented with (SourceAccuracyBoost*10,
   IntSampleBoost*10, tol1 scalings, qmin hacks) are ALL commented out — current AxiECAMB
   runs with stock Nov13 sampling except the tau_osc TimeSteps refinement (modules.f90)
   and the transfer outputs. Do not port any commented experiment.


Original-code analysis: recombination and reionization (recfast_axion.F90, reionization.f90)
---------------------------------------------------------------------------------------------

Sources analyzed:

- Diff: ``/Users/vivianmiranda/data/research/WayneHu/rayne/.port_analysis/diffs/recfast_axion.diff`` (OLDCAMB/recfast.f90 -> AxiECAMB/recfast_axion.F90)
- Diff: ``/Users/vivianmiranda/data/research/WayneHu/rayne/.port_analysis/diffs/reionization.f90.diff``
- Originals: ``/Users/vivianmiranda/data/research/WayneHu/rayne/AxiECAMB/recfast_axion.F90``, ``/Users/vivianmiranda/data/research/WayneHu/rayne/AxiECAMB/reionization.f90``, ``/Users/vivianmiranda/data/research/WayneHu/rayne/OLDCAMB/recfast.f90``, ``/Users/vivianmiranda/data/research/WayneHu/rayne/OLDCAMB/reionization.f90``
- Modern targets: ``/Users/vivianmiranda/data/research/WayneHu/rayne/CAMB/fortran/recfast.f90``, ``/Users/vivianmiranda/data/research/WayneHu/rayne/CAMB/fortran/reionization.f90``, ``/Users/vivianmiranda/data/research/WayneHu/rayne/CAMB/fortran/results.f90``
- Support: ``/Users/vivianmiranda/data/research/WayneHu/rayne/AxiECAMB/equations_ppf.f90`` (grhoax_frac, dtauda), ``/Users/vivianmiranda/data/research/WayneHu/rayne/AxiECAMB/modules.f90`` (call sites), ``/Users/vivianmiranda/data/research/WayneHu/rayne/AxiECAMB/axion_background.F90`` (omegar definition)

**Headline:** The 2929-line recfast diff is ~90% whitespace/re-indentation plus a verbatim private copy of the DVERK integrator (``recdverk``) added solely to thread extra arguments. The genuine physics content is small: one axion term in ``dHdz`` (used only in the Tmat tight-coupling smoothing term), a corrected ``OmegaK``, a replaced ``z_eq``, a loosened ODE tolerance, and integral-splitting at the KG->EFA switch scale factor ``a_osc`` in reionization timing integrals.

1. WHY recfast needed modification: what Nov13 hardcoded vs what AxiECAMB changed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Important nuance: Nov13 recfast did NOT hardcode the main H(z).** The Hubble rate used in *all* recombination rate equations (Peebles K coefficient, f(1), f(2), f(3) denominators, Sobolev optical depths) was already obtained from the full CAMB background via the external ``dtauda``:

OLDCAMB/recfast.f90:828 (unchanged at AxiECAMB/recfast_axion.F90:877):

.. code-block:: fortran

    Hz = 1/dtauda(1/(1._dl+z))*(1._dl+z)**2/MPC_in_sec

Since AxiECAMB's ``dtauda`` (equations_ppf.f90:299-335) adds the axion density,

.. code-block:: fortran

    grhoa2 = grhoa2 + grhoax_frac(a)*grhom*(a2**2._dl)
    dtauda=sqrt(3._dl/grhoa2)

``Hz`` inside recfast is automatically axion-aware with **no change to recfast itself**.

What Nov13 *did* hardcode analytically (all "only used for approximations where small effect", per RECDATA comment):

1. OLDCAMB/recfast.f90:528-529:

   .. code-block:: fortran

       OmegaT=OmegaC+OmegaB            !total dark matter + baryons
       OmegaK=1.d0-OmegaT-OmegaV       !curvature

   (ignores radiation, neutrinos — comment "DM11: notes that the neutrinos are not included here. Gives wrong curvature.")

2. OLDCAMB/recfast.f90:548-549 (analytic z_eq for 3 massless neutrinos):

   .. code-block:: fortran

       z_eq = (3.d0*(HO*C)**2/(8.d0*Pi*G*a_rad*(1.d0+fnu)*Tnow**4))*(OmegaB+OmegaC)
       z_eq = z_eq - 1.d0

   with ``fnu = (21.d0/8.d0)*(4.d0/11.d0)**(4.d0/3.d0)``.

3. OLDCAMB/recfast.f90:976-977 — the analytic dH/dz used ONLY in the "additional term to smooth transition to Tmat evolution (suggested by Adam Moss)" inside the tightly-coupled branch ``timeTh < H_frac*timeH`` of subroutine ION:

   .. code-block:: fortran

       dHdz = (HO**2/2.d0/Hz)*(4.d0*(1.d0+z)**3/(1.d0+z_eq)*OmegaT &
        + 3.d0*OmegaT*(1.d0+z)**2 + 2.d0*OmegaK*(1.d0+z) )

   (matter + radiation-via-z_eq + curvature only; no Lambda, no exotic components).

4. OLDCAMB/recfast.f90:926 (EdS Hubble-time estimate, only used for the Compton/Hubble timescale switch; **unchanged in AxiECAMB**, recfast_axion.F90:975):

   .. code-block:: fortran

       timeH=2./(3.*HO*(1._dl+z)**1.5)      !Hubble time

AxiECAMB replacements (the actual axion physics)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**[PHYSICS] dHdz with axion density derivative** — AxiECAMB/recfast_axion.F90:1024-1100 (inside ``ION``). The old expression is kept as a comment (lines 1029-1031). New code:

.. code-block:: fortran

    sfac=1.0d0/(1.0d0+z)
    ...
    dorp = grhoax_frac(sfac) !RL
    ...
    deriv_eps=1.d-3*real(sfac)
    ...
    if (sfac .ge. aosc) then !RL 121524 forward stepping if the present scale factor is already after the switch
       dorpa = grhoax_frac(sfac+deriv_eps)
       dorpa=(dorpa-dorp)/(deriv_eps)
    else
       dorpa = grhoax_frac(sfac-deriv_eps) !RL 121524 - backstepping if the present scale factor is before the switch, to avoid cases where the derivative straddles the switch. This is not the best fix; the best fix is to remove dH/dz altogether since it's in fact redundant, though doing so is more integrated
       dorpa=(dorp-dorpa)/(deriv_eps)
    end if

    dHdz = (HO**2/2.d0/Hz)*(4.d0*((1.d0+z)**3)*omegar &
         + 3.d0*OmegaT*(1.d0+z)**2 + 2.d0*OmegaK*(1.d0+z)&
         &-dorpa/((1.0d0+z)**2.0d0))

Interpretation (units): ``grhoax_frac(a) = rho_ax(a)/rho_crit,0`` (dimensionless; comment at recfast_axion.F90:1055-1057: "this is grhoax_table_internal=rho_ax/rho_crit^today so already in the form needed for recfast's normalization conventions"). ``dorpa = d(rho_ax/rho_crit,0)/da`` by one-sided finite difference with step ``deriv_eps = 1e-3*a``, one-sided *away from* ``a_osc`` so the stencil never straddles the KG->EFA switch. The ``-dorpa/(1+z)^2`` term is exactly ``d[f_ax(a(z))]/dz = (df_ax/da)(da/dz) = -a^2 df_ax/da`` inserted into ``dHdz = HO^2/(2Hz) * d/dz[ E^2(z) ]``. The radiation term uses ``omegar`` directly instead of ``OmegaT/(1+z_eq)``.

``dHdz`` is used (unchanged form) in the tightly-coupled Tmat derivative, recfast_axion.F90:1104-1107:

.. code-block:: fortran

    epsilon = Hz*(1.d0+x+fHe)/(CT*Trad**3*x)
    f(3) = Tnow &
         + epsilon*((1.d0+fHe)/(1.d0+fHe+x))*((f(1)+fHe*f(2))/x) &
         - epsilon* dHdz/Hz + 3.0d0*epsilon/(1.d0+z)

The supporting background function (defined in AxiECAMB/equations_ppf.f90:257-296, used by recfast via ``external grhoax_frac``):

.. code-block:: fortran

    if (a .lt. CP%a_osc) then
        ...
        call spline_out(loga_table,rhoaxh2ovrhom_logtable,rhoaxh2ovrhom_logtable_buff,ntable,dlog10(a),grhoaxh2_ov_grhom)
        ...
       grhoax_frac = (10._dl**grhoaxh2_ov_grhom)/(CP%H0**2.0d0/1.0d4)
    else
       wcorr_coeff = CP%ahosc_ETA*CP%a_osc/((CP%ma/CP%H0_eV)*(CP%H0/100.0d0)) !RL082924
       grhoax_frac=(CP%rhorefp_ovh2)*((CP%a_osc/a)**3.0d0)*dexp((wcorr_coeff**2.0d0)*3.0d0*CP%wEFA_c*(1.0d0/(a2**2.0d0) &
            &- 1.0d0/(CP%a_osc**4.0d0))/4.0d0)
    endif

i.e. KG-table spline of log10(rho_ax h^2 / rho_crit) before the switch; a^-3 scaling with an EFA equation-of-state correction exponential after the switch. (This function belongs to the background/equations report; recfast just consumes it.)

**[PHYSICS/ACCURACY] Self-consistent OmegaK** — AxiECAMB/recfast_axion.F90:562-563:

.. code-block:: fortran

    OmegaT=OmegaC+OmegaB        !total dark matter + baryons + axions
    OmegaK=1.d0-OmegaT-OmegaAx-OmegaV-omegar-Omegan      !curvature

(Note: despite the comment, ``OmegaT`` itself remains CDM+baryons only; the axion enters via the separate ``-dorpa`` term in dHdz and via OmegaK.) Comment block: "DG15: Have added massive and massless neutrino+photon contribution. Does not make a huge difference."

**[PHYSICS, but vestigial] z_eq from axion background solver** — AxiECAMB/recfast_axion.F90:583-586:

.. code-block:: fortran

    !z_eq = !(3.d0*(HO*C)**2/(8.d0*Pi*G*a_rad*(1.d0+fnu)*Tnow**4))*(OmegaB+OmegaC)
    z_eq=aeq**(-1.0d0)
    z_eq = z_eq - 1.d0

``aeq`` is the matter-radiation equality scale factor computed by ``axion_background.F90`` (counts the ULA as matter when appropriate). **However**, since AxiECAMB's new ``dHdz`` uses ``omegar`` directly instead of ``OmegaT/(1+z_eq)``, ``z_eq`` is no longer referenced anywhere in recfast_axion's ION — this assignment is effectively dead code inside recfast (it only fills the RECDATA module variable).

2. Every other change in recfast_axion.F90
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[PLUMBING] Signature and argument-threading changes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``Recombination_init`` signature (recfast_axion.F90:470-471):

  .. code-block:: fortran

      subroutine Recombination_init(Recomb, OmegaC, OmegaB, Omegan, Omegav, &
           h0inp,tcmb,yp,OmegaAx,omegar,aeq, aosc)

  vs OLD ``(Recomb, OmegaC, OmegaB, Omegan, Omegav, h0inp,tcmb,yp, nnu)``. Drops the (unused) optional ``nnu``; adds:

  - ``OmegaAx`` = Omega_axion today (used only in OmegaK),
  - ``omegar`` = Omega_photons + Omega_massless-neutrinos (set as ``Params%omegar=Params%omegah2_rad/hsq`` at axion_background.F90:1168; ``omegah2_rad`` = photon + massless-nu energy density /h^2, axion_background.F90:293-308),
  - ``aeq`` = equality scale factor from the axion background solver (used only for the vestigial z_eq),
  - ``aosc`` = scale factor where m = dfac*H (KG->EFA switch) — used only to pick the finite-difference direction in ION.

  Call site: AxiECAMB/modules.f90:2500-2501:

  .. code-block:: fortran

      call Recombination_Init(CP%Recomb, CP%omegac,CP%omegab,CP%Omegan,&
           CP%Omegav,CP%h0,CP%tcmb,CP%yhe,CP%omegaax,CP%omegar,CP%aeq, CP%a_osc)

- ``ION`` signature (recfast_axion.F90:819): ``subroutine ION(Recomb,OmegaAx,omegar,aeq,aosc,Ndim,z,Y,f)``. Inside ION only ``omegar`` (dHdz) and ``aosc`` (FD direction) are actually used; ``OmegaAx`` and ``aeq`` are dead parameters there.

- ``recdverk`` (recfast_axion.F90:1434-end): a full verbatim copy of the DVERK Runge-Kutta integrator whose only difference from ``dverk`` is that ``fcn`` is invoked as ``call fcn(EV,OmegaAx,omegar,aeq,aosc,n, x, y, w(1,1))`` so the axion scalars reach ION. Comment: "recombination tailored verk that properly passes around splines / put in to avoid a make catastrophe". Replaces the two ``call DVERK(Recomb,3,ION,...)`` / ``call DVERK(Recomb,nw,ION,...)`` calls with ``call recdverk(Recomb,OmegaAx,omegar,aeq,aosc,3,ION,zstart,y,zend,tol,ind,cw,nw,w)`` (recfast_axion.F90:707, 713).

  - **Latent bug (do not replicate):** the 7 inner-stage calls (recfast_axion.F90:2000, 2006, 2013, 2021, 2030, 2039, 2049) pass ``real(OmegaAx)`` — a default-precision (single) conversion — instead of ``OmegaAx``. Harmless only because ``OmegaAx`` is never used inside ION.

[ACCURACY] ODE tolerance loosened 100x
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

recfast_axion.F90:515:

.. code-block:: fortran

    real(dl), parameter :: tol=1.D-3                !Tolerance for R-K

vs OLDCAMB/recfast.f90:487 ``tol=1.D-5``. **Undocumented.** Most plausibly a workaround for RK step-failures when integrating across the kink at ``a_osc`` (the same kink that forced the FD-direction switch and the reionization rombint splitting). Modern CAMB recfast uses ``tol=1.D-5`` (CAMB/fortran/recfast.f90:574). Port decision: try keeping 1.D-5; only loosen (or better, make tolerance adaptive/split at the switch) if the EFA kink actually causes failures. Flag as an accuracy regression risk if blindly ported.

[COSMETIC]/dead-code items
^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``Recombination_Name`` changed ``'Recfast_1.5.2'`` -> ``'Recfast_1.5'`` (string only).
- Commented alternative grids ``!!integer, parameter :: Nz=5000`` / ``!!Nz=20000 !RL 112123`` next to the unchanged ``Nz=10000`` (zinitial=1e4, zfinal=0, delta_z unchanged). **No actual z-grid change.**
- Commented linear-interpolation experiment in ``Recombination_xe`` ("``!!Recombination_xe=az*xrec(ilo)+bz*xrec(ihi) !RL trying linear interpolation``") — spline interpolation retained.
- Numerous commented debug writes (``write(11162302,...)``, ``write(11162303,...)``, ``write(12212311,...)``), commented spline-table plumbing from an earlier implementation (loga_table/grhoax_table/rhorefp_hsq, replaced by ``grhoax_frac`` on 07/31/2023), and a commented earlier w-correction formula ``dorp=rhorefp_hsq*((aosc/sfac)**3.0d0)*dexp((wcorr_coeff**2.0d0)*9.0d0*(1.0d0/(sfac**4.0d0) - 1.0d0/(aosc**4.0d0))/8.0d0)``.
- Whole-file re-indentation (the bulk of the 2929-line diff).

Explicitly UNCHANGED physics (verified against the diff; nothing to port beyond stock recfast 1.5.2)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- All atomic data, PPB/VF/triplet rate fits, He_Boltz overflow guard.
- All fudge factors and defaults: ``RECFAST_fudge_default = 1.14``, ``RECFAST_fudge_default2 = 1.105d0 + 0.02d0``, ``RECFAST_fudge_He_default = 0.86``, ``RECFAST_Heswitch_default = 6``, ``RECFAST_Hswitch_default = .true.``, Gaussian fit parameters (AGauss1=-0.14, AGauss2=0.079, zGauss1=7.28, zGauss2=6.73, wGauss1=0.18, wGauss2=0.33).
- Saha-regime switch redshifts (z>8000, >5000, >3500, x_He>0.99, x_H>0.99/0.985) and ``H_frac = 1D-3`` (recfast_axion.F90:612).
- Peebles K, K_He, Sobolev/continuum-opacity He treatment, triplet channel, f(1), f(2) equations, the non-tightly-coupled f(3) branch, 21cm functions (``kappa_HH_21cm``, ``kappa_eH_21cm``, ``kappa_pH_21cm``), ``dDeltaxe_dtau`` (still uses plain ``K = CK/Hz`` with Hz from dtauda).
- Caching check ``dtauda(0.2352375823_dl) == Last_dtauda`` (works for axions because dtauda changes when axion params change).

3. reionization.f90 changes (AxiECAMB vs OLDCAMB)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[ACCURACY/PHYSICS-tweak] Helium full-reionization start redshift
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

AxiECAMB/reionization.f90:32-33:

.. code-block:: fortran

    !!real(dl) :: helium_fullreion_redshiftstart  = 5._dl !original
    real(dl) :: helium_fullreion_redshiftstart  = 7._dl  !RL tweaking 011023

This widens the redshift below which the second HeII reionization tanh contribution is added to xe (``if (include_helium_fullreion .and. a > (1/(1+ helium_fullreion_redshiftstart)))``, line 89). Hardcoded tweak, not axion physics per se (helium_fullreion_redshift=3.5 and deltaredshift=0.5 unchanged). Modern CAMB 1.6.7 default is ``helium_redshiftstart = 5.5_dl`` (CAMB/fortran/reionization.f90:44) and it is a user-settable parameter — **port as a parameter default override only if needed to reproduce AxiECAMB outputs**; do not hardcode.

[PLUMBING/ACCURACY] tau_start / tau_complete integrals split at a_osc
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``Reionization_Init`` gains an ``aosc`` argument (AxiECAMB/reionization.f90:149: "RL 050225 added aosc since there's no clean way to use DeltaTime here"; called from modules.f90:440 with ``CP%a_osc``). The two conformal-time rombint integrals over ``dtauda`` are split at ``a_osc`` (lines 197-219):

.. code-block:: fortran

    atol = 1d-3
    !RL 050225 reuse deltatime_reion for both time regimes
    if (0._dl .lt. aosc .and. astart .ge. aosc) then
       deltatime_reion=rombint(dtauda,0._dl,aosc*(1._dl-max(atol/100.0_dl,1.d-15)),atol) + rombint(dtauda, aosc, astart, atol)
    else
       deltatime_reion=rombint(dtauda,0._dl,astart,atol)
    end if
    ReionHist%tau_start = max(0.05_dl, deltatime_reion)
    ...
    aend = 1.d0/(1.d0+max(0.d0,Reion%redshift-Reion%delta_redshift*8))
    if (astart .lt. aosc .and. aend .ge. aosc) then
       deltatime_reion=rombint(dtauda,astart,aosc*(1._dl-max(atol/100.0_dl,1.d-15)),atol) + rombint(dtauda, aosc, aend, atol)
    else
       deltatime_reion=rombint(dtauda,astart, aend,atol)
    end if
    ReionHist%tau_complete = min(tau0, ReionHist%tau_start + deltatime_reion)

Reason: ``dtauda`` has a kink/derivative discontinuity at ``a_osc`` (KG -> effective-fluid switch), which degrades/breaks Romberg convergence when the interval straddles the switch. The upper limit of the first sub-integral is nudged just below ``a_osc`` (``aosc*(1 - max(atol/100, 1e-15))``) so neither piece touches the kink. (An attempted cleaner route via a ``DeltaTime_external`` declared at line 57 is commented "RL 050225 failed" — dead code.)

**Note for low-mass axions:** if ``m_ax`` is light enough that ``a_osc`` falls inside/after reionization, both branches matter; for heavy axions ``a_osc << astart`` and only the tau_start branch splits.

Optical depth integral: UNCHANGED (automatically axion-aware)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``Reionization_GetOptDepth`` and its integrand are untouched:

.. code-block:: fortran

    Reionization_doptdepth_dz = Reionization_xe(a)*ThisReionHist%akthom*dtauda(a)

(OLDCAMB/reionization.f90:261; identical in AxiECAMB). The axion background enters only through ``dtauda``. The tanh window (``Reionization_xe`` with ``WindowVarMid``, ``Rionization_zexp = 1.5``) is untouched. ``Reionization_SetParamsForZre``, ``Reionization_SetFromOptDepth``, default fraction, ``Reionization_maxz = 50`` — all untouched.

[PLUMBING — DO NOT PORT] zreFromOptDepth error handling silenced
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

AxiECAMB/reionization.f90:320 and 325-328:

.. code-block:: fortran

    !       if (i>100) call mpiStop('Reionization_zreFromOptDepth: failed to converge')
    ...
    if (abs(tau - Reion%optical_depth) > 0.002) then
     !write (*,*) 'Reionization_zreFromOptDepth: Did not converge to optical depth'
     !write (*,*) 'tau =',tau, 'optical_depth = ', Reion%optical_depth
     !write (*,*) try_t, try_b
     !call mpiStop()
    end if

Hard aborts on bisection non-convergence were commented out (presumably to keep MCMC chains alive in pathological corners of axion parameter space). This silently accepts a wrong zre/tau mapping. Modern CAMB raises a *recoverable* error instead (``call GlobalError('...failed to converge',error_reionization)``, CAMB/fortran/reionization.f90:242,250), which is strictly better — **mark [OBSOLETE]; do not port the silencing.**

4. CRITICAL ASSESSMENT: what survives into CAMB 1.6.7 and what must be ported
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Modern recfast.f90 architecture facts (verified)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Main expansion rate (CAMB/fortran/recfast.f90:882):

  .. code-block:: fortran

      Hz = ainv**2/dtauda(Recomb%State,1/ainv)/MPC_in_sec

  i.e. the **full cosmology-state background**. ``dtauda`` is the state-bound function; once the axion density is added to the modern background (results.f90 ``grho`` machinery / the new axion component), every Hz use in modern recfast is automatically axion-aware. ION already has ``Recomb%State`` via the bound calculator type — no argument threading needed.

- BUT the smoothing-term derivative is **still analytic** (CAMB/fortran/recfast.f90:1030-1031):

  .. code-block:: fortran

      dHdz = (Recomb%HO**2/2.d0/Hz)*(4.d0*ainv**3/(1.d0+Recomb%z_eq)*Recomb%OmegaT &
          + 3.d0*Recomb%OmegaT*ainv**2 + 2.d0*Recomb%OmegaK*ainv )

  with ``Calc%OmegaT=(State%CP%omch2+State%CP%ombh2)/H**2`` (line 611), ``Calc%OmegaK=State%CP%omk`` (line 612), ``Calc%z_eq = State%z_eq`` (line 628), and ``State%z_eq = (grhob+grhoc)/(grhog+grhornomass+sum(grhormass(...))) - 1`` (results.f90:476). **No dark-energy or exotic-component term — the axion term still needs porting.**

- Modern ION evolves ``y(3) = a*Tmat`` (not Tmat); dHdz enters as ``- epsilon*dHdz/(Hz*ainv)`` in ``daTmat_dz`` (recfast.f90:1034-1036). The ported axion term slots into the same ``dHdz`` expression unchanged.

- Modern tol is ``1.D-5`` (recfast.f90:574); Nz default 10000 but now run-time settable (``RECFAST_nz``, recfast.f90:276,369).

Classification of each AxiECAMB recfast change for the port
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Change (AxiECAMB location)
     - Class
     - Port verdict
   * - ``Hz`` from dtauda (recfast_axion.F90:877)
     - —
     - Nothing to do: identical mechanism already in 1.6.7 (recfast.f90:882), generic via State. Becomes axion-aware automatically once the background component exists.
   * - ``dHdz`` axion term ``-dorpa/((1.0d0+z)**2.0d0)`` with one-sided FD of rho_ax/rho_crit,0 away from a_osc (recfast_axion.F90:1058-1100)
     - **[PHYSICS]**
     - **MUST PORT.** Modern dHdz (recfast.f90:1030) has no exotic term. Two options: (a) literal port — add ``- (drho_ax_frac/da)/(1+z)^2`` using the modern axion component's density (and keep the one-sided FD away from a_osc, or use an analytic derivative which the EFA form admits in closed form); (b) cleaner — replace the whole analytic dHdz with a numerical d/dz of ``Hz`` from ``dtauda`` (RL's own comment at recfast_axion.F90:1089 says dH/dz "is in fact redundant"). Effect is confined to the tightly-coupled Tmat smoothing term (small but the term scales with the DM-like axion fraction at recombination: after a_osc the axion contributes ~\ ``3*f_ax*(1+z)^2``-like a matter term).
   * - Radiation term ``4*(1+z)^3*omegar`` replacing ``4*(1+z)^3*OmegaT/(1+z_eq)`` (recfast_axion.F90:1098)
     - [OBSOLETE]
     - Algebraically already handled: with modern ``z_eq=(grhob+grhoc)/grho_rad - 1`` (results.f90:476), ``OmegaT/(1+z_eq) == Omega_rad`` *identically* (OmegaT ∝ grhob+grhoc cancels), independent of axions. Minor bookkeeping difference: modern counts massive nu as radiation in z_eq's denominator; AxiECAMB's omegar excludes massive nu (it moved Omegan into OmegaK). Negligible for this approximation term; no port needed if z_eq form retained.
   * - ``z_eq = aeq^{-1} - 1`` from axion solver (recfast_axion.F90:585-586)
     - [OBSOLETE]
     - Dead inside recfast even in AxiECAMB (z_eq no longer referenced after the omegar substitution). Modern recfast's z_eq use is the identity above. Do not port into recfast. (Whether ``State%z_eq`` *as a derived parameter* should count ULAs as matter is a separate, results.f90-level question for the background report.)
   * - ``OmegaK = 1 - OmegaT - OmegaAx - OmegaV - omegar - Omegan`` (recfast_axion.F90:563)
     - [OBSOLETE]
     - Modern uses ``Calc%OmegaK = State%CP%omk`` — the true input curvature, exact for any composition. No port needed.
   * - Signature/threading: Recombination_init extra args, ION extra args, ``recdverk`` clone (recfast_axion.F90:470-471, 707, 713, 819, 1434-2161)
     - [OBSOLETE]/[PLUMBING]
     - Entirely superseded: modern ION reaches the state via ``Recomb%State``; modern Recombination_init takes ``(this, State, WantTSpin)``. Intent (give ION access to rho_ax(a), a_osc) is satisfied by giving the recfast module access to the axion component through State. **Do not port recdverk** (and do not replicate its ``real(OmegaAx)`` single-precision conversion bug at recfast_axion.F90:2000-2049).
   * - RK ``tol=1.D-3`` (recfast_axion.F90:515) vs 1.D-5
     - **[ACCURACY]**
     - Hardcoded loosening, 100x, undocumented; likely a workaround for the dtauda/grhoax kink at a_osc. Recommend keeping modern 1.D-5 and verifying; only act if DVERK fails near a_osc. Flag in regression tests (xe(z) differences at the 1e-3-tol level would matter for high-ell TT/EE).
   * - Nz comments, debug writes, name string, re-indentation
     - [COSMETIC]
     - Ignore.
   * - timeH EdS estimate, all fudges, all rates, Saha switches
     - —
     - Unchanged in AxiECAMB; stock in 1.6.7. Nothing to do.

Classification of reionization changes for the port
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Change (AxiECAMB location)
     - Class
     - Port verdict
   * - ``helium_fullreion_redshiftstart`` 5 -> 7 (reionization.f90:33)
     - [ACCURACY]
     - Hardcoded tweak (modern default 5.5, user-settable as ``helium_redshiftstart``). Set via parameter default in the axion model setup only if needed to reproduce AxiECAMB; document either way.
   * - tau_start/tau_complete rombint split at a_osc (reionization.f90:149, 197-219)
     - **[PLUMBING with physics motivation — must re-derive]**
     - In 1.6.7 these moved to results.f90:1817-1820: ``State%reion_tau_start = max(0.05_dl, State%TimeOfZ(reion_z_start, 1d-3))`` and ``State%DeltaTime(...)``, where ``CAMBdata_DeltaTime = Integrate_Romberg(this, dtauda,a1,a2,atol)`` (results.f90:610-619). The *same Romberg-across-the-kink* problem will recur there — and in every other ``DeltaTime``/``TimeOfZ``/sound-horizon integral in results.f90, not just reionization. **Re-derive generically:** either (a) make ``CAMBdata_DeltaTime`` (or the axion component's dtauda) split integrals at a_osc when ``a1 < a_osc < a2``, or (b) smooth rho_ax(a) across the switch so dtauda is C1. Option (a) at the DeltaTime level fixes all call sites at once; the AxiECAMB patch fixed only reionization.
   * - Optical depth integrand uses dtauda (unchanged)
     - —
     - Modern ``GetReionizationOptDepth`` (results.f90:1243-1253) integrates ``reion_doptdepth_dz`` with the state background — automatically axion-aware. Nothing to port. (Same caveat: ``Integrate_Romberg`` over 0..zstart crosses a_osc for light axions; covered by the generic DeltaTime/Romberg fix.)
   * - Silenced non-convergence aborts in zreFromOptDepth (reionization.f90:320, 325-328)
     - [OBSOLETE]
     - Modern raises recoverable ``GlobalError(error_reionization)`` — superior. Do not port.
   * - ``DeltaTime_external`` declaration + failed attempts (reionization.f90:57, 198, 212)
     - [COSMETIC]
     - Dead code; ignore.

5. Surprises / risks for the port designer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **The recfast port is far smaller than the diff suggests.** Main H(z) was never hardcoded in Nov13 — it already came from ``dtauda``. The only genuine recombination physics to port is the axion term in the analytic ``dHdz`` of the Tmat smoothing term (plus deciding whether to instead replace dHdz with a numerical derivative, which the AxiECAMB author himself recommends in a comment at recfast_axion.F90:1089).
2. **Accuracy regression risk:** RK tolerance silently loosened 1e-5 -> 1e-3 in AxiECAMB. Do not carry over without evidence it is needed; if needed, it signals the dtauda kink problem, which should be fixed at the source.
3. **The a_osc kink is a cross-cutting issue.** AxiECAMB patched it in three uncoordinated places visible in this analysis alone: one-sided FD in ION, split rombint in Reionization_Init, and (probably) the loosened RK tol. In 1.6.7 the kink will hit ``CAMBdata_DeltaTime``, ``TimeOfZ``, ``sound_horizon``, ``GetReionizationOptDepth``, etc. A single generic treatment (split Romberg at a_osc, or smooth the EFA matching) is strongly preferred.
4. **Latent precision bug not to replicate:** ``recdverk`` passes ``real(OmegaAx)`` (single precision) in 7 of 8 fcn calls (recfast_axion.F90:2000-2049); benign only because OmegaAx is unused in ION.
5. **Error-handling was disabled** in ``Reionization_zreFromOptDepth`` — silent acceptance of unconverged zre(tau). Modern recoverable-error mechanism supersedes this; porting the silencing would mask real failures in axion parameter scans.
6. **helium_fullreion_redshiftstart=7** is an unexplained hardcoded deviation (modern default 5.5, Nov13 default 5.0); it changes xe(z) at 3.5<z<7 and hence tau(zre) mapping at the ~1e-3 level. Must be set consciously to reproduce AxiECAMB spectra.
7. Modern recfast also changed variables (y(3)=a*Tmat, z_scale Tcmb rescaling of the Gaussian fudge) — the ported axion dHdz term must go into the ``daTmat_dz`` form (recfast.f90:1034-1036), not the old f(3) form.


Original-code analysis: driver, build system and halofit (inidriver_axion.F90, halofit_ppf.f90)
------------------------------------------------------------------------------------------------

Sources analyzed:

- ``/Users/vivianmiranda/data/research/WayneHu/rayne/.port_analysis/diffs/inidriver_axion.diff`` (OLDCAMB/inidriver.F90 -> AxiECAMB/inidriver_axion.F90)
- ``/Users/vivianmiranda/data/research/WayneHu/rayne/.port_analysis/diffs/Makefile.diff``
- ``/Users/vivianmiranda/data/research/WayneHu/rayne/.port_analysis/diffs/Makefile_main.diff``
- ``/Users/vivianmiranda/data/research/WayneHu/rayne/.port_analysis/diffs/halofit_ppf.diff`` (OLDCAMB/halofit.f90 -> AxiECAMB/halofit_ppf.f90)
- ``/Users/vivianmiranda/data/research/WayneHu/rayne/.port_analysis/diffs/params.ini.diff``
- Originals: ``/Users/vivianmiranda/data/research/WayneHu/rayne/AxiECAMB/inidriver_axion.F90`` (758 lines), ``AxiECAMB/halofit_ppf.f90``, ``AxiECAMB/Makefile``, ``AxiECAMB/Makefile_main``, ``AxiECAMB/params.ini``, plus cross-references into ``AxiECAMB/axion_background.F90``, ``AxiECAMB/modules.f90``, ``AxiECAMB/constants.f90``, ``AxiECAMB/equations_ppf.f90``.

All file:line references below are to the AxiECAMB files unless noted.

1. inidriver_axion.F90 — new ini parameters and validation/transformation logic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1.1 Module usage and declarations (top of file)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

[PLUMBING] New ``use`` statements (inidriver_axion.F90:14-22):

- ``use Precision``
- ``use NonLinear ! RH made this change for axions`` (needed for ``NonLinear_ReadParams`` / ``halofit_version``)
- ``use axion_background`` (the KG background solver module; provides ``w_evolve``, ``get_phase_info``, and the lookup tables)

[PLUMBING] New driver-local declarations (inidriver_axion.F90:36-64):

.. code-block:: fortran

    real(dl) output_factor, nmassive,omnuh2,nu_massless_degeneracy,fractional_number
    real(dl) actual_massless,neff_i
    type (CAMBdata)  :: AxionIsoData ! Adding this for the iso stuff
    type (CAMBdata)  :: AxionAdiData ! Adding this for the iso stuff
    real(dl) hnot
    integer iter_dfac, iter_dfacETA
    real(dl) aeq,omegar,phiinit,omegah2_rad,rhocrit, nnu, rh_num_nu_massless
    real clock_totstart, clock_totstop ! RH timing
    integer reni ! RH
    integer badflag
    real(dl), allocatable :: RHCl_temp(:,:,:), RHCl_temp_lensed(:,:,:), RHCl_temp_tensor(:,:,:)
    real(dl) twobeta_tgt, beta_coeff, y_phase, movHETA_beta, movHETA_new
    real(dl) twobeta_new, twobeta_old1, twobeta_old2, hETA_beta, hETA_new
    real(dl) hosc_new, hosc_old1, hosc_old2, hETA_old1, hETA_old2, beta_tol !RL 030624

[COSMETIC/OBSOLETE] ``AxionIsoData``, ``AxionAdiData``, ``hnot``, ``aeq``, ``omegar``, ``phiinit``, ``nnu``, ``reni``, the clock variables, and ``omnuh2`` are declared but never used in the driver (leftovers). ``omegah2_rad``, ``rhocrit``, ``rh_num_nu_massless``, ``badflag``, the ``iter_*``/``twobeta_*``/``hosc_*``/``hETA_*`` variables and ``RHCl_temp*`` arrays ARE used.

1.2 do_nonlinear warning
^^^^^^^^^^^^^^^^^^^^^^^^^

[PLUMBING] inidriver_axion.F90:98-103:

.. code-block:: fortran

    P%NonLinear = Ini_Read_Int('do_nonlinear',NonLinear_none)
    if (P%NonLinear == 1 .or. P%NonLinear == 2) then
       write(*, *) 'Warning: do_nonlinear options are not well-tested in this version'
    end if

(A commented-out earlier version forcibly reset ``P%NonLinear = 0``.) Note: warning checks only ==1 and ==2, not ==3 (NonLinear_both). Port intent: warn users halofit was not validated for ULA models.

1.3 H0-derived quantities (new CAMBparams fields)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

[PHYSICS]/[PLUMBING] inidriver_axion.F90:140-142, immediately after reading hubble:

.. code-block:: fortran

    P%H0     = Ini_Read_Double('hubble')
    P%H0_in_Mpc_inv = dble(P%H0)/dble(c/1.0d3) !RL
    P%H0_eV = h_P*P%H0_in_Mpc_inv/(Mpc_in_sec*2._dl*const_pi*elecV) !RL

- ``P%H0_in_Mpc_inv`` [Mpc^-1]: H0 (km/s/Mpc) / (c in km/s).
- ``P%H0_eV`` [eV]: hbar*H0 in eV, i.e. ``h_P * H0_in_Mpc_inv / (Mpc_in_sec * 2*pi * elecV)``. Constants from ``constants.f90``: ``h_P = 6.62606896e-34`` (J s), ``elecV = 1.60217646e-19`` (J/eV), ``Mpc_in_sec = Mpc/c``. These are new ``CAMBparams`` fields used everywhere ``m/H0`` is needed.

1.4 Early read of omk and omegan (reordering)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

[PLUMBING] inidriver_axion.F90:144-150. ``P%omegak = Ini_Read_Double('omk')`` and ``P%omegan`` (from ``omnuh2``/(h^2) or ``omega_neutrino``) are read BEFORE the main density block, with comment:

.. code-block:: fortran

    !RL 120124: I need radiation fraction for assigning the DE fraction of axions, which in turn needs the neutrino fraction. Hence I have to take out this assignment separately

Note ``omegak`` is now a stored ``CAMBparams`` field set directly (OLDCAMB driver only used omk to derive omegav).

1.5 The m/H switch ratio dfac is internal (not an ini parameter)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

[PHYSICS]+[ACCURACY] inidriver_axion.F90:152-154:

.. code-block:: fortran

    !!P%dfac = Ini_Read_Double('movH_switch') !RL 092623 switch time
    P%dfac = 10._dl !RL 121924 making movH internal
    ntable = nint(P%dfac*100) + 1 !RL 111123

- ``P%dfac`` = the KG->effective-fluid switch threshold, m/H at switch. Hardcoded default 10 (formerly the ini key ``movH_switch``, now commented out — do NOT expose it).
- ``ntable`` is a module variable in modules.f90:249 (``integer :: ntable = 5000 !RL 111123. ntable should be properly set in inidriver_axion.``) giving the number of rows of the axion background lookup tables; rule: ``ntable = nint(P%dfac*100) + 1`` (100 points per unit of m/H, +1). It is recomputed every time ``dfac`` changes during the driver-level iterations (see §2).

1.6 Radiation density computed in the driver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

[PLUMBING] (physics-required input to the axion budget, but in CAMB 1.6.7 these quantities exist natively). inidriver_axion.F90:160-163:

.. code-block:: fortran

    rhocrit=(8.0d0*const_pi*G*1.d3/(3.0d0*((1.d7/(MPC_in_sec*c*1.d2))**(2.0d0))))**(-1.0d0)
    omegah2_rad=((P%tcmb**4.0d0)/(rhocrit))/(c**2.0d0) !RL replaced the COBE temperature
    omegah2_rad=omegah2_rad*a_rad*1.d1/(1.d4)

i.e. omega_gamma h^2 = a_rad T_cmb^4 / (rho_crit c^2) in CGS-ish unit gymnastics (rhocrit here is the critical density for H0=100 km/s/Mpc divided by h^2... the 1e1/1e4 factors are unit conversions). Then photon+massless-neutrino:
inidriver_axion.F90:240-241:

.. code-block:: fortran

    grhog= ((kappa/(c**2.0d0)*4.0d0*sigma_boltz)/(c**3.0d0))*(P%tcmb**4.0d0)*(Mpc**2.0d0) !RL replaced the COBE temperature
    P%grhor = (7.0d0/8.0d0)*((4.0d0/11.0d0)**(4.0d0/3.0d0))*grhog

(``P%grhor`` is a NEW CAMBparams field — per-neutrino-species 8*pi*G*rho*a^4 in Mpc^-2, same formula CAMBparams_Set uses internally.)
inidriver_axion.F90:269-270:

.. code-block:: fortran

    omegah2_rad= omegah2_rad+(rh_Num_Nu_massless*P%grhor*(c**2.0d0)/((1.d5**2.0d0)))/3.0d0
    P%omegah2_rad = omegah2_rad

``P%omegah2_rad`` [dimensionless, omega_r h^2 including photons + massless neutrinos] is a new CAMBparams field used by the axion background solver and the DE-budget formulas below.

1.7 Neutrino bookkeeping moved/duplicated in the driver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

[PLUMBING] (re-derive: modern CAMB already does all of this in ``CAMBparams`` / ``ThermalNuBackground``). Two blocks:

**(a)** inidriver_axion.F90:199-238 — the standard ``CAMBparams_Set`` massive-neutrino fixups copied into the driver so that N_eff entering ``omegah2_rad`` is correct ("4/8 DG Error in original AxionCAMB: massless neutrino contribution wrong when neutrinos are massive"):

.. code-block:: fortran

    if (P%Omegan == 0 .and. P%Num_Nu_Massive /=0) then
       print*, 'omeganuh2 is set to 0 but we still have massive neutrinos, resetting'
       if (P%share_delta_neff) then
          P%Num_Nu_Massless = P%Num_Nu_Massless + P%Num_Nu_Massive
       else
          P%Num_Nu_Massless = P%Num_Nu_Massless + sum(P%Nu_mass_degeneracies(1:P%Nu_mass_eigenstates))
       end if
       P%Num_Nu_Massive  = 0
       P%Nu_mass_numbers = 0
    end if
    P%Nu_massless_degeneracy = P%Num_Nu_massless !N_eff for massless neutrinos !RL 061425

then (215-238) for ``Num_nu_massive > 0`` with ``share_delta_neff``:

.. code-block:: fortran

    fractional_number = P%Num_Nu_massless + P%Num_Nu_massive
    actual_massless = int(P%Num_Nu_massless + 1e-6_dl)
    neff_i = fractional_number/(actual_massless + P%Num_Nu_massive)
    nu_massless_degeneracy = neff_i*actual_massless
    P%Nu_massless_degeneracy=nu_massless_degeneracy
    P%Nu_mass_degeneracies(1:P%Nu_mass_eigenstates) = P%Nu_mass_numbers(1:P%Nu_mass_eigenstates)*neff_i

``P%Nu_massless_degeneracy`` is a NEW CAMBparams field (in OLDCAMB this was a local inside CAMBparams_Set).

**(b)** inidriver_axion.F90:244-265 — ``rh_num_nu_massless`` (driver local) = effective number of massless species used for the radiation density:

- if ``omegan==0 .and. Num_Nu_Massive/=0``: ``rh_num_Nu_Massless = Num_Nu_Massless + Num_Nu_Massive`` (share_delta_neff) or ``+ sum(Nu_mass_degeneracies)``;
- if ``omegan==0 .and. Num_Nu_Massive==0``: ``rh_num_nu_massless = P%Num_Nu_Massless``;
- if ``omegan>0 .and. Num_nu_massive>0``: ``rh_num_nu_massless = P%Num_Nu_Massless ! only using the massless neutrinos, and adding in omnuh2 later``.

Massive-neutrino energy is therefore NOT in ``P%omegah2_rad``; it enters the axion background solver separately via ``omegan``. In modern CAMB use ``CP%grhog``, ``CP%grhornomass``, etc. directly.

1.8 NEW ini parameters (complete list read by the driver)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All inside the ``use_physical = T`` branch (inidriver_axion.F90:272-325) unless noted:

.. list-table::
   :header-rows: 1
   :widths: auto

   * - ini key
     - type
     - default
     - sets
     - line
   * - ``use_axfrac``
     - logical
     - ``.false.``
     - ``P%use_axfrac``
     - 275
   * - ``m_ax``
     - double
     - none (required)
     - ``P%ma`` [eV]
     - 277 (also 333 in use_physical=F branch)
   * - ``omdah2``
     - double
     - none (req. if use_axfrac)
     - ``P%omegada = omdah2/(H0/100)^2``
     - 283
   * - ``axfrac``
     - double
     - none (req. if use_axfrac)
     - ``P%axfrac``
     - 285
   * - ``omch2``
     - double
     - none (req. if .not.use_axfrac)
     - ``P%omegac = omch2/(H0/100)^2``
     - 300
   * - ``omaxh2``
     - double
     - none (req. if .not.use_axfrac)
     - ``P%omegaax = omaxh2/(H0/100)^2``
     - 301
   * - ``Hinf``
     - double
     - none (required)
     - ``P%Hinf`` (read as log10 of H_inflation in GeV)
     - 316
   * - ``axion_isocurvature``
     - logical
     - ``.true.``
     - ``P%axion_isocurvature`` (then forced to F)
     - 319
   * - ``omega_axion``
     - double
     - none
     - ``P%omegaax = omega_axion/(H0/100)^2`` (use_physical=F branch; NOTE: divides by h^2 even though the key name suggests an Omega — looks like a bug/inconsistency in the non-physical branch)
     - 332
   * - ``halofit_version``
     - int
     - ``halofit_default`` (=1)
     - module var ``halofit_version`` via ``NonLinear_ReadParams``
     - driver 354 -> halofit_ppf.f90:49

NOT read by the driver despite appearing in params.ini: ``alpha_ax`` (dead ini key; ``P%alpha_ax`` is COMPUTED in axion_background.F90:1162, see §1.10) and ``tens_ratio`` (dead key, grep finds no reader anywhere).

Removed relative to OLDCAMB params.ini (still readable by ``DarkEnergy_ReadParams`` in equations_ppf, but commented from the ini): ``wa``, ``use_tabulated_w``, ``wafile``. ``w`` and ``cs2_lam`` remain.

1.9 Transformation/validation logic for the axion parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

[PHYSICS] m_ax log convention, inidriver_axion.F90:277-279:

.. code-block:: fortran

    P%ma     = Ini_Read_Double('m_ax') !! RH axion mass
    if (P%ma < 0) P%ma = 10**P%ma ! RH making this exponential from the inidriver
    P%m_ovH0 = P%ma/P%H0_eV !RL 050724

- ``m_ax`` in eV if positive; if negative it is interpreted as log10(m_ax/eV).
- ``P%m_ovH0`` (new CAMBparams field, dimensionless) = m/H0 with H0 in eV; this is the discriminator between "dark-matter-like" (m/H0 >= 10) and "dark-energy-like" (m/H0 < 10) axions used in the driver, halofit and modules.

[PHYSICS] Density assignment, ``use_axfrac = T`` (inidriver_axion.F90:281-296):

.. code-block:: fortran

    P%omegada = Ini_Read_Double('omdah2')/(P%H0/100)**2
    P%axfrac = Ini_Read_Double('axfrac')
    if (P%m_ovH0 .ge. 10._dl) then !RL 120124 - DM case
       P%omegaax = P%axfrac*P%omegada
       P%omegac = (1-P%axfrac)*P%omegada
    else !RL 120124 - DE case
       write(*, *) 'Note: m/H0 < 10, axfrac is ULA fraction in dark energy'
       P%omegac = P%omegada
       P%omegaax = P%axfrac*(1._dl-P%omegab-P%omegac - P%omegan -P%omegak - P%omegah2_rad/((P%H0/1.d2)**2.0d0))
    end if

i.e. DM case: axfrac splits the total dark matter ``omdah2`` between axion and CDM. DE case (m/H0 < 10): ALL of omdah2 goes to CDM (``omch2 = omdah2``), and axfrac is reinterpreted as the axion fraction of the dark-energy budget Omega_DE = 1 - Omega_b - Omega_c - Omega_nu - Omega_k - Omega_r.

[PHYSICS] ``use_axfrac = F`` (inidriver_axion.F90:298-308):

.. code-block:: fortran

    P%omegac = Ini_Read_Double('omch2')/(P%H0/100)**2
    P%omegaax = Ini_Read_Double('omaxh2')/(P%H0/100)**2
    if (P%m_ovH0 .ge. 10._dl) then !RL 120124 - DM case
       P%axfrac = P%omegaax/(P%omegac+P%omegaax)
    else !RL 120124 - DE case
       P%axfrac = P%omegaax/(1.0d0-P%omegab-P%omegac - P%omegan -P%omegak - P%omegah2_rad/((P%H0/1.d2)**2.0d0))
    end if

``P%axfrac`` is a stored CAMBparams field (used later, e.g. in the ETA-phase iteration condition and in isocurvature normalization).

[PHYSICS] Lambda closure including radiation (inidriver_axion.F90:311-312):

.. code-block:: fortran

    !Compute value of cosmological constant including curvature and radiation (photons + massless neutrinos) self consistently
    P%omegav = 1._dl-P%omegab-P%omegac - P%omegan -P%omegak-P%omegaax - P%omegah2_rad/((P%H0/1.d2)**2.0d0)

Difference from OLDCAMB: OLDCAMB driver had ``P%omegav = 1- omk - omegab - omegac - omegan`` (no radiation, no axion). [OBSOLETE in 1.6.7]: modern CAMB computes omega_de internally from ``omk`` and the physical densities including radiation, so only "subtract ``omegaax`` from the budget" must be ported.

[PHYSICS] Hinf transformation (inidriver_axion.F90:316-317):

.. code-block:: fortran

    P%Hinf = Ini_Read_Double('Hinf') ! H inflation in GeV
    P%Hinf = (10**P%Hinf)/mplanck ! computing the ratio of Hinflation to Mplanck

``mplanck = 2.435e18`` GeV (reduced Planck mass, constants.f90:44). So stored ``P%Hinf`` is dimensionless H_inf/M_pl; the ini value is log10(H_inf/GeV) (params.ini default 13.7).

[PLUMBING] Isocurvature disabled in this release (inidriver_axion.F90:319-324):

.. code-block:: fortran

    P%axion_isocurvature = Ini_Read_Logical('axion_isocurvature', .true.)
    !RL 121924 disable isocurvature
    if (P%axion_isocurvature .eqv. .true.) then
       write(*, *) 'WARNING: axion isocurvature disabled in this release, proceeding without'
       P%axion_isocurvature = .false.
    end if

[PLUMBING] use_physical=F branch (inidriver_axion.F90:327-335) is mostly broken/legacy: reads ``omega_baryon``, ``omega_cdm``, ``omega_lambda``, ``omega_axion`` (divided by h^2 — inconsistent), ``m_ax``; never sets ``m_ovH0``, ``axfrac``, ``Hinf``, ``axion_isocurvature``. Recommend supporting only physical densities in the port.

1.10 Isocurvature amplitude quantities (computed, not read)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

[PHYSICS] Computed inside ``w_evolve`` (axion_background.F90:1151, 1158-1165) once the initial field value is known:

.. code-block:: fortran

    Params%phiinit=vtwiddle_init*sqrt(6.0d0)        ! initial field in reduced-Planck units
    if (Params%axion_isocurvature) then
       Params%amp_i = Params%Hinf**2/(pi**2*Params%phiinit**2)
       Params%r_val  = 2*(Params%Hinf**2/(pi**2.*Params%InitPower%ScalarPowerAmp(1)))
       Params%alpha_ax = Params%amp_i/Params%InitPower%ScalarPowerAmp(1)
    end if

- ``amp_i`` = isocurvature power amplitude (Hinf/(pi*phi_i))^2 with both in reduced Planck units.
- ``r_val`` = tensor-to-scalar ratio 2 Hinf^2/(pi^2 As).
- ``alpha_ax`` = amp_i/As, the isocurvature fraction; it is a CAMBparams field (modules.f90:124) consumed elsewhere (derived-parameter output in modules; NOT read from ini even though ``alpha_ax = 0`` appears in params.ini).

These feed the second ``CAMB_GetResults`` call (see §2.4).

1.11 Other small driver changes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- [PLUMBING] inidriver_axion.F90:354: ``if (P%NonLinear/= NonLinear_none) call NonLinear_ReadParams(DefIni) ! RH axions`` — reads ``halofit_version``.
- [PLUMBING] inidriver_axion.F90:417-420: initial-vector read extended to 6 modes:

  .. code-block:: fortran

      read (numstr,*) P%InitialConditionVector(1:initial_iso_axion)

  with ``initial_iso_axion=6`` and ``initial_nummodes = initial_iso_axion`` defined in equations_ppf.f90:374-375. (OLDCAMB read ``1:initial_iso_neutrino_vel`` = 5.)
- [PLUMBING] ``derived_parameters`` ini read is unchanged (``P%DerivedParameters = Ini_Read_Logical('derived_parameters',.true.)``); derived-parameter output itself lives in camb.f90/modules.f90 (other reports). The key was deleted from params.ini, so the default ``.true.`` applies.
- [COSMETIC] Whole file re-indented (2-space); many ``! RH`` / ``!RL`` chatter comments; commented-out cpu_time timing calls.

2. Driver-level iteration / shooting (between ``Ini_Close`` and ``CAMB_GetResults``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the biggest structural addition: the axion background solver is run (and re-run) at the DRIVER level, before ``CAMB_ValidateParams``/``CAMB_GetResults``. All of this must find a home in modern CAMB (natural place: inside ``CAMBparams_SetParams``/``SetBackgroundOutputs`` or the dark-energy/axion component's ``Init``), since the Python wrapper never goes through inidriver.

2.1 Initial background solve
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

[PHYSICS]/[PLUMBING] inidriver_axion.F90:519-523:

.. code-block:: fortran

    call init_massive_nu(P%omegan /=0) !RL added 07/10/23
    P%a_skip = 1._dl/(800._dl + 1._dl) !RL for skipping
    P%a_skipst = 1._dl/(1300._dl + 1._dl) !lower threshold of recombination skip
    P%dfac_skip = 0._dl !RL initializing P%dfac_skip just to make sure it doesn't get assigned to some random numbers
    call w_evolve(P, badflag)

- ``init_massive_nu`` must be called BEFORE the axion solver because ``w_evolve``/``auxiIC`` need the massive-neutrino background (``Nu_background``). In OLDCAMB this was only called inside ``CAMBparams_Set``.
- ``P%a_skip = 1/801`` (z=800) and ``P%a_skipst = 1/1301`` (z=1300): the "recombination skip window" — if the KG->fluid switch (a_osc) lands inside [a_skipst, a_skip), the switch is pushed past the window (see §2.3). New CAMBparams fields: ``a_skip``, ``a_skipst``, ``dfac_skip`` (the dfac value that ``w_evolve`` reports would move a_osc past a_skip).
- ``w_evolve(P, badflag)`` (axion_background.F90:35) solves the KG background, shoots for ``phiinit`` to match ``omegaax``, fills module tables ``loga_table, phinorm_table, phidotnorm_table, (+_ddlga splines), rhoaxh2ovrhom_logtable(_buff)``, and sets outputs on P: ``a_osc``, ``ah_osc``, ``ahosc_ETA``, ``aeq``, ``phiinit``, ``omegar``, and the isocurvature amplitudes (§1.10). It also sets the module variable ``aeq_LCDM`` (modules.f90:250, "added for photon oscillation skipping").

2.2 ETA-phase shooting/bisection on dfac (m/H at switch)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

[PHYSICS] (accuracy-motivated physics: aligns the switch with a fixed photon-driving phase to suppress the spurious feature from switching mid-oscillation before recombination). Trigger condition, inidriver_axion.F90:525-528:

.. code-block:: fortran

    !First check if we're in the window before recombination where the photon ETA is an issue
    if (P%dfac .lt. 23._dl .and. P%m_ovH0 .ge. 10._dl .and. P%ma .lt. 1.e-25_dl &
         &.and. P%a_osc*(P%omegaax/(P%omegac+P%omegaax))/aeq_LCDM .gt. 0.03_dl &
         &.and. P%a_osc .lt. P%a_skipst) then !RL 022624 - is an empirically tuned number 0.03_dl

So this runs only when: dfac<23, DM-like axion (m/H0>=10), m < 1e-25 eV, the axion is dynamically important before equality (``a_osc * f_ax / aeq_LCDM > 0.03``, empirically tuned), and a_osc occurs before z=1300.

Target phase and first guess, inidriver_axion.F90:530-535:

.. code-block:: fortran

    twobeta_tgt = 7.08_dl*const_pi!22.25_dl!10.5_dl*const_pi!
    !P%dfac = twobeta_tgt + 0.75_dl*const_pi !First guess using radiation domination
    P%dfac = twobeta_tgt + 0.75_dl*const_pi - twobeta_tgt**2/&
         &(4._dl*(twobeta_tgt + 2._dl*P%m_ovH0*(aeq_LCDM**1.5_dl)/&
         &sqrt(2._dl*(P%omegac + P%omegab + P%omegan + P%omegaax)))) !First guess considering matter-radiation equality in LCDM
    ntable = nint(P%dfac*100) + 1
    call w_evolve(P, badflag)
    call get_phase_info(P, y_phase, beta_coeff, movHETA_new, twobeta_new)

``get_phase_info`` (axion_background.F90:1207-1218):

.. code-block:: fortran

    y_beta = Params%a_osc/aeq_LCDM
    movHETA = Params%dfac*Params%ah_osc/Params%ahosc_ETA
    beta_coeff = (4._dl*(y_beta**2 - y_beta - 2.0_dl + 2.0_dl*sqrt(1.0_dl + y_beta)))/(3._dl*(y_beta**2))
    beta2x = movHETA*beta_coeff - const_pi*3._dl*(1.0_dl + y_beta)/(4.0_dl + 3.0_dl*y_beta)

(``beta2x`` = "2*beta" phase of the field oscillation at the switch in the ETA (effective-time-averaged) Hubble; ``ah_osc`` = (aH) at a_osc instantaneous, ``ahosc_ETA`` = (aH) at a_osc with time-averaged H — both CAMBparams outputs of w_evolve.)

Bracket guess, inidriver_axion.F90:539-542:

.. code-block:: fortran

    movHETA_beta = (twobeta_tgt + const_pi*3._dl*(1.0_dl + y_phase)/(4.0_dl + 3.0_dl*y_phase))/beta_coeff
    hETA_beta = (P%dfac/movHETA_beta) * (P%ah_osc/P%a_osc)
    beta_tol = 2.e-2_dl*const_pi

If ``abs(twobeta_new - twobeta_tgt) .gt. beta_tol``, a two-stage root find on the variable ``hosc = ah_osc/a_osc`` (i.e. H at switch) is performed, where each new trial updates ``P%dfac = P%dfac * ((P%ah_osc/P%a_osc)/hosc_new)``, ``ntable = nint(P%dfac*100)+1``, then ``call w_evolve(P, badflag)`` and ``call get_phase_info(...)`` again:

1. Bracketing loop (inidriver_axion.F90:546-578, ``iter_dfacETA`` up to 500): step ``hosc_new = 1._dl*(hETA_beta-hETA_old1)+hosc_old2`` (first step uses factor 2: ``hosc_new = 2._dl*(hETA_beta-hETA_old1)+hosc_old1``), until ``(twobeta_old2-twobeta_tgt)*(twobeta_old1-twobeta_tgt) < 0``.
2. Bisection loop (inidriver_axion.F90:580-616, up to 500 iters): ``hosc_new = (hosc_old1 + hosc_old2) / 2._dl``, replacing whichever bracket end has the same sign of ``(twobeta - twobeta_tgt)``, until ``abs(twobeta_new - twobeta_tgt) < beta_tol``.

[ACCURACY] hardcoded numbers worth keeping verbatim: ``twobeta_tgt = 7.08*pi``, ``beta_tol = 0.02*pi``, trigger numbers ``23``, ``1e-25 eV``, ``0.03``, max 500 iterations.

2.3 Recombination-window skip loop
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

[PHYSICS] inidriver_axion.F90:622-636:

.. code-block:: fortran

    do iter_dfac = 1, 500
       if (P%a_osc .lt. P%a_skip*(1._dl - 1.e-2_dl) .and. P%a_osc .ge. P%a_skipst) then
          !RL 032024: 1e-2 is the tolerence to eliminate additional loops if we don't skip to exactly after a_skip due to numerical factors
          P%dfac = P%dfac_skip
          ntable = nint(P%dfac*100) + 1
          call w_evolve(P, badflag)
       else
          exit
       end if
    end do
    if (iter_dfac .gt. 500 .and. P%a_osc .lt. P%a_skip) then
       print*, 'Warning: maximum iteration reached, but aosc still not skipped sufficiently: ...'
    end if

If the switch epoch a_osc falls inside the z∈(800,1300) window, re-run the background with ``dfac = P%dfac_skip`` (a larger switch threshold computed inside ``w_evolve`` that pushes a_osc past ``a_skip``), iterating until it leaves the window. Tolerance factor ``(1 - 1e-2)``.

Only AFTER all this: ``if (.not. CAMB_ValidateParams(P)) stop ...`` then ``CAMB_GetResults(P)`` (inidriver_axion.F90:639, 652).

NOTE on "shooting for H0/densities": the relic-abundance shooting (adjusting the initial field value ``phiinit`` to hit ``omegaax``) happens INSIDE ``w_evolve`` (axion_background.F90), not in the driver. The driver-level iterations above only tune ``dfac`` (the switch time). There is no H0 shooting at the driver level.

2.4 Axion isocurvature two-pass Cl computation (currently disabled)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

[PHYSICS] (port if isocurvature is to be revived; gated by ``P%axion_isocurvature`` which is forced false — see §1.9). inidriver_axion.F90:650-693:

1. First ``CAMB_GetResults(P)`` = adiabatic run; copy ``Cl_scalar(lmin:P%Max_l,1,C_Temp:C_last)``, ``Cl_lensed(...,C_Temp:C_Cross)``, ``Cl_tensor(...)`` into ``RHCl_temp``, ``RHCl_temp_lensed``, ``RHCl_temp_tensor``.
2. Reconfigure for the isocurvature pass:

   .. code-block:: fortran

       P%Scalar_initial_condition = 6
       P%InitPower%rat(1) =  0
       P%InitPower%ant(1) = 0
       P%InitPower%ScalarPowerAmp(1) = P%amp_i
       P%InitPower%an(1)= 1-P%r_val/8.d0
       call CAMB_GetResults(P)

   (i.e. run mode-6 axion isocurvature with amplitude ``amp_i = (Hinf/(pi*phiinit))^2``, tilt ``n_iso = 1 - r/8`` with ``r = 2 Hinf^2/(pi^2 As)``, no tensors in the second pass.)
3. Add the two spectra: ``Cl_scalar = Cl_scalar + RHCl_temp`` (same for lensed over ``lmin:lmax_lensed`` and tensor over ``lmin:P%Max_l_tensor``).

[RISK] The totally-correlated ``initial_vector`` mechanism is NOT used; isocurvature is added in power (uncorrelated). Tensor ``r`` is forced to 0 in the iso pass via ``rat(1)=0``.

2.5 Cleanup
^^^^^^^^^^^^

[PLUMBING] inidriver_axion.F90:725-735: explicit deallocation of the module-level axion tables before ``CAMB_cleanup``:
``loga_table, phinorm_table, phidotnorm_table, phinorm_table_ddlga, phidotnorm_table_ddlga, rhoaxh2ovrhom_logtable, rhoaxh2ovrhom_logtable_buff`` plus the three ``RHCl_temp*`` arrays. In a class-based port these become components of the results/params object (no global state).

3. Build system (Makefile, Makefile_main)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

3.1 Makefile
^^^^^^^^^^^^^

All [PLUMBING]/[OBSOLETE] — modern CAMB 1.6.7 has its own Makefile system; nothing here needs literal porting:

- Default compiler switched from ``ifort`` (``FFLAGS = -openmp -fast -W0 -WB -fpp2 -vec_report0``, ``F90CRLINK = -cxxlib``) to:

  .. code-block:: make

      F90C     = gfortran
      FFLAGS = -O3 -fopenmp -ffpe-summary=none

  [ACCURACY-adjacent note]: ``-ffpe-summary=none`` only suppresses FPE summaries at exit (cosmetic); no precision-affecting flag changes. A commented ifort line suggests ``-O3 -qopenmp`` was used in testing.
- FITS/HEALPIX (``FITSDIR``, ``FITSLIB``, ``HEALPIXDIR``) and the ``FISHER``/``Matrix_utils.o``/``-mkl`` blocks are commented out -> ``EXTCAMBFILES`` is undefined (empty). camb_fits / Fisher bispectrum support effectively dropped.

3.2 Makefile_main
^^^^^^^^^^^^^^^^^^

Module selection (AxiECAMB/Makefile_main:1-15):

.. code-block:: make

    EQUATIONS     ?= equations_ppf        # was: equations
    POWERSPECTRUM ?= power_tilt           # unchanged
    REIONIZATION ?= reionization          # unchanged
    RECOMBINATION ?= recfast_axion        # was: recfast
    NONLINEAR     ?= halofit_ppf          # was: halofit
    DRIVER        ?= inidriver_axion.F90  # was: inidriver.F90

[PLUMBING] So the binary is built from: ``equations_ppf.f90`` (the heavily modified PPF+axion equations), ``recfast_axion.F90``, ``halofit_ppf.f90``, ``inidriver_axion.F90``. The shipped ``halofit.f90``, ``cmbmainOMP.f90``, ``cosmorec.F90``, ``hyrec.F90``, ``lenspen.f90``, ``sigma8.f90``, ``Matrix_utils.F90`` are NOT compiled into the default target.

- **halofit.f90 (shipped) is byte-identical to OLDCAMB/halofit.f90** (verified with ``diff -q``) — it is dead code; only ``halofit_ppf.f90`` is used.
- **cmbmainOMP.f90 is NOT in the build** — ``CAMBOBJ`` references ``cmbmain.o`` only and the dependency line ``cmbmain.o: lensing.o $(NONLINEAR).o $(EQUATIONS).o $(BISPECTRUM).o`` is unchanged. (There is a large diff cmbmainOMP_vs_cmbmain.diff; whatever it contains is an unused experimental variant from the build's perspective.)

Object list and order (AxiECAMB/Makefile_main:34-35):

.. code-block:: make

    CAMBOBJ = constants.o utils.o $(EXTCAMBFILES) subroutines.o inifile.o $(POWERSPECTRUM).o $(RECOMBINATION).o \
        $(REIONIZATION).o modules.o bessels.o $(EQUATIONS).o $(NONLINEAR).o lensing.o $(BISPECTRUM).o cmbmain.o camb.o axion_background.o

[PLUMBING] ``axion_background.o`` appended (link order is irrelevant; compile order is driven by the explicit rule). New dependency rule (Makefile_main:75-77):

.. code-block:: make

    #RL added dependency fix of axion_background on modules
    axion_background.o: axion_background.F90 modules.o
                $(F90C) $(F90FLAGS) -c axion_background.F90

So the module DAG is: constants/utils/subroutines/inifile -> power_tilt/recfast_axion/reionization -> modules -> {bessels, axion_background, equations_ppf} -> halofit_ppf -> lensing -> SeparableBispectrum -> cmbmain -> camb -> driver. Note ``axion_background`` depends only on ``modules`` (uses ModelParams, constants, Precision, MassiveNu), while the driver and equations_ppf use ``axion_background``. (The driver is compiled together with the link step: ``camb: $(CAMBOBJ) $(DRIVER)``.) For CAMB 1.6.7: add ``axion_background`` (or its class equivalent) to ``fortran/Makefile``'s object list between ``results``-level modules and the equations module, mirroring this dependency.

- ``F90CRLINK ?= -lstdc++`` commented out; large commented MKL link-flag graveyard added near ``F90FLAGS`` [COSMETIC].
- ``clean:`` now also removes the ``camb`` binary: ``-rm -f *.o *.a *.d core *.mod camb`` [COSMETIC].
- ``CAMBLIB = libcamb_$(RECOMBINATION).a`` unchanged (so the lib would be ``libcamb_recfast_axion.a``).

4. halofit_ppf.f90 vs OLDCAMB halofit.f90
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Baseline note: OLDCAMB Nov13 ``halofit.f90`` is the Takahashi-only variant (with Bird-style fnu terms ``(2.080-12.39*(omm0-0.3))/(1+1.201e-03*y**3)`` and ``26.29*rk**2``, constant ``w_lam``). AxiECAMB ``halofit_ppf.f90`` is NOT a hand-edit of that file: it is the later upstream ``halofit_ppf.f90`` (AL Sept 14 version with ``halofit_version`` selector, ``wa_ppf`` support, JD variable-w ``omega_m/omega_v``, Bird 2012 fnu coefficients ``0.977`` / ``47.48`` / ``-6.4868+1.4373*rn**2``) with three axion-specific modifications on top. CAMB 1.6.7's ``fortran/halofit.f90`` already contains everything in the upstream base, so only the axion deltas below need porting.

4.1 Axion-specific changes (the only [PHYSICS] deltas)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**(a)** Version default forced to original Smith et al. (halofit_ppf.f90:38-40):

.. code-block:: fortran

    integer, parameter :: halofit_default = halofit_original ! DM15: Takahashi is not stable for axion models. Other versions agree well and are sensible to percent level for lensing \ell<4000.
    integer :: halofit_version = halofit_default

[PHYSICS]/[ACCURACY] In stock CAMB the default is Takahashi (and 1.6.7 defaults to ``mead2020``). The AxiECAMB authors found Takahashi unstable for axion models; the port should default the ULA model to ``halofit_original`` (or at minimum document/validate the choice).

**(b)** Axion included in/excluded from the halofit matter budget by mass (halofit_ppf.f90:64-85):

.. code-block:: fortran

    !! DM16: modification to treat axions in non-linear lensing S4 fiducial model.
    ! Axion mass kluge: include in computation on non-linear ratio only
    ! for masses that are non-linear at z>2.
    ! Otherwise they are treated as quintessence, i.e. ignored in this.
    ! Boundary found by hand for fiducial "test field" approximation,
    ! valid for low axion density.
    !Make sure to set the same boundary mass, m/H0 = 10, in subroutine outtransf.
    if (CP%m_ovH0.ge.10._dl) then !RL 120224
       omm0 = CP%omegac+CP%omegab+CP%omegan+CP%omegaax
       omv0_axion = 0._dl !RL 090925 DM cases, not adding DE
    else
       omm0 = CP%omegac+CP%omegab+CP%omegan
       omv0_axion = CP%omegaax !RL 090925 DE cases
    end if
    fnu = CP%omegan/omm0

(new module variable ``real(dl):: omv0_axion !RL 090925`` at line 37; original code was just ``omm0 = CP%omegac+CP%omegab+CP%omegan``).

- DM-like axion (m/H0 >= 10): axion counts as matter in ``omm0`` (so it enters om_m(z), om_v(z) and the fnu normalization).
- DE-like axion (m/H0 < 10): axion excluded from matter and added to the dark-energy density:

  .. code-block:: fortran

      om_m = omega_m(a, omm0, CP%omegav + omv0_axion, w_lam,wa_ppf) !RL 090925 added DE-like axion
      om_v = omega_v(a, omm0, CP%omegav + omv0_axion, w_lam,wa_ppf) !RL 090925 added DE-like axion

  (halofit_ppf.f90:93-94; stock code passes ``CP%omegav`` only). The DE-like axion is treated with the SAME (w_lam, wa_ppf) as the PPF dark energy inside omega_m/omega_v — an approximation.
- The original mass boundary was ``CP%ma.ge.1.e-25`` (left commented at line 73: ``!if (CP%ma.ge.1.e-25) then``); RL changed it to ``m_ovH0 >= 10`` (120224) to match the driver's DM/DE split.
- [RISK/cross-file consistency] The comment requires the same boundary in ``subroutine outtransf`` (modules.f90), which builds the transfer-function/power-spectrum outputs: whether delta_ax is included in delta_tot (the spectrum halofit rescales and that sigma8 is computed from) is decided in modules.f90, not here. The nonlinear INPUT spectrum is whatever ``MatterPowerData`` (via ``transfer_power_var``, default ``transfer_tot``) contains — see the modules.f90 report; in the DM-like case the axion is included there, consistent with ``omm0`` here.

**(c)** Error trap for runaway nonlinearity (halofit_ppf.f90:115-120), inside the rknl bisection:

.. code-block:: fortran

    else if (xlogr1>3.4999) then
       ! Totally crazy non-linear
       global_error_flag=349
          write(*,*) 'Error in halofit'
          goto 101
       end if

[PLUMBING] This also exists in modern CAMB (which calls ``GlobalError('Error in halofit')``); no port needed beyond keeping the behavior.

**(d)** ``acur`` saved for Takahashi wa terms (halofit_ppf.f90:35 ``real(dl):: ... acur``, set at line 95 ``acur = a``, used in ``0.1749*om_v*(1.+w_lam+wa_ppf*(1-acur))`` and ``0.2279*om_v*(1.+w_lam+wa_ppf*(1-acur))`` at lines 195/198). This is stock upstream halofit_ppf behavior, not axion-specific. [OBSOLETE] (already in 1.6.7).

4.2 NonLinear_ReadParams / halofit_version ini plumbing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

[PLUMBING] halofit_ppf.f90:45-51:

.. code-block:: fortran

    subroutine NonLinear_ReadParams(Ini)
    use IniFile
    Type(TIniFile) :: Ini
    halofit_version = Ini_Read_Int_File(Ini, 'halofit_version', halofit_default)
    end subroutine NonLinear_ReadParams

Called from the driver only when ``P%NonLinear /= NonLinear_none`` (inidriver_axion.F90:354). In 1.6.7, ``halofit_version`` is already a settable property of ``THalofit``; only the default (=1/original for axion runs) matters.

4.3 sigma8
^^^^^^^^^^^

NO changes to sigma8 computation in halofit_ppf.f90 (it contains none). ``Transfer_output_sig8(MT)`` is called unchanged by the driver (inidriver_axion.F90:707); any sigma8 changes live in modules.f90 (``Transfer_Get_sigma8``/``outtransf``) — covered by the modules report. The ``wint`` integration in halofit_ppf is byte-equivalent to stock (nint=3000 etc.).

4.4 Everything else in the halofit_ppf diff
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The remaining diff bulk (halofit subroutine restructured with version dispatch, Bird/Peacock branches, ``omega_m/omega_v`` gaining ``waval`` argument with ``Qa2 = aa**(-1.0-3.0*(wval+waval))*dexp(-3.0*(1-aa)*waval)``, NonLinear_GetRatios_all stop message, reformatting) = upstream halofit_ppf inheritance, all [OBSOLETE] w.r.t. CAMB 1.6.7 which already has equal-or-newer versions of these.

5. params.ini changes (complete list)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

New axion/user-facing parameters (params.ini):

.. list-table::
   :header-rows: 1
   :widths: auto

   * - key
     - default
     - meaning
     - actually read?
   * - ``m_ax``
     - ``1.e-27``
     - axion mass in eV; negative value = log10(m/eV)
     - yes (driver:277)
   * - ``use_axfrac``
     - ``T``
     - if T use omdah2+axfrac; if F use omch2+omaxh2
     - yes (driver:275)
   * - ``omaxh2``
     - ``0.1200``
     - omega_ax h^2 (used if use_axfrac=F)
     - yes (driver:301)
   * - ``omdah2``
     - ``0.1200``
     - total dark "matter" h^2 (used if use_axfrac=T)
     - yes (driver:283)
   * - ``axfrac``
     - ``1.0000``
     - axion fraction of DM (m/H0>=10) or of DE (m/H0<10; then omch2=omdah2) — comment in ini states exactly this
     - yes (driver:285)
   * - ``axion_isocurvature``
     - ``F``
     - enable iso mode (force-disabled in code)
     - yes (driver:319)
   * - ``alpha_ax``
     - ``0``
     - DEAD KEY — never read; alpha_ax computed internally
     - no
   * - ``Hinf``
     - ``13.7``
     - log10(H_inflation/GeV)
     - yes (driver:316)
   * - ``halofit_version``
     - ``1``
     - 1=original (default), 2=bird, 3=peacock, 4=takahashi
     - yes (halofit_ppf:49 via driver:354)
   * - ``omega_axion``
     - (commented ``0``)
     - use_physical=F branch only
     - yes (driver:332)
   * - ``tens_ratio``
     - ``0``
     - DEAD KEY — no reader anywhere ("DM: gave this twice because of problem using ini_driver to read initial_ratio")
     - no

Changed defaults (none of these are code requirements — they retune the sample file to a Planck-2018-like ULA fiducial) [ACCURACY/COSMETIC]:

- ``get_transfer F -> T``; ``l_max_scalar 2200 -> 2700``; ``k_eta_max_scalar`` (commented 4000) -> ``6000``
- ``ombh2 0.0226 -> 0.0224``; ``omch2 0.112 -> 0.1200``; ``omnuh2 0.00064 -> 0.6451439e-3``; ``hubble 70 -> 67.36``
- ``scalar_amp(1) 2.1e-9 -> 2.196e-9``; ``scalar_spectral_index(1) 0.96 -> 0.9655``; ``initial_ratio(1) 1 -> 0``
- ``re_optical_depth 0.09 -> 0.05``; ``re_delta_redshift 1.5 -> 0.5``
- ``transfer_kmax 2 -> 5``; ``transfer_k_per_logint 0 -> 8``
- ``accurate_BB F -> T``; ``bispectrum_nfields 1 -> 2``
- ``massless_neutrinos = 2.046``, ``massive_neutrinos = 1``, ``share_delta_neff = T`` (unchanged values, reordered)

Removed keys: ``wa``, ``use_tabulated_w``, ``wafile`` (commented PPF extras), ``derived_parameters`` (driver default T applies), ``nu_mass_degeneracies`` (blank line removed), ``bispectrum_export_alpha_beta``, cosmorec comment block.

Comment changes: ``do_nonlinear`` help no longer lists option 3; ``initial_condition`` help adds "axion isocurvature = 6"; ``initial_vector`` example changed ``-1 0 0 0 0`` -> ``1 0 0 0 0 1`` (6 entries to match ``initial_iso_axion=6``).

6. Port-design notes, surprises, risks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Driver logic must move into the library for 1.6.7.** Everything in §1.5-§1.9 and §2 happens before ``CAMB_GetResults`` in the ini driver, so the Python wrapper path (``camb.set_params -> CAMBdata.calc_*``) would bypass it entirely. Natural target: a ULA component class whose ``Init`` (called from ``CAMBdata%SetParams``) runs ``w_evolve`` + the dfac iterations, plus ``CAMBparams`` validation/derivation in ``set_cosmology``-equivalent code for the omdah2/axfrac logic.
2. **The DM/DE discriminator is** ``m_ovH0 = m_ax/H0_eV >= 10``, identically used in driver (density assignment), halofit (matter vs DE budget) and reportedly ``outtransf`` in modules.f90. Keep ONE definition.
3. **dfac (switch at m/H=10) is internal and is MUTATED by two driver-level iteration schemes** (ETA-phase shooting with target ``2beta = 7.08*pi``, tol ``0.02*pi``; and recombination-window skip z∈(800,1300) via ``dfac_skip``), each re-running the full background solver with ``ntable = nint(dfac*100)+1``. These loops are pure physics/accuracy logic and must be ported as-is (they protect Cl accuracy for 1e-27..1e-25 eV masses).
4. **P%omegav closure includes radiation** (``P%omegah2_rad`` = photons + massless nu only); the driver duplicates CAMBparams_Set neutrino fixups to get N_eff right before that. In 1.6.7 use the native ``grhog/grhornomass`` machinery instead — port the intent (``Omega_de = 1 - .. - Omega_ax``), not the code.
5. **Halofit:** only three real deltas — default version = ``halofit_original`` (Takahashi declared unstable for axions, DM15 comment), ``omm0 += omegaax`` when m/H0>=10, and ``om_m/om_v`` called with ``CP%omegav + CP%omegaax`` when m/H0<10 (new ``omv0_axion``). Everything else in the file is stock Sept-2014 halofit_ppf already superseded by 1.6.7. The shipped ``AxiECAMB/halofit.f90`` is byte-identical to OLDCAMB's and unused (Makefile_main: ``NONLINEAR ?= halofit_ppf``).
6. **Isocurvature is dead code in this release** (forced ``axion_isocurvature=.false.`` at driver:321-324) but the full two-pass machinery exists: mode 6 ICs, ``amp_i = Hinf**2/(pi**2*phiinit**2)``, ``n_iso = 1 - r_val/8`` with ``r_val = 2*Hinf**2/(pi**2*As)``, spectra added in power. ``alpha_ax`` in params.ini is a dead key (computed, never read). Decide early whether to port it or strip it.
7. **Dead keys/quirks:** ``tens_ratio`` never read; ``use_physical=F`` branch divides ``omega_axion`` by h^2 (inconsistent with the key name) and skips setting ``m_ovH0``/``axfrac``/``Hinf`` — recommend requiring physical densities; ``do_nonlinear==3`` escapes the "not well-tested" warning.
8. **Build:** the only structural build change to port is adding the axion background module with a dependency on the parameters module (axion_background.o: modules.o), compiled before equations. ``cmbmainOMP.f90`` is shipped but never compiled — confirm with the cmbmain analyst that nothing in it is load-bearing.
9. **init_massive_nu(P%omegan/=0) must run before the axion background solver** (driver:519) because ``w_evolve`` evaluates massive-nu background densities; in 1.6.7 ensure the nu interpolation tables are initialized before the ULA ``Init``.

