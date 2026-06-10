    !!!!!!!!!!!!!!!! AxionBackground !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Ultralight-axion (ULA) background solver for AxiECAMB ported to modern CAMB.
    !
    ! This module computes the evolution of the background in the presence of an axion
    ! field with a quadratic potential. The exact Klein-Gordon (KG) equation is solved
    ! from deep radiation domination until m = dfac*H (dfac ~ 10, retuned internally),
    ! after which the axion is treated as an effective fluid (EFA) that tracks the
    ! cycle average of the fast field oscillations, following
    !   Hlozek et al 2014 (arXiv:1410.2896)  [axionCAMB heritage]
    !   Passaglia & Hu 2022 (arXiv:2201.10238)
    !   Liu, Hu et al 2024 (arXiv:2412.15192)  [AxiECAMB - cite this when using]
    !
    ! KG equation in conformal time: phidot_dot + 2 H phidot + m^2 a^2 phi = 0
    ! transformed with phi = sqrt(3/(4 pi G)) v1, dphi/dtau = H0 sqrt(3/(4 pi G)) v2.
    ! The dimensionless mass is m_ovH0 = m/H0; "lh" returns conformal aH in units of
    ! 100 km/s/Mpc (so lh/(a*hnot) = H/H0).
    !
    ! Ported from AxiECAMB axion_background.F90 (CAMB Nov13 base): the module globals
    ! became members of TAxionModel, the driver-level dfac orchestration from
    ! inidriver_axion.F90 (oscillation-phase targeting and recombination skip) moved
    ! into TAxionModel_Solve, and the background pieces (radiation, neutrinos, Lambda)
    ! are passed in from CAMBdata so that H(a) here is identical to the rest of CAMB.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    module AxionBackground
    use precision
    use constants
    use classes
    use config
    use MassiveNu
    use Interpolation, only : spline, SPLINE_DANGLE
    implicit none
    private

    !maximum number of shooting/bisection iterations for the initial field value
    integer, parameter :: nphi = 150
    !8*pi*G*rho/c^2 per unit Omega h^2, in Mpc^-2 (= grhocrit/h^2 = 3*(100 km/s/Mpc)^2/c^2)
    real(dl), parameter :: grhom_ovh2 = 3.0d0*1.d10/(c**2.0d0)
    !reduced Planck mass in GeV, used to normalize H_inflation (AxiECAMB constants.f90 mplanck)
    real(dl), parameter :: mplanck_GeV = 2.435d18

    !Background pieces from CAMBdata, so the Friedmann equation used in the KG solve
    !is identical to the one CAMB uses everywhere else (via grho_no_de/dtauda).
    type TAxionBgInputs
        real(dl) :: hnot = 0._dl                 !H0/100
        real(dl) :: omegah2_b = 0._dl            !Omega_b h^2
        real(dl) :: omegah2_dm = 0._dl           !standard CDM Omega_c h^2 (after axion split)
        real(dl) :: omegah2_lambda = 0._dl       !dark energy treated as Lambda in the solver
        real(dl) :: omegah2_nu = 0._dl           !massive neutrinos Omega_nu h^2
        real(dl) :: omk = 0._dl
        real(dl) :: omegah2_rad = 0._dl          !photons + massless neutrinos
        real(dl) :: omegah2_gamma = 0._dl        !photons only (for aeq_LCDM)
        real(dl) :: nu_massless_degeneracy = 0._dl
        real(dl) :: sum_nu_mass_degeneracies = 0._dl
        integer :: Nu_mass_eigenstates = 0
        real(dl), allocatable :: nu_masses(:)       !CAMBdata nu_masses (m c^2/k_B T_nu0)
        real(dl), allocatable :: lhsqcont_massive(:)!grhormass(i)/grhocrit*h^2: Omega h^2 per unit nu_rho
        real(dl) :: ScalarPowerAmp = 0._dl       !A_s, for isocurvature normalization only
        integer :: FeedbackLevel = 0
    end type TAxionBgInputs

    type, extends(TCambComponent) :: TAxionModel
        !---- input parameters ----
        logical  :: active = .false.       !axion component switched on
        real(dl) :: m_ax = 0._dl           !axion mass in eV
        real(dl) :: omaxh2 = 0._dl         !Omega_ax h^2 (input when use_axfrac=F; effective value stored in omaxh2_eff)
        logical  :: use_axfrac = .false.   !use (omdah2, axfrac) instead of omaxh2
        real(dl) :: omdah2 = 0._dl         !total dark matter Omega h^2 (when use_axfrac=T and m/H0>=10)
        real(dl) :: axfrac = 0._dl         !axion fraction of DM (m/H0>=10) or of DE (m/H0<10)
        logical  :: axion_isocurvature = .false. !accepted but force-disabled (AxiECAMB v1.0 parity)
        real(dl) :: Hinf = 13.7_dl         !log10 of H_inflation in GeV (converted to H_inf/M_pl at use)
        real(dl) :: dfac = 10._dl          !input KG->EFA switch threshold m = dfac*H
        !---- derived parameters / state (treat as read-only) ----
        real(dl) :: H0_eV = 0._dl          !H0 in eV
        real(dl) :: m_ovH0 = 0._dl         !m/H0
        logical  :: is_de_like = .false.   !m/H0 < 10: KG solved to a=1, no fluid switch
        logical  :: has_switch = .false.   !a KG->EFA switch exists at a_osc (<canonically <=1)
        real(dl) :: dfac_used = 10._dl     !dfac after internal retuning (phase targeting/recombination skip)
        real(dl) :: a_osc = 1._dl          !scale factor of the KG->EFA switch
        real(dl) :: tau_osc = 0._dl        !conformal time of the switch (set by CAMBdata after Solve)
        real(dl) :: aeq = 0._dl            !matter-radiation equality incl. exact axion density
        real(dl) :: aeq_LCDM = 0._dl       !analytic equality assuming axions scale as matter
        real(dl) :: phiinit = 0._dl        !initial field value in reduced Planck units
        real(dl) :: ah_osc = 0._dl         !instantaneous conformal aH at a_osc (units 100 km/s/Mpc)
        real(dl) :: ahosc_ETA = 0._dl      !cycle-averaged <aH> at a_osc from the EFA matching
        real(dl) :: A_coeff = 0._dl        !EFA matching coefficient (see arXiv:2412.15192)
        real(dl) :: A_coeff_alt = 0._dl    !A_coeff + 2 aH/(a m), used by the perturbations
        real(dl) :: tvarphi_c = 0._dl      !WKB cos/sin projections of (phi, phi') at the switch
        real(dl) :: tvarphi_s = 0._dl
        real(dl) :: tvarphi_cp = 0._dl
        real(dl) :: tvarphi_sp = 0._dl
        real(dl) :: wEFA_c = 9._dl/8._dl   !coefficient of w_ax = wEFA_c (H/m)^2 (iterated to consistency)
        real(dl) :: rhorefp_ovh2 = 0._dl   !EFA cycle-averaged Omega_ax(a_osc) (no h^2)
        real(dl) :: Prefp = 0._dl          !EFA cycle-averaged pressure at a_osc, Omega(a) h^2 units
        real(dl) :: wcorr_coeff = 0._dl    !<aH>_osc a_osc/(m_ovH0 hnot): (H/m) = wcorr_coeff/a^2 after switch
        real(dl) :: dfac_skip = 0._dl      !dfac that relocates the switch to past z=800
        real(dl) :: a_skip = 0._dl         !1/801: recombination-skip window upper edge
        real(dl) :: a_skipst = 0._dl       !1/1301: recombination-skip window lower edge
        real(dl) :: opac_tauosc = 0._dl    !opacity at tau_osc (filled by Thermo_Init)
        real(dl) :: expmmu_tauosc = 0._dl  !exp(-kappa) at tau_osc (filled by Thermo_Init)
        real(dl) :: amp_i = 0._dl          !isocurvature amplitude (inactive in this release)
        real(dl) :: r_val = 0._dl
        real(dl) :: alpha_ax = 0._dl
        real(dl) :: omaxh2_eff = 0._dl     !Omega_ax h^2 actually used after axfrac logic
        real(dl) :: omegah2_rad = 0._dl    !radiation Omega h^2 used in the solve (diagnostic)
        real(dl) :: hsq = 0._dl            !(H0/100)^2 cached from Solve
        integer  :: ntable = 1001          !background table size = nint(dfac_used*100)+1
        real(dl) :: a_table_min = 0._dl    !first scale factor in the tables
        !---- background tables (abscissa loga_table is log10(a)) ----
        real(dl), allocatable :: loga_table(:)
        real(dl), allocatable :: phinorm_table(:)            !v1 = phi sqrt(4 pi G/3)
        real(dl), allocatable :: phidotnorm_table(:)         !v2 = (dphi/dtau) sqrt(4 pi G/3)/H0
        real(dl), allocatable :: phinorm_table_ddlga(:)      !spline second derivatives
        real(dl), allocatable :: phidotnorm_table_ddlga(:)
        real(dl), allocatable :: rhoaxh2ovrhom_logtable(:)   !log10(Omega_ax(a) h^2), instantaneous KG
        real(dl), allocatable :: rhoaxh2ovrhom_logtable_buff(:)
    contains
    procedure, nopass :: PythonClass => TAxionModel_PythonClass
    procedure, nopass :: SelfPointer => TAxionModel_SelfPointer
    procedure :: ReadParams => TAxionModel_ReadParams
    procedure :: Validate => TAxionModel_Validate
    procedure :: Solve => TAxionModel_Solve
    procedure :: RhoaxH2AtA => TAxionModel_RhoaxH2AtA
    procedure :: GrhoAx => TAxionModel_GrhoAx
    procedure :: FieldValsAta => TAxionModel_FieldValsAta
    procedure, private :: w_evolve => TAxionModel_w_evolve
    procedure, private :: auxiIC => TAxionModel_auxiIC
    procedure :: get_phase_info => TAxionModel_get_phase_info
    end type TAxionModel

    public TAxionModel, TAxionBgInputs, mplanck_GeV, grhom_ovh2
    contains

    function TAxionModel_PythonClass()
    character(LEN=:), allocatable :: TAxionModel_PythonClass
    TAxionModel_PythonClass = 'AxionModel'
    end function TAxionModel_PythonClass

    subroutine TAxionModel_SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type(TAxionModel), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType
    end subroutine TAxionModel_SelfPointer

    subroutine TAxionModel_ReadParams(this, Ini)
    use IniObjects
    class(TAxionModel) :: this
    class(TIniFile), intent(in) :: Ini

    this%m_ax = Ini%Read_Double('m_ax', 0._dl)
    if (this%m_ax < 0) this%m_ax = 10**this%m_ax !negative input interpreted as log10(m_ax/eV)
    this%use_axfrac = Ini%Read_Logical('use_axfrac', .false.)
    if (this%use_axfrac) then
        this%omdah2 = Ini%Read_Double('omdah2')
        this%axfrac = Ini%Read_Double('axfrac')
    else
        this%omaxh2 = Ini%Read_Double('omaxh2', 0._dl)
    end if
    this%dfac = Ini%Read_Double('axion_dfac', 10._dl)
    this%Hinf = Ini%Read_Double('Hinf', 13.7_dl)  !log10 of H_inflation in GeV
    this%axion_isocurvature = Ini%Read_Logical('axion_isocurvature', .false.)
    if (this%axion_isocurvature) then
        write(*,*) 'WARNING: axion isocurvature disabled in this release, proceeding without'
        this%axion_isocurvature = .false.
    end if
    this%active = this%m_ax > 0 .and. (this%use_axfrac .and. this%axfrac > 0 &
        .or. .not. this%use_axfrac .and. this%omaxh2 > 0)
    end subroutine TAxionModel_ReadParams

    subroutine TAxionModel_Validate(this, OK)
    class(TAxionModel), intent(in) :: this
    logical, intent(inout) :: OK

    if (this%active) then
        if (this%m_ax <= 0) then
            write(*,*) 'TAxionModel: m_ax must be positive (eV)'
            OK = .false.
        end if
        if (this%use_axfrac) then
            if (this%axfrac < 0 .or. this%axfrac > 1 .or. this%omdah2 < 0) then
                write(*,*) 'TAxionModel: need 0 <= axfrac <= 1 and omdah2 >= 0'
                OK = .false.
            end if
        else if (this%omaxh2 < 0) then
            write(*,*) 'TAxionModel: omaxh2 must be >= 0'
            OK = .false.
        end if
        if (this%dfac <= 0) then
            write(*,*) 'TAxionModel: axion_dfac must be positive'
            OK = .false.
        end if
    end if
    end subroutine TAxionModel_Validate

    !> Omega_ax(a) h^2: instantaneous KG density from the table for a <= a_osc,
    !> analytic EFA cycle-averaged density beyond the switch (the table is never used
    !> past a_osc; in the oscillating case it only extends to ~1.1 a_osc).
    function TAxionModel_RhoaxH2AtA(this, a) result(rhoaxh2)
    class(TAxionModel), intent(in) :: this
    real(dl), intent(in) :: a
    real(dl) rhoaxh2, lgrho

    if (this%has_switch .and. a > this%a_osc) then
        rhoaxh2 = this%rhorefp_ovh2*this%hsq*((this%a_osc/a)**3.0d0)* &
            dexp((this%wcorr_coeff**2.0d0)*3.0d0*this%wEFA_c* &
            (1.0d0/(a**4.0d0) - 1.0d0/(this%a_osc**4.0d0))/4.0d0)
    else if (a < this%a_table_min) then
        !field frozen well before the table starts: density constant
        rhoaxh2 = 10._dl**this%rhoaxh2ovrhom_logtable(1)
    else
        call spline_out(this%loga_table, this%rhoaxh2ovrhom_logtable, &
            this%rhoaxh2ovrhom_logtable_buff, this%ntable, dlog10(a), lgrho)
        rhoaxh2 = 10._dl**lgrho
    end if
    end function TAxionModel_RhoaxH2AtA

    !> 8*pi*G*rho_ax*a^4 in Mpc^-2 (same convention as CAMBdata grho_no_de terms)
    function TAxionModel_GrhoAx(this, a) result(grhoax_a4)
    class(TAxionModel), intent(in) :: this
    real(dl), intent(in) :: a
    real(dl) grhoax_a4

    if (.not. this%active .or. a <= 0._dl) then
        grhoax_a4 = 0._dl
    else
        grhoax_a4 = grhom_ovh2*this%RhoaxH2AtA(a)*a**4.0d0
    end if
    end function TAxionModel_GrhoAx

    !> Background field variables v1, v2 at scale factor a (clamped to the table range).
    subroutine TAxionModel_FieldValsAta(this, a, v1, v2)
    class(TAxionModel), intent(in) :: this
    real(dl), intent(in) :: a
    real(dl), intent(out) :: v1, v2
    real(dl) lga

    lga = dlog10(max(a, this%a_table_min))
    lga = min(lga, this%loga_table(this%ntable))
    call spline_out(this%loga_table, this%phinorm_table, this%phinorm_table_ddlga, &
        this%ntable, lga, v1)
    call spline_out(this%loga_table, this%phidotnorm_table, this%phidotnorm_table_ddlga, &
        this%ntable, lga, v2)
    end subroutine TAxionModel_FieldValsAta

    !> Analytic WKB estimate of the field oscillation phase at the switch, used for
    !> the dfac phase targeting (AxiECAMB get_phase_info, verbatim).
    subroutine TAxionModel_get_phase_info(this, y_beta, beta_coeff, movHETA, beta2x)
    class(TAxionModel), intent(in) :: this
    real(dl), intent(out) :: y_beta, beta_coeff, movHETA, beta2x

    y_beta = this%a_osc/this%aeq_LCDM
    movHETA = this%dfac_used*this%ah_osc/this%ahosc_ETA
    beta_coeff = (4._dl*(y_beta**2 - y_beta - 2.0_dl + 2.0_dl*sqrt(1.0_dl + y_beta)))/(3._dl*(y_beta**2))
    beta2x = movHETA*beta_coeff - const_pi*3._dl*(1.0_dl + y_beta)/(4.0_dl + 3.0_dl*y_beta)
    end subroutine TAxionModel_get_phase_info

    !> Master driver: solves the KG background (shooting for the initial field value to
    !> match omaxh2_eff), then retunes dfac: (a) oscillation-phase targeting at
    !> 2 beta = 7.08 pi for DM-like axions with an early switch, (b) recombination skip
    !> pushing the switch out of z in (800, 1300). Ported from AxiECAMB
    !> inidriver_axion.F90:508-636.
    subroutine TAxionModel_Solve(this, BG, badflag)
    class(TAxionModel), intent(inout) :: this
    type(TAxionBgInputs), intent(in) :: BG
    integer, intent(out) :: badflag
    real(dl) twobeta_tgt, twobeta_new, twobeta_old1, twobeta_old2, beta_tol
    real(dl) y_phase, beta_coeff, movHETA_new, movHETA_beta, hETA_beta
    real(dl) hosc_old1, hosc_old2, hosc_new, hETA_old1, hETA_old2
    integer iter_dfacETA, iter_dfac

    badflag = 0
    this%hsq = BG%hnot**2.0d0
    this%omegah2_rad = BG%omegah2_rad
    this%a_skip = 1._dl/(800._dl + 1._dl)
    this%a_skipst = 1._dl/(1300._dl + 1._dl)
    this%dfac_skip = 0._dl
    this%dfac_used = this%dfac          !reset retuning so repeated SetParams calls are idempotent
    this%ntable = nint(this%dfac_used*100) + 1

    call this%w_evolve(BG, badflag)
    if (badflag /= 0 .or. global_error_flag /= 0) return

    !Oscillation-phase targeting: only when the switch is early enough to matter for the
    !photon eta at recombination (empirically tuned 0.03 threshold), the axion is DM-like
    !and light, and dfac has not already been raised
    if (this%dfac_used < 23._dl .and. this%m_ovH0 >= 10._dl .and. this%m_ax < 1.e-25_dl &
        .and. this%a_osc*(this%omaxh2_eff/(BG%omegah2_dm + this%omaxh2_eff))/this%aeq_LCDM > 0.03_dl &
        .and. this%a_osc < this%a_skipst) then
        twobeta_tgt = 7.08_dl*const_pi
        !First guess considering matter-radiation equality in LCDM
        this%dfac_used = twobeta_tgt + 0.75_dl*const_pi - twobeta_tgt**2/ &
            (4._dl*(twobeta_tgt + 2._dl*this%m_ovH0*(this%aeq_LCDM**1.5_dl)/ &
            sqrt(2._dl*(BG%omegah2_dm + BG%omegah2_b + BG%omegah2_nu + this%omaxh2_eff)/this%hsq)))
        this%ntable = nint(this%dfac_used*100) + 1
        call this%w_evolve(BG, badflag)
        if (badflag /= 0 .or. global_error_flag /= 0) return
        call this%get_phase_info(y_phase, beta_coeff, movHETA_new, twobeta_new)

        !Rough guess of target ETA values using beta, useful in bracket finding
        movHETA_beta = (twobeta_tgt + const_pi*3._dl*(1.0_dl + y_phase)/(4.0_dl + 3.0_dl*y_phase))/beta_coeff
        hETA_beta = (this%dfac_used/movHETA_beta)*(this%ah_osc/this%a_osc)
        beta_tol = 2.e-2_dl*const_pi
        if (abs(twobeta_new - twobeta_tgt) > beta_tol) then
            !Find the bracket by shooting on hosc = ah_osc/a_osc
            hosc_old1 = this%ah_osc/this%a_osc
            hETA_old1 = this%ahosc_ETA/this%a_osc
            twobeta_old1 = twobeta_new
            hosc_new = 2._dl*(hETA_beta - hETA_old1) + hosc_old1
            this%dfac_used = this%dfac_used*((this%ah_osc/this%a_osc)/hosc_new)
            this%ntable = nint(this%dfac_used*100) + 1
            call this%w_evolve(BG, badflag)
            if (badflag /= 0 .or. global_error_flag /= 0) return
            call this%get_phase_info(y_phase, beta_coeff, movHETA_new, twobeta_new)
            hosc_old2 = this%ah_osc/this%a_osc
            hETA_old2 = this%ahosc_ETA/this%a_osc
            twobeta_old2 = twobeta_new

            iter_dfacETA = 1
            do while (iter_dfacETA < 500)
                if (abs(twobeta_new - twobeta_tgt) < beta_tol) then
                    iter_dfacETA = -1
                    exit
                else if ((twobeta_old2 - twobeta_tgt)*(twobeta_old1 - twobeta_tgt) < 0._dl) then
                    exit
                else
                    !Still use hETA_beta-hETA_old1 to keep the shooting step from diminishing
                    hosc_new = 1._dl*(hETA_beta - hETA_old1) + hosc_old2
                    this%dfac_used = this%dfac_used*((this%ah_osc/this%a_osc)/hosc_new)
                    this%ntable = nint(this%dfac_used*100) + 1
                    call this%w_evolve(BG, badflag)
                    if (badflag /= 0 .or. global_error_flag /= 0) return
                    call this%get_phase_info(y_phase, beta_coeff, movHETA_new, twobeta_new)
                    hosc_old2 = this%ah_osc/this%a_osc
                    hETA_old2 = this%ahosc_ETA/this%a_osc
                    twobeta_old2 = twobeta_new
                end if
            end do

            if (iter_dfacETA /= -1) then
                !After the bracket is found, bisect on hosc
                hosc_new = (hosc_old1 + hosc_old2)/2._dl
                this%dfac_used = this%dfac_used*((this%ah_osc/this%a_osc)/hosc_new)
                this%ntable = nint(this%dfac_used*100) + 1
                call this%w_evolve(BG, badflag)
                if (badflag /= 0 .or. global_error_flag /= 0) return
                call this%get_phase_info(y_phase, beta_coeff, movHETA_new, twobeta_new)

                iter_dfacETA = 1
                do while (iter_dfacETA < 500)
                    if (abs(twobeta_new - twobeta_tgt) < beta_tol) then
                        exit
                    else
                        if ((twobeta_new - twobeta_tgt)*(twobeta_old1 - twobeta_tgt) < 0._dl) then
                            hosc_old2 = this%ah_osc/this%a_osc
                            hETA_old2 = this%ahosc_ETA/this%a_osc
                            twobeta_old2 = twobeta_new
                        else
                            hosc_old1 = this%ah_osc/this%a_osc
                            hETA_old1 = this%ahosc_ETA/this%a_osc
                            twobeta_old1 = twobeta_new
                        end if
                        hosc_new = (hosc_old1 + hosc_old2)/2._dl
                        this%dfac_used = this%dfac_used*((this%ah_osc/this%a_osc)/hosc_new)
                        this%ntable = nint(this%dfac_used*100) + 1
                        call this%w_evolve(BG, badflag)
                        if (badflag /= 0 .or. global_error_flag /= 0) return
                        call this%get_phase_info(y_phase, beta_coeff, movHETA_new, twobeta_new)
                    end if
                end do
            end if
        end if
    end if

    !Recombination skip: push the switch past z=800 if it would land in z in (800,1300)
    do iter_dfac = 1, 500
        if (this%a_osc < this%a_skip*(1._dl - 1.e-2_dl) .and. this%a_osc >= this%a_skipst) then
            this%dfac_used = this%dfac_skip
            this%ntable = nint(this%dfac_used*100) + 1
            call this%w_evolve(BG, badflag)
            if (badflag /= 0 .or. global_error_flag /= 0) return
        else
            exit
        end if
    end do
    if (iter_dfac > 500 .and. this%a_osc < this%a_skip) then
        write(*,*) 'Warning: maximum iteration reached, but aosc still not skipped sufficiently: ', &
            'a_osc, a_skip', this%a_osc, this%a_skip
    end if

    this%is_de_like = this%m_ovH0 < 10._dl
    this%has_switch = .not. this%is_de_like
    if (this%a_osc >= 1.0d0) this%a_osc = 1.0d0   !cap the sentinel (AxiECAMB modules.f90:417)
    this%wcorr_coeff = this%ahosc_ETA*this%a_osc/(this%m_ovH0*BG%hnot)
    this%a_table_min = 10._dl**this%loga_table(1)

    if (BG%FeedbackLevel > 0 .and. this%has_switch) then
        write(*,'(" Axion KG->EFA switch: a_osc = ",E13.5," (z = ",F10.2,"), dfac = ",F8.3)') &
            this%a_osc, 1/this%a_osc - 1, this%dfac_used
    end if
    end subroutine TAxionModel_Solve

    !> KG background solve + shooting for the initial field value (AxiECAMB w_evolve).
    subroutine TAxionModel_w_evolve(this, BG, badflag)
    class(TAxionModel), intent(inout) :: this
    type(TAxionBgInputs), intent(in) :: BG
    integer, intent(out) :: badflag
    integer i, j, ntable
    !contribution to H^2/(100 km/s/Mpc)^2 of massive and massless neutrinos
    real(dl) lhsqcont_massive(max(BG%Nu_mass_eigenstates,1)), lhsqcont_massless
    real(dl) omegah2_m, omegah2_lambda, omnuh2
    real(dl) omegah2_regm, omk
    real(dl) omegah2_b, omegah2_dm, omegah2_ax, maxion_twiddle, hnot
    real(dl) hsq, regzeq
    real(dl) dfac
    real(dl), allocatable :: a_arr(:), v_vec(:,:), littlehfunc(:), diagnostic(:)
    real(dl), allocatable :: rhoaxh2_ov_rhom(:)
    real(dl), allocatable :: f_arr(:), eq_arr(:)
    real(dl), allocatable :: v_buff(:), abuff(:), eq_arr_buff(:)
    real(dl) phiosc, phidosc, laosc
    real(dl) rhorefp, Prefp
    real(dl) v1_ref, v2_ref
    real(dl) littlehauxi, wcorr_coeff, omaxh2_wcorr, A_coeff
    real(dl) tvarphi_c, tvarphi_cp, tvarphi_s, tvarphi_sp
    real(dl) lh_skip
    real(dl) v1_initguess(3), omaxh2_guess(3), aosc_guess(3)
    integer iter_c
    logical bisec_bracketed
    real(dl) vtwiddle_init
    real(dl) d1, d2
    real(dl) a_init, a_m, a_lambda, a_rel, as_scalar, as_rad, as_matt, a_final
    real(dl) dloga, log_a_final, log_a_init
    real(dl), allocatable :: cmat(:,:), kvec(:,:), kfinal(:), svec(:), avec(:)

    badflag = 0
    ntable = this%ntable
    if (allocated(this%loga_table)) deallocate(this%loga_table)
    if (allocated(this%phinorm_table)) deallocate(this%phinorm_table)
    if (allocated(this%phidotnorm_table)) deallocate(this%phidotnorm_table)
    if (allocated(this%phinorm_table_ddlga)) deallocate(this%phinorm_table_ddlga)
    if (allocated(this%phidotnorm_table_ddlga)) deallocate(this%phidotnorm_table_ddlga)
    if (allocated(this%rhoaxh2ovrhom_logtable)) deallocate(this%rhoaxh2ovrhom_logtable)
    if (allocated(this%rhoaxh2ovrhom_logtable_buff)) deallocate(this%rhoaxh2ovrhom_logtable_buff)

    hnot = BG%hnot
    hsq = hnot**2.0d0
    omegah2_dm = BG%omegah2_dm
    omegah2_b = BG%omegah2_b
    omnuh2 = BG%omegah2_nu
    omegah2_lambda = BG%omegah2_lambda
    omegah2_ax = this%omaxh2_eff
    maxion_twiddle = this%m_ovH0
    omegah2_regm = omegah2_dm + omegah2_b
    omegah2_m = omegah2_regm + omegah2_ax

    !Neutrino/radiation pieces come from CAMBdata so H(a) here matches grho_no_de exactly
    lhsqcont_massless = BG%omegah2_rad - BG%omegah2_gamma
    lhsqcont_massive = 0._dl
    if (BG%Nu_mass_eigenstates > 0) &
        lhsqcont_massive(1:BG%Nu_mass_eigenstates) = BG%lhsqcont_massive(1:BG%Nu_mass_eigenstates)

    dfac = this%dfac_used
    this%wEFA_c = 9._dl/8._dl   !analytic RD value; iterated to self-consistency in auxiIC

    !16-stage Fehlberg classical 8th-order Runge-Kutta tableau (fixed step in ln a):
    !page 75 of E Fehlberg, NASA Huntsville 1968, http://hdl.handle.net/2060/19680027281
    !A reasonably accurate integrator is required to accurately obtain the adiabatic
    !sound speed at early times; this integrator with ~100 points per dfac unit was
    !needed to avoid exciting a non-physically large low-l ISW effect.
    allocate(cmat(16,16), kvec(2,16))
    allocate(kfinal(2), svec(16), avec(16))
    avec = 0.0d0
    kvec = 0.0d0
    avec(1)=0.4436894037649818d0
    avec(2)=0.6655341056474727d0
    avec(3)=0.9983011584712091d0
    avec(4)=0.31550d0
    avec(5)=0.5054410094816906d0
    avec(6)=0.1714285714285714d0
    avec(7)=0.8285714285714285d0
    avec(8)=0.6654396612101156d0
    avec(9)=0.2487831796806265d0
    avec(10)=0.1090d0
    avec(11)=0.8910d0
    avec(12)=0.3995d0
    avec(13)=0.6005d0
    avec(14)=1.0d0
    avec(15)=0.0d0
    avec(16)=1.0d0
    cmat = 0.0d0
    cmat(1,1)=avec(1)
    cmat(2,1)=0.1663835264118681d0
    cmat(2,2)=0.49915057d0
    cmat(3,1)=0.24957528d0
    cmat(3,3)=0.74872586d0
    cmat(4,1)=0.20661891d0
    cmat(4,3)=0.17707880d0
    cmat(4,4)=-0.68197715d-1
    cmat(5,1)=0.10927823d0
    cmat(5,4)=0.40215962d-2
    cmat(5,5)=0.39214118d0
    cmat(6,1)=0.98899281d-1
    cmat(6,4)=0.35138370d-1
    cmat(6,5)=0.12476099d0
    cmat(6,6)=-0.55745546d-1
    cmat(7,1)=-0.36806865d0
    cmat(7,5)=-0.22273897d1
    cmat(7,6)=0.13742908d1
    cmat(7,7)=0.20497390d1
    cmat(8,1)=0.45467962d-1
    cmat(8,6)=0.32542131d0
    cmat(8,7)=0.28476660d0
    cmat(8,8)=0.97837801d-2
    cmat(9,1)=0.60842071d-1
    cmat(9,6)=-0.21184565d-1
    cmat(9,7)=0.19596557d0
    cmat(9,8)=-0.42742640d-2
    cmat(9,9)=0.17434365d-1
    cmat(10,1)=0.54059783d-1
    cmat(10,7)=.11029325d0
    cmat(10,8)=-.12565008d-2
    cmat(10,9)=0.36790043d-2
    cmat(10,10)=-.57780542d-1
    cmat(11,1)=.12732477d0
    cmat(11,8)=0.11448805d0
    cmat(11,9)=0.28773020d0
    cmat(11,10)=0.50945379d0
    cmat(11,11)=-0.14799682d0
    cmat(12,1)=-0.36526793d-2
    cmat(12,6)=0.81629896d-1
    cmat(12,7)=-0.38607735d0
    cmat(12,8)=0.30862242d-1
    cmat(12,9)=-0.58077254d-1
    cmat(12,10)=0.33598659d0
    cmat(12,11)=0.41066880d0
    cmat(12,12)=-0.11840245d-1
    cmat(13,1)=-0.12375357d1
    cmat(13,6)=-0.24430768d2
    cmat(13,7)=0.54779568d0
    cmat(13,8)=-0.44413863d1
    cmat(13,9)=0.10013104d2
    cmat(13,10)=-0.14995773d2
    cmat(13,11)=0.58946948d1
    cmat(13,12)=0.17380377d1
    cmat(13,13)=0.27512330d2
    cmat(14,1)=-0.35260859d0
    cmat(14,6)=-0.18396103d0
    cmat(14,7)=-0.65570189d0
    cmat(14,8)=-.39086144d0
    cmat(14,9)=0.26794646d0
    cmat(14,10)=-0.10383022d1
    cmat(14,11)=0.16672327d1
    cmat(14,12)=0.49551925d0
    cmat(14,13)=.11394001d1
    cmat(14,14)=0.51336696d-1
    cmat(15,1)=0.10464847d-2
    cmat(15,9)=-0.67163886d-2
    cmat(15,10)=0.81828762d-2
    cmat(15,11)=-0.42640342d-2
    cmat(15,12)=0.280090294741d-3
    cmat(15,13)=-0.87835333d-2
    cmat(15,14)=0.10254505d-1
    cmat(16,1)=-0.13536550d1
    cmat(16,6)=-0.18396103d0
    cmat(16,7)=-0.65570189d0
    cmat(16,8)=-0.39086144d0
    cmat(16,9)=0.27466285d0
    cmat(16,10)=-0.10464851d1
    cmat(16,11)=0.16714967d1
    cmat(16,12)=0.49523916d0
    cmat(16,13)=0.11481836d1
    cmat(16,14)=0.41082191d-1
    cmat(16,16)=1.0d0

    svec = 0.0d0
    svec(1)=0.32256083d-1
    svec(9)=0.25983725d0
    svec(10)=0.92847805d-1
    svec(11)=.16452330d0
    svec(12)=0.176659510d0
    svec(13)=0.23920102d0
    svec(14)=0.39484274d-2
    svec(15)=0.3072649547580d-1

    !Initial aosc sentinel: 15 = "oscillation not (yet) found"
    aosc_guess = (/15.0d0, 15.0d0, 15.0d0/)
    omk = 1.0d0 - (omegah2_m + BG%omegah2_rad + omegah2_lambda + omnuh2)/hsq

    !Determine a range of epochs to be well before in starting the field evolution
    as_matt = (omegah2_regm/(maxion_twiddle**2.0d0))**(1.0d0/3.0d0)
    as_rad = (BG%omegah2_rad/(maxion_twiddle**2.0d0))**(1.0d0/4.0d0)
    a_m = (BG%omegah2_rad/(omegah2_regm))
    a_rel = 10.0d0
    a_lambda = (BG%omegah2_rad/omegah2_lambda)**(0.25d0)
    as_scalar = (omegah2_ax/(maxion_twiddle**2.0d0))**(1.0d0/3.0d0)
    a_init = min(a_rel, a_lambda, a_m, as_matt, as_rad, as_scalar)*1.d-8
    a_final = 1.0d0
    log_a_init = dlog(a_init)
    log_a_final = dlog(a_final)
    dloga = (log_a_final - log_a_init)/(dble(ntable - 1))
    allocate(a_arr(ntable))
    a_arr(1) = dexp(log_a_init)
    allocate(this%loga_table(ntable))
    !NOTE: loga_table holds ln(a) during the integration and is converted to log10(a)
    !below (AxiECAMB heritage) - all spline lookups against the final tables use log10.
    this%loga_table(1) = log_a_init
    do i = 2, ntable
        this%loga_table(i) = log_a_init + dloga*dble(i-1)
        a_arr(i) = dexp(this%loga_table(i))
    end do
    bisec_bracketed = .false.

    !Analytic first guess of the initial v1 using m/H_* = 3 asymptotics
    if (maxion_twiddle < 3.0d0) then !axions never oscillate (DE-like)
        v1_initguess(2) = dsqrt(omegah2_ax)/maxion_twiddle
    else if ((maxion_twiddle**2.0d0)/9.0d0 < &
        (omegah2_m**4.0d0)/(BG%omegah2_rad**3.0d0) + (BG%omegah2_rad**3.0d0)/(omegah2_m**4.0d0)) then
        !oscillation starts approximately in matter domination
        v1_initguess(2) = dsqrt(omegah2_ax)/(3.0d0*dsqrt(omegah2_m)/hnot)
    else !oscillation starts approximately in radiation domination
        v1_initguess(2) = dsqrt(omegah2_ax)/(((9.0d0*BG%omegah2_rad/hsq)**0.375d0)*(maxion_twiddle**0.25d0))
    end if
    v1_initguess(1) = v1_initguess(2)/2.0d0
    v1_initguess(3) = v1_initguess(2)*2.0d0
    v1_initguess(2) = (v1_initguess(1) + v1_initguess(3))/2.0d0
    !42 = "not computed yet" sentinel; any value > 1 triggers recomputation
    omaxh2_guess = (/42.0d0, 42.0d0, 42.0d0/)

    iter_c = 0
    allocate(v_vec(2,ntable))
    allocate(littlehfunc(ntable), diagnostic(ntable), f_arr(ntable), v_buff(ntable), abuff(ntable))
    do while (iter_c < nphi)
        do j = 1, 3
            if (omaxh2_guess(j) > 1.0d0) then !needs (re)computation
                aosc_guess(j) = 15.0d0
                vtwiddle_init = v1_initguess(j)

                !start the field at rest on the hill...
                v_vec(1,1) = vtwiddle_init
                v_vec(2,1) = 0.0d0
                call lh(omegah2_regm, BG%omegah2_rad, omegah2_lambda, omk, hsq, &
                    maxion_twiddle, a_arr(1), v_vec(1:2,1), littlehfunc(1), badflag, &
                    lhsqcont_massless, lhsqcont_massive, BG%Nu_mass_eigenstates, BG%nu_masses)
                !...then put v2 on the early-time attractor (RD slow roll, the 1/5 factor)
                v_vec(2,1) = -vtwiddle_init*(a_arr(1)**2.0d0)*(maxion_twiddle**2.0d0)*hnot/(5.0d0*littlehfunc(1))

                kvec = 0.0d0
                kfinal = 0.0d0
                call next_step(a_arr(1), v_vec(1:2,1), kvec(1:2,1:16), kfinal(1:2), avec(1:16), &
                    omegah2_regm, BG%omegah2_rad, omegah2_lambda, omk, hsq, &
                    maxion_twiddle, badflag, dloga, 16, cmat(1:16,1:16), lhsqcont_massless, &
                    lhsqcont_massive, BG%Nu_mass_eigenstates, BG%nu_masses)
                diagnostic(1) = dfac*littlehfunc(1)/(a_arr(1)*hnot)

                do i = 2, ntable
                    !8th-order Runge-Kutta update: only stages 1, 9-15 carry weight
                    v_vec(:,i) = v_vec(:,i-1) + (svec(1)*kvec(:,1) + svec(9)*kvec(:,9) + svec(10)*kvec(:,10) &
                        + svec(11)*kvec(:,11) + svec(12)*kvec(:,12) + svec(13)*kvec(:,13) &
                        + svec(14)*kvec(:,14) + svec(15)*kvec(:,15))
                    call lh(omegah2_regm, BG%omegah2_rad, omegah2_lambda, omk, hsq, &
                        maxion_twiddle, a_arr(i), v_vec(1:2,i), littlehfunc(i), badflag, &
                        lhsqcont_massless, lhsqcont_massive, BG%Nu_mass_eigenstates, BG%nu_masses)
                    kvec = 0.0d0
                    kfinal = 0.0d0
                    call next_step(a_arr(i), v_vec(1:2,i), kvec(1:2,1:16), kfinal(1:2), &
                        avec(1:16), omegah2_regm, BG%omegah2_rad, omegah2_lambda, omk, hsq, &
                        maxion_twiddle, badflag, dloga, 16, cmat(1:16,1:16), &
                        lhsqcont_massless, lhsqcont_massive, BG%Nu_mass_eigenstates, BG%nu_masses)
                    !diagnostic m/(dfac*H): first grid crossing locates the switch
                    diagnostic(i) = dfac*littlehfunc(i)/(a_arr(i)*hnot)
                    if (aosc_guess(j) == 15.0d0) then
                        if (maxion_twiddle < diagnostic(i-1)) then
                            if (maxion_twiddle >= diagnostic(i)) then
                                aosc_guess(j) = a_arr(i)
                            end if
                        end if
                    end if
                end do
                if (badflag /= 0) then
                    call GlobalError('TAxionModel: collapsing or NaN background history', &
                        error_unsupported_params)
                    return
                end if

                f_arr(1:ntable) = maxion_twiddle/(diagnostic(1:ntable))
                f_arr = dlog(f_arr)

                if (maxion_twiddle < 10._dl) then
                    !no switch for m/H0 < 10: use the present-day field density for shooting
                    omaxh2_guess(j) = (v_vec(2,ntable)/a_arr(ntable))**2.0d0 + (maxion_twiddle*v_vec(1,ntable))**2.0d0
                else
                    if (aosc_guess(j) /= 15._dl) then
                        !refine a_osc: natural-spline root of ln(m/(dfac*H)) = 0
                        d1 = SPLINE_DANGLE
                        d2 = SPLINE_DANGLE
                        call spline(f_arr(1:ntable), this%loga_table(1:ntable), ntable, d1, d2, abuff(1:ntable))
                        call spline_out(f_arr(1:ntable), this%loga_table(1:ntable), abuff(1:ntable), &
                            ntable, 0.0d0, laosc)
                    else
                        !switch would be after the present day
                        aosc_guess(j) = 1.0_dl - 1.e-3_dl
                        laosc = dlog(aosc_guess(j))
                    end if

                    !field value and derivative at the switch (analytic spline boundary
                    !derivatives from the exact KG equations; abscissa still ln a here)
                    d1 = v_vec(2,1)*hnot/littlehfunc(1)
                    d2 = v_vec(2,ntable)*hnot/littlehfunc(ntable)
                    call spline(this%loga_table(1:ntable), v_vec(1,1:ntable), ntable, d1, d2, v_buff(1:ntable))
                    call spline_out(this%loga_table(1:ntable), v_vec(1,1:ntable), v_buff(1:ntable), &
                        ntable, laosc, phiosc)
                    d1 = -(2.0d0*v_vec(2,1) + &
                        ((maxion_twiddle*dexp(this%loga_table(1)))**2.0d0)*v_vec(1,1)*hnot/littlehfunc(1))
                    d2 = -(2.0d0*v_vec(2,ntable) + &
                        ((maxion_twiddle*dexp(this%loga_table(ntable)))**2.0d0)*v_vec(1,ntable)*hnot/littlehfunc(ntable))
                    call spline(this%loga_table(1:ntable), v_vec(2,1:ntable), ntable, d1, d2, v_buff(1:ntable))
                    call spline_out(this%loga_table(1:ntable), v_vec(2,1:ntable), v_buff(1:ntable), &
                        ntable, laosc, phidosc)
                    aosc_guess(j) = dexp(laosc)

                    !EFA matching at the switch (iterates <H> and wEFA_c to consistency)
                    call this%auxiIC(omegah2_regm, BG%omegah2_rad, omegah2_lambda, omk, hnot, &
                        maxion_twiddle, aosc_guess(j), (/phiosc, phidosc/), badflag, lhsqcont_massless, &
                        lhsqcont_massive, BG%Nu_mass_eigenstates, BG%nu_masses, littlehauxi, &
                        this%ahosc_ETA, A_coeff, tvarphi_c, tvarphi_cp, tvarphi_s, tvarphi_sp, rhorefp, Prefp)

                    !present-day density: a^-3 dilution with the w = wEFA_c (H/m)^2 correction
                    wcorr_coeff = this%ahosc_ETA*aosc_guess(j)/(maxion_twiddle*hnot)
                    omaxh2_wcorr = rhorefp*(aosc_guess(j)**3.0d0)*dexp((wcorr_coeff**2.0d0)*3.0d0* &
                        this%wEFA_c*(1.0d0 - 1.0d0/(aosc_guess(j)**4.0d0))/4.0d0)
                    omaxh2_guess(j) = omaxh2_wcorr

                    !dfac that would relocate the switch to z=800 (for the recombination skip)
                    if (aosc_guess(j) < this%a_skip .and. aosc_guess(j) >= this%a_skipst) then
                        call lh(omegah2_regm, BG%omegah2_rad, omegah2_lambda, omk, hsq, &
                            maxion_twiddle, this%a_skip, v_vec(1:2,1), lh_skip, badflag, &
                            lhsqcont_massless, lhsqcont_massive, BG%Nu_mass_eigenstates, BG%nu_masses, &
                            rhorefp*((aosc_guess(j)/this%a_skip)**3.0d0)*dexp((wcorr_coeff**2.0d0)*3.0d0* &
                            this%wEFA_c*(1.0d0/(this%a_skip**4.0d0) - 1.0d0/(aosc_guess(j)**4.0d0))/4.0d0))
                        this%dfac_skip = min(littlehauxi/this%ahosc_ETA, 1._dl)*lh_skip
                        this%dfac_skip = (maxion_twiddle*this%a_skip*hnot/this%dfac_skip)
                    end if
                end if
            end if

            if (abs(omaxh2_guess(j)/omegah2_ax - 1.0d0) < 1.0d-6) then
                vtwiddle_init = v1_initguess(j)
                this%a_osc = aosc_guess(j)
                iter_c = -1
                exit
            end if
        end do

        if (iter_c == -1) exit
        if (iter_c == nphi - 1) then
            write(*,*) 'Warning: exceeding the maximum number of iterations for bisection: ', &
                'iterations:', iter_c, 'omaxh2 result:', omaxh2_guess(2), 'omaxh2 input:', &
                omegah2_ax, 'fractional error:', omaxh2_guess(2)/omegah2_ax - 1.0d0, &
                'aosc:', aosc_guess(2), '. The code will proceed - possibly wEF was not set to ', &
                'converge accurately enough which is ok, but please do check for unreasonable inputs.'
            vtwiddle_init = v1_initguess(2)
            this%a_osc = aosc_guess(2)
            if (this%a_osc >= 1.0d0) this%a_osc = 1.0d0
            exit
        end if

        if ((omaxh2_guess(1) - omegah2_ax)*(omaxh2_guess(3) - omegah2_ax) > 0 .and. .not. bisec_bracketed) then
            !still bracketing: expand the failing end (monotonicity assumed)
            if (omaxh2_guess(3) < omegah2_ax) then
                v1_initguess(3) = v1_initguess(3)*2.0d0
                omaxh2_guess(3) = 42.0d0
            else if (omaxh2_guess(1) > omegah2_ax) then
                v1_initguess(1) = v1_initguess(1)/2.0d0
                omaxh2_guess(1) = 42.0d0
            end if
        else if ((omaxh2_guess(1) - omegah2_ax)*(omaxh2_guess(3) - omegah2_ax) < 0 .and. .not. bisec_bracketed) then
            if ((omaxh2_guess(1) - omegah2_ax)*(omaxh2_guess(2) - omegah2_ax) < 0) then
                v1_initguess(3) = v1_initguess(2)
                omaxh2_guess(3) = omaxh2_guess(2)
            else if ((omaxh2_guess(3) - omegah2_ax)*(omaxh2_guess(2) - omegah2_ax) < 0) then
                v1_initguess(1) = v1_initguess(2)
                omaxh2_guess(1) = omaxh2_guess(2)
            end if
            bisec_bracketed = .true.
        else if ((omaxh2_guess(1) - omegah2_ax)*(omaxh2_guess(3) - omegah2_ax) < 0 .and. bisec_bracketed) then
            if ((omaxh2_guess(1) - omegah2_ax)*(omaxh2_guess(2) - omegah2_ax) < 0) then
                v1_initguess(3) = v1_initguess(2)
                omaxh2_guess(3) = omaxh2_guess(2)
                aosc_guess(3) = aosc_guess(2)
            else if ((omaxh2_guess(3) - omegah2_ax)*(omaxh2_guess(2) - omegah2_ax) < 0) then
                v1_initguess(1) = v1_initguess(2)
                omaxh2_guess(1) = omaxh2_guess(2)
                aosc_guess(1) = aosc_guess(2)
            end if
        else if ((omaxh2_guess(1) - omegah2_ax)*(omaxh2_guess(3) - omegah2_ax) > 0 .and. bisec_bracketed) then
            call GlobalError('TAxionModel: bisection already started, bracket lost', &
                error_unsupported_params)
            badflag = 1
            return
        end if

        v1_initguess(2) = (v1_initguess(1) + v1_initguess(3))/2.0_dl
        omaxh2_guess(2) = 42.0d0
        aosc_guess(2) = 15.0d0

        !once close (within ~10%), truncate the grid to a_final = 1.1 a_osc - in the
        !oscillating case the final tables do NOT reach a=1; consumers must use the
        !analytic EFA density beyond a_osc
        if (omaxh2_guess(3)/omegah2_ax < 1.1d0 .and. a_final == 1.0d0 .and. aosc_guess(3)*1.1d0 < 1.0d0) then
            a_final = aosc_guess(3)*1.1d0
            log_a_final = dlog(a_final)
            dloga = (log_a_final - log_a_init)/(dble(ntable - 1))
            this%loga_table(1) = log_a_init
            do i = 2, ntable
                this%loga_table(i) = log_a_init + dloga*dble(i-1)
                a_arr(i) = dexp(this%loga_table(i))
            end do
        end if
        iter_c = iter_c + 1
    end do
    deallocate(f_arr, v_buff, abuff)
    deallocate(cmat, kvec, kfinal, svec, avec)

    !Build the output tables from the converged history
    allocate(rhoaxh2_ov_rhom(ntable))
    allocate(this%phinorm_table(ntable), this%phidotnorm_table(ntable))
    allocate(this%phinorm_table_ddlga(ntable), this%phidotnorm_table_ddlga(ntable))
    do i = 1, ntable
        this%phinorm_table(i) = v_vec(1,i)
        this%phidotnorm_table(i) = v_vec(2,i)
        !instantaneous Omega_ax(a) h^2 = rho_ax(a) h^2/rho_crit,0
        rhoaxh2_ov_rhom(i) = (v_vec(2,i)/a_arr(i))**2.0d0 + (maxion_twiddle*v_vec(1,i))**2.0d0
    end do
    deallocate(diagnostic)

    !Change of base: loga_table becomes log10(a) from here on (note the log(10)
    !Jacobians in the analytic spline boundary derivatives below)
    this%loga_table = dlog10(dexp(this%loga_table))

    d1 = dlog(10._dl)*this%phidotnorm_table(1)*hnot/littlehfunc(1)
    d2 = dlog(10._dl)*this%phidotnorm_table(ntable)*hnot/littlehfunc(ntable)
    call spline(this%loga_table(1:ntable), this%phinorm_table, ntable, d1, d2, this%phinorm_table_ddlga)
    d1 = -dlog(10._dl)*(2._dl*this%phidotnorm_table(1) + &
        ((maxion_twiddle*(10._dl**(this%loga_table(1))))**2._dl)*this%phinorm_table(1)*hnot/littlehfunc(1))
    d2 = -dlog(10._dl)*(2._dl*this%phidotnorm_table(ntable) + &
        ((maxion_twiddle*(10._dl**(this%loga_table(ntable))))**2._dl)*this%phinorm_table(ntable)*hnot/littlehfunc(ntable))
    call spline(this%loga_table(1:ntable), this%phidotnorm_table, ntable, d1, d2, this%phidotnorm_table_ddlga)

    !Re-evaluate the field at the converged a_osc and do the final EFA matching
    call spline_out(this%loga_table, this%phinorm_table, this%phinorm_table_ddlga, ntable, &
        dlog10(this%a_osc), v1_ref)
    call spline_out(this%loga_table, this%phidotnorm_table, this%phidotnorm_table_ddlga, ntable, &
        dlog10(this%a_osc), v2_ref)
    call this%auxiIC(omegah2_regm, BG%omegah2_rad, omegah2_lambda, omk, hnot, &
        maxion_twiddle, this%a_osc, (/v1_ref, v2_ref/), badflag, lhsqcont_massless, &
        lhsqcont_massive, BG%Nu_mass_eigenstates, BG%nu_masses, littlehauxi, &
        this%ahosc_ETA, A_coeff, tvarphi_c, tvarphi_cp, tvarphi_s, tvarphi_sp, rhorefp, Prefp)

    this%ah_osc = littlehauxi
    this%A_coeff = A_coeff
    this%A_coeff_alt = A_coeff + 2.0d0*littlehauxi/(this%a_osc*hnot*maxion_twiddle)
    this%tvarphi_c = tvarphi_c
    this%tvarphi_s = tvarphi_s
    this%tvarphi_cp = tvarphi_cp
    this%tvarphi_sp = tvarphi_sp
    this%rhorefp_ovh2 = rhorefp/hsq
    this%Prefp = Prefp   !left in Omega(a) h^2 units (not divided by h^2)

    !Density log-table: purely the KG solution; EFA extrapolation beyond a_osc is
    !always done analytically by RhoaxH2AtA, never from this table
    allocate(this%rhoaxh2ovrhom_logtable(ntable), this%rhoaxh2ovrhom_logtable_buff(ntable))
    this%rhoaxh2ovrhom_logtable = dlog10(rhoaxh2_ov_rhom)
    !dlog10(rho)/dlog10(a) = -6 KE/rho, exact for KG
    d1 = -6.0_dl*(this%phidotnorm_table(1)/(10._dl**(this%loga_table(1))))**2._dl/rhoaxh2_ov_rhom(1)
    d2 = -6.0_dl*(this%phidotnorm_table(ntable)/(10._dl**(this%loga_table(ntable))))**2._dl/rhoaxh2_ov_rhom(ntable)
    call spline(this%loga_table, this%rhoaxh2ovrhom_logtable, ntable, d1, d2, this%rhoaxh2ovrhom_logtable_buff)

    !Matter-radiation equality with the exact axion density counted as matter
    !(massive neutrinos counted as radiation); crude finite-difference BCs kept from
    !AxiECAMB for output parity
    allocate(eq_arr(ntable), eq_arr_buff(ntable))
    do i = 1, ntable
        eq_arr(i) = ((omegah2_regm/(a_arr(i)**3.0d0) + rhoaxh2_ov_rhom(i))/(BG%omegah2_rad/(a_arr(i)**4.0d0)))
    end do
    eq_arr = dlog(eq_arr)
    d1 = (this%loga_table(2) - this%loga_table(1))/(eq_arr(2) - eq_arr(1))
    d2 = (this%loga_table(ntable) - this%loga_table(ntable-1))/(eq_arr(ntable) - eq_arr(ntable-1))
    call spline(eq_arr(1:ntable), this%loga_table(1:ntable), ntable, d1, d2, eq_arr_buff(1:ntable))
    call spline_out(eq_arr(1:ntable), this%loga_table(1:ntable), eq_arr_buff(1:ntable), &
        ntable, 0.0d0, this%aeq)
    this%aeq = 10._dl**(this%aeq)
    !fallback if the spline breaks when a_osc < equality
    regzeq = (BG%omegah2_rad + sum(lhsqcont_massive(1:max(BG%Nu_mass_eigenstates,1))))/ &
        (omegah2_b + omegah2_dm + omegah2_ax)
    !pure-LCDM analytic equality assuming axions scale as matter (for phase targeting)
    this%aeq_LCDM = BG%omegah2_gamma*(1._dl + (BG%nu_massless_degeneracy + &
        BG%sum_nu_mass_degeneracies)*(7._dl/8._dl)*((4._dl/11._dl)**(4._dl/3._dl)))/ &
        (omegah2_b + omegah2_dm + omegah2_ax)
    if (this%a_osc < regzeq) this%aeq = regzeq
    deallocate(eq_arr, eq_arr_buff)

    !initial field value in reduced Planck units
    this%phiinit = vtwiddle_init*sqrt(6.0d0)
    if (this%axion_isocurvature) then
        block
            real(dl) Hinf_ovMpl  !Hinf is stored as log10(GeV); convert at use
            Hinf_ovMpl = (10**this%Hinf)/mplanck_GeV
            this%amp_i = Hinf_ovMpl**2/(const_pi**2*this%phiinit**2)
            this%r_val = 2*(Hinf_ovMpl**2/(const_pi**2.*BG%ScalarPowerAmp))
            this%alpha_ax = this%amp_i/BG%ScalarPowerAmp
        end block
    end if

    deallocate(a_arr, v_vec, littlehfunc, rhoaxh2_ov_rhom)
    end subroutine TAxionModel_w_evolve

    !> EFA matching at the switch: iteratively determines <H> and wEFA_c, projects
    !> (phi, phi') onto the cos/sin WKB basis and returns cycle-averaged rho and P.
    !> Mutates this%wEFA_c (stateful across calls, as in AxiECAMB).
    subroutine TAxionModel_auxiIC(this, omegah2_regm, omegah2_rad, omegah2_lambda, omk, hnot, &
        maxion_twiddle, a, v, badflag, lhsqcont_massless, lhsqcont_massive, &
        Nu_mass_eigenstates, Numasses, littlehauxi, lhETA, A_coeff, &
        tvarphi_c, tvarphi_cp, tvarphi_s, tvarphi_sp, rhorefp, Prefp)
    class(TAxionModel), intent(inout) :: this
    integer badflag, Nu_mass_eigenstates, i, iter_EFA, maxiter
    real(dl) omegah2_regm, omegah2_rad, omegah2_lambda, omk, hnot, maxion_twiddle, a, a2
    real(dl) v(1:2), littlehauxi, lhETA, lhETA_upd, lhsqcont_massless
    real(dl) lhsqcont_massive(max(Nu_mass_eigenstates,1))
    real(dl) mass_correctors(max(Nu_mass_eigenstates,1)), w_nu(max(Nu_mass_eigenstates,1)), rhonu, pnu
    real(dl) numasses(max(Nu_mass_eigenstates,1))
    real(dl) w_ax, wEFA_c_upd
    real(dl) Hovm_ins, Hovm_ETA, dHsqdmt_term, A_coeff, A_denom
    real(dl) tvarphi_c, tvarphi_cp, tvarphi_s, tvarphi_sp, rhorefp, Prefp
    real(dl) tol_EFA

    !instantaneous conformal aH at the switch (full Friedmann with field energy)
    call lh(omegah2_regm, omegah2_rad, omegah2_lambda, omk, hnot**2.0d0, maxion_twiddle, a, v, &
        littlehauxi, badflag, lhsqcont_massless, lhsqcont_massive, Nu_mass_eigenstates, Numasses)
    lhETA = littlehauxi
    mass_correctors = 0._dl
    w_nu = 0._dl
    do i = 1, Nu_mass_eigenstates
        call ThermalNuBack%rho_P(a*numasses(i), rhonu, pnu)
        mass_correctors(i) = rhonu
        w_nu(i) = pnu/rhonu
    end do
    a2 = a**2.0d0
    !start from the instantaneous energy density ratio
    rhorefp = (v(2)/a)**2.0d0 + (maxion_twiddle*v(1))**2.0d0
    tol_EFA = 1.e-7_dl
    maxiter = 30
    Hovm_ins = littlehauxi/(a*hnot*maxion_twiddle)
    do iter_EFA = 1, maxiter
        Hovm_ETA = lhETA/(a*hnot*maxion_twiddle)
        w_ax = (Hovm_ETA**2.0d0)*this%wEFA_c
        !dH^2/d(mt) assembled from the continuity equations of all components
        !(dimensional bookkeeping: a^2 factors per species cancel against (H0/H)^2)
        dHsqdmt_term = omegah2_regm*(-3.0d0)/a + &
            omegah2_rad*(-4.0d0)/(a2) + &
            sum(lhsqcont_massive(1:max(Nu_mass_eigenstates,1))*mass_correctors(1:max(Nu_mass_eigenstates,1)) &
            *(-3.0d0*(1.0d0 + w_nu(1:max(Nu_mass_eigenstates,1)))))/(a2) + &
            rhorefp*(-3.0d0*(w_ax + 1.0d0))*(a2) + &
            omk*(hnot**2.0d0)*(-2.0d0)
        dHsqdmt_term = dHsqdmt_term/(lhETA**2.0d0)
        A_coeff = (-Hovm_ETA/2.0d0)*(3.0d0 - dHsqdmt_term)
        A_denom = A_coeff**2.0d0 + 3.0d0*A_coeff*Hovm_ins + 4.0d0
        !the four WKB projections (normalization sqrt(4 pi G/3) h)
        tvarphi_c = v(1)
        tvarphi_cp = -3.0d0*Hovm_ins*(2.0d0*v(1) + (A_coeff + 3.0d0*Hovm_ins)*v(2)/(a*maxion_twiddle))/A_denom
        tvarphi_s = v(2)/(a*maxion_twiddle) - tvarphi_cp
        tvarphi_sp = 3.0d0*Hovm_ins*(A_coeff*v(1) - 2.0d0*v(2)/(a*maxion_twiddle))/A_denom
        !cycle-averaged EFA density and pressure
        rhorefp = (maxion_twiddle**2.0d0)*(tvarphi_c**2.0d0 + tvarphi_s**2.0d0 + &
            (tvarphi_cp**2.0d0 + tvarphi_sp**2.0d0)/2.0d0 - tvarphi_c*tvarphi_sp + tvarphi_s*tvarphi_cp)
        Prefp = (maxion_twiddle**2.0d0)*(tvarphi_cp**2.0d0/2.0d0 + tvarphi_sp**2.0d0/2.0d0 - &
            tvarphi_c*tvarphi_sp + tvarphi_s*tvarphi_cp)

        wEFA_c_upd = (Prefp/rhorefp)/((lhETA/(maxion_twiddle*hnot*a))**2._dl)
        call lh(omegah2_regm, omegah2_rad, omegah2_lambda, omk, hnot**2.0d0, maxion_twiddle, &
            a, v, lhETA_upd, badflag, lhsqcont_massless, lhsqcont_massive, Nu_mass_eigenstates, &
            Numasses, rhorefp)
        if (abs(wEFA_c_upd/this%wEFA_c - 1.0_dl) < tol_EFA) then
            exit
        else
            this%wEFA_c = wEFA_c_upd
            lhETA = lhETA_upd
        end if
    end do
    if (iter_EFA == maxiter) then
        print '(A, I0, A, E13.5E3, A, E13.5E3)', &
            'Warning: exceeding the maximum number of iterations (', iter_EFA, &
            ') for the EFA matching: fractional error in w = ', &
            wEFA_c_upd/this%wEFA_c - 1.0_dl, ', fractional error in H_ETA = ', lhETA_upd/lhETA - 1.0_dl
    end if
    end subroutine TAxionModel_auxiIC

    !> RHS of the KG system in d/d(ln a)
    subroutine derivs_bg(a, v, dvt_dloga, omegah2_regm, omegah2_rad, omegah2_lambda, omk, hsq, &
        maxion_twiddle, badflag, lhsqcont_massless, lhsqcont_massive, Nu_mass_eigenstates, Nu_masses)
    integer badflag, Nu_mass_eigenstates
    real(dl) a
    real(dl) dvt_dloga(1:2), dvt_da(1:2), lhr
    real(dl) v(1:2)
    real(dl) omegah2_regm, omegah2_rad, omegah2_lambda
    real(dl) maxion_twiddle, Nu_masses(max(Nu_mass_eigenstates,1))
    real(dl) omk, hsq, lhsqcont_massless, lhsqcont_massive(max(Nu_mass_eigenstates,1))

    call lh(omegah2_regm, omegah2_rad, omegah2_lambda, omk, hsq, maxion_twiddle, a, v, lhr, badflag, &
        lhsqcont_massless, lhsqcont_massive, Nu_mass_eigenstates, Nu_masses)
    !dsqrt(hsq) = hnot converts lhr (100 km/s/Mpc units) to aH/H0
    dvt_da(1) = v(2)*dsqrt(hsq)/(a*(lhr))
    dvt_da(2) = -2.0d0*v(2)/(a) - (maxion_twiddle**2.0d0)*a*v(1)*dsqrt(hsq)/(lhr)
    dvt_dloga(1:2) = a*dvt_da(1:2)
    end subroutine derivs_bg

    !> Conformal aH in units of 100 km/s/Mpc (despite the "littleh" heritage name).
    !> If rho_f is present it replaces the field energy with the given EFA density.
    subroutine lh(omegah2_regm, omegah2_rad, omegah2_lambda, omk, hsq, maxion_twiddle, a, v, &
        littlehfunc, badflag, lhsqcont_massless, lhsqcont_massive, Nu_mass_eigenstates, nu_masses, rho_f)
    integer badflag, Nu_mass_eigenstates, i
    real(dl) omegah2_regm, omegah2_rad, omegah2_lambda, omk, hsq, maxion_twiddle, a
    real(dl) v(1:2), littlehfunc, lhsqcont_massless, lhsqcont_massive(max(Nu_mass_eigenstates,1))
    real(dl) mass_correctors(max(Nu_mass_eigenstates,1)), rhonu
    real(dl) nu_masses(max(Nu_mass_eigenstates,1))
    real(dl), optional :: rho_f

    mass_correctors = 0._dl
    do i = 1, Nu_mass_eigenstates
        call ThermalNuBack%rho(a*nu_masses(i), rhonu)
        mass_correctors(i) = rhonu
    end do
    !lhsqcont_massless is already folded into omegah2_rad (kept for interface parity)
    littlehfunc = (omegah2_regm/(a**3.0d0) + omegah2_rad/(a**4.0d0)) + &
        sum(lhsqcont_massive(1:max(Nu_mass_eigenstates,1))*mass_correctors(1:max(Nu_mass_eigenstates,1)))/(a**4.0d0)
    littlehfunc = littlehfunc + omegah2_lambda
    if (present(rho_f)) then
        littlehfunc = littlehfunc + rho_f
    else
        littlehfunc = littlehfunc + (maxion_twiddle*v(1))**2.0d0 + ((v(2)/a)**2.0d0)
    end if
    littlehfunc = littlehfunc*(a**2.0d0) + omk*hsq  !Hubble -> conformal Hubble
    !flag collapsing universes and NaN histories
    if (littlehfunc <= 0.0d0) then
        badflag = 1
    end if
    if (.not. littlehfunc == littlehfunc) then
        badflag = 1
    end if
    littlehfunc = dsqrt(littlehfunc)
    end subroutine lh

    !> One 16-stage Runge-Kutta evaluation: fills the stage derivatives kvec
    !> (each pre-multiplied by dloga)
    subroutine next_step(a, v, kvec, kfinal, avec, omegah2_regm, omegah2_rad, omegah2_lambda, omk, &
        hsq, maxion_twiddle, badflag, dloga, nstep, cmat, &
        lhsqcont_massless, lhsqcont_massive, Nu_mass_eigenstates, Nu_masses)
    integer nstep, cp, m, badflag
    real(dl) hsq, a, v(1:2), kvec(1:2,1:nstep), omegah2_regm, omegah2_rad, omegah2_lambda
    real(dl) maxion_twiddle, dloga, omk
    real(dl) vfeed(1:2), cmat(1:nstep,1:nstep), kfinal(1:2), avec(1:nstep)
    integer Nu_mass_eigenstates
    real(dl) lhsqcont_massless, lhsqcont_massive(max(Nu_mass_eigenstates,1))
    real(dl) Nu_masses(max(Nu_mass_eigenstates,1))

    kvec = 0.0d0
    kfinal = 0.0d0
    call derivs_bg(a, v(1:2), kvec(1:2,1), omegah2_regm, omegah2_rad, omegah2_lambda, omk, hsq, &
        maxion_twiddle, badflag, lhsqcont_massless, lhsqcont_massive, Nu_mass_eigenstates, Nu_masses)
    kvec(1:2,1) = kvec(1:2,1)*dloga
    do m = 1, nstep
        do cp = 1, 2
            vfeed(cp) = dot_product(cmat(m, 1:m), kvec(cp, 1:m))
        end do
        if (m <= (nstep-1)) then
            call derivs_bg(a*dexp(dloga*avec(m)), v(1:2) + vfeed(1:2), &
                kvec(1:2,m+1), omegah2_regm, omegah2_rad, omegah2_lambda, omk, hsq, &
                maxion_twiddle, badflag, lhsqcont_massless, lhsqcont_massive, Nu_mass_eigenstates, Nu_masses)
            kvec(1:2,m+1) = kvec(1:2,m+1)*dloga
        else
            call derivs_bg(a*dexp(dloga*avec(m)), v(1:2) + vfeed(1:2), &
                kfinal(1:2), omegah2_regm, omegah2_rad, omegah2_lambda, omk, hsq, &
                maxion_twiddle, badflag, lhsqcont_massless, lhsqcont_massive, Nu_mass_eigenstates, Nu_masses)
            kfinal(1:2) = kfinal(1:2)*dloga
        end if
    end do
    end subroutine next_step

    !> Cubic-spline evaluation with bisection lookup (NR splint; double precision
    !> implementation from AxiECAMB subroutines.f90, originally by Dan Grin)
    subroutine spline_out(xarr, yarr, yarr_buff, n, x, y)
    integer n, llo_out, lhi_out, midp
    real(dl) xarr(n), yarr(n), yarr_buff(n)
    real(dl) x, y, a0_out, b0_out, ho_out

    llo_out = 1
    lhi_out = n
    do while ((lhi_out - llo_out) > 1)
        midp = (llo_out + lhi_out)/2
        if (xarr(midp) > x) then
            lhi_out = midp
        else
            llo_out = midp
        end if
    end do
    ho_out = xarr(lhi_out) - xarr(llo_out)
    a0_out = (xarr(lhi_out) - x)/ho_out
    b0_out = (x - xarr(llo_out))/ho_out
    y = a0_out*yarr(llo_out) + b0_out*yarr(lhi_out) + ((a0_out**3.0d0 - a0_out)* &
        yarr_buff(llo_out) + (b0_out**3.0d0 - b0_out)* &
        yarr_buff(lhi_out))*ho_out*ho_out/6.d0
    end subroutine spline_out

    end module AxionBackground
