from ctypes import c_bool, c_double, c_int

from .baseconfig import AllocatableArrayDouble, F2003Class, fortran_class


@fortran_class
class AxionModel(F2003Class):
    """
    AxiECAMB ultralight-axion (ULA) component (arXiv:2412.15192).

    The exact Klein-Gordon equation is solved for the axion background and
    perturbations until m = dfac*H (dfac ~ 10, retuned internally), after which the
    axion is treated as an effective fluid tracking the cycle average of the fast
    field oscillations. For m/H0 < 10 the axion is dark-energy-like and the field
    equations are used to the present day.

    Set the input fields (m_ax and either omaxh2, or use_axfrac with omdah2/axfrac)
    before calling calc_background/get_results; everything else is derived output.
    """

    _fields_ = [
        ("active", c_bool, "axion component switched on"),
        ("m_ax", c_double, "axion mass in eV"),
        ("omaxh2", c_double, "Omega_ax h^2 (used when use_axfrac=False)"),
        ("use_axfrac", c_bool, "use (omdah2, axfrac) instead of omaxh2"),
        ("omdah2", c_double, "total dark matter Omega h^2 (when use_axfrac and m/H0>=10)"),
        ("axfrac", c_double, "axion fraction of DM (m/H0>=10) or of DE (m/H0<10)"),
        ("axion_isocurvature", c_bool, "accepted but force-disabled in this release"),
        ("Hinf", c_double, "log10 of H_inflation in GeV (only used by disabled isocurvature)"),
        ("dfac", c_double, "input KG->EFA switch threshold m = dfac*H (default 10)"),
        ("H0_eV", c_double, "H0 in eV (derived)"),
        ("m_ovH0", c_double, "m/H0 (derived)"),
        ("is_de_like", c_bool, "m/H0 < 10: KG solved to a=1, no fluid switch"),
        ("has_switch", c_bool, "a KG->EFA switch exists"),
        ("dfac_used", c_double, "dfac after internal retuning"),
        ("a_osc", c_double, "scale factor of the KG->EFA switch"),
        ("tau_osc", c_double, "conformal time of the switch"),
        ("aeq", c_double, "matter-radiation equality incl. exact axion density"),
        ("aeq_LCDM", c_double, "analytic equality assuming axions scale as matter"),
        ("phiinit", c_double, "initial field value in reduced Planck units"),
        ("ah_osc", c_double, "instantaneous conformal aH at a_osc (100 km/s/Mpc units)"),
        ("ahosc_ETA", c_double, "cycle-averaged <aH> at a_osc"),
        ("A_coeff", c_double, "EFA matching coefficient"),
        ("A_coeff_alt", c_double, "A_coeff + 2 aH/(a m)"),
        ("tvarphi_c", c_double, "WKB cos projection of the field at the switch"),
        ("tvarphi_s", c_double, "WKB sin projection"),
        ("tvarphi_cp", c_double, "WKB cos-derivative projection"),
        ("tvarphi_sp", c_double, "WKB sin-derivative projection"),
        ("wEFA_c", c_double, "coefficient of w_ax = wEFA_c (H/m)^2"),
        ("rhorefp_ovh2", c_double, "EFA cycle-averaged Omega_ax(a_osc), no h^2"),
        ("Prefp", c_double, "EFA cycle-averaged pressure at a_osc, Omega(a) h^2 units"),
        ("wcorr_coeff", c_double, "<aH>_osc a_osc/(m/H0 h)"),
        ("dfac_skip", c_double, "dfac that relocates the switch past z=800"),
        ("a_skip", c_double, "recombination-skip window upper edge (1/801)"),
        ("a_skipst", c_double, "recombination-skip window lower edge (1/1301)"),
        ("opac_tauosc", c_double, "opacity at tau_osc"),
        ("expmmu_tauosc", c_double, "exp(-kappa) at tau_osc"),
        ("amp_i", c_double, "isocurvature amplitude (inactive)"),
        ("r_val", c_double, "isocurvature tensor ratio bookkeeping (inactive)"),
        ("alpha_ax", c_double, "isocurvature fraction (inactive)"),
        ("omaxh2_eff", c_double, "Omega_ax h^2 actually used after axfrac logic"),
        ("omegah2_rad", c_double, "radiation Omega h^2 used in the solve"),
        ("hsq", c_double, "(H0/100)^2 cached from the solve"),
        ("ntable", c_int, "background table size"),
        ("a_table_min", c_double, "first scale factor in the tables"),
        ("loga_table", AllocatableArrayDouble),
        ("phinorm_table", AllocatableArrayDouble),
        ("phidotnorm_table", AllocatableArrayDouble),
        ("phinorm_table_ddlga", AllocatableArrayDouble),
        ("phidotnorm_table_ddlga", AllocatableArrayDouble),
        ("rhoaxh2ovrhom_logtable", AllocatableArrayDouble),
        ("rhoaxh2ovrhom_logtable_buff", AllocatableArrayDouble),
    ]

    _fortran_class_module_ = "AxionBackground"
    _fortran_class_name_ = "TAxionModel"

    def set_params(self, m_ax, omaxh2=None, omdah2=None, axfrac=None, dfac=10.0):
        """
        Configure the ultralight-axion component.

        :param m_ax: axion mass in eV (if negative, interpreted as log10(m_ax/eV))
        :param omaxh2: Omega_ax h^2 (mutually exclusive with omdah2/axfrac)
        :param omdah2: total dark matter density Omega h^2 (with axfrac)
        :param axfrac: axion fraction of the dark matter (m/H0>=10) or of the
            dark energy (m/H0<10)
        :param dfac: KG->EFA switch threshold m = dfac*H (default 10; retuned internally)
        """
        if m_ax < 0:
            m_ax = 10**m_ax
        self.m_ax = m_ax
        self.dfac = dfac
        if omaxh2 is not None:
            if omdah2 is not None or axfrac is not None:
                raise ValueError("set either omaxh2 or (omdah2, axfrac), not both")
            self.use_axfrac = False
            self.omaxh2 = omaxh2
            self.active = m_ax > 0 and omaxh2 > 0
        elif omdah2 is not None and axfrac is not None:
            self.use_axfrac = True
            self.omdah2 = omdah2
            self.axfrac = axfrac
            self.active = m_ax > 0 and axfrac > 0
        else:
            raise ValueError("give either omaxh2 or both omdah2 and axfrac")
        return self
