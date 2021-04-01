""""Panels for DDB files."""

import param
import panel as pn
import panel.widgets as pnw
import bokeh.models.widgets as bkw

from abipy.panels.core import AbipyParameterized


class DdbFilePanel(AbipyParameterized):
    """
    A panel to analyze a |DdbFile|.
    Provides widgets to invoke anaddb and visualize the results.
    """
    verbose = param.Integer(0, bounds=(0, None), doc="Verbosity Level")
    mpi_procs = param.Integer(1, bounds=(1, None), doc="Number of MPI processes used in anaddb")

    nqsmall = param.Integer(10, bounds=(1, None), doc="Number of divisions for smallest vector to generate Q-mesh")
    ndivsm = param.Integer(5, bounds=(1, None), doc="Number of divisions for smallest segment in q-path")
    lo_to_splitting = param.ObjectSelector(default="automatic", objects=["automatic", True, False])
    chneut = param.ObjectSelector(default=1, objects=[0, 1, 2], doc="Abinit variable")
    dipdip = param.ObjectSelector(default=1, objects=[0, 1], doc="Abinit variable")
    asr = param.ObjectSelector(default=2, objects=[0, 1, 2], doc="Abinit variable")
    units = param.ObjectSelector(default="eV", objects=["eV", "meV", "Ha", "cm-1", "Thz"], doc="Energy units")

    dos_method = param.ObjectSelector(default="tetra", objects=["tetra", "gaussian"], doc="Integration method for DOS")
    temp_range = pnw.RangeSlider(name="T-range", start=0.0, end=1000, value=(0.0, 300.0), step=20)

    gamma_ev = param.Number(1e-4, bounds=(1e-20, None), doc="Phonon linewidth in eV")
    w_range = pnw.RangeSlider(name="Frequency range (eV)", start=0.0, end=1.0,
                              value=(0.0, 0.1), step=0.001)

    get_epsinf_btn = pnw.Button(name="Compute", button_type='primary')
    plot_phbands_btn = pnw.Button(name="Plot Bands and DOS", button_type='primary')
    plot_eps0w_btn = pnw.Button(name="Plot eps0(omega)", button_type='primary')

    plot_vsound_btn = pnw.Button(name="Calculate speed of sound", button_type='primary')
    plot_check_asr_dipdip_btn = pnw.Button(name="Compute phonons with/wo ASR and DIPDIP", button_type='primary')
    plot_ifc_btn = pnw.Button(name="Compute IFC(R)", button_type='primary')

    def __init__(self, ddb, **params):
        super().__init__(**params)
        self.ddb = ddb

    @param.depends('get_epsinf_btn.clicks')
    def get_epsinf(self):
        """Compute eps_infinity and Born effective charges from DDB."""
        if self.get_epsinf_btn.clicks == 0: return

        epsinf, becs = self.ddb.anaget_epsinf_and_becs(chneut=self.chneut,
                                                       mpi_procs=self.mpi_procs, verbose=self.verbose)

        gen, inp = self.ddb.anaget_dielectric_tensor_generator(asr=self.asr, chneut=self.chneut, dipdip=self.dipdip,
                                                               mpi_procs=self.mpi_procs, verbose=self.verbose,
                                                               return_input=True)

        eps0 = gen.tensor_at_frequency(w=0, gamma_ev=self.gamma_ev)
        #eps0 = pnw.DataFrame(eps0.get_dataframe())
        return pn.Column(epsinf, eps0, becs, pn.pane.HTML(inp._repr_html_()))

    @param.depends('plot_eps0w_btn.clicks')
    def plot_eps0w(self):
        """Compute eps0(omega) from DDB and plot the results."""
        if self.plot_eps0w_btn.clicks == 0: return
        gen, inp = self.ddb.anaget_dielectric_tensor_generator(asr=self.asr, chneut=self.chneut, dipdip=self.dipdip,
                                                               mpi_procs=self.mpi_procs, verbose=self.verbose,
                                                               return_input=True)
        ws = self.w_range.value
        w_max = ws[1]
        if w_max == 1.0: w_max = None # Will compute w_max in plot routine from ph freqs.

        def p(component, reim):
            return gen.plot(w_min=ws[0], w_max=w_max, gamma_ev=self.gamma_ev, num=500, component=component,
                            reim=reim, units=self.units, **self.fig_kwargs)

        # Build grid
        gspec = pn.GridSpec(sizing_mode='scale_width')
        gspec[0, 0] = p("diag", "re")
        gspec[0, 1] = p("diag", "im")
        gspec[1, 0] = p("offdiag", "re")
        gspec[1, 1] = p("offdiag", "im")
        gspec[2, :] = gen.get_oscillator_dataframe(reim="all", tol=1e-6)
        # Add HTML pane with input.
        gspec[3, 0] = pn.pane.HTML(inp._repr_html_())

        return gspec

    @param.depends('plot_phbands_btn.clicks')
    def plot_phbands_and_phdos(self, event=None):
        """Compute phonon bands and ph-DOSes from DDB and plot the results."""
        if self.plot_phbands_btn.clicks == 0: return
        #self.plot_phbands_btn.button_type = "warning"

        #print("Computing phbands")
        with self.ddb.anaget_phbst_and_phdos_files(
                nqsmall=self.nqsmall, qppa=None, ndivsm=self.ndivsm,
                line_density=None, asr=self.asr, chneut=self.chneut, dipdip=self.dipdip,
                dos_method=self.dos_method, lo_to_splitting=self.lo_to_splitting,
                verbose=self.verbose, mpi_procs=self.mpi_procs, return_input=True) as g:

            phbst_file, phdos_file = g
            phbands, phdos = phbst_file.phbands, phdos_file.phdos
            #print("Computing phbands completed")

            # Build grid
            gspec = pn.GridSpec(sizing_mode='scale_width')
            #gspec[0, 0] = phbands.plot_with_phdos(phdos, units=self.units, **self.fig_kwargs)
            gspec[0, 0] = phbands.plotly_with_phdos(phdos, units=self.units, show=False)
            gspec[0, 1] = phdos_file.plot_pjdos_type(units=self.units, exchange_xy=True, **self.fig_kwargs)
            gspec[1, 0] = phdos_file.msqd_dos.plot(units=self.units, **self.fig_kwargs)
            temps = self.temp_range.value
            #gspec[1, 1] = phdos.plot_harmonic_thermo(tstart=temps[0], tstop=temps[1], num=50, **self.fig_kwargs)
            gspec[1, 1] = phdos.plotly_harmonic_thermo(tstart=temps[0], tstop=temps[1], num=50, show=False)
            #msqd_dos.plot_tensor(**self.fig_kwargs)
            #self.plot_phbands_btn.button_type = "primary"

            # Add HTML pane with input
            gspec[2,:] = pn.pane.HTML(g.input._repr_html_())

        return gspec

    @param.depends('plot_vsound_btn.clicks')
    def plot_vsound(self):
        """
        Compute the speed of sound by fitting phonon frequencies
        along selected directions by linear least-squares fit.
        """
        if self.plot_vsound_btn.clicks == 0: return
        from abipy.dfpt.vsound import SoundVelocity
        sv = SoundVelocity.from_ddb(self.ddb.filepath, num_points=20, qpt_norm=0.1,
                                    ignore_neg_freqs=True, asr=self.asr, chneut=self.chneut, dipdip=self.dipdip,
                                    verbose=self.verbose, mpi_procs=self.mpi_procs)

        # Insert results in grid.
        gspec = pn.GridSpec(sizing_mode='scale_width')
        gspec[0, :1] = sv.get_dataframe()
        gspec[1, :1] = sv.plot(**self.fig_kwargs)

        return gspec

    @param.depends('plot_check_asr_dipdip_btn.clicks')
    def plot_without_asr_dipdip(self):
        if self.plot_check_asr_dipdip_btn.clicks == 0: return

        # Insert results in grid.
        gspec = pn.GridSpec(sizing_mode='scale_width')

        asr_plotter = self.ddb.anacompare_asr(asr_list=(0, 2), chneut_list=(1, ), dipdip=1,
                                              lo_to_splitting=self.lo_to_splitting,
                                              nqsmall=self.nqsmall, ndivsm=self.ndivsm,
                                              dos_method=self.dos_method, ngqpt=None,
                                              verbose=self.verbose, mpi_procs=self.mpi_procs)
        gspec[0, :1] = asr_plotter.plot(**self.fig_kwargs)

        dipdip_plotter = self.ddb.anacompare_dipdip(chneut_list=(1,), asr=2, lo_to_splitting=self.lo_to_splitting,
                                                    nqsmall=self.nqsmall, ndivsm=self.ndivsm,
                                                    dos_method=self.dos_method, ngqpt=None,
                                                    verbose=self.verbose, mpi_procs=self.mpi_procs)
        gspec[1, :1] = dipdip_plotter.plot(**self.fig_kwargs)
        return gspec

    @param.depends('plot_ifc_btn.clicks')
    def plot_ifc(self):
        if self.plot_ifc_btn.clicks == 0: return

        ifc = self.ddb.anaget_ifc(asr=self.asr, chneut=self.chneut, dipdip=self.dipdip)

        # Insert results in grid.
        gspec = pn.GridSpec(sizing_mode='scale_width')
        gspec[0, :1] = ifc.plot_longitudinal_ifc(title="Longitudinal IFCs", show=False)
        gspec[1, :1] = ifc.plot_longitudinal_ifc_short_range(title="Longitudinal IFCs short range", show=False)
        gspec[2, :1] = ifc.plot_longitudinal_ifc_ewald(title="Longitudinal IFCs Ewald", show=False)

        return gspec

    def get_panel(self):
        """Return tabs with widgets to interact with the DDB file."""
        tabs = pn.Tabs(); app = tabs.append
        row = pn.Row(bkw.PreText(text=self.ddb.to_string(verbose=self.verbose), sizing_mode="scale_both"))
        app(("Summary", row))
        app(("Ph-bands", pn.Row(
            pn.Column("# PH-bands options",
                      *[self.param[k] for k in ("nqsmall", "ndivsm", "asr", "chneut", "dipdip",
                                                "lo_to_splitting", "dos_method")],
                      self.temp_range, self.plot_phbands_btn),
            self.plot_phbands_and_phdos)
        ))
        app(("BECs", pn.Row(
            pn.Column("# Born effective charges options",
                     *[self.param[k] for k in ("asr", "chneut", "dipdip", "gamma_ev")], self.get_epsinf_btn),
            self.get_epsinf)
        ))
        app(("eps0", pn.Row(
            pn.Column("# epsilon_0",
                      *[self.param[k] for k in ("asr", "chneut", "dipdip", "gamma_ev")],
                      self.w_range, self.plot_eps0w_btn),
            self.plot_eps0w)
        ))
        app(("Speed of sound", pn.Row(
            pn.Column("# Speed of sound options",
                      *[self.param[k] for k in ("asr", "chneut", "dipdip")],
                      self.plot_vsound_btn),
            self.plot_vsound)
        ))
        app(("Check ASR and DIPDIP", pn.Row(
            pn.Column("# Options",
                      *[self.param[k] for k in ("nqsmall", "ndivsm", "dos_method")],
                      self.plot_check_asr_dipdip_btn),
            self.plot_without_asr_dipdip)
        ))

        app(("IFCs", pn.Row(
            pn.Column("# Options",
                      *[self.param[k] for k in ("asr", "dipdip", "chneut")],
                      self.plot_ifc_btn),
            self.plot_ifc)
        ))

        app(("Global Opts", pn.Column("# Global Options", *[self.param[k] for k in ("units", "mpi_procs", "verbose")])))

        return tabs
