""""Panels for DDB files."""

import sys
import param
import panel as pn
import panel.widgets as pnw
import bokeh.models.widgets as bkw

from abipy.panels.core import AbipyParameterized, HasStructureParams, ButtonContext, mpl, ply, dfc
from abipy.dfpt.ddb import PhononBandsPlotter


class HasAnaddbParams(AbipyParameterized):
    """
    Base class for panel classes requiring widgets to invoke Anaddb via AbiPy.
    Used, for instance, by DdbFilePanel and DdbRobotPanel so that we don't have to
    repeat the same widget stuff over and over agains.
    """
    verbose = param.Integer(0, bounds=(0, None), doc="Verbosity level")
    mpi_procs = param.Integer(1, bounds=(1, None), doc="Number of MPI processes used in anaddb")

    nqsmall = param.Integer(10, bounds=(1, None), doc="Number of divisions for smallest vector to generate Q-mesh")
    ndivsm = param.Integer(5, bounds=(None, None), doc="Number of divisions for smallest segment in q-path")
    lo_to_splitting = param.ObjectSelector(default="automatic", objects=["automatic", True, False])
    chneut = param.ObjectSelector(default=1, objects=[0, 1, 2], doc="Abinit variable")
    dipdip = param.ObjectSelector(default=1, objects=[0, 1, -1], doc="Abinit variable")
    asr = param.ObjectSelector(default=2, objects=[0, 1, 2], doc="Abinit variable")
    units = param.ObjectSelector(default="eV", objects=["eV", "meV", "Ha", "cm-1", "Thz"], doc="Energy units")

    dos_method = param.ObjectSelector(default="tetra", objects=["tetra", "gaussian"], doc="Integration method for DOS")
    temp_range = pnw.RangeSlider(name="T-range", start=0.0, end=1000, value=(0.0, 300.0), step=20)

    gamma_ev = param.Number(1e-4, bounds=(1e-20, None), doc="Phonon linewidth in eV")
    w_range = pnw.RangeSlider(name="Frequency range (eV)", start=0.0, end=1.0,
                              value=(0.0, 0.1), step=0.001)

    nqsmall_list = pnw.LiteralInput(name='nsmalls (python list)', value=[10, 20, 30], type=list)
    #nqqpt = pnw.LiteralInput(name='nsmalls (list)', value=[10, 20, 30], type=list)

    warning = pn.pane.Markdown(
"""
Refresh the page if the plotly figure is not shown.

Note that widgets for input variables such as *asr*, *chneut*, *dipdip*, *dos_method*, *etc.*
are **shared by the different tabs**.

This means that if you change the value of one of these variables in the active tab,
the same value will **automagically** appear in the other tabs yet the results/figures
are not automatically recomputed when you change the value.

In other words, if you change some variable in the active tab and then you move to another tab,
the results/figures are stil computed with the **old input variables** and you will need to
recompute the new results by clicking the button.
""", name="warning")


    # Base buttons
    plot_check_asr_dipdip_btn = pnw.Button(name="Compute phonons with/wo ASR and DIPDIP", button_type='primary')

    def kwargs_for_anaget_phbst_and_phdos_files(self, **extra_kwargs):
        """
        Return the parameters require to invoke anaget_phbst_and_phdos_files
        Additional kwargs can be specified if needed.
        """
        d = dict(nqsmall=self.nqsmall, qppa=None, ndivsm=self.ndivsm,
                 line_density=None, asr=self.asr, chneut=self.chneut, dipdip=self.dipdip,
                 dos_method=self.dos_method, lo_to_splitting=self.lo_to_splitting,
                 verbose=self.verbose, mpi_procs=self.mpi_procs)

        if extra_kwargs: d.update(extra_kwargs)

        return d


class DdbFilePanel(HasStructureParams, HasAnaddbParams):
    """
    A panel to analyze a |DdbFile|.
    Provides widgets to invoke anaddb and visualize the results.
    """

    # Buttons
    get_epsinf_btn = pnw.Button(name="Compute", button_type='primary')
    plot_phbands_btn = pnw.Button(name="Plot Bands and DOS", button_type='primary')
    plot_eps0w_btn = pnw.Button(name="Plot eps0(omega)", button_type='primary')
    plot_vsound_btn = pnw.Button(name="Calculate speed of sound", button_type='primary')

    plot_ifc_btn = pnw.Button(name="Compute IFC(R)", button_type='primary')
    plot_phbands_quad_btn = pnw.Button(name="Plot PHbands with/without quadrupoles", button_type='primary')
    plot_dos_vs_qmesh_btn = pnw.Button(name="Plot PHDos vs Qmesh", button_type='primary')

    def __init__(self, ddb, **params):
        super().__init__(**params)
        self.ddb = ddb

    @property
    def structure(self):
        """Structure object provided by subclass."""
        return self.ddb.structure

    @param.depends('get_epsinf_btn.clicks')
    def get_epsinf(self):
        """Compute eps_infinity and Born effective charges from DDB."""
        if self.get_epsinf_btn.clicks == 0: return

        with ButtonContext(self.get_epsinf_btn):
            epsinf, becs = self.ddb.anaget_epsinf_and_becs(chneut=self.chneut,
                                                           mpi_procs=self.mpi_procs, verbose=self.verbose)

            gen, inp = self.ddb.anaget_dielectric_tensor_generator(asr=self.asr, chneut=self.chneut, dipdip=self.dipdip,
                                                                   mpi_procs=self.mpi_procs, verbose=self.verbose,
                                                                   return_input=True)

            # Fill column
            col = pn.Column(sizing_mode='stretch_width'); ca = col.append
            #df_kwargs = dict(auto_edit=False, autosize_mode="fit_viewport")
            df_kwargs = {}

            eps0 = gen.tensor_at_frequency(w=0, gamma_ev=self.gamma_ev)
            ca(r"## $\epsilon^0$ in Cart. coords (computed with Gamma_eV):")
            ca(dfc(eps0.get_dataframe(cmode="real"), **df_kwargs))
            ca(r"## $\epsilon^\infty$ in Cart. coords:")
            ca(dfc(epsinf.get_dataframe(), **df_kwargs))
            ca("## Born effective charges in Cart. coords:")
            ca(dfc(becs.get_voigt_dataframe(), **df_kwargs))
            ca("## Anaddb input file.")
            ca(pn.pane.HTML(inp._repr_html_()))

            return col

    @param.depends('plot_eps0w_btn.clicks')
    def plot_eps0w(self):
        """Compute eps0(omega) from DDB and plot the results."""
        if self.plot_eps0w_btn.clicks == 0: return

        with ButtonContext(self.plot_eps0w_btn):
            gen, inp = self.ddb.anaget_dielectric_tensor_generator(asr=self.asr, chneut=self.chneut, dipdip=self.dipdip,
                                                                   mpi_procs=self.mpi_procs, verbose=self.verbose,
                                                                   return_input=True)
            ws = self.w_range.value
            w_max = ws[1]
            if w_max == 1.0: w_max = None # Will compute w_max in plot routine from ph freqs.

            def p(component, reim):
                return gen.plot(w_min=ws[0], w_max=w_max, gamma_ev=self.gamma_ev, num=500, component=component,
                                reim=reim, units=self.units, **self.mpl_kwargs)

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
        """Compute phonon bands and DOSes from DDB and plot the results."""
        if self.plot_phbands_btn.clicks == 0: return

        with ButtonContext(self.plot_phbands_btn):
            # Computing phbands
            kwargs = self.kwargs_for_anaget_phbst_and_phdos_files(return_input=True)

            with self.ddb.anaget_phbst_and_phdos_files(**kwargs) as g:
                phbst_file, phdos_file = g
                phbands, phdos = phbst_file.phbands, phdos_file.phdos

            # Fill column
            col = pn.Column(sizing_mode='stretch_width'); ca = col.append

            ca("## Phonon band structure and DOS:")
            ca(ply(phbands.plotly_with_phdos(phdos, units=self.units, show=False)))
            #ca(mpl(phbands.plot_with_phdos(phdos, units=self.units, **self.mpl_kwargs)))
            #ca(mpl(phdos_file.plot_pjdos_type(units=self.units, exchange_xy=True, **self.mpl_kwargs)))
            #ca(mpl(phdos_file.msqd_dos.plot(units=self.units, **self.mpl_kwargs)))
            temps = self.temp_range.value
            ca("## Thermodynamic properties in the harmonic approximation:")
            #ca(phdos.plot_harmonic_thermo(tstart=temps[0], tstop=temps[1], num=50, **self.mpl_kwargs))
            ca(ply(phdos.plotly_harmonic_thermo(tstart=temps[0], tstop=temps[1], num=50, show=False)))
            #msqd_dos.plot_tensor(**self.mpl_kwargs)

            # Add Anaddb input file
            ca("## Anaddb input file:")
            ca(self.html_with_clipboard_btn(g.input._repr_html_()))

            return col

    @param.depends('plot_vsound_btn.clicks')
    def plot_vsound(self):
        """
        Compute the speed of sound by fitting phonon frequencies
        along selected directions by linear least-squares fit.
        """
        if self.plot_vsound_btn.clicks == 0: return

        with ButtonContext(self.plot_vsound_btn):
            from abipy.dfpt.vsound import SoundVelocity
            sv = SoundVelocity.from_ddb(self.ddb.filepath, num_points=20, qpt_norm=0.1,
                                        ignore_neg_freqs=True, asr=self.asr, chneut=self.chneut, dipdip=self.dipdip,
                                        verbose=self.verbose, mpi_procs=self.mpi_procs)

            # Insert results in grid.
            gspec = pn.GridSpec(sizing_mode='scale_width')
            gspec[0, :1] = sv.get_dataframe()
            gspec[1, :1] = sv.plot(**self.mpl_kwargs)

            return gspec

    @param.depends('plot_check_asr_dipdip_btn.clicks')
    def plot_without_asr_dipdip(self):
        """
        Compare phonon bands and DOSes computed with/without the acoustic sum rule
        and the treatment of the dipole-dipole interaction in the dynamical matrix.
        Requires DDB file with eps_inf, BECS.
        """
        if self.plot_check_asr_dipdip_btn.clicks == 0: return

        with ButtonContext(self.plot_check_asr_dipdip_btn):
            asr_plotter = self.ddb.anacompare_asr(asr_list=(0, 2), chneut_list=(1, ), dipdip=1,
                                                  lo_to_splitting=self.lo_to_splitting,
                                                  nqsmall=self.nqsmall, ndivsm=self.ndivsm,
                                                  dos_method=self.dos_method, ngqpt=None,
                                                  verbose=self.verbose, mpi_procs=self.mpi_procs)

            dipdip_plotter = self.ddb.anacompare_dipdip(chneut_list=(1,), asr=2, lo_to_splitting=self.lo_to_splitting,
                                                        nqsmall=self.nqsmall, ndivsm=self.ndivsm,
                                                        dos_method=self.dos_method, ngqpt=None,
                                                        verbose=self.verbose, mpi_procs=self.mpi_procs)

            # Fill column
            col = pn.Column(sizing_mode='stretch_width'); ca = col.append

            ca("## Phonon bands and DOS with/wo acoustic sum rule:")
            #ca(mpl(asr_plotter.plot(**self.mpl_kwargs)))
            ca(ply(asr_plotter.combiplotly(show=False)))
            ca("## Phonon bands and DOS with/without the treatment of the dipole-dipole interaction:")
            #ca(mpl(dipdip_plotter.plot(**self.mpl_kwargs)))
            ca(ply(dipdip_plotter.combiplotly(show=False)))

            return col

    @param.depends('plot_dos_vs_qmesh_btn.clicks')
    def plot_dos_vs_qmesh(self):
        """
        Compare phonon DOSes computed with/without the inclusion
        of the dipole-quadrupole and quadrupole-quadrupole terms in the dynamical matrix.
        Requires DDB file with eps_inf, BECS and dynamical quadrupoles.
        """
        if self.plot_dos_vs_qmesh_btn.clicks == 0: return

        with ButtonContext(self.plot_dos_vs_qmesh_btn):
            num_cpus = 1
            #print(self.nqsmall_list.value)
            r = self.ddb.anacompare_phdos(self.nqsmall_list.value, asr=self.asr, chneut=self.chneut, dipdip=self.dipdip,
                                          dos_method=self.dos_method, ngqpt=None,
                                          verbose=self.verbose, num_cpus=num_cpus, stream=sys.stdout)

            #r.phdoses: List of |PhononDos| objects

            # Fill column
            col = pn.Column(sizing_mode='stretch_width'); ca = col.append
            ca("## Phonon DOSes obtained with different q-meshes:")
            ca(ply(r.plotter.combiplotly(show=False)))

            ca("## Convergence of termodynamic properties.")
            temps = self.temp_range.value
            ca(mpl(r.plotter.plot_harmonic_thermo(tstart=temps[0], tstop=temps[1], num=50,
                                                  units=self.units, **self.mpl_kwargs)))

            return col

    @param.depends('plot_phbands_quad_btn.clicks')
    def plot_phbands_quad(self):
        """
        Compare phonon bands and DOSes computed with/without the inclusion
        of the dipole-quadrupole and quadrupole-quadrupole terms in the dynamical matrix.
        Requires DDB file with eps_inf, BECS and dynamical quadrupoles.
        """
        if self.plot_phbands_quad_btn.clicks == 0: return

        with ButtonContext(self.plot_dos_vs_qmesh_btn):
            plotter = self.ddb.anacompare_quad(asr=self.asr, chneut=self.chneut, dipdip=self.dipdip,
                                               lo_to_splitting=self.lo_to_splitting,
                                               nqsmall=0, ndivsm=self.ndivsm, dos_method=self.dos_method, ngqpt=None,
                                               verbose=self.verbose, mpi_procs=self.mpi_procs)

            # Fill column
            col = pn.Column(sizing_mode='stretch_width'); ca = col.append
            ca("## Phonon Bands obtained with different q-meshes:")
            ca(ply(plotter.combiplotly(show=False)))

            return col

    @param.depends('plot_ifc_btn.clicks')
    def plot_ifc(self):
        if self.plot_ifc_btn.clicks == 0: return

        ifc = self.ddb.anaget_ifc(asr=self.asr, chneut=self.chneut, dipdip=self.dipdip)

        # Fill column
        col = pn.Column(sizing_mode='stretch_width'); ca = col.append
        ca(mpl(ifc.plot_longitudinal_ifc(title="Longitudinal IFCs", **self.mpl_kwargs)))
        ca(mpl(ifc.plot_longitudinal_ifc_short_range(title="Longitudinal IFCs short range", **self.mpl_kwargs)))
        ca(mpl(ifc.plot_longitudinal_ifc_ewald(title="Longitudinal IFCs Ewald", **self.mpl_kwargs)))

        return col

    def get_panel(self):
        """Return tabs with widgets to interact with the DDB file."""

        tabs = pn.Tabs(); app = tabs.append

        app(("Summary", pn.Row(
            bkw.PreText(text=self.ddb.to_string(verbose=self.verbose), sizing_mode="scale_both"))
        ))
        app(("PH-bands", pn.Row(
            pn.Column("# PH-bands options",
                      *self.pws("nqsmall", "ndivsm", "asr", "chneut", "dipdip",
                                "lo_to_splitting", "dos_method", "temp_range", "plot_phbands_btn",
                                self.helpc("plot_phbands_and_phdos")),
                      ),
            self.plot_phbands_and_phdos)
        ))
        app(("BECs", pn.Row(
            pn.Column("# Born effective charges options",
                      *self.pws("asr", "chneut", "dipdip", "gamma_ev", "get_epsinf_btn", self.helpc("get_epsinf")),
                     ),
            self.get_epsinf)
        ))
        app(("eps0", pn.Row(
            pn.Column("# epsilon_0",
                      *self.pws("asr", "chneut", "dipdip", "gamma_ev", "w_range", "plot_eps0w_btn",
                                self.helpc("plot_eps0w")),
                      ),
            self.plot_eps0w)
        ))
        app(("Speed of sound", pn.Row(
            pn.Column("# Speed of sound options",
                      *self.pws("asr", "chneut", "dipdip", "plot_vsound_btn", self.helpc("plot_vsound")),
                      ),
            self.plot_vsound)
        ))
        app(("ASR & DIPDIP", pn.Row(
            pn.Column("# ASR & DIPDIP options",
                      *self.pws("nqsmall", "ndivsm", "dos_method", "plot_check_asr_dipdip_btn",
                                self.helpc("plot_without_asr_dipdip")),
                      ),
            self.plot_without_asr_dipdip)
        ))
        app(("DOS vs q-mesh", pn.Row(
            pn.Column("# DOS vs q-mesh options",
                      *self.pws("asr", "chneut", "dipdip", "dos_method", "nqsmall_list",
                                "temp_range", "plot_dos_vs_qmesh_btn", self.helpc("plot_dos_vs_qmesh")),
                      ),
            self.plot_dos_vs_qmesh)
        ))
        app(("Quadrupoles", pn.Row(
            pn.Column("# Quadrupoles options",
                      *self.pws("asr", "chneut", "dipdip", "lo_to_splitting", "ndivsm", "dos_method",
                                "plot_phbands_quad_btn", self.helpc("plot_phbands_quad")),
                      ),
            self.plot_phbands_quad)
        ))
        app(("IFCs", pn.Row(
            pn.Column("# IFCs options",
                      *self.pws("asr", "dipdip", "chneut", "plot_ifc_btn", self.helpc("plot_ifc")),
                      ),
            self.plot_ifc)
        ))
        app(self.get_struct_view_tab_entry())
        app(("Global options",
            pn.Column("# Global options",
                *self.pws("units", "mpi_procs", "verbose"),
        )))

        return tabs


class DdbRobotPanel(HasAnaddbParams):
    """
    A panel to analyze multiple |DdbFile| via the low-level API provided by DdbRobot.
    Provides widgets to invoke anaddb and visualize the results.
    """

    # Buttons
    plot_combiplot_btn = pnw.Button(name="Compute", button_type='primary')
    combiplot_check_btn = pnw.CheckButtonGroup(name='Check Button Group',
                                               value=['combiplot'], options=['combiplot', 'gridplot'])

    def __init__(self, robot, **params):
        super().__init__(**params)
        self.robot = robot

    def kwargs_for_anaget_phbst_and_phdos_files(self, **extra_kwargs):
        """Extend method of base class to handle lo_to_splitting"""
        kwargs = super().kwargs_for_anaget_phbst_and_phdos_files(**extra_kwargs)

        if kwargs["lo_to_splitting"] == "automatic":
            if any(not ddb.has_lo_to_data() for ddb in self.robot.abifiles):
                kwargs["lo_to_splitting"] = False
                if self.verbose:
                    print("Setting lo_to_splitting to False since at least one DDB file does not have LO-TO data.")

        return kwargs

    @param.depends('plot_combiplot_btn.clicks')
    def plot_combiplot(self):
        """Plot phonon band structures."""
        if self.plot_combiplot_btn.clicks == 0: return

        kwargs = self.kwargs_for_anaget_phbst_and_phdos_files()

        with ButtonContext(self.plot_combiplot_btn):
            #TODO: Recheck lo-to automatic.
            r = self.robot.anaget_phonon_plotters(**kwargs)
            #r = self.robot.anaget_phonon_plotters()

            # Fill column
            col = pn.Column(sizing_mode='stretch_width'); ca = col.append
            if "combiplot" in self.combiplot_check_btn.value:
                ca("## Phonon band structure and DOS:")
                ca(ply(r.phbands_plotter.combiplotly(units=self.units, show=False)))

            if "gridplot" in self.combiplot_check_btn.value:
                ca("## Phonon band structure and DOS:")
                ca(ply(r.phbands_plotter.gridplotly(units=self.units, show=False)))

            #if "temp_range" in self.combiplot_check_btn.value:
            #temps = self.temp_range.value
            #ca("## Thermodynamic properties in the harmonic approximation:")
            ##ca(phdos.plot_harmonic_thermo(tstart=temps[0], tstop=temps[1], num=50, **self.mpl_kwargs))
            #ca(ply(phdos.plotly_harmonic_thermo(tstart=temps[0], tstop=temps[1], num=50, show=False)))


            return col

    #@param.depends('get_epsinf_btn.clicks')
    #def get_epsinf(self):
    #    """Compute eps_infinity and Born effective charges from DDB."""
    #    if self.get_epsinf_btn.clicks == 0: return

    #    with ButtonContext(self.get_epsinf_btn):
    #        epsinf, becs = self.ddb.anaget_epsinf_and_becs(chneut=self.chneut,
    #                                                       mpi_procs=self.mpi_procs, verbose=self.verbose)

    #        gen, inp = self.ddb.anaget_dielectric_tensor_generator(asr=self.asr, chneut=self.chneut, dipdip=self.dipdip,
    #                                                               mpi_procs=self.mpi_procs, verbose=self.verbose,
    #                                                               return_input=True)

    #        # Fill column
    #        col = pn.Column(sizing_mode='stretch_width'); ca = col.append
    #        df_kwargs = dict(auto_edit=False, autosize_mode="fit_viewport")
    #        #l = pn.pane.LaTeX

    #        eps0 = gen.tensor_at_frequency(w=0, gamma_ev=self.gamma_ev)
    #        ca(r"## $\epsilon^0$ in Cart. coords (computed with Gamma_eV):")
    #        ca(dfc(eps0.get_dataframe(cmode="real"), **df_kwargs))
    #        ca(r"## $\epsilon^\infty$ in Cart. coords:")
    #        ca(dfc(epsinf.get_dataframe(), **df_kwargs))
    #        ca("## Born effective charges in Cart. coords:")
    #        ca(dfc(becs.get_voigt_dataframe(), **df_kwargs))
    #        ca("## Anaddb input file.")
    #        ca(pn.pane.HTML(inp._repr_html_()))

    #        return col

    #@param.depends('plot_eps0w_btn.clicks')
    #def plot_eps0w(self):
    #    """Compute eps0(omega) from DDB and plot the results."""
    #    if self.plot_eps0w_btn.clicks == 0: return

    #    with ButtonContext(self.plot_eps0w_btn):
    #        gen, inp = self.ddb.anaget_dielectric_tensor_generator(asr=self.asr, chneut=self.chneut, dipdip=self.dipdip,
    #                                                               mpi_procs=self.mpi_procs, verbose=self.verbose,
    #                                                               return_input=True)
    #        ws = self.w_range.value
    #        w_max = ws[1]
    #        if w_max == 1.0: w_max = None # Will compute w_max in plot routine from ph freqs.

    #        def p(component, reim):
    #            return gen.plot(w_min=ws[0], w_max=w_max, gamma_ev=self.gamma_ev, num=500, component=component,
    #                            reim=reim, units=self.units, **self.mpl_kwargs)

    #        # Build grid
    #        gspec = pn.GridSpec(sizing_mode='scale_width')
    #        gspec[0, 0] = p("diag", "re")
    #        gspec[0, 1] = p("diag", "im")
    #        gspec[1, 0] = p("offdiag", "re")
    #        gspec[1, 1] = p("offdiag", "im")
    #        gspec[2, :] = gen.get_oscillator_dataframe(reim="all", tol=1e-6)
    #        # Add HTML pane with input.
    #        gspec[3, 0] = pn.pane.HTML(inp._repr_html_())

    #        return gspec

    #@param.depends('plot_phbands_btn.clicks')
    #def plot_phbands_and_phdos(self, event=None):
    #    """Compute phonon bands and DOSes from DDB and plot the results."""
    #    if self.plot_phbands_btn.clicks == 0: return

    #    with ButtonContext(self.plot_phbands_btn):
    #        # Computing phbands
    #        with self.ddb.anaget_phbst_and_phdos_files(
    #                nqsmall=self.nqsmall, qppa=None, ndivsm=self.ndivsm,
    #                line_density=None, asr=self.asr, chneut=self.chneut, dipdip=self.dipdip,
    #                dos_method=self.dos_method, lo_to_splitting=self.lo_to_splitting,
    #                verbose=self.verbose, mpi_procs=self.mpi_procs, return_input=True) as g:

    #            phbst_file, phdos_file = g
    #            phbands, phdos = phbst_file.phbands, phdos_file.phdos

    #        # Fill column
    #        col = pn.Column(sizing_mode='stretch_width'); ca = col.append

    #        ca("## Phonon band structure and DOS:")
    #        ca(ply(phbands.plotly_with_phdos(phdos, units=self.units, show=False)))
    #        #ca(mpl(phbands.plot_with_phdos(phdos, units=self.units, **self.mpl_kwargs)))
    #        #ca(mpl(phdos_file.plot_pjdos_type(units=self.units, exchange_xy=True, **self.mpl_kwargs)))
    #        #ca(mpl(phdos_file.msqd_dos.plot(units=self.units, **self.mpl_kwargs)))
    #        temps = self.temp_range.value
    #        ca("## Thermodynamic properties in the harmonic approximation:")
    #        #ca(phdos.plot_harmonic_thermo(tstart=temps[0], tstop=temps[1], num=50, **self.mpl_kwargs))
    #        ca(ply(phdos.plotly_harmonic_thermo(tstart=temps[0], tstop=temps[1], num=50, show=False)))
    #        #msqd_dos.plot_tensor(**self.mpl_kwargs)
    #        #self.plot_phbands_btn.button_type = "primary"

    #        # Add HTML pane with input
    #        ca("## Anaddb input file:")
    #        ca(pn.pane.HTML(g.input._repr_html_()))

    #        return col

    #@param.depends('plot_vsound_btn.clicks')
    #def plot_vsound(self):
    #    """
    #    Compute the speed of sound by fitting phonon frequencies
    #    along selected directions by linear least-squares fit.
    #    """
    #    if self.plot_vsound_btn.clicks == 0: return

    #    with ButtonContext(self.plot_vsound_btn):
    #        from abipy.dfpt.vsound import SoundVelocity
    #        sv = SoundVelocity.from_ddb(self.ddb.filepath, num_points=20, qpt_norm=0.1,
    #                                    ignore_neg_freqs=True, asr=self.asr, chneut=self.chneut, dipdip=self.dipdip,
    #                                    verbose=self.verbose, mpi_procs=self.mpi_procs)

    #        # Insert results in grid.
    #        gspec = pn.GridSpec(sizing_mode='scale_width')
    #        gspec[0, :1] = sv.get_dataframe()
    #        gspec[1, :1] = sv.plot(**self.mpl_kwargs)

    #        return gspec

    # THIS OK but I don't think it's very useful
    @param.depends('plot_check_asr_dipdip_btn.clicks')
    def plot_without_asr_dipdip(self):
        """
        Compare phonon bands and DOSes computed with/without the acoustic sum rule
        and the treatment of the dipole-dipole interaction in the dynamical matrix.
        Requires DDB file with eps_inf, BECS.
        """
        if self.plot_check_asr_dipdip_btn.clicks == 0: return

        asr_plotter = PhononBandsPlotter()
        dipdip_plotter = PhononBandsPlotter()

        with ButtonContext(self.plot_check_asr_dipdip_btn):
            for label, ddb in self.robot.items():
                asr_p = ddb.anacompare_asr(asr_list=(0, 2), chneut_list=(1, ), dipdip=1,
                                           lo_to_splitting=self.lo_to_splitting,
                                           nqsmall=self.nqsmall, ndivsm=self.ndivsm,
                                           dos_method=self.dos_method, ngqpt=None,
                                           verbose=self.verbose, mpi_procs=self.mpi_procs,
                                           pre_label=label)

                asr_plotter.append_plotter(asr_p)

                dipdip_p = ddb.anacompare_dipdip(chneut_list=(1,), asr=2, lo_to_splitting=self.lo_to_splitting,
                                                 nqsmall=self.nqsmall, ndivsm=self.ndivsm,
                                                 dos_method=self.dos_method, ngqpt=None,
                                                 verbose=self.verbose, mpi_procs=self.mpi_procs,
                                                 pre_label=label)

                dipdip_plotter.append_plotter(dipdip_p)

            # Fill column
            col = pn.Column(sizing_mode='stretch_width'); ca = col.append

            ca("## Phonon bands and DOS with/wo acoustic sum rule:")
            ca(ply(asr_plotter.combiplotly(show=False)))
            ca("## Phonon bands and DOS with/without the treatment of the dipole-dipole interaction:")
            ca(ply(dipdip_plotter.combiplotly(show=False)))

            return col

    def get_panel(self):
        """Return tabs with widgets to interact with the DDB file."""

        robot = self.robot

        tabs = pn.Tabs(); app = tabs.append

        app(("Summary", pn.Row(
            bkw.PreText(text=robot.to_string(verbose=self.verbose), sizing_mode="scale_both"))
        ))

        dfs = robot.get_structure_dataframes()
        app(("Structures",
            pn.Column("# Lattice daframe", dfc(dfs.lattice),
                      "# Atomic positions", dfc(dfs.coords),)
        ))

        app(("Combiplot", pn.Row(
            pn.Column("# PH-bands options",
                      *self.pws("nqsmall", "ndivsm", "asr", "chneut", "dipdip",
                                "lo_to_splitting", "dos_method", "temp_range",
                                "combiplot_check_btn", "plot_combiplot_btn",
                                self.helpc("plot_combiplot")),
                      ),
            self.plot_combiplot)
        ))
        #app(("PH-bands", pn.Row(
        #    pn.Column("# PH-bands options",
        #              *self.pws("nqsmall", "ndivsm", "asr", "chneut", "dipdip",
        #                        "lo_to_splitting", "dos_method", "temp_range", "plot_phbands_btn",
        #                        self.helpc("plot_phbands_and_phdos")),
        #              ),
        #    self.plot_phbands_and_phdos)
        #))
        #app(("BECs", pn.Row(
        #    pn.Column("# Born effective charges options",
        #              *self.pws("asr", "chneut", "dipdip", "gamma_ev", "get_epsinf_btn",
        #                        self.helpc("get_epsinf")),
        #             ),
        #    self.get_epsinf)
        #))
        #app(("eps0", pn.Row(
        #    pn.Column("# epsilon_0",
        #              *self.pws("asr", "chneut", "dipdip", "gamma_ev", "w_range", "plot_eps0w_btn",
        #                        self.helpc("plot_eps0w")),
        #              ),
        #    self.plot_eps0w)
        #))
        #app(("Speed of sound", pn.Row(
        #    pn.Column("# Speed of sound options",
        #              *self.pws("asr", "chneut", "dipdip", "plot_vsound_btn",
        #                        self.helpc("plot_vsound")),
        #              ),
        #    self.plot_vsound)
        #))
        app(("ASR & DIPDIP", pn.Row(
            pn.Column("# ASR & DIPDIP options",
                      *self.pws("nqsmall", "ndivsm", "dos_method", "plot_check_asr_dipdip_btn",
                                self.helpc("plot_without_asr_dipdip")),
                      ),
            self.plot_without_asr_dipdip)
        ))
        #app(("DOS vs q-mesh", pn.Row(
        #    pn.Column("# DOS vs q-mesh options",
        #              *self.pws("asr", "chneut", "dipdip", "dos_method", "nqsmall_list", "plot_dos_vs_qmesh_btn",
        #                        self.helpc("plot_dos_vs_qmesh")),
        #              ),
        #    self.plot_dos_vs_qmesh)
        #))
        #app(("Quadrupoles", pn.Row(
        #    pn.Column("# Quadrupoles options",
        #              *self.pws("asr", "chneut", "dipdip", "lo_to_splitting", "ndivsm", "dos_method", "plot_phbands_quad_btn",
        #                        self.helpc("plot_phbands_quad")),
        #              ),
        #    self.plot_phbands_quad)
        #))
        #app(("IFCs", pn.Row(
        #    pn.Column("# IFCs options",
        #              *self.pws("asr", "dipdip", "chneut", "plot_ifc_btn",
        #                        self.helpc("plot_ifc")),
        #              ),
        #    self.plot_ifc)
        #))
        app(("Global options",
            pn.Column("# Global options",
                *self.pws("units", "mpi_procs", "verbose"),
        )))

        return tabs
