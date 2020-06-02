""""Base classes and mixins for AbiPy panels."""

#import abc
import param
import panel as pn
import panel.widgets as pnw
import bokeh.models.widgets as bkw

from monty.functools import lazy_property


def gen_id(n=1, pre="uuid-"):
    """
    Generate ``n`` universally unique identifiers prepended with ``pre`` string.
    Return string if n == 1 or list of strings if n > 1
    """
    # The HTML4 spec says:
    # ID and NAME tokens must begin with a letter ([A-Za-z]) and may be followed by any number of letters,
    # digits ([0-9]), hyphens ("-"), underscores ("_"), colons (":"), and periods (".").
    import uuid
    if n == 1:
        return pre + str(uuid.uuid4())
    elif n > 1:
        return [pre + str(uuid.uuid4()) for i in range(n)]
    else:
        raise ValueError("n must be > 0 but got %s" % str(n))


def sizing_mode_select(name="sizing_mode", value="scale_both"):
    """
    Widget to select the value of sizing_mode. See https://panel.holoviz.org/user_guide/Customization.html
    """
    return pn.widgets.Select(name=name, value=value, options=["fixed",
                 "stretch_width", "stretch_height", "stretch_both",
                 "scale_height", "scale_width", "scale_both"])


class AbipyParameterized(param.Parameterized):

    verbose = param.Integer(0, bounds=(0, None), doc="Verbosity Level")
    mpi_procs = param.Integer(1, bounds=(1, None), doc="Number of MPI processes used when runnin Fortran code")
    #fontsize =

    #def get_global_widgets(self):
    #    return pn.Column(self.verbose, self.mpi_procs)

    @lazy_property
    def fig_kwargs(self):
        """Default arguments passed to AbiPy plot methods."""
        return dict(show=False, fig_close=True)

    @staticmethod
    def _mp(fig):
        return pn.pane.Matplotlib(fig, sizing_mode="scale_width")

    @staticmethod
    def _df(df, disabled=True, sizing_mode="scale_width"):
        return pn.widgets.DataFrame(df, disabled=disabled, sizing_mode=sizing_mode)


#class PanelWithNcFile(AbipyParameterized): #, metaclass=abc.ABCMeta):
#    """
#    This frame allows the user to inspect the dimensions and the variables reported in  a netcdf file.
#    Tab showing information on the netcdf file.
#    """


class PanelWithElectronBands(AbipyParameterized): #, metaclass=abc.ABCMeta):

    # Bands plot
    with_gaps = pnw.Checkbox(name='Show gaps')
    #ebands_ylims
    #ebands_e0
    # e0: Option used to define the zero of energy in the band structure plot. Possible values:
    #     - `fermie`: shift all eigenvalues to have zero energy at the Fermi energy (`self.fermie`).
    #     -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
    #     -  None: Don't shift energies, equivalent to e0=0
    set_fermie_to_vbm = pnw.Checkbox(name="Set Fermie to VBM")

    plot_ebands_btn = pnw.Button(name="Plot e-bands", button_type='primary')

    # DOS plot.
    edos_method = pnw.Select(name="e-DOS method", options=["gaussian", "tetra"])
    edos_step = pnw.Spinner(name='e-DOS step (eV)', value=0.1, step=0.05, start=1e-6, end=None)
    edos_width = pnw.Spinner(name='e-DOS Gaussian broadening (eV)', value=0.2, step=0.05, start=1e-6, end=None)
    plot_edos_btn = pnw.Button(name="Plot e-DOS", button_type='primary')

    # Fermi surface plot.
    fs_viewer = pnw.Select(name="FS viewer", options=["matplotlib", "xcrysden"])
    plot_fermi_surface_btn = pnw.Button(name="Plot Fermi surface", button_type='primary')

    #@abc.abstractproperty
    #def ebands(self):
    #    """Returns the |ElectronBands| object."""

    def get_plot_ebands_widgets(self):
        """Widgets to plot ebands."""
        return pn.Column(self.with_gaps, self.set_fermie_to_vbm, self.plot_ebands_btn)

    @param.depends('plot_ebands_btn.clicks')
    def on_plot_ebands_btn(self):
        """Button triggering ebands plot."""
        if self.plot_ebands_btn.clicks == 0: return
        if self.set_fermie_to_vbm.value:
            self.ebands.set_fermie_to_vbm()

        fig1 = self.ebands.plot(e0="fermie", ylims=None,
            with_gaps=self.with_gaps.value, max_phfreq=None, fontsize=8, **self.fig_kwargs)

        fig2 = self.ebands.kpoints.plot(**self.fig_kwargs)
        row = pn.Row(self._mp(fig1), self._mp(fig2)) #, sizing_mode='scale_width')
        text = bkw.PreText(text=self.ebands.to_string(verbose=self.verbose))
        return pn.Column(row, text, sizing_mode='scale_width')

    def get_plot_edos_widgets(self):
        """Widgets to compute e-DOS."""
        return pn.Column(self.edos_method, self.edos_step, self.edos_width, self.plot_edos_btn)

    @param.depends('plot_edos_btn.clicks')
    def on_plot_edos_btn(self):
        """Button triggering edos plot."""
        if self.plot_edos_btn.clicks == 0: return
        edos = self.ebands.get_edos(method=self.edos_method.value, step=self.edos_step.value, width=self.edos_width.value)
        fig = edos.plot(**self.fig_kwargs)
        return pn.Row(self._mp(fig), sizing_mode='scale_width')

    def get_plot_fermi_surface_widgets(self):
        """Widgets to compute e-DOS."""
        return pn.Column(self.fs_viewer, self.plot_fermi_surface_btn)

    @param.depends('plot_fermi_surface_btn.clicks')
    def on_plot_fermi_surface_btn(self):
        if self.plot_fermi_surface_btn.clicks == 0: return

        # Cache eb3d
        if hasattr(self, "_eb3d"):
            eb3d = self._eb3d
        else:
            # Build ebands in full BZ.
            eb3d = self._eb3d = self.ebands.get_ebands3d()

        if self.fs_viewer.value == "matplotlib":
            # Use matplotlib to plot isosurfaces corresponding to the Fermi level (default)
            # Warning: requires skimage package, rendering could be slow.
            fig = eb3d.plot_isosurfaces(e0="fermie", cmap=None, **self.fig_kwargs)
            return pn.Row(self._mp(fig), sizing_mode='scale_width')

        else:
            raise ValueError("Invalid choice: %s" % self.fs_viewer.value)

        #elif self.fs_viewer.value == "xcrysden":
            # Alternatively, it's possible to export the data in xcrysden format
            # and then use `xcrysden --bxsf mgb2.bxsf`
            #eb3d.to_bxsf("mgb2.bxsf")
            # If you have mayavi installed, try:
            #eb3d.mvplot_isosurfaces()


class BaseRobotPanel(AbipyParameterized):
    """pass"""


class PanelWithEbandsRobot(BaseRobotPanel): #, metaclass=abc.ABCMeta):
    """
    Mixin class for panels with a robot that owns a list of of |ElectronBands|
    """

    # Widgets to plot ebands.
    ebands_plotter_mode = pnw.Select(name="Plot Mode", value="gridplot",
        options=["gridplot", "combiplot", "boxplot", "combiboxplot"]) # "animate",
    ebands_plotter_btn = pnw.Button(name="Plot", button_type='primary')
    ebands_df_checkbox = pnw.Checkbox(name='With Ebands DataFrame', value=False)

    # Widgets to plot edos.
    edos_plotter_mode = pnw.Select(name="Plot Mode", value="gridplot",
        options=["gridplot", "combiplot"])
    edos_plotter_btn = pnw.Button(name="Plot", button_type='primary')

    def get_ebands_plotter_widgets(self):
        return pn.Column(self.ebands_plotter_mode, self.ebands_df_checkbox, self.ebands_plotter_btn)

    @param.depends("ebands_plotter_btn.clicks")
    def on_ebands_plotter_btn(self):
        if self.ebands_plotter_btn.clicks == 0: return
        ebands_plotter = self.robot.get_ebands_plotter()
        plot_mode = self.ebands_plotter_mode.value
        plotfunc = getattr(ebands_plotter, plot_mode, None)
        if plotfunc is None:
            raise ValueError("Don't know how to handle plot_mode: %s" % plot_mode)

        fig = plotfunc(**self.fig_kwargs)
        col = pn.Column(self._mp(fig), sizing_mode='scale_width')
        if self.ebands_df_checkbox.value:
            df = ebands_plotter.get_ebands_frame(with_spglib=True)
            col.append(self._df(df))

        return pn.Row(col, sizing_mode='scale_width')

    def get_edos_plotter_widgets(self):
        return pn.Column(self.edos_plotter_mode, self.edos_plotter_btn)

    @param.depends("edos_plotter_btn.clicks")
    def on_edos_plotter_btn(self):
        if self.edos_plotter_btn.clicks == 0: return
        edos_plotter = self.robot.get_edos_plotter()
        plot_mode = self.edos_plotter_mode.value
        plotfunc = getattr(edos_plotter, plot_mode, None)
        if plotfunc is None:
            raise ValueError("Don't know how to handle plot_mode: %s" % plot_mode)

        fig = plotfunc(**self.fig_kwargs)
        return pn.Row(pn.Column(self._mp(fig)), sizing_mode='scale_width')
