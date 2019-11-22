""""Base classes and mixins for AbiPy panels."""
#import abc
import param
import panel as pn
import bokeh.models.widgets as bw

from monty.functools import lazy_property


def _mp(fig):
    return pn.pane.Matplotlib(fig)

def _df(df):
    return pn.widgets.DataFrame(df, disabled=True)


class AbipyParameterized(param.Parameterized):

    verbose = param.Integer(0, bounds=(0, None), doc="Verbosity Level")
    mpi_procs = param.Integer(1, bounds=(1, None), doc="Number of MPI processes used when runnin Fortran code")
    #fontsize =

    #def get_global_widgets(self):
    #    return pn.Column(self.verbose, self.mpi_procs)

    @lazy_property
    def fig_kwargs(self):
        return dict(show=False)
        #TODO requires new pymatgen
        #return dict(show=False, fig_close=True)


#class PanelWithNcFile(AbipyParameterized): #, metaclass=abc.ABCMeta):
#    """
#    This frame allows the user to inspect the dimensions and the variables reported in  a netcdf file.
#    Tab showing information on the netcdf file.
#    """


class PanelWithElectronBands(AbipyParameterized): #, metaclass=abc.ABCMeta):

    # Bands plot
    with_gaps = pn.widgets.Checkbox(name='Show gaps')
    #ebands_ylims
    #ebands_e0
    # e0: Option used to define the zero of energy in the band structure plot. Possible values:
    #     - `fermie`: shift all eigenvalues to have zero energy at the Fermi energy (`self.fermie`).
    #     -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
    #     -  None: Don't shift energies, equivalent to e0=0

    #set_fermie_to_vbm
    plot_ebands_btn = pn.widgets.Button(name="Plot e-bands", button_type='primary')

    # DOS plot.
    edos_method = pn.widgets.Select(name="e-DOS method", options=["gaussian", "tetra"])
    edos_step = pn.widgets.Spinner(name='e-DOS step (eV)', value=0.1, step=0.05, start=1e-6, end=None)
    edos_width = pn.widgets.Spinner(name='e-DOS Gaussian broadening (eV)', value=0.2, step=0.05, start=1e-6, end=None)
    plot_edos_btn = pn.widgets.Button(name="Plot e-DOS", button_type='primary')

    #@abc.abstractproperty
    #def ebands(self):
    #    """Returns the |ElectronBands| object."""

    def get_plot_ebands_widgets(self):
        """Widgets to plot ebands."""
        return pn.Column(self.with_gaps,
                         self.plot_ebands_btn)

    @param.depends('plot_ebands_btn.clicks')
    def on_plot_ebands_btn(self):
        """Button triggering ebands plot."""
        if self.plot_ebands_btn.clicks == 0: return
        fig1 = self.ebands.plot(e0="fermie", ylims=None,
            with_gaps=self.with_gaps.value, max_phfreq=None, fontsize=8, **self.fig_kwargs)

        fig2 = self.ebands.kpoints.plot(**self.fig_kwargs)
        row = pn.Row(_mp(fig1), _mp(fig2)) #, sizing_mode='scale_width')
        text = bw.PreText(text=self.ebands.to_string(verbose=self.verbose))
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
        #print(edos)
        return pn.Row(_mp(fig), sizing_mode='scale_width')


class BaseRobotPanel(AbipyParameterized):
    """pass"""



class PanelWithEbandsRobot(BaseRobotPanel): #, metaclass=abc.ABCMeta):
    """
    Mixin class for panels with a robot that owns a list of of |ElectronBands|
    """

    # Widgets to plot ebands.
    ebands_plotter_mode = pn.widgets.Select(name="Plot Mode", value="gridplot",
        options=["gridplot", "combiplot", "boxplot", "combiboxplot"]) # "animate",
    ebands_plotter_btn = pn.widgets.Button(name="Plot", button_type='primary')

    # Widgets to plot edos.
    edos_plotter_mode = pn.widgets.Select(name="Plot Mode", value="gridplot",
        options=["gridplot", "combiplot"])
    edos_plotter_btn = pn.widgets.Button(name="Plot", button_type='primary')

    def get_ebands_plotter_widgets(self):
        return pn.Column(self.ebands_plotter_mode, self.ebands_plotter_btn)

    @param.depends("ebands_plotter_btn.clicks")
    def on_ebands_plotter_btn(self):
        if self.ebands_plotter_btn.clicks == 0: return
        ebands_plotter = self.robot.get_ebands_plotter()
        plot_mode = self.ebands_plotter_mode.value
        plotfunc = getattr(ebands_plotter, plot_mode, None)
        if plotfunc is None:
            raise ValueError("Don't know how to handle plot_mode: %s" % plot_mode)

        fig = plotfunc(**self.fig_kwargs)
        df = ebands_plotter.get_ebands_frame(with_spglib=True)
        return pn.Row(pn.Column(_mp(fig), df),
                      sizing_mode='scale_width')

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
        return pn.Row(pn.Column(_mp(fig)), sizing_mode='scale_width')
