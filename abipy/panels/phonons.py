""""Panels for phonon-related objects."""
import param
import panel as pn

from abipy.panels.core import AbipyParameterized

def _mp(fig):
    return pn.pane.Matplotlib(fig)


class PhononBandsPlotterPanel(AbipyParameterized):

    phbands_plotter_mode = pn.widgets.Select(name="Plot Mode", value="gridplot",
        options=["gridplot", "combiplot", "boxplot", "combiboxplot"]) # "animate",
    phbands_plotter_units = pn.widgets.Select(name="Units", value="eV", options=["eV", "meV", "Ha", "cm-1", "Thz"])
    phbands_plotter_btn = pn.widgets.Button(name="Plot", button_type='primary')

    def __init__(self, plotter, **params):
        super().__init__(**params)
        self.plotter = plotter

    @param.depends("phbands_plotter_btn.clicks")
    def on_phbands_plot_btn(self):
        if self.phbands_plotter_btn.clicks == 0: return
        plot_mode = self.phbands_plotter_mode.value
        plotfunc = getattr(self.plotter, plot_mode, None)
        if plotfunc is None:
            raise ValueError("Don't know how to handle plot_mode: %s" % plot_mode)

        fig = plotfunc(units=self.phbands_plotter_units.value, **self.fig_kwargs)
        df = self.plotter.get_phbands_frame(with_spglib=True)
        return pn.Row(pn.Column(_mp(fig), df), sizing_mode='scale_width')

    def get_panel(self):
        """Return tabs with widgets to interact with the |PhononBandsPlotter|."""
        tabs = pn.Tabs()
        widgets = pn.Column(self.phbands_plotter_mode, self.phbands_plotter_units, self.phbands_plotter_btn)
        tabs.append(("PhbandsPlotter", pn.Row(widgets, self.on_phbands_plot_btn, sizing_mode='scale_width')))

        return tabs
