""""Panels for phonon-related objects."""
import param
import panel as pn
import panel.widgets as pnw

from abipy.panels.core import AbipyParameterized


class PhononBandsPlotterPanel(AbipyParameterized):

    phbands_plotter_mode = pnw.Select(name="Plot Mode", value="gridplot",
                                      options=["gridplot", "combiplot", "boxplot", "combiboxplot"]) # "animate",
    phbands_plotter_units = pnw.Select(name="Units", value="eV",
                                       options=["eV", "meV", "Ha", "cm-1", "Thz"])
    phbands_plotter_btn = pnw.Button(name="Plot", button_type='primary')

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
        return pn.Row(pn.Column(self._mp(fig), self._df(df)), sizing_mode='scale_width')

    def get_panel(self):
        """Return tabs with widgets to interact with the |PhononBandsPlotter|."""
        tabs = pn.Tabs()
        ws = pn.Column(self.phbands_plotter_mode, self.phbands_plotter_units, self.phbands_plotter_btn)
        tabs.append(("PhbandsPlotter", pn.Row(ws, self.on_phbands_plot_btn, sizing_mode='scale_width')))

        return tabs
