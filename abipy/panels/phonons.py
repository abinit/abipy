""""Panels for phonon-related objects."""
import param
import panel as pn
import panel.widgets as pnw

from abipy.panels.core import AbipyParameterized, mpl, ply, dfc, depends_on_btn_click


class PhononBandsPlotterPanel(AbipyParameterized):

    def __init__(self, plotter, **params):

        self.phbands_plotter_mode = pnw.Select(name="Plot Mode", value="gridplot",
                                          options=["gridplot", "combiplot", "boxplot", "combiboxplot"]) # "animate",
        self.phbands_plotter_units = pnw.Select(name="Units", value="eV",
                                           options=["eV", "meV", "Ha", "cm-1", "Thz"])
        self.phbands_plotter_btn = pnw.Button(name="Plot", button_type='primary')

        self.plotter = plotter
        super().__init__(**params)

    @depends_on_btn_click('phbands_blotter_btn')
    def on_phbands_plot_btn(self):
        plot_mode = self.phbands_plotter_mode.value
        plotfunc = getattr(self.plotter, plot_mode, None)
        if plotfunc is None:
            raise ValueError("Don't know how to handle plot_mode: %s" % plot_mode)

        fig = plotfunc(units=self.phbands_plotter_units.value, **self.mpl_kwargs)
        df = self.plotter.get_phbands_frame(with_spglib=True)
        return pn.Row(pn.Column(mpl(fig), dfc(df)), sizing_mode='scale_width')

    def get_panel(self, as_dict=False, **kwargs):
        """Return tabs with widgets to interact with the |PhononBandsPlotter|."""
        d = {}

        ws = pn.Column(self.phbands_plotter_mode, self.phbands_plotter_units, self.phbands_plotter_btn)
        d["PhbandsPlotter"] = pn.Row(ws, self.on_phbands_plot_btn, sizing_mode='scale_width')

        if as_dict: return d

        tabs = pn.Tabs(*d.items())
        return self.get_template_from_tabs(tabs, template=kwargs.get("template", None))
