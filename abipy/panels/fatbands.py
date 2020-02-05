"""Panels for interacting with FATBANDS files."""
import param
import panel as pn
import panel.widgets as pnw
import bokeh.models.widgets as bkw

from .core import PanelWithElectronBands, PanelWithEbandsRobot


class FatBandsFilePanel(PanelWithElectronBands):
    """
    Panel with widgets to interact with a |FatBandsFile|.
    """
    def __init__(self, ncfile, **params):
        super().__init__(**params)
        self.ncfile = ncfile

    @property
    def ebands(self):
        """|ElectronBands|."""
        return self.ncfile.ebands

    def get_panel(self):
        """Return tabs with widgets to interact with the DDB file."""
        tabs = pn.Tabs(); app = tabs.append
        app(("Summary", pn.Row(bkw.PreText(text=self.ncfile.to_string(verbose=self.verbose),
                               sizing_mode="scale_both"))))
        app(("e-Bands", pn.Row(self.get_plot_ebands_widgets(), self.on_plot_ebands_btn)))

        # Add DOS tab only if k-sampling.
        if self.ncfile.ebands.kpoints.is_ibz:
            app(("e-DOS", pn.Row(self.get_plot_edos_widgets(), self.on_plot_edos_btn)))

            if self.ncfile.ebands.supports_fermi_surface:
                # Fermi surface requires a gamma-centered k-mesh
                app(("Fermi Surface", pn.Row(self.get_plot_fermi_surface_widgets(), self.on_plot_fermi_surface_btn)))

        return tabs


#class GsrRobotPanel(PanelWithEbandsRobot):
#    """
#    A Panel to interoperate with multiple GSR files.
#    """
#
#    gsr_dataframe_btn = pnw.Button(name="Compute", button_type='primary')
#
#    def __init__(self, robot, **params):
#        super().__init__(**params)
#        self.robot = robot
#
#    @param.depends("gsr_dataframe_btn.clicks")
#    def on_gsr_dataframe_btn(self):
#        if self.gsr_dataframe_btn.clicks == 0: return
#        df = self.robot.get_dataframe(with_geo=True)
#        return pn.Column(self._df(df), sizing_mode='stretch_width')
#
#    def get_panel(self):
#        """Return tabs with widgets to interact with the |GsrRobot|."""
#        tabs = pn.Tabs(); app = tabs.append
#        app(("Summary", pn.Row(bkw.PreText(text=self.robot.to_string(verbose=self.verbose),
#                               sizing_mode="scale_both"))))
#        app(("e-Bands", pn.Row(self.get_ebands_plotter_widgets(), self.on_ebands_plotter_btn)))
#
#        # Add e-DOS tab only if all ebands have k-sampling.
#        if all(abifile.ebands.kpoints.is_ibz for abifile in self.robot.abifiles):
#            app(("e-DOS", pn.Row(self.get_edos_plotter_widgets(), self.on_edos_plotter_btn)))
#
#        app(("GSR-DataFrame", pn.Row(self.gsr_dataframe_btn, self.on_gsr_dataframe_btn)))
#
#        return tabs
