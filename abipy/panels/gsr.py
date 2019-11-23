"""Panels for GSR files."""
import param
import panel as pn
import panel.widgets as pnw
import bokeh.models.widgets as bkw

from .core import PanelWithElectronBands, PanelWithEbandsRobot

def _df(df):
    return pnw.DataFrame(df, disabled=True)


class GsrFilePanel(PanelWithElectronBands):
    """
    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GsrFilePanel
    """
    def __init__(self, gsr, **params):
        super().__init__(**params)
        self.gsr = gsr

    @property
    def ebands(self):
        return self.gsr.ebands

    def get_panel(self):
        """Return tabs with widgets to interact with the DDB file."""
        tabs = pn.Tabs()
        tabs.append(("Summary", pn.Row(bkw.PreText(text=self.gsr.to_string(verbose=self.verbose),
                     sizing_mode="scale_both"))))
        tabs.append(("e-Bands", pn.Row(self.get_plot_ebands_widgets(), self.on_plot_ebands_btn)))
        # Add DOS tab only if k-sampling.
        if self.gsr.ebands.kpoints.is_ibz:
            tabs.append(("e-DOS", pn.Row(self.get_plot_edos_widgets(), self.on_plot_edos_btn)))

        return tabs


class GsrRobotPanel(PanelWithEbandsRobot):
    """

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GsrRobotPanel
    """

    gsr_dataframe_btn = pnw.Button(name="Compute", button_type='primary')

    def __init__(self, robot, **params):
        super().__init__(**params)
        self.robot = robot

    @param.depends("gsr_dataframe_btn.clicks")
    def on_gsr_dataframe_btn(self):
        if self.gsr_dataframe_btn.clicks == 0: return
        df = self.robot.get_dataframe(with_geo=True)
        return pn.Column(_df(df), sizing_mode='stretch_width')

    def get_panel(self):
        """Return tabs with widgets to interact with the |GsrRobot|."""
        tabs = pn.Tabs()
        tabs.append(("Summary", pn.Row(bkw.PreText(text=self.robot.to_string(verbose=self.verbose),
                     sizing_mode="scale_both"))))
        tabs.append(("e-Bands", pn.Row(self.get_ebands_plotter_widgets(), self.on_ebands_plotter_btn)))
        # Add DOS tab only if all ebands have k-sampling.
        if all(abifile.ebands.kpoints.is_ibz for abifile in self.robot.abifiles):
            tabs.append(("e-DOS", pn.Row(self.get_edos_plotter_widgets(), self.on_edos_plotter_btn)))

        tabs.append(("GSR-DataFrame", pn.Row(self.gsr_dataframe_btn, self.on_gsr_dataframe_btn)))

        return tabs
