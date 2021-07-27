"""Panels to interact with GSR files."""

import param
import panel as pn
import panel.widgets as pnw
import bokeh.models.widgets as bkw

from .core import (PanelWithElectronBands,
  PanelWithEbandsRobot, ButtonContext, ply, mpl, dfc, depends_on_btn_click)


class GsrFilePanel(PanelWithElectronBands):
    """
    Panel with widgets to interact with a |GsrFile|.
    """
    def __init__(self, gsr, **params):
        PanelWithElectronBands.__init__(self, ebands=gsr.ebands, **params)
        self.gsr = gsr

    def get_panel(self, as_dict=False, **kwargs):
        """Return tabs with widgets to interact with the GSR file."""
        d = {}

        d["Summary"] = pn.Row(
            bkw.PreText(text=self.gsr.to_string(verbose=self.verbose), sizing_mode="scale_both")
        )
        d["e-Bands"] = self.get_plot_ebands_view()

        kpoints = self.gsr.ebands.kpoints
        if kpoints.is_ibz:
            # Add DOS tab but only if k-sampling.
            d["e-DOS"] = self.get_plot_edos_view()
            d["SKW"] = self.get_skw_view()

            if not self.gsr.ebands.isnot_ibz_sampling():
                d["ifermi"] = self.get_ifermi_view()
                #d["fsviewer"] = self.get_fsviewer_view()

        d["Structure"] = self.get_structure_view()
        d["NcFile"] = self.gsr.get_ncfile_view()

        # TODO
        #d["Global"] = pn.Row(
        #    pn.Column("# Global options",
        #              *self.pws("units", "mpi_procs", "verbose"),
        #              ),
        #    self.get_software_stack())
        #))

        if as_dict: return d
        return self.get_template_from_tabs(d, template=kwargs.get("template", None))


class GsrRobotPanel(PanelWithEbandsRobot):
    """
    A Panel to interact with multiple GSR files.
    """

    def __init__(self, robot, **params):
        PanelWithEbandsRobot.__init__(self, robot=robot, **params)

        self.gsr_dataframe_btn = pnw.Button(name="Compute", button_type='primary')
        self.transpose_gsr_dataframe = pnw.Checkbox(name='Transpose GSR dataframe')

    @depends_on_btn_click('gsr_dataframe_btn')
    def on_gsr_dataframe_btn(self):
        df = self.robot.get_dataframe(with_geo=True)
        transpose = self.transpose_gsr_dataframe.value

        return pn.Column(dfc(df, transpose=transpose), sizing_mode='stretch_width')

    def get_panel(self, as_dict=False, **kwargs):
        """Return tabs with widgets to interact with the |GsrRobot|."""
        d = {}

        d["Summary"] = pn.Row(bkw.PreText(text=self.robot.to_string(verbose=self.verbose),
                               sizing_mode="scale_both"))
        d["Plot-eBands"] = pn.Row(self.get_ebands_plotter_widgets(), self.on_ebands_plotter_btn)

        # Add e-DOS tab but only if all ebands have k-sampling.
        if all(abifile.ebands.kpoints.is_ibz for abifile in self.robot.abifiles):
            d["Plot-eDOS"] = pn.Row(self.get_edos_plotter_widgets(), self.on_edos_plotter_btn)

        d["Dataframe"] = pn.Row(
            pn.Column(self.transpose_gsr_dataframe, self.gsr_dataframe_btn),
            self.on_gsr_dataframe_btn)

        if as_dict: return d

        return self.get_template_from_tabs(d, template=kwargs.get("template", None))
