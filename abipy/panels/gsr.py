"""Panels to interact with GSR files."""

import param
import panel as pn
import panel.widgets as pnw
import bokeh.models.widgets as bkw

from .core import (PanelWithElectronBands, HasStructureParams, PanelWithNcFile,
  PanelWithEbandsRobot, ButtonContext, ply, mpl, dfc, depends_on_btn_click)


class GsrFilePanel(PanelWithElectronBands, HasStructureParams, PanelWithNcFile):
    """
    Panel with widgets to interact with a |GsrFile|.
    """
    def __init__(self, gsr, **params):
        self.gsr = gsr
        super().__init__(**params)

    @property
    def structure(self):
        """|Structure| object"""
        return self.gsr.structure

    @property
    def ebands(self):
        """|ElectronBands| object."""
        return self.gsr.ebands

    @property
    def ncfile(self):
        """This for for the PanelWithNcFile mixin"""
        return self.gsr

    def get_panel(self, as_dict=False, **kwargs):
        """Return tabs with widgets to interact with the GSR file."""
        d = {}

        d["Summary"] = pn.Row(
            bkw.PreText(text=self.gsr.to_string(verbose=self.verbose),  sizing_mode="scale_both")
        )
        d["e-Bands"] = pn.Row(
            pn.Column("### e-Bands Options", self.get_plot_ebands_widgets(), self.helpc("on_plot_ebands_btn"),
            self.on_plot_ebands_btn),
        )
        # Add DOS tab but only if k-sampling.
        kpoints = self.gsr.ebands.kpoints
        if kpoints.is_ibz:
            d["e-DOS"] = pn.Row(
                pn.Column("### Options", self.get_plot_edos_widgets(), self.helpc("on_plot_edos_btn")),
                self.on_plot_edos_btn
            )

            d["SKW"] = self.get_plot_skw_widgets()

            #if self.gsr.ebands.supports_fermi_surface:
            #    # Fermi surface requires Gamma-centered k-mesh
            #    app(("Fermi Surface", pn.Row(
            #        pn.Column("# Options",
            #            self.get_plot_fermi_surface_widgets(),
            #            self.helpc("on_plot_fermi_surface_btn"),
            #    ),
            #    self.on_plot_fermi_surface_btn)
            #    ))

        d["Structure"] = self.get_struct_view_tab_entry()
        # TODO
        #app(("NcFile", self.get_ncfile_panel()))

        #app(("Global", pn.Row(
        #    pn.Column("# Global options",
        #              *self.pws("units", "mpi_procs", "verbose"),
        #              ),
        #    self.get_software_stack())
        #))

        if as_dict: return d
        tabs = pn.Tabs(*d.items())
        return self.get_template_from_tabs(tabs, template=kwargs.get("template", None))


class GsrRobotPanel(PanelWithEbandsRobot):
    """
    A Panel to interoperate with multiple GSR files.
    """
    gsr_dataframe_btn = pnw.Button(name="Compute", button_type='primary')

    transpose_gsr_dataframe = pnw.Checkbox(name='Transpose GSR dataframe')

    def __init__(self, robot, **params):
        super().__init__(**params)
        self.robot = robot

    @depends_on_btn_click('gsr_dataframe_btn')
    def on_gsr_dataframe_btn(self):
        df = self.robot.get_dataframe(with_geo=True)
        transpose = self.transpose_gsr_dataframe.value
        return pn.Column(dfc(df, transpose=transpose), sizing_mode='stretch_width')

    def get_panel(self, **kwargs):
        """Return tabs with widgets to interact with the |GsrRobot|."""
        d = {}
        d["Summary"] = pn.Row(bkw.PreText(text=self.robot.to_string(verbose=self.verbose),
                               sizing_mode="scale_both"))
        d["e-Bands"] = pn.Row(self.get_ebands_plotter_widgets(), self.on_ebands_plotter_btn)

        # Add e-DOS tab but only if all ebands have k-sampling.
        if all(abifile.ebands.kpoints.is_ibz for abifile in self.robot.abifiles):
            d["e-DOS"] = pn.Row(self.get_edos_plotter_widgets(), self.on_edos_plotter_btn)

        d["GSR-dataframe"] = pn.Row(
            pn.Column(self.transpose_gsr_dataframe, self.gsr_dataframe_btn),
            self.on_gsr_dataframe_btn)

        tabs = pn.Tabs(*d.items())
        return self.get_template_from_tabs(tabs, template=kwargs.get("template", None))
