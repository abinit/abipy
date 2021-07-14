"""Panels for interacting with FATBANDS.nc files."""

#import param
import panel as pn
import panel.widgets as pnw
import bokeh.models.widgets as bkw

from .core import PanelWithElectronBands, ply, mpl, dfc, depends_on_btn_click #, PanelWithEbandsRobot


class FatBandsFilePanel(PanelWithElectronBands):
    """
    Panel with widgets to interact with a |FatBandsFile|.
    """

    def __init__(self, ncfile, **params):
        PanelWithElectronBands.__init__(self, ebands=ncfile.ebands, **params)
        self._ncfile = ncfile

        # Create buttons
        self.plot_fatbands_btn = pnw.Button(name="Plot fatbands", button_type='primary')
        self.plot_fatdos_btn = pnw.Button(name="Plot fatdos", button_type='primary')

    @property
    def ncfile(self):
        return self._ncfile

    @depends_on_btn_click('plot_fatbands_btn')
    def on_plot_fatbands_btn(self):

        sz_mode = "stretch_width"
        col = pn.Column(sizing_mode=sz_mode); ca = col.append

        # Plot the electronic fatbands grouped by atomic type.
        ca("## Electronic fatbands grouped by atomic type:")
        fig = self.ncfile.plotly_fatbands_typeview(e0="fermie", fact=1.0, lmax=None, fig=None, ylims=None,
                                                   blist=None, fontsize=12, band_and_dos=0, show=False)
        ca(ply(fig))

        # Plot the electronic fatbands grouped by L.
        ca("## Electronic fatbands grouped by L:")
        fig = self.ncfile.plotly_fatbands_lview(e0="fermie", fact=1.0, lmax=None, fig=None, ylims=None,
                                                blist=None, fontsize=12, band_and_dos=0, show=False)
        ca(ply(fig))

        return col

    @depends_on_btn_click('plot_fatdos_btn')
    def on_plot_fatdos_btn(self):

        sz_mode = "stretch_width"
        col = pn.Column(sizing_mode=sz_mode); ca = col.append

        # Plot the L-PJDOS grouped by atomic type.
        lmax = 2
        ca("## Electronic fatdos grouped by atomic type:")
        fig = self.ncfile.plotly_pjdos_typeview(lmax=lmax, show=False)
        ca(ply(fig))

        # Plot the L-PJDOS grouped by L.
        ca("## Electronic fatdos grouped by L:")
        fig = self.ncfile.plotly_pjdos_lview(lmax=lmax, show=False)
        ca(ply(fig))

        return col

    def get_panel(self, as_dict=False, **kwargs):
        """Return tabs with widgets to interact with the FATBANDS.nc file."""
        d = {}

        d["Summary"] = pn.Row(bkw.PreText(text=self.ncfile.to_string(verbose=self.verbose), sizing_mode="scale_both"))
        d["e-Bands"] = pn.Row(self.get_plot_ebands_widgets(), self.on_plot_ebands_btn)

        if self.ncfile.ebands.kpoints.is_ibz:
            # Add DOS tab but only if we have a k-sampling.
            d["e-DOS"] = pn.Row(self.get_plot_edos_widgets(), self.on_plot_edos_btn)

            d["FatDos"] = pn.Row(
                self.pws_col(["### Fatdos", "plot_fatdos_btn", self.helpc("on_plot_fatdos_btn")]),
                self.on_plot_fatdos_btn
            )

            # Fermi surface requires a gamma-centered k-mesh
            #if self.ncfile.ebands.supports_fermi_surface:
            if self.ncfile.ebands.kpoints.is_ibz:
                d["Fermi Surface"] = self.get_ifermi_view()

        elif self.ncfile.ebands.kpoints.is_path:
            # NC files have contributions up to L=4 (g channel)
            # but here we are intererested in s,p,d terms only so
            # we use the optional argument lmax
            lmax = 2

            d["FatBands"] = pn.Row(
                self.pws_col(["### Fatbands", "plot_fatbands_btn", self.helpc("on_plot_fatbands_btn")]),
                self.on_plot_fatbands_btn
            )

        else:
            raise ValueError("Neither a IBZ nor k-path!")

        if as_dict: return d

        return self.get_template_from_tabs(d, template=kwargs.get("template", None))


#class FatbandsRobotPanel(PanelWithEbandsRobot):
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
#        return pn.Column(dfc(df), sizing_mode='stretch_width')
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
