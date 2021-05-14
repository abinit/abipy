"""Panels for interacting with FATBANDS.nc files."""

#import param
import panel as pn
#import panel.widgets as pnw
import bokeh.models.widgets as bkw

from .core import PanelWithElectronBands, ply, mpl, dfc #, PanelWithEbandsRobot


class FatBandsFilePanel(PanelWithElectronBands):
    """
    Panel with widgets to interact with a |FatBandsFile|.
    """
    def __init__(self, ncfile, **params):
        self._ncfile = ncfile
        super().__init__(**params)

    @property
    def ncfile(self):
        return self._ncfile

    @property
    def ebands(self):
        """|ElectronBands|."""
        return self._ncfile.ebands

    def get_panel(self, as_dict=False, **kwargs):
        """Return tabs with widgets to interact with the FATBANDS.nc file."""
        d = {}

        d["Summary"] = pn.Row(bkw.PreText(text=self.ncfile.to_string(verbose=self.verbose), sizing_mode="scale_both"))
        d["e-Bands"] = pn.Row(self.get_plot_ebands_widgets(), self.on_plot_ebands_btn)

        if self.ncfile.ebands.kpoints.is_ibz:
            # Add DOS tab but only if k-sampling.
            d["e-DOS"] = pn.Row(self.get_plot_edos_widgets(), self.on_plot_edos_btn)

            # Plot the L-PJDOS grouped by atomic type.
            #self.ncfile.plot_pjdos_typeview(lmax=lmax, **self.mpl_kwargs)
            # Plot the L-PJDOS grouped by L.
            #self.ncfile.plot_pjdos_lview(lmax=lmax, **self.mpl_kwargs)

            # Fermi surface requires a gamma-centered k-mesh
            if self.ncfile.ebands.supports_fermi_surface:
                d["Fermi Surface"] = pn.Row(self.get_plot_fermi_surface_widgets(), self.on_plot_fermi_surface_btn)

        elif self.ncfile.ebands.kpoints.is_path:
            # NC files have contributions up to L=4 (g channel)
            # but here we are intererested in s,p,d terms only so
            # we use the optional argument lmax
            lmax = 2

            # Plot the electronic fatbands grouped by atomic type.
            #self.ncfile.plot_fatbands_typeview(lmax=lmax, **self.mpl_kwargs)
            # Plot the electronic fatbands grouped by L.
            #self.ncfile.plot_fatbands_lview(lmax=lmax, **self.mpl_kwargs)

        else:
            raise ValueError("Neither a IBZ nor k-path!")

        if as_dict: return d
        tabs = pn.Tabs(*d.items())
        return self.get_template_from_tabs(tabs, template=kwargs.get("template", None))


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
