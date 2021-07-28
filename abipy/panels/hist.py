""""Panels for HIST files."""

#import param
import panel as pn
import panel.widgets as pnw
import bokeh.models.widgets as bkw

from abipy.panels.core import AbipyParameterized, mpl, ply, depends_on_btn_click


class HistFilePanel(AbipyParameterized):
    """
    Panel with widgets to interact with a |HistFile|.
    """

    def __init__(self, hist, **params):
        self.hist = hist

        _what_list = ["abc", "angles", "energy", "volume", "pressure", "forces"]
        self.what_list = pnw.CheckBoxGroup(name="Select", value=_what_list, options=_what_list, inline=False)
        self.plot_relax_btn = pnw.Button(name="Plot relaxation", button_type="primary")

        self.appname = pnw.Select(name="Viewer", value="ovito", options=["ovito", "mayavi", "vtk"])
        self.to_unit_cell = pnw.Checkbox(name="To unit cell")
        self.view_relax_btn = pnw.Button(name="View relaxation", button_type="primary")

        super().__init__(**params)

    @depends_on_btn_click('plot_relax_btn')
    def on_plot_relax_btn(self):
        """
        Plot the evolution of structural parameters (lattice lengths, angles and volume)
        as well as pressure, info on forces and total energy.
        """
        col = pn.Column(sizing_mode="stretch_width"); ca = col.append
        for what in self.what_list.value:
            #ca(f"## {what}")
            ca(ply(self.hist.plotly(what, title=what, show=False)))

        return col

    @depends_on_btn_click('view_relax_btn')
    def on_view_relax_btn(self):
        """
        Visalize the structural relaxation with an external application.
        """
        return self.hist.visualize(appname=self.appname.value, to_unit_cell=self.to_unit_cell.value)

    def get_panel(self, as_dict=False, **kwargs):
        """Return tabs with widgets to interact with the HIST.nc file."""
        d = {}

        d["Summary"] = self.get_summary_view_for_abiobj(self.hist)
        d["Plot"] = pn.Row(
                self.pws_col(["## Plot Options", "what_list", "plot_relax_btn"]),
                self.on_plot_relax_btn
        )

        if not self.has_remote_server:
            # As we don't have visualizers that can work in remote server mode,
            # this tab should not be created.
            d["Visualize"] = pn.Row(
                    pn.Column(self.appname, self.to_unit_cell, self.view_relax_btn),
                    self.on_view_relax_btn
            )

        if as_dict: return d

        return self.get_template_from_tabs(d, template=kwargs.get("template", None))
