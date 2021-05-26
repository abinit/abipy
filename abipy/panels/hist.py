""""Panels for HIST files."""
import param
import panel as pn
import panel.widgets as pnw
import bokeh.models.widgets as bkw

from abipy.panels.core import AbipyParameterized, mpl, ply, depends_on_btn_click


_what_list = ["pressure", "forces", "energy", "abc", "angles", "volume"]


class HistFilePanel(AbipyParameterized):
    """
    Panel with widgets to interact with a |HistFile|.
    """

    def __init__(self, hist, **params):
        self.hist = hist
        self.what_list = pnw.CheckBoxGroup(name="Select", value=_what_list, options=_what_list, inline=False)
        self.plot_relax_btn = pnw.Button(name="Show relaxation", button_type="primary")

        self.appname = pnw.Select(name="Viewer", value="ovito", options=["ovito", "mayavi", "vtk"])
        self.to_unit_cell = pnw.Checkbox(name="To unit cell")
        self.view_relax_btn = pnw.Button(name="View relaxation", button_type="primary")

        super().__init__(**params)

    def get_plot_relax_widgets(self):
        """Widgets to visualize the structure relaxation."""
        return pn.Column(self.what_list, self.plot_relax_btn, sizing_mode="stretch_width")

    @depends_on_btn_click('plot_relax_btn')
    def on_plot_relax_btn(self):
        """
        Plot the evolution of structural parameters (lattice lengths, angles and volume)
        as well as pressure, info on forces and total energy.
        """
        num_plots, nrows, ncols = len(self.what_list.value), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots // ncols) + (num_plots % ncols)

        box = pn.GridBox(nrows=nrows, ncols=ncols, sizing_mode=self.sizing_mode.value) #'scale_width')
        for i, what in enumerate(self.what_list.value):
            irow, icol = divmod(i, ncols)
            box.append(mpl(self.hist.plot(what, title=what, **self.mpl_kwargs)))

        return box
        #return pn.Column(box, box.controls(jslink=True))

    @depends_on_btn_click('view_relax_btn')
    def on_view_relax_btn(self):
        return self.hist.visualize(appname=self.appname.value, to_unit_cell=self.to_unit_cell.value)

    def get_panel(self, as_dict=False, **kwargs):
        """Return tabs with widgets to interact with the HIST.nc file."""
        d = {}

        d["Summary"] = pn.Row(bkw.PreText(text=self.hist.to_string(verbose=self.verbose),
                              sizing_mode="scale_both"))

        d["Relaxation"] = pn.Row(self.get_plot_relax_widgets(), self.on_plot_relax_btn)

        d["Visualize"] = pn.Row(pn.Column(self.appname, self.to_unit_cell, self.view_relax_btn),
                                 self.on_view_relax_btn)

        if as_dict: return d
        tabs = pn.Tabs(*d.items())
        return self.get_template_from_tabs(tabs, template=kwargs.get("template", None))
