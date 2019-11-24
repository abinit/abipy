""""Panels for HIST files."""
import param
import panel as pn
import bokeh.models.widgets as bkw

from abipy.panels.core import AbipyParameterized


_what_list = ["pressure", "forces", "energy", "abc", "angles", "volume"]


class HistFilePanel(AbipyParameterized):
    """
    """

    plot_relax_btn = pn.widgets.Button(name="Show relaxation", button_type='primary')

    what_list = pn.widgets.CheckBoxGroup(name='Select', value=_what_list, options=_what_list, inline=False)

    def __init__(self, hist, **params):
        super().__init__(**params)
        self.hist = hist

    def get_plot_relax_widgets(self):
        """Widgets to visualize relaxation."""
        return pn.Column(self.what_list, self.plot_relax_btn)

    @param.depends('plot_relax_btn.clicks')
    def on_plot_relax_btn(self):
        """
        Plot the evolution of structural parameters (lattice lengths, angles and volume)
        as well as pressure, info on forces and total energy.
        """
        if self.plot_relax_btn.clicks == 0: return

        num_plots, nrows, ncols = len(self.what_list.value), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots // ncols) + (num_plots % ncols)

        box = pn.GridBox(nrows=nrows, ncols=ncols, sizing_mode='scale_width')
        for i, what in enumerate(self.what_list.value):
            irow, icol = divmod(i, ncols)
            box.append(self._mp(self.hist.plot(what, title=what, **self.fig_kwargs)))

        return box

    def get_panel(self):
        """Return tabs with widgets to interact with the DDB file."""
        tabs = pn.Tabs()
        tabs.append(("Summary", pn.Row(bkw.PreText(text=self.hist.to_string(verbose=self.verbose),
                     sizing_mode="scale_both"))))
        tabs.append(("Relax", pn.Row(self.get_plot_relax_widgets(), self.on_plot_relax_btn)))

        return tabs
