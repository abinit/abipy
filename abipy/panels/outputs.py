"""Panels for interacting with output files in text format."""

import param
import panel as pn
#import panel.widgets as pnw
import bokeh.models.widgets as bkw
from abipy.panels.core import AbipyParameterized


class AbinitOutputFilePanel(AbipyParameterized):
    """
    Panel with widgets to interact with an Abinit output file.
    """
    def __init__(self, outfile, **params):
        super().__init__(**params)
        self.outfile = outfile

    def _get_gridbox(self, what):
        """Return GridBox with matplotlib for GS/DFPT SCF cycles."""
        if what == "GS":
            cycles = self.outfile.get_all_gs_scf_cycles()
        elif what == "DFPT":
            cycles = self.outfile.get_all_d2de_scf_cycles()
        else:
            raise ValueError("Invalid value for what: %s" % what)

        if not cycles: return None

        num_plots, nrows, ncols = len(cycles), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots // ncols) + (num_plots % ncols)

        box = pn.GridBox(nrows=nrows, ncols=ncols) #, sizing_mode='scale_both')
        for icycle, cycle in enumerate(cycles):
            box.append(self._mp(cycle.plot(title="%s cycle #%d" % (what, icycle), **self.fig_kwargs)))

        return box

    def get_panel(self):
        """Return tabs with widgets to interact with the Abinit output file."""
        tabs = pn.Tabs(); app = tabs.append
        app(("Summary", pn.Row(bkw.PreText(text=self.outfile.to_string(verbose=self.verbose),
                               sizing_mode="scale_both"))))
        df = self.outfile.get_dims_spginfo_dataframe().transpose()
        df.index.name = "Dataset"
        app(("Dims", self._df(df)))

        # Add tabs with plots for the GS/DFPT SCF cycles.
        for what in ("GS", "DFPT"):
            box = self._get_gridbox(what)
            if box is not None:
                app(("%s Cycles" % what, box))

        #timer = self.get_timer()
        #timer.plot_all(**self.fig_kwargs)

        return tabs
