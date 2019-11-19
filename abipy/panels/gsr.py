"""Panels for GSR files."""
import panel as pn
import bokeh.models.widgets as bw

from .core import PanelWithElectronBands


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
        tabs.append(("Summary", pn.Row(bw.PreText(text=self.gsr.to_string(verbose=self.verbose),
                     sizing_mode="scale_both"))))
        tabs.append(("e-Bands", pn.Row(self.get_plot_ebands_widgets(), self.on_plot_ebands_btn)))
        # Add DOS tab only if k-sampling.
        if self.gsr.ebands.kpoints.is_ibz:
            tabs.append(("e-DOS", pn.Row(self.get_plot_edos_widgets(), self.on_plot_edos_btn)))

        return tabs
