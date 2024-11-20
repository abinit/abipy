"""Panels to interact with A2F files."""
from __future__ import annotations

#import param
import panel as pn
import panel.widgets as pnw

from .core import (PanelWithElectronBands, ply, mpl, dfc, depends_on_btn_click)


class A2fFilePanel(PanelWithElectronBands):
    """
    Panel with widgets to interact with a |GsrFile|.
    """
    def __init__(self, ncfile, **params):
        PanelWithElectronBands.__init__(self, ebands=ncfile.ebands, **params)
        self.ncfile = ncfile

        self.a2f_view_btn = pnw.Button(name="Plot a2F", button_type='primary')

    def get_a2f_view(self) -> pn.Row:
        """
        Return Row with widgets to visualize the structure.
        """
        return pn.Row(
            self.pws_col(["## Visualize a2F",
                          "a2f_view_btn",
                          ]),
            self.on_view_a2f,
        )

    @depends_on_btn_click('a2f_view_btn')
    def on_view_a2f(self):
        """Visualize a2f function."""
        col = pn.Column(sizing_mode='stretch_width'); ca = col.append

        ncfile = self.ncfile

        a2f = ncfile.a2f_qcoarse
        ca(f"## A2F")
        fig = a2f.plot_with_lambda(units="eV", **self.mpl_kwargs)
        ca(mpl(fig))

        ca(f"## A2F + nu")
        fig = a2f.plot_nuterms(units="eV", ax_mat=None, with_lambda=True, fontsize=12,
                               xlims=None, ylims=None, label=None, **self.mpl_kwargs)
        ca(mpl(fig))

        for what in ("gamma", "lambda"):
            ca(f"## Phbands with {what}")
            fig = ncfile.plot(what=what, units="eV", scale=None, alpha=0.6, ylims=None,
                              ax=None, colormap="jet", **self.mpl_kwargs)
            ca(mpl(fig))

        return col

    def get_panel(self, as_dict=False, **kwargs):
        """
        Return tabs with widgets to interact with the A2F file.
        """
        d = {}

        #d["Summary"] = self.get_summary_view_for_abiobj(self.gsr)
        d["a2F"] = self.get_a2f_view()
        d["e-Bands"] = self.get_plot_ebands_view()

        #tabs_dict = self.ncfile.phbands.get_panel(as_dict=True)
        #for k, v in tabs_dict.items():
        #    d[k] = v

        kpoints = self.ncfile.ebands.kpoints
        if kpoints.is_ibz:
            # Add DOS tab but only if k-sampling.
            d["e-DOS"] = self.get_plot_edos_view()
            d["SKW"] = self.get_skw_view()

            if not self.ncfile.ebands.isnot_ibz_sampling():
                d["Ifermi"] = self.get_ifermi_view()
                #d["fsviewer"] = self.get_fsviewer_view()

        d["Structure"] = self.get_structure_view()
        d["NcFile"] = self.ncfile.get_ncfile_view()

        # TODO
        #d["Global"] = pn.Row(
        #    pn.Column("# Global options",
        #              *self.pws("units", "mpi_procs", "verbose"),
        #              ),
        #    self.get_software_stack())
        #))

        if as_dict: return d
        return self.get_template_from_tabs(d, template=kwargs.get("template", None))
