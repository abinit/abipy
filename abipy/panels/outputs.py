"""Panels for interacting with output files in text format."""

import panel as pn
import panel.widgets as pnw
import bokeh.models.widgets as bkw

from abipy.panels.core import AbipyParameterized, Loading, mpl, ply, dfc


class AbinitOutputFilePanel(AbipyParameterized):
    """
    Panel with widgets to interact with the main Abinit output file.
    """
    def __init__(self, outfile, **params):
        super().__init__(**params)
        self.outfile = outfile

    def get_cycles_view(self, what):
        """
        Return GridBox with matplotlib plots for the GS/DFPT SCF cycles
        or None if no cycle is found.
        """
        if what == "GS_SCF":
            cycles = self.outfile.get_all_gs_scf_cycles()
        elif what == "DFPT_D2DE":
            cycles = self.outfile.get_all_d2de_scf_cycles()
        else:
            raise ValueError(f"Invalid value for what: `{what}`")

        if not cycles: return None

        num_plots, nrows, ncols = len(cycles), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots // ncols) + (num_plots % ncols)

        box = pn.GridBox(nrows=nrows, ncols=ncols, sizing_mode='stretch_width')
        #box = pn.Column(sizing_mode="stretch_width")

        for icycle, cycle in enumerate(cycles):
            #box.append(mpl(cycle.plot(title="%s cycle #%d" % (what, icycle), **self.mpl_kwargs)))
            f = ply(cycle.plotly(title="%s cycle #%d" % (what, icycle + 1), show=False))
            box.append(f)

        return box

    def get_outfile_view(self):
        col = pn.Column(sizing_mode="stretch_width")
        ca = col.append; cext = col.extend

        filepath = self.outfile.filepath
        with open(filepath, "rt") as fh:
            text = fh.read()

        ace = pnw.Ace(value=text, language='text', readonly=True,
                      sizing_mode='stretch_width', height=1200)
                      #sizing_mode='stretch_width', width=900)
        cext([f"## Output <small>{filepath}</small>", ace, pn.layout.Divider()])

        return col

    def get_panel(self, as_dict=False, **kwargs):
        """Return tabs with widgets to interact with the Abinit output file."""
        d = {}

        d["Summary"] = self.get_summary_view_for_abiobj(self.outfile)
        #d["Output"] = self.get_outfile_view()

        df = self.outfile.get_dims_spginfo_dataframe()
        df.index.name = "Dataset"
        d["Dims"] = dfc(df, transpose=True)

        # Add tabs with plots for the GS/DFPT SCF cycles.
        for what in ("GS_SCF", "DFPT_D2DE"):
            box = self.get_cycles_view(what)
            if box is not None:
                d[f"{what} cycles"] = box

        #timer = self.get_timer()
        #timer.plot_all(**self.mpl_kwargs)

        if as_dict: return d

        return self.get_template_from_tabs(d, template=kwargs.get("template", None))


class AbinitOutputFilePanelWithFileInput(AbipyParameterized):

    info_str = """
This application allows users to analyze the Abinit main output file
"""

    def __init__(self, **params):
        super().__init__(**params)

        help_md = pn.pane.Markdown(f"""
## Description

{self.info_str}
""")

        self.main_area = pn.Column(help_md, sizing_mode="stretch_width")
        self.abifile = None

        self.file_input = pnw.FileInput(height=60, css_classes=["pnx-file-upload-area"])
        self.file_input.param.watch(self.on_file_input, "value")

    def on_file_input(self, event):

        with Loading(self.main_area):
            new_abifile = self.get_abifile_from_file_input(self.file_input)

            if self.abifile is not None:
                self.abifile.remove()

            self.abifile = new_abifile
            self.main_area.objects = [self.abifile.get_panel()]

    def get_panel(self, **kwargs):

        col = pn.Column("## Upload (or drag & drop) an *.abo* file (main ABINIT output file):",
                        self.get_fileinput_section(self.file_input),
                        sizing_mode="stretch_width")

        main = pn.Column(col, self.main_area, sizing_mode="stretch_width")

        cls, cls_kwds = self.get_abinit_template_cls_kwds()

        return cls(main=main, title="Abo Analyzer", **cls_kwds)
