"""Panels for interacting with output files in text format."""

import panel as pn
import panel.widgets as pnw
import bokeh.models.widgets as bkw

from abipy.panels.core import AbipyParameterized, Loading, mpl, ply, dfc


class AbinitOutputFilePanel(AbipyParameterized):
    """
    Panel with widgets to interact with the Abinit main output file.
    """
    def __init__(self, outfile, **params):
        super().__init__(**params)
        self.outfile = outfile

    def _get_gridbox(self, what):
        """Return GridBox with matplotlib for the GS/DFPT SCF cycles."""
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
            box.append(mpl(cycle.plot(title="%s cycle #%d" % (what, icycle), **self.mpl_kwargs)))

        return box

    def get_outfile_view(self):
        col = pn.Column(sizing_mode="stretch_width")
        ca = col.append; cext = col.extend

        filepath = self.outfile.filepath
        with open(filepath) as fh:
            text = fh.read()

        ace = pnw.Ace(value=text, language='text', readonly=True,
                      sizing_mode='stretch_width', height=1200)
                      #sizing_mode='stretch_width', width=900)
        cext([f"## Output <small>{filepath}</small>", ace, pn.layout.Divider()])

        return col

    def get_panel(self, as_dict=False, **kwargs):
        """Return tabs with widgets to interact with the Abinit output file."""
        d = {}

        d["Summary"] = pn.Row(
            bkw.PreText(text=self.outfile.to_string(verbose=self.verbose), sizing_mode="scale_both")
        )

        d["Output"] = self.get_outfile_view()

        df = self.outfile.get_dims_spginfo_dataframe().transpose()
        df.index.name = "Dataset"
        d["Dims"] = dfc(df)

        # Add tabs with plots for the GS/DFPT SCF cycles.
        for what in ("GS", "DFPT"):
            box = self._get_gridbox(what)
            if box is not None:
                d["%s cycles" % what] = box

        #timer = self.get_timer()
        #timer.plot_all(**self.mpl_kwargs)

        if as_dict: return d

        tabs = pn.Tabs(*d.items())
        return self.get_template_from_tabs(tabs, template=kwargs.get("template", None))


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

        cls, kwds = self.get_abinit_template_cls_kwds()

        return cls(main=main, title="Abo Analyzer", **kwds)
