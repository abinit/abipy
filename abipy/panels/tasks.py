""""Panels to interact with the AbiPy tasks."""
import param
import panel as pn
import panel.widgets as pnw
import bokeh.models.widgets as bkw

from io import StringIO
from abipy.panels.core import AbipyParameterized, mpl, ply, dfc, depends_on_btn_click
from abipy.panels.nodes import NodeParameterized


class AbinitTaskPanel(NodeParameterized):
    """
    Provides widgets to interact with an AbiPy Task.
    """

    def __init__(self, task, **params):
        NodeParameterized.__init__(self, node=task, **params)

        self.task = task

        #self.structures_btn = pnw.Button(name="Show Structures", button_type='primary')
        #self.structures_io_checkbox = pnw.CheckBoxGroup(
        #    name='Input/Output Structure', value=['output'], options=['input', 'output'], inline=Tru

    def on_input_files_tab(self):
        """
        Show the input files of the task: submission script and Abinit input file.
        """
        file = self.task.job_file
        text = file.read() if file.exists else "Cannot find job_file!"
        job_file = pn.pane.Markdown(f"```shell\n{text}\n```")

        return pn.Column(
            "## Submission script:",
            job_file,
            pn.layout.Divider(),
            "## Input file:",
            self.html_with_clipboard_btn(self.task.input),
            sizing_mode="stretch_width",
        )

    def on_err_files_tab(self):
        """
        Show the error files of the task
        """
        col = pn.Column(sizing_mode="stretch_width"); cext = col.extend

        count = 0
        for fname in ["stderr_file", "mpiabort_file"]: # "log_file",
            file = getattr(self.task, fname)
            if file.exists:
                text = file.read().strip()
                if text:
                    cext([f"## {fname}", pn.pane.Markdown(f"```shell\n{text}\n```"), pn.layout.Divider()])
                    count += 1

        return col

    @depends_on_btn_click("structures_btn")
    def on_structures_btn(self):
        what = ""
        if "input" in self.structures_io_checkbox.value: what += "i"
        if "output" in self.structures_io_checkbox.value: what += "o"
        dfs = self.flow.compare_structures(nids=self.nids,
                                           what=what,
                                           verbose=self.verbose, with_spglib=False, printout=False,
                                           with_colors=False)

        return pn.Row(dfc(dfs.lattice), sizing_mode="scale_width")

    def get_panel(self, as_dict=False, **kwargs):
        """Return tabs with widgets to interact with the flow."""

        d =  super().get_panel(as_dict=True)
        #return d
        #d =  {}

        #d["Summary"] = pn.Row(bkw.PreText(text=str(self.task))) # .to_string(verbose=self.verbose)))
        ##d["Structures"] = pn.Row(pn.Column(self.structures_io_checkbox, self.structures_btn), self
        ###ws = pn.Column(self.ebands_plotter_mode, self.ebands_ksamp_checkbox, self.ebands_df_check
        ###d["Ebands"] = pn.Row(ws, self.on_ebands_btn)

        # This stuff is computed lazyly when the tab is activated.
        d["Input"] = pn.param.ParamMethod(self.on_input_files_tab, lazy=True)
        d["Err"] = pn.param.ParamMethod(self.on_err_files_tab, lazy=True)
        #d["Outdir"] = self.get_dirview("outdir")

        if as_dict: return d

        return self.get_template_from_tabs(d, template=kwargs.get("template", None),
                                           closable=False, dynamic=True)
