"""Panels to interact with the AbiPy tasks."""
from __future__ import annotations

import panel as pn
#import panel.widgets as pnw

from abipy.panels.viewers import AceViewer
#from abipy.panels.core import mpl, ply, dfc, depends_on_btn_click, Loading
from abipy.panels.nodes import NodeParameterized
from abipy.flowtk.tasks import AbinitTask


class TaskPanel(NodeParameterized):
    """
    Provides widgets to interact with an AbiPy Task.
    """

    def __init__(self, task: AbinitTask, **params):
        NodeParameterized.__init__(self, node=task, **params)
        self.task = task

        #self.structures_btn = pnw.Button(name="Show Structures", button_type='primary')
        #self.structures_io_checkbox = pnw.CheckBoxGroup(
        #    name='Input/Output Structure', value=['output'], options=['input', 'output'], inline=Tru

    def get_inputs_view(self):
        """
        Show the input files of the task: input file, submission script and TaskManager
        """
        file = self.task.job_file
        text = file.read() if file.exists else "Cannot find job_file!"
        job_file = pn.pane.Markdown(f"```shell\n{text}\n```")

        from .viewers import JSONViewer
        json_view = JSONViewer(self.task.manager.as_dict())

        def card(title, *items, collapsed=True):
            return pn.Card(*items,
                           title=title,
                           collapsed=collapsed,
                           sizing_mode="stretch_width",
                           header_color="blue",
                           #header_background="blue",
                           )

        return pn.Column(
            f"## Input files of `{repr(self.task)}`",
            card("Input file", self.html_with_clipboard_btn(self.task.input), collapsed=False),
            pn.layout.Divider(),
            card("Submission script", job_file),
            pn.layout.Divider(),
            card("TaskManager", json_view),
            pn.layout.Divider(),
            sizing_mode="stretch_width",
        )

    def get_errs_view(self):
        """
        Show the error files of the task
        Return None if no error is found so that we don't show this view in the GUI.
        """
        col = pn.Column(sizing_mode="stretch_width"); cext = col.extend

        count = 0
        for fname in ("stderr_file", "mpiabort_file", "qerr_file", "qout_file"):
            file = getattr(self.task, fname)
            if file.exists:
                text = file.read().strip()
                if text:
                    cext([f"## {fname}",
                         pn.pane.Markdown(f"```shell\n{text}\n```"),
                         pn.layout.Divider()
                         ])
                    count += 1

        return col if count > 0 else None

    def get_main_text_outs_view(self):
        """
        Show the main text output files of the task.
        """
        col = pn.Column(f"## Main output and log file of `{repr(self.task)}`",
                        sizing_mode="stretch_width")

        for fname in ("output_file", "log_file"):
            file = getattr(self.task, fname)
            col.append(AceViewer(file.path))

        return col

    #@depends_on_btn_click("structures_btn")
    #def on_structures_btn(self):
    #    what = ""
    #    if "input" in self.structures_io_checkbox.value: what += "i"
    #    if "output" in self.structures_io_checkbox.value: what += "o"
    #    dfs = self.flow.compare_structures(nids=self.nids,
    #                                       what=what,
    #                                       verbose=self.verbose, with_spglib=False, printout=False,
    #                                       with_colors=False)

    #    return pn.Row(dfc(dfs.lattice), sizing_mode="scale_width")

    def get_panel(self, as_dict=False, **kwargs):
        """Return tabs with widgets to interact with the flow."""
        d = dict()

        # This stuff is computed lazyly when the tab is activated.
        d["Inputs"] = self.get_inputs_view()
        d["Output"] = self.get_main_text_outs_view()
        view = self.get_errs_view()
        if view is not None: d["ErrFiles"] = view

        super_d =  super().get_panel(as_dict=True)
        d.update(super_d)
        ##d["Structures"] = pn.Row(pn.Column(self.structures_io_checkbox, self.structures_btn), self

        if as_dict: return d

        return self.get_template_from_tabs(d, template=kwargs.get("template", None))