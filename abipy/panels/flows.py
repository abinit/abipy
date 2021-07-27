""""Panels to interact with AbiPy flows."""
import param
import panel as pn
import panel.widgets as pnw
import bokeh.models.widgets as bkw

from panel.viewable import Viewer
from abipy.panels.core import AbipyParameterized, mpl, ply, dfc, depends_on_btn_click
from abipy.panels.nodes import NodeParameterized
from abipy import flowtk


class TaskSelector(Viewer):

    value = param.ClassSelector(class_=flowtk.AbinitTask, doc="Task object")

    def __init__(self, flow, **params):
        self._wstr2work = {f"w{i}": work for (i, work) in enumerate(flow.works)}

        self._work_select = pnw.Select(name="Select a Work",
                                       value="w0", options=list(self._wstr2work.keys()))
        self._task_select = pnw.Select(name="Select a Task of the Work",
                                       value="t0", options=[f"t{i}" for i in range(len(flow[0]))])

        super().__init__(**params)
        #self._layout = pn.Row(self._work_select, self._task_select)
        self._layout = pn.WidgetBox(self._work_select, self._task_select)
        self._sync_widgets()

    def __panel__(self):
        return self._layout

    @param.depends('_work_select.value', watch=True)
    def _sync_widgets(self):
        work = self._wstr2work[self._work_select.value]
        self._task_select.options = [f"t{i}" for i in range(len(work))]
        self.value = work[0]
        print("_sync_widgets", self.value)

    @param.depends('_task_select.value', watch=True)
    def _update_task(self):
        work = self._wstr2work[self._work_select.value]
        task_idx = int(self._task_select.value[1:])
        print("task_idx:", task_idx)
        self.value = work[task_idx]
        print("_update_task", self.value)


class FlowPanel(NodeParameterized):
    """
    Provides widgets and callbacks to interact with an AbiPy Flow.
    """

    #param.FileSelector()

    def __init__(self, flow, **params):
        NodeParameterized.__init__(self, node=flow, **params)

        self.structures_btn = pnw.Button(name="Show Structures", button_type='primary')
        self.structures_io_checkbox = pnw.CheckBoxGroup(
            name='Input/Output Structure', value=['output'], options=['input', 'output'], inline=True)

        self.task_selector = TaskSelector(flow)
        self.task_btn = pnw.Button(name="Show Task", button_type='primary')

    def get_task_view(self):
        return pn.Column(
            "## CLick the button to show the status of the task",
            #self.pws_row(["task_selector", "task_btn"]),
            self.task_selector,
            self.task_btn,
            pn.layout.Divider(),
            self.on_task_btn,
            sizing_mode='stretch_width',
            )

    @depends_on_btn_click("task_btn")
    def on_task_btn(self):
        """
        Return panel associated to the selected task.
        """
        task = self.task_selector.value
        return pn.Column(
                f"## {repr(task)}",
                task.get_panel(),
        )

    @depends_on_btn_click("structures_btn")
    def on_structures_btn(self):
        what = ""
        if "input" in self.structures_io_checkbox.value: what += "i"
        if "output" in self.structures_io_checkbox.value: what += "o"
        dfs = self.flow.compare_structures(nids=None, # select_nids(flow, options),
                                           what=what,
                                           verbose=self.verbose, with_spglib=False, printout=False,
                                           with_colors=False)

        return pn.Row(dfc(dfs.lattice), sizing_mode="scale_width")

    def get_panel(self, as_dict=False, **kwargs):
        """Return tabs with widgets to interact with the flow."""

        d = super().get_panel(as_dict=True)

        #row = pn.Row(bkw.PreText(text=self.ddb.to_string(verbose=self.verbose), sizing_mode="scale_both"))
        d["Task"] = self.get_task_view()
        #d["Task"] = self.get_work_view()
        #d["Task"] = self.get_work_view()
        #d["Structures"] = pn.Row(pn.Column(self.structures_io_checkbox, self.structures_btn), self.on_structures_btn)
        ###ws = pn.Column(self.ebands_plotter_mode, self.ebands_ksamp_checkbox, self.ebands_df_checkbox, self.ebands_plotter_btn)
        ###d["Ebands"] = pn.Row(ws, self.on_ebands_btn)
        d["Browse"] = self.get_workdir_view()

        if as_dict: return d

        return self.get_template_from_tabs(d, template=kwargs.get("template", None), closable=False)
