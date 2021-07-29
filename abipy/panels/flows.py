""""Panels to interact with AbiPy flows."""
import param
import panel as pn
import panel.widgets as pnw
import bokeh.models.widgets as bkw

from panel.viewable import Viewer
from abipy.panels.core import mpl, ply, dfc, depends_on_btn_click
from abipy.panels.nodes import NodeParameterized
from abipy import flowtk


class WorkTaskSelector(Viewer):

    task = param.ClassSelector(class_=flowtk.AbinitTask, doc="Task object")

    def __init__(self, flow, **params):
        self._wstr2work = {f"w{i} ({work.__class__.__name__})": work for (i, work) in enumerate(flow.works)}
        options = list(self._wstr2work.keys())
        self.work_select = pnw.Select(name="Select a Work", value=options[0], options=options)

        options = [f"t{i} ({task.__class__.__name__})" for (i, task) in enumerate(flow[0])]
        self.task_select = pnw.Select(name="Select a Task in the Work", value=options[0], options=options)

        super().__init__(**params)

        self.layout = pn.Column(self.work_select, self.task_select)
        #self.outarea = pn.Column("## Hello", sizing_mode="stretch_width")
        #self.layout = pn.Row(pn.WidgetBox(self.work_select, self.task_select), self.outarea)
        self.sync_widgets()

    def __panel__(self):
        return self.layout

    @pn.depends('work_select.value', watch=True)
    def update_work(self):
        work = self._wstr2work[self.work_select.value]
        self.task_select.options = [f"t{i} ({task.__class__.__name__})" for (i, task) in enumerate(work)]
        self.task = work[0]

    @pn.depends('task_select.value', watch=True)
    def sync_widgets(self):
        work = self._wstr2work[self.work_select.value]
        task_idx = int(self.task_select.value[1:].split()[0])
        self.task = work[task_idx]
        #print(repr(self.task))
        #self.outarea.objects[0] = repr(self.task)
        #self.outarea.append(pn.pane.Markdown(repr(self.task)))



class FlowPanel(NodeParameterized):
    """
    Provides widgets and callbacks to interact with an AbiPy Flow.
    """

    def __init__(self, flow, **params):
        NodeParameterized.__init__(self, node=flow, **params)

        self.structures_btn = pnw.Button(name="Show Structures", button_type='primary')
        self.structures_io_checkbox = pnw.CheckBoxGroup(
            name='Input/Output Structure', value=['output'], options=['input', 'output'], inline=True)

        self.wt_selector = WorkTaskSelector(flow)
        self.task_btn = pnw.Button(name="Analyze Task", button_type='primary')

    def get_task_view(self):

        wbox = pn.WidgetBox

        return pn.Column(
            wbox("## Select Work and Task",
                 self.wt_selector,
                 self.task_btn,
            ),
            pn.layout.Divider(),
            self.on_task_btn,
            sizing_mode='stretch_width',
        )

    @depends_on_btn_click("task_btn")
    def on_task_btn(self):
        """
        Return panel associated to the selected task.
        """
        task = self.wt_selector.task
        return pn.Column(
            f"## {repr(task)}",
            task.get_panel(),
            sizing_mode="stretch_width",
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
        #d["Browse"] = self.get_workdir_view()

        if as_dict: return d

        return self.get_template_from_tabs(d, template=kwargs.get("template", None), closable=False)
