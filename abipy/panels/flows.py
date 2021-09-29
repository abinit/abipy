""""Panels to interact with AbiPy flows."""
from __future__ import annotations

import param
import panel as pn
import panel.widgets as pnw
#import bokeh.models.widgets as bkw

from panel.viewable import Viewer
from abipy.panels.core import mpl, ply, dfc, depends_on_btn_click
from abipy.panels.nodes import NodeParameterized
from abipy.flowtk.tasks import AbinitTask
from abipy.flowtk.flows import Flow


class WorkTaskSelector(Viewer):

    task = param.ClassSelector(class_=AbinitTask, doc="Task object")

    def __init__(self, flow: Flow, **params):
        self._wstr2work = {f"w{i} ({work.__class__.__name__}, len: {len(work)})": work
                           for (i, work) in enumerate(flow.works)}
        options = list(self._wstr2work.keys())
        self.work_select = pnw.Select(name="Select a Work", value=options[0], options=options)

        options = [f"t{i} ({task.__class__.__name__}, {str(task.status)})"
                   for (i, task) in enumerate(flow[0])]
        self.task_select = pnw.Select(name="Select a Task in the Work", value=options[0], options=options)

        super().__init__(**params)

        self.layout = pn.Column(self.work_select, self.task_select)
        self.sync_widgets()

    def __panel__(self):
        return self.layout

    @pn.depends('work_select.value', watch=True)
    def update_work(self):
        self.work = self._wstr2work[self.work_select.value]
        self.task_select.options = [f"t{i} ({task.__class__.__name__}, {str(task.status)})"
                                    for (i, task) in enumerate(self.work)]
        self.task = self.work[0]

    @pn.depends('task_select.value', watch=True)
    def sync_widgets(self):
        self.work = self._wstr2work[self.work_select.value]
        task_idx = int(self.task_select.value[1:].split()[0])
        self.task = self.work[task_idx]


class FlowPanel(NodeParameterized):
    """
    Provides widgets and callbacks to interact with an AbiPy Flow.
    """

    def __init__(self, flow: Flow, **params):
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
        #d["Work"] = self.get_work_view()
        #d["Structures"] = pn.Row(pn.Column(self.structures_io_checkbox, self.structures_btn), self.on_structures_btn)
        ###ws = pn.Column(self.ebands_plotter_mode, self.ebands_ksamp_checkbox, self.ebands_df_checkbox, self.ebands_plotter_btn)
        ###d["Ebands"] = pn.Row(ws, self.on_ebands_btn)
        #d["Browse"] = self.get_workdir_view()

        if as_dict: return d

        return self.get_template_from_tabs(d, template=kwargs.get("template", None), closable=False)


class JsPane(pn.pane.HTML):
    """
    Based on: https://discourse.holoviz.org/t/how-to-make-a-dynamic-link-in-panel/2137
    """
    def __init__(self):
        super().__init__(width=0, height=0, margin=0, sizing_mode="fixed")

    def execute_js(self, script: str):
        script = f'<script type="text/javascript">{script}</script>'
        self.object = script
        self.object = ""


class FlowMultiPageApp():

    def __init__(self, flow: Flow, template, spectator_mode=True, **kwargs):

        if spectator_mode:
            # We are in read-only mode so we have to disable signals to avoid side effects and callbacks
            flow.set_spectator_mode()

        self.flow = flow
        self.template = template

        self.wt_selector = WorkTaskSelector(flow)
        goto_work_btn = pnw.Button(name="Go to Work", button_type='primary')
        goto_work_btn.on_click(self.on_goto_work_bnt)
        goto_task_btn = pnw.Button(name="Go to Task", button_type='primary')
        goto_task_btn.on_click(self.on_goto_task_bnt)
        self.new_tab = pnw.Checkbox(value=True, name="Open in new Tab")
        self.js_panel = JsPane()

        self.sidebar = pn.WidgetBox(
            self.wt_selector,
            self.new_tab,
            pn.layout.Divider(),
            goto_work_btn,
            goto_task_btn,
            self.js_panel,
        )

        def handle_home():
            app = FlowPanel(self.flow).get_panel(template=template)
            if hasattr(app, "sidebar"):
                app.sidebar.append(self.sidebar)
                #app.header.append(self.sidebar)

            return app

        # url --> handler
        self.routes = {
            "/": handle_home,
            r"/w\d+/?": self.handle_wt,
            r"/w\d+/t\d+": self.handle_wt,
        }

    def on_goto_work_bnt(self, event):
        work = self.wt_selector.work
        url = "/w%d" % work.pos
        code = f"window.open('{url}')" if self.new_tab.value else f"window.location.href='{url}'"
        self.js_panel.execute_js(code)

    def on_goto_task_bnt(self, event):
        task = self.wt_selector.task
        url = "/w%d/t%d" % (task.pos[0], task.pos[1])
        code = f"window.open('{url}')" if self.new_tab.value else f"window.location.href='{url}'"
        self.js_panel.execute_js(code)

    def handle_wt(self):
        # URL example: /w1/t5/w\d+/t\d+
        #print("in handle_wt with pn.state.app_url:", pn.state.app_url)
        tokens = pn.state.app_url.split("/")
        work_idx = int(tokens[1][1:])
        task_idx = None
        if tokens[2].startswith("t"):
            task_idx = int(tokens[2][1:])

        #print("got request with work_idx:", work_idx, "task_idx:", task_idx)
        #from abipy.panels.core import abipanel
        #abipanel()

        if task_idx is None:
            work = self.flow[work_idx]
            app = work.get_panel(template=self.template, title=repr(work))
        else:
            task = self.flow[work_idx][task_idx]
            app = task.get_panel(template=self.template, title=repr(task))

        if hasattr(app, "sidebar"):
            app.sidebar.append(self.sidebar)
            #app.header.append(self.sidebar)

        #app.title = title

        return app

    def serve(self, **serve_kwargs):
        return pn.serve(self.routes, **serve_kwargs)