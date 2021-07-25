""""Panels to interact with AbiPy flows."""
import param
import panel as pn
import panel.widgets as pnw
import bokeh.models.widgets as bkw

#from io import StringIO
from panel.viewable import Viewer
from abipy.panels.core import AbipyParameterized, mpl, ply, dfc, depends_on_btn_click
from abipy import flowtk


class TaskSelect(Viewer):

    value = param.ClassSelector(class_=flowtk.AbinitTask, doc="Task object")

    def __init__(self, flow, **params):
        self.wstr2work = {f"w{i}": work for (i, work) in enumerate(flow.works)}

        self.work_select = pnw.Select(value="w0", options=list(self.wstr2work.keys()))
        self.task_select = pnw.Select(value="t0", options=[f"t{i}" for i in range(len(flow[0]))])

        self.value = flow[0][0]

        super().__init__(**params)
        self._layout = pn.Row(self.work_select, self.task_select)
        #self._sync_widgets()

    def __panel__(self):
        return self._layout

    @param.depends('work_select.value', watch=True)
    def _sync_widgets(self):
        work = self.wstr2work[self.work_select.value]
        self.task_select.options = [f"t{i}" for i in range(len(work))]
        self.value = work[0]
        print("_sync_widgets", self.value)

    @param.depends('task_select.value', watch=True)
    def _update_task(self):
        work = self.wstr2work[self.work_select.value]
        task_idx = int(self.task_select.value[1:])
        print("task_idx:", task_idx)
        self.value = work[task_idx]
        print("_update_task", self.value)


class FlowPanel(AbipyParameterized):
    """
    Provides widgets and callbacks to interact with an AbiPy Flow.
    """

    #param.FileSelector()

    def __init__(self, flow, **params):
        self.flow = flow

        self.engine = pnw.Select(value="fdp",
                        options=['dot', 'neato', 'twopi', 'circo', 'fdp', 'sfdp', 'patchwork', 'osage'])
        self.dirtree = pnw.Checkbox(name='Dirtree', value=False)
        self.graphviz_btn = pnw.Button(name="Show graph", button_type='primary')

        self.status_btn = pnw.Button(name="Show status", button_type='primary')
        self.task_btn = pnw.Button(name="Show Task", button_type='primary')
        #self.work_btn = pnw.Button(name="Show Task", button_type='primary')
        self.history_btn = pnw.Button(name="Show history", button_type='primary')
        self.debug_btn = pnw.Button(name="Debug", button_type='primary')
        self.events_btn = pnw.Button(name="Events", button_type='primary')
        self.corrections_btn = pnw.Button(name="Corrections", button_type='primary')
        self.handlers_btn = pnw.Button(name="Handlers", button_type='primary')

        self.vars_text = pnw.TextInput(name='Abivars',
                placeholder='Enter list of variables separated by comma e.g. `ecut, natom`')
        self.vars_btn = pnw.Button(name="Show Variables", button_type='primary')

        self.dims_btn = pnw.Button(name="Show Dimensions", button_type='primary')

        self.structures_btn = pnw.Button(name="Show Structures", button_type='primary')
        self.structures_io_checkbox = pnw.CheckBoxGroup(
            name='Input/Output Structure', value=['output'], options=['input', 'output'], inline=True)

        # Widgets to plot ebands.
        #self.ebands_btn = pnw.Button(name="Show Ebands", button_type='primary')
        #self.ebands_plotter_mode = pnw.Select(name="Plot Mode", value="gridplot",
        #    options=["gridplot", "combiplot", "boxplot", "combiboxplot"]) # "animate",
        #self.ebands_plotter_btn = pnw.Button(name="Plot", button_type='primary')
        #self.ebands_df_checkbox = pnw.Checkbox(name='With Ebands DataFrame', value=False)
        #self.ebands_ksamp_checkbox = pnw.CheckBoxGroup(
        #    name='Input/Output Structure', value=["with_path", "with_ibz"], options=['with_path', 'with_ibz'], inline=True)

        #TODO: Implement widget for selected_nids(flow, options),
        #radio_group = pn.widgets.RadioButtonGroup(
        #   name='Radio Button Group', options=['Biology', 'Chemistry', 'Physics'], button_type='success')

        self.workdir_selector = pn.widgets.FileSelector(flow.workdir)
        #outdir_selector = pn.widgets.FileSelector('~')
        #indir_selector = pn.widgets.FileSelector('~')

        self.task_select = TaskSelect(flow)

        super().__init__(**params)

    def get_status_view(self):
        return pn.Column(
            "## CLick the button to show the status of the flow",
            self.pws_row(["status_btn", "verbose"]),
            pn.layout.Divider(),
            self.on_status_btn,
            sizing_mode='stretch_width',
            )

    @depends_on_btn_click("status_btn")
    def on_status_btn(self):
        """
        Show the status of the flow.
        """
        stream = pnw.Terminal(output="\n\n",
            height=1200, # Need this one else the terminal is not show properly
            sizing_mode='stretch_width',
        )

        self.flow.show_status(stream=stream, verbose=self.verbose)
        #return pn.Row(stream, sizing_mode="stretch_both")
        return stream

    def get_task_view(self):
        return pn.Column(
            "## CLick the button to show the status of the flow",
            self.pws_row(["task_select", "task_btn"]),
            pn.layout.Divider(),
            self.on_task_btn,
            sizing_mode='stretch_width',
            )

    @depends_on_btn_click("task_btn")
    def on_task_btn(self):
        """
        Return panel associated to the selected task.
        """
        #task = self.flow[0][0]
        task = self.task_select.value
        return task.get_panel()

    def get_history_view(self):
        return pn.Column(
            "## CLick the button to show the history of the flow",
            self.pws_row(["history_btn", "verbose"]),
            self.on_history_btn,
            pn.layout.Divider(),
            sizing_mode='stretch_width',
            )

    @depends_on_btn_click("history_btn")
    def on_history_btn(self):
        """
        Show the history of the flow.
        """
        stream = pnw.Terminal(output="\n\n",
            height=1200, # Need this one else the terminal is not show properly
            sizing_mode='stretch_width',
        )
        #flow.show_history(status=options.task_status, nids=selected_nids(flow, options),
        #                  full_history=options.full_history, metadata=options.metadata)
        self.flow.show_history(stream=stream)
        return stream

    @depends_on_btn_click("graphviz_btn")
    def on_graphviz_btn(self):
        """
        Visualize the flow with graphviz.
        """
        node = self.flow
        if self.dirtree.value:
            graph = node.get_graphviz_dirtree(engine=self.engine.value)
        else:
            graph = node.get_graphviz(engine=self.engine.value)
        return pn.Column(graph)

    def get_debug_view(self):
        return pn.Column(
            "## CLick the button to debug the flow",
            self.pws_row(["history_btn", "verbose"]),
            self.on_debug_btn,
            pn.layout.Divider(),
            sizing_mode='stretch_width',
            )

    @depends_on_btn_click("debug_btn")
    def on_debug_btn(self):
        stream = pnw.Terminal(output="\n\n",
            height=1200, # Need this one else the terminal is not show properly
            sizing_mode='stretch_width',
        )
        #flow.debug(status=options.task_status, nids=selected_nids(flow, options))
        self.flow.debug(stream=stream)
        return stream

    def get_events_view(self):
        return pn.Column(
            "## CLick the button to show the events of the flow",
            self.pws_row(["history_btn", "verbose"]),
            self.on_events_btn,
            pn.layout.Divider(),
            sizing_mode='stretch_width',
            )

    @depends_on_btn_click("events_btn")
    def on_events_btn(self):
        stream = pnw.Terminal(output="\n\n",
            height=1200, # Need this one else the terminal is not show properly
            sizing_mode='stretch_width',
        )
        self.flow.show_events(stream=stream)
        #flow.show_events(status=options.task_status, nids=selected_nids(flow, options))
        return stream

    @depends_on_btn_click("corrections_btn")
    def on_corrections_btn(self):
        stream = pnw.Terminal(output="\n\n",
            height=1200, # Need this one else the terminal is not show properly
            sizing_mode='stretch_width',
        )
        self.flow.show_corrections(stream=stream)
        return stream

    @depends_on_btn_click("handlers_btn")
    def on_handlers_btn(self):
        stream = pnw.Terminal(output="\n\n",
            height=1200, # Need this one else the terminal is not show properly
            sizing_mode='stretch_width',
        )
        #if options.doc:
        #    flowtk.autodoc_event_handlers()
        #else:
        #show_events(self, status=None, nids=None, stream=sys.stdout):
        self.flow.show_event_handlers(verbose=self.verbose, stream=stream)
        return stream

    @depends_on_btn_click("vars_btn")
    def on_vars_btn(self):
        if not self.vars_text.value: return
        varnames = [s.strip() for s in self.vars_text.value.split(",")]
        df = self.flow.compare_abivars(varnames=varnames, # nids=selected_nids(flow, options),
                                       printout=False, with_colors=False)
        return pn.Row(dfc(df))

    @depends_on_btn_click("dims_btn")
    def on_dims_btn(self):
        df = self.flow.get_dims_dataframe(# nids=selected_nids(flow, options),
                                          printout=False, with_colors=False)
        return pn.Row(dfc(df), sizing_mode="scale_width")

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

    def get_workdir_view(self):

        col = pn.Column(sizing_mode="stretch_width"); ca = col.append
        ca(self.workdir_selector)
        self.workdir_selector_btn = pnw.Button(name="Show file", button_type='primary')
        self.workdir_selector_btn.on_click(self.on_workdir_selector_btn)
        ca(self.workdir_selector_btn)
        self.workdir_selector_output_area = pn.Column(sizing_mode="stretch_width")
        ca(self.workdir_selector_output_area)

        return pn.Column(col, self.workdir_selector_output_area, sizing_mode="stretch_width")

    def on_workdir_selector_btn(self, event):
        """hello word"""
        filepaths = self.workdir_selector.value
        if not filepaths:
            objects = [pn.pane.Alert("## No file selected", altert_type="warning")]

        else:
            from abipy.abilab import abiopen
            objects = []
            for path in filepaths:
                try:
                    with abiopen(path) as abifile:
                        pn_obj = abifile.get_panel()
                except Exception:
                    pn_obj = pn.pane.Alert(str(exc), altert_type="warning")

                objects.append(pn_obj)

        self.workdir_selector_output_area.objects = objects

    def get_panel(self, as_dict=False, **kwargs):
        """Return tabs with widgets to interact with the flow."""

        d = {}

        #row = pn.Row(bkw.PreText(text=self.ddb.to_string(verbose=self.verbose), sizing_mode="scale_both"))
        d["Status"] = self.get_status_view()
        #d["Task"] = pn.Column(self.task_btn, self.on_task_btn)
        d["Task"] = self.get_task_view()
        d["History"] = self.get_history_view()
        d["Events"] = self.get_events_view()
        d["Corrections"] = pn.Row(self.corrections_btn, self.on_corrections_btn)
        d["Handlers"] = pn.Row(self.handlers_btn, self.on_handlers_btn)
        ##d["Structures"] = pn.Row(pn.Column(self.structures_io_checkbox, self.structures_btn), self.on_structures_btn)
        ###ws = pn.Column(self.ebands_plotter_mode, self.ebands_ksamp_checkbox, self.ebands_df_checkbox, self.ebands_plotter_btn)
        ###d["Ebands"] = pn.Row(ws, self.on_ebands_btn)
        d["Abivars"] = pn.Row(pn.Column(self.vars_text, self.vars_btn), self.on_vars_btn)
        ###d["Dims"] = pn.Row(pn.Column(self.dims_btn), self.on_dims_btn)
        d["Debug"] = self.get_debug_view()
        d["Browse"] = self.get_workdir_view()
        d["Graphviz"] = pn.Row(pn.Column(self.engine, self.dirtree, self.graphviz_btn), self.on_graphviz_btn)

        if as_dict: return d

        return self.get_template_from_tabs(d, template=kwargs.get("template", None), closable=False)
