""""Panels to interact with the AbiPy tasks."""
import param
import panel as pn
import panel.widgets as pnw
import bokeh.models.widgets as bkw

from io import StringIO
from abipy.panels.core import AbipyParameterized, mpl, ply, dfc, depends_on_btn_click


class AbinitTaskPanel(AbipyParameterized):
    """
    Panel to interact with an AbiPy Task
    """

    def __init__(self, task, **params):
        self.task = task
        self.flow = task.flow
        self.nids = self.task.node_id

        self.engine = pnw.Select(value="fdp",
                      options=['dot', 'neato', 'twopi', 'circo', 'fdp', 'sfdp', 'patchwork', 'osage'])
        self.dirtree = pnw.Checkbox(name='Dirtree', value=False)
        self.graphviz_btn = pnw.Button(name="Show graph", button_type='primary')

        self.status_btn = pnw.Button(name="Show status", button_type='primary')
        #self.task_btn = pnw.Button(name="Show Task", button_type='primary')
        ##self.work_btn = pnw.Button(name="Show Task", button_type='primary')
        self.history_btn = pnw.Button(name="Show history", button_type='primary')
        self.debug_btn = pnw.Button(name="Debug", button_type='primary')
        self.events_btn = pnw.Button(name="Events", button_type='primary')
        #self.corrections_btn = pnw.Button(name="Corrections", button_type='primary')
        self.handlers_btn = pnw.Button(name="Handlers", button_type='primary')

        self.files_btn = pnw.Button(name="Show Input", button_type='primary')

        #self.dims_btn = pnw.Button(name="Show Dimensions", button_type='primary')

        #self.structures_btn = pnw.Button(name="Show Structures", button_type='primary')
        #self.structures_io_checkbox = pnw.CheckBoxGroup(
        #    name='Input/Output Structure', value=['output'], options=['input', 'output'], inline=Tru

        self.workdir_selector = pn.widgets.FileSelector(flow.workdir)

        super().__init__(**params)

    @depends_on_btn_click("status_btn")
    def on_status_btn(self):
        stream = StringIO()
        self.flow.show_status(stream=stream, nids=self.nids, verbose=self.verbose)
        return pn.Row(bkw.PreText(text=stream.getvalue()))

    @depends_on_btn_click("history_btn")
    def on_history_btn(self):
        stream = StringIO()
        self.flow.show_history(nids=self.nids, stream=stream)
        return pn.Row(bkw.PreText(text=stream.getvalue()))

    @depends_on_btn_click("graphviz_btn")
    def on_graphviz_btn(self):
        """
        Visualize the flow with graphviz.
        """
        node = self.task
        if self.dirtree.value:
            graph = node.get_graphviz_dirtree(engine=self.engine.value)
        else:
            graph = node.get_graphviz(engine=self.engine.value)

        return pn.Column(graph)

    @depends_on_btn_click("debug_btn")
    def on_debug_btn(self):
        #TODO https://github.com/ralphbean/ansi2html ?
        stream = StringIO()
        #flow.debug(status=options.task_status, nids=selected_nids(flow, options))
        self.flow.debug(stream=stream, nids=self.nids)
        return pn.Row(bkw.PreText(text=stream.getvalue()))

    @depends_on_btn_click("events_btn")
    def on_events_btn(self):
        stream = StringIO()
        self.flow.show_events(nids=self.nids, stream=stream)
        return pn.Row(bkw.PreText(text=stream.getvalue()))

    @depends_on_btn_click("corrections_btn")
    def on_corrections_btn(self):
        stream = StringIO()
        self.flow.show_corrections(stream=stream, nids=self.nids)
        #flow.show_corrections(status=options.task_status, nids=selected_nids(flow, options))
        return pn.Row(bkw.PreText(text=stream.getvalue()))

    @depends_on_btn_click("handlers_btn")
    def on_handlers_btn(self):
        stream = StringIO()
        #if options.doc:
        #    flowtk.autodoc_event_handlers()
        #else:
        #show_events(self, status=None, nids=None, stream=sys.stdout):
        self.flow.show_event_handlers(verbose=self.verbose, nids=self.nids, stream=stream)
        return pn.Row(bkw.PreText(text=stream.getvalue()))

    @depends_on_btn_click("files_btn")
    def on_files_btn(self):
        """Show the input files of the task."""

        col = pn.Column(sizing_mode="stretch_width")
        ca = col.append; cext = col.extend

        #input_html = self.html_with_clipboard_btn(self.task.input._repr_html_())
        input_html = self.html_with_clipboard_btn(self.task.input)
        cext(["## Input", input_html, pn.layout.Divider()])

        #for fname in ["output_file", "job_file", "stderr_file", "mpiabort_file"]: # "log_file",
        #    file = getattr(self.task, fname)
        #    text = file.read() if file.exists else f"{fname} does not exist"
        #    #ace = pnw.Ace(value=text, language='text', readonly=True,
        #    #              sizing_mode='stretch_width', width=900)
        #    #              #sizing_mode='stretch_width', width=700, height=600)
        #    #ca(pn.Accordion((f"{fname}", ace), sizing_mode='stretch_both'))
        #    cext([f"## {fname}", entry, pn.layout.Divider()])

        return col

    @depends_on_btn_click("dims_btn")
    def on_dims_btn(self):
        df = self.flow.get_dims_dataframe(nids=self.nids,
                                          printout=False, with_colors=False)
        return pn.Row(dfc(df), sizing_mode="scale_width")

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

        d = {}

        #row = pn.Row(bkw.PreText(text=self.ddb.to_string(verbose=self.verbose), sizing_mode="scale_
        d["Status"] = pn.Row(self.status_btn, self.on_status_btn)
        d["History"] = pn.Row(self.history_btn, self.on_history_btn)
        d["Events"] = pn.Row(self.events_btn, self.on_events_btn)
        ##d["Corrections"] = pn.Row(self.corrections_btn, self.on_corrections_btn)
        ##d["Handlers"] = pn.Row(self.handlers_btn, self.on_handlers_btn)
        ##d["Structures"] = pn.Row(pn.Column(self.structures_io_checkbox, self.structures_btn), self
        ###ws = pn.Column(self.ebands_plotter_mode, self.ebands_ksamp_checkbox, self.ebands_df_check
        ###d["Ebands"] = pn.Row(ws, self.on_ebands_btn)
        d["Input"] = pn.Row(pn.Column("## Input button", self.files_btn), self.on_files_btn)
        ###d["Dims"] = pn.Row(pn.Column(self.dims_btn), self.on_dims_btn)
        ##d["Debug"] = pn.Row(self.debug_btn, self.on_debug_btn)
        d["Graphviz"] = pn.Row(pn.Column(self.engine, self.dirtree, self.graphviz_btn), self.on_graphviz_btn)

        if as_dict: return d

        return self.get_template_from_tabs(d, template=kwargs.get("template", None), closable=False)
