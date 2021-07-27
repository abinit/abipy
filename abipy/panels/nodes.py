""""Panels to interact with AbiPy flows."""
import textwrap
import traceback
import param
import bokeh.models.widgets as bkw
import panel as pn
import panel.widgets as pnw

#from panel.viewable import Viewer
from abipy.panels.core import AbipyParameterized, Loading, ButtonContext, depends_on_btn_click, dfc #, #, mpl, ply,
#from abipy import flowtk


class FilePathSelect(pnw.Select):

    @classmethod
    def from_filepaths(cls, filepaths, filter_files=True, **kwargs):
        import os
        #items = [(os.path.basename(p), p) for p in filepaths ]
        #if filter_files:
        #    def filter_basename(name):
        #    items = [t for t in items if not filte
        #
        #base2path = dict(*items)

        base2path = {os.path.basename(p): p for p in filepaths}
        new = cls(options=list(base2path.keys()), **kwargs)
        new._base2path = base2path
        return new

    @property
    def filepath(self):
        return self._base2path[self.value]

    def __bool__(self):
        return bool(self._base2path)



class NodeParameterized(AbipyParameterized):
    """

    """

    def __init__(self, node, **params):
        super().__init__(**params)
        self.node = node

        if node.is_flow:
            self.flow = node
            self.nids = []
            for work in node:
                self.nids.append(work.node_id)
                self.nids.extend([task.node_id for task in work])

        elif node.is_work:
            self.flow = self.node.flow
            self.nids = [task.node_id for task in node]

        elif node.is_task:
            self.flow = self.node.flow
            self.nids = self.node.node_id

        else:
            raise ValueError(f"Don't know how to hande type: `{type(node)}`")

        self.engine = pnw.Select(value="fdp",
                        options=['dot', 'neato', 'twopi', 'circo', 'fdp', 'sfdp', 'patchwork', 'osage'])
        self.dirtree = pnw.Checkbox(name='Dirtree', value=False)
        self.graphviz_btn = pnw.Button(name="Show Graph", button_type='primary')

        self.status_btn = pnw.Button(name="Show Status", button_type='primary')

        self.history_btn = pnw.Button(name="Show history", button_type='primary')
        self.debug_btn = pnw.Button(name="Debug", button_type='primary')
        self.events_btn = pnw.Button(name="Show Events", button_type='primary')
        self.corrections_btn = pnw.Button(name="Show Corrections", button_type='primary')
        self.handlers_btn = pnw.Button(name="Show Handlers", button_type='primary')
        self.vars_text = pnw.TextInput(name='Abivars',
                placeholder='Enter list of variables separated by comma e.g. `ecut, natom`')
        self.vars_btn = pnw.Button(name="Show Variables", button_type='primary')
        #self.dims_btn = pnw.Button(name="Show Dimensions", button_type='primary')

        self.workdir_fileselector = pn.widgets.FileSelector(node.workdir, only_files=True)
        self.outdir_fileselector = pn.widgets.FileSelector(node.outdir.path)
        self.indir_fileselector = pn.widgets.FileSelector(node.indir.path)

        # Create select widgets with the files in indir/outdir/workdir
        # Use basenames as items but remember that we need to abspath when opening the file.
        from abipy.flowtk.utils import Directory
        self.filepath_select_dir = {}
        for where in ("indir", "outdir", "workdir"):
            directory = Directory(self.node.workdir) if where == "workdir" else getattr(self.node, where)
            filepaths = directory.list_filepaths()
            self.filepath_select_dir[where] = FilePathSelect.from_filepaths(filepaths) #, name=f"Files in {where}")

    def get_status_view(self):
        return pn.Column(
            f"## Show the status of: `{repr(self.node)}`",
            pn.Row(
                self.pws_col(["verbose", "status_btn"]),
                bkw.PreText(text=self.node.str_deps()),
            ),
            pn.layout.Divider(),
            self.on_status_btn,
            sizing_mode='stretch_width',
            )

    @depends_on_btn_click("status_btn")
    def on_status_btn(self):
        """
        Show the status of the node.
        """
        stream = pnw.Terminal(output="\n\n",
            height=1200, # Need this one else the terminal is not show properly
            sizing_mode='stretch_width',
        )

        self.flow.show_status(nids=self.nids, stream=stream, verbose=self.verbose)

        if self.node.is_flow:
            stream.write("\n")
            self.flow.show_summary(stream=stream)

        return stream

    def get_history_view(self):
        return pn.Column(
            f"## Show the history of: `{repr(self.node)}`",
            self.pws_col(["verbose", "history_btn"]),
            self.on_history_btn,
            pn.layout.Divider(),
            sizing_mode='stretch_width',
            )

    @depends_on_btn_click("history_btn")
    def on_history_btn(self):
        """
        Show the history of the node.
        """
        stream = pnw.Terminal(output="\n\n",
            height=1200, # Need this one else the terminal is not show properly
            sizing_mode='stretch_width',
        )

        self.flow.show_history(nids=self.nids,
                               stream=stream,
                               #status=options.task_status,
                               #full_history=options.full_history, metadata=options.metadata
                               )
        #self.flow.show_history(stream=stream)
        return stream

    @depends_on_btn_click("graphviz_btn")
    def on_graphviz_btn(self):
        """
        Visualize the node dependencies with graphviz.
        """
        if self.dirtree.value:
            graph = self.node.get_graphviz_dirtree(engine=self.engine.value)
        else:
            graph = self.node.get_graphviz(engine=self.engine.value)

        #self.flow.plot_networkx(mode="network", with_edge_labels=False, ax=None, arrows=False,
        #                        node_size="num_cores", node_label="name_class", layout_type="spring", **kwargs):

        return pn.Column(graph)

    def get_debug_view(self):
        return pn.Column(
            f"## Debug node:`{repr(self.node)}`",
            self.pws_col(["verbose", "debug_btn"]),
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
        self.flow.debug(stream=stream, nids=self.nids) # status=options.task_status,
        return stream

    def get_events_view(self):
        return pn.Column(
            f"## Show the events of node: `{repr(self.node)}`",
            self.pws_col(["verbose", "events_btn"]),
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
        self.flow.show_events(stream=stream, nids=self.nids) # status=options.task_status,
        return stream

    @depends_on_btn_click("corrections_btn")
    def on_corrections_btn(self):
        stream = pnw.Terminal(output="\n\n",
            height=1200, # Need this one else the terminal is not show properly
            sizing_mode='stretch_width',
        )
        self.flow.show_corrections(stream=stream, nids=self.nids)
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
        #show_events(self, stream=sys.stdout):
        self.flow.show_event_handlers(stream=stream, verbose=self.verbose) #, nids=self.nids,  status=None,
        return stream

    @depends_on_btn_click("vars_btn")
    def on_vars_btn(self):
        if not self.vars_text.value: return
        varnames = [s.strip() for s in self.vars_text.value.split(",")]
        df = self.flow.compare_abivars(varnames=varnames, nids=self.nids,
                                       printout=False, with_colors=False)
        return pn.Row(dfc(df))

    @depends_on_btn_click("dims_btn")
    def on_dims_btn(self):
        df = self.flow.get_dims_dataframe(nids=self.nids, printout=False, with_colors=False)
        return pn.Row(dfc(df), sizing_mode="scale_width")

    #def get_workdir_view(self):

    #    col = pn.Column(sizing_mode="stretch_width"); ca = col.append
    #    ca(self.workdir_fileselector)
    #    self.workdir_selector_btn = pnw.Button(name="Show file", button_type='primary')
    #    self.workdir_selector_btn.on_click(self.on_workdir_selector_btn)
    #    ca(self.workdir_selector_btn)
    #    self.workdir_selector_output_area = pn.Column(sizing_mode="stretch_width")
    #    ca(self.workdir_selector_output_area)

    #    return pn.Column(col, self.workdir_selector_output_area, sizing_mode="stretch_width")

    def get_files_in_dir_view(self, where):
        select = self.filepath_select_dir[where]
        btn = pnw.Button(name="Analyze", button_type='primary')
        output_area = pn.Column(sizing_mode="stretch_width")

        from abipy.abilab import abiopen
        def update_output_area(event):
            with ButtonContext(btn), Loading(output_area):
                try:
                    # Cannot close the file at this level because it may be needed by the new app.
                    abifile = abiopen(select.filepath)
                    output_area.objects = [abifile.get_panel()]
                except Exception as exc:
                    #print(exc)
                    obj = pn.pane.Markdown("```shell\n%s\n```" % traceback.format_exc())
                    output_area.objects = [obj]

        btn.on_click(update_output_area)

        return pn.Column(
                "## Select the file and click the button to analyze the data",
                pn.Row(select, btn),
                pn.layout.Divider(),
                output_area,
                sizing_mode="stretch_width"
        )

    def on_workdir_selector_btn(self, event):
        """hello word"""
        filepaths = self.workdir_fileselector.value
        if not filepaths:
            objects = [pn.pane.Alert("## No file selected", altert_type="warning")]

        else:
            from abipy.abilab import abiopen
            objects = []
            for path in filepaths:
                try:
                    abifile = abiopen(path)
                    pn_obj = abifile.get_panel()
                except Exception:
                    pn_obj = pn.pane.Alert(str(exc), altert_type="warning")

                objects.append(pn_obj)

        self.workdir_selector_output_area.objects = objects

    def get_panel(self, as_dict=False, **kwargs):
        """Return tabs with widgets to interact with the flow."""

        d = {}

        d["Status"] = self.get_status_view()
        d["History"] = self.get_history_view()
        d["Events"] = self.get_events_view()
        d["Corrections"] = pn.Row(self.corrections_btn, self.on_corrections_btn)
        d["Handlers"] = pn.Row(self.handlers_btn, self.on_handlers_btn)
        d["Abivars"] = pn.Row(pn.Column(self.vars_text, self.vars_btn), self.on_vars_btn)
        #d["Dims"] = pn.Row(pn.Column(self.dims_btn), self.on_dims_btn)
        d["Debug"] = self.get_debug_view()
        #d["Browse"] = self.get_workdir_view()
        d["Graphviz"] = pn.Row(pn.Column(self.engine, self.dirtree, self.graphviz_btn), self.on_graphviz_btn)

        for where in ("workdir", "outdir", "indir"):
            if self.filepath_select_dir[where]:
                d[where.capitalize()] = self.get_files_in_dir_view(where)

        if as_dict: return d

        return self.get_template_from_tabs(d, template=kwargs.get("template", None), closable=False)
