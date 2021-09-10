""""Panels to interact with AbiPy flows."""
from __future__ import annotations

import textwrap
import traceback

import pandas as pd
import param
import bokeh.models.widgets as bkw
import panel as pn
import panel.widgets as pnw

from abipy.panels.core import AbipyParameterized, Loading, ButtonContext, depends_on_btn_click, dfc, ply
from abipy.flowtk.nodes import Node
#from abipy import flowtk


class FilePathSelect(pnw.Select):

    @classmethod
    def from_filepaths(cls, filepaths, filter_files=True, **kwargs):
        import os
        items = [(os.path.basename(p), p) for p in filepaths ]

        if filter_files:
            def filter_basename(name):
                if name.startswith(".") or name.endswith(".pickle"):
                    return False
                return True

            items = [t for t in items if filter_basename(t[0])]

        base2path = dict(items)
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

    def __init__(self, node: Node, **params):
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
            raise ValueError(f"Don't know how to handle type: `{type(node)}`")

        self.engine = pnw.Select(value="fdp", name="engine",
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

        self.workdir_fileselector = pnw.FileSelector(node.workdir, only_files=True)
        self.outdir_fileselector = pnw.FileSelector(node.outdir.path)
        self.indir_fileselector = pnw.FileSelector(node.indir.path)

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
                self.wdg_box(["verbose", "status_btn"]),
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
        term = pnw.Terminal(
            output="\n\n",
            height=1200, # Need this one else the terminal is not shown properly
            sizing_mode='stretch_width',
        )
        term.write("\n")

        if self.node.is_task:
            return term

        # Here it's important to enforce verbose 1 else show_status
        # does not analyze the tasks that are completed.
        df = self.flow.show_status(nids=self.nids, stream=term, verbose=1, return_df=True) #self.verbose)

        return pn.Column(
                StatusCards(df),
                sizing_mode="stretch_width",
        )

        # TODO: Finalize the implementation.
        # Generate heatmap with plotly
        max_num_tasks = max(len(work) for work in self.flow)
        y = [f"w{i}" for i in range(len(self.flow))]
        x = [f"t{i}" for i in range(max_num_tasks)]
        z = []
        for work in self.flow:
            row = [None for i in range(max_num_tasks)]
            for i, task in enumerate(work):
                #row[i] = task.status
                row[i] = task.mpi_procs
            z.append(row)

        import plotly.graph_objects as go
        fig = go.Figure(data=go.Heatmap(
                        x=x, y=y, z=z,
                        hoverongaps=False, transpose=False, colorscale="Viridis",
                        ))
        fig.update_layout(title_text="Number of MPI procs", title_x=0.5)

        from abipy.tools.plotting import add_colorscale_dropwdowns
        add_colorscale_dropwdowns(fig)

        return ply(fig)

    def get_history_view(self):
        return pn.Column(
            f"## Show the history of: `{repr(self.node)}`",
            self.wdg_box(["verbose", "history_btn"]),
            pn.layout.Divider(),
            self.on_history_btn,
            sizing_mode='stretch_width',
       )

    @depends_on_btn_click("history_btn")
    def on_history_btn(self):
        """
        Show the history of the node.
        """
        term = pnw.Terminal(
            output="\n\n",
            height=1200, # Need this one else the terminal is not show properly
            sizing_mode='stretch_width',
        )

        self.flow.show_history(nids=self.nids,
                               stream=term,
                               #status=options.task_status,
                               #full_history=options.full_history,
                               #metadata=options.metadata
                               )
        return term

    def get_graphviz_view(self):
        return pn.Column(
                f"## Graphviz options for node: `{repr(self.node)}`",
                pn.WidgetBox(self.engine, self.dirtree, self.graphviz_btn),
                pn.layout.Divider(),
                self.on_graphviz_btn,
                sizing_mode="stretch_width"
        )

    @depends_on_btn_click("graphviz_btn")
    def on_graphviz_btn(self):
        """
        Visualize node dependencies with [graphviz package](https://graphviz.readthedocs.io/en/stable/index.html)
        """
        if self.dirtree.value:
            graph = self.node.get_graphviz_dirtree(engine=self.engine.value)
        else:
            graph = self.node.get_graphviz(engine=self.engine.value)

        #self.flow.plot_networkx(mode="network", with_edge_labels=False, ax=None, arrows=False,
        #                        node_size="num_cores", node_label="name_class", layout_type="spring", **kwargs):

        return pn.Column(
            "## Dependency Graph:",
            pn.pane.SVG(graph),
            sizing_mode="stretch_width"
        )

    def get_debug_view(self):
        return pn.Column(
            f"## Debug node:`{repr(self.node)}`",
            self.pws_col(["verbose", "debug_btn"]),
            self.on_debug_btn,
            pn.layout.Divider(),
            sizing_mode='stretch_width',
        )

        #d["Corrections"] = pn.Row(self.corrections_btn, self.on_corrections_btn)
        #d["Handlers"] = pn.Row(self.handlers_btn, self.on_handlers_btn)

    @depends_on_btn_click("debug_btn")
    def on_debug_btn(self):
        term = pnw.Terminal(output="\n\n",
            height=1200, # Need this one else the terminal is not show properly
            sizing_mode='stretch_width',
        )
        self.flow.debug(stream=term, nids=self.nids) # status=options.task_status,
        return term

    def get_events_view(self):
        return pn.Column(
            f"## Show the events of: `{repr(self.node)}`",
            self.pws_col(["verbose", "events_btn"]),
            self.on_events_btn,
            pn.layout.Divider(),
            sizing_mode='stretch_width',
        )

    @depends_on_btn_click("events_btn")
    def on_events_btn(self):
        term = pnw.Terminal(
            output="\n\n",
            height=1200, # Need this one else the terminal is not show properly
            sizing_mode='stretch_width',
        )
        self.flow.show_events(stream=term, nids=self.nids) # status=options.task_status,
        return term

    @depends_on_btn_click("corrections_btn")
    def on_corrections_btn(self):
        term = pnw.Terminal(
            output="\n\n",
            height=1200, # Need this one else the terminal is not show properly
            sizing_mode='stretch_width',
        )
        self.flow.show_corrections(stream=term, nids=self.nids)
        return term

    @depends_on_btn_click("handlers_btn")
    def on_handlers_btn(self):
        term = pnw.Terminal(
            output="\n\n",
            height=1200, # Need this one else the terminal is not show properly
            sizing_mode='stretch_width',
        )
        self.flow.show_event_handlers(stream=term, verbose=self.verbose) #, nids=self.nids,  status=None,
        return term

    def get_dims_and_vars_view(self):
        row = pn.Row(pn.Column(self.vars_text, self.vars_btn), self.on_vars_btn)
        return row
        #d["Dims"] = pn.Row(pn.Column(self.dims_btn), self.on_dims_btn)

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

    def get_files_in_dir_view(self, where):
        """
        Return None if no file is found
        """
        select = self.filepath_select_dir[where]
        if not select: return None

        btn = pnw.Button(name="Analyze", button_type='primary')
        output_area = pn.Column(sizing_mode="stretch_width")

        from abipy.abilab import abiopen
        #from .core import NcFileViewer

        def update_output_area(event):
            with ButtonContext(btn), Loading(output_area):
                try:
                    # Cannot close the file at this level because it may be needed by the new app.
                    abifile = abiopen(select.filepath)
                    output_area.objects = [abifile.get_panel()]
                except Exception as exc:
                    #print(exc)
                    #if select.filepath.endswith(".nc"):
                    #    # We have a nc file but it's not supported by abiopen.
                    #    # Let's create a minimalistic view of the netcdf dims/vars
                    #    #abifile = AbinitNcFile(select.filepath)
                    #    NcFileViewer(self).get_ncfile_view(**kwargs)
                    #    output_area.objects = [abifile.get_ncfile_view()]
                    #else:
                    obj = pn.pane.Markdown("```shell\n%s\n```" % traceback.format_exc())
                    output_area.objects = [obj]

        btn.on_click(update_output_area)

        return pn.Column(
            "## Select a file and click the button to analyze the data",
            pn.WidgetBox(select, btn),
            pn.layout.Divider(),
            output_area,
            sizing_mode="stretch_width"
        )

    #def on_workdir_selector_btn(self, event):
    #    """hello word"""
    #    filepaths = self.workdir_fileselector.value
    #    if not filepaths:
    #        objects = [pn.pane.Alert("## No file selected", altert_type="warning")]

    #    else:
    #        from abipy.abilab import abiopen
    #        objects = []
    #        for path in filepaths:
    #            try:
    #                abifile = abiopen(path)
    #                pn_obj = abifile.get_panel()
    #            except Exception:
    #                pn_obj = pn.pane.Alert(str(exc), altert_type="warning")

    #            objects.append(pn_obj)

    #    self.workdir_selector_output_area.objects = objects

    def get_panel(self, as_dict=False, **kwargs):
        """
        Return tabs with widgets to interact with the flow.
        """
        d = {}

        d["Status"] = self.get_status_view()
        d["History"] = self.get_history_view()
        d["Events"] = self.get_events_view()
        for where in ("workdir", "outdir", "indir"):
            if self.filepath_select_dir[where]:
                view = self.get_files_in_dir_view(where)
                if view is not None: d[where.capitalize()] = view

        d["Debug"] = self.get_debug_view()

        if not self.node.is_task:
            d["Dims & Vars"] = self.get_dims_and_vars_view()
        d["Graphviz"] = self.get_graphviz_view()

        if as_dict: return d

        return self.get_template_from_tabs(d, template=kwargs.get("template", None), closable=False)


class StatusCards(param.Parameterized):

    def __init__(self, df: pd.DataFrame, **params):
        self.df = df
        super().__init__(**params)

        # For each work_idx, compute the min/max index
        # we need this dict in add_vrect to shadow the region associated to the work.
        self.w_start_stop = {}
        for w_idx in df["work_idx"].unique():
            lst = sorted(df.index[(df['work_idx'] == w_idx)].tolist())
            self.w_start_stop[w_idx] = (lst[0], lst[-1])

        header_list = [
            "## DataFrame",
            "## Task status histogram",
            "## Task Timeline",
            "## Task class histogram",
            "## Runtime in seconds for each task in the flow (-1 if task is not running)",
            "## Number of WARNINGs found in the log file of the Task",
            "## Barplot with number of MPI procs",
        ]

        #if len(df['status'].unique()) != 1:
        # Show histogram with task status only if we have different status values.

        self.cards, self.done = {}, {}
        for header in header_list:
            card = pn.layout.Card(None, header=header, collapsed=True, sizing_mode="stretch_width")
            # Compute stuff only when user opens the card.
            card.param.watch(self.update_card, ['collapsed'])
            self.cards[header] = card
            self.done[header] = False

        # Compute this card
        self.cards["## Task status histogram"].collapsed = False

        open_btn = pnw.Button(name="Open all cards", button_type='primary')
        open_btn.on_click(self.open_all_cards)
        close_btn = pnw.Button(name="Close all cards", button_type='primary')
        close_btn.on_click(self.close_all_cards)

        self.layout = pn.Column(pn.Row(open_btn, close_btn),
                                *list(self.cards.values()),
                                sizing_mode="stretch_width"
                                )

    def __panel__(self):
        return self.layout

    def open_all_cards(self, event):
        for card in self.cards.values():
            card.collapsed = False

    def close_all_cards(self, event):
        for card in self.cards.values():
            card.collapsed = True

    def add_vrect_to_fig(self, fig):
        """
        Add vertical rectangles to the plotly fig in order to group tasks belonging to the same Work.
        """
        for w_idx, (x0, x1) in self.w_start_stop.items():
            fig.add_vrect(x0=x0, x1=x1,
                          annotation_text=f"w{w_idx}", annotation_position="top left",
                          fillcolor="grey", opacity=0.1, line_width=0)
        return fig

    def update_card(self, event):
        """
        This callback is triggered when the user opens/closes the card.
        Here we compute and display the plot/table the first time the card is opened.
        """
        header = event.obj.header
        if self.done[header]: return

        card = self.cards[header]
        has_pane = False
        import plotly.express as px
        df = self.df

        with Loading(card):
            if header == "## Task status histogram":
                fig = px.histogram(df, x="status")

            elif header == "## Task class histogram":
                fig = px.histogram(df, x="task_class", color="status")

            elif header == "## Runtime in seconds for each task in the flow (-1 if task is not running)":

                fig = px.scatter(df, x=df.index, y="task_runtime_s", color="status", #, symbol="work_idx", #size=
                                 hover_data =["num_warnings", "num_comments", "work_idx", "task_class"],
                                 hover_name="name")
                self.add_vrect_to_fig(fig)

            elif header == "## Number of WARNINGs found in the log file of the Task":

                fig = px.scatter(df, x=df.index, y="num_warnings", color="status", #, symbol="work_idx", #size=
                                 hover_data=["task_runtime_s", "num_comments", "work_idx", "task_class"],
                                 hover_name="name")
                self.add_vrect_to_fig(fig)

            elif header == "## Barplot with number of MPI procs":

                fig = px.bar(df, x="work_idx", y="mpi_procs", color="task_widx",
                             #barmode="group",
                             #color="status",
                             #pattern_shape="task_class",
                             #pattern_shape="status",
                             #pattern_shape_sequence=[".", "x", "+"],
                             hover_data =["num_warnings", "num_comments", "task_class"],
                             hover_name="name")

            elif header == "## Task Timeline":
                fig = px.timeline(df, x_start="start_datetime", x_end="end_datetime", y="name",
                                  color="task_class", hover_name="name")

            elif header == "## DataFrame":
                # Remove some columns as well as the index.
                simple_df = df.drop(columns=["node_id", "queue_id", "qname",
                                             "task_queue_time_s", "submission_datetime",
                                             "start_datetime", "end_datetime", "task_widx"]) # "work_idx",
                simple_df.set_index('name', inplace=True)

                pane = pnw.Tabulator(simple_df, groupby=['work_idx']) #, height=240)
                has_pane = True

            else:
                raise ValueError(f"No action registered for: `{header}`")

        if not has_pane:
            pane = ply(fig, with_divider=False)

        card.objects = [pane]
        self.done[header] = True
