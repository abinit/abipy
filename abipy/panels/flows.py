""""Panels for AbiPy flows."""
import param
import panel as pn
import panel.widgets as pnw
import bokeh.models.widgets as bkw

from io import StringIO
from abipy.panels.core import AbipyParameterized


class FlowPanel(AbipyParameterized):
    """
    """
    verbose = pn.widgets.IntSlider(start=0, end=10, step=1, value=0)

    engine = pn.widgets.Select(value="fdp",
        options=['dot', 'neato', 'twopi', 'circo', 'fdp', 'sfdp', 'patchwork', 'osage'])
    dirtree = pn.widgets.Checkbox(name='Dirtree', value=False)
    graphviz_btn = pn.widgets.Button(name="Show graph", button_type='primary')

    status_btn = pn.widgets.Button(name="Show status", button_type='primary')
    history_btn = pn.widgets.Button(name="Show history", button_type='primary')
    debug_btn = pn.widgets.Button(name="Debug", button_type='primary')
    events_btn = pn.widgets.Button(name="Events", button_type='primary')
    corrections_btn = pn.widgets.Button(name="Corrections", button_type='primary')
    handlers_btn = pn.widgets.Button(name="Handlers", button_type='primary')

    vars_text = pn.widgets.TextInput(name='Abivars', placeholder='Enter list of variables separated by comma')
    vars_btn = pn.widgets.Button(name="Show Variables", button_type='primary')

    dims_btn = pn.widgets.Button(name="Show Dimensions", button_type='primary')

    structures_btn = pn.widgets.Button(name="Show Structures", button_type='primary')
    structures_io_checkbox = pn.widgets.CheckBoxGroup(
        name='Input/Output Structure', value=['output'], options=['input', 'output'], inline=True)

    # Widgets to plot ebands.
    ebands_btn = pn.widgets.Button(name="Show Ebands", button_type='primary')
    ebands_plotter_mode = pnw.Select(name="Plot Mode", value="gridplot",
        options=["gridplot", "combiplot", "boxplot", "combiboxplot"]) # "animate",
    ebands_plotter_btn = pnw.Button(name="Plot", button_type='primary')
    ebands_df_checkbox = pnw.Checkbox(name='With Ebands DataFrame', value=False)
    ebands_ksamp_checkbox = pn.widgets.CheckBoxGroup(
        name='Input/Output Structure', value=["with_path", "with_ibz"], options=['with_path', 'with_ibz'], inline=True)

    #TODO: Implement widget for selected_nids(flow, options),
    #radio_group = pn.widgets.RadioButtonGroup(
    #   name='Radio Button Group', options=['Biology', 'Chemistry', 'Physics'], button_type='success')

    def __init__(self, flow, **params):
        super().__init__(**params)
        self.flow = flow

    @param.depends('status_btn.clicks')
    def on_status_btn(self):
        if self.status_btn.clicks == 0: return
        stream = StringIO()
        self.flow.show_status(stream=stream, verbose=self.verbose.value)
        return pn.Row(bkw.PreText(text=stream.getvalue()))

    @param.depends('history_btn.clicks')
    def on_history_btn(self):
        if self.history_btn.clicks == 0: return
        stream = StringIO()
        #flow.show_history(status=options.task_status, nids=selected_nids(flow, options),
        #                  full_history=options.full_history, metadata=options.metadata)
        self.flow.show_history(stream=stream)
        return pn.Row(bkw.PreText(text=stream.getvalue()))

    @param.depends('graphviz_btn.clicks')
    def on_graphviz_btn(self):
        """
        """
        if self.graphviz_btn.clicks == 0: return
        node = self.flow
        if self.dirtree.value:
            graph = node.get_graphviz_dirtree(engine=self.engine.value)
        else:
            graph = node.get_graphviz(engine=self.engine.value)
        return pn.Column(graph)

    @param.depends('debug_btn.clicks')
    def on_debug_btn(self):
        if self.debug_btn.clicks == 0: return
        #TODO https://github.com/ralphbean/ansi2html ?
        stream = StringIO()
        #flow.debug(status=options.task_status, nids=selected_nids(flow, options))
        self.flow.debug(stream=stream)
        return pn.Row(bkw.PreText(text=stream.getvalue()))

    @param.depends('events_btn.clicks')
    def on_events_btn(self):
        if self.events_btn.clicks == 0: return
        stream = StringIO()
        self.flow.show_events(stream=stream)
        #flow.show_events(status=options.task_status, nids=selected_nids(flow, options))
        return pn.Row(bkw.PreText(text=stream.getvalue()))

    @param.depends('corrections_btn.clicks')
    def on_corrections_btn(self):
        if self.corrections_btn.clicks == 0: return
        stream = StringIO()
        self.flow.show_corrections(stream=stream)
        #flow.show_corrections(status=options.task_status, nids=selected_nids(flow, options))
        return pn.Row(bkw.PreText(text=stream.getvalue()))

    @param.depends('handlers_btn.clicks')
    def on_handlers_btn(self):
        #if self.handlers_btn.clicks == 0: return
        stream = StringIO()
        #if options.doc:
        #    flowtk.autodoc_event_handlers()
        #else:
        #show_events(self, status=None, nids=None, stream=sys.stdout):
        self.flow.show_event_handlers(verbose=self.verbose.value, stream=stream)
        return pn.Row(bkw.PreText(text=stream.getvalue()))

    @param.depends('vars_btn.clicks')
    def on_vars_btn(self):
        if self.vars_btn.clicks == 0: return
        if not self.vars_text.value: return
        varnames = [s.strip() for s in self.vars_text.value.split(",")]
        df = self.flow.compare_abivars(varnames=varnames, # nids=selected_nids(flow, options),
                                       printout=False, with_colors=False)
        return pn.Row(self._df(df))

    @param.depends('dims_btn.clicks')
    def on_dims_btn(self):
        if self.dims_btn.clicks == 0: return
        df = self.flow.get_dims_dataframe(# nids=selected_nids(flow, options),
                                          printout=False, with_colors=False)
        return pn.Row(self._df(df), sizing_mode="scale_width")

    @param.depends('structures_btn.clicks')
    def on_structures_btn(self):
        if self.structures_btn.clicks == 0: return
        what = ""
        if "input" in self.structures_io_checkbox.value: what += "i"
        if "output" in self.structures_io_checkbox.value: what += "o"
        dfs = self.flow.compare_structures(nids=None, # select_nids(flow, options),
                                           what=what,
                                           verbose=self.verbose.value, with_spglib=False, printout=False,
                                           with_colors=False)

        return pn.Row(self._df(dfs.lattice), sizing_mode="scale_width")

    @param.depends('ebands_plotter_btn.clicks')
    def on_ebands_btn(self):
        if self.ebands_plotter_btn.clicks == 0: return

        df, ebands_plotter = self.flow.compare_ebands(
                                nids=None, # select_nids(flow, options),
                                with_path="with_path" in self.ebands_ksamp_checkbox.value,
                                with_ibz="with_ibz" in self.ebands_ksamp_checkbox.value,
                                verbose=self.verbose.value,
                                with_spglib=False
                                )

        if ebands_plotter is None:
            return

        plot_mode = self.ebands_plotter_mode.value
        plotfunc = getattr(ebands_plotter, plot_mode, None)
        if plotfunc is None:
            raise ValueError("Don't know how to handle plot_mode: %s" % plot_mode)

        fig = plotfunc(**self.fig_kwargs)
        col = pn.Column(self._mp(fig))
        if self.ebands_df_checkbox.value:
            col.append(self._df(df))

        return pn.Row(col) #, sizing_mode='scale_width')

    def get_panel(self):
        """Return tabs with widgets to interact with the flow."""
        tabs = pn.Tabs(); app = tabs.append
        #row = pn.Row(bkw.PreText(text=self.ddb.to_string(verbose=self.verbose.value), sizing_mode="scale_both"))
        app(("Status", pn.Row(self.status_btn, self.on_status_btn)))
        app(("History", pn.Row(self.history_btn, self.on_history_btn)))
        app(("Events", pn.Row(self.events_btn, self.on_events_btn)))
        app(("Corrections", pn.Row(self.corrections_btn, self.on_corrections_btn)))
        app(("Handlers", pn.Row(self.handlers_btn, self.on_handlers_btn)))
        app(("Structures", pn.Row(pn.Column(self.structures_io_checkbox, self.structures_btn), self.on_structures_btn)))
        ws = pn.Column(self.ebands_plotter_mode, self.ebands_ksamp_checkbox, self.ebands_df_checkbox, self.ebands_plotter_btn)
        app(("Ebands", pn.Row(ws, self.on_ebands_btn)))
        app(("Abivars", pn.Row(pn.Column(self.vars_text, self.vars_btn), self.on_vars_btn)))
        app(("Dims", pn.Row(pn.Column(self.dims_btn), self.on_dims_btn)))
        app(("Debug", pn.Row(self.debug_btn, self.on_debug_btn)))
        app(("Graphviz", pn.Row(pn.Column(self.engine, self.dirtree, self.graphviz_btn),
                                self.on_graphviz_btn)))
        return tabs
