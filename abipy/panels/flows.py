""""Panels for AbiPy flows."""
import param
import panel as pn
import bokeh.models.widgets as bw

from io import StringIO
from abipy.dynamics.hist import HistFile
from abipy.panels.core import AbipyParameterized


#def _mp(fig):
#    return pn.pane.Matplotlib(fig)


class FlowPanel(AbipyParameterized):
    """

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: FlowPanel
    """
    verbose = pn.widgets.IntSlider(start=0, end=10, step=1, value=0)

    engine = pn.widgets.Select(value="dot",
        options=['dot', 'neato', 'twopi', 'circo', 'fdp', 'sfdp', 'patchwork', 'osage'])
    dirtree = pn.widgets.Checkbox(name='Dirtree', value=False)
    graphviz_btn = pn.widgets.Button(name="Show graph", button_type='primary')

    status_btn = pn.widgets.Button(name="Show status", button_type='primary')
    history_btn = pn.widgets.Button(name="Show history", button_type='primary')
    debug_btn = pn.widgets.Button(name="Debug", button_type='primary')
    events_btn = pn.widgets.Button(name="Events", button_type='primary')
    corrections_btn = pn.widgets.Button(name="Corrections", button_type='primary')
    handlers_btn = pn.widgets.Button(name="Handlers", button_type='primary')

    #what_list = pn.widgets.CheckBoxGroup(name='Select', value=_what_list, options=_what_list, inline=False)

    #TODO: Implement widget for selected_nids(flow, options),
    #radio_group = pn.widgets.RadioButtonGroup(
    #   name='Radio Button Group', options=['Biology', 'Chemistry', 'Physics'], button_type='success')

    def __init__(self, flow, **params):
        super().__init__(**params)
        self.flow = flow

    @param.depends('status_btn.clicks')
    def on_status_btn(self):
        stream = StringIO()
        self.flow.show_status(stream=stream, verbose=self.verbose.value)
        return pn.Row(bw.PreText(text=stream.getvalue()))

    @param.depends('history_btn.clicks')
    def on_history_btn(self):
        stream = StringIO()
        #flow.show_history(status=options.task_status, nids=selected_nids(flow, options),
        #                  full_history=options.full_history, metadata=options.metadata)
        self.flow.show_history(stream=stream)
        return pn.Row(bw.PreText(text=stream.getvalue()))

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
        #TODO https://github.com/ralphbean/ansi2html ?
        stream = StringIO()
        #flow.debug(status=options.task_status, nids=selected_nids(flow, options))
        self.flow.debug(stream=stream)
        return pn.Row(bw.PreText(text=stream.getvalue()))

    @param.depends('events_btn.clicks')
    def on_events_btn(self):
        stream = StringIO()
        self.flow.show_events(stream=stream)
        #flow.show_events(status=options.task_status, nids=selected_nids(flow, options))
        return pn.Row(bw.PreText(text=stream.getvalue()))

    @param.depends('corrections_btn.clicks')
    def on_corrections_btn(self):
        stream = StringIO()
        self.flow.show_corrections(stream=stream)
        #flow.show_corrections(status=options.task_status, nids=selected_nids(flow, options))
        return pn.Row(bw.PreText(text=stream.getvalue()))

    @param.depends('handlers_btn.clicks')
    def on_handlers_btn(self):
        stream = StringIO()
        #if options.doc:
        #    flowtk.autodoc_event_handlers()
        #else:
        #show_events(self, status=None, nids=None, stream=sys.stdout):
        self.flow.show_event_handlers(verbose=self.verbose.value, stream=stream)
        return pn.Row(bw.PreText(text=stream.getvalue()))

    def get_panel(self):
        """Return tabs with widgets to interact with the flow."""
        tabs = pn.Tabs()
        #row = pn.Row(bw.PreText(text=self.ddb.to_string(verbose=self.verbose.value), sizing_mode="scale_both"))
        tabs.append(("Status", pn.Row(self.status_btn, self.on_status_btn)))
        tabs.append(("History", pn.Row(self.history_btn, self.on_history_btn)))
        tabs.append(("Events", pn.Row(self.events_btn, self.on_events_btn)))
        tabs.append(("Corrections", pn.Row(self.corrections_btn, self.on_corrections_btn)))
        tabs.append(("Handlers", pn.Row(self.handlers_btn, self.on_handlers_btn)))
        tabs.append(("Debug", pn.Row(self.debug_btn, self.on_debug_btn)))
        tabs.append(("Graphviz", pn.Row(pn.Column(self.engine, self.dirtree, self.graphviz_btn),
                                        self.on_graphviz_btn)))
        return tabs
