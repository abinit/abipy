""""MongoGUI."""
from __future__ import annotations

import json
#import param
import functools
import panel as pn
import panel.widgets as pnw

from typing import Type  #  List
from abipy.panels.core import AbipyParameterized, Loading  # depends_on_btn_click, ply,
from abipy.htc.base_models import QueryResults, MongoConnector
from abipy.htc.flow_models import FlowModel


def display_qr(method):
    @functools.wraps(method)
    def decorated(*args, **kwargs):
        self = args[0]
        with Loading(self.out_area):
            qr = method(*args, **kwargs)
            self._display_qr(qr)
    return decorated


class MongoGui(AbipyParameterized):
    """
    A MongoGUI is connected to MongoDB collection FlowModel documents (actually, a concrete
    subclass of FlowModel) via a MongoConnector.
    The GUI uses the abstract protocol of the FlowModel to provide widgets to filter
    analyze, aggregate and visualize the results.
    The abstract protocol should be already enough to cover the most important cases.
    In needed, one can define a subclass of MongoGUI providing extra tools to interact with that
    particular FlowModel.
    """

    def __init__(self, mng_connector: MongoConnector, flow_model_cls: Type[FlowModel]):
        super().__init__()
        self.mng_connector = mng_connector
        self.flow_model_cls = flow_model_cls

        self.build_ui()

    #def collection(self):
    #    return self.mng_connector.get_collection()

    def build_ui(self):
        # This is the part that should be abstracted out
        # so that one can implement customized subclasses
        # that rely on model.get_panel_view

        self.limit = 0

        self.filter_query = pnw.LiteralInput(name='MongoDB filter (python dict)',
                                             placeholder="e.g. {'key': value}", type=dict,
                                             value={'kppa': 300})

        self.filter_btn = pnw.Button(name="Run query", button_type='primary')
        self.filter_btn.on_click(self.on_filter_btn)

        collection = self.mng_connector.get_collection()

        # Get all spg_numbers in the collection.
        spg_numbers = self.flow_model_cls.mng_get_spg_numbers_incoll(collection)
        self.spg_number_wdg = pnw.Select(name='Filter by space group number', options=spg_numbers)
        #self.spg_number_wdg = pnw.IntInput(name='Filter by space group number', value=None, step=1, start=1, end=230)

        self.filter_by_spg_number_btn = pnw.Button(name="Run query", button_type='primary')
        self.filter_by_spg_number_btn.on_click(self.on_filter_by_spg_number_btn)

        self.formula_wdg = pnw.TextInput(name='Filter by reduced formula', placeholder='Enter e.g. Si or MgO...')
        self.filter_by_formula_btn = pnw.Button(name="Run query", button_type='primary')
        self.filter_by_formula_btn.on_click(self.on_filter_by_formula_btn)

        crystal_systems = self.flow_model_cls.mng_get_crystal_systems_incoll(collection)
        self.crystal_system_wdg = pnw.Select(name='Filter by crystal system', options=crystal_systems)
        self.filter_by_crystal_system_btn = pnw.Button(name="Run query", button_type='primary')
        self.filter_by_crystal_system_btn.on_click(self.on_filter_by_crystal_system_btn)

        self.preset_queries_btn = pnw.Button(name="Run query", button_type='primary')
        self.preset_queries_btn.on_click(self.on_preset_queries_btn)
        options = []
        # TODO: Fixme
        #if hasattr(self.flow_model_cls, "get_preset_queries"):
        #    options = [json.dumps(q) for q in self.flow_model_cls.get_preset_queries()]
        self.preset_queries_wdg = pnw.Select(name="Pommon MongoDB queries", options=options)

        md = pn.pane.Markdown("## Perform MongoDB queries using one of the widgets below:")

        widget_box = pn.WidgetBox(
            md,
            self.filter_query, self.filter_btn,
            pn.Row(self.spg_number_wdg, self.formula_wdg),
            pn.Row(self.filter_by_spg_number_btn, self.filter_by_formula_btn),
            self.crystal_system_wdg, self.filter_by_crystal_system_btn,
            self.preset_queries_wdg, self.preset_queries_btn,
            #sizing_mode="stretch_width",
        )

        self.out_area = pn.Column(sizing_mode="stretch_width")

        self.main = pn.Column(self.mng_connector,
                              widget_box,
                              #pn.layout.Divider(),
                              #"## Results:",
                              self.out_area,
                              sizing_mode="stretch_width")

    def get_app(self):
        """Build and return the application."""
        template = pn.template.FastListTemplate
        #template = pn.template.BootstrapTemplate
        #template = pn.template.ReactTemplate
        return template(main=self.main, title="AbiPyFlowModel MongoDB GUI")

    #@depends_on_btn_click('filter_btn', show_shared_wdg_warning=False)
    def on_filter_btn(self, _):
        with Loading(self.out_area):
            collection = self.mng_connector.get_collection()
            query = self.filter_query.value
            qr = self.flow_model_cls.mng_find(query, collection, limit=self.limit)
            self._display_qr(qr)

    @display_qr
    def on_filter_by_spg_number_btn(self, _):
        """Filter documents by space group number."""
        spg_number = self.spg_number_wdg.value
        collection = self.mng_connector.get_collection()
        qr = self.flow_model_cls.mng_find_by_spg_number(spg_number, collection, limit=self.limit)
        return qr

    @display_qr
    def on_filter_by_formula_btn(self, _):
        """Filter documents by formula."""
        reduced_formula = self.formula_wdg.value
        collection = self.mng_connector.get_collection()
        qr = self.flow_model_cls.mng_find_by_formula(reduced_formula, collection, limit=self.limit)
        return qr

    @display_qr
    def on_filter_by_crystal_system_btn(self, _):
        """Filter documents by crystal system."""
        crystal_system = self.crystal_system_wdg.value
        collection = self.mng_connector.get_collection()
        qr = self.flow_model_cls.mng_find_by_crystal_system(crystal_system, collection, limit=self.limit)
        return qr

    @display_qr
    def on_preset_queries_btn(self, _):
        """
        Show the list of preset queries if not index is provided by the user
        else perform the query and return the results.
        """
        query = self.preset_queries_wdg.value
        collection = self.mng_connector.get_collection()
        if query is None:
            return QueryResults.empty_from_query({})

        query = json.loads(query)
        qr = self.flow_model_cls.mng_find(query, collection, limit=self.limit)
        return qr

    def _display_qr(self, qr):
        """
        Helper function that builds the view and update `self.out_area`.
        """
        objects = []
        if qr:
            for model in qr.models:
                objects.append(model.get_panel_view(self.mng_connector))
        else:
            alert = pn.pane.Alert(f"No model found in collection for query: {qr.query}.",
                                  alert_type="danger")
            objects.append(alert)

        self.out_area.objects = objects
