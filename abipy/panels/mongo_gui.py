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

    def __init__(self, mongo_connector: MongoConnector, flow_model_cls: Type[FlowModel]):
        super().__init__()
        self.mongo_connector = mongo_connector
        self.flow_model_cls = flow_model_cls

        self.limit = 0
        self.build_ui()

    #def collection(self):
    #    return self.mongo_connector.get_collection()

    def build_ui(self):
        # This is the part that should be abstracted out
        # so that one can implement customized subclasses
        # that rely on model.get_panel_view

        self.filter_query = pnw.LiteralInput(name='MongoDB filter (python dict)',
                                             placeholder="e.g. {'key': value}", type=dict,
                                             value={'kppa': 300})

        self.filter_btn = pnw.Button(name="Run query", button_type='primary')
        self.filter_btn.on_click(self.on_filter_btn)

        collection = self.mongo_connector.get_collection()

        # Get all spg_numbers in the collection.
        spg_numbers = self.flow_model_cls.mongo_get_spg_numbers_incoll(collection)
        self.spg_number_wdg = pnw.Select(name='Filter by space group number', options=spg_numbers)
        #self.spg_number_wdg = pnw.IntInput(name='Filter by space group number', value=None, step=1, start=1, end=230)

        self.filter_by_spg_number_btn = pnw.Button(name="Run query", button_type='primary')
        self.filter_by_spg_number_btn.on_click(self.on_filter_by_spg_number_btn)

        self.formula_wdg = pnw.TextInput(name='Filter by reduced formula', placeholder='Enter e.g. Si or MgO...')
        self.filter_by_formula_btn = pnw.Button(name="Run query", button_type='primary')
        self.filter_by_formula_btn.on_click(self.on_filter_by_formula_btn)

        crystal_systems = self.flow_model_cls.mongo_get_crystal_systems_incoll(collection)
        self.crystal_system_wdg = pnw.Select(name='Filter by crystal system', options=crystal_systems)
        self.filter_by_crystal_system_btn = pnw.Button(name="Run query", button_type='primary')
        self.filter_by_crystal_system_btn.on_click(self.on_filter_by_crystal_system_btn)

        self.common_queries_btn = pnw.Button(name="Run query", button_type='primary')
        self.common_queries_btn.on_click(self.on_common_queries_btn)
        options = []
        if hasattr(self.flow_model_cls, "get_common_queries"):
            options = [json.dumps(q) for q in self.flow_model_cls.get_common_queries()]
        self.common_queries_wdg = pnw.Select(name="Common MongoDB queries", options=options)

        md = pn.pane.Markdown("## Perform MongoDB queries using one of the widgets below:")

        widget_box = pn.WidgetBox(
            md,
            self.filter_query, self.filter_btn,
            pn.Row(self.spg_number_wdg, self.formula_wdg),
            pn.Row(self.filter_by_spg_number_btn, self.filter_by_formula_btn),
            self.crystal_system_wdg, self.filter_by_crystal_system_btn,
            self.common_queries_wdg, self.common_queries_btn,
            #sizing_mode="stretch_width",
        )

        self.out_area = pn.Column(sizing_mode="stretch_width")

        self.main = pn.Column(self.mongo_connector,
                              widget_box,
                              #pn.layout.Divider(),
                              #"## Results:",
                              self.out_area,
                              sizing_mode="stretch_width")

    def get_app(self):
        template = pn.template.FastListTemplate
        #template = pn.template.BootstrapTemplate
        #template = pn.template.ReactTemplate
        return template(main=self.main, title="AbiPyFlowModel MongoDB GUI")

    #@depends_on_btn_click('filter_btn', show_shared_wdg_warning=False)
    def on_filter_btn(self, _):
        with Loading(self.out_area):
            collection = self.mongo_connector.get_collection()
            query = self.filter_query.value
            qr = self.flow_model_cls.mongo_find(query, collection, limit=self.limit)
            self._display_qr(qr)

    @display_qr
    def on_filter_by_spg_number_btn(self, _):
        spg_number = self.spg_number_wdg.value
        collection = self.mongo_connector.get_collection()
        qr = self.flow_model_cls.mongo_find_by_spg_number(spg_number, collection, limit=self.limit)
        return qr

    @display_qr
    def on_filter_by_formula_btn(self, _):
        reduced_formula = self.formula_wdg.value
        collection = self.mongo_connector.get_collection()
        qr = self.flow_model_cls.mongo_find_by_formula(reduced_formula, collection, limit=self.limit)
        return qr

    @display_qr
    def on_filter_by_crystal_system_btn(self, _):
        crystal_system = self.crystal_system_wdg.value
        collection = self.mongo_connector.get_collection()
        qr = self.flow_model_cls.mongo_find_by_crystal_system(crystal_system, collection, limit=self.limit)
        return qr

    @display_qr
    def on_common_queries_btn(self, _):
        query = self.common_queries_wdg.value
        collection = self.mongo_connector.get_collection()
        if query is None:
            return QueryResults.empty_from_query({})

        query = json.loads(query)
        qr = self.flow_model_cls.mongo_find(query, collection, limit=self.limit)
        return qr

    def _display_qr(self, qr):
        objects = []
        if qr:
            for model in qr.models:
                objects.append(model.get_panel_view())
        else:
            alert = pn.pane.Alert(f"No model found in collection for query: {qr.query}.",
                                  alert_type="danger")
            objects.append(alert)

        self.out_area.objects = objects


#if __name__ == "__main__":
