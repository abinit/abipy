import numpy as np

from pprint import pprint
from typing import List
from pydantic.main import ModelMetaclass
from pymatgen.core.structure import Structure as pmg_Structure
from abipy.core.testing import AbipyTest
from abipy.core.structure import Structure as abi_Structure
from abipy.htc.base_models import AbipyBaseModel, MongoConnector


class TestAbipyBaseModels(AbipyTest):

    def test_base_models(self):

        class SubModel(AbipyBaseModel):
            abi_structure: abi_Structure

        class TopModel(AbipyBaseModel):
            pmg_structure: pmg_Structure
            submodel: SubModel
            class_type: ModelMetaclass = SubModel
            items: List[int] = [1, 2, 3]
            #items: np.ndarray = np.array([1, 2, 3])

        # This is a pymatgen structure.
        pmg_structure = self.get_structure("Si")
        # This is an Abipy structure.
        abi_structure = abi_Structure.as_structure(self.get_structure("Si"))

        sub_model = SubModel(abi_structure=abi_structure)
        top_model = TopModel(pmg_structure=pmg_structure, submodel=sub_model)

        # Test pydantic API
        pydantic_json_string = top_model.json()
        assert "TopModel" not in pydantic_json_string

        monty_json_string = top_model.to_json()
        print(monty_json_string)
        assert "TopModel" in monty_json_string
        #assert 0

        monty_dict = top_model.as_dict()
        pprint(monty_dict)
        same_top_model = TopModel.from_dict(monty_dict)

        assert same_top_model.class_type is SubModel
        assert isinstance(pmg_structure, pmg_Structure)
        assert same_top_model.pmg_structure == pmg_structure
        assert isinstance(same_top_model.submodel.abi_structure, abi_Structure)
        assert same_top_model.submodel.abi_structure == pmg_structure
        self.assertMSONable(top_model, test_if_subclass=True)
        #assert 0

    def test_mongo_connector(self):
        c = MongoConnector(host="localhost", port=27017, collection_name="collection_name")
        assert c.host == "localhost"
        assert str(c.port) in c._repr_markdown_()

        c = MongoConnector.for_localhost(collection_name="foobar")
        assert "foobar" in c._repr_markdown_()
        self.assertMSONable(c, test_if_subclass=True)

        #collection = mongo_connector.get_collection()
        #mongo_connector.open_mongoflow_gui(**serve_kwargs)
        #client = c.get_client()
        #client = c.get_collection()
