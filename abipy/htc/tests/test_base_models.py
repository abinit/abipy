"""
Tests for base_models module.
"""
#from pprint import pprint
import os
import mongomock

from typing import List
from pydantic.main import ModelMetaclass
from pymatgen.core.structure import Structure as pmg_Structure
from abipy.core.testing import AbipyTest
from abipy.core.structure import Structure as abi_Structure
from abipy.abio.inputs import AbinitInput
from abipy.htc.base_models import (AbipyModel, MongoConnector, MockedMongoConnector, QueryResults, GridFsDesc,
        MongoModel, mongo_insert_models)


class SubModel(AbipyModel):
    abi_structure: abi_Structure


class TopModel(MongoModel):
    pmg_structure: pmg_Structure
    submodel: SubModel
    class_type: ModelMetaclass = SubModel
    items: List[int] = [1, 2, 3]
    abinit_input: AbinitInput


class TestAbipyBaseModels(AbipyTest):

    def test_base_models(self):

        pmg_structure = self.get_structure("Si")  # This is a pymatgen structure.
        abi_structure = abi_Structure.as_structure(self.get_structure("Si"))  # This is an Abipy structure.

        sub_model = SubModel(abi_structure=abi_structure)
        abinit_input = self.get_gsinput_si()
        top_model = TopModel(pmg_structure=pmg_structure, submodel=sub_model,
                             abinit_input=abinit_input, items=[4.2, 5.3, 6.0])


        # Note how pydantic has converted items to int.
        assert top_model.items == [4, 5, 6]

        # Test pydantic API
        pydantic_json_string = top_model.json()
        assert "TopModel" not in pydantic_json_string

        monty_json_string = top_model.to_json()
        assert "TopModel" in monty_json_string

        monty_dict = top_model.as_dict()
        same_top_model = TopModel.from_dict(monty_dict)

        assert same_top_model.class_type is SubModel

        assert isinstance(pmg_structure, pmg_Structure)
        assert same_top_model.pmg_structure == pmg_structure
        assert isinstance(same_top_model.submodel.abi_structure, abi_Structure)
        assert same_top_model.submodel.abi_structure == pmg_structure
        #assert same_top_model.abinit_input == abinit_input
        #self.assertMSONable(top_model, test_if_subclass=True)
        #assert 0

        collection = mongomock.MongoClient().db.collection
        oid = top_model.mongo_insert(collection)
        assert oid
        same_model = TopModel.from_mongo_oid(oid, collection)
        assert same_model.pmg_structure == pmg_structure
        same_model.mongo_full_update_oid(oid, collection)

        models = [top_model, same_model]
        oids = mongo_insert_models(models, collection)
        assert len(oids) == 2

        query = {}
        qr = TopModel.mongo_find(query, collection)
        assert bool(qr) and len(qr) == 3
        assert len(qr.models) == 3
        assert set(qr.oids) == set(oids + [oid])

    def test_mongo_connector_base_api(self):
        connector = MongoConnector(host="example.com", port=27017, collection_name="collection_name")
        assert connector.host == "example.com"
        assert str(connector.port) in str(connector)

        connector = MongoConnector.for_localhost(collection_name="foobar")
        assert "foobar" in str(connector)
        self.assertMSONable(connector, test_if_subclass=True)

        # Username requires password.
        with self.assertRaises(ValueError):
            MongoConnector(host="example.com", username="John", collection_name="collection_name")

        connector = MongoConnector(host="example.com", username="John",
                                   password="password", collection_name="collection_name")
        assert repr(connector)
        assert str(connector)

    def test_query_results_api(self):
        query = {"_id": "foo"}
        qr = QueryResults.empty_from_query(query)
        assert qr.query == query
        assert not qr
        assert len(qr) == 0

    def test_mocked_mongo_connector(self):
        """Testing MockedMongoConnector"""
        collection_name = "dummy_collection_name"
        mocked_connector = MockedMongoConnector(host="example.com", port=27017, collection_name=collection_name)
        assert str(mocked_connector)
        assert repr(mocked_connector)
        collection = mocked_connector.get_collection()
        assert collection
        #print(collection)
        assert collection.insert_one({"foo": 1, "bar": "foo"})

        fs, fs_collname = mocked_connector.get_gridfs_and_name()

        assert fs_collname == f"gridfs_{collection_name}"
        buf = b"hello world"
        oid = fs.put(buf)
        assert fs.get(oid).read() == buf

        # It seems that mongomock doesn't support db.list_collection_names
        #print("names:", mocked_connector.list_collection_names())
        #assert collection_name in mocked_connector.list_collection_names()
        #assert 0

        # Testing GridFsDesc
        import abipy.data as abidata
        test_dir = os.path.join(os.path.dirname(__file__), "..", "..", 'test_files')
        ddb_filepath = os.path.join(test_dir, "AlAs_444_nobecs_DDB")
        assert ddb_filepath

        # Insert the DDB in GridFs.
        gfsd = mocked_connector.gridfs_put_filepath(ddb_filepath)
        assert gfsd.oid is not None
        assert gfsd.filepath == ddb_filepath
        assert gfsd.collection_name == fs_collname
        assert str(gfsd)
        assert gfsd.json()
        assert fs.exists({"filename": os.path.basename(ddb_filepath)})
        assert fs.exists({"filepath": ddb_filepath})
        assert fs.exists({"parent_collection": collection_name})

        # Extract the DDB from GridFs and open it with abiopen.
        with mocked_connector.abiopen_gfsd(gfsd) as ddb:
            assert ddb.structure

        gsr_filepath = abidata.ref_file("si_scf_GSR.nc")
        gfsd = mocked_connector.gridfs_put_filepath(gsr_filepath)
        with mocked_connector.abiopen_gfsd(gfsd) as gsr:
            assert gsr.ebands.structure
