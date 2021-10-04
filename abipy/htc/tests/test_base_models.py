"""
Tests for base_models module.
"""
#from pprint import pprint
import os
import zlib
import mongomock

from typing import List
from pydantic.main import ModelMetaclass
from pymatgen.core.structure import Structure as pmg_Structure
from abipy.core.testing import AbipyTest
from abipy.core.structure import Structure as abi_Structure
from abipy.abio.inputs import AbinitInput
from abipy.htc.base_models import (AbipyModel, MongoConnector, MockedMongoConnector, QueryResults, GfsFileDesc,
        TopLevelModel, mng_insert_models)


class SubModel(AbipyModel):
    abi_structure: abi_Structure


class Top(TopLevelModel):
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
        top_model = Top(pmg_structure=pmg_structure, submodel=sub_model,
                             abinit_input=abinit_input, items=[4.2, 5.3, 6.0])


        # Note how pydantic has converted items to int.
        assert top_model.items == [4, 5, 6]

        # Test pydantic API
        pydantic_json_string = top_model.json()
        assert "Top" not in pydantic_json_string

        monty_json_string = top_model.to_json()
        assert "Top" in monty_json_string

        monty_dict = top_model.as_dict()
        same_top_model = Top.from_dict(monty_dict)

        assert same_top_model.class_type is SubModel

        assert isinstance(pmg_structure, pmg_Structure)
        assert same_top_model.pmg_structure == pmg_structure
        assert isinstance(same_top_model.submodel.abi_structure, abi_Structure)
        assert same_top_model.submodel.abi_structure == pmg_structure
        #assert same_top_model.abinit_input == abinit_input
        #self.assertMSONable(top_model, test_if_subclass=True)
        #assert 0

        collection = mongomock.MongoClient().db.collection
        oid = top_model.mng_insert(collection)
        assert oid
        same_model = Top.from_oid(oid, collection)
        assert same_model.pmg_structure == pmg_structure
        same_model.mng_full_update_oid(oid, collection)

        models = [top_model, same_model]
        oids = mng_insert_models(models, collection)
        assert len(oids) == 2

        query = {}
        qr = Top.mng_find(query, collection)
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

        new_conn = connector.new_with_collection_name("foobar")
        assert new_conn.host == connector.host
        assert new_conn.collection_name == "foobar"

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
        assert collection.insert_one({"foo": 1, "bar": "foo"})

        fs, fs_collname = mocked_connector.get_gridfs_and_name()

        assert fs_collname == f"{collection_name}_gridfs"
        buf = b"hello world"
        oid = fs.put(buf)
        assert fs.get(oid).read() == buf

        # It seems that mongomock doesn't support db.list_collection_names
        #print("names:", mocked_connector.list_collection_names())
        #assert collection_name in mocked_connector.list_collection_names()
        #assert 0

        # Testing GfsFileDesc
        import abipy.data as abidata
        test_dir = os.path.join(os.path.dirname(__file__), "..", "..", 'test_files')
        ddb_filepath = os.path.join(test_dir, "AlAs_444_nobecs_DDB")
        assert ddb_filepath

        # Insert the DDB in GridFs.
        for zlib_level in [zlib.Z_DEFAULT_COMPRESSION, zlib.Z_NO_COMPRESSION]:
            gfsd = mocked_connector.gfs_put_filepath(ddb_filepath, zlib_level=zlib_level)
            assert gfsd.oid is not None
            assert gfsd.filepath == ddb_filepath
            assert gfsd.collection_name == fs_collname
            assert gfsd.zlib_level == zlib_level
            assert str(gfsd)
            assert gfsd.json()
            assert fs.exists({"filename": os.path.basename(ddb_filepath)})
            assert fs.exists({"filepath": ddb_filepath})

            # Extract the DDB from GridFs and open it with abiopen.
            with mocked_connector.abiopen_gfsd(gfsd) as ddb:
                assert ddb.structure

        # Insert a netcdf file.
        gsr_filepath = abidata.ref_file("si_scf_GSR.nc")
        for zlib_level in [zlib.Z_DEFAULT_COMPRESSION, zlib.Z_NO_COMPRESSION]:
            gfsd = mocked_connector.gfs_put_filepath(gsr_filepath, zlib_level=zlib_level)
            assert gfsd.oid is not None
            assert gfsd.zlib_level == zlib_level
            with mocked_connector.abiopen_gfsd(gfsd) as gsr:
                assert gsr.ebands.structure

        # Insert a MSONable object
        for zlib_level in [zlib.Z_DEFAULT_COMPRESSION, zlib.Z_NO_COMPRESSION]:
            gfsd = mocked_connector.gfs_put_mson_obj(gsr.ebands, zlib_level=zlib_level)
            assert gfsd.zlib_level == zlib_level
            saved_ebands = mocked_connector.gfs_get_mson_obj(gfsd)
            assert saved_ebands.nsppol == gsr.ebands.nsppol
            assert saved_ebands.structure.formula == gsr.ebands.structure.formula

        #assert 0

        # Test drop_collection
        from unittest.mock import patch
        with patch('abipy.tools.iotools.get_input', return_value='no'):
            assert not mocked_connector.drop_collection(ask_for_confirmation=True)

        with patch('abipy.tools.iotools.get_input', return_value='yes'):
            assert mocked_connector.drop_collection(ask_for_confirmation=True)

        assert mocked_connector.drop_collection(ask_for_confirmation=False)
