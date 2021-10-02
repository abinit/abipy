"""Tests for test_flow_models.py module"""
import mongomock
import abipy.data as abidata

from pymongo.collection import Collection
from abipy.core.testing import AbipyTest
from abipy.core.structure import Structure
from abipy.electrons.ebands import ElectronBands
from abipy.flowtk import TaskManager, Flow
from abipy.htc.base_models import mongo_insert_models, MockedMongoConnector
from abipy.htc.structure_models import StructureData
from abipy.htc.pseudos_models import PseudoSpecs
from abipy.htc.gs_models import GsData, NscfData
from abipy.htc.flow_models import FlowModel, ExecStatus
from abipy.htc.gs_flow_models import EbandsFlowModelWithParams, EbandsFlowModelWithInputs


class TestFlowModels(AbipyTest):

    def test_ebands_flow_model(self):
        """Testing EbandsFlowModelWithParams."""
        collection_name = "dummy_collection_name"
        mocked_connector = MockedMongoConnector(host="example.com", port=27017, collection_name=collection_name)
        pseudos_specs = PseudoSpecs.from_repo_table_name("ONCVPSP-PBEsol-SR-PDv0.4", "standard")

        si = Structure.from_file(abidata.ref_file("refs/si_ebands/run.abi"))
        collection = mocked_connector.get_collection()

        # This should raise as the collection does not contain the magic document with the FlowModel class.
        with self.assertRaises(RuntimeError):
            EbandsFlowModelWithParams.get_subclass_from_collection(collection)

        assert collection.count_documents({}) == 0

        oid = EbandsFlowModelWithParams.init_collection(collection)
        assert EbandsFlowModelWithParams.init_collection(collection) == oid
        assert collection.count_documents({}) == 1
        sub_class = EbandsFlowModelWithParams.get_subclass_from_collection(collection)
        assert sub_class is EbandsFlowModelWithParams

        with self.assertRaises(TypeError):
            # This should raise as we are trying to register another model in the same collection.
            FlowModel.init_collection(collection)

        with self.assertRaises(TypeError):

            class DifferentFlowModel(FlowModel):
                """This should raise as the class in the collection is not the same as the base FlowModel class."""

            DifferentFlowModel.get_subclass_from_collection(collection)

        # but one can still use the base class.
        assert FlowModel.get_subclass_from_collection(collection) is EbandsFlowModelWithParams

        items = EbandsFlowModelWithParams.find_runnable_oid_models(collection)
        assert not items

        model_list = []
        structures = [si]
        for structure in structures:
            in_structure_data = StructureData.from_structure(structure)
            kppa = 300
            model = EbandsFlowModelWithParams(in_structure_data=in_structure_data,
                                              pseudos_specs=pseudos_specs, kppa=kppa)
            model_list.append(model)
            assert model.is_init
            assert not model.is_built
            assert not model.is_completed
            assert not model.is_errored
            assert model.abipy_version and model.pymatgen_version
            assert model.kppa == kppa
            assert model.scf_data is None
            assert model.nscf_kpath_data is None

            oid = model.mongo_insert(collection)
            assert oid

            workdir = self.mkdtemp()
            manager = TaskManager.from_user_config()
            flow = model.build_flow(workdir, manager)
            #flow = model.build_flow_and_update_collection(workdir, oid, collection, worker)
            #assert model.is_built
            assert isinstance(flow, Flow)
            assert flow.workdir == workdir
            assert flow[0][0].manager is not None

            collection_name = "dummy_collection_name"
            mocked_connector = MockedMongoConnector(host="example.com", port=27017, collection_name=collection_name)

            with_gsr = True
            #with_gsr = False

            model.scf_data = GsData.from_gsr_filepath(abidata.ref_file("si_scf_GSR.nc"), mocked_connector, with_gsr)
            #print(model.scf_data.json())
            #assert 0
            model.nscf_kpath_data = NscfData.from_gsr_filepath(abidata.ref_file("si_nscf_GSR.nc"),
                                                               mocked_connector, with_gsr)

            # Save the updated model
            model.mongo_full_update_oid(oid, collection)

            # Now retrieve the same model from the collection.
            same_model = EbandsFlowModelWithParams.from_mongo_oid(oid, collection)
            assert isinstance(same_model.scf_data, GsData)
            assert model.scf_data.pressure_gpa == float(same_model.scf_data.pressure_gpa)
            assert isinstance(same_model.scf_data.ebands, ElectronBands)
            assert isinstance(same_model.nscf_kpath_data.ebands, ElectronBands)

            #d = model.as_dict()
            #print(d, type(d["structure"]))
            ##print(monty_json_dumps(model))
            ##new_model = EbandsFlowModelWithParams.from_dict(d)
            #print("same_model", type(same_model)) #, model)
            #print("same_model.structure", type(same_model.structure)) #, model)
            #assert same_model.flow_status == 0
            #assert same_model.structure == structure
            #print(same_model.scf_data.ebands.structure)

        systems = EbandsFlowModelWithParams.mongo_get_crystal_systems_incoll(collection)
        assert systems == ["Cubic"]
        spg_numbers = EbandsFlowModelWithParams.mongo_get_spg_numbers_incoll(collection)
        assert spg_numbers == [227]

        qr = EbandsFlowModelWithParams.mongo_find_by_formula("Si", collection)
        assert qr and qr.models[0].in_structure_data.structure.composition.reduced_formula == "Si"

        qr = EbandsFlowModelWithParams.mongo_find_by_spg_number(1, collection)
        assert not qr
        qr = EbandsFlowModelWithParams.mongo_find_by_spg_number(227, collection)
        assert qr and qr.models[0].in_structure_data.spg_number == 227

        qr = EbandsFlowModelWithParams.mongo_find_by_crystal_system("Cubic", collection)
        assert qr and qr.models[0].in_structure_data.crystal_system == "Cubic"

        oid_models = EbandsFlowModelWithParams.find_runnable_oid_models(collection, limit=1)
        assert len(oid_models) == 1
        #runnable_oid, runnable_model = oid_models[0]
        #flow = runnable_model.build_flow()
        #assert runnable_model.exec_status == ""

        status2oids = EbandsFlowModelWithParams.mongo_get_status2oids(collection)
        assert not status2oids[ExecStatus.errored]
        assert not status2oids[ExecStatus.completed]
        assert len(status2oids[ExecStatus.init]) == collection.count_documents({}) - 1

        oids = mongo_insert_models(model_list, collection, verbose=1)
        assert len(oids) == len(model_list)

        #import tempfile
        #import os
        #with tempfile.TemporaryDirectory() as tmp:
        #    filepath = os.path.join(tmp, "structure_data.json")
        #    structure_data.json_write(filepath, indent=4)
        #    same_data = StructureData.from_json_file(filepath)
        #    assert type(structure_data) is type(same_data)
        #    assert same_data.structure == structure_data.structure

    def test_ebands_flow_model_with_inputs(self):
        """Testing EbandsFlowModeWithInputs"""

        scf_input = nscf_input = self.get_gsinput_si()
        in_structure_data = StructureData.from_structure(scf_input.structure)
        pseudos_specs = PseudoSpecs.from_repo_table_name("ONCVPSP-PBEsol-SR-PDv0.4", "standard")

        with self.assertRaises(ValueError):
            EbandsFlowModelWithInputs(in_structure_data=in_structure_data, pseudos_specs=pseudos_specs)

        model = EbandsFlowModelWithInputs(scf_input=scf_input, nscf_input=nscf_input,
                                          in_structure_data=in_structure_data, pseudos_specs=pseudos_specs)

        collection: Collection = mongomock.MongoClient().db.collection
        oid = model.mongo_insert(collection)
        assert oid
