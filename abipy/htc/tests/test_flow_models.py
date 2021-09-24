"""Tests for test_flow_models.py module"""
import mongomock
import abipy.data as abidata

from pymongo.collection import Collection
from abipy.core.testing import AbipyTest
from abipy.core.structure import Structure
from abipy.electrons.ebands import ElectronBands
from abipy.flowtk import TaskManager, Flow
from abipy.htc.base_models import mongo_insert_models
from abipy.htc.structure_models import StructureData
from abipy.htc.pseudos_models import PseudoSpecs
from abipy.htc.gs_models import GsData
from abipy.htc.flow_models import EbandsFlowModel


class TestFlowModels(AbipyTest):

    def test_ebands_flow_model(self):
        """Testing EbandsFlowModel."""
        si = Structure.from_file(abidata.ref_file("refs/si_ebands/run.abi"))
        pseudos_specs = PseudoSpecs.from_repo_name("ONCVPSP-PBEsol-SR-PDv0.4")
        collection: Collection = mongomock.MongoClient().db.collection

        with self.assertRaises(RuntimeError):
            EbandsFlowModel.get_subclass_from_collection(collection)

        assert collection.count_documents({}) == 0

        EbandsFlowModel.init_collection(collection)
        assert collection.count_documents({}) == 1
        sub_class = EbandsFlowModel.get_subclass_from_collection(collection)
        assert sub_class is EbandsFlowModel

        items = EbandsFlowModel.find_runnable_oid_models(collection)
        assert not items

        model_list = []
        structures = [si]
        for structure in structures:
            input_structure_data = StructureData.from_structure(structure)
            kppa = 300
            model = EbandsFlowModel(input_structure_data=input_structure_data,
                                    pseudos_specs=pseudos_specs, kppa=kppa)
            model_list.append(model)
            #assert model.flow_data.exec_status  == "Initialized"
            assert model.abipy_version and model.pymatgen_version
            assert model.kppa == kppa
            assert model.scf_data is None
            assert model.nscf_kpath_data is None

            oid = model.mongo_insert(collection)
            assert oid

            workdir = self.mkdtemp()
            manager = TaskManager.from_user_config()
            flow = model.build_flow(workdir, manager)
            assert isinstance(flow, Flow)
            assert flow.workdir == workdir
            assert flow[0][0].manager is not None

            #d = model.dict()
            #for k, v in d.items(): print(f"{type(k)} --> {type(v)}")
            model.scf_data = GsData.from_gsr_filepath(abidata.ref_file("si_scf_GSR.nc"))
            model.nscf_kpath_data = GsData.from_gsr_filepath(abidata.ref_file("si_nscf_GSR.nc"))

            # Save the updated model
            model.mongo_full_update_oid(oid, collection)

            # Now retrieve the same model from the collection.
            same_model = EbandsFlowModel.from_mongo_oid(oid, collection)
            assert isinstance(same_model.scf_data, GsData)
            assert model.scf_data.pressure_gpa == float(same_model.scf_data.pressure_gpa)
            assert isinstance(same_model.scf_data.ebands, ElectronBands)
            assert isinstance(same_model.nscf_kpath_data.ebands, ElectronBands)

            #d = model.as_dict()
            #print(d, type(d["structure"]))
            ##print(monty_json_dumps(model))
            ##new_model = EbandsFlowModel.from_dict(d)
            #print("same_model", type(same_model)) #, model)
            #print("same_model.structure", type(same_model.structure)) #, model)
            #assert same_model.flow_status == 0
            #assert same_model.structure == structure
            #print(same_model.scf_data.ebands.structure)

        #qr = EbandsFlowModel.mongo_find_by_formula("Si", collection)
        #assert qr
        #qr = EbandsFlowModel.mongo_find_by_spg_number(1, collection)
        #assert not qr
        #qr = EbandsFlowModel.mongo_find_by_spg_number(32, collection)
        #assert not qr

        oid_models = EbandsFlowModel.find_runnable_oid_models(collection, limit=1)
        assert len(oid_models) == 1
        #runnable_oid, runnable_model = oid_models[0]
        #flow = runnable_model.build_flow()
        #assert runnable_model.exec_status == ""

        oids = mongo_insert_models(model_list, collection)
        assert len(oids) == len(model_list)

        #import tempfile
        #import os
        #with tempfile.TemporaryDirectory() as tmp:
        #    filepath = os.path.join(tmp, "structure_data.json")
        #    structure_data.json_write(filepath, indent=4)
        #    same_data = StructureData.from_json_file(filepath)
        #    assert type(structure_data) is type(same_data)
        #    assert same_data.structure == structure_data.structure
