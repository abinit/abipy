import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.core.structure import Structure
from abipy.flowtk.tasks import TaskManager
from abipy.htc.structure_models import StructureData
from abipy.htc.pseudos_models import PseudoDojoSpecs
from abipy.htc.gs_models import GsData
from abipy.htc.flow_models import EbandsFlowModel


class TestAbipyBaseModels(AbipyTest):

    def test_ebands_flow_model(self):
        #mongo_connector = MongoConnector.for_localhost(collection_name="foobar")
        #collection = mongo_connector.get_collection()
        #collection.drop()
        #print("Building collection....")

        si = Structure.from_file(abidata.ref_file("refs/si_ebands/run.abi"))
        pseudos_specs = PseudoDojoSpecs.from_table_name("Foo")

        structures = [si]
        models = []
        for structure in structures:
            input_structure_data = StructureData.from_structure(structure)
            kppa = 300
            model = EbandsFlowModel(input_structure_data=input_structure_data,
                                    pseudos_specs=pseudos_specs, kppa=kppa)
            models.append(model)
            assert model.flow_data.exec_status == "Initialized"
            assert model.abipy_version and model.pymatgen_version
            assert model.kppa == kppa

            workdir = self.mkdtemp()
            manager = TaskManager.from_user_config()
            flow = model.build_flow(workdir, manager)

            #oid = model.mongo_insert(collection)
            #assert oid

            model.scf_data = GsData.from_gsr_filepath(abidata.ref_file("si_scf_GSR.nc"))
            model.nscf_data = GsData.from_gsr_filepath(abidata.ref_file("si_nscf_GSR.nc"))
            #d = model.dict()
            #for k, v in d.items(): print(f"{type(k)} --> {type(v)}")

            #oid = model.mongo_insert(collection)
            #model.mongo_full_update_oid(oid, collection)

            d = model.as_dict()
            #print(d, type(d["structure"]))
            ##print(monty_json_dumps(model))
            ##new_model = EbandsFlowModel.from_dict(d)
            ##sys.exit(1)

            #same_model = EbandsFlowModel.from_mongo_oid(oid, collection)
            #print("same_model", type(same_model)) #, model)
            #print("same_model.structure", type(same_model.structure)) #, model)
            #assert same_model.flow_status == 0
            #assert same_model.structure == structure
            #assert isinstance(same_model.scf_data.ebands, ElectronBands)
            #assert isinstance(same_model.nscf_data.ebands, ElectronBands)
            #print(same_model.scf_data.ebands.structure)

        ##mongo_insert_models(models, collection)
        #print("OK")

        ##items = EbandsFlowModel.find_runnable_oid_models(collection, limit=1)
        ##if items:
        ##    for item in items:
        ##        print(item)

        ##for model_cls in [MongoConnector, PseudoDojoSpecs,]: # EbandsFlowModel.
        ##    print(model_cls.schema_json(indent=2))
        ##sys.exit(1)
