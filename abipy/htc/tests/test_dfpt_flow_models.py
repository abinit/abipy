from abipy.core.testing import AbipyTest
from abipy.htc.base_models import MongoConnector
from abipy.htc.pseudos_models import PseudoSpecs
from abipy.htc.structure_models import StructureData
from abipy.htc.dfpt_flow_models import PhononFlowModel #, PhononFlowModelWithInput


class TestPhononFlowModels(AbipyTest):

    def test_phonon_flow_model_with_input(self):
        #PhononFlowModelWithInput

        #mongo_connector = MongoConnector.for_localhost(collection_name="ebands")

        # Take options from ~/.abinit/abipy/config.yml configuration file.
        mongo_connector = MongoConnector.from_abipy_config(collection_name="ebands")

        collection = mongo_connector.get_collection()

         # Take options from ~/.abinit/abipy/config.yml configuration file.
        mongo_connector = MongoConnector.from_abipy_config(collection_name="ebands")
        #mongo_connector = MongoConnector.for_localhost(collection_name="ebands")

        # Pseudopotential specifications.
        # Note that we still need to specify the accuracy.
        pseudos_specs = PseudoSpecs.from_repo_name("ONCVPSP-PBE-SR-PDv0.4", table_name="standard")

        # Generate list of structures.
        from abipy.data.ucells import structure_from_ucell
        structures = [structure_from_ucell(name) for name in ("Si", "Si-shifted")]

        # Get pseudopotential tables with hints.
        #accuracy = "standard"
        pseudos = pseudos_specs.get_pseudos()

        models = []
        for i, structure in enumerate(structures):
            in_structure_data = StructureData.from_structure(structure)

            if i == 0: PhononFlowModel.init_collection(collection)
            #scf_input = make_scf_input(structure, pseudos, accuracy)
            model = PhononFlowModel(in_structure_data=in_structure_data, scf_input=scf_input,
                                    pseudos_specs=pseudos_specs,
                                    with_becs=False, with_quad=False)

            #print(model)
            models.append(model)

        #mongo_insert_models(models, collection, verbose=1)
        #mongo_connector_insert_models(models)
