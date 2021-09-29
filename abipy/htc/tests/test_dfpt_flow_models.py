from abipy.core.testing import AbipyTest
from abipy.htc.base_models import MockedMongoConnector
from abipy.htc.pseudos_models import PseudoSpecs
from abipy.htc.structure_models import StructureData
from abipy.htc.dfpt_flow_models import PhononFlowModelWithParams, PhononFlowModelWithInput



def make_scf_input(structure, pseudos, accuracy, paral_kgb=0):
    """
    This function constructs the input file for the GS calculation:
    """
    # Crystalline AlAs: computation of the second derivative of the total energy
    from abipy import abilab
    #structure = abidata.structure_from_ucell("AlAs")
    gs_inp = abilab.AbinitInput(structure, pseudos=pseudos)

    #num_valence_electrons = gs_inp.structure.num_valence_electrons(pseudos)
    num_valence_electrons = gs_inp.num_valence_electrons
    nsppol = 2
    nband = num_valence_electrons // nsppol + 4
    nband += nband % 2

    gs_inp.set_vars(
        nband=nband,
        #nband=4,
        #ecut=2.0,
        ngkpt=[4, 4, 4],
        #nshiftk=4,
        #shiftk=[0.0, 0.0, 0.5,   # This gives the usual fcc Monkhorst-Pack grid
        #        0.0, 0.5, 0.0,
        #        0.5, 0.0, 0.0,
        #        0.5, 0.5, 0.5],
        shiftk=[0, 0, 0],
        paral_kgb=paral_kgb,
        tolvrs=1.0e-8,
        diemac=9.0,
    )

    gs_inp.set_cutoffs_for_accuracy(accuracy)
    #gs_inp.set_auto_scf_nband(nsppol=1, nspinor=1, nspden=1, occopt, tsmear)
    return gs_inp



class TestPhononFlowModels(AbipyTest):

    def test_phonon_flow_model_with_input(self):

        mongo_connector = MockedMongoConnector.for_localhost(collection_name="phbands")
        collection = mongo_connector.get_collection()

        # Pseudopotential specifications.
        pseudos_specs = PseudoSpecs.from_repo_table_name("ONCVPSP-PBE-SR-PDv0.4", "standard")
        # Get pseudopotential table with hints.
        pseudos = pseudos_specs.get_pseudos()

        # Generate list of structures.
        from abipy.data.ucells import structure_from_ucell
        structures = [structure_from_ucell(name) for name in ("Si", "Si-shifted")]

        models = []
        for i, structure in enumerate(structures):
            in_structure_data = StructureData.from_structure(structure)

            if i == 0: PhononFlowModelWithParams.init_collection(collection)
            scf_input = make_scf_input(structure, pseudos, accuracy="normal")

            #model = PhononFlowModelWithParams(in_structure_data=in_structure_data,
            #                        scf_input=scf_input,
            #                        pseudos_specs=pseudos_specs,
            #                        with_becs=False, with_quad=False)

            #PhononFlowModelWithInput

            #print(model)
            #models.append(model)

        #mongo_insert_models(models, collection, verbose=1)
        #mongo_connector_insert_models(models)
