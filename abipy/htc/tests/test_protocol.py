from abipy.core.testing import AbipyTest
from abipy.htc.protocol import GsScfSpecs, Protocol


class TestProtocol(AbipyTest):

    def test_specs(self):
        with self.assertRaises(ValueError):
            # Invalid abinit variable
            GsScfSpecs(meta_params={}, abivars={"ecut_foo": 10})

        with self.assertRaises(ValueError):
            # kppa_foo is not supported by GsScfSpecs
            GsScfSpecs(meta_params={"kppa_foo": 1000}, abivars={"ecut": 10})

        with self.assertRaises(ValueError):
            # kppa and ngkpt are mutually exclusive
            GsScfSpecs(meta_params={"kppa": 1000}, abivars={"ngkpt": [10, 10, 10]})

        specs = GsScfSpecs(meta_params={"kppa": 1000}, abivars={"ecut": 10})
        assert specs.abivars["ecut"] == 10

    def test_protocol_api(self):
        """Testing Protocol API."""

        all_protocols = Protocol.get_all_abipy_protocols()
        assert len(all_protocols) > 0
        proto = Protocol.from_name("NC-PBE-SR-PDv0.4.yml")
        assert str(proto)

        yaml_string = """
# Variables common to all the protocols
global_abivars:
  fband: 2.00  # increase the number of bands to > fband * natoms
  nstep: 100
  paral_kgb: 1
  rmm_diis: 1
  expert_user: 1

pseudo_specs:
  repo_name: "ONCVPSP-PBEsol-SR-PDv0.4"
  table_name: "standard"

info: "Fast protocol"
#cutoff_stringency: "low"

gs_scf_specs:
  meta_params:
    kppa: 500
  abivars:
    nshiftk: 1
    shiftk: [[0.0, 0.0, 0.0]]
    tolvrs: 1.0e-7
    tsmear: 0.008 # Ha

#moderate:
#    info: "Moderate protocol"
#    cutoff_stringency: "normal"
#
#    gs_scf_specs:
#      meta_params:
#        kppa: 1000
#      abivars:
#        #kpoints_distance: 0.20
#        nshiftk: 1
#        shiftk: [[0.0, 0.0, 0.0]]
#        tolvrs: 1.0e-9
#        tsmear: 0.008  # Ha
"""

        proto = Protocol.from_yaml_string(yaml_string)

        assert proto.info == "Fast protocol"
        #assert proto.cutoff_stringency == "low"
        specs = proto.gs_scf_specs
        assert specs.abivars["rmm_diis"] == 1
        assert specs.abivars["tolvrs"] == 1.0e-7
        assert specs.meta_params["kppa"] == 500

        #moderate = proto.protocols["moderate"]
        #print(moderate)
        #assert moderate.cutoff_stringency == "normal"
        #assert moderate.info == "Moderate protocol"
        #gs_scf = moderate.gs_scf_specs
        #assert gs_scf.abivars["rmm_diis"] == 1
        #assert gs_scf.abivars["tolvrs"] == 1.0e-9
        #assert gs_scf.meta_params["kppa"] == 1000

        #from abipy.data.ucells import structure_from_ucell
        #structures = [structure_from_ucell(name) for name in ("Si",)] # "Si-shifted")]
        #scf_inp = proto.get_gs_scf_input(structure)
        #self.abivalidate_inp(scf_inp)

        #with self.assertRaises(ValueError)
        #proto.get_ebands_input(structure)
