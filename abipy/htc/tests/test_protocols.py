from abipy.core.testing import AbipyTest
from abipy.htc.protocols import ProtocolParser  # Protocol,


class TestProtocol(AbipyTest):

    def test_base_api(self):

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

fast:
    description: "Fast protocol"
    cutoff_stringency: "low"

    gs_scf_specs:
      meta_params:
        kppa: 500
      abivars:
        nshiftk: 1
        shiftk: [[0.0, 0.0, 0.0]]
        tolvrs: 1.0e-7
        tsmear: 0.008 # Ha

moderate:
    description: "Moderate protocol"
    cutoff_stringency: "normal"

    gs_scf_specs:
      meta_params:
        kppa: 1000
      abivars:
        #kpoints_distance: 0.20
        nshiftk: 1
        shiftk: [[0.0, 0.0, 0.0]]
        tolvrs: 1.0e-9
        tsmear: 0.008  # Ha
"""

        p = ProtocolParser.from_yaml_string(yaml_string)

        fast = p.protocols["fast"]
        print(fast)
        assert fast.description == "Fast protocol"
        assert fast.cutoff_stringency == "low"
        gs_scf = fast.gs_scf_specs
        assert gs_scf.abivars["rmm_diis"] == 1
        assert gs_scf.abivars["tolvrs"] == 1.0e-7
        assert gs_scf.meta_params["kppa"] == 500

        moderate = p.protocols["moderate"]
        print(moderate)
        assert moderate.cutoff_stringency == "normal"
        assert moderate.description == "Moderate protocol"
        gs_scf = moderate.gs_scf_specs
        assert gs_scf.abivars["rmm_diis"] == 1
        assert gs_scf.abivars["tolvrs"] == 1.0e-9
        assert gs_scf.meta_params["kppa"] == 1000
