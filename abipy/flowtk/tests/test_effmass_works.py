"""Tests for abiphonopy module"""
import abipy.data as abidata
import abipy.flowtk as flowtk
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.abio.factories import gs_input
from abipy.flowtk.effmass_works import EffMassLineWork, EffMassDFPTWork, EffMassAutoDFPTWork, FrohlichZPRFlow


class TestEffMassWorks(AbipyTest):

    def test_effmassline_work(self):
        """Testing EffMassLineWork."""
        si_structure = abidata.structure_from_cif("si.cif")

        scf_input = gs_input(si_structure, pseudos=abidata.pseudos("14si.pspnc"), ecut=4, spin_mode="unpolarized")
        flow = flowtk.Flow.temporary_flow()

        with self.assertRaises(ValueError):
            work = EffMassLineWork.from_scf_input(scf_input, k0_list=(0, 0, 0), step=0.01, npts=2)

        # Run SCF from scratch.
        work = EffMassLineWork.from_scf_input(scf_input, k0_list=(0, 0, 0), step=0.01, npts=9)

        t0, t1 = work[0], work[1]
        assert t1.depends_on(t0)
        assert t1.input["kptopt"] == 0
        assert t1.input["iscf"] == -2

        flow.register_work(work)
        flow.allocate()
        flow.check_status()
        isok, checks = flow.abivalidate_inputs()
        assert isok

        # From DEN file
        flow = flowtk.Flow.temporary_flow()
        work = EffMassLineWork.from_scf_input(scf_input, k0_list=(0, 0, 0), step=0.01, npts=10,
                                              red_dirs=(1, 0, 0),
                                              cart_dirs=[(1, 0, 0), (1, 1, 0)],
                                              den_node=abidata.ref_file("si_DEN.nc"))
        flow.register_work(work)
        flow.allocate()
        flow.check_status()
        isok, checks = flow.abivalidate_inputs()
        assert isok

    def test_effmass_dfpt_work(self):
        """Testing EffMassDFPTWork."""
        si_structure = abidata.structure_from_cif("si.cif")

        scf_input = gs_input(si_structure, pseudos=abidata.pseudos("14si.pspnc"), ecut=4, spin_mode="unpolarized")
        flow = flowtk.Flow.temporary_flow()

        # Run SCF from scratch.
        work = EffMassDFPTWork.from_scf_input(scf_input, k0_list=(0, 0, 0), effmass_bands_f90=[2, 4])

        assert len(work) == 3
        t0, t1, t2 = work
        assert t1.input["iscf"] == -2
        assert t2.input["kptopt"] == 0
        assert t2.input["efmas"] == 1 and t2.input["rfelfd"] == 2
        assert t2.depends_on(t1)

        flow.register_work(work)
        flow.allocate()
        flow.check_status()
        isok, checks = flow.abivalidate_inputs()
        assert isok

    def test_effmass_autodfpt_work(self):
        """Testing EffMassAutoDFPTWork."""
        si_structure = abidata.structure_from_cif("si.cif")

        scf_input = gs_input(si_structure, pseudos=abidata.pseudos("14si.pspnc"), ecut=4, spin_mode="unpolarized")
        flow = flowtk.Flow.temporary_flow()

        # Run SCF from scratch.
        work = EffMassAutoDFPTWork.from_scf_input(scf_input)

        assert len(work) == 2
        assert work[1].depends_on(work[0])

        flow.register_work(work)
        flow.allocate()
        flow.check_status()
        isok, checks = flow.abivalidate_inputs()
        assert isok

    def test_frohlich_zpr_flow(self):
        """Testing FrohlichZPRFlow"""
        # Read structure from DDB file.
        from abipy import abilab
        ddb_path = abidata.ref_file("refs/mgo_v8t57/mgo_zpr_t57o_DS3_DDB")
        with abilab.abiopen(ddb_path) as ddb:
            structure = ddb.structure

        pseudos = abidata.pseudos("Ca.psp8", "O.psp8")
        scf_input = abilab.AbinitInput(structure=structure, pseudos=pseudos)
        # Set other input variables. These quantities are system-depedent.
        # Here we use parameters similar to https://docs.abinit.org/tests/v8/Input/t57.in
        scf_input.set_vars(
            nband=12,
            nbdbuf=2,
            diemac=6,
            ecut=30,                # Underconverged ecut.
            #ecut=15,
            nstep=100,
            tolvrs=1e-16,
            kptrlatt=[-2,  2,  2,   # In cartesian coordinates, this grid is simple cubic
                       2, -2,  2,
                       2,  2, -2],
        )

        workdir = self.mkdtemp()
        flow = FrohlichZPRFlow.from_scf_input(workdir, scf_input, ddb_node=ddb_path, ndivsm=2, tolwfr=1e-20,
                                              metadata={"mp_id": "mp-149"})
        flow.allocate()
        flow.check_status()
        isok, checks = flow.abivalidate_inputs()
        assert isok
