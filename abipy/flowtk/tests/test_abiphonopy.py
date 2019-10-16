"""Tests for abiphonopy module"""
import abipy.data as abidata
import abipy.flowtk as flowtk
import abipy.flowtk.abiphonopy as abiph

from abipy.abio.factories import gs_input
from abipy.core.testing import AbipyTest


class TestAbiPhonopy(AbipyTest):
    """Unit tests for abiphonopy module."""

    def test_phonopy_work(self):
        """Testing PhononWork."""
        self.skip_if_not_phonopy()
        si_structure = abidata.structure_from_cif("si.cif")
        same_structure = abiph.structure_from_atoms(abiph.atoms_from_structure(si_structure))
        assert si_structure == same_structure

        # TODO: Spin
        gsinp = gs_input(si_structure, pseudos=abidata.pseudos("14si.pspnc"), ecut=4, spin_mode="unpolarized")
        #gsinp = gs_input(si_structure, pseudos=abidata.pseudos("14si.pspnc"), ecut=4, spin_mode="polarized")
        flow = flowtk.Flow.temporary_flow()
        scdims = [2, 2, 2]
        phpy_work = abiph.PhonopyWork.from_gs_input(gsinp, scdims=scdims, phonopy_kwargs=None, displ_kwargs=None)
        flow.register_work(phpy_work)

        self.assert_equal(scdims, phpy_work.scdims)
        assert hasattr(phpy_work, "phonon")
        assert len(phpy_work.phonopy_tasks) == len(phpy_work)
        assert len(phpy_work.phonopy_tasks) == 1
        assert len(phpy_work.bec_tasks) == 0

        # Gruneisen with phonopy
        grun_work = abiph.PhonopyGruneisenWork.from_gs_input(gsinp, voldelta=0.1, scdims=scdims,
                                                             phonopy_kwargs=None, displ_kwargs=None)
        flow.register_work(grun_work)

        self.assert_equal(scdims, grun_work.scdims)

        flow.allocate()
        flow.check_status()
        isok, checks = flow.abivalidate_inputs()
        assert isok
