"""Tests for abiphonopy module"""
import abipy.data as abidata
import abipy.flowtk as flowtk

from abipy.core.testing import AbipyTest
from abipy.abio.factories import gs_input
from abipy.flowtk.eph_flows import GkqPathFlow


class TestEphFlows(AbipyTest):
    """Unit tests for eph_flows module."""

    def test_gkqpathflows(self):
        """Testing GkqPathFlow."""
        ngkpt = [4, 4, 4]

        structure = abidata.structure_from_ucell("AlAs")
        pseudos = abidata.pseudos("13al.981214.fhi", "33as.pspnc")
        from abipy.abio.inputs import AbinitInput
        gs_inp = AbinitInput(structure, pseudos=pseudos)

        gs_inp.set_vars(
            nband=5,
            ecut=8.0,
            ngkpt=ngkpt,
            nshiftk=1,
            shiftk=[0, 0, 0],
            tolvrs=1.0e-6,
            diemac=12.0,
        )

        ngqpt = (2, 2, 2)
        qbounds = [[0.01, 0, 0], [0.02, 0, 0], [0.04, 0, 0], [0.05, 0, 0]]
        qbounds = [[0, 0, 0], [0.5, 0, 0]]

        workdir = self.mkdtemp()
        flow = GkqPathFlow.from_scf_input(workdir, gs_inp, ngqpt, qbounds, ndivsm=2, with_becs=True, ddk_tolerance=None,
                                          test_ft_interpolation=True, prepgkk=0)

        #assert len(work.relax_tasks) == 3
        #assert all(t.input["optcell"] == 3 for t in work)
        #assert all(t.input["ionmov"] == 3 for t in work)
        #assert all(t.input["ecutsm"] == 0.5 for t in work)
        #assert all(t.input["dilatmx"] == 1.05 for t in work)

        flow.allocate()
        flow.check_status()
        isok, checks = flow.abivalidate_inputs()
        assert isok

        #with self.assertRaises(ValueError):
        #    GruneisenWork.from_gs_input(gs_inp, voldelta, ngqpt=[3, 3, 3], with_becs=False)
