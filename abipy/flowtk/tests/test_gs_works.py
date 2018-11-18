"""Tests for gs_works module"""
from __future__ import print_function, division, unicode_literals, absolute_import

import abipy.data as abidata
import abipy.flowtk as flowtk

from abipy.core.testing import AbipyTest
from abipy.flowtk import gs_works


class TestGsWorks(AbipyTest):

    def test_eos_work(self):
        """Testing EosWork."""
        scf_input = self.get_gsinput_si()
        work = gs_works.EosWork.from_scf_input(scf_input, npoints=4, deltap_vol=0.25, ecutsm=2.0, move_atoms=True)
        assert len(work) == 9
        assert all(isinstance(task, flowtk.RelaxTask) for task in work)
        assert all(task.input["ecutsm"] == 2.0 for task in work)
        assert all(task.input["ionmov"] == 2 for task in work)

        work = gs_works.EosWork.from_scf_input(scf_input, npoints=3, deltap_vol=0.25, ecutsm=0.5, move_atoms=False)
        assert len(work) == 7
        assert all(isinstance(task, flowtk.ScfTask) for task in work)
        self.abivalidate_work(work)
