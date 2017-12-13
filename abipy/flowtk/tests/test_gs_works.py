"""Tests for gs_works module"""
from __future__ import print_function, division, unicode_literals, absolute_import

import abipy.data as abidata
import abipy.flowtk as flowtk

from abipy.core.testing import AbipyTest
#from abipy.abio.factories import gs_input
from abipy.flowtk import gs_works
#from abipy.flowtk import mocks


class TestGsWorks(AbipyTest):

    def test_eoswork(self):
        """Testing EosWork."""
        return
        work = gs_works.EosWork.from_scf_input(scf_input, npoints=4, deltap_vol=0.25, ecutsm=0.5, move_atoms=True)

        #mocks.change_task_start(gstask, mocked_status="Error")
        #assert gstask.start() == 1 and gstask.status == gstask.S_ERROR
