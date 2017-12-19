"""Tests for dfpt_works module"""
from __future__ import print_function, division, unicode_literals, absolute_import

import abipy.data as abidata
import abipy.flowtk as flowtk
#import abipy.flowtk

from abipy.core.testing import AbipyTest
from abipy.flowtk import dfpt_works


class TestDfptWorks(AbipyTest):

    def test_nscfddkswork(self):
        """Testing NscfDdksWork."""
        scf_task = self.get_gsinput_si(as_task=True)
        work = dfpt_works.NscfDdksWork.from_scf_task(scf_task, ddk_ngkpt=[8, 8, 8],
            ddk_shiftk=[0, 0, 0], ddk_nband=10)
        assert len(work) == 4
