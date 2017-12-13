"""Tests for dfpt_works module"""
from __future__ import print_function, division, unicode_literals, absolute_import

import abipy.data as abidata
import abipy.flowtk as flowtk
#import abipy.flowtk

from abipy.core.testing import AbipyTest
#from abipy.abio.factories import gs_input
#from abipy.flowtk import mocks
from abipy.flowtk import dfpt_works


class TestDfptWorks(AbipyTest):

    def test_nscfddkswork(self):
        """Testing NscfDdksWork."""
        return
        work = dfpt_works.NscfDdksWork.from_scf_task(cls, scf_task, ddk_ngkpt, ddk_shiftk, ddk_nband, manager=None)
