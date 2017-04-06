"""Tests for phonons"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import numpy as np

from abipy.core.testing import *
from abipy.dfpt.anaddbnc import AnaddbNcFile
import abipy.data as abidata


class AnaddbNcFileTest(AbipyTest):

    def test_base(self):
        """Base tests for AnaddbNcFile"""
        anaddbnc_fname = abidata.ref_file("AlAs_nl_dte_anaddb.nc")

        with AnaddbNcFile(anaddbnc_fname) as anc:
            self.assertIsNotNone(anc.becs)
            self.assertIsNotNone(anc.emacro)
            self.assertIsNotNone(anc.emacro_rlx)
            self.assertIsNotNone(anc.dchidt)
            self.assertIsNotNone(anc.dchide)
            self.assertIsNotNone(anc.oscillator_strength)
