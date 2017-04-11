"""Tests for tensors module"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import numpy as np

from abipy.core.testing import AbipyTest
from abipy.dfpt.tensors import DielectricTensor, NLOpticalSusceptibilityTensor
import abipy.data as abidata


class NLOpticalSusceptibilityTensorTest(AbipyTest):

    def test_base(self):
        """Base tests for NLOpticalSusceptibilityTensor"""
        anaddbnc_fname = abidata.ref_file("AlAs_nl_dte_anaddb.nc")

        tensor = NLOpticalSusceptibilityTensor.from_file(anaddbnc_fname)
        repr(tensor)
        str(tensor)


class DielectricTensorTest(AbipyTest):

    def test_base(self):
        """Base tests for DielectricTensor"""
        dt = DielectricTensor(np.diag([1,2,3]))
        repr(dt)
        str(dt)

        self.assertArrayAlmostEqual(dt.reflectivity(), [0., 0.029437251522859434, 0.071796769724490825])