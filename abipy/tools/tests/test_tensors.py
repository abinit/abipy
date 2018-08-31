"""Tests for tensors module"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import numpy as np

from abipy.core.testing import AbipyTest
from abipy.tools.tensors import DielectricTensor, NLOpticalSusceptibilityTensor
import abipy.data as abidata


class DielectricTensorTest(AbipyTest):

    def test_base(self):
        """Base tests for DielectricTensor"""
        eps = DielectricTensor(np.diag([1, 2, 3]))
        repr(eps); str(eps)
        assert eps._repr_html_()
        assert len(eps.get_dataframe()) == 3
        assert len(eps.get_voigt_dataframe().keys()) == 6
        self.assertArrayAlmostEqual(eps.reflectivity(), [0., 0.029437251522859434, 0.071796769724490825])


class NLOpticalSusceptibilityTensorTest(AbipyTest):

    def test_base(self):
        """Base tests for NLOpticalSusceptibilityTensor"""
        anaddbnc_fname = abidata.ref_file("AlAs_nl_dte_anaddb.nc")

        tensor = NLOpticalSusceptibilityTensor.from_file(anaddbnc_fname)
        repr(tensor); str(tensor)