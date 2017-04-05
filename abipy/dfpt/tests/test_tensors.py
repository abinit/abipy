"""Tests for phonons"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import numpy as np

from abipy.core.testing import *
from abipy.dfpt.tensors import DielectricTensor, NLOpticalSusceptibilityTensor
from abipy.dfpt.anaddbnc import AnaddbNcFile
from abipy.dfpt.phonons import PhononBands
import abipy.data as abidata


class NLOpticalSusceptibilityTensorTest(AbipyTest):

    def test_base(self):
        """Base tests for NLOpticalSusceptibilityTensor"""
        anaddbnc_fname = abidata.ref_file("AlAs_nl_dte_anaddb.nc")

        NLOpticalSusceptibilityTensor.from_file(anaddbnc_fname)


class DielectricTensorTest(AbipyTest):

    def test_base(self):
        """Base tests for DielectricTensor"""
        anaddbnc_fname = abidata.ref_file("AlAs_nl_dte_anaddb.nc")
        phbstnc_fname = abidata.ref_file("AlAs_nl_dte_PHBST.nc")

        d = DielectricTensor.from_files(phbstnc_fname, anaddbnc_fname)

        self.assertAlmostEqual(d.tensor_at_frequency(0.001)[0,0], 11.917178775812721)

        d = DielectricTensor.from_objects(PhononBands.from_file(phbstnc_fname), AnaddbNcFile.from_file(anaddbnc_fname))

        self.assertAlmostEqual(d.tensor_at_frequency(0.001)[0,0], 11.917178775812721)