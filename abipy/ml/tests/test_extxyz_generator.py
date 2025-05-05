"""Tests for extxyz_generator module"""
import os
import numpy as np
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.ml.extxyz_generator import ExtxyzIOWriter


class AbimlTest(AbipyTest):

    def test_from_hist(self):
        """Testing ExtxyzIOWriter from ABINIT HIST.nc file."""
        filepaths_set = [
            [abidata.ref_file("si_scf_GSR.nc")],
            [abidata.ref_file("sic_relax_HIST.nc")],
            # TODO: Vasprun and HIST.nc
        ]

        for filepaths in filepaths_set:
            writer = ExtxyzIOWriter(filepaths)
            assert writer.to_string(verbose=1)
            print(writer)
            xyz_filepath = self.tmpfileindir("foo.xyz")
            writer.write(xyz_filepath)
            os.remove(xyz_filepath)
