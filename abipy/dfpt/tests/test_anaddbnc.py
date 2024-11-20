"""Tests for phonons"""
import os
import numpy as np
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.dfpt.anaddbnc import AnaddbNcFile


class AnaddbNcFileTest(AbipyTest):

    def test_base(self):
        """Testing AnaddbNcFile API"""
        anaddbnc_fname = abidata.ref_file("AlAs_nl_dte_anaddb.nc")

        with AnaddbNcFile(anaddbnc_fname) as anc:
            repr(anc); str(anc)
            anc.to_string(verbose=2)
            assert anc.structure.formula == "Al1 As1"
            assert not anc.has_elastic_data
            assert not anc.has_piezoelectric_data
            assert anc.becs is not None
            assert anc.epsinf is not None
            assert anc.eps0 is not None
            assert anc.dchidt is not None
            assert anc.dchide is not None
            assert anc.oscillator_strength is not None
            assert anc.ifc is None
            assert anc.params["chneut"] == 1
            assert anc.amu[13] == 26.981539

            if self.has_nbformat():
                assert anc.write_notebook(nbpath=self.get_tmpname(text=True))


class AnaddbNcRobotTest(AbipyTest):

    def test_base(self):
        """Testing AnaddbNcRobot API"""
        from abipy import abilab
        assert abilab.AnaddbNcRobot.class_handles_filename("anaddb.nc")
        assert abilab.AnaddbNcRobot.class_handles_filename("/foo/bar/anaddb.nc")
        robot = abilab.AnaddbNcRobot()
        assert robot.EXT == "anaddb"
        repr(robot); str(robot)
        assert robot.to_string(verbose=2)

        if self.has_nbformat():
            assert robot.write_notebook(nbpath=self.get_tmpname(text=True))
