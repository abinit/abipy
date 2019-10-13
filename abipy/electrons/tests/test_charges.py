# coding: utf-8
"""Tests for charges."""
import numpy as np
import os
import abipy.data as abidata

from abipy import abilab
from abipy.core.testing import AbipyTest
from abipy.electrons.denpot import DensityNcFile
from abipy.electrons.charges import HirshfeldCharges, BaderCharges


class HirshfeldChargesTest(AbipyTest):

    def test_hirshfeld_charges(self):
        si_structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))
        hc = HirshfeldCharges([0,0], si_structure)
        hc.structure


class BaderChargesTest(AbipyTest):

    def test_bader_charges(self):
        si_structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))
        bc = BaderCharges([-4,-4], si_structure, [-4,-4])
        bc.structure
        self.assertArrayAlmostEqual(bc.net_charges, [0,0])

    def test_bader_generation(self):
        if not self.which("bader"):
            raise self.SkipTest("bader missing")

        # pseudodojo required to get the core charges
        self.skip_if_not_pseudodojo()

        import pseudo_dojo
        si_pseudo_path = os.path.join(pseudo_dojo.dojotable_absdir("ONCVPSP-PBE-PDv0.4"),"Si", "Si.psp8")
        den_path = abidata.ref_file("si_DEN.nc")

        bc = BaderCharges.from_files(den_path, {"Si": si_pseudo_path},
                                     maxr=1., small_dist_mesh=(5, 5, 5), small_dist_factor=1.)

        # check with very low accuracy, since a different version of Bader or of the DEN file may give different results
        self.assertArrayAlmostEqual(bc.net_charges, [0.,0.], decimal=0)

