from __future__ import division, print_function, unicode_literals

__author__ = 'setten'

import unittest
import os
import subprocess

import tempfile
from pymatgen.util.testing import PymatgenTest
from pymatgen.core.structure import Structure
from pymatgen.io.vaspio.GWvaspinputsets import GWDFTDiagVaspInputSet, GWG0W0VaspInputSet, GWscDFTPrepVaspInputSet
from pymatgen.io.vaspio.GWvaspinputsets import SingleVaspGWWork
from abipy.gw.datastructures import GWSpecs, GWConvergenceData, get_spec
from abipy.gw.codeinterfaces import AbinitInterface, VaspInterface, get_code_interface
from abipy.gw.tests.test_helpers import structure

#test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",'test_files')

from pymatgen.io.vaspio.vasp_input import get_potcar_dir
POTCAR_DIR = get_potcar_dir()


class GWSetupTest(PymatgenTest):
    def test_setup(self):
        """
        testing the main functions called in the abiGWsetup script
        """

        spec_in = get_spec('GW')
        self.assertIsInstance(spec_in, GWSpecs)

        spec_in.test()

        wdir = tempfile.mkstemp()
        base = os.getcwd()

        os.chdir(wdir)

        #copy refference spec.in

        #test the main function used in abiGWsetup

        spec_in.write_to_file('spec.in')
        spec_in.loop_structures('i')

        os.chdir(base)
        os.remove(wdir)

        # test the folder structure

        # read some of the input files
