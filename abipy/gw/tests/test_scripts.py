from __future__ import division, print_function, unicode_literals

__author__ = 'setten'

import unittest
import os
import shutil
import subprocess

import tempfile
import abipy.data as abidata
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

        self.assert_equal(spec_in.test(), 0)

        wdir = tempfile.mkdtemp()
        base = os.getcwd()
        print(wdir)

        os.chdir(wdir)

        shutil.copyfile(abidata.cif_file("si.cif"), os.path.join(wdir, 'si.cif'))
        shutil.copyfile(abidata.pseudo("14si.pspnc").path, os.path.join(wdir, 'Si.pspnc'))
        shutil.copyfile(os.path.join(abidata.dirpath, 'managers', 'shell_manager.yml'), os.path.join(wdir, 'manager.yml'))
        shutil.copyfile(os.path.join(abidata.dirpath, 'managers', 'scheduler.yml'), os.path.join(wdir, 'scheduler.yml'))
        spec_in.data['source'] = 'cif'
        print(abidata.dirpath)

        spec_in.write_to_file('spec.in')
        spec_in.loop_structures('i')

        os.chdir(base)
        os.remove(wdir)

        # test the folder structure

        # read some of the input files
