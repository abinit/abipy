from __future__ import division, print_function, unicode_literals

import unittest
try:
    raise ImportError("No module named sets_deprecated")
except ImportError:
    raise unittest.SkipTest("Skipping all tests in test_classes due to sets_deprecated")

import os
import shutil
import tempfile
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.gw.datastructures import GWSpecs, get_spec  # , GWConvergenceData
from abipy.gw.tests.test_helpers import structure


#from pymatgen import SETTINGS
#POTCAR_DIR = SETTINGS.get("VASP_PSP_DIR")


__author__ = 'setten'


class GWSetupTest(AbipyTest):
    def test_setup(self):
        """
        Testing the main functions called in the abiGWsetup script
        """

        spec_in = get_spec('GW')
        self.assertIsInstance(spec_in, GWSpecs)

        self.assert_equal(spec_in.test(), 0)
        self.assert_equal(len(spec_in.to_dict()), 10)
        self.assert_equal(spec_in.get_code(), 'ABINIT')


        spec_in.data['source'] = 'cif'

        self.assert_equal(spec_in.hash_str(), "code:ABINIT;source:cifjobs:[u'prep', u'G0W0'];mode:ceci;functional:PBE;kp_grid_dens:500prec:m;tol:0.0001;test:False;converge:False")

        wdir = tempfile.mkdtemp()
        base = os.getcwd()
        print('wdir', wdir)

        os.chdir(wdir)

        shutil.copyfile(abidata.cif_file("si.cif"), os.path.join(wdir, 'si.cif'))
        shutil.copyfile(abidata.pseudo("14si.pspnc").path, os.path.join(wdir, 'Si.pspnc'))
        shutil.copyfile(os.path.join(abidata.dirpath, 'managers', 'shell_manager.yml'), os.path.join(wdir, 'manager.yml'))
        shutil.copyfile(os.path.join(abidata.dirpath, 'managers', 'scheduler.yml'), os.path.join(wdir, 'scheduler.yml'))

        try:
            temp_ABINIT_PS_EXT = os.environ['ABINIT_PS_EXT']
            temp_ABINIT_PS = os.environ['ABINIT_PS']
        except KeyError:
            temp_ABINIT_PS_EXT = None
            temp_ABINIT_PS = None

        os.environ['ABINIT_PS_EXT'] = '.pspnc'
        os.environ['ABINIT_PS'] = wdir

        spec_in.data['source'] = 'cif'
        print('dirpath', abidata.dirpath)

        spec_in.write_to_file('spec.in')
        self.assertTrue(os.path.isfile(os.path.join(wdir, 'spec.in')))


        # broken due to strategy refactoring
        # spec_in.loop_structures('i')

        os.chdir(base)

        shutil.rmtree(wdir)

        if temp_ABINIT_PS is not None:
            os.environ['ABINIT_PS_EXT'] = temp_ABINIT_PS_EXT
            os.environ['ABINIT_PS'] = temp_ABINIT_PS


        # test the folder structure

        # read some of the input files
