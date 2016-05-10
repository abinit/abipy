from __future__ import division, print_function, unicode_literals
import os
import shutil
import tempfile
import abipy.data as abidata
from pymatgen.util.testing import PymatgenTest
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.GWvaspinputsets import GWDFTDiagVaspInputSet, GWG0W0VaspInputSet, GWscDFTPrepVaspInputSet
from pymatgen.io.vasp.GWvaspinputsets import SingleVaspGWWork
from abipy.gw.datastructures import GWSpecs, get_spec  # , GWConvergenceData
from abipy.gw.codeinterfaces import AbinitInterface, VaspInterface, get_code_interface
from abipy.gw.tests.test_helpers import structure
from pymatgen.io.vasp.inputs import get_potcar_dir
POTCAR_DIR = get_potcar_dir()


__author__ = 'setten'


class GWSetupTest(PymatgenTest):
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
