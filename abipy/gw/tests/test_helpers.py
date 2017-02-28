from __future__ import division, print_function, unicode_literals

import unittest
try:
    raise ImportError("No module named sets_deprecated")
except ImportError:
    raise unittest.SkipTest("Skipping all tests in test_classes due to sets_deprecated")



import os
import tempfile
import shutil

from abipy.abilab import Structure as AbiStructure
from pymatgen.util.testing import PymatgenTest
from pymatgen.core.structure import Structure
from pymatgen.io.abinit.helpers import clean, read_extra_abivars, expand
import abipy.data as abidata
from abipy.gw.datastructures import get_spec
from abipy.gw.GWworks import GWG0W0VaspInputSet, SingleAbinitGWWork

__author__ = 'setten'

have_abinit_ps_ext = True
try:
    os.environ['ABINIT_PS_EXT']
except KeyError:
    have_abinit_ps_ext = False


structure_dict = {'lattice': {'a': 3.866974622849504,
                              'gamma': 60.0000000241681,
                              'c': 3.86697462,
                              'b': 3.866974623775052,
                              'matrix': [[3.34889826, 0.0, 1.93348731], [1.11629942, 3.15737156, 1.93348731], [0.0, 0.0, 3.86697462]],
                              'volume': 40.888291888494884,
                              'alpha': 60.000000032293386,
                              'beta': 60.00000002437586},
                  'sites': [{'properties': {u'coordination_no': 5, u'forces': [0.0, 0.0, 0.0]},
                             'abc': [0.875, 0.875, 0.875],
                             'xyz': [3.9070479700000003, 2.762700115, 6.767205585],
                             'species': [{'occu': 1, 'element': 'Si'}], 'label': 'Si'},
                            {'properties': {u'coordination_no': 5, u'forces': [0.0, 0.0, 0.0]},
                             'abc': [0.125, 0.125, 0.125], 'xyz': [0.55814971, 0.394671445, 0.966743655],
                             'species': [{'occu': 1, 'element': 'Si'}], 'label': 'Si'}],
                  '@class': 'Structure', '@module': 'pymatgen.core.structure'}
structure = Structure.from_dict(structure_dict)


class GWTestHelpers(PymatgenTest):

    def test_clean(self):
        """
        Testing helper function fo string cleaning
        """
        string = 'MddmmdDD  '
        self.assertEqual(clean(string), 'mddmmddd')

    def test_read_extra_abivars(self):
        """
        Testing helper function to read extra variables
        """
        vars_out = {'ecut': 40}
        f = open('extra_abivars', 'w')
        f.write(str(vars_out))
        f.close()
        vars_in = read_extra_abivars()
        self.assertEqual(vars_out, vars_in)
        os.remove('extra_abivars')

    def test_expand(self):
        """
        Testing helper function to extend the convergence grid
        """
        self.maxDiff = None
        spec = get_spec('GW')
        struc = AbiStructure.from_file(abidata.cif_file("si.cif"))
        struc.item = 'test'

        wdir = tempfile.mkdtemp()
        print('wdir', wdir)

        os.chdir(wdir)
        shutil.copyfile(abidata.cif_file("si.cif"), os.path.join(wdir, 'si.cif'))
        shutil.copyfile(abidata.pseudo("14si.pspnc").path, os.path.join(wdir, 'Si.pspnc'))
        shutil.copyfile(os.path.join(abidata.dirpath, 'managers', 'shell_manager.yml'),
                        os.path.join(wdir, 'manager.yml'))
        shutil.copyfile(os.path.join(abidata.dirpath, 'managers', 'scheduler.yml'), os.path.join(wdir, 'scheduler.yml'))

        try:
            temp_ABINIT_PS_EXT = os.environ['ABINIT_PS_EXT']
            temp_ABINIT_PS = os.environ['ABINIT_PS']
        except KeyError:
            temp_ABINIT_PS_EXT = None
            temp_ABINIT_PS = None

        os.environ['ABINIT_PS_EXT'] = '.pspnc'
        os.environ['ABINIT_PS'] = wdir

        tests = SingleAbinitGWWork(struc, spec).convs
        tests_out = {'nscf_nbands': {'test_range': (40,),
                                     'control': 'gap', 'method': 'set_bands', 'level': 'nscf'},
                     'ecut': {'test_range': (50, 48, 46, 44),
                              'control': 'e_ks_max', 'method': 'direct', 'level': 'scf'},
                     'ecuteps': {'test_range': (4, 8, 12, 16, 20),
                                 'control': 'gap', 'method': 'direct', 'level': 'sigma'}}
        self.assertEqual(expand(tests, 1), tests_out)
        spec.data['code'] = 'VASP'

#        if "VASP_PSP_DIR" in os.environ:
#            spec.update_code_interface()
#            tests = GWG0W0VaspInputSet(struc, spec).convs
#            tests_out = {'ENCUTGW': {'test_range': (200, 400, 600, 800), 'control': 'gap', 'method': 'incar_settings'}}
#            self.assertEqual(expand(tests, 1), tests_out)

        if temp_ABINIT_PS is not None:
            os.environ['ABINIT_PS_EXT'] = temp_ABINIT_PS_EXT
            os.environ['ABINIT_PS'] = temp_ABINIT_PS
