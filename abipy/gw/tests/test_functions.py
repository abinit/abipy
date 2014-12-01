from __future__ import division, print_function, unicode_literals

__author__ = 'setten'

from pymatgen.util.testing import PymatgenTest
from abipy.gw.codeinterfaces import AbinitInterface, NewCodeInterface, VaspInterface
from abipy.gw.codeinterfaces import get_all_ecuteps, get_all_nbands, CODE_CLASSES


class GWFunctionsTest(PymatgenTest):

    def test_get_all_ecuteps(self):
        self.assertEqual(set(get_all_ecuteps()), set(['ecuteps', 'ENCUTGW', 'new_code_nbands']))

    def test_get_all_nbands(self):
        self.assertEqual(set(get_all_nbands()), set(['nscf_nbands', 'NBANDS', 'new_code_nbands']))


class GWConstantsTest(PymatgenTest):

    def test_code_classes(self):
        self.assertEqual(CODE_CLASSES, {'VASP': VaspInterface, 'ABINIT': AbinitInterface, 'NEW_CODE': NewCodeInterface})
