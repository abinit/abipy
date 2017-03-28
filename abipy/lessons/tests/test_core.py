from __future__ import unicode_literals, division, print_function

import abipy.data as abidata
import abipy.abilab as abilab
from abipy.core.testing import AbipyTest
from abipy.lessons.core import get_pseudos
from pymatgen.io.abinit.pseudos import NcAbinitPseudo


class CoreTest(AbipyTest):
    """
    Testing helper functions in core
    """

    def test_get_pseudo(self):
        """
        Testing the get_pseudo method
        """
        structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))
        pseudos = get_pseudos(structure)
        print(type(pseudos))
        self.assertIsInstance(pseudos, list)
        self.assertIsInstance(pseudos[0], NcAbinitPseudo)
        self.assertEqual(len(pseudos), 1)
