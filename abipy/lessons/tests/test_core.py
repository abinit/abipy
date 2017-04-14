from __future__ import unicode_literals, division, print_function

import os
import copy
import abipy.data as abidata
import abipy.abilab as abilab
from abipy.core.testing import AbipyTest
from abipy.lessons.core import get_pseudos, BaseLesson
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

    def test_base_lesson(self):

        a, c, p = "a", "c", "p.py"

        class Lesson(BaseLesson):
            @property
            def abipy_string(self):
                return copy.copy(a)

            @property
            def comline_string(self):
                return copy.copy(c)

            @property
            def pyfile(self):
                return copy.copy(p)

        lesson = Lesson()
        lesson.abidata.cif_file("si.cif")

        self.assertEqual(lesson.abipy_string, a)
        self.assertEqual(lesson.comline_string, c)
        self.assertEqual(lesson.pyfile, p)
        self.assertEqual(lesson.manpath.replace("u'", "'"), p.replace('py','man'))
        self.assertTrue(len(str(lesson.docvar('ecut'))) > 4)
        lesson._gen_manfile()
        os.remove(p.replace('py', 'man'))
