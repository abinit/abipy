# coding: utf-8
"""Test manager files."""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import abipy.data as abidata
from abipy import abilab

from abipy.core.testing import AbipyTest


class AbilabTest(AbipyTest):

    def test_abilab(self):
        """Testing abilab"""
        abilab.abiopen_ext2class_table()
        assert abilab.abifile_subclass_from_filename("GSR.nc") is not None

        assert len(abilab.dir2abifiles(top=abidata.dirpath))
        with self.assertRaises(ValueError):
            abilab.abifile_subclass_from_filename("foobar")

        assert abilab.isabifile("foo_GSR.nc")
        assert not abilab.isabifile("foobar")

        d = abilab.software_stack()
        assert d

        assert not abilab.abicheck(verbose=1)

        abilab.abipy_logo1()
        abilab.abipy_logo2()
        abilab.abipy_logo3()