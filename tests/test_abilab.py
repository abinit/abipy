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
        assert len(abilab.dir2abifiles(top=abidata.dirpath, recurse=False)) == 1
        with self.assertRaises(ValueError):
            abilab.abifile_subclass_from_filename("foobar")

        assert abilab.isabifile("foo_GSR.nc")
        assert not abilab.isabifile("foobar")

        import pandas
        df = pandas.DataFrame({"a": [1, 2], "b": [3, 4]})
        abilab.print_dataframe(df, title="foo")

        d = abilab.software_stack()
        assert d

        assert not abilab.abicheck(verbose=1)

        abilab.abipy_logo1()
        abilab.abipy_logo2()
        abilab.abipy_logo3()

        assert not abilab.in_notebook()
        abilab.enable_notebook(with_seaborn=True)
        assert abilab.in_notebook()
        abilab.disable_notebook()
        assert not abilab.in_notebook()