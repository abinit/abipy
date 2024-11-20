# coding: utf-8
"""Test abilab module."""
import os
import abipy.data as abidata
import json

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
        assert d and "pymatgen" in d
        df = abilab.software_stack(as_dataframe=True)
        assert df is not None

        filepath = self.get_tmpname(text=True, suffix=".json")
        data = {"foo": "bar"}
        abilab.mjson_write(data, filepath, indent=4)
        data = abilab.mjson_load(filepath)
        assert data["foo"] == "bar"

        same_data = abilab.mjson_loads(json.dumps(data))
        assert len(same_data) == len(data)
        assert same_data["foo"] == data["foo"]

        assert not abilab.abicheck(verbose=1)

        abilab.abipy_logo1()
        abilab.abipy_logo2()
        abilab.abipy_logo3()

        assert not abilab.in_notebook()
        abilab.enable_notebook(with_seaborn=True)
        assert abilab.in_notebook()
        abilab.disable_notebook()
        assert not abilab.in_notebook()

        assert abilab.install_config_files(workdir=self.mkdtemp()) == 0

    def extscls_supporint_panel(self):
        table = abilab.extcls_supporting_panel(as_table=True)
        exscls = abilab.extcls_supporting_panel(as_table=False)
        exts = [item[0] for item in extcls]
        assert "DDB" in exts
        assert ".cube" not in exts
