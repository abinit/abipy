"""Tests for htc.FilesFile."""
from __future__ import print_function, division, unicode_literals

from abipy.core.testing import AbipyTest
from abipy.abio.abivars_db import get_abinit_variables, abinit_help, docvar


class AbinitVariableDatabaseTest(AbipyTest):

    def test_database(self):
        """Testing database of ABINIT variables."""
        database = get_abinit_variables()
        assert database is get_abinit_variables()

        # The text of this variable contaings greek symbols in HTML.
        var = database["cd_frqim_method"]
        print(var)

        # Print all variables in the database.
        for name, var in database.items():
            #print("testing variable: ", name)
            assert var.name == name
            repr(var)
            str(var)
            str(var.info)
            str(var._repr_html_())

        # Database methods.
        database.apropos("ecut")
        assert len(database.json_dumps_varnames())

        print("vargeo section:\n", database.vars_with_section("vargeo"))
        for section in database.sections:
            assert len(database.vars_with_section(section))

        for charact in database.characteristics:
            print("character:", charact)
            assert len(database.vars_with_char(charact))

        name2section = database.name2section
        assert name2section["ecut"] == "varbas" and name2section["ionmov"] == "varrlx"

        assert database.group_by_section("ecut") == {'varbas': ['ecut']}

        assert not database["ecut"].isarray
        assert database["spinat"].isarray

        abinit_help("ecut", info=True)
        # Should not raise
        abinit_help("foobar", info=True)

        ecut = docvar("ecut")
        assert ecut.name == "ecut"
