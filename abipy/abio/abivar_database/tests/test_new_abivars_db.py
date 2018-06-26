"""Tests for htc.FilesFile."""
from __future__ import print_function, division, unicode_literals

from abipy.core.testing import AbipyTest

from abipy.abio.abivar_database.variables import get_codevars


class AbinitVariableDatabaseTest(AbipyTest):

    def test_database(self):
        """Testing database of ABINIT variables."""
        varscode = get_codevars()
        assert varscode is get_codevars()
        database = varscode["abinit"]

        # The text of this variable contaings greek symbols in HTML.
        var = database["cd_frqim_method"]
        repr(var); str(var)

        # Print all variables in the database.
        for name, var in database.items():
            #print("testing variable: ", name)
            assert var.name == name
            repr(var); str(var)
            str(var.info)
            # FIXME
            #str(var._repr_html_())

        # Database methods.
        database.apropos("ecut")
        # FIXME
        #assert len(database.json_dumps_varnames())

        for setname in [
            "basic", "rlx", "gstate", "eph", "ffield", "paral", "gw", "dfpt", 
            "geo", "bse", "dev", "paw", "dmft", "files", "internal", "w90"]:
            assert database.vars_with_varset(setname)

        for section in database.my_varset_list:
            assert len(database.vars_with_varset(section))

        for charact in database.my_characteristics:
            print("character:", charact)
            assert len(database.vars_with_char(charact))

        name2varset = database.name2varset
        assert name2varset["ecut"] == "basic" and name2varset["ionmov"] == "rlx"

        assert database.group_by_varset("ecut") == {'basic': ['ecut']}

        natom_var = database["natom"]
        ecut_var = database["ecut"]
        assert ecut_var.name == "ecut"
        assert not ecut_var.isarray
        assert not ecut_var.depends_on_dimension("natom")
        assert not ecut_var.depends_on_dimension(natom_var)
        # FIXME
        #assert ecut_var.url
        #assert ecut_var.html_link() and ecut_var.html_link(label="foo")

        spinat_var = database["spinat"]
        assert spinat_var.isarray
        assert spinat_var.depends_on_dimension("natom")
        assert spinat_var.depends_on_dimension(natom_var)
        assert not spinat_var.depends_on_dimension("ntypat")

        # FIXME
        #abinit_help("ecut", info=True)
        # Should not raise
        #abinit_help("foobar", info=True)

        #ecut_var = docvar("ecut")
        #assert ecut_var.name == "ecut"
