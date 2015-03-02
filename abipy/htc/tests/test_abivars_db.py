"""Tests for htc.FilesFile."""
from __future__ import print_function, division, unicode_literals

from abipy.htc.abivars_db import get_abinit_variables

def test_database():
    """Testing database of ABINIT variables."""
    database = get_abinit_variables()

    # The text of this variable contaings greek symbols in HTML.
    var = database["cd_frqim_method"]
    print(var)

    for name, var in database.items():
        print("testing variable: ", name)
        print(repr(var))
        print(var)
        print(var.info)
        print(var._repr_html_())

    # Database methods.
    database.apropos("ecut")
    print(database.json_dumps_varnames())

    for section in database.sections:
        database.vars_with_section(section)
