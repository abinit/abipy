"""This module contains lookup table with the name of the ABINIT variables."""
from __future__ import print_function, division

import os
import json

__all__ = [
    "is_anaddb_var",
    "is_abitoken",
    "is_abivar",
    "is_abiunit",
    "has_abiop",
]


with open(os.path.join(os.path.dirname(__file__), "anaddb_vars.json")) as fh:
    #print(fh.read())
    _anaddb_varnames = set(json.load(fh))


def is_anaddb_var(varname):
    """True if varname is an Anaddb variable."""
    return varname in _anaddb_varnames


def is_abitoken(s):
    """
    True if s is one of the token supported by the ABINIT parser
    i.e. variable name or unit name.
    """
    return s in ABI_ALLTOKENS


def is_abivar(s):
    """
    True if s is an ABINIT variable. 
    """
    return s in ABI_VARNAMES


def is_abiunit(s):
    """
    True if s is one of the units supported by the ABINIT parser
    """
    return s in ABI_UNITS


def has_abiop(s):
    """True is string contains one of the operators supported by the ABINIT parser."""
    return any(op in s for op in ABI_OPS)


def extract_abivars(f90file, pretty_print=True):
    """This routine extract the variable names from src/57_iovars/chkvars.F90."""
    # !<ABINIT_VARS>
    # list_vars=                 ' accesswff accuracy acell algalch amu angdeg atvshift autoparal awtr'
    # list_vars=trim(list_vars)//' fband fermie_nest'
    # list_logicals=' SpinPolarized '
    # !</ABINIT_VARS>

    #!Extra token, also admitted :
    #!<ABINIT_UNITS>
    # list_vars=trim(list_vars)//' au Angstr Angstrom Angstroms Bohr Bohrs eV Ha'
    # list_vars=trim(list_vars)//' Hartree Hartrees K Ry Rydberg Rydbergs T Tesla'
    #!</ABINIT_UNITS>

    #!<ABINIT_OPERATORS>
    #list_vars=trim(list_vars)//' sqrt end'
    #!</ABINIT_OPERATORS>

    with open(f90file, "r") as fh:
        lines = [l.strip() for l in fh]

    start = lines.index("!<ABINIT_VARS>")
    stop  = lines.index("!</ABINIT_VARS>", start)
    var_field = [l for l in lines[start+1:stop] if l and not l.startswith("!")]
    #print(var_field)

    start = lines.index("!<ABINIT_UNITS>")
    stop  = lines.index("!</ABINIT_UNITS>", start)
    unit_field = [l for l in lines[start+1:stop] if l and not l.startswith("!")]

    start = lines.index("!<ABINIT_OPERATORS>")
    stop  = lines.index("!</ABINIT_OPERATORS>", start)
    op_field = [l for l in lines[start+1:stop] if l and not l.startswith("!")]

    def extract_names(field):
        names = []
        for line in field:
            i = line.find("'")
            if i == -1: i = line.find('"')
            if i == -1:
                raise ValueError("Cannot find string in line %s" % line)

            s = line[i+1:]
            s = s.replace("'", "").replace('"', "")
            tokens = [t for t in s.split() if not t.startswith("!")]
            names.extend(tokens)

        return names

    var_names = extract_names(var_field)
    unit_names = extract_names(unit_field)
    op_names = extract_names(op_field)

    if pretty_print:
        from pprint import pprint
        print("# " + 80*"-")
        print("# Begin computer generated code")
        print()
        print("# Variable names.")
        print("ABI_VARNAMES = ", end="")
        pprint(var_names)
        print("")
        print("# Unit names.")
        print("ABI_UNITS = ", end="")
        pprint(unit_names)
        print("")
        print("# Operators.")
        print("ABI_OPS = ", end="")
        pprint(op_names + ["*", "/"])
        print("")
        print("# End computer generated code")
        print("# " + 80*"-")

    return var_names

# The variables below have been extracted from the file src/57_iovars/chkvars.F90. 
# See the function extract_abivars defined in this module.
# Tokens are divides in 3 classes: variable names, unit names, operators
from .abivars_db import ABI_VARNAMES, ABI_UNITS, ABI_OPS

# All tokens supported by the abinit parser.
ABI_ALLTOKENS = ABI_VARNAMES + ABI_OPS + ABI_UNITS 

# Build sets to speedup search.
ABI_ALLTOKENS = set(ABI_ALLTOKENS)
ABI_VARNAMES = set(ABI_VARNAMES)
ABI_OPS = set(ABI_OPS)
ABI_UNITS = set(ABI_UNITS) 


if __name__ == "__main__":
    import sys
    f90file = sys.argv[1]
    varnames = extract_abivars(f90file)
    sys.exit(0)
