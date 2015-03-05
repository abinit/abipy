"""This module contains lookup table with the name of the ABINIT variables."""
from __future__ import print_function, division, unicode_literals

import os
import json

__all__ = [
    "is_anaddb_var",
    "is_abitoken",
    "is_abivar",
    "is_abiunit",
    "has_abiop",
]

_anaddb_varnames = None

def _get_anaddb_varnames():
    global _anaddb_varnames 
    if _anaddb_varnames is not None: 
        return _anaddb_varnames

    with open(os.path.join(os.path.dirname(__file__), "anaddb_vars.json")) as fh:
        _anaddb_varnames = set(json.load(fh))
        return _anaddb_varnames


def is_anaddb_var(varname):
    """True if varname is an Anaddb variable."""
    return varname in _get_anaddb_varnames()


def is_abitoken(s):
    """
    True if s is one of the token supported by the ABINIT parser
    i.e. variable name or unit name.
    """
    return s in ABI_ALLTOKENS


def is_abivar(s):
    """True if s is an ABINIT variable."""
    return s in ABI_VARNAMES


def is_abiunit(s):
    """True if s is one of the units supported by the ABINIT parser"""
    return s in ABI_UNITS


def has_abiop(s):
    """True is string contains one of the operators supported by the ABINIT parser."""
    return any(op in s for op in ABI_OPS)


# Tokens are divides in 3 classes: variable names, unit names, operators
from .abivars_db import ABI_VARNAMES, ABI_UNITS, ABI_OPS

# All tokens supported by the abinit parser.
ABI_ALLTOKENS = ABI_VARNAMES + ABI_OPS + ABI_UNITS 

# Build sets to speedup search.
ABI_ALLTOKENS = set(ABI_ALLTOKENS)
ABI_VARNAMES = set(ABI_VARNAMES)
ABI_OPS = set(ABI_OPS)
ABI_UNITS = set(ABI_UNITS) 
