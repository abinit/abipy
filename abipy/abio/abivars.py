"""This module contains lookup table with the name of the ABINIT variables."""
from __future__ import division, print_function, unicode_literals, absolute_import

import json
import numpy as np

from monty.string import is_string
from monty.functools import lazy_property
from pymatgen.core.units import bohr_to_ang
from abipy.core.structure import Structure, frame_from_structures
from abipy.core.mixins import Has_Structure

import logging
logger = logging.getLogger(__name__)

__all__ = [
    "AbinitInputFile",
    "is_anaddb_var",
    "is_abivar",
    "is_abiunit",
]

_anaddb_varnames = None

def _get_anaddb_varnames():
    global _anaddb_varnames
    if _anaddb_varnames is not None:
        return _anaddb_varnames

    from abipy import data as abidata
    with open(abidata.var_file("anaddb_vars.json")) as fh:
        _anaddb_varnames = set(json.load(fh))
        return _anaddb_varnames


def is_anaddb_var(varname):
    """True if varname is a valid Anaddb variable."""
    return varname in _get_anaddb_varnames()


ABI_VARNAMES = None

def is_abivar(s):
    """True if s is an ABINIT variable."""
    global ABI_VARNAMES
    if ABI_VARNAMES is None:
        from abipy import data as abidata
        with open(abidata.var_file("abinit_vars.json")) as fh:
            ABI_VARNAMES = set(json.load(fh))

    return s in ABI_VARNAMES


ABI_UNIT_NAMES = {s.lower() for s in (
"au",
"Angstr", "Angstrom", "Angstroms", "Bohr", "Bohrs",
"eV", "Ha", "Hartree", "Hartrees", "K", "Ry", "Rydberg", "Rydbergs",
"T", "Tesla")}


def is_abiunit(s):
    """True if string is one of the units supported by the ABINIT parser"""
    if not is_string(s): return False
    global ABI_UNIT_NAMES
    return s.lower() in ABI_UNIT_NAMES


ABI_OPERATORS = set(["sqrt", ])

def is_abioperator(s):
    """True is string contains one of the operators supported by the ABINIT parser."""
    global ABI_UNIT_NAMES
    return s in ABI_OPERATORS

ABI_MULTI_TAGS = set(["+", ":", "*", "?"])


def analyze_token(tok):
    """
    Analyze a string.
    Return (True, dataset_index) or (None, None).
    """
    name = None
    isvar = False
    isunit = False
    dataset = None
    postfix = None

    if tok[0].isalpha():
        # We have a new variable or a string giving the unit.
        dtidx = None

        #unit = None
        #if is_abiunit(tok):
        #    unit = tok

        if tok[-1].isdigit() and "?" not in tok:
            # Handle dataset index.
            l = []
            for i, c in enumerate(tok[::-1]):
                if c.isalpha(): break
                l.append(c)
            else:
                raise ValueError("Cannot find dataset index in %s" % tok)

            l.reverse()
            dataset = int("".join(l))
            tok = tok[:len(tok)-i-1]

    #return dict2namedtuple(name, isvar, isunit, dataset, postfix)


def expand_star_syntax(s):
    """
    Evaluate star syntax. Return new string

    >>> assert expand_star_syntax("3*2") == '2 2 2'
    >>> assert expand_star_syntax("2 *1") == '1 1'
    >>> assert expand_star_syntax("1 2*2") == '1 2 2'
    >>> assert expand_star_syntax("*2") == '* 2'
    """
    if "*" not in s: return s
    s = s.replace("*", " * ").strip()
    tokens = s.split()

    # Handle "*2" case i.e. return "* 2"
    if len(tokens) == 2 and tokens[0] == "*":
        assert tokens[1] != "*"
        return " ".join(tokens)

    #print(s, tokens)
    l = []
    while tokens:
        c = tokens.pop(0)
        if c == "*":
            num = int(l.pop(-1))
            val = tokens.pop(0)
            l.extend(num * [val])
        else:
            l.append(c)

    return " ".join(l)


def str2array_bohr(obj):
    if not is_string(obj): return np.asarray(obj)

    # Treat e.g. acell 3 * 1.0
    obj = expand_star_syntax(obj)
    tokens = obj.split()

    if not tokens[-1].isalpha():
        # No unit
        return np.fromstring(obj, sep=" ")

    unit = tokens[-1].lower()
    if unit in ("angstr", "angstrom", "angstroms"):
        return np.fromstring(" ".join(tokens[:-1]), sep=" ") / bohr_to_ang
    elif unit in ("bohr", "bohrs"):
        return np.fromstring(" ".join(tokens[:-1]), sep=" ")
    else:
        raise ValueError("Don't know how to handle unit: %s" % unit)


def str2array(obj):
    if not is_string(obj): return np.asarray(obj)
    return np.fromstring(expand_star_syntax(obj), sep=" ")


def varname_dtindex(tok):
    """
    >>> assert varname_dtindex("acell1") == ("acell", 1)
    >>> assert varname_dtindex("fa1k2") == ("fa1k", 2)
    """
    l = []
    for i, c in enumerate(tok[::-1]):
        if c.isalpha(): break
        l.append(c)
    else:
        raise ValueError("Cannot find dataset index in %s" % tok)

    assert i > 0
    l.reverse()
    dtidx = int("".join(l))
    tok = tok[:len(tok)-i]

    return tok, dtidx


def eval_abinit_operators(tokens):
    """
    Receive a list of strings, find the occurences of operators supported
    in the input file (e.g. sqrt), evalute the expression and return new list of strings.
    """
    import math
    import re
    re_sqrt = re.compile("[+|-]?sqrt\((.+)\)")

    values = []
    for tok in tokens:
        m = re_sqrt.match(tok)
        if m:
            #print("in sqrt with token:", tok)
            tok = tok.replace("sqrt", "math.sqrt")
            tok = str(eval(tok))
        values.append(tok)
    return values


class Dataset(dict, Has_Structure):

    @lazy_property
    def structure(self):
        kwargs = {}
        if "angdeg" in self:
            assert "rprim" not in self
            raise NotImplementedError("angdeg")
            #kwargs["angdeg"] =
        else:
            # Handle structure specified with rprim.
            kwargs["rprim"] = str2array_bohr(self.get("rprim", "1.0 0 0 0 1 0 0 0 1"))

        # Default value for acell
        acell = str2array_bohr(self.get("acell", "1.0 1.0 1.0"))

        check = {k: 1 for k in ("xred", "xcart", "xangst") if k in self}
        if len(check) != 1:
            raise ValueError("Atomic positions are not specified correctly:\n%s" % str(check))

        if "xred" in self:
            kwargs["xred"] = np.fromstring(self["xred"], sep=" ")
        elif "xcart" in self:
            kwargs["xcart"] = str2array_bohr(self["xcart"])
        elif "xangst" in self:
            kwargs["xangst"] = np.fromstring(self["xangst"], sep=" ")

        #print(kwargs["rprim"])
        return Structure.from_abivars(
            acell=acell,
            znucl=str2array(self["znucl"]),
            typat=str2array(self["typat"]),
            **kwargs
        )


class AbinitInputFile(Has_Structure):
    """
    This object parses the Abinit input file, stores the variables in
    dict-like objects (Datasets) and build `Structure` objects from
    the input variables. Mainly used for inspecting the structure
    declared in the Abinit input file.
    """
    @classmethod
    def from_file(cls, path):
        """Build the object from file."""
        with open(path, "rt") as fh:
            return cls.from_string(fh.read())

    @classmethod
    def from_string(cls, string):
        """Build the object from string."""
        d = parse_abinit_string(string)
        return cls(d, string)

    def __init__(self, dvars, string):
        """
        Args:
            dvars: Dictionary {varname: value} returned by `parse_abinit_string`.
            string: String with the Abinit input (used in __str__)
        """
        self.string = string
        self.ndtset = int(dvars.pop("ndtset", 1))
        self.udtset = dvars.pop("udtset", None)
        self.jdtset = dvars.pop("jdtset", None)
        if self.udtset is not None or self.jdtset is not None:
            raise ValueError("udtset and jdtset are not supported")

        self.dtsets = [Dataset() for i in range(self.ndtset)]

        # Treat all variables without a dataset index
        kv_list = list(dvars.items())
        for k, v in kv_list:
            if k[-1].isdigit() or any(c in k for c in ("?", ":", "+", "*")): continue
            for d in self.dtsets: d[k] = v
            dvars.pop(k)

        # Treat all variables with a dataset index except those with "?", ":", "+"
        kv_list = list(dvars.items())
        for k, v in kv_list:
            if any(c in k for c in ("?", ":", "+", "*")): continue
            varname, idt = varname_dtindex(k)
            dvars.pop(k)

            if idt > self.ndtset:
                #print("Ignoring key: %s because ndtset: %d" % (k, self.ndtset))
                continue
            self.dtsets[idt-1][varname] = v

        # Now treat series e.g. ecut: 10 ecut+ 5 (? is not treated)
        kv_list = list(dvars.items())
        for k, v in kv_list:
            if "?" in k: continue
            if ":" not in k: continue
            # TODO units
            vname = k[:-1]
            start = str2array(dvars.pop(k))

            incr = dvars.pop(vname + "+", None)
            if incr is not None:
                incr = str2array(incr)
                for dt in self.dtsets:
                    dt[vname] = start.copy()
                    start += incr

            else:
                mult = dvars.pop(vname + "*")
                mult = str2array(mult)
                for dt in self.dtsets:
                    dt[vname] = start.copy()
                    start *= mult

        wrong = []
        for i, dt in enumerate(self.dtsets):
            wlist = [k for k in dt if not is_abivar(k)]
            if wlist:
                wrong.extend(("dataset %d" % i, wlist))
                #print("Wrong variables in dataset %d. wrong_vars = %s" % (i, str(wlist)))
        if wrong:
            raise ValueError("Found wrong variables.\n%s" % str(wrong))

        if dvars:
            raise ValueError("Don't know how handle variables in %s" % str(dvars))

    def __str__(self):
        """String representation."""
        lines = []
        app = lines.append
        header = 10 * "=" + " Input File " + 10 * "="
        app(header)
        app(self.string)
        app(len(header) * "=" + "\n")

        # Print info on structure(s).
        if self.structure is not None:
            app(self.structure.spglib_summary())
        else:
            structures = [dt.structure for dt in self.dtsets]
            app("Input file contains %d structures:" % len(structures))
            for i, struct in enumerate(structures):
                app("Dataset %d" % (i+1))
                app(structure.spglib_summary())
            app(str(self.frame_from_structures()))

        return "\n".join(lines)

    @lazy_property
    def structure(self):
        """
        The structure defined in the input file.
        If the input file contains multiple datasets **AND** the datasets have different structures.
        this property returns None. In this case, one has to access the structure of the individual datasets.
        For example:

            input.dtsets[0].structure

        gives the structure of the first dataset.
        """
        for dt in self.dtsets[1:]:
            if dt.structure != self.dtsets[0].structure:
                logger.info("Datasets have different structures. Returning None. Use input.dtsets[i].structure")
                return None
        return self.dtsets[0].structure

    def frame_from_structures(self, **kwargs):
        """Wraps frame_from_structures."""
        structures = [dt.structure for dt in self.dtsets]
        return frame_from_structures(structures, **kwargs)


def parse_abinit_string(s):
    """
    This function receives a string `s` with the Abinit input and return
    a dictionary {var_name: var_value} when var_value is still a string.
    """
    # Remove comments from lines.
    lines = []
    for line in s.splitlines():
        line.strip()
        i = line.find("#")
        if i != -1: line = line[:i]
        i = line.find("!")
        if i != -1: line = line[:i]
        if line: lines.append(line)

    # 1) Build string of the form "var1 value1 var2 value2"
    # 2) split string in tokens.
    # 3) Evaluate star syntax i.e. "3*2" ==> '2 2 2'
    # 4) Evaluate operators e.g. sqrt(0.75)
    #
    # Step 3 is needed because we are gonna use python to evaluate the operators and
    # in abinit `2*sqrt(0.75)` means `sqrt(0.75) sqrt(0.75)` and not math multiplication!
    #tokens = []
    #for t in " ".join(lines).split():
    #    tokens.extend(expand_star_syntax(t).split())
    #tokens = eval_abinit_operators(tokens)
    # The code above does not work everywhere e.g. typat 2 * 1
    # I think this section should be performed afterwards when we have already constructed the dictionary
    tokens = " ".join(lines).split()

    # Get value of ndtset.
    #try:
    #    i = tokens.index("ndtset")
    #except ValueError:
    #    i = None
    #ndtset = 1 if i is None else int(tokens[i+1])

    varpos = []
    for pos, tok in enumerate(tokens):
        #t = analyze_token(tok)
        #if t.isvar: elif t.isunit:

        if tok[0].isalpha():
            # Either new variable, string defining the unit or operator e.g. sqrt
            if is_abiunit(tok) or is_abioperator(tok): continue

            # new variable
            if tok[-1].isdigit() and "?" not in tok:
                # Handle dataset index.
                l = []
                for i, c in enumerate(tok[::-1]):
                    if c.isalpha(): break
                    l.append(c)
                else:
                    raise ValueError("Cannot find dataset index in %s" % tok)
                l.reverse()

            varpos.append(pos)

    varpos.append(len(tokens) + 1)

    d = {}
    for i, pos in enumerate(varpos[:-1]):
        varname = tokens[pos]
        d[varname] = " ".join(tokens[pos+1: varpos[i+1]])

    err_lines = []
    for k, v in d.items():
        if not v:
            err_lines.append("key %s was not parsed correctly (empty value)" % k)
    if err_lines:
        raise RuntimeError("\n".join(err_lines))

    return d
