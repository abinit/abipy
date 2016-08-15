"""This module contains lookup table with the name of the ABINIT variables."""
from __future__ import division, print_function, unicode_literals, absolute_import

import json
import numpy as np

from monty.string import is_string, boxed
from monty.functools import lazy_property
from pymatgen.core.units import bohr_to_ang
from abipy.core.structure import Structure, frame_from_structures
from abipy.core.mixins import Has_Structure

import logging
logger = logging.getLogger(__name__)

__all__ = [
    "is_abivar",
    "is_abiunit",
    "AbinitInputFile",
    "AbinitInputParser",
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
            # Add include statement
            ABI_VARNAMES.add("include")

    return s in ABI_VARNAMES


ABI_OPERATORS = set(["sqrt", ])
#ABI_MULTI_TAGS = set(["+", ":", "*", "?"])

ABI_UNIT_NAMES = {
    s.lower() for s in (
        "au",
        "Angstr", "Angstrom", "Angstroms", "Bohr", "Bohrs",
        "eV", "Ha", "Hartree", "Hartrees", "K", "Ry", "Rydberg", "Rydbergs",
        "T", "Tesla",)
}

def is_abiunit(s):
    """
    True if string is one of the units supported by the ABINIT parser
    """
    if not is_string(s): return False
    return s.lower() in ABI_UNIT_NAMES


def expand_star_syntax(s):
    """
    Evaluate star syntax. Return new string
    Remember that Abinit does not accept white spaces.
    For example `typat 2 * 1` is not valid.

    >>> assert expand_star_syntax("3*2") == '2 2 2'
    >>> assert expand_star_syntax("2 *1") == '1 1'
    >>> assert expand_star_syntax("1 2*2") == '1 2 2'
    >>> assert expand_star_syntax("*2") == '* 2'
    """
    s = s.strip()
    if "*" not in s:
        return s
    else:
        # Handle e.g `pawecutdg*`
        if s[0].isalpha() and s[-1] == "*": return s

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
    if obj.startswith("*"):
        raise ValueError("This case should be treated by the caller: %s" % str(obj))
    return np.fromstring(expand_star_syntax(obj), sep=" ")
    #return np.fromstring(obj, sep=" ")


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

        if "xred" in self:
            # v3/Input/t05.in
            #kwargs["xred"] = np.fromstring(self["xred"], sep=" ")
            #print(self["xred"])
            kwargs["xred"] = str2array(self["xred"])
            #print(kwargs["xred"])
        elif "xcart" in self:
            kwargs["xcart"] = str2array_bohr(self["xcart"])
        elif "xangst" in self:
            kwargs["xangst"] = np.fromstring(self["xangst"], sep=" ")
        else:
            raise ValueError("xred|xcart|xangst must be given in input")

        # Get important dimensions.
        ntypat = int(self.get("ntypat", 1))
        natom = int(self.get("natom", 1))

        # znucl(npsp)
        znucl = self["znucl"]
        if znucl.startswith("*"):
            i = znucl.find("*")
            znucl_size = natom if "npsp" not in self else int(self["npsp"])
            znucl = zcnul_size * [float(znucl[i+1:])]
        else:
            znucl = str2array(self["znucl"])

        # v67mbpt/Input/t12.in
        typat = self["typat"]
        if typat.startswith("*"):
            i = typat.find("*")
            typat = ntypat * [int(typat[i+1:])]
        else:
            typat = str2array(self["typat"])

        #print(kwargs["rprim"])
        return Structure.from_abivars(
            acell=acell,
            znucl=znucl,
            typat=typat,
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
        return cls(string)

    def __init__(self, string):
        """
        Args:
            dvars: Dictionary {varname: value}.
            string: String with the Abinit input (used in __str__)
        """
        self.string = string
        self.dtsets = AbinitInputParser().parse(string)
        self.ndtset = len(self.dtsets)

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
            for i, structure in enumerate(structures):
                app(boxed("Dataset: %d" % (i+1)))
                app(structure.spglib_summary())
                app("")

            app(boxed("Tabular view (each row corresponds to a dataset structure)"))
            app(str(frame_from_structures(structures, index=[i+1 for i in range(self.ndtset)])))

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


class AbinitInputParser(object):
    verbose = 0

    def parse(self, s):
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
        tokens = " ".join(lines).split()
        # Step 3 is needed because we are gonna use python to evaluate the operators and
        # in abinit `2*sqrt(0.75)` means `sqrt(0.75) sqrt(0.75)` and not math multiplication!
        # /Users/gmatteo/git/abinit/tests/v7/Input/t03.in
        print("tokens", tokens)
        new_tokens = []
        for t in tokens:
            l = expand_star_syntax(t).split()
            print("t", t, "l", l)
            new_tokens.extend(l)
        tokens = new_tokens
        print("new_tokens", new_tokens)

        tokens = self.eval_abinit_operators(tokens)
        # The code above does not work everywhere e.g. typat 2 * 1
        # I think this section should be performed afterwards when we have already constructed the dictionary
        #tokens = " ".join(lines).split()

        varpos = []
        for pos, tok in enumerate(tokens):

            if tok[0].isalpha():
                # Either new variable, string defining the unit or operator e.g. sqrt
                if is_abiunit(tok) or tok in ABI_OPERATORS: continue

                # Have new variable
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

        dvars = {}
        for i, pos in enumerate(varpos[:-1]):
            varname = tokens[pos]
            dvars[varname] = " ".join(tokens[pos+1: varpos[i+1]])

        err_lines = []
        for k, v in dvars.items():
            if not v:
                err_lines.append("key %s was not parsed correctly (empty value)" % k)
        if err_lines:
            raise RuntimeError("\n".join(err_lines))

        # Get value of ndtset.
        ndtset = int(dvars.pop("ndtset", 1))
        udtset = dvars.pop("udtset", None)
        jdtset = dvars.pop("jdtset", None)
        if udtset is not None or jdtset is not None:
            raise NotImplementedError("udtset and jdtset are not supported")

        if "spgroup" in dvars or "nobj" in dvars:
            raise NotImplementedError(
                "Abinit spgroup builder is not supported. Structure must be given explicitly!")

        dtsets = [Dataset() for i in range(ndtset)]

        # Treat all variables without a dataset index
        kv_list = list(dvars.items())
        for k, v in kv_list:
            if k[-1].isdigit() or any(c in k for c in ("?", ":", "+", "*")): continue
            for d in dtsets: d[k] = v
            dvars.pop(k)

        # Treat all variables with a dataset index except those with "?", ":", "+"
        kv_list = list(dvars.items())
        for k, v in kv_list:
            if any(c in k for c in ("?", ":", "+", "*")): continue
            varname, idt = self.varname_dtindex(k)
            dvars.pop(k)

            if idt > ndtset:
                if self.verbose: print("Ignoring key: %s because ndtset: %d" % (k, ndtset))
                continue
            dtsets[idt-1][varname] = v

        # Now treat series e.g. ecut: 10 ecut+ 5 (? is not treated here)
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
                for dt in dtsets:
                    dt[vname] = start.copy()
                    start += incr

            else:
                mult = dvars.pop(vname + "*")
                mult = str2array(mult)
                for dt in dtsets:
                    dt[vname] = start.copy()
                    start *= mult

        wrong = []
        for i, dt in enumerate(dtsets):
            wlist = [k for k in dt if not is_abivar(k)]
            if wlist:
                wrong.extend(("dataset %d" % i, wlist))
        if wrong:
            raise ValueError("Found wrong variables.\n%s" % str(wrong))

        if dvars:
            raise ValueError("Don't know how handle variables in %s" % str(dvars))

        return dtsets

    @staticmethod
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
                tok = eval(tok)
            if "/" in tok:
                tok = eval(tok)
            values.append(str(tok))
        return values

    @staticmethod
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
