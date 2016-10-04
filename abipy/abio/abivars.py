"""This module contains lookup table with the name of the ABINIT variables."""
from __future__ import division, print_function, unicode_literals, absolute_import

import json
import numpy as np

from pprint import pformat
from monty.string import is_string, boxed
from monty.functools import lazy_property
from pymatgen.core.units import bohr_to_ang
from abipy.core.structure import Structure, frames_from_structures
from abipy.core.mixins import Has_Structure, TextFile

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
            # FIXME: This should be added to the database.
            ABI_VARNAMES.add("include")
            ABI_VARNAMES.add("xyzfile")

    return s in ABI_VARNAMES


ABI_OPERATORS = set(["sqrt", ])

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
    if not is_string(obj):
        return np.asarray(obj)

    # Treat e.g. acell 3 * 1.0
    obj = expand_star_syntax(obj)
    # Numpy does not understand "0.00d0 0.00d0"
    obj = obj.lower().replace("d", "e")

    tokens = obj.split()
    if not tokens[-1].isalpha():
        # No unit
        return np.fromstring(obj, sep=" ")

    unit = tokens[-1]
    if unit in ("angstr", "angstrom", "angstroms"):
        return np.fromstring(" ".join(tokens[:-1]), sep=" ") / bohr_to_ang
    elif unit in ("bohr", "bohrs", "au"):
        return np.fromstring(" ".join(tokens[:-1]), sep=" ")
    else:
        raise ValueError("Don't know how to handle unit: %s" % unit)


def str2array(obj, dtype=float):
    if not is_string(obj): return np.asarray(obj)
    if obj.startswith("*"):
        raise ValueError("This case should be treated by the caller: %s" % str(obj))
    s = expand_star_syntax(obj)
    # Numpy does not understand "0.00d0 0.00d0"
    s = s.lower().replace("d", "e")
    return np.fromstring(s, sep=" ", dtype=dtype)


class Dataset(dict, Has_Structure):

    @lazy_property
    def structure(self):
        # Get lattice.
        kwargs = {}
        if "angdeg" in self:
            if "rprim" in self:
                raise ValueError("rprim and angdeg cannot be used together!")
            angdeg = str2array(self["angdeg"])
            angdeg.shape = (3)
            kwargs["angdeg"] = angdeg
        else:
            # Handle structure specified with rprim.
            kwargs["rprim"] = str2array_bohr(self.get("rprim", "1.0 0 0 0 1 0 0 0 1"))

        # Default value for acell.
        acell = str2array_bohr(self.get("acell", "1.0 1.0 1.0"))

        # Get important dimensions.
        ntypat = int(self.get("ntypat", 1))
        natom = int(self.get("natom", 1))

        # znucl(npsp)
        znucl = self["znucl"]
        if znucl.startswith("*"):
            i = znucl.find("*")
            znucl_size = natom if "npsp" not in self else int(self["npsp"])
            znucl = znucl_size * [float(znucl[i+1:])]
        else:
            znucl = str2array(self["znucl"])

        # v67mbpt/Input/t12.in
        typat = self["typat"]
        if typat.startswith("*"):
            i = typat.find("*")
            typat = np.array(natom * [int(typat[i+1:])], dtype=int)
        else:
            typat = str2array(self["typat"], dtype=int)

        # Extract atomic positions.
        # Select first natom entries (needed if multidatasets with different natom)
        #    # v3/Input/t05.in
        typat = typat[:natom]
        for k in ("xred", "xcart", "xangst"):
            toarray = str2array_bohr if k == "xcart" else str2array
            if k in self:
                arr = np.reshape(toarray(self[k]), (-1, 3))
                kwargs[k] = arr[:natom]
                break
        else:
            raise ValueError("xred|xcart|xangst must be given in input")

        try:
            return Structure.from_abivars(acell=acell, znucl=znucl, typat=typat, **kwargs)
        except Exception as exc:
            print("Wrong inputs passed to Structure.from_abivars:")
            print("  acell", acell)
            print("  znucl", znucl)
            print("  typat", typat)
            print("  kwargs", kwargs)
            raise exc


class AbinitInputFile(TextFile, Has_Structure):
    """
    This object parses the Abinit input file, stores the variables in
    dict-like objects (Datasets) and build `Structure` objects from
    the input variables. Mainly used for inspecting the structure
    declared in the Abinit input file.
    """
    @classmethod
    def from_string(cls, string):
        """Build the object from string."""
        import tempfile
        _, filename = tempfile.mkstemp(suffix=".abi", text=True)
        with open(filename, "wt") as fh:
            fh.write(string)
        return cls(filename)

    def __init__(self, filepath):
        super(AbinitInputFile, self).__init__(filepath)

        with open(filepath, "rt") as fh:
            string = fh.read()

        self.string = string
        self.datasets = AbinitInputParser().parse(string)
        self.ndtset = len(self.datasets)

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
            structures = [dt.structure for dt in self.datasets]
            app("Input file contains %d structures:" % len(structures))
            for i, structure in enumerate(structures):
                app(boxed("Dataset: %d" % (i+1)))
                app(structure.spglib_summary())
                app("")

            dfs = frames_from_structures(structures, index=[i+1 for i in range(self.ndtset)])
            app(boxed("Tabular view (each row corresponds to a dataset structure)"))
            app("")
            app("Lattice parameters:")
            app(str(dfs.lattice))
            app("")
            app("Atomic positions:")
            app(str(dfs.coords))

        return "\n".join(lines)

    def close(self):
        """NOP, required by ABC."""

    @lazy_property
    def structure(self):
        """
        The structure defined in the input file.

        If the input file contains multiple datasets **AND** the datasets
        have different structures, this property returns None.
        In this case, one has to access the structure of the individual datasets.
        For example:

            input.datasets[0].structure

        gives the structure of the first dataset.
        """
        for dt in self.datasets[1:]:
            if dt.structure != self.datasets[0].structure:
                logger.info("Datasets have different structures. Returning None. Use input.datasets[i].structure")
                return None
        return self.datasets[0].structure


class AbinitInputParser(object):
    verbose = 0

    def parse(self, s):
        """
        This function receives a string `s` with the Abinit input and return
        a list of :class:`Dataset` objects.
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

        tokens = " ".join(lines).split()
        # Step 3 is needed because we are gonna use python to evaluate the operators and
        # in abinit `2*sqrt(0.75)` means `sqrt(0.75) sqrt(0.75)` and not math multiplication!
        if self.verbose: print("tokens", tokens)
        new_tokens = []
        for t in tokens:
            l = expand_star_syntax(t).split()
            #print("t", t, "l", l)
            new_tokens.extend(l)
        tokens = new_tokens
        if self.verbose: print("new_tokens", new_tokens)

        tokens = self.eval_abinit_operators(tokens)
        #print(tokens)

        varpos = []
        for pos, tok in enumerate(tokens):
            #if not isnewvar(ok): continue

            if tok[0].isalpha():
                # Either new variable, string defining the unit or operator e.g. sqrt
                if is_abiunit(tok) or tok in ABI_OPERATORS or "?" in tok:
                    continue

                # Have new variable
                if tok[-1].isdigit(): # and "?" not in tok:
                    # Handle dataset index.
                    l = []
                    for i, c in enumerate(tok[::-1]):
                        if c.isalpha(): break
                        l.append(c)
                    else:
                        raise ValueError("Cannot find dataset index in token: %s" % tok)
                    l.reverse()
                    #if not is_abivar(tok):
                        #continue
                        #raise ValueError("Expecting variable but got: %s" % tok)

                #print("new var", tok, pos)
                varpos.append(pos)

        varpos.append(len(tokens))

	# Build dict {varname --> value_string}
        dvars = {}
        for i, pos in enumerate(varpos[:-1]):
            varname = tokens[pos]
            if pos + 2 == len(tokens):
                dvars[varname] = tokens[-1]
            else:
                dvars[varname] = " ".join(tokens[pos+1: varpos[i+1]])

        #print(dvars)
        err_lines = []
        for k, v in dvars.items():
            if not v:
                err_lines.append("key `%s` was not parsed correctly (empty value)" % k)
        if err_lines:
            raise RuntimeError("\n".join(err_lines))

        # Get value of ndtset.
        ndtset = int(dvars.pop("ndtset", 1))
        udtset = dvars.pop("udtset", None)
        jdtset = dvars.pop("jdtset", None)
        if udtset is not None:
            raise NotImplementedError("udtset is not supported")

	# Build list of datasets.
        datasets = [Dataset() for i in range(ndtset)]

        # Treat all variables without a dataset index
        kv_list = list(dvars.items())
        for k, v in kv_list:
            if k[-1].isdigit() or any(c in k for c in ("?", ":", "+", "*")): continue
            for d in datasets: d[k] = v
            dvars.pop(k)

        # Treat all variables with a dataset index except those with "?", ":", "+"
        kv_list = list(dvars.items())
        for k, v in kv_list:
            if any(c in k for c in ("?", ":", "+", "*")): continue
            varname, idt = self.varname_dtindex(k)
            dvars.pop(k)
            #if varname == "angdeg": raise ValueError("got angdeg")
            if idt > ndtset:
                if self.verbose: print("Ignoring key: %s because ndtset: %d" % (k, ndtset))
                continue
            datasets[idt-1][varname] = v

        # Now treat series e.g. ecut: 10 ecut+ 5 (NB: ? is not treated here)
        kv_list = list(dvars.items())
        for k, v in kv_list:
            if "?" in k: continue
            if ":" not in k: continue
            # TODO units
            vname = k[:-1]
            start = str2array(dvars.pop(k))

	    # Handle ecut+ or ecut*
            incr = dvars.pop(vname + "+", None)
            if incr is not None:
                incr = str2array(incr)
                for dt in datasets:
                    dt[vname] = start.copy()
                    start += incr

            else:
                mult = dvars.pop(vname + "*")
                mult = str2array(mult)
                for dt in datasets:
                    dt[vname] = start.copy()
                    start *= mult

	# Consistency check
	# 1) dvars should be empty
        if dvars:
            raise ValueError("Don't know how handle variables in:\n%s" % pformat(dvars), indent=4)

	# 2) Keys in datasets should be valid Abinit input variables.
        wrong = []
        for i, dt in enumerate(datasets):
            wlist = [k for k in dt if not is_abivar(k)]
            if wlist:
                wrong.extend(("dataset %d" % i, wlist))
        if wrong:
            raise ValueError("Found variables that are not registered in the abipy database:\n%s" % pformat(wrong, indent=4))

	# 3) We don't support spg builder: dataset.structure will fail or, even worse,
        #    spglib will segfault so it's better to raise here!
        for dt in datasets:
            if "spgroup" in dt or "nobj" in dt:
                raise NotImplementedError(
                    "Abinit spgroup builder is not supported. Structure must be given explicitly!")

        if jdtset is not None:
            # Return the datasets selected by jdtset.
            datasets = [datasets[i-1] for i in np.fromstring(jdtset, sep=" ", dtype=int)]

        return datasets

    @staticmethod
    def eval_abinit_operators(tokens):
        """
        Receive a list of strings, find the occurences of operators supported
        in the input file (e.g. sqrt), evalute the expression and return new list of strings.

	.. note:

	    This function is not recursive hence expr like sqrt(1/2) are not supported
        """
        import math
        import re
        re_sqrt = re.compile("[+|-]?sqrt\((.+)\)")

        values = []
        for tok in tokens:
            m = re_sqrt.match(tok)
            if m:
                tok = tok.replace("sqrt", "math.sqrt")
                tok = str(eval(tok))
            if "/" in tok: # Note true_division from __future__
                tok = str(eval(tok))
            values.append(tok)
        return values

    @staticmethod
    def varname_dtindex(tok):
        """
        >>> p = AbinitInputParser()
        >>> assert p.varname_dtindex("acell1") == ("acell", 1)
        >>> assert p.varname_dtindex("fa1k2") == ("fa1k", 2)
        """
        l = []
        for i, c in enumerate(tok[::-1]):
            if c.isalpha(): break
            l.append(c)
        else:
            raise ValueError("Cannot find dataset index in: %s" % tok)

        assert i > 0
        l.reverse()
        dtidx = int("".join(l))
        varname = tok[:len(tok)-i]

        return varname, dtidx


def validate_input_parser():
    """
    Script to validate/test AbinitInput parser.

    Usage:
	$ test_input_parser DIRECTORY
	$ test_input_parser FILES

    Return: Exit code.
    """
    import sys
    import os
    from abipy import abilab

    try:
        args = sys.argv[1:]
    except:
        print(validate_input_parser.__doc__)
        return 1

    if "-h" in args or "--help" in args:
        print(validate_input_parser.__doc__)
        return 0

    def is_abinit_input(path):
        """
        True if path is one of the input files used in the Abinit Test suite.
        """
        if path.endswith(".abi"): return True
        if not path.endswith(".in"): return False

        with open(path, "rt") as fh:
            for line in fh:
                if "executable" in line and "abinit" in line: return True
            return False

    # Collect files
    paths = []

    if len(args) == 1 and os.path.isdir(args[0]):
        print("Analyzing directory %s for input files" % args[0])
        for dirpath, dirnames, filenames in os.walk(args[0]):
            for fname in filenames:
                path = os.path.join(dirpath, fname)
                if is_abinit_input(path): paths.append(path)
    else:
        for arg in args:
            if is_abinit_input(arg):
                paths.append(arg)

    nfiles = len(paths)
    if nfiles == 0:
        print("Empty list of input files.")
        return 0

    print("Found %d Abinit input files" % len(paths))
    errpaths = []
    for path in paths:
        print("About to analyze:", path)
        try:
            inp = abilab.AbinitInputFile.from_file(path)
            s = str(inp)
            print(s)
        except Exception as exc:
            if not isinstance(exc, NotImplementedError):
                errpaths.append(path)
                import traceback
                print(traceback.format_exc())
                #print("[%s]: Exception:\n%s" % (path, str(exc)))

                #with open(path, "rt") as fh:
                #    print(10*"=" + "Input File" + 10*"=")
                #    print(fh.read())
                #    print()

    print("failed: %d/%d [%.1f%%]" % (len(errpaths), nfiles, 100 * len(errpaths)/nfiles))
    if errpaths:
        for i, epath in enumerate(errpaths):
            print("[%d] %s" % (i, epath))

    return len(errpaths)
