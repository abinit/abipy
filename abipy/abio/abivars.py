"""This module contains lookup table with the name of the ABINIT variables."""
import os
import warnings
import numpy as np

from pprint import pformat #, pprint
from monty.string import is_string, boxed
from monty.functools import lazy_property
from monty.termcolor import cprint
from pymatgen.core.units import bohr_to_ang
from abipy.core.structure import Structure, dataframes_from_structures
from abipy.core.mixins import Has_Structure, TextFile, NotebookWriter
from abipy.abio.abivar_database.variables import get_codevars

__all__ = [
    "is_abivar",
    "is_abiunit",
    "AbinitInputFile",
    "AbinitInputParser",
    "structure_from_abistruct_fmt",
]


def is_anaddb_var(varname):
    """True if varname is a valid Anaddb variable."""
    return varname in get_codevars()["anaddb"]


def is_abivar(varname):
    """True if s is an ABINIT variable."""
    # Add include statement
    # FIXME: These variables should be added to the database.
    extra = ["include", "xyzfile"]
    return varname in get_codevars()["abinit"] or varname in extra


# TODO: Move to new directory
ABI_OPERATORS = set(["sqrt", ])

ABI_UNIT_NAMES = {
    s.lower() for s in (
        "au", "nm",
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
    >>> assert expand_star_syntax("*2") == '*2'
    """
    s = s.strip()
    if "*" not in s:
        return s
    else:
        # Handle e.g `pawecutdg*`
        if s[0].isalpha() and s[-1] == "*": return s

    s = s.replace("*", " * ").strip()
    tokens = s.split()
    #tokens = [c.rstrip().lstrip() for c in s.split()]

    # Handle "*2" case i.e. return "*2"
    if len(tokens) == 2 and tokens[0] == "*":
        assert tokens[1] != "*"
        return "".join(tokens)

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
        raise ValueError("Don't know how to handle unit: %s" % str(unit))


def str2array(obj, dtype=float):
    if not is_string(obj): return np.asarray(obj)
    if obj.startswith("*"):
        raise ValueError("This case should be treated by the caller: %s" % str(obj))
    s = expand_star_syntax(obj)
    # Numpy does not understand "0.00d0 0.00d0"
    s = s.lower().replace("d", "e")
    return np.fromstring(s, sep=" ", dtype=dtype)


def eval_abinit_operators(tokens):
    """
    Receive a list of strings, find the occurences of operators supported
    in the input file (e.g. sqrt), evalute the expression and return new list of strings.

    .. note:

        This function is not recursive hence expr like sqrt(1/2) are not supported
    """
    import math # noqa: F401
    import re
    re_sqrt = re.compile(r"[+|-]?sqrt\((.+)\)")

    values = []
    for tok in tokens:
        if tok.startswith("'") or tok.startswith('"'):
            values.append(tok)
        else:
            m = re_sqrt.match(tok)
            if m:
                tok = tok.replace("sqrt", "math.sqrt")
                tok = str(eval(tok))
            if "/" in tok:
                tok = str(eval(tok))
            values.append(tok)

    return values


def abi_tokenize(string):

    string = string.strip()

    # stack = [(char, [start, end]), ...] where char is " or '
    stack = []
    for i, c in enumerate(string):
        if c in ("'", '"'):
            if stack and stack[-1][1][1] == -1:
                # Close quotation mark. Note that the algo is not very robust as nested quotation marks are not supported!
                entry = stack[-1]
                if c != entry[0]:
                    raise ValueError("Found unclosed quotation marks in: %s" % string)
                entry = stack[-1][1][1] = i
            else:
                # First hit.
                stack.append((c, [i, -1]))

    if not stack:
        # No ' or " found in string.
        return string.split()

    # Unroll stack and buid tokens so that e.g. "foo bar" is treated as a single token.
    tokens = []
    base = 0
    for item in stack:
        start, stop = item[1][0], item[1][1]
        # Consistency check:
        if stop == -1:
            raise ValueError("Found unclosed quotation mark in: %s" % string)
        if start > 0:
            tokens.extend(string[base:start-1].strip().split())
            base = stop + 1
        tokens.append(string[start:stop + 1])

    if base < len(string):
        tokens.append(string[base:])

    return tokens


class Dataset(dict, Has_Structure):

    @lazy_property
    def structure(self):
        """
        The initial structure associated to the dataset.
        """

        # First of all check whether the structure is defined through an external file.
        if "structure" in self:
            s = self["structure"].replace('"', "")
            filetype, path = s.split(":")
            #print("filetype:", filetype)
            from abipy import abilab
            if filetype == "poscar":
                return abilab.Structure.from_file(path)

            elif filetype == "abivars":
                with open(path, "rt") as fh:
                    return structure_from_abistruct_fmt(fh.read())

            else:
                try:
                    with abilab.abiopen(path) as abifile:
                        if hasattr(abifile, "final_structure"):
                            return abifile.final_structure
                        else:
                            return abifile.structure
                except Exception as exc:
                    raise RuntimeError("Error while opening file: `%s`:\n%s" % (path, exc))

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
            print("acell:", acell, "znucl:", znucl, "typat:", typat, "kwargs:", kwargs, sep="\n")
            raise exc

    def get_vars(self):
        """
        Return dictionary with variables. The variables describing the crystalline structure
        are removed from the output dictionary.
        """
        geovars = {"acell", "angdeg", "rprim", "ntypat", "natom", "znucl", "typat", "xred", "xcart", "xangst"}
        return {k: self[k] for k in self if k not in geovars}

    def __str__(self):
        return self.to_string()

    def to_string(self, post=None, mode="text", verbose=0):
        """
        String representation.

        Args:
            post: String that will be appended to the name of the variables
            mode: Either `text` or `html` if HTML output with links is wanted.
            verbose: Verbosity level.
        """
        post = post if post is not None else ""
        if mode == "html":
            from abipy.abio.abivars_db import get_abinit_variables
            var_database = get_abinit_variables()

        lines = []
        app = lines.append
        for k in sorted(list(self.keys())):
            vname = k + post
            if mode == "html": vname = var_database[k].html_link(label=vname)
            app("%s %s" % (vname, str(self[k])))

        return "\n".join(lines) if mode == "text" else "\n".join(lines).replace("\n", "<br>")

    def _repr_html_(self):
        """Integration with jupyter_ notebooks."""
        return self.to_string(mode="html")


class AbinitInputFile(TextFile, Has_Structure, NotebookWriter):
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
        super().__init__(filepath)

        with open(filepath, "rt") as fh:
            self.string = fh.read()

        self.datasets = AbinitInputParser().parse(self.string)
        self.ndtset = len(self.datasets)

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        lines = []
        app = lines.append
        header = 10 * "=" + " Input File " + 10 * "="
        app(header)
        app(self.string)
        app(len(header) * "=" + "\n")

        # Print info on structure(s).
        if self.structure is not None:
            app(self.structure.spget_summary())
        else:
            structures = [dt.structure for dt in self.datasets]
            app("Input file contains %d structures:" % len(structures))
            for i, structure in enumerate(structures):
                app(boxed("Dataset: %d" % (i+1)))
                app(structure.spget_summary())
                app("")

            dfs = dataframes_from_structures(structures, index=[i+1 for i in range(self.ndtset)])
            app(boxed("Tabular view (each row corresponds to a dataset structure)"))
            app("")
            app("Lattice parameters:")
            app(str(dfs.lattice))
            app("")
            app("Atomic positions:")
            app(str(dfs.coords))

        return "\n".join(lines)

    @lazy_property
    def has_multi_structures(self):
        """True if input defines multiple structures."""
        return self.structure is None

    def _repr_html_(self):
        """Integration with jupyter notebooks."""
        from abipy.abio.abivars_db import repr_html_from_abinit_string
        return repr_html_from_abinit_string(self.string)
        #return self.to_string(mode="html"))

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
                warnings.warn("Datasets have different structures. Returning None. Use input.datasets[i].structure")
                return None

        return self.datasets[0].structure

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        if not self.has_multi_structures:
            yield self.structure.plot(show=False)
            yield self.structure.plot_bz(show=False)
        else:
            for dt in self.datasets:
                yield dt.structure.plot(show=False)
                yield dt.structure.plot_bz(show=False)

    def write_notebook(self, nbpath=None):
        """
        Write an ipython notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("abinp = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(abinp)"),
        ])

        if self.has_multi_structures:
            nb.cells.extend([
                nbv.new_code_cell("""\
for dataset in inp.datasets:
    print(dataset.structure)"""),
            ])

        if self.ndtset > 1:
            nb.cells.extend([
                nbv.new_code_cell("""\
for dataset in abinp.datasets:
    print(dataset)"""),
            ])

        return self._write_nb_nbpath(nb, nbpath)


class AbinitInputParser(object):

    verbose = 0

    def parse(self, s):
        """
        This function receives a string `s` with the Abinit input and return
        a list of :class:`Dataset` objects.
        """
        # Remove "=" and comments from lines.
        lines = []
        s = s.replace("=", " ")
        for line in s.splitlines():
            line.strip()
            i = line.find("#")
            if i != -1: line = line[:i]
            i = line.find("!")
            if i != -1: line = line[:i]
            if line: lines.append(line)

        # 1) Build string of the form `var1 value1 var2 value2`
        # 2) Split string in tokens.
        # 3) Evaluate star syntax i.e. "3*2" ==> '2 2 2'
        # 4) Evaluate operators e.g. sqrt(0.75)
        #
        # NB: Step 3 is needed because we are gonna use python to evaluate the operators and
        #     in abinit `2*sqrt(0.75)` means `sqrt(0.75) sqrt(0.75)` and not math multiplication!
        #
        #tokens = " ".join(lines).split()
        #print("tokens:\n", pformat(tokens))

        tokens = abi_tokenize(" ".join(lines))  #print("new tokens:\n", pformat(tokens))
        if self.verbose: print("tokens", tokens)

        # Step 3-4
        new_tokens = []
        for t in tokens:
            if t.startswith("'") or t.startswith('"'):
                new_tokens.append(t)
            else:
                l = expand_star_syntax(t).split() # ;print(l)
                new_tokens.extend(l)

        tokens = new_tokens
        if self.verbose: print("new_tokens", new_tokens)

        tokens = self.eval_abinit_operators(tokens) # ; print(tokens)

        varpos = []
        for pos, tok in enumerate(tokens):
            #if not isnewvar(ok): continue
            #print("token:", tok)

            if tok[0].isalpha():
                #print("tok[0]", tok[0])
                # Either new variable or string defining the unit or operator e.g. sqrt
                if is_abiunit(tok) or tok in ABI_OPERATORS or "?" in tok:
                    continue

                # Got new variable
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

                #print("new varname `", tok, "` at positions", pos)
                varpos.append(pos)

        varpos.append(len(tokens))

        # Build dict {varname --> value_string}
        dvars = {}
        for i, pos in enumerate(varpos[:-1]):
            varname = tokens[pos]  # ; print("varname:", varname)
            if pos + 2 == len(tokens):
                dvars[varname] = tokens[-1]
            else:
                dvars[varname] = " ".join(tokens[pos+1: varpos[i+1]])
        # print("dvars:\n", pformat(dvars))

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
        return eval_abinit_operators(tokens)

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


def structure_from_abistruct_fmt(string):
    """
    Parse geometrical information given in the structure:abivars format return Structure object

    A typical input file in "structure:abivars" format looks like::

        # MgB2 lattice structure.
        natom   3
        acell   2*3.086  3.523 Angstrom
        rprim   0.866025403784439  0.5  0.0
               -0.866025403784439  0.5  0.0
                0.0                0.0  1.0

        # Atomic positions
        xred_symbols
            0.0  0.0  0.0 Mg
            1/3  2/3  0.5 B
            2/3  1/3  0.5 B

    i.e. the zcucl, typat and ntypat vars are implicitly defined via the xred_symbols table.
    """
    # Algorithm is as follows:
    #     1) parse the xred_symbols section using a record-based parser.
    #     2) Compute xred, znucl, ntypat and typat from the xred_symbols section
    #     3) Add these values in string format to the head to get a "standard" Abinit input file.
    #     4) Init structure from the "standard" Abinit input file.
    found = 0
    head = []
    xred, symbols = [], []

    for line in string.split("\n"):
        line = line.strip()

        if line.lower() == "xred_symbols":
            found += 1
            continue
        if not found:
            head.append(line)
        if found and line:
            tokens = line.split() #; print(line)
            xred.append([float(t) for t in eval_abinit_operators(tokens[:3])])
            symbols.append(tokens[3])

    if not xred:
        raise ValueError("Cannot find xred_symbols section in string:\n %s" % string)

    xred = np.reshape(xred, (-1, 3))

    # Set small values to zero. This usually happens when the CIF file
    # does not give structure parameters with enough digits.
    xred = np.where(np.abs(xred) > 1e-8, xred, 0.0)

    # Compute ntypat, typat, znucl from xred and symbols.
    natom = len(xred)
    ntypat = 1
    typat = [0]
    for iatom in range(1, natom):
        found = -1
        for jj in range(ntypat):
            if symbols[iatom] == symbols[typat[jj]]:
                found = jj
                break
        if found == -1:
            ntypat += 1
            typat.append(ntypat - 1)
        else:
            typat.append(found)

    from pymatgen.core.periodic_table import Element
    znucl = np.empty(ntypat)
    for iatom in range(natom):
        itypat = typat[iatom]
        znucl[itypat] = Element(symbols[iatom]).Z

    # Convert C indexing to Fortran.
    typat = np.array(typat) + 1

    # Build standard Abinit input string.
    lines = []
    app = lines.append
    app("ntypat %d" % ntypat)
    app("typat " + " ".join("%d" % t for t in typat))
    app("znucl " + " ".join("%f" % z for z in znucl))
    app("xred ")
    for iatom in range(natom):
        v = xred[iatom]
        lines.append("%.16f %.16f %.16f" % (v[0], v[1], v[2]))

    s = "\n".join(head) + "\n" + "\n".join(lines) #; print(s)

    return AbinitInputFile.from_string(s).structure


def validate_input_parser(abitests_dir=None, input_files=None):
    """
    validate/test AbinitInput parser.

    Args:
        dirpath: Abinit tests directory.
        input_files: List of Abinit input files.

    Return: Exit code.
    """
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

    # Files are collected in paths.
    paths = []

    if abitests_dir is not None:
        print("Analyzing directory %s for input files" % abitests_dir)

        for dirpath, dirnames, filenames in os.walk(abitests_dir):
            for fname in filenames:
                path = os.path.join(dirpath, fname)
                if is_abinit_input(path): paths.append(path)

            #import ast
            #init_path = os.path.join(dirpath, "__init__.py")
            #with open(init_path, "rt") as f:
            #    source = f.read()
            #    start = source.find("inp_files = [")
            #    if start == -1:
            #        print("ignoring ", init_path)
            #        continue
            #    stop = source.find("]", start)
            #    if stop == -1:
            #        raise ValueError("Invalid code in %s" % init_path)
            #    print(init_path)
            #    inp_basenames = ast.literal_eval(source[start:stop+1])
            #    print(int_basenames)

    if input_files is not None:
        print("Analyzing files ", str(input_files))
        for arg in input_files:
            if is_abinit_input(arg): paths.append(arg)

    nfiles = len(paths)
    if nfiles == 0:
        cprint("Empty list of input files.", "red")
        return 0

    print("Found %d Abinit input files" % len(paths))
    errpaths = []
    for path in paths:
        print(path + ": ", end="")
        try:
            inp = AbinitInputFile.from_file(path)
            s = str(inp)
            cprint("OK", "green")
        except Exception as exc:
            if not isinstance(exc, NotImplementedError):
                cprint("FAILED", "red")
                errpaths.append(path)
                import traceback
                print(traceback.format_exc())
                #print("[%s]: Exception:\n%s" % (path, str(exc)))
                #with open(path, "rt") as fh:
                #    print(10*"=" + "Input File" + 10*"=")
                #    print(fh.read())
                #    print()
            else:
                cprint("NOTIMPLEMENTED", "magenta")

    if errpaths:
        cprint("failed: %d/%d [%.1f%%]" % (len(errpaths), nfiles, 100 * len(errpaths)/nfiles), "red")
        for i, epath in enumerate(errpaths):
            cprint("[%d] %s" % (i, epath), "red")
    else:
        cprint("All input files successfully parsed!", "green")

    return len(errpaths)
