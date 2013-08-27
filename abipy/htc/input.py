from __future__ import print_function, division

import os
import collections
import warnings
import itertools
import numpy as np

from pymatgen.io.abinitio.pseudos import PseudoTable
from abipy.core import Structure

__all__ = ['AbiInput']

# variables that must have a unique value throughout all the datasets.
_NO_MULTI = [
    "jdtset",
    "ndtset",
    "ntypat",
    "znucl",
    "cpus",
    "cpum",
    "cpuh",
    "npsp",
    "timopt",
    "iatcon",
    "natcon",
    "nconeq",
    "prtxtypat",
    "wtatcon",
]


class AbinitVariable(object):

    def __init__(self, name, value):
        self.name = name
        self.value = value

    def __repr__(self):
        return "<%s, name=%s, value=%s>" % (self.__class__.__name__, self.name, self.value)

    def __str__(self):
        return format_var(self.name, self.value)


class AbinitInputError(Exception):
    """Base error class for exceptions raised by `AbiInput`"""


class AbiInput(object):

    Error = AbinitInputError

    def __init__(self, pseudos, pseudo_dir="", ndtset=1, comment=""):
        """
        Args:
            pseudos:
                String or list of string with the name of the pseudopotential files.
            pseudo_dir:
                Name of the directory where pseudo files are located.
            ndtset:
                Number of datasets.
            comment:
                Optional string with a comment that will be placed at the beginning of the file.
        """
        # Dataset[0] contains the global variables common to the different datasets
        # Dataset[1:ndtset+1] stores the variables specific to the different datasets.
        self._ndtset = ndtset

        self._datasets = []
        for i in range(ndtset+1):
            dt0 = None
            if i > 0: dt0 = self._datasets[0]
            self._datasets.append(Dataset(index=i, dt0=dt0))

        self._datasets[0]["ndtset"] = ndtset

        # Setup of the pseudopotential files.
        self._pseudo_dir = os.path.abspath(pseudo_dir)
        if isinstance(pseudos, str): pseudos = [pseudos]
        pseudo_paths = [os.path.join(self._pseudo_dir, p) for p in pseudos]

        missing = [p for p in pseudo_paths if not os.path.exists(p)]
        if missing:
            raise self.Error("Cannot find the following pseudopotential files:\n%s" % str(missing)) 

        self._pseudos = PseudoTable(pseudo_paths)

        if comment:
            self.set_comment(comment)

    @property
    def ndtset(self):
        """Number of datasets."""
        return self[0]["ndtset"]

    @property
    def global_vars(self):
        """Dictionary with global variables (common to the different datasets)."""
        return self[0]

    def __str__(self):
        s = ""
        for (i, dataset) in enumerate(self):
            header = "*** DATASET %d ***" % i
            if i == 0: 
                header = "*** GLOBAL VARIABLES ***"

            str_dataset = str(dataset)
            if str_dataset:
                header = len(header) * "*" + "\n" + header + "\n" + len(header) * "*" + "\n"
                s += "\n" + header + str(dataset) + "\n"

        return s

    def __getitem__(self, key):
        return self._datasets[key]

    def __setattr__(self, varname, value):

        if varname.startswith("ndtset"):
            raise self.Error("Cannot change read-only variable ndtset")

        if varname.startswith("_"):
            # Standard behaviour for "private" variables.
            super(AbiInput, self).__setattr__(varname, value)

        else:
            if not varname[-1].isdigit():
                # Global variable.
                self[0].set_variable(varname, value)

            else:
                # Find the index of the dataset.
                idt = ""
                for i, c in enumerate(reversed(varname)):
                    if c.isdigit():
                        idt += c
                    else:
                        break

                idt = int("".join(reversed(idt)))
                if not (self.ndtset >= idt  >= 1):
                    raise self.Error("Invalid dataset index %d, ndtset %d " % (idt, self.ndtset))

                # Strip the numeric index.
                varname = varname[:len(varname)-i]

                self[idt].set_variable(varname, value)

                # Handle no_multi variables.
                if varname in _NO_MULTI:

                    if varname in self[0]:
                        glob_value = np.array(self[0][varname])
                        isok  = np.all(glob_value == np.array(value))
                        if not isok:
                            err_msg = "NO_MULTI variable: dataset 0: %s, dataset %d: %s" % (str(glob_value), idt, str(value))
                            raise self.Error(err_msg)
                    else:
                        self[0].set_variable(varname, value)

    def __getattr__(self, varname):

        if not varname[-1].isdigit() and varname in self.global_vars:
            # Global variable.
            value = self[0][varname]
            return AbinitVariable(varname, value)
                                                                        
        elif varname[-1].isdigit():
            # Find the index of the dataset.
            idt = ""
            for i, c in enumerate(reversed(varname)):
                if c.isdigit():
                    idt += c
                else:
                    break

            idt = int("".join(reversed(idt)))
            if not (self.ndtset >= idt  >= 1):
                raise self.Error("Invalid dataset index %d, ndtset %d " % (idt, self.ndtset))

            # Strip the numeric index.
            var = varname[:len(varname)-i]

            try:
                value = self[idt][var]
                return AbinitVariable(var, value)

            except:
                raise AttributeError("cannot find attribute %s" % varname)

        else:
            raise AttributeError("cannot find attribute %s" % varname)

    @property
    def pseudo_dir(self):
        """The directory where pseudos are located."""
        return self._pseudo_dir

    @property
    def pseudos(self):
        """List of pseudopotentials."""
        return self._pseudos

    @staticmethod
    def _dtset2range(dtset):
        """
        Helper function that coverts dtset into a range. Accepts int, slice objects or iterable.
        """
        if isinstance(dtset, int):
            return [dtset]

        elif isinstance(dtset, slice): 
            start = dtset.start if dtset.start is not None else 0
            stop = dtset.stop
            step = dtset.step if dtset.step is not None else 1
            return range(start, stop, step)

        elif isinstance(dtset, collections.Iterable):
            return dtset

        raise self.Error("Don't know how to convert %s to a range-like object" % str(dtset))

    def set_variables(self, dtset=0, **vars):
        """
        Set the value of a set of input variables.

        Args:
            dtset:
                Int with the index of the dataset, slice object of iterable 
            vars:
                Dictionary with the variables.
        """
        for idt in self._dtset2range(dtset):
            self[idt].set_variables(**vars)

    def list_variable(self, varname):
        # Global value
        glob = self[0].get(varname, None)
        # Values defined in the datasets.
        vals = [self[idt].get(varname, None) for idt in range(1,self.ndtset+1)]
        
        # Replace None with glob.
        values = []
        for v in vals:
            if v is not None:
                values.append(v)
            else:
                values.append(glob)

        return values

    def set_comment(self, comment, dtset=0):
        """Set the main comment in the input file."""
        for idt in self._dtset2range(dtset):
            self[idt].set_comment(comment)

    def linspace(self, varname, start, stop, endpoint=True):
        """
        Returns `ndtset` evenly spaced samples, calculated over the interval [`start`, `stop` ].

        The endpoint of the interval can optionally be excluded.

        Args:
            start: 
                The starting value of the sequence.
            stop:
                The end value of the sequence, unless `endpoint` is set to False.
                In that case, the sequence consists of all but the last of ``ndtset + 1``
                evenly spaced samples, so that `stop` is excluded.  Note that the step
                size changes when `endpoint` is False.
        """
        if varname in self[0]:
            raise self.Error("varname %s is already defined in globals" % varname)
            
        lspace = np.linspace(start, stop, num=self.ndtset, endpoint=endpoint, retstep=False)

        for (dataset, value) in zip(self[1:], lspace):
            dataset.set_variable(varname, value)

    def arange(self, varname, start, stop, step):
        """
        Return evenly spaced values within a given interval.

        Values are generated within the half-open interval ``[start, stop)``
        (in other words, the interval including `start` but excluding `stop`).

        When using a non-integer step, such as 0.1, the results will often not
        be consistent.  It is better to use ``linspace`` for these cases.

        Args:
            start: 
                Start of interval.  The interval includes this value. The default start value is 0.
            stop:
                End of interval.  The interval does not include this value, except
                in some cases where `step` is not an integer and floating point
            step: 
                Spacing between values.  For any output `out`, this is the distance
                between two adjacent values, ``out[i+1] - out[i]``.  The default
                step size is 1.  If `step` is specified, `start` must also be given.
        """
        if varname in self[0]:
            raise self.Error("varname %s is already defined in globals" % varname)

        arang = np.arange(start=start, stop=stop, step=step)
        if len(arang) != self.ndtset:
            raise self.Error("len(arang) %d != ndtset %d" % (len(arang), self.ndtset))

        for (dataset, value) in zip(self[1:], arang):
            dataset.set_variable(varname, value)

    def product(self, *items):
        """
        Cartesian product of input iterables.  Equivalent to nested for-loops.

        .. example:

            inp.product("tsmear", "ngkpt", [[2,2,2], [4,4,4]], [0.1, 0.2, 0.3])
        """
        # Split items into varnames and values
        for i, item in enumerate(items):
            if not isinstance(item, str):
                break

        varnames, values = items[:i], items[i:]
        if len(varnames) != len(values):
            raise self.Error("The number of variables must equal the number of lists")

        varnames = [ [varnames[i]] * len(values[i]) for i in range(len(values))]
        varnames = itertools.product(*varnames)
        values = itertools.product(*values)
        #print([v for v in varnames])
        #print([v for v in values])

        idt = 0
        for names, values in zip(varnames, values):
            idt += 1
            self[idt].set_variables(**{k:v for k, v in zip(names, values)})
        
        if idt != self.ndtset:
            raise self.Error("The number of configurations must equal ndtset while %d != %d" % (idt, self.ndtset))

    def set_structure(self, structure, dtset=0):
        """Set the `Structure` object for the specified dtset."""
        for idt in self._dtset2range(dtset):
            self[idt].set_structure(structure)

    def set_structure_from_file(self, filepath, dtset=0):
        """Set the `Structure` object for the specified dtset (data is read from filepath)."""
        structure = Structure.from_file(filepath)
        self.set_structure(structure, dtset=dtset)
        return structure

    def set_kmesh(self, ngkpt, shiftk, kptopt=1, dtset=0):
        """
        Set the variables defining the k-point sampling for the specified dtset.
        """
        for idt in self._dtset2range(dtset):
            self[idt].set_kmesh(ngkpt, shiftk, kptopt=kptopt)

    def set_kpath(self, ndivsm, kptbounds=None, dtset=0):
        """
        Set the variables defining the k-path for the specified dtset.

        The list of K-points is taken from the pymatgen database if kptbounds is None

        Args:
            ndivsm:
                Number of divisions for the smallest segment.
            kptbounds:
                k-points defining the path in k-space.
            dtset:
                Index of the dataset. 0 for global variables.
        """
        for idt in self._dtset2range(dtset):
            self[idt].set_kpath(ndivsm, kptbounds=kptbounds)

    def set_kptgw(self, kptgw, bdgw, dtset=0):
        """
        Set the variables defining the k point for the GW corrections
        """
        for idt in self._dtset2range(dtset):
            self[idt].set_kptgw(kptgw, bdgw)

    #def make_filesfiles(self):
    #    lines = []
    #    app = lines.append

    #    app(".abi")
    #    app(".abo")
    #    app("i")
    #    app("o")
    #    app("t")
    #    for pseudo in self.pseudos:
    #        app(pseudo.filepath)
    #    return "\n".join(lines)

    def write(self, filepath):
        """
        Write the input file to file filepath.
        """
        dirname = os.path.dirname(filepath)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
                                                                                      
        # Write input file.
        with open(filepath, "w") as fh:
           fh.write(str(self))

        # TODO filesfile.
        #if filesfile is not None
        #    with open(filesfiles.filepath, "w") as fh:
        #        fh.write(self.make_filesfiles(filesfile))

_UNITS = {
        'bohr' : 1.0,
        'angstrom' : 1.8897261328856432,
        'hartree' : 1.0,
        'Ha' : 1.0,
        'eV' : 0.03674932539796232,
        }


class Dataset(collections.Mapping):
    """
    This object stores the variables for a single dataset.
    """
    Error = AbinitInputError

    def __init__(self, index, dt0, *args, **kwargs):
        self.vars = collections.OrderedDict(*args, **kwargs)
        self._index = index
        self._dt0 = dt0

    def __repr__(self):
        "<%s at %s>" % self.__class__.__name__, id(self)

    def __str__(self):
        return self.to_string()

    def __len__(self):
        return self.vars

    def __iter__(self):
        return self.vars.__iter__()

    def __getitem__(self, key):
        return self.vars[key]

    def __setitem__(self, key, value):
        self.vars[key] = value

    def __contains__(self, key):
        return key in self.vars

    def keys(self):
        return self.vars.keys()

    def items(self):
        return self.vars.items()

    def get(self, value, default=None):
        try:
            return self[value]
        except KeyError:
            return default

    @property
    def index(self):
        """The index of the dataset within the input."""
        return self._index

    @property
    def dt0(self):
        return self._dt0

    @property
    def global_vars(self):
        """Copy of the dictionary with the global variables."""
        return self.dt0.vars.copy()

    @property
    def allvars(self):
        """
        Dictionay with the variables of the dataset + the global variables.
        Variables local to the dataset overwrite the global values (if any).
        """
        all_vars = self.global_vars
        all_vars.update(self.vars)
        return all_vars

    def to_string(self, sortmode=None):
        """String representation."""
        lines = []
        app = lines.append

        if self.comment:
            app("# " + self.comment.replace("\n", "\n#"))

        post = "" if self.index == 0 else str(self.index)

        if sortmode is None:
            # no sorting.
            keys = self.keys()
        elif sortmode is "a":
            # alphabetical order.
            keys = sorted(self.keys())
        #elif sortmode is "t": TODO
        # group variables by topic
        else:
            raise ValueError("Unsupported value for sortmode %s" % str(sortmode))

        for var in keys:
            if self.index != 0 and var in _NO_MULTI:
                continue
            varname = var + post
            value = self[var]
            app("%s %s" % (varname, format_var(var, value)))

        return "\n".join(lines)

    @property
    def comment(self):
        try:
            return self._comment

        except AttributeError:
            return None

    def set_comment(self, comment):
        """Set a comment to be included at the top of the file."""
        self._comment = comment

    def set_variable(self, varname, value):
        """Set a single variable."""
        if varname in self:
            warnings.warn("Variable %s is already defined with value:\n%s" % (varname, str(self[varname])))

        self[varname] = value

        # Handle no_multi variables.
        if varname in _NO_MULTI and self.index != 0:
                                                                                                                      
            if varname in self.dt0:
                glob_value = np.array(self.dt0[varname])
                isok  = np.all(glob_value == np.array(value))
                if not isok:
                    err_msg = "NO_MULTI variable: dataset 0: %s, dataset %d: %s" % (str(glob_value), self.index, str(value))
                    raise self.Error(err_msg)
            else:
                self.dt0.set_variable(varname, value)

    def set_variables(self, **vars):
        """Sets variables by providing a dictionary"""
        for (varname, varvalue) in vars.items():
            self.set_variable(varname, varvalue)

    @property
    def structure(self):
        try:
            return self._structure

        except AttributeError:
            try:
                structure = Structure.from_abivars(self.allvars)
            except Exception as exc:
                raise
                print(str(exc))
                structure = None

            self.set_structure(structure)
            return self._structure

    def set_structure(self, structure):
        if hasattr(self, "_structure") and self._structure != structure:
            raise self.Error("Dataset object already has a structure object, cannot overwrite")

        self._structure = structure
        if structure is None:
            return

        self.set_variables(**structure.to_abivars())

        #if self.index == 0:
        #    self.set_variables(**structure.to_abivars())

        #if self.index != 0 and self.dt0.structure is not None and structure != self.dt0.structure:
        #    self.set_variables(**structure.to_abivars())

    # Helper functions to facilitate the specification of several variables.
    def set_kmesh(self, ngkpt, shiftk, kptopt=1):
        """
        Set the variables for the sampling of the BZ.

        Args:
            ngkpt:
                Monkhorst-Pack divisions
            shiftk:
                List of shifts.
            kptopt:
                Option for the generation of the mesh.
        """
        shiftk = np.reshape(shiftk, (-1,3))
        
        self.set_variables(ngkpt=ngkpt,
                           kptopt=kptopt,
                           nshiftk=len(shiftk),
                           shiftk=shiftk,
                           )

    def set_kpath(self, ndivsm, kptbounds=None):
        """
        Set the variables for the computation of the band structure.

        Args:
            ndivsm:
                Number of divisions for the smallest segment.
            kptbounds
                k-points defining the path in k-space.
        """
        if kptbounds is None:
            kptbounds = [k.frac_coords for k in self.structure.hsym_kpoints]

        kptbounds = np.reshape(kptbounds, (-1,3))

        self.set_variables(kptbounds=kptbounds,
                           kptopt=-(len(kptbounds)-1),
                           ndivsm=ndivsm,
                           )

    def set_kptgw(self, kptgw, bdgw):
        """
        Set the variables (k-points, bands) for the computation of the GW corrections.

        Args
            kptgw:
                List of k-points in reduced coordinates.
            bdgw:
                Specifies the range of bands for the GW corrections.
                Accepts iterable that be reshaped to (nkptgw, 2) 
                or a tuple of two integers if the extrema are the same for each k-point.
        """
        kptgw = np.reshape(kptgw, (-1,3))
        nkptgw = len(kptgw)

        if len(bdgw) == 2:
            bdgw = len(kptgw) * bdgw

        self.set_variables(kptgw=kptgw,
                           nkptgw=nkptgw,
                           bdgw=np.reshape(bdgw, (nkptgw, 2)),
                           )


def format_var(varname, value):
    """
    Return a string for a single input variable to be written in the input file.
    Always ends with a carriage return.
    """
    if value is None or not str(value):
        return ""

    # By default, do not impose a number of decimal points ...but not for those
    digits = {
        ('xred','xcart','qpt','kpt') : 10,
        ('ngkpt','kptrlatt', 'ngqpt'): 0,
    }

    precision = {}
    for t, numd in digits.items():
        for varname in t:
            precision[varname] = numd 

    floatdecimal = precision.get(varname, 0)

    line = " "

    if isinstance(value, np.ndarray):
        n = 1
        for i in np.shape(value):
            n *= i
        value = list(np.reshape(value, n))

    if isinstance(value, (list, tuple)):
        # values in lists
        # Reshape a list of lists into a single list
        if all(isinstance(v, list) for v in value):
            line += format_list_of_list(value, floatdecimal)

        else:
            # Maximum number of values per line.
            valperline = 3
            if any(inp in varname for inp in ['bdgw']):
                valperline = 2

            line += format_list(value, valperline, floatdecimal)

    else:
        # scalar values
        line += ' ' + str(value)

    return line


def format_scalar(val, floatdecimal=0):
    """
    Format a single numerical value into a string
    with the appropriate number of decimals.
    """
    sval = str(val)
    if sval.isdigit() and floatdecimal == 0:
        return sval

    try:
        fval = float(val)
    except:
        return sval

    if fval == 0 or (abs(fval) > 1e-3 and abs(fval) < 1e4):
        form = 'f'
        addlen = 5
    else:
        form = 'e'
        addlen = 8

    ndec = max(len(str(fval-int(fval)))-2, floatdecimal)
    ndec = min(ndec, 10)

    sval = '{v:>{l}.{p}{f}}'.format(v=fval, l=ndec+addlen, p=ndec, f=form)

    return sval.replace('e', 'd')


def format_list_of_list(values, floatdecimal=0):
    """Format a list of lists."""

    lvals = flattens(values)

    # Determine the representation
    if all(isinstance(v, int) for v in lvals):
        type_all = int
    else:
        try:
            for v in lvals:
                float(v)
            type_all = float
        except:
            type_all = str

    # Determine the format
    width = max(len(str(s)) for s in lvals)
    if type_all == int:
        formatspec = '>{}d'.format(width)

    elif type_all == str:
        formatspec = '>{}'.format(width)

    else:
        # Number of decimal
        maxdec = max(len(str(f-int(f)))-2 for f in lvals)
        ndec = min(max(maxdec, floatdecimal), 10)

        if all(f == 0 or (abs(f) > 1e-3 and abs(f) < 1e4) for f in lvals):
            formatspec = '>{w}.{p}f'.format(w=ndec+5, p=ndec)
        else:
            formatspec = '>{w}.{p}e'.format(w=ndec+8, p=ndec)

    line = '\n'
    for L in values:
        for val in L:
            line += ' {v:{f}}'.format(v=val, f=formatspec)
        line += '\n'

    return line


def format_list(values, valperline=3, floatdecimal=0):
    """
    Format a list of values into a string.  The result might be spread among
    several lines, each line starting with a blank space.
    'valperline' specifies the number of values per line.
    """
    line = ""
    # Format the line declaring the value
    for i, val in enumerate(values):
        line += ' ' + format_scalar(val, floatdecimal)
        if (i+1) % valperline == 0:
            line += '\n'

    # Add a carriage return in case of several lines
    if '\n' in line.rstrip('\n'):
        line = '\n' + line

    # Add a carriage return at the end
    if not line.endswith('\n'):
        line += '\n'

    return line


def string_to_value(sval):
    """
    Interpret a string variable and attempt to return a value of the appropriate type.  
    If all else fails, return the initial string.
    """
    value = None
    try:
        for part in sval.split():
            if '*' in part:
                if part[0] == '*':
                    # cases like istwfk *1
                    value = None
                    break

                else:
                    # cases like acell 3*3.52
                    n = int(part.split('*')[0])
                    f = convert_number(part.split('*')[1])
                    if value is None:
                        value = []
                    value += n * [f]
                    continue

            if '/' in part:
                # Fractions
                (num, den) = (float(part.split('/')[i]) for i in range(2))
                part = num / den

            if part in _UNITS.keys():
                # Unit
                if value is None:
                    msg = "Could not apply the unit tolken '%s'." %(part)
                    warnings.warn(msg)
                elif isinstance(value, list):
                    value.append(part)
                else:
                    value = [value, part]

                # Convert
                if False:
                    if isinstance(value, list):
                        for i in range(len(value)):
                            value[i] *= _UNITS[part]
                    elif isinstance(value, str):
                        value = None
                        break
                    else:
                        value *= _UNITS[part]

                continue

            # Convert
            try:
                val = convert_number(part)
            except:
                val = part

            if value is None:
                value = val
            elif isinstance(value, list):
                value.append(val)
            else:
                value = [value, val]
    except:
        value = None

    if value is None:
        value = sval

    return value


def format_variable_dict(variables):
    """Format a dictionary of input variables into a multi-line string."""
    s = ""
    for (var, val) in variables.items():
        s += format_var(var, val)

    return s


def flattens(lists):
    """Transform a list of lists into a single list."""
    if not isinstance(lists[0], list):
        return lists
    lst = []
    for l in lists:
        lst += l

    return lst


def listify(obj):
    """Transform any object, iterable or not, to a list."""
    if isinstance(obj, collections.Iterable):
        return obj

    elif isinstance(obj, collections.Iterator):
        return list(obj)

    else:
        return [obj]


def is_number(s):
    """Returns True if the argument can be converted to float."""
    try:
        float(s)
        return True
    except:
        return False


def convert_number(value):
    """
    Converts some object to a float or a string.
    If the argument is an integer or a float, returns the same object.
    If the argument is a string, tries to convert to an integer, then to a float.
    The string '1.0d-03' will be treated the same as '1.0e-03'
    """
    if isinstance(value (float, int)):
        return value

    elif isinstance(value, str):

        if is_number(value):
            try:
                return int(value)
            except ValueError:
                return float(value)

        else:
            val = value.replace('d', 'e')
            if is_number(val):
                return float(val)
            else:
                raise ValueError("convert_number failed")

    else:
        raise ValueError("convert_number failed")
