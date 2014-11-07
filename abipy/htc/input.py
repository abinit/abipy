"""
This module defines objects that faciliate the creation of the 
ABINIT input files. The syntax is similar to the one used 
in ABINIT with small differences. 
"""
from __future__ import print_function, division, unicode_literals

import os
import collections
import warnings
import itertools
import copy
import six
import abc
import numpy as np
import pymatgen.core.units as units
import abipy.tools.mixins as mixins

from collections import OrderedDict
from monty.string import is_string, list_strings
from pymatgen.io.abinitio.pseudos import PseudoTable, Pseudo
from pymatgen.core.units import Energy
from abipy.core import Structure
from abipy.htc.variable import InputVariable
from abipy.htc.abivars import is_abivar, is_anaddb_var

import logging
logger = logging.getLogger(__file__)


__all__ = [
    "AbinitInputError",
    "AbiInput",
    "LdauParams",
    "LexxParams",
    "input_gen",
    "AnaddbInput",
]

# Variables that must have a unique value throughout all the datasets.
_ABINIT_NO_MULTI = [
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


def straceback():
    """Returns a string with the traceback."""
    import traceback
    return traceback.format_exc()


@six.add_metaclass(abc.ABCMeta)
class Input(object): 
    """
    Base class foor Input objects.

    An input object must define have a make_input method 
    that returns the string representation used by the client code.
    """
    def copy(self):
        """Shallow copy of the input."""
        return copy.copy(self)
                                   
    def deepcopy(self):
        """Deep copy of the input."""
        return copy.deepcopy(self)

    def write(self, filepath):
        """
        Write the input file to file. Returns a string with the input.
        """
        dirname = os.path.dirname(filepath)
        if not os.path.exists(dirname): os.makedirs(dirname)
                                                                                      
        # Write the input file.
        input_string = str(self)
        with open(filepath, "w") as fh:
            fh.write(input_string)

        return input_string

    #@abc.abstractmethod
    #def make_input(self)

    #@abc.abstractproperty
    #def structure(self):

    #def set_structure(self, structure):

    #def set_variables(self, dtset=0, **vars):
    #    """
    #    Set the value of a set of input variables.

    #    Args:
    #        dtset:
    #            Int with the index of the dataset, slice object of iterable 
    #        vars:
    #            Dictionary with the variables.
    #    """
    #    for idt in self._dtset2range(dtset):
    #        self[idt].set_variables(**vars)

    #def remove_variables(self, keys, dtset=0):
    #    """
    #    Remove the variable listed in keys
    #                                                                         
    #    Args:
    #        dtset:
    #            Int with the index of the dataset, slice object of iterable 
    #        keys:
    #            List of variables to remove.
    #    """
    #    for idt in self._dtset2range(dtset):
    #        self[idt].remove_variables(keys)


class AbinitInputError(Exception):
    """Base error class for exceptions raised by `AbiInput`"""


class AbiInput(Input):
    """
    This object represents an ABINIT input file. It supports multi-datasets a
    and provides an easy-to-use interface for defining the variables of the calculation.

    Usage example:

    .. code-block:: python

        inp = AbiInput(pseudos="si.pspnc", pseudo_dir="directory_with_pseudos")

        inp.set_structure_from_file("si.cif")
        inp.set_variables(ecut=10, nband=3)
        # Add other variables. See the other methods provided by this object.

        print(inp)
    """
    Error = AbinitInputError

    def __init__(self, pseudos, pseudo_dir="", ndtset=1, comment=""):
        """
        Args:
            pseudos:
                String or list of string with the name of the pseudopotential files.
            pseudo_dir:
                Name of the directory where the pseudopotential files are located.
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
        if isinstance(pseudos, PseudoTable):
            self._pseudos = pseudos

        elif all(isinstance(p, Pseudo) for p in pseudos):
            self._pseudos = PseudoTable(pseudos)

        else:
            # String(s)
            pseudo_dir = os.path.abspath(pseudo_dir)

            pseudo_paths = [os.path.join(pseudo_dir, p) for p in list_strings(pseudos)]

            missing = [p for p in pseudo_paths if not os.path.exists(p)]
            if missing:
                raise self.Error("Cannot find the following pseudopotential files:\n%s" % str(missing)) 

            try:
                self._pseudos = PseudoTable(pseudo_paths)

            except Exception as exc:
                msg = "\nIgnoring error raised while parsing pseudopotential files:\n Backtrace:" + straceback()
                warnings.warn(msg)
                self._pseudos = []

        if comment:
            self.set_comment(comment)

    #def make_input(self):
    #    return str(self)

    def __str__(self):
        """String representation i.e. the input file read by Abinit."""
        if self.ndtset > 1:
            s = ""
            for (i, dataset) in enumerate(self):
                header = "### DATASET %d ###" % i
                if i == 0: 
                    header = "### GLOBAL VARIABLES ###"

                str_dataset = str(dataset)
                if str_dataset:
                    header = len(header) * "#" + "\n" + header + "\n" + len(header) * "#" + "\n"
                    s += "\n" + header + str(dataset) + "\n"
        else:
            # single datasets ==> don't append the dataset index to the variables in self[1]
            # this trick is needed because Abinit complains if ndtset is not specified 
            # and we have variables that end with the dataset index e.g. acell1
            # We don't want to specify ndtset here since abinit will start to add DS# to 
            # the input and output files thus complicating the algorithms we have to use
            # to locate the files.
            d = self[0].deepcopy()
            d.update(self[1])
            s = d.to_string(post="")

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
                if not (self.ndtset >= idt >= 1):
                    raise self.Error("Invalid dataset index %d, ndtset %d " % (idt, self.ndtset))

                # Strip the numeric index.
                varname = varname[:len(varname)-i]

                self[idt].set_variable(varname, value)

                # Handle no_multi variables.
                if varname in _ABINIT_NO_MULTI:

                    if varname in self[0]:
                        glob_value = np.array(self[0][varname])
                        isok = np.all(glob_value == np.array(value))
                        if not isok:
                            err_msg = "NO_MULTI variable: dataset 0: %s, dataset %d: %s" % (str(glob_value), idt, str(value))
                            raise self.Error(err_msg)
                    else:
                        self[0].set_variable(varname, value)

    @property
    def ndtset(self):
        """Number of datasets."""
        return self[0]["ndtset"]

    @property
    def pseudos(self):
        """List of `Pseudo` objecst."""
        return self._pseudos

    @property
    def ispaw(self):
        """True if we have a PAW calculation."""
        return all(p.ispaw for p in self.pseudos)

    @property
    def isnc(self):
        """True if we have a norm-conserving calculation."""
        return all(p.isnc for p in self.pseudos)

    @property
    def structure(self):
        """
        Returns the `Structure` associated to the ABINIT input.

        Raises:
            ValueError if we have multi datasets with different 
            structure so that users will learn that multiple datasets are bad.
        """
        structures = [dt.structure for dt in self[1:]]

        if all(s == structures[0] for s in structures):
            return structures[0]

        raise ValueError("Cannot extract a unique structure from an input file with multiple datasets!\n" + 
                         "Please DON'T use multiple datasets with different unit cells!")

    def is_abivar(self, varname):
        """True if varname is a valid Abinit variable."""
        return is_abivar(varname)

    def _dtset2range(self, dtset):
        """
        Helper function to convert dtset into a range. Accepts int, slice objects or iterable.
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

    def split_datasets(self):
        """
        Split an input file with multiple datasets into a  list of `ndtset` distinct input files.
        """
        # Propagate subclasses (if any)
        cls = self.__class__ 
        news = []

        for i in range(self.ndtset):
            my_vars = self[i+1].allvars
            my_vars.pop("ndtset", None)

            # Cannot use get*, ird* variables since links must be explicit.
            for varname in my_vars:
                if varname.startswith("get") or varname.startswith("ird"):
                    err_msg = ("get* or ird* variables should not be present in the input when you split it into datasets")
                    raise self.Error(err_msg)

            new = cls(pseudos=self.pseudos, ndtset=1)
            new.set_variables(**my_vars)
            #new.set_geoformat(self.geoformat)
            news.append(new)
    
        return news

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

    def remove_variables(self, keys, dtset=0):
        """
        Remove the variable listed in keys
                                                                             
        Args:
            dtset:
                Int with the index of the dataset, slice object of iterable 
            keys:
                List of variables to remove.
        """
        for idt in self._dtset2range(dtset):
            self[idt].remove_variables(keys)

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

    #@property
    #def geoformat(self):
    #    """
    #    angdeg if the crystalline structure should be specified with angdeg and acell, 
    #    rprim otherwise (default format)
    #    """
    #    try:
    #        return self._geoformat
    #    except AttributeError
    #        return "rprim" # default

    #def set_geoformat(self, format):
    #    assert format in ["angdeg", "rprim"]
    #    self._geoformat = format

    def linspace(self, varname, start, stop, endpoint=True):
        """
        Returns `ndtset` evenly spaced samples, calculated over the interval [`start`, `stop`].

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

        .. code-block:: python

            inp.product("tsmear", "ngkpt", [[2,2,2], [4,4,4]], [0.1, 0.2, 0.3])
        """
        # Split items into varnames and values
        for i, item in enumerate(items):
            if not is_string(item):
                break

        varnames, values = items[:i], items[i:]
        if len(varnames) != len(values):
            raise self.Error("The number of variables must equal the number of lists")

        varnames = [ [varnames[i]] * len(values[i]) for i in range(len(values))]
        varnames = itertools.product(*varnames)
        values = itertools.product(*values)

        idt = 0
        for names, values in zip(varnames, values):
            idt += 1
            self[idt].set_variables(**{k: v for k, v in zip(names, values)})
        
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

    def set_autokmesh(self, nksmall, kptopt=1, dtset=0):
        """
        Set the variables (ngkpt, shift, kptopt) for the sampling of the BZ.
                                                       
        Args:
            nksmall:
                Number of k-points used to sample the smallest lattice vector.
            kptopt:
                Option for the generation of the mesh.
        """
        for idt in self._dtset2range(dtset):
            self[idt].set_autokmesh(nksmall, kptopt=kptopt)

    def set_autokpath(self, ndivsm, dtset=0):
        for idt in self._dtset2range(dtset):
            self[idt].set_kpath(ndivsm, kptbounds=None)

    def set_kpath(self, ndivsm, kptbounds=None, iscf=-2, dtset=0):
        """
        Set the variables defining the k-path for the specified dtset.

        The list of K-points is taken from the pymatgen database if kptbounds is None

        Args:
            ndivsm:
                Number of divisions for the smallest segment.
            kptbounds:
                k-points defining the path in k-space.
            iscf:
                iscf variable.
            dtset:
                Index of the dataset. 0 for global variables.
        """
        for idt in self._dtset2range(dtset):
            self[idt].set_kpath(ndivsm, kptbounds=kptbounds, iscf=iscf)

    def set_kptgw(self, kptgw, bdgw, dtset=0):
        """
        Set the variables defining the k point for the GW corrections
        """
        for idt in self._dtset2range(dtset):
            self[idt].set_kptgw(kptgw, bdgw)


class Dataset(mixins.MappingMixin):
    """
    This object stores the ABINIT variables for a single dataset.
    """
    # TODO this should be the "actual" input file
    Error = AbinitInputError

    def __init__(self, index, dt0, *args, **kwargs):
        self._mapping_mixin_ = collections.OrderedDict(*args, **kwargs)
        self._index = index
        self._dt0 = dt0

    def __repr__(self):
        "<%s at %s>" % (self.__class__.__name__, id(self))

    def __str__(self):
        return self.to_string()

    @property
    def vars(self):
        return self._mapping_mixin_

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
        return all_vars.copy()

    def deepcopy(self):
        """Deepcopy of the `Dataset`"""
        return copy.deepcopy(self)

    #@property
    #def geoformat(self):
    #    """
    #    angdeg if the crystalline structure should be specified with angdeg and acell, 
    #    rprim otherwise (default format)
    #    """
    #    try:
    #        return self._geoformat
    #    except AttributeError
    #        return "rprim" # default
    #                                                                                    
    #def set_geoformat(self, format):
    #    assert format in ["angdeg", "rprim"]
    #    self._geoformat = format

    def to_string(self, sortmode=None, post=None):
        """
        String representation.

        Args:
            sortmode:
                "a" for alphabetical order, None if no sorting is wanterd
            post:
                String that will be appended to the name of the variables
                Note that post is usually autodetected when we have multiple datatasets
                It is mainly used when we have an input file with a single dataset
                so that we can prevent the code from adding "1" to the name of the variables 
                (In this case, indeed, Abinit complains if ndtset=1 is not specified 
                and we don't want ndtset=1 simply because the code will start to add 
                _DS1_ to all the input and output files.
        """
        lines = []
        app = lines.append

        if self.comment:
            app("# " + self.comment.replace("\n", "\n#"))

        if post is None:
            post = "" if self.index == 0 else str(self.index)
        else:
            post = post

        if sortmode is None:
            # no sorting.
            keys = self.keys()
        elif sortmode == "a":
            # alphabetical order.
            keys = sorted(self.keys())
        else:
            raise ValueError("Unsupported value for sortmode %s" % str(sortmode))

        for var in keys:
            value = self[var]
            # Do not print NO_MULTI variables except for dataset 0.
            if self.index != 0 and var in _ABINIT_NO_MULTI:
                continue

            # Print ndtset only if we really have datasets. 
            # This trick prevents abinit from adding DS# to the output files.
            # thus complicating the creation of workflows made of separated calculations.
            if var == "ndtset" and value == 1:
                continue

            varname = var + post
            variable = InputVariable(varname, value)
            app(str(variable))

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
        """Set a the value of a variable."""
        # Check if varname is in the internal database.
        if not is_abivar(varname):
            raise self.Error("%s is not a valid ABINIT variable." % varname)

        if varname in self:
            try: 
                iseq = (self[varname] == value)
                iseq = np.all(iseq)
            except ValueError:
                # array like.
                iseq = np.allclose(self[varname], value)
            except:
                iseq = False

            if not iseq:
                msg = "%s is already defined with a different value:\nOLD:\n %s,\nNEW\n %s" % (
                    varname, str(self[varname]), str(value))
                warnings.warn(msg)

        self[varname] = value

        # Handle no_multi variables.
        if varname in _ABINIT_NO_MULTI and self.index != 0:
                                                                                                                      
            if varname in self.dt0:
                glob_value = np.array(self.dt0[varname])
                isok = np.all(glob_value == np.array(value))
                if not isok:
                    err_msg = "NO_MULTI variable: dataset 0: %s, dataset %d: %s" % (
                        str(glob_value), self.index, str(value))
                    raise self.Error(err_msg)
            else:
                self.dt0.set_variable(varname, value)

    def set_variables(self, **vars):
        """Set the value of the variables provied in the dictionary **vars"""
        for varname, varvalue in vars.items():
            self.set_variable(varname, varvalue)

    def remove_variables(self, keys):
        """Remove the variables listed in keys."""
        for key in list_strings(keys):
            if key not in self:
                raise KeyError("key: %s not in self:\n %s" % (key, list(self.keys())))
            #self.pop(key, None)
            self.pop(key)

    @property
    def structure(self):
        """Returns the `Structure` associated to this dataset."""
        try:
            return self._structure

        except AttributeError:
            structure = Structure.from_abivars(self.allvars)

            self.set_structure(structure)
            return self._structure

    def set_structure(self, structure):
        #if hasattr(self, "_structure") and self._structure != structure:
        #    raise self.Error("Dataset object already has a structure object, cannot overwrite")

        self._structure = structure
        if structure is None:
            return

        geoformat = "rprim"
        self.set_variables(**structure.to_abivars())

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
                           shiftk=shiftk)

    def set_autokmesh(self, nksmall, kptopt=1):
        """
        Set the variables (ngkpt, shift, kptopt) for the sampling of the BZ.
                                                       
        Args:
            nksmall:
                Number of k-points used to sample the smallest lattice vector.
            kptopt:
                Option for the generation of the mesh.
        """
        shiftk = self.structure.calc_shiftk()
        ngkpt = self.structure.calc_ngkpt(nksmall)
        
        self.set_variables(ngkpt=ngkpt,
                           kptopt=kptopt,
                           nshiftk=len(shiftk),
                           shiftk=shiftk)

    def set_kpath(self, ndivsm, kptbounds=None, iscf=-2):
        """
        Set the variables for the computation of the band structure.

        Args:
            ndivsm:
                Number of divisions for the smallest segment.
            kptbounds
                k-points defining the path in k-space.
                If None, we use the default high-symmetry k-path 
                defined in the pymatgen database.
        """
        if kptbounds is None:
            kptbounds = self.structure.calc_kptbounds()

        kptbounds = np.reshape(kptbounds, (-1,3))

        self.set_variables(kptbounds=kptbounds,
                           kptopt=-(len(kptbounds)-1),
                           ndivsm=ndivsm,
                           iscf=iscf)

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
                           bdgw=np.reshape(bdgw, (nkptgw, 2)))


class LujForSpecie(collections.namedtuple("LdauForSpecie", "l u j unit")):
    """
    This object stores the value of l, u, j used for a single atomic specie.
    """
    def __new__(cls, l, u, j, unit):
        """
        Args:
            l: 
                Angular momentum (int or string).
            u:
                U value
            j:
                J Value
            unit:
                Energy unit for u and j.
        """
        l = l
        u = Energy(u, unit)
        j = Energy(j, unit)
        return super(cls, LujForSpecie).__new__(cls, l, u, j, unit)


class LdauParams(object):
    """
    This object stores the parameters for LDA+U calculations with the PAW method
    It facilitates the specification of the U-J parameters in the Abinit input file.
    (see `to_abivars`). The U-J operator will be applied only on  the atomic species 
    that have been selected by calling `lui_for_symbol`.

    To setup the Abinit variables for a LDA+U calculation in NiO with a 
    U value of 5 eV applied on the nickel atoms:

    .. code-block:: python

        luj_params = LdauParams(usepawu=1, structure=nio_structure)
        # Apply U-J on Ni only.
        u = 5.0 
        luj_params.luj_for_symbol("Ni", l=2, u=u, j=0.1*u, unit="eV")

        print(luj_params.to_abivars())
    """
    def __init__(self, usepawu, structure):
        """
        Arg:
            usepawu:
                Abinit variable `usepawu` defining the LDA+U method.
            structure:
                `Structure` object.
        """
        self.usepawu = usepawu
        self.structure = structure
        self._params = {} 

    @property
    def symbols_by_typat(self):
        return [specie.symbol for specie in self.structure.types_of_specie]

    def luj_for_symbol(self, symbol, l, u, j, unit="eV"):
        """
        Args:
            symbol:
                Chemical symbol of the atoms on which LDA+U should be applied.
            l:
                Angular momentum.
            u: 
                Value of U.
            j:
                Value of J.
            unit:
                Energy unit of U and J.
        """
        if symbol not in self.symbols_by_typat:
            err_msg = "Symbol %s not in symbols_by_typat:\n%s" % (symbol, self.symbols_by_typat)
            raise ValueError(err_msg)

        if symbol in self._params:
            err_msg = "Symbol %s is already present in LdauParams! Cannot overwrite:\n" % symbol
            raise ValueError(err_msg)

        self._params[symbol] = LujForSpecie(l=l, u=u, j=j, unit=unit)

    def to_abivars(self):
        """
        Returns a dict with the Abinit variables.
        """
        lpawu, upawu, jpawu = [], [], []

        for symbol in self.symbols_by_typat:
            p = self._params.get(symbol, None)

            if p is not None:
                l, u, j = p.l, p.u.to("eV"), p.j.to("eV")
            else:
                l, u, j = -1, 0, 0

            lpawu.append(int(l)) 
            upawu.append(float(u))
            jpawu.append(float(j))

        # convert upawu and jpaw to string so that we can use 
        # eV unit in the Abinit input file (much more readable).
        return dict(
            usepawu=self.usepawu,
            lpawu=" ".join(map(str, lpawu)),
            upawu=" ".join(map(str, upawu)) + " eV",
            jpawu=" ".join(map(str, jpawu)) + " eV",
        )


class LexxParams(object):
    """
    This object stores the parameters for local exact exchange calculations with the PAW method
    It facilitates the specification of the LEXX parameters in the Abinit input file.
    (see `to_abivars`). The LEXX operator will be applied only on the atomic species 
    that have been selected by calling `lexx_for_symbol`.

    To perform a LEXX calculation for NiO in which the LEXX is compute only for the l=2
    channel of the nickel atoms:
                                                                         
    .. code-block:: python

        lexx_params = LexxParams(nio_structure)
        lexx_params.lexx_for_symbol("Ni", l=2)                                                                         

        print(lexc_params.to_abivars())
    """
    def __init__(self, structure):
        """
        Arg:
            structure:
                `Structure` object.
        """
        self.structure = structure
        self._lexx_for_symbol = {} 

    @property
    def symbols_by_typat(self):
        return [specie.symbol for specie in self.structure.types_of_specie]

    def lexx_for_symbol(self, symbol, l):
        """
        Enable LEXX for the given chemical symbol and the angular momentum l 

        Args:
            symbol:
                Chemical symbol of the atoms on which LEXX should be applied.
            l:
                Angular momentum.
        """
        if symbol not in self.symbols_by_typat:
            err_msg = "Symbol %s not in symbols_by_typat:\n%s" % (symbol, self.symbols_by_typat)
            raise ValueError(err_msg)

        if symbol in self._lexx_for_symbol:
            err_msg = "Symbol %s is already present in LdauParams! Cannot overwrite:\n" % symbol
            raise ValueError(err_msg)

        self._lexx_for_symbol[symbol] = l

    def to_abivars(self):
        """
        Returns a dict with the Abinit variables.
        """
        lexx_typat = []

        for symbol in self.symbols_by_typat:
            l = self._lexx_for_symbol.get(symbol, -1)
            lexx_typat.append(int(l)) 

        return dict(useexexch=1, lexexch=" ".join(map(str, lexx_typat)))


def input_gen(inp, **kwargs):
    """
    This function receives an `AbiInput` and generates
    new inputs by replacing the variables specified in kwargs.

    Args:
        inp:
            `AbiInput` file.
        kwargs:
            keyword arguments with the values used for each variable.

    .. code-block:: python

        gs_inp = call_function_to_generate_initial_template()

        # To generate two input files with different values of ecut:
        for inp_ecut in input_gen(gs_inp, ecut=[10, 20]):
            print("do something with inp_ecut %s" % inp_ecut)

        # To generate four input files with all the possible combinations of ecut and nsppol:
        for inp_ecut in input_gen(gs_inp, ecut=[10, 20], nsppol=[1, 2]):
            print("do something with inp_ecut %s" % inp_ecut)
    """
    for new_vars in product_dict(kwargs):
        new_inp = inp.deepcopy()
        # Remove the variable names to avoid annoying warnings.
        # if the variable is overwritten.
        new_inp.remove_variables(new_vars.keys())
        new_inp.set_variables(**new_vars)

        yield new_inp


def product_dict(d):
    """
    This function receives a dictionary d where each key defines a list of items or a simple scalar.
    It constructs the Cartesian product of the values (equivalent to nested for-loops),
    and returns a list of dictionaries with the values that would be used inside the loop.

    >>> d = OrderedDict([("foo", [2, 4]), ("bar", 1)])
    >>> product_dict(d)
    [OrderedDict([('foo', 2), ('bar', 1)]), OrderedDict([('foo', 4), ('bar', 1)])]
    >>> d = OrderedDict([("bar", [1,2]), ('foo', [3,4])]) 
    >>> product_dict(d) == [{'bar': 1, 'foo': 3},
    ... {'bar': 1, 'foo': 4},
    ... {'bar': 2, 'foo': 3},
    ... {'bar': 2, 'foo': 4}]
    True

    .. warning:
        Dictionaries are not ordered, therefore one cannot assume that 
        the order of the keys in the output equals the one used to loop.
        If the order is important, one should pass a `OrderedDict` in input
    """
    keys, vals = d.keys(), d.values()

    # Each item in vals must be iterable.
    values = []

    for v in vals:
        if not isinstance(v, collections.Iterable): v = [v]
        values.append(v)

    # Build list of dictionaries. Use ordered dicts so that 
    # we preserve the order when d is an OrderedDict.
    vars_prod = [] 

    for prod_values in itertools.product(*values):
        dprod = OrderedDict(zip(keys, prod_values))
        vars_prod.append(dprod)

    return vars_prod


class AnaddbInputError(Exception):
    """Exceptions raised by AnaddbInput."""


class AnaddbInput(mixins.MappingMixin):
    #TODO: Abstract interface to that we can provide
    # tools for AbinitInput and AnaddbInput
    #deepcopy
    #removevariable
    Error = AnaddbInputError

    def __init__(self, structure, comment="", **kwargs):
        """
        Args:
            structure:
                Crystalline structure.
            comment:
                Optional string with a comment that will be placed at the beginning of the file.
        """
        self._structure = structure
        self.comment = comment

        for k in kwargs:
            if not self.is_anaddb_var(k):
                raise self.Error("%s is not a registered Anaddb variable" % k)

        self._mapping_mixin_ = collections.OrderedDict(**kwargs)

    @property
    def vars(self):
        return self._mapping_mixin_

    @classmethod
    def phbands_and_dos(cls, structure, ngqpt, ndivsm, nqsmall, q1shft=(0,0,0),
                        qptbounds=None, asr=1, chneut=0, dipdip=1, dos_method="tetra", **kwargs):
        """
        Build an anaddb input file for the computation of phonon bands and phonon DOS.

        Args:
            Structure:
                Structure object
            ngqpt:
                Monkhorst-Pack divisions for the phonon Q-mesh (coarse one)
            ndivsm:
                Used to generate a normalized path for the phonon bands.
                If gives the number of divisions for the smallest segment of the path.
            nqsmall:
                Used to generate the (dense) mesh for the DOS.
                It defines the number of q-points used to sample the smallest lattice vector.
            q1shft:
                Shifts used for the coarse Q-mesh
            qptbounds
                Boundaries of the path. If None, the path is generated from an internal database
                depending on the input structure.
            asr, chneut, dipdp:
                Anaddb input variable. See official documentation.
            dos_method:
                Possible choices: "tetra", "gaussian" or "gaussian:0.001 eV".
                In the later case, the value 0.001 eV is used as gaussian broadening
        """
        dosdeltae, dossmear = None, None

        if dos_method == "tetra":
            prtdos = 2
        elif "gaussian" in dos_method:
            prtdos = 1
            i = dos_method.find(":")
            if i != -1:
                value, eunit = dos_method[i+1:].split()
                dossmear = units.Energy(float(value), eunit).to("Ha")
        else:
            raise cls.Error("Wrong value for dos_method: %s" % dos_method)

        new = cls(structure, comment="ANADB input for phonon bands and DOS", **kwargs)

        new.set_qpath(ndivsm, qptbounds=qptbounds)
        new.set_autoqmesh(nqsmall)

        q1shft = np.reshape(q1shft, (-1, 3))

        new.set_variables(
            ifcflag=1,
            ngqpt=np.array(ngqpt),
            q1shft=q1shft,
            nqshft=len(q1shft),
            asr=asr,
            chneut=chneut,
            dipdip=dipdip,
            prtdos=prtdos,
            dosdeltae=dosdeltae,
            dossmear=dossmear,
        )

        return new

    @classmethod
    def thermo(cls, structure, ngqpt, nqsmall, q1shft=(0, 0, 0), nchan=1250, nwchan=5, thmtol=0.5,
               ntemper = 79, temperinc = 5, tempermin = 5., asr=1, chneut=1, dipdip=1, ngrids=10, **kwargs):
        """
        Build an anaddb input file for the computation of phonon bands and phonon DOS.

        Args:
            Structure:
                Structure object
            ngqpt:
                Monkhorst-Pack divisions for the phonon Q-mesh (coarse one)
            nqsmall:
                Used to generate the (dense) mesh for the DOS.
                It defines the number of q-points used to sample the smallest lattice vector.
            q1shft:
                Shifts used for the coarse Q-mesh
            qptbounds
                Boundaries of the path. If None, the path is generated from an internal database
                depending on the input structure.
            asr, chneut, dipdp:
                Anaddb input variable. See official documentation.
            dos_method:
                Possible choices: "tetra", "gaussian" or "gaussian:0.001 eV".
                In the later case, the value 0.001 eV is used as gaussian broadening

            #!Flags
             # ifcflag   1     ! Interatomic force constant flag
             # thmflag   1     ! Thermodynamical properties flag
            #!Wavevector grid number 1 (coarse grid, from DDB)
            #  brav    2      ! Bravais Lattice : 1-S.C., 2-F.C., 3-B.C., 4-Hex.)
             #  ngqpt   4  4  4   ! Monkhorst-Pack indices
             #  nqshft  1         ! number of q-points in repeated basic q-cell
             #  q1shft  3*0.0
            #!Effective charges
             #     asr   1     ! Acoustic Sum Rule. 1 => imposed asymetrically
             #  chneut   1     ! Charge neutrality requirement for effective charges.
            #!Interatomic force constant info
             #  dipdip  1      ! Dipole-dipole interaction treatment
            #!Wavevector grid number 2 (series of fine grids, extrapolated from interat forces)
            #  ng2qpt   20 20 20  ! sample the BZ up to ngqpt2
            #  ngrids   5         ! number of grids of increasing size#  q2shft   3*0.0
            #!Thermal information
             #  nchan   1250   ! # of channels for the DOS with channel width 1 cm-1
             #  nwchan  5      ! # of different channel widths from this integer down to 1 cm-1
             #  thmtol  0.120  ! Tolerance on thermodynamical function fluctuations
             #  ntemper 10     ! Number of temperatures
             #  temperinc 20.  ! Increment of temperature in K for temperature dependency
             #  tempermin 20.  ! Minimal temperature in Kelvin
            # This line added when defaults were changed (v5.3) to keep the previous, old behaviour
            #  symdynmat 0

        """

        new = cls(structure, comment="ANADB input for themodynamics", **kwargs)

        #new.set_qpath(ndivsm, qptbounds=qptbounds)
        new.set_autoqmesh(nqsmall)

        q1shft = np.reshape(q1shft, (-1, 3))

        new.set_variables(
            ifcflag=1,
            thmflag=1,
            ngqpt=np.array(ngqpt),
            ngrids=ngrids,
            q1shft=q1shft,
            nqshft=len(q1shft),
            asr=asr,
            chneut=chneut,
            dipdip=dipdip,
            nchan=nchan,
            nwchan=nwchan,
            thmtol=thmtol,
            ntemper=ntemper,
            temperinc=temperinc,
            tempermin=tempermin,
        )

        return new

    @classmethod
    def modes(cls, structure, nqsmall, q1shft=(0, 0, 0), enunit=2, asr=1, chneut=1, **kwargs):
        """
        Build an anaddb input file for the computation of phonon modes.

        Args:
            Structure:
                Structure object
            ngqpt:
                Monkhorst-Pack divisions for the phonon Q-mesh (coarse one)
            nqsmall:
                Used to generate the (dense) mesh for the DOS.
                It defines the number of q-points used to sample the smallest lattice vector.
            q1shft:
                Shifts used for the coarse Q-mesh
            qptbounds
                Boundaries of the path. If None, the path is generated from an internal database
                depending on the input structure.
            asr, chneut, dipdp:
                Anaddb input variable. See official documentation.

        #!General information
         #enunit    2
         #eivec     1
        #!Flags
         #dieflag   1
         #ifcflag   1
         #ngqpt     1 1 1
        #!Effective charges
         #asr       2
         #chneut    2
        #!Wavevector list number 1
         #nph1l     1
         #qph1l   0.0  0.0  0.0    1.0   ! (Gamma point)
        #!Wavevector list number 2
         #nph2l     3      ! number of phonons in list 1
         #qph2l   1.0  0.0  0.0    0.0
         #        0.0  1.0  0.0    0.0
         #        0.0  0.0  1.0    0.0
        """

        new = cls(structure, comment="ANADB input for modes", **kwargs)

        new.set_variables(
            enunit=enunit,
            eivec=1,
            ifcflag=1,
            dieflag=1,
            ngqpt=[1.0, 1.0, 1.0],
            asr=asr,
            chneut=chneut,
            nph1l=1,
            qph1l=[0.0, 0.0, 0.0, 1.0],
            nph2l=3,
            qph2l=[[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]
        )

        return new



    @property
    def structure(self):
        return self._structure

    def __repr__(self):
        return "<%s at %s>" % (self.__class__.__name__, id(self))

    def __str__(self):
        return self.to_string()

    def make_input(self):
        return self.to_string()

    def to_string(self, sortmode=None):
        """
        String representation.

        Args:
            sortmode:
                "a" for alphabetical order, None if no sorting is wanterd
        """
        lines = []
        app = lines.append

        if self.comment:
            app("# " + self.comment.replace("\n", "\n#"))

        if sortmode is None:
            # no sorting.
            keys = self.keys()
        elif sortmode == "a":
            # alphabetical order.
            keys = sorted(self.keys())
        else:
            raise ValueError("Unsupported value for sortmode %s" % str(sortmode))

        for varname in keys:
            value = self[varname]
            app(str(InputVariable(varname, value)))

        return "\n".join(lines)

    def set_variable(self, varname, value):
        """Set a single variable."""
        if varname in self:
            try:
                iseq = (self[varname] == value)
                iseq = np.all(iseq)
            except ValueError:
                # array like.
                iseq = np.allclose(self[varname], value)
            else:
                iseq = False

            if not iseq:
                msg = "%s is already defined with a different value:\nOLD:\n %s,\nNEW\n %s" % (
                    varname, str(self[varname]), str(value))
                warnings.warn(msg)

        if not self.is_anaddb_var(varname):
            raise self.Error("%s is not a valid ANADDB variable." % varname)

        self[varname] = value

    def set_variables(self, **vars):
        """Set the value of the variables provided in the dictionary vars"""
        for varname, varvalue in vars.items():
            self.set_variable(varname, varvalue)

    def add_extra_abivars(self, abivars):
        """
        This method is needed not to break the API used for strategies

        Connection is explicit via the input file
        since we can pass the paths of the output files
        produced by the previous runs.
        """

    def is_anaddb_var(self, varname):
        """"True if varname is a valid anaddb variable."""
        return is_anaddb_var(varname)

    def set_qpath(self, ndivsm, qptbounds=None):
        """
        Set the variables for the computation of the phonon band structure.

        Args:
            ndivsm:
                Number of divisions for the smallest segment.
            qptbounds
                q-points defining the path in k-space.
                If None, we use the default high-symmetry k-path
                defined in the pymatgen database.
        """
        if qptbounds is None:
            qptbounds = self.structure.calc_kptbounds()

        qptbounds = np.reshape(qptbounds, (-1, 3))

        # TODO: Add support for ndivsm in anaddb
        self.set_variables(
            #ndivsm=ndivsm,
            nqpath=len(qptbounds),
            qpath=qptbounds,
        )

    def set_autoqmesh(self, nqsmall):
        """
        Set the variabls nqpt for the sampling of the BZ.

        Args:
            nqsmall:
                Number of q-points used to sample the smallest lattice vector.
        """
        self.set_variables(ng2qpt=self.structure.calc_ngkpt(nqsmall))
