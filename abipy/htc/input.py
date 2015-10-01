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
import json
import numpy as np

from collections import OrderedDict
from monty.dev import deprecated
from monty.collections import dict2namedtuple
from monty.string import is_string, list_strings
from monty.os.path import which
from pymatgen.core.units import Energy
from pymatgen.serializers.json_coders import PMGSONable, pmg_serialize
from pymatgen.io.abinit.pseudos import PseudoTable, Pseudo
from pymatgen.io.abinit.tasks import TaskManager, AbinitTask
from pymatgen.io.abinit.netcdf import NetcdfReader
from pymatgen.io.abinit.abiinspect import yaml_read_irred_perts
from abipy.core.structure import Structure
from abipy.core.mixins import Has_Structure
from .variable import InputVariable
from abipy.abio.abivars import is_abivar, is_anaddb_var
from abipy.abio.abivars_db import get_abinit_variables

import logging
logger = logging.getLogger(__file__)


__all__ = [
    "AbiInput",
    "LdauParams",
    "LexxParams",
    "input_gen",
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


def _idt_varname(varname):
    """
    Find the dataset index from the name of the variable.
    Return the dataset index and the name of the variable with the numeric index removed.
    dataset index is set to 0 for global variable.

    >>> assert _idt_varname("foo") == (0, "foo")
    >>> assert _idt_varname("foo1") == (1, "foo")
    """
    if not varname[-1].isdigit(): 
        return 0, varname

    # Find the index of the dataset.
    idt = ""
    for i, c in enumerate(reversed(varname)):
        if c.isdigit():
            idt += c
        else:
            break

    # Strip the numeric index in varname
    return int("".join(reversed(idt))), varname[:len(varname)-i]


@six.add_metaclass(abc.ABCMeta)
class Input(six.with_metaclass(abc.ABCMeta, PMGSONable, object)):
    """
    Base class for Input objects.

    An input object must define have a make_input method 
    that returns the string representation used by the client code.
    """
    def deepcopy(self):
        """Deep copy of the input."""
        return copy.deepcopy(self)

    def write(self, filepath):
        """
        Write the input file to file to `filepath`. Returns a string with the input.
        """
        dirname = os.path.dirname(filepath)
        if not os.path.exists(dirname): os.makedirs(dirname)
                                                                                      
        # Write the input file.
        input_string = str(self)
        with open(filepath, "wt") as fh:
            fh.write(input_string)

        return input_string

    #@abc.abstractmethod
    #def make_input(self):

    #@abc.abstractproperty
    #def structure(self):

    #def set_structure(self, structure):


class AbinitInputError(Exception):
    """Base error class for exceptions raised by `AbiInput`"""

#from monty.dev import deprecated
#from abipy.abio.inputs import AbinitInput
#@deprecated(replacement=AbinitInput,
#            message="This class is deprecated and will be removed in abipy0.2.")
class AbiInput(Input, Has_Structure):
    """
    This object represents an ABINIT input file. It supports multi-datasets a
    and provides an easy-to-use interface for defining the variables of the calculation.

    Usage example:

    .. code-block:: python

        inp = AbiInput(pseudos="si.pspnc", pseudo_dir="directory_with_pseudos")
        inp.set_structure("si.cif")
        inp.set_vars(ecut=10, nband=3)

        # Print it to get the final input.
        print(inp)
    """
    Error = AbinitInputError

    def __init__(self, pseudos, pseudo_dir="", structure=None, ndtset=1, comment="", decorators=None):
        """
        Args:
            pseudos: String or list of string with the name of the pseudopotential files.
            pseudo_dir: Name of the directory where the pseudopotential files are located.
            structure: file with the structure, :class:`Structure` object or dictionary with ABINIT geo variable
            ndtset: Number of datasets.
            comment: Optional string with a comment that will be placed at the beginning of the file.
            decorators: List of `AbinitInputDecorator` objects.
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

            self._pseudos = PseudoTable(pseudo_paths)

        if structure is not None: self.set_structure(structure)
        if comment is not None: self.set_comment(comment)

        self._decorators = [] if not decorators else decorators

    def __str__(self):
        """String representation i.e. the input file read by Abinit."""
        if self.ndtset > 1:
            s = ""
            for i, dataset in enumerate(self):
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

        # Add JSON section with pseudo potentials.
        ppinfo = ["\n\n\n#<JSON>"]
        d = {"pseudos": [p.as_dict() for p in self.pseudos]}
        ppinfo.extend(json.dumps(d, indent=4).splitlines())
        ppinfo.append("</JSON>")

        return s + "\n#".join(ppinfo)

    def __getitem__(self, key):
        return self._datasets[key]

    def __getattr__(self, varname):
        try:
            return super(AbiInput, self).__geattr__(self, varname)
        except AttributeError:
            try:
                idt, varname = _idt_varname(varname)
                return self._datasets[idt][varname]
            except Exception as exc:
                raise AttributeError(str(exc))

    def __setattr__(self, varname, value):
        if varname.startswith("ndtset"):
            raise self.Error("Cannot change read-only variable ndtset")

        if varname.startswith("_"):
            # Standard behaviour for "private" variables.
            super(AbiInput, self).__setattr__(varname, value)
        else:
            idt, varname = _idt_varname(varname)

            if idt == 0:
                # Global variable.
                self[0].set_var(varname, value)
            else:
                if not (self.ndtset >= idt >= 1):
                    raise self.Error("Invalid dataset index %d, ndtset %d " % (idt, self.ndtset))

                self[idt].set_var(varname, value)

                # Handle no_multi variables.
                if varname in _ABINIT_NO_MULTI:

                    if varname in self[0]:
                        glob_value = np.array(self[0][varname])
                        isok = np.all(glob_value == np.array(value))
                        if not isok:
                            err_msg = "NO_MULTI variable: dataset 0: %s, dataset %d: %s" % (str(glob_value), idt, str(value))
                            raise self.Error(err_msg)
                    else:
                        self[0].set_var(varname, value)

    @property
    def ndtset(self):
        """Number of datasets."""
        return self[0]["ndtset"]

    @property
    def pseudos(self):
        """List of :class:`Pseudo` objects."""
        return self._pseudos

    @property
    def ispaw(self):
        """True if PAW calculation."""
        return all(p.ispaw for p in self.pseudos)

    @property
    def isnc(self):
        """True if norm-conserving calculation."""
        return all(p.isnc for p in self.pseudos)

    @property
    def num_valence_electrons(self):
        """Number of valence electrons computed from the pseudos and the structure."""
        return self.structure.num_valence_electrons(self.pseudos)

    @property
    def structure(self):
        """
        Returns the :class:`Structure` associated to the ABINIT input.

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

    def new_with_vars(self, *args, **kwargs):
        """Generate a new input (deep copy) with these variables"""
        new = self.deepcopy()
        new.set_vars(*args, **kwargs)
        return new

    @property
    def decorators(self):
        return self._decorators

    def register_decorator(self, decorator):
        """Register a :class:`AbinitInputDecorator`."""
        self._decorators.append(decorator.as_dict())

    def generate(self, **kwargs):
        """
        Generate new inputs by replacing the variables specified in kwargs.

        .. code-block:: python

            # To generate two input files with different values of ecut:
            for inp_ecut in gs_inp.generate(ecut=[10, 20])):
                print("do something with inp_ecut %s" % inp_ecut)

            # To generate four input files with all the possible combinations of ecut and nsppol:
            for inp_ecut in gs_inp.generate(ecut=[10, 20], nsppol=[1, 2]):
                print("do something with inp_ecut %s" % inp_ecut)
        """
        if self.ndtset != 1: raise ValueError("Cannot use generate methods when ndtset != 1")
        return input_gen(self, **kwargs)

    def _dtset2range(self, dtset):
        """Helper function to convert dtset into a range. Accepts int, slice objects or iterable."""
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
        news = []
        for i in range(self.ndtset):
            my_vars = self[i+1].allvars
            my_vars.pop("ndtset", None)

            # Cannot use get*, ird* variables since links must be explicit.
            for varname in my_vars:
                if varname.startswith("get") or varname.startswith("ird"):
                    err_msg = "get* or ird variables should not be present in the input when you split it into datasets"
                    raise self.Error(err_msg)

            new = self.__class__(pseudos=self.pseudos, ndtset=1, decorators=self._decorators)
            new.set_vars(**my_vars)
            new.set_mnemonics(self.mnemonics)
            news.append(new)
    
        return news

    def set_vars(self, *args, **kwargs):
        """
        Set the value of a set of input variables.

        Args:
            dtset: Int with the index of the dataset, slice object of iterable 
            kwargs: Dictionary with the variables.

        Returns:
            self
        """
        dtset = kwargs.pop("dtset", 0)
        kwargs.update(dict(*args))
        for idt in self._dtset2range(dtset):
            self[idt].set_vars(**kwargs)

        return self

    # Alias
    set_variables = set_vars
    set_variables = deprecated(replacement=set_vars)(set_variables)

    def remove_vars(self, keys, dtset=0):
        """
        Remove the variables listed in keys
                                                                             
        Args:
            dtset: Int with the index of the dataset, slice object of iterable 
            keys: List of variables to remove.
        """
        for idt in self._dtset2range(dtset):
            self[idt].remove_vars(keys)

    # Alias
    remove_variables = remove_vars
    remove_variables = deprecated(replacement=remove_vars)(remove_variables)

    def set_comment(self, comment, dtset=0):
        """Set the main comment in the input file."""
        for idt in self._dtset2range(dtset):
            self[idt].set_comment(comment)

    def set_mnemonics(self, boolean):
        """True if mnemonics should be printed"""
        self._mnemonics = bool(boolean)
        for idt in range(self.ndtset):
            self[idt].set_mnemonics(boolean)

    @property
    def mnemonics(self):
        """Return True if mnemonics should be printed"""
        try:
            return self._mnemonics
        except AttributeError:
            return False
        
    def get_ibz(self, ngkpt=None, shiftk=None, kptopt=None, qpoint=None, workdir=None, manager=None):
        """
        This function, computes the list of points in the IBZ and the corresponding weights.
        It should be called with an input file that contains all the mandatory variables required by ABINIT.

        Args:
            ngkpt: Number of divisions for the k-mesh (default None i.e. use ngkpt from self)
            shiftk: Shiftks (default None i.e. use shiftk from self)
            qpoint: qpoint in reduced coordinates. Used to shift the k-mesh (default None i.e no shift)
            workdir: Working directory of the fake task used to compute the ibz. Use None for temporary dir.
            manager: :class:`TaskManager` of the task. If None, the manager is initialized from the config file.

        Returns:
            `namedtuple` with attributes:
                points: `ndarray` with points in the IBZ in reduced coordinates.
                weights: `ndarray` with weights of the points.

        .. warning::

            Multiple datasets are ignored. Only the list of k-points for dataset 1 are returned.
        """
        if self.ndtset != 1:
            raise RuntimeError("get_ibz cannot be used if the input contains more than one dataset")

        # Avoid modifications in self.
        inp = self.split_datasets()[0].deepcopy()

        # The magic value that makes ABINIT print the ibz and then stop.
        inp.prtkpt = -2

        if ngkpt is not None: inp.ngkpt = ngkpt
        if shiftk is not None:
            inp.shiftk = np.reshape(shiftk, (-1,3))
            inp.nshiftk = len(inp.shiftk)

        if kptopt is not None:
            inp.kptopt = kptopt

        if qpoint is not None:
            inp.qptn, inp.nqpt = qpoint, 1

        # Build a Task to run Abinit in a shell subprocess
        task = AbinitTask.temp_shell_task(inp, workdir=workdir, manager=manager)
        task.start_and_wait(autoparal=False)

        # Read the list of k-points from the netcdf file.
        try:
            with NetcdfReader(os.path.join(task.workdir, "kpts.nc")) as r:
                ibz = collections.namedtuple("ibz", "points weights")
                return ibz(points=r.read_value("reduced_coordinates_of_kpoints"),
                           weights=r.read_value("kpoint_weights"))

        except Exception as exc:
            # Try to understand if it's a problem with the Abinit input.
            report = task.get_event_report()
            if report.errors: raise self.Error(str(report))
            raise exc

    def get_irred_perts(self, ngkpt=None, shiftk=None, kptopt=None, qpoint=None, workdir=None, manager=None):
        """
        This function, computes the list of irreducible perturbations for DFPT.
        It should be called with an input file that contains all the mandatory variables required by ABINIT.

        Args:
            ngkpt: Number of divisions for the k-mesh (default None i.e. use ngkpt from self)
            shiftk: Shiftks (default None i.e. use shiftk from self)
            qpoint: qpoint in reduced coordinates. Used to shift the k-mesh (default None i.e no shift)
            workdir: Working directory of the fake task used to compute the ibz. Use None for temporary dir.
            manager: :class:`TaskManager` of the task. If None, the manager is initialized from the config file.

        Returns:
            List of dictionaries with the Abinit variables defining the irreducible perturbation
            Example:

                [{'idir': 1, 'ipert': 1, 'qpt': [0.25, 0.0, 0.0]},
                 {'idir': 2, 'ipert': 1, 'qpt': [0.25, 0.0, 0.0]}]

        .. warning::

            Multiple datasets are ignored. Only the list of k-points for dataset 1 are returned.
        """
        if self.ndtset != 1:
            raise RuntimeError("get_irred_perts cannot be used if the input contains more than one dataset")

        warnings.warn("get_irred_perts is still under development.")
        # Avoid modifications in self.
        inp = self.split_datasets()[0].deepcopy()

        # Use the magic value paral_rf = -1 to get the list of irreducible perturbations for this q-point.
        d = dict(
            paral_rf=-1,
            rfatpol=[1, len(inp.structure)],  # Set of atoms to displace.
            rfdir=[1, 1, 1],                  # Along this set of reduced coordinate axis.
        )
        inp.set_vars(d)

        # Build a Task to run Abinit in a shell subprocess
        task = AbinitTask.temp_shell_task(inp, workdir=workdir, manager=manager)
        task.start_and_wait(autoparal=False)

        # Parse the file to get the perturbations.
        try:
            return yaml_read_irred_perts(task.log_file.path)
        except Exception as exc:
            # Try to understand if it's a problem with the Abinit input.
            report = task.get_event_report()
            if report.errors: raise self.Error(str(report))
            raise exc

    def linspace(self, varname, start, stop, endpoint=True):
        """
        Returns `ndtset` evenly spaced samples, calculated over the interval [`start`, `stop`].

        The endpoint of the interval can optionally be excluded.

        Args:
            start: The starting value of the sequence.
            stop: The end value of the sequence, unless `endpoint` is set to False.
                In that case, the sequence consists of all but the last of ``ndtset + 1``
                evenly spaced samples, so that `stop` is excluded.  Note that the step
                size changes when `endpoint` is False.
        """
        if varname in self[0]:
            raise self.Error("varname %s is already defined in globals" % varname)
            
        lspace = np.linspace(start, stop, num=self.ndtset, endpoint=endpoint, retstep=False)

        for dataset, value in zip(self[1:], lspace):
            dataset.set_var(varname, value)

    def arange(self, varname, start, stop, step):
        """
        Return evenly spaced values within a given interval.

        Values are generated within the half-open interval ``[start, stop)``
        (in other words, the interval including `start` but excluding `stop`).

        When using a non-integer step, such as 0.1, the results will often not
        be consistent.  It is better to use ``linspace`` for these cases.

        Args:
            start:  Start of interval. The interval includes this value. The default start value is 0.
            stop: End of interval.  The interval does not include this value, except
                in some cases where `step` is not an integer and floating point
            step: Spacing between values.  For any output `out`, this is the distance
                between two adjacent values, ``out[i+1] - out[i]``.  The default
                step size is 1.  If `step` is specified, `start` must also be given.
        """
        if varname in self[0]:
            raise self.Error("varname %s is already defined in globals" % varname)

        arang = np.arange(start=start, stop=stop, step=step)
        if len(arang) != self.ndtset:
            raise self.Error("len(arang) %d != ndtset %d" % (len(arang), self.ndtset))

        for (dataset, value) in zip(self[1:], arang):
            dataset.set_var(varname, value)

    def product(self, *items):
        """
        Cartesian product of input iterables. Equivalent to nested for-loops.

        .. code-block:: python

            inp.product("tsmear", "ngkpt", [[2,2,2], [4,4,4]], [0.1, 0.2, 0.3])
        """
        # Split items into varnames and values
        for i, item in enumerate(items):
            if not is_string(item): break

        varnames, values = items[:i], items[i:]
        if len(varnames) != len(values):
            raise self.Error("The number of variables must equal the number of lists")

        varnames = [ [varnames[i]] * len(values[i]) for i in range(len(values))]
        varnames = itertools.product(*varnames)
        values = itertools.product(*values)

        idt = 0
        for names, values in zip(varnames, values):
            idt += 1
            self[idt].set_vars(**{k: v for k, v in zip(names, values)})
        
        if idt != self.ndtset:
            raise self.Error("The number of configurations must equal ndtset while %d != %d" % (idt, self.ndtset))

    def set_structure(self, structure, dtset=0):
        """Set the :class:`Structure` object for the specified dtset."""
        structure = Structure.as_structure(structure)

        if dtset is None:
            dtset = slice(self.ndtset+1)
        for idt in self._dtset2range(dtset):
            self[idt].set_structure(structure)

    @deprecated(set_structure)
    def set_structure_from_file(self, filepath, dtset=0, cls=Structure):
        """
        Set the :class:`Structure` object for the specified dtset. 
        data is read from filepath. cls specifies the class to instantiate.
        """
        structure = cls.from_file(filepath)
        self.set_structure(structure, dtset=dtset)
        return structure

    def set_kmesh(self, ngkpt, shiftk, kptopt=1, dtset=0):
        """Set the variables defining the k-point sampling for the specified dtset."""
        for idt in self._dtset2range(dtset):
            self[idt].set_kmesh(ngkpt, shiftk, kptopt=kptopt)

    def set_autokmesh(self, nksmall, kptopt=1, dtset=0):
        """
        Set the variables (ngkpt, shift, kptopt) for the sampling of the BZ.
                                                       
        Args:
            nksmall: Number of k-points used to sample the smallest lattice vector.
            kptopt: Option for the generation of the mesh (ABINIT variable).
        """
        for idt in self._dtset2range(dtset):
            self[idt].set_autokmesh(nksmall, kptopt=kptopt)

    def set_autokpath(self, ndivsm, dtset=0):
        """
        Set automatically the k-path from the lattice 
        and the number of divisions for the smallest segment (ndivsm)
        """
        for idt in self._dtset2range(dtset):
            self[idt].set_kpath(ndivsm, kptbounds=None)

    def set_kpath(self, ndivsm, kptbounds=None, iscf=-2, dtset=0):
        """
        Set the variables defining the k-path for the specified dtset.
        The list of K-points is taken from the pymatgen database if `kptbounds` is None.

        Args:
            ndivsm: Number of divisions for the smallest segment.
            kptbounds: k-points defining the path in k-space.
            iscf: iscf variable.
            dtset: Index of the dataset. 0 for global variables.
        """
        for idt in self._dtset2range(dtset):
            self[idt].set_kpath(ndivsm, kptbounds=kptbounds, iscf=iscf)

    def set_kptgw(self, kptgw, bdgw, dtset=0):
        """Set the variables defining the k point for the GW corrections"""
        for idt in self._dtset2range(dtset):
            self[idt].set_kptgw(kptgw, bdgw)

    @pmg_serialize
    def as_dict(self):
        dtsets = []
        for ds in self:
            ds_copy = ds.deepcopy()
            for key, value in ds_copy.items():
                if isinstance(value, np.ndarray):
                    ds_copy[key] = value.tolist()
            dtsets.append(dict(ds_copy))

        return {'pseudos': [p.as_dict() for p in self.pseudos], 
                'datasets': dtsets,
                "decorators": [dec.as_dict() for dec in self._decorators],
                }

    @classmethod
    def from_dict(cls, d):
        pseudos = []
        for p in d['pseudos']:
            pseudos.append(Pseudo.from_file(p['filepath']))

        dtsets = d['datasets']
        abiinput = cls(pseudos, ndtset=dtsets[0]['ndtset'], decorators=d["decorators"])

        for n, ds in enumerate(dtsets):
            abiinput.set_vars(dtset=n, **ds)

        return abiinput

    def new_from_decorators(self, decorators):
        """
        This function receives a list of :class:`AbinitInputDecorator` objects or just a single object,
        applyes the decorators to the input and returns a new :class:`AbiInput` object.
        self is not changed.
        """
        if not isinstance(decorators, (list, tuple)): decorators = [decorators]

        # Deepcopy only at the first step to improve performance.
        inp = self
        for i, dec in enumerate(decorators):
            inp = dec(inp, deepcopy=(i == 0))

        return inp

    def validate(self):
        """
        Run ABINIT in dry mode to validate the input file.

        Return:
            `namedtuple` with the following attributes:

                retcode: Return code. 0 if OK.
                log_file:  log file of the Abinit run, use log_file.read() to access its content.
                stderr_file: stderr file of the Abinit run. use stderr_file.read() to access its content.

        Raises:
            `RuntimeError` if executable is not in $PATH.
        """
        task = AbinitTask.temp_shell_task(inp=self) 
        retcode = task.start_and_wait(autoparal=False, exec_args=["--dry-run"])
        return dict2namedtuple(retcode=retcode, log_file=task.log_file, stderr_file=task.stderr_file)



class MappingMixin(collections.Mapping):
    """
    Mixin class implementing the mapping protocol. Useful to avoid boilerplate code if you want
    to define a object that behaves as a Mapping but without inheriting from dict or similar classes
    because you don't want to expose/support all the methods of dict.

    Client code must initialize a Mapping object either in new or in init and bound it to _mapping_mixin_
    The implemention of the Mapping interface is delegated to _mapping_mixin_

    .. Example:

    >>> class Foo(MappingMixin):
    ...     def __init__(self, attr, **kwargs):
    ...         self._mapping_mixin_ = kwargs
    ...         self.attr = attr
    >>> obj = Foo(attr=1, spam=2)
    >>> obj.attr, obj["spam"]
    (1, 2)
    >>> obj.pop("spam")
    2
    >>> len(obj), "spam" in obj
    (0, False)
    """
    def __len__(self):
        return len(self._mapping_mixin_)

    def __iter__(self):
        return self._mapping_mixin_.__iter__()

    def __getitem__(self, key):
        return self._mapping_mixin_[key]

    def __setitem__(self, key, value):
        self._mapping_mixin_[key] = value

    def __contains__(self, key):
        return key in self._mapping_mixin_

    def keys(self):
        return self._mapping_mixin_.keys()

    def items(self):
        return self._mapping_mixin_.items()

    def get(self, value, default=None):
        try:
            return self[value]
        except KeyError:
            return default

    def pop(self, k, *d):
        """
        D.pop(k[,d]) -> v, remove specified key and return the corresponding value.
        If key is not found, d is returned if given, otherwise KeyError is raised
        """
        if d:
            return self._mapping_mixin_.pop(k, d[0])
        else:
            return self._mapping_mixin_.pop(k)

    def update(self, e, **f):
        """
        D.update([E, ]**F) -> None.  Update D from dict/iterable E and F.
        If E present and has a .keys() method, does:     for k in E: D[k] = E[k]
        If E present and lacks .keys() method, does:     for (k, v) in E: D[k] = v
        In either case, this is followed by: for k in F: D[k] = F[k]
        """
        self._mapping_mixin_.update(e, **f)



class Dataset(MappingMixin, Has_Structure):
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
        return "<%s at %s>" % (self.__class__.__name__, id(self))

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

    #@abc.property
    #def runlevel(self):
    #    """String defining the Runlevel. See _runl2optdriver."""
    # Mapping runlevel --> optdriver variable
    #_runl2optdriver = {
    #    "scf": 0,
    #    "nscf": 0,
    #    "relax": 0,
    #    "dfpt": 1,
    #    "screening": 3,
    #    "sigma": 4,
    #    "bse": 99,
    #}
    #    # Find the value of optdriver (firt in self, then in globals finally use default value.
    #    optdriver = self.get("optdriver")
    #    if optdriver is None: optdriver = self.dt0.get("optdriver")
    #    if optdriver is None: optdriver = 0

    #    # At this point we have to understand the type of calculation.

    def set_mnemonics(self, boolean):
        """True if mnemonics should be printed"""
        self._mnemonics = bool(boolean)

    @property
    def mnemonics(self):
        """Return True if mnemonics should be printed"""
        try:
            return self._mnemonics
        except AttributeError:
            return False

    def to_string(self, sortmode=None, post=None):
        """
        String representation.

        Args:
            sortmode: "a" for alphabetical order, None if no sorting is wanted
            post: String that will be appended to the name of the variables
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

        with_mnemonics = self.mnemonics
        if with_mnemonics:
            var_database = get_abinit_variables()

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

            if with_mnemonics:
                v = var_database[var]
                app("# <" + v.definition + ">")

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

    def set_var(self, varname, value):
        """Set a the value of a variable."""
        # Check if varname is in the internal database.
        if not is_abivar(varname):
            raise self.Error("varname %s is not a valid ABINIT variable\n."
                             "Modify abipy/htc/abinit_vars.json" % varname)

        # Debugging section.
        if False and varname in self:
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
                logger.debug(msg)

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
                self.dt0.set_var(varname, value)

    # Alias
    set_variable = set_var
    set_variable = deprecated(replacement=set_var)(set_variable)

    def set_vars(self, *args, **kwargs):
        """Set the value of the variables provied in the dictionary **kwargs"""
        kwargs.update(dict(*args))
        for varname, varvalue in kwargs.items():
            self.set_var(varname, varvalue)

    # Alias
    set_variables = set_vars
    set_variables = deprecated(replacement=set_vars)(set_variables)

    def remove_vars(self, keys):
        """Remove the variables listed in keys."""
        for key in list_strings(keys):
            if key not in self:
                raise KeyError("key: %s not in self:\n %s" % (key, list(self.keys())))
            self.pop(key)

    # Alias
    remove_variables = remove_vars
    remove_variables = deprecated(replacement=remove_vars)(remove_variables)

    @property
    def structure(self):
        """Returns the :class:`Structure` associated to this dataset."""
        # TODO: Avoid calling Structure.from_abivars, find a way to cache the object and invalidate it.
        return Structure.from_abivars(self.allvars)
        # Cannot use lazy properties here because we may want to change the structure

        #if hasattr(self, "_structure") and self._structure is None:
        #    return self._structure
        #try:
        #    return self._structure
        #except AttributeError:
        #    structure = Structure.from_abivars(self.allvars)
        #    self.set_structure(structure)
        #    return self._structure

    def set_structure(self, structure):
        structure = Structure.as_structure(structure)

        self._structure = structure
        if structure is None: return

        self.set_vars(**structure.to_abivars())

    # Helper functions to facilitate the specification of several variables.
    def set_kmesh(self, ngkpt, shiftk, kptopt=1):
        """
        Set the variables for the sampling of the BZ.

        Args:
            ngkpt: Monkhorst-Pack divisions
            shiftk: List of shifts.
            kptopt: Option for the generation of the mesh.
        """
        shiftk = np.reshape(shiftk, (-1,3))
        
        self.set_vars(ngkpt=ngkpt, kptopt=kptopt, nshiftk=len(shiftk), shiftk=shiftk)

    def set_autokmesh(self, nksmall, kptopt=1):
        """
        Set the variables (ngkpt, shift, kptopt) for the sampling of the BZ.
                                                       
        Args:
            nksmall: Number of k-points used to sample the smallest lattice vector.
            kptopt: Option for the generation of the mesh.
        """
        shiftk = self.structure.calc_shiftk()
        
        self.set_vars(ngkpt=self.structure.calc_ngkpt(nksmall), kptopt=kptopt, 
                      nshiftk=len(shiftk), shiftk=shiftk)

    def set_kpath(self, ndivsm, kptbounds=None, iscf=-2):
        """
        Set the variables for the computation of the band structure.

        Args:
            ndivsm: Number of divisions for the smallest segment.
            kptbounds: k-points defining the path in k-space.
                If None, we use the default high-symmetry k-path defined in the pymatgen database.
        """
        if kptbounds is None: kptbounds = self.structure.calc_kptbounds()
        kptbounds = np.reshape(kptbounds, (-1,3))

        self.set_vars(kptbounds=kptbounds, kptopt=-(len(kptbounds)-1), ndivsm=ndivsm, iscf=iscf)

    def set_kptgw(self, kptgw, bdgw):
        """
        Set the variables (k-points, bands) for the computation of the GW corrections.

        Args
            kptgw: List of k-points in reduced coordinates.
            bdgw: Specifies the range of bands for the GW corrections.
                Accepts iterable that be reshaped to (nkptgw, 2) 
                or a tuple of two integers if the extrema are the same for each k-point.
        """
        kptgw = np.reshape(kptgw, (-1,3))
        nkptgw = len(kptgw)
        if len(bdgw) == 2: bdgw = len(kptgw) * bdgw

        self.set_vars(kptgw=kptgw, nkptgw=nkptgw, bdgw=np.reshape(bdgw, (nkptgw, 2)))


class LujForSpecie(collections.namedtuple("LdauForSpecie", "l u j unit")):
    """
    This object stores the value of l, u, j used for a single atomic specie.
    """
    def __new__(cls, l, u, j, unit):
        """
        Args:
            l: Angular momentum (int or string).
            u: U value
            j: J Value
            unit: Energy unit for u and j.
        """
        l = l
        u, j = Energy(u, unit), Energy(j, unit)
        return super(cls, LujForSpecie).__new__(cls, l, u, j, unit)


class LdauParams(object):
    """
    This object stores the parameters for LDA+U calculations with the PAW method
    It facilitates the specification of the U-J parameters in the Abinit input file.
    (see `to_abivars`). The U-J operator will be applied only on the atomic species 
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
            usepawu: ABINIT variable `usepawu` defining the LDA+U method.
            structure: :class:`Structure` object.
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
            symbol: Chemical symbol of the atoms on which LDA+U should be applied.
            l: Angular momentum.
            u: Value of U.
            j: Value of J.
            unit: Energy unit of U and J.
        """
        if symbol not in self.symbols_by_typat:
            raise ValueError("Symbol %s not in symbols_by_typat:\n%s" % (symbol, self.symbols_by_typat))

        if symbol in self._params:
            raise ValueError("Symbol %s is already present in LdauParams! Cannot overwrite:\n" % symbol)

        self._params[symbol] = LujForSpecie(l=l, u=u, j=j, unit=unit)

    def to_abivars(self):
        """Returns a dict with the Abinit variables."""
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
            jpawu=" ".join(map(str, jpawu)) + " eV")


class LexxParams(object):
    """
    This object stores the parameters for local exact exchange calculations with the PAW method
    It facilitates the specification of the LEXX parameters in the Abinit input file.
    (see `to_abivars`). The LEXX operator will be applied only on the atomic species 
    that have been selected by calling `lexx_for_symbol`.

    To perform a LEXX calculation for NiO in which the LEXX is computed only for the l=2
    channel of the nickel atoms:
                                                                         
    .. code-block:: python

        lexx_params = LexxParams(nio_structure)
        lexx_params.lexx_for_symbol("Ni", l=2)

        print(lexc_params.to_abivars())
    """
    def __init__(self, structure):
        """
        Arg:
            structure: :class:`Structure` object.
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
            symbol: Chemical symbol of the atoms on which LEXX should be applied.
            l: Angular momentum.
        """
        if symbol not in self.symbols_by_typat:
            err_msg = "Symbol %s not in symbols_by_typat:\n%s" % (symbol, self.symbols_by_typat)
            raise ValueError(err_msg)

        if symbol in self._lexx_for_symbol:
            raise ValueError("Symbol %s is already present in LdauParams! Cannot overwrite:" % symbol)

        self._lexx_for_symbol[symbol] = l

    def to_abivars(self):
        """Returns a dict with the Abinit variables."""
        lexx_typat = []

        for symbol in self.symbols_by_typat:
            l = self._lexx_for_symbol.get(symbol, -1)
            lexx_typat.append(int(l)) 

        return dict(useexexch=1, lexexch=" ".join(map(str, lexx_typat)))


def input_gen(inp, **kwargs):
    """
    This function receives an :class:`AbiInput` and generates
    new inputs by replacing the variables specified in kwargs.

    Args:
        inp: :class:`AbiInput` file.
        kwargs: keyword arguments with the values used for each variable.

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
        # Remove the variable names to avoid annoying warnings if the variable is overwritten.
        new_inp.remove_vars(new_vars.keys())
        new_inp.set_vars(**new_vars)

        yield new_inp


def product_dict(d):
    """
    This function receives a dictionary d where each key defines a list of items or a simple scalar.
    It constructs the Cartesian product of the values (equivalent to nested for-loops),
    and returns a list of dictionaries with the values that would be used inside the loop.

    >>> d = OrderedDict([("foo", [2, 4]), ("bar", 1)])
    >>> product_dict(d) == [OrderedDict([('foo', 2), ('bar', 1)]), OrderedDict([('foo', 4), ('bar', 1)])]
    True
    >>> d = OrderedDict([("bar", [1,2]), ('foo', [3,4])]) 
    >>> product_dict(d) == [{'bar': 1, 'foo': 3},
    ... {'bar': 1, 'foo': 4},
    ... {'bar': 2, 'foo': 3},
    ... {'bar': 2, 'foo': 4}]
    True

    .. warning:

        Dictionaries are not ordered, therefore one cannot assume that 
        the order of the keys in the output equals the one used to loop.
        If the order is important, one should pass a :class:`OrderedDict` in input.
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
