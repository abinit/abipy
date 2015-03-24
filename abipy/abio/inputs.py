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

from collections import OrderedDict, MutableMapping
from monty.collections import dict2namedtuple
from monty.string import is_string, list_strings
from pymatgen.serializers.json_coders import PMGSONable, pmg_serialize
from pymatgen.io.abinitio.pseudos import PseudoTable, Pseudo
from pymatgen.io.abinitio.tasks import AbinitTask, ParalHintsParser
from pymatgen.io.abinitio.netcdf import NetcdfReader
from pymatgen.io.abinitio.abiinspect import yaml_read_irred_perts
from abipy.core.structure import Structure
from abipy.core.mixins import Has_Structure
from abipy.htc.variable import InputVariable
from abipy.htc.abivars import is_abivar, is_anaddb_var
from abipy.htc.abivars_db import get_abinit_variables

import logging
logger = logging.getLogger(__file__)


class AbinitInputError(Exception):
    """Base error class for exceptions raised by `AbiInput`"""


class AbinitInput(six.with_metaclass(abc.ABCMeta, MutableMapping, PMGSONable, Has_Structure, object)):
    """
    This object stores the ABINIT variables for a single dataset.
    """
    Error = AbinitInputError

    def __init__(self, pseudos, pseudo_dir=None, structure=None, comment=None, decorators=None, abivars=None):
        """
        Args:
            pseudos: String or list of strings with the name of the pseudopotential files.
            pseudo_dir: Name of the directory where the pseudopotential files are located.
            structure: file with the structure, :class:`Structure` object or dictionary with ABINIT geo variable.
            ndtset: Number of datasets.
            comment: Optional string with a comment that will be placed at the beginning of the file.
            decorators: List of `AbinitInputDecorator` objects.
            abivars: Dictionary with the initial set of variables. Default: Empty
        """
        if pseudo_dir is not None:
            pseudo_dir = os.path.abspath(pseudo_dir)
            if not os.path.exists(pseudo_dir): raise self.Error("Directory  %s does not exist")
            pseudos = [os.path.join(pseudo_dir, p) for p in list_strings(pseudos)]

        self._pseudos = PseudoTable.as_table(pseudos)
        abivars = {} if abivars is None else abivars
        self._abivars = OrderedDict(**abivars)

        if structure is not None: self.set_structure(structure)
        if comment is not None: self.set_comment(comment)

        self._decorators = [] if not decorators else decorators

    @pmg_serialize
    def as_dict(self):
        abivars = {}
        for key, value in self.items():
            if isinstance(value, np.ndarray): value = value.tolist()
            abivars[key] = value

        return dict(pseudos=[p.as_dict() for p in self.pseudos], 
                    comment=self.comment,
                    decorators=[dec.as_dict() for dec in self.decorators],
                    abivars=abivars)

    @classmethod
    def from_dict(cls, d):
        pseudos = [Pseudo.from_file(p['filepath']) for p in d['pseudos']]
        return cls(pseudos, decorators=d["decorators"], comment=d["comment"], abivars=d["abivars"])

    # ABC protocol: __delitem__, __getitem__, __iter__, __len__, __setitem__
    def __delitem__(self, key):
        return self._abivars.__delitem__(key)
        
    def __getitem__(self, key):
        return self._abivars.__getitem__(key)

    def __iter__(self):
        return self._abivars.__iter__()

    def __len__(self):
        return len(self._abivars)

    def __setitem__(self, key, value):
        if not is_abivar(key):
            raise self.Error("%s is not a valid ABINIT variable.\n"
                             "If you are sure the name is correct, please contact the abipy developers\n" 
                             "or modify the JSON file abipy/htc/abinit_vars.json" % key)

        return self._abivars.__setitem__(key, value)

    def __repr__(self):
        return "<%s at %s>" % (self.__class__.__name__, id(self))

    def __str__(self):
        return self.to_string()

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

    def deepcopy(self):
        """Deep copy of the input."""
        return copy.deepcopy(self)

    @property
    def decorators(self):
        return self._decorators

    def register_decorator(self, decorator):
        """Register a :class:`AbinitInputDecorator`."""
        self._decorators.append(decorator.as_dict())

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

    def to_string(self, sortmode=None, post=None, with_pseudos=True):
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
            with_pseudos: False if JSON section with pseudo data should not be added.
        """
        lines = []
        app = lines.append

        if self.comment: app("# " + self.comment.replace("\n", "\n#"))

        post = post if post is not None else ""

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
            #if self.index != 0 and var in _ABINIT_NO_MULTI:
            #    continue

            # Print ndtset only if we really have datasets. 
            # This trick prevents abinit from adding DS# to the output files.
            # thus complicating the creation of workflows made of separated calculations.
            #if var == "ndtset" and value == 1:
            #    continue

            if with_mnemonics:
                v = var_database[var]
                app("# <" + v.definition + ">")

            varname = var + post
            variable = InputVariable(varname, value)
            app(str(variable))

        s = "\n".join(lines)

        if not with_pseudos: return s 

        # Add JSON section with pseudo potentials.
        ppinfo = ["\n\n\n#<JSON>"]
        d = {"pseudos": [p.as_dict() for p in self.pseudos]}
        ppinfo.extend(json.dumps(d, indent=4).splitlines())
        ppinfo.append("</JSON>")
                                                             
        return s + "\n#".join(ppinfo)

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

    @property
    def comment(self):
        try:
            return self._comment
        except AttributeError:
            return None

    def set_comment(self, comment):
        """Set a comment to be included at the top of the file."""
        self._comment = comment

    def set_vars(self, *args, **kwargs):
        """Set the value of the variables provied in the dictionary **kwargs"""
        kwargs.update(dict(*args))
        for varname, varvalue in kwargs.items():
            self[varname] = varvalue
        return kwargs

    def remove_vars(self, keys):
        """Remove the variables listed in keys."""
        values = []
        for key in list_strings(keys):
            if key not in self:
                raise KeyError("key: %s not in self:\n %s" % (key, list(self.keys())))
            values.append(self.pop(key))
        return values

    @property
    def structure(self):
        """Returns the :class:`Structure` associated to this dataset."""
        # TODO: Avoid calling Structure.from_abivars, find a way to cache the object and invalidate it.
        # Copy a structure to avoid references?
        return Structure.from_abivars(self._abivars)
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
        return structure

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
        return self.set_vars(ngkpt=ngkpt, kptopt=kptopt, nshiftk=len(shiftk), shiftk=shiftk)

    def set_autokmesh(self, nksmall, kptopt=1):
        """
        Set the variables (ngkpt, shift, kptopt) for the sampling of the BZ.
                                                       
        Args:
            nksmall: Number of k-points used to sample the smallest lattice vector.
            kptopt: Option for the generation of the mesh.
        """
        shiftk = self.structure.calc_shiftk()
        return self.set_vars(ngkpt=self.structure.calc_ngkpt(nksmall), kptopt=kptopt, 
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

        return self.set_vars(kptbounds=kptbounds, kptopt=-(len(kptbounds)-1), ndivsm=ndivsm, iscf=iscf)

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

        return self.set_vars(kptgw=kptgw, nkptgw=nkptgw, bdgw=np.reshape(bdgw, (nkptgw, 2)))

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

    def linspace(self, varname, start, stop, num=50, endpoint=True):
        """
        Returns `num` evenly spaced samples, calculated over the interval [`start`, `stop`].

        The endpoint of the interval can optionally be excluded.

        Args:
            start: The starting value of the sequence.
            stop: The end value of the sequence, unless `endpoint` is set to False.
                In that case, the sequence consists of all but the last of ``ndtset + 1``
                evenly spaced samples, so that `stop` is excluded.  Note that the step
                size changes when `endpoint` is False.
            num : int, optional
                Number of samples to generate. Default is 50.
            endpoint : bool, optional
                If True, `stop` is the last sample. Otherwise, it is not included.
                Default is True.
        """
        inps = []
        for value in np.linspace(start, stop, num=num, endpoint=endpoint, retstep=False):
            inp = self.deepcopy()
            inp[varname] = value
            inps.append(inp)
        return inps

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
        inps = []
        for value in np.arange(start=start, stop=stop, step=step):
            inp = self.deepcopy()
            inp[varname] = value
            inps.append(inp)
        return inps

    def product(self, *items):
        """
        Cartesian product of input iterables. Equivalent to nested for-loops.

        .. code-block:: python

            inp.product("ngkpt", "tsmear", [[2,2,2], [4,4,4]], [0.1, 0.2, 0.3])
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

        inps = []
        for names, values in zip(varnames, values):
            inp = self.deepcopy()
            inp.set_vars(**{k: v for k, v in zip(names, values)})
            inps.append(inp)
        return inps

    def new_with_decorators(self, decorators):
        """
        This function receives a list of :class:`AbinitInputDecorator` objects or just a single object,
        applyes the decorators to the input and returns a new :class:`AbinitInput` object.
        self is not changed.
        """
        if not isinstance(decorators, (list, tuple)): decorators = [decorators]

        # Deepcopy only at the first step to improve performance.
        inp = self
        for i, dec in enumerate(decorators):
            inp = dec(inp, deepcopy=(i == 0))

        return inp
        
    def abivalidate(self):
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

    def abiget_ibz(self, ngkpt=None, shiftk=None, kptopt=None, workdir=None, manager=None):
        """
        This function, computes the list of points in the IBZ and the corresponding weights.
        It should be called with an input file that contains all the mandatory variables required by ABINIT.

        Args:
            ngkpt: Number of divisions for the k-mesh (default None i.e. use ngkpt from self)
            shiftk: Shiftks (default None i.e. use shiftk from self)
            kptopt: Option for k-point generation. If None, the value in self is used.
            workdir: Working directory of the fake task used to compute the ibz. Use None for temporary dir.
            manager: :class:`TaskManager` of the task. If None, the manager is initialized from the config file.

        Returns:
            `namedtuple` with attributes:
                points: `ndarray` with points in the IBZ in reduced coordinates.
                weights: `ndarray` with weights of the points.

        .. warning::

            Multiple datasets are ignored. Only the list of k-points for dataset 1 are returned.
        """
        # Avoid modifications in self.
        inp = self.deepcopy()

        # The magic value that makes ABINIT print the ibz and then stop.
        inp["prtkpt"] = -2

        if ngkpt is not None: inp["ngkpt"] = ngkpt
        if shiftk is not None:
            shiftk = np.reshape(shiftk, (-1,3))
            inp.set_vars(shiftk=shiftk, nshiftk=len(inp.shiftk))

        if kptopt is not None: inp["kptopt"] = kptopt

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
            if report and report.errors: raise self.Error(str(report))
            raise self.Error("Problem in temp Task executed in %s\n%s" % (task.workdir, exc))

    def abiget_irred_phperts(self, qpt=None, ngkpt=None, shiftk=None, kptopt=None, workdir=None, manager=None):
        """
        This function, computes the list of irreducible perturbations for DFPT.
        It should be called with an input file that contains all the mandatory variables required by ABINIT.

        Args:
            qpt: qpoint of the phonon in reduced coordinates. Used to shift the k-mesh 
                if qpt is not passed, self must already contain "qpt" otherwise an exception is raised.
            ngkpt: Number of divisions for the k-mesh (default None i.e. use ngkpt from self)
            shiftk: Shiftks (default None i.e. use shiftk from self)
            kptopt: Option for k-point generation. If None, the value in self is used.
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
        warnings.warn("get_irred_perts is still under development.")

        # Avoid modifications in self.
        inp = self.deepcopy()

        qpt = inp.get("qpt") if qpt is None else qpt
        if qpt is None:
            raise ValueError("qpt is not in the input and therefore it must be passed explicitly")

        if ngkpt is not None: inp["ngkpt"] = ngkpt
        if shiftk is not None:
            shiftk = np.reshape(shiftk, (-1,3))
            inp.set_vars(shiftk=shiftk, nshiftk=len(inp.shiftk))
                                                                 
        if kptopt is not None: inp["kptopt"] = kptopt

        inp.set_vars(
            rfphon=1,                         # Will consider phonon-type perturbation
            nqpt=1,                           # One wavevector is to be considered
            qpt=qpt,                          # q-wavevector.
            rfatpol=[1, len(inp.structure)],  # Set of atoms to displace.
            rfdir=[1, 1, 1],                  # Along this set of reduced coordinate axis.
            paral_rf=-1,                      # Magic value to get the list of irreducible perturbations for this q-point.
        )

        # Build a Task to run Abinit in a shell subprocess
        task = AbinitTask.temp_shell_task(inp, workdir=workdir, manager=manager)
        task.start_and_wait(autoparal=False)

        # Parse the file to get the perturbations.
        try:
            return yaml_read_irred_perts(task.log_file.path)
        except Exception as exc:
            # Try to understand if it's a problem with the Abinit input.
            report = task.get_event_report()
            if report and report.errors: raise self.Error(str(report))
            raise self.Error("Problem in temp Task executed in %s\n%s" % (task.workdir, exc))

    def get_autoparal_pconfs(self, max_ncpus, autoparal=1, workdir=None, manager=None):
        """Get all the possible configurations up to max_ncpus"""
        inp = self.deepcopy()
        inp.set_vars(autoparal=autoparal, max_ncpus=max_ncpus)

        # Run the job in a shell subprocess with mpi_procs = 1
        # Return code is always != 0 
        task = AbinitTask.temp_shell_task(inp, workdir=workdir, manager=manager)
        task.start_and_wait(autoparal=False)

        ##############################################################
        # Parse the autoparal configurations from the main output file
        ##############################################################
        parser = ParalHintsParser()
        try:
            pconfs = parser.parse(task.output_file.path)
            return pconfs
        except parser.Error:
            # Try to understand if it's a problem with the Abinit input.
            report = task.get_event_report()
            if report and report.errors: raise self.Error(str(report))
            raise self.Error("Problem in temp Task executed in %s\n%s" % (task.workdir, exc))


class MultiDataset(object):

    def __init__(self, pseudos, pseudo_dir="", structure=None, ndtset=1):
        """
        Args:
            pseudos: String or list of string with the name of the pseudopotential files.
            pseudo_dir: Name of the directory where the pseudopotential files are located.
            structure: file with the structure, :class:`Structure` object or dictionary with ABINIT geo variable
            ndtset: Number of datasets.
        """
        # Setup of the pseudopotential files.
        if isinstance(pseudos, PseudoTable):
            pseudos = pseudos

        elif all(isinstance(p, Pseudo) for p in pseudos):
            pseudos = PseudoTable(pseudos)

        else:
            # String(s)
            pseudo_dir = os.path.abspath(pseudo_dir)
            pseudo_paths = [os.path.join(pseudo_dir, p) for p in list_strings(pseudos)]

            missing = [p for p in pseudo_paths if not os.path.exists(p)]
            if missing:
                raise self.Error("Cannot find the following pseudopotential files:\n%s" % str(missing)) 

            pseudos = PseudoTable(pseudo_paths)

        self._inputs = [AbinitInput(pseudos, structure=structure) for i in range(ndtset)]

    @property
    def ndtset(self):
        return len(self)

    def __len__(self):
        return len(self._inputs)

    def __getitem__(self, key):
        return self._inputs[key]

    def __iter__(self):
        return self._inputs.__iter__()

    def __getattr__(self, name):
        #print("in getname with name: %s" % name)
        m = getattr(self._inputs[0], name)
        if m is None:
            raise AttributeError("Cannot find attribute %s in AbinitInput object" % name)
        isattr = not callable(m)

        def on_all(*args, **kwargs):
            results = []
            for obj in self._inputs:
                a = getattr(obj, name)
                #print("name", name, ", type:", type(a), "callable: ",callable(a))
                if callable(a):
                    results.append(a(*args, **kwargs))
                else:
                    results.append(a)

            return results

        if isattr: on_all = on_all()
        return on_all

    def append_input(self, abinit_input):
        self._inputs.append(abinit_input)

    def addnew_from(self, dtindex):
        self.append_input(self[dtindex].deepcopy())

    def split_datasets(self):
        return self._inputs

    @property
    def has_same_structures(self):
        return all(self[0].structure == inp.structure for inp in self)

    def __str__(self):
        """String representation i.e. the input file read by Abinit."""

        if self.ndtset > 1:
            # Multi dataset mode.
            lines = ["ndtset %d" % self.ndtset]

            #same_structures = self.has_same_structures
            for i, inp in enumerate(self):
                header = "### DATASET %d ###" % (i + 1)
                is_last = (i==self.ndtset - 1)

                s = inp.to_string(post=str(i+1), with_pseudos=is_last)
                if s:
                    header = len(header) * "#" + "\n" + header + "\n" + len(header) * "#" + "\n"
                    s = "\n" + header + s + "\n"

                lines.append(s)

            return "\n".join(lines)

        else:
            # single datasets ==> don't append the dataset index to the variables.
            # this trick is needed because Abinit complains if ndtset is not specified 
            # and we have variables that end with the dataset index e.g. acell1
            # We don't want to specify ndtset here since abinit will start to add DS# to 
            # the input and output files thus complicating the algorithms we have to use to locate the files.
            return self[0].to_string()

    #def __dir__(self):
    #    """Interactive prompt"""
    #    #return dir(self) + dir(self._inputs[0])
    #    return dir(self._inputs[0])
