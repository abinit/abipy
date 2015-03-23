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
import abipy.tools.mixins as mixins

from collections import OrderedDict, MutableMapping
from monty.dev import deprecated
from monty.collections import dict2namedtuple
from monty.string import is_string, list_strings
from monty.os.path import which
from pymatgen.core.units import Energy
from pymatgen.serializers.json_coders import PMGSONable, pmg_serialize
from pymatgen.io.abinitio.pseudos import PseudoTable, Pseudo
from pymatgen.io.abinitio.tasks import TaskManager, AbinitTask
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


#class AbinitInput(Input, Has_Structure):
#@six.add_metaclass(abc.ABCMeta)
#class Input(six.with_metaclass(abc.ABCMeta, MutableMapping, PMGSONable, object)):
class AbinitInput(MutableMapping, Has_Structure):
    """
    This object stores the ABINIT variables for a single dataset.
    """
    Error = AbinitInputError

    def __init__(self, pseudos, pseudo_dir=None, structure=None, comment=None, decorators=None):
        """
        Args:
            pseudos: String or list of strings with the name of the pseudopotential files.
            pseudo_dir: Name of the directory where the pseudopotential files are located.
            structure: file with the structure, :class:`Structure` object or dictionary with ABINIT geo variable.
            ndtset: Number of datasets.
            comment: Optional string with a comment that will be placed at the beginning of the file.
            decorators: List of `AbinitInputDecorator` objects.
        """
        if pseudo_dir is not None:
            pseudo_dir = os.path.abspath(pseudo_dir)
            if not os.path.exists(pseudo_dir): raise self.Error("Directory  %s does not exist")
            pseudos = [os.path.join(pseudo_dir, p) for p in list_strings(pseudos)]

        self._pseudos = PseudoTable.as_table(pseudos)
        self._vars = OrderedDict()

        if structure is not None: self.set_structure(structure)
        if comment is not None: self.set_comment(comment)

        self._decorators = [] if not decorators else decorators

    #@classmethod
    #def from_dict(cls, d):

    #@pmg_serialize
    #def as_dict(self):

    # ABC protocol: __delitem__, __getitem__, __iter__, __len__, __setitem__
    def __delitem__(self, key):
        return self._vars.__delitem__(key)
        
    def __getitem__(self, key):
        return self._vars.__getitem__(key)

    def __iter__(self):
        return self._vars.__iter__()

    def __len__(self):
        return len(self._vars)

    def __setitem__(self, key, value):
        if not is_abivar(key):
            raise self.Error("%s is not a valid ABINIT variable\n."
                             "Modify abipy/htc/abinit_vars.json" % key)

        return self._vars.__setitem__(key, value)

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

    def set_vars(self, *args, **kwargs):
        """Set the value of the variables provied in the dictionary **kwargs"""
        kwargs.update(dict(*args))
        for varname, varvalue in kwargs.items():
            self[varname] = varvalue

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

        if self.comment: app("# " + self.comment.replace("\n", "\n#"))

        post = ""
        #if post is None:
        #    post = "" if self.index == 0 else str(self.index)
        #else:
        #    post = post

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

        return "\n".join(lines)

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
        return Structure.from_abivars(self._vars)
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
        # Avoid modifications in self.
        inp = self.deepcopy()

        # The magic value that makes ABINIT print the ibz and then stop.
        inp["prtkpt"] = -2

        # TODO: Redesign API
        #if ngkpt is not None: inp.ngkpt = ngkpt
        #if shiftk is not None:
        #    inp.shiftk = np.reshape(shiftk, (-1,3))
        #    inp.nshiftk = len(inp.shiftk)

        #if kptopt is not None:
        #    inp.kptopt = kptopt

        #if qpoint is not None:
        #    inp.qptn, inp.nqpt = qpoint, 1

        # Build a Task to run Abinit in a shell subprocess
        task = AbinitTask.temp_shell_task(inp, workdir=workdir, manager=manager)
        task.start_and_wait(autoparal=False)
        #print(task.workdir)

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
            raise self.Error("Problem in temp Task executed in %s\n%s" % (task.workdir, exc))

    def get_irred_phperts(self, qpt=None, ngkpt=None, shiftk=None, kptopt=None, workdir=None, manager=None):
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
        # TODO: Redesign API
        warnings.warn("get_irred_perts is still under development.")

        # Avoid modifications in self.
        inp = self.deepcopy()

        if qpt is None:
            qpt = inp.get("qpt")
            if qpt is None: raise ValueError("qpt must be specified")

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
        print(task.workdir)

        # Parse the file to get the perturbations.
        try:
            return yaml_read_irred_perts(task.log_file.path)
        except Exception as exc:
            # Try to understand if it's a problem with the Abinit input.
            report = task.get_event_report()
            if report.errors: raise self.Error(str(report))
            raise self.Error("Problem in temp Task executed in %s\n%s" % (task.workdir, exc))
