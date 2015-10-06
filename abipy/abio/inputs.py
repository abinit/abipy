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
from monty.json import MontyEncoder, MontyDecoder
from pymatgen.core.units import Energy
from pymatgen.serializers.json_coders import PMGSONable, pmg_serialize
from pymatgen.io.abinit.pseudos import PseudoTable, Pseudo
from pymatgen.io.abinit.tasks import AbinitTask, ParalHintsParser
from pymatgen.io.abinit.netcdf import NetcdfReader
from pymatgen.io.abinit.abiinspect import yaml_read_irred_perts
from abipy.core.structure import Structure
from abipy.core.mixins import Has_Structure
from abipy.htc.variable import InputVariable
from abipy.abio.abivars import is_abivar, is_anaddb_var
from abipy.abio.abivars_db import get_abinit_variables
from abipy.abio.input_tags import *

import logging
logger = logging.getLogger(__file__)


# List of Abinit variables used to specify the structure.
# This variables should not be passed to set_vars since 
# they will be generated with structure.to_abivars()
_GEOVARS = set([
    "acell",
    "rprim",
    "rprimd"
    "angdeg",
    "xred",
    "xcart",
    "xangst",
    "znucl",
    "typat",
    "ntypat",
    "natom",
])

# Variables defining tolerances (used in pop_tolerances)
_TOLVARS = set([
    'toldfe', 
    'tolvrs', 
    'tolwfr', 
    'tolrff',
    "toldff",  
    "tolimg", # ?
    "tolmxf", 
    "tolrde",  
])

# Variables defining tolerances for the SCF cycle that are mutally exclusive
_TOLVARS_SCF = set([
    'toldfe',
    'tolvrs',
    'tolwfr',
    'tolrff',
    "toldff",
])

# Variables determining if data files should be read in input
_IRDVARS = set([
    "irdbseig",
    "irdbsreso",
    "irdhaydock",
    "irdddk",
    "irdden",
    "ird1den",
    "irdqps",
    "irdkss",
    "irdscr",
    "irdsuscep",
    "irdvdw",
    "irdwfk",
    "irdwfkfine",
    "irdwfq",
    "ird1wf",
])


class AbstractInput(six.with_metaclass(abc.ABCMeta, MutableMapping, object)):

    # ABC protocol: __delitem__, __getitem__, __iter__, __len__, __setitem__
    def __delitem__(self, key):
        return self.vars.__delitem__(key)
        
    def __getitem__(self, key):
        return self.vars.__getitem__(key)
                                                                             
    def __iter__(self):
        return self.vars.__iter__()
                                                                             
    def __len__(self):
        return len(self.vars)
                                                                             
    def __setitem__(self, key, value):
        self._check_varname(key)
        return self.vars.__setitem__(key, value)

    def __repr__(self):
        return "<%s at %s>" % (self.__class__.__name__, id(self))

    def __str__(self):
        return self.to_string()

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

    def deepcopy(self):
        """Deep copy of the input."""
        return copy.deepcopy(self)

    def set_vars(self, *args, **kwargs):
        """Set the value of the variables"""
        kwargs.update(dict(*args))
        for varname, varvalue in kwargs.items():
            self[varname] = varvalue
        return kwargs

    def add_abiobjects(self, *abi_objects):
        """
        This function receive a list of `AbiVarable` objects and add 
        the corresponding variables to the input.
        """
        d = {}
        for aobj in abi_objects:
            if not hasattr(aobj, "to_abivars"):
                raise ValueError("type %s: %s does not have `to_abivars` method" % (type(aobj), repr(aobj)))
            d.update(self.set_vars(aobj.to_abivars()))
        return d

    def remove_vars(self, keys, strict=True):
        """
        Remove the variables listed in keys.
        Return dictionary with the variables that have been removed.
        """
        removed = {}
        for key in list_strings(keys):
            if strict and key not in self:
                raise KeyError("key: %s not in self:\n %s" % (key, list(self.keys())))
            if key in self:
                removed[key] = self.pop(key)

        return removed

    @abc.abstractproperty
    def vars(self):
        """Dictionary with the input variables. Used to implement dict-like interface."""

    @abc.abstractmethod
    def _check_varname(self, key):
        """Check if key is a valid name. Raise self.Error if not valid."""

    @abc.abstractmethod
    def to_string(self):
        """Returns a string with the input."""



class AbinitInputError(Exception):
    """Base error class for exceptions raised by `AbiInput`"""


# TODO: API to understand if one can use time-reversal symmetry and/or spatial symmetries
#       Very important especially when we have to select the value of kptopt

class AbinitInput(six.with_metaclass(abc.ABCMeta, AbstractInput, PMGSONable, Has_Structure, object)):
    """
    This object stores the ABINIT variables for a single dataset.
    """
    Error = AbinitInputError

    def __init__(self, structure, pseudos, pseudo_dir=None, comment=None, decorators=None, abi_args=None,
                 abi_kwargs=None, tags=None):
        """
        Args:
            structure: Parameters defining the crystalline structure. Accepts :class:`Structure` object 
            file with structure (CIF, netcdf file, ...) or dictionary with ABINIT geo variables.
            pseudos: Pseudopotentials to be used for the calculation. Accepts: string or list of strings with the name 
                of the pseudopotential files, list of :class:`Pseudo` objects or :class:`PseudoTable` object.
            pseudo_dir: Name of the directory where the pseudopotential files are located.
            ndtset: Number of datasets.
            comment: Optional string with a comment that will be placed at the beginning of the file.
            decorators: List of `AbinitInputDecorator` objects.
            abi_args: list of tuples (key, value) with the initial set of variables. Default: Empty
            abi_kwargs: Dictionary with the initial set of variables. Default: Empty
            tags: list/set of tags describing the input
        """
        # Internal dict with variables. we use an ordered dict so that 
        # variables will be likely grouped by `topics` when we fill the input.
        abi_args = [] if abi_args is None else abi_args
        for key, value in abi_args:
            self._check_varname(key)

        abi_kwargs = {} if abi_kwargs is None else abi_kwargs
        for key in abi_kwargs:
            self._check_varname(key)

        args = list(abi_args)[:]
        args.extend(list(abi_kwargs.items()))
        #print(args)

        self._vars = OrderedDict(args)

        self.set_structure(structure)

        if pseudo_dir is not None:
            pseudo_dir = os.path.abspath(pseudo_dir)
            if not os.path.exists(pseudo_dir): raise self.Error("Directory  %s does not exist")
            pseudos = [os.path.join(pseudo_dir, p) for p in list_strings(pseudos)]

        try:
            self._pseudos = PseudoTable.as_table(pseudos).get_pseudos_for_structure(self.structure)
        except ValueError as exc:
            raise self.Error(str(exc))

        if comment is not None: self.set_comment(comment)

        self._decorators = [] if not decorators else decorators

        self.tags = set() if not tags else set(tags)

    @pmg_serialize
    def as_dict(self):
        #vars = OrderedDict()
        # Use a list of (key, value) to serialize the OrderedDict
        abi_args = []
        for key, value in self.items():
            if isinstance(value, np.ndarray): value = value.tolist()
            #vars[key] = value
            abi_args.append((key, value))

        return dict(structure=self.structure.as_dict(),
                    pseudos=[p.as_dict() for p in self.pseudos], 
                    comment=self.comment,
                    decorators=[dec.as_dict() for dec in self.decorators],
                    abi_args=abi_args,
                    tags=list(self.tags))

    @property
    def vars(self):
        return self._vars

    @classmethod
    def from_dict(cls, d):
        pseudos = [Pseudo.from_file(p['filepath']) for p in d['pseudos']]
        dec = MontyDecoder()
        return cls(d["structure"], pseudos, decorators=dec.process_decoded(d["decorators"]),
                   comment=d["comment"], abi_args=d["abi_args"], tags=d["tags"])

    def __setitem__(self, key, value):
        if key in _TOLVARS_SCF and hasattr(self, '_vars') and any(t in self._vars and t != key for t in _TOLVARS_SCF):
            logger.info("Replacing previously set tolerance variable: {}."
                        .format(self.remove_vars(_TOLVARS_SCF, strict=False)))

        return super(AbinitInput, self).__setitem__(key, value)

    def _check_varname(self, key):
        if not is_abivar(key):
            raise self.Error("%s is not a valid ABINIT variable.\n"
                             "If you are sure the name is correct, please contact the abipy developers\n" 
                             "or modify the JSON file abipy/abio/abinit_vars.json" % key)

        if key in _GEOVARS:
            raise self.Error("You cannot set the value of a variable associated to the structure. Use set_structure")

    #def __eq__(self, other)
    #def __ne__(self, other)
    #    return not self.__eq__(other)

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

    @property
    def runlevel(self):
        """
        A Set of strings defining the type of run of the current input.
        """
        optdriver = self.get("optdriver")

        if optdriver is None:
            if self.get("rfddk") or self.get("rfelfd") or self.get("rfphon") or self.get("rfstrs"):
                optdriver = 1
            else:
                optdriver = 0

        runlevel = set()
        if optdriver == 0:
            runlevel.add(GROUND_STATE)
            iscf = self.get("iscf", 17 if self.pseudos[0].ispaw else 7)
            ionmov = self.get("ionmov", 0)
            optcell = self.get("optcell", 0)
            if ionmov == 0:
                if iscf < -1:
                    runlevel.add(NSCF)
                    if self.get("kptbounds") is not None:
                        runlevel.add(BANDS)
                else:
                    runlevel.add(SCF)
            elif ionmov in (2, 3, 4, 5, 7, 10, 11, 20):
                runlevel.add(RELAX)
                if optcell == 0:
                    runlevel.add(ION_RELAX)
                else:
                    runlevel.add(IONCELL_RELAX)
            elif ionmov in [1, 6, 8, 9, 12, 13, 14, 23]:
                runlevel.add(MOLECULAR_DYNACMICS)
        elif optdriver == 1:
            runlevel.add(DFPT)
            rfelfd = self.get("rfelfd")
            rfphon = self.get("rfphon")
            if self.get("rfddk") == 1 or rfelfd == 2:
                runlevel.add(DDK)
            elif rfelfd == 3:
                if rfphon == 1:
                    runlevel.add(BEC)
                else:
                    runlevel.add(DDE)
            elif rfphon == 1:
                runlevel.add(PH_Q_PERT)
            elif self.get("rfstrs ") > 0:
                runlevel.add(STRAIN)
        elif optdriver == 3:
            runlevel.update([MANY_BODY, SCREENING])
        elif optdriver == 4:
            runlevel.update([MANY_BODY, SIGMA])
            gwcalctyp = self.get("gwcalctyp")
            if gwcalctyp > 100:
                runlevel.add(HYBRID)
        elif optdriver == 99:
            runlevel.update([MANY_BODY, BSE])

        return runlevel

    @property
    def decorators(self):
        return self._decorators

    def register_decorator(self, decorator):
        """Register a :class:`AbinitInputDecorator`."""
        self._decorators.append(decorator)

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

    def to_string(self, sortmode="section", post=None, with_mnemonics=False, with_structure=True, with_pseudos=True):
        """
        String representation.

        Args:
            sortmode: "a" for alphabetical order, None if no sorting is wanted
            with_mnemonics: True if mnemonics should be added.
            post: String that will be appended to the name of the variables
                Note that post is usually autodetected when we have multiple datatasets
                It is mainly used when we have an input file with a single dataset
                so that we can prevent the code from adding "1" to the name of the variables 
                (In this case, indeed, Abinit complains if ndtset=1 is not specified 
                and we don't want ndtset=1 simply because the code will start to add 
                _DS1_ to all the input and output files.
            with_structure: False if section with structure variables should not be printed.
            with_pseudos: False if JSON section with pseudo data should not be added.
        """
        lines = []
        app = lines.append

        if self.comment: app("# " + self.comment.replace("\n", "\n#"))

        post = post if post is not None else ""

        mnemonics = self.mnemonics
        if with_mnemonics: mnemonics = with_mnemonics

        if mnemonics or sortmode == "section": 
            var_database = get_abinit_variables()

        if sortmode in (None, "a"):
            # Default is no sorting else alphabetical order.
            keys = list(self.keys())
            if sortmode == "a": keys = sorted(keys)

            # Extract the items from the dict and add the geo variables at the end
            items = list(self.items())
            if with_structure:
                items.extend(list(self.structure.to_abivars().items()))

            for name, value in items:
                if mnemonics:
                    app("# <" + var_database[name].definition + ">")

                # Build variable, convert to string and append it
                app(str(InputVariable(name + post, value)))

        elif sortmode == "section":
            # Group variables by section.
            # Get dict mapping section_name --> list of variable names belonging to the section.
            sec2names = var_database.group_by_section(list(self.keys()))
            w = 92

            for sec, names in sec2names.items():
                app(w * "#")
                app("#" + ("SECTION: %s" % sec).center(w - 1))
                app(w * "#")
                for name in names:
                    value = self[name]
                    if mnemonics:
                        app("# <" + var_database[name].definition + ">")

                    # Build variable, convert to string and append it
                    app(str(InputVariable(name + post, value)))

            if with_structure:
                app(w * "#")
                app("#" + ("STRUCTURE").center(w - 1))
                app(w * "#")
                for name, value in self.structure.to_abivars().items():
                    if mnemonics:
                        app("# <" + var_database[name].definition + ">")
                    app(str(InputVariable(name + post, value)))

        else:
            raise ValueError("Unsupported value for sortmode %s" % str(sortmode))

        s = "\n".join(lines)

        if not with_pseudos: return s 

        # Add JSON section with pseudo potentials.
        ppinfo = ["\n\n\n#<JSON>"]
        d = {"pseudos": [p.as_dict() for p in self.pseudos]}
        ppinfo.extend(json.dumps(d, indent=4).splitlines())
        ppinfo.append("</JSON>")
                                                             
        return s + "\n#".join(ppinfo)

    @property
    def comment(self):
        try:
            return self._comment
        except AttributeError:
            return None

    def set_comment(self, comment):
        """Set a comment to be included at the top of the file."""
        self._comment = comment

    @property
    def structure(self):
        """The :class:`Structure` associated to this input."""
        return self._structure

    def set_structure(self, structure):
        self._structure = Structure.as_structure(structure)

        # Check volume
        m = self.structure.lattice.matrix
        if np.dot(np.cross(m[0], m[1]), m[2]) <= 0:
            raise self.Error("The triple product of the lattice vector is negative. Use structure abi_sanitize.")

        return self._structure

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

    def set_autospinat(self, default=0.6):
        """
        Set the variable spinat for collineat calculation in the format (0, 0, m) with the value of m determined
        with the following order of preference:

        1. If the site of the structure has a magmom setting, that is used.
        2. If the species on the site has a spin setting, that is used.
        3. If the species itself has a particular setting in the config file, that
           is used, e.g., Mn3+ may have a different magmom than Mn4+.
        4. The element symbol itself is checked in the config file.
        5. If there are no settings, the default value is used.

        """
        #TODO copy the values in abipy instead?
        # We rely on the values from Materials Project for defaults of the magnetic moments
        import pymatgen.io.vasp as vasp
        from monty.serialization import loadfn
        magmom_mp_conf = loadfn(os.path.join(os.path.dirname(vasp.__file__), "MPVaspInputSet.yaml"))['INCAR']['MAGMOM']

        spinat = []
        for site in self.structure:
            if hasattr(site, 'magmom'):
                spinat.append((0., 0., site.magmom))
            elif hasattr(site.specie, 'spin'):
                spinat.append((0., 0., site.specie.spin))
            elif str(site.specie) in magmom_mp_conf:
                spinat.append((0., 0., magmom_mp_conf.get(str(site.specie))))
            else:
                spinat.append((0., 0., magmom_mp_conf.get(site.specie.symbol, default)))

        return self.set_vars(spinat=spinat)

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
            print(varnames, values)
            raise self.Error("The number of variables must equal the number of lists")

        # TODO: group varnames and varvalues!
        #varnames = [t[0] for t in items]
        #values = [t[1] for t in items]

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

    def pop_tolerances(self):
        """
        Remove all the tolerance variables present in self. 
        Return dictionary with the variables that have been removed."""
        return self.remove_vars(_TOLVARS, strict=False)

    def pop_irdvars(self):
        """
        Remove all the ird variables present in self.
        Return dictionary with the variables that have been removed."""
        return self.remove_vars(_IRDVARS, strict=False)

    @property
    def scf_tolvar(self):
        """
        Returns the tolerance variable and value relative to the scf convergence.
        If more than one is present raise an error
        """

        tolvar = None
        value = None
        for t in _TOLVARS_SCF:
            if t in self and self[t]:
                if tolvar:
                    raise self.Error('More than one tolerance set.')
                tolvar = t
                value = self[t]

        return tolvar, value

    def make_ph_inputs_qpoint(self, qpt, tolerance=None):
        """
        This functions should be called with an input the represents a GS run.
        """
        if tolerance is None:
            tolerance = {"tolvrs": 1.0e-10}

        if len(tolerance) != 1 or any(k not in _TOLVARS for k in tolerance):
            raise self.Error("Invalid tolerance: %s" % tolerance)

        # Call Abinit to get the list of irred perts.
        perts = self.abiget_irred_phperts(qpt=qpt)

        # Build list of datasets (one input per perturbation)
        ph_inputs = MultiDataset.replicate_input(input=self, ndtset=len(perts))

        # TODO: Should propagate info on symmetries.
        # use time-reversal if Gamma
        kptopt = 3
        if np.allclose(qpt, 0): kptopt = 2

        for pert, ph_input in zip(perts, ph_inputs):
            # TODO this will work for phonons, but not for the other types of perturbations.
            rfdir = 3 * [0]
            rfdir[pert.idir -1] = 1

            ph_input.set_vars(
                rfphon=1,                           # Will consider phonon-type perturbation
                nqpt=1,                             # One wavevector is to be considered
                qpt=pert.qpt,                       # q-wavevector.
                rfatpol=[pert.ipert, pert.ipert], 
                rfdir=rfdir,
                kptopt=kptopt,                           
            )
            ph_input.pop_tolerances()
            ph_input.set_vars(tolerance)

        return ph_inputs

    def make_ddk_inputs(self, tolerance=None):
        """
        Return inputs for the DDK calculation.
        This functions should be called with an input the represents a GS run.
        """
        if tolerance is None:
            tolerance = {"tolwfr": 1.0e-22}

        if len(tolerance) != 1 or any(k not in _TOLVARS for k in tolerance):
            raise self.Error("Invalid tolerance: %s" % tolerance)

        if "tolvrs" in tolerance:
            raise self.Error("tolvrs should not be used in a DDK calculation")

        # Call Abinit to get the list of irred perts.
        #perts = self.abiget_irred_phperts(qpt=qpt)
        # TODO Add symmetries
        ddk_rfdirs = [(1,0,0), (0,1,0), (0,0,1)]

        # Build list of datasets (one input per perturbation)
        ddk_inputs = MultiDataset.replicate_input(input=self, ndtset=len(ddk_rfdirs))

        # See tutorespfn/Input/trf1_5.in
        for rfdir, ddk_input in zip(ddk_rfdirs, ddk_inputs):
            ddk_input.set_vars(
                rfelfd=2,             # Activate the calculation of the d/dk perturbation
                rfdir=rfdir,          # Direction of the per ddk.
                nqpt=1,               # One wavevector is to be considered
                qpt=(0, 0, 0),        # q-wavevector.
                kptopt=2,             # Take into account time-reversal symmetry.           
                iscf=-3,              # The d/dk perturbation must be treated in a non-self-consistent way
            )

            ddk_input.pop_tolerances()
            ddk_input.set_vars(tolerance)

        return ddk_inputs

    def make_dde_inputs(self, tolerance=None):
        """
        Return inputs for the calculation of the electric field perturbations.

        This functions should be called with an input the represents a gs run.
        """
        if tolerance is None:
            tolerance = {"tolvrs": 1.0e-10}

        if len(tolerance) != 1 or any(k not in _TOLVARS for k in tolerance):
            raise self.Error("Invalid tolerance: %s" % tolerance)

        # Call Abinit to get the list of irred perts.
        perts = self.abiget_irred_ddeperts()

        # Build list of datasets (one input per perturbation)
        multi = MultiDataset.replicate_input(input=self, ndtset=len(perts))

        # See tutorespfn/Input/trf1_5.in dataset 3
        for pert, inp in zip(perts, multi):
            rfdir = 3 * [0]
            rfdir[pert.idir -1] = 1

            inp.set_vars(
                rfelfd=3,             # Activate the calculation of the electric field perturbation
                rfdir=rfdir,          # Direction of the dde perturbation.
                nqpt=1,               # One wavevector is to be considered
                qpt=(0, 0, 0),        # q-wavevector.
                kptopt=2,             # Take into account time-reversal symmetry.
            )

            inp.pop_tolerances()
            inp.set_vars(tolerance)

        return multi

    def make_bec_inputs(self, tolerance=None):
        """
        Return inputs for the calculation of the Born effective charges.

        This functions should be called with an input the represents a gs run.
        """
        if tolerance is None:
            tolerance = {"tolvrs": 1.0e-10}

        if len(tolerance) != 1 or any(k not in _TOLVARS for k in tolerance):
            raise self.Error("Invalid tolerance: %s" % tolerance)

        # Call Abinit to get the list of irred perts.
        # TODO:
        # Check that one can use the same list of irred perts as in phonons
        perts = self.abiget_irred_phperts(qpt=(0, 0, 0))

        # Build list of datasets (one input per perturbation)
        multi = MultiDataset.replicate_input(input=self, ndtset=len(perts))

        # See tutorespfn/Input/trf1_5.in dataset 3
        for pert, inp in zip(perts, multi):
            rfdir = 3 * [0]
            rfdir[pert.idir -1] = 1

            inp.set_vars(
                rfphon=1,             # Activate the calculation of the atomic dispacement perturbations
                rfatpol=[pert.ipert, pert.ipert], 
                rfdir=rfdir,
                rfelfd=3,             # Activate the calculation of the electric field perturbation
                nqpt=1,               # One wavevector is to be considered
                qpt=(0, 0, 0),        # q-wavevector.
                kptopt=2,             # Take into account time-reversal symmetry.           
            )

            inp.pop_tolerances()
            inp.set_vars(tolerance)

        return multi

    def pycheck(self):
        errors = []
        eapp = errors.append

        m = self.structure.lattice.matrix
        volume = np.dot(np.cross(m[0], m[1]), m[2])
        if volume < 0:
            eapp("The triple product of the lattice vector is negative. Use structure abi_sanitize.")

        #if sel.ispaw and "pawecutdg not in self 
        #if errors: raise self.Error("\n".join(errors))

        return dict2namedtuple(errors=errors, warnings=warnings)
        
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
        This function computes the list of points in the IBZ and the corresponding weights.
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
        """
        # Avoid modifications in self.
        inp = self.deepcopy()

        # The magic value that makes ABINIT print the ibz and then stop.
        inp["prtkpt"] = -2

        if ngkpt is not None: inp["ngkpt"] = ngkpt
        if shiftk is not None:
            shiftk = np.reshape(shiftk, (-1,3))
            inp.set_vars(shiftk=shiftk, nshiftk=len(shiftk))

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

    def _abiget_irred_perts(self, perts_vars, qpt=None, ngkpt=None, shiftk=None, kptopt=None, workdir=None, manager=None):
        """
        This function, computes the list of irreducible perturbations for DFPT.
        It should be called with an input file that contains all the mandatory variables required by ABINIT.

        Args:
            perts_vars: list of variables to be added to get the appropriate perturbation
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

        """
        # Avoid modifications in self.
        inp = self.deepcopy()

        qpt = inp.get("qpt") if qpt is None else qpt
        if qpt is None:
            raise ValueError("qpt is not in the input and therefore it must be passed explicitly")

        if ngkpt is not None: inp["ngkpt"] = ngkpt
        if shiftk is not None:
            shiftk = np.reshape(shiftk, (-1,3))
            inp.set_vars(shiftk=shiftk, nshiftk=len(inp.shiftk))

        inp.set_vars(
            nqpt=1,                           # One wavevector is to be considered
            qpt=qpt,                          # q-wavevector.
            paral_rf=-1,                      # Magic value to get the list of irreducible perturbations for this q-point.
            **perts_vars
        )

        if kptopt is not None: inp["kptopt"] = kptopt

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

        """
        phperts_vars = dict(rfphon=1,                         # Will consider phonon-type perturbation
                            rfatpol=[1, len(self.structure)],  # Set of atoms to displace.
                            rfdir=[1, 1, 1],                  # Along this set of reduced coordinate axis.
                            )

        return self._abiget_irred_perts(phperts_vars, qpt=qpt, ngkpt=ngkpt, shiftk=shiftk, kptopt=kptopt,
                                        workdir=workdir, manager=manager)

    def abiget_irred_ddeperts(self, ngkpt=None, shiftk=None, kptopt=None, workdir=None, manager=None):
        """
        This function, computes the list of irreducible perturbations for DFPT.
        It should be called with an input file that contains all the mandatory variables required by ABINIT.

        Args:
            ngkpt: Number of divisions for the k-mesh (default None i.e. use ngkpt from self)
            shiftk: Shiftks (default None i.e. use shiftk from self)
            kptopt: Option for k-point generation. If None, the value in self is used.
            workdir: Working directory of the fake task used to compute the ibz. Use None for temporary dir.
            manager: :class:`TaskManager` of the task. If None, the manager is initialized from the config file.

        Returns:
            List of dictionaries with the Abinit variables defining the irreducible perturbation
            Example:

                [{'idir': 1, 'ipert': 4, 'qpt': [0.0, 0.0, 0.0]},
                 {'idir': 2, 'ipert': 4, 'qpt': [0.0, 0.0, 0.0]}]

        """
        ddeperts_vars = dict(rfphon=0,           # No phonon-type perturbation
                             rfelfd=3,           # Electric field
                             kptopt=2,           # kpt time reversal symmetry
                             )

        return self._abiget_irred_perts(ddeperts_vars, qpt=(0, 0, 0), ngkpt=ngkpt, shiftk=shiftk, kptopt=kptopt,
                                        workdir=workdir, manager=manager)

    def abiget_autoparal_pconfs(self, max_ncpus, autoparal=1, workdir=None, manager=None):
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
        except parser.Error as exc:
            # Try to understand if it's a problem with the Abinit input.
            report = task.get_event_report()
            if report and report.errors: raise self.Error(str(report))
            raise self.Error("Problem in temp Task executed in %s\n%s" % (task.workdir, exc))

    def add_tags(self, tags):
        """
        Add tags to the input
        Args:
            tags: A single tag or list/tuple/set of tags
        """
        if isinstance(tags, (list, tuple, set)):
            self.tags.update(tags)
        else:
            self.tags.add(tags)

    def remove_tags(self, tags):
        """
        Remove tags from the input
        Args:
            tags: A single tag or list/tuple/set of tags
        """
        if isinstance(tags, (list, tuple, set)):
            self.tags.difference_update(tags)
        else:
            self.tags.discard(tags)


class MultiDataset(object):
    """
    This object is essentially a list of :class:`AbinitInput` objects.
    that provides an easy-to-use interface to apply global changes to the 
    the inputs stored in the objects.

    Let's assume for example that multi contains two `AbinitInput` objects and we
    want to set `ecut` to 1 in both dictionaries. The direct approach would be:

        for inp in multi:
            inp.set_vars(ecut=1)

    or alternatively:

        for i in range(multi.ndtset):
            multi[i].set_vars(ecut=1)


    MultiDataset provides its own implementaion of __getattr__ so that one can simply use:

         multi.set_vars(ecut=1)

    .. warning::

        MultiDataset does not support calculations done with different sets of pseudopotentials.
        The inputs can have different crystalline structures (as long as the atom types are equal)
        but each input in MultiDataset must have the same set of pseudopotentials.
    """
    Error = AbinitInputError

    @classmethod
    def from_inputs(cls, inputs):
        for inp in inputs:
            if any(p1 != p2 for p1, p2 in zip(inputs[0].pseudos, inp.pseudos)):
                raise ValueError("Pseudos must be consistent when from_inputs is invoked.")

        # Build MultiDataset from input structures and pseudos and add inputs.
        multi = cls(structure=[inp.structure for inp in inputs], pseudos=inputs[0].pseudos, ndtset=len(inputs))

        for inp, new_inp in zip(inputs, multi):
            new_inp.set_vars(**inp)
            new_inp._decorators = inp.decorators
            new_inp.tags = set(inp.tags)

        return multi

    @classmethod
    def replicate_input(cls, input, ndtset):
        """Constructur a multidataset with ndtset from the :class:`AbinitInput` input."""
        multi = cls(input.structure, input.pseudos, ndtset=ndtset)

        for inp in multi: 
            inp.set_vars({k: v for k, v in input.items()})
            inp.tags = set(input.tags)

        return multi

    def __init__(self, structure, pseudos, pseudo_dir="", ndtset=1):
        """
        Args:
            structure: file with the structure, :class:`Structure` object or dictionary with ABINIT geo variable
                Accepts also list of objects that can be converted to Structure object.
                In this case, however, ndtset must be equal to the length of the list.
            pseudos: String or list of string with the name of the pseudopotential files.
            pseudo_dir: Name of the directory where the pseudopotential files are located.
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

        # Build the list of AbinitInput objects.
        if ndtset <= 0:
            raise ValueError("ndtset %d cannot be <=0" % ndtset)

        if not isinstance(structure, (list, tuple)):
            self._inputs = [AbinitInput(structure=structure, pseudos=pseudos) for i in range(ndtset)]
        else:
            assert len(structure) == ndtset
            self._inputs = [AbinitInput(structure=s, pseudos=pseudos) for s in structure]

        # Check pseudos
        #for i in range(self.ndtset):
        #    if any(p1 != p2 for p1, p2 in zip(self[0].pseudos, self[i].pseudos)):
        #        raise selfError("Pseudos must be consistent when from_inputs is invoked.")

    @property
    def ndtset(self):
        """Number of inputs in self."""
        return len(self)

    @property
    def pseudos(self):
        return self[0].pseudos

    @property
    def ispaw(self):
        """True if PAW calculation."""
        return all(p.ispaw for p in self.pseudos)

    @property
    def isnc(self):
        """True if norm-conserving calculation."""
        return all(p.isnc for p in self.pseudos)

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
            raise AttributeError("Cannot find attribute %s. Tried in %s and then in AbinitInput object" 
                                 % (self.__class__.__name__, name))
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

    def __add__(self, other):
        if isinstance(other, AbinitInput):
            new_mds = MultiDataset.from_inputs(self)
            new_mds.append(other)
            return new_mds
        elif isinstance(other, MultiDataset):
            new_mds = MultiDataset.from_inputs(self)
            new_mds.extend(other)
            return new_mds
        else:
            return NotImplemented

    def __radd__(self, other):
        if isinstance(other, AbinitInput):
            new_mds = MultiDataset.from_inputs([other])
            new_mds.extend(self)
        elif isinstance(other, MultiDataset):
            new_mds = MultiDataset.from_inputs(other)
            new_mds.extend(self)
        else:
            return NotImplemented

    def append(self, abinit_input):
        """Add a :class:`AbinitInput` to the list."""
        assert isinstance(abinit_input, AbinitInput)
        if any(p1 != p2 for p1, p2 in zip(abinit_input.pseudos, abinit_input.pseudos)):
            raise ValueError("Pseudos must be consistent when from_inputs is invoked.")
        self._inputs.append(abinit_input)

    def extend(self, abinit_inputs):
        """Extends self with a list of :class:`AbinitInput` objects."""
        assert all(isinstance(inp, AbinitInput) for inp in abinit_inputs)
        for inp in abinit_inputs:
            if any(p1 != p2 for p1, p2 in zip(self[0].pseudos, inp.pseudos)):
                raise ValueError("Pseudos must be consistent when from_inputs is invoked.")
        self._inputs.extend(abinit_inputs)

    def addnew_from(self, dtindex):
        self.append(self[dtindex].deepcopy())

    def split_datasets(self):
        return self._inputs

    def deepcopy(self):
        """Deep copy of the object."""
        return copy.deepcopy(self)

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
                #with_structure = True 
                #if same_structure and not is_last: with_structure = False

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

    def filter_by_tags(self, tags):
        """
        Filters the input according to the tags
        Args:
            tags: A single tag or list/tuple/set of tags
        Returns:
            A multidata containing the inputs containing all the requested tags
        """
        if isinstance(tags, (list, tuple, set)):
            tags = set(tags)
        elif not isinstance(tags, set):
            tags = {tags}

        inputs = [i for i in self if tags.issubset(i.tags)]

        return MultiDataset.from_inputs(inputs) if inputs else None

    def add_tags(self, tags, dtindeces=None):
        """
        Add tags to the selected inputs
        Args:
            tags: A single tag or list/tuple/set of tags
            dtindeces: a list of indices to which the tags will be added. None=all the inputs.
        """
        for i in dtindeces if dtindeces else range(len(self)):
            self[i].add_tags(tags)

    def remove_tags(self, tags, dtindeces=None):
        """
        Remove tags from the selected inputs
        Args:
            tags: A single tag or list/tuple/set of tags
            dtindeces: a list of indices from which the tags will be removed. None=all the inputs.
        """
        for i in dtindeces if dtindeces else range(len(self)):
            self[i].remove_tags(tags)

    def filter_by_runlevel(self, runlevel):
        if isinstance(runlevel, (list, tuple, set)):
            runlevel = set(runlevel)
        elif not isinstance(runlevel, set):
            runlevel = {runlevel}

        inputs = [i for i in self if runlevel.issubset(i.runlevel)]

        return MultiDataset.from_inputs(inputs) if inputs else None


class AnaddbInputError(Exception):
    """Base error class for exceptions raised by `AnaddbInput`"""


class AnaddbInput(AbstractInput, Has_Structure):

    Error = AnaddbInputError

    def __init__(self, structure, comment="", anaddb_args=None, anaddb_kwargs=None):
        """
        Args:
            structure: :class:`Structure` object 
            comment: Optional string with a comment that will be placed at the beginning of the file.
            anaddb_args: List of tuples (key, value) with Anaddb input variables (default: empty)
            anaddb_kwargs: Dictionary with Anaddb input variables (default: empty)
        """
        self._structure = structure
        self.comment = comment

        anaddb_args = [] if anaddb_args is None else anaddb_args
        for key, value in anaddb_args:
            self._check_varname(key)

        anaddb_kwargs = {} if anaddb_kwargs is None else anaddb_kwargs
        for key in anaddb_kwargs:
            self._check_varname(key)

        args = list(anaddb_args)[:]
        args.extend(list(anaddb_kwargs.items()))

        self._vars = OrderedDict(args)

    @property
    def vars(self):
        return self._vars

    def _check_varname(self, key):
        if not is_anaddb_var(key):
            raise self.Error("%s is not a registered Anaddb variable\n"
                             "If you are sure the name is correct, please contact the abipy developers\n" 
                             "or modify the JSON file abipy/abio/anaddb_vars.json" % key)

    @classmethod
    def modes_at_qpoint(cls, structure, qpoint, asr=2, chneut=1, dipdip=1, 
                        anaddb_args=None, anaddb_kwargs=None):
        """
        Input file for the calculation of the phonon frequencies at a given q-point.

        Args:
            Structure: :class:`Structure` object
            qpoint: Reduced coordinates of the q-point where phonon frequencies and modes are wanted
            asr, chneut, dipdp: Anaddb input variable. See official documentation.
            anaddb_args: List of tuples (key, value) with Anaddb input variables (default: empty)
            anaddb_kwargs: Dictionary with Anaddb input variables (default: empty)
        """
        new = cls(structure, comment="ANADB input for phonon frequencies at one q-point", 
                  anaddb_args=anaddb_args, anaddb_kwargs=anaddb_kwargs)

        # We need a numpy array.
        if hasattr(qpoint, "frac_coords"):
            qpoint = qpoint.frac_coords
        else:
            qpoint = np.array(qpoint)

        if len(qpoint) != 3:
            #print(type(qpoint), qpoint.shape)
            raise ValueError("Wrong q-point %s" % qpoint)

        new.set_vars(
            ifcflag=1,        # Interatomic force constant flag
            asr=asr,          # Acoustic Sum Rule
            chneut=chneut,    # Charge neutrality requirement for effective charges.
            dipdip=dipdip,    # Dipole-dipole interaction treatment
            # This part is fixed
            ngqpt=(1, 1, 1), 
            nqshft=1,         
            q1shft=qpoint,
            nqpath=2,
            # FIXME
            # ndivsm requires at least two q-points
            qpath=np.array([qpoint, qpoint + 1]).ravel(),
            ndivsm=1
        )

        return new

    @classmethod
    def piezo_elastic(cls, structure, anaddb_args=None, anaddb_kwargs=None):
        new = cls(structure, comment="ANADB input for piezoelectric and elastic tensor calculation",
                  anaddb_args=anaddb_args, anaddb_kwargs=anaddb_kwargs)

        new.set_vars(
            elaflag=3,
            piezoflag=3,
            instrflag=1,
            chneut=1,
            asr=0,
            symdynmat=1
        )

        return new

    #@classmethod
    #def phbands(cls, structure, ngqpt, nqsmall, q1shft=(0,0,0), asr=2, chneut=0, dipdip=1, 
    #           anaddb_args=None, anaddb_kwargs=None):
    #    """
    #    Build an anaddb input file for the computation of phonon band structure.
    #    """
    #    return self.phbands_and_dos(structure, ngqpt, nqsmall, ndivsm=20, q1shft=(0,0,0),
    #                                qptbounds=None, asr=2, chneut=0, dipdip=1, dos_method="tetra", 
    #                                anaddb_args=anaddb_args, anaddb_kwargs=anaddb_kwargs)

    #@classmethod
    #def phdos(cls, structure, ngqpt, nqsmall, q1shft=(0,0,0), asr=2, chneut=0, dipdip=1, dos_method="tetra", 
    #           anaddb_args=None, anaddb_kwargs=None):
    #    """
    #    Build an anaddb input file for the computation of phonon DOS.
    #    """
    #    return self.phbands_and_dos(structure, ngqpt, nqsmall, ndivsm=20, q1shft=(0,0,0),
    #                                qptbounds=None, asr=2, chneut=0, dipdip=1, dos_method="tetra", 
    #                                anaddb_args=anaddb_args, anaddb_kwargs=anaddb_kwargs)

    @classmethod
    def phbands_and_dos(cls, structure, ngqpt, nqsmall, ndivsm=20, q1shft=(0,0,0),
                        qptbounds=None, asr=2, chneut=0, dipdip=1, dos_method="tetra", lo_to_splitting=False,
                        anaddb_args=None, anaddb_kwargs=None):
        """
        Build an anaddb input file for the computation of phonon bands and phonon DOS.

        Args:
            structure: :class:`Structure` object
            ngqpt: Monkhorst-Pack divisions for the phonon Q-mesh (coarse one)
            nqsmall: Used to generate the (dense) mesh for the DOS.
                It defines the number of q-points used to sample the smallest lattice vector.
            ndivsm: Used to generate a normalized path for the phonon bands.
                If gives the number of divisions for the smallest segment of the path.
            q1shft: Shifts used for the coarse Q-mesh
            qptbounds Boundaries of the path. If None, the path is generated from an internal database
                depending on the input structure.
            asr, chneut, dipdp: Anaddb input variable. See official documentation.
            dos_method: Possible choices: "tetra", "gaussian" or "gaussian:0.001 eV".
                In the later case, the value 0.001 eV is used as gaussian broadening
            lo_to_splitting: if True calculation of the LO-TO splitting will be included
            anaddb_args: List of tuples (key, value) with Anaddb input variables (default: empty)
            anaddb_kwargs: Dictionary with Anaddb input variables (default: empty)
        """
        dosdeltae, dossmear = None, None

        if dos_method == "tetra":
            prtdos = 2
        elif "gaussian" in dos_method:
            prtdos = 1
            i = dos_method.find(":")
            if i != -1:
                value, eunit = dos_method[i+1:].split()
                dossmear = Energy(float(value), eunit).to("Ha")
        else:
            raise cls.Error("Wrong value for dos_method: %s" % dos_method)

        new = cls(structure, comment="ANADB input for phonon bands and DOS", 
                  anaddb_args=anaddb_args, anaddb_kwargs=anaddb_kwargs)

        # Parameters for the dos.
        new.set_autoqmesh(nqsmall)
        new.set_vars(prtdos=prtdos, dosdeltae=dosdeltae, dossmear=dossmear)

        new.set_qpath(ndivsm, qptbounds=qptbounds)
        qptbounds = new['qpath']
        q1shft = np.reshape(q1shft, (-1, 3))

        new.set_vars(
            ifcflag=1,
            ngqpt=np.array(ngqpt),
            q1shft=q1shft,
            nqshft=len(q1shft),
            asr=asr,
            chneut=chneut,
            dipdip=dipdip,
        )

        if lo_to_splitting:
            directions = []
            for i, qpt in enumerate(qptbounds):
                if np.array_equal(qpt, (0, 0, 0)):
                    # anaddb expects cartesian coordinates for the qph2l list
                    if i>0:
                        directions.extend(structure.lattice.reciprocal_lattice_crystallographic.get_cartesian_coords(qptbounds[i-1]))
                        directions.append(0)
                    if i<len(qptbounds)-1:
                        directions.extend(structure.lattice.reciprocal_lattice_crystallographic.get_cartesian_coords(qptbounds[i+1]))
                        directions.append(0)

            if directions:
                directions = np.reshape(directions, (-1, 4))
                new.set_vars(
                    nph2l=len(directions),
                    qph2l=directions
                )

        return new

    @classmethod
    def thermo(cls, structure, ngqpt, nqsmall, q1shft=(0, 0, 0), nchan=1250, nwchan=5, thmtol=0.5,
               ntemper=199, temperinc=5, tempermin=5., asr=2, chneut=1, dipdip=1, ngrids=10, 
               anaddb_args=None, anaddb_kwargs=None):

        """
        Build an anaddb input file for the computation of phonon bands and phonon DOS.

        Args:
            structure: :class:`Structure` object
            ngqpt: Monkhorst-Pack divisions for the phonon Q-mesh (coarse one)
            nqsmall: Used to generate the (dense) mesh for the DOS.
                It defines the number of q-points used to sample the smallest lattice vector.
            q1shft: Shifts used for the coarse Q-mesh
            nchan:
            nwchan:
            thmtol:
            ntemper:
            temperinc:
            tempermin:
            asr, chneut, dipdp: Anaddb input variable. See official documentation.
            ngrids:
            anaddb_args: List of tuples (key, value) with Anaddb input variables (default: empty)
            anaddb_kwargs: Dictionary with Anaddb input variables (default: empty)

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
        new = cls(structure, comment="ANADB input for thermodynamics", 
                  anaddb_args=anaddb_args, anaddb_kwargs=anaddb_kwargs)
        new.set_autoqmesh(nqsmall)

        q1shft = np.reshape(q1shft, (-1, 3))

        new.set_vars(
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
    def modes(cls, structure, enunit=2, asr=2, chneut=1, anaddb_args=None, anaddb_kwargs=None):
        """
        Build an anaddb input file for the computation of phonon modes.

        Args:
            Structure: :class:`Structure` object
            ngqpt: Monkhorst-Pack divisions for the phonon Q-mesh (coarse one)
            nqsmall: Used to generate the (dense) mesh for the DOS.
                It defines the number of q-points used to sample the smallest lattice vector.
            q1shft: Shifts used for the coarse Q-mesh
            qptbounds Boundaries of the path. If None, the path is generated from an internal database
                depending on the input structure.
            asr, chneut, dipdp: Anaddb input variable. See official documentation.
            anaddb_args: List of tuples (key, value) with Anaddb input variables (default: empty)
            anaddb_kwargs: Dictionary with Anaddb input variables (default: empty)

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
        # Wavevector list number 1
        #nph1l     1
        #qph1l   0.0  0.0  0.0    1.0   ! (Gamma point)
        #!Wavevector list number 2
        #nph2l     3      ! number of phonons in list 1
        #qph2l   1.0  0.0  0.0    0.0
        #        0.0  1.0  0.0    0.0
        #        0.0  0.0  1.0    0.0
        """
        new = cls(structure, comment="ANADB input for modes", anaddb_args=anaddb_args, anaddb_kwargs=anaddb_kwargs)

        new.set_vars(
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

    def to_string(self, sortmode=None):
        """
        String representation.

        Args:
            sortmode: "a" for alphabetical order, None if no sorting is wanted
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

    def set_qpath(self, ndivsm, qptbounds=None):
        """
        Set the variables for the computation of the phonon band structure.

        Args:
            ndivsm: Number of divisions for the smallest segment.
            qptbounds: q-points defining the path in k-space.
                If None, we use the default high-symmetry k-path defined in the pymatgen database.
        """
        if qptbounds is None: qptbounds = self.structure.calc_kptbounds()
        qptbounds = np.reshape(qptbounds, (-1, 3))

        return self.set_vars(ndivsm=ndivsm, nqpath=len(qptbounds), qpath=qptbounds)

    def set_autoqmesh(self, nqsmall):
        """
        Set the variable nqpt for the sampling of the BZ.

        Args:
            nqsmall: Number of divisions used to sample the smallest lattice vector.
        """
        return self.set_vars(ng2qpt=self.structure.calc_ngkpt(nqsmall))


class OpticVar(collections.namedtuple("OpticVar", "name default group help")):
    def __str__(self):
        sval = str(self.default)
        return (4*" ").join(sval, "!" + self.help)


class OpticError(Exception):
    """Error class raised by OpticInput."""


# TODO: OpticInput should implement AbstractInput!
#class OpticInput(collections.MutableMapping):
class OpticInput(AbstractInput):
    """
    abo_1WF7      ! Name of the first d/dk response wavefunction file, produced by abinit
    abo_1WF8      ! Name of the second d/dk response wavefunction file, produced by abinit
    abo_1WF9      ! Name of the third d/dk response wavefunction file, produced by abinit
    abo_WFK       ! Name of the ground-state wavefunction file, produced by abinit
    0.01          ! Value of the *smearing factor*, in Hartree
    0.010   1     ! frequency *step* and *maximum* frequency (Ha)
    0.000         ! *Scissor* shift if needed, in Hartree
    0.001         ! *Tolerance* on closeness of singularities (in Hartree)
    3             ! *Number of components* of linear optic tensor to be computed
    11 33 23      ! Linear *coefficients* to be computed (x=1, y=2, z=3)
    2             ! Number of components of nonlinear optic tensor to be computed
    123 222       ! Non-linear coefficients to be computed
    """

    """
    &FILES
     ddkfile_1 = 'abo_1WF7',
     ddkfile_2 = 'abo_1WF8',
     ddkfile_3 = 'abo_1WF9',
     wfkfile = 'abo_WFK'
    /
    &PARAMETERS
     broadening = 0.002,
     domega = 0.0003,
     maxomega = 0.3,
     scissor = 0.000,
     tolerance = 0.002
    /
    &COMPUTATIONS
     num_lin_comp = 1,
     lin_comp = 11,
     num_nonlin_comp = 2,
     nonlin_comp = 123,222,
     num_linel_comp = 0,
     num_nonlin2_comp = 0,
    /
    """


    Error = OpticError

    # variable name --> default value.
    _VARIABLES = [
        #OpticVar(name="ddkfile_x", default=None, help="Name of the first d/dk response wavefunction file"),
        #OpticVar(name="ddkfile_y", default=None, help="Name of the second d/dk response wavefunction file"),
        #OpticVar(name="ddkfile_z", default=None, help="Name of the third d/dk response wavefunction file"),
        #OpticVar(name="wfkfile",   default=None, help="Name of the ground-state wavefunction file"),
        OpticVar(name="broadening",      default=0.01, group='PARAMETERS', help="Value of the *smearing factor*, in Hartree"),
        OpticVar(name="domega",     default=0.010, group='PARAMETERS', help="Frequency *step* (Ha)"),
        OpticVar(name="maxomega",     default=1, group='PARAMETERS', help="Maximum frequency (Ha)"),
        OpticVar(name="scissor",   default=0.000, group='PARAMETERS', help="*Scissor* shift if needed, in Hartree"),
        OpticVar(name="tolerance",  default=0.001, group='PARAMETERS', help="*Tolerance* on closeness of singularities (in Hartree)"),
        OpticVar(name="num_lin_comp", default=0, group='COMPUTATIONS', help="*Number of components* of linear optic tensor to be computed"),
        OpticVar(name="lin_comp",     default=0, group='COMPUTATIONS', help="Linear *coefficients* to be computed (x=1, y=2, z=3)"),
        OpticVar(name="num_nonlin_comp", default=0, group='COMPUTATIONS', help="Number of components of nonlinear optic tensor to be computed"),
        OpticVar(name="nonlin_comp", default=0, group='COMPUTATIONS', help="Non-linear coefficients to be computed"),
        OpticVar(name="num_linel_comp", default=0, group='COMPUTATIONS', help="Number of components of linear electro-optic tensor to be computed"),
        OpticVar(name="linel_comp", default=0, group='COMPUTATIONS', help="Linear electro-optic coefficients to be computed"),
        OpticVar(name="num_nonlin2_comp", default=0, group='COMPUTATIONS', help="Number of components of nonlinear optic tensor v2 to be computed"),
        OpticVar(name="nonlin2_comp", default=0, group='COMPUTATIONS', help="Non-linear coefficients v2 to be computed"),
    ]
 
    _GROUPS = ['PARAMETERS','COMPUTATIONS']

    # Variable names supported
    _VARNAMES = [v.name for v in _VARIABLES]

    # Mapping name --> var object.
    _NAME2VAR = {v.name: v for v in _VARIABLES}

    def __init__(self, **kwargs):
        # Init with default values.
        self._vars = collections.OrderedDict((v.name, v.default) for v in self._VARIABLES)

        # Update the variables with the values passed by the user
        for k, v in kwargs.items():
            if k not in self._VARNAMES:
                raise self.Error("varname %s not in %s" % (k, str(self._VARNAMES)))
            self[k] = v

    def __str__(self):
        return self.to_string()

    @property
    def vars(self):
        return self._vars

    def _check_varname(self, key):
        if key not in self._VARNAMES:
            raise self.Error("%s is not a valid optic variable.\n"
                             "If you are sure the name is correct, please change the _VARIABLES list in:\n%s"  % 
                             (key, __file__))

    def get_default(self, key):
        for var in self._VARIABLES:
            if var.name == key: return var.default
        raise KeyError("Cannot find %s in _VARIABLES" % key)

    def as_dict(self):
        my_dict = OrderedDict()
        for grp in self._GROUPS:
            my_dict[grp] = OrderedDict()
        for name in self._VARNAMES:
            value = self.vars.get(name)
            if value is None: value = self.get_default(name)
            if value is None:
                raise self.Error("Variable %s is missing" % name)
           
            var = self._NAME2VAR[name]
            grp = var.group
            my_dict[grp].update({name : value})

        return my_dict

    def to_string(self):
        table = []
        app = table.append

        for name in self._VARNAMES:
            value = self.vars.get(name)
            if value is None: value = self.get_default(name)
            if value is None:
                raise self.Error("Variable %s is missing" % name)

            # One line per variable --> valperline set to None
            variable = InputVariable("", value, valperline=None)
            app([str(variable).strip(), "! " + self._NAME2VAR[name].help])

        # Align
        width = max(len(row[0]) for row in table)
        lines = []
        for row in table:
            s = row[0].ljust(width) + "\t" + row[1]
            lines.append(s)

        return "\n".join(lines)
