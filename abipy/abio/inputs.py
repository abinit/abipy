"""
This module defines objects to facilitate the creation of ABINIT input files.
The syntax is similar to the one used in ABINIT with small differences.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import collections
import warnings
import itertools
import copy
import time
import six
import abc
import json
import numpy as np

from collections import OrderedDict
try: # py3k
    from collections.abc import MutableMapping
except ImportError:
    from collections import MutableMapping
from monty.collections import dict2namedtuple
from monty.string import is_string, list_strings
from monty.json import MontyDecoder, MSONable
from pymatgen.core.units import Energy
try:
    from pymatgen.util.serialization import pmg_serialize
except ImportError:
    from pymatgen.serializers.json_coders import pmg_serialize

from abipy.tools.numtools import is_diagonal
from abipy.core.structure import Structure
from abipy.core.mixins import Has_Structure
from abipy.core.kpoints import has_timrev_from_kptopt
from abipy.abio.variable import InputVariable
from abipy.abio.abivars import is_abivar, is_anaddb_var
from abipy.abio.abivars_db import get_abinit_variables
from abipy.abio.input_tags import *
from abipy.flowtk import PseudoTable, Pseudo, AbinitTask, AnaddbTask, ParalHintsParser, NetcdfReader
from abipy.flowtk.abiinspect import yaml_read_irred_perts
from abipy.flowtk import abiobjects as aobj
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.io.abinit.abiobjects import KSampling

import logging
logger = logging.getLogger(__file__)


# List of Abinit variables used to specify the structure.
# This variables should not be passed to set_vars since
# they will be generated with structure.to_abivars()
GEOVARS = set([
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

# FIXME __mul__ operator in pymatgen should allow for grouping atoms by individual cells
# The present version group by image.
#def _repeat_array(name, values, from_natom, numcells):
#    if name in {"spinat",}:
#        # [:, natom]
#        values = np.reshape(values, (-1, from_natom))
#        return np.repeat(values, numcells, axis=0)
#    else:
#        raise ValueError("Don't know how to reallocate variable %s" % str(name))


class AbstractInput(six.with_metaclass(abc.ABCMeta, MutableMapping, object)):
    """
    Abstract class defining the methods that must be implemented by Input objects.
    """

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

    def write(self, filepath="run.abi"):
        """
        Write the input file to file to ``filepath``.
        """
        dirname = os.path.dirname(os.path.abspath(filepath))
        if not os.path.exists(dirname): os.makedirs(dirname)

        # Write the input file.
        with open(filepath, "wt") as fh:
            fh.write(str(self))

    def deepcopy(self):
        """Deep copy of the input."""
        return copy.deepcopy(self)

    def set_vars(self, *args, **kwargs):
        """
        Set the value of the variables.
        Return dict with the variables added to the input.

        Example:

            input.set_vars(ecut=10, ionmov=3)
        """
        kwargs.update(dict(*args))
        for varname, varvalue in kwargs.items():
            self[varname] = varvalue
        return kwargs

    def set_vars_ifnotin(self, *args, **kwargs):
        """
        Set the value of the variables but only if the variable is not already present.
        Return dict with the variables added to the input.

        Example:

            input.set_vars(ecut=10, ionmov=3)
        """
        kwargs.update(dict(*args))
        added = {}
        for varname, varvalue in kwargs.items():
            if varname not in self:
                self[varname] = varvalue
                added[varname] = varvalue
        return added

    def pop_vars(self, keys):
        """
        Remove the variables listed in keys.
        Return dictionary with the variables that have been removed.
        Unlike remove_vars, no exception is raised if the variables are not in the input.

        Args:
            keys: string or list of strings with variable names.

        Example:
            inp.pop_vars(["ionmov", "optcell", "ntime", "dilatmx"])
        """
        return self.remove_vars(keys, strict=False)

    def remove_vars(self, keys, strict=True):
        """
        Remove the variables listed in keys.
        Return dictionary with the variables that have been removed.

        Args:
            keys: string or list of strings with variable names.
            strict: If True, KeyError is raised if at least one variable is not present.
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

    def generate(self, **kwargs):
        """
        This function generates new inputs by replacing the variables specified in kwargs.

        Args:
            kwargs: keyword arguments with the values used for each variable.

        .. code-block:: python

            gs_inp = call_function_to_generate_initial_template()

            # To generate two input files with different values of ecut:
            for inp_ecut in gs_inp.generate(ecut=[10, 20]):
                print("do something with inp_ecut %s" % inp_ecut)

            # To generate four input files with all the possible combinations of ecut and nsppol:
            for inp_ecut in gs_inpt.generate(ecut=[10, 20], nsppol=[1, 2]):
                print("do something with inp_ecut %s" % inp_ecut)
        """
        for new_vars in product_dict(kwargs):
            new_inp = self.deepcopy()
            # Remove the variable names to avoid annoying warnings if the variable is overwritten.
            new_inp.remove_vars(new_vars.keys())
            new_inp.set_vars(**new_vars)
            yield new_inp


class AbiAbstractInput(AbstractInput):
    """
    Abstract class defining the methods that must be implemented by Input objects.
    associated to Abinit executables.
    """

    def add_abiobjects(self, *abi_objects):
        """
        This function receive a list of ``AbiVarable`` objects and add
        the corresponding variables to the input.
        """
        d = {}
        for aobj in abi_objects:
            if not hasattr(aobj, "to_abivars"):
                raise TypeError("type %s: %s does not have `to_abivars` method" % (type(aobj), repr(aobj)))
            d.update(self.set_vars(aobj.to_abivars()))
        return d

    @abc.abstractmethod
    def abivalidate(self, workdir=None, manager=None):
        """
        This method should invoke the executable associated to the input object.
        to test whether the input variables are correct and consistent.
        The executable is supposed to implemente some sort of `--dry-run` option
        that invokes the parser to validate the input and exits.

        Args:
            workdir: Working directory of the fake task used to compute the ibz. Use None for temporary dir.
            manager: |TaskManager| of the task. If None, the manager is initialized from the config file.

        Return:
            ``namedtuple`` with the following attributes:

                retcode: Return code. 0 if OK.
                output_file: output file of the run.
                log_file:  log file of the Abinit run, use log_file.read() to access its content.
                stderr_file: stderr file of the Abinit run. use stderr_file.read() to access its content.
                task: Task object
        """


class AbinitInputError(Exception):
    """Base error class for exceptions raised by ``AbinitInput``."""


class AbinitInput(six.with_metaclass(abc.ABCMeta, AbiAbstractInput, MSONable, Has_Structure, object)):
    """
    This object stores the ABINIT variables for a single dataset.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: AbinitInput
    """
    Error = AbinitInputError

    def __init__(self, structure, pseudos, pseudo_dir=None, comment=None, decorators=None, abi_args=None,
                 abi_kwargs=None, tags=None):
        """
        Args:
            structure: Parameters defining the crystalline structure. Accepts |Structure| object
            file with structure (CIF, netcdf file, ...) or dictionary with ABINIT geo variables.
            pseudos: Pseudopotentials to be used for the calculation. Accepts: string or list of strings
                with the name of the pseudopotential files, list of |Pseudo| objects
                or |PseudoTable| object.
            pseudo_dir: Name of the directory where the pseudopotential files are located.
            ndtset: Number of datasets.
            comment: Optional string with a comment that will be placed at the beginning of the file.
            decorators: List of `AbinitInputDecorator` objects.
            abi_args: list of tuples (key, value) with the initial set of variables. Default: Empty
            abi_kwargs: Dictionary with the initial set of variables. Default: Empty
            tags: list/set of tags describing the input
        """
        self._spell_check = True

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
            if not os.path.exists(pseudo_dir): raise self.Error("Directory %s does not exist" % pseudo_dir)
            pseudos = [os.path.join(pseudo_dir, p) for p in list_strings(pseudos)]

        try:
            self._pseudos = PseudoTable.as_table(pseudos).get_pseudos_for_structure(self.structure)
        except ValueError as exc:
            raise self.Error(str(exc))

        if comment is not None: self.set_comment(comment)

        self._decorators = [] if not decorators else decorators[:]
        self.tags = set() if not tags else set(tags)

    def variable_checksum(self):
        """
        Return string with sha1 value in hexadecimal format.
        This method is mainly used in unit tests to check the invariance
        of the input objects. Note, indeed, that AbintInput is mutable and therefore
        should not be used as keyword in dictionaries.
        """
        # Use sha1 from hashlib because python builtin hash is not deterministic
        # (hash is version- and machine-dependent)
        import hashlib
        sha1 = hashlib.sha1()

        try:
            tos = unicode
        except NameError:
            # Py3K
            def tos(s):
                return str(s).encode(encoding="utf-8")

        # Add key, values to sha1
        # (not sure this is code is portable: roundoff errors and conversion to string)
        # We could just compute the hash from the keys (hash equality does not necessarily imply __eq__!)
        for key in sorted(self.keys()):
            value = self[key]
            if isinstance(value, np.ndarray): value = value.tolist()
            sha1.update(tos(key))
            sha1.update(tos(value))

        # Use string representation to compute hash
        # Not perfect but it supposed to be better than the version above
        # Use alphabetical sorting, don't write pseudos (treated below).
        #s = self.to_string(sortmode="a", with_mnemonics=False, with_structure=True, with_pseudos=False)
        #sha1.update(tos(s))

        sha1.update(tos(self.comment))
        # add pseudos (this is easy because we have md5)
        sha1.update(tos([p.md5 for p in self.pseudos]))
        # add the decorators, do we need to add them ?
        sha1.update(tos([dec.__class__.__name__ for dec in self.decorators]))

        return sha1.hexdigest()

    @pmg_serialize
    def as_dict(self):
        """
        JSON interface used in pymatgen for easier serialization.
        """
        #vars = OrderedDict()
        # Use a list of (key, value) to serialize the OrderedDict
        abi_args = []
        for key, value in self.items():
            if isinstance(value, np.ndarray): value = value.tolist()
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
        """
        JSON interface used in pymatgen for easier serialization.
        """
        pseudos = [Pseudo.from_file(p['filepath']) for p in d['pseudos']]
        dec = MontyDecoder()
        return cls(d["structure"], pseudos, decorators=dec.process_decoded(d["decorators"]),
                   comment=d["comment"], abi_args=d["abi_args"], tags=d["tags"])

    def __setitem__(self, key, value):
        if key in _TOLVARS_SCF and hasattr(self, '_vars') and any(t in self._vars and t != key for t in _TOLVARS_SCF):
            logger.info("Replacing previously set tolerance variable: {0}."
                        .format(self.remove_vars(_TOLVARS_SCF, strict=False)))

        return super(AbinitInput, self).__setitem__(key, value)

    def _check_varname(self, key):
        if key in GEOVARS:
            raise self.Error("You cannot set the value of a variable associated to the structure.\n"
                             "Use Structure objects to prepare the input file.")

        if self.spell_check and not is_abivar(key):
            raise self.Error("""
Cannot find variable `%s` in internal database.
If you think this is not a typo, use:

    input.set_spell_check(False)

to disable spell checking. Perhaps the internal database is not in synch
with the Abinit version you are using. Please contact the AbiPy developers.""" % key)

    #def __eq__(self, other)
    #def __ne__(self, other)
    #    return not self.__eq__(other)

    @property
    def runlevel(self):
        """
        A set of strings defining the type of run of the current input.
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
            gwcalctyp = self.get("gwcalctyp")
            if int(gwcalctyp) > 100:
                runlevel.add(HYBRID)
            else:
                runlevel.update([MANY_BODY, SIGMA])
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

    @property
    def uses_ktimereversal(self):
        """
        True if time-reversal symmetry is used to generate k-points in the IBZ.
        """
        return has_timrev_from_kptopt(self.get("kptopt", 1))

    def set_spell_check(self, false_or_true):
        """Activate/Deactivate spell-checking"""
        self._spell_check = bool(false_or_true)

    @property
    def spell_check(self):
        """True if spell checking is activated."""
        try:
            return self._spell_check
        except AttributeError:  # TODO: This is to maintain compatibility with pickle
            return False

    def to_string(self, sortmode="section", post=None, with_mnemonics=False, mode="text",
                  with_structure=True, with_pseudos=True, exclude=None, verbose=0):
        r"""
        String representation.

        Args:
            sortmode: "section" if variables should be grouped by sections.
                "a" for alphabetical order, None if no sorting is wanted.
            with_mnemonics: True if mnemonics should be added.
            mode: Either `text` or `html` if HTML output with links is wanted.
            post: String that will be appended to the name of the variables
                Note that post is usually autodetected when we have multiple datatasets
                It is mainly used when we have an input file with a single dataset
                so that we can prevent the code from adding "1" to the name of the variables
                (In this case, indeed, Abinit complains if ndtset=1 is not specified
                and we don't want ndtset=1 simply because the code will start to add
                _DS1_ to all the input and output files.
            with_structure: False if section with structure variables should not be printed.
            with_pseudos: False if JSON section with pseudo data should not be added.
            exclude: List of variable names that should be ignored.
        """
        if mode == "html":
            if six.PY2:
                import cgi
                def escape(text):
                    return cgi.escape(text, quote=True)
            else:
                import html
                def escape(text):
                    return html.escape(text, quote=True)
        else:
            def escape(text):
                return text

        lines = []
        app = lines.append

        if self.comment: app("# " + self.comment.replace("\n", "\n#"))

        post = post if post is not None else ""
        mnemonics = self.mnemonics
        if with_mnemonics: mnemonics = with_mnemonics
        exclude = set(exclude) if exclude is not None else set()

        # If spell checking is deactivated, we cannot use mmemonics or sormode == "section"
        if not self.spell_check:
            mnemonics = False
            sortmode = "a"

        if mnemonics or sortmode == "section":
            var_database = get_abinit_variables()

        if sortmode in (None, "a"):
            # Default is no sorting else alphabetical order.
            keys = [k for k, v in self.items() if k not in exclude and v is not None]
            if sortmode == "a": keys = sorted(keys)

            # Extract the items from the dict and add the geo variables at the end
            items = [(k, self[k]) for k in keys]
            if with_structure:
                items.extend(list(self.structure.to_abivars().items()))

            for name, value in items:
                if mnemonics and value is not None:
                    app("# <" + var_database[name].mnemonics + ">")

                # Build variable, convert to string and append it
                vname = name + post
                if mode == "html": vname = var_database[name].html_link(label=vname)
                app(str(InputVariable(vname, value)))

        elif sortmode == "section":
            # Group variables by section.
            # Get dict mapping section_name --> list of variable names belonging to the section.
            keys = [k for (k, v) in self.items() if k not in exclude and v is not None]
            sec2names = var_database.group_by_varset(keys)
            w = 92

            for sec, names in sec2names.items():
                app(w * "#")
                app("#" + ("SECTION: %s" % sec).center(w - 1))
                app(w * "#")
                for name in names:
                    value = self[name]
                    if mnemonics and value is not None:
                        app(escape("# <" + var_database[name].mnemonics + ">"))

                    # Build variable, convert to string and append it
                    vname = name + post
                    if mode == "html": vname = var_database[name].html_link(label=vname)

                    app(str(InputVariable(vname, value)))

            if with_structure:
                app(w * "#")
                app("#" + ("STRUCTURE").center(w - 1))
                app(w * "#")
                for name, value in self.structure.to_abivars().items():
                    if mnemonics and value is not None:
                        app(escape("# <" + var_database[name].mnemonics + ">"))
                    vname = name + post
                    if mode == "html": vname = var_database[name].html_link(label=vname)
                    app(str(InputVariable(vname, value)))

        else:
            raise ValueError("Unsupported value for sortmode %s" % str(sortmode))

        s = "\n".join(lines)
        if not with_pseudos:
            return s if mode != "html" else s.replace("\n", "<br>")

        # Add JSON section with pseudo potentials.
        ppinfo = ["\n\n\n#<JSON>"]
        d = {"pseudos": [p.as_dict() for p in self.pseudos]}
        ppinfo.extend(json.dumps(d, indent=4).splitlines())
        ppinfo.append("</JSON>")

        s += escape("\n#".join(ppinfo))
        if mode == "html": s = s.replace("\n", "<br>")
        return s

    def _repr_html_(self):
        """Integration with jupyter notebooks."""
        return self.to_string(sortmode="section", with_mnemonics=False, mode="html",
                              with_structure=True, with_pseudos=False)

    @property
    def comment(self):
        """Optional string with comment. None if comment is not set."""
        try:
            return self._comment
        except AttributeError:
            return None

    def set_comment(self, comment):
        """Set a comment to be included at the top of the file."""
        self._comment = comment

    @property
    def structure(self):
        """The |Structure| object associated to this input."""
        return self._structure

    def set_structure(self, structure):
        """Set structure."""
        self._structure = Structure.as_structure(structure)

        # Check volume
        m = self.structure.lattice.matrix
        if np.dot(np.cross(m[0], m[1]), m[2]) <= 0:
            raise self.Error("The triple product of the lattice vector is negative. Use structure.abi_sanitize.")

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
        shiftk = np.reshape(shiftk, (-1, 3))
        return self.set_vars(ngkpt=ngkpt, kptopt=kptopt, nshiftk=len(shiftk), shiftk=shiftk)

    def set_gamma_sampling(self):
        """Gamma-only sampling of the BZ."""
        return self.set_kmesh(ngkpt=(1, 1, 1), shiftk=(0, 0, 0))

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

    def get_ngkpt_shiftk(self):
        """
        Return info on the k-point sampling from the input file,
        more specifically a tuple with nkgpt and shift.
        ngkpt is set to None if the BZ sampling cannot be described in terms of
        three divisions + one shift.
        """
        # ngkpt and kptrlatt are mutually exclusive.
        nshiftk = self.get("nshiftk", 1)
        shiftk = np.reshape(self.get("shiftk", [0.5, 0.5, 0.5]), (-1, 3))[:nshiftk, :]

        kptrlatt = self.get("kptrlatt", None)
        if kptrlatt is not None:
            kptrlatt = np.reshape(kptrlatt, (3, 3))
            # Check whether is diagonal with one shift.
            ngkpt = None
            if nshiftk == 1 and is_diagonal(kptrlatt):
                ngkpt = np.diag(kptrlatt)
        else:
            ngkpt = np.array(self.get("ngkpt"))

        return ngkpt, shiftk

    def set_phdos_qmesh(self, nqsmall, method="tetra", ph_qshift=(0, 0, 0)):
        """
        Set the variables (ngkpt, shift, kptopt) for the computation of the Phonon DOS in Abinit.
        Remember that the Phdos is computed via Fourier interpolation so there's no constraint
        of the q-mesh.

        Args:
            nqsmall: Number of k-points used to sample the smallest lattice vector.
            method: gaussian or tetra.
            ph_qshift: Shift for the mesh.
        """
        # q-mesh for Fourier interpolatation of IFC and a2F(w)
        ph_ngqpt = self.structure.calc_ngkpt(nqsmall)
        ph_qshift = np.reshape(ph_qshift, (-1, 3))

        # TODO: Test default values of wstep and smear
        ph_intmeth = {"gaussian": 1, "tetra": 2}[method]
        ph_smear = "0.001 eV" if method == "gaussian" else None

        return self.set_vars(
            ph_intmeth=ph_intmeth,
            ph_smear=ph_smear,
            ph_wstep="0.0001 eV",
            ph_ngqpt=ph_ngqpt,
            ph_qshift=ph_qshift,
            ph_nqshift=len(ph_qshift),
        )

    def set_kpath(self, ndivsm, kptbounds=None, iscf=-2):
        """
        Set the variables for the computation of the electronic band structure.

        Args:
            ndivsm: Number of divisions for the smallest segment.
            kptbounds: k-points defining the path in k-space.
                If None, we use the default high-symmetry k-path defined in the pymatgen database.
        """
        if kptbounds is None: kptbounds = self.structure.calc_kptbounds()
        kptbounds = np.reshape(kptbounds, (-1,3))
        #self.pop_vars(["ngkpt", "shiftk"]) ??

        return self.set_vars(kptbounds=kptbounds, kptopt=-(len(kptbounds)-1), ndivsm=ndivsm, iscf=iscf)

    def set_qpath(self, ndivsm, qptbounds=None):
        """
        Set the variables for the computation of the phonon band structure
        and phonon linewidths.

        Args:
            ndivsm: Number of divisions for the smallest segment.
            qptbounds: q-points defining the path in q-space.
                If None, we use the default high-symmetry q-path defined in the pymatgen database.
        """
        if qptbounds is None: qptbounds = self.structure.calc_kptbounds()
        qptbounds = np.reshape(qptbounds, (-1, 3))

        return self.set_vars(ph_ndivsm=ndivsm, ph_nqpath=len(qptbounds), ph_qpath=qptbounds)

    def set_kptgw(self, kptgw, bdgw):
        """
        Set the variables (k-points, bands) for the computation of GW corrections.

        Args:
            kptgw: List of k-points in reduced coordinates.
            bdgw: Specifies the range of bands for the GW corrections.
                Accepts iterable that be reshaped to (nkptgw, 2)
                or a tuple of two integers if the extrema are the same for each k-point.
        """
        kptgw = np.reshape(kptgw, (-1,3))
        nkptgw = len(kptgw)
        if len(bdgw) == 2: bdgw = len(kptgw) * bdgw

        return self.set_vars(kptgw=kptgw, nkptgw=nkptgw, bdgw=np.reshape(bdgw, (nkptgw, 2)))

    def set_spin_mode(self, spin_mode):
        """
        Set the variables used to the treat the spin degree of freedom.
        Return dictionary with the variables that have been removed.

        Args:
            spin_mode: :class:`SpinMode` object or string. Possible values for string are:

            - polarized
            - unpolarized
            - afm (anti-ferromagnetic)
            - spinor (non-collinear magnetism)
            - spinor_nomag (non-collinear, no magnetism)
        """
        # Remove all variables used to treat spin
        old_vars = self.pop_vars(["nsppol", "nspden", "nspinor"])
        self.add_abiobjects(aobj.SpinMode.as_spinmode(spin_mode))
        return old_vars

    def set_autospinat(self, default=0.6):
        """
        Set the variable spinat for collinear calculation in the format (0, 0, m) with the value of m determined
        with the following order of preference:

        1. If the site of the structure has a magmom setting, that is used.
        2. If the species on the site has a spin setting, that is used.
        3. If the species itself has a particular setting in the config file, that
           is used, e.g., Mn3+ may have a different magmom than Mn4+.
        4. The element symbol itself is checked in the config file.
        5. If there are no settings, the default value is used.
        """
        # These magnetic moments are from the Materials Project
        # (MPVaspInputSet.yaml, short_sha1 = a63bcdf)

        magmom_mp_conf = {
            "Co": 5,
            "Co3+": 0.6,
            "Co4+": 1,
            "Cr": 5,
            "Fe": 5,
            "Mn": 5,
            "Mn3+": 4,
            "Mn4+": 3,
            "Mo": 5,
            "Ni": 5,
            "V": 5,
            "W": 5,
            "Ce": 5,
            "Eu": 10,
        }

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
        """List of |Pseudo| objects."""
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
    def valence_electrons_per_atom(self):
        """Number of valence electrons for each atom in the structure."""
        return self.structure.valence_electrons_per_atom(self.pseudos)

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
            num (int): Number of samples to generate. Default is 50.
            endpoint (bool): optional. If True, `stop` is the last sample. Otherwise, it is not included.
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
            raise self.Error("The number of variables must equal the number of lists\n"
                              "varnames: %s\nvalues %s" % (str(varnames), str(values)))

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

    def new_with_vars(self, *args, **kwargs):
        """
        Return a new input with the given variables.

        Example:
            new = input.new_with_vars(ecut=20)
        """
        # Avoid modifications in self.
        new = self.deepcopy()
        new.set_vars(*args, **kwargs)
        return new

    #def new_with_supercell(self, scdims):
    #    sucell = self.structure * scdims
    #    return new_with_structure(sucell, scdims=scdims)

    def new_with_structure(self, new_structure, scdims=None, verbose=1):
        """
        Return a new |AbinitInput| with different structure.
        See notes below for the constraints that must be fulfilled by the new structure

        Args:
            new_structure: Parameters defining the crystalline structure. Accepts |Structure| object
                file with structure (CIF, netcdf file, ...) or dictionary with ABINIT geo variables.
            scdims: 3 integer giving with the number of cells in the supercell along the three reduced directions.
                Must be used when `new_structure` represents a supercell of the initial structure defined
                in the input file.
            verbose: Verbosity level.

        .. warning::

            If ``scdims`` is None (i.e. no supercell), the two structures must have the same value of
            `natom` and `typat`, they can only differ at the level of the lattice and of the atomic positions.
            When structure represents a supercell, `scdims` must be coherent with the `new_structure` passed
            as argument.
        """
        # Check structure
        if scdims is None:
            # Assume same value of natom and typat
            if len(self.structure) != len(new_structure):
                raise ValueError("Structures must have same value of natom")
            errors = []
            for i, (site1, site2) in enumerate(zip(self.structure, new_structure)):
                if site1.specie.symbol != site2.specie.symbol:
                    errors.append("[%d] %s != %s" % (i, site1.specie.symbol, site2.specie.symbol))
            if errors:
                raise ValueError("Structures must have same order of atomic types:\n" + "\n".join(errors))

        else:
            scdims = np.array(scdims)
            if scdims.shape != (3,):
                raise ValueError("Expecting 3 int in scdims but got %s" % str(scdims))

            numcells = np.product(scdims)
            if len(new_structure) != numcells * len(self.structure):
                errmsg = "Number of atoms in the input structure should be %d * %d but found %d" % (
                    numcells, len(self.structure), len(new_structure))
                raise ValueError(errmsg)

            expected_symbols = numcells * [site.specie.symbol for site in self.structure]
            supcell_symbols = [site.specie.symbol for site in new_structure]
            if not np.array_equal(expected_symbols, supcell_symbols):
                msg = ("Wrong supercell. The routine assumes the atoms in the other cells have the\n"
                       "same ordering as the atoms in the original cell.\n"
                       "expected_symbols: %s\nsupercell_symbols: %s" % (str(expected_symbols), str(supcell_symbols))
                       )
                raise ValueError(msg)
            # TODO Check angles and lengths

        # Build new input
        new = AbinitInput(new_structure, self.pseudos, abi_args=list(self.items()),
                          decorators=self.decorators, tags=self.tags)

        if scdims is not None:
            # This is the tricky part because variables whose shape depends on natom
            # must be changed in order to be consistent with the supercell.
            # Here we use the database of abinit variables to find the variables whose shape depends on `natom`.
            # The method raises ValueError if an array that depends on `natom` is found and no handler is implemented.
            # It's better to raise an exception here than having a error when Abinit parses the input file!

            errors = []
            var_database = get_abinit_variables()
            for name in new:
                var = var_database[name]
                #if var.depends_on_dimension("ntypat"):
                #    errors.append("Found variable %s with ntypat in dimensions %s" % (name, str(var.dimensions)))

                if var.depends_on_dimension("natom"):
                    errors.append("Found variable %s with natom in dimensions %s" % (name, str(var.dimensions)))
                    #new[name] = _repeat_array(name, new[name], len(self.structure), numcells)

            if errors:
                errmsg = ("\n".join(errors) +
                          "\nThe present version of new_with_structure is not able to handle this case.")
                raise ValueError(errmsg)

            # Rescale nband and k-point sampling
            iscale = int(np.ceil(len(new.structure) / len(self.structure)))
            if "nband" in new:
                new["nband"] = int(self["nband"] * iscale)
                if verbose: print("self['nband']", self["nband"], "new['nband']", new["nband"])

            if "ngkpt" in new:
                new["ngkpt"] = (np.rint(np.array(new["ngkpt"]) / scdims)).astype(int)
                if verbose: print("new ngkpt:", new["ngkpt"])

            # TODO
            elif "kptrlatt" in new:
                raise NotImplementedError("kptrlatt in new_with_structure")
                #new["kptrlatt"] = (np.rint(np.array(new["kptrlatt"]) / iscale)).astype(int)
            else:
               # Single k-point
               pass

            # Add chkprim if not yet done.
            new.set_vars_ifnotin(chkprim=0)

        return new

    def new_with_decorators(self, decorators):
        """
        This function receives a list of :class:`AbinitInputDecorator` objects or just a single object,
        applies the decorators to the input and returns a new |AbinitInput| object. self is not changed.
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
        Return dictionary with the variables that have been removed.
        """
        return self.remove_vars(_TOLVARS, strict=False)

    def pop_irdvars(self):
        """
        Remove all the `ird*` variables present in self.
        Return dictionary with the variables that have been removed.
        """
        return self.remove_vars(_IRDVARS, strict=False)

    #def pop_relax_vars(self):
    #    """
    #    Remove all the relax variables present in self.
    #    Return dictionary with the variables that have been removed.
    #    """
    #    return self.pop_vars(["ionmov", "optcell", "ntime", "dilatmx"])

    @property
    def scf_tolvar(self):
        """
        Returns the tolerance variable and value relative to the SCF convergence.
        If more than one is present raise an error
        """
        tolvar, value = None, None
        for t in _TOLVARS_SCF:
            if t in self and self[t]:
                if tolvar:
                    raise self.Error('More than one tolerance set.')
                tolvar = t
                value = self[t]

        return tolvar, value

    def make_ph_inputs_qpoint(self, qpt, tolerance=None, prtwf=-1, prepgkk=0, manager=None):
        """
        Builds and returns a |MultiDataset| list of input files
        for the calculation of phonons at the given q-point `qpt`.
        It should be called with an input the represents a GS run.

        Args:
            qpt: q-point in reduced coordinates.
            tolerance: dict {varname: value} with the tolerance to be used in the DFPT run.
                Defaults to {"tolvrs": 1.0e-10}.
            prtwf: 1WFs are only needed for restarting or non-linear response.
                Since these files are huge, we use prtwf -1 so that the 1WF file is produced
                only if the calculation is not converged so that AbiPy can restart it.
            manager: |TaskManager| of the task. If None, the manager is initialized from the config file.

        .. WARNING::

            The routine assumes the q-point is such that k + q belongs to the initial GS mesh.
            so that the DFPT run can be started from the WFK file directly without having
            to generate WFQ files.
        """
        if tolerance is None: tolerance = {"tolvrs": 1.0e-10}

        if len(tolerance) != 1 or any(k not in _TOLVARS for k in tolerance):
            raise self.Error("Invalid tolerance: %s" % str(tolerance))

        # Call Abinit to get the list of irred perts.
        perts = self.abiget_irred_phperts(qpt=qpt, manager=manager, prepgkk=prepgkk)

        # Build list of datasets (one input per perturbation)
        # Remove iscf if any (required if we pass an input for NSCF calculation)
        ph_inputs = MultiDataset.replicate_input(input=self, ndtset=len(perts))
        ph_inputs.pop_vars("iscf")

        # Set kptopt depending on the q-points i.e use time-reversal if Gamma
        kptopt = 3
        if np.allclose(qpt, 0): kptopt = 2

        # Note: this will work for phonons, but not for the other types of perturbations.
        for pert, ph_input in zip(perts, ph_inputs):
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
            #if "prtwf" not in  ph_input: ph_input["prtwf"] = prtwf
            ph_input.pop_tolerances()
            ph_input.set_vars(tolerance)

        return ph_inputs

    def make_ddk_inputs(self, tolerance=None, kptopt=2, only_vk=False, manager=None):
        """
        Return inputs for performing DDK calculations.
        This functions should be called with an input the represents a GS run.

        Args:
            tolerance: dict {varname: value} with the tolerance to be used in the DFPT run.
                Defaults to {"tolwfr": 1.0e-22}.
            kptopt: 2 to take into account time-reversal symmetry. Note that kptopt 1 is not available.
            only_vk: If only matrix elements of the velocity operator are needed.
                First-order wavefunctions won't be converged --> not usable for other DFPT calculations.
            manager: |TaskManager| of the task. If None, the manager is initialized from the config file.

        Return:
            List of |AbinitInput| objects for DFPT runs.
        """
        if tolerance is None: tolerance = {"tolwfr": 1.0e-22}

        if len(tolerance) != 1 or any(k not in _TOLVARS for k in tolerance):
            raise self.Error("Invalid tolerance: %s" % str(tolerance))

        if "tolvrs" in tolerance:
            raise self.Error("tolvrs should not be used in a DDK calculation")

        # Call Abinit to get the list of irred perts.
        #perts = self.abiget_irred_phperts(qpt=qpt)
        # TODO Add symmetries when implemented.
        ddk_rfdirs = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]

        # Build list of datasets (one input per perturbation)
        ddk_inputs = MultiDataset.replicate_input(input=self, ndtset=len(ddk_rfdirs))

        # See tutorespfn/Input/trf1_5.in
        for rfdir, ddk_input in zip(ddk_rfdirs, ddk_inputs):
            ddk_input.set_vars(
                rfelfd=2,             # Activate the calculation of the d/dk perturbation
                                      # only the derivative of ground-state wavefunctions with respect to k
                rfdir=rfdir,          # Direction of the per ddk.
                nqpt=1,               # One wavevector is to be considered
                qpt=(0, 0, 0),        # q-wavevector.
                kptopt=kptopt,        # 2 to take into account time-reversal symmetry.
                iscf=-3,              # The d/dk perturbation must be treated in a non-self-consistent way
            )

            ddk_input.pop_tolerances()
            ddk_input.set_vars(tolerance)

        if only_vk:
            ddk_inputs.set_vars(nstep=1, nline=1)

        return ddk_inputs

    def make_dde_inputs(self, tolerance=None, use_symmetries=True, manager=None):
        """
        Return |MultiDataset| inputs for the calculation of the electric field perturbations.
        This functions should be called with an input the represents a gs run.

        Args:
            tolerance: dict {varname: value} with the tolerance to be used in the DFPT run.
                Defaults to {"tolvrs": 1.0e-22}.
            use_symmetries: boolean that computes the irreducible components of the perturbation.
                Default to True. Should be set to False for nonlinear coefficients calculation.
            manager: |TaskManager| of the task. If None, the manager is initialized from the config file.
        """
        if tolerance is None:
            tolerance = {"tolvrs": 1.0e-22}

        if len(tolerance) != 1 or any(k not in _TOLVARS for k in tolerance):
            raise self.Error("Invalid tolerance: %s" % str(tolerance))

        if use_symmetries:
            # Call Abinit to get the list of irred perts.
            perts = self.abiget_irred_ddeperts(manager=manager)

            # Build list of datasets (one input per irreducible perturbation)
            multi = MultiDataset.replicate_input(input=self, ndtset=len(perts))

            # See tutorespfn/Input/trf1_5.in dataset 3
            for pert, inp in zip(perts, multi):
                rfdir = 3 * [0]
                rfdir[pert.idir - 1] = 1

                inp.set_vars(
                    rfdir=rfdir,  # Direction of the dde perturbation.
                )

        else:
            # Compute all the directions of the perturbation
            dde_rfdirs = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]

            # Build list of datasets (one input per perturbation)
            multi = MultiDataset.replicate_input(input=self, ndtset=len(dde_rfdirs))

            # See tutorespfn/Input/tnlo_2.in dataset 4
            for rfdir, inp in zip(dde_rfdirs, multi):
                inp.set_vars(
                    rfdir=rfdir,  # Direction of the per ddk.
                    prepanl=1,    # Prepare Non-linear RF calculations.
                )

        multi.set_vars(
            rfelfd=3,       # Activate the calculation of the electric field perturbation
                            # Assuming the data on derivative of ground-state wavefunction with respect
                            # to k (DDK) is available on disk and will be read with getddk/irddk
            nqpt=1,         # One wavevector is to be considered
            qpt=(0, 0, 0),  # q-wavevector.
            kptopt=2,       # Take into account time-reversal symmetry.
        )

        multi.pop_tolerances()
        multi.set_vars(tolerance)

        return multi

    def make_dte_inputs(self, phonon_pert=False, skip_permutations=False, manager=None):
        """
        Return |MultiDataset| inputs for DTE calculation.
        This functions should be called with an input that represents a GS run.

        Args:
            phonon_pert: is True also the phonon perturbations will be considered. Default False.
            skip_permutations: Since the current version of abinit always performs all the permutations
                of the perturbations, even if only one is asked, if True avoids the creation of inputs that
                will produce duplicated outputs.
            manager: |TaskManager| of the task. If None, the manager is initialized from the config file.
        """
        # Call Abinit to get the list of irred perts.
        perts = self.abiget_irred_dteperts(phonon_pert=phonon_pert, manager=manager)

        if skip_permutations:
            perts_to_skip = []
            reduced_perts = []
            for pert in perts:
                p = ((pert.i1pert, pert.i1dir), (pert.i2pert, pert.i2dir), (pert.i3pert, pert.i3dir))
                if p not in perts_to_skip:
                    reduced_perts.append(pert)
                    perts_to_skip.extend(itertools.permutations(p))

            perts = reduced_perts

        # Build list of datasets (one input per perturbation)
        multi = MultiDataset.replicate_input(input=self, ndtset=len(perts))

        # See tutorespfn/Input/tnlo_2.in
        na = len(self.structure)

        for pert, inp in zip(perts, multi):
            rfdir1 = 3 * [0]
            rfdir1[pert.i1dir - 1] = 1
            rfdir2 = 3 * [0]
            rfdir2[pert.i2dir - 1] = 1
            rfdir3 = 3 * [0]
            rfdir3[pert.i3dir - 1] = 1

            # atpol if needed. Since there can be only one spatial perturbation
            m = min(pert.i1pert, pert.i2pert, pert.i3pert)
            atpol = [m, m] if m <= na else None

            inp.set_vars(
                # Activate the calculation of the electric field perturbation
                d3e_pert1_elfd=1 if pert.i1pert == na+2 else 0,
                d3e_pert2_elfd=1 if pert.i2pert == na+2 else 0,
                d3e_pert3_elfd=1 if pert.i3pert == na+2 else 0,
                d3e_pert1_dir=rfdir1,  # Direction of the dte perturbation.
                d3e_pert2_dir=rfdir2,
                d3e_pert3_dir=rfdir3,
                d3e_pert1_phon = 1 if pert.i1pert <= na else 0,
                d3e_pert2_phon = 1 if pert.i2pert <= na else 0,
                d3e_pert3_phon = 1 if pert.i3pert <= na else 0,
                d3e_pert1_atpol = atpol,
                nqpt=1,         # One wavevector is to be considered
                qpt=(0, 0, 0),  # q-wavevector.
                optdriver=5,    # non-linear response functions, using the 2n+1 theorem.
                kptopt=2,       # Take into account time-reversal symmetry.
            )

            inp.pop_tolerances()

        return multi

    def make_bec_inputs(self, tolerance=None, manager=None):
        """
        Return |MultiDataset| inputs for the calculation of the Born effective charges.
        This functions should be called with an input that represents a GS run.
        """
        if tolerance is None: tolerance = {"tolvrs": 1.0e-10}

        if len(tolerance) != 1 or any(k not in _TOLVARS for k in tolerance):
            raise self.Error("Invalid tolerance: %s" % str(tolerance))

        # Call Abinit to get the list of irred perts.
        # TODO:
        # Check that one can use the same list of irred perts as in phonons
        perts = self.abiget_irred_phperts(qpt=(0, 0, 0), manager=manager)

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

    def make_strain_perts_inputs(self, tolerance=None, phonon_pert=True, kptopt=2, manager=None):
        """
        Return |MultiDataset| inputs for strain perturbation calculation.
        This functions should be called with an input that represents a GS run.

        Args:
            tolerance: dict {varname: value} with the tolerance to be used in the DFPT run.
                Defaults to {"tolvrs": 1.0e-12}.
            manager: |TaskManager| of the task. If None, the manager is initialized from the config file.
            phonon_pert: is True also the phonon perturbations will be considered. Default False.
            kptopt: 2 to take into account time-reversal symmetry.
        """
        if tolerance is None:
            tolerance = {"tolvrs": 1.0e-12}

        if len(tolerance) != 1 or any(k not in _TOLVARS for k in tolerance):
            raise self.Error("Invalid tolerance: {}".format(str(tolerance)))

        perts = self.abiget_irred_strainperts(kptopt=kptopt, manager=manager, phonon_pert=phonon_pert)
        #print("Stress perts:", perts)

        # Build list of datasets (one input per perturbation)
        multi = MultiDataset.replicate_input(input=self, ndtset=len(perts))

        for pert, inp in zip(perts, multi):
            rfdir = 3 * [0]
            rfdir[pert.idir -1] = 1

            if pert.ipert <= len(self.structure):
                inp.set_vars(rfphon=1,             # Activate the calculation of the atomic dispacement perturbations
                             rfatpol=[pert.ipert, pert.ipert],
                             rfdir=rfdir,
                             nqpt=1,               # One wavevector is to be considered
                             qpt=(0, 0, 0),        # q-wavevector.
                             kptopt=kptopt,        # No symmetries
                             #iscf=7,
                             paral_kgb=0
                             )

            elif pert.ipert == len(self.structure) + 3:
                inp.set_vars(rfstrs=1,             # Activate the calculation of the strain perturbations (uniaxial)
                             rfdir=rfdir,
                             nqpt=1,               # One wavevector is to be considered
                             qpt=(0, 0, 0),        # q-wavevector.
                             kptopt=kptopt,        # No symmetries
                             #iscf=7,
                             paral_kgb=0
                             )

            elif pert.ipert == len(self.structure) + 4:
                inp.set_vars(rfstrs=2,             # Activate the calculation of the strain perturbations (shear)
                             rfdir=rfdir,
                             nqpt=1,               # One wavevector is to be considered
                             qpt=(0, 0, 0),        # q-wavevector.
                             kptopt=kptopt,        # No symmetries
                             #iscf=7,
                             paral_kgb=0
                             )

            inp.pop_tolerances()
            inp.set_vars(tolerance)

        return multi

    def abivalidate(self, workdir=None, manager=None):
        """
        Run ABINIT in dry-run mode to validate the input file.

        Args:
            workdir: Working directory of the fake task used to compute the ibz. Use None for temporary dir.
            manager: |TaskManager| of the task. If None, the manager is initialized from the config file.

        Return:
            `namedtuple` with the following attributes:

                retcode: Return code. 0 if OK.
                output_file: output file of the run.
                log_file:  log file of the Abinit run, use log_file.read() to access its content.
                stderr_file: stderr file of the Abinit run. use stderr_file.read() to access its content.
                task: Task object
        """
        task = AbinitTask.temp_shell_task(inp=self, workdir=workdir, manager=manager)
        retcode = task.start_and_wait(autoparal=False, exec_args=["--dry-run"])
        return dict2namedtuple(retcode=retcode, output_file=task.output_file, log_file=task.log_file,
                               stderr_file=task.stderr_file, task=task)

    def abiget_spacegroup(self, tolsym=None, retdict=False, workdir=None, manager=None, verbose=0):
        """
        This function invokes Abinit to get the space group (as detected by Abinit, not by spglib)
        It should be called with an input file that contains all the mandatory variables required by ABINIT.

        Args:
            tolsym: Abinit tolsym input variable. None correspondes to the default value.
            retdict: True to return dictionary with space group information instead of Structure.
            workdir: Working directory of the fake task used to compute the ibz. Use None for temporary dir.
            manager: |TaskManager| of the task. If None, the manager is initialized from the config file.
	    verbose: Verbosity level.

        Return:
            |Structure| object with AbinitSpaceGroup obtained from the main output file if retdict is False
	    else dict with e.g. {'bravais': 'Bravais cF (face-center cubic)', 'spg_number': 227, 'spg_symbol': 'Fd-3m'}.
        """
        # Avoid modifications in self.
        inp = self.deepcopy()
        if tolsym is not None: inp["tolsym"] = float(tolsym)

        # Bypass Abinit check as we always want to return results.
        inp["chksymbreak"] = 0
        # Disable memory check.
        inp["mem_test"] = 0

        # Build a Task to run Abinit in --dry-run mode.
        task = AbinitTask.temp_shell_task(inp, workdir=workdir, manager=manager)
        task.start_and_wait(autoparal=False, exec_args=["--dry-run"])

        # Parse the output file and return structure extracted from run.abo
        from abipy.abio.outputs import AbinitOutputFile
        try:
            with AbinitOutputFile(task.output_file.path) as out:
                if not retdict:
                    return out.initial_structure
                else:
                    dims_dataset, spginfo_dataset = out.get_dims_spginfo_dataset(verbose=verbose)
                    return spginfo_dataset[1]

        except Exception as exc:
            self._handle_task_exception(task, exc)

    def abiget_ibz(self, ngkpt=None, shiftk=None, kptopt=None, workdir=None, manager=None, verbose=0):
        """
        This function computes the list of points in the IBZ and the corresponding weights.
        It should be called with an input file that contains all the mandatory variables required by ABINIT.

        Args:
            ngkpt: Number of divisions for the k-mesh (default None i.e. use ngkpt from self)
            shiftk: List of shifts (default None i.e. use shiftk from self)
            kptopt: Option for k-point generation. If None, the value in self is used.
            workdir: Working directory of the fake task used to compute the ibz. Use None for temporary dir.
            manager: |TaskManager| of the task. If None, the manager is initialized from the config file.
            verbose: verbosity level.

        Returns:
            `namedtuple` with attributes:
                points: |numpy-array| with points in the IBZ in reduced coordinates.
                weights: |numpy-array| with weights of the points.
        """
        # Avoid modifications in self.
        inp = self.deepcopy()

        # The magic value that makes ABINIT print the ibz and then stop.
        inp["prtkpt"] = -2
        # Bypass Abinit check as we always want to return results.
        inp["chksymbreak"] = 0
        # Disable memory check.
        inp["mem_test"] = 0

        if ngkpt is not None: inp["ngkpt"] = ngkpt
        if shiftk is not None:
            shiftk = np.reshape(shiftk, (-1, 3))
            inp.set_vars(shiftk=shiftk, nshiftk=len(shiftk))

        if kptopt is not None: inp["kptopt"] = kptopt
        if verbose:
            print("Computing ibz with input:\n", str(inp))

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
            self._handle_task_exception(task, exc)

    def _handle_task_exception(self, task, prev_exc):
        """
        This method is called when we have executed a temporary task but we encounter
        an exception when we try to extract data from the output results produced by Abinit
        It tries to extract information about the error and finally raises self.Error.

        .. example::

            try:

                do_something_with_the_output_files_produced_by_the_task

            except Exception as exc:

                self._handle_task_exception(task, exc)
        """
        # Check if there are errors in the log file.
        report = task.get_event_report()
        if report and report.errors:
            raise self.Error(str(report))

        # Weird condition. Possible explanations:
        # 1) Abinit cannot be executed or runtime errors due e.g to libraries
        # 2) IO buffering (Abinit called MPI_ABORT but files are not flushed before aborting.
        # Try to return as much iformation as possible to aid debugging
        errors = ["Problem in temp Task executed in %s" %  task.workdir,
                  "Previous exception %s" % prev_exc]

        try:
            errors.append("Last 50 line from %s:" % str(task.log_file.path))
            log_lines = task.log_file.readlines()
            i = len(log_lines) - 50 if len(log_lines) >= 50 else 0
            errors.extend(s.strip() for s in log_lines[i:])
        except Exception as exc:
            errors.append(str(exc))

        emsg = "\n".join(errors)

        try:
            # TODO: in principle task.debug() but I have to change pymatgen.io.abinit
            task.flow.debug()
        finally:
            raise self.Error(emsg)

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
            manager: |TaskManager| of the task. If None, the manager is initialized from the config file.

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

        # Bypass Abinit check as we always want to return results.
        inp["chksymbreak"] = 0

        if ngkpt is not None: inp["ngkpt"] = ngkpt
        if shiftk is not None:
            shiftk = np.reshape(shiftk, (-1,3))
            inp.set_vars(shiftk=shiftk, nshiftk=len(inp['shiftk']))

        inp.set_vars(
            nqpt=1,       # One wavevector is to be considered
            qpt=qpt,      # q-wavevector.
            paral_rf=-1,  # Magic value to get the list of irreducible perturbations for this q-point.
            **perts_vars
        )

        if kptopt is not None: inp["kptopt"] = kptopt
        #print("Computing irred_perts with input:\n", str(inp))

        # Build a Task to run Abinit in a shell subprocess
        task = AbinitTask.temp_shell_task(inp, workdir=workdir, manager=manager)
        task.start_and_wait(autoparal=False)

        # Parse the file to get the perturbations.
        try:
            return yaml_read_irred_perts(task.log_file.path)
        except Exception as exc:
            # Sometimes the previous call raises: Cannot find next YAML document in /tmp/tmpskvdr_bo/run.log
            # perhaps because the log file is still being written (?) so let's wait a bit.
            time.sleep(5.0)
            try:
                return yaml_read_irred_perts(task.log_file.path)
            except Exception as exc:
                self._handle_task_exception(task, exc)

    def abiget_irred_phperts(self, qpt=None, ngkpt=None, shiftk=None, kptopt=None, prepgkk=0, workdir=None, manager=None):
        """
        This function, computes the list of irreducible perturbations for DFPT.
        It should be called with an input file that contains all the mandatory variables required by ABINIT.

        Args:
            qpt: qpoint of the phonon in reduced coordinates. Used to shift the k-mesh
                if qpt is not passed, self must already contain "qpt" otherwise an exception is raised.
            ngkpt: Number of divisions for the k-mesh (default None i.e. use ngkpt from self)
            shiftk: Shiftks (default None i.e. use shiftk from self)
            kptopt: Option for k-point generation. If None, the value in self is used.
            prepgkk: 1 to activate computation of all 3*natom perts (debugging option).
            workdir: Working directory of the fake task used to compute the ibz. Use None for temporary dir.
            manager: |TaskManager| of the task. If None, the manager is initialized from the config file.

        Returns:
            List of dictionaries with the Abinit variables defining the irreducible perturbation

        Example:

                [{'idir': 1, 'ipert': 1, 'qpt': [0.25, 0.0, 0.0]},
                 {'idir': 2, 'ipert': 1, 'qpt': [0.25, 0.0, 0.0]}]
        """
        phperts_vars = dict(rfphon=1,                         # Will consider phonon-type perturbation
                            rfatpol=[1, len(self.structure)], # Set of atoms to displace.
                            rfdir=[1, 1, 1],                  # Along this set of reduced coordinate axis.
                            prepgkk=prepgkk,
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
            manager: |TaskManager| of the task. If None, the manager is initialized from the config file.

        Returns:
            List of dictionaries with the Abinit variables defining the irreducible perturbation

        Example:

            [{'idir': 1, 'ipert': 4, 'qpt': [0.0, 0.0, 0.0]},
             {'idir': 2, 'ipert': 4, 'qpt': [0.0, 0.0, 0.0]}]

        """
        ddeperts_vars = dict(rfphon=0,  # No phonon-type perturbation
                             rfelfd=3,  # Electric field
                             kptopt=2,  # kpt time reversal symmetry
                             )

        return self._abiget_irred_perts(ddeperts_vars, qpt=(0, 0, 0), ngkpt=ngkpt, shiftk=shiftk, kptopt=kptopt,
                                        workdir=workdir, manager=manager)

    def abiget_irred_dteperts(self, ngkpt=None, shiftk=None, kptopt=None, workdir=None, manager=None,
                              phonon_pert=False):
        """
        This function, computes the list of irreducible perturbations for DFPT.
        It should be called with an input file that contains all the mandatory variables required by ABINIT.

        Args:
            ngkpt: Number of divisions for the k-mesh (default None i.e. use ngkpt from self)
            shiftk: Shiftks (default None i.e. use shiftk from self)
            kptopt: Option for k-point generation. If None, the value in self is used.
            workdir: Working directory of the fake task used to compute the ibz. Use None for temporary dir.
            manager: |TaskManager| of the task. If None, the manager is initialized from the config file.
            phonon_pert: if True also the phonon perturbations will be considered. Default False.

        Returns:
            List of dictionaries with the Abinit variables defining the irreducible perturbation

        Example:

            [{'idir': 1, 'ipert': 4, 'qpt': [0.0, 0.0, 0.0]},
             {'idir': 2, 'ipert': 4, 'qpt': [0.0, 0.0, 0.0]}]

        """
        dteperts_vars = dict(d3e_pert1_phon=1 if phonon_pert else 0,           # phonon-type perturbation
                             d3e_pert2_phon=0,
                             d3e_pert3_phon=0,
                             d3e_pert1_atpol=[1, len(self.structure)] if phonon_pert else None,
                             d3e_pert1_elfd=1,           # Electric field perturbation
                             d3e_pert2_elfd=1,
                             d3e_pert3_elfd=1,
                             d3e_pert1_dir=[1,1,1],
                             d3e_pert2_dir=[1,1,1],
                             d3e_pert3_dir=[1,1,1],
                             optdriver=5,                # non-linear response functions , using the 2n+1 theorem
                             kptopt=2,                   # kpt time reversal symmetry
                             )

        return self._abiget_irred_perts(dteperts_vars, qpt=(0, 0, 0), ngkpt=ngkpt, shiftk=shiftk, kptopt=kptopt,
                                        workdir=workdir, manager=manager)

    def abiget_irred_strainperts(self, ngkpt=None, shiftk=None, kptopt=None, workdir=None, manager=None,
                                 phonon_pert=True):
        """
        This function, computes the list of irreducible perturbations for strain perturbations in DFPT.
        It should be called with an input file that contains all the mandatory variables required by ABINIT.

        Args:
            ngkpt: Number of divisions for the k-mesh (default None i.e. use ngkpt from self)
            shiftk: Shiftks (default None i.e. use shiftk from self)
            kptopt: Option for k-point generation. If None, the value in self is used.
            workdir: Working directory of the fake task used to compute the ibz. Use None for temporary dir.
            manager: |TaskManager| of the task. If None, the manager is initialized from the config file.
            phonon_pert: if True the phonon perturbation at gamma will be included.

        Returns:
            List of dictionaries with the Abinit variables defining the irreducible perturbation

        Example:

            [{'idir': 1, 'ipert': 4, 'qpt': [0.0, 0.0, 0.0]},
             {'idir': 2, 'ipert': 4, 'qpt': [0.0, 0.0, 0.0]}]
        """
        strainperts_vars = dict(rfstrs=3,                        # Do the strain perturbations
                                rfdir=(1,1,1),                   # All directions
                                # nqpt=1,                        # One wavevector is to be considered
                                # qpt=(0, 0, 0),                 # q-wavevector.
                                kptopt=kptopt,                   # Take into account time-reversal symmetry.
                                iscf=7                           # Just so that it works with PAW ... #TODO: check this
                             )

        if phonon_pert:
            strainperts_vars['rfphon'] = 1                        # No phonon-type perturbation
            strainperts_vars['rfatpol'] = (1,len(self.structure)) # Perturbation of all atoms

        return self._abiget_irred_perts(strainperts_vars, qpt=(0, 0, 0), ngkpt=ngkpt, shiftk=shiftk,
                                        kptopt=kptopt,
                                        workdir=workdir, manager=manager)

    def pop_par_vars(self, all=False):
        """
        Remove all the variables associated to parallelism from the input file.
        Useful in case of a restart when we need to remove the parallel variables before rerunning autoparal
        """
        parvars = ['npkpt', 'npfft', 'npband', 'npspinor', 'npimage']
        if all:
            parvars.append('gwpara')
        popped = {}
        for var in parvars:
            popped[var] = self.pop(var, None)

        return popped

    def abiget_autoparal_pconfs(self, max_ncpus, autoparal=1, workdir=None, manager=None, verbose=0):
        """
        Get all the possible configurations up to ``max_ncpus``.
        Return list of parallel configurations.
        """
        inp = self.deepcopy()
        inp.set_vars(autoparal=autoparal, max_ncpus=max_ncpus)

        # Bypass Abinit check as we always want to return results.
        inp["chksymbreak"] = 0
        # Disable memory check.
        inp["mem_test"] = 0

        # Run the job in a shell subprocess with mpi_procs = 1
        # Return code is always != 0
        task = AbinitTask.temp_shell_task(inp, workdir=workdir, manager=manager)
        if verbose: print("Running in:", task.workdir)
        task.start_and_wait(autoparal=False)

        ##############################################################
        # Parse the autoparal configurations from the main output file
        ##############################################################
        parser = ParalHintsParser()
        try:
            pconfs = parser.parse(task.output_file.path)
            return pconfs
        except parser.Error as exc:
            self._handle_task_exception(task, exc)

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
    This object is essentially a list of |AbinitInput| objects.
    that provides an easy-to-use interface to apply global changes to the
    the inputs stored in the objects.

    Let's assume for example that multi contains two ``AbinitInput`` objects and we
    want to set `ecut` to 1 in both dictionaries. The direct approach would be:

        for inp in multi:
            inp.set_vars(ecut=1)

    or alternatively:

        for i in range(multi.ndtset):
            multi[i].set_vars(ecut=1)

    MultiDataset provides its own implementaion of __getattr__ so that one can simply use:

         multi.set_vars(ecut=1)

        multi.get("ecut") returns a list of values. It's equivalent to:

            [inp["ecut"] for inp in multi]

        Note that if "ecut" is not present in one of the input of multi, the corresponding entry is set to None.
        A default value can be specified with:

            multi.get("paral_kgb", 0)

    .. warning::

        MultiDataset does not support calculations done with different sets of pseudopotentials.
        The inputs can have different crystalline structures (as long as the atom types are equal)
        but each input in MultiDataset must have the same set of pseudopotentials.
    """
    Error = AbinitInputError

    @classmethod
    def from_inputs(cls, inputs):
        """Build object from a list of |AbinitInput| objects."""
        for inp in inputs:
            if any(p1 != p2 for p1, p2 in zip(inputs[0].pseudos, inp.pseudos)):
                raise ValueError("Pseudos must be consistent when from_inputs is invoked.")

        # Build MultiDataset from input structures and pseudos and add inputs.
        multi = cls(structure=[inp.structure for inp in inputs], pseudos=inputs[0].pseudos, ndtset=len(inputs))

        # Add variables, decorators and tags.
        for inp, new_inp in zip(inputs, multi):
            new_inp.set_vars(**inp)
            new_inp._decorators = inp.decorators
            new_inp.tags = set(inp.tags)

        return multi

    @classmethod
    def replicate_input(cls, input, ndtset):
        """Construct a multidataset with ndtset from the |AbinitInput| input."""
        multi = cls(input.structure, input.pseudos, ndtset=ndtset)

        for inp in multi:
            inp.set_vars({k: v for k, v in input.items()})
            inp.tags = set(input.tags)

        return multi

    def __init__(self, structure, pseudos, pseudo_dir="", ndtset=1):
        """
        Args:
            structure: file with the structure, |Structure| object or dictionary with ABINIT geo variable
                Accepts also list of objects that can be converted to Structure object.
                In this case, however, ndtset must be equal to the length of the list.
            pseudos: String or list of string with the name of the pseudopotential files.
            pseudo_dir: Name of the directory where the pseudopotential files are located.
            ndtset: Number of datasets.
        """
        # Setup of the pseudopotential files.
        if isinstance(pseudos, Pseudo):
            pseudos = [pseudos]

        elif isinstance(pseudos, PseudoTable):
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
        """Pseudopotential objects."""
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
        #m = getattr(self._inputs[0], name)
        _inputs = object.__getattribute__(self, "_inputs")
        m = getattr(_inputs[0], name)
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
        """self + other"""
        if isinstance(other, AbinitInput):
            new_mds = MultiDataset.from_inputs(self)
            new_mds.append(other)
            return new_mds
        elif isinstance(other, MultiDataset):
            new_mds = MultiDataset.from_inputs(self)
            new_mds.extend(other)
            return new_mds
        else:
            raise NotImplementedError("Operation not supported")

    def __radd__(self, other):
        if isinstance(other, AbinitInput):
            new_mds = MultiDataset.from_inputs([other])
            new_mds.extend(self)
        elif isinstance(other, MultiDataset):
            new_mds = MultiDataset.from_inputs(other)
            new_mds.extend(self)
        else:
            raise NotImplementedError("Operation not supported")

    def append(self, abinit_input):
        """Add a |AbinitInput| to the list."""
        assert isinstance(abinit_input, AbinitInput)
        if any(p1 != p2 for p1, p2 in zip(abinit_input.pseudos, abinit_input.pseudos)):
            raise ValueError("Pseudos must be consistent when from_inputs is invoked.")
        self._inputs.append(abinit_input)

    def extend(self, abinit_inputs):
        """Extends self with a list of |AbinitInput| objects."""
        assert all(isinstance(inp, AbinitInput) for inp in abinit_inputs)
        for inp in abinit_inputs:
            if any(p1 != p2 for p1, p2 in zip(self[0].pseudos, inp.pseudos)):
                raise ValueError("Pseudos must be consistent when from_inputs is invoked.")
        self._inputs.extend(abinit_inputs)

    def addnew_from(self, dtindex):
        """Add a new entry in the multidataset by copying the input with index ``dtindex``."""
        self.append(self[dtindex].deepcopy())

    def split_datasets(self):
        """Return list of |AbinitInput| objects.."""
        return self._inputs

    def deepcopy(self):
        """Deep copy of the MultiDataset."""
        return copy.deepcopy(self)

    @property
    def has_same_structures(self):
        """True if all inputs in MultiDataset are equal."""
        return all(self[0].structure == inp.structure for inp in self)

    def __str__(self):
        return self.to_string()

    def to_string(self, mode="text", verbose=0, with_pseudos=True):
        """
        String representation i.e. the input file read by Abinit.

        Args:
            mode: Either ``text`` or ``html`` if HTML output with links is wanted.
            with_pseudos: False if JSON section with pseudo data should not be added.
        """
        if mode == "html":
            var_database = get_abinit_variables()

        if self.ndtset > 1:
            # Multi dataset mode.
            lines = ["ndtset %d" % self.ndtset]

            def has_same_variable(kref, vref, other_inp):
                """True if variable kref is present in other_inp with the same value."""
                if kref not in other_inp: return False
                otherv = other_inp[kref]
                return np.array_equal(vref, otherv)

            # Don't repeat variable that are common to the different datasets.
            # Put them in the `Global Variables` section and exclude these variables in inp.to_string
            global_vars = set()
            for k0, v0 in self[0].items():
                isame = True
                for i in range(1, self.ndtset):
                    isame = has_same_variable(k0, v0, self[i])
                    if not isame:
                        break
                if isame:
                    global_vars.add(k0)
            #print("global_vars vars", global_vars)

            w = 92
            if global_vars:
                lines.append(w * "#")
                lines.append("### Global Variables.")
                lines.append(w * "#")
                for key in global_vars:
                    vname = key if mode == "text" else var_database[key].html_link(label=key)
                    lines.append(str(InputVariable(vname, self[0][key])))

            has_same_structures = self.has_same_structures
            if has_same_structures:
                # Write structure here and disable structure output in input.to_string
                lines.append(w * "#")
                lines.append("#" + ("STRUCTURE").center(w - 1))
                lines.append(w * "#")
                for key, value in self[0].structure.to_abivars().items():
                    vname = key if mode == "text" else var_database[key].html_link(label=key)
                    lines.append(str(InputVariable(vname, value)))

            for i, inp in enumerate(self):
                header = "### DATASET %d ###" % (i + 1)
                is_last = (i==self.ndtset - 1)
                s = inp.to_string(post=str(i + 1), with_pseudos=is_last and with_pseudos, mode=mode,
                                  with_structure=not has_same_structures, exclude=global_vars)
                if s:
                    header = len(header) * "#" + "\n" + header + "\n" + len(header) * "#" + "\n"
                    s = "\n" + header + s + "\n"

                lines.append(s)

            return "\n".join(lines) if mode=="text" else "\n".join(lines).replace("\n", "<br>")

        else:
            # single datasets ==> don't append the dataset index to the variables.
            # this trick is needed because Abinit complains if ndtset is not specified
            # and we have variables that end with the dataset index e.g. acell1
            # We don't want to specify ndtset here since abinit will start to add DS# to
            # the input and output files thus complicating the algorithms we have to use to locate the files.
            return self[0].to_string(mode=mode, with_pseudos=with_pseudos)

    def _repr_html_(self):
        """Integration with jupyter_ notebooks."""
        return self.to_string(mode="html")

    def get_vars_dataframe(self, *varnames):
        """
        Return pandas DataFrame with the value of the variables specified in `varnames`.

        .. example:

            df = multi.get_vars_dataframe("ecut", "ngkpt")
        """
        import pandas as pd
        frames = []
        for i, inp in enumerate(self):
            df = pd.DataFrame([{v: inp.get(v, None) for v in varnames}],
                              index=["dataset %d" % i], columns=varnames)
            frames.append(df)
        return pd.concat(frames)

    def filter_by_tags(self, tags=None, exclude_tags=None):
        """
        Filters the input according to the tags

        Args:
            tags: A single tag or list/tuple/set of tags
            exclude_tags: A single tag or list/tuple/set of tags that should be excluded

        Returns:
            A |MultiDataset| containing the inputs containing all the requested tags.
        """
        if isinstance(tags, (list, tuple, set)):
            tags = set(tags)
        elif not isinstance(tags, set) and tags is not None:
            tags = {tags}

        if isinstance(exclude_tags, (list, tuple, set)):
            exclude_tags = set(exclude_tags)
        elif not isinstance(exclude_tags, set) and exclude_tags is not None:
            exclude_tags = {exclude_tags}

        if tags is not None:
            inputs = [i for i in self if tags.issubset(i.tags)]
        else:
            inputs = self

        if exclude_tags is not None:
            inputs = [i for i in inputs if not exclude_tags.intersection(i.tags)]

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
        """
        Return new |MultiDataset| object in which only the inputs with the given runlevel are selected.
        """
        if isinstance(runlevel, (list, tuple, set)):
            runlevel = set(runlevel)
        elif not isinstance(runlevel, set):
            runlevel = {runlevel}

        inputs = [i for i in self if runlevel.issubset(i.runlevel)]

        return MultiDataset.from_inputs(inputs) if inputs else None

    def write(self, filepath="run.abi"):
        """
        Write ``ndset`` input files to disk. The name of the file
        is constructed from the dataset index e.g. run0.abi
        """
        root, ext = os.path.splitext(filepath)
        for i, inp in enumerate(self):
            p = root + "DS%d" % i + ext
            inp.write(filepath=p)


class AnaddbInputError(Exception):
    """Base error class for exceptions raised by `AnaddbInput`"""


class AnaddbInput(AbiAbstractInput, Has_Structure):
    """
    This object stores the anaddb variables.


    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: AnaddbInput
    """

    Error = AnaddbInputError

    def __init__(self, structure, comment="", anaddb_args=None, anaddb_kwargs=None, spell_check=True):

        """
        Args:
            structure: |Structure| object
            comment: Optional string with a comment that will be placed at the beginning of the file.
            anaddb_args: List of tuples (key, value) with Anaddb input variables (default: empty)
            anaddb_kwargs: Dictionary with Anaddb input variables (default: empty)
            spell_check: False to disable spell checking for input variables.
        """
        self.set_spell_check(spell_check)
        self._structure = structure
        self.comment = "" if comment is None else str(comment)

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

    def set_spell_check(self, false_or_true):
        """Activate/Deactivate spell-checking"""
        self._spell_check = bool(false_or_true)

    @property
    def spell_check(self):
        """True if spell checking is activated."""
        try:
            return self._spell_check
        except AttributeError: # This is to maintain compatibility with pickle
            return False

    def _check_varname(self, key):
        if self.spell_check and not is_anaddb_var(key):
            raise self.Error("""
Cannot find variable `%s` in internal database.
If you think this is not a typo, use:

    input.set_spell_check(False)

to disable spell checking. Perhaps the internal database is not in synch
with the Abinit version you are using. Please contact the AbiPy developers.""" % key)

    @classmethod
    def modes_at_qpoint(cls, structure, qpoint, asr=2, chneut=1, dipdip=1, ifcflag=0, lo_to_splitting=False,
                        directions=None, anaddb_args=None, anaddb_kwargs=None, spell_check=False):
        """
        Build an |AnaddbInput| for the calculation of the phonon frequencies at a given q-point.

        Args:
            structure: |Structure| object
            qpoint: Reduced coordinates of the q-point where phonon frequencies and modes are wanted
            asr, chneut, dipdp, ifcflag: Anaddb input variable. See official documentation.
            lo_to_splitting: if True calculation of the LO-TO splitting will be included if qpoint==Gamma
            directions: list of 3D directions along which the LO-TO splitting will be calculated. If None the three
                cartesian direction will be used
            anaddb_args: List of tuples (key, value) with Anaddb input variables (default: empty)
            anaddb_kwargs: Dictionary with Anaddb input variables (default: empty)
            spell_check: False to disable spell checking for input variables.
        """
        new = cls(structure, comment="ANADDB input for phonon frequencies at one q-point",
                  anaddb_args=anaddb_args, anaddb_kwargs=anaddb_kwargs)

        # We need a numpy array.
        qpoint = qpoint.frac_coords if hasattr(qpoint, "frac_coords") else np.array(qpoint)

        if len(qpoint) != 3:
            raise ValueError("Wrong q-point %s" % qpoint)

        new.set_vars(
            ifcflag=ifcflag,        # Interatomic force constant flag
            asr=asr,                # Acoustic Sum Rule
            chneut=chneut,          # Charge neutrality requirement for effective charges.
            dipdip=dipdip,          # Dipole-dipole interaction treatment
            # This part is fixed
            nph1l=1,
            qph1l=np.append(qpoint, 1)
        )

        if lo_to_splitting and np.allclose(qpoint, [0, 0, 0]):
            if directions is None:
                directions = [1, 0, 0, 0, 1, 0, 0, 0, 1]
            directions = np.reshape(directions, (-1, 3))
            # append 0 to specify that these are directions,
            directions = np.c_[directions, np.zeros(len(directions))]
            # add
            new.set_vars(
                nph2l=len(directions),
                qph2l=directions
            )

        return new

    @classmethod
    def piezo_elastic(cls, structure, relaxed_ion=True, stress_correction=False,
                      asr=2, chneut=1, dipdip=1, anaddb_args=None, anaddb_kwargs=None):
        """
        Build an |AnaddbInput| for the calculation of piezoelectric and elastic tensor calculations.

        Args:
            asr, chneut, dipdp: Anaddb input variable. See official documentation.
        """
        comment = "ANADDB input for piezoelectric and elastic tensor calculation"

        new = cls.dfpt(structure, relaxed_ion=relaxed_ion, piezo=True, dde=False, strain=True, dte=False,
                       stress_correction=stress_correction, asr=asr, chneut=chneut, dipdip=dipdip,
                       anaddb_args=anaddb_args, anaddb_kwargs=anaddb_kwargs, comment=comment)
        return new

    @classmethod
    def phbands_and_dos(cls, structure, ngqpt, nqsmall, qppa=None, ndivsm=20, line_density=None, q1shft=(0, 0, 0),
                        qptbounds=None, asr=2, chneut=0, dipdip=1, dos_method="tetra", lo_to_splitting=False,
                        anaddb_args=None, anaddb_kwargs=None, spell_check=False, comment=None):
        """
        Build an |AnaddbInput| for the computation of phonon bands and phonon DOS.

        Args:
            structure: |Structure| object
            ngqpt: Monkhorst-Pack divisions for the phonon Q-mesh (coarse one)
            nqsmall: Used to generate the (dense) mesh for the DOS.
                It defines the number of q-points used to sample the smallest lattice vector.
            qppa: Defines the homogeneous q-mesh used for the DOS in units of q-points per reciproval atom.
                Overrides nqsmall.
            line_density: Defines the a density of k-points per reciprocal atom to plot the phonon dispersion.
                Overrides ndivsm.
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
            spell_check: False to disable spell checking for input variables.
            comment: Optional string with a comment that will be placed at the beginning of the file.
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
            raise NotImplementedError("Wrong value for dos_method: %s" % str(dos_method))

        new = cls(structure, comment="ANADDB input for phonon bands and DOS" if not comment else comment,
                  anaddb_args=anaddb_args, anaddb_kwargs=anaddb_kwargs, spell_check=spell_check)

        # Parameters for the DOS
        if qppa:
            ng2qpt = KSampling.automatic_density(structure, kppa=qppa).kpts[0]
            # Set new variables
            new.set_vars(ng2qpt=ng2qpt,prtdos=prtdos,dossmear=dossmear)
        else:
            new.set_autoqmesh(nqsmall)
            new.set_vars(prtdos=prtdos, dosdeltae=dosdeltae, dossmear=dossmear)
            if nqsmall == 0:
                new["prtdos"] = 0

        # Parameters for the Bandstructure.
        if line_density:
            hs = HighSymmKpath(structure, symprec=1e-2)
            qpts, labels_list = hs.get_kpoints(line_density=line_density, coords_are_cartesian=False)
            # remove repeated q-points since those do
            # interfere with _set_split_vals in the phonon plotter
            qph1l = [qpts[0]]
            for qpt in qpts[1:]:
                if not np.array_equal(qpt,qph1l[-1]):
                    qph1l.append(qpt)
            #set new variables
            new['qph1l'] = [q.tolist()+[1] for q in qph1l]
            new['nph1l'] = len(qph1l)
            qptbounds = structure.calc_kptbounds()
        else:
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
                    if i > 0:
                        directions.extend(structure.lattice.reciprocal_lattice_crystallographic.get_cartesian_coords(qptbounds[i-1]))
                        directions.append(0)

                    if i < len(qptbounds) - 1:
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
    def modes(cls, structure, enunit=2, asr=2, chneut=1, anaddb_args=None, anaddb_kwargs=None):
        """
        Build an |AnaddbInput| for the computation of phonon modes.

        Args:
            Structure: |Structure| object
            ngqpt: Monkhorst-Pack divisions for the phonon Q-mesh (coarse one)
            nqsmall: Used to generate the (dense) mesh for the DOS.
                It defines the number of q-points used to sample the smallest lattice vector.
            q1shft: Shifts used for the coarse Q-mesh
            qptbounds Boundaries of the path. If None, the path is generated from an internal database
                depending on the input structure.
            asr, chneut, dipdp: Anaddb input variable. See official documentation.
            anaddb_args: List of tuples (key, value) with Anaddb input variables (default: empty)
            anaddb_kwargs: Dictionary with Anaddb input variables (default: empty)
        """
        new = cls(structure, comment="ANADDB input for modes", anaddb_args=anaddb_args, anaddb_kwargs=anaddb_kwargs)

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

    @classmethod
    def ifc(cls, structure, ngqpt, ifcout=None, q1shft=(0, 0, 0), asr=2, chneut=1, dipdip=1, anaddb_args=None,
            anaddb_kwargs=None):
        """
        Build an |AnaddbInput| for the computation of interatomic force constants.

        Args:
            structure: |Structure| object
            ngqpt: Monkhorst-Pack divisions for the phonon Q-mesh (coarse one)
            ifcout: Number of neighbouring atoms for which the ifc's will be output. If None all the atoms in the big box.
            q1shft: Shifts used for the coarse Q-mesh
            asr, chneut, dipdip: Anaddb input variable. See official documentation.
            anaddb_args: List of tuples (key, value) with Anaddb input variables (default: empty)
            anaddb_kwargs: Dictionary with Anaddb input variables (default: empty)
        """
        new = cls(structure, comment="ANADDB input for IFC",
                  anaddb_args=anaddb_args, anaddb_kwargs=anaddb_kwargs)

        q1shft = np.reshape(q1shft, (-1, 3))

        #TODO add in anaddb an option to get all the atoms if ifcout<0
        # Huge number abinit will limit to the big box
        ifcout = ifcout or 10000000

        new.set_vars(
            ifcflag=1,
            ngqpt=np.array(ngqpt),
            q1shft=q1shft,
            nqshft=len(q1shft),
            asr=asr,
            chneut=chneut,
            dipdip=dipdip,
            ifcout=ifcout,
            natifc=len(structure),
            atifc=list(range(1, len(structure)+1)),
            ifcana=1,
            prt_ifc=1
        )

        return new

    @classmethod
    def dfpt(cls, structure, ngqpt=None, relaxed_ion=False, piezo=False, dde=False, strain=False, dte=False,
             stress_correction=False, nqsmall=None, qppa=None, ndivsm=20, line_density=None, q1shft=(0, 0, 0),
             qptbounds=None, asr=2, chneut=1, dipdip=1, dos_method="tetra", anaddb_args=None, anaddb_kwargs=None, comment=None):
        """
        Builds an |AnaddbInput| to post-process a generic DFPT calculation.

        Args:
            structure: |Structure| object.
            ngqpt: Monkhorst-Pack divisions for the phonon Q-mesh (coarse one)
            stress_correction: True to activate computation of  stress correction in elastic tensor.
                Requires DDB with stress entries.
            relaxed_ion: True to activate computation of relaxed-ion elastic and piezoelectric tensors.
                (assume the DDB has atomic perturbations at Gamma)
            piezo: if True the piezoelectric tensor are calculated (requires piezoelectric perturbations)
            dde: if True dielectric tensors will be calculated. If phonon band
                structure is calculated will also enable the calculation of the lo_to splitting
                (requires the DDE perturbations)
            strain: if True the elastic tensors will be calculated (requires the strain perturbations)
            dte: if True properties related to the nonlinear tensors will be calculated
                (requires third orders perturbations)
            nqsmall: Used to generate the (dense) mesh for the DOS.
                It defines the number of q-points used to sample the smallest lattice vector.
            qppa: Defines the homogeneous q-mesh used for the DOS in units of q-points per reciproval atom.
                Overrides nqsmall.
            line_density: Defines the a density of k-points per reciprocal atom to plot the phonon dispersion.
                Overrides ndivsm.
            ndivsm: Used to generate a normalized path for the phonon bands.
                If gives the number of divisions for the smallest segment of the path.
            q1shft: Shifts used for the coarse Q-mesh
            qptbounds Boundaries of the path. If None, the path is generated from an internal database
                depending on the input structure.
            asr, chneut, dipdp: Anaddb input variable. See official documentation.
            dos_method: Possible choices: "tetra", "gaussian" or "gaussian:0.001 eV".
                In the later case, the value 0.001 eV is used as gaussian broadening
            anaddb_args: List of tuples (key, value) with Anaddb input variables (default: empty)
            anaddb_kwargs: Dictionary with Anaddb input variables (default: empty)
            comment: Optional string with a comment that will be placed at the beginning of the file.
        """
        # use the phonon BS and DOS input as starting point is required, otherwise
        if ngqpt:
            anaddb_input = cls.phbands_and_dos(structure=structure, ngqpt=ngqpt, ndivsm=ndivsm, nqsmall=nqsmall,
                                               qppa=qppa, line_density=line_density, asr=asr, chneut=chneut,
                                               dipdip=dipdip, qptbounds=qptbounds, dos_method=dos_method,
                                               lo_to_splitting=dde, q1shft=q1shft, comment=comment)
        else:
            anaddb_input = AnaddbInput(structure, comment=comment)
            anaddb_input.set_vars(asr=asr, chneut=chneut)

        dieflag = 0
        if dde:
            dieflag = 3 if (relaxed_ion and strain) else 2

        elaflag = 0
        if strain:
            if not relaxed_ion:
                elaflag = 1
            elif stress_correction:
                elaflag = 5
            else:
                elaflag = 3

        piezoflag = 0
        if piezo:
            if not relaxed_ion:
                piezoflag = 1
            elif dde and strain:
                piezoflag = 7
            else:
                piezoflag = 3

        anaddb_input.set_vars(dieflag=dieflag, elaflag=elaflag, piezoflag=piezoflag)

        if dieflag == 3 and 'nph2l' not in anaddb_input:
            anaddb_input['nph2l'] = 1

        if elaflag > 1:
            anaddb_input["instrflag"] = 1

        if dte:
            anaddb_input.set_vars(nlflag=1,
                                  ramansr=1,
                                  alphon=1,
                                  prtmbm=1)

        anaddb_args = [] if anaddb_args is None else anaddb_args
        anaddb_kwargs = {} if anaddb_kwargs is None else anaddb_kwargs
        args = list(anaddb_args)[:]
        args.extend(list(anaddb_kwargs.items()))
        anaddb_input.set_vars(args)

        return anaddb_input

    @property
    def structure(self):
        """|Structure| object."""
        return self._structure

    def to_string(self, sortmode=None, mode="text", verbose=0):
        """
        String representation.

        Args:
            sortmode: "a" for alphabetical order, None if no sorting is wanted
            mode: Either `text` or `html` if HTML output with links is wanted.
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

        root = "https://docs.abinit.org/variables/anaddb/"
        for varname in keys:
            value = self[varname]
            if mode == "html": varname = root + "#%s" % varname
            app(str(InputVariable(varname, value)))

        return "\n".join(lines) if mode == "text" else "\n".join(lines).replace("\n", "<br>")

    def _repr_html_(self):
        """Integration with jupyter_ notebooks."""
        return self.to_string(mode="html")

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

    def abivalidate(self, workdir=None, manager=None):
        """
        Run ANADDB in dry-run mode to validate the input file.

        Args:
            workdir: Working directory of the fake task used to compute the ibz. Use None for temporary dir.
            manager: |TaskManager| of the task. If None, the manager is initialized from the config file.

        Return:
            `namedtuple` with the following attributes:

                retcode: Return code. 0 if OK.
                output_file: output file of the run.
                log_file:  log file of the Abinit run, use log_file.read() to access its content.
                stderr_file: stderr file of the Abinit run. use stderr_file.read() to access its content.
                task: Task object
        """
        task = AnaddbTask.temp_shell_task(self, ddb_node="fake_DDB", workdir=workdir, manager=manager)
        # TODO: Anaddb does not support --dry-run
        #retcode = task.start_and_wait(autoparal=False, exec_args=["--dry-run"])
        return dict2namedtuple(retcode=0, output_file=task.output_file, log_file=task.log_file,
                               stderr_file=task.stderr_file, task=task)


class OpticVar(collections.namedtuple("OpticVar", "name default group help")):

    def __str__(self):
        sval = str(self.default)
        return (4*" ").join([sval, "!" + self.help])

    @property
    def url(self):
        """The url associated to the variable."""
        root = "https://docs.abinit.org/variables/optic/"
        return root + "#%s" % self.name

    def html_link(self, label=None):
        """String with the URL of the web page."""
        return '<a href="%s" target="_blank">%s</a>' % (self.url, self.name if label is None else label)



class OpticError(Exception):
    """Error class raised by OpticInput."""


class OpticInput(AbiAbstractInput, MSONable):
    """
    Input file for optic executable
    """
    Error = OpticError

    # variable name --> default value.
    _VARIABLES = [
        #OpticVar(name="ddkfile_x", default=None, help="Name of the first d/dk response wavefunction file"),
        #OpticVar(name="ddkfile_y", default=None, help="Name of the second d/dk response wavefunction file"),
        #OpticVar(name="ddkfile_z", default=None, help="Name of the third d/dk response wavefunction file"),
        #OpticVar(name="wfkfile",   default=None, help="Name of the ground-state wavefunction file"),

        # PARAMETERS section:
        OpticVar(name="broadening", default=0.01, group='PARAMETERS',
                 help="Value of the smearing factor, in Hartree"),
        OpticVar(name="domega", default=0.010, group='PARAMETERS',
                 help="Frequency step (Ha)"),
        OpticVar(name="maxomega", default=1, group='PARAMETERS',
                 help="Maximum frequency (Ha)"),
        OpticVar(name="scissor", default=0.000, group='PARAMETERS',
                 help="Scissor shift if needed, in Hartree"),
        OpticVar(name="tolerance", default=0.001, group='PARAMETERS',
                 help="Tolerance on closeness of singularities (in Hartree)"),
        OpticVar(name="autoparal", default=0, group='PARAMETERS',
                 help="Autoparal option"),
        OpticVar(name="max_ncpus", default=0, group='PARAMETERS',
                 help="Max number of CPUs considered in autoparal mode"),

        # COMPUTATIONS section:
        OpticVar(name="num_lin_comp", default=0, group='COMPUTATIONS',
                 help="Number of components of linear optic tensor to be computed"),
        OpticVar(name="lin_comp", default=0, group='COMPUTATIONS',
                 help="Linear coefficients to be computed (x=1, y=2, z=3)"),
        OpticVar(name="num_nonlin_comp", default=0, group='COMPUTATIONS',
                 help="Number of components of nonlinear optic tensor to be computed"),
        OpticVar(name="nonlin_comp", default=0, group='COMPUTATIONS',
                 help="Non-linear coefficients to be computed"),
        OpticVar(name="num_linel_comp", default=0, group='COMPUTATIONS',
                 help="Number of components of linear electro-optic tensor to be computed"),
        OpticVar(name="linel_comp", default=0, group='COMPUTATIONS',
                 help="Linear electro-optic coefficients to be computed"),
        OpticVar(name="num_nonlin2_comp", default=0, group='COMPUTATIONS',
                 help="Number of components of nonlinear optic tensor v2 to be computed"),
        OpticVar(name="nonlin2_comp", default=0, group='COMPUTATIONS',
                 help="Non-linear coefficients v2 to be computed"),
    ]

    _GROUPS = ['PARAMETERS','COMPUTATIONS']

    # Variable names supported
    _VARNAMES = [v.name for v in _VARIABLES]

    # Mapping name --> var object.
    _NAME2VAR = {v.name: v for v in _VARIABLES}

    def __init__(self, **kwargs):
        # Initalize with default values.
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
        """Return the default value of variable `key`."""
        for var in self._VARIABLES:
            if var.name == key: return var.default
        raise self.Error("Cannot find %s in _VARIABLES" % str(key))

    @classmethod
    def from_dict(cls, d):
        """
        JSON interface used in pymatgen for easier serialization.
        """
        kwargs = {}
        for grp, section in d.items():
            if grp in ("@module", "@class"): continue
            kwargs.update(**section)
        return cls(**kwargs)

    # TODO
    #@pmg_serialize
    def as_dict(self):
        """
        JSON interface used in pymatgen for easier serialization.
        """
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
            my_dict[grp].update({name: value})

        return my_dict

    def to_string(self, verbose=0):
        """String representation."""
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

    #def _repr_html_(self):
    #    """Integration with jupyter notebooks."""
    #    return self.to_string(mode="html")

    def only_independent_chi_components(self, structure, assume_symmetric_tensor=False,
                                        symprec=1e-3, angle_tolerance=5):
        """
        Use the crystal system returned by spglib to find the independent components
        of the linear susceptibility tensor and set the appropriate variables.

        Args:
            structure: Crystalline structure
            assume_symmetric_tensor: True if tensor can be assumed symmetric.
                Note that the tensor is symmetric only for a lossless and non-optically active material.
            symprec, angle_tolerance: Parameters passed to spglib.

        Return:
            Set internal variables and return list of components to compute.
        """
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        spgan = SpacegroupAnalyzer(structure, symprec=symprec, angle_tolerance=angle_tolerance)
        system = spgan.get_crystal_system()

        # Table 1.5.1 of https://booksite.elsevier.com/samplechapters/9780123694706/Sample_Chapters/02~Chapter_1.pdf.
        # Note that the tensor is symmetric only for a lossless and non-optically active material.
        components_for_system = {
            "triclinic": "xx yy zz xy yx xz zx yz zy",
            "monoclinic": "xx yy zz xz zx",
            "orthorhombic": "xx yy zz",
            "tetragonal": "xx zz",
            "cubic": "xx",
        }

        if assume_symmetric_tensor:
            components_for_system["triclinic"] = "xx yy zz xy xz yz"
            components_for_system["monoclinic"] = "xx yy zz xz"

        components_for_system["trigonal"] = components_for_system["tetragonal"]
        components_for_system["hexagonal"] = components_for_system["tetragonal"]

        for k, v in components_for_system.items():
            components_for_system[k] = v.split()

        ind_comps = components_for_system[system]
        d = {"x": 1, "y": 2, "z": 3}
        self["num_lin_comp"] = len(ind_comps)
        self["lin_comp"] = [10 * d[comp[0]] + d[comp[1]] for comp in ind_comps]

        return ind_comps

    def abivalidate(self, workdir=None, manager=None):
        """
        Run OPTIC in dry-run mode to validate the input file.
        Note: This method is a stub, it always return retcode 0

        Args:
            workdir: Working directory of the fake task used to compute the ibz. Use None for temporary dir.
            manager: |TaskManager| of the task. If None, the manager is initialized from the config file.

        Return:
            `namedtuple` with the following attributes:

                retcode: Return code. 0 if OK.
                output_file: output file of the run.
                log_file:  log file of the Abinit run, use log_file.read() to access its content.
                stderr_file: stderr file of the Abinit run. use stderr_file.read() to access its content.
                task: Task object
        """
        # TODO: Optic does not support --dry-run
        #task = OpticTask.temp_shell_task(inp=self, workdir=workdir, manager=manager)
        #retcode = task.start_and_wait(autoparal=False, exec_args=["--dry-run"])
        return dict2namedtuple(retcode=0, output_file=None, log_file=None,
                               stderr_file=None, task=None)


class Cut3DInput(MSONable, object):
    """
    This object stores the options to run a single cut3d analysis.

    .. warning::

        Converters with nspden > 1 won't work since cut3d asks for the ispden index.
    """
    def __init__(self, infile_path=None, output_filepath=None, options=None):
        """
        Args:
            infile_path: absolute or relative path to the input file produced by abinit (e.g. DEN, WFK, ...). Can be
                None to be defined at a later time.
            output_filepath: path to the file that should be produced by cut3D, if required. At this stage it would be
                safer to use just the file name, as using an absolute or relative path may fail depending on
                the compiler.
            options: a list of strings that defines the options to be passed to cut3d
        """
        self.infile_path = infile_path
        self.output_filepath = output_filepath
        self.options = [str(o) for o in options]

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """Returns a string with the input."""
        lines = [self.infile_path]
        lines.extend(self.options)
        return '\n'.join(lines)

    def write(self, filepath):
        """Writes the input to a file."""
        if self.infile_path is None or self.options is None:
            raise ValueError("Infile path and options should be provided")

        with open(filepath, 'wt') as f:
            f.write(self.to_string())

    @classmethod
    def _convert(cls, infile_path, output_filepath, out_option):
        """
        Generic function used to generate the input for convertions using cut3d

        Args:
            infile_path: absolute or relative path to the input file produced by abinit (e.g. DEN, WFK, ...). Can be
                None to be defined at a later time.
            output_filepath: path to the file that should be produced by cut3D, if required. At this stage it would be
                safer to use just the file name, as using an absolute or relative path may fail depending on
                the compiler.
            out_option: a number corresponding to the required converting option in cut3d
        """
        options = [str(out_option)]  # Option to convert a _DEN file
        options.append(output_filepath)  # Name of the output file
        options.append('0')  # No more analysis
        return cls(infile_path=infile_path, output_filepath=output_filepath, options=options)

    @classmethod
    def den_to_3d_formatted(cls, density_filepath, output_filepath):
        """
        Generates a cut3d input for the conversion to a 3D formatted format.

        Args:
            density_filepath: absolute or relative path to the input density produced by abinit.
                Can be None to be defined at a later time.
            output_filepath: path to the file that should be produced by cut3D, if required. At this stage it would be
                safer to use just the file name, as using an absolute or relative path may fail depending on
                the compiler.
        """
        return cls._convert(density_filepath, output_filepath, 5)

    @classmethod
    def den_to_3d_indexed(cls, density_filepath, output_filepath):
        """
        Generates a cut3d input for the conversion to a 3D indexed format.

        Args:
            density_filepath: absolute or relative path to the input density produced by abinit. Can be None to be
                defined at a later time.
            output_filepath: path to the file that should be produced by cut3D, if required. At this stage it would be
                safer to use just the file name, as using an absolute or relative path may fail depending on
                the compiler.
        """
        return cls._convert(density_filepath, output_filepath, 6)

    @classmethod
    def den_to_molekel(cls, density_filepath, output_filepath):
        """
        Generates a cut3d input for the conversion to a Molekel format.

        Args:
            density_filepath: absolute or relative path to the input density produced by abinit. Can be None to be
                defined at a later time.
            output_filepath: path to the file that should be produced by cut3D, if required. At this stage it would be
                safer to use just the file name, as using an absolute or relative path may fail depending on
                the compiler.
        """
        return cls._convert(density_filepath, output_filepath, 7)

    @classmethod
    def den_to_tecplot(cls, density_filepath, output_filepath):
        """
        Generates a cut3d input for the conversion to a Tecplot format.

        Args:
            density_filepath: absolute or relative path to the input density produced by abinit. Can be None to be
                defined at a later time.
            output_filepath: path to the file that should be produced by cut3D, if required. At this stage it would be
                safer to use just the file name, as using an absolute or relative path may fail depending on
                the compiler.
        """
        return cls._convert(density_filepath, output_filepath, 8)

    @classmethod
    def den_to_xsf(cls, density_filepath, output_filepath, shift=None):
        """
        Generates a cut3d input for the conversion to an xsf format.

        Args:
            density_filepath: absolute or relative path to the input density produced by abinit. Can be None to be
                defined at a later time.
            output_filepath: path to the file that should be produced by cut3D, if required. At this stage it would be
                safer to use just the file name, as using an absolute or relative path may fail depending on
                the compiler.
            shift: a list of three integers defining the shift along the x, y, z axis.
                None if no shift is required.
        """
        options = ['9']  # Option to convert a _DEN file to an .xsf file
        options.append(output_filepath)  # Name of the output .xsf file
        if shift is not None:
            options.append('y')
            options.append("{} {} {} ".format(*shift))
        else:
            options.append('n')
        options.append('0')  # No more analysis
        return cls(infile_path=density_filepath, output_filepath=output_filepath, options=options)

    @classmethod
    def den_to_cube(cls, density_filepath, output_filepath):
        """
        Generates a cut3d input for the conversion to a cube format.

        Args:
            density_filepath: absolute or relative path to the input density produced by abinit. Can be None to be
                defined at a later time.
            output_filepath: path to the file that should be produced by cut3D, if required. At this stage it would be
                safer to use just the file name, as using an absolute or relative path may fail depending on
                the compiler.
        """
        return cls._convert(density_filepath, output_filepath, 14)

    @classmethod
    def hirshfeld(cls, density_filepath, all_el_dens_paths):
        """
        Generates a cut3d input for the calculation of the Hirshfeld charges from the density.

        Args:
            density_filepath: absolute or relative path to the input density produced by abinit. Can be None to be
                defined at a later time.
            all_el_dens_paths: a list of paths to the all-electron density files corresponding to the elements defined
                in the abinit input. See https://www.abinit.org/downloads/all_core_electron for files.
        """
        options = ['11']  # Option to convert _DEN file to a .cube file
        for p in all_el_dens_paths:
            options.append(p)
        options.append('0')

        return cls(infile_path=density_filepath, options=options)

    @classmethod
    def hirshfeld_from_fhi_path(cls, density_filepath, structure, fhi_all_el_path):
        """
        Generates a cut3d input for the calculation of the Hirshfeld charges from the density. Automatically
        selects the all-electron density files from a folder containing the fhi all-electron density files:
        https://www.abinit.org/downloads/all_core_electron

        This will work only if the input has been generated with AbinitInput and the Structure object is the same
        provided to AbinitInput.

        Args:
            density_filepath: absolute or relative path to the input density produced by abinit. Can be None to be
                defined at a later time.
            structure: the structure used for the ground state calculation. Used to determine the elements
            fhi_all_el_path: path to the folder containing the fhi all-electron density files
        """
        all_el_dens_paths = []
        # This relies on AbinitInput using Structure.types_of_specie to define znucl
        for e in structure.types_of_specie:
            all_el_dens_paths.append(os.path.join(fhi_all_el_path, "0.{:02}-{}.8.density.AE".format(e.number, e.name)))

        return cls.hirshfeld(density_filepath, all_el_dens_paths)

    @pmg_serialize
    def as_dict(self):
        """
        JSON interface used in pymatgen for easier serialization.
        """
        return dict(infile_path=self.infile_path, output_filepath=self.output_filepath, options=self.options)

    @classmethod
    def from_dict(cls, d):
        """
        JSON interface used in pymatgen for easier serialization.
        """
        return cls(infile_path=d.get('infile_path', None), output_filepath=d.get('output_filepath', None),
                   options=d.get('options', None))


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
