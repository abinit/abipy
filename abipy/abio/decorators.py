# coding: utf-8
"""Decorators for AbinitInput or MultiDataset objects."""
from __future__ import print_function, division, unicode_literals, absolute_import

import six
import abc
import pymatgen.io.abinit.abiobjects as aobj

from monty.inspect import initializer
from monty.json import MSONable
try:
    from pymatgen.util.serialization import pmg_serialize
except ImportError:
    from pymatgen.serializers.json_coders import pmg_serialize

from abipy.flowtk.abiobjects import LdauParams, LexxParams
from .inputs import AbinitInput, MultiDataset

import logging
logger = logging.getLogger(__file__)


class InputDecoratorError(Exception):
    """Error class raised by :class:`AbinitInputDecorator`."""


class AbinitInputDecorator(six.with_metaclass(abc.ABCMeta, MSONable)):
    """
    An `AbinitInputDecorator` adds new options to an existing :class:`AbinitInput`
    or an existing :class:`MultiDataset` without altering its structure. This is an abstract Base class.

    Example:

        decorator = MyDecorator(arguments)

        new_abinit_input = decorator(abinit_input)
        new_multidataset = decorator(multidataset)

    Note that a decorator does not modify the object on which it acts.

    .. warning::

        Please avoid introducing decorators acting on the structure (in particular the lattice)
        since the initial input may use the initial structure to  compute important variables.
        For instance, the list of k-points for band structure calculation depend on the bravais lattice
        and a decorator that changes it should recompute the path.
        This should  not represent a serious limitation because it's always possible to change the structure
        with its methods and then call the factory function without having to decorate an already existing object.
    """
    Error = InputDecoratorError

    def __str__(self):
        return str(self.as_dict())

    def __call__(self, obj, deepcopy=True):
        """
        Decorate an `AbinitInput` or a `MultiDataset` object
        This is the public API that calls the concrete implementation of the subclass

        Returns:
            New `AbinitInput` or new `MultiDataset` depending on obj.
        """

        if isinstance(obj, AbinitInput):
            new_inp = self._decorate(obj, deepcopy=deepcopy)
            # Log the decoration in new_inp.
            new_inp.register_decorator(self)
            return new_inp

        elif isinstance(obj, MultiDataset):
            new_inputs = []
            for inp in obj:
                new_inp = self._decorate(inp, deepcopy=deepcopy)
                # Log the decoration in new_inp.
                new_inp.register_decorator(self)
                new_inputs.append(new_inp)

            return MultiDataset.from_inputs(new_inputs)

        else:
            raise TypeError("Don't know how to decorate type %s" % type(obj))

    @abc.abstractmethod
    def _decorate(self, inp, deepcopy=True):
        """
        Abstract method that must be implemented by the concrete classes.
        It receives a :class:`AbinitInput` object, applies the decoration and returns a new `AbinitInput`.

        Args:
            inp: :class:`AbinitInput` object.
            deepcopy: True if a deepcopy of inp should be performed before changing the object.

        Returns:
            decorated :class:`AbinitInput` object (new object)
        """


class SpinDecorator(AbinitInputDecorator):
    """This decorator changes the spin polarization."""
    def __init__(self, spinmode, kptopt_ifspinor=4):
        self.spinmode = aobj.SpinMode.as_spinmode(spinmode)
        self.kptopt_ifspinor = kptopt_ifspinor

    @pmg_serialize
    def as_dict(self):
        return dict(spinmode=self.spinmode.as_dict(), kptopt_ifspinor=self.kptopt_ifspinor)

    @classmethod
    def from_dict(cls, d):
        return cls(aobj.SpinMode.from_dict(d["spinmode"]), kptopt_ifspinor=d["kptopt_ifspinor"])

    def _decorate(self, inp, deepcopy=True):
        if deepcopy: inp = inp.deepcopy()

        inp.set_vars(self.spinmode.to_abivars())

        # in version 7.11.5
        # When non-collinear magnetism is activated (nspden=4),
        # time-reversal symmetry cannot be used in the present
        # state of the code (to be checked and validated).
        # Action: choose kptopt different from 1 or 2.

        # Here we set kptopt to 4 (spatial symmetries, no time-reversal)
        # unless we already have a dataset with kptopt == 3 (no tr, no spatial)
        # This case is needed for DFPT at q != 0.
        if self.spinmode.nspinor == 2 and inp.get("kptopt") != 3:
            inp.set_vars(kptopt=self.kptopt_ifspinor)

        return inp


class SmearingDecorator(AbinitInputDecorator):
    """This decorator changes the electronic smearing."""
    def __init__(self, smearing):
        self.smearing = aobj.Smearing.as_smearing(smearing)

    @pmg_serialize
    def as_dict(self):
        return {"smearing": self.smearing.as_dict()}

    @classmethod
    def from_dict(cls, d):
        return cls(aobj.Smearing.from_dict(d["smearing"]))

    def _decorate(self, inp, deepcopy=True):
        if deepcopy: inp = inp.deepcopy()
        inp.set_vars(self.smearing.to_abivars())
        return inp


class XcDecorator(AbinitInputDecorator):
    """Change the exchange-correlation functional."""
    def __init__(self, ixc):
        """
        Args:
            ixc: Abinit input variable
        """
        self.ixc = ixc

    @pmg_serialize
    def as_dict(self):
        return {"ixc": self.ixc}

    @classmethod
    def from_dict(cls, d):
        return cls(d["ixc"])

    def _decorate(self, inp, deepcopy=True):
        if deepcopy: inp = inp.deepcopy()
        # TODO: Don't understand why abinit does not enable usekden if MGGA!
        usekden = None
        #usekden = 1 if ixc.ismgga() else None
        inp.set_vars(ixc=self.ixc, usekden=usekden)

        return inp


class LdaUDecorator(AbinitInputDecorator):
    """This decorator adds LDA+U parameters to an :class:`AbinitInput` object."""
    def __init__(self, symbols_luj, usepawu=1, unit="eV"):
        """
        Args:
            symbols_luj: dictionary mapping chemical symbols to another dict with (l, u, j) values
            usepawu: Abinit input variable.
            unit: Energy unit for U and J
        """
        self.symbols_luj, self.usepawu, self.unit = symbols_luj, usepawu, unit

    @pmg_serialize
    def as_dict(self):
        return dict(symbols_luj=self.symbols_luj, usepawu=self.usepawu, unit=self.unit)

    @classmethod
    def from_dict(cls, d):
        return cls(**{k: v for k, v in d.items() if not k.startswith("@")})

    def _decorate(self, inp, deepcopy=True):
        if not inp.ispaw: raise self.Error("LDA+U requires PAW!")
        if deepcopy: inp = inp.deepcopy()
        luj_params = LdauParams(usepawu=self.usepawu, structure=inp.structure)

        # Apply UJ on all the symbols present in symbols_lui.
        for symbol in inp.structure.symbol_set:
            if symbol not in self.symbols_luj: continue
            args = self.symbols_luj[symbol]
            luj_params.luj_for_symbol(symbol, l=args["l"], u=args["u"], j=args["j"], unit=self.unit)
            #luj_params.luj_for_symbol("Ni", l=2, u=u, j=0.1*u, unit=self.unit)

        inp.set_vars(luj_params.to_abivars())

        return inp


class LexxDecorator(AbinitInputDecorator):
    """This decorator add local exact exchange to an :class:`AbinitInput` object."""
    def __init__(self, symbols_lexx, exchmix=None):
        """
        Args:
            symbols_lexx: dictionary mapping chemical symbols to the angular momentum l on which lexx is applied.
            exchmix: ratio of exact exchange when useexexch is used. The default value of 0.25 corresponds to PBE0.

            Example. To perform a LEXX calculation for NiO in which the LEXX is computed only for the l=2
            channel of the nickel atoms:

                {"Ni": 2}
        """
        self.symbols_lexx, self.exchmix = symbols_lexx, exchmix

    @classmethod
    def from_dict(cls, d):
        return cls(**{k:v for k, v in d.items() if not k.startswith("@")})

    @pmg_serialize
    def as_dict(self):
        return {"symbols_lexx": self.symbols_lexx, "exchmix": self.exchmix}

    def _decorate(self, inp, deepcopy=True):
        if not inp.ispaw: raise self.Error("LEXX requires PAW!")
        if deepcopy: inp = inp.deepcopy()

        lexx_params = LexxParams(inp.structure)
        for symbol in inp.structure.symbol_set:
            if symbol not in self.symbols_lexx: continue
            lexx_params.lexx_for_symbol(symbol, l=self.symbols_lexx[symbol])

        # Context : the value of the variable useexexch is   1.
        # The value of the input variable ixc is    7, while it must be
        # equal to one of the following:  11  23
        # Action : you should change the input variables ixc or useexexch.
        inp.set_vars(lexx_params.to_abivars())
        dt_ixc = inp.get("ixc")
        if dt_ixc is None or ixc not in [11, 23]: inp.set_vars(ixc=11)
        if self.exchmix is not None: inp.set_vars(exchmix=self.exchmix)

        return inp


# Stubs
#class ScfMixingDecorator(AbinitInputDecorator):


#class MagneticMomentDecorator(AbinitInputDecorator):
#    """Add reasoanble guesses for the initial magnetic moments."""


#class SpinOrbitDecorator(AbinitInputDecorator):
#    """Enable spin-orbit in the input."""
#     def __init__(self, no_spatial_symmetries=True, no_time_reversal=False, spnorbscl=None):
#        self.use_spatial_symmetries = use_spati
#        self.use_spatial_symmetries
#
#    def _decorate(self, inp, deepcopy=True)
#        if deepcopy: inp = inp.deepcopy()
#        kptopt =
#        if inp.ispaw:
#            for dt in inp.datasets:
#               dt.set_vars(pawspnorb=1, kptopt=kptopt)
#        return inp


#class PerformanceDecorator(AbinitInputDecorator):
#    """Change the variables in order to speedup the calculation."""
#    fftgw
#    boxcutmin
#    fft
#    def __init__(self, accuracy):
#        self.accuracy = accuracy
#
#    def _decorate(self, inp, deepcopy=True)
#        if deepcopy: inp = inp.deepcopy()
#         for dt in inp[1:]:
#            runlevel = dt.runlevel
#        return inp


#class DmftDecorator(AbinitInputDecorator):
#    """Add DMFT variables."""
