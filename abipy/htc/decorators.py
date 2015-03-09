# coding: utf-8
"""Decorators for AbiInput objects."""
from __future__ import print_function, division, unicode_literals

import six
import abc

from pymatgen.serializers.json_coders import PMGSONable


import logging
logger = logging.getLogger(__file__)


class DecoratorError(Exception):
    """Error class raised by :class:`InputDecorator`."""


class InputDecorator(six.with_metaclass(abc.ABCMeta, PMGSONable)):
    """
    An `InputDecorator` adds new options to an existing :class:`AbiInput` without altering its structure. 
    This is an abstract Base class.

    .. warning::

        Please avoid introducing decorators acting on the structure (in particular the lattice) 
        since the initial input may use the initial structure to  compute important variables. 
        For instance, the list of k-points for band structure calculation depend on the bravais lattice 
        and a decorator that changes it should recompute the path. 
        This should  not represent a serious limitation because it's always possible to change the structure 
        with its methods and then call the factory function without having to decorate an already existing object.
    """
    Error = DecoratorError

    #@abc.abstractmethod
    #def as_dict(self):
    #    """Return a dict with the PMGSON representation."""

    def decorate(self, inp, deepcopy=True):
        new_inp = self._decorate(inp, deepcopy=deepcopy)
        # Log the decoration in new_inp.
        new_inp._decorators.append(self.as_dict())
        return new_inp

    @abc.abstractmethod
    def _decorate(self, inp, deepcopy=True):
        """
        Abstract method that must be implemented by the concrete classes.
        It receives a :class:`AbiInput` object, applies the decoration and returns a new `AbiInput`.

        Args:
            inp: :class:`AbiInput` object.
            deepcopy: True if a deepcopy of inp should be performed before changing the object.

        Returns:
            decorated :class:`AbiInput` object (new object)
        """

# Stubs
#class SpinDecorator(InputDecorator):
#    def __init__(self, spinmode):
#        """Change the spin polarization."""
#        self.spinmode = aobj.SpinMpde.as_spinmode(spin_mode)
#
#    def _decorate(self, inp, deepcopy=True)
#        if deepcopy: inp = inp.deepcopy()
#        inp.set_vars(self.spinmode.to_abivars())
#        return inp
#
#
#class SmearingDecorator(InputDecorator):
#    """Change the electronic smearing."""
#    def __init__(self, spinmode):
#        self.smearing = aobj.Smearing.as_smearing(smearing)
#
#    def _decorate(self, inp, deepcopy=True)
#        if deepcopy: inp = inp.deepcopy()
#        inp.set_vars(self.smearing.to_abivars())
#        return inp


#class ScfMixingDecorator(InputDecorator):

#class MagneticMomentDecorator(InputDecorator):

#class AccuracyDecorator(InputDecorator):
#    """Change the electronic smearing."""
#    def __init__(self, accuracy):
#        self.accuracy = accuracy
#
#    def _decorate(self, inp, deepcopy=True)
#        if deepcopy: inp = inp.deepcopy()
#         for dt in inp[1:]:
#            runlevel = dt.runlevel 
#        return inp
#
#
#class UJDecorator(InputDecorator):
#    """Add LDA+U to an :class:`AbiInput` object."""
#    def __init__(self, luj_for_symbol, usepawu=1):
#        self.usepawu = usepawu
#        self.luj_for_symbol = luj_for_symbol
#
#    def _decorate(self, inp, deepcopy=True)
#        if not inp.ispaw: raise self.Error("LDA+U requires PAW!")
#        if deepcopy: inp = inp.deepcopy()
#
#        inp.set_vars(usepawu=self.usepawu)
#        return inp
#
#
#class LexxDecorator(InputDecorator):
#    """Add LDA+U to an :class:`AbiInput` object."""
#    def __init__(self, luj_for_symbol, usepawu=1):
#        self.usepawu = usepawu
#        self.luj_for_symbol = luj_for_symbol
#
#    def _decorate(self, inp, deepcopy=True)
#        if not inp.ispaw: raise self.Error("LDA+U requires PAW!")
#        if deepcopy: inp = inp.deepcopy()
#
#        inp.set_vars(usepawu=self.usepawu)
#        return inp
