from __future__ import print_function, division, unicode_literals, absolute_import

import collections

from pymatgen.core.units import Energy
from pymatgen.io.abinit.abiobjects import *


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
            structure: |Structure| object.
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
            structure: |Structure| object.
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
