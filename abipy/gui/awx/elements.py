# -*- coding: utf-8 -*-
# elements.py

# Copyright (c) 2005-2014, Christoph Gohlke
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# * Neither the name of the copyright holders nor the names of any
#   contributors may be used to endorse or promote products derived
#   from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

"""Properties of the chemical elements.

Each chemical element is represented as an object instance. Physicochemical
and descriptive properties of the elements are stored as instance attributes.

:Author: `Christoph Gohlke <http://www.lfd.uci.edu/~gohlke/>`_

:Version: 2013.03.18

Requirements
------------
* `CPython 2.7 or 3.3 <http://www.python.org>`_

References
----------
(1) http://physics.nist.gov/PhysRefData/Compositions/
(2) http://physics.nist.gov/PhysRefData/IonEnergy/tblNew.html
(3) http://en.wikipedia.org/wiki/%(element.name)s
(4) http://www.miranda.org/~jkominek/elements/elements.db

Examples
--------
>>> from elements import ELEMENTS
>>> len(ELEMENTS)
109
>>> str(ELEMENTS[109])
'Meitnerium'
>>> ele = ELEMENTS['C']
>>> ele.number, ele.symbol, ele.name, ele.eleconfig
(6, 'C', 'Carbon', '[He] 2s2 2p2')
>>> ele.eleconfig_dict
{(1, 's'): 2, (2, 'p'): 2, (2, 's'): 2}
>>> sum(ele.mass for ele in ELEMENTS)
14659.1115599
>>> for ele in ELEMENTS:
...     ele.validate()
...     ele = eval(repr(ele))

"""

from __future__ import division, print_function

__version__ = '2013.03.18'
__docformat__ = 'restructuredtext en'
__all__ = ['ELEMENTS']


class lazyattr(object):
    """Lazy object attribute whose value is computed on first access."""
    __slots__ = ['func']

    def __init__(self, func):
        self.func = func

    def __get__(self, instance, owner):
        result = self.func(instance)
        if result is NotImplemented:
            return getattr(super(owner, instance), self.func.__name__)
        setattr(instance, self.func.__name__, result)
        return result


class Element(object):
    """Chemical element.

    Attributes
    ----------
    number : int
        Atomic number
    symbol : str of length 1 or 2
        Chemical symbol
    name : str
        Name in english
    group : int
        Group in periodic table
    period : int
        Period in periodic table
    block : int
        Block in periodic table
    series : int
        Index to chemical series
    protons : int
        Number of protons
    neutrons : int
        Number of neutrons in the most abundant naturally occurring stable
        isotope
    nominalmass : int
        Mass number of the most abundant naturally occurring stable isotope
    electrons : int
        Number of electrons
    mass : float
        Relative atomic mass. Ratio of the average mass of atoms
        of the element to 1/12 of the mass of an atom of 12C
    exactmass : float
        Relative atomic mass calculated from the isotopic composition
    eleneg : float
        Electronegativity (Pauling scale)
    covrad : float
        Covalent radius in Angstrom
    atmrad :
        Atomic radius in Angstrom
    vdwrad : float
        Van der Waals radius in Angstrom
    tboil : float
        Boiling temperature in K
    tmelt : float
        Melting temperature in K
    density : float
        Density at 295K in g/cm3 respectively g/L
    oxistates : str
        Oxidation states
    eleaffin : float
        Electron affinity in eV
    eleconfig : str
        Ground state electron configuration
    eleconfig_dict : dict
        Ground state electron configuration (shell, subshell): electrons
    eleshells : int
        Number of electrons per shell
    ionenergy : tuple
        Ionization energies in eV
    isotopes : dict
        Isotopic composition.
        keys: isotope mass number
        values: Isotope(relative atomic mass, abundance)

    """
    def __init__(self, number, symbol, name, **kwargs):
        self.number = number
        self.symbol = symbol
        self.name = name
        self.electrons = number
        self.protons = number
        self.__dict__.update(kwargs)

    def __str__(self):
        return self.name

    def __repr__(self):
        ionenergy = []
        for i, j in enumerate(self.ionenergy):
            if i and (i % 5 == 0):
                ionenergy.append("\n" + " " * 15)
            ionenergy.append("%s, " % j)
        ionenergy = "".join(ionenergy)

        isotopes = []
        for massnum in sorted(self.isotopes):
            iso = self.isotopes[massnum]
            isotopes.append("%i: Isotope(%s, %s, %i)" % (
                massnum, iso.mass, iso.abundance, massnum))
        isotopes = ",\n              ".join(isotopes)

        description = word_wrap(self.description, linelen=66, indent=0,
                                joinstr=""" "\n        \"""")
        description = """    e['%s'].description = (\n        "%s\")""" % (
            self.symbol, description)
        # return description

        result = [
            "Element(\n    %i, '%s', '%s'" % (
                self.number, self.symbol, self.name),
            "group=%s, period=%s, block='%s', series=%i" % (
                self.group, self.period, self.block, self.series),
            "mass=%s, eleneg=%s, eleaffin=%s" % (
                self.mass, self.eleneg, self.eleaffin),
            "covrad=%s, atmrad=%s, vdwrad=%s" % (
                self.covrad, self.atmrad, self.vdwrad),
            "tboil=%s, tmelt=%s, density=%s" % (
                self.tboil, self.tmelt, self.density),
            "eleconfig='%s'" % self.eleconfig,
            "oxistates='%s'" % self.oxistates,
            "ionenergy=(%s)" % ionenergy,
            "isotopes={%s})" % isotopes
        ]
        return ",\n    ".join(result)

    @lazyattr
    def nominalmass(self):
        """Return mass number of most abundant natural stable isotope."""
        nominalmass = 0
        maxabundance = 0
        for massnum, iso in self.isotopes.items():
            if iso.abundance > maxabundance:
                maxabundance = iso.abundance
                nominalmass = massnum
        return nominalmass

    @lazyattr
    def neutrons(self):
        """Return number neutrons in most abundant natural stable isotope."""
        return self.nominalmass - self.protons

    @lazyattr
    def exactmass(self):
        """Return relative atomic mass calculated from isotopic composition."""
        return sum(iso.mass * iso.abundance for iso in self.isotopes.values())

    @lazyattr
    def eleconfig_dict(self):
        """Return electron configuration as dict."""
        adict = {}
        if self.eleconfig.startswith('['):
            base = self.eleconfig.split(' ', 1)[0][1:-1]
            adict.update(ELEMENTS[base].eleconfig_dict)
        for e in self.eleconfig.split()[bool(adict):]:
            adict[(int(e[0]), e[1])] = int(e[2:]) if len(e) > 2 else 1
        return adict

    @lazyattr
    def eleshells(self):
        """Return number of electrons in shell as tuple."""
        eleshells = [0, 0, 0, 0, 0, 0, 0]
        for key, val in self.eleconfig_dict.items():
            eleshells[key[0] - 1] += val
        return tuple(ele for ele in eleshells if ele)

    @lazyattr
    def description(self):
        """Return text description of element."""
        return _descriptions(self.symbol)

    def validate(self):
        """Check consistency of data. Raise Error on failure."""
        assert self.period in PERIODS
        assert self.group in GROUPS
        assert self.block in BLOCKS
        assert self.series in SERIES

        if self.number != self.protons:
            raise ValueError(
                "%s - atomic number must equal proton number" % self.symbol)
        if self.protons != sum(self.eleshells):
            raise ValueError(
                "%s - number of protons must equal electrons" % self.symbol)

        mass = 0.0
        frac = 0.0
        for iso in self.isotopes.values():
            mass += iso.abundance * iso.mass
            frac += iso.abundance
        if abs(mass - self.mass) > 0.03:
            raise ValueError(
                "%s - average of isotope masses (%.4f) != mass (%.4f)" % (
                    self.symbol, mass, self.mass))
        if abs(frac - 1.0) > 1e-9:
            raise ValueError(
                "%s - sum of isotope abundances != 1.0" % self.symbol)


class Isotope(object):
    """Isotope massnumber, relative atomic mass, and abundance."""
    __slots__ = ['massnumber', 'mass', 'abundance']

    def __init__(self, mass=0.0, abundance=1.0, massnumber=0):
        self.mass = mass
        self.abundance = abundance
        self.massnumber = massnumber

    def __str__(self):
        return "%i, %.4f, %.6f%%" % (self.massnumber, self.mass,
                                     self.abundance * 100)

    def __repr__(self):
        return "Isotope(%s, %s, %s)" % (
            repr(self.mass), repr(self.abundance), repr(self.massnumber))


class ElementsDict(object):
    """Ordered dict of Elements with lookup by number, symbol, and name."""
    def __init__(self, *elements):
        self._list = []
        self._dict = {}
        for element in elements:
            if element.number > len(self._list) + 1:
                raise ValueError("Elements must be added in order")
            if element.number <= len(self._list):
                self._list[element.number - 1] = element
            else:
                self._list.append(element)
            self._dict[element.number] = element
            self._dict[element.symbol] = element
            self._dict[element.name] = element

    def __str__(self):
        return "[%s]" % ", ".join(ele.symbol for ele in self._list)

    def __contains__(self, item):
        return item in self._dict

    def __iter__(self):
        return iter(self._list)

    def __len__(self):
        return len(self._list)

    def __getitem__(self, key):
        try:
            return self._dict[key]
        except KeyError:
            try:
                start, stop, step = key.indices(len(self._list))
                return self._list[slice(start - 1, stop - 1, step)]
            except:
                raise KeyError("key %s not found" % key)


ELEMENTS = ElementsDict(
    Element(
        1, 'H', 'Hydrogen',
        group=1, period=1, block='s', series=1,
        mass=1.00794, eleneg=2.2, eleaffin=0.75420375,
        covrad=0.32, atmrad=0.79, vdwrad=1.2,
        tboil=20.28, tmelt=13.81, density=0.084,
        eleconfig='1s',
        oxistates='1*, -1',
        ionenergy=(13.5984, ),
        isotopes={1: Isotope(1.0078250321, 0.999885, 1),
                  2: Isotope(2.014101778, 0.000115, 2)}),
    Element(
        2, 'He', 'Helium',
        group=18, period=1, block='s', series=2,
        mass=4.002602, eleneg=0.0, eleaffin=0.0,
        covrad=0.93, atmrad=0.49, vdwrad=1.4,
        tboil=4.216, tmelt=0.95, density=0.1785,
        eleconfig='1s2',
        oxistates='*',
        ionenergy=(24.5874, 54.416, ),
        isotopes={3: Isotope(3.0160293097, 1.37e-06, 3),
                  4: Isotope(4.0026032497, 0.99999863, 4)}),
    Element(
        3, 'Li', 'Lithium',
        group=1, period=2, block='s', series=3,
        mass=6.941, eleneg=0.98, eleaffin=0.618049,
        covrad=1.23, atmrad=2.05, vdwrad=1.82,
        tboil=1615.0, tmelt=453.7, density=0.53,
        eleconfig='[He] 2s',
        oxistates='1*',
        ionenergy=(5.3917, 75.638, 122.451, ),
        isotopes={6: Isotope(6.0151223, 0.0759, 6),
                  7: Isotope(7.016004, 0.9241, 7)}),
    Element(
        4, 'Be', 'Beryllium',
        group=2, period=2, block='s', series=4,
        mass=9.012182, eleneg=1.57, eleaffin=0.0,
        covrad=0.9, atmrad=1.4, vdwrad=0.0,
        tboil=3243.0, tmelt=1560.0, density=1.85,
        eleconfig='[He] 2s2',
        oxistates='2*',
        ionenergy=(9.3227, 18.211, 153.893, 217.713, ),
        isotopes={9: Isotope(9.0121821, 1.0, 9)}),
    Element(
        5, 'B', 'Boron',
        group=13, period=2, block='p', series=5,
        mass=10.811, eleneg=2.04, eleaffin=0.279723,
        covrad=0.82, atmrad=1.17, vdwrad=0.0,
        tboil=4275.0, tmelt=2365.0, density=2.46,
        eleconfig='[He] 2s2 2p',
        oxistates='3*',
        ionenergy=(8.298, 25.154, 37.93, 59.368, 340.217, ),
        isotopes={10: Isotope(10.012937, 0.199, 10),
                  11: Isotope(11.0093055, 0.801, 11)}),
    Element(
        6, 'C', 'Carbon',
        group=14, period=2, block='p', series=1,
        mass=12.0107, eleneg=2.55, eleaffin=1.262118,
        covrad=0.77, atmrad=0.91, vdwrad=1.7,
        tboil=5100.0, tmelt=3825.0, density=3.51,
        eleconfig='[He] 2s2 2p2',
        oxistates='4*, 2, -4*',
        ionenergy=(11.2603, 24.383, 47.877, 64.492, 392.077,
                   489.981, ),
        isotopes={12: Isotope(12.0, 0.9893, 12),
                  13: Isotope(13.0033548378, 0.0107, 13)}),
    Element(
        7, 'N', 'Nitrogen',
        group=15, period=2, block='p', series=1,
        mass=14.0067, eleneg=3.04, eleaffin=-0.07,
        covrad=0.75, atmrad=0.75, vdwrad=1.55,
        tboil=77.344, tmelt=63.15, density=1.17,
        eleconfig='[He] 2s2 2p3',
        oxistates='5, 4, 3, 2, -3*',
        ionenergy=(14.5341, 39.601, 47.488, 77.472, 97.888,
                   522.057, 667.029, ),
        isotopes={14: Isotope(14.0030740052, 0.99632, 14),
                  15: Isotope(15.0001088984, 0.00368, 15)}),
    Element(
        8, 'O', 'Oxygen',
        group=16, period=2, block='p', series=1,
        mass=15.9994, eleneg=3.44, eleaffin=1.461112,
        covrad=0.73, atmrad=0.65, vdwrad=1.52,
        tboil=90.188, tmelt=54.8, density=1.33,
        eleconfig='[He] 2s2 2p4',
        oxistates='-2*, -1',
        ionenergy=(13.6181, 35.116, 54.934, 54.934, 77.412,
                   113.896, 138.116, 739.315, 871.387, ),
        isotopes={16: Isotope(15.9949146221, 0.99757, 16),
                  17: Isotope(16.9991315, 0.00038, 17),
                  18: Isotope(17.9991604, 0.00205, 18)}),
    Element(
        9, 'F', 'Fluorine',
        group=17, period=2, block='p', series=6,
        mass=18.9984032, eleneg=3.98, eleaffin=3.4011887,
        covrad=0.72, atmrad=0.57, vdwrad=1.47,
        tboil=85.0, tmelt=53.55, density=1.58,
        eleconfig='[He] 2s2 2p5',
        oxistates='-1*',
        ionenergy=(17.4228, 34.97, 62.707, 87.138, 114.24,
                   157.161, 185.182, 953.886, 1103.089, ),
        isotopes={19: Isotope(18.9984032, 1.0, 19)}),
    Element(
        10, 'Ne', 'Neon',
        group=18, period=2, block='p', series=2,
        mass=20.1797, eleneg=0.0, eleaffin=0.0,
        covrad=0.71, atmrad=0.51, vdwrad=1.54,
        tboil=27.1, tmelt=24.55, density=0.8999,
        eleconfig='[He] 2s2 2p6',
        oxistates='*',
        ionenergy=(21.5645, 40.962, 63.45, 97.11, 126.21,
                   157.93, 207.27, 239.09, 1195.797, 1362.164, ),
        isotopes={20: Isotope(19.9924401759, 0.9048, 20),
                  21: Isotope(20.99384674, 0.0027, 21),
                  22: Isotope(21.99138551, 0.0925, 22)}),
    Element(
        11, 'Na', 'Sodium',
        group=1, period=3, block='s', series=3,
        mass=22.98977, eleneg=0.93, eleaffin=0.547926,
        covrad=1.54, atmrad=2.23, vdwrad=2.27,
        tboil=1156.0, tmelt=371.0, density=0.97,
        eleconfig='[Ne] 3s',
        oxistates='1*',
        ionenergy=(5.1391, 47.286, 71.64, 98.91, 138.39,
                   172.15, 208.47, 264.18, 299.87, 1465.091,
                   1648.659, ),
        isotopes={23: Isotope(22.98976967, 1.0, 23)}),
    Element(
        12, 'Mg', 'Magnesium',
        group=2, period=3, block='s', series=4,
        mass=24.305, eleneg=1.31, eleaffin=0.0,
        covrad=1.36, atmrad=1.72, vdwrad=1.73,
        tboil=1380.0, tmelt=922.0, density=1.74,
        eleconfig='[Ne] 3s2',
        oxistates='2*',
        ionenergy=(7.6462, 15.035, 80.143, 109.24, 141.26,
                   186.5, 224.94, 265.9, 327.95, 367.53,
                   1761.802, 1962.613, ),
        isotopes={24: Isotope(23.9850419, 0.7899, 24),
                  25: Isotope(24.98583702, 0.1, 25),
                  26: Isotope(25.98259304, 0.1101, 26)}),
    Element(
        13, 'Al', 'Aluminium',
        group=13, period=3, block='p', series=7,
        mass=26.981538, eleneg=1.61, eleaffin=0.43283,
        covrad=1.18, atmrad=1.82, vdwrad=0.0,
        tboil=2740.0, tmelt=933.5, density=2.7,
        eleconfig='[Ne] 3s2 3p',
        oxistates='3*',
        ionenergy=(5.9858, 18.828, 28.447, 119.99, 153.71,
                   190.47, 241.43, 284.59, 330.21, 398.57,
                   442.07, 2085.983, 2304.08, ),
        isotopes={27: Isotope(26.98153844, 1.0, 27)}),
    Element(
        14, 'Si', 'Silicon',
        group=14, period=3, block='p', series=5,
        mass=28.0855, eleneg=1.9, eleaffin=1.389521,
        covrad=1.11, atmrad=1.46, vdwrad=2.1,
        tboil=2630.0, tmelt=1683.0, density=2.33,
        eleconfig='[Ne] 3s2 3p2',
        oxistates='4*, -4',
        ionenergy=(8.1517, 16.345, 33.492, 45.141, 166.77,
                   205.05, 246.52, 303.17, 351.1, 401.43,
                   476.06, 523.5, 2437.676, 2673.108, ),
        isotopes={28: Isotope(27.9769265327, 0.922297, 28),
                  29: Isotope(28.97649472, 0.046832, 29),
                  30: Isotope(29.97377022, 0.030871, 30)}),
    Element(
        15, 'P', 'Phosphorus',
        group=15, period=3, block='p', series=1,
        mass=30.973761, eleneg=2.19, eleaffin=0.7465,
        covrad=1.06, atmrad=1.23, vdwrad=1.8,
        tboil=553.0, tmelt=317.3, density=1.82,
        eleconfig='[Ne] 3s2 3p3',
        oxistates='5*, 3, -3',
        ionenergy=(10.4867, 19.725, 30.18, 51.37, 65.023,
                   220.43, 263.22, 309.41, 371.73, 424.5,
                   479.57, 560.41, 611.85, 2816.943, 3069.762, ),
        isotopes={31: Isotope(30.97376151, 1.0, 31)}),
    Element(
        16, 'S', 'Sulfur',
        group=16, period=3, block='p', series=1,
        mass=32.065, eleneg=2.58, eleaffin=2.0771029,
        covrad=1.02, atmrad=1.09, vdwrad=1.8,
        tboil=717.82, tmelt=392.2, density=2.06,
        eleconfig='[Ne] 3s2 3p4',
        oxistates='6*, 4, 2, -2',
        ionenergy=(10.36, 23.33, 34.83, 47.3, 72.68,
                   88.049, 280.93, 328.23, 379.1, 447.09,
                   504.78, 564.65, 651.63, 707.14, 3223.836,
                   3494.099, ),
        isotopes={32: Isotope(31.97207069, 0.9493, 32),
                  33: Isotope(32.9714585, 0.0076, 33),
                  34: Isotope(33.96786683, 0.0429, 34),
                  36: Isotope(35.96708088, 0.0002, 36)}),
    Element(
        17, 'Cl', 'Chlorine',
        group=17, period=3, block='p', series=6,
        mass=35.453, eleneg=3.16, eleaffin=3.612724,
        covrad=0.99, atmrad=0.97, vdwrad=1.75,
        tboil=239.18, tmelt=172.17, density=2.95,
        eleconfig='[Ne] 3s2 3p5',
        oxistates='7, 5, 3, 1, -1*',
        ionenergy=(12.9676, 23.81, 39.61, 53.46, 67.8,
                   98.03, 114.193, 348.28, 400.05, 455.62,
                   529.97, 591.97, 656.69, 749.75, 809.39,
                   3658.425, 3946.193, ),
        isotopes={35: Isotope(34.96885271, 0.7578, 35),
                  37: Isotope(36.9659026, 0.2422, 37)}),
    Element(
        18, 'Ar', 'Argon',
        group=18, period=3, block='p', series=2,
        mass=39.948, eleneg=0.0, eleaffin=0.0,
        covrad=0.98, atmrad=0.88, vdwrad=1.88,
        tboil=87.45, tmelt=83.95, density=1.66,
        eleconfig='[Ne] 3s2 3p6',
        oxistates='*',
        ionenergy=(15.7596, 27.629, 40.74, 59.81, 75.02,
                   91.007, 124.319, 143.456, 422.44, 478.68,
                   538.95, 618.24, 686.09, 755.73, 854.75,
                   918.0, 4120.778, 4426.114, ),
        isotopes={36: Isotope(35.96754628, 0.003365, 36),
                  38: Isotope(37.9627322, 0.000632, 38),
                  40: Isotope(39.962383123, 0.996003, 40)}),
    Element(
        19, 'K', 'Potassium',
        group=1, period=4, block='s', series=3,
        mass=39.0983, eleneg=0.82, eleaffin=0.501459,
        covrad=2.03, atmrad=2.77, vdwrad=2.75,
        tboil=1033.0, tmelt=336.8, density=0.86,
        eleconfig='[Ar] 4s',
        oxistates='1*',
        ionenergy=(4.3407, 31.625, 45.72, 60.91, 82.66,
                   100.0, 117.56, 154.86, 175.814, 503.44,
                   564.13, 629.09, 714.02, 787.13, 861.77,
                   968.0, 1034.0, 4610.955, 4933.931, ),
        isotopes={39: Isotope(38.9637069, 0.932581, 39),
                  40: Isotope(39.96399867, 0.000117, 40),
                  41: Isotope(40.96182597, 0.067302, 41)}),
    Element(
        20, 'Ca', 'Calcium',
        group=2, period=4, block='s', series=4,
        mass=40.078, eleneg=1.0, eleaffin=0.02455,
        covrad=1.74, atmrad=2.23, vdwrad=0.0,
        tboil=1757.0, tmelt=1112.0, density=1.54,
        eleconfig='[Ar] 4s2',
        oxistates='2*',
        ionenergy=(6.1132, 11.71, 50.908, 67.1, 84.41,
                   108.78, 127.7, 147.24, 188.54, 211.27,
                   591.25, 656.39, 726.03, 816.61, 895.12,
                   974.0, 1087.0, 1157.0, 5129.045, 5469.738, ),
        isotopes={40: Isotope(39.9625912, 0.96941, 40),
                  42: Isotope(41.9586183, 0.00647, 42),
                  43: Isotope(42.9587668, 0.00135, 43),
                  44: Isotope(43.9554811, 0.02086, 44),
                  46: Isotope(45.9536928, 4e-05, 46),
                  48: Isotope(47.952534, 0.00187, 48)}),
    Element(
        21, 'Sc', 'Scandium',
        group=3, period=4, block='d', series=8,
        mass=44.95591, eleneg=1.36, eleaffin=0.188,
        covrad=1.44, atmrad=2.09, vdwrad=0.0,
        tboil=3109.0, tmelt=1814.0, density=2.99,
        eleconfig='[Ar] 3d 4s2',
        oxistates='3*',
        ionenergy=(6.5615, 12.8, 24.76, 73.47, 91.66,
                   11.1, 138.0, 158.7, 180.02, 225.32,
                   225.32, 685.89, 755.47, 829.79, 926.0, ),
        isotopes={45: Isotope(44.9559102, 1.0, 45)}),
    Element(
        22, 'Ti', 'Titanium',
        group=4, period=4, block='d', series=8,
        mass=47.867, eleneg=1.54, eleaffin=0.084,
        covrad=1.32, atmrad=2.0, vdwrad=0.0,
        tboil=3560.0, tmelt=1935.0, density=4.51,
        eleconfig='[Ar] 3d2 4s2',
        oxistates='4*, 3',
        ionenergy=(6.8281, 13.58, 27.491, 43.266, 99.22,
                   119.36, 140.8, 168.5, 193.5, 193.2,
                   215.91, 265.23, 291.497, 787.33, 861.33, ),
        isotopes={46: Isotope(45.9526295, 0.0825, 46),
                  47: Isotope(46.9517638, 0.0744, 47),
                  48: Isotope(47.9479471, 0.7372, 48),
                  49: Isotope(48.9478708, 0.0541, 49),
                  50: Isotope(49.9447921, 0.0518, 50)}),
    Element(
        23, 'V', 'Vanadium',
        group=5, period=4, block='d', series=8,
        mass=50.9415, eleneg=1.63, eleaffin=0.525,
        covrad=1.22, atmrad=1.92, vdwrad=0.0,
        tboil=3650.0, tmelt=2163.0, density=6.09,
        eleconfig='[Ar] 3d3 4s2',
        oxistates='5*, 4, 3, 2, 0',
        ionenergy=(6.7462, 14.65, 29.31, 46.707, 65.23,
                   128.12, 150.17, 173.7, 205.8, 230.5,
                   255.04, 308.25, 336.267, 895.58, 974.02, ),
        isotopes={50: Isotope(49.9471628, 0.0025, 50),
                  51: Isotope(50.9439637, 0.9975, 51)}),
    Element(
        24, 'Cr', 'Chromium',
        group=6, period=4, block='d', series=8,
        mass=51.9961, eleneg=1.66, eleaffin=0.67584,
        covrad=1.18, atmrad=1.85, vdwrad=0.0,
        tboil=2945.0, tmelt=2130.0, density=7.14,
        eleconfig='[Ar] 3d5 4s',
        oxistates='6, 3*, 2, 0',
        ionenergy=(6.7665, 16.5, 30.96, 49.1, 69.3,
                   90.56, 161.1, 184.7, 209.3, 244.4,
                   270.8, 298.0, 355.0, 384.3, 1010.64, ),
        isotopes={50: Isotope(49.9460496, 0.04345, 50),
                  52: Isotope(51.9405119, 0.83789, 52),
                  53: Isotope(52.9406538, 0.09501, 53),
                  54: Isotope(53.9388849, 0.02365, 54)}),
    Element(
        25, 'Mn', 'Manganese',
        group=7, period=4, block='d', series=8,
        mass=54.938049, eleneg=1.55, eleaffin=0.0,
        covrad=1.17, atmrad=1.79, vdwrad=0.0,
        tboil=2235.0, tmelt=1518.0, density=7.44,
        eleconfig='[Ar] 3d5 4s2',
        oxistates='7, 6, 4, 3, 2*, 0, -1',
        ionenergy=(7.434, 15.64, 33.667, 51.2, 72.4,
                   95.0, 119.27, 196.46, 221.8, 248.3,
                   286.0, 314.4, 343.6, 404.0, 435.3,
                   1136.2, ),
        isotopes={55: Isotope(54.9380496, 1.0, 55)}),
    Element(
        26, 'Fe', 'Iron',
        group=8, period=4, block='d', series=8,
        mass=55.845, eleneg=1.83, eleaffin=0.151,
        covrad=1.17, atmrad=1.72, vdwrad=0.0,
        tboil=3023.0, tmelt=1808.0, density=7.874,
        eleconfig='[Ar] 3d6 4s2',
        oxistates='6, 3*, 2, 0, -2',
        ionenergy=(7.9024, 16.18, 30.651, 54.8, 75.0,
                   99.0, 125.0, 151.06, 235.04, 262.1,
                   290.4, 330.8, 361.0, 392.2, 457.0,
                   485.5, 1266.1, ),
        isotopes={54: Isotope(53.9396148, 0.05845, 54),
                  56: Isotope(55.9349421, 0.91754, 56),
                  57: Isotope(56.9353987, 0.02119, 57),
                  58: Isotope(57.9332805, 0.00282, 58)}),
    Element(
        27, 'Co', 'Cobalt',
        group=9, period=4, block='d', series=8,
        mass=58.9332, eleneg=1.88, eleaffin=0.6633,
        covrad=1.16, atmrad=1.67, vdwrad=0.0,
        tboil=3143.0, tmelt=1768.0, density=8.89,
        eleconfig='[Ar] 3d7 4s2',
        oxistates='3, 2*, 0, -1',
        ionenergy=(7.881, 17.06, 33.5, 51.3, 79.5,
                   102.0, 129.0, 157.0, 186.13, 276.0,
                   305.0, 336.0, 376.0, 411.0, 444.0,
                   512.0, 546.8, 1403.0, ),
        isotopes={59: Isotope(58.9332002, 1.0, 59)}),
    Element(
        28, 'Ni', 'Nickel',
        group=10, period=4, block='d', series=8,
        mass=58.6934, eleneg=1.91, eleaffin=1.15716,
        covrad=1.15, atmrad=1.62, vdwrad=1.63,
        tboil=3005.0, tmelt=1726.0, density=8.91,
        eleconfig='[Ar] 3d8 4s2',
        oxistates='3, 2*, 0',
        ionenergy=(7.6398, 18.168, 35.17, 54.9, 75.5,
                   108.0, 133.0, 162.0, 193.0, 224.5,
                   321.2, 352.0, 384.0, 430.0, 464.0,
                   499.0, 571.0, 607.2, 1547.0, ),
        isotopes={58: Isotope(57.9353479, 0.680769, 58),
                  60: Isotope(59.9307906, 0.262231, 60),
                  61: Isotope(60.9310604, 0.011399, 61),
                  62: Isotope(61.9283488, 0.036345, 62),
                  64: Isotope(63.9279696, 0.009256, 64)}),
    Element(
        29, 'Cu', 'Copper',
        group=11, period=4, block='d', series=8,
        mass=63.546, eleneg=1.9, eleaffin=1.23578,
        covrad=1.17, atmrad=1.57, vdwrad=1.4,
        tboil=2840.0, tmelt=1356.6, density=8.92,
        eleconfig='[Ar] 3d10 4s',
        oxistates='2*, 1',
        ionenergy=(7.7264, 20.292, 26.83, 55.2, 79.9,
                   103.0, 139.0, 166.0, 199.0, 232.0,
                   266.0, 368.8, 401.0, 435.0, 484.0,
                   520.0, 557.0, 633.0, 671.0, 1698.0, ),
        isotopes={63: Isotope(62.9296011, 0.6917, 63),
                  65: Isotope(64.9277937, 0.3083, 65)}),
    Element(
        30, 'Zn', 'Zinc',
        group=12, period=4, block='d', series=8,
        mass=65.409, eleneg=1.65, eleaffin=0.0,
        covrad=1.25, atmrad=1.53, vdwrad=1.39,
        tboil=1180.0, tmelt=692.73, density=7.14,
        eleconfig='[Ar] 3d10 4s2',
        oxistates='2*',
        ionenergy=(9.3942, 17.964, 39.722, 59.4, 82.6,
                   108.0, 134.0, 174.0, 203.0, 238.0,
                   274.0, 310.8, 419.7, 454.0, 490.0,
                   542.0, 579.0, 619.0, 698.8, 738.0,
                   1856.0, ),
        isotopes={64: Isotope(63.9291466, 0.4863, 64),
                  66: Isotope(65.9260368, 0.279, 66),
                  67: Isotope(66.9271309, 0.041, 67),
                  68: Isotope(67.9248476, 0.1875, 68),
                  70: Isotope(69.925325, 0.0062, 70)}),
    Element(
        31, 'Ga', 'Gallium',
        group=13, period=4, block='p', series=7,
        mass=69.723, eleneg=1.81, eleaffin=0.41,
        covrad=1.26, atmrad=1.81, vdwrad=1.87,
        tboil=2478.0, tmelt=302.92, density=5.91,
        eleconfig='[Ar] 3d10 4s2 4p',
        oxistates='3*',
        ionenergy=(5.9993, 20.51, 30.71, 64.0, ),
        isotopes={69: Isotope(68.925581, 0.60108, 69),
                  71: Isotope(70.924705, 0.39892, 71)}),
    Element(
        32, 'Ge', 'Germanium',
        group=14, period=4, block='p', series=5,
        mass=72.64, eleneg=2.01, eleaffin=1.232712,
        covrad=1.22, atmrad=1.52, vdwrad=0.0,
        tboil=3107.0, tmelt=1211.5, density=5.32,
        eleconfig='[Ar] 3d10 4s2 4p2',
        oxistates='4*',
        ionenergy=(7.8994, 15.934, 34.22, 45.71, 93.5, ),
        isotopes={70: Isotope(69.9242504, 0.2084, 70),
                  72: Isotope(71.9220762, 0.2754, 72),
                  73: Isotope(72.9234594, 0.0773, 73),
                  74: Isotope(73.9211782, 0.3628, 74),
                  76: Isotope(75.9214027, 0.0761, 76)}),
    Element(
        33, 'As', 'Arsenic',
        group=15, period=4, block='p', series=5,
        mass=74.9216, eleneg=2.18, eleaffin=0.814,
        covrad=1.2, atmrad=1.33, vdwrad=1.85,
        tboil=876.0, tmelt=1090.0, density=5.72,
        eleconfig='[Ar] 3d10 4s2 4p3',
        oxistates='5, 3*, -3',
        ionenergy=(9.7886, 18.633, 28.351, 50.13, 62.63,
                   127.6, ),
        isotopes={75: Isotope(74.9215964, 1.0, 75)}),
    Element(
        34, 'Se', 'Selenium',
        group=16, period=4, block='p', series=1,
        mass=78.96, eleneg=2.55, eleaffin=2.02067,
        covrad=1.16, atmrad=1.22, vdwrad=1.9,
        tboil=958.0, tmelt=494.0, density=4.82,
        eleconfig='[Ar] 3d10 4s2 4p4',
        oxistates='6, 4*, -2',
        ionenergy=(9.7524, 21.9, 30.82, 42.944, 68.3,
                   81.7, 155.4, ),
        isotopes={74: Isotope(73.9224766, 0.0089, 74),
                  76: Isotope(75.9192141, 0.0937, 76),
                  77: Isotope(76.9199146, 0.0763, 77),
                  78: Isotope(77.9173095, 0.2377, 78),
                  80: Isotope(79.9165218, 0.4961, 80),
                  82: Isotope(81.9167, 0.0873, 82)}),
    Element(
        35, 'Br', 'Bromine',
        group=17, period=4, block='p', series=6,
        mass=79.904, eleneg=2.96, eleaffin=3.363588,
        covrad=1.14, atmrad=1.12, vdwrad=1.85,
        tboil=331.85, tmelt=265.95, density=3.14,
        eleconfig='[Ar] 3d10 4s2 4p5',
        oxistates='7, 5, 3, 1, -1*',
        ionenergy=(11.8138, 21.8, 36.0, 47.3, 59.7,
                   88.6, 103.0, 192.8, ),
        isotopes={79: Isotope(78.9183376, 0.5069, 79),
                  81: Isotope(80.916291, 0.4931, 81)}),
    Element(
        36, 'Kr', 'Krypton',
        group=18, period=4, block='p', series=2,
        mass=83.798, eleneg=0.0, eleaffin=0.0,
        covrad=1.12, atmrad=1.03, vdwrad=2.02,
        tboil=120.85, tmelt=116.0, density=4.48,
        eleconfig='[Ar] 3d10 4s2 4p6',
        oxistates='2*',
        ionenergy=(13.9996, 24.359, 36.95, 52.5, 64.7,
                   78.5, 110.0, 126.0, 230.39, ),
        isotopes={78: Isotope(77.920386, 0.0035, 78),
                  80: Isotope(79.916378, 0.0228, 80),
                  82: Isotope(81.9134846, 0.1158, 82),
                  83: Isotope(82.914136, 0.1149, 83),
                  84: Isotope(83.911507, 0.57, 84),
                  86: Isotope(85.9106103, 0.173, 86)}),
    Element(
        37, 'Rb', 'Rubidium',
        group=1, period=5, block='s', series=3,
        mass=85.4678, eleneg=0.82, eleaffin=0.485916,
        covrad=2.16, atmrad=2.98, vdwrad=0.0,
        tboil=961.0, tmelt=312.63, density=1.53,
        eleconfig='[Kr] 5s',
        oxistates='1*',
        ionenergy=(4.1771, 27.28, 40.0, 52.6, 71.0,
                   84.4, 99.2, 136.0, 150.0, 277.1, ),
        isotopes={85: Isotope(84.9117893, 0.7217, 85),
                  87: Isotope(86.9091835, 0.2783, 87)}),
    Element(
        38, 'Sr', 'Strontium',
        group=2, period=5, block='s', series=4,
        mass=87.62, eleneg=0.95, eleaffin=0.05206,
        covrad=1.91, atmrad=2.45, vdwrad=0.0,
        tboil=1655.0, tmelt=1042.0, density=2.63,
        eleconfig='[Kr] 5s2',
        oxistates='2*',
        ionenergy=(5.6949, 11.03, 43.6, 57.0, 71.6,
                   90.8, 106.0, 122.3, 162.0, 177.0,
                   324.1, ),
        isotopes={84: Isotope(83.913425, 0.0056, 84),
                  86: Isotope(85.9092624, 0.0986, 86),
                  87: Isotope(86.9088793, 0.07, 87),
                  88: Isotope(87.9056143, 0.8258, 88)}),
    Element(
        39, 'Y', 'Yttrium',
        group=3, period=5, block='d', series=8,
        mass=88.90585, eleneg=1.22, eleaffin=0.307,
        covrad=1.62, atmrad=2.27, vdwrad=0.0,
        tboil=3611.0, tmelt=1795.0, density=4.47,
        eleconfig='[Kr] 4d 5s2',
        oxistates='3*',
        ionenergy=(6.2173, 12.24, 20.52, 61.8, 77.0,
                   93.0, 116.0, 129.0, 146.52, 191.0,
                   206.0, 374.0, ),
        isotopes={89: Isotope(88.9058479, 1.0, 89)}),
    Element(
        40, 'Zr', 'Zirconium',
        group=4, period=5, block='d', series=8,
        mass=91.224, eleneg=1.33, eleaffin=0.426,
        covrad=1.45, atmrad=2.16, vdwrad=0.0,
        tboil=4682.0, tmelt=2128.0, density=6.51,
        eleconfig='[Kr] 4d2 5s2',
        oxistates='4*',
        ionenergy=(6.6339, 13.13, 22.99, 34.34, 81.5, ),
        isotopes={90: Isotope(89.9047037, 0.5145, 90),
                  91: Isotope(90.905645, 0.1122, 91),
                  92: Isotope(91.9050401, 0.1715, 92),
                  94: Isotope(93.9063158, 0.1738, 94),
                  96: Isotope(95.908276, 0.028, 96)}),
    Element(
        41, 'Nb', 'Niobium',
        group=5, period=5, block='d', series=8,
        mass=92.90638, eleneg=1.6, eleaffin=0.893,
        covrad=1.34, atmrad=2.08, vdwrad=0.0,
        tboil=5015.0, tmelt=2742.0, density=8.58,
        eleconfig='[Kr] 4d4 5s',
        oxistates='5*, 3',
        ionenergy=(6.7589, 14.32, 25.04, 38.3, 50.55,
                   102.6, 125.0, ),
        isotopes={93: Isotope(92.9063775, 1.0, 93)}),
    Element(
        42, 'Mo', 'Molybdenum',
        group=6, period=5, block='d', series=8,
        mass=95.94, eleneg=2.16, eleaffin=0.7472,
        covrad=1.3, atmrad=2.01, vdwrad=0.0,
        tboil=4912.0, tmelt=2896.0, density=10.28,
        eleconfig='[Kr] 4d5 5s',
        oxistates='6*, 5, 4, 3, 2, 0',
        ionenergy=(7.0924, 16.15, 27.16, 46.4, 61.2,
                   68.0, 126.8, 153.0, ),
        isotopes={92: Isotope(91.90681, 0.1484, 92),
                  94: Isotope(93.9050876, 0.0925, 94),
                  95: Isotope(94.9058415, 0.1592, 95),
                  96: Isotope(95.9046789, 0.1668, 96),
                  97: Isotope(96.906021, 0.0955, 97),
                  98: Isotope(97.9054078, 0.2413, 98),
                  100: Isotope(99.907477, 0.0963, 100)}),
    Element(
        43, 'Tc', 'Technetium',
        group=7, period=5, block='d', series=8,
        mass=97.907216, eleneg=1.9, eleaffin=0.55,
        covrad=1.27, atmrad=1.95, vdwrad=0.0,
        tboil=4538.0, tmelt=2477.0, density=11.49,
        eleconfig='[Kr] 4d5 5s2',
        oxistates='7*',
        ionenergy=(7.28, 15.26, 29.54, ),
        isotopes={98: Isotope(97.907216, 1.0, 98)}),
    Element(
        44, 'Ru', 'Ruthenium',
        group=8, period=5, block='d', series=8,
        mass=101.07, eleneg=2.2, eleaffin=1.04638,
        covrad=1.25, atmrad=1.89, vdwrad=0.0,
        tboil=4425.0, tmelt=2610.0, density=12.45,
        eleconfig='[Kr] 4d7 5s',
        oxistates='8, 6, 4*, 3*, 2, 0, -2',
        ionenergy=(7.3605, 16.76, 28.47, ),
        isotopes={96: Isotope(95.907598, 0.0554, 96),
                  98: Isotope(97.905287, 0.0187, 98),
                  99: Isotope(98.9059393, 0.1276, 99),
                  100: Isotope(99.9042197, 0.126, 100),
                  101: Isotope(100.9055822, 0.1706, 101),
                  102: Isotope(101.9043495, 0.3155, 102),
                  104: Isotope(103.90543, 0.1862, 104)}),
    Element(
        45, 'Rh', 'Rhodium',
        group=9, period=5, block='d', series=8,
        mass=102.9055, eleneg=2.28, eleaffin=1.14289,
        covrad=1.25, atmrad=1.83, vdwrad=0.0,
        tboil=3970.0, tmelt=2236.0, density=12.41,
        eleconfig='[Kr] 4d8 5s',
        oxistates='5, 4, 3*, 1*, 2, 0',
        ionenergy=(7.4589, 18.08, 31.06, ),
        isotopes={103: Isotope(102.905504, 1.0, 103)}),
    Element(
        46, 'Pd', 'Palladium',
        group=10, period=5, block='d', series=8,
        mass=106.42, eleneg=2.2, eleaffin=0.56214,
        covrad=1.28, atmrad=1.79, vdwrad=1.63,
        tboil=3240.0, tmelt=1825.0, density=12.02,
        eleconfig='[Kr] 4d10',
        oxistates='4, 2*, 0',
        ionenergy=(8.3369, 19.43, 32.93, ),
        isotopes={102: Isotope(101.905608, 0.0102, 102),
                  104: Isotope(103.904035, 0.1114, 104),
                  105: Isotope(104.905084, 0.2233, 105),
                  106: Isotope(105.903483, 0.2733, 106),
                  108: Isotope(107.903894, 0.2646, 108),
                  110: Isotope(109.905152, 0.1172, 110)}),
    Element(
        47, 'Ag', 'Silver',
        group=11, period=5, block='d', series=8,
        mass=107.8682, eleneg=1.93, eleaffin=1.30447,
        covrad=1.34, atmrad=1.75, vdwrad=1.72,
        tboil=2436.0, tmelt=1235.1, density=10.49,
        eleconfig='[Kr] 4d10 5s',
        oxistates='2, 1*',
        ionenergy=(7.5762, 21.49, 34.83, ),
        isotopes={107: Isotope(106.905093, 0.51839, 107),
                  109: Isotope(108.904756, 0.48161, 109)}),
    Element(
        48, 'Cd', 'Cadmium',
        group=12, period=5, block='d', series=8,
        mass=112.411, eleneg=1.69, eleaffin=0.0,
        covrad=1.48, atmrad=1.71, vdwrad=1.58,
        tboil=1040.0, tmelt=594.26, density=8.64,
        eleconfig='[Kr] 4d10 5s2',
        oxistates='2*',
        ionenergy=(8.9938, 16.908, 37.48, ),
        isotopes={106: Isotope(105.906458, 0.0125, 106),
                  108: Isotope(107.904183, 0.0089, 108),
                  110: Isotope(109.903006, 0.1249, 110),
                  111: Isotope(110.904182, 0.128, 111),
                  112: Isotope(111.9027572, 0.2413, 112),
                  113: Isotope(112.9044009, 0.1222, 113),
                  114: Isotope(113.9033581, 0.2873, 114),
                  116: Isotope(115.904755, 0.0749, 116)}),
    Element(
        49, 'In', 'Indium',
        group=13, period=5, block='p', series=7,
        mass=114.818, eleneg=1.78, eleaffin=0.404,
        covrad=1.44, atmrad=2.0, vdwrad=1.93,
        tboil=2350.0, tmelt=429.78, density=7.31,
        eleconfig='[Kr] 4d10 5s2 5p',
        oxistates='3*',
        ionenergy=(5.7864, 18.869, 28.03, 28.03, ),
        isotopes={113: Isotope(112.904061, 0.0429, 113),
                  115: Isotope(114.903878, 0.9571, 115)}),
    Element(
        50, 'Sn', 'Tin',
        group=14, period=5, block='p', series=7,
        mass=118.71, eleneg=1.96, eleaffin=1.112066,
        covrad=1.41, atmrad=1.72, vdwrad=2.17,
        tboil=2876.0, tmelt=505.12, density=7.29,
        eleconfig='[Kr] 4d10 5s2 5p2',
        oxistates='4*, 2*',
        ionenergy=(7.3439, 14.632, 30.502, 40.734, 72.28, ),
        isotopes={112: Isotope(111.904821, 0.0097, 112),
                  114: Isotope(113.902782, 0.0066, 114),
                  115: Isotope(114.903346, 0.0034, 115),
                  116: Isotope(115.901744, 0.1454, 116),
                  117: Isotope(116.902954, 0.0768, 117),
                  118: Isotope(117.901606, 0.2422, 118),
                  119: Isotope(118.903309, 0.0859, 119),
                  120: Isotope(119.9021966, 0.3258, 120),
                  122: Isotope(121.9034401, 0.0463, 122),
                  124: Isotope(123.9052746, 0.0579, 124)}),
    Element(
        51, 'Sb', 'Antimony',
        group=15, period=5, block='p', series=5,
        mass=121.76, eleneg=2.05, eleaffin=1.047401,
        covrad=1.4, atmrad=1.53, vdwrad=0.0,
        tboil=1860.0, tmelt=903.91, density=6.69,
        eleconfig='[Kr] 4d10 5s2 5p3',
        oxistates='5, 3*, -3',
        ionenergy=(8.6084, 16.53, 25.3, 44.2, 56.0,
                   108.0, ),
        isotopes={121: Isotope(120.903818, 0.5721, 121),
                  123: Isotope(122.9042157, 0.4279, 123)}),
    Element(
        52, 'Te', 'Tellurium',
        group=16, period=5, block='p', series=5,
        mass=127.6, eleneg=2.1, eleaffin=1.970875,
        covrad=1.36, atmrad=1.42, vdwrad=2.06,
        tboil=1261.0, tmelt=722.72, density=6.25,
        eleconfig='[Kr] 4d10 5s2 5p4',
        oxistates='6, 4*, -2',
        ionenergy=(9.0096, 18.6, 27.96, 37.41, 58.75,
                   70.7, 137.0, ),
        isotopes={120: Isotope(119.90402, 0.0009, 120),
                  122: Isotope(121.9030471, 0.0255, 122),
                  123: Isotope(122.904273, 0.0089, 123),
                  124: Isotope(123.9028195, 0.0474, 124),
                  125: Isotope(124.9044247, 0.0707, 125),
                  126: Isotope(125.9033055, 0.1884, 126),
                  128: Isotope(127.9044614, 0.3174, 128),
                  130: Isotope(129.9062228, 0.3408, 130)}),
    Element(
        53, 'I', 'Iodine',
        group=17, period=5, block='p', series=6,
        mass=126.90447, eleneg=2.66, eleaffin=3.059038,
        covrad=1.33, atmrad=1.32, vdwrad=1.98,
        tboil=457.5, tmelt=386.7, density=4.94,
        eleconfig='[Kr] 4d10 5s2 5p5',
        oxistates='7, 5, 1, -1*',
        ionenergy=(10.4513, 19.131, 33.0, ),
        isotopes={127: Isotope(126.904468, 1.0, 127)}),
    Element(
        54, 'Xe', 'Xenon',
        group=18, period=5, block='p', series=2,
        mass=131.293, eleneg=0.0, eleaffin=0.0,
        covrad=1.31, atmrad=1.24, vdwrad=2.16,
        tboil=165.1, tmelt=161.39, density=4.49,
        eleconfig='[Kr] 4d10 5s2 5p6',
        oxistates='2, 4, 6',
        ionenergy=(12.1298, 21.21, 32.1, ),
        isotopes={124: Isotope(123.9058958, 0.0009, 124),
                  126: Isotope(125.904269, 0.0009, 126),
                  128: Isotope(127.9035304, 0.0192, 128),
                  129: Isotope(128.9047795, 0.2644, 129),
                  130: Isotope(129.9035079, 0.0408, 130),
                  131: Isotope(130.9050819, 0.2118, 131),
                  132: Isotope(131.9041545, 0.2689, 132),
                  134: Isotope(133.9053945, 0.1044, 134),
                  136: Isotope(135.90722, 0.0887, 136)}),
    Element(
        55, 'Cs', 'Caesium',
        group=1, period=6, block='s', series=3,
        mass=132.90545, eleneg=0.79, eleaffin=0.471626,
        covrad=2.35, atmrad=3.34, vdwrad=0.0,
        tboil=944.0, tmelt=301.54, density=1.9,
        eleconfig='[Xe] 6s',
        oxistates='1*',
        ionenergy=(3.8939, 25.1, ),
        isotopes={133: Isotope(132.905447, 1.0, 133)}),
    Element(
        56, 'Ba', 'Barium',
        group=2, period=6, block='s', series=4,
        mass=137.327, eleneg=0.89, eleaffin=0.14462,
        covrad=1.98, atmrad=2.78, vdwrad=0.0,
        tboil=2078.0, tmelt=1002.0, density=3.65,
        eleconfig='[Xe] 6s2',
        oxistates='2*',
        ionenergy=(5.2117, 100.004, ),
        isotopes={130: Isotope(129.90631, 0.00106, 130),
                  132: Isotope(131.905056, 0.00101, 132),
                  134: Isotope(133.904503, 0.02417, 134),
                  135: Isotope(134.905683, 0.06592, 135),
                  136: Isotope(135.90457, 0.07854, 136),
                  137: Isotope(136.905821, 0.11232, 137),
                  138: Isotope(137.905241, 0.71698, 138)}),
    Element(
        57, 'La', 'Lanthanum',
        group=3, period=6, block='f', series=9,
        mass=138.9055, eleneg=1.1, eleaffin=0.47,
        covrad=1.69, atmrad=2.74, vdwrad=0.0,
        tboil=3737.0, tmelt=1191.0, density=6.16,
        eleconfig='[Xe] 5d 6s2',
        oxistates='3*',
        ionenergy=(5.5769, 11.06, 19.175, ),
        isotopes={138: Isotope(137.907107, 0.0009, 138),
                  139: Isotope(138.906348, 0.9991, 139)}),
    Element(
        58, 'Ce', 'Cerium',
        group=3, period=6, block='f', series=9,
        mass=140.116, eleneg=1.12, eleaffin=0.5,
        covrad=1.65, atmrad=2.7, vdwrad=0.0,
        tboil=3715.0, tmelt=1071.0, density=6.77,
        eleconfig='[Xe] 4f 5d 6s2',
        oxistates='4, 3*',
        ionenergy=(5.5387, 10.85, 20.2, 36.72, ),
        isotopes={136: Isotope(135.90714, 0.00185, 136),
                  138: Isotope(137.905986, 0.00251, 138),
                  140: Isotope(139.905434, 0.8845, 140),
                  142: Isotope(141.90924, 0.11114, 142)}),
    Element(
        59, 'Pr', 'Praseodymium',
        group=3, period=6, block='f', series=9,
        mass=140.90765, eleneg=1.13, eleaffin=0.5,
        covrad=1.65, atmrad=2.67, vdwrad=0.0,
        tboil=3785.0, tmelt=1204.0, density=6.48,
        eleconfig='[Xe] 4f3 6s2',
        oxistates='4, 3*',
        ionenergy=(5.473, 10.55, 21.62, 38.95, 57.45, ),
        isotopes={141: Isotope(140.907648, 1.0, 141)}),
    Element(
        60, 'Nd', 'Neodymium',
        group=3, period=6, block='f', series=9,
        mass=144.24, eleneg=1.14, eleaffin=0.5,
        covrad=1.64, atmrad=2.64, vdwrad=0.0,
        tboil=3347.0, tmelt=1294.0, density=7.0,
        eleconfig='[Xe] 4f4 6s2',
        oxistates='3*',
        ionenergy=(5.525, 10.72, ),
        isotopes={142: Isotope(141.907719, 0.272, 142),
                  143: Isotope(142.90981, 0.122, 143),
                  144: Isotope(143.910083, 0.238, 144),
                  145: Isotope(144.912569, 0.083, 145),
                  146: Isotope(145.913112, 0.172, 146),
                  148: Isotope(147.916889, 0.057, 148),
                  150: Isotope(149.920887, 0.056, 150)}),
    Element(
        61, 'Pm', 'Promethium',
        group=3, period=6, block='f', series=9,
        mass=144.912744, eleneg=1.13, eleaffin=0.5,
        covrad=1.63, atmrad=2.62, vdwrad=0.0,
        tboil=3273.0, tmelt=1315.0, density=7.22,
        eleconfig='[Xe] 4f5 6s2',
        oxistates='3*',
        ionenergy=(5.582, 10.9, ),
        isotopes={145: Isotope(144.912744, 1.0, 145)}),
    Element(
        62, 'Sm', 'Samarium',
        group=3, period=6, block='f', series=9,
        mass=150.36, eleneg=1.17, eleaffin=0.5,
        covrad=1.62, atmrad=2.59, vdwrad=0.0,
        tboil=2067.0, tmelt=1347.0, density=7.54,
        eleconfig='[Xe] 4f6 6s2',
        oxistates='3*, 2',
        ionenergy=(5.6437, 11.07, ),
        isotopes={144: Isotope(143.911995, 0.0307, 144),
                  147: Isotope(146.914893, 0.1499, 147),
                  148: Isotope(147.914818, 0.1124, 148),
                  149: Isotope(148.91718, 0.1382, 149),
                  150: Isotope(149.917271, 0.0738, 150),
                  152: Isotope(151.919728, 0.2675, 152),
                  154: Isotope(153.922205, 0.2275, 154)}),
    Element(
        63, 'Eu', 'Europium',
        group=3, period=6, block='f', series=9,
        mass=151.964, eleneg=1.2, eleaffin=0.5,
        covrad=1.85, atmrad=2.56, vdwrad=0.0,
        tboil=1800.0, tmelt=1095.0, density=5.25,
        eleconfig='[Xe] 4f7 6s2',
        oxistates='3*, 2',
        ionenergy=(5.6704, 11.25, ),
        isotopes={151: Isotope(150.919846, 0.4781, 151),
                  153: Isotope(152.921226, 0.5219, 153)}),
    Element(
        64, 'Gd', 'Gadolinium',
        group=3, period=6, block='f', series=9,
        mass=157.25, eleneg=1.2, eleaffin=0.5,
        covrad=1.61, atmrad=2.54, vdwrad=0.0,
        tboil=3545.0, tmelt=1585.0, density=7.89,
        eleconfig='[Xe] 4f7 5d 6s2',
        oxistates='3*',
        ionenergy=(6.1498, 12.1, ),
        isotopes={152: Isotope(151.919788, 0.002, 152),
                  154: Isotope(153.920862, 0.0218, 154),
                  155: Isotope(154.922619, 0.148, 155),
                  156: Isotope(155.92212, 0.2047, 156),
                  157: Isotope(156.923957, 0.1565, 157),
                  158: Isotope(157.924101, 0.2484, 158),
                  160: Isotope(159.927051, 0.2186, 160)}),
    Element(
        65, 'Tb', 'Terbium',
        group=3, period=6, block='f', series=9,
        mass=158.92534, eleneg=1.2, eleaffin=0.5,
        covrad=1.59, atmrad=2.51, vdwrad=0.0,
        tboil=3500.0, tmelt=1629.0, density=8.25,
        eleconfig='[Xe] 4f9 6s2',
        oxistates='4, 3*',
        ionenergy=(5.8638, 11.52, ),
        isotopes={159: Isotope(158.925343, 1.0, 159)}),
    Element(
        66, 'Dy', 'Dysprosium',
        group=3, period=6, block='f', series=9,
        mass=162.5, eleneg=1.22, eleaffin=0.5,
        covrad=1.59, atmrad=2.49, vdwrad=0.0,
        tboil=2840.0, tmelt=1685.0, density=8.56,
        eleconfig='[Xe] 4f10 6s2',
        oxistates='3*',
        ionenergy=(5.9389, 11.67, ),
        isotopes={156: Isotope(155.924278, 0.0006, 156),
                  158: Isotope(157.924405, 0.001, 158),
                  160: Isotope(159.925194, 0.0234, 160),
                  161: Isotope(160.92693, 0.1891, 161),
                  162: Isotope(161.926795, 0.2551, 162),
                  163: Isotope(162.928728, 0.249, 163),
                  164: Isotope(163.929171, 0.2818, 164)}),
    Element(
        67, 'Ho', 'Holmium',
        group=3, period=6, block='f', series=9,
        mass=164.93032, eleneg=1.23, eleaffin=0.5,
        covrad=1.58, atmrad=2.47, vdwrad=0.0,
        tboil=2968.0, tmelt=1747.0, density=8.78,
        eleconfig='[Xe] 4f11 6s2',
        oxistates='3*',
        ionenergy=(6.0215, 11.8, ),
        isotopes={165: Isotope(164.930319, 1.0, 165)}),
    Element(
        68, 'Er', 'Erbium',
        group=3, period=6, block='f', series=9,
        mass=167.259, eleneg=1.24, eleaffin=0.5,
        covrad=1.57, atmrad=2.45, vdwrad=0.0,
        tboil=3140.0, tmelt=1802.0, density=9.05,
        eleconfig='[Xe] 4f12 6s2',
        oxistates='3*',
        ionenergy=(6.1077, 11.93, ),
        isotopes={162: Isotope(161.928775, 0.0014, 162),
                  164: Isotope(163.929197, 0.0161, 164),
                  166: Isotope(165.93029, 0.3361, 166),
                  167: Isotope(166.932045, 0.2293, 167),
                  168: Isotope(167.932368, 0.2678, 168),
                  170: Isotope(169.93546, 0.1493, 170)}),
    Element(
        69, 'Tm', 'Thulium',
        group=3, period=6, block='f', series=9,
        mass=168.93421, eleneg=1.25, eleaffin=0.5,
        covrad=1.56, atmrad=2.42, vdwrad=0.0,
        tboil=2223.0, tmelt=1818.0, density=9.32,
        eleconfig='[Xe] 4f13 6s2',
        oxistates='3*, 2',
        ionenergy=(6.1843, 12.05, 23.71, ),
        isotopes={169: Isotope(168.934211, 1.0, 169)}),
    Element(
        70, 'Yb', 'Ytterbium',
        group=3, period=6, block='f', series=9,
        mass=173.04, eleneg=1.1, eleaffin=0.5,
        covrad=1.74, atmrad=2.4, vdwrad=0.0,
        tboil=1469.0, tmelt=1092.0, density=9.32,
        eleconfig='[Xe] 4f14 6s2',
        oxistates='3*, 2',
        ionenergy=(6.2542, 12.17, 25.2, ),
        isotopes={168: Isotope(167.933894, 0.0013, 168),
                  170: Isotope(169.934759, 0.0304, 170),
                  171: Isotope(170.936322, 0.1428, 171),
                  172: Isotope(171.9363777, 0.2183, 172),
                  173: Isotope(172.9382068, 0.1613, 173),
                  174: Isotope(173.9388581, 0.3183, 174),
                  176: Isotope(175.942568, 0.1276, 176)}),
    Element(
        71, 'Lu', 'Lutetium',
        group=3, period=6, block='d', series=9,
        mass=174.967, eleneg=1.27, eleaffin=0.5,
        covrad=1.56, atmrad=2.25, vdwrad=0.0,
        tboil=3668.0, tmelt=1936.0, density=9.84,
        eleconfig='[Xe] 4f14 5d 6s2',
        oxistates='3*',
        ionenergy=(5.4259, 13.9, ),
        isotopes={175: Isotope(174.9407679, 0.9741, 175),
                  176: Isotope(175.9426824, 0.0259, 176)}),
    Element(
        72, 'Hf', 'Hafnium',
        group=4, period=6, block='d', series=8,
        mass=178.49, eleneg=1.3, eleaffin=0.0,
        covrad=1.44, atmrad=2.16, vdwrad=0.0,
        tboil=4875.0, tmelt=2504.0, density=13.31,
        eleconfig='[Xe] 4f14 5d2 6s2',
        oxistates='4*',
        ionenergy=(6.8251, 14.9, 23.3, 33.3, ),
        isotopes={174: Isotope(173.94004, 0.0016, 174),
                  176: Isotope(175.9414018, 0.0526, 176),
                  177: Isotope(176.94322, 0.186, 177),
                  178: Isotope(177.9436977, 0.2728, 178),
                  179: Isotope(178.9458151, 0.1362, 179),
                  180: Isotope(179.9465488, 0.3508, 180)}),
    Element(
        73, 'Ta', 'Tantalum',
        group=5, period=6, block='d', series=8,
        mass=180.9479, eleneg=1.5, eleaffin=0.322,
        covrad=1.34, atmrad=2.09, vdwrad=0.0,
        tboil=5730.0, tmelt=3293.0, density=16.68,
        eleconfig='[Xe] 4f14 5d3 6s2',
        oxistates='5*',
        ionenergy=(7.5496, ),
        isotopes={180: Isotope(179.947466, 0.00012, 180),
                  181: Isotope(180.947996, 0.99988, 181)}),
    Element(
        74, 'W', 'Tungsten',
        group=6, period=6, block='d', series=8,
        mass=183.84, eleneg=2.36, eleaffin=0.815,
        covrad=1.3, atmrad=2.02, vdwrad=0.0,
        tboil=5825.0, tmelt=3695.0, density=19.26,
        eleconfig='[Xe] 4f14 5d4 6s2',
        oxistates='6*, 5, 4, 3, 2, 0',
        ionenergy=(7.864, ),
        isotopes={180: Isotope(179.946706, 0.0012, 180),
                  182: Isotope(181.948206, 0.265, 182),
                  183: Isotope(182.9502245, 0.1431, 183),
                  184: Isotope(183.9509326, 0.3064, 184),
                  186: Isotope(185.954362, 0.2843, 186)}),
    Element(
        75, 'Re', 'Rhenium',
        group=7, period=6, block='d', series=8,
        mass=186.207, eleneg=1.9, eleaffin=0.15,
        covrad=1.28, atmrad=1.97, vdwrad=0.0,
        tboil=5870.0, tmelt=3455.0, density=21.03,
        eleconfig='[Xe] 4f14 5d5 6s2',
        oxistates='7, 6, 4, 2, -1',
        ionenergy=(7.8335, ),
        isotopes={185: Isotope(184.9529557, 0.374, 185),
                  187: Isotope(186.9557508, 0.626, 187)}),
    Element(
        76, 'Os', 'Osmium',
        group=8, period=6, block='d', series=8,
        mass=190.23, eleneg=2.2, eleaffin=1.0778,
        covrad=1.26, atmrad=1.92, vdwrad=0.0,
        tboil=5300.0, tmelt=3300.0, density=22.61,
        eleconfig='[Xe] 4f14 5d6 6s2',
        oxistates='8, 6, 4*, 3, 2, 0, -2',
        ionenergy=(8.4382, ),
        isotopes={184: Isotope(183.952491, 0.0002, 184),
                  186: Isotope(185.953838, 0.0159, 186),
                  187: Isotope(186.9557479, 0.0196, 187),
                  188: Isotope(187.955836, 0.1324, 188),
                  189: Isotope(188.9581449, 0.1615, 189),
                  190: Isotope(189.958445, 0.2626, 190),
                  192: Isotope(191.961479, 0.4078, 192)}),
    Element(
        77, 'Ir', 'Iridium',
        group=9, period=6, block='d', series=8,
        mass=192.217, eleneg=2.2, eleaffin=1.56436,
        covrad=1.27, atmrad=1.87, vdwrad=0.0,
        tboil=4700.0, tmelt=2720.0, density=22.65,
        eleconfig='[Xe] 4f14 5d7 6s2',
        oxistates='6, 4*, 3, 2, 1*, 0, -1',
        ionenergy=(8.967, ),
        isotopes={191: Isotope(190.960591, 0.373, 191),
                  193: Isotope(192.962924, 0.627, 193)}),
    Element(
        78, 'Pt', 'Platinum',
        group=10, period=6, block='d', series=8,
        mass=195.078, eleneg=2.28, eleaffin=2.1251,
        covrad=1.3, atmrad=1.83, vdwrad=1.75,
        tboil=4100.0, tmelt=2042.1, density=21.45,
        eleconfig='[Xe] 4f14 5d9 6s',
        oxistates='4*, 2*, 0',
        ionenergy=(8.9588, 18.563, ),
        isotopes={190: Isotope(189.95993, 0.00014, 190),
                  192: Isotope(191.961035, 0.00782, 192),
                  194: Isotope(193.962664, 0.32967, 194),
                  195: Isotope(194.964774, 0.33832, 195),
                  196: Isotope(195.964935, 0.25242, 196),
                  198: Isotope(197.967876, 0.07163, 198)}),
    Element(
        79, 'Au', 'Gold',
        group=11, period=6, block='d', series=8,
        mass=196.96655, eleneg=2.54, eleaffin=2.30861,
        covrad=1.34, atmrad=1.79, vdwrad=1.66,
        tboil=3130.0, tmelt=1337.58, density=19.32,
        eleconfig='[Xe] 4f14 5d10 6s',
        oxistates='3*, 1',
        ionenergy=(9.2255, 20.5, ),
        isotopes={197: Isotope(196.966552, 1.0, 197)}),
    Element(
        80, 'Hg', 'Mercury',
        group=12, period=6, block='d', series=8,
        mass=200.59, eleneg=2.0, eleaffin=0.0,
        covrad=1.49, atmrad=1.76, vdwrad=0.0,
        tboil=629.88, tmelt=234.31, density=13.55,
        eleconfig='[Xe] 4f14 5d10 6s2',
        oxistates='2*, 1',
        ionenergy=(10.4375, 18.756, 34.2, ),
        isotopes={196: Isotope(195.965815, 0.0015, 196),
                  198: Isotope(197.966752, 0.0997, 198),
                  199: Isotope(198.968262, 0.1687, 199),
                  200: Isotope(199.968309, 0.231, 200),
                  201: Isotope(200.970285, 0.1318, 201),
                  202: Isotope(201.970626, 0.2986, 202),
                  204: Isotope(203.973476, 0.0687, 204)}),
    Element(
        81, 'Tl', 'Thallium',
        group=13, period=6, block='p', series=7,
        mass=204.3833, eleneg=2.04, eleaffin=0.377,
        covrad=1.48, atmrad=2.08, vdwrad=1.96,
        tboil=1746.0, tmelt=577.0, density=11.85,
        eleconfig='[Xe] 4f14 5d10 6s2 6p',
        oxistates='3, 1*',
        ionenergy=(6.1082, 20.428, 29.83, ),
        isotopes={203: Isotope(202.972329, 0.29524, 203),
                  205: Isotope(204.974412, 0.70476, 205)}),
    Element(
        82, 'Pb', 'Lead',
        group=14, period=6, block='p', series=7,
        mass=207.2, eleneg=2.33, eleaffin=0.364,
        covrad=1.47, atmrad=1.81, vdwrad=2.02,
        tboil=2023.0, tmelt=600.65, density=11.34,
        eleconfig='[Xe] 4f14 5d10 6s2 6p2',
        oxistates='4, 2*',
        ionenergy=(7.4167, 15.032, 31.937, 42.32, 68.8, ),
        isotopes={204: Isotope(203.973029, 0.014, 204),
                  206: Isotope(205.974449, 0.241, 206),
                  207: Isotope(206.975881, 0.221, 207),
                  208: Isotope(207.976636, 0.524, 208)}),
    Element(
        83, 'Bi', 'Bismuth',
        group=15, period=6, block='p', series=7,
        mass=208.98038, eleneg=2.02, eleaffin=0.942363,
        covrad=1.46, atmrad=1.63, vdwrad=0.0,
        tboil=1837.0, tmelt=544.59, density=9.8,
        eleconfig='[Xe] 4f14 5d10 6s2 6p3',
        oxistates='5, 3*',
        ionenergy=(7.2855, 16.69, 25.56, 45.3, 56.0,
                   88.3, ),
        isotopes={209: Isotope(208.980383, 1.0, 209)}),
    Element(
        84, 'Po', 'Polonium',
        group=16, period=6, block='p', series=5,
        mass=208.982416, eleneg=2.0, eleaffin=1.9,
        covrad=1.46, atmrad=1.53, vdwrad=0.0,
        tboil=0.0, tmelt=527.0, density=9.2,
        eleconfig='[Xe] 4f14 5d10 6s2 6p4',
        oxistates='6, 4*, 2',
        ionenergy=(8.414, ),
        isotopes={209: Isotope(208.982416, 1.0, 209)}),
    Element(
        85, 'At', 'Astatine',
        group=17, period=6, block='p', series=6,
        mass=209.9871, eleneg=2.2, eleaffin=2.8,
        covrad=1.45, atmrad=1.43, vdwrad=0.0,
        tboil=610.0, tmelt=575.0, density=0.0,
        eleconfig='[Xe] 4f14 5d10 6s2 6p5',
        oxistates='7, 5, 3, 1, -1*',
        ionenergy=(),
        isotopes={210: Isotope(209.987131, 1.0, 210)}),
    Element(
        86, 'Rn', 'Radon',
        group=18, period=6, block='p', series=2,
        mass=222.0176, eleneg=0.0, eleaffin=0.0,
        covrad=0.0, atmrad=1.34, vdwrad=0.0,
        tboil=211.4, tmelt=202.0, density=9.23,
        eleconfig='[Xe] 4f14 5d10 6s2 6p6',
        oxistates='2*',
        ionenergy=(10.7485, ),
        isotopes={222: Isotope(222.0175705, 1.0, 222)}),
    Element(
        87, 'Fr', 'Francium',
        group=1, period=7, block='s', series=3,
        mass=223.0197307, eleneg=0.7, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=950.0, tmelt=300.0, density=0.0,
        eleconfig='[Rn] 7s',
        oxistates='1*',
        ionenergy=(4.0727, ),
        isotopes={223: Isotope(223.0197307, 1.0, 223)}),
    Element(
        88, 'Ra', 'Radium',
        group=2, period=7, block='s', series=4,
        mass=226.025403, eleneg=0.9, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=1413.0, tmelt=973.0, density=5.5,
        eleconfig='[Rn] 7s2',
        oxistates='2*',
        ionenergy=(5.2784, 10.147, ),
        isotopes={226: Isotope(226.0254026, 1.0, 226)}),
    Element(
        89, 'Ac', 'Actinium',
        group=3, period=7, block='f', series=10,
        mass=227.027747, eleneg=1.1, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=3470.0, tmelt=1324.0, density=10.07,
        eleconfig='[Rn] 6d 7s2',
        oxistates='3*',
        ionenergy=(5.17, 12.1, ),
        isotopes={227: Isotope(227.027747, 1.0, 227)}),
    Element(
        90, 'Th', 'Thorium',
        group=3, period=7, block='f', series=10,
        mass=232.0381, eleneg=1.3, eleaffin=0.0,
        covrad=1.65, atmrad=0.0, vdwrad=0.0,
        tboil=5060.0, tmelt=2028.0, density=11.72,
        eleconfig='[Rn] 6d2 7s2',
        oxistates='4*',
        ionenergy=(6.3067, 11.5, 20.0, 28.8, ),
        isotopes={232: Isotope(232.0380504, 1.0, 232)}),
    Element(
        91, 'Pa', 'Protactinium',
        group=3, period=7, block='f', series=10,
        mass=231.03588, eleneg=1.5, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=4300.0, tmelt=1845.0, density=15.37,
        eleconfig='[Rn] 5f2 6d 7s2',
        oxistates='5*, 4',
        ionenergy=(5.89, ),
        isotopes={231: Isotope(231.0358789, 1.0, 231)}),
    Element(
        92, 'U', 'Uranium',
        group=3, period=7, block='f', series=10,
        mass=238.02891, eleneg=1.38, eleaffin=0.0,
        covrad=1.42, atmrad=0.0, vdwrad=1.86,
        tboil=4407.0, tmelt=1408.0, density=18.97,
        eleconfig='[Rn] 5f3 6d 7s2',
        oxistates='6*, 5, 4, 3',
        ionenergy=(6.1941, ),
        isotopes={234: Isotope(234.0409456, 5.5e-05, 234),
                  235: Isotope(235.0439231, 0.0072, 235),
                  238: Isotope(238.0507826, 0.992745, 238)}),
    Element(
        93, 'Np', 'Neptunium',
        group=3, period=7, block='f', series=10,
        mass=237.048167, eleneg=1.36, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=4175.0, tmelt=912.0, density=20.48,
        eleconfig='[Rn] 5f4 6d 7s2',
        oxistates='6, 5*, 4, 3',
        ionenergy=(6.2657, ),
        isotopes={237: Isotope(237.0481673, 1.0, 237)}),
    Element(
        94, 'Pu', 'Plutonium',
        group=3, period=7, block='f', series=10,
        mass=244.064198, eleneg=1.28, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=3505.0, tmelt=913.0, density=19.74,
        eleconfig='[Rn] 5f6 7s2',
        oxistates='6, 5, 4*, 3',
        ionenergy=(6.026, ),
        isotopes={244: Isotope(244.064198, 1.0, 244)}),
    Element(
        95, 'Am', 'Americium',
        group=3, period=7, block='f', series=10,
        mass=243.061373, eleneg=1.3, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=2880.0, tmelt=1449.0, density=13.67,
        eleconfig='[Rn] 5f7 7s2',
        oxistates='6, 5, 4, 3*',
        ionenergy=(5.9738, ),
        isotopes={243: Isotope(243.0613727, 1.0, 243)}),
    Element(
        96, 'Cm', 'Curium',
        group=3, period=7, block='f', series=10,
        mass=247.070347, eleneg=1.3, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=0.0, tmelt=1620.0, density=13.51,
        eleconfig='[Rn] 5f7 6d 7s2',
        oxistates='4, 3*',
        ionenergy=(5.9914, ),
        isotopes={247: Isotope(247.070347, 1.0, 247)}),
    Element(
        97, 'Bk', 'Berkelium',
        group=3, period=7, block='f', series=10,
        mass=247.070299, eleneg=1.3, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=0.0, tmelt=1258.0, density=13.25,
        eleconfig='[Rn] 5f9 7s2',
        oxistates='4, 3*',
        ionenergy=(6.1979, ),
        isotopes={247: Isotope(247.070299, 1.0, 247)}),
    Element(
        98, 'Cf', 'Californium',
        group=3, period=7, block='f', series=10,
        mass=251.07958, eleneg=1.3, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=0.0, tmelt=1172.0, density=15.1,
        eleconfig='[Rn] 5f10 7s2',
        oxistates='4, 3*',
        ionenergy=(6.2817, ),
        isotopes={251: Isotope(251.07958, 1.0, 251)}),
    Element(
        99, 'Es', 'Einsteinium',
        group=3, period=7, block='f', series=10,
        mass=252.08297, eleneg=1.3, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=0.0, tmelt=1130.0, density=0.0,
        eleconfig='[Rn] 5f11 7s2',
        oxistates='3*',
        ionenergy=(6.42, ),
        isotopes={252: Isotope(252.08297, 1.0, 252)}),
    Element(
        100, 'Fm', 'Fermium',
        group=3, period=7, block='f', series=10,
        mass=257.095099, eleneg=1.3, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=0.0, tmelt=1800.0, density=0.0,
        eleconfig='[Rn] 5f12 7s2',
        oxistates='3*',
        ionenergy=(6.5, ),
        isotopes={257: Isotope(257.095099, 1.0, 257)}),
    Element(
        101, 'Md', 'Mendelevium',
        group=3, period=7, block='f', series=10,
        mass=258.098425, eleneg=1.3, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=0.0, tmelt=1100.0, density=0.0,
        eleconfig='[Rn] 5f13 7s2',
        oxistates='3*',
        ionenergy=(6.58, ),
        isotopes={258: Isotope(258.098425, 1.0, 258)}),
    Element(
        102, 'No', 'Nobelium',
        group=3, period=7, block='f', series=10,
        mass=259.10102, eleneg=1.3, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=0.0, tmelt=1100.0, density=0.0,
        eleconfig='[Rn] 5f14 7s2',
        oxistates='3, 2*',
        ionenergy=(6.65, ),
        isotopes={259: Isotope(259.10102, 1.0, 259)}),
    Element(
        103, 'Lr', 'Lawrencium',
        group=3, period=7, block='d', series=10,
        mass=262.10969, eleneg=1.3, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=0.0, tmelt=1900.0, density=0.0,
        eleconfig='[Rn] 5f14 6d 7s2',
        oxistates='3*',
        ionenergy=(4.9, ),
        isotopes={262: Isotope(262.10969, 1.0, 262)}),
    Element(
        104, 'Rf', 'Rutherfordium',
        group=4, period=7, block='d', series=8,
        mass=261.10875, eleneg=0.0, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=0.0, tmelt=0.0, density=0.0,
        eleconfig='[Rn] 5f14 6d2 7s2',
        oxistates='*',
        ionenergy=(6.0, ),
        isotopes={261: Isotope(261.10875, 1.0, 261)}),
    Element(
        105, 'Db', 'Dubnium',
        group=5, period=7, block='d', series=8,
        mass=262.11415, eleneg=0.0, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=0.0, tmelt=0.0, density=0.0,
        eleconfig='[Rn] 5f14 6d3 7s2',
        oxistates='*',
        ionenergy=(),
        isotopes={262: Isotope(262.11415, 1.0, 262)}),
    Element(
        106, 'Sg', 'Seaborgium',
        group=6, period=7, block='d', series=8,
        mass=266.12193, eleneg=0.0, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=0.0, tmelt=0.0, density=0.0,
        eleconfig='[Rn] 5f14 6d4 7s2',
        oxistates='*',
        ionenergy=(),
        isotopes={266: Isotope(266.12193, 1.0, 266)}),
    Element(
        107, 'Bh', 'Bohrium',
        group=7, period=7, block='d', series=8,
        mass=264.12473, eleneg=0.0, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=0.0, tmelt=0.0, density=0.0,
        eleconfig='[Rn] 5f14 6d5 7s2',
        oxistates='*',
        ionenergy=(),
        isotopes={264: Isotope(264.12473, 1.0, 264)}),
    Element(
        108, 'Hs', 'Hassium',
        group=8, period=7, block='d', series=8,
        mass=269.13411, eleneg=0.0, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=0.0, tmelt=0.0, density=0.0,
        eleconfig='[Rn] 5f14 6d6 7s2',
        oxistates='*',
        ionenergy=(),
        isotopes={269: Isotope(269.13411, 1.0, 269)}),
    Element(
        109, 'Mt', 'Meitnerium',
        group=9, period=7, block='d', series=8,
        mass=268.13882, eleneg=0.0, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=0.0, tmelt=0.0, density=0.0,
        eleconfig='[Rn] 5f14 6d7 7s2',
        oxistates='*',
        ionenergy=(),
        isotopes={268: Isotope(268.13882, 1.0, 268)}))


PERIODS = {1: 'K', 2: 'L', 3: 'M', 4: 'N', 5: 'O', 6: 'P', 7: 'Q'}

BLOCKS = {'s': '', 'g': '', 'f': '', 'd': '', 'p': ''}

GROUPS = {
    1: ('IA', 'Alkali metals'),
    2: ('IIA', 'Alkaline earths'),
    3: ('IIIB', ''),
    4: ('IVB', ''),
    5: ('VB', ''),
    6: ('VIB', ''),
    7: ('VIIB', ''),
    8: ('VIIIB', ''),
    9: ('VIIIB', ''),
    10: ('VIIIB', ''),
    11: ('IB', 'Coinage metals'),
    12: ('IIB', ''),
    13: ('IIIA', 'Boron group'),
    14: ('IVA', 'Carbon group'),
    15: ('VA', 'Pnictogens'),
    16: ('VIA', 'Chalcogens'),
    17: ('VIIA', 'Halogens'),
    18: ('VIIIA', 'Noble gases')}

SERIES = {
    1: 'Nonmetals',
    2: 'Noble gases',
    3: 'Alkali metals',
    4: 'Alkaline earth metals',
    5: 'Metalloids',
    6: 'Halogens',
    7: 'Poor metals',
    8: 'Transition metals',
    9: 'Lanthanides',
    10: 'Actinides'}


def _descriptions(symbol):
    """Delay load descriptions."""
    e = ELEMENTS
    e['H'].description = (
        "Colourless, odourless gaseous chemical element. Lightest and "
        "most abundant element in the universe. Present in water and in "
        "all organic compounds. Chemically reacts with most elements. "
        "Discovered by Henry Cavendish in 1776.")
    e['He'].description = (
        "Colourless, odourless gaseous nonmetallic element. Belongs to "
        "group 18 of the periodic table. Lowest boiling point of all "
        "elements and can only be solidified under pressure. Chemically "
        "inert, no known compounds. Discovered in the solar spectrum in "
        "1868 by Lockyer.")
    e['Li'].description = (
        "Socket silvery metal. First member of group 1 of the periodic "
        "table. Lithium salts are used in psychomedicine.")
    e['Be'].description = (
        "Grey metallic element of group 2 of the periodic table. Is toxic "
        "and can cause severe lung diseases and dermatitis. Shows high "
        "covalent character. It was isolated independently by F. Wohler "
        "and A.A. Bussy in 1828.")
    e['B'].description = (
        "An element of group 13 of the periodic table. There are two "
        "allotropes, amorphous boron is a brown power, but metallic boron "
        "is black. The metallic form is hard (9.3 on Mohs' scale) and a "
        "bad conductor in room temperatures. It is never found free in "
        "nature. Boron-10 is used in nuclear reactor control rods and "
        "shields. It was discovered in 1808 by Sir Humphry Davy and by "
        "J.L. Gay-Lussac and L.J. Thenard.")
    e['C'].description = (
        "Carbon is a member of group 14 of the periodic table. It has "
        "three allotropic forms of it, diamonds, graphite and fullerite. "
        "Carbon-14 is commonly used in radioactive dating. Carbon occurs "
        "in all organic life and is the basis of organic chemistry. Carbon "
        "has the interesting chemical property of being able to bond with "
        "itself, and a wide variety of other elements.")
    e['N'].description = (
        "Colourless, gaseous element which belongs to group 15 of the "
        "periodic table. Constitutes ~78% of the atmosphere and is an "
        "essential part of the ecosystem. Nitrogen for industrial purposes "
        "is acquired by the fractional distillation of liquid air. "
        "Chemically inactive, reactive generally only at high temperatures "
        "or in electrical discharges. It was discovered in 1772 by D. "
        "Rutherford.")
    e['O'].description = (
        "A colourless, odourless gaseous element belonging to group 16 of "
        "the periodic table. It is the most abundant element present in "
        "the earth's crust. It also makes up 20.8% of the Earth's "
        "atmosphere. For industrial purposes, it is separated from liquid "
        "air by fractional distillation. It is used in high temperature "
        "welding, and in breathing. It commonly comes in the form of "
        "Oxygen, but is found as Ozone in the upper atmosphere. It was "
        "discovered by Priestley in 1774.")
    e['F'].description = (
        "A poisonous pale yellow gaseous element belonging to group 17 of "
        "the periodic table (The halogens). It is the most chemically "
        "reactive and electronegative element. It is highly dangerous, "
        "causing severe chemical burns on contact with flesh. Fluorine was "
        "identified by Scheele in 1771 and first isolated by Moissan in "
        "1886.")
    e['Ne'].description = (
        "Colourless gaseous element of group 18 on the periodic table "
        "(noble gases). Neon occurs in the atmosphere, and comprises "
        "0.0018% of the volume of the atmosphere. It has a distinct "
        "reddish glow when used in discharge tubes and neon based lamps. "
        "It forms almost no chemical compounds. Neon was discovered in "
        "1898 by Sir William Ramsey and M.W. Travers.")
    e['Na'].description = (
        "Soft silvery reactive element belonging to group 1 of the "
        "periodic table (alkali metals). It is highly reactive, oxidizing "
        "in air and reacting violently with water, forcing it to be kept "
        "under oil. It was first isolated by Humphrey Davy in 1807.")
    e['Mg'].description = (
        "Silvery metallic element belonging to group 2 of the periodic "
        "table (alkaline-earth metals). It is essential for living "
        "organisms, and is used in a number of light alloys. Chemically "
        "very reactive, it forms a protective oxide coating when exposed "
        "to air and burns with an intense white flame. It also reacts with "
        "sulphur, nitrogen and the halogens. First isolated by Bussy in "
        "1828.")
    e['Al'].description = (
        "Silvery-white lustrous metallic element of group 3 of the "
        "periodic table. Highly reactive but protected by a thin "
        "transparent layer of the oxide which quickly forms in air. There "
        "are many alloys of aluminum, as well as a good number of "
        "industrial uses. Makes up 8.1% of the Earth's crust, by weight. "
        "Isolated in 1825 by H.C. Oersted.")
    e['Si'].description = (
        "Metalloid element belonging to group 14 of the periodic table. "
        "It is the second most abundant element in the Earth's crust, "
        "making up 25.7% of it by weight. Chemically less reactive than "
        "carbon. First identified by Lavoisier in 1787 and first isolated "
        "in 1823 by Berzelius.")
    e['P'].description = (
        "Non-metallic element belonging to group 15 of the periodic "
        "table. Has a multiple allotropic forms. Essential element for "
        "living organisms. It was discovered by Brandt in 1669.")
    e['S'].description = (
        "Yellow, nonmetallic element belonging to group 16 of the "
        "periodic table. It is an essential element in living organisms, "
        "needed in the amino acids cysteine and methionine, and hence in "
        "many proteins. Absorbed by plants from the soil as sulphate ion.")
    e['Cl'].description = (
        "Halogen element. Poisonous greenish-yellow gas. Occurs widely in "
        "nature as sodium chloride in seawater. Reacts directly with many "
        "elements and compounds, strong oxidizing agent. Discovered by "
        "Karl Scheele in 1774. Humphrey David confirmed it as an element "
        "in 1810.")
    e['Ar'].description = (
        "Monatomic noble gas. Makes up 0.93% of the air. Colourless, "
        "odorless. Is inert and has no true compounds. Lord Rayleigh and "
        "Sir william Ramsey identified argon in 1894.")
    e['K'].description = (
        "Soft silvery metallic element belonging to group 1 of the "
        "periodic table (alkali metals). Occurs naturally in seawater and "
        "a many minerals. Highly reactive, chemically, it resembles sodium "
        "in its behavior and compounds. Discovered by Sir Humphry Davy in "
        "1807.")
    e['Ca'].description = (
        "Soft grey metallic element belonging to group 2 of the periodic "
        "table. Used a reducing agent in the extraction of thorium, "
        "zirconium and uranium. Essential element for living organisms.")
    e['Sc'].description = (
        "Rare soft silvery metallic element belonging to group 3 of the "
        "periodic table. There are ten isotopes, nine of which are "
        "radioactive and have short half-lives. Predicted in 1869 by "
        "Mendeleev, isolated by Nilson in 1879.")
    e['Ti'].description = (
        "White metallic transition element. Occurs in numerous minerals. "
        "Used in strong, light corrosion-resistant alloys. Forms a passive "
        "oxide coating when exposed to air. First discovered by Gregor in "
        "1789.")
    e['V'].description = (
        "Soft and ductile, bright white metal. Good resistance to "
        "corrosion by alkalis, sulphuric and hydrochloric acid. It "
        "oxidizes readily about 933K. There are two naturally occurring "
        "isotopes of vanadium, and 5 radioisotopes, V-49 having the "
        "longest half-life at 337 days. Vanadium has nuclear applications, "
        "the foil is used in cladding titanium to steel, and "
        "vanadium-gallium tape is used to produce a superconductive "
        "magnet. Originally discovered by Andres Manuel del Rio of Mexico "
        "City in 1801. His discovery went unheeded, however, and in 1820, "
        "Nils Gabriel Sefstron of Sweden rediscovered it. Metallic "
        "vanadium was isolated by Henry Enfield Roscoe in 1867. The name "
        "vanadium comes from Vanadis, a goddess of Scandinavian mythology. "
        "Silvery-white metallic transition element. Vanadium is essential "
        "to Ascidians. Rats and chickens are also known to require it. "
        "Metal powder is a fire hazard, and vanadium compounds should be "
        "considered highly toxic. May cause lung cancer if inhaled.")
    e['Cr'].description = (
        "Hard silvery transition element. Used in decorative "
        "electroplating. Discovered in 1797 by Vauquelin.")
    e['Mn'].description = (
        "Grey brittle metallic transition element. Rather "
        "electropositive, combines with some non-metals when heated. "
        "Discovered in 1774 by Scheele.")
    e['Fe'].description = (
        "Silvery malleable and ductile metallic transition element. Has "
        "nine isotopes and is the fourth most abundant element in the "
        "earth's crust. Required by living organisms as a trace element "
        "(used in hemoglobin in humans.) Quite reactive, oxidizes in moist "
        "air, displaces hydrogen from dilute acids and combines with "
        "nonmetallic elements.")
    e['Co'].description = (
        "Light grey transition element. Some meteorites contain small "
        "amounts of metallic cobalt. Generally alloyed for use. Mammals "
        "require small amounts of cobalt salts. Cobalt-60, an artificially "
        "produced radioactive isotope of Cobalt is an important "
        "radioactive tracer and cancer-treatment agent. Discovered by G. "
        "Brandt in 1737.")
    e['Ni'].description = (
        "Malleable ductile silvery metallic transition element. "
        "Discovered by A.F. Cronstedt in 1751.")
    e['Cu'].description = (
        "Red-brown transition element. Known by the Romans as 'cuprum.' "
        "Extracted and used for thousands of years. Malleable, ductile and "
        "an excellent conductor of heat and electricity. When in moist "
        "conditions, a greenish layer forms on the outside.")
    e['Zn'].description = (
        "Blue-white metallic element. Occurs in multiple compounds "
        "naturally. Five stable isotopes are six radioactive isotopes have "
        "been found. Chemically a reactive metal, combines with oxygen and "
        "other non-metals, reacts with dilute acids to release hydrogen.")
    e['Ga'].description = (
        "Soft silvery metallic element, belongs to group 13 of the "
        "periodic table. The two stable isotopes are Ga-69 and Ga-71. "
        "Eight radioactive isotopes are known, all having short "
        "half-lives. Gallium Arsenide is used as a semiconductor. Corrodes "
        "most other metals by diffusing into their lattice. First "
        "identified by Francois Lecoq de Boisbaudran in 1875.")
    e['Ge'].description = (
        "Lustrous hard metalloid element, belongs to group 14 of the "
        "periodic table. Forms a large number of organometallic compounds. "
        "Predicted by Mendeleev in 1871, it was actually found in 1886 by "
        "Winkler.")
    e['As'].description = (
        "Metalloid element of group 15. There are three allotropes, "
        "yellow, black, and grey. Reacts with halogens, concentrated "
        "oxidizing acids and hot alkalis. Albertus Magnus is believed to "
        "have been the first to isolate the element in 1250.")
    e['Se'].description = (
        "Metalloid element, belongs to group 16 of the periodic table. "
        "Multiple allotropic forms exist. Chemically resembles sulphur. "
        "Discovered in 1817 by Jons J. Berzelius.")
    e['Br'].description = (
        "Halogen element. Red volatile liquid at room temperature. Its "
        "reactivity is somewhere between chlorine and iodine. Harmful to "
        "human tissue in a liquid state, the vapour irritates eyes and "
        "throat. Discovered in 1826 by Antoine Balard.")
    e['Kr'].description = (
        "Colorless gaseous element, belongs to the noble gases. Occurs in "
        "the air, 0.0001% by volume. It can be extracted from liquid air "
        "by fractional distillation. Generally not isolated, but used with "
        "other inert gases in fluorescent lamps. Five natural isotopes, "
        "and five radioactive isotopes. Kr-85, the most stable radioactive "
        "isotope, has a half-life of 10.76 years and is produced in "
        "fission reactors. Practically inert, though known to form "
        "compounds with Fluorine.")
    e['Rb'].description = (
        "Soft silvery metallic element, belongs to group 1 of the "
        "periodic table. Rb-97, the naturally occurring isotope, is "
        "radioactive. It is highly reactive, with properties similar to "
        "other elements in group 1, like igniting spontaneously in air. "
        "Discovered spectroscopically in 1861 by W. Bunsen and G.R. "
        "Kirchoff.")
    e['Sr'].description = (
        "Soft yellowish metallic element, belongs to group 2 of the "
        "periodic table. Highly reactive chemically. Sr-90 is present in "
        "radioactive fallout and has a half-life of 28 years. Discovered "
        "in 1798 by Klaproth and Hope, isolated in 1808 by Humphry Davy.")
    e['Y'].description = (
        "Silvery-grey metallic element of group 3 on the periodic table. "
        "Found in uranium ores. The only natural isotope is Y-89, there "
        "are 14 other artificial isotopes. Chemically resembles the "
        "lanthanoids. Stable in the air below 400 degrees, celsius. "
        "Discovered in 1828 by Friedrich Wohler.")
    e['Zr'].description = (
        "Grey-white metallic transition element. Five natural isotopes "
        "and six radioactive isotopes are known. Used in nuclear reactors "
        "for a Neutron absorber. Discovered in 1789 by Martin Klaproth, "
        "isolated in 1824 by Berzelius.")
    e['Nb'].description = (
        "Soft, ductile grey-blue metallic transition element. Used in "
        "special steels and in welded joints to increase strength. "
        "Combines with halogens and oxidizes in air at 200 degrees "
        "celsius. Discovered by Charles Hatchett in 1801 and isolated by "
        "Blomstrand in 1864. Called Columbium originally.")
    e['Mo'].description = (
        "Silvery-white, hard metallic transition element. It is "
        "chemically unreactive and is not affected by most acids. It "
        "oxidizes at high temperatures. There are seven natural isotopes, "
        "and four radioisotopes, Mo-93 being the most stable with a "
        "half-life of 3500 years. Molybdenum is used in almost all "
        "high-strength steels, it has nuclear applications, and is a "
        "catalyst in petroleum refining. Discovered in 1778 by Carl "
        "Welhelm Scheele of Sweden. Impure metal was prepared in 1782 by "
        "Peter Jacob Hjelm. The name comes from the Greek word molybdos "
        "which means lead. Trace amounts of molybdenum are required for "
        "all known forms of life. All molybdenum compounds should be "
        "considered highly toxic, and will also cause severe birth "
        "defects.")
    e['Tc'].description = (
        "Radioactive metallic transition element. Can be detected in some "
        "stars and the fission products of uranium. First made by Perrier "
        "and Segre by bombarding molybdenum with deutrons, giving them "
        "Tc-97. Tc-99 is the most stable isotope with a half-life of "
        "2.6*10^6 years. Sixteen isotopes are known. Organic technetium "
        "compounds are used in bone imaging. Chemical properties are "
        "intermediate between rhenium and manganese.")
    e['Ru'].description = (
        "Hard white metallic transition element. Found with platinum, "
        "used as a catalyst in some platinum alloys. Dissolves in fused "
        "alkalis, and is not attacked by acids. Reacts with halogens and "
        "oxygen at high temperatures. Isolated in 1844 by K.K. Klaus.")
    e['Rh'].description = (
        "Silvery white metallic transition element. Found with platinum "
        "and used in some platinum alloys. Not attacked by acids, "
        "dissolves only in aqua regia. Discovered in 1803 by W.H. "
        "Wollaston.")
    e['Pd'].description = (
        "Soft white ductile transition element. Found with some copper "
        "and nickel ores. Does not react with oxygen at normal "
        "temperatures. Dissolves slowly in hydrochloric acid. Discovered "
        "in 1803 by W.H. Wollaston.")
    e['Ag'].description = (
        "White lustrous soft metallic transition element. Found in both "
        "its elemental form and in minerals. Used in jewellery, tableware "
        "and so on. Less reactive than silver, chemically.")
    e['Cd'].description = (
        "Soft bluish metal belonging to group 12 of the periodic table. "
        "Extremely toxic even in low concentrations. Chemically similar to "
        "zinc, but lends itself to more complex compounds. Discovered in "
        "1817 by F. Stromeyer.")
    e['In'].description = (
        "Soft silvery element belonging to group 13 of the periodic "
        "table. The most common natural isotope is In-115, which has a "
        "half-life of 6*10^4 years. Five other radioisotopes exist. "
        "Discovered in 1863 by Reich and Richter.")
    e['Sn'].description = (
        "Silvery malleable metallic element belonging to group 14 of the "
        "periodic table. Twenty-six isotopes are known, five of which are "
        "radioactive. Chemically reactive. Combines directly with chlorine "
        "and oxygen and displaces hydrogen from dilute acids.")
    e['Sb'].description = (
        "Element of group 15. Multiple allotropic forms. The stable form "
        "of antimony is a blue-white metal. Yellow and black antimony are "
        "unstable non-metals. Used in flame-proofing, paints, ceramics, "
        "enamels, and rubber. Attacked by oxidizing acids and halogens. "
        "First reported by Tholden in 1450.")
    e['Te'].description = (
        "Silvery metalloid element of group 16. Eight natural isotopes, "
        "nine radioactive isotopes. Used in semiconductors and to a degree "
        "in some steels. Chemistry is similar to Sulphur. Discovered in "
        "1782 by Franz Miller.")
    e['I'].description = (
        "Dark violet nonmetallic element, belongs to group 17 of the "
        "periodic table. Insoluble in water. Required as a trace element "
        "for living organisms. One stable isotope, I-127 exists, in "
        "addition to fourteen radioactive isotopes. Chemically the least "
        "reactive of the halogens, and the most electropositive metallic "
        "halogen. Discovered in 1812 by Courtois.")
    e['Xe'].description = (
        "Colourless, odourless gas belonging to group 18 on the periodic "
        "table (the noble gases.) Nine natural isotopes and seven "
        "radioactive isotopes are known. Xenon was part of the first "
        "noble-gas compound synthesized. Several others involving Xenon "
        "have been found since then. Xenon was discovered by Ramsey and "
        "Travers in 1898.")
    e['Cs'].description = (
        "Soft silvery-white metallic element belonging to group 1 of the "
        "periodic table. One of the three metals which are liquid at room "
        "temperature. Cs-133 is the natural, and only stable, isotope. "
        "Fifteen other radioisotopes exist. Caesium reacts explosively "
        "with cold water, and ice at temperatures above 157K. Caesium "
        "hydroxide is the strongest base known. Caesium is the most "
        "electropositive, most alkaline and has the least ionization "
        "potential of all the elements. Known uses include the basis of "
        "atomic clocks, catalyst for the hydrogenation of some organic "
        "compounds, and in photoelectric cells. Caesium was discovered by "
        "Gustav Kirchoff and Robert Bunsen in Germany in 1860 "
        "spectroscopically. Its identification was based upon the bright "
        "blue lines in its spectrum. The name comes from the latin word "
        "caesius, which means sky blue. Caesium should be considered "
        "highly toxic. Some of the radioisotopes are even more toxic.")
    e['Ba'].description = (
        "Silvery-white reactive element, belonging to group 2 of the "
        "periodic table. Soluble barium compounds are extremely poisonous. "
        "Identified in 1774 by Karl Scheele and extracted in 1808 by "
        "Humphry Davy.")
    e['La'].description = (
        "(From the Greek word lanthanein, to line hidden) Silvery "
        "metallic element belonging to group 3 of the periodic table and "
        "oft considered to be one of the lanthanoids. Found in some "
        "rare-earth minerals. Twenty-five natural isotopes exist. La-139 "
        "which is stable, and La-138 which has a half-life of 10^10 to "
        "10^15 years. The other twenty-three isotopes are radioactive. It "
        "resembles the lanthanoids chemically. Lanthanum has a low to "
        "moderate level of toxicity, and should be handled with care. "
        "Discovered in 1839 by C.G. Mosander.")
    e['Ce'].description = (
        "Silvery metallic element, belongs to the lanthanoids. Four "
        "natural isotopes exist, and fifteen radioactive isotopes have "
        "been identified. Used in some rare-earth alloys. The oxidized "
        "form is used in the glass industry. Discovered by Martin .H. "
        "Klaproth in 1803.")
    e['Pr'].description = (
        "Soft silvery metallic element, belongs to the lanthanoids. Only "
        "natural isotope is Pr-141 which is not radioactive. Fourteen "
        "radioactive isotopes have been artificially produced. Used in "
        "rare-earth alloys. Discovered in 1885 by C.A. von Welsbach.")
    e['Nd'].description = (
        "Soft bright silvery metallic element, belongs to the "
        "lanthanoids. Seven natural isotopes, Nd-144 being the only "
        "radioactive one with a half-life of 10^10 to 10^15 years. Six "
        "artificial radioisotopes have been produced. The metal is used in "
        "glass works to color class a shade of violet-purple and make it "
        "dichroic. One of the more reactive rare-earth metals, quickly "
        "reacts with air. Used in some rare-earth alloys. Neodymium is "
        "used to color the glass used in welder's glasses. Neodymium is "
        "also used in very powerful, permanent magnets (Nd2Fe14B). "
        "Discovered by Carl F. Auer von Welsbach in Austria in 1885 by "
        "separating didymium into its elemental components Praseodymium "
        "and neodymium. The name comes from the Greek words 'neos didymos' "
        "which means 'new twin'. Neodymium should be considered highly "
        "toxic, however evidence would seem to show that it acts as little "
        "more than a skin and eye irritant. The dust however, presents a "
        "fire and explosion hazard.")
    e['Pm'].description = (
        "Soft silvery metallic element, belongs to the lanthanoids. "
        "Pm-147, the only natural isotope, is radioactive and has a "
        "half-life of 252 years. Eighteen radioisotopes have been "
        "produced, but all have very short half-lives. Found only in "
        "nuclear decay waste. Pm-147 is of interest as a beta-decay "
        "source, however Pm-146 and Pm-148 have to be removed from it "
        "first, as they generate gamma radiation. Discovered by J.A. "
        "Marinsky, L.E. Glendenin and C.D. Coryell in 1947.")
    e['Sm'].description = (
        "Soft silvery metallic element, belongs to the lanthanoids. Seven "
        "natural isotopes, Sm-147 is the only radioisotope, and has a "
        "half-life of 2.5*10^11 years. Used for making special alloys "
        "needed in the production of nuclear reactors. Also used as a "
        "neutron absorber. Small quantities of samarium oxide is used in "
        "special optical glasses. The largest use of the element is its "
        "ferromagnetic alloy which produces permanent magnets that are "
        "five times stronger than magnets produced by any other material. "
        "Discovered by Francois Lecoq de Boisbaudran in 1879.")
    e['Eu'].description = (
        "Soft silvery metallic element belonging to the lanthanoids. "
        "Eu-151 and Eu-153 are the only two stable isotopes, both of which "
        "are Neutron absorbers. Discovered in 1889 by Sir William Crookes.")
    e['Gd'].description = (
        "Soft silvery metallic element belonging to the lanthanoids. "
        "Seven natural, stable isotopes are known in addition to eleven "
        "artificial isotopes. Gd-155 and Gd-157 and the best neutron "
        "absorbers of all elements. Gadolinium compounds are used in "
        "electronics. Discovered by J.C.G Marignac in 1880.")
    e['Tb'].description = (
        "Silvery metallic element belonging to the lanthanoids. Tb-159 is "
        "the only stable isotope, there are seventeen artificial isotopes. "
        "Discovered by G.G. Mosander in 1843.")
    e['Dy'].description = (
        "Metallic with a bright silvery-white lustre. Dysprosium belongs "
        "to the lanthanoids. It is relatively stable in air at room "
        "temperatures, it will however dissolve in mineral acids, evolving "
        "hydrogen. It is found in from rare-earth minerals. There are "
        "seven natural isotopes of dysprosium, and eight radioisotopes, "
        "Dy-154 being the most stable with a half-life of 3*10^6 years. "
        "Dysprosium is used as a neutron absorber in nuclear fission "
        "reactions, and in compact disks. It was discovered by Paul Emile "
        "Lecoq de Boisbaudran in 1886 in France. Its name comes from the "
        "Greek word dysprositos, which means hard to obtain.")
    e['Ho'].description = (
        "Relatively soft and malleable silvery-white metallic element, "
        "which is stable in dry air at room temperature. It oxidizes in "
        "moist air and at high temperatures. It belongs to the "
        "lanthanoids. A rare-earth metal, it is found in the minerals "
        "monazite and gadolinite. It possesses unusual magnetic "
        "properties. One natural isotope, Ho-165 exists, six radioisotopes "
        "exist, the most stable being Ho-163 with a half-life of 4570 "
        "years. Holmium is used in some metal alloys, it is also said to "
        "stimulate the metabolism. Discovered by Per Theodor Cleve and "
        "J.L. Soret in Switzerland in 1879. The name homium comes from the "
        "Greek word Holmia which means Sweden. While all holmium compounds "
        "should be considered highly toxic, initial evidence seems to "
        "indicate that they do not pose much danger. The metal's dust "
        "however, is a fire hazard.")
    e['Er'].description = (
        "Soft silvery metallic element which belongs to the lanthanoids. "
        "Six natural isotopes that are stable. Twelve artificial isotopes "
        "are known. Used in nuclear technology as a neutron absorber. It "
        "is being investigated for other possible uses. Discovered by Carl "
        "G. Mosander in 1843.")
    e['Tm'].description = (
        "Soft grey metallic element that belongs to the lanthanoids. One "
        "natural isotope exists, Tm-169, and seventeen artificial isotopes "
        "have been produced. No known uses for the element. Discovered in "
        "1879 by Per Theodor Cleve.")
    e['Yb'].description = (
        "Silvery metallic element of the lanthanoids. Seven natural "
        "isotopes and ten artificial isotopes are known. Used in certain "
        "steels. Discovered by J.D.G. Marignac in 1878.")
    e['Lu'].description = (
        "Silvery-white rare-earth metal which is relatively stable in "
        "air. It happens to be the most expensive rare-earth metal. Its "
        "found with almost all rare-earth metals, but is very difficult to "
        "separate from other elements. Least abundant of all natural "
        "elements. Used in metal alloys, and as a catalyst in various "
        "processes. There are two natural, stable isotopes, and seven "
        "radioisotopes, the most stable being Lu-174 with a half-life of "
        "3.3 years. The separation of lutetium from Ytterbium was "
        "described by Georges Urbain in 1907. It was discovered at "
        "approximately the same time by Carl Auer von Welsbach. The name "
        "comes from the Greek word lutetia which means Paris.")
    e['Hf'].description = (
        "Silvery lustrous metallic transition element. Used in tungsten "
        "alloys in filaments and electrodes, also acts as a neutron "
        "absorber. First reported by Urbain in 1911, existence was finally "
        "established in 1923 by D. Coster, G.C. de Hevesy in 1923.")
    e['Ta'].description = (
        "Heavy blue-grey metallic transition element. Ta-181 is a stable "
        "isotope, and Ta-180 is a radioactive isotope, with a half-life in "
        "excess of 10^7 years. Used in surgery as it is unreactive. Forms "
        "a passive oxide layer in air. Identified in 1802 by Ekeberg and "
        "isolated in 1820 by Jons J. Berzelius.")
    e['W'].description = (
        "White or grey metallic transition element,formerly called "
        "Wolfram. Forms a protective oxide in air and can be oxidized at "
        "high temperature. First isolated by Jose and Fausto de Elhuyer in "
        "1783.")
    e['Re'].description = (
        "Silvery-white metallic transition element. Obtained as a "
        "by-product of molybdenum refinement. Rhenium-molybdenum alloys "
        "are superconducting.")
    e['Os'].description = (
        "Hard blue-white metallic transition element. Found with platinum "
        "and used in some alloys with platinum and iridium.")
    e['Ir'].description = (
        "Very hard and brittle, silvery metallic transition element. It "
        "has a yellowish cast to it. Salts of iridium are highly colored. "
        "It is the most corrosion resistant metal known, not attacked by "
        "any acid, but is attacked by molten salts. There are two natural "
        "isotopes of iridium, and 4 radioisotopes, the most stable being "
        "Ir-192 with a half-life of 73.83 days. Ir-192 decays into "
        "Platinum, while the other radioisotopes decay into Osmium. "
        "Iridium is used in high temperature apparatus, electrical "
        "contacts, and as a hardening agent for platinumpy. Discovered in "
        "1803 by Smithson Tennant in England. The name comes from the "
        "Greek word iris, which means rainbow. Iridium metal is generally "
        "non-toxic due to its relative unreactivity, but iridium compounds "
        "should be considered highly toxic.")
    e['Pt'].description = (
        "Attractive greyish-white metal. When pure, it is malleable and "
        "ductile. Does not oxidize in air, insoluble in hydrochloric and "
        "nitric acid. Corroded by halogens, cyandies, sulphur and alkalis. "
        "Hydrogen and Oxygen react explosively in the presence of "
        "platinumpy. There are six stable isotopes and three "
        "radioisotopes, the most stable being Pt-193 with a half-life of "
        "60 years. Platinum is used in jewelry, laboratory equipment, "
        "electrical contacts, dentistry, and anti-pollution devices in "
        "cars. PtCl2(NH3)2 is used to treat some forms of cancer. "
        "Platinum-Cobalt alloys have magnetic properties. It is also used "
        "in the definition of the Standard Hydrogen Electrode. Discovered "
        "by Antonio de Ulloa in South America in 1735. The name comes from "
        "the Spanish word platina which means silver. Platinum metal is "
        "generally not a health concern due to its unreactivity, however "
        "platinum compounds should be considered highly toxic.")
    e['Au'].description = (
        "Gold is gold colored. It is the most malleable and ductile metal "
        "known. There is only one stable isotope of gold, and five "
        "radioisotopes of gold, Au-195 being the most stable with a "
        "half-life of 186 days. Gold is used as a monetary standard, in "
        "jewelry, dentistry, electronics. Au-198 is used in treating "
        "cancer and some other medical conditions. Gold has been known to "
        "exist as far back as 2600 BC. Gold comes from the Anglo-Saxon "
        "word gold. Its symbol, Au, comes from the Latin word aurum, which "
        "means gold. Gold is not particularly toxic, however it is known "
        "to cause damage to the liver and kidneys in some.")
    e['Hg'].description = (
        "Heavy silvery liquid metallic element, belongs to the zinc "
        "group. Used in thermometers, barometers and other scientific "
        "apparatus. Less reactive than zinc and cadmium, does not displace "
        "hydrogen from acids. Forms a number of complexes and "
        "organomercury compounds.")
    e['Tl'].description = (
        "Pure, unreacted thallium appears silvery-white and exhibits a "
        "metallic lustre. Upon reacting with air, it begins to turn "
        "bluish-grey and looks like lead. It is very malleable, and can be "
        "cut with a knife. There are two stable isotopes, and four "
        "radioisotopes, Tl-204 being the most stable with a half-life of "
        "3.78 years. Thallium sulphate was used as a rodenticide. Thallium "
        "sulphine's conductivity changes with exposure to infrared light, "
        "this gives it a use in infrared detectors. Discovered by Sir "
        "William Crookes via spectroscopy. Its name comes from the Greek "
        "word thallos, which means green twig. Thallium and its compounds "
        "are toxic and can cause cancer.")
    e['Pb'].description = (
        "Heavy dull grey ductile metallic element, belongs to group 14. "
        "Used in building construction, lead-place accumulators, bullets "
        "and shot, and is part of solder, pewter, bearing metals, type "
        "metals and fusible alloys.")
    e['Bi'].description = (
        "White crystalline metal with a pink tinge, belongs to group 15. "
        "Most diamagnetic of all metals and has the lowest thermal "
        "conductivity of all the elements except mercury. Lead-free "
        "bismuth compounds are used in cosmetics and medical procedures. "
        "Burns in the air and produces a blue flame. In 1753, C.G. Junine "
        "first demonstrated that it was different from lead.")
    e['Po'].description = (
        "Rare radioactive metallic element, belongs to group 16 of the "
        "periodic table. Over 30 known isotopes exist, the most of all "
        "elements. Po-209 has a half-life of 103 years. Possible uses in "
        "heating spacecraft. Discovered by Marie Curie in 1898 in a sample "
        "of pitchblende.")
    e['At'].description = (
        "Radioactive halogen element. Occurs naturally from uranium and "
        "thorium decay. At least 20 known isotopes. At-210, the most "
        "stable, has a half-life of 8.3 hours. Synthesized by nuclear "
        "bombardment in 1940 by D.R. Corson, K.R. MacKenzie and E. Segre "
        "at the University of California.")
    e['Rn'].description = (
        "Colorless radioactive gaseous element, belongs to the noble "
        "gases. Of the twenty known isotopes, the most stable is Rn-222 "
        "with a half-life of 3.8 days. Formed by the radioactive decay of "
        "Radium-226. Radon itself decays into Polonium. Used in "
        "radiotherapy. As a noble gas, it is effectively inert, though "
        "radon fluoride has been synthesized. First isolated in 1908 by "
        "Ramsey and Gray.")
    e['Fr'].description = (
        "Radioactive element, belongs to group 1 of the periodic table. "
        "Found in uranium and thorium ores. The 22 known isotopes are all "
        "radioactive, with the most stable being Fr-223. Its existence was "
        "confirmed in 1939 by Marguerite Perey.")
    e['Ra'].description = (
        "Radioactive metallic transuranic element, belongs to group 2 of "
        "the periodic table. Most stable isotope, Ra-226 has a half-life "
        "of 1602 years, which decays into radon. Isolated from pitchblende "
        "in 1898 Marie and Pierre Curie.")
    e['Ac'].description = (
        "Silvery radioactive metallic element, belongs to group 3 of the "
        "periodic table. The most stable isotope, Ac-227, has a half-life "
        "of 217 years. Ac-228 (half-life of 6.13 hours) also occurs in "
        "nature. There are 22 other artificial isotopes, all radioactive "
        "and having very short half-lives. Chemistry similar to "
        "lanthanumpy. Used as a source of alpha particles. Discovered by "
        "A. Debierne in 1899.")
    e['Th'].description = (
        "Grey radioactive metallic element. Belongs to actinoids. Found "
        "in monazite sand in Brazil, India and the US. Thorium-232 has a "
        "half-life of 1.39x10^10 years. Can be used as a nuclear fuel for "
        "breeder reactors. Thorium-232 captures slow Neutrons and breeds "
        "uranium-233. Discovered by Jons J. Berzelius in 1829.")
    e['Pa'].description = (
        "Radioactive metallic element, belongs to the actinoids. The most "
        "stable isotope, Pa-231 has a half-life of 2.43*10^4 years. At "
        "least 10 other radioactive isotopes are known. No practical "
        "applications are known. Discovered in 1917 by Lise Meitner and "
        "Otto Hahn.")
    e['U'].description = (
        "White radioactive metallic element belonging to the actinoids. "
        "Three natural isotopes, U-238, U-235 and U-234. Uranium-235 is "
        "used as the fuel for nuclear reactors and weapons. Discovered by "
        "Martin H. Klaproth in 1789.")
    e['Np'].description = (
        "Radioactive metallic transuranic element, belongs to the "
        "actinoids. Np-237, the most stable isotope, has a half-life of "
        "2.2*10^6 years and is a by product of nuclear reactors. The other "
        "known isotopes have mass numbers 229 through 236, and 238 through "
        "241. Np-236 has a half-life of 5*10^3 years. First produced by "
        "Edwin M. McMillan and P.H. Abelson in 1940.")
    e['Pu'].description = (
        "Dense silvery radioactive metallic transuranic element, belongs "
        "to the actinoids. Pu-244 is the most stable isotope with a "
        "half-life of 7.6*10^7 years. Thirteen isotopes are known. Pu-239 "
        "is the most important, it undergoes nuclear fission with slow "
        "neutrons and is hence important to nuclear weapons and reactors. "
        "Plutonium production is monitored down to the gram to prevent "
        "military misuse. First produced by Gleen T. Seaborg, Edwin M. "
        "McMillan, J.W. Kennedy and A.C. Wahl in 1940.")
    e['Am'].description = (
        "Radioactive metallic transuranic element, belongs to the "
        "actinoids. Ten known isotopes. Am-243 is the most stable isotope, "
        "with a half-life of 7.95*10^3 years. Discovered by Glenn T. "
        "Seaborg and associates in 1945, it was obtained by bombarding "
        "Uranium-238 with alpha particles.")
    e['Cm'].description = (
        "Radioactive metallic transuranic element. Belongs to actinoid "
        "series. Nine known isotopes, Cm-247 has a half-life of 1.64*10^7 "
        "years. First identified by Glenn T. Seaborg and associates in "
        "1944, first produced by L.B. Werner and I. Perlman in 1947 by "
        "bombarding americium-241 with Neutrons. Named for Marie Curie.")
    e['Bk'].description = (
        "Radioactive metallic transuranic element. Belongs to actinoid "
        "series. Eight known isotopes, the most common Bk-247, has a "
        "half-life of 1.4*10^3 years. First produced by Glenn T. Seaborg "
        "and associates in 1949 by bombarding americium-241 with alpha "
        "particles.")
    e['Cf'].description = (
        "Radioactive metallic transuranic element. Belongs to actinoid "
        "series. Cf-251 has a half life of about 700 years. Nine isotopes "
        "are known. Cf-252 is an intense Neutron source, which makes it an "
        "intense Neutron source and gives it a use in Neutron activation "
        "analysis and a possible use as a radiation source in medicine. "
        "First produced by Glenn T. Seaborg and associates in 1950.")
    e['Es'].description = (
        "Appearance is unknown, however it is most probably metallic and "
        "silver or gray in color. Radioactive metallic transuranic element "
        "belonging to the actinoids. Es-254 has the longest half-life of "
        "the eleven known isotopes at 270 days. First identified by Albert "
        "Ghiorso and associates in the debris of the 1952 hydrogen bomb "
        "explosion. In 1961 the first microgram quantities of Es-232 were "
        "separated. While einsteinium never exists naturally, if a "
        "sufficient amount was assembled, it would pose a radiation "
        "hazard.")
    e['Fm'].description = (
        "Radioactive metallic transuranic element, belongs to the "
        "actinoids. Ten known isotopes, most stable is Fm-257 with a "
        "half-life of 10 days. First identified by Albert Ghiorso and "
        "associates in the debris of the first hydrogen-bomb explosion in "
        "1952.")
    e['Md'].description = (
        "Radioactive metallic transuranic element. Belongs to the "
        "actinoid series. Only known isotope, Md-256 has a half-life of "
        "1.3 hours. First identified by Glenn T. Seaborg, Albert Ghiorso "
        "and associates in 1955. Alternative name Unnilunium has been "
        "proposed. Named after the 'inventor' of the periodic table, "
        "Dmitri Mendeleev.")
    e['No'].description = (
        "Radioactive metallic transuranic element, belongs to the "
        "actinoids. Seven known isotopes exist, the most stable being "
        "No-254 with a half-life of 255 seconds. First identified with "
        "certainty by Albert Ghiorso and Glenn T. Seaborg in 1966. "
        "Unnilbium has been proposed as an alternative name.")
    e['Lr'].description = (
        "Appearance unknown, however it is most likely silvery-white or "
        "grey and metallic. Lawrencium is a synthetic rare-earth metal. "
        "There are eight known radioisotopes, the most stable being Lr-262 "
        "with a half-life of 3.6 hours. Due to the short half-life of "
        "lawrencium, and its radioactivity, there are no known uses for "
        "it. Identified by Albert Ghiorso in 1961 at Berkeley. It was "
        "produced by bombarding californium with boron ions. The name is "
        "temporary IUPAC nomenclature, the origin of the name comes from "
        "Ernest O. Lawrence, the inventor of the cyclotron. If sufficient "
        "amounts of lawrencium were produced, it would pose a radiation "
        "hazard.")
    e['Rf'].description = (
        "Radioactive transactinide element. Expected to have similar "
        "chemical properties to those displayed by hafnium. Rf-260 was "
        "discovered by the Joint Nuclear Research Institute at Dubna "
        "(U.S.S.R.) in 1964. Researchers at Berkeley discovered Unq-257 "
        "and Unq-258 in 1964.")
    e['Db'].description = (
        "Also known as Hahnium, Ha. Radioactive transactinide element. "
        "Half-life of 1.6s. Discovered in 1970 by Berkeley researchers. So "
        "far, seven isotopes have been discovered.")
    e['Sg'].description = (
        "Half-life of 0.9 +/- 0.2 s. Discovered by the Joint Institute "
        "for Nuclear Research at Dubna (U.S.S.R.) in June of 1974. Its "
        "existence was confirmed by the Lawrence Berkeley Laboratory and "
        "Livermore National Laboratory in September of 1974.")
    e['Bh'].description = (
        "Radioactive transition metal. Half-life of approximately 1/500 "
        "s. Discovered by the Joint Institute for Nuclear Research at "
        "Dubna (U.S.S.R.) in 1976. Confirmed by West German physicists at "
        "the Heavy Ion Research Laboratory at Darmstadt.")
    e['Hs'].description = (
        "Radioactive transition metal first synthesized in 1984 by a "
        "German research team led by Peter Armbruster and Gottfried "
        "Muenzenberg at the Institute for Heavy Ion Research at Darmstadt.")
    e['Mt'].description = (
        "Half-life of approximately 5 ms. The creation of this element "
        "demonstrated that fusion techniques could indeed be used to make "
        "new, heavy nuclei. Made and identified by physicists of the Heavy "
        "Ion Research Laboratory, Darmstadt, West Germany in 1982. Named "
        "in honor of Lise Meitner, the Austrian physicist.")
    return e[symbol].description


def sqlite_script():
    """Return SQL script to create sqlite database of elements.

    Examples
    --------
    >>> import sqlite3
    >>> con = sqlite3.connect(':memory:')
    >>> cur = con.executescript(sqlite_script())
    >>> con.commit()
    >>> for r in cur.execute("SELECT name FROM element WHERE number=6"):
    ...     str(r[0])
    'Carbon'
    >>> con.close()

    """
    sql = ["""
        CREATE TABLE "period" (
            "number" TINYINT NOT NULL PRIMARY KEY,
            "label" CHAR NOT NULL UNIQUE,
            "description" VARCHAR(64)
        );
        CREATE TABLE "group" (
            "number" TINYINT NOT NULL PRIMARY KEY,
            "label" VARCHAR(8) NOT NULL,
            "description" VARCHAR(64)
        );
        CREATE TABLE "block" (
            "label" CHAR NOT NULL PRIMARY KEY,
            "description" VARCHAR(64)
        );
        CREATE TABLE "series" (
            "id" TINYINT NOT NULL PRIMARY KEY,
            "label"  VARCHAR(32) NOT NULL,
            "description" VARCHAR(256)
        );
        CREATE TABLE "element" (
            "number" TINYINT NOT NULL PRIMARY KEY,
            "symbol" VARCHAR(2) UNIQUE NOT NULL,
            "name" VARCHAR(16) UNIQUE NOT NULL,
            "period" TINYINT NOT NULL,
            --FOREIGN KEY("period") REFERENCES "period"(number),
            "group" TINYINT NOT NULL,
            --FOREIGN KEY("group") REFERENCES "group"(number),
            "block" CHAR NOT NULL,
            --FOREIGN KEY("block") REFERENCES "block"(label),
            "series" TINYINT NOT NULL,
            --FOREIGN KEY("series") REFERENCES "series"(id),
            "mass" REAL NOT NULL,
            "eleneg" REAL,
            "covrad" REAL,
            "atmrad" REAL,
            "vdwrad" REAL,
            "tboil" REAL,
            "tmelt" REAL,
            "density" REAL,
            "eleaffin" REAL,
            "eleconfig" VARCHAR(32),
            "oxistates" VARCHAR(32),
            "description" VARCHAR(2048)
        );
        CREATE TABLE "isotope" (
            "element" TINYINT NOT NULL,
            --FOREIGN KEY ("element") REFERENCES "element"("number"),
            "massnum" TINYINT NOT NULL,
            "mass" REAL NOT NULL,
            "abundance" REAL NOT NULL,
            PRIMARY KEY ("element", "massnum")
        );
        CREATE TABLE "eleconfig" (
            "element" TINYINT NOT NULL,
            --FOREIGN KEY ("element") REFERENCES "element"("number"),
            "shell" TINYINT NOT NULL,
            --FOREIGN KEY ("shell") REFERENCES "period"("number"),
            "subshell" CHAR NOT NULL,
            --FOREIGN KEY ("subshell") REFERENCES "block"("label"),
            "count" TINYINT,
            PRIMARY KEY ("element", "shell", "subshell")
        );
        CREATE TABLE "ionenergy" (
            "element" TINYINT NOT NULL,
            --FOREIGN KEY ("element") REFERENCES "element"("number"),
            "number" TINYINT NOT NULL,
            "energy" REAL NOT NULL,
            PRIMARY KEY ("element", "number")
        );
    """]

    for key, label in PERIODS.items():
        sql.append("""INSERT INTO "period" VALUES (%i, '%s', NULL);""" % (
            key, label))

    for key, (label, descr) in GROUPS.items():
        sql.append("""INSERT INTO "group" VALUES (%i, '%s', '%s');""" % (
            key, label, descr))

    for data in BLOCKS.items():
        sql.append("""INSERT INTO "block" VALUES ('%s', '%s');""" % data)

    for series in sorted(SERIES):
        sql.append("""INSERT INTO "series" VALUES (%i, '%s', '');""" % (
            series, SERIES[series]))

    for ele in ELEMENTS:
        sql.append("""
        INSERT INTO "element" VALUES (%i, '%s', '%s', %i, %i, '%s', %i,
            %.10f, %.4f, %.4f, %.4f, %.4f,
            %.4f, %.4f, %.4f, %.8f,
            '%s', '%s',
            '%s'
        );""" % (
            ele.number, ele.symbol, ele.name, ele.period, ele.group,
            ele.block, ele.series, ele.mass, ele.eleneg,
            ele.covrad, ele.atmrad, ele.vdwrad, ele.tboil, ele.tmelt,
            ele.density, ele.eleaffin, ele.eleconfig, ele.oxistates,
            word_wrap(
                ele.description.replace("'", "\'\'").replace("\"", "\"\""),
                linelen=74, indent=0, joinstr="\n ")))

    for ele in ELEMENTS:
        for iso in ele.isotopes.values():
            sql.append(
                """INSERT INTO "isotope" VALUES (%i, %i, %.10f, %.8f);""" % (
                    ele.number, iso.massnumber, iso.mass, iso.abundance))

    for ele in ELEMENTS:
        for (shell, subshell), count in ele.eleconfig_dict.items():
            sql.append(
                """INSERT INTO "eleconfig" VALUES (%i, %i, '%s', %i);""" % (
                    ele.number, shell, subshell, count))

    for ele in ELEMENTS:
        for i, ionenergy in enumerate(ele.ionenergy):
            sql.append("""INSERT INTO "ionenergy" VALUES (%i, %i, %.4f);""" % (
                ele.number, i + 1, ionenergy))

    return '\n'.join(sql).replace("        ", "")


def word_wrap(text, linelen=80, indent=0, joinstr="\n"):
    """Return string, word wrapped at linelen."""
    if len(text) < linelen:
        return text
    result = []
    line = []
    llen = -indent
    for word in text.split():
        llen += len(word) + 1
        if llen < linelen:
            line.append(word)
        else:
            result.append(" ".join(line))
            line = [word]
            llen = len(word)
    if line:
        result.append(" ".join(line))
    return joinstr.join(result)


if __name__ == "__main__":
    for ele in ELEMENTS:
        print(repr(ele), '\n')
    import doctest
    doctest.testmod(verbose=False)

