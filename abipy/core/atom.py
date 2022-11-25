# coding: utf-8
"""This module provides objects and helper functions for atomic calculations."""
from __future__ import annotations

import collections
import numpy as np

from io import StringIO
from typing import Any, List, Union, Optional, Iterable, Tuple
from abipy.data import nist_database
from scipy.interpolate import UnivariateSpline
from scipy.integrate import cumtrapz

from monty.functools import lazy_property

__version__ = "0.1"
__author__ = "Matteo Giantomassi"
__maintainer__ = "Matteo Giantomassi"

__all__ = [
    "NlkState",
    "QState",
    "AtomicConfiguration",
    "RadialFunction",
    "RadialWaveFunction",
]

char2l = {
    "s": 0,
    "p": 1,
    "d": 2,
    "f": 3,
    "g": 4,
    "h": 5,
    "i": 6,
}

l2char = {
    0: "s",
    1: "p",
    2: "d",
    3: "f",
    4: "g",
    5: "h",
    6: "i",
}

# Accept strings as keys as well.
l2char.update({str(l): l2char[l] for l in l2char})


def _asl(obj: Any) -> int:
    try:
        return char2l[obj]
    except KeyError:
        return int(obj)


def states_from_string(confstr: str) -> List[QState]:
    """
    Parse a string with an atomic configuration and build a list of `QState` instance.
    """
    states = []
    tokens = confstr.split()

    if tokens[0].startswith("["):
        noble_gas = AtomicConfiguration.neutral_from_symbol(tokens.pop(0)[1:-1])
        states = noble_gas.states

    states.extend(parse_orbtoken(t) for t in tokens)
    return states


def parse_orbtoken(orbtoken: str) -> QState:
    import re
    m = re.match(r"(\d+)([spdfghi]+)(\d+)", orbtoken.strip())
    if m:
        return QState(n=m.group(1), l=m.group(2), occ=m.group(3))

    raise ValueError("Don't know how to interpret %s" % str(orbtoken))


class NlkState(collections.namedtuple("NlkState", "n, l, k")):
    """
    Named tuple storing (n, l) or (n, l, k) if relativistic pseudos.
    k is set to None if we are in non-relativistic mode.
    """

    # The relativistic code utilizes a new main program, oncvpsp_r, rather than
    # complicate oncvpsp with lots of branching.  Many arrays have acquired an
    # additional final index, and most loops over angular momenta include inner
    # loops over this "ikap=1,2" index referring to the additional quantum number of
    # the radial Dirac equations kappa (kap) =l, -(l+1) for j=l -/+ 1/2.

    def __new__(cls, n: int, l: int, k: Optional[int] = None):
        """Extends super.__new__ adding type conversion and default values."""
        if k is not None:
            if k not in (1, 2):
                raise ValueError(f"Invalid k index: {k} for l: {l}")
            if l == 0 and k != 1:
                raise ValueError(f"When l is 0, k must be 1 while it is: {k}")

        #if n <= 0:
        #    raise ValueError(f"Invalid value for n: {n}")
        if l < 0:
            raise ValueError(f"Invalid value for l: {l}")

        return super().__new__(cls, n, l, k)

    @classmethod
    def from_nlkap(cls, n: int, l: int, kap: Union[int, None]) -> NlkState:
        k = None
        if kap is not None:
            #if(ikap==1) kap=-(ll+1)
            #if(ikap==2) kap=  ll
            k = 1
            if l != 0:
                k = {-(l + 1): 1, l:2}[kap]

        return cls(n=n, l=l, k=k)

    def __repr__(self) -> str:
        return f"n={self.n}, l={self.l}" if self.k is None else f"n={self.n}, l={self.l}, k={self.k}"

    def __str__(self) -> str:
        lc = l2char[self.l]
        if self.k is None:
            return f"{self.n}{lc}"  # e.g. 2s
        else:
            return f"{self.n}{lc}{self.ksign}"  # e.g. 2s+

    @lazy_property
    def latex(self):
        lc = l2char[self.l]
        if self.k is None:
            return f"${self.n}{lc}$"  # e.g. 2s
        else:
            return f"${self.n}{lc}^{self.ksign}$"  # e.g. 2s^+

    @lazy_property
    def latex_l(self):
        lc = l2char[self.l]
        if self.k is None:
            return f"${lc}$"  # e.g. s
        else:
            return f"${lc}^{self.ksign}$"  # e.g. s^+

    @lazy_property
    def ksign(self) -> str:
        return {1: "+", 2: "-"}[self.k]

    @lazy_property
    def j(self) -> int:
        """Total angular momentum"""
        l = self.l
        if self.k is None: return l
        return l - 1/2 if self.k == l else l + 1/2

    #@property
    #def to_dict(self) -> dict:
    #    return {"n": self.n, "l": self.l, "k": self.k}


class QState(collections.namedtuple("QState", "n, l, occ, eig, j, s")):
    """
    This object collects the set of quantum numbers and the physical properties
    associated to a single electron in a spherically symmetric atom.

    .. attributes:

        n: Principal quantum number.
        l: Angular momentum.
        occ: Occupancy of the atomic orbital.
        eig: Eigenvalue of the atomic orbital.
        j: J quantum number. None if spin is a good quantum number.
        s: Spin polarization. None if spin is not taken into account.
    """
    # TODO
    # Spin +1, -1 or 1,2 or 0,1?

    def __new__(cls, n: int, l: int, occ: float,
                eig: Optional[float] = None,
                j: Optional[int] = None,
                s: Optional[int] = None):
        """
        Extends super.__new__ adding type conversion and default values.
        """
        eig = float(eig) if eig is not None else eig
        j = int(j) if j is not None else j
        s = int(s) if s is not None else s
        return super().__new__(cls, int(n), _asl(l), float(occ), eig=eig, j=j, s=s)

    @property
    def has_j(self) -> bool:
        return self.j is not None

    @property
    def has_s(self) -> bool:
        return self.s is not None


class AtomicConfiguration:
    """
    Atomic configuration of an all-electron atom.
    """
    def __init__(self, Z: int, states: List[QState]) -> None:
        """
        Args:
            Z: Atomic number.
            states: List of :class:`QState` instances.
        """
        self.Z = Z
        self.states = states

    @classmethod
    def from_string(cls, Z: int, string: str,
                    has_s: bool = False, has_j: bool = False) -> AtomicConfiguration:
        if not has_s and not has_j:
            # Ex: [He] 2s2 2p3
            states = states_from_string(string)
        else:
            raise NotImplementedError("")

        return cls(Z, states)

    @classmethod
    def neutral_from_symbol(cls, symbol: Union[str, int]) -> AtomicConfiguration:
        """
        symbol: str or int
            Can be a chemical symbol (str) or an atomic number (int).
        """
        entry = nist_database.get_neutral_entry(symbol)
        states = [QState(n=s[0], l=s[1], occ=s[2]) for s in entry.states]
        return cls(entry.Z, states)

    def __str__(self) -> str:
        lines = ["%s: " % self.Z]
        lines += [str(state) for state in self]
        return "\n".join(lines)

    def __len__(self) -> int:
        return len(self.states)

    def __iter__(self) -> Iterable:
        return self.states.__iter__()

    def __eq__(self, other: AtomicConfiguration) -> bool:
        if len(self.states) != len(other.states):
            return False

        return (self.Z == other.Z and
                all(s1 == s2 for s1, s2 in zip(self.states, other.states)))

    def __ne__(self, other: AtomicConfiguration) -> bool:
        return not self == other

    def copy(self) -> AtomicConfiguration:
        """Shallow copy of self."""
        return AtomicConfiguration(self.Z, [s for s in self.states])

    @property
    def symbol(self) -> str:
        """Chemical symbol"""
        return nist_database.symbol_from_Z(self.Z)

    @property
    def spin_mode(self) -> str:
        """
        unpolarized: Spin-unpolarized calculation.
        polarized: Spin-polarized calculation.
        """
        for state in self:
            # FIXME
            if state.s is not None and state.s == 2:
                return "polarized"
        return "unpolarized"

    @property
    def echarge(self) -> float:
        """Electronic charge (float < 0 )."""
        return -sum(state.occ for state in self)

    @property
    def isneutral(self) -> bool:
        """True if self is a neutral configuration."""
        return abs(self.echarge + self.Z) < 1.e-8

    def add_state(self, **qnumbers):
        """Add a list of :class:`QState` instances to self."""
        self._push(QState(**qnumbers))

    def remove_state(self, **qnumbers):
        """Remove a quantum state from self."""
        self._pop(QState(**qnumbers))

    def _push(self, state):
        # TODO check that ordering in the input does not matter!
        if state in self.states:
            raise ValueError("state %s is already in self" % str(state))
        self.states.append(state)

    def _pop(self, state):
        try:
            self.states.remove(state)
        except ValueError:
            raise


class RadialFunction:
    """
    A RadialFunction has a name, a radial mesh and values defined on this mesh.
    """

    def __init__(self, name: str, rmesh, values):
        """
        Args:
            name: Name of the function (string).
            rmesh: Iterable with the points of the radial mesh.
            values: Iterable with the values of the function on the radial mesh.
        """
        self.name = name
        self.rmesh = np.ascontiguousarray(rmesh)
        self.values = np.ascontiguousarray(values)
        assert len(self.rmesh) == len(self.values)

    @classmethod
    def from_filename(cls, filename: str, rfunc_name=None, cols=(0, 1)) -> RadialFunction:
        """
        Initialize the object reading data from filename (txt format).

        Args:
            filename: Path to the file containing data.
            rfunc_name: Optional name for the RadialFunction (defaults to filename)
            cols: List with the index of the columns containing the radial mesh and the values.
        """
        data = np.loadtxt(filename)
        rmesh, values = data[:,cols[0]], data[:,cols[1]]
        name = filename if rfunc_name is None else rfunc_name
        return cls(name, rmesh, values)

    def __len__(self) -> int:
        return len(self.values)

    def __iter__(self) -> Iterable:
        """Iterate over (rpoint, value)."""
        return iter(zip(self.rmesh, self.values))

    def __getitem__(self, rslice):
        return self.rmesh[rslice], self.values[rslice]

    def __repr__(self) -> str:
        return "<%s, name=%s at %s>" % (self.__class__.__name__, self.name, id(self))

    def __str__(self) -> str:
        stream = StringIO()
        self.pprint(stream=stream)
        return stream.getvalue()

    #def __add__(self, other):
    #def __sub__(self, other):
    #def __mul__(self, other):

    def __abs__(self) -> RadialFunction:
        return self.__class__(self.rmesh, np.abs(self.values))

    @property
    def to_dict(self) -> dict:
        return dict(
            name=str(self.name),
            rmesh=list(self.rmesh),
            values=list(self.values),
        )

    def pprint(self, what: str = "rmesh+values", stream=None) -> None:
        """pprint method (useful for debugging)"""
        from pprint import pprint
        if "rmesh" in what:
            pprint("rmesh:", stream=stream)
            pprint(self.rmesh, stream=stream)

        if "values" in what:
            pprint("values:", stream=stream)
            pprint(self.values, stream=stream)

    @property
    def rmax(self) -> float:
        """Outermost point of the radial mesh."""
        return self.rmesh[-1]

    @property
    def rsize(self) -> float:
        """Size of the radial mesh."""
        return len(self.rmesh)

    @property
    def minmax_ridx(self) -> Tuple[int, int]:
        """
        Returns the indices of the values in a list with the maximum and minimum value.
        """
        minimum = min(enumerate(self.values), key=lambda s: s[1])
        maximum = max(enumerate(self.values), key=lambda s: s[1])
        return minimum[0], maximum[0]

    @property
    def inodes(self) -> List[int]:
        """"
        List with the index of the nodes of the radial function.
        """
        inodes = []
        for i in range(len(self.values)-1):
            if self.values[i] * self.values[i+1] <= 0:
                inodes.append(i)
        return inodes

    @property
    def spline(self):
        """Cubic spline."""
        try:
            return self._spline
        except AttributeError:
            self._spline = UnivariateSpline(self.rmesh, self.values, s=0)
            return self._spline

    @property
    def roots(self):
        """Return the zeros of the spline."""
        return self.spline.roots()

    def derivatives(self, r):
        """Return all derivatives of the spline at the point r."""
        return self.spline.derivatives(r)

    def integral(self) -> RadialFunction:
        r"""
        Cumulatively integrate y(x) using the composite trapezoidal rule.

        Returns:
            `Function1d` with :math:`\int y(x) dx`
        """
        integ = cumtrapz(self.values, x=self.rmesh)
        pad_intg = np.zeros(len(self.values))
        pad_intg[1:] = integ

        return self.__class__(self.rmesh, pad_intg)

    def integral3d(self, a=None, b=None):
        """
        Return definite integral of the spline of (r**2 values**2) between two given points a and b

        Args:
            a: First point. rmesh[0] if a is None
            b: Last point. rmesh[-1] if a is None
        """
        a = self.rmesh[0] if a is None else a
        b = self.rmesh[-1] if b is None else b
        r2v2_spline = UnivariateSpline(self.rmesh, (self.rmesh * self.values) ** 2, s=0)

        return r2v2_spline.integral(a, b)

    def ifromr(self, rpoint):
        """
        The index of the point in the radial mesh.
        """
        for (i, r) in enumerate(self.rmesh):
            if r > rpoint:
                return i - 1

        if rpoint == self.rmesh[-1]:
            return len(self.rmesh)
        else:
            raise ValueError("Cannot find %s in rmesh" % rpoint)

    def ir_small(self, abs_tol: float = 0.01) -> int:
        """
        Returns the rightmost index where the abs value of the wf becomes greater than abs_tol

        Args:
            abs_tol: Absolute tolerance.

        .. warning::

            Assumes that self.values are tending to zero for r --> infinity.
        """
        for i in range(len(self.rmesh)-1, -1, -1):
            if abs(self.values[i]) > abs_tol:
                break
        return i

    def r2f2_integral(self):
        """
        Cumulatively integrate r**2 f**2(r) using the composite trapezoidal rule.
        """
        integ = cumtrapz(self.rmesh**2 * self.values**2, x=self.rmesh)
        pad_intg = np.zeros(len(self))
        pad_intg[1:] = integ

        return pad_intg

    def r2f_integral(self):
        """
        Cumulatively integrate r**2 f(r) using the composite trapezoidal rule.
        """
        integ = cumtrapz(self.rmesh**2 * self.values, x=self.rmesh)
        pad_intg = np.empty(len(self))
        pad_intg[1:] = integ
        return pad_intg

    def get_intr2j0(self, ecut: float, numq: float = 3001):
        r"""
        Compute 4\pi\int[(\frac{\sin(2\pi q r)}{2\pi q r})(r^2 n(r))dr].
        """
        qmax = np.sqrt(ecut / 2) / np.pi
        qmesh = np.linspace(0, qmax, num=numq, endpoint=True)
        outs = np.empty(len(qmesh))

        # Treat q == 0. Note that rmesh[0] > 0
        f = 4 * np.pi * self.rmesh**2 * self.values
        outs[0] = cumtrapz(f, x=self.rmesh)[-1]

        for i, q in enumerate(qmesh[1:]):
            twopiqr = 2 * np.pi * q * self.rmesh
            f = 4 * np.pi * (np.sin(twopiqr) / twopiqr) * self.rmesh**2 * self.values
            outs[i+1] = cumtrapz(f, x=self.rmesh)[-1]

        from abipy.core.func1d import Function1D
        ecuts = 2 * np.pi**2 * qmesh**2
        return Function1D(ecuts, outs)


class RadialWaveFunction(RadialFunction):
    """
    Extends :class:`RadialFunction` adding info on the set of quantum numbers.
    and methods specialized for electronic wavefunctions.
    """
    TOL_BOUND = 1.e-10

    def __init__(self, nlk: NlkState, name: str, rmesh, values):
        super().__init__(name, rmesh, values)
        self.nlk = nlk

    @lazy_property
    def isbound(self) -> bool:
        """True if self is a bound state."""
        back = min(10, len(self))
        return np.all(np.abs(self.values[-back:]) < self.TOL_BOUND)

    #@property
    #def to_dict(self) -> dict:
    #    d = super().to_dict
    #    d.update(self.nlk.to_dict)
    #    return d

