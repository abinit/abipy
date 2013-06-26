"""Tools for k-points."""
from __future__ import print_function, division

from numpy import allclose,  where, around, asarray #, int, all, abs,

__all__ = [
    "issamek",
    "wrap_to_ws",
    "wrap_to_bz",
]

def isinteger(x, atol=1e-08):
    """
    True if all x is integer within the absolute tolerance atol.

    >>> isinteger([1., 2.])
    True
    >>> isinteger(1.01, atol=0.011)
    True
    >>> isinteger([1.01, 2])
    False
    """
    int_x = around(x)
    return allclose(int_x, x, atol=atol)


def issamek(k1, k2, atol=1e-08):
    """
    True if k1 and k2 are equal modulo a lattice vector.

    Examples

    >>> issamek([1,1,1], [0,0,0])
    True
    >>> issamek([1.1,1,1], [0,0,0], atol=0.1)
    True
    >>> issamek(0.00003, 1)
    False
    """
    kdiff = asarray(k1) - asarray(k2)
    return isinteger(kdiff, atol=atol)


def wrap_to_ws(x):
    """
    Transforms x in its corresponding reduced number in the interval ]-1/2,1/2].
    """
    w = x % 1
    return where(w > 0.5, w-1.0, w)


def wrap_to_bz(x):
    """
    Transforms x in its corresponding reduced number in the interval [0,1[."
    """
    return x % 1

#########################################################################################

if __name__ == "__main__":
    import doctest
    doctest.testmod()
