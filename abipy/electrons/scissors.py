from __future__ import print_function, division

import numpy as np

from abipy.core.constants import Ha_eV

__all__ = [
    "ScissorsOperator",
]


class ScissorsOperator(object):
    """
    This object represents an energy-dependent scissors operator.

    The standard way to instanciate this object is via the methods
    provided by the factory class ...

    Once the instance has been create, one can correct the band structure by 
    calling the apply method.
    """

    def __init__(self, func_list, domains, bounds=None):
        """
        Args:
            func_list:
                List of functions
            domains:
                List of domains of each function.
            bounds:
                Boundaries

        .. note:

            - Domains should not overlap, cover e0mesh, and given in increasing order.

            - Holes are permitted but the interpolation will raise an exception if the
              eigenvalue falls inside the hole.
        """
        # TODO Add consistency check.
        self.funclist = func_list
        self.domains = np.atleast_2d(domains)
        assert len(self.funclist) == len(self.domains)

        # Treat the out-of-boundary conditions.
        blow, bhigh  = "c", "c"
        if bounds is not None:
            blow  = bounds[0][0]
            bhigh = bounds[0][1]

        if blow.lower() == "c":
            try:
                self.func_low = lambda x: float(bounds[0][1])
            except:
                x_low = self.domains[0,0]
                fx_low = func_list[0](x_low)
                self.func_low = lambda x: fx_low
        else:
            raise NotImplementedError("Only constant boundaries are implemented")

        if bhigh.lower() == "c":
            try:
                self.func_high = lambda x: float(bounds[1][1])
            except:
                x_high  = self.domains[1, -1]
                fx_high = func_list[-1](x_high)
                self.func_high = lambda x: fx_high
        else:
            raise NotImplementedError("Only constant boundaries are implemented")

        # Init the internal counter.
        self.out_bounds = np.zeros(3, np.int)

    def apply(self, eig):
        """Correct the eigenvalue eig (eV units)."""
        domains = self.domains

        if eig < domains[0,0]:
            print("left ", eig, domains[0,0])
            self.out_bounds[0] += 1
            return self.func_low(eig)

        if eig > domains[-1,1]:
            print("right ", eig, domains[-1,1])
            self.out_bounds[1] += 1
            return self.func_high(eig)

        for idx, dms in enumerate(domains):
            if dms[1] >= eig >= dms[0]:
                return self.funclist[idx](eig)

        self.out_bounds[2] += 1
        raise ValueError("Cannot bracket eig %s" % eig)

#########################################################################################

