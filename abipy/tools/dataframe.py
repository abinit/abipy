# coding: utf-8
""""""
from __future__ import division, print_function, absolute_import

from collections import namedtuple

import numpy as np
import pandas as pd


class PrettyTable(pd.DataFrame):
    """
    Extends :class:`DataFrame` with methods for numerical analysis.
    """
    def quadfit(self, xname=None, yname=None):
        """
        Quadratic fit. Return namedtuple with the parameters of the fix a*x**2 + b*x + c,
        the position of the minimum in x0 and the value of the minimum in y0
        """
        xvals, yvals = self[xname], self[yname]

        a, b, c = np.polyfit(xvals, yvals, 2)
        x0 = -b/(2*a)
        y0 = a*x0**2 + b*x0 + c

        return namedtuple("quadfit_results", "a, b, c, x0, y0")(a=a, b=b, c=c, x0=x0, y0=y0)
