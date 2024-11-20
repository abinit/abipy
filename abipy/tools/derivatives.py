# coding: utf-8
"""Tools for computing derivatives by finite differences."""
import numpy as np

from monty.collections import dict2namedtuple

__all__ = [
    "finite_diff",
]


def rearr(array):
    return np.array(array, dtype=float)

# This table contains the coefficients of the central differences, for several order of accuracy: [1]
# See http://en.wikipedia.org/wiki/Finite_difference_coefficients
# Derivative Accuracy -4 -3 -2 -1 0 1 2 3 4


central_fdiff_weights = {
1: {
    2: rearr([-1/2, 0, 1/2]),
    4: rearr([1/12, -2/3, 0, 2/3, -1/12]),
    6: rearr([-1/60, 3/20, -3/4, 0, 3/4, -3/20, 1/60]),
    8: rearr([1/280, -4/105, 1/5, -4/5, 0, 4/5, -1/5, 4/105, -1/280]),
    },
2: {
    2: rearr([1, -2, 1]),
    4: rearr([-1/12, 4/3, -5/2, 4/3, -1/12]),
    6: rearr([1/90, -3/20, 3/2, -49/18, 3/2, -3/20, 1/90]),
    8: rearr([-1/560, 8/315, -1/5, 8/5, -205/72, 8/5, -1/5, 8/315, -1/560]),
    },
3: {
    2: rearr([-1/2, 1, 0, -1, 1/2]),
    4: rearr([1/8, -1, 13/8, 0, -13/8, 1, -1/8]),
    6: rearr([-7/240, 3/10, -169/120, 61/30, 0, -61/30, 169/120, -3/10, 7/240]),
    },
4: {
    2: rearr([1, -4, 6, -4, 1]),
    4: rearr([-1/6, 2, -13/2, 28/3, -13/2, 2, -1/6]),
    6: rearr([7/240, -2/5, 169/60, -122/15, 91/8, -122/15, 169/60, -2/5, 7/240]),
    },
}


#Derivative  Accuracy   0 1 2 3 4 5 6 7 8
forward_fdiff_weights = {
1: {
    1: rearr([-1, 1]),
    2: rearr([-3/2, 2, -1/2]),
    3: rearr([-11/6, 3, -3/2, 1/3]),
    4: rearr([-25/12, 4, -3, 4/3, -1/4]),
    5: rearr([-137/60, 5, -5, 10/3, -5/4, 1/5]),
    6: rearr([-49/20, 6, -15/2, 20/3, -15/4, 6/5, -1/6]),
},
2: {
    1: rearr([1, -2, 1]),
    2: rearr([2, -5, 4, -1]),
    3: rearr([35/12, -26/3, 19/2, -14/3, 11/12]),
    4: rearr([15/4,  -77/6, 107/6, -13, 61/12, -5/6]),
    5: rearr([203/45, -87/5, 117/4, -254/9, 33/2, -27/5, 137/180]),
    6: rearr([469/90, -223/10, 879/20, -949/18, 41, -201/10, 1019/180, -7/10]),
},
3: {
    1: rearr([-1, 3, -3, 1]),
    2: rearr([-5/2, 9, -12, 7, -3/2]),
    3: rearr([-17/4, 71/4, -59/2, 49/2, -41/4, 7/4]),
    4: rearr([-49/8, 29, -461/8, 62, -307/8, 13, -15/8]),
    5: rearr([-967/120, 638/15, -3929/40, 389/3, -2545/24, 268/5, -1849/120, 29/15]),
    6: rearr([-801/80, 349/6, -18353/120, 2391/10, -1457/6, 4891/30, -561/8, 527/30, -469/240]),
},
4: {
    1: rearr([1, -4, 6, -4, 1]),
    2: rearr([3, -14, 26, -24, 11, -2]),
    3: rearr([35/6, -31, 137/2, -242/3, 107/2, -19, 17/6]),
    4: rearr([28/3, -111/2, 142, -1219/6, 176, -185/2,  82/3, -7/2]),
    5: rearr([1069/80, -1316/15, 15289/60, -2144/5, 10993/24, -4772/15, 2803/20, -536/15, 967/240]),
},
}
del rearr

# To get the coefficients of the backward approximations,
# give all odd derivatives listed in the table the opposite sign,
# whereas for even derivatives the signs stay the same.
backward_fdiff_weights = d = {}

for ord, v in forward_fdiff_weights.items():
    d[ord] = {}
    for accuracy, weights in v.items():
        d[ord][accuracy] = ((-1)**ord) * weights[-1::-1]


def finite_diff(arr, h, order=1, acc=4, index=None):
    """
    Compute the derivative of order `order` by finite difference.
    For each point in arr, the function tries to use central differences
    and fallbacks to forward/backward approximations for points that are close to the extrema.
    Note that high accuracy levels can fail and raise `ValueError` if not enough points are available in `arr`.

    Args:
        arr: Input array with y-values.
        h: Spacing along x
        order: Derivative order
        acc: accuracy level.
        index: If not None, gives the index of the single element in arr where the derivative is wanted.
            In this case a namedtuple with the derivative, the number of points used and the mode is returned

    Return:
        numpy array or (value, npts, mode) if index is not None .
    """
    arr = np.asarray(arr)

    if np.iscomplexobj(arr):
        raise ValueError("Derivatives of complex functions are not supported!")

    # Retrieve weights.
    try:
        centr_ws = central_fdiff_weights[order][acc]
    except KeyError:
        raise ValueError("Centeral diff weights for order: %s, and accuracy: %s are missing!" % (order, acc))

    npsum = np.sum
    ders = np.empty(arr.shape)
    n = len(arr)
    cpad = len(centr_ws) // 2

    for i in range(n):
        if index is not None and i != index: continue
        start = i - cpad
        stop = i + cpad + 1

        if start >= 0 and stop <= n:
            # Can do central difference.
            ders[i] = npsum(centr_ws * arr[start:stop])
            npts = len(centr_ws)
            mode = "central"

        elif start < 0:
            # Try forward.
            forw_ws = forward_fdiff_weights[order][acc]
            stop = i + len(forw_ws)
            if stop > n:
                raise ValueError(
                        ("\n\tDon't have enough points for index: %s in array of lenght: %s\n" +
                         "\tto compute forward finite difference with order: %s, and acc: %s (num_weights: %s)\n" +
                         "\tDecrease acc or increase the number of sampling points.") % (i, n, order, acc, len(forw_ws)))
            ders[i] = npsum(forw_ws * arr[i:stop])
            npts = len(forw_ws)
            mode = "forward"

        elif stop > n:
            # Try backward.
            back_ws = backward_fdiff_weights[order][acc]
            start = i - len(back_ws) + 1
            if start < 0:
                raise ValueError(
                    ("\n\tDon't have enough points for index: %s in array of length: %s\n" +
                    "\tto compute backward finite difference with order: %s, and acc: %s (num_weights: %s)\n" +
                    "\tDecrease acc or increase the number of sampling points.") % (i, n, order, acc, len(back_ws)))
            ders[i] = npsum(back_ws * arr[start:i+1])
            npts = len(back_ws)
            mode = "backward"

    if index is None:
        return ders / (h ** order)
    else:
        return dict2namedtuple(value=ders[index] / (h ** order), npts=npts, mode=mode)
