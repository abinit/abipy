# coding: utf-8
"""Numeric tools."""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
import bisect as bs

from monty.collections import dict2namedtuple
from abipy.tools import duck

#########################################################################################
# Array tools
#########################################################################################

def transpose_last3dims(arr):
    """
    Transpose the last three dimensions of arr: (...,x,y,z) --> (...,z,y,x).
    """
    axes = np.arange(arr.ndim)
    axes[-3:] = axes[::-1][:3]

    view = np.transpose(arr, axes=axes)
    return np.ascontiguousarray(view)


def add_periodic_replicas(arr):
    """
    Returns a new array of shape=(..., nx+1,ny+1,nz+1) with redundant data points.

    Periodicity in enforced only on the last three dimensions.
    """
    ishape, ndim = np.array(arr.shape), arr.ndim
    oshape = ishape.copy()

    if ndim == 1:
        oarr = np.empty(ishape+1, dtype=arr.dtype)
        oarr[:-1] = arr
        oarr[-1] = arr[0]

    elif ndim == 2:
        oarr = np.empty(oshape+1, dtype=arr.dtype)
        oarr[:-1,:-1] = arr
        oarr[-1,:-1] = arr[0,:]
        oarr[:-1,-1] = arr[:,0]
        oarr[-1,-1] = arr[0,0]

    else:
        # Add periodic replica along the last three directions.
        oshape[-3:] = oshape[-3:] + 1
        oarr = np.empty(oshape, dtype=arr.dtype)

        for x in range(ishape[-3]):
            for y in range(ishape[-2]):
                oarr[..., x, y, :-1] = arr[..., x, y, :]
                oarr[..., x, y, -1] = arr[..., x, y, 0]

            oarr[..., x, y + 1, :-1] = arr[..., x, 0, :]
            oarr[..., x, y + 1, -1] = arr[..., x, 0, 0]

        oarr[..., x + 1, :, :] = oarr[..., 0, :, :]

    return oarr


def data_from_cplx_mode(cplx_mode, arr, tol=None):
    """
    Extract the data from the numpy array ``arr`` depending on the values of ``cplx_mode``.

    Args:
        cplx_mode: Possible values in ("re", "im", "abs", "angle")
            "re" for the real part,
            "im" for the imaginary part.
            "all" for both re and im.
            "abs" means that the absolute value of the complex number is shown.
            "angle" will display the phase of the complex number in radians.
	tol: If not None, values below tol are set to zero. Cannot be used with "angle"
    """
    if cplx_mode == "re":
        val = arr.real
    elif cplx_mode == "im":
        val = arr.imag
    elif cplx_mode == "all":
        val = arr
    elif cplx_mode == "abs":
        val = np.abs(arr)
    elif cplx_mode == "angle":
        val = np.angle(arr, deg=False)
        if tol is not None:
            raise ValueError("Tol cannot be used with cplx_mode = angle")
    else:
        raise ValueError("Unsupported mode `%s`" % str(cplx_mode))

    return val if tol is None else np.where(np.abs(val) > tol, val, 0)


def is_diagonal(matrix, atol=1e-12):
    """
    Return True if matrix is diagonal.
    """
    m = matrix.copy()
    np.fill_diagonal(m, 0)

    if issubclass(matrix.dtype.type, np.integer):
        return np.all(m == 0)
    else:
        return np.all(np.abs(m) <= atol)


#########################################################################################
# Tools to facilitate iterations
#########################################################################################


def alternate(*iterables):
    """
    [a[0], b[0], ... , a[1], b[1], ..., a[n], b[n] ...]
    >>> alternate([1,4], [2,5], [3,6])
    [1, 2, 3, 4, 5, 6]
    """
    items = []
    for tup in zip(*iterables):
        items.extend([item for item in tup])
    return items


def iflat(iterables):
    """
    Iterator over all elements of a nested iterable. It's recursive!

    >>> list(iflat([[0], [1,2, [3,4]]]))
    [0, 1, 2, 3, 4]
    """
    for item in iterables:
        if not hasattr(item, "__iter__"):
            yield item
        else:
            # iterable object.
            for it in iflat(item):
                yield it


def grouper(n, iterable, fillvalue=None):
    """
    >>> assert grouper(3, "ABCDEFG", "x") == [('A', 'B', 'C'), ('D', 'E', 'F'), ('G', 'x', 'x')]
    >>> assert grouper(3, [1, 2, 3, 4]) == [(1, 2, 3), (4, None, None)]
    """
    # https://stackoverflow.com/questions/434287/what-is-the-most-pythonic-way-to-iterate-over-a-list-in-chunks/434411#434411
    try:
        from itertools import zip_longest
    except ImportError:
        from itertools import izip_longest as zip_longest

    args = [iter(iterable)] * n
    return list(zip_longest(fillvalue=fillvalue, *args))


def sort_and_groupby(items, key=None, reverse=False, ret_lists=False):
    """
    Sort ``items`` using ``key`` function and invoke itertools.groupby to group items.
    If ret_lists is True, a tuple of lists (keys, groups) is returned else iterator.
    See itertools.groupby for further info.

    >>> sort_and_groupby([1, 2, 1], ret_lists=True)
    ([1, 2], [[1, 1], [2]])
    """
    from itertools import groupby
    if not ret_lists:
        return groupby(sorted(items, key=key, reverse=reverse), key=key)
    else:
        keys, groups = [], []
        for hvalue, grp in groupby(sorted(items, key=key, reverse=reverse), key=key):
            keys.append(hvalue)
            groups.append(list(grp))
        return keys, groups


#########################################################################################
# Sorting and ordering
#########################################################################################

def prune_ord(alist):
    """
    Return new list where all duplicated items in alist are removed

    1) The order of items in alist is preserved.
    2) items in alist MUST be hashable.

    Taken from http://code.activestate.com/recipes/52560/
    >>> prune_ord([1, 1, 2, 3, 3])
    [1, 2, 3]
    """
    mset = {}
    return [mset.setdefault(e, e) for e in alist if e not in mset]

#########################################################################################
# Special functions
#########################################################################################

def gaussian(x, width, center=0.0, height=None):
    """
    Returns the values of gaussian(x) where x is array-like.

    Args:
        x: Input array.
        width: Width of the gaussian.
        center: Center of the gaussian.
        height: height of the gaussian. If height is None, a normalized gaussian is returned.
    """
    x = np.asarray(x)
    if height is None: height = 1.0 / (width * np.sqrt(2 * np.pi))

    return height * np.exp(-((x - center) / width) ** 2 / 2.)


def lorentzian(x, width, center=0.0, height=None):
    """
    Returns the values of gaussian(x) where x is array-like.

    Args:
        x: Input array.
        width: Width of the Lorentzian (half-width at half-maximum)
        center: Center of the Lorentzian.
        height: height of the Lorentzian. If height is None, a normalized Lorentzian is returned.
    """
    x = np.asarray(x)
    if height is None: height = 1.0 / (width * np.pi)

    return height * width**2 / ((x - center) ** 2 + width ** 2)

#=====================================
# === Data Interpolation/Smoothing ===
#=====================================

def smooth(x, window_len=11, window='hanning'):
    """
    smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    Taken from http://www.scipy.org/Cookbook/SignalSmooth

    Args:
        x:
            the input signal
        window_len:
            the dimension of the smoothing window. it should be an odd integer
        window:
            the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'.
            'flat' window will produce a moving average smoothing.

    Returns:
        the smoothed signal.

    example::

        t = linspace(-2,2,0.1)
        x = sin(t)+randn(len(t))*0.1
        y = smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    """
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

    if window_len % 2 == 0:
        raise ValueError("window_len should be odd.")

    windows = ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']

    if window not in windows:
        raise ValueError("window must be in: " + str(windows))

    s = np.r_[x[window_len - 1:0:-1], x, x[-1:-window_len:-1]]

    if window == 'flat': # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.' + window + '(window_len)')

    y = np.convolve(w / w.sum(), s, mode='valid')

    s = window_len // 2
    e = s + len(x)
    return y[s:e]


def find_convindex(values, tol, min_numpts=1, mode="abs", vinf=None):
    """
    Given a list of values and a tolerance tol, returns the leftmost index for which

        abs(value[i] - vinf) < tol if mode == "abs"

    or
        abs(value[i] - vinf) / vinf < tol if mode == "rel"

    Args:
        tol: Tolerance
        min_numpts: Minimum number of points that must be converged.
        mode: "abs" for absolute convergence, "rel" for relative convergence.
        vinf: Used to specify an alternative value instead of values[-1].
            By default, vinf = values[-1]

    Return:
        -1 if convergence is not achieved else the index in values.
    """
    vinf = values[-1] if vinf is None else vinf

    if mode == "abs":
        vdiff = [abs(v - vinf) for v in values]
    elif mode == "rel":
        vdiff = [abs(v - vinf) / vinf for v in values]
    else:
        raise ValueError("Wrong mode %s" % mode)

    numpts, i = len(vdiff), -2
    if numpts > min_numpts and vdiff[-2] < tol:
        for i in range(numpts-1, -1, -1):
            if vdiff[i] > tol:
                break
        if (numpts - i - 1) < min_numpts: i = -2

    return i + 1


class BlochRegularGridInterpolator(object):
    """
    This object interpolates the periodic part of a Bloch state in real space.
    """

    def __init__(self, structure, datar, add_replicas=True):
        """
        Args:
            structure: :class:`Structure` object.
            datar: [ndt, nx, ny, nz] array.
            add_replicas: If True, data is padded with redundant data points.
                in order to have a periodic 3D array of shape=[ndt, nx+1, ny+1, nz+1].
        """
        from scipy.interpolate import RegularGridInterpolator
        self.structure = structure

        if add_replicas:
            datar = add_periodic_replicas(datar)

        self.dtype = datar.dtype
        # We want a 4d array (ndt arrays of shape (nx, ny, nz)
        nx, ny, nz = datar.shape[-3:]
        datar = np.reshape(datar, (-1,) + (nx, ny, nz))
        self.ndt = len(datar)
        x = np.linspace(0, 1, num=nx)
        y = np.linspace(0, 1, num=ny)
        z = np.linspace(0, 1, num=nz)

        # Build `ndt` interpolators. Note that RegularGridInterpolator supports
        # [nx, ny, nz, ...] arrays but then each call operates on the full set of
        # ndt components and this complicates the declation of callbacks
        # operating on a single component.
        self._interpolators = [None] * self.ndt
        for i in range(self.ndt):
            self._interpolators[i] = RegularGridInterpolator((x, y, z), datar[i])

    def eval_line(self, point1, point2, num=200, cartesian=False, kpoint=None):
        """
        Interpolate values along a line.

        Args:
            point1: First point of the line. Accepts 3d vector or integer.
                The vector is in reduced coordinates unless `cartesian == True`.
                If integer, the first point of the line is given by the i-th site of the structure
                e.g. `point1=0, point2=1` gives the line passing through the first two atoms.
            point2: Second point of the line. Same API as `point1`.
            num: Number of points sampled along the line.
            cartesian: By default, `point1` and `point1` are interpreted as points in fractional
                coordinates (if not integers). Use True to pass points in cartesian coordinates.
            kpoint: k-point in reduced coordinates. If not None, the phase-factor e^{ikr} is included.

        Return: named tuple with
            site1, site2: None if the points do not represent atomic sites.
            points: Points in fractional coords.
            dist: the distance of points along the line in Ang.
            values: numpy array of shape [ndt, num] with interpolated values.
        """
        site1 = None
        if duck.is_intlike(point1):
            if point1 > len(self.structure):
                raise ValueError("point1: %s > natom: %s" % (point1, len(self.structure)))
            site1 = self.structure[point1]
            point1 = site1.coords if cartesian else site1.frac_coords

        site2 = None
        if duck.is_intlike(point2):
            if point2 > len(self.structure):
                raise ValueError("point2: %s > natom: %s" % (point2, len(self.structure)))
            site2 = self.structure[point2]
            point2 = site2.coords if cartesian else site2.frac_coords

        point1 = np.reshape(point1, (3,))
        point2 = np.reshape(point2, (3,))
        if cartesian:
            red_from_cart = self.structure.lattice.inv_matrix.T
            point1 = np.dot(red_from_cart, point1)
            point2 = np.dot(red_from_cart, point2)

        p21 = point2 - point1
        line_points = np.reshape([alpha * p21 for alpha in np.linspace(0, 1, num=num)], (-1, 3))
        dist = self.structure.lattice.norm(line_points)
        line_points += point1

        return dict2namedtuple(site1=site1, site2=site2, points=line_points, dist=dist,
                               values=self.eval_points(line_points, kpoint=kpoint))

    def eval_points(self, frac_coords, idt=None, cartesian=False, kpoint=None):
        """
        Interpolate values on an arbitrary list of points.

        Args:
            frac_coords: List of points in reduced coordinates unless `cartesian`.
            idt: Index of the sub-array to interpolate. If None, all sub-arrays are interpolated.
            cartesian: True if points are in cartesian coordinates.
            kpoint: k-point in reduced coordinates. If not None, the phase-factor e^{ikr} is included.

        Return:
            [ndt, npoints] array or [1, npoints] if idt is not None
        """
        frac_coords = np.reshape(frac_coords, (-1, 3))
        if cartesian:
            red_from_cart = self.structure.lattice.inv_matrix.T
            frac_coords = [np.dot(red_from, v) for v in frac_coords]

        uc_coords = np.reshape(frac_coords, (-1, 3)) % 1

        if idt is None:
            values = np.empty((self.ndt, len(uc_coords)), dtype=self.dtype)
            for idt in range(self.ndt):
                values[idt] = self._interpolators[idt](uc_coords)
        else:
            values = self._interpolators[idt](uc_coords)

        if kpoint is not None:
            if hasattr(kpoint, "frac_coords"): kpoint = kpoint.frac_coords
            kpoint = np.reshape(kpoint, (3,))
            values *= np.exp(2j * np.pi * np.dot(frac_coords, kpoint))

        return values
