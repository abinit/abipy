# coding: utf-8
"""Numeric tools."""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
import bisect as bs

from abipy.tools import duck

#########################################################################################
# Array tools
#########################################################################################


def transpose_last3dims(arr):
    """Transpose the last three dimensions of arr: (...,x,y,z) --> (...,z,y,x)."""
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

#########################################################################################
# Bisection algorithms
#########################################################################################


def minloc(iterable):
    """Return the min value and its position."""
    min_val, min_idx = iterable[0], 0

    for idx, item in enumerate(iterable[1:]):
        if item < min_val:
            min_val, min_idx = item, idx

    return min_val, min_idx


def maxloc(iterable):
    """Return the max value and its position."""
    max_val, max_idx = iterable[0], 0

    for idx, item in enumerate(iterable[1:]):
        if item > max_val:
            max_val, max_idx = item, idx

    return max_val, max_idx

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


#########################################################################################
# Sorting and ordering
#########################################################################################

def prune_ord(alist):
    """
    Return new list where all duplicated items in alist are removed

    1) The order of items in alist is preserved.
    2) items in alist MUST be hashable.

    Taken from http://code.activestate.com/recipes/52560/
    >>> prune_ord([1,1,2,3,3])
    [1, 2, 3]
    """
    mset = {}
    return [mset.setdefault(e, e) for e in alist if e not in mset]

#########################################################################################
# Special functions
#########################################################################################


def gauss_ufunc(width, center=0.0, height=None):
    """
    Returns a gaussian (u)function with the given parameters.

    If height is None, a normalized gaussian is returned.
    """
    if height is None:
        height = 1.0 / (width * np.sqrt(2 * np.pi))

    return lambda x: height * np.exp(-((x - center) / width) ** 2 / 2.)


def gaussian(x, width, center=0.0, height=None):
    """
    Returns the values of gaussian(x) where x is array-like.

    If height is None, a normalized gaussian is returned.
    """
    x = np.asarray(x)

    if height is None:
        height = 1.0 / (width * np.sqrt(2 * np.pi))

    return height * np.exp(-((x - center) / width) ** 2 / 2.)

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
    #return y
    return y[s:e]


def integrator_linspace(y, dx=1, method="simps", axis=-1):
    """
    Integrate y using samples along the given axis with spacing dx.
    method is in ['simps', 'trapz', 'romb']. Note that romb requires len(y) = 2**k + 1.

    Example:
    >>> nx = 2**8+1
    >>> x = np.linspace(0, np.pi, nx)
    >>> dx = np.pi/(nx-1)
    >>> y = np.sin(x)
    >>> for method in ["simps", "trapz", "romb"]: print(integrator_linspace(y, dx=dx, method=method))
    2.00000000025
    1.99997490024
    2.0
    """
    from scipy.integrate import simps, trapz, romb

    if method == "simps":
        return simps(y, x=None, dx=dx, axis=axis, even='avg')
    elif method == "trapz":
        return trapz(y, x=None, dx=dx, axis=axis)
    elif method == "romb":
        return romb(y, dx=dx, axis=axis, show=False)
    else:
        raise ValueError("Wrong method: %s" % method)


class Interpol3D(object):
    """Container to store and apply interpolating schemes."""

    _allowed_mesh_types = [
        'cubic',
        'tetrahedral',
    ]

    _allowed_interpol_schemes = dict({
        'cubic': [
            'trilin',
            'tricub',
            'tripen',
            'tricos',
        ],
        'tetrahedral': [
        ],
    })

    def __init__(self, mesh_type=None, interpol_scheme=None):
        if mesh_type:
            self.set_mesh_type(mesh_type)
            if interpol_scheme:
                self.set_interpol_scheme(interpol_scheme)
        else:
            self.mesh_type = None
            self.interpol_scheme = None

    def set_mesh_type(self, mesh_type):
        if not mesh_type in self._allowed_mesh_types:
            raise ValueError("Wrong mesh type: " + str(mesh_type))
        self.mesh_type = mesh_type

    def set_interpol_scheme(self, interpol_scheme):
        if not self.mesh_type:
            raise ValueError("Mesh type not specified")
        if not interpol_scheme in self._allowed_interpol_schemes[self.mesh_type]:
            raise ValueError(
                "Wrong interpolation scheme \"%s\" for mesh type \"%s\"" % (str(interpol_scheme), str(self.mesh_type)))
        self.interpol_scheme = interpol_scheme

    def set_values(self, values):
        if not self.mesh_type:
            raise ValueError("Mesh type not specified")
        if not self.interpol_scheme:
            raise ValueError("Interpolation scheme not specified")
        if self.mesh_type == 'cubic':
            if len(values) != 8:
                raise ValueError("Wrong length of array \"values\" in set_values for mesh type %s" % self.mesh_type)
        elif self.mesh_type == 'tetrahedral':
            raise NotImplementedError("")
            if len(values) != 4:
                raise ValueError("Wrong length of array \"values\" in set_values for mesh type %s" % self.mesh_type)
        self.values = values
        self.initialized = True

    def interp(self, x, y, z):
        if not self.initialized:
            raise ValueError("Interpolation not fully initialized")
        if self.mesh_type == 'cubic':
            if self.interpol_scheme == 'trilin':
                return self.interp_cubic_trilin(x, y, z)
            elif self.interpol_scheme == 'tricub':
                return self.interp_cubic_tricub(x, y, z)
            elif self.interpol_scheme == 'tripen':
                return self.interp_cubic_tripen(x, y, z)
            elif self.interpol_scheme == 'tricos':
                return self.interp_cubic_tricos(x, y, z)

    def interp_cubic_trilin(self, x, y, z):
        fx00 = self.values[0] + x * (self.values[1] - self.values[0])
        fx10 = self.values[2] + x * (self.values[3] - self.values[2])
        fx01 = self.values[4] + x * (self.values[5] - self.values[4])
        fx11 = self.values[6] + x * (self.values[7] - self.values[6])
        fxy0 = fx00 + y * (fx10 - fx00)
        fxy1 = fx01 + y * (fx11 - fx01)
        return fxy0 + z * (fxy1 - fxy0)

    def interp_cubic_tricub(self, x, y, z):
        px = -2 * x * x * x + 3 * x * x
        fx00 = self.values[0] + px * (self.values[1] - self.values[0])
        fx10 = self.values[2] + px * (self.values[3] - self.values[2])
        fx01 = self.values[4] + px * (self.values[5] - self.values[4])
        fx11 = self.values[6] + px * (self.values[7] - self.values[6])
        py = -2 * y * y * y + 3 * y * y
        fxy0 = fx00 + py * (fx10 - fx00)
        fxy1 = fx01 + py * (fx11 - fx01)
        return fxy0 + (-2 * z * z * z + 3 * z * z) * (fxy1 - fxy0)

    def interp_cubic_tripen(self, x, y, z, k=None):
        if k is None:
            a = -4
            b = 10
            c = -10
            d = 5
        else:
            a = -(24 + 16 * k)
            b = 60 + 40 * k
            c = -50 - 32 * k
            d = 15 + 8 * k
        px = a * x * x * x * x * x + b * x * x * x * x + c * x * x * x + d * x * x
        fx00 = self.values[0] + px * (self.values[1] - self.values[0])
        fx10 = self.values[2] + px * (self.values[3] - self.values[2])
        fx01 = self.values[4] + px * (self.values[5] - self.values[4])
        fx11 = self.values[6] + px * (self.values[7] - self.values[6])
        py = a * y * y * y * y * y + b * y * y * y * y + c * y * y * y + d * y * y
        fxy0 = fx00 + py * (fx10 - fx00)
        fxy1 = fx01 + py * (fx11 - fx01)
        return fxy0 + (a * z * z * z * z * z + b * z * z * z * z + c * z * z * z + d * z * z) * (fxy1 - fxy0)

    def interp_cubic_tricos(self, x, y, z):
        px = 0.5 - 0.5 * np.cos(np.pi * x)
        fx00 = self.values[0] + px * (self.values[1] - self.values[0])
        fx10 = self.values[2] + px * (self.values[3] - self.values[2])
        fx01 = self.values[4] + px * (self.values[5] - self.values[4])
        fx11 = self.values[6] + px * (self.values[7] - self.values[6])
        py = 0.5 - 0.5 * np.cos(np.pi * y)
        fxy0 = fx00 + py * (fx10 - fx00)
        fxy1 = fx01 + py * (fx11 - fx01)
        return fxy0 + (0.5 - 0.5 * np.cos(np.pi * z)) * (fxy1 - fxy0)


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

    def __init__(self, structure, datar, kpoint=None, add_replicas=True):
        """
        Args:
            structure: :class:`Structure` object.
            datar: [ndt, nx, ny, nz] array
            kpoint: k-point in reduced coordinates. If not None, the phase-factor e^{ikr} is included.
            add_replicas: If True, data is padded with redundant data points.
                in order to have a periodic 3D array of shape=[ndt, nx+1, ny+1, nz+1].
        """
        from scipy.interpolate import RegularGridInterpolator
        self.structure = structure

        self.kpoint = None
        if kpoint is not None:
            if hasattr(kpoint, "frac_coords"): kpoint = kpoint.frac_coords
            self.kpoint = np.reshape(kpoint, (3,))

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

    def eval_line(self, point1, point2, num=200, cartesian=False):
        """
        Interpolate values along a line.
        """
        if duck.is_intlike(point1):
            site1 = self.structure[point1]
            point1 = site1.coords if cartesian else site1.frac_coords

        if duck.is_intlike(point2):
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

        return dist, self.eval_points(line_points)

    def eval_points(self, frac_coords, idt=None, cartesian=False):
        """
        Interpolate values on an arbitrary list of points.
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

        if self.kpoint is not None:
            values *= np.exp(2j * np.pi * np.dot(frac_coords, self.kpoint))

        return values
