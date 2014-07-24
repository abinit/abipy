"""
This module defines the object Function1D that described a functions.
of a single variables and provides simple interfaces for performing
common tasks such as algebraic operations, integrations, differentiations, plots ...
"""
from __future__ import print_function, division

import itertools
import cStringIO as StringIO
import numpy as np

from abipy.tools.derivatives import finite_diff

__all__ = [
    "Function1D",
]


class Function1D(object):
    """A (real|complex) function of a real variable."""
    def __init__(self, mesh, values):
        """
        Args:
            mesh:
                array-like object with the real points of the grid.
            values:
                array-like object with the values of the function (can be complex-valued)
        """
        self._mesh = np.ascontiguousarray(mesh)
        self._values = np.ascontiguousarray(values)
        assert len(self.mesh) == len(self.values)

    @property
    def mesh(self):
        """Array with the mesh points"""
        return self._mesh

    @property
    def values(self):
        """Values of the functions."""
        return self._values

    def __len__(self):
        return len(self.mesh)

    def __iter__(self):
        return itertools.izip(self.mesh, self.values)

    def __getitem__(self, slice):
        return self.mesh[slice], self.values[slice]

    def __eq__(self, other):
        if other is None: return False
        return (self.has_same_mesh(other) and 
                np.allclose(self.values, other.values))

    def __ne__(self, other):
        return not (self == other)

    def __neg__(self):
        return self.__class__(self.mesh, -self.values)

    def __pos__(self):
        return self

    def __abs__(self):
        return self.__class__(self.mesh, np.abs(self.values))

    def __add__(self, other):
        cls = self.__class__
        if isinstance(other, cls):
            assert self.has_same_mesh(other)
            return cls(self.mesh, self.values+other.values)
        else:
            return cls(self.mesh, self.values+np.array(other))
    __radd__ = __add__

    def __sub__(self, other):
        cls = self.__class__
        if isinstance(other, cls):
            assert self.has_same_mesh(other)
            return cls(self.mesh, self.values-other.values)
        else:
            return cls(self.mesh, self.values-np.array(other))

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        cls = self.__class__
        if isinstance(other, cls):
            assert self.has_same_mesh(other)
            return cls(self.mesh, self.values*other.values)
        else:
            return cls(self.mesh, self.values*other)
    __rmul__ = __mul__

    def __truediv__(self, other):
        cls = self.__class__
        if isinstance(other, cls):
            assert self.has_same_mesh(other)
            return cls(self.mesh, self.values/other.values)
        else:
            return cls(self.mesh, self.values/other)
    __rtruediv__ = __truediv__

    def __pow__(self, other):
        return self.__class__(self.mesh, self.values**other)

    @property
    def real(self):
        """Return new `Function1D` with the real part of self."""
        return self.__class__(self.mesh, self.values.real)

    @property
    def imag(self):
        """Return new `Function1D` with the imaginary part of self."""
        return self.__class__(self.mesh, self.values.real)

    def conjugate(self):
        """Return new `Function1D` with the complex conjugate."""
        return self.__class__(self.mesh, self.values.conjugate)

    @classmethod
    def from_func(cls, func, mesh):
        """
        Initialize the object from a callable.

        .. example:

            Function1D.from_func(np.sin, mesh=np.arange(1,10))
        """
        return cls(mesh, np.vectorize(func)(mesh))

    @classmethod
    def from_file(cls, path, comments="#", delimiter=None, usecols=(0, 1)):
        """
        Initialize an instance by reading data from path (txt format)
        see also :func:`np.loadtxt`

        Args:
            path:
                Path to the file containing data.
            comments:
                The character used to indicate the start of a comment; default: '#'.
            delimiter: str, optional
                The string used to separate values. By default, this is any whitespace.
            usecols:
                sequence, optional
                Which columns to read, with 0 being the first.
                For example, usecols = (1,4) will extract data from the 2nd, and 5th columns.
        """
        mesh, values = np.loadtxt(path, comments=comments, delimiter=delimiter,
                                  usecols=usecols, unpack=True)
        return cls(mesh, values)

    def to_file(self, path, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# '):
        """
        Save self to a text file. See :func:`np.savetext` for the description of the variables
        """
        data = zip(self.mesh, self.values)
        np.savetxt(path, data, fmt=fmt, delimiter=delimiter, newline=newline,
                   header=header, footer=footer, comments=comments)

    def __repr__(self):
        return "%s at %s, size = %d" % (self.__class__.__name__, id(self), len(self))

    def __str__(self):
        stream = StringIO.StringIO()
        self.to_file(stream)
        return "\n".join(stream.getvalue())

    def has_same_mesh(self, other):
        """True if self and other have the same mesh."""
        return np.allclose(self.mesh, other.mesh)

    def bma(self):
        """Return b-a. f(x) is defined in [a,b]"""
        return self.mesh[-1] - self.mesh[0]

    @property
    def max(self):
        """Max of f(x) if f is real, max of :math:`|f(x)|` if complex."""
        if not self.iscomplexobj:
            return self.values.max()
        else:
            # Max of |z|
            return np.max(np.abs(self.values))

    @property
    def min(self):
        """Min of f(x) if f is real, min of :math:`|f(x)|` if complex."""
        if not self.iscomplexobj:
            return self.values.min()
        else:
            return np.max(np.abs(self.values))

    @property
    def iscomplexobj(self):
        """
        Check if values is array of complex numbers.
        The type of the input is checked, not the value. Even if the input
        has an imaginary part equal to zero, `iscomplexobj` evaluates to True.
        """
        return np.iscomplexobj(self.values)

    @property
    def h(self):
        """The spacing of the mesh. None if mesh is not homogeneous."""
        try:
            return self._h
        except AttributeError:
            return self.dx[0] if np.allclose(self.dx[0], self.dx) else None

    @property
    def dx(self):
        """
        ndarray of len(self)-1 elements giving the distance between two
        consecutive points of the mesh, i.e. dx[i] = ||x[i+1] - x[i]||.
        """
        try:
            return self._dx
        except AttributeError:
            self._dx = dx = np.zeros(len(self)-1)
            for (i, x) in enumerate(self.mesh[:-1]):
                dx[i] = self.mesh[i+1] - x
            return self._dx

    def find_mesh_index(self, value):
        """
        Return the index of the first point in the mesh whose value is >= value
        -1 if not found
        """
        for (i, x) in enumerate(self.mesh):
            if x >= value:
                return i
        return -1

    def finite_diff(self, order=1, acc=4):
        """
        Compute the derivatives by finite differences.

        Args:
            order:
                Order of the derivative.
            acc:
                Accuracy. 4 is fine in many cases.

        Returns:
            new `Function1d` instance with the derivative.
        """
        if self.h is None:
            raise ValueError("Finite differences with inhomogeneous meshes are not supported")

        return self.__class__(self.mesh, finite_diff(self.values, self.h, order=order, acc=acc))

    def integral(self):
        """
        Cumulatively integrate y(x) using the composite trapezoidal rule.

        Returns:
            `Function1d` with :math:`\int y(x) dx`
        """
        from scipy.integrate import cumtrapz
        integ = cumtrapz(self.values, x=self.mesh)
        pad_intg = np.zeros(len(self.values))
        pad_intg[1:] = integ

        return self.__class__(self.mesh, pad_intg)

    @property
    def spline(self):
        """Cubic spline with s=0"""
        try:
            return self._spline
        except AttributeError:
            from scipy.interpolate import UnivariateSpline
            self._spline = UnivariateSpline(self.mesh, self.values, s=0)
            return self._spline

    @property
    def spline_roots(self):
        """Zeros of the spline."""
        return self.spline.roots()

    def spline_derivatives(self, x):
        """Returns all derivatives of the spline at point x."""
        return self.spline.derivatives(x)

    def spline_integral(self, a=None, b=None):
        """
        Returns the definite integral of the spline of f between two given points a and b

        Args:
            a:
                First point. mesh[0] if a is None
            b:
                Last point. mesh[-1] if a is None
        """
        a = self.mesh[0] if a is None else a
        b = self.mesh[-1] if b is None else b
        return self.spline.integral(a, b)

    @property
    def l1_norm(self):
        """Compute :math:`\int |f(x)| dx`."""
        return abs(self).integral()[-1][1]

    @property
    def l2_norm(self):
        """Compute :math:`\int |f(x)|^2 dx`."""
        return (abs(self)**2).integral()[-1][1]

    def fft(self):
        """Compute the FFT transform (negative sign)."""
        # Compute FFT and frequencies.
        from scipy import fftpack
        n, d = len(self), self.h
        fft_vals = fftpack.fft(self.values, n=n)
        freqs = fftpack.fftfreq(n, d=d)

        # Shift the zero-frequency component to the center of the spectrum.
        fft_vals = fftpack.fftshift(fft_vals)
        freqs = fftpack.fftshift(freqs)

        return self.__class__(freqs, fft_vals)

    def ifft(self, x0=None):
        """Compute the FFT transform :math:`\int e+i`"""
        # Rearrange values in the standard order then perform IFFT.
        from scipy import fftpack
        n, d = len(self), self.h
        fft_vals = fftpack.ifftshift(self.values)
        fft_vals = fftpack.ifft(fft_vals, n=n)

        # Compute the mesh of the IFFT output.
        x0 = 0.0 if x0 is None else x0
        stop = 1./d - 1./(n*d) + x0
        mesh = np.linspace(x0, stop, num=n)

        return self.__class__(mesh, fft_vals)

    #def convolve(self, other):
    #    ""Convolution with FFT."""
    #    assert self.has_same_mesh(other)
    #    from scipy.signal import convolve
    #    conv = convolve(self.values, other.values, mode="same") * self.h
    #    return self.__class__(self.mesh, conv)

    #def gauss_convolve(self, width):
    #   """Convolve self with a gaussian of standard deviation width."""
    #    from abipy.tools import gauss_ufunc
    #    gauss = gauss_ufunc(width, center=0.0, height=None)
    #    gdata = Function1D(self.mesh, gauss(self.mesh))
    #    return self.convolve(gdata)

    #def lorentz_convolve(self, gamma):
    #   """Convolve self with a Lorentzian."""
    #    lorentz = lorentz_ufunc(gamma, center=0.0, height=None)
    #    return self.convolve(Function1D.from_func(self.mesh, lorentz))

    #def smooth(self, window_len=21, window="hanning"):
    #    from abipy.tools import smooth
    #    smooth_vals = smooth(self.values, window_len=window_len, window=window)
    #    return self.__class__(self.mesh, smooth_vals)

    def plot_ax(self, ax, exchange_xy=False, *args, **kwargs):
        """
        Helper function to plot self on axis ax.

        Args:
            ax:
                `matplotlib` axis.
            exchange_xy:
                True to exchange the axis in the plot
            args:
                Positional arguments passed to ax.plot
            kwargs:
                Keyword arguments passed to `matplotlib`.
                Accepts also:

                cplx_mode:
                    string defining the data to print. Possible choices are (case-insensitive):

                        - "re"  for real part.
                        - "im" for imaginary part.
                        - "abs" for the absolute value

                    Options can be concatenated with "-"

        Returns:
            List of lines added.
        """
        xx, yy = self.mesh, self.values
        if exchange_xy:
            xx, yy = yy, xx

        cplx_mode = kwargs.pop("cplx_mode", "re-im").lower().split("-")

        if np.iscomplexobj(yy):
            lines = []
            if "re" in cplx_mode:
                lines.extend(ax.plot(xx, yy.real, *args, **kwargs))

            if "im" in cplx_mode:
                lines.extend(ax.plot(xx, yy.imag, *args, **kwargs))

            if "abs" in cplx_mode:
                lines.extend(ax.plot(xx, np.abs(yy), *args, **kwargs))

        else:
            lines = ax.plot(xx, yy, *args, **kwargs)

        return lines

    def plot(self, **kwargs):
        """
        Args:
            args:
                Positional arguments passed to `matplotlib`.

        ==============  ==============================================================
        kwargs          Meaning
        ==============  ==============================================================
        title           Title of the plot (Default: None).
        show            True to show the figure (Default: True).
        savefig:        'abc.png' or 'abc.eps'* to save the figure to a file.
        ==============  ==============================================================

        Returns:
            `matplotlib` figure.
        """
        title = kwargs.pop("title", None)
        show = kwargs.pop("show", True)
        savefig = kwargs.pop("savefig", None)

        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        ax.grid(True)
        if title:
            ax.set_title(title)

        exchange_xy = kwargs.pop("exchange_xy", False)
        self.plot_ax(ax, exchange_xy=exchange_xy, **kwargs)

        if show:
            plt.show()

        if savefig is not None:
            fig.savefig(savefig)

        return fig
