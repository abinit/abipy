# coding: utf-8
"""This module ..."""
from __future__ import print_function, division, unicode_literals, absolute_import

import abc
import os
import six
import collections
import tempfile
import pickle

from time import ctime
import numpy as np
from monty.os.path import which
from monty.termcolor import cprint
from monty.string import is_string, list_strings
from monty.collections import dict2namedtuple
from monty.functools import lazy_property
from abipy.flowtk.netcdf import NetcdfReader, NO_DEFAULT


__all__ = [
    "AbinitNcFile",
    "Has_Structure",
    "Has_ElectronBands",
    "Has_PhononBands",
    "NotebookWriter",
    "Has_Header",
]

@six.add_metaclass(abc.ABCMeta)
class BaseFile(object):
    """
    Abstract base class defining the methods that must be implemented
    by the concrete classes representing the different files produced by ABINIT.
    """
    def __init__(self, filepath):
        self._filepath = os.path.abspath(filepath)

        # Save stat values
        stat = os.stat(filepath)
        self._last_atime = stat.st_atime
        self._last_mtime = stat.st_mtime
        self._last_ctime = stat.st_ctime

    def __repr__(self):
        return "<%s, %s>" % (self.__class__.__name__, self.relpath)

    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a string."""
        if isinstance(filepath, cls):
            return filepath

        #print("Perhaps the subclass", cls, "must redefine the classmethod from_file.")
        return cls(filepath)

    @property
    def filepath(self):
        """Absolute path of the file."""
        return self._filepath

    @property
    def relpath(self):
        """Relative path."""
        try:
            return os.path.relpath(self.filepath)
        except OSError:
            # current working directory may not be defined!
            return self.filepath

    @property
    def basename(self):
        """Basename of the file."""
        return os.path.basename(self.filepath)

    @property
    def filetype(self):
        """String defining the filetype."""
        return self.__class__.__name__

    def filestat(self, as_string=False):
        """
        Dictionary with file metadata
        if ``as_string`` is True, a string is returned.
        """
        d = get_filestat(self.filepath)
        if not as_string: return d
        return "\n".join("%s: %s" % (k, v) for k, v in d.items())

    @abc.abstractmethod
    def close(self):
        """Close file."""

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Activated at the end of the with statement. It automatically closes the file."""
        self.close()

    #def __del__(self):
    #    """
    #    Called when the instance is about to be destroyed.
    #    """
    #    try:
    #        self.close()
    #    finally:
    #        super(BaseFile, self).__close__(self)


class TextFile(BaseFile):

    #@classmethood
    #def from_string(cls, s):
    #    return cls.from_file(filepath)

    def __enter__(self):
        # Open the file
        self._file
        return self

    def __iter__(self):
        return iter(self._file)

    @lazy_property
    def _file(self):
        """File object open in read-only mode."""
        return open(self.filepath, mode="rt")

    def close(self):
        """Close the file."""
        try:
            self._file.close()
        except Exception:
            pass

    def seek(self, offset, whence=0):
        """Set the file's current position, like stdio's fseek()."""
        self._file.seek(offset, whence)


@six.add_metaclass(abc.ABCMeta)
class AbinitNcFile(BaseFile):
    """
    Abstract class representing a Netcdf file with data saved
    according to the ETSF-IO specifications (when available).
    An AbinitNcFile has a netcdf reader to read data from file and build objects.
    """
    def ncdump(self, *nc_args, **nc_kwargs):
        """Returns a string with the output of ncdump."""
        return NcDumper(*nc_args, **nc_kwargs).dump(self.filepath)

    @lazy_property
    def abinit_version(self):
        """String with abinit version: three digits separated by comma."""
        return self.reader.rootgrp.getncattr("abinit_version")

    @abc.abstractproperty
    def params(self):
        """
        :class:`OrderedDict` with the convergence parameters
        Used to construct |pandas-DataFrames|.
        """


@six.add_metaclass(abc.ABCMeta)
class AbinitFortranFile(BaseFile):
    """
    Abstract class representing a fortran file containing output data from abinit.
    """
    def close(self):
        pass


class CubeFile(BaseFile):
    """

    .. attribute:: structure

        |Structure| object

    .. attribute:: mesh

        |Mesh3d| object with information on the uniform 3d mesh.

    .. attribute:: data

        |numpy-array| of shape [nx, ny, nz] with numerical values on the real-space mesh.
    """
    def __init__(self, filepath):
        from abipy.iotools.cube import cube_read_structure_mesh_data
        super(CubeFile, self).__init__(filepath)
        self.structure, self.mesh, self.data = cube_read_structure_mesh_data(self.filepath)

    def close(self):
        """nop, just to fulfill the abstract interface."""

    #@classmethod
    #def write_structure_mesh_data(cls, path, structure, mesh, data):
    #    with open(path, "wt") as fh:
    #        cube_write_structure_mesh(fh, structure, mesh)
    #        cube_write_data(fh, data, mesh):


@six.add_metaclass(abc.ABCMeta)
class Has_Structure(object):
    """Mixin class for :class:`AbinitNcFile` containing crystallographic data."""

    @abc.abstractproperty
    def structure(self):
        """Returns the |Structure| object."""

    def plot_bz(self, **kwargs):
        """
        Gives the plot (as a matplotlib object) of the symmetry line path in the Brillouin Zone.
        """
        return self.structure.plot_bz(**kwargs)

    # To maintain backward compatbility
    show_bz = plot_bz

    def export_structure(self, filepath):
        """
        Export the structure on file.

        returns: |Visualizer| instance.
        """
        return self.structure.export(filepath)

    def visualize_structure_with(self, appname):
        """
        Visualize the crystalline structure with the specified visualizer.

        See |Visualizer| for the list of applications and formats supported.
        """
        from abipy.iotools.visualizer import Visualizer
        visu = Visualizer.from_name(appname)

        for ext in visu.supported_extensions():
            ext = "." + ext
            try:
                return self.export_structure(ext)
            except visu.Error:
                pass
        else:
            raise visu.Error("Don't know how to export data for appname %s" % appname)

    def _get_atomview(self, view, select_symbols=None, verbose=0):
        """
        Helper function used to select (inequivalent||all) atoms depending on view.
        Uses spglib to find inequivalent sites.

        Args:
            view: "inequivalent" to show only inequivalent atoms. "all" for all sites.
            select_symbols: String or list of strings with chemical symbols.
                Used to select only atoms of this type.

        Return named tuple with:

                * iatom_list: list of site index.
                * wyckoffs: Wyckoff letters
                * site_labels: Labels for each site in `iatom_list` e.g Si2a
        """
        natom = len(self.structure)
        if natom == 1: verbose = False
        if verbose:
            print("Calling spglib to find inequivalent sites. Magnetic symmetries (if any) are not taken into account.")

        ea = self.structure.spget_equivalent_atoms(printout=verbose > 0)

        # Define iatom_list depending on view
        if view == "all":
            iatom_list = np.arange(natom)
        elif view == "inequivalent":
            iatom_list = ea.irred_pos
        else:
            raise ValueError("Wrong value for view: %s" % str(view))

        # Filter by element symbol.
        if select_symbols is not None:
            select_symbols = set(list_strings(select_symbols))
            iatom_list = [i for i in iatom_list if self.structure[i].specie.symbol in select_symbols]
            iatom_list = np.array(iatom_list, dtype=np.int)

        # Slice full arrays.
        wyckoffs = ea.wyckoffs[iatom_list]
        wyck_labels = ea.wyck_labels[iatom_list]
        site_labels = ea.site_labels[iatom_list]

        return dict2namedtuple(iatom_list=iatom_list, wyckoffs=wyckoffs, wyck_labels=wyck_labels, site_labels=site_labels)

    def yield_structure_figs(self, **kwargs):
        """*Generates* a predefined list of matplotlib figures with minimal input from the user."""
        yield self.structure.plot(show=False)


@six.add_metaclass(abc.ABCMeta)
class Has_ElectronBands(object):
    """Mixin class for :class:`AbinitNcFile` containing electron data."""

    @abc.abstractproperty
    def ebands(self):
        """Returns the |ElectronBands| object."""

    @property
    def nsppol(self):
        """Number of spin polarizations"""
        return self.ebands.nsppol

    @property
    def nspinor(self):
        """Number of spinors"""
        return self.ebands.nspinor

    @property
    def nspden(self):
        """Number of indepedendent spin-density components."""
        return self.ebands.nspden

    @property
    def mband(self):
        """Maximum number of bands."""
        return self.ebands.mband

    @property
    def nband(self):
        """Maximum number of bands."""
        return self.ebands.nband

    @property
    def nelect(self):
        """Number of electrons per unit cell"""
        return self.ebands.nelect

    @property
    def nkpt(self):
        """Number of k-points."""
        return self.ebands.nkpt

    @property
    def kpoints(self):
        """Iterable with the Kpoints."""
        return self.ebands.kpoints

    @lazy_property
    def tsmear(self):
        return self.ebands.smearing.tsmear_ev.to("Ha")

    def get_ebands_params(self):
        """:class:`OrderedDict` with the convergence parameters."""
        return collections.OrderedDict([
            ("nsppol", self.nsppol),
            ("nspinor", self.nspinor),
            ("nspden", self.nspden),
            ("nband", self.nband),
            ("nkpt", self.nkpt),
        ])

    def plot_ebands(self, **kwargs):
        """Plot the electron energy bands. See the :func:`ElectronBands.plot` for the signature."""
        return self.ebands.plot(**kwargs)

    def plot_ebands_with_edos(self, edos, **kwargs):
        """Plot the electron energy bands with DOS. See the :func:`ElectronBands.plot_with_edos` for the signature."""
        return self.ebands.plot_with_edos(edos, **kwargs)

    def get_edos(self, **kwargs):
        """Compute the electronic DOS on a linear mesh. Wraps ebands.get_edos."""
        return self.ebands.get_edos(**kwargs)

    def yield_ebands_figs(self, **kwargs):
        """*Generates* a predefined list of matplotlib figures with minimal input from the user."""
        with_gaps = not self.ebands.has_metallic_scheme
        if self.ebands.kpoints.is_path:
            yield self.ebands.plot(with_gaps=with_gaps, show=False)
            yield self.ebands.kpoints.plot(show=False)
        else:
            edos = self.ebands.get_edos()
            yield self.ebands.plot_with_edos(edos, with_gaps=with_gaps, show=False)
            yield edos.plot(show=False)

    def expose_ebands(self, slide_mode=False, slide_timeout=None, **kwargs):
        """
        Shows a predefined list of matplotlib figures for electron bands with minimal input from the user.
        """
        from abipy.tools.plotting import MplExpose
        with MplExpose(slide_mode=slide_mode, slide_timeout=slide_mode, verbose=1) as e:
            e(self.yield_ebands_figs(**kwargs))


@six.add_metaclass(abc.ABCMeta)
class Has_PhononBands(object):
    """
    Mixin class for :class:`AbinitNcFile` containing phonon data.
    """

    @abc.abstractproperty
    def phbands(self):
        """Returns the |PhononBands| object."""

    def get_phbands_params(self):
        """:class:`OrderedDict` with the convergence parameters."""
        return collections.OrderedDict([
            ("nqpt", len(self.phbands.qpoints)),
        ])

    def plot_phbands(self, **kwargs):
        """
        Plot the electron energy bands. See the :func:`PhononBands.plot` for the signature.""
        """
        return self.phbands.plot(**kwargs)

    #def plot_phbands_with_phdos(self, phdos, **kwargs):
    #    return self.phbands.plot_with_phdos(phdos, **kwargs)

    def yield_phbands_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """
        units = kwargs.get("units", "mev")
        yield self.phbands.qpoints.plot(show=False)
        yield self.phbands.plot(units=units, show=False)
        yield self.phbands.plot_colored_matched(units=units, show=False)

    def expose_phbands(self, slide_mode=False, slide_timeout=None, **kwargs):
        """
        Shows a predefined list of matplotlib figures for phonon bands with minimal input from the user.
        """
        from abipy.tools.plotting import MplExpose
        with MplExpose(slide_mode=slide_mode, slide_timeout=slide_mode, verbose=1) as e:
            e(self.yield_phbands_figs(**kwargs))


class NcDumper(object):
    """Wrapper object for the ncdump tool."""

    def __init__(self, *nc_args, **nc_kwargs):
        """
        Args:
            nc_args: Arguments passed to ncdump.
            nc_kwargs: Keyword arguments passed to ncdump
        """
        self.nc_args = nc_args
        self.nc_kwargs = nc_kwargs
        self.ncdump = which("ncdump")

    def dump(self, filepath):
        """Returns a string with the output of ncdump."""
        if self.ncdump is None:
            return "Cannot find ncdump tool in $PATH"
        else:
            from subprocess import check_output
            return check_output(["ncdump", filepath])


_ABBREVS = [
    (1 << 50, 'Pb'),
    (1 << 40, 'Tb'),
    (1 << 30, 'Gb'),
    (1 << 20, 'Mb'),
    (1 << 10, 'kb'),
    (1, 'b'),
]


def size2str(size):
    """Convert size to string with units."""
    for factor, suffix in _ABBREVS:
        if size > factor:
            break
    return "%.2f " % (size / factor) + suffix


def get_filestat(filepath):
    stat = os.stat(filepath)
    return collections.OrderedDict([
        ("Name", os.path.basename(filepath)),
        ("Directory", os.path.dirname(filepath)),
        ("Size", size2str(stat.st_size)),
        ("Access Time", ctime(stat.st_atime)),
        ("Modification Time", ctime(stat.st_mtime)),
        ("Change Time", ctime(stat.st_ctime)),
    ])


@six.add_metaclass(abc.ABCMeta)
class NotebookWriter(object):
    """
    Mixin class for objects that are able to generate jupyter_ notebooks.
    Subclasses must provide a concrete implementation of `write_notebook`.
    """
    def make_and_open_notebook(self, nbpath=None, foreground=False):  # pragma: no cover
        """
        Generate an jupyter_ notebook and open it in the browser.

        Args:
            nbpath: If nbpath is None, a temporay file is created.
            foreground: By default, jupyter is executed in background and stdout, stderr are redirected
            to devnull. Use foreground to run the process in foreground

        Return:
            system exit code.

        Raise:
            `RuntimeError` if jupyter_ is not in $PATH
        """
        nbpath = self.write_notebook(nbpath=nbpath)

        if which("jupyter") is None:
            raise RuntimeError("Cannot find jupyter in $PATH. Install it with `conda install jupyter or `pip install jupyter`")

        # Use jupyter-lab instead of classic notebook if possible.
        has_jupyterlab = which("jupyter-lab") is not None
        #has_jupyterlab = False
        appname = "jupyter-lab" if has_jupyterlab else "jupyter notebook"

        if foreground:
            return os.system("%s %s" % (appname, nbpath))
        else:
            fd, tmpname = tempfile.mkstemp(text=True)
            print(tmpname)
            cmd = "%s %s" % (appname, nbpath)
            print("Executing:", cmd)
            print("stdout and stderr redirected to %s" % tmpname)
            import subprocess
            process = subprocess.Popen(cmd.split(), shell=False, stdout=fd, stderr=fd)
            cprint("pid: %s" % str(process.pid), "yellow")
            return 0

    @staticmethod
    def get_nbformat_nbv():
        """Return nbformat module, notebook version module"""
        import nbformat
        nbv = nbformat.v4
        return nbformat, nbv

    def get_nbformat_nbv_nb(self, title=None):
        """
        Return ``nbformat`` module, notebook version module
        and new notebook with title and import section
        """
        nbformat, nbv = self.get_nbformat_nbv()
        nb = nbv.new_notebook()

        if title is not None:
            nb.cells.append(nbv.new_markdown_cell("## %s" % title))

        nb.cells.extend([
            nbv.new_code_cell("""\
from __future__ import print_function, division, unicode_literals, absolute_import

import sys, os
import numpy as np

%matplotlib notebook
from IPython.display import display

# This to render pandas DataFrames with https://github.com/quantopian/qgrid
#import qgrid
#qgrid.nbinstall(overwrite=True)  # copies javascript dependencies to your /nbextensions folder

# This to view Mayavi visualizations. See http://docs.enthought.com/mayavi/mayavi/tips.html
#from mayavi import mlab; mlab.init_notebook(backend='x3d', width=None, height=None, local=True)

from abipy import abilab

# Tell AbiPy we are inside a notebook and use seaborn settings for plots.
# See https://seaborn.pydata.org/generated/seaborn.set.html#seaborn.set
abilab.enable_notebook(with_seaborn=True)

# AbiPy widgets for pandas and seaborn plot APIs
#import abipy.display.seabornw import snw
#import abipy.display.pandasw import pdw""")
        ])

        return nbformat, nbv, nb

    @abc.abstractmethod
    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to nbpath. If nbpath is None, a temporay file is created.
        Return path to the notebook. A typical template:

        .. code-block:: python

            # Preable.
            nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

            #####################
            # Put your code here
            nb.cells.extend([
                nbv.new_markdown_cell("# This is a markdown cell"),
                nbv.new_code_cell("a = 1"),
            ])
            #####################

            # Call _write_nb_nbpath
            return self._write_nb_nbpath(nb, nbpath)
        """

    @staticmethod
    def _write_nb_nbpath(nb, nbpath):
        """
        This method must be called at the end of ``write_notebook``.
        nb is the jupyter notebook and nbpath the argument passed to ``write_notebook``.
        """
        import io, os, tempfile
        if nbpath is None:
            _, nbpath = tempfile.mkstemp(prefix="abinb_", suffix='.ipynb', dir=os.getcwd(), text=True)

        # Write notebook
        import nbformat
        with io.open(nbpath, 'wt', encoding="utf8") as fh:
            nbformat.write(nb, fh)
            return nbpath

    @classmethod
    def pickle_load(cls, filepath):
        """
        Loads the object from a pickle file.
        """
        with open(filepath, "rb") as fh:
            new = pickle.load(fh)
            #assert cls is new.__class__
            return new

    def pickle_dump(self, filepath=None):
        """
        Save the status of the object in pickle format.
        If filepath is None, a temporary file is created.

        Return:
            name of the pickle file.
        """
        if filepath is None:
            _, filepath = tempfile.mkstemp(suffix='.pickle')

        with open(filepath, "wb") as fh:
            pickle.dump(self, fh)
            return filepath

    # TODO: Activate this
    @abc.abstractmethod
    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """

    def expose(self, slide_mode=False, slide_timeout=None, **kwargs):
        """
        Shows a predefined list of matplotlib figures with minimal input from the user.
        """
        from abipy.tools.plotting import MplExpose
        with MplExpose(slide_mode=slide_mode, slide_timeout=slide_mode, verbose=1) as e:
            e(self.yield_figs(**kwargs))


class Has_Header(object):
    """Mixin class for netcdf_ files containing the Abinit header."""

    @lazy_property
    def hdr(self):
        """|AttrDict| with the Abinit header e.g. hdr.ecut."""
        return self.reader.read_abinit_hdr()

    #def get_hdr_params(self):
    #    """:class:`OrderedDict` with the convergence parameters."""
    #    return collections.OrderedDict([

    #def compare_hdr(self, other_hdr):
