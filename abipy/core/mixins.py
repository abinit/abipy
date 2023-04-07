# coding: utf-8
"""This module provides mixin classes"""
from __future__ import annotations

import abc
import os
import collections
import tempfile
import pickle
import numpy as np
import pandas as pd

from time import ctime
from monty.os.path import which
from monty.termcolor import cprint
from monty.string import list_strings
from monty.collections import dict2namedtuple
from monty.functools import lazy_property


__all__ = [
    "AbinitNcFile",
    "Has_Structure",
    "Has_ElectronBands",
    "Has_PhononBands",
    "NotebookWriter",
    "Has_Header",
    "SlotPickleMixin"
]


class BaseFile(metaclass=abc.ABCMeta):
    """
    Abstract base class defining the methods that must be implemented
    by the concrete classes representing the different files produced by ABINIT.
    """
    def __init__(self, filepath: str):
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
        if isinstance(filepath, cls): return filepath

        #print("Perhaps the subclass", cls, "must redefine the classmethod from_file.")
        return cls(filepath)

    @property
    def filepath(self) -> str:
        """Absolute path of the file."""
        return self._filepath

    @property
    def relpath(self) -> str:
        """Relative path."""
        try:
            return os.path.relpath(self.filepath)
        except OSError:
            # current working directory may not be defined!
            return self.filepath

    @property
    def basename(self) -> str:
        """Basename of the file."""
        return os.path.basename(self.filepath)

    @property
    def filetype(self) -> str:
        """String defining the filetype."""
        return self.__class__.__name__

    def filestat(self, as_string: bool = False) -> dict:
        """
        Dictionary with file metadata, if ``as_string`` is True, a string is returned.
        """
        d = get_filestat(self.filepath)
        if not as_string: return d
        return "\n".join("%s: %s" % (k, v) for k, v in d.items())

    @abc.abstractmethod
    def close(self) -> None:
        """Close the file."""

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Activated at the end of the with statement. It automatically closes the file."""
        self.close()

    def remove(self) -> None:
        """Close the file handle, remove the file from disk."""
        self.close()
        try:
            os.remove(self.filepath)
        except FileNotFoundError:
            pass


class TextFile(BaseFile):

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

    def close(self) -> None:
        """Close the file."""
        try:
            self._file.close()
        except Exception:
            pass

    def seek(self, offset, whence: int = 0) -> None:
        """Set the file's current position, like stdio's fseek()."""
        self._file.seek(offset, whence)

    def get_panel(self, **kwargs):
        import panel as pn
        import panel.widgets as pnw

        root, ext = os.path.splitext(self.basename)
        text = open(self.filepath, "rt").read()
        if not text: text = "This file is empty!"

        # Use Markdown for selected extensions else Ace editor.
        if ext and len(ext) > 1: ext = ext[1:]
        ext2format = dict(sh="shell", py="python", stdin="shell", stdout="shell", stderr="shell")

        if ext in ext2format:
            fmt = ext2format[ext]
            obj = pn.pane.Markdown(f"```{fmt}\n{text}\n```", sizing_mode="stretch_both")
        else:
            obj = pnw.Ace(value=text, language='text', readonly=True,
                          sizing_mode='stretch_width', height=1200)

        return pn.Column(f"## File: {self.filepath}",
                         obj,
                         pn.layout.Divider(),
                         sizing_mode="stretch_width")


class JsonFile(TextFile):
    """
    A TextFile containing JSON data.
    Provides get_panel method so that we can visualize the file with `abiopen.py FILE --panel`
    """

    def get_panel(self, **kwargs):
        import json
        from abipy.panels.viewers import JSONViewer
        with self:
            return JSONViewer(json.load(self._file))


class AbinitNcFile(BaseFile):
    """
    Abstract class representing a Netcdf file with data saved
    according to the ETSF-IO specifications (when available).
    An AbinitNcFile has a netcdf reader to read data from file and build objects.
    """

    @classmethod
    def from_binary_string(cls, bstring) -> AbinitNcFile:
        """
        Build object from a binary string with the netcdf data.
        Useful for implementing GUIs in which widgets returns binary data.
        """
        workdir = tempfile.mkdtemp()
        fd, tmp_path = tempfile.mkstemp(suffix=".nc")
        with open(tmp_path, "wb") as fh:
            fh.write(bstring)
            return cls.from_file(tmp_path)

    def ncdump(self, *nc_args, **nc_kwargs) -> str:
        """Returns a string with the output of ncdump."""
        return NcDumper(*nc_args, **nc_kwargs).dump(self.filepath)

    @lazy_property
    def abinit_version(self) -> str:
        """String with the abinit version: three digits separated by comma."""
        return self.reader.rootgrp.getncattr("abinit_version")

    @abc.abstractproperty
    def params(self) -> dict:
        """
        :class:`OrderedDict` with the convergence parameters
        Used to construct |pandas-DataFrames|.
        """

    def get_dims_dataframe(self, as_dict=False, path="/") -> pd.DataFrame:
        """
        Return: |pandas-Dataframe| with the dimensions defined in the `path` group.
            or dict if as_dict is True.
        """
        grp = self.reader.rootgrp if path == "/" else self.path2group[path]
        d = {k: len(v) for k, v in grp.dimensions.items()}

        if as_dict: return d

        # Since this is a Series but we want a dataframe to facilitate interoperability.
        # we have to call init with additional kwargs.
        return pd.DataFrame.from_dict(d, orient='index', columns=['value'])

    def get_input_string(self) -> str:
        """
        Read and return input string stored in the netcdf.
        Only nc files generared by Abinit9 have this variable.
        """
        if "input_string" in self.reader.rootgrp.variables:
            return self.reader.read_string("input_string")
        else:
            return "Nc file does not contain `input_string`"

    def get_ncfile_view(self, **kwargs):
        """
        Return panel Parameterized object with widgets to visualize
        netcdf dimensions and variables.
        """
        from abipy.panels.core import NcFileViewer
        return NcFileViewer(self).get_ncfile_view(**kwargs)


class AbinitFortranFile(BaseFile):
    """
    Abstract class representing a Fortran file containing output data from abinit.
    """
    def close(self) -> None:
        """nop, just to fulfill the abstract interface."""


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
        super().__init__(filepath)
        self.structure, self.mesh, self.data = cube_read_structure_mesh_data(self.filepath)

    def close(self) -> None:
        """nop, just to fulfill the abstract interface."""

    #@classmethod
    #def write_structure_mesh_data(cls, path, structure, mesh, data):
    #    with open(path, "wt") as fh:
    #        cube_write_structure_mesh(fh, structure, mesh)
    #        cube_write_data(fh, data, mesh):


class Has_Structure(metaclass=abc.ABCMeta):
    """
    Mixin class for |AbinitNcFile| containing crystallographic data.
    """

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

    def export_structure(self, filepath: str):
        """
        Export the structure on file.

        returns: |Visualizer| instance.
        """
        return self.structure.export(filepath)

    def visualize_structure_with(self, appname: str):
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
            iatom_list = np.array(iatom_list, dtype=int)

        # Slice full arrays.
        wyckoffs = ea.wyckoffs[iatom_list]
        wyck_labels = ea.wyck_labels[iatom_list]
        site_labels = ea.site_labels[iatom_list]

        return dict2namedtuple(iatom_list=iatom_list, wyckoffs=wyckoffs, wyck_labels=wyck_labels, site_labels=site_labels)

    def yield_structure_figs(self, **kwargs):
        """*Generates* a predefined list of matplotlib figures with minimal input from the user."""
        yield self.structure.plot(show=False)

    def yield_structure_plotly_figs(self, **kwargs):
        """*Generates* a predefined list of plotly figures with minimal input from the user."""
        yield self.structure.plotly(show=False)


class Has_ElectronBands(metaclass=abc.ABCMeta):
    """Mixin class for |AbinitNcFile| containing electron data."""

    @abc.abstractproperty
    def ebands(self):
        """Returns the |ElectronBands| object."""

    @property
    def nsppol(self) -> int:
        """Number of spin polarizations"""
        return self.ebands.nsppol

    @property
    def nspinor(self) -> int:
        """Number of spinors"""
        return self.ebands.nspinor

    @property
    def nspden(self) -> int:
        """Number of indepedendent spin-density components."""
        return self.ebands.nspden

    @property
    def mband(self) -> int:
        """Maximum number of bands."""
        return self.ebands.mband

    @property
    def nband(self) -> int:
        """Maximum number of bands."""
        return self.ebands.nband

    @property
    def nelect(self) -> float:
        """Number of electrons per unit cell"""
        return self.ebands.nelect

    @property
    def nkpt(self) -> int:
        """Number of k-points."""
        return self.ebands.nkpt

    @property
    def kpoints(self):
        """Iterable with the Kpoints."""
        return self.ebands.kpoints

    @lazy_property
    def tsmear(self):
        return self.ebands.smearing.tsmear_ev.to("Ha")

    def get_ebands_params(self) -> dict:
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
        with_gaps = not self.ebands.is_metal
        if self.ebands.kpoints.is_path:
            yield self.ebands.plot(with_gaps=with_gaps, show=False)
            if len(self.ebands.kpoints) > 1:
                yield self.ebands.kpoints.plot(show=False)
        else:
            edos = self.ebands.get_edos()
            yield self.ebands.plot_with_edos(edos, with_gaps=with_gaps, show=False)
            yield edos.plot(show=False)

    def yield_ebands_plotly_figs(self, **kwargs):
        """*Generates* a predefined list of plotly figures with minimal input from the user."""
        with_gaps = not self.ebands.is_metal

        if self.ebands.kpoints.is_path:
            yield self.ebands.plotly(with_gaps=with_gaps, show=False)
            if len(self.ebands.kpoints) > 1:
                yield self.ebands.kpoints.plotly(show=False)
        else:
            edos = self.ebands.get_edos()
            yield self.ebands.plotly_with_edos(edos, with_gaps=with_gaps, show=False)
            yield edos.plotly(show=False)

    def expose_ebands(self, slide_mode=False, slide_timeout=None, expose_web=False, **kwargs):
        """
        Shows a predefined list of matplotlib figures for electron bands with minimal input from the user.
        """
        from abipy.tools.plotting import MplExpose, PanelExpose

        if expose_web:
            e = PanelExpose(title=f"e-Bands of {self.structure.formula}")
        else:
            e = MplExpose(slide_mode=slide_mode, slide_timeout=slide_mode, verbose=1)

        with e:
            e(self.yield_ebands_figs(**kwargs))

    #def plotly_expose_ebands(self, **kwargs):
    #    """
    #    This function *generates* a predefined list of plotly figures with minimal input from the user.
    #    """
    #    chart_studio = kwargs.pop("chart_studio", None)
    #    verbose = kwargs.pop("verbose", 0)
    #    kwargs.update(dict(
    #        renderer="chart_studio" if chart_studio else None,
    #        title=f"Band structure of {self.ebands.structure.formula}",
    #        with_gaps = not self.ebands.has_metallic_scheme,
    #    ))
    #    self.ebands.plotly(**kwargs)


class Has_PhononBands(metaclass=abc.ABCMeta):
    """
    Mixin class for |AbinitNcFile| containing phonon data.
    """

    @abc.abstractproperty
    def phbands(self):
        """Returns the |PhononBands| object."""

    def get_phbands_params(self) -> dict:
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

    def yield_phbands_plotly_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of plotly figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """
        units = kwargs.get("units", "mev")
        yield self.phbands.qpoints.plotly(show=False)
        yield self.phbands.plotly(units=units, show=False)
        yield self.phbands.plot_colored_matched(units=units, show=False)

    def expose_phbands(self, slide_mode=False, slide_timeout=None, **kwargs):
        """
        Shows a predefined list of matplotlib figures for phonon bands with minimal input from the user.
        """
        from abipy.tools.plotting import MplExpose
        with MplExpose(slide_mode=slide_mode, slide_timeout=slide_mode, verbose=1) as e:
            e(self.yield_phbands_figs(**kwargs))


class NcDumper(object):
    """
    Wrapper object for the ncdump tool.
    """

    def __init__(self, *nc_args, **nc_kwargs):
        """
        Args:
            nc_args: Arguments passed to ncdump.
            nc_kwargs: Keyword arguments passed to ncdump
        """
        self.nc_args = nc_args
        self.nc_kwargs = nc_kwargs
        self.ncdump = which("ncdump")

    def dump(self, filepath: str) -> str:
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


class HasNotebookTools(object):

    def has_panel(self):
        """
        Return panel module (that evaluates to True) if panel is installed else False.
        """
        try:
            import panel as pn
            return pn
        except ImportError:
            return False

    def make_and_open_notebook(self, nbpath=None, foreground=False,
                               classic_notebook=False, no_browser=False):  # pragma: no cover
        """
        Generate an jupyter_ notebook and open it in the browser.

        Args:
            nbpath: If nbpath is None, a temporay file is created.
                foreground: By default, jupyter is executed in background and stdout, stderr are redirected.
                to devnull. Use foreground to run the process in foreground
            classic_notebook: True to use the classic notebook instead of jupyter-lab (default)
            no_browser: Start the jupyter server to serve the notebook but don't open the notebook in the browser.
                        Use this option to connect remotely from localhost to the machine running the kernel

        Return: system exit code.

        Raise: `RuntimeError` if jupyter executable is not in $PATH
        """
        nbpath = self.write_notebook(nbpath=nbpath)

        if not classic_notebook:
            # Use jupyter-lab.
            app_path = which("jupyter-lab")
            if app_path is None:
                raise RuntimeError("""
Cannot find jupyter-lab application in $PATH. Install it with:

    conda install -c conda-forge jupyterlab

or:

    pip install jupyterlab

See also https://jupyterlab.readthedocs.io/
""")

        else:
            # Use classic notebook
            app_path = which("jupyter")
            if app_path is None:
                raise RuntimeError("""
Cannot find jupyter application in $PATH. Install it with:

    conda install -c conda-forge jupyter

or:

    pip install jupyterlab

See also https://jupyter.readthedocs.io/en/latest/install.html
""")
            app_path = app_path + " notebook "

        if not no_browser:
            if foreground:
                return os.system("%s %s" % (app_path, nbpath))
            else:
                fd, tmpname = tempfile.mkstemp(text=True)
                print(tmpname)
                cmd = "%s %s" % (app_path, nbpath)
                print("Executing:", cmd, "\nstdout and stderr redirected to %s" % tmpname)
                import subprocess
                process = subprocess.Popen(cmd.split(), shell=False, stdout=fd, stderr=fd)
                cprint("pid: %s" % str(process.pid), "yellow")
                return 0

        else:
            # Based on https://github.com/arose/nglview/blob/master/nglview/scripts/nglview.py
            notebook_name = os.path.basename(nbpath)
            dirname = os.path.dirname(nbpath)
            print("nbpath:", nbpath)

            import socket

            def find_free_port():
                """https://stackoverflow.com/questions/1365265/on-localhost-how-do-i-pick-a-free-port-number"""
                from contextlib import closing
                with closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as s:
                    s.bind(('', 0))
                    s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
                    return s.getsockname()[1]

            username = os.getlogin()
            hostname = socket.gethostname()
            port = find_free_port()

            client_cmd = "ssh -NL localhost:{port}:localhost:{port} {username}@{hostname}".format(
                username=username, hostname=hostname, port=port)

            print(f"""
Using port: {port}

\033[32m In your local machine, run: \033[0m

                {client_cmd}

\033[32m NOTE: you might want to replace {hostname} by full hostname with domain name \033[0m
\033[32m Then open your web browser, copy and paste the URL: \033[0m

http://localhost:{port}/notebooks/{notebook_name}
""")
            if not classic_notebook:
                cmd = f'{app_path} {notebook_name} --no-browser --port {port} --notebook-dir {dirname}'
            else:
                cmd = f'{app_path} notebook {notebook_name} --no-browser --port {port} --notebook-dir {dirname}'

            print("Executing:", cmd)
            print('NOTE: make sure to open `{}` in your local machine\n'.format(notebook_name))

            return os.system(cmd)

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
import sys, os
import numpy as np

%matplotlib notebook

# Use this magic for jupyterlab.
# For installation instructions, see https://github.com/matplotlib/jupyter-matplotlib
#%matplotlib widget

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
""")
        ])

        return nbformat, nbv, nb

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


class NotebookWriter(HasNotebookTools, metaclass=abc.ABCMeta):
    """
    Mixin class for objects that are able to generate jupyter_ notebooks.
    Subclasses must provide a concrete implementation of `write_notebook`.
    """

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

    @classmethod
    def pickle_load(cls, filepath: str):
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

        Return: The name of the pickle file.
        """
        if filepath is None:
            _, filepath = tempfile.mkstemp(suffix='.pickle')

        with open(filepath, "wb") as fh:
            pickle.dump(self, fh)
            return filepath

    @abc.abstractmethod
    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """

    #@abc.abstractmethod
    #def yield_plotly_figs(self, **kwargs):  # pragma: no cover
    #    """
    #    This function *generates* a predefined list of matplotlib figures with minimal input from the user.
    #    Used in abiview.py to get a quick look at the results.
    #    """

    def _get_panel_and_template(self):
        # Create panel template with matplotlib figures and show them in the browser.
        import panel as pn
        pn.config.sizing_mode = 'stretch_width'
        from abipy.panels.core import get_template_cls_from_name
        cls = get_template_cls_from_name("FastGridTemplate")

        title = self.__class__.__name__
        if hasattr(self, "structure"): title = f"{title} <small>({self.structure.formula})</small>"
        template = cls(
            title=title,
            header_background="#ff8c00 ", # Dark orange
        )

        return pn, template

    def expose(self, slide_mode=False, slide_timeout=None, use_web=False, **kwargs):
        """
        Shows a predefined list of matplotlib figures with minimal input from the user.
        Relies on the ``yield_fig``s methods implemented by the subclass to generate matplotlib figures.

        Args:
            use_web: True to show all figures inside a panel template executed in the local browser.
               False to show figures in different GUIs
        """
        if not use_web:
            # Produce all matplotlib figures and show them with the X-server.
            from abipy.tools.plotting import MplExpose
            with MplExpose(slide_mode=slide_mode, slide_timeout=slide_mode, verbose=1) as e:
                e(self.yield_figs(**kwargs))

        else:
            # Create panel template with matplotlib figures and show them in the browser.
            pn, template = self._get_panel_and_template()
            pn.config.sizing_mode = 'stretch_width'
            from abipy.panels.core import mpl
            for i, fig in enumerate(self.yield_figs()):
                row, col = divmod(i, 2)
                p = mpl(fig, with_divider=False, dpi=82)
                if hasattr(template.main, "append"):
                    template.main.append(p)
                else:
                    # Assume .main area acts like a GridSpec
                    row_slice = slice(3 * row, 3 * (row + 1))
                    if col == 0: template.main[row_slice, :6] = p
                    if col == 1: template.main[row_slice, 6:] = p

            return template.show()

    def plotly_expose(self, **kwargs):
        """
        This function *generates* a predefined list of plotly figures with minimal input from the user.
        Relies on the yield_plotly_figs method implemented by the subclass in order to generate the figures.
        """
        #print("in plotly expose")
        pn, template = self._get_panel_and_template()
        pn.config.sizing_mode = 'stretch_width'
        from abipy.panels.core import mpl, ply

        # Insert figure in template.main.
        from abipy.tools.plotting import is_mpl_figure, is_plotly_figure
        for i, fig in enumerate(self.yield_plotly_figs()):
            row, col = divmod(i, 2)
            # Handle both matplotlib and plotly figures since we dont' support plotly everywhere.
            if is_plotly_figure(fig):
                p = ply(fig, with_divider=False)
            elif is_mpl_figure(fig):
                p = mpl(fig, with_divider=False)
            else:
                raise TypeError(f"Don't know how to handle type: `{type(fig)}`")

            if hasattr(template.main, "append"):
                template.main.append(p)
            else:
                # Assume main area acts like a panel GridSpec
                row_slice = slice(3 * row, 3 * (row + 1))
                if col == 0: template.main[row_slice, :6] = p
                if col == 1: template.main[row_slice, 6:] = p

        return template.show()


class Has_Header(object):
    """Mixin class for netcdf files containing the Abinit header."""

    @lazy_property
    def hdr(self):
        """|AttrDict| with the Abinit header e.g. hdr.ecut."""
        return self.reader.read_abinit_hdr()

    #def compare_hdr(self, other_hdr):


class SlotPickleMixin:
    """
    This mixin makes it possible to pickle/unpickle objects with __slots__
    defined.
    """

    def __getstate__(self):
        return {slot: getattr(self, slot) for slot in self.__slots__ if hasattr(self, slot)}

    def __setstate__(self, state):
        for slot, value in state.items():
            setattr(self, slot, value)
