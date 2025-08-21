# coding: utf-8
"""
AnaddbNcFile provides a high-level interface to the data stored in the anaddb.nc file.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from collections import OrderedDict
from functools import cached_property
from monty.string import marquee
from monty.termcolor import cprint
from abipy.core.mixins import AbinitNcFile, Has_Structure, NotebookWriter
from abipy.abio.robots import Robot
from abipy.iotools import ETSF_Reader
from abipy.tools.typing import Figure
from abipy.tools.plotting import add_fig_kwargs, get_axarray_fig_plt, rotate_ticklabels
from abipy.tools.tensors import Tensor, DielectricTensor, NLOpticalSusceptibilityTensor
from abipy.dfpt.ifc import InteratomicForceConstants
from abipy.dfpt.ddb import Zeffs, DynQuad
from abipy.dfpt.elastic import ElasticData


class AnaddbNcFile(AbinitNcFile, Has_Structure, NotebookWriter):
    """
    AnaddbNcFile provides a high-level interface to the data stored in the anaddb.nc file.
    This object is usually instantiated with `abiopen("anaddb.nc")`.

    .. attribute:: structure

        |Structure| object.

    .. attribute:: epsinf

        Macroscopic dielectric tensor. None if the file does not contain this information.

    .. attribute:: becs

        Born effective charges. None if the file does not contain this inf

    .. attribute:: ifc

        :class:`InteratomicForceConstants` object with the interatomic force constants calculated by anaddb.
        None, if the netcdf file does not contain the IFCs.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: AnaddbNcFile
    """

    @classmethod
    def from_file(cls, filepath: str) -> AnaddbNcFile:
        """Initialize the object from file."""
        return cls(filepath)

    def __init__(self, filepath: str):
        super().__init__(filepath)
        self.r = ETSF_Reader(filepath)

    def close(self):
        self.r.close()

    @cached_property
    def structure(self):
        return self.r.read_structure()

    @cached_property
    def params(self):
        # -666 to support old anaddb.nc files without metadata
        return OrderedDict([
            ("asr", int(self.r.read_value("asr", default=-666))),
            ("chneut", int(self.r.read_value("chneut", default=-666))),
            ("dipdip", int(self.r.read_value("dipdip", default=-666))),
            ("symdynmat", int(self.r.read_value("symdynmat", default=-666))),
        ])

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        """
        String representation

        Args:
            verbose: verbosity level.
        """
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))

        import json
        app(marquee("Parameters", mark="="))
        app(json.dumps(self.params, indent=2, sort_keys=True))
        app("")

        if self.has_elastic_data:
            app("")
            df = self.elastic_data.get_average_elastic_dataframe(tensor="elastic_relaxed")
            if not df.empty:
                app(marquee("Averaged elastic properties (relaxed ions)", mark="="))
                app(df.to_string(index=False))
                app("")
            df = self.elastic_data.get_elast_properties_dataframe()
            if not df.empty:
                app(marquee("Averaged elastic properties (relaxed ions)", mark="="))
                app(df.T.to_string(index=True))

            if verbose:
                df = self.elastic_data.get_voigt_dataframe()
                app(df.T.to_string())

        tol = 1e-3
        if self.epsinf is not None:
            app("Electronic dielectric tensor (eps_inf) in Cartesian coordinates. Set to zero below %.2e." % tol)
            app(self.epsinf.get_dataframe(tol=tol).to_string())
            app("")

        if self.eps0 is not None:
            app("Zero-frequency dielectric tensor (eps_zero) in Cartesian coordinates. Set to zero below %.2e." % tol)
            app(self.eps0.get_dataframe(tol=tol).to_string())
            app("")

        if self.dchide is not None:
            app("Non-linear optical susceptibility tensor.")
            app(str(self.dchide))
            app("")

        if self.dchidt is not None:
            app("First-order change in the linear dielectric susceptibility.")
            app(str(self.dchidt))
            app("")

        #if self.has_piezoelectric_data:
        #    df = self.elastic_data.get_piezoelectric_dataframe()

        return "\n".join(lines)

    @cached_property
    def epsinf(self) -> DielectricTensor:
        """
        Macroscopic electronic |DielectricTensor| in Cartesian coordinates (a.k.a. epsilon_infinity)
        None if the file does not contain this information.
        """
        try:
            return DielectricTensor(self.r.read_value("emacro_cart").T.copy())
        except Exception as exc:
            #print(exc, "Returning None", sep="\n")
            return None

    @cached_property
    def eps0(self) -> DielectricTensor:
        """
        Relaxed ion macroscopic |DielectricTensor| in Cartesian coordinates (a.k.a. epsilon_zero)
        None if the file does not contain this information.
        """
        try:
            return DielectricTensor(self.r.read_value("emacro_cart_rlx").T.copy())
        except Exception as exc:
            #print(exc, "Requires dieflag > 0", "Returning None", sep="\n")
            return None

    @cached_property
    def becs(self) -> Zeffs | None:
        """
        Born effective charges. None if the file does not contain this information.
        """
        try:
            return Zeffs("Ze", self.r.read_value("becs_cart"), self.structure, self.params)
        except Exception as exc:
            return None

    @cached_property
    def dyn_quad(self) -> DynQuad | None:
        try:
            # First two dimension are associated to the q-point, then atomic perturbation in Cart coords.
            # nctkarr_t('quadrupoles_cart', "dp", 'three, three, three, number_of_atoms')
            quad_cart = self.r.read_value("quadrupoles_cart")
            quad_cart = np.swapaxes(quad_cart, -1, -2)
            return DynQuad(quad_cart, self.structure, self.params)
        except Exception as exc:
            #print(exc)
            return None

    @cached_property
    def ifc(self) -> InteratomicForceConstants:
        """
        The interatomic force constants calculated by anaddb.
        The following anaddb variables should be used in the run: ``ifcflag``, ``natifc``, ``atifc``, ``ifcout``.
        Return None, if the netcdf file does not contain the IFCs,
        """
        try:
            return InteratomicForceConstants.from_file(self.filepath)
        except Exception as exc:
            cprint("Interatomic force constants have not been calculated. Returning None", "red")
            return None

    @cached_property
    def dchide(self) -> NLOpticalSusceptibilityTensor:
        """
        Non-linear optical susceptibility tensor.
        Returns a :class:`NLOpticalSusceptibilityTensor` or None if the file does not contain this information.
        """
        try:
            return NLOpticalSusceptibilityTensor(self.r.read_value("dchide"))
        except Exception as exc:
            #print(exc, "Requires nlflag > 0", "Returning None", sep="\n")
            return None

    @cached_property
    def dchidt(self):
        """
        First-order change in the linear dielectric susceptibility.
        Returns a list of lists of 3x3 Tensor object with shape (number of atoms, 3).
        The [i][j] element of the list contains the Tensor representing the change due to the
        displacement of the ith atom in the jth direction.
        None if the file does not contain this information.
        """
        try:
            a = self.r.read_value("dchidt").T.copy()
        except Exception as exc:
            #print(exc, "Requires 0 < nlflag < 3", "Returning None", sep="\n")
            return None

        dchidt = []
        for i in a:
            d = []
            for j in i:
                d.append(Tensor(j))
            dchidt.append(d)

        return dchidt

    @cached_property
    def oscillator_strength(self) -> np.ndarray:
        """
        A complex |numpy-array| containing the oscillator strengths with shape [number of phonon modes, 3, 3],
        in a.u. (1 a.u.=253.2638413 m3/s2).
        None if the file does not contain this information.
        """
        try:
            carr = self.r.read_value("oscillator_strength", cmode="c")
            carr = carr.transpose((0, 2, 1)).copy()
            return carr
        except Exception as exc:
            #print(exc, "Oscillator strengths require dieflag == 1, 3 or 4", "Returning None", sep="\n")
            return None

    @cached_property
    def has_elastic_data(self) -> bool:
        """True if elastic tensors have been computed."""
        return self.r.read_value("elaflag", default=0) != 0

    @cached_property
    def has_piezoelectric_data(self) -> bool:
        """True if piezoelectric tensors have been computed."""
        return self.r.read_value("piezoflag", default=0) != 0

    @cached_property
    def elastic_data(self) -> ElasticData:
        """
        Container with the different (piezo)elastic tensors computed by anaddb.
        stored in pymatgen tensor objects.
        """
        return ElasticData.from_ncreader(self.r)

    @cached_property
    def amu(self) -> dict:
        """
        Dictionary with atomic_number as keys and the atomic massu units as values.
        """
        amu_list = self.r.read_value("atomic_mass_units", default=None)
        if amu_list is not None:
            atomic_numbers = self.r.read_value("atomic_numbers")
            amu = {at: a for at, a in zip(atomic_numbers, amu_list)}
        else:
            amu = None

        return amu

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """
        yield None

    def write_notebook(self, nbpath=None) -> str:
        """
        Write an jupyter_ notebook to nbpath. If ``nbpath`` is None, a temporary file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("ananc = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(ananc)"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class AnaddbNcRobot(Robot):
    """
    This robot analyzes the results contained in multiple anaddb.nc files.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: AnaddbNcRobot
    """
    EXT = "anaddb"

    @property
    def has_elastic_data(self) -> bool:
        return all(ncfile.has_elastic_data for ncfile in self.abifiles)

    def get_dataframe(self):
        if self.has_elastic_data:
            return self.get_elastic_dataframe()
        return None

    def get_elastic_dataframe(self, with_geo=True, abspath=False, with_params=False,
                              funcs=None, **kwargs) -> pd.DataFrame:
        """
        Return a |pandas-DataFrame| with properties derived from the elastic tensor
        and an associated structure. Filename is used as index.

        Args:
            with_geo: True if structure info should be added to the dataframe
            abspath: True if paths in index should be absolute. Default: Relative to getcwd().
            with_params: False to exclude calculation parameters from the dataframe.

        kwargs:
            attrs:
                List of additional attributes of the |GsrFile| to add to the DataFrame.
            funcs: Function or list of functions to execute to add more data to the DataFrame.
                Each function receives a |GsrFile| object and returns a tuple (key, value)
                where key is a string with the name of column and value is the value to be inserted.
        """
        # Add attributes specified by the users
        attrs = [
            #"energy", "pressure", "max_force",
            #"nsppol", "nspinor", "nspden",
        ] + kwargs.pop("attrs", [])

        rows, index = [], []
        for label, ncfile in self.items():
            index.append(label)
            d = {}

            # Add info on structure.
            if with_geo:
                d.update(ncfile.structure.get_dict4pandas(with_spglib=True))

            if with_params:
                d.update(self.params)

            # Execute functions
            if funcs is not None: d.update(self._exec_funcs(funcs, ncfile))

            df = ncfile.elastic_data.get_elast_properties_dataframe(etypes="elastic_relaxed")
            d.update(df.to_dict("records")[0])

            rows.append(d)

        return pd.DataFrame(rows, index=index, columns=list(rows[0].keys() if rows else None))

    @add_fig_kwargs
    def plot_elastic_properties(self, fontsize=8, **kwargs) -> Figure:
        """
        Args:
            fontsize: legend and label fontsize.

        Returns: |matplotlib-Figure|
        """
        df = self.get_elastic_dataframe(with_geo=False, abspath=False, with_params=False)
        from pandas.api.types import is_numeric_dtype
        keys = [k for k in df.keys() if is_numeric_dtype(df[k])]
        i = keys.index("fitted_to_structure")
        if i != -1:
            keys.pop(i)

        num_plots, ncols, nrows = len(keys), 1, 1
        if num_plots > 1:
            ncols = 3
            nrows = (num_plots // ncols) + (num_plots % ncols)

        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=False, sharey=False, squeeze=False)
        ax_list = ax_list.ravel()

        for ix, (key, ax) in enumerate(zip(keys, ax_list)):
            irow, icol = divmod(ix, ncols)
            xn = range(len(df.index))
            ax.plot(xn, df[key], marker="o")
            ax.grid(True)
            ax.set_xticks(xn)
            ax.set_ylabel(key, fontsize=fontsize)
            ax.set_xticklabels([])

        ax.set_xticklabels(self.keys(), fontsize=fontsize)
        rotate_ticklabels(ax, 15)

        if ix != len(ax_list) - 1:
            for ix in range(ix + 1, len(ax_list)):
                ax_list[ix].axis('off')

        return fig

    #def get_voigt_dataframe(self, tensor_names):
    #    ncfile.get_voigt_dataframe(self, tensor_names):

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """
        if self.has_elastic_data:
            yield self.plot_elastic_properties(show=False)
        else:
            yield None

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporary file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        args = [(l, f.filepath) for l, f in self.items()]
        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("robot = abilab.AnaddbNcRobot(*%s)\nrobot.trim_paths()\nrobot" % str(args)),
            #nbv.new_code_cell("df = ebands_plotter.get_ebands_frame()\ndisplay(df)"),
        ])

        # Mixins
        nb.cells.extend(self.get_baserobot_code_cells())

        return self._write_nb_nbpath(nb, nbpath)
