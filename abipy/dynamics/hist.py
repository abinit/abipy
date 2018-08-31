# coding: utf-8
"""History file with structural relaxation results."""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import numpy as np
import pymatgen.core.units as units

from collections import OrderedDict
from monty.functools import lazy_property
from monty.collections import AttrDict
from monty.string import marquee, list_strings
from pymatgen.core.periodic_table import Element
from pymatgen.analysis.structure_analyzer import RelaxationAnalyzer
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, set_visible
from abipy.core.structure import Structure
from abipy.core.mixins import AbinitNcFile, NotebookWriter
from abipy.abio.robots import Robot
from abipy.iotools import ETSF_Reader
import abipy.core.abinit_units as abu


class HistFile(AbinitNcFile, NotebookWriter):
    """
    File with the history of a structural relaxation or molecular dynamics calculation.

    Usage example:

    .. code-block:: python

        with HistFile("foo_HIST") as hist:
            hist.plot()


    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: HistFile
    """
    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a netcdf_ file"""
        return cls(filepath)

    def __init__(self, filepath):
        super(HistFile, self).__init__(filepath)
        self.reader = HistReader(filepath)

    def close(self):
        """Close the file."""
        self.reader.close()

    @lazy_property
    def params(self):
        """:class:`OrderedDict` with parameters that might be subject to convergence studies."""
        return {}

    def __str__(self):
        return self.to_string()

    # TODO: Add more metadata.
    #@lazy_property
    #def nsppol(self):
    #    """Number of independent spins."""
    #    return self.reader.read_dimvalue("nsppol")

    #@lazy_property
    #def nspden(self):
    #    """Number of independent spin densities."""
    #    return self.reader.read_dimvalue("nspden")


    #@lazy_property
    #def nspinor(self):
    #    """Number of spinor components."""
    #    return self.reader.read_dimvalue("nspinor")

    @lazy_property
    def final_energy(self):
        """Total energy in eV of the last iteration."""
        return self.etotals[-1]

    @lazy_property
    def final_pressure(self):
        """Final pressure in Gpa."""
        cart_stress_tensors, pressures = self.reader.read_cart_stress_tensors()
        return pressures[-1]

    #@lazy_property
    #def final_max_force(self):

    def get_fstats_dict(self, step):
        """
        Return |AttrDict| with stats on the forces at the given ``step``.
        """
        # [time, natom, 3]
        var = self.reader.read_variable("fcart")
        forces = units.ArrayWithUnit(var[step], "Ha bohr^-1").to("eV ang^-1")
        fmods = np.array([np.linalg.norm(force) for force in forces])

        return AttrDict(
            fmin=fmods.min(),
            fmax=fmods.max(),
            fmean=fmods.mean(),
            fstd=fmods.std(),
            drift=np.linalg.norm(forces.sum(axis=0)),
        )

    def to_string(self, verbose=0, title=None):
        """String representation."""
        lines = []; app = lines.append
        if title is not None: app(marquee(title, mark="="))

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.initial_structure.to_string(verbose=verbose, title="Initial Structure"))
        app("")
        app("Number of relaxation steps performed: %d" % self.num_steps)
        app(self.final_structure.to_string(verbose=verbose, title="Final structure"))
        app("")

        an = self.get_relaxation_analyzer()
        app("Volume change in percentage: %.2f%%" % (an.get_percentage_volume_change() * 100))
        d = an.get_percentage_lattice_parameter_changes()
        vals = tuple(d[k] * 100 for k in ("a", "b", "c"))
        app("Percentage lattice parameter changes:\n\ta: %.2f%%, b: %.2f%%, c: %2.f%%" % vals)
        #an.get_percentage_bond_dist_changes(max_radius=3.0)
        app("")

        cart_stress_tensors, pressures = self.reader.read_cart_stress_tensors()
        app("Stress tensor (Cartesian coordinates in GPa):\n%s" % cart_stress_tensors[-1])
        app("Pressure: %.3f [GPa]" % pressures[-1])

        return "\n".join(lines)

    @property
    def num_steps(self):
        """Number of iterations performed."""
        return self.reader.num_steps

    @lazy_property
    def steps(self):
        """Step indices."""
        return list(range(self.num_steps))

    @property
    def initial_structure(self):
        """The initial |Structure|."""
        return self.structures[0]

    @property
    def final_structure(self):
        """The |Structure| of the last iteration."""
        return self.structures[-1]

    @lazy_property
    def structures(self):
        """List of |Structure| objects at the different steps."""
        return self.reader.read_all_structures()

    @lazy_property
    def etotals(self):
        """|numpy-array| with total energies in eV at the different steps."""
        return self.reader.read_eterms().etotals

    def get_relaxation_analyzer(self):
        """
        Return a pymatgen :class:`RelaxationAnalyzer` object to analyze the relaxation in a calculation.
        """
        return RelaxationAnalyzer(self.initial_structure, self.final_structure)

    def to_xdatcar(self, filepath=None, groupby_type=True, to_unit_cell=False, **kwargs):
        """
        Return Xdatcar pymatgen object. See write_xdatcar for the meaning of arguments.

        Args:
            to_unit_cell (bool): Whether to translate sites into the unit cell.
            kwargs: keywords arguments passed to Xdatcar constructor.
        """
        filepath = self.write_xdatcar(filepath=filepath, groupby_type=groupby_type,
                                      to_unit_cell=to_unit_cell, overwrite=True)
        from pymatgen.io.vasp.outputs import Xdatcar
        return Xdatcar(filepath, **kwargs)

    def write_xdatcar(self, filepath="XDATCAR", groupby_type=True, overwrite=False, to_unit_cell=False):
        """
        Write Xdatcar file with unit cell and atomic positions to file ``filepath``.

        Args:
            filepath: Xdatcar filename. If None, a temporary file is created.
            groupby_type: If True, atoms are grouped by type. Note that this option
                may change the order of the atoms. This option is needed because
                there are post-processing tools (e.g. ovito) that do not work as expected
                if the atoms in the structure are not grouped by type.
            overwrite: raise RuntimeError, if False and filepath exists.
            to_unit_cell (bool): Whether to translate sites into the unit cell.

        Return:
            path to Xdatcar file.
        """
        if filepath is not None and os.path.exists(filepath) and not overwrite:
            raise RuntimeError("Cannot overwrite pre-existing file `%s`" % filepath)
        if filepath is None:
            import tempfile
            fd, filepath = tempfile.mkstemp(text=True, suffix="_XDATCAR")

        # int typat[natom], double znucl[npsp]
        # NB: typat is double in the HIST.nc file
        typat = self.reader.read_value("typat").astype(int)
        znucl = self.reader.read_value("znucl")
        ntypat = self.reader.read_dimvalue("ntypat")
        num_pseudos = self.reader.read_dimvalue("npsp")
        if num_pseudos != ntypat:
            raise NotImplementedError("Alchemical mixing is not supported, num_pseudos != ntypat")
        #print("znucl:", znucl, "\ntypat:", typat)

        symb2pos = OrderedDict()
        symbols_atom = []
        for iatom, itype in enumerate(typat):
            itype = itype - 1
            symbol = Element.from_Z(int(znucl[itype])).symbol
            if symbol not in symb2pos: symb2pos[symbol] = []
            symb2pos[symbol].append(iatom)
            symbols_atom.append(symbol)

        if not groupby_type:
            group_ids = np.arange(self.reader.natom)
        else:
            group_ids = []
            for pos_list in symb2pos.values():
                group_ids.extend(pos_list)
            group_ids = np.array(group_ids, dtype=np.int)

        comment = " %s\n" % self.initial_structure.formula
        with open(filepath, "wt") as fh:
            # comment line  + scaling factor set to 1.0
            fh.write(comment)
            fh.write("1.0\n")
            for vec in self.initial_structure.lattice.matrix:
                fh.write("%.12f %.12f %.12f\n" % (vec[0], vec[1], vec[2]))
            if not groupby_type:
                fh.write(" ".join(symbols_atom) + "\n")
                fh.write("1 " * len(symbols_atom) + "\n")
            else:
                fh.write(" ".join(symb2pos.keys()) + "\n")
                fh.write(" ".join(str(len(p)) for p in symb2pos.values()) + "\n")

            # Write atomic positions in reduced coordinates.
            xred_list = self.reader.read_value("xred")
            if to_unit_cell:
                xred_list = xred_list % 1

            for step in range(self.num_steps):
                fh.write("Direct configuration= %d\n" % (step + 1))
                frac_coords = xred_list[step, group_ids]
                for fs in frac_coords:
                    fh.write("%.12f %.12f %.12f\n" % (fs[0], fs[1], fs[2]))

        return filepath

    def visualize(self, appname="ovito", to_unit_cell=False):  # pragma: no cover
        """
        Visualize the crystalline structure with visualizer.
        See :class:`Visualizer` for the list of applications and formats supported.

        Args:
            to_unit_cell (bool): Whether to translate sites into the unit cell.
        """
        if appname == "mayavi": return self.mayaview()

        # Get the Visualizer subclass from the string.
        from abipy.iotools import Visualizer
        visu = Visualizer.from_name(appname)
        if visu.name != "ovito":
            raise NotImplementedError("visualizer: %s" % visu.name)

        filepath = self.write_xdatcar(filepath=None, groupby_type=True, to_unit_cell=to_unit_cell)

        return visu(filepath)()
        #if options.trajectories:
        #    hist.mvplot_trajectories()
        #else:
        #    hist.mvanimate()

    def plot_ax(self, ax, what, fontsize=12, **kwargs):
        """
        Helper function to plot quantity ``what`` on axis ``ax``.

        Args:
            fontsize: fontsize for legend
            kwargs are passed to matplotlib plot method
        """
        label = None
        if what == "energy":
            # Total energy in eV.
            marker = kwargs.pop("marker", "o")
            label = kwargs.pop("label", "Energy")
            ax.plot(self.steps, self.etotals, label=label, marker=marker, **kwargs)
            ax.set_ylabel('Energy (eV)')

        elif what == "abc":
            # Lattice parameters.
            mark = kwargs.pop("marker", None)
            markers = ["o", "^", "v"] if mark is None else 3 * [mark]
            for i, label in enumerate(["a", "b", "c"]):
                ax.plot(self.steps, [s.lattice.abc[i] for s in self.structures], label=label,
                        marker=markers[i], **kwargs)
            ax.set_ylabel("abc (A)")

        elif what in ("a", "b", "c"):
            i =  ("a", "b", "c").index(what)
            marker = kwargs.pop("marker", None)
            if marker is None:
                marker = {"a": "o", "b": "^", "c": "v"}[what]
            label = kwargs.pop("label", what)
            ax.plot(self.steps, [s.lattice.abc[i] for s in self.structures], label=label,
                    marker=marker, **kwargs)
            ax.set_ylabel('%s (A)' % what)

        elif what == "angles":
            # Lattice Angles
            mark = kwargs.pop("marker", None)
            markers = ["o", "^", "v"] if mark is None else 3 * [mark]
            for i, label in enumerate(["alpha", "beta", "gamma"]):
                ax.plot(self.steps, [s.lattice.angles[i] for s in self.structures], label=label,
                        marker=markers[i], **kwargs)
            ax.set_ylabel(r"$\alpha\beta\gamma$ (degree)")

        elif what in ("alpha", "beta", "gamma"):
            i =  ("alpha", "beta", "gamma").index(what)
            marker = kwargs.pop("marker", None)
            if marker is None:
                marker = {"alpha": "o", "beta": "^", "gamma": "v"}[what]

            label = kwargs.pop("label", what)
            ax.plot(self.steps, [s.lattice.angles[i] for s in self.structures], label=label,
                    marker=marker, **kwargs)
            ax.set_ylabel(r"$\%s$ (degree)" % what)

        elif what == "volume":
            marker = kwargs.pop("marker", "o")
            ax.plot(self.steps, [s.lattice.volume for s in self.structures], marker=marker, **kwargs)
            ax.set_ylabel(r'$V\, (A^3)$')

        elif what == "pressure":
            stress_cart_tensors, pressures = self.reader.read_cart_stress_tensors()
            marker = kwargs.pop("marker", "o")
            label = kwargs.pop("label", "P")
            ax.plot(self.steps, pressures, label=label, marker=marker, **kwargs)
            ax.set_ylabel('P (GPa)')

        elif what == "forces":
            forces_hist = self.reader.read_cart_forces()
            fmin_steps, fmax_steps, fmean_steps, fstd_steps = [], [], [], []
            for step in range(self.num_steps):
                forces = forces_hist[step]
                fmods = np.sqrt([np.dot(force, force) for force in forces])
                fmean_steps.append(fmods.mean())
                fstd_steps.append(fmods.std())
                fmin_steps.append(fmods.min())
                fmax_steps.append(fmods.max())

            mark = kwargs.pop("marker", None)
            markers = ["o", "^", "v", "X"] if mark is None else 4 * [mark]
            ax.plot(self.steps, fmin_steps, label="min |F|", marker=markers[0], **kwargs)
            ax.plot(self.steps, fmax_steps, label="max |F|", marker=markers[1], **kwargs)
            ax.plot(self.steps, fmean_steps, label="mean |F|", marker=markers[2], **kwargs)
            ax.plot(self.steps, fstd_steps, label="std |F|", marker=markers[3], **kwargs)
            label = "std |F"
            ax.set_ylabel('F stats (eV/A)')

        else:
            raise ValueError("Invalid value for what: `%s`" % str(what))

        ax.set_xlabel('Step')
        ax.grid(True)
        if label is not None:
            ax.legend(loc='best', fontsize=fontsize, shadow=True)


    @add_fig_kwargs
    def plot(self, ax_list=None, fontsize=8, **kwargs):
        """
        Plot the evolution of structural parameters (lattice lengths, angles and volume)
        as well as pressure, info on forces and total energy.

        Args:
            ax_list: List of |matplotlib-Axes|. If None, a new figure is created.
            fontsize: fontsize for legend

        Returns: |matplotlib-Figure|
        """
        what_list = ["abc", "angles", "volume", "pressure", "forces", "energy"]
        nrows, ncols = 3, 2
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=False, squeeze=False)
        ax_list = ax_list.ravel()
        assert len(ax_list) == len(what_list)

        for what, ax in zip(what_list, ax_list):
            self.plot_ax(ax, what, fontsize=fontsize, marker="o")

        return fig

    @add_fig_kwargs
    def plot_energies(self, ax=None, fontsize=12, **kwargs):
        """
        Plot the total energies as function of the iteration step.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: Legend and title fontsize.

        Returns: |matplotlib-Figure|
        """
        # TODO max force and pressure
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        terms = self.reader.read_eterms()
        for key, values in terms.items():
            if np.all(values == 0.0): continue
            ax.plot(self.steps, values, marker="o", label=key)

        ax.set_xlabel('Step')
        ax.set_ylabel('Energies (eV)')
        ax.grid(True)
        ax.legend(loc='best', fontsize=fontsize, shadow=True)

        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        yield self.plot(show=False)
        yield self.plot_energies(show=False)

    def mvplot_trajectories(self, colormap="hot", sampling=1, figure=None, show=True,
                            with_forces=True, **kwargs):  # pragma: no cover
        """
        Call mayavi_ to plot atomic trajectories and the variation of the unit cell.
        """
        from abipy.display import mvtk
        figure, mlab = mvtk.get_fig_mlab(figure=figure)
        style = "labels"
        line_width = 100
        mvtk.plot_structure(self.initial_structure, style=style, unit_cell_color=(1, 0, 0), figure=figure)
        mvtk.plot_structure(self.final_structure, style=style, unit_cell_color=(0, 0, 0), figure=figure)

        steps = np.arange(start=0, stop=self.num_steps, step=sampling)
        xcart_list = self.reader.read_value("xcart") * units.bohr_to_ang
        for iatom in range(self.reader.natom):
            x, y, z = xcart_list[::sampling, iatom, :].T
            #for i in zip(x, y, z): print(i)
            trajectory = mlab.plot3d(x, y, z, steps, colormap=colormap, tube_radius=None,
                                    line_width=line_width, figure=figure)
            mlab.colorbar(trajectory, title='Iteration', orientation='vertical')

        if with_forces:
            fcart_list = self.reader.read_cart_forces(unit="eV ang^-1")
            for iatom in range(self.reader.natom):
                x, y, z = xcart_list[::sampling, iatom, :].T
                u, v, w = fcart_list[::sampling, iatom, :].T
                q = mlab.quiver3d(x, y, z, u, v, w, figure=figure, colormap=colormap,
                                  line_width=line_width, scale_factor=10)
                #mlab.colorbar(q, title='Forces [eV/Ang]', orientation='vertical')

        if show: mlab.show()
        return figure

    def mvanimate(self, delay=500):  # pragma: no cover
        from abipy.display import mvtk
        figure, mlab = mvtk.get_fig_mlab(figure=None)
        style = "points"
        #mvtk.plot_structure(self.initial_structure, style=style, figure=figure)
        #mvtk.plot_structure(self.final_structure, style=style, figure=figure)

        xcart_list = self.reader.read_value("xcart") * units.bohr_to_ang
        #t = np.arange(self.num_steps)
        #line_width = 2
        #for iatom in range(self.reader.natom):
        #    x, y, z = xcart_list[:, iatom, :].T
        #    trajectory = mlab.plot3d(x, y, z, t, colormap=colormap, tube_radius=None, line_width=line_width, figure=figure)
        #mlab.colorbar(trajectory, title='Iteration', orientation='vertical')

        #x, y, z = xcart_list[0, :, :].T
        #nodes = mlab.points3d(x, y, z)
        #nodes.glyph.scale_mode = 'scale_by_vector'
        #this sets the vectors to be a 3x5000 vector showing some random scalars
        #nodes.mlab_source.dataset.point_data.vectors = np.tile( np.random.random((5000,)), (3,1))
        #nodes.mlab_source.dataset.point_data.scalars = np.random.random((5000,))

        @mlab.show
        @mlab.animate(delay=delay, ui=True)
        def anim():
            """Animate."""
            for it, structure in enumerate(self.structures):
            #for it in range(self.num_steps):
                print('Updating scene for iteration:', it)
                #mlab.clf(figure=figure)
                mvtk.plot_structure(structure, style=style, figure=figure)
                #x, y, z = xcart_list[it, :, :].T
                #nodes.mlab_source.set(x=x, y=y, z=z)
                #figure.scene.render()
                mlab.draw(figure=figure)
                yield

        anim()

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("hist = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(hist)"),
            nbv.new_code_cell("hist.plot_energies();"),
            nbv.new_code_cell("hist.plot();"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class HistRobot(Robot):
    """
    This robot analyzes the results contained in multiple HIST.nc_ files.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: HistRobot
    """
    EXT = "HIST"

    def to_string(self, verbose=0):
        """String representation with verbosity level ``verbose``."""
        s = ""
        if verbose:
            s = super(HistRobot, self).to_string(verbose=0)
        df = self.get_dataframe()
        s_df = "Table with final structures, pressures in GPa and force stats in eV/Ang:\n\n%s" % str(df)
        if s:
            return "\n".join([s, str(s_df)])
        else:
            return str(s_df)

    def get_dataframe(self, with_geo=True, index=None, abspath=False, with_spglib=True, funcs=None, **kwargs):
        """
        Return a |pandas-DataFrame| with the most important final results and the filenames as index.

        Args:
            with_geo: True if structure info should be added to the dataframe
            abspath: True if paths in index should be absolute. Default: Relative to getcwd().
            index: Index of the dataframe, if None, robot labels are used
            with_spglib: If True, spglib_ is invoked to get the space group symbol and number

        kwargs:
            attrs:
                List of additional attributes of the |GsrFile| to add to the |pandas-DataFrame|.
            funcs: Function or list of functions to execute to add more data to the DataFrame.
                Each function receives a |GsrFile| object and returns a tuple (key, value)
                where key is a string with the name of column and value is the value to be inserted.
        """
        # Add attributes specified by the users
        attrs = [
            "num_steps", "final_energy", "final_pressure",
            "final_fmin", "final_fmax", "final_fmean", "final_fstd", "final_drift",
            "initial_fmin", "initial_fmax", "initial_fmean", "initial_fstd", "initial_drift",
            # TODO add more columns but must update HIST file
            #"nsppol", "nspinor", "nspden",
            #"ecut", "pawecutdg", "tsmear", "nkpt",
        ] + kwargs.pop("attrs", [])

        rows, row_names = [], []
        for label, hist in self.items():
            row_names.append(label)
            d = OrderedDict()

            initial_fstas_dict = hist.get_fstats_dict(step=0)
            final_fstas_dict = hist.get_fstats_dict(step=-1)

            # Add info on structure.
            if with_geo:
                d.update(hist.final_structure.get_dict4pandas(with_spglib=with_spglib))

            for aname in attrs:
                if aname in ("final_fmin", "final_fmax", "final_fmean", "final_fstd", "final_drift",):
                    value = final_fstas_dict[aname.replace("final_", "")]
                elif aname in ("initial_fmin", "initial_fmax", "initial_fmean", "initial_fstd", "initial_drift"):
                    value = initial_fstas_dict[aname.replace("initial_", "")]
                else:
                    value = getattr(hist, aname, None)
                d[aname] = value

            # Execute functions
            if funcs is not None: d.update(self._exec_funcs(funcs, hist))
            rows.append(d)

        import pandas as pd
        row_names = row_names if not abspath else self._to_relpaths(row_names)
        index = row_names if index is None else index
        return pd.DataFrame(rows, index=index, columns=list(rows[0].keys()))

    @property
    def what_list(self):
        """List with all quantities that can be plotted (what_list)."""
        return ["energy", "abc", "angles", "volume", "pressure", "forces"]

    @add_fig_kwargs
    def gridplot(self, what_list=None, sharex="row", sharey="row", fontsize=8, **kwargs):
        """
        Plot the ``what`` value extracted from multiple HIST.nc_ files on a grid.

        Args:
            what_list: List of quantities to plot.
                Must be in ["energy", "abc", "angles", "volume", "pressure", "forces"]
            sharex: True if xaxis should be shared.
            sharey: True if yaxis should be shared.
            fontsize: fontsize for legend.

        Returns: |matplotlib-Figure|
        """
        what_list = list_strings(what_list) if what_list is not None else self.what_list

        # Build grid of plots.
        nrows, ncols = len(what_list), len(self)

        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=sharex, sharey=sharey, squeeze=False)
        ax_mat = np.reshape(ax_mat, (nrows, ncols))

        for irow, what in enumerate(what_list):
            for icol, hist in enumerate(self.abifiles):
                ax = ax_mat[irow, icol]
                ax.grid(True)
                hist.plot_ax(ax_mat[irow, icol], what, fontsize=fontsize, marker="o")

                if irow == 0:
                    ax.set_title(hist.relpath, fontsize=fontsize)
                if irow != nrows - 1:
                    set_visible(ax, False, "xlabel")
                if icol != 0:
                    set_visible(ax, False, "ylabel")

        return fig

    @add_fig_kwargs
    def combiplot(self, what_list=None, colormap="jet", fontsize=6, **kwargs):
        """
        Plot multiple HIST.nc_ files on a grid. One plot for each ``what`` value.

        Args:
            what_list: List of strings with the quantities to plot. If None, all quanties are plotted.
            colormap: matplotlib color map.
            fontsize: fontisize for legend.

        Returns: |matplotlib-Figure|.
        """
        what_list = (list_strings(what_list) if what_list is not None
            else ["energy", "a", "b", "c", "alpha", "beta", "gamma", "volume", "pressure"])

        num_plots, ncols, nrows = len(what_list), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots // ncols) + (num_plots % ncols)

        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=False, squeeze=False)
        ax_list = ax_list.ravel()
        cmap = plt.get_cmap(colormap)

        for i, (ax, what) in enumerate(zip(ax_list, what_list)):
            for ih, hist in enumerate(self.abifiles):
                label= None if i != 0 else hist.relpath
                hist.plot_ax(ax, what, color=cmap(ih / len(self)), label=label, fontsize=fontsize)

            if label is not None:
                ax.legend(loc="best", fontsize=fontsize, shadow=True)

            if i == len(ax_list) - 1:
                ax.set_xlabel("Step")
            else:
                ax.set_xlabel("")

        # Get around a bug in matplotlib.
        if num_plots % ncols != 0: ax_list[-1].axis('off')

        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        yield self.gridplot(show=False)
        yield self.combiplot(show=False)

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        args = [(l, f.filepath) for l, f in self.items()]
        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("robot = abilab.HistRobot(*%s)\nrobot.trim_paths()\nrobot" % str(args)),
            nbv.new_code_cell("robot.get_dataframe()"),
            nbv.new_code_cell("for what in robot.what_list: robot.gridplot(what=what, tight_layout=True);"),
        ])

        # Mixins
        #nb.cells.extend(self.get_baserobot_code_cells())

        return self._write_nb_nbpath(nb, nbpath)


class HistReader(ETSF_Reader):
    """
    This object reads data from the HIST file.


    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: HistReader
    """

    @lazy_property
    def num_steps(self):
        """Number of iterations present in the HIST.nc_ file."""
        return self.read_dimvalue("time")

    @lazy_property
    def natom(self):
        """Number of atoms un the unit cell."""
        return self.read_dimvalue("natom")

    def read_all_structures(self):
        """Return the list of structures at the different iteration steps."""
        rprimd_list = self.read_value("rprimd")
        xred_list = self.read_value("xred")

        # Alchemical mixing is not supported.
        num_pseudos = self.read_dimvalue("npsp")
        ntypat = self.read_dimvalue("ntypat")
        if num_pseudos != ntypat:
            raise NotImplementedError("Alchemical mixing is not supported, num_pseudos != ntypat")

        znucl, typat = self.read_value("znucl"), self.read_value("typat").astype(int)
        #print(znucl.dtype, typat)
        cart_forces_step = self.read_cart_forces(unit="eV ang^-1")

        structures = []
        #print("typat", type(typat))
        for step in range(self.num_steps):
            s = Structure.from_abivars(
                xred=xred_list[step],
                rprim=rprimd_list[step],
                acell=3 * [1.0],
                # FIXME ntypat, typat, znucl are missing!
                znucl=znucl,
                typat=typat,
            )
            s.add_site_property("cartesian_forces", cart_forces_step[step])
            structures.append(s)

        return structures

    def read_eterms(self, unit="eV"):
        """|AttrDict| with the decomposition of the total energy in units ``unit``"""
        return AttrDict(
            etotals=units.EnergyArray(self.read_value("etotal"), "Ha").to(unit),
            kinetic_terms=units.EnergyArray(self.read_value("ekin"), "Ha").to(unit),
            entropies=units.EnergyArray(self.read_value("entropy"), "Ha").to(unit),
        )

    def read_cart_forces(self, unit="eV ang^-1"):
        """
        Read and return a |numpy-array| with the cartesian forces in unit ``unit``.
        Shape (num_steps, natom, 3)
        """
        return units.ArrayWithUnit(self.read_value("fcart"), "Ha bohr^-1").to(unit)

    def read_reduced_forces(self):
        """
        Read and return a |numpy-array| with the forces in reduced coordinates
        Shape (num_steps, natom, 3)
        """
        return self.read_value("fred")

    def read_cart_stress_tensors(self):
        """
        Return the stress tensors (nstep x 3 x 3) in cartesian coordinates (GPa)
        and the list of pressures in GPa unit.
        """
        # Abinit stores 6 unique components of this symmetric 3x3 tensor:
        # Given in order (1,1), (2,2), (3,3), (3,2), (3,1), (2,1).
        c = self.read_value("strten")
        tensors = np.empty((self.num_steps, 3, 3), dtype=np.float)

        for step in range(self.num_steps):
            for i in range(3): tensors[step, i,i] = c[step, i]
            for p, (i, j) in enumerate(((2,1), (2,0), (1,0))):
                tensors[step, i,j] = c[step, 3+p]
                tensors[step, j,i] = c[step, 3+p]

        tensors *= abu.HaBohr3_GPa
        pressures = np.empty(self.num_steps)
        for step, tensor in enumerate(tensors):
            pressures[step] = - tensor.trace() / 3

        return tensors, pressures
