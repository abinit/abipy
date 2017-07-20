# coding: utf-8
"""Robots."""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import numpy as np
import pandas as pd

from collections import OrderedDict, deque
from monty.string import is_string, list_strings
from monty.termcolor import cprint
from monty.collections import dict2namedtuple
#from monty.functools import lazy_property
from pymatgen.analysis.eos import EOS
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt
from abipy.flowtk import Flow
from abipy.core.mixins import NotebookWriter


#__all__ = [
#    "abirobot",
#]


def abirobot(obj, ext, nids=None):
    """
    Factory function that builds and return the :class:`Robot` subclass from the file
    extension `ext`. `obj` can be a directory path, or a :class:`Flow` instance.
    `nids` is an optional list of node identifiers used to filter the tasks in the flow.

    Usage example:

    .. code-block:: python

        with abirobot(flow, "GSR") as robot:
            # do something with robot and close the GSR files when done.

        with abirobot("dirpath", "SIGRES") as robot:
            # do something with robot and close the SIGRES files when done.
    """
    for cls in Robot.__subclasses__():
        if cls.EXT in (ext, ext.upper()):
            return cls.open(obj, nids=nids)

    raise ValueError("Cannot find Robot subclass associated to extension %s\n" % ext +
                     "The list of supported extensions is:\n%s" %
                     [cls.EXT for cls in Robot.__subclasses__()])


class Robot(object):
    """
    The main function of a `Robot` is facilitating the extraction of the output data produced by
    multiple tasks in a :class:`Flow`. This is the base class from which all Robot subclasses should derive.
    A Robot supports the `with` context manager:

    Usage example:

    .. code-block:: python

        with Robot([("label1", "file1"), (label2, "file2")]) as robot:
            # Do something with robot. files are automatically closed when we exit.
    """
    # TODO
    # 1) Abstract interface from collections
    # 2) should __iter__  return (label, ncfile) or ncfile (not __getitem__ returns ncfiles.__getitem__ !!!
    # 3) replace ncfiles with files just to be consistent since we have DdbRobot!

    # filepaths are relative to `start`. None for asbolute paths. This flag is set in trim_paths
    start = None

    def __init__(self, *args):
        """
        Args:
            args is a list of tuples (label, filepath)
        """
        self._ncfiles, self._do_close = OrderedDict(), OrderedDict()
        self._exceptions = deque(maxlen=100)

        for label, ncfile in args:
            self.add_file(label, ncfile)

    @classmethod
    def class_for_ext(cls, ext):
        """Return the Robot subclass associated to the given extension."""
        for subcls in cls.__subclasses__():
            if subcls.EXT in (ext, ext.upper()):
                return subcls

        raise ValueError("Cannot find Robot subclass associated to extension %s\n" % ext +
                         "The list of supported extensions is:\n%s" %
                         [cls.EXT for cls in Robot.__subclasses__()])

    # Deprecated. Use class_for_ext
    for_ext = class_for_ext

    @classmethod
    def from_dir(cls, top, walk=True):
        """
        This class method builds a robot by scanning all files located within directory `top`.
        This method should be invoked with a concrete robot class, for example:

            robot = GsrRobot.from_dir(".")

        Args:
            top (str): Root directory
	    walk: if True, directories inside `top` are included as well.
        """
        return cls(*cls._open_files_in_dir(top, walk))

    @classmethod
    def _open_files_in_dir(cls, top, walk):
        """Open files in directory tree starting from `top`. Return list of Abinit files."""
        from abipy.abilab import abiopen
        items = []
        if walk:
            for dirpath, dirnames, filenames in os.walk(top):
                filenames = [f for f in filenames if cls.class_handles_filename(f)]
                for f in filenames:
                    ncfile = abiopen(os.path.join(dirpath, f))
                    if ncfile is not None: items.append((ncfile.filepath, ncfile))
        else:
            filenames = [f for f in os.listdir(top) if cls.class_handles_filename(f)]
            for f in filenames:
                ncfile = abiopen(os.path.join(top, f))
                if ncfile is not None: items.append((ncfile.filepath, ncfile))

        return items

    @classmethod
    def class_handles_filename(cls, filename):
        """True if robot class handles filename."""
        return filename.endswith("_" + cls.EXT + ".nc")

    @classmethod
    def from_files(cls, filenames):
        """Build a Robot from a list of `filenames`."""
        filenames = list_strings(filenames)
        from abipy.abilab import abiopen
        filenames = [f for f in filenames if cls.class_handles_filename(f)]
        items = []
        for f in filenames:
            try:
                ncfile = abiopen(f)
            except Exception:
                ncfile = None
            if ncfile is not None: items.append((ncfile.filepath, ncfile))
        return cls(*items)

    @classmethod
    def from_flow(cls, flow, outdirs="all", nids=None):
        """
        Build a robot from a Flow.

        Args:
            flow: :class:`Flow` object
            outdirs: String used to select/ignore the files in the output directory of flow, works and tasks
                outdirs="work" selects only the outdir of the Works,
                outdirs="flow+task" selects the outdir of the Flow and the outdirs of the tasks
                outdirs="-work" excludes the outdir of the Works.
                Cannot use `+` and `-` flags in the same string.
                Default: `all` that is equivalent to "flow+work+task"
            nids: List of node identifiers used to select particular nodes. Not used if None

        Returns:
            `Robot` subclass.
        """
        robot = cls()
        all_opts = ("flow", "work", "task")

        if outdirs == "all":
            tokens = all_opts
        elif "+" in outdirs:
            assert "-" not in outdirs
            tokens = outdirs.split("+")
        elif "-" in outdirs:
            assert "+" not in outdirs
            tokens = [s for s in all if s not in outdirs.split("-")]
        else:
            tokens = list_strings(outdirs)

        if not all(t in all_opts for t in tokens):
            raise ValueError("Wrong outdirs string %s" % outdirs)

        if "flow" in tokens:
            robot.add_extfile_of_node(flow, nids=nids)

        if "work" in tokens:
            for work in flow:
                robot.add_extfile_of_node(work, nids=nids)

        if "task" in tokens:
            for task in flow.iflat_tasks():
                #print("task %s, nids %s" %  (task, nids))
                robot.add_extfile_of_node(task, nids=nids)

        return robot

    def add_extfile_of_node(self, node, nids=None):
        """
        Add the file produced by this node to the robot.
        """
        if nids and node.node_id not in nids: return
        filepath = node.outdir.has_abiext(self.EXT)
        if filepath:
            try:
                label = os.path.relpath(filepath)
            except OSError:
                # current working directory may not be defined!
                label = filepath

            self.add_file(label, filepath)

    def scan_dir(self, top, walk=True):
        """
        Scan directory tree starting from `top`. Add files to the robot instance.

        Args:
            top (str): Root directory
            walk: if True, directories inside `top` are included as well.

        Return:
            Number of files found.
	"""
        count = 0
        for filepath, ncfile in self.__class__._open_files_in_dir(top, walk):
            count += 1
            self.add_file(filepath, ncfile)

        return count

    def add_file(self, label, ncfile):
        """
        Add a file to the robot with the given label.

        Args:
            label: String used to identify the file (must be unique, ax exceptions is
                raised if label is already present.
            ncfile: Specify the file to be added. Accepts strings (filepath) or abipy file-like objects.
        """
        if is_string(ncfile):
            from abipy.abilab import abiopen
            ncfile = abiopen(ncfile)
            # Open file here --> have to close it.
            self._do_close[ncfile.filepath] = True

        if label in self._ncfiles:
            raise ValueError("label %s is already present!")

        self._ncfiles[label] = ncfile

    #def pop_filepath(self, filepath):
    #    """
    #    Remove the file with the given `filepath` and close it.
    #    """
    #    if label, ncfile in self._ncfiles.items():
    #        if ncfile.filepath != filepath: continue
    #        self._ncfiles.pop(label)
    #        if self._do_close.pop(ncfile.filepath, False):
    #            try:
    #                ncfile.close()
    #            except Exception as exc:
    #                print("Exception while closing: ", ncfile.filepath)
    #                print(exc)
    #                #raise

    def pop_label(self, label):
        """
        Remove file with the given `label` and close it.
        """
        if label in self._ncfiles:
            ncfile = self._ncfiles.pop(label)
            if self._do_close.pop(ncfile.filepath, False):
                try:
                    ncfile.close()
                except Exception as exc:
                    print("Exception while closing: ", ncfile.filepath)
                    print(exc)
                    #raise

    def change_labels(self, new_labels, dryrun=False):
        """
        Change labels of the files. Return mapping new --> old.
        """
        if len(new_labels) != len(self):
            raise ValueError("Robot has %d files while len(new_labels) = %d" % (len(new_labels), len(self)))

        old_labels = list(self._ncfiles.keys())
        old_ncfiles = self._ncfiles
        self._ncfiles = OrderedDict()
        new2old = OrderedDict()
        for old, new in zip(old_labels, new_labels):
            print("old [%s] --> new [%s]" % (old, new))
            new2old[new] = old
            if not dryrun:
                self._ncfiles[new] = old_ncfiles[old]

        return new2old

    def trim_paths(self, start=None):
        """
        Replace absolute filepaths in the robot with relative paths wrt to `start` directory.
        If start is None, os.getcwd() is used. Set `self.start` attribute, return self.start.
        """
        self.start = os.getcwd() if start is None else start
        old_paths = list(self._ncfiles.keys())
        old_new_paths = [(p, os.path.relpath(os.path.abspath(p), start=self.start)) for p in old_paths]

        old_ncfiles = self._ncfiles
        self._ncfiles = OrderedDict()
        for old, new in old_new_paths:
            self._ncfiles[new] = old_ncfiles[old]

        return self.start

    @property
    def exceptions(self):
        """List of exceptions."""
        return self._exceptions

    def __len__(self):
        return len(self._ncfiles)

    def __iter__(self):
        return iter(self._ncfiles.items())

    def __getitem__(self, key):
        # self[key]
        return self.ncfiles.__getitem__(key)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Activated at the end of the with statement."""
        self.close()

    def items(self):
        return self._ncfiles.items()

    def show_files(self, stream=sys.stdout):
        s = "\n".join(["%s --> %s" % (label, ncfile.filepath) for label, ncfile in self])
        stream.write(s)

    def __repr__(self):
        """Invoked by repr(obj)"""
        return self.to_string(func=repr)

    def __str__(self):
        """Invoked by str."""
        return self.to_string(func=str)

    def to_string(self, func=str, verbose=0):
        """String representation."""
        lines = ["%s with %d files in memory" % (self.__class__.__name__, len(self.ncfiles))]
        app = lines.append
        for i, f in enumerate(self.ncfiles):
            app(" ")
            app(func(f))
            app(" ")

        return "\n".join(lines)

    @property
    def ncfiles(self):
        """List of netcdf files."""
        return list(self._ncfiles.values())

    def close(self):
        """
        Close all files that have been opened by the Robot
        """
        for ncfile in self.ncfiles:
            if self._do_close.pop(ncfile.filepath, False):
                try:
                    ncfile.close()
                except:
                    print("Exception while closing: ", ncfile.filepath)
                    print(exc)
                    #raise

    @classmethod
    def open(cls, obj, nids=None, **kwargs):
        """
        Flexible constructor. obj can be a :class:`Flow` or a string with the directory containing the Flow.
        `nids` is an optional list of :class:`Node` identifiers used to filter the set of :class:`Task` in the Flow.
        """
        has_dirpath = False
        if is_string(obj):
            try:
                obj = Flow.pickle_load(obj)
            except:
                has_dirpath = True

        if not has_dirpath:
            # We have a Flow. smeth is the name of the Task method used to open the file.
            items = []
            smeth = "open_" + cls.EXT.lower()
            for task in obj.iflat_tasks(nids=nids): #, status=obj.S_OK):
                open_method = getattr(task, smeth, None)
                if open_method is None: continue
                ncfile = open_method()
                if ncfile is not None: items.append((task.pos_str, ncfile))
            return cls(*items)

        else:
            # directory --> search for files with the appropriate extension and open it with abiopen.
            if nids is not None: raise ValueError("nids cannot be used when obj is a directory.")
            return cls.from_dir(obj)

    #def get_attributes(self, attr_name, obj=None, retdict=False):
    #    od = OrderedDict()
    #    for label, ncfile in self:
    #        obj = ncfile if obj is None else getattr(ncfile, obj)
    #        od[label] = getattr(obj, attr_name)

    #    if retdict:
    #        return od
    #    else:
    #        return list(od.values())

    def _exec_funcs(self, funcs, arg):
        """
        Execute list of callable functions. Each function receives arg as argument.
        """
        if not isinstance(funcs, (list, tuple)): funcs = [funcs]
        d = {}
        for func in funcs:
            try:
                key, value = func(arg)
                d[key] = value
            except Exception as exc:
                self._exceptions.append(str(exc))
        return d

    #def get_structure_dataframe(self):
    #    from abipy.core.structure import frames_from_structures
    #    index = list(self.keys())
    #    frames_from_structures(self.ncfiles, index=None, with_spglib=True, cart_coords=False)

    #def get_ebands_dataframe(self):
    #    from abipy.electrons.ebands import frames_from_ebands
    #    for abifile in self.ncfiles:
    #    index = list(self.keys())
    #    frame_from_ebands(self.ncfiles, index=None, with_spglib=True)

    #def ncread_and_plot_variables(self, varname_x, varname_y, hspan=None, **kwargs):
    #    """
    #    Ex:
    #        plot_variables("ecut", "etotal")
    #    """
    #    title = kwargs.pop("title", None)
    #    show = kwargs.pop("show", True)
    #    savefig = kwargs.pop("savefig", None)

    #    # Read the value of varname from the files.
    #    xx, yy = [], []
    #    for filepath in self.filepaths:
    #        with GsrReader(filepath) as r:
    #            xx.append(r.read_value(varname_x))
    #            yy.append(r.read_value(varname_y))

    #    import matplotlib.pyplot as plt
    #    fig = plt.figure()
    #    ax = fig.add_subplot(1,1,1)

    #    ax.plot(xx, yy, "o-", **kwargs)

    #    if hspan is not None:
    #        last = yy[-1]
    #        ax.axhspan(last - hspan, last + hspan, facecolor='0.5', alpha=0.5)

    #    if title is not None: fig.suptitle(title)
    #    if show: plt.show()
    #    if savefig is not None: fig.savefig(savefig)

    #    return fig


class GsrRobot(Robot, NotebookWriter):
    """
    This robot analyzes the results contained in multiple GSR files.
    """
    EXT = "GSR"

    def get_dataframe(self, with_geo=True, **kwargs):
        """
        Return a pandas DataFrame with the most important GS results.

        Args:
            with_geo: True if structure info should be added to the dataframe

        kwargs:
            attrs:
                List of additional attributes of the :class:`GsrFile` to add to
                the pandas :class:`DataFrame`
            funcs:
                Function or list of functions to execute to add more data to the DataFrame.
                Each function receives a :class:`GsrFile` object and returns a tuple (key, value)
                where key is a string with the name of column and value is the value to be inserted.
        """
        # Add attributes specified by the users
        # TODO add more columns
        attrs = [
            "energy", "pressure", "max_force",
            "ecut", "pawecutdg",
            "tsmear", "nkpt",
            "nsppol", "nspinor", "nspden",
        ] + kwargs.pop("attrs", [])

        rows, row_names = [], []
        for label, gsr in self:
            row_names.append(label)
            d = OrderedDict()
            for aname in attrs:
                if aname == "nkpt":
                    value = len(gsr.ebands.kpoints)
                else:
                    value = getattr(gsr, aname, None)
                    if value is None: value = getattr(gsr.ebands, aname, None)
                d[aname] = value

            # Add info on structure.
            if with_geo:
                d.update(gsr.structure.get_dict4frame(with_spglib=True))

            # Execute functions
            d.update(self._exec_funcs(kwargs.get("funcs", []), gsr))
            rows.append(d)

        return pd.DataFrame(rows, index=row_names, columns=list(rows[0].keys()))

    def get_ebands_plotter(self, cls=None):
        """
        Build and return an instance of `ElectronBandsPlotter` or a subclass is cls is not None.

        Args:
            cls: subclass of `ElectronBandsPlotter`
        """
        from abipy.electrons.ebands import ElectronBandsPlotter
        plotter = ElectronBandsPlotter() if cls is None else cls()

        for label, gsr in self:
            plotter.add_ebands(label, gsr.ebands)

        return plotter

    def get_edos_plotter(self, cls=None, **kwargs):
        """
        Build and return an instance of `ElectronDosPlotter` or a subclass is cls is not None.

        Args:
            cls: subclass of `ElectronDosPlotter`
            kwargs: Arguments passed to ebands.get_edos
        """
        from abipy.electrons.ebands import ElectronDosPlotter
        plotter = ElectronDosPlotter() if cls is None else cls()

        for label, gsr in self:
            if not gsr.ebands.kpoints.is_ibz:
                cprint("Skipping %s because kpoint sampling not IBZ" % gsr.filepath, "magenta")
                continue
            plotter.add_edos(label, gsr.ebands.get_edos(**kwargs))

        return plotter

    # FIXME: EOS has been changed in pymatgen.
    def eos_fit(self, eos_name="murnaghan"):
        """
        Fit energy as function of volume to get the equation of state, equilibrium volume,
        bulk modulus and its derivative wrt to pressure.

        Args:
            eos_name:
                For the list of available models, see pymatgen.analysis.eos.

        Return
        """
        # Read volumes and energies from the GSR files.
        energies, volumes = [], []
        for label, gsr in self:
            energies.append(gsr.energy)
            volumes.append(gsr.structure.volume)

        # Note that eos.fit expects lengths in Angstrom, and energies in eV.
        # I'm also monkey-patching the plot method.
        if eos_name != "all":
            fit = EOS(eos_name=eos_name).fit(volumes, energies)
            fit.plot = my_fit_plot
            return fit
        else:
            # Use all the available models.
            fits, rows = [], []
            models = list(EOS.MODELS.keys())
            for eos_name in models:
                fit = EOS(eos_name=eos_name).fit(volumes, energies)
                fit.plot = my_fit_plot
                fits.append(fit)
                rows.append(OrderedDict([(aname, getattr(fit, aname)) for aname in
                    ("v0", "e0", "b0_GPa", "b1")]))

            frame = pd.DataFrame(rows, index=models, columns=list(rows[0].keys()))
            return fits, frame
            #return dict2namedtuple(fits=fits, frame=frame)

    #def get_phasediagram_results(self):
    #    from abipy.core.restapi import PhaseDiagramResults
    #    entries = []
    #    for label, gsr in self:
    #        entries.append(gsr.get_computed_entry(inc_structure=True, parameters=None, data=None))
    #    return PhaseDiagramResults(entries)

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        args = [(l, f.filepath) for l, f in self.items()]
        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("robot = abilab.GsrRobot(*%s)\nrobot.trim_paths()\nprint(robot)" % str(args)),
            nbv.new_code_cell("ebands_plotter = robot.get_ebands_plotter()"),
            nbv.new_code_cell("df = ebands_plotter.get_ebands_frame()\ndisplay(df)"),
            nbv.new_code_cell("ebands_plotter.ipw_select_plot()"),
            nbv.new_code_cell("#anim = ebands_plotter.animate()"),
            nbv.new_code_cell("edos_plotter = robot.get_edos_plotter()"),
            nbv.new_code_cell("edos_plotter.ipw_select_plot()"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


@add_fig_kwargs
def my_fit_plot(self, ax=None, **kwargs):
    """
    Plot the equation of state.

    Args:

    Returns:
        Matplotlib figure object.
    """
    ax, fig, plt = get_ax_fig_plt(ax=ax)

    color = kwargs.get("color", "r")
    label = kwargs.get("label", "{} fit".format(self.__class__.__name__))
    lines = ["Equation of State: %s" % self.__class__.__name__,
             "Minimum energy = %1.2f eV" % self.e0,
             "Minimum or reference volume = %1.2f Ang^3" % self.v0,
             "Bulk modulus = %1.2f eV/Ang^3 = %1.2f GPa" %
             (self.b0, self.b0_GPa),
             "Derivative of bulk modulus wrt pressure = %1.2f" % self.b1]
    text = "\n".join(lines)
    text = kwargs.get("text", text)

    # Plot input data.
    ax.plot(self.volumes, self.energies, linestyle="None", marker="o", color=color)

    # Plot eos fit.
    vmin, vmax = min(self.volumes), max(self.volumes)
    vmin, vmax = (vmin - 0.01 * abs(vmin), vmax + 0.01 * abs(vmax))
    vfit = np.linspace(vmin, vmax, 100)

    ax.plot(vfit, self.func(vfit), linestyle="dashed", color=color, label=label)

    ax.grid(True)
    ax.xlabel("Volume $\\AA^3$")
    ax.ylabel("Energy (eV)")
    ax.legend(loc="best", shadow=True)
    # Add text with fit parameters.
    ax.text(0.4, 0.5, text, transform=ax.transAxes)

    return fig


class SigresRobot(Robot, NotebookWriter):
    """
    This robot analyzes the results contained in multiple SIGRES files.
    """
    EXT = "SIGRES"

    def merge_dataframes_sk(self, spin, kpoint, **kwargs):
        for i, (label, sigr) in enumerate(self):
            frame = sigr.get_dataframe_sk(spin, kpoint, index=label)
            if i == 0:
                table = frame
            else:
                table = table.append(frame)

        return table

    def get_qpgaps_dataframe(self, spin=None, kpoint=None, with_geo=False, **kwargs):
        """
        Return a pandas DataFrame with the most important results for the given (spin, kpoint).

        Args:
            spin: Spin index.
            kpoint
            with_geo: True if structure info should be added to the dataframe
        """
        # TODO: Ideally one should select the k-point for which we have the fundamental gap for the given spin
        # TODO: In principle the SIGRES might have different k-points
        if spin is None: spin = 0
        if kpoint is None: kpoint = 0

        attrs = [
            "nsppol",
            #"nspinor", "nspden", #"ecut", "pawecutdg",
            #"tsmear", "nkibz",
        ] + kwargs.pop("attrs", [])

        rows, row_names = [], []
        for label, sigres in self:
            row_names.append(label)
            d = OrderedDict()
            for aname in attrs:
                d[aname] = getattr(sigres, aname, None)

            qpgap = sigres.get_qpgap(spin, kpoint)
            d.update({"qpgap": qpgap})

            # Add convergence parameters
            d.update(sigres.params)

            # Add info on structure.
            if with_geo:
                d.update(sigres.structure.get_dict4frame(with_spglib=True))

            # Execute functions.
            d.update(self._exec_funcs(kwargs.get("funcs", []), sigres))
            rows.append(d)

        return pd.DataFrame(rows, index=row_names, columns=list(rows[0].keys()))

    def plot_conv_qpgap(self, x_vars, show=True, **kwargs):
        """
        Plot the convergence of the Quasi-particle gap.
        kwargs are passed to :class:`seaborn.PairGrid`.
        """
        import matplotlib.pyplot as plt
        import seaborn.apionly as sns

        data = self.get_qpgaps_dataframe()
        print(list(data.keys()))
        grid = sns.PairGrid(data, x_vars=x_vars, y_vars="qpgap", **kwargs)
        grid.map(plt.plot, marker="o")
        grid.add_legend()
        if show:
            plt.show()

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        args = [(l, f.filepath) for l, f in self.items()]
        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("robot = abilab.SigresRobot(*%s)\nrobot.trim_paths()\nprint(robot)" % str(args)),
            #nbv.new_code_cell("df = robot.get_qpgaps_dataframe(spin=None, kpoint=None, with_geo=False, **kwargs)"),
            #nbv.new_code_cell("plotter = robot.get_ebands_plotter()"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class MdfRobot(Robot, NotebookWriter):
    """
    This robot analyzes the results contained in multiple MDF files.
    """
    EXT = "MDF"

    def get_multimdf_plotter(self, cls=None):
        """
        Return an instance of MultipleMdfPlotter to compare multiple dielectric functions.
        """
        from abipy.electrons.bse import MultipleMdfPlotter
        plotter = MultipleMdfPlotter() if cls is None else cls()

        for label, mdf in self:
            plotter.add_mdf_file(label, mdf)

        return plotter

    def get_dataframe(self, with_geo=False, **kwargs):
        """
        Build and return Pandas dataframe with the most import BSE results.

        Args:
            with_geo: True if structure info should be added to the dataframe
            funcs:

        Return:
            pandas DataFrame
        """
        rows, row_names = [], []
        for i, (label, mdf) in enumerate(self):
            row_names.append(label)
            d = OrderedDict([
                ("exc_mdf", mdf.exc_mdf),
                ("rpa_mdf", mdf.rpanlf_mdf),
                ("gwrpa_mdf", mdf.gwnlf_mdf),
            ])
            #d = {aname: getattr(mdf, aname) for aname in attrs}
            #d.update({"qpgap": mdf.get_qpgap(spin, kpoint)})

            # Add convergence parameters
            d.update(mdf.params)

            # Add info on structure.
            if with_geo:
                d.update(mdf.structure.get_dict4frame(with_spglib=True))

            # Execute functions.
            d.update(self._exec_funcs(kwargs.get("funcs", []), mdf))
            rows.append(d)

        return pd.DataFrame(rows, index=row_names, columns=list(rows[0].keys()))

    #@add_fig_kwargs
    #def plot_conv_mdf(self, hue, mdf_type="exc_mdf", **kwargs):
    #    import matplotlib.pyplot as plt
    #    frame = self.get_dataframe()
    #    grouped = frame.groupby(hue)

    #    fig, ax_list = plt.subplots(nrows=len(grouped), ncols=1, sharex=True, sharey=True, squeeze=True)

    #    for i, (hue_val, group) in enumerate(grouped):
    #        #print(group)
    #        mdfs = group[mdf_type]
    #        ax = ax_list[i]
    #        ax.set_title("%s = %s" % (hue, hue_val))
    #        for mdf in mdfs:
    #            mdf.plot_ax(ax)

    #    return fig

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        args = [(l, f.filepath) for l, f in self.items()]
        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("robot = abilab.MdfRobot(*%s)\nrobot.trim_paths()\nprint(robot)" % str(args)),
            nbv.new_code_cell("#df = robot.get_dataframe(with_geo=False"),
            nbv.new_code_cell("plotter = robot.get_multimdf_plotter()"),
            nbv.new_code_cell('fig = plotter.plot(mdf_type="exc", qview="avg", xlim=None, ylim=None)'),
            #nbv.new_code_cell("fig = plotter.combiboxplot()"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class DdbRobot(Robot, NotebookWriter):
    """
    This robot analyzes the results contained in multiple DDB files.
    """
    EXT = "DDB"

    @classmethod
    def class_handles_filename(cls, filename):
        """Exclude DDB.nc files."""
        return filename.endswith("_" + cls.EXT)

    #def get_qpoints_union(self):
    #    """
    #    Return numpy array with the q-points in reduced coordinates found in the DDB files.
    #    """
    #    qpoints = []
    #    for label, ddb in enumerate(self):
    #        qpoints.extend(q.frac_coords for q in ddb.qpoints if q not in qpoints)

    #    return np.array(qpoints)

    #def get_qpoints_intersection(self):
    #    """Return numpy array with the q-points in reduced coordinates found in the DDB files."""
    #    qpoints = []
    #    for label, ddb in enumerate(self):
    #        qpoints.extend(q.frac_coords for q in ddb.qpoints if q not in qpoints)
    #
    #    return np.array(qpoints)

    def get_dataframe_at_qpoint(self, qpoint=None, asr=2, chneut=1, dipdip=1, with_geo=True, **kwargs):
        """
	Call anaddb to compute the phonon frequencies at a single q-point using the DDB files treated
	by the robot and the given anaddb input arguments. Build and return a pandas dataframe with results

        Args:
            qpoint: Reduced coordinates of the qpoint where phonon modes are computed
            asr, chneut, dipdp: Anaddb input variable. See official documentation.
            with_geo: True if structure info should be added to the dataframe

        Return:
            pandas DataFrame
        """
        # If qpoint is None, all the DDB must contain have the same q-point .
        if qpoint is None:
            if not all(len(ddb.qpoints) == 1 for ddb in self.ncfiles):
                raise ValueError("Found more than one q-point in the DDB file. qpoint must be specified")

            qpoint = self[0].qpoints[0]
            if any(np.any(ddb.qpoints[0] != qpoint) for ddb in self.ncfiles):
                raise ValueError("All the q-points in the DDB files must be equal")

        rows, row_names = [], []
        for i, (label, ddb) in enumerate(self):
            row_names.append(label)
            d = OrderedDict()
            #d = {aname: getattr(ddb, aname) for aname in attrs}
            #d.update({"qpgap": mdf.get_qpgap(spin, kpoint)})

            # Call anaddb to get the phonon frequencies.
            phbands = ddb.anaget_phmodes_at_qpoint(qpoint=qpoint, asr=asr, chneut=chneut, dipdip=dipdip)
            freqs = phbands.phfreqs[0, :]  # [nq, nmodes]

            d.update({"mode" + str(i): freqs[i] for i in range(len(freqs))})

            # Add convergence parameters
            d.update(ddb.params)

            # Add info on structure.
            if with_geo:
                d.update(phbands.structure.get_dict4frame(with_spglib=True))

            # Execute functions.
            d.update(self._exec_funcs(kwargs.get("funcs", []), ddb))

            rows.append(d)

        return pd.DataFrame(rows, index=row_names, columns=list(rows[0].keys()))

    def plot_conv_phfreqs_qpoint(self, x_vars, qpoint=None, **kwargs):
        """
        Plot the convergence of the phonon frequencies.
        kwargs are passed to :class:`seaborn.PairGrid`.
        """
        import matplotlib.pyplot as plt
        import seaborn.apionly as sns

        # Get the dataframe for this q-point.
        data = self.get_dataframe_at_qpoint(qpoint=qpoint)

        y_vars = sorted([k for k in data if k.startswith("mode")])
        #print(y_vars)

        # Call seaborn.
        grid = sns.PairGrid(data, x_vars=x_vars, y_vars=y_vars, **kwargs)
        grid.map(plt.plot, marker="o")
        grid.add_legend()
        plt.show()

    # TODO Test
    def get_phonon_plotters(self, **kwargs):
        """
        Invoke anaddb to compute phonon bands and DOS using the arguments passed via **kwargs.
        Collect results and return `namedtuple` with the following attributes:

            phbands_plotter: `PhononBandsPlotter` object.
            phdos_plotter: `PhononDosPlotter` object.
        """
        if "workdir" in kwargs:
            raise ValueError("Cannot specify `workdir` when multiple DDB file are executed.")

        from abipy.dfpt.phonons import PhononBandsPlotter, PhononDosPlotter
        phbands_plotter, phdos_plotter = PhononBandsPlotter(), PhononDosPlotter()

        for label, ddb in self:
            # Invoke anaddb to get phonon bands and DOS.
            phbst_file, phdos_file = ddb.anaget_phbst_and_phdos_files(**kwargs)

            # Phonon frequencies with non analytical contributions, if calculated, are saved in anaddb.nc
            # Those results should be fetched from there and added to the phonon bands.
            if kwargs.get("lo_to_splitting", False):
                anaddb_path = os.path.join(os.path.dirname(phbst_file.filepath), "anaddb.nc")
                phbst_file.phbands.read_non_anal_from_file(anaddb_path)

            phbands_plotter.add_phbands(label, phbst_file, phdos=phdos_file)
            phdos_plotter.add_phdos(label, phdos=phdos_file.phdos)
            phbst_file.close()
            phdos_file.close()

        return dict2namedtuple(phbands_plotter=phbands_plotter, phdos_plotter=phdos_plotter)

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        get_phonon_plotters_kwargs = ( "\n"
            '\tnqsmall=10, ndivsm=20, asr=2, chneut=1, dipdip=1, dos_method="tetra",\n'
            '\tlo_to_splitting=False, ngqpt=None, qptbounds=None,\n'
            '\tanaddb_kwargs=None, verbose=0')

        args = [(l, f.filepath) for l, f in self.items()]
        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("robot = abilab.DdbRobot(*%s)\nrobot.trim_paths()\nprint(robot)" % str(args)),
            nbv.new_code_cell("#dfq = robot.get_dataframe_at_qpoint(qpoint=None)"),
            nbv.new_code_cell("r = robot.get_phonon_plotters(%s)" % get_phonon_plotters_kwargs),
            nbv.new_code_cell("r.phbands_plotter.get_phbands_frame()"),
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("r.phbands_plotter.ipw_select_plot()"),
            nbv.new_code_cell("r.phdos_plotter.ipw_select_plot()"),
            nbv.new_code_cell("r.phdos_plotter.ipw_harmonic_thermo()"),
        ])

        return self._write_nb_nbpath(nb, nbpath)