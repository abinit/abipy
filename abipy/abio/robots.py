# coding: utf-8
"""
This module defines the Robot BaseClass. Robots operates on multiple files and provide helper
functions to plot the data e.g. convergence studies and to build pandas dataframes from the output files.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import inspect

from collections import OrderedDict, deque
from functools import wraps
from monty.string import is_string, list_strings
from monty.termcolor import cprint
from abipy.core.mixins import NotebookWriter
from abipy.tools.plotting import plot_xy_with_hue


class Robot(NotebookWriter):
    """
    This is the base class from which all Robot subclasses should derive.
    A Robot supports the `with` context manager:

    Usage example:

    .. code-block:: python

        with Robot([("label1", "file1"), (label2, "file2")]) as robot:
            # Do something with robot. files are automatically closed when we exit.

    .. note::

        __iter__  returns (label, abifile) so:

        for label, abifile in self:
            print(label)
    """
    # filepaths are relative to `start`. None for asbolute paths. This flag is set in trim_paths
    start = None

    def __init__(self, *args):
        """
        Args:
            args is a list of tuples (label, filepath)
        """
        self._abifiles, self._do_close = OrderedDict(), OrderedDict()
        self._exceptions = deque(maxlen=100)

        for label, ncfile in args:
            self.add_file(label, ncfile)

    @classmethod
    def get_supported_extensions(self):
        """List of strings with extensions supported by Robot subclasses."""
        # This is needed to have all subclasses.
        from abipy.abilab import Robot
        return sorted([cls.EXT for cls in Robot.__subclasses__()])

    @classmethod
    def class_for_ext(cls, ext):
        """Return the Robot subclass associated to the given extension."""
        for subcls in cls.__subclasses__():
            if subcls.EXT in (ext, ext.upper()):
                return subcls

        raise ValueError("Cannot find Robot subclass associated to extension %s\n" % ext +
                         "The list of supported extensions (case insensitive) is:\n%s" %
                         str(cls.get_supported_extensions()))

    @classmethod
    def from_dir(cls, top, walk=True, abspath=False):
        """
        This class method builds a robot by scanning all files located within directory `top`.
        This method should be invoked with a concrete robot class, for example:

            robot = GsrRobot.from_dir(".")

        Args:
            top (str): Root directory
	    walk: if True, directories inside `top` are included as well.
            abspath: True if paths in index should be absolute. Default: Relative to `top`.
        """
        new = cls(*cls._open_files_in_dir(top, walk))
        if not abspath: new.trim_paths(start=top)
        return new

    @classmethod
    def from_dirs(cls, dirpaths, walk=True, abspath=False):
        """
        Similar to `from_dir` but accepts a list of directories instead of a single directory.

        Args:
	    walk: if True, directories inside `top` are included as well.
            abspath: True if paths in index should be absolute. Default: Relative to `top`.
        """
        items = []
        for top in list_strings(dirpaths):
            items.extend(cls._open_files_in_dir(top, walk))
        new = cls(*items)
        if not abspath: new.trim_paths(start=os.getcwd())
        return new

    @classmethod
    def from_dir_glob(cls, pattern, walk=True, abspath=False):
        """
        This class method builds a robot by scanning all files located within the directories
        matching `pattern` as implemented by glob.glob
        This method should be invoked with a concrete robot class, for example:

            robot = GsrRobot.from_dir_glob("flow_dir/w*/outdata/")

        Args:
            pattern: Pattern string
	    walk: if True, directories inside `top` are included as well.
            abspath: True if paths in index should be absolute. Default: Relative to getcwd().
        """
        import glob
        items = []
        for top in filter(os.path.isdir, glob.iglob(pattern)):
            items += cls._open_files_in_dir(top, walk=walk)
        new = cls(*items)
        if not abspath: new.trim_paths(start=os.getcwd())
        return new

    @classmethod
    def _open_files_in_dir(cls, top, walk):
        """Open files in directory tree starting from `top`. Return list of Abinit files."""
        if not os.path.isdir(top):
            raise ValueError("%s: no such directory" % str(top))
        from abipy.abilab import abiopen
        items = []
        if walk:
            for dirpath, dirnames, filenames in os.walk(top):
                filenames = [f for f in filenames if cls.class_handles_filename(f)]
                for f in filenames:
                    abifile = abiopen(os.path.join(dirpath, f))
                    if abifile is not None: items.append((abifile.filepath, abifile))
        else:
            filenames = [f for f in os.listdir(top) if cls.class_handles_filename(f)]
            for f in filenames:
                abifile = abiopen(os.path.join(top, f))
                if abifile is not None: items.append((abifile.filepath, abifile))

        return items

    @classmethod
    def class_handles_filename(cls, filename):
        """True if robot class handles filename."""
        return (filename.endswith("_" + cls.EXT + ".nc") or
                filename.endswith("." + cls.EXT))  # This for .abo

    @classmethod
    def from_files(cls, filenames, labels=None):
        """
        Build a Robot from a list of `filenames`.
        if labels is None, labels are automatically generated from absolute paths.
        """
        filenames = list_strings(filenames)
        from abipy.abilab import abiopen
        filenames = [f for f in filenames if cls.class_handles_filename(f)]
        items = []
        for i, f in enumerate(filenames):
            try:
                abifile = abiopen(f)
            except Exception:
                abifile = None

            if abifile is not None:
                label = abifile.filepath if labels is None else labels[i]
                items.append((label, abifile))

        return cls(*items)

    @classmethod
    def from_flow(cls, flow, outdirs="all", nids=None, ext=None, task_class=None):
        """
        Build a robot from a |Flow| object.

        Args:
            flow: |Flow| object
            outdirs: String used to select/ignore the files in the output directory of flow, works and tasks
                outdirs="work" selects only the outdir of the Works,
                outdirs="flow+task" selects the outdir of the Flow and the outdirs of the tasks
                outdirs="-work" excludes the outdir of the Works.
                Cannot use ``+`` and ``-`` flags in the same string.
                Default: `all` that is equivalent to "flow+work+task"
            nids: List of node identifiers used to select particular nodes. Not used if None
            ext: File extension associated to the robot. Mainly used if method is invoked with the BaseClass
            task_class: Task class or string with the class name used to select the tasks in the flow.
                None implies no filtering.

        Usage example:

        .. code-block:: python

            with abilab.GsrRobot.from_flow(flow) as robot:
                print(robot)

            # That is equivalent to:
            with Robot.from_flow(flow, ext="GSR") as robot:
                print(robot)

        Returns:
            ``Robot`` subclass.
        """
        robot = cls() if ext is None else cls.class_for_ext(ext)()
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
            robot.add_extfile_of_node(flow, nids=nids, task_class=task_class)

        if "work" in tokens:
            for work in flow:
                robot.add_extfile_of_node(work, nids=nids, task_class=task_class)

        if "task" in tokens:
            for task in flow.iflat_tasks():
                robot.add_extfile_of_node(task, nids=nids, task_class=task_class)

        return robot

    def add_extfile_of_node(self, node, nids=None, task_class=None):
        """
        Add the file produced by this node to the robot.

        Args:
            node: |Flow| or |Work| or |Task| object.
            nids: List of node identifiers used to select particular nodes. Not used if None
            task_class: Task class or string with class name used to select the tasks in the flow.
                None implies no filtering.
        """
        if nids and node.node_id not in nids: return
        filepath = node.outdir.has_abiext(self.EXT)
        if not filepath:
            # Look in run.abi directory.
            filepath = node.wdir.has_abiext(self.EXT)

        if filepath:
            try:
                label = os.path.relpath(filepath)
            except OSError:
                # current working directory may not be defined!
                label = filepath

            # Filter by task_class (class or string with class name)
            if task_class is not None and not node.isinstance(task_class):
                return None

            self.add_file(label, filepath)

    def scan_dir(self, top, walk=True):
        """
        Scan directory tree starting from ``top``. Add files to the robot instance.

        Args:
            top (str): Root directory
            walk: if True, directories inside ``top`` are included as well.

        Return:
            Number of files found.
	"""
        count = 0
        for filepath, abifile in self.__class__._open_files_in_dir(top, walk):
            count += 1
            self.add_file(filepath, abifile)

        return count

    def add_file(self, label, abifile, filter_abifile=None):
        """
        Add a file to the robot with the given label.

        Args:
            label: String used to identify the file (must be unique, ax exceptions is
                raised if label is already present.
            abifile: Specify the file to be added. Accepts strings (filepath) or abipy file-like objects.
            filter_abifile: Function that receives an ``abifile`` object and returns
                True if the file should be added to the plotter.
        """
        if is_string(abifile):
            from abipy.abilab import abiopen
            abifile = abiopen(abifile)
            if filter_abifile is not None and not filter_abifile(abifile):
                abifile.close()
                return

            # Open file here --> have to close it.
            self._do_close[abifile.filepath] = True

        if label in self._abifiles:
            raise ValueError("label %s is already present!" % label)

        self._abifiles[label] = abifile

    #def pop_filepath(self, filepath):
    #    """
    #    Remove the file with the given `filepath` and close it.
    #    """
    #    if label, abifile in self._abifiles.items():
    #        if abifile.filepath != filepath: continue
    #        self._abifiles.pop(label)
    #        if self._do_close.pop(abifile.filepath, False):
    #            try:
    #                abifile.close()
    #            except Exception as exc:
    #                print("Exception while closing: ", abifile.filepath)
    #                print(exc)
    #                #raise

    @staticmethod
    def ordered_intersection(list_1, list_2):
        """Return ordered intersection of two lists. Items must be hashable."""
        set_2 = frozenset(list_2)
        return [x for x in list_1 if x in set_2]

    #def _get_ointersection_i(self, iattrname):
    #    if len(self.abifiles) == 0: return []
    #    values = list(range(getattr(self.abifiles[0], iattrname)))
    #    if len(self.abifiles) == 1: return values
    #    for abifile in self.abifiles[1:]:
    #        values = self.ordered_intersection(values, range(getattr(abifile, iattrname)))
    #    return values

    @staticmethod
    def _to_relpaths(paths):
        """Convert a list of absolute paths to relative paths."""
        root = os.getcwd()
        return [os.path.relpath(p, root) for p in paths]

    def pop_label(self, label):
        """
        Remove file with the given ``label`` and close it.
        """
        if label in self._abifiles:
            abifile = self._abifiles.pop(label)
            if self._do_close.pop(abifile.filepath, False):
                try:
                    abifile.close()
                except Exception as exc:
                    print("Exception while closing: ", abifile.filepath)
                    print(exc)
                    #raise

    def change_labels(self, new_labels, dryrun=False):
        """
        Change labels of the files. Return mapping new --> old.
        """
        if len(new_labels) != len(self):
            raise ValueError("Robot has %d files while len(new_labels) = %d" % (len(new_labels), len(self)))

        old_labels = list(self._abifiles.keys())
        old_abifiles = self._abifiles
        self._abifiles = OrderedDict()
        new2old = OrderedDict()
        for old, new in zip(old_labels, new_labels):
            print("old [%s] --> new [%s]" % (old, new))
            new2old[new] = old
            if not dryrun:
                self._abifiles[new] = old_abifiles[old]

        return new2old

    def trim_paths(self, start=None):
        """
        Replace absolute filepaths in the robot with relative paths wrt to ``start`` directory.
        If start is None, os.getcwd() is used. Set ``self.start`` attribute, return ``self.start``.
        """
        self.start = os.getcwd() if start is None else start
        old_paths = list(self._abifiles.keys())
        old_new_paths = [(p, os.path.relpath(os.path.abspath(p), start=self.start)) for p in old_paths]

        old_abifiles = self._abifiles
        self._abifiles = OrderedDict()
        for old, new in old_new_paths:
            self._abifiles[new] = old_abifiles[old]

        return self.start

    @property
    def exceptions(self):
        """List of exceptions."""
        return self._exceptions

    def __len__(self):
        return len(self._abifiles)

    def __iter__(self):
        return iter(self._abifiles.items())

    def __getitem__(self, key):
        # self[key]
        return self._abifiles.__getitem__(key)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Activated at the end of the with statement."""
        self.close()

    def items(self):
        return self._abifiles.items()

    @property
    def labels(self):
        return list(self._abifiles.keys())

    def get_label_files_str(self):
        """Return string with [label, filepath]."""
        from tabulate import tabulate
        return tabulate([(label, abifile.relpath) for label, abifile in self], headers=["Label", "Relpath"]) + "\n"

    def show_files(self, stream=sys.stdout):
        """Show label --> file path"""
        stream.write(self.get_label_files_str())

    def __str__(self):
        """Invoked by str."""
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        lines = ["%s with %d files in memory:\n" % (self.__class__.__name__, len(self.abifiles))]
        app = lines.append
        for i, f in enumerate(self.abifiles):
            app(f.to_string(verbose=verbose))
            app("\n")

        return "\n".join(lines)

    def _repr_html_(self):
        """Integration with jupyter_ notebooks."""
        return "<ol>\n{}\n</ol>".format("\n".join("<li>%s</li>" % label for label, abifile in self))

    @property
    def abifiles(self):
        """List of netcdf files."""
        return list(self._abifiles.values())

    def is_sortable(self, aname, raise_exc=False):
        """
        Return True if ``aname`` is an attribute of the netcdf file
        If raise_exc is True, AttributeError with an explicit message is raised.
        """
        try:
            obj = getattr(self.abifiles[0], aname)
            #return not (callable(obj) or hasattr(obj, "__len__"))
            return True
        except AttributeError:
            if not raise_exc: return False
            attrs = []
            for key, obj in inspect.getmembers(self.abifiles[0]):
                # Ignores anything starting with underscore
                if key.startswith('_') or callable(obj) or hasattr(obj, "__len__"): continue
                attrs.append(key)

            raise AttributeError("""\
`%s` object has no attribute `%s`. Choose among:

    %s

Note that this list is automatically generated.
Not all entries are sortable (Please select number-like quantities)""" % (self.__class__.__name__, aname, str(attrs)))

    def sortby(self, func_or_string, reverse=False):
        """
        Sort files in the robot by ``func_or_string``.
        Return list of (label, abifile, param) tuples where param is obtained via ``func_or_string``.

        Args:
            func_or_string: Either None, string, callable defining the quantity to be used for sorting.
                If string, it's assumed that the abifile has an attribute with the same name and getattr is invoked.
                If callable, the output of callable(abifile) is used.
                If None, no sorting is performed.
            reverse: If set to True, then the list elements are sorted as if each comparison were reversed.
        """
        if not func_or_string:
            return [(label, abifile, label) for (label, abifile) in self]
        elif callable(func_or_string):
            items = [(label, abifile, func_or_string(abifile)) for (label, abifile) in self]
        else:
            # Assume string and attribute with the same name.
            self.is_sortable(func_or_string, raise_exc=True)
            items = [(label, abifile, getattr(abifile, func_or_string)) for (label, abifile) in self]

        return sorted(items, key=lambda t: t[2], reverse=reverse)

    def close(self):
        """
        Close all files that have been opened by the Robot.
        """
        for abifile in self.abifiles:
            if self._do_close.pop(abifile.filepath, False):
                try:
                    abifile.close()
                except:
                    print("Exception while closing: ", abifile.filepath)
                    print(exc)

    #@classmethod
    #def open(cls, obj, nids=None, **kwargs):
    #    """
    #    Flexible constructor. obj can be a :class:`Flow` or a string with the directory containing the Flow.
    #    `nids` is an optional list of :class:`Node` identifiers used to filter the set of :class:`Task` in the Flow.
    #    """
    #    has_dirpath = False
    #    if is_string(obj):
    #        try:
    #            from abipy.flowtk import Flow
    #            obj = Flow.pickle_load(obj)
    #        except:
    #            has_dirpath = True

    #    if not has_dirpath:
    #        # We have a Flow. smeth is the name of the Task method used to open the file.
    #        items = []
    #        smeth = "open_" + cls.EXT.lower()
    #        for task in obj.iflat_tasks(nids=nids): #, status=obj.S_OK):
    #            open_method = getattr(task, smeth, None)
    #            if open_method is None: continue
    #            abifile = open_method()
    #            if abifile is not None: items.append((task.pos_str, abifile))
    #        return cls(*items)

    #    else:
    #        # directory --> search for files with the appropriate extension and open it with abiopen.
    #        if nids is not None: raise ValueError("nids cannot be used when obj is a directory.")
    #        return cls.from_dir(obj)

    #def get_attributes(self, attr_name, obj=None, retdict=False):
    #    od = OrderedDict()
    #    for label, abifile in self:
    #        obj = abifile if obj is None else getattr(abifile, obj)
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
                cprint("Exception: %s" % str(exc), "red")
                self._exceptions.append(str(exc))
        return d

    @staticmethod
    def sortby_label(sortby, param):
        """Return the label to be used when files are sorted with ``sortby``."""
        return "%s %s" % (sortby, param) if not callable(sortby) else str(param)

    def get_structure_dataframes(self, abspath=False, filter_abifile=None, **kwargs):
        """
        Wrap dataframes_from_structures function.

        Args:
            abspath: True if paths in index should be absolute. Default: Relative to getcwd().
            filter_abifile: Function that receives an ``abifile`` object and returns
                True if the file should be added to the plotter.
        """
        from abipy.core.structure import dataframes_from_structures
        if "index" not in kwargs:
            index = list(self._abifiles.keys())
            if not abspath: index = self._to_relpaths(index)
            kwargs["index"] = index

        abifiles = self.abifiles if filter_abifile is not None else list(filter(filter_abifile, self.abifiles))
        return dataframes_from_structures(struct_objects=abifiles, **kwargs)

    def get_lattice_dataframe(self, **kwargs):
        """Return |pandas-DataFrame| with lattice parameters."""
        dfs = self.get_structure_dataframes(**kwargs)
        return dfs.lattice

    def get_coords_dataframe(self, **kwargs):
        """Return |pandas-DataFrame| with atomic positions."""
        dfs = self.get_structure_dataframes(**kwargs)
        return dfs.coords

    ##############################################
    # Helper functions to plot pandas dataframes #
    ##############################################

    @staticmethod
    @wraps(plot_xy_with_hue)
    def plot_xy_with_hue(*args, **kwargs):
        return plot_xy_with_hue(*args, **kwargs)

    def get_baserobot_code_cells(self, title=None):
        """
        Return list of notebook_ cells with calls to methods provides by the baseclass.
        """
        # Try not pollute namespace with lots of variables.
        nbformat, nbv = self.get_nbformat_nbv()
        title = "## Code to compare multiple Structure objects" if title is None else str(title)
        return [
            nbv.new_markdown_cell(title),
            nbv.new_code_cell("robot.get_lattice_dataframe()"),
            nbv.new_code_cell("#robot.get_coords_dataframe()"),
        ]