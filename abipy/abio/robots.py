# coding: utf-8
"""
This module defines the Robot BaseClass. Robots operates on multiple files and provide helper
functions to plot the data e.g. convergence studies and to build pandas dataframes from the output files.
"""
from __future__ import annotations

import sys
import os
import inspect
import itertools
import json
import numpy as np
import pandas as pd

from collections import OrderedDict, deque
from typing import Callable, Union, Any
from functools import wraps
from monty.string import is_string, list_strings
from monty.termcolor import cprint
from monty.json import MontyEncoder
from abipy.tools.serialization import pmg_serialize
from abipy.tools.iotools import make_executable
from abipy.core.structure import Structure
from abipy.core.mixins import NotebookWriter
from abipy.tools.numtools import sort_and_groupby
from abipy.tools import duck
from abipy.tools.typing import Figure
from abipy.tools.plotting import (plot_xy_with_hue, add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt,
    rotate_ticklabels, set_visible, ConvergenceAnalyzer)


class Robot(NotebookWriter):
    """
    This is the base class from which all Robot subclasses should derive.
    A Robot supports the `with` context manager:

    Usage example:

    .. code-block:: python

        with Robot([("label1", "file1"), (label2, "file2")]) as robot:
            # Do something with robot. files are automatically closed when we exit.
            for label, abifile in self.items():
                print(label)
    """
    # filepaths are relative to `start`. None for asbolute paths. This flag is set in trim_paths
    start = None

    # Used in iter_lineopt to generate matplotlib linestyles.
    _LINE_COLORS = ["b", "r", "g", "m", "y", "k", "c"]
    _LINE_STYLES = ["-", ":", "--", "-.",]
    _LINE_WIDTHS = [2, ]

    # matplotlib option to fill convergence window.
    HATCH = "/"

    def __init__(self, *args):
        """
        Args:
            args is a list of tuples (label, filepath)
        """
        self._abifiles, self._do_close = OrderedDict(), OrderedDict()
        self._exceptions = deque(maxlen=100)

        for label, abifile in args:
            self.add_file(label, abifile)

    @classmethod
    def get_supported_extensions(self) -> list[str]:
        """List of strings with extensions supported by Robot subclasses."""
        # This is needed to have all subclasses.
        from abipy.abilab import Robot
        return sorted([cls.EXT for cls in Robot.__subclasses__()])

    @classmethod
    def class_for_ext(cls, ext: str):
        """Return the Robot subclass associated to the given extension."""
        for subcls in cls.__subclasses__():
            if subcls.EXT in (ext, ext.upper()):
                return subcls

        # anaddb.nc does not follow the extension rule...
        if ext.lower() == "anaddb":
            from abipy.dfpt.anaddbnc import AnaddbNcRobot as subcls
            return subcls

        raise ValueError("Cannot find Robot subclass associated to extension %s\n" % ext +
                         "The list of supported extensions (case insensitive) is:\n%s" %
                         str(cls.get_supported_extensions()))

    @classmethod
    def from_dir(cls, top: str, walk: bool = True, abspath: bool = False) -> Robot:
        """
        Build a robot by scanning all files located within directory `top`.
        This method should be invoked with a concrete robot class, for example:

            robot = GsrRobot.from_dir(".")

        Args:
            top: Root directory
            walk: if True, directories inside `top` are included as well.
            abspath: True if paths in index should be absolute. Default: Relative to `top`.
        """
        new = cls(*cls._open_files_in_dir(top, walk))
        if not abspath: new.trim_paths(start=top)
        return new

    @classmethod
    def from_dirs(cls, dirpaths: list[str], walk: bool = True, abspath: bool = False) -> Robot:
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
    def from_dir_glob(cls, pattern: str, walk: bool = True, abspath: bool = False) -> Robot:
        """
        This class method builds a robot by scanning all files located within the directories
        matching `pattern` as implemented by glob.glob
        This method should be invoked with a concrete robot class, for example:

            robot = GsrRobot.from_dir_glob("flow_dir/w*/outdata/")

        Args:
            pattern: Pattern string.
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
    def _open_files_in_dir(cls, top: str, walk: bool) -> list:
        """
        Open files in directory tree starting from `top`. Return list of Abinit files.
        """
        if not os.path.isdir(top):
            raise ValueError("%s: no such directory" % str(top))
        from abipy.abilab import abiopen
        items = []
        if walk:
            for dirpath, dirnames, filenames in os.walk(top):
                filenames = sorted([f for f in filenames if cls.class_handles_filename(f)])
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
    def class_handles_filename(cls, filename: str) -> bool:
        """
        True if robot class handles filename.
        """
        # Special treatment of AnaddbNcRobot
        if cls.EXT == "anaddb" and os.path.basename(filename).lower() == "anaddb.nc":
            return True

        return (filename.endswith("_" + cls.EXT + ".nc") or
                filename.endswith("." + cls.EXT))  # This for .abo

    @classmethod
    def from_files(cls, filenames, labels=None, abspath=False) -> Robot:
        """
        Build a Robot from a list of `filenames`.
        If labels is None, labels are automatically generated from absolute paths.

        Args:
            abspath: True if paths in index should be absolute. Default: Relative to `top`.
        """
        filenames = list_strings(filenames)
        from abipy.abilab import abiopen
        filenames = [f for f in filenames if cls.class_handles_filename(f)]
        items = []
        for i, f in enumerate(filenames):
            try:
                abifile = abiopen(f)
            except Exception as exc:
                cprint("Exception while opening file: `%s`" % str(f), "red")
                cprint(exc, "red")
                abifile = None

            if abifile is not None:
                label = abifile.filepath if labels is None else labels[i]
                items.append((label, abifile))

        new = cls(*items)
        if labels is None and not abspath: new.trim_paths(start=None)
        return new

    @classmethod
    def from_json_file(cls, filepath: str):
        """
        Build Robot from a file in json format.
        """
        from abipy.tools.serialization import mjson_load
        new = mjson_load(filepath)
        return new

    @classmethod
    def from_top_and_json_basename(cls, top: str, json_basename: str) -> Robot:
        """
        Build a robot by scanning all json files located within directory `top`.
        and matching json_basename.

        Example:

            gwr_robot = Robot.from_top_and_json_basename(".", "gwr_robot.json")
        """
        json_paths = []
        for root, dirs, files in os.walk(top):
            for name in files:
                if name == json_basename:
                    json_paths.append(os.path.join(root, name))

        return cls.from_json_files(json_paths)

    @classmethod
    def from_json_files(cls, json_paths: list[str]):
        """
        Build Robot from a list of json files.
        Each json file should have a list of filepaths and @module and @class as required by msonable.
        """
        # Merge filepaths and chech that module and class are equal across files.
        robot_filepaths = []
        _module, _class = None, None
        for path in json_paths:
            with open(path, "rt") as fh:
                d = json.load(fh)

            if "filepaths" not in d:
                raise ValueError(f"{path=} does not provide `filepaths` entry")

            # Merge filepaths
            robot_filepaths.extend(d["filepaths"])

            if _module is None:
                _module, _class = d["@module"], d["@class"]
            else:
                if _module != d["@module"]:
                    raise ValueError(f'{_module=} != {d["@module"]=}')
                if _class != d["@class"]:
                    raise ValueError(f'{_class=} != {d["@class"]=}')

        # Build Robot instance from string in json format.
        #print("Merged filepaths:", robot_filepaths)
        from abipy.tools.serialization import mjson_loads
        data = {"filepaths": robot_filepaths, "@class": _class, "@module": _module}
        return mjson_loads(json.dumps(data))

    @classmethod
    def from_work(cls, work, outdirs="all", nids=None, ext=None, task_class=None) -> Robot:
        """
        Build a robot from a |Work| object.
        """
        nids = [work.node_id] + [task.node_id for task in work]
        return cls.from_flow(work.flow, nids=nids)

    @classmethod
    def from_flow(cls, flow, outdirs="all", nids=None, ext=None, task_class=None) -> Robot:
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


    def __len__(self):
        return len(self._abifiles)

    #def __iter__(self):
    #    return iter(self._abifiles)

    def __getitem__(self, key):
        # self[key]
        return self._abifiles.__getitem__(key)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Activated at the end of the with statement."""
        self.close()

    def keys(self):
        return self._abifiles.keys()

    def items(self):
        return self._abifiles.items()

    def add_extfile_of_node(self, node, nids=None, task_class=None) -> None:
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

        # This to ignore DDB.nc files (only text DDB are supported)
        if filepath and filepath.endswith("_DDB.nc"):
            return

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

    def scan_dir(self, top: str, walk: bool = True) -> int:
        """
        Scan directory tree starting from ``top``.
        Add files to the robot instance.
        Return: Number of files found.

        Args:
            top: Root directory
            walk: if True, directories inside ``top`` are included as well.
        """
        count = 0
        for filepath, abifile in self.__class__._open_files_in_dir(top, walk):
            count += 1
            self.add_file(filepath, abifile)

        return count

    def add_file(self, label, abifile, filter_abifile=None) -> None:
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

    @pmg_serialize
    def as_dict(self) -> dict:
        """Return dict with filepaths that can be used to reconstruct the Robot."""
        return dict(filepaths=[abifile.filepath for abifile in self.abifiles])

    @classmethod
    def from_dict(cls, d: dict):
        """Recontruct object from dictionary with filepaths."""
        return cls.from_files(d["filepaths"])

    def to_json(self) -> str:
        """
        Returns a JSON string representation of the object.
        """
        return json.dumps(self.as_dict(), cls=MontyEncoder)

    def get_pyscript(self, filepath: str) -> RobotPythonScript:
        """Return RobotPythonScript to br used as context manager."""
        return RobotPythonScript(self, filepath)

    #def pop_filepath(self, filepath: str) -> None:
    #    """
    #    Remove the file with the given `filepath`.
    #    """
    #    if label, abifile in self._abifiles.items():
    #        if abifile.filepath != filepath: continue
    #        self._abifiles.pop(label)

    def iter_lineopt(self):
        """Generates matplotlib linestyles."""
        for o in itertools.product( self._LINE_WIDTHS,  self._LINE_STYLES, self._LINE_COLORS):
            yield {"linewidth": o[0], "linestyle": o[1], "color": o[2]}

    @staticmethod
    def ordered_intersection(list_1, list_2) -> list:
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
    def _to_relpaths(paths) -> list[str]:
        """Convert a list of absolute paths to relative paths."""
        root = os.getcwd()
        return [os.path.relpath(p, root) for p in paths]

    def remove(self) -> None:
        """Close the file handle, remove the file from disk for each file in the robot."""
        for abifile in self.abifiles:
            abifile.remove()

    def pop_label(self, label: str) -> None:
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

    def change_labels(self, new_labels: list[str], dryrun: bool = False) -> dict:
        """
        Change labels of the files.

        Args:
            new_labels: List of strings (same length as self.abifiles)
            dryrun: True to activate dryrun mode.

        Return:
            mapping new_label --> old_label.
        """
        if len(new_labels) != len(self):
            raise ValueError("Robot has %d files while len(new_labels) = %d" % (len(new_labels), len(self)))

        old_labels = list(self._abifiles.keys())
        if not dryrun:
            old_abifiles, self._abifiles = self._abifiles, OrderedDict()
        new2old = OrderedDict()
        for old, new in zip(old_labels, new_labels):
            new2old[new] = old
            if not dryrun:
                self._abifiles[new] = old_abifiles[old]
            else:
                print("old [%s] --> new [%s]" % (old, new))

        return new2old

    def remap_labels(self, function: Callable, dryrun: bool = False) -> dict:
        """
        Change labels of the files by executing ``function``

        Args:
            function: Callable object e.g. lambda function. The output of function(abifile) is used as
                new label. Note that the function shall not return duplicated labels when applied to self.abifiles.
            dryrun: True to activate dryrun mode.

        Return:
            mapping new_label --> old_label.
        """
        new_labels = [function(afile) for afile in self.abifiles]
        # Labels must be unique and hashable.
        if len(set(new_labels)) != len(new_labels):
            raise ValueError("Duplicated labels are not allowed. Change input function.\nnew_labels %s" % str(new_labels))

        return self.change_labels(new_labels, dryrun=dryrun)

    def trim_paths(self, start=None) -> str:
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
    def exceptions(self) -> list:
        """List of exceptions."""
        return self._exceptions

    @property
    def labels(self) -> list[str]:
        """
        List of strings used to create labels in matplotlib figures when plotting results
        taked from multiple files. By default, labels is initialized with the path of the files in the robot.
        Use change_labels to change the list.
        """
        return list(self._abifiles.keys())

    def get_label_files_str(self):
        """Return string with [label, filepath]."""
        from tabulate import tabulate
        return tabulate([(label, abifile.relpath) for label, abifile in self.items()],
                        headers=["Label", "Relpath"]) + "\n"

    def show_files(self, stream=sys.stdout) -> None:
        """Show label --> file path"""
        stream.write(self.get_label_files_str())

    def __repr__(self):
        """Invoked by repr."""
        return self.get_label_files_str()

    def __str__(self):
        """Invoked by str."""
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        """String representation."""
        lines = ["%s with %d files in memory:\n" % (self.__class__.__name__, len(self.abifiles))]
        app = lines.append
        for i, f in enumerate(self.abifiles):
            app(f.to_string(verbose=verbose))
            app("\n")

        return "\n".join(lines)

    def _repr_html_(self) -> str:
        """Integration with jupyter_ notebooks."""
        return '<ol start="0">\n{}\n</ol>'.format("\n".join("<li>%s</li>" % label for label, abifile in self.items()))

    def getattr_alleq(self, aname : str):
        """
        Return the value of attribute aname.
        Raises ValueError if value is not the same across all the files in the robot.
        """
        val1 = getattr(self.abifiles[0], aname)

        for abifile in self.abifiles[1:]:
            val2 = getattr(abifile, aname)
            if isinstance(val1, (str, int, float)):
                eq = val1 == val2
            elif isinstance(val1, np.ndarray):
                eq = np.allclose(val1, val2)
            if not eq:
                raise ValueError(f"Different values of {aname=}, {val1=}, {val2=}")

        return val1

    @property
    def abifiles(self) -> list:
        """List of netcdf files."""
        return list(self._abifiles.values())

    def has_different_structures(self, rtol=1e-05, atol=1e-08) -> str:
        """
        Check if structures are equivalent, return string with info about differences (if any).
        """
        if len(self) <= 1: return ""
        formulas = set([af.structure.composition.formula for af in self.abifiles])
        if len(formulas) != 1:
            return "Found structures with different full formulas: %s" % str(formulas)

        lines = []
        s0 = self.abifiles[0].structure
        for abifile in self.abifiles[1:]:
            s1 = abifile.structure
            if not np.allclose(s0.lattice.matrix, s1.lattice.matrix, rtol=rtol, atol=atol):
                lines.append("Structures have different lattice:")
            if not np.allclose(s0.frac_coords, s1.frac_coords, rtol=rtol, atol=atol):
                lines.append("Structures have different atomic positions:")

        return "\n".join(lines)

    def _get_ref_abifile_from_basename(self, ref_basename: str | None):
        """
        Find reference abifile. If None, the first file in the robot is used.
        """
        ref_file = self.abifiles[0]
        if ref_basename is None:
            return ref_file

        for i, abifile in enumerate(self.abifiles):
            if abifile.basename == ref_basename:
                return abifile

        raise ValueError(f"Cannot find {ref_basename=}")

    @staticmethod
    def _compare_attr_name(aname: str, ref_abifile, other_abifile) -> None:
        """
        Compare the value of attribute `aname` in two files.
        """
        # Get attributes in abifile first, then in abifile.r, else raise.
        if hasattr(ref_abifile, aname):
            val1, val2 = getattr(ref_abifile, aname), getattr(other_abifile, aname)

        elif hasattr(ref_abifile , "r") and hasattr(ref_abifile.r, aname):
            val1, val2 = getattr(ref_abifile.r, aname), getattr(other_abifile.r, aname)

        else:
            raise AttributeError(f"Cannot find attribute `{aname =}`")

        # Now compare val1 and val2 taking into account the type.
        if isinstance(val1, (str, int, float, Structure)):
            eq = val1 == val2
        elif isinstance(val1, np.ndarray):
            eq = np.allclose(val1, val2)
        else:
            raise TypeError(f"Don't know how to handle comparison for type: {type(val1)}")

        if not eq:
            raise ValueError(f"Different values of {aname=}, {val1=}, {val2=}")

    def is_sortable(self, aname: str, raise_exc: bool = False) -> bool:
        """
        Return True if ``aname`` is an attribute of the netcdf file
        If raise_exc is True, AttributeError with an explicit message is raised.
        """
        try:
            obj = None
            try:
                # abifile.foo.bar ?
                obj = duck.getattrd(self.abifiles[0], aname)
            except AttributeError:
                # abifile.params[aname] ?
                if hasattr(self.abifiles[0], "params") and aname in self.abifiles[0].params:
                    obj = self.abifiles[0].params[aname]

            # Let's try to convert obj to scalar.
            float(obj)
            return True

        except Exception:
            if not raise_exc: return False
            attrs = []
            for key, obj in inspect.getmembers(self.abifiles[0]):
                # Ignores anything starting with underscore
                if key.startswith('_') or callable(obj) or hasattr(obj, "__len__"): continue
                attrs.append(key)

            # Add all keys in params dict.
            if hasattr(self.abifiles[0], "params") and hasattr(self.abifiles[0].params, "keys"):
                attrs.extend(self.abifiles[0].params.keys())

            raise AttributeError("""\
`%s` object has no attribute `%s`. Choose among:

    %s

Note that this list is automatically generated.
Not all entries are sortable (Please select number-like quantities)""" % (self.__class__.__name__, aname, str(attrs)))

    def _sortby_labelfile_list(self, labelfile_list, func_or_string, reverse=False, unpack=False):
        """
        Return: list of (label, abifile, xs) tuples where xs is obtained via ``func_or_string``.
            or labels, abifiles, xs if ``unpack``
        """
        if not func_or_string:
            # Catch None or empty
            items = [(label, abifile, label) for (label, abifile) in labelfile_list]
            if not unpack:
                return items
            else:
                return [t[0] for t in items], [t[1] for t in items], [t[2] for t in items]

        elif callable(func_or_string):
            items = [(label, abifile, func_or_string(abifile)) for (label, abifile) in labelfile_list]

        else:
            # Assume string and attribute with the same name.
            # try in abifile.params if not hasattrd(abifile, func_or_string)
            self.is_sortable(func_or_string, raise_exc=True)
            if duck.hasattrd(self.abifiles[0], func_or_string):
                items = [(label, abifile, duck.getattrd(abifile, func_or_string)) for (label, abifile) in labelfile_list]

            #elif hasattr(self.abifiles[0], "reader") and duck.hasattrd(self.abifiles[0].reader, func_or_string):
            #    items = [(label, abifile, duck.getattrd(abifile.reader, func_or_string)) for (label, abifile) in labelfile_list]

            #elif hasattr(self.abifiles[0], "r") and duck.hasattrd(self.abifiles[0].r, func_or_string):
            #elif hasattr(self.abifiles[0], "r"):
            #    items = [(label, abifile, duck.getattrd(abifile.r, func_or_string)) for (label, abifile) in labelfile_list]

            else:
                #print("Cannot find", func_or_string, " in self.abifiles[0]")
                items = [(label, abifile, abifile.params[func_or_string]) for (label, abifile) in labelfile_list]

        items = sorted(items, key=lambda t: t[2], reverse=reverse)
        if not unpack:
            return items
        else:
            return [t[0] for t in items], [t[1] for t in items], [t[2] for t in items]

    def sortby(self, func_or_string: Union[Callable, str, None],
               reverse: bool = False, unpack: bool = False) -> list[tuple]:
        """
        Sort files in the robot by ``func_or_string``.

        Args:
            func_or_string: Can be None, string or callable defining the quantity to be used for sorting.
                If string, it's assumed that the abifile has an attribute with the same name and getattr is invoked.
                If callable, the output of func_or_string(abifile) is used.
                If None, no sorting is performed.
            reverse: If set to True, then the list elements are sorted as if each comparison were reversed.
            unpack: Return (labels, abifiles, xs) if True

        Return: list of (label, abifile, xs) tuples where xs is obtained via ``func_or_string``.
            or labels, abifiles, xs if ``unpack``
        """
        labelfile_list = list(self.items())
        return self._sortby_labelfile_list(labelfile_list, func_or_string, reverse=reverse, unpack=unpack)

    def group_and_sortby(self,
                         hue: Union[Callable, str],
                         func_or_string: Union[Callable, str, None]) -> list[HueGroup]:
        """
        Group files by ``hue`` and, inside each group` sort items by ``func_or_string``.

        Args:
            hue: Variable that defines subsets of the data, which will be drawn on separate lines.
                Accepts callable or string
                If string, it's assumed that the abifile has an attribute with the same name and getattr is invoked.
                Dot notation is also supported e.g. hue="structure.formula" --> abifile.structure.formula
                If callable, the output of hue(abifile) is used.
            func_or_string: Either None, string, callable defining the quantity to be used for sorting.
                If string, it's assumed that the abifile has an attribute with the same name and getattr is invoked.
                If callable, the output of func_or_string(abifile) is used.
                If None, no sorting is performed.

        Return: List of :class:`HueGroup` instance.
        """
        # Group by hue.
        # This is the section in which we support: callable, abifile.attr.name syntax or abifile.params["key"]
        items = list(self.items())

        if callable(hue):
            key = lambda t: hue(t[1])
        else:
            # Assume string.
            if duck.hasattrd(self.abifiles[0], hue):
                key = lambda t: duck.getattrd(t[1], hue)
            else:
                # Try in abifile.params
                if hasattr(self.abifiles[0], "params") and hue in self.abifiles[0].params:
                    key = lambda t: t[1].params[hue]
                else:
                    raise TypeError("""\
Cannot interpret hue argument of type `%s` and value `%s`.
Expecting callable or attribute name or key in abifile.params""" % (type(hue), str(hue)))

        groups = []
        for hvalue, labelfile_list in sort_and_groupby(items, key=key):
            # Use func_or_string to sort each group
            labels, abifiles, xvalues = self._sortby_labelfile_list(labelfile_list, func_or_string, unpack=True)
            groups.append(HueGroup(hvalue, xvalues, abifiles, labels))

        return groups

    def close(self) -> None:
        """
        Close all files that have been opened by the Robot.
        """
        for abifile in self.abifiles:
            if self._do_close.pop(abifile.filepath, False):
                try:
                    abifile.close()
                except Exception as exc:
                    print("Exception while closing: ", abifile.filepath)
                    print(exc)

    def _exec_funcs(self, funcs, arg) -> dict:
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
    def sortby_label(sortby, param) -> str:
        """Return the label to be used when files are sorted with ``sortby``."""
        return "%s %s" % (sortby, param) if not (callable(sortby) or sortby is None) else str(param)

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

    def get_lattice_dataframe(self, **kwargs) -> pd.DataFrame:
        """Return |pandas-DataFrame| with lattice parameters."""
        dfs = self.get_structure_dataframes(**kwargs)
        return dfs.lattice

    def get_coords_dataframe(self, **kwargs) -> pd.DataFrame:
        """Return |pandas-DataFrame| with atomic positions."""
        dfs = self.get_structure_dataframes(**kwargs)
        return dfs.coords

    def get_params_dataframe(self, abspath: bool = False) -> pd.DataFrame:
        """
        Return |pandas-DataFrame| with the most important parameters.
        that are usually subject to convergence studies.

        Args:
            abspath: True if paths in index should be absolute. Default: Relative to `top`.
        """
        rows, row_names = [], []
        for label, abifile in self.items():
            if not hasattr(abifile, "params"):
                import warnings
                warnings.warn("%s does not have `params` attribute" % type(abifile))
                break
            rows.append(abifile.params)
            row_names.append(label)

        row_names = row_names if abspath else self._to_relpaths(row_names)
        return pd.DataFrame(rows, index=row_names, columns=list(rows[0].keys()))

    ##############################################
    # Helper functions to plot pandas dataframes #
    ##############################################

    @staticmethod
    @wraps(plot_xy_with_hue)
    def plot_xy_with_hue(*args, **kwargs) -> Figure:
        return plot_xy_with_hue(*args, **kwargs)

    @staticmethod
    def _get_label(func_or_string) -> str:
        """
        Return label associated to ``func_or_string``.
        If callable, docstring __doc__ is used.
        """
        if func_or_string is None:
            return ""
        elif callable(func_or_string):
            if getattr(func_or_string, "__doc__", ""):
                return func_or_string.__doc__.strip()
            else:
                return func_or_string.__name__
        else:
            return str(func_or_string)

    @add_fig_kwargs
    def plot_convergence(self, item: Union[str, Callable],
                         sortby=None, hue=None, abs_conv=None,
                         ax=None, fontsize=8, **kwargs) -> Figure:
        """
        Plot the convergence of ``item`` wrt the ``sortby`` parameter.
        Values can optionally be grouped by ``hue``.

        Args:
            item: Define the quantity to plot. Accepts callable or string
                If string, it's assumed that the abifile has an attribute with the same name and `getattr` is invoked.
                Dot notation is also supported e.g. hue="structure.formula" --> `abifile.structure.formula`.
                If callable, the output of `item(abifile)` is used.
            sortby: Define the convergence parameter, sort files and produce plot labels.
                Can be None, string or function. If None, no sorting is performed.
                If string and not empty it's assumed that the abifile has an attribute
                with the same name and `getattr` is invoked.
                If callable, the output of sortby(abifile) is used.
            hue: Variable that define subsets of the data, which will be drawn on separate lines.
                Accepts callable or string.
                If string, it's assumed that the abifile has an attribute with the same name and getattr is invoked.
                If callable, the output of hue(abifile) is used.
            abs_conv: If not None, plot f(x) and abs(f(x) - f(x_inf)) in log scale.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: legend and label fontsize.
            kwargs: keyword arguments passed to matplotlib plot method.

        Returns: |matplotlib-Figure|

        Example:

             robot.plot_convergence("energy")
             robot.plot_convergence("energy", sortby="nkpt")
             robot.plot_convergence("pressure", sortby="nkpt", hue="tsmear")
             robot.plot_convergence("pressure", sortby="nkpt", hue="tsmear", abs_conv=1e-3)
        """
        if "marker" not in kwargs:
            kwargs["marker"] = "o"

        if abs_conv is None:
            # Plot f(x) only
            ax1, fig, plt = get_ax_fig_plt(ax=ax)
            ax2 = None
        else:
            # Plot f(x) and |f(x) - f(x_inf)| in log scale.
            (ax1, ax2), fig, plt = get_axarray_fig_plt(ax, nrows=2, ncols=1,
                                                       sharex=True, sharey=False, squeeze=True)
            ax2.grid(True)

        if hue is None:
            labels, abifiles, xs = self.sortby(sortby, unpack=True)
            yvals = self.get_yvals_item_abifiles(item, abifiles)
            ax1.plot(xs, yvals, **kwargs)

            if ax2:
                self.plot_abs_conv(ax1, ax2, xs, yvals, abs_conv, self._get_label(sortby),
                                   fontsize, self.HATCH, **kwargs)

        else:
            groups = self.group_and_sortby(hue, sortby)
            for g in groups:
                yvals = self.get_yvals_item_abifiles(item, g.abifiles)
                label = "%s: %s" % (self._get_label(hue), g.hvalue)
                ax1.plot(g.xvalues, yvals, label=label, **kwargs)

                if ax2:
                    self.plot_abs_conv(ax1, ax2, g.xvalues, yvals, abs_conv, self._get_label(sortby),
                                       fontsize, self.HATCH, **kwargs)

        ax1.grid(True)
        ax1.set_xlabel("%s" % self._get_label(sortby))
        if sortby is None: rotate_ticklabels(ax1, 15)
        ax1.set_ylabel("%s" % self._get_label(item))

        if hue is not None:
            ax1.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    @add_fig_kwargs
    def plot_convergence_items(self, items: list[Union[str, Callable]],
                               sortby=None, hue=None, abs_conv=None,
                               fontsize=8, **kwargs) -> Figure:
        """
        Plot the convergence of a list of ``items`` wrt to the ``sortby`` parameter.
        Values can optionally be grouped by ``hue``.

        Args:
            items: List of attributes (or callables) to be analyzed.
            sortby: Define the convergence parameter, sort files and produce plot labels.
                Can be None, string or function. If None, no sorting is performed.
                If string and not empty it's assumed that the abifile has an attribute
                with the same name and `getattr` is invoked.
                If callable, the output of sortby(abifile) is used.
            hue: Variable that define subsets of the data, which will be drawn on separate lines.
                Accepts callable or string
                If string, it's assumed that the abifile has an attribute with the same name and getattr is invoked.
                Dot notation is also supported e.g. hue="structure.formula" --> abifile.structure.formula
                If callable, the output of hue(abifile) is used.
            abs_conv: If not None, plot f(x) and abs(f(x) - f(x_inf)) in log scale.
                Since we are plotting multiple quantities, abs_conv is a dict mapping the name of the item to
                to the convergence.
            fontsize: legend and label fontsize.
            kwargs: keyword arguments are passed to ax.plot

        Returns: |matplotlib-Figure|
        """
        # Note: in principle one could call plot_convergence inside a loop but
        # this one is faster as sorting is done only once.

        # Build grid plot.
        nrows, ncols = len(items), 1 if abs_conv is None else 2
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=False, squeeze=False)

        # Sort and group files if hue.
        if hue is None:
            labels, ncfiles, xs = self.sortby(sortby, unpack=True)
        else:
            groups = self.group_and_sortby(hue, sortby)

        if "marker" not in kwargs:
            kwargs["marker"] = "o"

        for i, (ax_row, item) in enumerate(zip(ax_mat, items)):
            ax1 = ax_row[0]
            ax2 = ax_row[1] if abs_conv is not None else None

            if hue is None:
                # Extract data.
                yvals = self.get_yvals_item_abifiles(item, self.abifiles)

                self.plot_xvals_or_xstr_ax(ax1, xs, yvals, fontsize, **kwargs)

                if ax2:
                    xlabel = self._get_label(sortby) if i == len(items) - 1 else None
                    self.plot_abs_conv(ax1, ax2, xs, yvals, abs_conv[item], xlabel, fontsize, self.HATCH,
                                       **kwargs)

            else:
                for g in groups:
                    # Extract data in group
                    yvals = self.get_yvals_item_abifiles(item, g.abifiles)

                    label = "%s: %s" % (self._get_label(hue), g.hvalue)
                    ax1.plot(g.xvalues, yvals, label=label, **kwargs)

                    if ax2:
                        self.plot_abs_conv(ax1, ax2, g.xvalues, yvals, abs_conv[item], None,
                                       fontsize, self.HATCH, **kwargs)

            ax1.grid(True)
            ax1.set_ylabel(self._get_label(item))
            if ax2:
                ax2.grid(True)

            if i == len(items) - 1:
                ax1.set_xlabel("%s" % self._get_label(sortby))
                if sortby is None: rotate_ticklabels(ax1, 15)

            if i == 0 and hue is not None:
                ax1.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    def get_convergence_analyzer(self, xname: str, ytols_dict: dict) -> ConvergenceAnalyzer:
        """
        The main difference is that ConvergenceAnalyze supports multiple convergence tolerances
        for a given y-value.

        Args:
            xname: Name of the x-variable.
            ytols_dict: dict mapping the name of the y-variable to the tolerance(s).

        Example:
        """
        df = self.get_dataframe()
        return ConvergenceAnalyzer.from_dataframe(df, xname, ytols_dict)

    @add_fig_kwargs
    def plot_lattice_convergence(self, what_list=None, sortby=None, hue=None, fontsize=8, **kwargs) -> Figure:
        """
        Plot the convergence of the lattice parameters (a, b, c, alpha, beta, gamma).
        wrt the``sortby`` parameter. Values can optionally be grouped by ``hue``.

        Args:
            what_list: List of strings with the quantities to plot e.g. ["a", "alpha", "beta"].
                None means all.
            item: Define the quantity to plot. Accepts callable or string
                If string, it's assumed that the abifile has an attribute
                with the same name and `getattr` is invoked.
                If callable, the output of item(abifile) is used.
            sortby: Define the convergence parameter, sort files and produce plot labels.
                Can be None, string or function. If None, no sorting is performed.
                If string and not empty it's assumed that the abifile has an attribute
                with the same name and `getattr` is invoked.
                If callable, the output of sortby(abifile) is used.
            hue: Variable that define subsets of the data, which will be drawn on separate lines.
                Accepts callable or string
                If string, it's assumed that the abifile has an attribute with the same name and getattr is invoked.
                Dot notation is also supported e.g. hue="structure.formula" --> abifile.structure.formula
                If callable, the output of hue(abifile) is used.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: legend and label fontsize.

        Returns: |matplotlib-Figure|

        Example:

             robot.plot_lattice_convergence()

             robot.plot_lattice_convergence(sortby="nkpt")

             robot.plot_lattice_convergence(sortby="nkpt", hue="tsmear")
        """
        if not self.abifiles: return None

        # The majority of AbiPy files have a structure object
        # whereas Hist.nc defines final_structure. Use geattr and key to extract structure object.
        key = "structure"
        if not hasattr(self.abifiles[0], "structure"):
            if hasattr(self.abifiles[0], "final_structure"):
                key = "final_structure"
            else:
                raise TypeError("Don't know how to extract structure from %s" % type(self.abifiles[0]))

        # Define callbacks. docstrings will be used as ylabels.
        def a(afile):
            "a (Ang)"
            return getattr(afile, key).lattice.a
        def b(afile):
            "b (Ang)"
            return getattr(afile, key).lattice.b
        def c(afile):
            "c (Ang)"
            return getattr(afile, key).lattice.c
        def volume(afile):
            r"$V$"
            return getattr(afile, key).lattice.volume
        def alpha(afile):
            r"$\alpha$"
            return getattr(afile, key).lattice.alpha
        def beta(afile):
            r"$\beta$"
            return getattr(afile, key).lattice.beta
        def gamma(afile):
            r"$\gamma$"
            return getattr(afile, key).lattice.gamma

        items = [a, b, c, volume, alpha, beta, gamma]
        if what_list is not None:
            locs = locals()
            items = [locs[what] for what in list_strings(what_list)]

        # Build plot grid.
        nrows, ncols = len(items), 1
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=False, squeeze=False)

        marker = kwargs.pop("marker", "o")
        for i, (ax, item) in enumerate(zip(ax_list.ravel(), items)):
            self.plot_convergence(item, sortby=sortby, hue=hue, ax=ax, fontsize=fontsize,
                                  marker=marker, show=False)
            if i != 0:
                set_visible(ax, False, "legend")
            if i != len(items) - 1:
                set_visible(ax, False, "xlabel")

        return fig

    def get_baserobot_code_cells(self, title=None) -> list:
        """
        Return list of jupyter_ cells with calls to methods provided by the base class.
        """
        # Try not pollute namespace with lots of variables.
        nbformat, nbv = self.get_nbformat_nbv()
        title = "## Code to compare multiple Structure objects" if title is None else str(title)
        return [
            nbv.new_markdown_cell(title),
            nbv.new_code_cell("robot.get_lattice_dataframe()"),
            nbv.new_code_cell("""# robot.plot_lattice_convergence(sortby="nkpt", hue="tsmear")"""),
            nbv.new_code_cell("#robot.get_coords_dataframe()"),
        ]

    ##########################################################
    # Helper functions used by Robot to extract and plot data
    ##########################################################

    @staticmethod
    def get_yvals_item_abifiles(item: Any, abifiles: list) -> np.ndarray:
        """Extract values for a list of Abinit files."""
        def _float(obj):
            if obj is None: return obj
            return float(obj)

        if callable(item):
            return np.array([_float(item(a)) for a in abifiles])
        else:
            return np.array([_float(duck.getattrd(a, item)) for a in abifiles])

    @staticmethod
    def plot_xvals_or_xstr_ax(ax, xs, yvals, fontsize, **kwargs) -> list:
        """Plot xs where xs can contain either numbers or strings."""
        if not is_string(xs[0]):
            lines = ax.plot(xs, yvals, **kwargs)
        else:
            # Must handle list of strings in a different way.
            xn = range(len(xs))
            lines = ax.plot(xn, yvals, **kwargs)
            ax.set_xticks(xn)
            ax.set_xticklabels(xs, fontsize=fontsize)

        return lines

    @staticmethod
    def plot_abs_conv(ax1, ax2, xs, yvals, abs_conv, xlabel, fontsize, hatch, **kwargs) -> None:
        """
        Plot |y - y_xmax| in log scale on ax2 and add hspan to ax1.
        """
        y_xmax = yvals[-1]
        span_style = dict(alpha=0.2, color="green", hatch=hatch)
        ax1.axhspan(y_xmax - abs_conv, y_xmax + abs_conv, label=r"$|y-y(x_{max})| \leq %s$" % abs_conv, **span_style)

        # Plot |y - y_xmax| in log scale on ax2.
        ax2.plot(xs, np.abs(yvals - y_xmax), **kwargs)
        ax2.set_yscale("log")
        ax2.set_ylabel(r"$|y-y_{x_{max}}|$", fontsize=fontsize)
        ax2.axhspan(0, abs_conv, label=r"$|y-y(x_{max})| \leq %s$" % abs_conv, **span_style)
        ax2.legend(loc="best", fontsize=fontsize, shadow=True)
        if xlabel:
            ax2.set_xlabel("%s" % xlabel)



class HueGroup:
    """
    This small object is used by ``group_and_sortby`` to store information about the group.
    """

    def __init__(self, hvalue, xvalues, abifiles, labels):
        """
        Args:
            hvalue: Hue value.
            xvalues: abifiles are sorted by ``func_or_string`` and these are the values
                associated to ``abifiles``.
            abifiles: List of file with this hue value.
            labels: List of labels associated to ``abifiles``.
        """
        self.hvalue = hvalue
        self.abifiles = abifiles
        self.labels = labels
        self.xvalues = xvalues
        assert len(abifiles) == len(labels)
        assert len(abifiles) == len(xvalues)

    def __len__(self) -> int:
        return len(self.abifiles)

    def __iter__(self):
        """Iterate over (label, abifile, xvalue)."""
        return zip(self.labels, self.abifiles, self.xvalues)



class RobotPythonScript:
    """
    Small object used to generate a python script that reconstructs the
    robot from a json file containing the list of files.
    Client code can then add additional logic to the script and write it to disk.

    This object is typically used in the `on_all_ok` method of Works
    to generate ready-to-use python scripts to post-process/visualize the results.

    Example:

        with gsr_robot.get_pyscript(work.outdir.path_in("gsr_robot.py")) as script:
            script.add_text("a = 1")
    """

    def __init__(self, robot: Robot, filepath_py: str):
        self.robot = robot
        self.filepath_py = filepath_py

        root, _ = os.path.splitext(filepath_py)
        self.filepath_json = root + ".json"
        basename = os.path.basename(filepath_py)

        self.pytext = f"""\
#!/usr/bin/env python
# This script has been automatically generated by AbiPy.
from __future__ import annotations

if False:
    import seaborn as sns
    sns.set(context="paper", style='darkgrid', palette='deep',
           font='sans-serif', font_scale=0.8, color_codes=False, rc=None)

from abipy.abio.robots import Robot
robot = Robot.from_json_file("{self.filepath_json}")
print(robot)

# Uncomment these two lines to produce an excel file
#df = robot.get_dataframe(with_geo=False)
#df.to_excel("{basename}.xls")
#df.to_csv("{basename}.csv")
"""

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Activated at the end of the with statement."""
        self.write()

    def add_text(self, text: str) -> None:
        """Add `text` to the script."""
        self.pytext += "\n" + text

    def write(self) -> None:
        """
        Write python script and json file with the list of files in the Robot.
        """
        with open(self.filepath_py, "wt") as fh:
            fh.write(self.pytext)
        make_executable(self.filepath_py)

        with open(self.filepath_json, "wt") as fh:
            fh.write(self.robot.to_json())
