# coding: utf-8
"""Define a class used to execute a visualizer within the Python interpreter."""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import abc
import six

from monty.os.path import which
from monty.termcolor import cprint
from monty.functools import lazy_property


__all__ = [
    "Visualizer",
]


def is_macosx():
    """True if we are running on Mac."""
    return "darwin" in sys.platform


def find_loc(app_name):
    """
    Returns the location of the application from its name. None if not found.
    """
    path = _find_loc(app_name)
    if path is not None: return path
    path = _find_loc(app_name.upper())
    return path


def _find_loc(app_name):  # pragma: no cover
    # Try command line version
    path = which(app_name)
    if path is not None: return path

    # Treat Mac OsX applications.
    if is_macosx():
        system_apps = [f.lower() for f in os.listdir("/Applications")]
        try:
            i = system_apps.index(app_name)
            return os.path.join("/Applications", system_apps[i])
        except ValueError:
            try:
                i = system_apps.index(app_name + ".app")
                return os.path.join("/Applications", system_apps[i] + ".app")
            except ValueError:
                pass

        try:
            user_dir = os.path.expanduser("~/Applications/")
            user_apps = [f.lower() for f in os.listdir(user_dir)]
            try:
                i = user_apps.index(app_name)
                return os.path.join(user_dir, user_apps[i])
            except ValueError:
                try:
                    i = user_apps.index(app_name + ".app")
                    return os.path.join("/Applications", user_apps[i] + ".app")
                except ValueError:
                    pass
        except Exception:
            pass

    return None


class VisualizerError(Exception):
    """Base class for Visualizer errors"""


@six.add_metaclass(abc.ABCMeta)
class Visualizer(object):
    """
    Handle the visualization of data.
    """
    # True if its a Mac OsX applications (default is unix executable).
    # If we have a Mac OsX application we have to run it with "open -a app_name --args"
    is_macosx_app = False

    Error = VisualizerError

    def __init__(self, filepath):
        """
        Args:
            filepath: Name of the file to visualize
        """
        self.filepath = os.path.abspath(filepath)

    def __str__(self):
        return "%s: %s, is_macosx_app %s, filepath: %s" % (
            self.__class__.__name__, self.binpath, self.is_macosx_app, self.filepath)

    def __call__(self):  # pragma: no cover
        """
        Call the visualizer in a subprocess.

        Returns: exit status of the subprocess.
        """
        from subprocess import call
        if not self.is_macosx_app:
            cprint("Executing: %s %s %s" % (self.binpath, self.cmdarg, self.filepath), "yellow")
            return call([self.binpath, self.cmdarg, self.filepath])

        else:
            # Mac-OSx applications can be launched with `open -a Vesta --args si.cif`
            cmd = "open -a %s --args %s %s" % (self.name, self.cmdarg, self.filepath)
            cprint("Executing MacOSx open command: %s" % cmd, "yellow")
            return call(cmd, shell=True)

    @property
    def cmdarg(self):
        """Arguments that must be used to visualize the file."""
        root, ext = os.path.splitext(self.filepath)
        ext = ext.replace(".", "")
        #print(root, ext)
        for e, args in self.EXTS:
            if e == ext:
                return args
        return " "

    @lazy_property
    def binpath(self):
        """Absolute path of the binary. None if app is not found"""
        return find_loc(self.name)

    @property
    def is_available(self):
        """True is the visualizer is available on the local machine."""
        return self.binpath is not None

    @classmethod
    def get_available(cls, ext=None):
        """
        List of visualizers available on the local host.
        If ext is not None, only the visualizers supporting this extension are returned.
        """
        visus = [v for v in cls.__subclasses__() if v.is_available]
        if ext is None: return visus
        return [v for v in visus if v.support_ext(ext)]

    @classmethod
    def support_ext(cls, ext):
        """True if visualizer supports the extension ext."""
        if ext.startswith("."): ext = ext[1:]
        return any(e == ext for e, _ in cls.EXTS)

    @classmethod
    def from_file(cls, filepath):
        """
        Initialize a subclass of :class:`Visualizer` from filepath, the application
        is chosen automatically depending on the file extension.

        Raise:
            `VisualizerError` if no visualizer is available for the given file.
        """
        # Get the file extension.
        root, ext = os.path.splitext(filepath)
        if not ext:
            raise ValueError("Cannot detect file extension in %s " % filepath)

        avail_visus = cls.get_available(ext=ext)

        if not avail_visus:
            raise cls.Error("Cannot find available visualizer for extension %s " % ext)

        return avail_visus[0](filepath)

    @classmethod
    def supported_extensions(cls):
        """List of file extensions supported by the visualizer."""
        return [e for (e, args) in cls.EXTS]

    @classmethod
    def from_name(cls, appname):
        """Return the visualizer class from the name of the application."""
        for visu in cls.__subclasses__():
            if visu.name == appname:
                return visu

        raise cls.Error("appname is not among the list of supported visualizers %s " % appname)

    @classmethod
    def all_visunames(cls):
        """List with the names of the visualizers supported."""
        return sorted([visu.name for visu in cls.__subclasses__()])


####################
# Concrete classes #
####################

class Xcrysden(Visualizer):
    name = "xcrysden"

    EXTS = [
        ("xsf", "--xsf"),
        ("bxsf", "--bxsf"),
    ]


class V_Sim(Visualizer):
    name = "v_sim"

    EXTS = [
        ("xsf", "--xsf"),
    ]


class Vesta(Visualizer):
    is_macosx_app = is_macosx()

    name = "vesta"

    EXTS = [
        ("xsf", ""),
    ]


class Ovito(Visualizer):
    is_macosx_app = is_macosx()

    name = "ovito"

    EXTS = [
        ("POSCAR", ""),
    ]


class Avogadro(Visualizer):
    is_macosx_app = is_macosx()

    name = "avogadro"

    EXTS = [
       ("cif", ""),
    ]
