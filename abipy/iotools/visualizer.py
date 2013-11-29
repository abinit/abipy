"Define a class used to execute a visualizer within the Python interpreter."
from __future__ import division, print_function

import os

from collections import namedtuple, OrderedDict
from abipy.core import AbipyException
from abipy.tools import ask_yes_no, which

__all__ = [
    "Visualizer",
    "supported_visunames",
]

# cmdarg is the command line option that has to be provided to visualize a file with the given extension.
App = namedtuple("App", "name cmdarg")

# One-to-many mapping file extension --> applications
_EXT2APPS = OrderedDict( {
    "xsf": [App("xcrysden", "--xsf"),
            App("v_sim", ""),
            App("VESTA", ""),
            ],

    "bxsf": [App("xcrysden", "--bxsf"),
             ],
} )

# NOTE: Mac-OSx applications can be launched with
#open -a Vesta --args /Users/gmatteo/Coding/abipy/abipy/data/cifs/si.cif

#: One-to-many mapping app_name --> file extensions supported.
appname2exts = {}

for (ext, applications) in _EXT2APPS.items():
    for app in applications:
        aname = app.name
        if aname not in appname2exts:
            appname2exts[aname] = [ext]
        else:
            appname2exts[aname].append(ext)

for aname in appname2exts: # Remove duplicated entries
    appname2exts[aname] = set(appname2exts[aname])


def supported_visunames():
    """List of strings with the name of the supported visualizers."""
    return list(appname2exts.keys())


class VisualizerError(AbipyException):
    """Base class for Visualizer errors"""


class Visualizer(object):
    """
    Handle the visualization of data.
    """
    Error = VisualizerError

    @classmethod
    def from_file(cls, filename):
        """
        Initialize the instance from filename, application is chosen automatically 
        depending on the file extension.
        """
        root, ext = os.path.splitext(filename)

        if not ext:
            raise cls.Error("Cannot detect file extension in %s " % filename)

        try:
            app, executable = cls.appath_from_ext(ext)

        except Exception as exc:
            raise cls.Error(str(exc))

        return Visualizer(filename, executable, app.cmdarg)

    @staticmethod
    def appath_from_ext(ext):
        """Return the absolute path of the first (available) application that supports extension ext"""
        if ext.startswith("."): ext = ext[1:]
        try:
            applications = _EXT2APPS[ext]
        except KeyError:
            raise self.Error("Don't know how to handle extension: %s" % ext)

        for app in applications:
            executable = which(app.name)
            if executable is not None:
                return app, executable
        else:
            raise self.Error("No executable found for file: %s" % filename)

    @staticmethod
    def exts_from_appname(app_name):
        """Return the set of extensions supported by app_name"""
        name = os.path.basename(app_name)
        try:
            return appname2exts[name]
        except KeyError:
            raise self.Error("application %s is not supported" % app_name)

    def __init__(self, filename, executable, cmdarg=""):
        """
        Args:
            filename: 
                Name of the file to visualize
            executable: 
                Name of the visualizer (absolute path or simple name).
            cmdarg: 
                command line option passed to the visualizer.
        """
        self.executable = which(executable)
        self.filename = os.path.abspath(filename)
        self.cmdarg = cmdarg

    def __str__(self):
        return "%s" % self.executable

    def __call__(self):
        return self.show()

    def show(self):
        """
        Call the visualizer in a subprocess to visualize the data.

        Returns: exit status of the subprocess.
        """
        #print("Executing: ", self.executable, self.cmdarg, self.filename)
        from subprocess import call
        return call([self.executable, self.cmdarg, self.filename])
