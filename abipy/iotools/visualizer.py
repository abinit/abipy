"Define a class used to execute a visualizer within the Python interpreter."
from __future__ import division, print_function

import os.path

from subprocess import call
from collections import namedtuple, OrderedDict

from abipy.core import AbipyException
from abipy.tools import ask_yes_no, which

__all__ = [
    "Visualizer",
]

# cmdarg is the command line option that has to be provided to visualize a file with the given extension.
Application = namedtuple("Application", "name cmdarg")

#: One-to-many mapping file extension --> applications
ext2apps = OrderedDict( {
    "xsf": [Application("xcrysden", "--xsf"),
            Application("vsim", None),
            ],

    "bxsf": [Application("xcrysden", "--bxsf"),
             ],
} )

#: One-to-many mapping application_name --> file extensions supported.
appname2exts = dict()

for (ext, applications) in ext2apps.items():
    for app in applications:
        aname = app.name
        if aname not in appname2exts:
            appname2exts[aname] = [ext]
        else:
            appname2exts[aname].append(ext)

for aname in appname2exts: # Remove duplicated entries
    appname2exts[aname] = set(appname2exts[aname])

##########################################################################################


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
            application, executable = cls.appath_from_ext(ext)
        except Exception as exc:
            raise cls.Error(exc)

        return Visualizer(filename, executable, application.cmdarg)

    @staticmethod
    def appath_from_ext(ext):
        "Return the absolute path of the first (available) application that supports extension ext"
        if ext.startswith("."): ext = ext[1:]
        try:
            applications = ext2apps[ext]
        except KeyError:
            raise KeyError("Don't know how to handle extension: %s" % ext)

        for app in applications:
            executable = which(app.name)
            if executable is not None:
                return app, executable
        else:
            raise ValueError("No executable found for file: %s" % filename)

    @staticmethod
    def exts_from_appname(application_name):
        "Return the set of extensions supported by application_name"
        aname = os.path.basename(application_name)
        try:
            return appname2exts[aname]
        except KeyError:
            err_msg = "application %s is not supported" % application_name
            raise KeyError(err_msg)

    def __init__(self, filename, executable, cmdarg=None, wait=False):
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
        self.wait = wait

    def __call__(self):
        return self.show()

    def __str__(self):
        return "%s" % self.executable

    def set_wait(self, wait=True):
        self.wait = wait

    def show(self):
        """Call the visualizer in a subprocess to visualize the data."""
        if self.cmdarg is None:
            raise NotImplementedError()

        print("Executing: ", self.executable, self.cmdarg, self.filename)
        retcode = call([self.executable, self.cmdarg, self.filename])

        if self.wait:
            ask_yes_no("Enter [y/n] to continue or to exit.", default="yes")

        return retcode
