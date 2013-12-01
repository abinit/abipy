"Define a class used to execute a visualizer within the Python interpreter."
from __future__ import division, print_function

import sys
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
_EXT2APPS = OrderedDict({
    "xsf": [App("xcrysden", "--xsf"),
            App("v_sim", ""),
            App("VESTA", ""),
            ],

    "bxsf": [App("xcrysden", "--bxsf"),
             ],
})

# One-to-many mapping: app_name --> file extensions supported.
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


def is_macosx():
    """True if we are running on Mac."""
    return "darwin" in sys.platform


def find_app(app_name):
    """Returns the location of the application from its name."""
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
            pass

        user_dir = os.path.expanduser("~/Applications/")
        user_apps = [f.lower() for f in os.listdir(user_dir)]
        try:
            i = user_apps.index(app_name)
            return os.path.join(user_dir, user_apps[i])
        except ValueError:
            pass

    return None


class VisualizerError(AbipyException):
    """Base class for Visualizer errors"""


class Visualizer(object):
    """
    Handle the visualization of data.
    """
    Error = VisualizerError

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
            app, executable = cls.app_path_from_ext(ext)

        except Exception as exc:
            raise cls.Error(str(exc))

        return Visualizer(filename, executable, app.cmdarg)

    @staticmethod
    def app_path_from_ext(ext):
        """Return the absolute path of the first (available) application that supports extension ext"""
        if ext.startswith("."): 
            ext = ext[1:]

        try:
            apps = _EXT2APPS[ext]
        except KeyError:
            raise self.Error("Don't know how to handle extension: %s" % ext)

        for app in apps:
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

    def show(self):
        """
        Call the visualizer in a subprocess to visualize the data.

        Returns: exit status of the subprocess.
        """
        #print("Executing: ", self.executable, self.cmdarg, self.filename)
        from subprocess import call
        return call([self.executable, self.cmdarg, self.filename])

        # NOTE: Mac-OSx applications can be launched with
        #open -a Vesta --args /Users/gmatteo/Coding/abipy/abipy/data/cifs/si.cif


class _Visualizer(object):
    """
    Handle the visualization of data.
    """
    # True if its a MacOsx applications (default is unix executable).
    is_macosx_app = False
    
    Error = VisualizerError

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
        return "%s: %s, is_macosx_app %s" % (self.__class__.__name__, self.executable, self.is_macosx_app)

    def __call__(self):
        """
        Call the visualizer in a subprocess to visualize the data.

        Returns: exit status of the subprocess.
        """
        #print("Executing: ", self.executable, self.cmdarg, self.filename)
        from subprocess import call

        if not self.is_macosx_app:
            return call([self.executable, self.cmdarg, self.filename])

        else:
            # NOTE: Mac-OSx applications can be launched with
            #open -a Vesta --args si.cif
            cmd = "open -a %s --args %s %s" % (self.executable, self.cmdarg, self.filename)
            return call(cmd)

    @classmethod
    def get_available(cls, ext=None):
        """
        List of visualizers available on the local host.
        ext is used to filter the visualizers that support the given extension.
        """
        visus = [v for v in cls.__subclasses__() if v.app_name is not None]
        if ext is None: return visus

        return [v for v in visus if v.support_ext(ext)]

    def support_ext(self, ext):
        """True if visualizer supports the extension ext."""
        return any(e == ext for e, _ in self.EXTS)
        #for e, _ in self.EXTS:
        #    if e == ext: return True
        #return False

    @classmethod
    def from_file(cls, filename):
        """
        Initialize the instance from filename, application 
        is chosen automatically depending on the file extension.
        """
        # Get the file extension.
        root, ext = os.path.splitext(filename)

        if not ext:
            raise cls.Error("Cannot detect file extension in %s " % filename)

        try:
            app, executable = cls.app_path_from_ext(ext)

        except Exception as exc:
            raise cls.Error(str(exc))

        return Visualizer(filename, executable, app.cmdarg)

    @staticmethod
    def app_path_from_ext(ext):
        """
        Return the absolute path of the first (available) 
        application that supports extension ext
        """
        if ext.startswith("."): 
            ext = ext[1:]

        try:
            apps = _EXT2APPS[ext]
        except KeyError:
            raise self.Error("Don't know how to handle extension: %s" % ext)

        for app in apps:
            executable = which(app.name)
            if executable is not None:
                return app, executable
        else:
            raise self.Error("No executable found for file: %s" % filename)


class Xcrysden(_Visualizer):
    app_name = find_app("xcrysden")

    EXTS = [
        ("xsf", "--xsf"),
        ("bxsf", "--bxsf")
    ]

class V_Sim(_Visualizer):
    app_name = find_app("v_sim")

    EXTS = [
        ("xsf", "--xsf"),
    ]

class Vesta(_Visualizer):
    is_macosx_app = is_macosx()

    app_name = find_app("vesta")
        
    EXTS = [
        ("xsf", ""),
    ]

class Avogadro(_Visualizer):
    is_macosx_app = is_macosx()
    app_name = find_app("avogadro")
        
    EXTS = [
        ("xsf", ""),
    ]

if __name__ == "__main__":
    #print("ismac", is_macosx())
    print("available:", _Visualizer.get_available())

    #Vesta
