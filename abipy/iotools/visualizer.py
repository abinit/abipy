"Define a class used to execute a visualizer within the Python interpreter."
from __future__ import division, print_function

import sys
import os

from collections import namedtuple, OrderedDict
from abipy.core import AbipyException
from abipy.tools import ask_yes_no, which

__all__ = [
    "Visualizer",
]

# cmdarg is the command line option that has to be provided to visualize a file with the given extension.
#App = namedtuple("App", "name cmdarg")
#
## One-to-many mapping file extension --> applications
#_EXT2APPS = OrderedDict({
#    "xsf": [App("xcrysden", "--xsf"),
#            App("v_sim", ""),
#            App("VESTA", ""),
#            ],
#
#    "bxsf": [App("xcrysden", "--bxsf"),
#             ],
#})
#
## One-to-many mapping: app_name --> file extensions supported.
#appname2exts = {}
#
#for (ext, applications) in _EXT2APPS.items():
#    for app in applications:
#        aname = app.name
#        if aname not in appname2exts:
#            appname2exts[aname] = [ext]
#        else:
#            appname2exts[aname].append(ext)
#
#
#for aname in appname2exts: # Remove duplicated entries
#    appname2exts[aname] = set(appname2exts[aname])
#
#
#def supported_visunames():
#    """List of strings with the name of the supported visualizers."""
#    return list(appname2exts.keys())


def is_macosx():
    """True if we are running on Mac."""
    return "darwin" in sys.platform


def find_loc(app_name):
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
            try:
                i = system_apps.index(app_name + ".app")
                return os.path.join("/Applications", system_apps[i] + ".app")
            except ValueError:
                pass

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

    return None


class VisualizerError(AbipyException):
    """Base class for Visualizer errors"""


#class Visualizer(object):
#    """
#    Handle the visualization of data.
#    """
#    Error = VisualizerError
#
#    def __init__(self, filename, executable, cmdarg=""):
#        """
#        Args:
#            filename: 
#                Name of the file to visualize
#            executable: 
#                Name of the visualizer (absolute path or simple name).
#            cmdarg: 
#                command line option passed to the visualizer.
#        """
#        self.executable = which(executable)
#        self.filename = os.path.abspath(filename)
#        self.cmdarg = cmdarg
#
#    def __str__(self):
#        return "%s" % self.executable
#
#    def __call__(self):
#        return self.show()
#
#    @classmethod
#    def from_file(cls, filename):
#        """
#        Initialize the instance from filename, application is chosen automatically 
#        depending on the file extension.
#        """
#        root, ext = os.path.splitext(filename)
#
#        if not ext:
#            raise cls.Error("Cannot detect file extension in %s " % filename)
#
#        try:
#            app, executable = cls.app_path_from_ext(ext)
#
#        except Exception as exc:
#            raise cls.Error(str(exc))
#
#        return Visualizer(filename, executable, app.cmdarg)
#
#    @staticmethod
#    def app_path_from_ext(ext):
#        """Return the absolute path of the first (available) application that supports extension ext"""
#        if ext.startswith("."): 
#            ext = ext[1:]
#
#        try:
#            apps = _EXT2APPS[ext]
#        except KeyError:
#            raise self.Error("Don't know how to handle extension: %s" % ext)
#
#        for app in apps:
#            executable = which(app.name)
#            if executable is not None:
#                return app, executable
#        else:
#            raise self.Error("No executable found for file: %s" % filename)
#
#    @staticmethod
#    def exts_from_appname(app_name):
#        """Return the set of extensions supported by app_name"""
#        name = os.path.basename(app_name)
#        try:
#            return appname2exts[name]
#        except KeyError:
#            raise self.Error("application %s is not supported" % app_name)
#
#    def show(self):
#        """
#        Call the visualizer in a subprocess to visualize the data.
#
#        Returns: exit status of the subprocess.
#        """
#        #print("Executing: ", self.executable, self.cmdarg, self.filename)
#        from subprocess import call
#        return call([self.executable, self.cmdarg, self.filename])
#
#        # NOTE: Mac-OSx applications can be launched with
#        #open -a Vesta --args /Users/gmatteo/Coding/abipy/abipy/data/cifs/si.cif



class MetaClass(type):
    def __str__(self):
        return "%s: bin: %s, macosx_app: %s" % (self.__class__.__name__, self.bin, self.is_macosx_app)


class Visualizer(object):
    """
    Handle the visualization of data.
    """
    __metaclass__ = MetaClass

    # True if its a MacOsx applications (default is unix executable).
    is_macosx_app = False
    
    Error = VisualizerError

    def __init__(self, filepath):
        """
        Args:
            filepath: 
                Name of the file to visualize
        """
        self.filepath = os.path.abspath(filepath)

    def __str__(self):
        return "%s: %s, is_macosx_app %s, filepath: %s" % (
            self.__class__.__name__, self.bin, self.is_macosx_app, self.filepath)

    def __call__(self):
        """
        Call the visualizer in a subprocess to visualize the data.

        Returns: exit status of the subprocess.
        """
        from subprocess import call

        if not self.is_macosx_app:
            print("Executing: ", self.bin, self.cmdarg, self.filepath)
            return call([self.bin, self.cmdarg, self.filepath])

        else:
            # NOTE: Mac-OSx applications can be launched with
            #open -a Vesta --args si.cif
            cmd = "open -a %s --args %s %s" % (self.name, self.cmdarg, self.filepath)
            print("Executing Mac open: %s" % cmd)
            return call(cmd, shell=True)

    @property
    def cmdarg(self):
        root, ext = os.path.splitext(self.filepath)
        ext = ext.replace(".", "")
        #print(root, ext)
        for e, args in self.EXTS:
            if e == ext:
                return args
        return " "

    @property
    def is_available(self):
        """True is the visualizer is available on the local machine."""
        return self.bin is not None

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
        Initialize a subclass of `Visualizer` from filepath, the application 
        is chosen automatically depending on the file extension.

        Raiseds:
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
        return [e for (e, args) in cls.EXTS]

    @classmethod
    def from_name(cls, visu_name):
        for visu in cls.__subclasses__(): 
            if visu.name == visu_name:
                return visu

        raise cls.Error("visu_name is not among the list of supported visualizers %s " % visu_name)

####################
# Concrete classes #
####################

class Xcrysden(Visualizer):
    name = "xcrysden"
    bin = find_loc(name)

    EXTS = [
        ("xsf", "--xsf"),
        ("bxsf", "--bxsf")
    ]

class V_Sim(Visualizer):
    name = "v_sim"
    bin = find_loc(name)

    EXTS = [
        ("xsf", "--xsf"),
    ]

class Vesta(Visualizer):
    is_macosx_app = is_macosx()

    name = "vesta"
    bin = find_loc(name)
        
    EXTS = [
        ("xsf", ""),
    ]

#class Avogadro(Visualizer):
#    is_macosx_app = is_macosx()
#
#    name = "avogadro"
#    bin = find_loc(name)
#        
#    EXTS = [
#    ]

if __name__ == "__main__":
    #print("ismac", is_macosx())
    print("available:") 
    for visu in _Visualizer.get_available():
        print(visu)
    import sys
    filename = sys.argv[1]
    v = Xcrysden(filename)
    v = Vesta(filename)
    v()

