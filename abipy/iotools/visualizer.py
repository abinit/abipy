# coding: utf-8
"""Define a class used to execute a visualizer within the Python interpreter."""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import abc
import six

from monty.os.path import which

__all__ = [
    "Visualizer",
]


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
        except:
            pass

    return None


class VisualizerError(Exception):
    """Base class for Visualizer errors"""


class MetaClass(type):
    """Provides str representation of the class."""
    def __str__(self):
        return "%s: bin: %s, macosx_app: %s" % (self.__class__.__name__, self.bin, self.is_macosx_app)


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
            self.__class__.__name__, self.bin, self.is_macosx_app, self.filepath)

    def __call__(self):
        """
        Call the visualizer in a subprocess.

        Returns: exit status of the subprocess.
        """
        from subprocess import call

        if not self.is_macosx_app:
            #print("Executing: ", self.bin, self.cmdarg, self.filepath)
            return call([self.bin, self.cmdarg, self.filepath])

        else:
            # NOTE: Mac-OSx applications can be launched with
            #open -a Vesta --args si.cif
            cmd = "open -a %s --args %s %s" % (self.name, self.cmdarg, self.filepath)
            #print("Executing Mac open: %s" % cmd)
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
    def from_name(cls, visu_name):
        """Return the visualizer class from the name of the application."""
        for visu in cls.__subclasses__(): 
            if visu.name == visu_name:
                return visu

        raise cls.Error("visu_name is not among the list of supported visualizers %s " % visu_name)

    @classmethod
    def all_visunames(cls):
        """List with the names of the visualizers supported."""
        return [visu.name for visu in cls.__subclasses__()]


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
    print("available visualizers:") 
    for visu in Visualizer.get_available(): print(visu)
    import sys
    filename = sys.argv[1]
    v = Xcrysden(filename)
    v = Vesta(filename)
    v()
