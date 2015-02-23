# coding: utf-8
"""Base classes and utils for lessons."""
from __future__ import print_function, division, unicode_literals

import sys
import os
import six
import abc
import shutil
import abipy.data as abidata

from monty.functools import lazy_property
from abipy.tools.notebooks import mpld3_enable_notebook


class BaseLesson(six.with_metaclass(abc.ABCMeta, object)):

    def __init__(self, **kwargs):
        mode = kwargs.get("mode") #, "ipython-shell")
        if mode == "mpld3": 
            mpld3_enable_notebook()

    @abc.abstractproperty
    def abipy_string(self):
        """the abipy lesson."""

    @abc.abstractproperty
    def comline_string(self):
        """the commandline lesson."""

    @abc.abstractproperty
    def pyfile(self):
        """Absolute Path of the python script."""

    def get_local_copy(self):
        """Copy this script to the current working dir to explore it and edit"""
        dst = os.path.basename(self.pyfile)
        if os.path.exists(dst):
            raise RuntimeError("file %s already exists. Remove it before calling get_local_copy" % dst)
        shutil.copyfile(self.pyfile, dst)
        self.get_local_copy()

    def __repr__(self):
        """String representation."""
        return self.abipy_string

    @property
    def manpath(self):
        return self.pyfile.replace(".py", ".man")

    #def manfile(self, what=None):
    #    """The path to the man file of the lesson. Use `%man %lesson.manfile` to open it in ipython"""
    #    _, ext = os.path.splitext(self.pyfile)
    #    man_path = self.pyfile.replace(ext, '.man')
    #    try:
    #        shutil.copy(os.path.join(os.path.dirname(os.path.realpath(__file__)), man_path), '.')
    #    except IOError:
    #        with open(man_path, "wt") as fh:
    #            fh.write(self._pandoc_convert(to="man", what=what, extra_args=("-s",)))

    def setup(self):
        lesson_name = os.path.basename(self.pyfile)[:-3]
        print("\n"
              "The lesson %s " % lesson_name + "is prepared.\n"
              "Use\n\n"
              "     man ./" + lesson_name + ".man\n\n"
              "to view the text of the lesson.\n")

        # Copy man file (except when we are running inside lessons in develop mode
        # Use. python __setup__.py to regenerate the man files.
        dst = os.path.join(os.path.abspath(os.getcwd()), os.path.basename(self.manpath))
        if not dst == self.manpath:
            shutil.copyfile(self.manpath, dst)

    def _pandoc_convert(self, to, what, extra_args=()):
        if what is None:
            what = self.abipy_string
        try:
            import pypandoc
            return pypandoc.convert(what, to, "rst", extra_args=extra_args)
        except (OSError, ImportError):
            return "pypandoc.convert failed. Please install pandoc and pypandoc"

    def _gen_manfile(self):
        _, ext = os.path.splitext(self.pyfile)
        man_path = self.pyfile.replace(ext, '.man')
        with open(man_path, "wt") as fh:
            fh.write(self._pandoc_convert(to="man", what=self.comline_string, extra_args=("-s",)))

    def _repr_html_(self):
        """Support for ipython notebooks."""
        from docutils.core import publish_string, publish_parts
        return publish_string(self.abipy_string, writer_name="html")

    @staticmethod
    def docvar(varname):
        from abipy.htc.abivars_db import get_abinit_variables
        if varname == "inputvariable":
            return ("inputvariable is a very complicated input variable. Better to ask for help immediately")
#            return ("When we said that you should pass `inputvariable` to docvar we meant that"
#                    "you should pass a string with the name of a valid ABINT variable e.g. `ecut`"
#                    "not `inputvariable` :)")
        return get_abinit_variables()[varname]

    @property
    def abidata(self):
        """Abipy data files."""
        return abidata


def get_pseudos(structure, extension='oncvpsp'):
    """
    returns a list of pseudos names for structure. This list should be fed to abidata.pseudos like
    abidata.pseudos(get_pseudos(structure))
    """
    pseudos = []
    for element in structure.composition.elements:
        pseudos.append(abidata.pseudo(str(element)+'.'+extension))
        #todo test if the pseudo potential file exists
    return pseudos
