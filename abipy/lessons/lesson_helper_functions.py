#!/usr/bin/env python
"""
reusable functions for lessons
"""

from __future__ import division, print_function

import sys
import os
import shutil
import yaml
import html2text
from variables import *
import abipy.abilab as abilab
import abipy.data as abidata


class VariableDatabase(object):

    all_vars = None

    def __init__(self):
        self.load_vars('abinit_vars.yml')

    def load_vars(self, file_yml):

        f_var = open(file_yml, 'r')
        variables = yaml.load(f_var)
        f_var.close()

        self.all_vars = dict()

        for var in variables:
            self.all_vars[var.varname] = var

    def get_var(self, variable):
        return self.all_vars[variable]

def help(stream=sys.stdout):
    """
    Display the tutorial text.
    """
    stream.write(__doc__)


def get_local_copy():
    """
    Copy this script to the current working dir to explore and edit
    """
    dst = os.path.basename(__file__[:-1])
    if os.path.exists(dst):
        raise RuntimeError("file %s already exists. Remove it before calling get_local_copy" % dst)
    shutil.copyfile(__file__[:-1], dst)


def abinit_help(inputvariable):
    """
    Print the abinit documentation on the abinit input variable 'inputvariable'
    """
    database = VariableDatabase()
    var = database.get_var(inputvariable)
    text = html2text.html2text("<h2>Default value : </h2>"+str(var.defaultval)+"<br /><h2>Description</h2>"+str(var.text))
    print(text.replace("[[", "\033[1m").replace("]]", "\033[0m"))


