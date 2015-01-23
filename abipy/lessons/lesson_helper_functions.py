#!/usr/bin/env python
"""
reusable functions for lessons
"""
from __future__ import division, print_function, unicode_literals, division

import sys
import os
import shutil
import yaml
import html2text
import abipy.abilab as abilab
import abipy.data as abidata

from variables import *

__VARS_DATABASE = None


def get_abinit_variables():
    """Returns the database with the description of the ABINIT variables."""
    global __VARS_DATABASE
    if __VARS_DATABASE is None: __VARS_DATABASE = VariableDatabase()
    return __VARS_DATABASE
        


class VariableDatabase(object):

    all_vars = None

    def __init__(self):
        self.load_vars(abidata.var_file('abinit_vars.yml'))

    def load_vars(self, file_yml):

        f_var = open(file_yml, 'r')
        variables = yaml.load(f_var)
        f_var.close()

        self.all_vars = dict()

        for var in variables:
            self.all_vars[var.varname] = var

    def get_var(self, variable):
        return self.all_vars[variable]


def get_pseudos(structure, extension='oncvpsp'):
    """
    returns a list of pseudos names for structure. This list should be fed to abidata.pseudos like
    abidata.pseudos(get_pseudos(structure))
    """
    pseudos = []
    for element in structure.composition.elements:
        pseudos.append(element+'.'+extension)
        #todo test if the pseudo potential file exists
    return pseudos


def abinit_help(inputvariable):
    """
    Print the abinit documentation on the abinit input variable 'inputvariable'
    """
    database = get_abinit_variables()
    var = database.get_var(inputvariable)
    text = html2text.html2text("<h2>Default value : </h2>"+str(var.defaultval)+"<br /><h2>Description</h2>"+str(var.text))
    print(text.replace("[[", "\033[1m").replace("]]", "\033[0m"))
