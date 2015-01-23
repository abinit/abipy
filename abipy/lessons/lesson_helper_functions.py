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
from pymatgen.matproj.rest import MPRester, MPRestError
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


class MPConnection(object):
    """
    connection to the materials project database, should maby be sitting elsewhere ..
    """
    def __init__(self, mp_key=None):
        if mp_key is None:
            try:
                mp_key = os.environ['MP_KEY']
            except OSError:
                mp_key = raw_input("there is no key for accesing the materials projects database\n"
                                   "please provide one. (if you don't have one visit the materials \n"
                                   "project website to generate one) :")
        if len(mp_key) > 4:
            self.mp_key = mp_key

    def structure_from_mp(self, mpid):
        """
        return a structure from the materials datebase
        """
        with MPRester(self.mp_key) as mp_database:
            return mp_database.get_structure_by_material_id(mpid, final=True)


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
