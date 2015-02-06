#!/usr/bin/env python
"""
reusable functions for lessons
"""
from __future__ import division, print_function, unicode_literals, division

import sys
import os
import shutil
import yaml

import abipy.abilab as abilab
import abipy.data as abidata

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

