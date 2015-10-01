#!/usr/bin/env python
# coding: utf-8
"""Factory functions for Abinit input files """
from __future__ import print_function, division, unicode_literals

import numpy as np
import pymatgen.io.abinit.abiobjects as aobj

from collections import namedtuple
from pymatgen.io.abinit.pseudos import PseudoTable
from abipy.core.structure import Structure
#from .input import AbiInput
from abipy.htc.factories import *

import abipy.abilab as abilab
import abipy.data as abidata

import logging
logger = logging.getLogger(__file__)

from pymatgen.io.abinit.launcher import BatchLauncher


@abilab.flow_main
def main(options):
    structure = abidata.structure_from_ucell("MgB2")
    pseudos = abilab.PseudoTable(abidata.pseudos("12mg.pspnc", "5b.pspnc"))

    #structure = abidata.structure_from_ucell("NiO")
    #pseudos = abidata.pseudos("28ni.paw", "8o.2.paw")

    structure.scale_lattice(structure.volume * 1.2)

    inp = ion_ioncell_relax_input(structure, pseudos, ecut=4, pawecutdg=8)

    inp.chkprim = 0
    inp.paral_kgb = 1
    ion_inp, ioncell_inp = inp.split_datasets()

    # Create the flow
    flow = abilab.Flow("flow_debug_relax")

    relax_work = abilab.RelaxWork(ion_inp, ioncell_inp)
    flow.register_work(relax_work)
    return flow


if __name__ == "__main__":
    import sys
    sys.exit(main())
