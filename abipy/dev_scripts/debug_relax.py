#!/usr/bin/env python
# coding: utf-8
"""Factory functions for Abinit input files """
from __future__ import print_function, division, unicode_literals

import numpy as np
import pymatgen.io.abinitio.abiobjects as aobj

from collections import namedtuple
from pymatgen.io.abinitio.pseudos import PseudoTable
from abipy.core.structure import Structure
#from .input import AbiInput
from abipy.htc.factories import *

import abipy.abilab as abilab
import abipy.data as abidata

import logging
logger = logging.getLogger(__file__)


@abilab.flow_main
def main(options):
    #structure = abidata.structure_from_ucell("MgB2")
    #pseudos = abilab.PseudoTable(abidata.pseudos("12mg.pspnc", "5b.pspnc"))

    structure = abidata.structure_from_ucell("NiO")
    pseudos = abidata.pseudos("28ni.paw", "8o.2.paw")


    inp = ion_ioncell_relax_input(structure, pseudos, ecut=8, pawecutdg=12)

    inp.chkprim = 0
    ion_inp, ioncell_inp = inp.split_datasets()

    # Create the flow
    flow = abilab.Flow("flow_debug_relax")

    relax_work = abilab.RelaxWork(ion_inp, ioncell_inp)
    flow.register_work(relax_work)
    return flow.allocate()


if __name__ == "__main__":
    import sys
    sys.exit(main())
