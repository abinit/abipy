# coding: utf-8
from __future__ import unicode_literals, division, print_function

import sys
import abipy.data as abidata
import abipy.abilab as abilab
import abipy.fworks.fw_tasks as fw_tasks

from abipy.abio.factories import *
from abipy.core.testing import AbipyTest



class AbiFireTaskTest(AbipyTest):

    def setUp(self):
        si = abilab.Structure.from_file(abidata.cif_file("si.cif"))
        self.si_scf_input = ebands_input(si, abidata.pseudos("14si.pspnc"), ecut=2, kppa=10).split_datasets()[0]

    def testAbiFireTask(self):
        task = fw_tasks.AbiFireTask(self.si_scf_input)
        print(task.to_dict())
        self.assertFwSerializable(task)
