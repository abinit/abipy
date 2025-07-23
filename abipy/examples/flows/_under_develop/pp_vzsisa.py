#!/usr/bin/env python
import numpy as np

from abipy.dfpt.vzsisa import Vzsisa
nqsmall = 4
qha = Vzsisa.from_json_file("flow_qha_vzsisa/outdata/vzsisa.json", nqsmall, verbose=0)

#print(qha)
#print(qha.bo_energies)
#tot_en = qha.bo_energies[np.newaxis, :].T
#print(tot_en)

