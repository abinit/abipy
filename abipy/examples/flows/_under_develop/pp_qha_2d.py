#!/usr/bin/env python

from abipy.dfpt.qha_2D import QHA_2D

nqsmall = 4
qha = QHA_2D.from_json_file("flow_qha_2d/outdata/qha_2d.json", nqsmall, verbose=0)
print(qha)
