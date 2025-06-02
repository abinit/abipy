#!/usr/bin/env python

from abipy.dfpt.qha_general_stress import QHA_ZSISA

nqsmall_or_qppa = 4
qha = QHA_ZSISA.from_json_file("flow_qha_zsisa/outdata/zsisa.json", nqsmall_or_qppa, verbose=1)
print(qha)
