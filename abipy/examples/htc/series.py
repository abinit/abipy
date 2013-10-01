#!/usr/bin/env python
from __future__ import division

import abipy.data as data
from abipy.htc.input import AbiInput

tsmear_list = [0.005, 0.01]
ngkpt_list = [[4,4,4], [8,8,8]]
occopt_list = [3, 4]

inp = AbiInput(pseudos=data.pseudos("14si.pspnc"), ndtset=len(tsmear_list))

inp.linspace("tsmear", start=tsmear_list[0], stop=tsmear_list[-1])
print(inp)

inp = AbiInput(pseudos=data.pseudos("14si.pspnc"), ndtset=len(tsmear_list) * len(ngkpt_list))

inp.product("tsmear", "ngkpt", tsmear_list, ngkpt_list)
print(inp)

# If you don't want to use multiple datasets in your calculation,
# you can split the initial input into ndtset different inputs.
separated_inps = inp.split_datasets()

for inp in separated_inps:
    print(inp)

# product accepts an arbitrary number of variables.
inp = AbiInput(pseudos=data.pseudos("14si.pspnc"), ndtset=len(tsmear_list) * len(ngkpt_list) * len(occopt_list))

inp.product("tsmear", "ngkpt", "occopt", tsmear_list, ngkpt_list, occopt_list)
print(inp)
