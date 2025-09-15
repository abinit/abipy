#!/usr/bin/env python

import sys
import os
#import abipy.data as abidata
#import abipy.abilab as abilab
import abipy.flowtk as flowtk

from abipy.dfpt.qha_general_stress import QHA_ZSISA
from abipy.flowtk.zsisa import ThermalRelaxWork

def build_flow(options):
    #nqsmall_or_qppa = 2
    #zsisa = QHA_ZSISA.from_json_file("flow_qha_zsisa/outdata/zsisa.json", nqsmall_or_qppa, verbose=1)
    #zsisa = QHA_ZSISA.pickle_load("/Users/giantomassi/git_repos/abipy/abipy/examples/flows/flow_qha_zsisa/outdata", basename="zsisa.json.pickle")
    #print(zsisa)

    temperatures = [10, 50]
    temperatures = [10, 50, 100, 200]
    pressures = [0]

    flow = flowtk.Flow(workdir="flow_thermal")
    work = ThermalRelaxWork.from_zsisa_flow("flow_qha_zsisa", temperatures, pressures)
    print(work)
    flow.register_work(work)
    return flow


@flowtk.flow_main
def main(options):
    """
    This is our main function that will be invoked by the script.
    flow_main is a decorator implementing the command line interface.
    Command line args are stored in `options`.
    """
    return build_flow(options)


if __name__ == "__main__":
    sys.exit(main())
