#!/usr/bin/env python
from __future__ import print_function
import sys

from abipy.data import AbinitFilesGenerator

class MyGenerator(AbinitFilesGenerator):
    """This class generates the output files used in the unit tests and in the examples."""
    # Subclasses must define the following class attributes:
    # List of pseudos in (basenames in abipy/data/pseudos
    pseudos = ["14si.pspnc"]

    # Mapping old_name --> new_name for the output files that must be preserved.
    files_to_save = {
        "out_DS1_DEN-etsf.nc": "si_DEN-etsf.nc",
        "out_DS2_GSR.nc": "si_nscf_GSR.nc",
        "out_DS2_WFK_0-etsf.nc": "si_nscf_WFK-etsf.nc",
        "out_DS1_GSR.nc": "si_scf_GSR.nc",
        "out_DS1_WFK_0-etsf.nc": "si_scf_WFK-etsf.nc",
    }


def main():
    return MyGenerator().run()

if __name__ == "__main__":
    sys.exit(main())
