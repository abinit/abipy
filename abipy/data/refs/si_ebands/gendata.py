#!/usr/bin/env python
from __future__ import print_function, unicode_literals, division, absolute_import
import sys

from abipy.data import AbinitFilesGenerator


class MyGenerator(AbinitFilesGenerator):
    """This class generates the output files used in the unit tests and in the examples."""
    # Subclasses must define the following class attributes:
    # List of pseudos (basenames) in abipy/data/pseudos
    pseudos = ["14si.pspnc"]

    # Mapping old_name --> new_name for the output files that must be preserved.
    files_to_save = {
        "out_DS1_DEN.nc": "si_DEN.nc",
        "out_DS2_GSR.nc": "si_nscf_GSR.nc",
        "out_DS2_WFK.nc": "si_nscf_WFK.nc",
        "out_DS1_GSR.nc": "si_scf_GSR.nc",
        "out_DS1_WFK.nc": "si_scf_WFK.nc",
    }


if __name__ == "__main__":
    sys.exit(MyGenerator().run())
