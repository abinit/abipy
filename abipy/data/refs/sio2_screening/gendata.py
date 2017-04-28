#!/usr/bin/env python
from __future__ import print_function, unicode_literals, division, absolute_import
import sys

from abipy.data import AbinitFilesGenerator


class MyGenerator(AbinitFilesGenerator):
    """This class generates the output files used in the unit tests and in the examples."""
    # Subclasses must define the following class attributes:
    # List of pseudos (basenames) in abipy/data/pseudos
    pseudos = ["14si.pspnc", "8o.pspnc"]

    # Mapping old_name --> new_name for the output files that must be preserved.
    files_to_save = {
        "out_DS1_DEN.nc": "sio2_DEN.nc",
        "out_DS2_GSR.nc": "sio2_kpath_GSR.nc",
        "out_DS3_SCR.nc": "sio2_SCR.nc",
        "out_DS3_EM1_NLF": "sio2_EM1_NLF",
        "out_DS3_EM1_LF": "sio2_EM1_LF",
        "out_DS3_EELF": "sio2_EELF",
    }


if __name__ == "__main__":
    sys.exit(MyGenerator().run())
