#!/usr/bin/env python
from __future__ import print_function, unicode_literals, division, absolute_import
import sys

from abipy.data import AbinitFilesGenerator


class MyGenerator(AbinitFilesGenerator):
    """This class generates the output files used in the unit tests and in the examples."""
    # Subclasses must define the following class attributes:
    # List of pseudos (basenames in abipy/data/pseudos)
    pseudos = ["14si.pspnc"]

    # Mapping old_name --> new_name for the output files that must be preserved.
    files_to_save = {
        "out_DS4_SIGRES.nc": "si_g0w0ppm_nband10_SIGRES.nc",
        "out_DS5_SIGRES.nc": "si_g0w0ppm_nband20_SIGRES.nc",
        "out_DS6_SIGRES.nc": "si_g0w0ppm_nband30_SIGRES.nc",
    }


if __name__ == "__main__":
    sys.exit(MyGenerator().run())
