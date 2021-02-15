#!/usr/bin/env python
from __future__ import print_function, unicode_literals, division, absolute_import
import sys

from abipy.data import AbinitFilesGenerator


class MyGenerator(AbinitFilesGenerator):
    """This class generates the output files used in the unit tests and in the examples."""
    # Subclasses must define the following class attributes:
    # List of pseudos (basenames) in abipy/data/pseudos
    pseudos = ["12mg.pspnc", "B.psp8"]
    #pseudos = ["B.psp8", "12mg.pspnc"]


    # Mapping old_name --> new_name for the output files that must be preserved.
    files_to_save = {
        "out_DS1_FATBANDS.nc": "mgb2_kmesh181818_FATBANDS.nc",
        "out_DS2_FATBANDS.nc": "mgb2_kpath_FATBANDS.nc",
    }


if __name__ == "__main__":
    sys.exit(MyGenerator().run())
