#!/usr/bin/env python
from __future__ import print_function, unicode_literals, division, absolute_import
import sys

from abipy.data import AbinitFilesGenerator


class MyGenerator(AbinitFilesGenerator):
    """This class generates the output files used in the unit tests and in the examples."""
    # Subclasses must define the following class attributes:
    # List of pseudos (basenames) in abipy/data/pseudos
    pseudos = ["Ni-sp.psp8"]

    # Mapping old_name --> new_name for the output files that must be preserved.
    files_to_save = {
        "run.abo": "run.abo",
        "out_DS1_DEN.nc": "ni_666k_DEN.nc",
        "out_DS1_GSR.nc": "ni_666k_GSR.nc",
        "out_DS1_POT.nc": "ni_666k_POT.nc",
        "out_DS1_VHA.nc": "ni_666k_VHA.nc",
        "out_DS1_VHXC.nc": "ni_666k_VHXC.nc",
        "out_DS1_VXC.nc": "ni_666k_VXC.nc",
        "out_DS1_FATBANDS.nc": "ni_666k_FATBANDS.nc",
        "out_DS2_GSR.nc": "ni_kpath_GSR.nc",
        "out_DS2_FATBANDS.nc": "ni_kpath_FATBANDS.nc",
    }


if __name__ == "__main__":
    sys.exit(MyGenerator().run())
