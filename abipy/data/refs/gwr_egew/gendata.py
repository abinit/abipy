#!/usr/bin/env python
import sys

from abipy.data import AbinitFilesGenerator


class MyGenerator(AbinitFilesGenerator):
    """This class generates the output files used in the unit tests and in the examples."""

    # Subclasses must define the following class attributes:
    # List of pseudos (basenames in abipy/data/pseudos)
    pseudos = ["Ga-low_r.psp8", "As_r.psp8"]

    # Mapping old_name --> new_name for the output files that must be preserved.
    files_to_save = {
        "out_DS3_GWR.nc": "out_DS3_GWR.nc",
        "out_DS4_GSR.nc": "out_DS4_GSR.nc",
    }

    mpiexec = "mpirun"
    mpinp = "8"


if __name__ == "__main__":
    sys.exit(MyGenerator().run())
