#!/usr/bin/env python
from __future__ import print_function
import sys

from abipy.data import AnaddbFilesGenerator

class MyGenerator(AbinitFilesGenerator):
    """This class generates the output files used in the unit tests and in the examples."""

    # Mapping old_name --> new_name for the output files that must be preserved.
    files_to_save = {
        "trf2_5.out_PHBST.nc": "trf2_5.out_PHBST.nc",
        "trf2_5.out_PHDOS.nc": "trf2_5.out_PHDOS.nc",
    }


def main():
    return MyGenerator().run()

if __name__ == "__main__":
    sys.exit(main())
