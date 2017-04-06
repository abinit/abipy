#!/usr/bin/env python
from __future__ import print_function, unicode_literals, division, absolute_import
import sys

from abipy.data import AnaddbFilesGenerator


class MyGenerator(AnaddbFilesGenerator):
    """This class generates the output files used in the unit tests and in the examples."""

    def __init__(self, **kwargs):
        super(MyGenerator, self).__init__(**kwargs)

        self.files_to_keep.add("ddb_notes")


    # Mapping old_name --> new_name for the output files that must be preserved.
    files_to_save = {
        "out_PHBST.nc": "ZnSe_hex_886.out_PHBST.nc",
        "out_PHDOS.nc": "ZnSe_hex_886.out_PHDOS.nc",
        "anaddb.nc": "ZnSe_hex_886.anaddb.nc",
    }

    in_ddb = "ZnSe_hex_qpt_DDB"


if __name__ == "__main__":
    sys.exit(MyGenerator().run())
