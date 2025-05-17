#!/usr/bin/env python
from __future__ import print_function, unicode_literals, division, absolute_import
import sys
import os

from abipy.data import AbinitFilesGenerator


class MyGenerator(AbinitFilesGenerator):
    """This class generates the output files used in the unit tests and in the examples."""
    # Subclasses must define the following class attributes:
    # List of pseudos (basenames in abipy/data/pseudos)
    pseudos = ["14si.pspnc"]

    # Mapping old_name --> new_name for the output files that must be preserved.

    files_to_save = {"tmp.abi": "tmp.abi"}

    def __init__(self, params=None):
        """Constructor.

        Args:
            params: Dictionary with the parameters to be used in the abi file.
        """

        if params is not None:
            inputfile = open("tmp.abi", "r")
            lines = inputfile.readlines()
            for i in range(len(lines)):
                for j in params.keys():
                    if lines[i].strip().startswith(j):
                        lines[i] = lines[i].replace(j, j+"  "+str(params[j]))
                        break
            inputfile.close()
            inputfile = open("run.abi", "w")
            for i in range(len(lines)):
                inputfile.write(lines[i])

        filename = "si_g0w0_ac_%s_SIGRES.nc" % ('_'.join([i+"_"+str(params[i]) for i in params.keys()]).replace('.', '').replace('-', ''))

        self.files_to_save["out_DS4_SIGRES.nc"] = filename

        MyGenerator.files_to_save[filename] = filename

        self.workdir = os.path.dirname(__file__)

        super().__init__()


if __name__ == "__main__":
    files_to_keep = ["tmp.abi"]
    for nomegasi in range(-10,-20,-2):
        MyGenerator(
            params={"nomegasi4": nomegasi}
        ).run()