#!/usr/bin/env python
import sys
from abipy.core import Structure
from abipy.lumi.deltaSCF import DeltaSCF
import matplotlib.pyplot as plt
import numpy as np

gs_file='gs_file/out_GSR.nc'
ex_file='ex_file/out_GSR.nc'

scf_object=DeltaSCF.from_relax_file(filepaths=(gs_file,ex_file))


def main():
    return scf_object.delta_r()

if __name__ == '__main__':
    sys.exit(main())
