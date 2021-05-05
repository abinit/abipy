#!/usr/bin/env python
import sys
from abipy.core import Structure
from abipy.lumi import DeltaSCF
import matplotlib.pyplot as plt
import numpy as np

gs_file='gs_file/out_GSR.nc'
ex_file='ex_file/out_GSR.nc'

scf_object=DeltaSCF.from_four_points_file(filepaths=(gs_file,ex_file))


#def main():
#    return scf_object.

#if __name__ == '__main__':
#    sys.exit(main())
