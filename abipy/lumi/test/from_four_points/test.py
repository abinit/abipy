#!/usr/bin/env python
import sys
from abipy.core import Structure
from abipy.lumi.deltaSCF import DeltaSCF
import matplotlib.pyplot as plt
import numpy as np


files=('Ag_file/out_GSR.nc',
       'Agstar_file/out_GSR.nc',
       'Aestar_file/out_GSR.nc',
       'Ae_file/out_GSR.nc',)

scf_object=DeltaSCF.from_four_points_file(filepaths=files)
df=scf_object.get_dataframe_element()
#scf_object.displacements_visu()
#plt.show()


def main():
    return df

if __name__ == '__main__':
    sys.exit(main())
