#!/usr/bin/env python
#
# This example shows how to compute and plot the derivatives of the
# KS eigenvalues along a high symmetry path in K-space.
from abipy.abilab import abiopen
import abipy.data as data

gsr_file = abiopen(data.ref_file("si_nscf_GSR.nc"))

kpath = gsr_file.kpoints

#print("ds", kpath.ds)
#for vers in kpath.versors:
#    print("versors", vers)

#print(kpath)
#print(kpath.lines)

#branch = bands.get_branch(spin=0, band=3)
#ders = kpath.finite_diff(branch, order=1)
#ders = kpath.finite_diff(branch, order=2)
#print("order2",ders)

bands = gsr_file.ebands

xys = bands.derivatives(spin=0,band=0,order=1, asmarker="DER1-band0")

xys = bands.derivatives(spin=0,band=1,order=1, asmarker="DER1-band1")

xys = bands.derivatives(spin=0,band=2,order=1, asmarker="DER1-band2")

xys = bands.derivatives(spin=0,band=3,order=1, asmarker="DER1-band3")

bands.plot(marker="DER1-band0:100")

emasses = bands.effective_masses(spin=0,band=0)
print("emasses", emasses)

emasses = bands.effective_masses(spin=0,band=3)
print("emasses", emasses)

emasses = bands.effective_masses(spin=0,band=4)
print("emasses", emasses)
