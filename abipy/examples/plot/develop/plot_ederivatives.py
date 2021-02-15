#!/usr/bin/env python
r"""
Effective masses
================

This example shows how to compute and plot the derivatives of the
KS eigenvalues along a high symmetry path in K-space.
"""
from abipy.abilab import abiopen
import abipy.data as abidata

gsr_file = abiopen(abidata.ref_file("si_nscf_GSR.nc"))

kpath = gsr_file.kpoints

#print("ds", kpath.ds)
#for vers in kpath.versors:
#    print("versors", vers)
print(kpath)
#print(kpath.lines)

#branch = bands.get_branch(spin=0, band=3)
#ders = kpath.finite_diff(branch, order=1)
#ders = kpath.finite_diff(branch, order=2)
#print("order2",ders)

ebands = gsr_file.ebands

xys = ebands.derivatives(spin=0, band=0, order=1)
print("xys\n", xys)
#xys = ebands.derivatives(spin=0, band=1, order=1)
#xys = ebands.derivatives(spin=0, band=2, order=1)
#xys = ebands.derivatives(spin=0, band=3, order=1)

#ebands.plot()

emasses = ebands.effective_masses(spin=0, band=0, acc=2)
print("emasses", emasses)

#emasses = ebands.effmass(spin=0, kpoint=[0,0,0], bands=0, acc=2)
#print("emasses", emasses)

#emasses = ebands.effective_masses(spin=0, band=3, acc=2)
#print("emasses", emasses)

#emasses = ebands.effective_masses(spin=0, band=4, acc=2)
#print("emasses", emasses)

gsr_file.close()
