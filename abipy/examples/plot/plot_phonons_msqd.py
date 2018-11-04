#!/usr/bin/env python
r"""
Projected phonon DOS
====================

This example shows how to plot the generalized phonon DOS with the mean square
displacement tensor in cartesian coords and how to calculate Debye Waller factors
as a function of temperature.
See :cite:`Lee1995` for the further details about the internal implementation and
:cite:`Trueblood1996` for the different conventions used by crystallographers.
"""
from abipy import abilab
import abipy.data as abidata

# Read the Phonon DOS from the netcdf file produced by anaddb with prtdos 2
# (alternatively one can use the shell and `abiopen.py OUT_PHDOS.nc -nb`
# to open the file in a jupyter notebook.
#with abiopen(abidata.ref_file("trf2_5.out_PHDOS.nc")) as phdos_file:
#    # Plot PJDOS.
#    phdos_file.plot_pjdos_type(units="cm-1", title="AlAs type-projected phonon DOS")
#    # To have the projection along the cartesian directions (summed over atomic types)
#    phdos_file.plot_pjdos_cartdirs_type(units="Thz", stacked=True,
#            title="Type-projected ph-DOS decomposed along the three Cartesian directions.")
#    # To plot the PJDOS for all the inequivalent sites.
#    phdos_file.plot_pjdos_cartdirs_site(view="inequivalent", stacked=True)

filepath = "~/git_repos/abipy/lifetimes/alphaSiO2/run.abo_PHDOS.nc"
phdos = abilab.abiopen(filepath)

msqd_dos = phdos.msqd_dos
print(msqd_dos)
#for fmt in ("cartesian", "cif", "ustar", "beta", "B"):
for fmt in ("cartesian", "cif"):
    print("Format:", fmt)
    #df = msqd_dos.get_dataframe(temp=100, fmt=fmt)
    df = msqd_dos.get_dataframe(temp=300, view="all", select_symbols="Si", fmt=fmt)
    abilab.print_dataframe(df)

# Plot generalized phonon DOS for each inequivalent atom in the unit cell.
msqd_dos.plot()

# Plot tensor(T) for each inequivalent atom.
msqd_dos.plot_tensor()

msqd_dos.plot_uiso()

#msqd_dos.write_cif_file("foo.cif", temp=300)
#msqd_dos.vesta_open(temp=300)

#msqd_dos.plot(view="all")
#msqd_dos.plot_tensor(view="all")
#msqd_dos.plot_uiso(view="all")
#msqd_dos.plot_uiso(view="all", what="vel")
