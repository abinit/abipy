#!/usr/bin/env python
from abipy import abilab

path = "/Users/gmatteo/git_repos/abinit_eph/build_gcc/tests/grunesein/run.abo_GRUNS.nc"
ncfile = abilab.abiopen(path)

ncfile.plot_doses(title="DOSes available in the GRUNS file.")

# Arrow up for positive values, down for negative values.
ncfile.plot_phbands_with_gruns(title="Phonon bands with markers proportional to Gruneisen parameters + DOSes")

plotter = ncfile.get_plotter()
plotter.combiboxplot()

ncfile.close()
