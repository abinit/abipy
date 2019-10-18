#!/usr/bin/env python
r"""
Infrared spectrum of AlAs
=========================

This example shows how to plot the infrared spectrum of a polar semiconductor (AlAs)
from the DDB file  See tutorial/lesson_rf2.html

For a command line interfase, use:

    abiview.py ddb_ir in_DDB
"""
import os
import abipy.data as abidata

from abipy import abilab

# Open DDB file for alpha-SiO2 taken from https://materialsproject.org/materials/mp-7000/
filepath = os.path.join(abidata.dirpath, "refs", "mp-7000_DDB.bz2")
ddb = abilab.abiopen(filepath)

# Invoke anaddb to compute dielectric tensor and oscillator strength.
tgen = ddb.anaget_dielectric_tensor_generator(asr=2, chneut=1, dipdip=1, verbose=1)
print(tgen)

# Set phonon damping factor in eV (full width).
gamma_ev = 1e-3

# Plot IR spectrum in Cartesian coordinates.
tgen.plot_all(gamma_ev=gamma_ev, title="Diagonal and off-diagonal components")

tgen.plot(component="diag", reim="re", gamma_ev=gamma_ev, title="Real part, diagonal components")

tgen.plot(component="diag", reim="im", gamma_ev=gamma_ev, title="Imaginary part, diagonal components")
