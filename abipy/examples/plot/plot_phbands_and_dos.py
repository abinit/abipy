#!/usr/bin/env python
#
# This example shows how to plot the phonon band structure of AlAs.
# See tutorial/lesson_rf2.html

# FIXME: LO-TO splitting and phonon displacements instead of eigenvectors.
from abipy.phonons import PhononBands, PHDOS_Reader, PHDOS_File
import abipy.data as data

# Path to the PHBST file produced by anaddb.
phbst_file = data.ref_file("trf2_5.out_PHBST.nc")

# Create the object from file.
phbands = PhononBands.from_file(phbst_file)

# Read the Phonon DOS from the netcd file produced by anaddb (prtdos 2)
phdos_file = data.ref_file("trf2_5.out_PHDOS.nc")

with PHDOS_Reader(phdos_file) as r:
    phdos = r.read_phdos()

# plot phonon bands and DOS.
phbands.plot_with_phdos(phdos, title="AlAs Phonon bands and DOS")
