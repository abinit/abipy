# This example shows how to plot the phonon band structure of AlAs.
# See tutorial/lesson_rf2.html

# FIXME: LO-TO splitting and phonon displacements instead of eigenvectors.
from abipy.tests import get_reference_file
from abipy.phonons import PhononBands, PHDOS_Reader, PHDOS_File

# Path to the PHBST file produced by anaddb.
phbst_file = get_reference_file("trf2_5.out_PHBST.nc")

# Create the object from file.
phbands = PhononBands.from_file(phbst_file)

# Read the Phonon DOS from the netcd file produced by anaddb (prtdos 2)
phdos_file = get_reference_file("trf2_5.out_PHDOS.nc")

with PHDOS_Reader(phdos_file) as ncdata:
    dos = ncdata.read_phdos()

# plot phonon bands and DOS.
phbands.plot_with_dos(dos, title="AlAs Phonon bands and DOS")
