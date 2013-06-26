# This example shows how to plot the projected phonon DOS of AlAs.
# See tutorial/lesson_rf2.html
from abipy.tests import get_reference_file
from abipy.phonons import PHDOS_File

# Read the Phonon DOS from the netcd file produced by anaddb (prtdos 2)
phdos_file = get_reference_file("trf2_5.out_PHDOS.nc")

phdos_file = PHDOS_File(phdos_file)

# Plot data.
phdos_file.plot_pjdos_type(title="AlAs type-projected phonon DOS")
