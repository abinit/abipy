# This example shows how to plot the phonon band structure of AlAs.
# See tutorial/lesson_rf2.html

from abipy.tests import get_reference_file
from abipy.phonons import PhononBands

# Path to the PHBST file produced by anaddb.
filename = get_reference_file("trf2_5.out_PHBST.nc")

# Create the object from file.
phbands = PhononBands.from_file(filename)

# Mapping reduced coordinates -> labels
qlabels = {
    (0,0,0): "$\Gamma$",
    (0.375, 0.375, 0.7500): "K",
    (0.5, 0.5, 1.0): "X",
    (0.5, 0.5, 0.5): "L",
    (0.5, 0.0, 0.5): "X",
    (0.5, 0.25, 0.75): "W",
}

# Plot the phonon band structure.
phbands.plot(title="AlAs without LO-TO splitting", qlabels=qlabels)
