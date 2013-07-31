# This example shows how to plot the phonon band structure with markers.
from abipy import *

# Create the object from file.
phbands = PhononBands.from_file(get_reference_file("trf2_5.out_PHBST.nc"))

# Create the marker. Here we just use the phonon frequency as size of the marker.
x = []
for q in range(phbands.num_qpoints):
    x.extend(phbands.num_branches * [q])

xys = [ x,
        phbands.phfreqs.ravel(),
        phbands.phfreqs.ravel(),
]

phbands.set_marker("fake marker", xys)

# Plot the phonon frequencies. Note that the labels for the q-points
# are found automatically by searching in an internal database.
phbands.plot(title="AlAs phonon band structure with omega(q, nu) as marker", marker="fake marker:10000")
