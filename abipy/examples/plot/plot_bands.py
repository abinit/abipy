# This example shows how to plot a band structure
# using the eigenvalues stored in a netCDF WFK file produced by abinit.
from abipy import *

# Here we use one of the WFK files shipped with abipy.
# Replace filename with the path to your WFK file.
filename = get_ncfile("si_nscf_WFK-etsf.nc")

# Open the WKF file and read data. 
wfk_file = WFK_File(filename)

# Extract the band structure. 
bands = wfk_file.get_bands()

# Plot the band energies. Note that the labels for the k-points
# are found automatically by searching in an internal database.
bands.plot(title="Silicon band structure")

# Alternatively you can use the optional argument klabels 
# that defines the mapping reduced_coordinates --> name of the k-point.
klabels = {
    (0.5, 0.0, 0.0) : "L",
    (0.0, 0.0, 0.0) : "$\Gamma$",
    (0.0, 0.5, 0.5) : "X",
}

# and pass it to the plot method:
#bands.plot(title="Silicon band structure", klabels=klabels)
