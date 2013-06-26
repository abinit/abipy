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

# Define the mapping reduced_coordinates -> label of the k-point.
klabels = {
    (0.5,  0.0,  0.0) : "L",
    (0.0,  0.0,  0.0) : "$\Gamma$",
    (0.0,  0.5,  0.5) : "X",
}

# Plot the band energies.
bands.plot(title="Silicon band structure", klabels = klabels)
