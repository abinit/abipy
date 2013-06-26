# This example shows how to plot the macroscopic 
# dielectric function (MDF) computed in the Bethe-Salpeter code
from abipy import *

# Path to the MDF file produced in the tutorial.
mdf_file = get_reference_file("tbs_4o_DS2_MDF.nc")

# Open the file and extract data.
with MDF_Reader(mdf_file) as r:

    # Build the plotter.
    plotter = MDF_Plotter()
    
    # Excitonic MDF.
    exc_mdf = r.read_exc_mdf()
    plotter.add_mdf("EXC", exc_mdf)
        
    # KS-RPA MDF
    rpanlf_mdf = r.read_rpanlf_mdf()
    plotter.add_mdf("KS-RPA", rpanlf_mdf)

    # GW-RPA MDF (obtained with the scissors operator).
    gwnlf_mdf = r.read_gwnlf_mdf()
    plotter.add_mdf("GW-RPA", gwnlf_mdf)

    # Plot spectra in the range 2-5 eV.
    title = "Si absorption spectrum: EXC vs RPA"
    plotter.plot(title=title, xlim=(2, 5))
