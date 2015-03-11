This document offers some suggestion on how to install the dependencies for abipy. We focus especially on the
dependencies that require external non-python libraries.

under construction


general
-------

If there are problems installing a new version of a python package, after a git pull, remover the build folder.
This may solve issues like:

ValueError: bad marshal data (unknown type code)




python
------

Indeed on a cluster where things are not working out of the box the easiest way to get it all going is to install python
from scratch in a place where you have full control. To find out if this in nesesary try the following imports:

import netCDF4
import _tkinter
import zlib

if these turn up errors it may be easier to start from scratch



netCDF4
-------

The python bindings need the netcdf4 and hdf5 libraries/headers installed. Get them from

http://www.unidata.ucar.edu/downloads/netcdf/index.jsp

and follow the install in structions in netcdf how to get zlib and hdf5 installed properly, and finally netcdf.

You cluster may provide modules for hdf5 and netcdf4.

Before installing the python bindings make sure the NETCDF4_DIR and HDF5_DIR environ mental variables point to the
folder that contains the lib and include folders.

Finally pip install netcdf4 should work.


matplotlib
----------

To make optimal use of abipy, plotting results directly, matplotlib needs a graphical backend. When installing
matplotlib it will search for various ones and displays which ones are found. If non, except of agg, are found you will
not be able to make plots. The simples on to install is tk/tcl.

get the tcl and tk tarballs from

www.tcl.tk/software/tcltk/downloads.html

configure with --prefix

run make again on python. This will produce the _tkinter module. Finally pip install matplotlib.

installing matplotlib is easier with gcc than icc....




