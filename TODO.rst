TODO List
=========

#. Fix problem with circular dependency due to abiopen 
   Rationalize the treatment of the Abinit files and 
   the extensions associated to them.

#. Rationalize the different readers defined in pymatgen and abipy.
   Try to integrate some method of abipy.core.Structure in pymatgen.core.Structure. 
   See the hack used in iotools.__init__.py

#. Make sure that the k-sampling description reported in the WFK file is equivalent
   to the one reported in the WFK file so that the two files can be used for plotting band structures 
   and other types of post-processing 

#. ecut is not reported in the GSR file. Similar problem for the k-sampling (see SIGRES.nc)

#. Move unit conversion to NetcdfReader e.g reader.read_value(varname, unit=None)

#. Write new unit tests for Xcrysden and the other visualizers (move these tools to pymatgen?)

#. Rewrite user profile from scratch. Look at matplotlib rc.
