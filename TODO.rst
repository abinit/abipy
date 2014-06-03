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

#. Add the fermi level to the DEN file (netcdf and fortran version) so that the NSCF run can read 
   it and can report this value in the final band structure.

#. Use different and cleaner rules for file extensions in ABINIT. Why _DEN12 and _1WF13 instead
   of the simpler syntax 12_1DEN, 13_1WF in which the extension is preserved?

#. split long lines in the abinit input (e.g. typat 1 1 1 2 --> typat 3*1 2)

#. Better support for PBS
