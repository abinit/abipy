TODO List
=========

#. Rationalize the different readers defined in pymatgen and abipy.
   Try to integrate some method of abipy.core.Structure in pymatgen.core.Structure. 
   See the hack used in iotools.__init__.py

#. Move unit conversion to NetcdfReader e.g reader.read_value(varname, unit=None)

#. Use different and cleaner rules for file extensions in ABINIT. Why _DEN12 and _1WF13 instead
   of the simpler syntax 12_1DEN, 13_1WF in which the extension is preserved?

#. split long lines in the abinit input (e.g. typat 1 1 1 2 --> typat 3*1 2)

#. Better support for PBS
