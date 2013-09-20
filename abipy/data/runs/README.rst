This directory contains the scripts used to generate most of the Netcdf files
employed in the unit tests and in the Abipy examples.
Using an automated procedure for the generation of the netcdf files has two distinct
advantages:

1) It is possible to check whether abipy is compatible with the different versions 
   of Abinit.

2) The automated generation allows one to test and validate the python interface 
   used for the generation of input files and high throughput calculations.

In order to facilitate the automatic execution and validation, the python scripts 
must satisfy some basic rules and conventions.

1) The name of the script must match the regular expression:
  run_[*].py
