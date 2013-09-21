This directory contains the scripts used to generate most of the Netcdf files
employed in the unit tests and in the Abipy examples.

Using an automated procedure for the generation of the netcdf files has two distinct
advantages:

1) It is possible to check whether abipy is compatible with the different versions 
   of Abinit.

2) The automated generation allows one to test and validate the python interface 
   used for the generation of the input files and the API for high throughput calculations.

In order to facilitate the automatic execution and validation, the python scripts 
must satisfy the following rules and conventions.

#. The name of the script must match the regular expression run_[*].py so that 
   we can run all the tests easily with run_all.py

#. The execution of the test should be managed by a `Tester` object
   The `Tester` is responsible for the definition of the working directory (constructed
   from the name of the script by just removing the prefix `run_`), the submission
   of the calculation (tester.set_work_and_run) and the analysis of the final results
   (tester.finalize).

#. The script should remove all the output files produced by the run that are not needed 
   for the automatic tests and/or the tutorials. Each file should have a unique (meaningful) name 
   so that we can easily access it with the syntax:

        import abipy.data as data
        path_to_reference_file = data.ref_file("basename_of_the_file")

   An exception is raised during the automatic tests if this rule is not respected.
