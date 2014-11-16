This directory contains the netcdf files used for the unit tests
and for the generation of the matplotlib plots.

Each directory contains:

    #. An input file `run.abi`
    #. A python script `gendata` that executes abinit to produce the output files.


Note that:
    #. The basename of the netcdf file must be unique so that we can get its 
       absolute path with abidata.ref_file(basename)
