.. _analyzing_results:

=======
abiopen
=======

This script opens Abinit outputs file (usually in `netcdf` format but other files are supported as well). 
By default the script starts an interactive `ipython` session so that the use can interact with the file 
(called `abifile` in the `ipython` terminal) and have access to its methods.
Alternatively, it is possible to generate automatically a `jupyter` notebook (``-nb`` option)
to execute code and plot the results.
`abiopen.py` uses the file extension to understand what to do with the file and the kind of python object
that should be instantiated.
The list of file extensions supported can be obtained with:

.. command-output:: abiopen.py --help

.. WARNING::

    AbiPy uses `.abi` for Abinit input files and `.abo` for output files.
    Try to follow this convention for your files as well to faciliate the integration with AbiPy

=========
abistruct
=========

Script to analyze/export/visualize the crystal structure saved in the `netcdf` files produced by ABINIT.

.. command-output:: abistruct.py --help

==========
abicomp.py
==========

Script to analyze/compare results stored in multiple netcdf files.
By default the script displays the results/plots in the shell.
Use `--ipython` to start an ipython terminal or `-nb` to generate a jupyter notebook.

.. command-output:: abicomp.py --help
