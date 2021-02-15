.. _abiopen.py:

^^^^^^^^^^^^^^
``abiopen.py``
^^^^^^^^^^^^^^

AbiPy provides python objects associated to several Abinit output files.
These objects implement methods to analyze and plot the results.
The examples in our :ref:`plot-gallery` use this API to plot data with matplotlib_.

The ``abiopen.py`` script provides a handy interface to the AbiPy objects.
It can be used to open Abinit files directly in the ipython_ shell or inside a jupyter_ 
notebook and interact with the associated object (called ``abifile`` in the ``ipython`` terminal).
The syntax of the script is::

    abiopen.py FILE [options]

where ``FILE`` is one of the files supported by AbiPy (usually in netcdf_ format but other 
files are supported as well e.g. Abinit input and output files in text format).
By default ``abiopen`` starts an ``ipython`` session and the user can interact with the ``abifile``
and invoke its methods.

Alternatively, it is possible to generate automatically a jupyter_ notebook with the ``-nb`` option e.g.::

    abiopen.py out_FATBANDS.nc -nb

will produce a notebook to visualize the electronic fatbands inside your default web browser.

Use the ``-p`` option if you just want to get information on the file without opening it, e.g.::

    abiopen.py out_GSR.nc -p

or the ``-e`` (``--expose``) to generate matplotlib plots automatically::

    abiopen.py out_GSR.nc -e -sns=talk

seaborn_ plot style and settings can be changed from the command line interface with the `-sns` option 

The script uses the file extension to decide what to do with the file and the type
of python object that should be instantiated.
The list of supported file extensions is obtained with:

.. command-output:: abiopen.py --help

.. WARNING::

    AbiPy uses the ``.abi`` extension for Abinit input files, ``.abo`` for output files and ``.log`` for log files.
    Please follow this convention to facilitate the integration with AbiPy.

Complete command line reference

.. argparse::
   :ref: abipy.scripts.abiopen.get_parser
   :prog: abiopen.py
