.. _abiopen:

^^^^^^^^^^^^^^
``abiopen.py``
^^^^^^^^^^^^^^

AbiPy provides python objects associated to several Abinit output files.
These objects implement methods to analyze and plot the results.
The examples in our :doc:`gallery ../gallery>` use this API to plot data with ``matplotlib``.

The ``abiopen.py`` script provides a handy interface to the AbiPy objects.
It can be used to open Abinit files directly in the ``ipython`` shell or in a ``jupyter`` notebook and interact with
the associated object (called ``abifile`` in the ``ipython`` terminal).
The syntax of the script is::

    $ abiopen.py FILE [options]

where ``FILE`` is one of the files supported by AbiPy (usually in ``netcdf`` format but other 
files are supported as well).
By default ``abiopen`` starts the ``ipython`` terminal.
Alternatively, it is possible to generate automatically a ``jupyter`` notebook with the ``-nb`` option e.g.::

    $ abiopen.py out_FATBANDS.nc -nb

which will produce a notebook to visualize the electronic fatbands produced with ``prtdos 3`` inside a web browser.

Use the ``-p`` option if you just want to get information on the file without opening it, e.g.::

    $ abiopen.py out_GSR.nc -p

``abiopen.py`` uses the file extension to decide what to do with the file and the type
of python object that should be instantiated.
The list of supported file extensions is obtained with:

.. command-output:: abiopen.py --help

.. WARNING::

    AbiPy uses the ``.abi`` extension for Abinit input files, ``.abo`` for output files and ``.log`` for log files.
    Try to follow this convention in your calculations to facilitate the integration with AbiPy.

.. argparse::
   :ref: abipy.scripts.abiopen.get_parser
   :prog: abiopen.py
