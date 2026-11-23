.. _abiopen-script:

^^^^^^^^^^^^^^
``abiopen.py``
^^^^^^^^^^^^^^

AbiPy provides Python objects associated with several Abinit output files.
These objects implement methods to analyze and plot the results.
The examples in our :ref:`plot-gallery` use this API to plot data with matplotlib_.

The ``abiopen.py`` script provides a handy interface to these objects.
It opens Abinit files directly in the ipython_ shell or in a jupyter_
notebook, where you can interact with the associated object (called ``abifile`` in the ``ipython`` terminal).
The syntax of the script is::

    abiopen.py FILE [options]

where ``FILE`` is one of the files supported by AbiPy (usually in netcdf_ format, although other
files are supported as well, e.g. Abinit input and output files in text format).
By default, ``abiopen`` starts an ``ipython`` session in which you can interact with the ``abifile``
and invoke its methods.

Alternatively, the ``-nb`` option automatically generates a jupyter_ notebook. For example::

    abiopen.py out_FATBANDS.nc -nb

produces a notebook to visualize the electronic fatbands in your default web browser.

Use the ``-p`` option if you just want to print information about the file without opening an ipython session, e.g.::

    abiopen.py out_GSR.nc -p

or the ``-e`` (``--expose``) option to generate matplotlib plots automatically::

    abiopen.py out_GSR.nc -e -sns=talk

The seaborn_ plot style and settings can be changed from the command line with the `-sns` option.

The script uses the file extension to decide what to do with the file and which
Python object to instantiate.
To get the list of supported file extensions, use:

.. command-output:: abiopen.py --help

.. WARNING::

    AbiPy uses the ``.abi`` extension for Abinit input files, ``.abo`` for output files and ``.log`` for log files.
    Please follow this convention to ease the integration with AbiPy.

Complete command line reference:

.. argparse::
   :ref: abipy.scripts.abiopen.get_parser
   :prog: abiopen.py
