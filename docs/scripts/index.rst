.. :Release: |version| :Date: |today|

.. _scripts-index:

=======
Scripts
=======

This page documents the AbiPy scripts,
their subcommands and the options supported by each subcommand.

To analyze the crystalline structure stored in FILE, use ``abistruct.py``.
To operate on a **single** FILE, use ``abiopen.py``.
To compare **multiple** FILES of the same type, use ``abicomp.py``.
If the analysis requires additional steps
(e.g. computing phonons with anaddb from a DDB file), use ``abiview.py``.
To generate a minimal Abinit input file, use ``abinp.py``.
For a command line interface to the Abinit documentation, use ``abidoc.py``.

Finally, use ``abicheck.py`` to validate your AbiPy + Abinit installation **before running** AbiPy flows,
and ``abirun.py`` to launch Abinit calculations.

.. important::

    Each script provides a ``--help`` option that documents all the available commands
    and gives a list of typical examples.
    To list the options supported by a **COMMAND**, use e.g. `abicomp.py COMMAND --help`.


.. toctree::
   :maxdepth: 3

   abistruct.rst
   abiopen.rst
   abicomp.rst
   abiview.rst
   abinp.rst
   abidoc.rst
   abicheck.rst
   abirun.rst
   abips.rst
   oncv.rst
