.. :Release: |version| :Date: |today|

.. _scripts-index:

=======
Scripts
=======

This page documents the usage of the AbiPy scripts,
the subcommands available and the options supported by each subcommand.

To analyze the cystalline structure stored in FILE, use ``abistruct.py``.
To operate on a **single** FILE, use ``abiopen.py``.
To compare **multiple** FILES of the same type, use the ``abicomp.py`` script.
If the analysis requires the execution of additional logic
(e.g. the computation of phonons with anaddb from the DDB file), use ``abiview.py``.
To generate a minimalist input file for Abinit, use ``abinp.py``.
For a command line interface to the Abinit documentation, use ``abidoc.py``.

Finally, use ``abicheck.py`` to validate your AbiPy + Abinit installation **before running** AbiPy flows
and use ``abirun.py`` to launch Abinit calculations.

.. important::

    Each script provides a ``--help`` option that documents all the commands available
    and provides a list of typical examples.
    To list of the options supported by **COMMAND** use e.g. `abicomp.py COMMAND --help`.


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
