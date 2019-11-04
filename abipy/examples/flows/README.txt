.. _flow-gallery:

Flow Gallery
============

This gallery contains python scripts to generate AbiPy flows from the command line.

Run the scripts to generate the directory with the flow, then use the :ref:`abirun.py` script to execute the flow.
Alternatively, one can use the ``-s`` option to generate the flow and run it immediately with the scheduler.
Use ``--help`` for further information on the available options.

.. warning::

    The following examples show how to use python and the AbiPy API to generate and run
    Abinit calculations in a semi-automatic way.
    These examples are not supposed to produce physically meaningful results
    as input parameters are usually underconverged.
    Note also that the figures can only show the initial configuration of the Flow.
    Additional Works generated at runtime won't be displayed.
    To visualize the entire Flow, you need to run the script and then use::

        abirun.py FLOWDIR graphviz

    where `FLOWDIR` is the directory of the Flow.
