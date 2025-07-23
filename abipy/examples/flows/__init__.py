r"""
AbiPy Flow Gallery
==================

This gallery contains python scripts to generate AbiPy flows from the command line.

Quick Start
===========

To quickly run a flow, execute the script with the following command::

    ./run_si_ebands.py -rs

The `-rs` options will:

- Remove the existing flow directory if it already exists.
- Generate a new directory for the flow.
- Immediately start the scheduler.

To access the full list of options, use::

    ./run_si_ebands.py --help

Production Runs
===============

For production calculations — which may require hours or even days to complete —
it is recommended to first generate the flow directory (FLOWDIR) without starting the scheduler.
This can be done by running the script without any arguments::

    ./run_si_ebands.py

Then, launch the scheduler using::

    nohup abirun.py FLOWDIR scheduler > log 2> err &

This approach allows better control and monitoring of long-running tasks.
"""
