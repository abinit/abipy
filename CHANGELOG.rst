.. _whats-new-0-3:

v0.3
----

* Add :ref:`abiview.py` script for a quick visualization of results.
* Fix bug in in Sigres.get_dataframe if QP states do not start with the same band index.
* Add ``mpirunner_options`` and ``shell_runner_options`` to TaskManager.
* Autodetect presence in the DDB_ of data required for the LO-TO splitting.
* DONE Solve problem with visualize in jupyter notebooks (files should be produced in workdir)
* DONE Change shifts default value in g0w0_with_ppmodel_inputs.
* DONE Finalize interface with phononwebsite_.
* Add robots for comparing/analyzing multiple files of the same type (DdbRobot, GsrRobot ...)
  Some of the robot capabilities are exposed via the :ref:`abicomp.py` and the :ref:`abirun.py` scripts.
* Add several new options to :ref:`abirun.py`, :ref:`abicomp.py`, :ref:`abistruct.py` scripts.
* Significant improvements to the documentation and the website: add :ref:`plot-gallery` with matplotlib plots
  and :ref:`flow-gallery` with AbiPy flows are now automatically generated.

.. _whats-new-0-2:

v0.2
----
Mar 10 2017

This is the first official release in which we have reached a relatively stable API
and a well-defined interface with the netcdf files produced by Abinit.
We recommend Abinit >= 8.0.8b, version 8.2.2 is required to analyze the electronic fatbands
saved in the FATBANDS.nc_ file.
