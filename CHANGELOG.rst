* :release:`0.3.0 <2017-12-26>`
* :feature:`0` Add :ref:`abiview.py` script for a quick visualization of results.
* :feature:`0` :ref:`abicheck.py` accepts ``-with-flow`` options 
* :feature:`0` Add AbinitInput.set_spell_check to activate/deactivate spell-checker
* :feature:`0` Improve coverage
* :bug:`0 major` Fix bug in in `SigresFile.get_dataframe` if QP states do not start with the same band index.
* :bug:`0 major` Fix bug in thermodinamical properties (zero-point energy was included twice)
* :feature:`0` Add ``mpirunner_options`` and ``shell_runner_options`` to |TaskManager|.
* :feature:`0` Autodetect presence in |DdbFile| of data required for the LO-TO splitting.
* :feature:`0` Solve problem with visualize in jupyter_ notebooks (files should be produced in workdir)
* :feature:`0` Change default value of ``shifts`` in :func:`abipy.abio.factories.g0w0_with_ppmodel_inputs`.
* :feature:`0` Add interface with phononwebsite_: `abiview.py phbands out_PHBST.nc -web`.
* :feature:`0` Add robots for comparing/analyzing multiple files of the same type (|DdbRobot|, |GsrRobot| ...)
  Some of the robot capabilities are exposed via the :ref:`abicomp.py` and the :ref:`abirun.py` scripts.
* :feature:`0` Add several new options to :ref:`abirun.py`, :ref:`abicomp.py` and :ref:`abistruct.py` scripts.
* :support:`0` Significant improvements to the documentation and the website: add :ref:`plot-gallery` with matplotlib plots
  and :ref:`flow-gallery` with AbiPy flows are now automatically generated.
* :feature:`0`: Add Shankland-Koelling-Wood Fourier interpolation scheme.
  For the theoretical background see :cite:`Euwema1969,Koelling1986,Pickett1988,Madsen2006`.
* :release:`0.2.0 <2017-03-10>`
* :feature:`0` This is the first official release in which we have reached a relatively stable API
  and a well-defined interface with the netcdf files produced by Abinit.
  We recommend Abinit >= 8.0.8b, version 8.2.2 is required to analyze the electronic fatbands
  saved in the FATBANDS.nc_ file.
