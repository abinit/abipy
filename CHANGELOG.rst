Release 0.8.0: xxxx-xx-xx

    * Add abicheck.py --create-config option to install predefined yaml configuration files
    * Add support for NSCF calculations with meta-GGA.
    * Preliminary support for panel dashboards exposed via `abiopen FILE --panel` and `abistruct panel FILE`.
      Note that not all Abinit files are supported at present.
    * Add examples and flows for effective mass calculations
    * Add examples for quasi-harmonic calculations and post-processing tools
    * Add support for JSON files (including MSONable format) to abiopen.py
      Supports `--notebook`, `--panel` options such as `abiopen.py FILE.json --panel`
    * Improved support for EPH calculations.
    * Add `primitive` command to `abistruct.py` to get primitive structure from spglib

Release 0.7.0: 2019-10-18

    * Remove support for py2. Now Abipy requires py >= 3.6 (3.8 is not yet supported)
    * AbiPy now requires pymatgen >= 2019.10.16
    * Move workflow code from pymatgen to abipy.flowtk
    * Improved support for EPH calculations.

Release:0.3.0 2017-12-26

    * Add ``abiview.py`` script for a quick visualization of results.
    * ``abicheck.py`` accepts ``-with-flow`` option
    * Add AbinitInput.set_spell_check to activate/deactivate spell-checker
    * Improve coverage
    * Fix bug in in ``SigresFile.get_dataframe`` if QP states do not start with the same band index.
    * Fix bug in thermodinamical properties (zero-point energy was included twice)
    * Add ``mpirunner_options`` and ``shell_runner_options`` to TaskManager.
    * Autodetect presence in DdbFile of data required for the LO-TO splitting.
    * Solve problem with visualize in jupyter notebooks (files should be produced in workdir)
    * Change default value of ``shifts`` in ``abipy.abio.factories.g0w0_with_ppmodel_inputs``.
    * Add interface with phononwebsite: ``abiview.py phbands out_PHBST.nc -web``.
    * Add robots for comparing/analyzing multiple files of the same type (``DdbRobot``, ``GsrRobot`` ...)
      Some of the robot capabilities are exposed via the ``abicomp.py`` and the ``abirun.py`` scripts.
    * Add several new options to ``abirun.py``, ``abicomp.py`` and ``abistruct.py`` scripts.
    * Significant improvements to the documentation and the website: add ``plot-gallery`` with matplotlib plots
      and ``flow-gallery`` with AbiPy flows are now automatically generated.
    * Add Shankland-Koelling-Wood Fourier interpolation scheme.

Release 0.2.0 <2017-03-10>

    This is the first official release in which we have reached a relatively stable API
    and a well-defined interface with the netcdf files produced by Abinit.
    We recommend Abinit >= 8.0.8b, version 8.2.2 is required to analyze the electronic fatbands
    saved in the FATBANDS.nc file.
