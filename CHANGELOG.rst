
Release 0.9.3:

    * Require pymatgen >= 2023.3.23
    * Add ``ndivsm`` argument to PhononWork and EffMassLineWork to activate band 
      structure computation along k-path.
    * Add new option to TaskManager (`limits_for_task_class`) to specify custom limits 
      depending on the name of the task class. See `abidoc.py manager` for syntax
    * Deprecate `single` and `rapid` commands of abirun.py.


Release 0.9.2:

   * Require pymatgen >= 2022.0.14
   * Add abiview.py ifermi_fs
   * G0W0WithQptdmFlow is deprecated and will be removed in v 0.9.3
   * Add abigui.py script
   * Add abilab.abirobot function to create a Robot from a list of filepaths.
   * Add plotly version of fatbands
   * Support for dipquad and quadquad

Release 0.9.1:

   * Add  "-ew", "--expose-web",
   * abiopen, abiview and abicomp now supports --plotly to produce plotly figures in the local browser
     and --chart-studio to push the figure to the chart studio service.
     Note that, at present, only DDB files support plotly.
   * AbinitInput set_kpath and make_ebands_input now support negative values of ndivsm that
     are interpreted as line_density following pymatgen conventions.
     This option is the recommended one if the k-path contains two consecutive high symmetry k-points
     that are very close as ndivsm > 0 may produce a very large number of wavevectors.
   * Preliminary support for plotly plots (phonons).
   * AbinitInputParser now can parse strings in the input file and read structure is the `structure:abivars`
     syntax in used.

Release 0.9.0:

    * Require pymatgen >= 2019.12.22
    * Integration with abinit 9.4.0
    * Support for py3.8 and py3.9
    * Converters for phonopy and tdep.
    * Preliminary support for panel dashboards.
    * Introduce ABIPY_TMPDIR environment variable to specify the workdir of the temporary tasks at runtime.
    * Use last version of apscheduler.
    * Minor bug fixes
    * New tools for Phonon and EPH calculations.
    * Note that this is the last AbiPy version supporting Abinit8.
      AbiPy version 1.0 will start to take advantage of features and ameliorations introduced in Abinit9
      We will also take the opportunity to refactor the code base so backward incompatibe changes in the API
      are expected in the next major version.

Release 0.8.0:

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

Release 0.2.0 2017-03-10

    This is the first official release in which we have reached a relatively stable API
    and a well-defined interface with the netcdf files produced by Abinit.
    We recommend Abinit >= 8.0.8b, version 8.2.2 is required to analyze the electronic fatbands
    saved in the FATBANDS.nc file.
