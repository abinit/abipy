TODO list:

## High priority

* DONE Get rid of readthedocs

* Reorganize modules in flowtk to prepare future migration. Modules with gs_works, dfpt_works ...
  qadapter package ... (postponed to v0.7)

* Use angdeg instead of rprimd in structure_to_abivars if hex or rhomboedral lattice 
  (tricky because input settings should be preserved)

* introduce new status for tasks that are removed at runtime e.g. S_CANCELLED
  and handle new case in flow machinery. Be careful with pickle, status comparison and ordering though.

* introduce new status WAITING_FOR_RESTART
  so that we don't have to restart task in callbacks

* Fix annoying warnings about k-point sampling.

* DONE Reintegrate AbiPy with new abivars (cleanup?)

* Check Positive gw_qprange in EPH (Fixed by Henrique)

* DONE abicomp should accept tolsym args

* Add support for PSML/UPF format

* Add iscf to GSR.nc so that we know if we have SCF|NSCF run.

* Improve exception handling in NetcdfReader

* Read forces in read_structure ?

* Automate CHANGELOG creation.

* Fix DFPT with iomode 3 (make_links logic)

* Refactor S_QCRITICAL logic (logic injected by user, since qcritical errors are cluster-specific)

* Refactor wrappers for mrgddb and mrgdvdb (problems with subprocess when
  merging large number of partial files (likely due to Popen with large stderr/stdout)

## Medium priority

* remove phononflow

* Add DOS to GSR file (useful if tetra)  Create Dosfile ? Fortran exec?

* videos in README (atom and hydrogen)

* ALMOST DONE: Fix travis warnings.

* Refactor/improve Visualizer

* Read LO-TO data from PHBST.nc instead of anaddb.nc (not easy as directions should be computed by AbiPy)

* add possibility of changing amu in anaddb/abinit and API to "mix" DDB files
  phonon group velocities (requires extension in netcdf files).

* DONE Solve problem with visualize in jupyter notebooks (files should be produced in workdir)

* Scheduler should report info on exceptions (especially if at the end when on_all_ok is invoked)

* ALMOST DONE: Replace core.tensor with pymatgen tensor
  DONE Use pmg tensor for stress as well.
  Check DielectricTensor in Anaddb from DDB.

* Add nsppol, nspinor, nspden to HIST file (and other stuff?)

* Fix bug with SCGW and SKW interpolation reported by Ahn.

* Optimize SKW (slow if dense IBZ). Add possibility of initializing SKW
  from nc file produced by Fortran version.

* Add integration test for dilatmx error handler

* Add ExpectedAbort to Abinit so that one can call the code to get data without triggering
  report.errors in task.get_event.report

* Add memory error to Abinit errors

* Investigate NaN issue in BECS reported by Ahn if tolvrs instead of tolwfr (tolwfr could activate nbdbuf)

* DONE Check infra-red dielectric function from DDB.

* Add input file to NC files (?)

* Add phonon plot with Longitudinal/transverse character and Z q 

## Low priority

* Rationalze wrappers for mrgdddb .... (raise exception in python if clear error, retcode 
  is the returncode of the script not necessarily the retcode of the exe, need to
  parse log file and make sure that all scripts write log files in "abinit" format
  that can be read with EventsParser.

* Refactor PyLauncher logic

* Add python API to support discontinuous paths (Abinit is not able to handle that
  but python code should be agnostic

* Finalize DDK.nc (EVK.nc)

* Fix issue with DOJO_REPORT and PAW XML files.

* DONE plot_networkx does not work with flows containing callbacks e.g. run_qptdm_flow
  FIXED with graphviz

* Check xsf_write_data and visualization of potentials.

* Add phbands.to_bxsf and histogram for phonon modes at a given q-point.
  overlap matrix for displacements?

* Add possibility of specifying the max number of CPUs that can be used  
  for a flow at the level of the scheduler.

* Fix problem with AbiniEvent format, src_file and scr_line (see src/67_common/scprqt.F90)
  Introduce an integer flag (msg_level) to be passed to msg_hndl

* ABINIT abort file should not be produced if the exit is expected otherwise we 
  can have IO race conditions and ABI_CRITICAL events!!!!!!!

* Add option max_num_launchers in scheduler.yml

* Add extra metadata to netcdf files (try to propagate info on space group from parser to crystal_t
  as well as Abinit input as string)

* Improvement in the dilatmx error handler:

        [30/03/15 15:15:58] guido petretto: cmq ci sarebbe un'altra cosa che non so se avevi già considerato. 
         Questo non porta a errori, ma non so se è il
        modo più corretto di gestire la cosa, anche perché non sono sicuro di cosa faaccia abinit esattamente
        [30/03/15 15:16:10] guido petretto: esempio:
        [30/03/15 15:17:12] guido petretto: - calcolo relax con ntime basso
        [30/03/15 15:18:23] guido petretto: - viene stoppato a un certo punto e riavviato -> vengono settate le variabili tipo irdden
        [30/03/15 15:18:54] guido petretto: - dopo il restart c'è un errore (dilatm?) e viene chiamato il reset_from_scraatch
        [30/03/15 15:19:48] guido petretto: a questo punto il job viene riavviato, ma in in c'è ancora la vecchia density e nell'input c'è irdden=1, ma la
        struttura è diversa

* FFTProf (use file extension and interface it with abiopen)

* Create new github package for benchmarks/Abinit tests + template for new python projects.

* Remove GUI code.

* nbjsmol (build system, refactor API?)

* fatbands with SOC (waiting for Matthieu's refactoring)

* integrate improvements in skw by Nicholas.
  Finalize baseclass for ElectronInterpolator

* ALMOST DONE lobster interface from Guido

* context manager to change variables (e.g. autoparal)

* Cleanup and refactoring in OpticTask

* Replace SIGRES with new fileformat based on SIGEPH (long-term)

* Update spack recipe, add support for EasyBuild, revamp homebrew (?)

* Classification of phonons/electrons

* Error handler for tolwfr to increase nband / nbdduf and resubmit
