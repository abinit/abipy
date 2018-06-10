TODO list:

## High priority

* Add ebands.method to fix the fermi level in semiconductors

* Send PR to phononwebsite with phonon band structures for tutorials.

* Use angdeg instead of rprimd in structure_to_abivars if hex or rhomboedral lattice 
  (tricky because input settings should be preserved)

* introduce new status for tasks that are removed at runtime e.g. S_CANCELLED
  and handle new case in flow machinery. Be careful with pickle, status comparison and ordering though.

* DONE Check PJDOS in abinit@gitlab

* DONE Add mpirun_args see e.g nic4 and mpirun --bind-to None

* DONE Re-implement max_njobs in the queue using a counter local to the Launcher.

* Fix annoying warnings about k-point sampling.

* Reorganize modules in flowtk to prepare future migration. Modules with gs_works, dfpt_works ...
  qadapter package ... (postponed to v0.4)

* DONE Reintegrate AbiPy with new abivars

* Almost DONE: Add support for https://mybinder.readthedocs.io/en/latest/sample_repos.html#conda-environment-with-environment-yml

* DONE Add https://github.com/mcmtroffaes/sphinxcontrib-bibtex

* DONE Add support for DVDV, DDB in plot_networkx (but graphviz is much better)

* Check Positive gw_qprange in EPH

* DONE Fix problem with get_edos if we don't have enough bands 

* abicomp should accept tolsym args

## Medium priority

* remove phononflow

* Add with_becs to PhononWork

* video with atom and hydrogen

* ALMOST DONE: Fix travis warnings.

* DONE: Fix sphinx warnings.

* Refactor/improve Visualizer

* Read LO-TO data from PHBST.nc instead of anaddb.nc (postponed to v0.4)

* add possibility of changing amu in anaddb/abinit and API to "mix" DDB files
  phonon group velocities (requires extension in netcdf files).

* DONE: Autodetect presence of data for lo_to_splitting in DDB.

* DONE Solve problem with visualize in jupyter notebooks (files should be produced in workdir)

* DONE Change shifts default value in g0w0_with_ppmodel_inputs

* Scheduler should report info on exceptions (especially if at the end when on_all_ok is invoked)

* Replace core.tensor with pymatgen tensor (postponed to v0.4)

* Add nsppol, nspinor, nspden to HIST file (and other stuff?)

* Fix bug with SCGW and SKW interpolation reported by Ahn.

* Add integration test for dilatmx error handler

* Add ExpectedAbort to Abinit so that one can call the code to get data without triggering
  report.errors in task.get_event.report

* Add memory error to Abinit errors

* Investigate NaN issue in BECS reported by Ahn if tolvrs instead of tolwfr (tolwfr could activate nbdbuf)

* DONE: Fix possible error reported by Henrique (NscfTask from file that tries to change ngfft)

* Add iscf to GSR.nc so that we know if we have SCF|NSCF run.

* Had to increase fermie again to get correct gap in diamond treated if FD smearing

* Create git repo for Abipy website to facilitate integration with binder + sphinx-gallery?

## Low priority

* Use parser subclass to avoid boiler plate code.

* Rationalze wrappers for mrgdddb .... (raise exception in python if clear error, retcode 
  is the returncode of the script not necessarily the retcode of the exe, need to
  parse log file and make sure that all scripts write log files in "abinit" format
  that can be read with EventsParser.

* Add support for PSML/UPF format

* Refactor PyLauncher logic

* Add python API to support discontinuous paths (Abinit is not able to handle that
  but python code should be agnostic

* Finalize DDK.nc 

* Remove abipy.core.mixis.AbinitOutNcFile (deprecated, will be removed in 0.4)

* Fix issue with DOJO_REPORT and PAW XML files.

* DONE plot_networkx does not work with flows containing callbacks e.g. run_qptdm_flow
  FIXED with graphviz

* Check xsf_write_data and visualization of potentials.

* Add phbands.to_bxsf and histogram for phonon modes at a given q-point.

* Add treatment of out-of-boundary conditions in scissors operator.

* Add possibility of specifying the max number of CPUs that can be used  
  for a flow at the level of the scheduler.

* Fix problem with AbiniEvent format, src_file and scr_line (see src/67_common/scprqt.F90)
  Introduce an integer flag (msg_level) to be passed to msg_hndl

* Try to generalize the (very nice) approach used by Guido to handle target_dilatmx
  WorkWithCondition(...) could be used to simulate a convergence study. The main 
  problem is how to include input_generators (__next__) while preserving data persistence
  with pickle! I've already done something related to this problem in FlowCallback...

* ABINIT abort file should not be produced if the exit is expected otherwise we 
  can have IO race conditions and ABI_CRITICAL events!!!!!!!

* Add option max_num_launchers in scheduler.yml

* Add extra metadata to netcdf files (try to propagate info on space group from parser to crystal_t
  as well as Abinit input as string)

* Initialize job.sh with max number of MPI procs?

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

* DONE Add the fermi level to the DEN file (netcdf and fortran version) so that the NSCF run can read 
  it and can report this value in the final band structure.

* DONE ecut is not reported in the GSR file. Similar problem for the k-sampling (see SIGRES.nc)

* FFTProf (use file extension and interface it with abiopen)

* Create new github package for benchmarks/Abinit tests + template for new python projects.

* _repr_html_ for structure and other basic objects (pymatgen/and abipy)
   matplotlib plot for structure?

* nbjsmol (build system, refactor API?)

* fatbands with SOC (waiting for Matthieu's refactoring)

* integrate improvements in skw by Nicholas.

* lobster interface from Guido

* interface abistruct with marylin

* Remove Prettytable

* context manager to change variables (e.g. autoparal)

* Cleanup and refactoring in OpticTask

* Replace SIGRES with new fileformat based on SIGEPH (long-term)

* [eV] --> (ev) in title
