TODO list:

## High priority

* Add new section to manager.yml that allows users to customize limits according
  to task.__class__.__name__.
  Possible Yaml syntax for `task_class_limits`:

      limits:
         min_cores: 1
         max_cores: 1000
         timelimit: 2:0:0

         task_class_limits:
            # TaskClassName --> dict with new limits
            # Accept absolute values or `scale_name` syntax to scale `name` (mutually exclusive)
            # If a new limit is not specified, the global value is used.
            #
            NscfTask: {scale_max_cores: 0.5, scale_timelimit: 0.2}
            KerangeTask: {max_cores: 2, timelimit: 0:5:0}

  You should be able to specify the Qadapter as well.

* Recheck autoparal code in abinit. For the DFPT part, the computation of the weights is clearly 
  wrong because with nkpt 128 Abinit suggests to use nprocs 130 !!!!!!

* Implement Task modifier i.e. operations that change the input file if some condition occurs.
  This extra logic is require to handle problematic cases in which for instance the ScfTask does not converge
  and modification in the input file are required
  For instance, one may need to increase nline and/or diemac before restarting.

* Use angdeg instead of rprimd in structure_to_abivars if hex or rhomboedral lattice
  (tricky because input settings should be preserved)

* introduce new status for tasks that are removed at runtime e.g. S_CANCELLED
  and handle new case in flow machinery. Be careful with pickle, status comparison and ordering though.

* introduce new status WAITING_FOR_RESTART
  so that we don't have to restart task in callbacks

* Check Positive gw_qprange in EPH (Fixed by Henrique)

* Add iscf to GSR.nc so that we know if we have SCF|NSCF run.

* Improve exception handling in NetcdfReader

* Read forces in read_structure ? Fix problem with  MSONable and ArrayWithUnit/complex numbers

* Automate CHANGELOG creation.

* Refactor S_QCRITICAL logic (logic injected by user, since qcritical errors are cluster-specific)

* Refactor wrappers for mrgddb and mrgdvdb (problems with subprocess when
  merging large number of partial files (likely due to Popen with large stderr/stdout)

* BECS: 3x3 Tensor is not symmetric. Remove get_voigt_dataframe

* Parse stderr to detect runtime errors such as

    forrtl: severe (24): end-of-file during read, unit 5, file /proc/59090/fd/0
    Image              PC                Routine            Line        Source
    abinit             0000000008914AC2  for__io_return        Unknown  Unknown
    abinit             000000000894378D  for_read_seq_fmt      Unknown  Unknown
    abinit             000000000194409E  Unknown               Unknown  Unknown
    abinit             000000000042C671  Unknown               Unknown  Unknown
    abinit             000000000042C30E  Unknown               Unknown  Unknown
    libc-2.17.so       00002AAAB5DD8505  __libc_start_main     Unknown  Unknown
    abinit             000000000042C229  Unknown               Unknown  Unknown


  and kill the scheduler else the code gets stuck here (issue reported on lemaitre3)

* Remove/check the usage of line_density. Use ndivsm < 0 to activate line_density a la pymatgen.

* Check whether it's possible to reduce the memory of the SIGRES file wihouth brealing backward compatibility.
  eigvec_qp = self.read_variable("eigvec_qp") scales quite badly with nands

* Migrate to pyproject.toml. Implement script to generate requirements.yml
  Very likely `conda install --file requirements-optional.txt` is now broken to the introduction of `-r`.


## Medium priority

* Add support for PSML/UPF format

* Add support for new Abinit9 interface (getden_path, getwfk_path, pp_dirpath and pseudos)
  but remember that strings in the input should not be too long.
  Use common root for pseudos, what about getwfk_path? Need to refactor treatment of string lengths in Abinit!

* Interface abitk with AbiPy to compute DOS with tetra.

* videos in README (atom and hydrogen) or screenshot based on jupyterlab

* Refactor/improve Visualizer. See also jsmol, nglview and crystaltoolkit

* add possibility of changing amu in anaddb/abinit and API to "mix" DDB files
  phonon group velocities (requires extension in netcdf files).

* Scheduler should report info on exceptions (especially if at the end when on_all_ok is invoked)

* Add nsppol, nspinor, nspden to HIST file (and other stuff?)

* Fix bug with SCGW and SKW interpolation reported by Ahn. Sort energies

* Optimize SKW (slow if dense IBZ). Add possibility of initializing SKW
  from nc file produced by Fortran version.

* Add integration test for dilatmx error handler

* Add ExpectedAbort to Abinit so that one can call the code to get data without triggering
  report.errors in task.get_event.report

* Add memory error to Abinit errors

* Investigate NaN issue in BECS reported by Ahn if tolvrs instead of tolwfr (tolwfr could activate nbdbuf)

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
  as well as Abinit input as string). Input file has been added in Abini9 (input_string)

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

* fatbands with SOC (waiting for Matthieu's refactoring)

* Improvements in SKW. Finalize baseclass for ElectronInterpolator
  Average degenerate states.

* context manager to change variables (e.g. autoparal)

* Cleanup and refactoring in OpticTask (Well, optic should be rewritten from scratch)

* Replace SIGRES with new fileformat based on SIGEPH (long-term project)

* Update spack recipe and EasyBuild

* Classification of phonons/electrons

* Error handler for tolwfr to increase nband / nbdduf and resubmit
