"""Wrappers for ABINIT main executables"""

from __future__ import annotations

import os
from io import StringIO

import numpy as np
from monty.string import list_strings

from abipy.core.globals import get_workdir
from abipy.tools.numtools import data_from_cplx_mode
from abipy.tools.plotting import add_fig_kwargs
from abipy.tools.typing import Figure

__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"
__email__ = "gmatteo at gmail.com"
__status__ = "Development"
__date__ = "$Feb 21, 2013M$"

__all__ = [
    "Abitk",
    "Cut3D",
    "Fold2Bloch",
    "Lruj",
    "Mrgddb",
    "Mrgdvdb",
    "Mrgscr",
]


class ExecError(Exception):
    """Error class raised by :class:`ExecWrapper`"""


class ExecWrapper:
    """
    Base class that runs an executable in a subprocess.
    """

    Error = ExecError

    def __init__(self, manager=None, executable=None, verbose=0):
        """
        Args:
            manager: :class:`TaskManager` object responsible for the submission of the jobs.
                if manager is None, the default manager is used.
            executable: path to the executable.
            verbose: Verbosity level.
        """
        from .tasks import TaskManager

        self.manager = manager if manager is not None else TaskManager.from_user_config()
        self.manager = self.manager.to_shell_manager(mpi_procs=1)

        self.executable = executable if executable is not None else self.name
        assert os.path.basename(self.executable) == self.name
        self.verbose = int(verbose)

    def __str__(self) -> str:
        return "%s" % self.executable

    @property
    def name(self) -> str:
        return self._name

    def execute(self, workdir, exec_args=None) -> int:
        # Try to execute binary without and with mpirun.
        try:
            return self._execute(workdir, with_mpirun=True, exec_args=exec_args)
        except self.Error:
            return self._execute(workdir, with_mpirun=False, exec_args=exec_args)

    def _execute(self, workdir, with_mpirun=False, exec_args=None) -> int:
        """
        Execute the executable in a subprocess inside workdir.

        Some executables fail if we try to launch them with mpirun.
        Use with_mpirun=False to run the binary without it.
        """
        qadapter = self.manager.qadapter
        if not with_mpirun:
            qadapter.name = None
        if self.verbose:
            print("Working in:", workdir)

        script = qadapter.get_script_str(
            job_name=self.name,
            launch_dir=workdir,
            executable=self.executable,
            qout_path="qout_file.path",
            qerr_path="qerr_file.path",
            stdin=self.stdin_fname,
            stdout=self.stdout_fname,
            stderr=self.stderr_fname,
            exec_args=exec_args,
        )

        # Write the script.
        script_file = os.path.join(workdir, "run_" + self.name + ".sh")
        with open(script_file, "w") as fh:
            fh.write(script)
            os.chmod(script_file, 0o740)

        qjob, process = qadapter.submit_to_queue(script_file)
        self.stdout_data, self.stderr_data = process.communicate()
        self.returncode = process.returncode

        return self.returncode


class Mrgscr(ExecWrapper):
    """
    Wraps the mrgddb Fortran executable.
    """

    _name = "mrgscr"

    def merge_qpoints(self, workdir: str, files_to_merge: list[str], out_prefix: str) -> None:
        """
        Execute mrgscr inside directory `workdir` to merge `files_to_merge` over q-points.
        Produce new file in workdir with prefix `out_prefix`.
        """
        # We work with absolute paths.
        files_to_merge = [os.path.abspath(s) for s in list_strings(files_to_merge)]
        nfiles = len(files_to_merge)

        if self.verbose:
            print("Will merge %d files with output_prefix %s" % (nfiles, out_prefix))
            for i, f in enumerate(files_to_merge):
                print(" [%d] %s" % (i, f))

        if nfiles == 1:
            raise self.Error("merge_qpoints does not support nfiles == 1")

        self.stdin_fname, self.stdout_fname, self.stderr_fname = map(
            os.path.join, 3 * [workdir], ["mrgscr.stdin", "mrgscr.stdout", "mrgscr.stderr"]
        )

        inp = StringIO()
        inp.write(str(nfiles) + "\n")  # Number of partial files to merge.
        inp.write(out_prefix + "\n")  # Prefix for the final output file (_SCR extension will be added)
        for filename in files_to_merge:
            inp.write(filename + "\n")  # List with the files to merge.
        inp.write("1\n")  # Option for merging q-points.

        self.stdin_data = [s for s in inp.getvalue()]

        with open(self.stdin_fname, "w") as fh:
            fh.writelines(self.stdin_data)
            # Force OS to write data to disk.
            fh.flush()
            os.fsync(fh.fileno())

        self.execute(workdir)

    def merge_omegas(self, workdir: str, files_to_merge: list[str], out_prefix: str) -> None:
        """
        Execute mrgscr inside directory `workdir` to merge `files_to_merge` over frequencies.
        Produce new file in workdir with prefix `out_prefix`.
        """
        # We work with absolute paths.
        files_to_merge = [os.path.abspath(s) for s in list_strings(files_to_merge)]
        nfiles = len(files_to_merge)

        if self.verbose:
            print("Will merge %d files with output_prefix %s" % (nfiles, out_prefix))
            for i, f in enumerate(files_to_merge):
                print(" [%d] %s" % (i, f))

        if nfiles == 1:
            raise self.Error("merge_omegas does not support nfiles == 1")

        self.stdin_fname, self.stdout_fname, self.stderr_fname = map(
            os.path.join, 3 * [workdir], ["mrgscr.stdin", "mrgscr.stdout", "mrgscr.stderr"]
        )

        inp = StringIO()
        inp.write(str(nfiles) + "\n")  # Number of partial SCR files to merge.
        inp.write(out_prefix + "\n")  # Prefix for the final output file (_SCR extension will be added)
        for filename in files_to_merge:
            inp.write(filename + "\n")  # List with the files to merge.
        inp.write("2\n")  # Option for merging frequencies.
        inp.write("0.0\n")  # To use all freqs found.

        self.stdin_data = [s for s in inp.getvalue()]

        with open(self.stdin_fname, "w") as fh:
            fh.writelines(self.stdin_data)
            # Force OS to write data to disk.
            fh.flush()
            os.fsync(fh.fileno())

        self.execute(workdir)


class Mrgddb(ExecWrapper):
    """
    Wraps the mrgddb Fortran executable.
    """

    _name = "mrgddb"

    def merge(self, workdir, ddb_files, out_ddb, description, delete_source_ddbs=True) -> str:
        """Merge DDB file, return the absolute path of the new database in workdir."""
        # We work with absolute paths.
        ddb_files = [os.path.abspath(s) for s in list_strings(ddb_files)]
        if not os.path.isabs(out_ddb):
            out_ddb = os.path.join(os.path.abspath(workdir), os.path.basename(out_ddb))

        if self.verbose:
            print("Will merge %d files into output DDB %s" % (len(ddb_files), out_ddb))
            for i, f in enumerate(ddb_files):
                print(" [%d] %s" % (i, f))

        # Handle the case of a single file since mrgddb uses 1 to denote GS files!
        if len(ddb_files) == 1:
            with open(ddb_files[0]) as in_fh, open(out_ddb, "w") as out:
                out.writelines(in_fh)
            return out_ddb

        self.stdin_fname, self.stdout_fname, self.stderr_fname = map(
            os.path.join, 3 * [os.path.abspath(workdir)], ["mrgddb.stdin", "mrgddb.stdout", "mrgddb.stderr"]
        )

        inp = StringIO()
        inp.write(out_ddb + "\n")  # Name of the output file.
        inp.write(str(description) + "\n")  # Description.
        inp.write(str(len(ddb_files)) + "\n")  # Number of input DDBs.

        # Names of the DDB files.
        for fname in ddb_files:
            inp.write(fname + "\n")

        self.stdin_data = [s for s in inp.getvalue()]

        with open(self.stdin_fname, "w") as fh:
            fh.writelines(self.stdin_data)
            # Force OS to write data to disk.
            fh.flush()
            os.fsync(fh.fileno())

        retcode = self.execute(workdir, exec_args=["--nostrict"])
        if retcode == 0 and delete_source_ddbs:
            # Remove ddb files.
            for f in ddb_files:
                try:
                    os.remove(f)
                except OSError:
                    pass

        return out_ddb


class Mrgdvdb(ExecWrapper):
    """
    Wraps the mrgdvdb Fortran executable.
    """

    _name = "mrgdv"

    def merge(self, workdir, pot_files, out_dvdb, delete_source=True) -> str:
        """
        Merge POT files containing 1st order DFPT potential
        return the absolute path of the new database in workdir.

        Args:
            delete_source: True if POT1 files should be removed after (successful) merge.
        """
        # We work with absolute paths.
        pot_files = [os.path.abspath(s) for s in list_strings(pot_files)]
        if not os.path.isabs(out_dvdb):
            out_dvdb = os.path.join(os.path.abspath(workdir), os.path.basename(out_dvdb))

        if self.verbose:
            print("Will merge %d files into output DVDB %s" % (len(pot_files), out_dvdb))
            for i, f in enumerate(pot_files):
                print(" [%d] %s" % (i, f))

        # Handle the case of a single file since mrgddb uses 1 to denote GS files!
        if len(pot_files) == 1:
            with open(pot_files[0]) as in_fh, open(out_dvdb, "w") as out:
                out.writelines(in_fh)
            return out_dvdb

        self.stdin_fname, self.stdout_fname, self.stderr_fname = map(
            os.path.join, 3 * [os.path.abspath(workdir)], ["mrgdvdb.stdin", "mrgdvdb.stdout", "mrgdvdb.stderr"]
        )

        inp = StringIO()
        inp.write(out_dvdb + "\n")  # Name of the output file.
        inp.write(str(len(pot_files)) + "\n")  # Number of input POT files.

        # Names of the POT files.
        for fname in pot_files:
            inp.write(fname + "\n")

        self.stdin_data = [s for s in inp.getvalue()]

        with open(self.stdin_fname, "w") as fh:
            fh.writelines(self.stdin_data)
            # Force OS to write data to disk.
            fh.flush()
            os.fsync(fh.fileno())

        retcode = self.execute(workdir)
        if retcode == 0 and delete_source:
            # Remove pot files.
            for f in pot_files:
                try:
                    os.remove(f)
                except OSError:
                    pass

        return out_dvdb

    def test_ftinterp(self, dvdb_path, ngqpt, workdir=None, coarse_ngqpt=None, ddb_path="",
                       dvdb_add_lr=1, rspace_cell=0, symv1scf=0, qdamp=0.1, potfile=None) -> str:
        """
        Execute ``mrgdv test_ftinterp`` to test the Fourier interpolation of the DFPT
        potentials stored in the DVDB file ``dvdb_path``, dumping the ab-initio and
        interpolated V1(r) to a netcdf file for further analysis (see also
        :meth:`plot_ftinterp_parity`).

        Args:
            dvdb_path: Path to the input DVDB file.
            ngqpt: [nx, ny, nz] divisions of the ab-initio q-mesh used to build the DVDB.
            workdir: Working directory. If None, a temporary directory is created.
            coarse_ngqpt: Optional [nx, ny, nz] coarser sub-mesh, commensurate with ``ngqpt``
                (its IBZ q-points must already be present in the DVDB). If given, also
                Fourier-interpolates from this coarser mesh and compares against the literal
                (ab-initio) values at every q-point of the dense ``ngqpt`` mesh -- this is the
                actual interpolation-accuracy test. If None, only the native-mesh
                self-consistency round-trip is performed.
            ddb_path: Optional path to a DDB file with Born effective charges/dielectric
                tensor, needed to activate the long-range term in the interpolation.
            dvdb_add_lr: 1 to include the long-range term (requires ``ddb_path``), 0 to disable it.
            rspace_cell: 0 for the default real-space box, 1 for the Wigner-Seitz cell construction.
            symv1scf: Symmetrize the interpolated potentials (0, 1 or 2).
            qdamp: Gaussian damping parameter for the long-range term.
            potfile: Path to the output netcdf file with ab-initio/interpolated V1(r).
                If None, a default name is used inside ``workdir``.

        Return: absolute path to the netcdf file with the ab-initio/interpolated potentials.
        """
        workdir = get_workdir(workdir)
        os.makedirs(workdir, exist_ok=True)

        dvdb_path = os.path.abspath(dvdb_path)
        if ddb_path:
            ddb_path = os.path.abspath(ddb_path)

        if potfile is None:
            potfile = os.path.basename(dvdb_path) + "_FTINTERP.nc"
        if not os.path.isabs(potfile):
            potfile = os.path.join(os.path.abspath(workdir), potfile)

        # --coarse-ngqpt must always be passed on the command line (even as 0 0 0 to mean
        # "disabled") due to how mrgdv's CLI parser handles `want_len` for this option.
        coarse_ngqpt = [0, 0, 0] if coarse_ngqpt is None else list(coarse_ngqpt)

        exec_args = [
            "test_ftinterp", dvdb_path,
            "--ngqpt", *(str(n) for n in ngqpt),
            "--coarse-ngqpt", *(str(n) for n in coarse_ngqpt),
            "--dvdb-add-lr", str(dvdb_add_lr),
            "--rspace_cell", str(rspace_cell),
            "--symv1scf", str(symv1scf),
            "--qdamp", str(qdamp),
            "--potfile", potfile,
        ]
        if ddb_path:
            exec_args += ["--ddb-path", ddb_path]

        self.stdin_fname = None
        self.stdout_fname, self.stderr_fname = map(
            os.path.join, 2 * [workdir], ["test_ftinterp.stdout", "test_ftinterp.stderr"])

        if retcode := self.execute(workdir, exec_args=exec_args):
            print("stdout:\n", self.stdout_data)
            print("stderr:\n", self.stderr_data)
            raise RuntimeError(f"Error while running mrgdv test_ftinterp in {workdir}")

        return potfile

    @add_fig_kwargs
    def plot_ftinterp_parity(self, dvdb_path, ngqpt, workdir=None, coarse_ngqpt=None, ddb_path="",
                              dvdb_add_lr=1, rspace_cell=0, symv1scf=0, qdamp=0.1, potfile=None,
                              group=None, min_mag_frac=1e-6, bins=200, **kwargs) -> Figure:
        """
        Run :meth:`test_ftinterp` and produce a two-panel parity plot comparing ab-initio and
        Fourier-interpolated V1(r): |V1(r)| (log-log) and arg(V1(r)) (radians), both ab-initio
        (x-axis) vs interpolated (y-axis), as density-colored hexbins with a y=x reference line.

        Args:
            dvdb_path, ngqpt, workdir, coarse_ngqpt, ddb_path, dvdb_add_lr, rspace_cell,
                symv1scf, qdamp, potfile: See :meth:`test_ftinterp`.
            group: "self" (native-mesh round-trip) or "coarse" (the real accuracy test,
                requires ``coarse_ngqpt``). If None, defaults to "coarse" if ``coarse_ngqpt``
                is given, else "self".
            min_mag_frac: Both panels exclude points with |V1_abinitio| below this fraction of
                max|V1_abinitio|: such points are dominated by FFT/roundoff noise rather than
                real signal, and phase is meaningless there.
            bins: hexbin grid resolution.

        Return: |matplotlib-Figure|
        """
        potfile = self.test_ftinterp(dvdb_path, ngqpt, workdir=workdir, coarse_ngqpt=coarse_ngqpt,
                                      ddb_path=ddb_path, dvdb_add_lr=dvdb_add_lr, rspace_cell=rspace_cell,
                                      symv1scf=symv1scf, qdamp=qdamp, potfile=potfile)

        if group is None:
            group = "coarse" if coarse_ngqpt is not None else "self"

        import matplotlib.pyplot as plt

        from abipy.iotools import ETSF_Reader

        with ETSF_Reader(potfile) as r:
            # (nqpt, natom3, nspden, nfft) -- keep the block structure for the vdiff stats below.
            ai_c = r.read_value(f"{group}_v1r_abinitio", cmode="c")
            it_c = r.read_value(f"{group}_v1r_interp", cmode="c")

        # Global: pool every point into one L1-relative-error norm over the whole dataset
        # (sum|f1-f2| / sum|f2|) -- complements Pearson r, which stays close to 1 even under
        # large relative errors because it is dominated by the largest-magnitude points.
        global_stats = _vdiff_stats(ai_c.ravel(), it_c.ravel())

        # Worst single (iqpt, iatom3, ispden) perturbation: reproduces exactly what `mrgdv
        # test_ftinterp` itself prints as "Max values over q-points and perturbations".
        nfft = ai_c.shape[-1]
        block_stats = _vdiff_stats(ai_c.reshape(-1, nfft), it_c.reshape(-1, nfft), axis=-1)
        worst_stats = {k: v.max() for k, v in block_stats.items()}

        ai_c, it_c = ai_c.ravel(), it_c.ravel()
        mag_ai = data_from_cplx_mode("abs", ai_c)
        mag_it = data_from_cplx_mode("abs", it_c)
        ang_ai = data_from_cplx_mode("angle", ai_c)
        ang_it = data_from_cplx_mode("angle", it_c)

        # Drop the FFT/roundoff noise floor -- see docstring of `min_mag_frac`.
        keep = mag_ai >= min_mag_frac * mag_ai.max()

        fig, axes = plt.subplots(1, 2, figsize=(11, 5.6))

        _hexbin_parity(axes[0], mag_ai[keep], mag_it[keep], bins, log=True,
                        label=f"|V1(r)|  ({group}, |V1|>{min_mag_frac:g}*max)")
        axes[0].set_xlabel("ab-initio  |V1(r)|  (Ha)")
        axes[0].set_ylabel("interpolated  |V1(r)|  (Ha)")

        _hexbin_parity(axes[1], ang_ai[keep], ang_it[keep], bins, log=False, wrap=True,
                        label=f"arg(V1(r))  ({group}, |V1|>{min_mag_frac:g}*max)")
        axes[1].set_xlabel("ab-initio  arg(V1(r))  (rad)")
        axes[1].set_ylabel("interpolated  arg(V1(r))  (rad, branch-aligned)")

        def fmt_row(row_label, s):
            return (f"{row_label}: L1_rerr={100*s['l1_rerr']:.2f}%  mean|diff|={s['mean_adiff']:.3e}  "
                    f"max|diff|={s['max_adiff']:.3e}  stdev|diff|={s['stdev_adiff']:.3e}")

        fig.text(0.5, 0.02,
                  "V1(r) complex-valued error (ABINIT's vdiff_t convention, L1_rerr = sum|f1-f2| / sum|f2|):\n"
                  + fmt_row("global, all points pooled       ", global_stats) + "\n"
                  + fmt_row("worst single perturbation (=mrgdv)", worst_stats),
                  ha="center", va="bottom", fontsize=8, family="monospace")

        fig.suptitle(os.path.basename(potfile))
        fig.tight_layout(rect=[0, 0.14, 1, 0.96])

        return fig

    def test_symcheck(self, dvdb_path, ngqpt, qpt, workdir=None, ddb_path="", dvdb_add_lr=1,
                       rspace_cell=0, symv1scf=0, qdamp=0.1, potfile=None) -> str:
        """
        Execute ``mrgdv test_symcheck`` to test CROSS-Q-POINT symmetry consistency of the
        Fourier interpolation: interpolate at ``qpt`` and, independently, at ``S.qpt`` for
        every symmetry ``S`` of the crystal (both directions of time reversal), and compare
        against the prediction obtained by rotating the ``qpt`` interpolation with
        ``v1phq_rotate`` -- the same formula used throughout ABINIT to expand an IBZ q-point
        to the full BZ. Unlike :meth:`test_ftinterp`, this never touches literal/ab-initio
        data: both sides being compared are themselves Fourier-interpolated, so it isolates
        whether the interpolation is internally consistent with the crystal's own symmetry,
        as opposed to :meth:`test_ftinterp`'s literal-ground-truth accuracy question. Dumps
        the two compared quantities to a netcdf file for further analysis (see also
        :meth:`plot_symcheck_parity`).

        Args:
            dvdb_path: Path to the input DVDB file.
            ngqpt: [nx, ny, nz] divisions of the ab-initio q-mesh used to build the DVDB.
            qpt: [qx, qy, qz] source q-point (reduced coordinates). Need not be on the
                ab-initio mesh -- interpolating it is itself a genuine off-grid test.
            workdir: Working directory. If None, a temporary directory is created.
            ddb_path, dvdb_add_lr, rspace_cell, symv1scf, qdamp: See :meth:`test_ftinterp`.
            potfile: Path to the output netcdf file. If None, a default name is used inside
                ``workdir``.

        Return: absolute path to the netcdf file with the target/predicted potentials.
        """
        workdir = get_workdir(workdir)
        os.makedirs(workdir, exist_ok=True)

        dvdb_path = os.path.abspath(dvdb_path)
        if ddb_path:
            ddb_path = os.path.abspath(ddb_path)

        if potfile is None:
            potfile = os.path.basename(dvdb_path) + "_SYMCHECK.nc"
        if not os.path.isabs(potfile):
            potfile = os.path.join(os.path.abspath(workdir), potfile)

        exec_args = [
            "test_symcheck", dvdb_path,
            "--ngqpt", *(str(n) for n in ngqpt),
            "--qpt", *(str(q) for q in qpt),
            "--dvdb-add-lr", str(dvdb_add_lr),
            "--rspace_cell", str(rspace_cell),
            "--symv1scf", str(symv1scf),
            "--qdamp", str(qdamp),
            "--potfile", potfile,
        ]
        if ddb_path:
            exec_args += ["--ddb-path", ddb_path]

        self.stdin_fname = None
        self.stdout_fname, self.stderr_fname = map(
            os.path.join, 2 * [workdir], ["test_symcheck.stdout", "test_symcheck.stderr"])

        if retcode := self.execute(workdir, exec_args=exec_args):
            print("stdout:\n", self.stdout_data)
            print("stderr:\n", self.stderr_data)
            raise RuntimeError(f"Error while running mrgdv test_symcheck in {workdir}")

        return potfile

    @add_fig_kwargs
    def plot_symcheck_parity(self, dvdb_path, ngqpt, qpt, workdir=None, ddb_path="", dvdb_add_lr=1,
                              rspace_cell=0, symv1scf=0, qdamp=0.1, potfile=None,
                              min_mag_frac=1e-6, bins=200, **kwargs) -> Figure:
        """
        Run :meth:`test_symcheck` and produce a two-panel parity plot: |V1(r)| (log-log) and
        arg(V1(r)) (radians), independently-interpolated target (x-axis) vs the
        `v1phq_rotate`-predicted value from the source q-point (y-axis), pooled over every
        symmetry operation and both signs of time reversal, as density-colored hexbins with
        a y=x reference line.

        Args:
            dvdb_path, ngqpt, qpt, workdir, ddb_path, dvdb_add_lr, rspace_cell, symv1scf,
                qdamp, potfile: See :meth:`test_symcheck`.
            min_mag_frac: Both panels exclude points with |V1_target| below this fraction of
                max|V1_target|: such points are dominated by FFT/roundoff noise rather than
                real signal, and phase is meaningless there.
            bins: hexbin grid resolution.

        Return: |matplotlib-Figure|
        """
        potfile = self.test_symcheck(dvdb_path, ngqpt, qpt, workdir=workdir, ddb_path=ddb_path,
                                      dvdb_add_lr=dvdb_add_lr, rspace_cell=rspace_cell,
                                      symv1scf=symv1scf, qdamp=qdamp, potfile=potfile)

        import matplotlib.pyplot as plt

        from abipy.iotools import ETSF_Reader

        with ETSF_Reader(potfile) as r:
            # (nsym, ntimrev, natom3, nspden, nfft) -- keep the block structure for vdiff stats.
            pred_c = r.read_value("v1r_predicted", cmode="c")
            tgt_c = r.read_value("v1r_target", cmode="c")

        # f1=predicted, f2=target, matching dvdb_test_symcheck's own `vd%eval` argument order
        # (Fortran's L1_rerr is normalized by f2) so these numbers match mrgdv's own stdout.
        global_stats = _vdiff_stats(pred_c.ravel(), tgt_c.ravel())
        nfft = pred_c.shape[-1]
        block_stats = _vdiff_stats(pred_c.reshape(-1, nfft), tgt_c.reshape(-1, nfft), axis=-1)
        worst_stats = {k: v.max() for k, v in block_stats.items()}

        pred_c, tgt_c = pred_c.ravel(), tgt_c.ravel()
        mag_pred, mag_tgt = data_from_cplx_mode("abs", pred_c), data_from_cplx_mode("abs", tgt_c)
        ang_pred, ang_tgt = data_from_cplx_mode("angle", pred_c), data_from_cplx_mode("angle", tgt_c)

        # Drop the FFT/roundoff noise floor -- see docstring of `min_mag_frac`.
        keep = mag_tgt >= min_mag_frac * mag_tgt.max()

        fig, axes = plt.subplots(1, 2, figsize=(11, 5.6))

        _hexbin_parity(axes[0], mag_tgt[keep], mag_pred[keep], bins, log=True,
                        label=f"|V1(r)|  (|V1|>{min_mag_frac:g}*max)")
        axes[0].set_xlabel("independent interp. (target)  |V1(r)|  (Ha)")
        axes[0].set_ylabel("v1phq_rotate prediction  |V1(r)|  (Ha)")

        _hexbin_parity(axes[1], ang_tgt[keep], ang_pred[keep], bins, log=False, wrap=True,
                        label=f"arg(V1(r))  (|V1|>{min_mag_frac:g}*max)")
        axes[1].set_xlabel("independent interp. (target)  arg(V1(r))  (rad)")
        axes[1].set_ylabel("v1phq_rotate prediction  arg(V1(r))  (rad, branch-aligned)")

        def fmt_row(row_label, s):
            return (f"{row_label}: L1_rerr={100*s['l1_rerr']:.2f}%  mean|diff|={s['mean_adiff']:.3e}  "
                    f"max|diff|={s['max_adiff']:.3e}  stdev|diff|={s['stdev_adiff']:.3e}")

        fig.text(0.5, 0.02,
                  "V1(r) complex-valued error (ABINIT's vdiff_t convention, L1_rerr = sum|predicted-target| / sum|target|):\n"
                  + fmt_row("global, all (isym, itimrev) pooled  ", global_stats) + "\n"
                  + fmt_row("worst single (isym, itimrev) (=mrgdv)", worst_stats),
                  ha="center", va="bottom", fontsize=8, family="monospace")

        fig.suptitle(os.path.basename(potfile))
        fig.tight_layout(rect=[0, 0.14, 1, 0.96])

        return fig


def _vdiff_stats(f1, f2, axis=-1) -> dict:
    """
    Reproduce ABINIT's own ``vdiff_t%eval`` (``m_numeric_tools.F90``) error estimators -- the
    same ones printed by ``mrgdv test_ftinterp`` itself -- from complex arrays ``f1``
    (ab-initio) and ``f2`` (interpolated), reducing over ``axis``. ``l1_rerr`` is normalized
    by ``|f2|`` (the INTERPOLATED array), matching the Fortran convention exactly (verified
    against ``mrgdv``'s own printed output: the block-wise max reproduces it to 5 sig figs).
    """
    diff = np.abs(f1 - f2)
    num = diff.sum(axis=axis)
    den = np.abs(f2).sum(axis=axis)
    l1_rerr = np.divide(num, den, out=np.zeros_like(num, dtype=float), where=den != 0)
    return dict(l1_rerr=l1_rerr, mean_adiff=diff.mean(axis=axis), max_adiff=diff.max(axis=axis),
                min_adiff=diff.min(axis=axis), stdev_adiff=diff.std(axis=axis))


def _hexbin_parity(ax, x, y, bins, log=False, wrap=False, label="") -> None:
    """
    Helper for :meth:`Mrgdvdb.plot_ftinterp_parity`: density-colored hexbin parity plot
    of ``y`` vs ``x`` with a y=x reference line, Pearson r / RMSE annotated.

    Args:
        wrap: True for angle data (radians), which is only defined modulo 2*pi: a
            physically tiny difference near the +-pi branch cut (e.g. x=+3.14, y=-3.14)
            would otherwise register as a difference of ~2*pi and dominate the Pearson r
            / RMSE statistics. Realigns y to the branch closest to x
            (y -> x + wrap_to_pi(y - x), which does not change y's physical value)
            before doing anything else, so both the plot and the stats reflect the true
            agreement.
    """
    import matplotlib.pyplot as plt

    if wrap:
        y = x + np.angle(np.exp(1j * (y - x)))

    if log:
        m = (x > 0) & (y > 0)
        x, y = x[m], y[m]

    lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
    kwargs = dict(xscale="log", yscale="log") if log else {}
    h = ax.hexbin(x, y, gridsize=bins, cmap="Blues", mincnt=1, bins="log", **kwargs)
    ax.plot([lo, hi], [lo, hi], "--", color="0.35", lw=1.2, zorder=3)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")

    cb = plt.colorbar(h, ax=ax, shrink=0.85)
    cb.set_label("point count (log)")

    ax.grid(True, which="both", ls=":", lw=0.5, color="0.85", zorder=0)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    r = np.corrcoef(x, y)[0, 1]
    rmse = np.sqrt(np.mean((x - y) ** 2))
    ax.text(0.03, 0.97, f"{label}\nN = {x.size:,}\nPearson r = {r:.5f}\nRMSE = {rmse:.3e}",
            transform=ax.transAxes, ha="left", va="top", fontsize=9,
            bbox=dict(boxstyle="round", fc="white", ec="0.8", alpha=0.85))


class Cut3D(ExecWrapper):
    """
    Wraps the cut3d Fortran executable.
    """

    _name = "cut3d"

    def cut3d(self, cut3d_input, workdir) -> tuple[str, str]:
        """
        Runs cut3d with a Cut3DInput

        Args:
            cut3d_input: a Cut3DInput object.
            workdir: directory where cut3d is executed.

        Returns:
            (string) absolute path to the standard output of the cut3d execution.
            (string) absolute path to the output filepath. None if output is required.
        """
        self.stdin_fname, self.stdout_fname, self.stderr_fname = map(
            os.path.join, 3 * [os.path.abspath(workdir)], ["cut3d.stdin", "cut3d.stdout", "cut3d.stderr"]
        )

        cut3d_input.write(self.stdin_fname)

        if retcode := self._execute(workdir, with_mpirun=False):
            stdout = os.path.join(workdir, "cut3d.stdout")
            stderr = os.path.join(workdir, "cut3d.stderr")
            if os.path.exists(stdout):
                with open(stdout) as fh:
                    print(fh.read())
            if os.path.exists(stderr):
                with open(stderr) as fh:
                    print(fh.read())

            raise RuntimeError("Error while running cut3d in %s" % workdir)

        output_filepath = cut3d_input.output_filepath

        if output_filepath is not None:
            if not os.path.isabs(output_filepath):
                output_filepath = os.path.abspath(os.path.join(workdir, output_filepath))

            if not os.path.isfile(output_filepath):
                raise RuntimeError("The file was not converted correctly in %s." % workdir)

        return self.stdout_fname, output_filepath


class Fold2Bloch(ExecWrapper):
    """
    Wraps the fold2Bloch Fortran executable.
    """

    _name = "fold2Bloch"

    def unfold(self, wfkpath, folds, workdir=None) -> str:
        """Unfold the wavefunctions from the supercell to the primitive cell."""
        workdir = get_workdir(workdir)

        self.stdin_fname = None
        self.stdout_fname, self.stderr_fname = map(
            os.path.join, 2 * [workdir], ["fold2bloch.stdout", "fold2bloch.stderr"]
        )

        folds = np.array(folds, dtype=int).flatten()
        if len(folds) not in (3, 9):
            raise ValueError("Expecting 3 ints or 3x3 matrix but got %s" % (str(folds)))
        fold_arg = ":".join(str(f) for f in folds)
        wfkpath = os.path.abspath(wfkpath)
        if not os.path.isfile(wfkpath):
            raise RuntimeError("WFK file `%s` does not exist in %s" % (wfkpath, workdir))

        # Usage: $ fold2Bloch file_WFK x:y:z (folds)
        if retcode := self.execute(workdir, exec_args=[wfkpath, fold_arg]):
            print("stdout:\n", self.stdout_data)
            print("stderr:\n", self.stderr_data)
            raise RuntimeError("fold2bloch returned %s in %s" % (retcode, workdir))

        filepaths = [f for f in os.listdir(workdir) if f.endswith("_FOLD2BLOCH.nc")]
        if len(filepaths) != 1:
            raise RuntimeError("Cannot find *_FOLD2BLOCH.nc file in: %s" % str(os.listdir(workdir)))

        return os.path.join(workdir, filepaths[0])


class Lruj(ExecWrapper):
    """
    Wraps the lruj Fortran executable.
    """

    _name = "lruj"

    def run(self, nc_paths: list[str], workdir=None) -> int:
        """
        Execute lruj inside directory `workdir` to analyze `nc_paths`.
        """
        workdir = get_workdir(workdir)

        self.stdin_fname = None
        self.stdout_fname, self.stderr_fname = map(os.path.join, 2 * [workdir], ["lruj.stdout", "lruj.stderr"])

        # We work with absolute paths.
        nc_paths = [os.path.abspath(s) for s in list_strings(nc_paths)]

        if retcode := self.execute(workdir, exec_args=nc_paths):
            print("stdout:\n", self.stdout_data)
            print("stderr:\n", self.stderr_data)
            raise RuntimeError(f"Error while running lruj in {workdir}")

        return retcode


class Abitk(ExecWrapper):
    """
    Wraps the abitk Fortran executable.
    """

    _name = "abitk"

    stdin_fname = None

    def run(self, exec_args: list, workdir=None) -> int:
        """
        Execute abitk inside directory `workdir`.
        """
        workdir = get_workdir(workdir)

        self.stdout_fname, self.stderr_fname = map(os.path.join, 2 * [workdir], ["abitk.stdout", "abitk.stderr"])

        if retcode := self.execute(workdir, exec_args=exec_args):
            print("stdout:\n", self.stdout_data)
            print("stderr:\n", self.stderr_data)
            raise RuntimeError(f"Error while running abitk in {workdir}")

        return retcode
