"""
This module contains objects for postprocessing polaron calculations
using the results stored in the VARPEQ.nc file.

For a theoretical introduction see ...
"""
from __future__ import annotations

import dataclasses
import numpy as np
#import pandas as pd
import abipy.core.abinit_units as abu

from collections import defaultdict
from monty.string import marquee #, list_strings
from monty.functools import lazy_property
from monty.termcolor import cprint
from abipy.core.func1d import Function1D
from abipy.core.structure import Structure
from abipy.core.kpoints import kpoints_indices
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_ElectronBands, Has_Header, NotebookWriter
from abipy.tools.typing import PathLike
from abipy.tools.plotting import (add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, set_axlims, set_visible,
    rotate_ticklabels, ax_append_title, set_ax_xylabels, linestyles, Marker)
#from abipy.tools import duck
from abipy.electrons.ebands import ElectronBands, RobotWithEbands
from abipy.dfpt.phonons import PhononBands
from abipy.tools.typing import Figure
from abipy.tools.numtools import BzRegularGridInterpolator, gaussian
from abipy.abio.robots import Robot
from abipy.eph.common import BaseEphReader

ITER_LABELS = [
    r'$E_{pol}$',
    r'$E_{el}$',
    r'$E_{ph}$',
    r'$E_{elph}$',
    r'$\varepsilon$'
]


class VarpeqFile(AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands, NotebookWriter):
    """
    This file stores the results of a VARPEQ calculations: SCF cycle, A_nk, B_qnu
    and provides methods to analyze and plot results.

    Usage example:

    .. code-block:: python

        with VarpeqFile("out_VARPEQ.nc") as varpeq:
            print(varpeq)
            varpeq.plot_scf_cycle()

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: VarpeqFile
    """

    @classmethod
    def from_file(cls, filepath: PathLike) -> VarpeqFile:
        """Initialize the object from a netcdf file."""
        return cls(filepath)

    def __init__(self, filepath: PathLike):
        super().__init__(filepath)
        self.r = VarpeqReader(filepath)

    @lazy_property
    def ebands(self) -> ElectronBands:
        """|ElectronBands| object."""
        return self.r.read_ebands()

    @property
    def structure(self) -> Structure:
        """|Structure| object."""
        return self.ebands.structure

    def close(self) -> None:
        """Close the file."""
        self.r.close()

    @lazy_property
    def polaron_spin(self) -> list[Polaron]:
        """List of polaron objects, one for each spin (if any)."""
        return [Polaron.from_varpeq(self, spin) for spin in range(self.r.nsppol)]

    @lazy_property
    def params(self) -> dict:
        """dict with the convergence parameters, e.g. ``nbsum``."""
        #od = OrderedDict([
        #    ("nbsum", self.nbsum),
        #    ("nqbz", self.r.nqbz),
        #    ("nqibz", self.r.nqibz),
        #])
        ## Add EPH parameters.
        #od.update(self.r.common_eph_params)

        od = {}
        return od

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose=0) -> str:
        """String representation with verbosiy level ``verbose``."""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))

        app("")
        app(self.ebands.to_string(with_structure=False, verbose=verbose, title="Electronic Bands"))
        #if verbose > 1:
        #    app("")
        #    app(self.hdr.to_string(verbose=verbose, title="Abinit Header"))

        app("")
        app("VARPEQ parameters:")
        app(f"varpeq_pkind: {self.r.varpeq_pkind}")
        #app(f"gstore_completed: {bool(self.r.completed)}")
        #app(f"gstore_cplex: {self.r.cplex}")
        #app(f"gstore_kzone: {self.r.kzone}")
        #app(f"gstore_kfilter: {self.r.kfilter}")
        #app(f"gstore_gmode: {self.r.gmode}")
        #app(f"gstore_qzone: {self.r.qzone}")
        #app(f"gstore_with_vk: {self.r.with_vk}")
        #app(f"gstore_kptopt: {self.r.kptopt}")
        #app(f"gstore_qptopt: {self.r.qptopt}")

        return "\n".join(lines)

    def get_last_iteration_dict_ev(self, spin: int) -> dict:
        """
        Return dictionary mapping the latex label to the value of the last iteration
        for the given spin index. All energies are in eV.
        """
        nstep2cv_spin = self.r.read_value('nstep2cv')
        iter_rec_spin = self.r.read_value('iter_rec')
        nstep2cv = nstep2cv_spin[spin]
        last_iteration = iter_rec_spin[spin, nstep2cv-1, :] * abu.Ha_eV

        return dict(zip(ITER_LABELS, last_iteration))

    @add_fig_kwargs
    def plot_scf_cycle(self, ax_mat=None, fontsize=12, **kwargs) -> Figure:
        """
        Plot the VARPEQ SCF cycle.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: fontsize for legends and titles
        """
        nsppol = self.r.nsppol
        nstep2cv_spin = self.r.read_value('nstep2cv')
        iter_rec_spin = self.r.read_value('iter_rec')

        # Build grid of plots.
        nrows, ncols = nsppol, 2
        ax_mat, fig, plt = get_axarray_fig_plt(ax_mat, nrows=nrows, ncols=ncols,
                                               sharex=False, sharey=False, squeeze=False)

        for spin in range(nsppol):
            nstep2cv = nstep2cv_spin[spin]
            iterations = iter_rec_spin[spin, :nstep2cv, :] * abu.Ha_eV
            xs = np.arange(1, nstep2cv + 1)

            for iax, ax in enumerate(ax_mat[spin]):
                for ilab, label in enumerate(ITER_LABELS):
                    ys = iterations[:,ilab]
                    if iax == 0:
                        # Plot energies in linear scale.
                        ax.plot(xs, ys, label=label)
                    else:
                        # Plot deltas in logscale.
                        ax.plot(xs, np.abs(ys - ys[-1]), label=label)
                        ax.set_yscale("log")

                ax.set_xlim(1, nstep2cv)
                ax.grid(True)
                ax.legend(loc="best", shadow=True, fontsize=fontsize)
                ax.set_xlabel("Iteration", fontsize=fontsize)
                ax.set_ylabel("Energy (eV)" if iax == 0 else r"$|\Delta|$ Energy (eV)", fontsize=fontsize)
                if nsppol == 2: ax.set_title(f"{spin=}", fontsize=fontsize)

        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        yield self.plot_scf_cycle(show=False)

    def write_notebook(self, nbpath=None) -> str:
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("varpeq = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(varpeq)"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


@dataclasses.dataclass(kw_only=True)
class Polaron:
    """
    This object stores the polaron coefficients A_kn, B_qnu for a given spin.
    Provides methods to plot |A_nk|^2 or |B_qnu|^2 together with band structures.
    """
    spin: int          # Spin index.
    nb: int            # Number of bands.
    nk: int            # Number of k-points (including filtering if any)
    nq: int            # Number of q-points (including filtering if any)
    bstart: int        # First band starts at bstart
    bstop: int

    kpoints: np.ndarray
    qpoints: np.ndarray
    a_kn: np.ndarray
    b_qnu: np.ndarray
    varpeq: VarpeqFile

    @classmethod
    def from_varpeq(cls, varpeq: VarpeqFile, spin: int) -> Polaron:
        """
        Build an istance from a VarpeqFile and the spin index.
        """
        r = varpeq.r
        nk, nq, nb = r.nk_spin[spin], r.nq_spin[spin], r.nb_spin[spin]
        bstart, bstop = r.brange_spin[spin]
        kpoints = r.read_value("kpts_spin")[spin, :nk]
        qpoints = r.read_value("qpts_spin")[spin, :nq]
        a_kn = r.read_value("a_spin", cmode="c")[spin, :nk, :nb]
        b_qnu = r.read_value("b_spin", cmode="c")[spin, :nq]

        data = locals()
        return cls(**{k: data[k] for k in [field.name for field in dataclasses.fields(Polaron)]})

    @property
    def structure(self):
        return self.varpeq.structure

    @property
    def ebands(self):
        return self.varpeq.ebands

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose=0) -> str:
        """String representation with verbosiy level ``verbose``."""
        lines = []; app = lines.append

        app(marquee(f"Ank for spin: {self.spin}", mark="="))
        app(f"nb: {self.nb}")
        app(f"nk: {self.nk}")
        app(f"nq: {self.nq}")
        app(f"bstart: {self.bstart}")
        app(f"bstop: {self.bstop}")
        ksampling = self.ebands.kpoints.ksampling
        ngkpt, shifts = ksampling.mpdivs, ksampling.shifts
        app(f"ksampling: {str(ksampling)}")
        ngqpt = self.varpeq.r.ngqpt
        app(f"q-mesh: {ngqpt}")

        #if verbose:
        norm = np.sum(np.abs(self.a_kn) ** 2) / self.nk
        app("1/N_k sum_{nk} |A_nk|^2: %f" % norm)
        #norm = np.sum(np.abs(self.b_qnu) ** 2)
        #app("sum_{qnu} |B_qnu|^2: %f" % norm)

        return "\n".join(lines)

    def get_a2_interpolator(self, method: str, check_mesh: int = 0) -> BzRegularGridInterpolator:
        """
        Build and return an interpolator for |A_nk|^2

        Args:
            method: String defining the interpolation method.
            check_mesh: Check whether k-points belong to the mesh ff !=0.
        """
        # Neeed to know the size of the k-mesh.
        ksampling = self.ebands.kpoints.ksampling
        ngkpt, shifts = ksampling.mpdivs, ksampling.shifts

        if ngkpt is None:
            raise ValueError("Non diagonal k-meshes are not supported")
        if len(shifts) > 1:
            raise ValueError("Multiple k-shifts are not supported!")

        k_indices = kpoints_indices(self.kpoints, ngkpt, check_mesh=check_mesh)

        nx, ny, nz = ngkpt
        a_data = np.empty((self.nb, nx, ny, nz), dtype=complex)
        for ib in range(self.nb):
            for a_cplx, k_inds in zip(self.a_kn[:,ib], k_indices):
                ix, iy, iz = k_inds
                a_data[ib, ix, iy, iz] = a_cplx

        return BzRegularGridInterpolator(self.structure, shifts, np.abs(a_data) ** 2,
                                         method=method)

    def get_b2_interpolator(self, method: str, check_mesh: int = 0) -> BzRegularGridInterpolator:
        """
        Build and return an interpolator for |B_qnu|^2.

        Args
            method: String defining the interpolation method.
            check_mesh: Check whether k-points belong to the mesh ff !=0.
        """
        # Neeed to know the size of the q-mesh (always Gamma-centered)
        ngqpt, shifts = self.varpeq.r.ngqpt, [0, 0, 0]
        q_indices = kpoints_indices(self.qpoints, ngqpt, check_mesh=check_mesh)

        natom3 = 3 * len(self.structure)
        nx, ny, nz = ngqpt
        b_data = np.empty((natom3, nx, ny, nz), dtype=complex)
        for nu in range(natom3):
            for b_cplx, q_inds in zip(self.b_qnu[:,nu], q_indices):
                ix, iy, iz = q_inds
                b_data[nu, ix, iy, iz] = b_cplx

        return BzRegularGridInterpolator(self.structure, shifts, np.abs(b_data) ** 2,
                                         method=method)

    @add_fig_kwargs
    def plot_bz_sampling(self, what="kpoints", fold=False,
                        ax=None, pmg_path=True, with_labels=True, **kwargs) -> Figure:
        """
        Plots a 3D representation of the Brillouin zone with the sampling.

        Args:
            what: "kpoints" or "qpoints"
            fold: whether the points should be folded inside the first Brillouin Zone.
                Defaults to False.
        """
        bz_points = dict(kpoints=self.kpoints, qpoints=self.qpoints)[what]
        kws = dict(ax=ax, pmg_path=pmg_path, with_labels=with_labels, fold=fold,
                   kpoints=bz_points)

        return self.structure.plot_bz(show=False, **kws)

    @add_fig_kwargs
    def plot_ank_with_ebands(self, ebands_kpath, ebands_kmesh=None,
                             ax=None, scale=50, fontsize=12, **kwargs) -> Figure:
        """
        Plot electronic energies with markers whose size is proportional to |A_nk|^2.

        Args:
            ebands: ElectronBands or Abipy file providing an electronic band structure along a path.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            scale: Scaling factor for |A_nk|^2.
        """
        ebands_kpath = ElectronBands.as_ebands(ebands_kpath)

        # Interpolate A_nk
        a2_interp = self.get_a2_interpolator("linear")

        # DEBUG SECTION
        #ref_kn = np.abs(self.a_kn) ** 2
        #for ik, kpoint in enumerate(self.kpoints):
        #    interp = a2_interp.eval_kpoint(kpoint)
        #    print("MAX (A2 ref - A2 interp) at qpoint", kpoint)
        #    print((np.abs(ref_kn[ik] - interp)).max())

        ymin, ymax = +np.inf, -np.inf
        x, y, s = [], [], []
        for ik, kpoint in enumerate(ebands_kpath.kpoints):
            enes_n = ebands_kpath.eigens[self.spin, ik, self.bstart:self.bstop]
            a2_n = a2_interp.eval_kpoint(kpoint.frac_coords)
            for e, a2 in zip(enes_n, a2_n):
                x.append(ik); y.append(e); s.append(scale * a2)
                ymin = min(ymin, e)
                ymax = max(ymax, e)

        points = Marker(x, y, s)

        nrows, ncols = 1, 2
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=False, sharey=True, squeeze=False)
        ax_list = ax_list.ravel()

        ax = ax_list[0]
        ebands_kpath.plot(ax=ax, points=points, show=False)

        vertices_names = [(k.frac_coords, k.name) for k in ebands_kpath.kpoints]
        lpratio = 5
        kmesh = np.array([12, 12, 12]) * 3
        step = 0.1
        width = 0.2

        if ebands_kmesh is None:
            # Compute ebands_kmesh with Star-function interpolation.
            r = self.ebands.interpolate(lpratio=lpratio, vertices_names=vertices_names, kmesh=kmesh)
            ebands_kmesh = r.ebands_kmesh

        # Get electronic DOS.
        edos = ebands_kmesh.get_edos(step=step, width=width)

        mesh = edos.spin_dos[self.spin].mesh
        ank_dos = np.zeros(len(mesh))
        e0 = self.ebands.fermie

        ymin -= 0.5 * abs(ymin)
        ymin -= e0
        ymax += 0.5 * abs(ymax)
        ymax -= e0
        #ymax = None

        #################
        # Compute Ank DOS
        #################
        # NB: This is just to sketch the ideas. I don't think the present version
        # is correct as only k --> -k symmetry can be used.

        for ik, kpoint in enumerate(ebands_kmesh.kpoints):
            weight = kpoint.weight
            enes_n = ebands_kmesh.eigens[self.spin, ik, self.bstart:self.bstop]
            a2_n = a2_interp.eval_kpoint(kpoint)
            for e, a2 in zip(enes_n, a2_n):
                ank_dos += weight * a2 * gaussian(mesh, width, center=e-e0)

        ax = ax_list[1]
        edos.plot_ax(ax, e0, spin=self.spin, exchange_xy=True, label="eDOS(E)")

        ank_dos = Function1D(mesh, ank_dos)
        ank_dos.plot_ax(ax, exchange_xy=True, label=r"$A^2$(E)")
        ax.grid(True)
        ax.legend(loc="best", shadow=True, fontsize=fontsize)
        print("A2(E) integrates to", ank_dos.integral_value)

        for ax in ax_list:
            ax.set_ylim(ymin, ymax)

        return fig

    @add_fig_kwargs
    def plot_bqnu_with_phbands(self, phbands_qpath, ax=None, scale=10, **kwargs) -> Figure:
        """
        Plot phonon energies with markers whose size is proportional to |B_qnu|^2.

        Args:
            phbands_qpath: PhononBands or Abipy file providing a phonon band structure.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            scale: Scaling factor for |B_qnu|^2.
        """
        phbands_qpath = PhononBands.as_phbands(phbands_qpath)
        b2_interp = self.get_b2_interpolator("linear")

        # DEBUG SECTION
        #ref_qnu = np.abs(self.b_qnu) ** 2
        #for iq, qpoint in enumerate(self.qpoints):
        #    print("MAX (B2 ref - B2 interp) at qpoint", qpoint)
        #    interp = b2_interp.eval_kpoint(qpoint)
        #    print((np.abs(ref_qnu[iq] - interp)).max())

        x, y, s = [], [], []
        for iq, qpoint in enumerate(phbands_qpath.qpoints):
            omegas_nu = phbands_qpath.phfreqs[iq,:]
            b2_nu = b2_interp.eval_kpoint(qpoint.frac_coords)
            assert len(omegas_nu) == len(b2_nu)
            for w, b2 in zip(omegas_nu, b2_nu):
                x.append(iq); y.append(w); s.append(scale * b2)
        points = Marker(x, y, s)

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        phbands_qpath.plot(ax=ax, points=points, show=False)

        return fig


class VarpeqReader(BaseEphReader):
    """
    Reads data from file and constructs objects.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: VarpeqReader
    """

    def __init__(self, filepath: PathLike):
        super().__init__(filepath)

        #char input_string(input_len) ;
        #int eph_task ;
        #int varpeq_nstep ;
        #int nkbz ;
        #int nqbz ;
        #double tolgrs ;
        #int nstep2cv(nsppol) ;
        #double iter_rec(nsppol, nstep, six) ;
        #int nk_spin(nsppol) ;
        #int nq_spin(nsppol) ;
        #int nb_spin(nsppol) ;
        #double kpts_spin(nsppol, max_nk, three) ;
        #double qpts_spin(nsppol, max_nq, three) ;
        #double a_spin(nsppol, max_nk, max_nb, two) ;
        #double b_spin(nsppol, max_nq, natom3, two) ;

        # Read important dimensions.
        self.nsppol = self.read_dimvalue("nsppol")
        self.nk_spin = self.read_value("nk_spin")
        self.nq_spin = self.read_value("nq_spin")
        #self.nb_spin = self.read_value("nb_spin")
        #self.nkbz = self.read_dimvalue("gstore_nkbz")
        #self.nkibz = self.read_dimvalue("gstore_nkibz")
        #self.nqbz = self.read_dimvalue("gstore_nqbz")
        #self.nqibz = self.read_dimvalue("gstore_nqibz")
        self.varpeq_pkind = self.read_string("varpeq_pkind")
        self.ngqpt = self.read_value("gstore_ngqpt")

        # Read important variables.
        #self.completed = self.read_value("gstore_completed")
        #self.done_spin_qbz = self.read_value("gstore_done_qbz_spin")
        #self.qptopt = self.read_value("gstore_qptopt")
        #self.kptopt = self.read_value("kptopt")
        #self.kzone = self.read_string("gstore_kzone")
        #self.qzone = self.read_string("gstore_qzone")
        #self.kfilter = self.read_string("gstore_kfilter")
        #self.gmode = self.read_string("gstore_gmode")

        # Note conversion Fortran --> C for the isym index.
        self.brange_spin = self.read_value("brange_spin")
        self.brange_spin[:,0] -= 1
        self.nb_spin = self.brange_spin[:,1] - self.brange_spin[:,0]
        #self.erange_spin = self.read_value("gstore_erange_spin")
        # Total number of k/q points for each spin after filtering (if any)
        #self.glob_spin_nq = self.read_value("gstore_glob_nq_spin")
        #self.glob_nk_spin = self.read_value("gstore_glob_nk_spin")

        # K-points and q-points in the IBZ
        #self.kibz = self.read_value("reduced_coordinates_of_kpoints")
        #self.qibz = self.read_value("gstore_qibz")

        # K-points and q-points in the BZ
        #self.kbz = self.read_value("gstore_kbz")
        #self.qbz = self.read_value("gstore_qbz")

        # Mapping BZ --> IBZ. Note conversion Fortran --> C for the isym index.
        # nctkarr_t("gstore_kbz2ibz", "i", "six, gstore_nkbz"), &
        # nctkarr_t("gstore_qbz2ibz", "i", "six, gstore_nqbz"), &
        #self.kbz2ibz = self.read_value("gstore_kbz2ibz")
        #self.kbz2ibz[:,0] -= 1

        #self.qbz2ibz = self.read_value("gstore_qbz2ibz")
        #self.qbz2ibz[:,0] -= 1

        # Mapping q/k points in gqk --> BZ. Note conversion Fortran --> C for indexing.
        # nctkarr_t("gstore_qglob2bz", "i", "gstore_max_nq, number_of_spins"), &
        # nctkarr_t("gstore_kglob2bz", "i", "gstore_max_nk, number_of_spins") &
        #self.qglob2bz = self.read_value("gstore_qglob2bz")
        #self.qglob2bz -= 1

        #self.kglob2bz = self.read_value("gstore_kglob2bz")
        #self.kglob2bz -= 1

    #def find_iq_glob_qpoint(self, qpoint, spin: int):
    #    """Find the internal index of the qpoint needed to access the gvals array."""
    #    qpoint = np.asarray(qpoint)
    #    for iq_g, iq_bz in enumerate(self.qglob2bz[spin]):
    #        if np.allclose(qpoint, self.qbz[iq_bz]):
    #            return iq_g, qpoint

    #    raise ValueError(f"Cannot find {qpoint=} in GSTORE.nc")

    #def find_ik_glob_kpoint(self, kpoint, spin: int):
    #    """Find the internal indices of the kpoint needed to access the gvals array."""
    #    kpoint = np.asarray(kpoint)
    #    for ik_g, ik_bz in enumerate(self.kglob2bz[spin]):
    #        if np.allclose(kpoint, self.kbz[ik_bz]):
    #            return ik_g, kpoint

    #    raise ValueError(f"Cannot find {kpoint=} in GSTORE.nc")


class VarpeqRobot(Robot, RobotWithEbands):
    """
    This robot analyzes the results contained in multiple VARPEQ.nc files.

    Usage example:

    .. code-block:: python

        robot = VarpeqRobot.from_files([
            "out1_VARPEQ.nc",
            "out2_VARPEQ.nc",
            ])

        robot.neq(verbose=1)

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: VarpeqRobot
    """

    EXT = "VARPEQ"

    def neq(self, ref_basename: str | None = None, verbose: int = 0) -> int:
        """
        Compare all VARPEQ.nc files stored in the VarpeqRobot
        """
        # Find reference abifile. By default the first file in the robot is used.
        ref_varpeq = self._get_ref_abifile_from_basename(ref_basename)

        exc_list = []
        ierr = 0
        for other_varpeq in self.abifiles:
            if ref_varpeq.filepath == other_varpeq.filepath:
                continue
            print("Comparing: ", ref_varpeq.basename, " with: ", other_varpeq.basename)
            try:
                ierr += self._neq_two_varpeqs(ref_varpeq, other_varpeq, verbose)
                cprint("EQUAL", color="green")
            except Exception as exc:
                exc_list.append(str(exc))

        for exc in exc_list:
            cprint(exc, color="red")

        return ierr

    def _neq_two_varpeqs(self, varpeq1: VarpeqFile, varpeq2: VarpeqFile, verbose: int) -> int:
        """
        Helper function to compare two VARPEQ files.
        """
        # These quantities must be the same to have a meaningfull comparison.
        aname_list = ["structure", "nsppol",
                      #"nkbz", "nkibz",
                      #"nqbz", "nqibz", "completed", "kzone", "qzone", "kfilter", "gmode",
                      #"brange_spin", "erange_spin", "glob_spin_nq", "glob_nk_spin",
                     ]

        for aname in aname_list:
            self._compare_attr_name(aname, varpeq1, varpeq2)

        # Now compare the gkq objects for each spin.
        ierr = 0
        #for spin in range(varpeq1.nsppol):
        #    gqk1, gqk2 = varpeq1.gqk_spin[spin], varpeq2.gqk_spin[spin]
        #    ierr += gqk1.neq(gqk2, verbose)
        return ierr

    def get_kdata_spin(self, spin: int) -> dict:
        """
        """
        # First of all sort the files in reverse order using the total number of k-points in the mesh.
        def sort_func(abifile):
            ksampling = abifile.ebands.kpoints.ksampling
            ngkpt, shifts = ksampling.mpdivs, ksampling.shifts
            return np.prod(ngkpt)

        labels, abifiles, nktot_list = self.sortby(sort_func, reverse=True, unpack=True)
        data = defaultdict(list)

        # Now loop over the sorted files and extract the results of the final iteration.
        for i, (label, abifile, nktot) in zip(labels, abifiles, nktot_list):
            for k, v in abifile.get_last_iteration_dict_ev(spin).items():
                data[k].append(v)

            ksampling = abifile.ebands.kpoints.ksampling
            ngkpt, shifts = ksampling.mpdivs, ksampling.shifts
            minibz_vol = abifile.structure.reciprocal_lattice.volume / np.prod(ngkpt)
            data["ngkpt"].append(ngkpt)
            data["minibz_vol"].append(minibz_vol)

        for k in data:
            data[k] = np.array(data[k])

        return data

    @add_fig_kwargs
    def plot_scf_cycle(self, **kwargs) -> Figure:
        """
        Plot the VARPEQ SCF cycle for all the files stored in the Robot.
        """
        nsppol = self.getattr_alleq("nsppol")

        # Build grid of plots.
        nrows, ncols = nsppol * len(self), 2
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=True, squeeze=False)

        for ifile, abifile in enumerate(self.abifiles):
            row_start = nsppol * ifile
            row_stop = row_start + nsppol
            this_ax_mat = ax_mat[row_start:row_stop]
            abifile.plot_scf_cycle(ax_mat=this_ax_mat, show=False)

        return fig

    @add_fig_kwargs
    def plot_kdata(self, fontsie=12, **kwargs) -> Figure:
        """
        """
        nsppol = self.getattr_alleq("nsppol")

        # Build grid of plots.
        nrows, ncols = len(ITER_LABELS), nsppol
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=True, squeeze=False)
        deg = 1
        for spin in range(nsppol):
            kdata = self.get_kdata_spin(spin)
            xs = kdata["minibz_vol"]
            xvals = np.linspace(0, 1.1 * xs.max(), 100)
            for ix, label in enumerate(ITER_LABELS):
                ys = kdata[label]
                p = np.poly1d(np.polyfit(xs, ys, deg))
                ax = ax_mat[ix,spin]
                ax.scatter(xs, ys, marker="o")
                ax.plot(xvals, p[xvals], style="k--")
                ax.grid(True)
                ax.legend(loc="best", shadow=True, fontsize=fontsize)
                #ax.set_xlabel("Iteration", fontsize=fontsize)
                #ax.set_ylabel("Energy (eV)" if iax == 0 else r"$|\Delta|$ Energy (eV)", fontsize=fontsize)

        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """
        yield self.plot_scf_cycle(show=False)

    def write_notebook(self, nbpath=None) -> str:
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporary file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        args = [(l, f.filepath) for l, f in self.items()]
        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("robot = abilab.VarpeqRobot(*%s)\nrobot.trim_paths()\nrobot" % str(args)),
            #nbv.new_code_cell("ebands_plotter = robot.get_ebands_plotter()"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


"""
def plot_data(kleninv, energy, p, str_label):
    xrange = np.linspace(0, 1/6, 100)
    plt.plot(kleninv, energy, 's-', label=str_label)
    plt.plot(xrange, p(xrange), 'k--')
    plt.xlim(0, 0.6)
    plt.xlabel('Inverse k-grid')

def interp_data(kleninv, energy, n):
    fit = np.polyfit(kleninv[-n:], energy[-n:], 1)
    p = np.poly1d(fit)
    return p

def transform_data(klen, enpol, eps):
    kleninv = 1/klen
    enpol = (enpol - cbm)*ha_ev
    eps = -(eps - cbm)*ha_ev
    return kleninv, enpol, eps

def analyze(filename, mode, label):
    klen, enpol, eps = get_data(filename)
    kleninv, enpol, eps = transform_data(klen, enpol, eps)
    p_enpol = interp_data(kleninv, enpol, 3)
    p_eps = interp_data(kleninv, eps, 3)

    if mode == 'enpol':
        plot_data(kleninv, enpol, p_enpol, label)
    elif mode == 'eps':
        plot_data(kleninv, eps, p_eps, label)

analyze('energy.dat', 'enpol', 'no sym')
analyze('energy_ksym.dat', 'enpol', r'$g(Sk,q) = g(k,S^{-1}q)$')
plt.ylabel('Hole polaron formation energy (eV)')
"""

