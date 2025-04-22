# coding: utf-8
"""Work subclasses related to GS calculations."""
from __future__ import annotations

import sys
import json
import itertools
#import pickle
import dataclasses
import numpy as np
import pandas as pd

#from monty.json import MSONable
from monty.string import list_strings #, marquee
from pymatgen.analysis.eos import EOS
from abipy.core.structure import Structure
from abipy.tools.numtools import build_mesh
from abipy.tools.derivatives import central_fdiff_weights, check_num_points_for_order # finite_diff
from abipy.tools.typing import Figure
from abipy.abio.inputs import AbinitInput
from abipy.tools.serialization import HasPickleIO
from abipy.tools.plotting import (set_axlims, add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt,
    get_ax3d_fig_plt, rotate_ticklabels, set_visible, set_ax_xylabels)

#from abipy.electrons.gsr import GsrRobot
#from abipy.tools.serialization import mjson_write, pmg_serialize
from .works import Work

def centered_indices(n):
    half = n // 2
    if n % 2 == 0:
        return list(range(-half, half)), half
    #else:
    #return list(range(-half, half + 1)),


class FiniteDiffForcesData(HasPickleIO):

    def __init__(self,
                 iatom,
                 step,
                 mesh_type,
                 deltas,
                 ix0,
                 energies_ev,
                 structures,
                 cart_forces_list,
                 cart_stress_tensor_list,
                 ):

        self.iatom = iatom
        self.step = step
        self.mesh_type = mesh_type
        self.deltas = deltas
        self.ix0 = ix0
        self.energies_ev = energies_ev
        self.structures = structures
        self.cart_forces_list = cart_forces_list
        self.cart_stress_tensor_list = cart_stress_tensor_list

    def get_results_nn_name(self, nn_name: str, with_delta: bool) -> dict:
        """
        Args:
            nn_names: String or list of strings defining the NN potential. See also CalcBuilder.
        """
        from abipy.ml.aseml import CalcBuilder
        calc = CalcBuilder(nn_name).get_calculator()
        energies_ev, forces_list, stress_list = [], [], []
        for structure in self.structures:
            atoms = structure.to_ase_atoms(calc=calc)
            energies_ev.append(float(atoms.get_potential_energy()))
            forces_list.append(atoms.get_forces())
            stress_list.append(atoms.get_stress(voigt=False))

        return dict(energies_ev=np.array(energies_ev),
                    forces_list=np.array(forces_list),
                    stress_list=np.array(stress_list),
                    )

    #def compare_with_nn_names(self, nn_names):
    #    for nn_name in list_strings(nn_names):
    #        d = self.get_results_nn_name(nn_name)



class FiniteDiffForcesWork(Work):
    """
    Work for the computation of forces with finite differences.
    """

    @classmethod
    def from_scf_input(cls, scf_input, iatom, direction,
                       frac_coords=True, num_points=5, step=0.01, manager=None):
        """
        Build the work from an AbinitInput representing a GS SCF calculation.

        Args:
            scf_input: AbinitInput for GS SCF used as template to generate the other inputs.
            iatom: Index of the atom to displace.
            frac_coords:
                Boolean stating whether the vector corresponds to fractional or
                Cartesian coordinates.
            num_points:
            step:
            manager: TaskManager instance. Use default manager if None.
        """
        work = cls(manager=manager)

        work.iatom = int(iatom)
        work.step = float(step)
        work.mesh_type = mesh_type
        work.deltas, work.ix0 = build_mesh(0.0, num_points, step, "=")

        structure = scf_input.structure
        norm = structure.lattice.norm(direction, frac_coords=frac_coords)
        versor = np.array(direction) / norm

        for delta in work.deltas:
            vector = delta * versor
            new_structure = structure.copy()
            new_structure.translate_sites([work.iatom], delta * versor, frac_coords=frac_coords, to_unit_cell=False)
            new_input = scf_input.new_with_structure(new_structure)
            work.register_scf_task(new_input)

        return work

    def get_data(self) -> FiniteDiffForcesData:
        """
	    Read data from the GSR files, and produce a JSON file
        in the outdata directory of the work.
        """
        energies_ev, cart_forces_list, cart_stress_tensor_list = [], [], []
        for task in self:
            with task.open_gsr() as gsr:
                energies_ev.append(float(gsr.energy))
                cart_forces_list.append(np.array(gsr.cart_forces))
                cart_stress_tensor_list.append(np.array(gsr.cart_stress_tensor))

        return FiniteDiffForcesData(
            iatom=self.iatom,
            step=self.step,
            mesh_type=self.mesh_type,
            deltas=self.deltas,
            ix0=self.ix0,
            energies_ev=energies_ev,
            structures=[task.input.structure for task in self],
            cart_forces_list=cart_forces_list,
            cart_stress_tensor_list=cart_stress_tensor_list,
        )

    def on_all_ok(self):
        """
        This method is called when all tasks reach S_OK.
        """
        data = self.get_data()
        data.pickle_dump(self.outdir.path)
        #mjson_write(data, self.outdir.path_in("forces.json"))
        return super().on_all_ok()


#class FiniteDiffStressWork(Work):
#    """
#    Work for the computation of the stress tensor with finite differences
#    """
#
#    @classmethod
#    def from_scf_input(cls, scf_input, delta=1e-4, ecutsm=0.5, manager=None):
#        """
#        Build a Work from an AbinitInput representing a GS SCF calculation.
#
#        Args:
#           scf_input: AbinitInput for GS SCF used as template to generate the other inputs.
#	        delta:
#           ecutsm: Value of ecutsm input variable. If `scf_input` does not provide ecutsm, this
#               value will be used else the vale in `scf_input`.
#           manager: TaskManager instance. Use default if None.
#        """
#	    scf_input = scf_input.deepcopy()
#
#        if "ecutsm" not in scf_input:
#            scf_input.set_vars(ecutsm=ecutsm)
#            print("Input does not define ecutsm input variable.\n",
#                  "A default value of %s will be added to all the EOS inputs" % ecutsm)
#
#        new_work = cls(manager=manager)
#        self.unstrained_task = new_work.register_scf_task(scf_input)
#
#        voigt_inds = [(0, 0)]
#
#        self.delta = float(delta)
#        self.tasks_voigt = {}
#        for voigt in voigt_inds:
#            self.tasks_voigt[voigt] = defaultdict(list)
#            for isign in (-1, +1):
#            # Apply strain to the lattice.
#                    strain = np.zeros((3, 3))
#                    strain[voigt] = float(isign) * delta
#                    new_structure = scf_input.structure.deepcopy()
#                    new_structure.apply_strain(strain)
#                    new_input = scf_input.new_with_structure(new_structure)
#
#                    # Perform GS calculations with strained cell.
#                    task = new_work.register_scf_task(new_input)
#                    self.tasks_voigt[voigt].append(task)
#
#            return new_work
#
#    def get_data(self):
#        """
#        This method is called when all tasks reach S_OK.
#	     It reads the energies and the volumes from the GSR file, computes the EOS and produce a
#        JSON file `eos_data.json` in outdata.
#        """
#        with self.unstrained_task.open_gsr() as gsr0:
#            e0, v0 = gsr0.energy, gsr0.structure.volume
#            cart_stress_tensor = gsr0.cart_stress_tensor
#
#        fd_tensor = np.empty((3, 3))
#        for voigt in voigt_inds:
#            tasks = self.tasks_voigt[voigt]
#            energies_ev, volumes = np.empty(len(tasks), np.empty(len(tasks))
#            for i, task in enumerate(task):
#                    with task.open_gsr() as gsr:
#                       energies_ev[i] = float(gsr.energy))
#
#            d = (e_plus - e_minus) / (2 * self.delta * v0)
#                fd_tensor[voigt] = d
#                fd_tensor[voigt[1], voigt[0]] = d
#
#        data = {
#                 "dfpt_cart_stress_tensor": cart_stress_tensor,
#                 "finite_diff_cart_stress_tensor": fd_cart_stress_tensor,
#            }
#
#            with open(self.outdir.path_in("stress.json"), "wt") as fh:
#                json.dump(data, fh, indent=4, sort_keys=True)
#
#            return data
#
#    def on_all_ok(self):
#        """
#        This method is called when all tasks reach S_OK. It reads the energies
#        and the volumes from the GSR file, computes the EOS and produce a
#        JSON file `eos_data.json` in outdata.
#        """
#        self.compute_and_write_stress()
#        return dict(returncode=0, message="EOS computed and file written")


class FdDynMagneticChargeWork(Work):
    r"""
    Work for the computation of the dynamical magnetic charges with finite differences.

    The dynamical magnetic charges are defined as:

        Z_jv^m=Ω_0 (∂M_v)/(∂u_j ) = (∂F_j)/(∂H_v ) = Ω_0 (∂^2 E)/(∂H_β ∂u_i).

    Here we compute them as derivatives of forces wrt to the Zeeman magnetic field.
    Non collinear calculations with SOC are required to have non-zero results.
    """

    @classmethod
    def from_scf_input(cls,
                       scf_input: AbinitInput,
                       num_points: int,
                       delta_h: float = 0.01,
                       relax: bool = True,
                       relax_opts: dict | None = None,
                       manager=None) -> FdDynMagneticChargeWork:
        """
        Build the work from an AbinitInput representing a GS SCF calculation.

        Args:
            scf_input: AbinitInput for GS SCF calculation used as template to generate the other inputs.
            num_points: Number of points for finite difference.
            delta_h: Finite difference step for the magnetic field in a.u.
            relax: False if the initial structural relaxation should not be performed.
            relax_opts: optional dictionary with relaxation options.
            manager: TaskManager instance. Use default manager if None.
        """
        work = cls(manager=manager)
        work.delta_h = delta_h
        work.h_values, work.ih0 = build_mesh(0.0, num_points, delta_h, "=")
        work.num_points = len(work.h_values)
        check_num_points_for_order(num_points=work.num_points, order=1, kind="=")

        nspinor = scf_input.get("nspinor", 1)
        #if nspinor != 2:
        #    raise ValueError("nspinor should be 2 to have non-zero dyn magnetic charges while it is {nspinor}")

        work.h_cart_dirs = np.array([
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1),
        ])
        work.scf_input_template = scf_input.deepcopy()

        work.relax = relax
        if work.relax:
            relax_input = scf_input.make_relax_input()
            work.initial_relax_task = work.register_relax_task(relax_input)
        else:
            work._add_tasks_with_zeemanfield(scf_input.structure)
            work.relaxed_structure = scf_input.structure

        return work

    def _add_tasks_with_zeemanfield(self, structure: Structure) -> None:
        """Build new GS tasks with zeemanfield."""
        scf_input = self.scf_input_template.new_with_structure(structure)

        self.tasks_hdir_h = np.empty((len(self.h_cart_dirs), len(self.h_values)), dtype=object)
        task_h0 = None
        for hdir, h_cart_dir in enumerate(self.h_cart_dirs):
            for ih, h_val in enumerate(self.h_values):
                is_h0 = abs(h_val) < 1e-16
                new_inp = scf_input.new_with_vars(zeemanfield=h_val * h_cart_dir)
                if is_h0:
                    # Avoid computing H=0 multiple times.
                    if task_h0 is None:
                        task_h0 = self.register_scf_task(new_inp)
                    self.tasks_hdir_h[hdir, ih] = task_h0
                else:
                    self.tasks_hdir_h[hdir, ih] = self.register_scf_task(new_inp)

    def on_ok(self, sender):
        """
        This method is called when one task reaches status `S_OK`.
        """
        if self.relax and sender == self.initial_relax_task:
            # Get relaxed structure from GSR file.
            self.relaxed_structure = sender.get_final_structure()
            self._add_tasks_with_zeemanfield(self.relaxed_structure)
            self.flow.allocate(build=True)

        return super().on_ok(sender)

    def on_all_ok(self):
        data = self.get_data()
        data.pickle_dump(self.outdir.path)
        return super().on_all_ok()

    def get_data(self) -> ZeemanData:
        """
	    Read data from the GSR files, and compute Zm with finite differences.
        """
        data = {
            "input_structure": self[0].input.structure,
            "relaxed_structure": self.relaxed_structure,
            "delta_h": self.delta_h,
            "h_values": self.h_values,
            "h_cart_dirs": self.h_cart_dirs,
        }

        # Read forces and stresses from the GSR files.
        natom = len(self[0].input.structure)
        nh_dirs, nh_vals = self.nh_dirs, len(self.h_values)

        etotals_dh = np.empty((nh_dirs, nh_vals)
        eterms_dh = np.empty((nh_dirs, nh_vals), dtype=object)
        cart_forces_dh = np.empty((nh_dirs, nh_vals, natom, 3))
        cart_stresses_dh np.empty((nh_dirs, nh_vals, 3, 3))

        for hdir, h_cart_dir in enumerate(self.h_cart_dirs):
            for ih, h_val in enumerate(self.h_values):
                with self.tasks_hdir_h[hdir, ih].open_gsr() as gsr:
                    etotals_dh[hdir, ih] = gsr.r.read_value("etotal")
                    eterms_dh[hdir, ih] = gsr.read_energy_terms(unit="Ha")
                    cart_forces_dh[hdir, ih] = gsr.r.read_value("cartesian_forces") # Ha/Bohr units.
                    cart_stresses_dh[hdir, ih] = gsr.r.read_cart_stress_tensor(units="au")

        data["etotals_dh"] = etotals_dh
        data["eterms_dh"] = eterms_dh
        data["cart_forces_dh"] = cart_forces_dh
        data["cart_stresses_dh"] = cart_stresses_dh

        # Finite difference: ∂F_i/∂H_β
        # Use all stencils compatible with input num_points so that we can monitor the convergence.
        data["zm_npts_adh"] = {}
        data["piezom_npts_hij"] = {}
        for acc, weights in central_fdiff_weights[1].items():
            if self.num_points < len(weights): continue
            nn = acc // 2
            zm_adh = np.zeros((natom, 3, 3))
            for iat, iat_dir, hdir in itertools.product(range(natom), range(3), range(3)):
                fvals_h = cart_forces_dh[hdir, :, iat, iat_dir]
                zm_adh[iat, iat_dir, hdir] = np.sum(fvals_h[self.ih0-nn:self.ih0+nn+1] * weights) / self.delta_h
            data["zm_npts_adh"][len(weights)] = zm_adh

            # Now the stresses compute the piezomagnetic tensor.
            # e_ijk=(∂M_i)/(∂ε_jk )=(∂^2 E)/(∂ε_jk ∂H_i )
            piezom_hij = np.zeros((nh_dirs, 3, 3))
            for ii, jj, hdir in itertools.product(range(3), range(3), range(nh_dirs)):
                svals_h = cart_stresses_dh[hdir, :, ii, jj] * self.relaxed_structure.volume # * abu.Angs2Bohr ** 3
                piezom_hij[hdir, ii, jj] = np.sum(svals_h[self.ih0-nn:self.ih0+nn+1] * weights) / self.delta_h
            data["piezom_npts_hij"][len(weights)] = piezom_hij

        return ZeemanData(**data)


def idir2s(idir: int):
    """Convert direction index to string."""
    return {0: "x", 1: "y", 2: "z"}[idir]


@dataclasses.dataclass(kw_only=True)
class ZeemanData(HasPickleIO):
    """
    This object stores the dynamical magnetic charges Zm computed with finite differences.
    All values are in a.u. and Zm are in Cartesian coordinates.
    """
    input_structure: Structure
    relaxed_structure: Structure
    delta_h: float
    h_values: np.ndarray
    h_cart_dirs: np.ndarray
    etotals_dh: np.ndarray
    eterms_dh: np.ndarray
    cart_forces_dh: np.ndarray
    cart_stresses_dh: np.ndarray
    zm_npts_adh: dict[np.array]      # Mapping npts -> Zm[iat, iat_dir, hdir] in Cart. coords.
    piezom_npts_hij: dict[np.array]  # Mapping npts -> Zm[iat, iat_dir, hdir] in Cart. coords.

    @propery
    def nh_dirs(self) -> int:
        """Number of field directions."""
        return len(self.h_cart_dirs)

    def get_zmdf_iatom(self, iatom: int) -> pd.Dataframe:
        """
        Return dataframe with Zm values for the given atom index.
        """
        components = "xx xy xz yx yy yz zx zy zz".split()
        site = self.relaxed_structure[iatom]
        rows = []
        for npts, zm_adh in self.zm_npts_adh.items():
            d = {"npts": npts}
            d.update({c: v for c, v in zip(components, zm_adh[iatom].flatten(), strict=True)})
            d["trace"] = np.trace(zm_adh[iatom])
            d["det"] = np.linalg.det(zm_adh[iatom])
            rows.append(d)

        return pd.DataFrame(rows)

    def print(self, elements: None | list[str] = None, file=sys.stdout) -> None:
        """
        Print Zm to `file`. Show only elements in `elements` if not None.
        """
        def _p(*args, **kwargs):
            print(*args, file=file, **kwargs)

        if elements is not None: elements = list_strings(elements)

        _p("Input structure:")
        _p(self.input_structure)
        _p("")
        _p("Relaxed structure:")
        _p(self.relaxed_structure)
        _p("")

        for iatom, site in enumerate(self.relaxed_structure):
            if elements is not None and site.species_string not in elements: continue
            df = self.get_zmdf_iatom(iatom)
            _p(f"Zm[atom_dir, H_dir] in Cart. coords for {iatom=}: element: {site.species_string}, frac_coords: {site.frac_coords}")
            _p(df)
            _p("")

    @add_fig_kwargs
    def plot_etotal(self, ax=None, fontsize=8, **kwargs) -> Figure:
        """Plot energies as a function of H."""
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        for hdir, h_cart_dir in enumerate(self.h_cart_dirs):
            ax.plot(self.delta_h, self.etotals[hdir], marker="o", label=f"H_dir: {idir2s(hdir)}")
        ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    @add_fig_kwargs
    def plot_forces(self, elements=None, fontsize=8, **kwargs) -> Figure:
        """Plot Cartesian forces as a function of H."""
        if elements is not None: elements = list_strings(elements)
        nrows, ncols = 3, self.nh_dirs
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=False, squeeze=False)
        for iat_dir in range(3):
            for hdir, h_cart_dir in enumerate(self.h_cart_dirs):
                ax = ax_mat[iat_dir, hdir]
                ax.set_title(f"H_dir: {idir2s(hdir)}, Atom_dir: {idir2s(iat_dir)}", fontsize=fontsize)
                for iat, site in enumerate(self.relaxed_structure):
                    if elements is not None and site.species_string not in elements: continue
                    ax.plot(self.h_values, self.cart_forces_dh[hdir, :, iat, iat_dir], marker="o",
                            label=site.species_string} + r"$_{\text{%s}}$" % iat,
                    )
                ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    @add_fig_kwargs
    def plot_stresses(self, fontsize=8, **kwargs) -> Figure:
        """Plot Cartesian stresses as a function of H."""
        nrows, ncols = self.nh_dirs, 1
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=False, squeeze=False)
        for hdir, h_cart_dir in enumerate(self.h_cart_dirs):
            ax = ax_mat[hdir]
            ax.set_title(f"H_dir: {idir2s(hdir)}", fontsize=fontsize)
            for ii, jj in itertools.product(range(3), range(3)):
                ax.plot(self.h_values, self.cart_stresses_dh[hdir, :, ii, jj], marker="o",
                        label="$\sigma_{%s}$" % (f"{ii}, {jj}"),
                        )
            ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig
