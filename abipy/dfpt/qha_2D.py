"""
Added to compute ZSISA-QHA for systems with two degrees of freedom (2DOF).
Capable of calculating anisotropic thermal expansion and lattice constants for uniaxial configurations.
Requires PHDOS.nc and DDB files for GSR calculations or _GSR.nc files.
If PHDOS.nc is available for all structures, normal interpolation for QHA will be applied.
Supports the use of six PHDOS.nc files for specific structures to employ the EinfVib2 approximation.
"""
from __future__ import annotations

import os
import abc
import numpy as np
import abipy.core.abinit_units as abu

from scipy.interpolate import RectBivariateSpline #, RegularGridInterpolator
from functools import cached_property
#from monty.collections import dict2namedtuple
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt
from abipy.tools.typing import Figure
from abipy.tools.serialization import HasPickleIO, mjson_load
from abipy.electrons.gsr import GsrFile
from abipy.dfpt.ddb import DdbFile
from abipy.dfpt.phonons import PhdosFile # PhononBandsPlotter, PhononDos,
from abipy.dfpt.vzsisa import anaget_phdoses_with_gauss



class QHA_2D(HasPickleIO):
    """
    Quasi-Harmonic Approximation (QHA) analysis in 2D.
    Provides methods for calculating and visualizing energy, free energy, and thermal expansion.
    """

    @classmethod
    def from_json_file(cls,
                       filepath: PathLike,
                       nqsmall_or_qppa: int,
                       anaget_kwargs: dict | None = None,
                       smearing_ev: float | None = None,
                       verbose: int = 0) -> Vzsisa:
        """
        Build an instance from a json file `filepath` typically produced by an Abipy flow.
        For the meaning of the other arguments see from_gsr_ddb_paths.
        """
        data = mjson_load(filepath)
        return cls.from_gsr_ddb_paths(nqsmall_or_qppa,
                                      data["gsr_relax_paths"], data["ddb_relax_paths"],
                                      data["bo_strains_ac"], data["phdos_strains_ac"],
                                      anaget_kwargs=anaget_kwargs, smearing_ev=smearing_ev, verbose=verbose)

    @classmethod
    def from_gsr_ddb_paths(cls,
                           nqsmall_or_qppa: int,
                           gsr_paths,
                           ddb_paths,
                           bo_strains_ac,
                           phdos_strains_ac,
                           anaget_kwargs: dict | None = None,
                           smearing_ev: float | None = None,
                           verbose: int = 0) -> QHA_2D:
        """
        Creates an instance from a list of GSR files and a list of DDB files.
        This is a simplified interface that computes the PHDOS.nc files automatically
        from the DDB files by invoking anaddb.

        Args:
            nqsmall_or_qppa: Define the q-mesh for the computation of the PHDOS.
                if > 0, it is interpreted as nqsmall.
                if < 0, it is interpreted as q-point per atom qppa.
            gsr_paths: list of paths to GSR files.
            ddb_paths: list of paths to DDB files.
            bo_strains_ac: List of strains for the a and the c lattice vector.
            phdos_strains_ac: List of strains for the a and the c lattice vector.
            anaget_kwargs: dict with extra arguments passed to anaget_phdoses_with_gauss.
            smearing_ev: Gaussian smearing in eV.
            verbose: Verbosity level.
        """
        phdos_paths, phbands_paths = anaget_phdoses_with_gauss(nqsmall_or_qppa, smearing_ev, ddb_paths, anaget_kwargs, verbose)

        new = cls.from_files(ddb_paths, phdos_paths_2D, bo_strains_ac, phdos_strains_ac, gsr_file="GSR.nc")
        #new.pickle_dump(workdir, basename=None)
        return new

    @classmethod
    def from_files(cls, gsr_paths_2D, phdos_paths_2D, bo_strains_ac, phdos_strains_ac, gsr_file="GSR.nc") -> QHA_2D:
        """
        Creates an instance of QHA from a 2D array of GSR and PHDOS files.

        Args:
            gsr_paths_2D: 2D list of paths to GSR files.
            phdos_paths_2D: 2D list of paths to PHDOS.nc files.
            bo_strains_ac: List of strains for the a and the c lattice vector.
            phdos_strains_ac: List of strains for the a and the c lattice vector.
        """
        energies, structures, phdoses , structures_from_phdos = [], [], [],[]

        #shape = (len(strains_a), len(strains_c))
        #gsr_paths_2d = np.reshape(gsr_paths_2D, shape)
        #phdos_paths_2d = np.reshape(phdos_paths_2D, shape)

        if gsr_file == "GSR.nc":
            # Process GSR files
            for row in gsr_paths_2D:
                row_energies, row_structures = [], []
                for gp in row:
                    if os.path.exists(gp):
                        with GsrFile.from_file(gp) as g:
                            row_energies.append(g.energy)
                            row_structures.append(g.structure)
                    else:
                        row_energies.append(None)
                        row_structures.append(None)
                energies.append(row_energies)
                structures.append(row_structures)

        elif gsr_file == "DDB":
            # Process DDB files
            for row in gsr_paths_2D:
                row_energies, row_structures = [], []
                for gp in row:
                    if os.path.exists(gp):
                        with DdbFile.from_file(gp) as g:
                            row_energies.append(g.total_energy)
                            row_structures.append(g.structure)
                    else:
                        row_energies.append(None)
                        row_structures.append(None)
                energies.append(row_energies)
                structures.append(row_structures)

        else:
            raise ValueError(f"Invalid {gsr_file=}")

        # Process PHDOS files
        for row in phdos_paths_2D:
            row_doses , row_structures = [],[]
            for path in row:
                if os.path.exists(path):
                    with PhdosFile(path) as p:
                        row_doses.append(p.phdos)
                        row_structures.append(p.structure)
                else:
                    row_doses.append(None)
                    row_structures.append(None)

            phdoses.append(row_doses)
            structures_from_phdos.append(row_structures)

        return cls(structures, phdoses, energies, structures_from_phdos, bo_strains_ac, phdos_strains_ac)

    def __init__(self, structures, phdoses, energies, structures_from_phdos,
                 bo_strains_ac, phdos_strains_ac,
                 eos_name: str='vinet', pressure: float=0.0):
        """
        Args:
            structures (list): List of structures at different volumes.
            phdoses: List of density of states (DOS) data for phonon calculations.
            energies (list): SCF energies for the structures in eV.
            bo_strains_ac: List of strains for the a and the c lattice vector.
            phdos_strains_ac: List of strains for the a and the c lattice vector.
            eos_name (str): Expression used to fit the energies (e.g., 'vinet').
            pressure (float): External pressure in GPa to include in p*V term.
        """
        self.phdoses = phdoses
        self.structures = structures
        self.structures_from_phdos = structures_from_phdos
        self.energies = np.array(energies, dtype=np.float64)

        self.bo_strains_ac = bo_strains_ac
        self.phdos_strains_ac = phdos_strains_ac

        self.eos_name = eos_name
        self.pressure = pressure
        self.volumes = np.array([[s.volume if s else np.nan for s in row] for row in structures])
        energies_array = np.array(energies)
        energies_array[energies_array == None] = np.nan

        # Find the indices of the minimum values
        self.ix0, self.iy0 = np.unravel_index(np.nanargmin(energies_array), energies_array.shape)

        # Extract lattice parameters and angles
        self.lattice_a = np.array([[s.lattice.abc[0] if s is not None else None for s in row] for row in structures])
        self.lattice_c = np.array([[s.lattice.abc[2] if s is not None else None for s in row] for row in structures])

        self.lattice_a_from_phdos = np.array([[s.lattice.abc[0] if s is not None else None for s in row] for row in structures_from_phdos])
        self.lattice_c_from_phdos = np.array([[s.lattice.abc[2] if s is not None else None for s in row] for row in structures_from_phdos])

        # Find index of minimum energy
        self.min_energy_idx = np.unravel_index(np.nanargmin(self.energies), self.energies.shape)

    @cached_property
    def use_qha(self) -> bool:
        """True if we are in full QHA_2D mode."""
        return len(self.lattice_a_from_phdos) == len(self.lattice_a) and len(self.lattice_c_from_phdos) == len(self.lattice_c)

    @cached_property
    def use_einfvib2(self) -> bool:
        return len(self.lattice_a_from_phdos) == 3 and len(self.lattice_c_from_phdos) == 3

    def get_initial_guess_ac(self) -> np.array:
        """Return the initial guess for (a, c):"""
        initial_guess = [1.005 * self.lattice_a[self.ix0, 0], 1.005 * self.lattice_c[0,self.iy0]]
        return np.array(initial_guess)

    @add_fig_kwargs
    def plot_energies(self, ax=None, **kwargs) -> Figure:
        """
        Plot BO energy surface and visualize minimum in a 3D plot.

        Args:
            ax: Matplotlib axis for the plot. If None, creates a new figure.
        """
        ax, fig, plt = get_ax_fig_plt(ax, figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')  # Create a 3D subplot

        a0 = self.lattice_a[:,0]
        c0 = self.lattice_c[0,:]

        X, Y = np.meshgrid(c0, a0)

        # Plot the surface
        ax.plot_wireframe(X, Y, self.energies, cmap='viridis')
        ax.scatter(self.lattice_c[0,self.iy0], self.lattice_a[self.ix0,0], self.energies[self.ix0, self.iy0], color='red', s=100)

        f_interp = RectBivariateSpline(a0, c0, self.energies, kx=4, ky=4)

        xy_init = self.get_initial_guess_ac()

        min_x0, min_y0, min_energy = self.find_minimum( f_interp, xy_init, tol=1e-6, max_iter=1000, step_size=0.01)

        x_new = np.linspace(min(self.lattice_a[:,0]), max(self.lattice_a[:,0]), 100)
        y_new = np.linspace(min(self.lattice_c[0,:]), max(self.lattice_c[0,:]), 100)
        x_grid, y_grid = np.meshgrid(y_new, x_new)

        energy_interp = f_interp(x_new, y_new)

        ax.plot_surface(x_grid, y_grid, energy_interp, cmap='viridis', alpha=0.6)

        # Set labels
        ax.set_xlabel('Lattice parameter C (Å)')
        ax.set_ylabel('Lattice parameter A (Å)')
        ax.set_zlabel('Energy (eV)')
        ax.set_title('BO Energy Surface in 3D')

        return fig

    def find_minimum(self, f_interp, xy_init, tol=1e-6, max_iter=1000, step_size=0.01) -> tuple:
        """
        Gradient descent to find the minimum of the interpolated BO energy surface.

        Args:
            f_interp: Interpolating function for energy.
            xy_init (list): Initial guess for [a, c].
            tol (float): Convergence tolerance for gradient norm.
            max_iter (int): Maximum number of iterations.
            step_size (float): Step size for gradient descent.

        Returns:
            tuple: Optimized [a, c] coordinates and minimum energy.
        """
        xy = np.array(xy_init)
        dx = dy = 0.001

        for it in range(max_iter):
            grad = [
                (f_interp(xy[0] + dx, xy[1]) - f_interp(xy[0] - dx, xy[1])) / (2 * dx),
                (f_interp(xy[0], xy[1] + dy) - f_interp(xy[0], xy[1] - dy)) / (2 * dy),
            ]
            xy -= step_size * np.ravel(grad)
            if np.linalg.norm(grad) < tol:
                #print(f"Converged after {it} iterations with {tol=}")
                break
        else:
            raise RuntimeError(f"Could not reach {tol=} after {max_iter=}")

        min_energy = f_interp(xy[0], xy[1])
        return xy[0], xy[1], min_energy

    @add_fig_kwargs
    def plot_free_energies(self, tstart=0, tstop=800, num=5, ax=None, **kwargs) -> Figure:
        """
        Plot free energy as a function of temperature in a 3D plot.

        Args:
            ax: Matplotlib axis for the plot.
        """
        ax, fig, plt = get_ax_fig_plt(ax, figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')  # Create a 3D subplot

        tmesh = np.linspace(tstart, tstop, num)
        ph_energies = self.get_vib_free_energies(tstart, tstop, num)
        a0 = self.lattice_a[:,0]
        c0 = self.lattice_c[0,:]

        if self.use_qha:
            tot_en = self.energies[np.newaxis, :].T + ph_energies + self.volumes[np.newaxis, :].T * self.pressure / abu.eVA3_GPa

            X, Y = np.meshgrid(self.lattice_c[0,:], self.lattice_a[:,0])
            for  e in ( tot_en.T ):
                ax.plot_surface(X, Y, e, cmap='viridis', alpha=0.7)
                ax.plot_wireframe(X, Y, e, cmap='viridis')

            xy_init = self.get_initial_guess_ac()
            min_x, min_y, min_tot_en = np.zeros(num), np.zeros(num), np.zeros(num)
            for j,e in enumerate(tot_en.T):
                f_interp = RectBivariateSpline(a0, c0, e, kx=4, ky=4)
                min_x[j], min_y[j], min_tot_en[j]= self.find_minimum(f_interp, xy_init, tol=1e-6, max_iter=1000, step_size=0.01)
                xy_init = min_x[j], min_y[j]

            ax.scatter(min_y, min_x, min_tot_en, color='c', s=100)
            ax.plot(min_y, min_x, min_tot_en, color='c')

        elif self.use_einfvib2:
            a0 = self.lattice_a[1,1]
            c0 = self.lattice_c[1,1]
            da = self.lattice_a[0,1]-self.lattice_a[1,1]
            dc = self.lattice_c[1,0]-self.lattice_c[1,1]
            dF_dA, dF_dC, d2F_dA2, d2F_dC2, d2F_dAdC = (np.zeros(num) for _ in range(5))

            for i, e in enumerate(ph_energies.T):
                dF_dA[i]=(e[0,1]-e[2,1])/(2*da)
                dF_dC[i]=(e[1,0]-e[1,2])/(2*dc)
                d2F_dA2[i]=(e[0,1]-2*e[1,1]+e[2,1])/(da)**2
                d2F_dC2[i]=(e[1,0]-2*e[1,1]+e[1,2])/(dc)**2
                d2F_dAdC[i] = (e[1,1] - e[1, 0] - e[0, 1] + e[0, 0]) / ( da * dc)

            tot_en2 = self.energies[np.newaxis, :].T + ph_energies[1,1] + self.volumes[np.newaxis, :].T * self.pressure / abu.eVA3_GPa
            tot_en2 = tot_en2+ (self.lattice_a[np.newaxis, :].T - a0)*dF_dA + 0.5*(self.lattice_a[np.newaxis, :].T - a0)**2*d2F_dA2
            tot_en2 = tot_en2+ (self.lattice_c[np.newaxis, :].T - c0)*dF_dC + 0.5*(self.lattice_c[np.newaxis, :].T - c0)**2*d2F_dC2
            tot_en2 = tot_en2+ (self.lattice_c[np.newaxis, :].T - c0)*(self.lattice_a[np.newaxis, :].T - a0)*d2F_dAdC

            a = self.lattice_a[:,0]
            c = self.lattice_c[0,:]
            a_phdos = self.lattice_a[:,0]
            c_phdos = self.lattice_c[0,:]

            xy_init = self.get_initial_guess_ac()
            min_x, min_y, min_tot_en2 = np.zeros(num), np.zeros(num), np.zeros(num)
            for j, e in enumerate(tot_en2.T):
                f_interp = RectBivariateSpline(a, c, e, kx=4, ky=4)
                min_x[j], min_y[j], min_tot_en2[j] = self.find_minimum(f_interp, xy_init, tol=1e-6, max_iter=1000, step_size=0.01)
                xy_init = min_x[j], min_y[j]

            X, Y = np.meshgrid(c, a)
            for e in tot_en2.T:
                ax.plot_wireframe(X, Y, e, cmap='viridis')
                ax.plot_surface(X, Y, e, cmap='viridis', alpha=0.7)

            ax.scatter(min_y, min_x, min_tot_en2, color='c', s=100)
            ax.plot(min_y, min_x, min_tot_en2, color='c')

        else:
            raise RuntimeError("Invalid branch")

        ax.scatter(self.lattice_c[0,self.iy0], self.lattice_a[self.ix0,0], self.energies[self.ix0, self.iy0], color='red', s=100)

        ax.set_xlabel('C')
        ax.set_ylabel('A')
        ax.set_zlabel('Free energy (eV)')
        #ax.set_title('Free energies as a 3D Plot')
        plt.savefig("energy.pdf", format="pdf", bbox_inches="tight")

        return fig

    @add_fig_kwargs
    def plot_thermal_expansion(self, tstart=0, tstop=800, num=81, ax=None, **kwargs) -> Figure:
        """
        Plots thermal expansion coefficients along the a-axis, c-axis, and volumetric alpha.
        Uses both QHA and a 9-point stencil for comparison.

        Args:
            tstart: Start temperature.
            tstop: Stop temperature.
            num: Number of temperature points.
            ax: Matplotlib axis object for plotting.
        """
        ax, fig, plt = get_ax_fig_plt(ax, figsize=(10, 8))  # Ensure a valid plot axis

        tmesh = np.linspace(tstart, tstop, num)
        ph_energies = self.get_vib_free_energies(tstart, tstop, num)
        min_x, min_y, min_tot_energy = np.zeros(num), np.zeros(num), np.zeros(num)

        if self.use_qha:
            tot_energies = self.energies[np.newaxis, :].T + ph_energies + self.volumes[np.newaxis, :].T * self.pressure / abu.eVA3_GPa

            # Initial guess for minimization
            xy_init = self.get_initial_guess_ac()

            # Perform minimization for each temperature
            for j, energy in enumerate(tot_energies.T):
                f_interp = RectBivariateSpline(self.lattice_a[:, 0], self.lattice_c[0, :], energy, kx=4, ky=4)
                min_x[j], min_y[j], min_tot_energy[j] = self.find_minimum( f_interp, xy_init, tol=1e-6, max_iter=1000, step_size=0.01)
                xy_init = min_x[j], min_y[j]

            # Calculate thermal expansion coefficients
            A0, C0 = self.lattice_a[self.ix0, self.iy0], self.lattice_c[self.ix0, self.iy0]
            scale = self.volumes[self.ix0, self.iy0] / A0**2 / C0
            min_volumes = min_x**2 * min_y * scale

            dt = tmesh[1] - tmesh[0]
            alpha_a = (min_x[2:] - min_x[:-2]) / (2 * dt) / min_x[1:-1]
            alpha_c = (min_y[2:] - min_y[:-2]) / (2 * dt) / min_y[1:-1]
            alpha_v = (min_volumes[2:] - min_volumes[:-2]) / (2 * dt) / min_volumes[1:-1]

            ax.plot(tmesh[1:-1], alpha_a, color='b', label=r"$\alpha_a$ (QHA)", linewidth=2)
            ax.plot(tmesh[1:-1], alpha_c, color='r', label=r"$\alpha_c$ (QHA)", linewidth=2)
            #ax.plot(tmesh[1:-1], alpha_v, color='purple', label=r"$\alpha_v$ (QHA)", linewidth=2)

        elif self.use_einfvib2:

            a0 = self.lattice_a_from_phdos[1,1]
            c0 = self.lattice_c_from_phdos[1,1]
            da = self.lattice_a_from_phdos[0,1]-self.lattice_a_from_phdos[1,1]
            dc = self.lattice_c_from_phdos[1,0]-self.lattice_c_from_phdos[1,1]

            dF_dA, dF_dC, d2F_dA2, d2F_dC2, d2F_dAdC = (np.zeros(num) for _ in range(5))
            for i, e in enumerate(ph_energies.T):
                dF_dA[i]=(e[0,1]-e[2,1])/(2*da)
                dF_dC[i]=(e[1,0]-e[1,2])/(2*dc)
                d2F_dA2[i]=(e[0,1]-2*e[1,1]+e[2,1])/(da)**2
                d2F_dC2[i]=(e[1,0]-2*e[1,1]+e[1,2])/(dc)**2
                d2F_dAdC[i] = (e[1,1] - e[1, 0] - e[0, 1] + e[0, 0]) / ( da * dc)

            tot_en2 = self.energies[np.newaxis, :].T + ph_energies[1,1] + self.volumes[np.newaxis, :].T * self.pressure / abu.eVA3_GPa
            tot_en2 = tot_en2+ (self.lattice_a[np.newaxis, :].T - a0)*dF_dA + 0.5*(self.lattice_a[np.newaxis, :].T - a0)**2*d2F_dA2
            tot_en2 = tot_en2+ (self.lattice_c[np.newaxis, :].T - c0)*dF_dC + 0.5*(self.lattice_c[np.newaxis, :].T - c0)**2*d2F_dC2
            tot_en2 = tot_en2+ (self.lattice_c[np.newaxis, :].T - c0)*(self.lattice_a[np.newaxis, :].T - a0)*d2F_dAdC

            gradient = np.zeros(2)

            # Initial guess for minimization
            xy_init = self.get_initial_guess_ac()

            for j, energy in enumerate(tot_en2.T):
                f_interp = RectBivariateSpline(self.lattice_a[:, 0], self.lattice_c[0, :], energy, kx=4, ky=4)
                min_x[j], min_y[j], min_tot_energy[j] = self.find_minimum(f_interp, xy_init, tol=1e-6, max_iter=1000, step_size=0.01)
                xy_init = min_x[j], min_y[j]

            A0 = self.lattice_a[self.ix0,self.iy0]
            C0 = self.lattice_c[self.ix0,self.iy0]
            scale = self.volumes[self.ix0,self.iy0]/A0**2/C0
            min_v = min_x**2*min_y*scale

            dt = tmesh[1] - tmesh[0]
            alpha_a = (min_x[2:] - min_x[:-2]) / (2 * dt) / min_x[1:-1]
            alpha_c = (min_y[2:] - min_y[:-2]) / (2 * dt) / min_y[1:-1]
            alpha_v = (min_v[2:] - min_v[:-2]) / (2 * dt) / min_v[1:-1]

            ax.plot(tmesh[1:-1], alpha_a, linestyle='--', color='gold', label=r"$\alpha_a$ E$\infty$Vib2")
            ax.plot(tmesh[1:-1], alpha_c, linestyle='--', color='teal', label=r"$\alpha_c$ E$\infty$Vib2")
            #ax.plot(tmesh[1:-1], alpha_v, linestyle='--', color='darkorange', label=r"$\alpha_v$ E$\infty$Vib2")

        else:
            raise RuntimeError("Invalid branch.")

        # Save the data
        data_to_save = np.column_stack((tmesh[1:-1], alpha_v, alpha_a, alpha_c))
        columns = ['#Tmesh', 'alpha_v', 'alpha_a', 'alpha_c']
        file_path = 'thermal-expansion_data.txt'
        print(f"Writing thermal expansion data to: {file_path}")
        np.savetxt(file_path, data_to_save, fmt='%4.6e', delimiter='\t\t',  header='\t\t\t'.join(columns), comments='')

        ax.grid(True)
        ax.legend(loc="best", shadow=True)
        ax.set_xlabel('Temperature (K)')
        ax.set_ylabel(r'Thermal Expansion Coefficients ($\alpha$)')
        plt.savefig("thermal_expansion.pdf", format="pdf", bbox_inches="tight")

        return fig

    @add_fig_kwargs
    def plot_lattice(self, tstart=0, tstop=800, num=81, ax=None, **kwargs) -> Figure:
        """
        Plots thermal expansion coefficients along the a-axis, c-axis, and volumetric alpha.
        Uses both QHA and a 9-point stencil for comparison.

        Args:
            tstart: Start temperature.
            tstop: Stop temperature.
            num: Number of temperature points.
            ax: Matplotlib axis object for plotting.
        """
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(1, 3, figsize=(18, 6), sharex=True)

        tmesh = np.linspace(tstart, tstop, num)
        ph_energies = self.get_vib_free_energies(tstart, tstop, num)
        min_x, min_y, min_tot_energy = np.zeros(num), np.zeros(num), np.zeros(num)

        if self.use_qha:
            tot_energies = self.energies[np.newaxis, :].T + ph_energies+ self.volumes[np.newaxis, :].T * self.pressure / abu.eVA3_GPa

            # Initial guess for minimization
            xy_init = self.get_initial_guess_ac()

            # Perform minimization for each temperature
            for j, energy in enumerate(tot_energies.T):
                f_interp = RectBivariateSpline(self.lattice_a[:, 0], self.lattice_c[0, :], energy, kx=4, ky=4)
                min_x[j], min_y[j], min_tot_energy[j] = self.find_minimum(f_interp, xy_init, tol=1e-6, max_iter=1000, step_size=0.01)
                xy_init = min_x[j], min_y[j]

            # Calculate thermal expansion coefficients
            A0, C0 = self.lattice_a[self.ix0, self.iy0], self.lattice_c[self.ix0, self.iy0]
            scale = self.volumes[self.ix0, self.iy0] / A0**2 / C0
            min_volumes = min_x**2 * min_y * scale

            # Plot min_x in the first subplot
            axs[0].plot(tmesh, min_x, color='c', label=r"$a$ (QHA)", linewidth=2)
            axs[1].set_title("Plots of a, c, and V (QHA)")
            axs[1].plot(tmesh, min_y, color='r', label=r"$c$ (QHA)", linewidth=2)
            axs[2].plot(tmesh, min_volumes, color='b', label=r"$V$ (QHA)", linewidth=2)

        elif self.use_einfvib2:

            a0 = self.lattice_a[1,1]
            c0 = self.lattice_c[1,1]
            da = self.lattice_a[0,1]-self.lattice_a[1,1]
            dc = self.lattice_c[1,0]-self.lattice_c[1,1]

            dF_dA, dF_dC, d2F_dA2, d2F_dC2, d2F_dAdC = (np.zeros(num) for _ in range(5))
            for i, e in enumerate(ph_energies.T):
                dF_dA[i]=(e[0,1]-e[2,1])/(2*da)
                dF_dC[i]=(e[1,0]-e[1,2])/(2*dc)
                d2F_dA2[i]=(e[0,1]-2*e[1,1]+e[2,1])/(da)**2
                d2F_dC2[i]=(e[1,0]-2*e[1,1]+e[1,2])/(dc)**2
                d2F_dAdC[i] = (e[1,1] - e[1, 0] - e[0, 1] + e[0, 0]) / ( da * dc)

            tot_en2 = self.energies[np.newaxis, :].T + ph_energies[1,1] + self.volumes[np.newaxis, :].T * self.pressure / abu.eVA3_GPa
            tot_en2 = tot_en2 + (self.lattice_a[np.newaxis, :].T - a0)*dF_dA + 0.5*(self.lattice_a[np.newaxis, :].T - a0)**2*d2F_dA2
            tot_en2 = tot_en2 + (self.lattice_c[np.newaxis, :].T - c0)*dF_dC + 0.5*(self.lattice_c[np.newaxis, :].T - c0)**2*d2F_dC2
            tot_en2 = tot_en2 + (self.lattice_c[np.newaxis, :].T - c0)*(self.lattice_a[np.newaxis, :].T - a0)*d2F_dAdC

            # Initial guess for minimization
            xy_init = self.get_initial_guess_ac()

            for j, energy in enumerate(tot_en2.T):
                f_interp = RectBivariateSpline(self.lattice_a[:, 0], self.lattice_c[0, :], energy, kx=4, ky=4)
                min_x[j], min_y[j], min_tot_energy[j] = self.find_minimum(f_interp, xy_init, tol=1e-6, max_iter=1000, step_size=0.01)
                xy_init = min_x[j], min_y[j]

            A0 = self.lattice_a[self.ix0, self.iy0]
            C0 = self.lattice_c[self.ix0, self.iy0]
            scale = self.volumes[self.ix0, self.iy0] / A0**2 / C0
            min_volumes = min_x**2 * min_y * scale

            axs[0].plot(tmesh, min_x, color='c', label=r"$a$ (E$\infty$Vib2)", linewidth=2)
            axs[1].set_title(r"Plots of a, c, and V (E$\infty$Vib2)")
            axs[1].plot(tmesh, min_y, color='r', label=r"$c$ (E$\infty$Vib2)", linewidth=2)
            axs[2].plot(tmesh, min_volumes, color='b', label=r"$V$ (E$\infty$Vib2)", linewidth=2)

        else:
            raise RuntimeError("Invalid branch.")

        axs[0].set_ylabel("a")
        axs[1].set_ylabel("c")
        axs[2].set_ylabel("Volume")

        for ax in axs:
            ax.legend(loc="best", shadow=True)
            ax.grid(True)
            ax.set_xlabel("Temperature (T)")

        # Adjust layout and show the figure
        plt.tight_layout()

        return fig

    def get_vib_free_energies(self, tstart: float, tstop: float, num: int) -> np.ndarray:
        """
        Computes the vibrational free energies from phonon density of states.

        Args:
            tstart: Start temperature.
            tstop: Stop temperature.
            num: Number of temperature points.

        Return: A 3D array of vibrational free energies of shape (num_c, num_a, num_temp)
        """
        f = np.zeros((len(self.lattice_c_from_phdos[0]), len(self.lattice_a_from_phdos[:, 0]), num))

        for i in range(len(self.lattice_a_from_phdos[:, 0])):
            for j in range(len(self.lattice_c_from_phdos[0])):
                if (phdos := self.phdoses[i][j]) is not None:
                    f[j, i] = phdos.get_free_energy(tstart, tstop, num).values
        return f
