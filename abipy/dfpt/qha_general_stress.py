"""
Computes thermal stress using EinfVib2 for specific configurations
across various crystallographic structures, from cubic to triclinic.
"""
from __future__ import annotations

import numpy as np
import sys
import os
import math
import dataclasses
import abipy.core.abinit_units as abu

from abipy.abio.enums import StrEnum
from abipy.tools.serialization import HasPickleIO, Serializable
from abipy.core.structure import Structure
from abipy.core.symmetries import AbinitSpaceGroup
from abipy.electrons.gsr import GsrFile
from abipy.dfpt.ddb import DdbFile
from abipy.dfpt.phonons import PhdosFile, PhononDosPlotter
from abipy.dfpt.vzsisa import anaget_phdoses_with_gauss


def spgnum_to_crystal_system(spgrp_number: int) -> str:
    if 1 <= spgrp_number <= 2:
        # triclinic crystal system
        sym = "triclinic"
    elif 3 <= spgrp_number <= 15:
        # monoclinic crystal system
        sym = "monoclinic"
    elif 16 <= spgrp_number <= 74:
        # orthorhombic crystal system
        sym = "orthorhombic"
    elif 75 <= spgrp_number <= 142:
        # Tetragonal crystal system
        sym = "tetragonal"
    elif 143 <= spgrp_number <= 167:
        # trigonal systems
        sym = "trigonal"
    elif 168 <= spgrp_number <= 194:
        # hexagonal crystal systems
        sym = "hexagonal"
    elif 195 <= spgrp_number <= 230:
        # cubic crystal system
        sym = "cubic"
    else:
        raise RuntimeError(f"Unknown {spgrp_number=}")

    return sym


class QhaModel(StrEnum):
    """Enumerator for the different kinds of QHA models."""
    zsisa = 'zsisa'
    v_zsisa = 'v_zsisa'
    zsisa_slab = 'zsisa_slab'


@dataclasses.dataclass(kw_only=True)
class ThermalData(Serializable):

    temperature: float           # Temperature in K
    pressure_gpa: float          # Pressure in GPa.
    converged: bool              # True if self-consistent condition in the relaxation has been reached
    dtol: np.ndarray
    stress_au: np.ndarray        # Stress in a.u. (Voigt form)
    elastic: np.ndarray | None   # Elastic constants.
    therm: np.ndarray | None     # Thermal expansion
    gibbs: float                 # Gibbs free energy in eV.

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        def format_array(arr):
            if arr is None:
                return "None"
            return np.array2string(arr, precision=10, suppress_small=True)

        return (
            f"ThermalData(\n"
            f"  temperature = {self.temperature:.3f} K,\n"
            f"  pressure_gpa = {self.pressure_gpa:.3f} GPa,\n"
            f"  converged = {self.converged},\n"
            f"  dtol = {format_array(self.dtol)},\n"
            f"  stress_au = {format_array(self.stress_au)},\n"
            #f"  elastic = {format_array(self.elastic)},\n"
            #f"  therm = {format_array(self.therm)}\n"
            f")"
        )


class QHA_ZSISA(HasPickleIO):
    """
    NB: The present implementation does not include electronic entropic contributions for metals.
    """

    @classmethod
    def from_json_file(cls,
                       json_filepath: str,
                       nqsmall_or_qppa: int,
                       anaget_kwargs: dict | None = None,
                       smearing_ev: float | None = None,
                       verbose: int = 0) -> QHA_ZSISA:
        """
        Build an instance from a json file produced by an AbiPy flow.

        Args:
            json_filepath: path to the json file.
            nqsmall_or_qppa:
            anaget_kwargs: kwargs passed to anaget_phdoses_with_gauss.
            smearing_ev: Gaussian smearing in eV. None to use Abinit default value.
            verbose: Verbosity level
        """
        json_filepath = os.path.abspath(json_filepath)
        from abipy.flowtk.zsisa import ZsisaResults
        data = ZsisaResults.json_load(json_filepath)

        phdos_paths, phbands_paths = anaget_phdoses_with_gauss(nqsmall_or_qppa, smearing_ev,
                                                               data.ddb_relax_paths, anaget_kwargs, verbose)

        # Create a 6D array initialized with None.
        phdos_paths_6d = np.full((3, 3, 3, 3, 3, 3), None, dtype=object)

        for phdos_path, ind in zip(phdos_paths, data.inds_6d, strict=True):
            phdos_paths_6d[tuple(ind)] = phdos_path

        new = cls.from_files(phdos_paths_6d, data.gsr_bo_path, qha_model=data.qha_model, verbose=verbose)

        workdir = os.path.dirname(json_filepath)
        new.pickle_dump(workdir, basename=os.path.basename(json_filepath) + ".pickle")

        return new

    @classmethod
    def from_files(cls,
                   phdos_paths_6D,
                   gsr_bo_path,
                   qha_model: str = 'zsisa',
                   verbose: int = 0,
                   ) -> QHA_ZSISA:
        """
        Creates an instance of QHA from a 6D array of PHDOS.nc files and a Born-Oppenheimer (BO) GSR file.

        Args:
            phdos_paths_6D: A 6D list of paths to PHDOS.nc files.
                The PHDOS files must be provided according to the deformations defined in Table IV of the paper.
                For the 'zsisa' qha_model:
                A 6D array of PHDOS files is required, but the array does not need to be completely filled.
                Only the necessary deformations from Table IV need to exist.
                For cubic cases, a 1D or 3D array is also accepted for the TEC case.
                For uniaxial cases (hexagonal, trigonal, and tetragonal), a 2D or 3D
                array is also accepted for the TEC case.
            gsr_bo_path: Path to the GSR file for the Born-Oppenheimer structure,
                or the reference structure used to build deformations.
                This is needed to reconstruct strains from Eqs. (24) and (25) in the paper.
                and find the crystallographic symmetry of the structure.
            qha_model:
                Specifies the QHA model type. Options are:
                'zsisa': Standard ZSISA model.
                'v_zsisa': v-ZSISA model.
                'zsisa_slab': ZSISA model adapted for slab geometries.
        """
        # If the BO GSR file exists, read the structure and stress tensor.
        if not os.path.exists(gsr_bo_path):
            raise FileNotFoundError(f"Error: Born-Oppenheimer GSR file at {gsr_bo_path=} does not exist.")

        if gsr_bo_path.endswith("DDB"):
            with DdbFile.from_file(gsr_bo_path) as g:
                structure_bo = g.structure

        elif gsr_bo_path.endswith("GSR.nc"):
            with GsrFile.from_file(gsr_bo_path) as g:
                structure_bo = g.structure
        else:
            raise TypeError(f"Unknown file type: {type(gsr_bo_path)=}")

        spgrp = AbinitSpaceGroup.from_structure(structure_bo)
        spgrp_number = spgrp.spgid
        sym = 'unknown'
        #print(spgrp)

        # Find the crystallographic symmetry from the BO structure.
        if qha_model == 'zsisa':
            sym = spgnum_to_crystal_system(spgrp_number)

        phdos_paths_6D = np.array(phdos_paths_6D)

        current_shape = phdos_paths_6D.shape
        dims_to_add = 6 - len(current_shape)

        # If the array has fewer than 6 dimensions, reshape it to have 6 dimensions
        if dims_to_add > 0:
            new_shape = current_shape + (1,) * dims_to_add
            phdos_paths_6D = phdos_paths_6D.reshape(new_shape)

        dim = phdos_paths_6D.shape

        structures = []
        phdoses = []

        # Loop through each dimension of the 6D matrix (phdos_paths_6D) to read PHDOS data and
        # structures from the given file paths, handling missing files by appending None if the file does not exist.
        for dim1_idx, dim1_list in enumerate(phdos_paths_6D):
            dim1_doses = []
            dim1_structures = []
            dim1_energies = []
            for dim2_idx, dim2_list in enumerate(dim1_list):
                dim2_doses = []
                dim2_structures = []
                dim2_energies = []
                for dim3_idx, dim3_list in enumerate(dim2_list):
                    dim3_doses = []
                    dim3_structures = []
                    dim3_energies = []
                    for dim4_idx, dim4_list in enumerate(dim3_list):
                        dim4_doses = []
                        dim4_structures = []
                        dim4_energies = []
                        for dim5_idx, dim5_list in enumerate(dim4_list):
                            dim5_doses = []
                            dim5_structures = []
                            dim5_energies = []
                            for path_idx, path in enumerate(dim5_list):
                                if path is not None and os.path.exists(path):
                                #if path is not None:
                                #    if not os.path.exists(path):
                                #        raise FileNotFoundError(f"Cannot find {path=}")
                                    with PhdosFile(path) as p:
                                        dim5_doses.append(p.phdos)
                                        dim5_structures.append(p.structure)
                                else:
                                    # Handle the case when the file does not exist
                                    dim5_doses.append(None)
                                    dim5_structures.append(None)

                            dim4_doses.append(dim5_doses)
                            dim4_structures.append(dim5_structures)
                        dim3_doses.append(dim4_doses)
                        dim3_structures.append(dim4_structures)
                    dim2_doses.append(dim3_doses)
                    dim2_structures.append(dim3_structures)
                dim1_doses.append(dim2_doses)
                dim1_structures.append(dim2_structures)
            phdoses.append(dim1_doses)
            structures.append(dim1_structures)

        # If the structure is uniaxial and the input PHDOS data is 2D, expand it to 3D format.
        if list(dim) == [3, 3, 1, 1, 1, 1] and qha_model == 'zsisa':
            if not (sym in ('hexagonal', 'trigonal', 'tetragonal')):
                raise RuntimeError("Only uniaxial structures (e.g., hexagonal, trigonal, tetragonal) are allowed to have 2D PHDOS data.")
            new_shape = (3, 3, 3, 1, 1, 1)
            dim = [3, 3, 3, 1, 1, 1]
            structures2 = np.empty(new_shape, dtype=object)
            phdoses2 = np.empty(new_shape, dtype=object)

            # Copy data: Fill the new dim3 by copying dim2 data
            for i in range(3):  # Loop over dim1
                for j in range(3):  # Loop over dim2
                    structures2[i][i][j][0][0][0] = structures[i][j][0][0][0][0]
                    phdoses2[i][i][j][0][0][0] = phdoses[i][j][0][0][0][0]

            structures = structures2
            phdoses = phdoses2

        # If the structure is cubic and the input PHDOS data is 1D, expand it to a 3D format.
        if list(dim) == [3, 1, 1, 1, 1, 1] and qha_model == 'zsisa':
            if sym != 'cubic' :
                raise RuntimeError("Only cubic structure is allowed to have 1D PHDOS data.")

            new_shape = (3, 3, 3, 1, 1, 1)
            dim = [3, 3, 3, 1, 1, 1]
            structures2 = np.empty(new_shape, dtype=object)  # Use dtype=object for lists
            phdoses2 = np.empty(new_shape, dtype=object)

            # Copy data: Fill the new dim3 by copying dim1 data
            for i in range(3):  # Loop over dim1
                structures2[i][i][i][0][0][0] = structures[i][0][0][0][0][0]
                phdoses2[i][i][i][0][0][0] = phdoses[i][0][0][0][0][0]
            structures = structures2
            phdoses = phdoses2

        return cls(structures, phdoses, dim, structure_bo, sym, verbose=verbose)

    def __init__(self,
                 structures: list[Structure],
                 phdoses,
                 dim,
                 structure_bo: Structure,
                 sym: str,
                 verbose: int = 0):
        """
        Args:
            structures: List of structures at different strains.
            phdoses: List of phonon density of states (PHDOS) objects computed at different strains.
            dim: Shape of the 6D dataset.
            structure_bo: Born-Oppenheimer reference structure.
            sym: Crystallographic symmetry of the reference structure.
        """
        self.structures = structures
        self.sym = sym
        #phdoses = np.reshape(phdoses, (3, 3, 3, 3, 3, 3))
        self.phdoses = phdoses
        self.dim = dim
        #print(f"{self.phdoses.shape=}, {self.dim=}")
        self.verbose = verbose

        self.dtol_tolerance = 1e-8

        def extract_attribute(structures, attribute_func) -> np.ndarray:
            return np.array([[[[[[attribute_func(s) if s is not None else None for s in col] for col in row] for row in d1]
                            for d1 in d2] for d2 in d3] for d3 in structures])

        # Extract lattice volume for each structure
        self.volumes = extract_attribute(structures, lambda s: s.volume)

        # Extract lattice parameters a, b, and c
        self.lattice_a = extract_attribute(structures, lambda s: s.lattice.abc[0])
        self.lattice_b = extract_attribute(structures, lambda s: s.lattice.abc[1])
        self.lattice_c = extract_attribute(structures, lambda s: s.lattice.abc[2])

        # Extract lattice angle
        self.alpha = extract_attribute(structures, lambda s: s.lattice.angles[0])
        self.beta = extract_attribute(structures, lambda s: s.lattice.angles[1])
        self.gamma = extract_attribute(structures, lambda s: s.lattice.angles[2])

        # Lattice vectors from the lattice matrix
        # Eq (17) from the paper:
        # R1 = (ax, ay, az) corresponds to (R1x, R1y, R1z)
        # R2 = (bx, by, bz) corresponds to (R2x, R2y, R2z)
        # R3 = (cx, cy, cz) corresponds to (R3x, R3y, R3z)
        self.ax = extract_attribute(structures, lambda s: s.lattice.matrix[0, 0])
        self.ay = extract_attribute(structures, lambda s: s.lattice.matrix[0, 1])
        self.az = extract_attribute(structures, lambda s: s.lattice.matrix[0, 2])
        self.bx = extract_attribute(structures, lambda s: s.lattice.matrix[1, 0])
        self.by = extract_attribute(structures, lambda s: s.lattice.matrix[1, 1])
        self.bz = extract_attribute(structures, lambda s: s.lattice.matrix[1, 2])
        self.cx = extract_attribute(structures, lambda s: s.lattice.matrix[2, 0])
        self.cy = extract_attribute(structures, lambda s: s.lattice.matrix[2, 1])
        self.cz = extract_attribute(structures, lambda s: s.lattice.matrix[2, 2])
        self.ave_x = np.zeros_like(self.volumes)
        self.ave_y = np.zeros_like(self.volumes)
        self.ave_z = np.zeros_like(self.volumes)

        # Create a mask where self.ax is not None
        mask = self.ax != None

        # Compute averages. Eq (43) from the paper:
        # ave_x corresponds to A_x in the paper
        # ave_y corresponds to B_y in the paper
        # ave_z corresponds to C_z in the paper
        self.ave_x[mask] = (abs(self.ax[mask]) + abs(self.bx[mask]) + abs(self.cx[mask]))
        self.ave_y[mask] = (abs(self.ay[mask]) + abs(self.by[mask]) + abs(self.cy[mask]))
        self.ave_z[mask] = (abs(self.az[mask]) + abs(self.bz[mask]) + abs(self.cz[mask]))

        # Store structure parameters for BO
        self.ax_bo = structure_bo.lattice.matrix[0,0]
        self.ay_bo = structure_bo.lattice.matrix[0,1]
        self.az_bo = structure_bo.lattice.matrix[0,2]
        self.bx_bo = structure_bo.lattice.matrix[1,0]
        self.by_bo = structure_bo.lattice.matrix[1,1]
        self.bz_bo = structure_bo.lattice.matrix[1,2]
        self.cx_bo = structure_bo.lattice.matrix[2,0]
        self.cy_bo = structure_bo.lattice.matrix[2,1]
        self.cz_bo = structure_bo.lattice.matrix[2,2]

        # Eq (43)
        self.ave_x_bo = (abs(self.ax_bo)+abs(self.bx_bo)+abs(self.cx_bo))
        self.ave_y_bo = (abs(self.ay_bo)+abs(self.by_bo)+abs(self.cy_bo))
        self.ave_z_bo = (abs(self.az_bo)+abs(self.bz_bo)+abs(self.cz_bo))

    def _set_structure_stress_guess(self, structure_guess: Structure, stress_guess, energy_guess) -> None:
        """
        Set structure parameters and stress for the initial guess.
        """
        if stress_guess.shape != (3, 3):
            raise ValueError(f"{stress_guess.shape=} != (3, 3)")

        self.structure_guess = structure_guess
        self.volume_guess = structure_guess.volume
        self.lattice_a_guess = structure_guess.lattice.abc[0]
        self.lattice_b_guess = structure_guess.lattice.abc[1]
        self.lattice_c_guess = structure_guess.lattice.abc[2]
        self.matrix_guess = structure_guess.lattice.matrix
        self.stress_guess = stress_guess
        self.energy_guess = energy_guess
        self.frac_coords_guess = structure_guess.frac_coords
        self.angles_guess = structure_guess.lattice.angles
        self.ax_guess = structure_guess.lattice.matrix[0,0]
        self.ay_guess = structure_guess.lattice.matrix[0,1]
        self.az_guess = structure_guess.lattice.matrix[0,2]
        self.bx_guess = structure_guess.lattice.matrix[1,0]
        self.by_guess = structure_guess.lattice.matrix[1,1]
        self.bz_guess = structure_guess.lattice.matrix[1,2]
        self.cx_guess = structure_guess.lattice.matrix[2,0]
        self.cy_guess = structure_guess.lattice.matrix[2,1]
        self.cz_guess = structure_guess.lattice.matrix[2,2]

        # Eq (43)
        self.ave_x_guess = (abs(self.ax_guess)+abs(self.bx_guess)+abs(self.cx_guess))
        self.ave_y_guess = (abs(self.ay_guess)+abs(self.by_guess)+abs(self.cy_guess))
        self.ave_z_guess = (abs(self.az_guess)+abs(self.bz_guess)+abs(self.cz_guess))

    def print_mat_entries(self, what: str, file=sys.stdout) -> None:
        mat6d = getattr(self, what, None)
        if mat6d is None:
            raise ValueError("Invalid value for {what=}")

        for index, value in np.ndenumerate(mat6d):
            if value is None: continue
            print(f"Index: {index}, {what}: {value}", file=file)

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        """
        String representation with verbosity level `verbose`.
        """
        from io import StringIO
        file = StringIO()
        def _p(*args, **kwargs):
            print(*args, file=file, **kwargs)

        _p("Structure:\n", self.structures)
        _p(f"sym={self.sym}, dym={self.dim}")

        what_list = ("ax", "ay", "az")
        for what in what_list:
            _p(f"Attribute: {what}")
            self.print_mat_entries(what, file=file)
            _p("")

        return file.getvalue()

    def stress_v_ZSISA(self, temp: float, pressure: float) -> tuple:

        e, S = self.get_vib_free_energies(temp)

        v = self.volume_guess
        dv = self.volumes[0,0,0,0,0,0] - self.volumes[1,0,0,0,0,0]

        if self.dim[0] == 3:
            v0 = self.volumes[1,0,0,0,0,0]
            dF_dV = (e[0,0,0,0,0,0]-e[2,0,0,0,0,0])/(2*dv)
            d2F_dV2 = (e[0,0,0,0,0,0]-2*e[1,0,0,0,0,0]+e[2,0,0,0,0,0])/(dv)**2
            dfdv = dF_dV + (v-v0)*d2F_dV2
            gibbs = self.energy_guess + e[1,0,0,0,0,0] + pressure*v/ abu.eVA3_HaBohr3 +\
                    (v-v0)*dF_dV + 0.5*(v-v0)**2*d2F_dV2

        elif self.dim[0] == 5:
            v0 = self.volumes[2,0,0,0,0,0]
            dF_dV = (-e[0,0,0,0,0,0]+ 8*e[1,0,0,0,0,0]-8*e[3,0,0,0,0,0]+e[4,0,0,0,0,0])/(12*dv)
            d2F_dV2 = (-e[0,0,0,0,0,0]+16*e[1,0,0,0,0,0]-30*e[2,0,0,0,0,0]+16*e[3,0,0,0,0,0]-e[4,0,0,0,0,0])/(12*dv**2)
            d3F_dV3 = (e[0,0,0,0,0,0]-2*e[1,0,0,0,0,0]+2*e[3,0,0,0,0,0]-e[4,0,0,0,0,0])/(2*dv**3)
            d4F_dV4 = (e[0,0,0,0,0,0]-4*e[1,0,0,0,0,0]+6*e[2,0,0,0,0,0]-4*e[3,0,0,0,0,0]+e[4,0,0,0,0,0])/(dv**4)
            dfdv = dF_dV + (v-v0)*d2F_dV2+ 0.5*(v-v0)**2*d3F_dV3+ 1/6.0*(v-v0)**3*d4F_dV4
            gibbs = self.energy_guess + e[1,0,0,0,0,0] + pressure*v/ abu.eVA3_HaBohr3+\
                    (v-v0)*dF_dV + 0.5*(v-v0)**2*d2F_dV2 + 1/6.0*(v-v0)**3*d3F_dV3+ 1/24.0*(v-v0)**4*d4F_dV4

        else:
            raise ValueError(f"Unsupported value of {self.dim[0]=}")

        stress_a = -dfdv * abu.eVA3_HaBohr3
        stress = np.zeros(6)
        stress[0] = stress_a - pressure
        stress[1] = stress_a - pressure
        stress[2] = stress_a - pressure

        dtol = np.zeros(6)
        dtol[0] = abs(stress[0] - self.stress_guess[0,0])
        dtol[1] = abs(stress[1] - self.stress_guess[1,1])
        dtol[2] = abs(stress[2] - self.stress_guess[2,2])

        return dtol, gibbs, stress

    def stress_ZSISA_1DOF(self, temp: float, pressure: float) -> tuple:
        """
        Compute thermal stress and thermal expansion for cubic structures,
        given temperature and pressure. If self-consistent convergence is achieved
        and BO-elastic constants exist, the thermal expansion is computed.
        """
        # Get vibrational free energy and entropy at a specific temperature
        # e = Vibrational free energy (F_vib)
        # S = Entropy (S)
        e, S = self.get_vib_free_energies(temp)

        XBO = self.ave_x_bo
        X0 = self.ave_x[0,0,0,0,0,0] # Slightly deformed structure: reference structure - exx0 -eyy0 -ezz0
        X1 = self.ave_x[1,1,1,0,0,0] # Reference structure

        # Compute strain-related quantities. Eq (50)
        dexx = (X0 - X1) / XBO  # Strain step
        exx0 = X1 / XBO - 1  # Reference strain
        #print(f"{dexx=}, {exx0=}")

        # Compute first and second derivatives of free energy w.r.t. strain
        dF_dX = (e[0,0,0,0,0,0]-e[2,2,2,0,0,0])/(2*dexx)
        d2F_dX2 = (e[0,0,0,0,0,0]-2*e[1,1,1,0,0,0]+e[2,2,2,0,0,0])/(dexx)**2

        # Compute first and second derivatives of entropy w.r.t. strain
        dS_dX = (S[0,0,0,0,0,0]-S[2,2,2,0,0,0])/(2*dexx)
        d2S_dX2 = (S[0,0,0,0,0,0]-2*S[1,1,1,0,0,0]+S[2,2,2,0,0,0])/(dexx)**2

        x = self.ave_x_guess
        v = self.volume_guess
        # compute BO strain at guess, Eq (50)
        exx_n = x/XBO-1

        # Free energy and entropy derivative at guess structure
        dfdx = dF_dX + (exx_n-exx0)*d2F_dX2
        dsdx = dS_dX + (exx_n-exx0)*d2S_dX2

        # Gibbs free energy
        gibbs = self.energy_guess + e[1,1,1,0,0,0] + (exx_n-exx0)*dF_dX + 0.5*(exx_n-exx0)**2*d2F_dX2 + pressure*v / abu.eVA3_HaBohr3

        # Compute thermal stress. Eq (51)
        stress_xx = -dfdx/v*(exx_n+1)/3.0 * abu.eVA3_HaBohr3
        if self.verbose:
            print("x/XBO: ", x/XBO)
            print(f"{stress_xx=}")

        # Apply external pressure
        stress = np.zeros(6)
        stress[0] = stress_xx - pressure
        stress[1] = stress_xx - pressure
        stress[2] = stress_xx - pressure

        # Calculate the absolute difference (tolerance) between the current stress and the guessed stress
        # This is used to check the convergence of the stress values during the optimization process
        dtol = np.zeros(6)
        dtol[0] = abs(stress[0]-self.stress_guess[0,0])
        dtol[1] = abs(stress[1]-self.stress_guess[1,1])
        dtol[2] = abs(stress[2]-self.stress_guess[2,2])

        therm = None
        # Check if the stress has converged (all tolerances below self.dtol_tolerance)
        if all(dtol[i] < self.dtol_tolerance for i in range(6)): # Check convergence
            if self.bo_elastic_voigt is not None:
                # Read elastic constants from the output of DFPT obtained by abiopen.py (ELASTIC_RELAXED)
                matrix_elastic = self.bo_elastic_voigt
                # Convert elastic constants to second derivative of BO energy (Eq. 39)
                matrix_elastic = matrix_elastic*v/abu.eVA3_GPa
                # Initialize second derivative matrix M using free energy second derivatives.
                M = np.array([[ d2F_dX2/3*(exx0+1)**2  , 0          , 0         ,0,0,0],
                              [    0       , d2F_dX2/3*(exx0+1)**2  , 0         ,0,0,0],
                              [    0       , 0          , d2F_dX2/3*(exx0+1)**2 ,0,0,0],
                              [0,0,0,0,0,0],
                              [0,0,0,0,0,0],
                              [0,0,0,0,0,0]])
                # Contribution of pressure to the second derivative matrix
                P = np.array([[0.0 ,pressure*v,pressure*v  ,0,0,0],
                              [pressure*v ,0.0,pressure*v  ,0,0,0],
                              [pressure*v ,pressure*v,0.0  ,0,0,0],
                              [0,0,0,0,0,0],
                              [0,0,0,0,0,0],
                              [0,0,0,0,0,0]])
                # Compute the final second derivative matrix M, including elastic and pressure contributions
                M = M + matrix_elastic + P / abu.eVA3_HaBohr3
                # Scale the first derivative of entropy
                S1 = dsdx*(exx_n+1)/3
                dSde = np.array([S1,S1,S1,0,0,0])
                # Compute thermal expansion using Eq. (37)
                dstrain_dt = np.linalg.inv(M) @ dSde
                # Scale the thermal expansion
                therm = [dstrain_dt[0]*(exx_n+1) , dstrain_dt[1]*(exx_n+1),dstrain_dt[2]*(exx_n+1),0,0,0]
                if self.verbose: print(f"{therm=}")

        return dtol, gibbs, stress, therm

    def stress_ZSISA_2DOF(self, temp: float, pressure: float) -> tuple:

        # Get vibrational free energy and entropy at a specific temperature
        # e = Vibrational free energy (F_vib)
        # S = Entropy (S)
        e, S = self.get_vib_free_energies(temp)

        # A_x and C_z from deformed structures based on table IV and Eq (52)
        X0 = self.ave_x[0,0,1,0,0,0] # reference structure - exx0 -eyy0
        Z0 = self.ave_z[1,1,0,0,0,0] # reference structure - ezz0

        # A_x and C_z from reference structure
        X1 = self.ave_x[1,1,1,0,0,0]
        Z1 = self.ave_z[1,1,1,0,0,0]

        XBO = self.ave_x_bo
        ZBO = self.ave_z_bo

        # Compute strain-related quantities. Eq (53)
        # Strain steps
        dexx = (X0-X1)/XBO
        dezz = (Z0-Z1)/ZBO
        # Reference strains
        exx0 = X1/XBO-1
        ezz0 = Z1/ZBO-1

        # Compute first and second derivatives of free energy w.r.t. strain
        dF_dX = (e[0,0,1,0,0,0]-e[2,2,1,0,0,0])/(2*dexx)
        dF_dZ = (e[1,1,0,0,0,0]-e[1,1,2,0,0,0])/(2*dezz)

        d2F_dX2 = (e[0,0,1,0,0,0]-2*e[1,1,1,0,0,0]+e[2,2,1,0,0,0])/(dexx)**2
        d2F_dZ2 = (e[1,1,0,0,0,0]-2*e[1,1,1,0,0,0]+e[1,1,2,0,0,0])/(dezz)**2
        d2F_dXdZ = (e[1,1,1,0,0,0] - e[0,0,1,0,0,0] - e[1,1,0,0,0,0] + e[0,0,0,0,0,0]) / (dexx *dezz)

        # Compute first and second derivatives of entropy w.r.t. strain
        dS_dX = (S[0,0,1,0,0,0]-S[2,2,1,0,0,0])/(2*dexx)
        dS_dZ = (S[1,1,0,0,0,0]-S[1,1,2,0,0,0])/(2*dezz)
        d2S_dX2 = (S[0,0,1,0,0,0]-2*S[1,1,1,0,0,0]+S[2,2,1,0,0,0])/(dexx)**2
        d2S_dZ2 = (S[1,1,0,0,0,0]-2*S[1,1,1,0,0,0]+S[1,1,2,0,0,0])/(dezz)**2
        d2S_dXdZ = (S[1,1,1,0,0,0] - S[0,0,1,0,0,0] - S[1,1,0,0,0,0] + S[0,0,0,0,0,0]) / (dexx *dezz)

        x = self.ave_x_guess
        z = self.ave_z_guess
        v = self.volume_guess
        # compute BO strain at guess, Eq (53)
        exx_n = x/XBO-1
        ezz_n = z/ZBO-1

        # Free energy and entropy derivative at guess structure
        dfdx = dF_dX + (exx_n-exx0)*d2F_dX2+(ezz_n-ezz0)*d2F_dXdZ
        dfdz = dF_dZ + (ezz_n-ezz0)*d2F_dZ2+(exx_n-exx0)*d2F_dXdZ

        dsdx = dS_dX + (exx_n-exx0)*d2S_dX2+(ezz_n-ezz0)*d2S_dXdZ
        dsdz = dS_dZ + (ezz_n-ezz0)*d2S_dZ2+(exx_n-exx0)*d2S_dXdZ

        # Gibbs free energy
        gibbs = self.energy_guess + e[1,1,1,0,0,0] + pressure*v/ abu.eVA3_HaBohr3 + \
                  (exx_n-exx0)*dF_dX + 0.5*(exx_n-exx0)**2*d2F_dX2 +\
                  (ezz_n-ezz0)*dF_dZ + 0.5*(ezz_n-ezz0)**2*d2F_dZ2 +\
                  (exx_n-exx0)*(ezz_n-ezz0)*d2F_dXdZ

        # Compute thermal stresses. Eq (54)
        stress_xx = -dfdx/v*(exx_n+1)*0.5 * abu.eVA3_HaBohr3
        stress_zz = -dfdz/v*(ezz_n+1)     * abu.eVA3_HaBohr3

        # Apply external pressure
        stress = np.zeros(6)
        stress[0] = stress_xx -pressure
        stress[1] = stress_xx -pressure
        stress[2] = stress_zz -pressure
        if self.verbose:
            print("x/XBO, z/ZBO: ", x/XBO, z/ZBO)
            print(f"{pressure=}")
            print(f"{stress[0]=}, {stress[2]=}")

        # Calculate the absolute difference (tolerance) between the current stress and the guessed stress
        # This is used to check the convergence of the stress values during the optimization process
        dtol = np.zeros(6)
        dtol[0] = abs(stress[0]-self.stress_guess[0,0])
        dtol[1] = abs(stress[1]-self.stress_guess[1,1])
        dtol[2] = abs(stress[2]-self.stress_guess[2,2])

        therm = None
        # Check if the stress has converged (all tolerances below self.dtol_tolerance)
        if all(dtol[i] < self.dtol_tolerance for i in range(6)):
            if self.bo_elastic_voigt is not None:
                # Read elastic constants from the output of DFPT obtained by abiopen.py (ELASTIC_RELAXED)
                matrix_elastic = self.bo_elastic_voigt
                # Convert elastic constants to second derivative of BO energy (Eq. 39)
                matrix_elastic = matrix_elastic*v/abu.eVA3_GPa
                # Initialize second derivative matrix M using free energy second derivatives.
                M = np.array([[ d2F_dX2/2 *(exx0+1)**2       , 0                            , d2F_dXdZ/2*(exx0+1)*(ezz0+1) ,0,0,0],
                              [    0                         , d2F_dX2/2 *(exx0+1)**2       , d2F_dXdZ/2*(exx0+1)*(ezz0+1) ,0,0,0],
                              [ d2F_dXdZ/2*(exx0+1)*(ezz0+1) , d2F_dXdZ/2*(exx0+1)*(ezz0+1) , d2F_dZ2*(ezz0+1)**2          ,0,0,0],
                              [0,0,0,0,0,0],
                              [0,0,0,0,0,0],
                              [0,0,0,0,0,0]])
                # Contribution of pressure to the second derivative matrix
                P = np.array([[0.0 ,pressure*v,pressure*v  ,0,0,0],
                              [pressure*v ,0.0,pressure*v  ,0,0,0],
                              [pressure*v ,pressure*v,0.0  ,0,0,0],
                              [0,0,0,0,0,0],
                              [0,0,0,0,0,0],
                              [0,0,0,0,0,0]])
                # Compute the final second derivative matrix M, including elastic and pressure contributions
                M = M + matrix_elastic + P / abu.eVA3_HaBohr3
                # Scale the first derivative of entropy
                S1 = dsdx*(exx_n+1)*0.5
                S3 = dsdz*(ezz_n+1)
                dSde = np.array([S1,S1,S3,0,0,0])
                # Compute thermal expansion using Eq. (37)
                dstrain_dt = np.linalg.inv(M) @ dSde
                # Scale the thermal expansion
                therm = [dstrain_dt[0]*(exx_n+1), dstrain_dt[1]*(exx_n+1),dstrain_dt[2]*(ezz_n+1),0,0,0]
                if self.verbose: print(f"{therm=}")

        return dtol, gibbs, stress, therm

    def stress_ZSISA_3DOF(self, temp: float, pressure: float, mode: str) -> tuple:

        # Get vibrational free energy and entropy at a specific temperature
        # e = Vibrational free energy (F_vib)
        # S = Entropy (S)
        e, S = self.get_vib_free_energies(temp)

        # A_x, B_y and C_z from deformed structures based on table IV and Eq (42)
        X0 = self.ave_x[0,1,1,0,0,0] # reference structure - exx0
        Y0 = self.ave_y[1,0,1,0,0,0] # reference structure - eyy0
        Z0 = self.ave_z[1,1,0,0,0,0] # reference structure - ezz0

        # A_x, B_y and C_z from reference structure
        X1 = self.ave_x[1,1,1,0,0,0]
        Y1 = self.ave_y[1,1,1,0,0,0]
        Z1 = self.ave_z[1,1,1,0,0,0]

        XBO = self.ave_x_bo
        YBO = self.ave_y_bo
        ZBO = self.ave_z_bo

        # Compute strain-related quantities. Eq (44)
        # Strain steps and reference strains
        dexx = (X0-X1)/XBO
        exx0 = X1/XBO-1

        # Apply symmetries for cubic and uniaxial cases.
        if self.sym == "cubic":
            deyy = dezz = dexx
            eyy0 = ezz0 = exx0

        elif self.sym in ("trigonal", "hexagonal", "tetragonal"):
            deyy = dexx
            eyy0 = exx0
            dezz = (Z0-Z1)/ZBO
            ezz0 = Z1/ZBO-1
        else:
            dexx = (X0-X1)/XBO
            dezz = (Z0-Z1)/ZBO
            deyy = (Y0-Y1)/YBO

            exx0 = X1/XBO-1
            eyy0 = Y1/YBO-1
            ezz0 = Z1/ZBO-1

        # Compute first and second derivatives of free energy w.r.t. strains (symmetries are applied)
        dF_dX = (e[0,1,1,0,0,0]-e[2,1,1,0,0,0])/(2*dexx)
        d2F_dX2 = (e[0,1,1,0,0,0]-2*e[1,1,1,0,0,0]+e[2,1,1,0,0,0])/(dexx)**2
        dS_dX = (S[0,1,1,0,0,0]-S[2,1,1,0,0,0])/(2*dexx)
        d2S_dX2 = (S[0,1,1,0,0,0]-2*S[1,1,1,0,0,0]+S[2,1,1,0,0,0])/(dexx)**2

        if self.sym in ("cubic", "trigonal", "hexagonal", "tetragonal"):
            d2F_dXdY = (e[1,1,1,0,0,0] - e[0,1,1,0,0,0] - e[0,1,1,0,0,0] + e[0,0,1,0,0,0]) / (dexx *deyy)
            d2S_dXdY = (S[1,1,1,0,0,0] - S[0,1,1,0,0,0] - S[0,1,1,0,0,0] + S[0,0,1,0,0,0]) / (dexx *deyy)

        if self.sym in ("trigonal", "hexagonal", "tetragonal", "orthorhombic"):
            dF_dZ = (e[1,1,0,0,0,0]-e[1,1,2,0,0,0])/(2*dezz)
            d2F_dZ2 = (e[1,1,0,0,0,0]-2*e[1,1,1,0,0,0]+e[1,1,2,0,0,0])/(dezz)**2
            d2F_dXdZ = (e[1,1,1,0,0,0] - e[0,1,1,0,0,0] - e[1,1,0,0,0,0] + e[0,1,0,0,0,0]) / (dexx *dezz)
            dS_dZ = (S[1,1,0,0,0,0]-S[1,1,2,0,0,0])/(2*dezz)
            d2S_dZ2 = (S[1,1,0,0,0,0]-2*S[1,1,1,0,0,0]+S[1,1,2,0,0,0])/(dezz)**2
            d2S_dXdZ = (S[1,1,1,0,0,0] - S[0,1,1,0,0,0] - S[1,1,0,0,0,0] + S[0,1,0,0,0,0]) / (dexx *dezz)

        if self.sym == "cubic":
            dF_dY = dF_dX
            dF_dZ = dF_dX
            d2F_dY2 = d2F_dX2
            d2F_dZ2 = d2F_dX2
            d2F_dXdZ = d2F_dXdY
            d2F_dYdZ = d2F_dXdY
            dS_dY = dS_dX
            dS_dZ = dS_dX
            d2S_dY2 = d2S_dX2
            d2S_dZ2 = d2S_dX2
            d2S_dXdZ = d2S_dXdY
            d2S_dYdZ = d2S_dXdY

        elif self.sym in ("trigonal", "hexagonal", "tetragonal"):
            dF_dY = dF_dX
            d2F_dY2 = d2F_dX2
            d2F_dYdZ = d2F_dXdZ
            dS_dY = dS_dX
            d2S_dY2 = d2S_dX2
            d2S_dYdZ = d2S_dXdZ

        elif self.sym == "orthorhombic":
            dF_dY = (e[1,0,1,0,0,0]-e[1,2,1,0,0,0])/(2*deyy)
            d2F_dY2 = (e[1,0,1,0,0,0]-2*e[1,1,1,0,0,0]+e[1,2,1,0,0,0])/(deyy)**2
            d2F_dXdY = (e[1,1,1,0,0,0] - e[0,1,1,0,0,0] - e[1,0,1,0,0,0] + e[0,0,1,0,0,0]) / (dexx *deyy)
            d2F_dYdZ = (e[1,1,1,0,0,0] - e[1,0,1,0,0,0] - e[1,1,0,0,0,0] + e[1,0,0,0,0,0]) / (deyy *dezz)

            dS_dY = (S[1,0,1,0,0,0]-S[1,2,1,0,0,0])/(2*deyy)
            d2S_dY2 = (S[1,0,1,0,0,0]-2*S[1,1,1,0,0,0]+S[1,2,1,0,0,0])/(deyy)**2
            d2S_dXdY = (S[1,1,1,0,0,0] - S[0,1,1,0,0,0] - S[1,0,1,0,0,0] + S[0,0,1,0,0,0]) / (dexx *deyy)
            d2S_dYdZ = (S[1,1,1,0,0,0] - S[1,0,1,0,0,0] - S[1,1,0,0,0,0] + S[1,0,0,0,0,0]) / (deyy *dezz)

        d2F_dXY2 = 0.0
        d2F_dXdYZ = 0.0
        d2F_dYZ2 = 0.0
        d2F_dXZ2 = 0.0

        # If elastic constants are requested, compute shear strain steps and their related second derivatives
        # by fitting a quadratic curve. Section G in APPENDIX.
        # The necessary derivatives are related to the symmetries. See table II.
        if mode == 'ECs':
            deyz= (self.ave_y[1,1,1,0,0,0] - self.ave_y[1,1,1,1,0,0])/ZBO
            #deyz= (self.ave_x[1,1,1,0,0,0] - self.ave_x[1,1,1,1,0,0])/ZBO
            Fyz = [e[1,1,1,2,0,0],e[1,1,1,1,0,0],e[1,1,1,1,0,0],e[1,1,1,2,0,0]]
            Dyz = [2*deyz,deyz,-deyz,-2*deyz]
            param = np.polyfit(Dyz,Fyz,2)
            d2F_dYZ2 = 2*param[0]

            if self.sym == "trigonal":
                d2F_dXdYZ = (e[1,1,1,1,0,0] - e[0,1,1,1,0,0] - e[1,1,1,0,0,0] + e[0,1,1,0,0,0]) / (dexx *deyy)

            if self.sym in ("tetragonal", "orthorhombic"):
                dexy = (self.ave_x[1,1,1,0,0,0] - self.ave_x[1,1,1,0,0,1])/YBO
                Fxy = [e[1,1,1,0,0,2],e[1,1,1,0,0,1],e[1,1,1,0,0,1],e[1,1,1,0,0,2]]
                Dxy = [2*dexy,dexy,-dexy,-2*dexy]
                param = np.polyfit(Dxy,Fxy,2)
                d2F_dXY2 = 2*param[0]

            if self.sym == "orthorhombic":
                dexz = (self.ave_x[1,1,1,0,0,0] - self.ave_x[1,1,1,0,1,0])/ZBO
                Fxz = [e[1,1,1,0,2,0],e[1,1,1,0,1,0],e[1,1,1,0,1,0],e[1,1,1,0,2,0]]
                Dxz = [2*dexz,dexz,-dexz,-2*dexz]
                param = np.polyfit(Dxz,Fxz,2)
                d2F_dXZ2 = 2*param[0]

        x = self.ave_x_guess
        y = self.ave_y_guess
        z = self.ave_z_guess
        v = self.volume_guess

        # Compute BO strain at guess, Eq (44)
        exx_n = x/XBO-1
        eyy_n = y/YBO-1
        ezz_n = z/ZBO-1

        # Free energy and entropy derivative at guess structure
        dfdx = dF_dX + (exx_n-exx0)*d2F_dX2+(eyy_n-eyy0)*d2F_dXdY+(ezz_n-ezz0)*d2F_dXdZ
        dfdy = dF_dY + (eyy_n-eyy0)*d2F_dY2+(exx_n-exx0)*d2F_dXdY+(ezz_n-ezz0)*d2F_dYdZ
        dfdz = dF_dZ + (ezz_n-ezz0)*d2F_dZ2+(exx_n-exx0)*d2F_dXdZ+(eyy_n-eyy0)*d2F_dYdZ

        dsdx = dS_dX + (exx_n-exx0)*d2S_dX2+(eyy_n-eyy0)*d2S_dXdY+(ezz_n-ezz0)*d2S_dXdZ
        dsdy = dS_dY + (eyy_n-eyy0)*d2S_dY2+(exx_n-exx0)*d2S_dXdY+(ezz_n-ezz0)*d2S_dYdZ
        dsdz = dS_dZ + (ezz_n-ezz0)*d2S_dZ2+(exx_n-exx0)*d2S_dXdZ+(eyy_n-eyy0)*d2S_dYdZ

        # Gibbs free energy
        gibbs = self.energy_guess + e[1,1,1,0,0,0] + pressure*v / abu.eVA3_HaBohr3 + \
                  (exx_n-exx0)*dF_dX + 0.5*(exx_n-exx0)**2*d2F_dX2 +\
                  (eyy_n-eyy0)*dF_dY + 0.5*(eyy_n-eyy0)**2*d2F_dY2 +\
                  (ezz_n-ezz0)*dF_dZ + 0.5*(ezz_n-ezz0)**2*d2F_dZ2 +\
                  (exx_n-exx0)*(ezz_n-ezz0)*d2F_dXdZ +\
                  (exx_n-exx0)*(ezz_n-eyy0)*d2F_dXdY +\
                  (eyy_n-eyy0)*(ezz_n-ezz0)*d2F_dYdZ

        # Compute thermal stresses. Eq (45)
        stress_xx = -dfdx/v*(exx_n+1) * abu.eVA3_HaBohr3
        stress_yy = -dfdy/v*(eyy_n+1) * abu.eVA3_HaBohr3
        stress_zz = -dfdz/v*(ezz_n+1) * abu.eVA3_HaBohr3

        stress = np.zeros(6)
        stress[0] = stress_xx - pressure
        stress[1] = stress_yy - pressure
        stress[2] = stress_zz - pressure
        if self.verbose:
            print("x/XBO, y/YBO, z/ZBO: ", x/XBO, y/YBO, z/ZBO)
            print(f"{stress[0]=}, {stress[1]=}, {stress[2]=}")

        dtol = np.zeros(6)
        dtol[0] = abs(stress[0]-self.stress_guess[0,0])
        dtol[1] = abs(stress[1]-self.stress_guess[1,1])
        dtol[2] = abs(stress[2]-self.stress_guess[2,2])

        # Check if the stress has converged (all tolerances below self.dtol_tolerance)
        therm = None
        elastic = None
        if all(dtol[i] < self.dtol_tolerance for i in range(6)):
            if self.bo_elastic_voigt is not None:
                # Read elastic constants from the output of DFPT obtained by abiopen.py (ELASTIC_RELAXED)
                matrix_elastic = self.bo_elastic_voigt
                # Convert elastic constants to second derivative of BO energy (Eq. 39)
                matrix_elastic = matrix_elastic*v/abu.eVA3_GPa
                scale_xx = (exx0+1)*(exx0+1)
                scale_yy = (eyy0+1)*(eyy0+1)
                scale_zz = (ezz0+1)*(ezz0+1)
                scale_xz = (exx0+1)*(ezz0+1)
                scale_yz = (eyy0+1)*(ezz0+1)
                scale_xy = (exx0+1)*(eyy0+1)
                # Initialize the second derivative matrix M using free energy second derivatives.
                # For 'TEC' mode or certain symmetries, the terms d2F_dXdYZ, d2F_dYZ2, d2F_dXZ2, d2F_dXY2 are zero.
                # These terms don't contribute to thermal expansion calculations because, due to the structure symmetries,
                # the first derivatives of free energy and entropy are zero for shear strains.
                M = np.array([[ d2F_dX2*scale_xx, d2F_dXdY*scale_xy ,d2F_dXdZ*scale_xz ,d2F_dXdYZ*scale_xz  ,0.0               ,0.0              ],
                              [d2F_dXdY*scale_xy, d2F_dY2 *scale_yy ,d2F_dYdZ*scale_yz ,0.0                 ,0.0               ,0.0              ],
                              [d2F_dXdZ*scale_xz, d2F_dYdZ*scale_yz ,d2F_dZ2 *scale_zz ,0.0                 ,0.0               ,0.0              ],
                              [d2F_dXdYZ*scale_xz,0.0               ,0.0               ,d2F_dYZ2*scale_zz   ,0.0               ,0.0              ],
                              [0.0              , 0.0               ,0.0               ,0.0                 ,d2F_dXZ2*scale_zz ,0.0              ],
                              [0.0              , 0.0               ,0.0               ,0.0                 ,0.0               ,d2F_dXY2*scale_yy]])
                # Contribution of pressure to the second derivative matrix
                P = np.array([[0.0 ,pressure*v,pressure*v  ,0,0,0],
                              [pressure*v ,0.0,pressure*v  ,0,0,0],
                              [pressure*v ,pressure*v,0.0  ,0,0,0],
                              [0,0,0,0,0,0],
                              [0,0,0,0,0,0],
                              [0,0,0,0,0,0]])
                # Compute the final second derivative matrix M, including elastic and pressure contributions
                M = M + matrix_elastic + P / abu.eVA3_HaBohr3
                # Scale the first derivative of entropy
                S1 = dsdx*(exx_n+1)
                S2 = dsdy*(eyy_n+1)
                S3 = dsdz*(ezz_n+1)
                dSde = np.array([S1,S2,S3,0,0,0])
                # Compute thermal expansion using Eq. (37)
                dstrain_dt = np.linalg.inv(M) @ dSde
                # Scale the thermal expansion
                therm = [dstrain_dt[0]*(exx_n+1) , dstrain_dt[1]*(eyy_n+1),dstrain_dt[2]*(ezz_n+1),0,0,0]
                #print(f"{therm=}")
                elastic = M / v * abu.eVA3_GPa

        return dtol, gibbs, stress, therm, elastic

    def stress_ZSISA_monoclinic(self, temp: float, pressure: float, mode: str) -> tuple:

        # Get vibrational free energy and entropy at a specific temperature
        # e = Vibrational free energy (F_vib)
        # S = Entropy (S)
        e, S = self.get_vib_free_energies(temp)
        AxBO = self.ax_bo
        ByBO = self.by_bo
        CzBO = self.cz_bo
        BxBO = self.bx_bo
        CxBO = self.cx_bo
        CyBO = self.cy_bo

        # Reference structure
        Ax1 = self.ax[1,1,1,1,1,1]
        Bx1 = self.bx[1,1,1,1,1,1]
        By1 = self.by[1,1,1,1,1,1]
        Cx1 = self.cx[1,1,1,1,1,1]
        Cy1 = self.cy[1,1,1,1,1,1]
        Cz1 = self.cz[1,1,1,1,1,1]

        # Adjust lattice variation from deformed structures based on table IV and eq(46)
        Ax0 = self.ax[0,1,1,1,1,1] # reference structure - exx0
        By0 = self.by[1,0,1,1,1,1] # reference structure - eyy0
        Cz0 = self.cz[1,1,0,1,1,1] # reference structure - ezz0
        Bx0 = self.bx[1,1,1,1,1,0] # reference structure - exy0
        Cx0 = self.cx[1,1,1,1,0,1] # reference structure - exz0
        Cy0 = self.cy[1,1,1,0,1,1] # reference structure - eyz0

        # Compute strain-related quantities. Eq (47)
        # Strain steps and
        dexx = (Ax0-Ax1)/AxBO
        deyy = (By0-By1)/ByBO
        dezz = (Cz0-Cz1)/CzBO
        #dexz = (AxBO*(Cx0-Cx1)-(Ax0-Ax1)*CxBO)/(AxBO*CzBO)
        dexz = (Cx0-Cx1)/CzBO

        # Reference strains
        exx0 = Ax1/AxBO-1
        eyy0 = By1/ByBO-1
        ezz0 = Cz1/CzBO-1
        exz0 = (AxBO*Cx1-Ax1*CxBO)/(AxBO*CzBO)
        if self.verbose:
            print(f"{dexx=}, {deyy=}, {dezz=}, {dexz=}")
            print(f"{exx0=}, {eyy0=}, {ezz0=}, {exz0=}")

        # Compute first and second derivatives of free energy w.r.t. strains
        dF_dA1 = (e[0,1,1,1,1,1]-e[2,1,1,1,1,1])/(2*dexx)
        dF_dB2 = (e[1,0,1,1,1,1]-e[1,2,1,1,1,1])/(2*deyy)
        dF_dC3 = (e[1,1,0,1,1,1]-e[1,1,2,1,1,1])/(2*dezz)
        dF_dC1 = (e[1,1,1,1,0,1]-e[1,1,1,1,2,1])/(2*dexz)

        d2F_dA12 = (e[0,1,1,1,1,1]-2*e[1,1,1,1,1,1]+e[2,1,1,1,1,1])/(dexx)**2
        d2F_dB22 = (e[1,0,1,1,1,1]-2*e[1,1,1,1,1,1]+e[1,2,1,1,1,1])/(deyy)**2
        d2F_dC32 = (e[1,1,0,1,1,1]-2*e[1,1,1,1,1,1]+e[1,1,2,1,1,1])/(dezz)**2
        d2F_dC12 = (e[1,1,1,1,0,1]-2*e[1,1,1,1,1,1]+e[1,1,1,1,2,1])/(dexz)**2

        d2F_dA1dB2 = (e[1,1,1,1,1,1] - e[0,1,1,1,1,1] - e[1,0,1,1,1,1] + e[0,0,1,1,1,1]) / (dexx *deyy)
        d2F_dA1dC3 = (e[1,1,1,1,1,1] - e[0,1,1,1,1,1] - e[1,1,0,1,1,1] + e[0,1,0,1,1,1]) / (dexx *dezz)
        d2F_dA1dC1 = (e[1,1,1,1,1,1] - e[0,1,1,1,1,1] - e[1,1,1,1,0,1] + e[0,1,1,1,0,1]) / (dexx *dexz)

        d2F_dB2dC3 = (e[1,1,1,1,1,1] - e[1,0,1,1,1,1] - e[1,1,0,1,1,1] + e[1,0,0,1,1,1]) / (deyy *dezz)
        d2F_dB2dC1 = (e[1,1,1,1,1,1] - e[1,0,1,1,1,1] - e[1,1,1,1,0,1] + e[1,0,1,1,0,1]) / (deyy *dexz)
        d2F_dC3dC1 = (e[1,1,1,1,1,1] - e[1,1,0,1,1,1] - e[1,1,1,1,0,1] + e[1,1,0,1,0,1]) / (dezz *dexz)

        dS_dA1 = (S[0,1,1,1,1,1]-S[2,1,1,1,1,1])/(2*dexx)
        dS_dB2 = (S[1,0,1,1,1,1]-S[1,2,1,1,1,1])/(2*deyy)
        dS_dC3 = (S[1,1,0,1,1,1]-S[1,1,2,1,1,1])/(2*dezz)
        dS_dC1 = (S[1,1,1,1,0,1]-S[1,1,1,1,2,1])/(2*dexz)

        d2S_dA12 = (S[0,1,1,1,1,1]-2*S[1,1,1,1,1,1]+S[2,1,1,1,1,1])/(dexx)**2
        d2S_dB22 = (S[1,0,1,1,1,1]-2*S[1,1,1,1,1,1]+S[1,2,1,1,1,1])/(deyy)**2
        d2S_dC32 = (S[1,1,0,1,1,1]-2*S[1,1,1,1,1,1]+S[1,1,2,1,1,1])/(dezz)**2
        d2S_dC12 = (S[1,1,1,1,0,1]-2*S[1,1,1,1,1,1]+S[1,1,1,1,2,1])/(dexz)**2

        d2S_dA1dB2 = (S[1,1,1,1,1,1] - S[0,1,1,1,1,1] - S[1,0,1,1,1,1] + S[0,0,1,1,1,1]) / (dexx *deyy)
        d2S_dA1dC3 = (S[1,1,1,1,1,1] - S[0,1,1,1,1,1] - S[1,1,0,1,1,1] + S[0,1,0,1,1,1]) / (dexx *dezz)
        d2S_dA1dC1 = (S[1,1,1,1,1,1] - S[0,1,1,1,1,1] - S[1,1,1,1,0,1] + S[0,1,1,1,0,1]) / (dexx *dexz)

        d2S_dB2dC3 = (S[1,1,1,1,1,1] - S[1,0,1,1,1,1] - S[1,1,0,1,1,1] + S[1,0,0,1,1,1]) / (deyy *dezz)
        d2S_dB2dC1 = (S[1,1,1,1,1,1] - S[1,0,1,1,1,1] - S[1,1,1,1,0,1] + S[1,0,1,1,0,1]) / (deyy *dexz)

        d2S_dC3dC1 = (S[1,1,1,1,1,1] - S[1,1,0,1,1,1] - S[1,1,1,1,0,1] + S[1,1,0,1,0,1]) / (dezz *dexz)

        a = self.lattice_a_guess
        b = self.lattice_b_guess
        c = self.lattice_c_guess
        v = self.volume_guess

        # Reconstruct the lattice vectors to find the standard lattice vector of monoclinic
        ax = a
        by = b
        cz = c * math.sin(math.pi*self.angles_guess[1]/180)
        cx = c * math.cos(math.pi*self.angles_guess[1]/180)

        # compute BO strain at guess, Eq (47)
        exx_n = ax/AxBO-1
        eyy_n = by/ByBO-1
        ezz_n = cz/CzBO-1
        exz_n = (AxBO*cx-ax*CxBO)/(AxBO*CzBO)

        # Free energy and entropy derivative at guess structure
        dfda1 = dF_dA1 + (exx_n-exx0)*d2F_dA12+(eyy_n-eyy0)*d2F_dA1dB2+(ezz_n-ezz0)*d2F_dA1dC3+(exz_n-exz0)*d2F_dA1dC1
        dfdb2 = dF_dB2 + (eyy_n-eyy0)*d2F_dB22+(exx_n-exx0)*d2F_dA1dB2+(ezz_n-ezz0)*d2F_dB2dC3+(exz_n-exz0)*d2F_dB2dC1
        dfdc3 = dF_dC3 + (ezz_n-ezz0)*d2F_dC32+(exx_n-exx0)*d2F_dA1dC3+(eyy_n-eyy0)*d2F_dB2dC3+(exz_n-exz0)*d2F_dC3dC1
        dfdc1 = dF_dC1 + (exz_n-exz0)*d2F_dC12+(exx_n-exx0)*d2F_dA1dC1+(eyy_n-eyy0)*d2F_dB2dC1+(ezz_n-ezz0)*d2F_dC3dC1

        dsda1 = dS_dA1 + (exx_n-exx0)*d2S_dA12+(eyy_n-eyy0)*d2S_dA1dB2+(ezz_n-ezz0)*d2S_dA1dC3+(exz_n-exz0)*d2S_dA1dC1
        dsdb2 = dS_dB2 + (eyy_n-eyy0)*d2S_dB22+(exx_n-exx0)*d2S_dA1dB2+(ezz_n-ezz0)*d2S_dB2dC3+(exz_n-exz0)*d2S_dB2dC1
        dsdc3 = dS_dC3 + (ezz_n-ezz0)*d2S_dC32+(exx_n-exx0)*d2S_dA1dC3+(eyy_n-eyy0)*d2S_dB2dC3+(exz_n-exz0)*d2S_dC3dC1
        dsdc1 = dS_dC1 + (exz_n-exz0)*d2S_dC12+(exx_n-exx0)*d2S_dA1dC1+(eyy_n-eyy0)*d2S_dB2dC1+(ezz_n-ezz0)*d2S_dC3dC1

        # Gibbs free energy
        gibbs = self.energy_guess + e[1,1,1,1,1,1] + pressure*v / abu.eVA3_HaBohr3 + \
                  (exx_n-exx0)*dF_dA1 + 0.5*(exx_n-exx0)**2*d2F_dA12 +\
                  (eyy_n-eyy0)*dF_dB2 + 0.5*(eyy_n-eyy0)**2*d2F_dB22 +\
                  (ezz_n-ezz0)*dF_dC3 + 0.5*(ezz_n-ezz0)**2*d2F_dC32 +\
                  (exz_n-exz0)*dF_dC1 + 0.5*(exz_n-exz0)**2*d2F_dC12 +\
                  (exx_n-exx0)*(eyy_n-eyy0)*d2F_dA1dB2 +\
                  (exx_n-exx0)*(ezz_n-ezz0)*d2F_dA1dC3 +\
                  (exx_n-exx0)*(exz_n-exz0)*d2F_dA1dC1 +\
                  (eyy_n-eyy0)*(ezz_n-ezz0)*d2F_dB2dC3 +\
                  (eyy_n-eyy0)*(exz_n-exz0)*d2F_dB2dC1 +\
                  (ezz_n-ezz0)*(exz_n-exz0)*d2F_dC3dC1

        # Compute thermal stresses. Eq (47)
        stress_a1 = -dfda1/v*(exx_n+1) * abu.eVA3_HaBohr3
        stress_b2 = -dfdb2/v*(eyy_n+1) * abu.eVA3_HaBohr3
        stress_c3 = -dfdc3/v*(ezz_n+1) * abu.eVA3_HaBohr3
        stress_c1 = -1.0/v*(dfdc1*(ezz_n+1)+dfda1*exz_n) * abu.eVA3_HaBohr3
        if self.verbose:
            print("ax/AxBO, by/ByBO, cz/CzBO: ", ax/AxBO, by/ByBO, cz/CzBO)
            print(f"{stress_a1=}, {stress_b2=}, {stress_c3=}, {stress_c1=}")

        therm = None

        stress = np.zeros(6)
        stress[0] = stress_a1 -pressure
        stress[1] = stress_b2 -pressure
        stress[2] = stress_c3 -pressure
        stress[4] = stress_c1

        dtol = np.zeros(6)
        dtol[0] = abs(stress[0] - self.stress_guess[0,0])
        dtol[1] = abs(stress[1] - self.stress_guess[1,1])
        dtol[2] = abs(stress[2] - self.stress_guess[2,2])
        dtol[4] = abs(stress[4] - self.stress_guess[2,0])

        # Check if the stress has converged (all tolerances below self.dtol_tolerance)
        elastic = None
        if all(dtol[i] < self.dtol_tolerance for i in range(6)):
            if self.bo_elastic_voigt is not None:
                # If elastic constants are requested, compute the second derivatives for C44, C66, and C46
                if mode == 'ECs':
                    dexy = (Bx0-Bx1)/(ByBO)
                    deyz = (Cy0-Cy1)/(CzBO)
                    d2F_dB12 = (e[1,1,1,1,1,0]-2*e[1,1,1,1,1,1]+e[1,1,1,1,1,0])/(deyz)**2
                    d2F_dC22 = (e[1,1,1,0,1,1]-2*e[1,1,1,1,1,1]+e[1,1,1,0,1,1])/(dexy)**2
                    d2F_dB1dC2 = (e[1,1,1,1,1,1] - e[1,1,1,0,1,1] - e[1,1,1,1,1,0] + e[1,1,1,0,1,0]) / (dexy *deyz)
                else:
                    d2F_dB12 =  0.0
                    d2F_dC22 =  0.0
                    d2F_dB1dC2 = 0.0

                # Read elastic constants from the output of DFPT obtained by abiopen.py (ELASTIC_RELAXED)
                matrix_elastic = self.bo_elastic_voigt
                # Convert elastic constants to second derivative of BO energy (Eq. 39)
                matrix_elastic = matrix_elastic*v/abu.eVA3_GPa
                scale_xx = (exx0+1)*(exx0+1)
                scale_yy = (eyy0+1)*(eyy0+1)
                scale_zz = (ezz0+1)*(ezz0+1)
                scale_xz = (exx0+1)*(ezz0+1)
                scale_yz = (eyy0+1)*(ezz0+1)
                scale_xy = (exx0+1)*(eyy0+1)
                # Initialize the second derivative matrix M using free energy second derivatives.
                # d2F_dC22, d2F_dB1dC2, d2F_dB12 are zero in 'TEC' mode and do not contribute to thermal expansion calculations
                # due to monoclinic symmetry, where the first derivatives of free energy and entropy in the yz and xy directions are zero.
                M = np.array([[ d2F_dA12 *scale_xx, d2F_dA1dB2*scale_xy ,d2F_dA1dC3*scale_xz ,0.0                 ,d2F_dA1dC1*scale_xz ,0.0                ],
                              [d2F_dA1dB2*scale_xy, d2F_dB22  *scale_yy ,d2F_dB2dC3*scale_yz ,0.0                 ,d2F_dB2dC1*scale_yz ,0.0                ],
                              [d2F_dA1dC3*scale_xz, d2F_dB2dC3*scale_yz ,d2F_dC32  *scale_zz ,0.0                 ,d2F_dC3dC1*scale_zz ,0.0                ],
                              [0.0                , 0.0                 ,0.0                 ,d2F_dC22  *scale_zz ,0.0                 ,d2F_dB1dC2*scale_yz],
                              [d2F_dA1dC1*scale_xz, d2F_dB2dC1*scale_yz ,d2F_dC3dC1*scale_zz ,0.0                 ,d2F_dC12  *scale_zz ,0.0                ],
                              [0.0                , 0.0                 ,0.0                 ,d2F_dB1dC2*scale_yz ,0.0                 ,d2F_dB12  *scale_yy]])
                # Contribution of pressure to the second derivative matrix
                P = np.array([[0.0 ,pressure*v,pressure*v  ,0,0,0],
                              [pressure*v ,0.0,pressure*v  ,0,0,0],
                              [pressure*v ,pressure*v,0.0  ,0,0,0],
                              [0,0,0,0,0,0],
                              [0,0,0,0,0,0],
                              [0,0,0,0,0,0]])
                # Compute the final second derivative matrix M, including elastic and pressure contributions
                M = M + matrix_elastic + P / abu.eVA3_HaBohr3
                # Scale the first derivative of entropy
                S1 = dsda1*(exx_n+1)
                S2 = dsdb2*(eyy_n+1)
                S3 = dsdc3*(ezz_n+1)
                S4 = 0.0
                S5 = dsdc1*(ezz_n+1)+dsda1*exz_n
                S6 = 0.0
                dSde = np.array([S1,S2,S3,S4,S5,S6])
                # Compute thermal expansion using Eq. (37)
                dstrain_dt = np.linalg.inv(M) @ dSde
                # Scale the thermal expansion (TOFIX)
                #therm=[dstrain_dt[0]*(exx_n+1) , dstrain_dt[1]*(eyy_n+1),dstrain_dt[2]*(ezz_n+1),0,dstrain_dt[4]*(ezz_n+1)+dstrain_dt[0]*exz_n,0]
                therm = [dstrain_dt[0]*(exx_n+1), dstrain_dt[1]*(eyy_n+1),dstrain_dt[2]*(ezz_n+1)-(dstrain_dt[4]*(ezz_n+1)+dstrain_dt[0]*exz_n)/cz,0,-dstrain_dt[4],0]
                elastic = M / v*abu.eVA3_GPa

        return dtol, gibbs, stress, therm, elastic

    def stress_ZSISA_triclinic(self, temp: float, pressure: float, mode: str) -> tuple:

        # Get vibrational free energy and entropy at a specific temperature
        # e = Vibrational free energy (F_vib)
        # S = Entropy (S)
        e, S = self.get_vib_free_energies(temp)

        AxBO = self.ax_bo
        ByBO = self.by_bo
        CzBO = self.cz_bo
        BxBO = self.bx_bo
        CxBO = self.cx_bo
        CyBO = self.cy_bo

        # Adjust lattice variation from deformed structures based on table IV and eq(55)
        Ax0 = self.ax[0,1,1,1,1,1] # reference structure - exx0
        Bx0 = self.bx[1,1,1,1,1,0] # reference structure - exy0
        By0 = self.by[1,0,1,1,1,1] # reference structure - eyy0
        Cx0 = self.cx[1,1,1,1,0,1] # reference structure - exz0
        Cy0 = self.cy[1,1,1,0,1,1] # reference structure - eyz0
        Cz0 = self.cz[1,1,0,1,1,1] # reference structure - ezz0

        # Reference structure
        Ax1 = self.ax[1,1,1,1,1,1]
        Bx1 = self.bx[1,1,1,1,1,1]
        By1 = self.by[1,1,1,1,1,1]
        Cx1 = self.cx[1,1,1,1,1,1]
        Cy1 = self.cy[1,1,1,1,1,1]
        Cz1 = self.cz[1,1,1,1,1,1]

        # Compute strain-related quantities. Eq (56)
        # Strain steps and
        dexx = (Ax0-Ax1)/AxBO
        deyy = (By0-By1)/ByBO
        dezz = (Cz0-Cz1)/CzBO
        dexy = (Bx0-Bx1)/ByBO
        deyz = (Cy0-Cy1)/CzBO
        dexz = (Cx0-Cx1)/CzBO

        # Reference strains
        exx0 = Ax1/AxBO-1
        eyy0 = By1/ByBO-1
        ezz0 = Cz1/CzBO-1
        exy0 = (AxBO*Bx1-Ax1*BxBO)/(AxBO*ByBO)
        eyz0 = (ByBO*Cy1-By1*CyBO)/(ByBO*CzBO)
        exz0 = (AxBO*(ByBO*Cx1-Bx1*CyBO)-Ax1*(ByBO*CxBO-BxBO*CyBO))/(AxBO*ByBO*CzBO)
        if self.verbose:
            print(f"{dexx=}, {deyy=}, {dezz=}, {dexz=}, {deyz=}, {dexz=}")
            print(f"{exx0=}, {eyy0=}, {ezz0=}, {exz0=}, {eyz0=}, {exz0=}")

        # Compute first and second derivatives of free energy w.r.t. strains
        dF_dA1 = (e[0,1,1,1,1,1]-e[2,1,1,1,1,1])/(2*dexx)
        dF_dB2 = (e[1,0,1,1,1,1]-e[1,2,1,1,1,1])/(2*deyy)
        dF_dC3 = (e[1,1,0,1,1,1]-e[1,1,2,1,1,1])/(2*dezz)
        dF_dC2 = (e[1,1,1,0,1,1]-e[1,1,1,2,1,1])/(2*deyz)
        dF_dC1 = (e[1,1,1,1,0,1]-e[1,1,1,1,2,1])/(2*dexz)
        dF_dB1 = (e[1,1,1,1,1,0]-e[1,1,1,1,1,2])/(2*dexy)

        d2F_dA12 = (e[0,1,1,1,1,1]-2*e[1,1,1,1,1,1]+e[2,1,1,1,1,1])/(dexx)**2
        d2F_dB22 = (e[1,0,1,1,1,1]-2*e[1,1,1,1,1,1]+e[1,2,1,1,1,1])/(deyy)**2
        d2F_dC32 = (e[1,1,0,1,1,1]-2*e[1,1,1,1,1,1]+e[1,1,2,1,1,1])/(dezz)**2
        d2F_dC22 = (e[1,1,1,0,1,1]-2*e[1,1,1,1,1,1]+e[1,1,1,2,1,1])/(deyz)**2
        d2F_dC12 = (e[1,1,1,1,0,1]-2*e[1,1,1,1,1,1]+e[1,1,1,1,2,1])/(dexz)**2
        d2F_dB12 = (e[1,1,1,1,1,0]-2*e[1,1,1,1,1,1]+e[1,1,1,1,1,2])/(dexy)**2

        d2F_dA1dB2 = (e[1,1,1,1,1,1] - e[0,1,1,1,1,1] - e[1,0,1,1,1,1] + e[0,0,1,1,1,1]) / (dexx *deyy)
        d2F_dA1dC3 = (e[1,1,1,1,1,1] - e[0,1,1,1,1,1] - e[1,1,0,1,1,1] + e[0,1,0,1,1,1]) / (dexx *dezz)
        d2F_dA1dC2 = (e[1,1,1,1,1,1] - e[0,1,1,1,1,1] - e[1,1,1,0,1,1] + e[0,1,1,0,1,1]) / (dexx *deyz)
        d2F_dA1dC1 = (e[1,1,1,1,1,1] - e[0,1,1,1,1,1] - e[1,1,1,1,0,1] + e[0,1,1,1,0,1]) / (dexx *dexz)
        d2F_dA1dB1 = (e[1,1,1,1,1,1] - e[0,1,1,1,1,1] - e[1,1,1,1,1,0] + e[0,1,1,1,1,0]) / (dexx *dexy)

        d2F_dB2dC3 = (e[1,1,1,1,1,1] - e[1,0,1,1,1,1] - e[1,1,0,1,1,1] + e[1,0,0,1,1,1]) / (deyy *dezz)
        d2F_dB2dC2 = (e[1,1,1,1,1,1] - e[1,0,1,1,1,1] - e[1,1,1,0,1,1] + e[1,0,1,0,1,1]) / (deyy *deyz)
        d2F_dB2dC1 = (e[1,1,1,1,1,1] - e[1,0,1,1,1,1] - e[1,1,1,1,0,1] + e[1,0,1,1,0,1]) / (deyy *dexz)
        d2F_dB2dB1 = (e[1,1,1,1,1,1] - e[1,0,1,1,1,1] - e[1,1,1,1,1,0] + e[1,0,1,1,1,0]) / (deyy *dexy)

        d2F_dC3dC2 = (e[1,1,1,1,1,1] - e[1,1,0,1,1,1] - e[1,1,1,0,1,1] + e[1,1,0,0,1,1]) / (dezz *deyz)
        d2F_dC3dC1 = (e[1,1,1,1,1,1] - e[1,1,0,1,1,1] - e[1,1,1,1,0,1] + e[1,1,0,1,0,1]) / (dezz *dexz)
        d2F_dC3dB1 = (e[1,1,1,1,1,1] - e[1,1,0,1,1,1] - e[1,1,1,1,1,0] + e[1,1,0,1,1,0]) / (dezz *dexy)

        d2F_dC1dC2 = (e[1,1,1,1,1,1] - e[1,1,1,0,1,1] - e[1,1,1,1,0,1] + e[1,1,1,0,0,1]) / (deyz *dexz)
        d2F_dB1dC2 = (e[1,1,1,1,1,1] - e[1,1,1,0,1,1] - e[1,1,1,1,1,0] + e[1,1,1,0,1,0]) / (deyz *dexy)

        d2F_dB1dC1 = (e[1,1,1,1,1,1] - e[1,1,1,1,0,1] - e[1,1,1,1,1,0] + e[1,1,1,1,0,0]) / (dexz *dexy)

        dS_dA1 = (S[0,1,1,1,1,1]-S[2,1,1,1,1,1])/(2*dexx)
        dS_dB2 = (S[1,0,1,1,1,1]-S[1,2,1,1,1,1])/(2*deyy)
        dS_dC3 = (S[1,1,0,1,1,1]-S[1,1,2,1,1,1])/(2*dezz)
        dS_dC2 = (S[1,1,1,0,1,1]-S[1,1,1,2,1,1])/(2*deyz)
        dS_dC1 = (S[1,1,1,1,0,1]-S[1,1,1,1,2,1])/(2*dexz)
        dS_dB1 = (S[1,1,1,1,1,0]-S[1,1,1,1,1,2])/(2*dexy)

        d2S_dA12 = (S[0,1,1,1,1,1]-2*S[1,1,1,1,1,1]+S[2,1,1,1,1,1])/(dexx)**2
        d2S_dB22 = (S[1,0,1,1,1,1]-2*S[1,1,1,1,1,1]+S[1,2,1,1,1,1])/(deyy)**2
        d2S_dC32 = (S[1,1,0,1,1,1]-2*S[1,1,1,1,1,1]+S[1,1,2,1,1,1])/(dezz)**2
        d2S_dC22 = (S[1,1,1,0,1,1]-2*S[1,1,1,1,1,1]+S[1,1,1,2,1,1])/(deyz)**2
        d2S_dC12 = (S[1,1,1,1,0,1]-2*S[1,1,1,1,1,1]+S[1,1,1,1,2,1])/(dexz)**2
        d2S_dB12 = (S[1,1,1,1,1,0]-2*S[1,1,1,1,1,1]+S[1,1,1,1,1,2])/(dexy)**2

        d2S_dA1dB2 = (S[1,1,1,1,1,1] - S[0,1,1,1,1,1] - S[1,0,1,1,1,1] + S[0,0,1,1,1,1]) / (dexx *deyy)
        d2S_dA1dC3 = (S[1,1,1,1,1,1] - S[0,1,1,1,1,1] - S[1,1,0,1,1,1] + S[0,1,0,1,1,1]) / (dexx *dezz)
        d2S_dA1dC2 = (S[1,1,1,1,1,1] - S[0,1,1,1,1,1] - S[1,1,1,0,1,1] + S[0,1,1,0,1,1]) / (dexx *deyz)
        d2S_dA1dC1 = (S[1,1,1,1,1,1] - S[0,1,1,1,1,1] - S[1,1,1,1,0,1] + S[0,1,1,1,0,1]) / (dexx *dexz)
        d2S_dA1dB1 = (S[1,1,1,1,1,1] - S[0,1,1,1,1,1] - S[1,1,1,1,1,0] + S[0,1,1,1,1,0]) / (dexx *dexy)

        d2S_dB2dC3 = (S[1,1,1,1,1,1] - S[1,0,1,1,1,1] - S[1,1,0,1,1,1] + S[1,0,0,1,1,1]) / (deyy *dezz)
        d2S_dB2dC2 = (S[1,1,1,1,1,1] - S[1,0,1,1,1,1] - S[1,1,1,0,1,1] + S[1,0,1,0,1,1]) / (deyy *deyz)
        d2S_dB2dC1 = (S[1,1,1,1,1,1] - S[1,0,1,1,1,1] - S[1,1,1,1,0,1] + S[1,0,1,1,0,1]) / (deyy *dexz)
        d2S_dB2dB1 = (S[1,1,1,1,1,1] - S[1,0,1,1,1,1] - S[1,1,1,1,1,0] + S[1,0,1,1,1,0]) / (deyy *dexy)

        d2S_dC3dC2 = (S[1,1,1,1,1,1] - S[1,1,0,1,1,1] - S[1,1,1,0,1,1] + S[1,1,0,0,1,1]) / (dezz *deyz)
        d2S_dC3dC1 = (S[1,1,1,1,1,1] - S[1,1,0,1,1,1] - S[1,1,1,1,0,1] + S[1,1,0,1,0,1]) / (dezz *dexz)
        d2S_dC3dB1 = (S[1,1,1,1,1,1] - S[1,1,0,1,1,1] - S[1,1,1,1,1,0] + S[1,1,0,1,1,0]) / (dezz *dexy)

        d2S_dC1dC2 = (S[1,1,1,1,1,1] - S[1,1,1,0,1,1] - S[1,1,1,1,0,1] + S[1,1,1,0,0,1]) / (deyz *dexz)
        d2S_dB1dC2 = (S[1,1,1,1,1,1] - S[1,1,1,0,1,1] - S[1,1,1,1,1,0] + S[1,1,1,0,1,0]) / (deyz *dexy)

        d2S_dB1dC1 = (S[1,1,1,1,1,1] - S[1,1,1,1,0,1] - S[1,1,1,1,1,0] + S[1,1,1,1,0,0]) / (dexz *dexy)

        # Reconstruct the lattice vectors to find the standard lattice vector of monoclinic
        a = self.lattice_a_guess
        b = self.lattice_b_guess
        c = self.lattice_c_guess

        cos_ab = math.cos(math.pi*self.angles_guess[2]/180)
        cos_ac = math.cos(math.pi*self.angles_guess[1]/180)
        cos_bc = math.cos(math.pi*self.angles_guess[0]/180)

        ax = 1.0
        bx = cos_ab
        by = np.sqrt(1-cos_ab**2)
        cx = cos_ac
        cy = (cos_bc-bx*cx)/by
        cz = np.sqrt(1.0-cx**2-cy**2)
        ax = ax*a
        bx = bx*b
        by = by*b
        cx = cx*c
        cy = cy*c
        cz = cz*c
        v = self.volume_guess

        # Compute BO strain at guess, Eq (56)
        exx_n = ax/Ax0-1
        eyy_n = by/By0-1
        ezz_n = cz/Cz0-1
        exy_n = (AxBO*bx-ax*BxBO)/(AxBO*ByBO)
        eyz_n = (ByBO*cy-by*CyBO)/(ByBO*CzBO)
        exz_n = (AxBO*(ByBO*cx-bx*CyBO)-ax*(ByBO*CxBO-BxBO*CyBO))/(AxBO*ByBO*CzBO)

        dfda1 = dF_dA1 + (exx_n-exx0)*d2F_dA12+(eyy_n-eyy0)*d2F_dA1dB2+(ezz_n-ezz0)*d2F_dA1dC3+(exy_n-exy0)*d2F_dA1dB1+(exz_n-exz0)*d2F_dA1dC1+(eyz_n-eyz0)*d2F_dA1dC2
        dfdb2 = dF_dB2 + (eyy_n-eyy0)*d2F_dB22+(exx_n-exx0)*d2F_dA1dB2+(ezz_n-ezz0)*d2F_dB2dC3+(exy_n-exy0)*d2F_dB2dB1+(exz_n-exz0)*d2F_dB2dC1+(eyz_n-eyz0)*d2F_dB2dC2
        dfdc3 = dF_dC3 + (ezz_n-ezz0)*d2F_dC32+(exx_n-exx0)*d2F_dA1dC3+(eyy_n-eyy0)*d2F_dB2dC3+(exy_n-exy0)*d2F_dC3dB1+(exz_n-exz0)*d2F_dC3dC1+(eyz_n-eyz0)*d2F_dC3dC2
        dfdb1 = dF_dB1 + (exy_n-exy0)*d2F_dB12+(exx_n-exx0)*d2F_dA1dB1+(eyy_n-eyy0)*d2F_dB2dB1+(ezz_n-ezz0)*d2F_dC3dB1+(exz_n-exz0)*d2F_dB1dC1+(eyz_n-eyz0)*d2F_dB1dC2
        dfdc1 = dF_dC1 + (exz_n-exz0)*d2F_dC12+(exx_n-exx0)*d2F_dA1dC1+(eyy_n-eyy0)*d2F_dB2dC1+(ezz_n-ezz0)*d2F_dC3dC1+(exy_n-exy0)*d2F_dB1dC1+(eyz_n-eyz0)*d2F_dC1dC2
        dfdc2 = dF_dC2 + (eyz_n-eyz0)*d2F_dC22+(exx_n-exx0)*d2F_dA1dC2+(eyy_n-eyy0)*d2F_dB2dC2+(ezz_n-ezz0)*d2F_dC3dC2+(exy_n-exy0)*d2F_dB1dC2+(exz_n-exz0)*d2F_dC1dC2

        dsda1 = dS_dA1 + (exx_n-exx0)*d2S_dA12+(eyy_n-eyy0)*d2S_dA1dB2+(ezz_n-ezz0)*d2S_dA1dC3+(exy_n-exy0)*d2S_dA1dB1+(exz_n-exz0)*d2S_dA1dC1+(eyz_n-eyz0)*d2S_dA1dC2
        dsdb2 = dS_dB2 + (eyy_n-eyy0)*d2S_dB22+(exx_n-exx0)*d2S_dA1dB2+(ezz_n-ezz0)*d2S_dB2dC3+(exy_n-exy0)*d2S_dB2dB1+(exz_n-exz0)*d2S_dB2dC1+(eyz_n-eyz0)*d2S_dB2dC2
        dsdc3 = dS_dC3 + (ezz_n-ezz0)*d2S_dC32+(exx_n-exx0)*d2S_dA1dC3+(eyy_n-eyy0)*d2S_dB2dC3+(exy_n-exy0)*d2S_dC3dB1+(exz_n-exz0)*d2S_dC3dC1+(eyz_n-eyz0)*d2S_dC3dC2
        dsdb1 = dS_dB1 + (exy_n-exy0)*d2S_dB12+(exx_n-exx0)*d2S_dA1dB1+(eyy_n-eyy0)*d2S_dB2dB1+(ezz_n-ezz0)*d2S_dC3dB1+(exz_n-exz0)*d2S_dB1dC1+(eyz_n-eyz0)*d2S_dB1dC2
        dsdc1 = dS_dC1 + (exz_n-exz0)*d2S_dC12+(exx_n-exx0)*d2S_dA1dC1+(eyy_n-eyy0)*d2S_dB2dC1+(ezz_n-ezz0)*d2S_dC3dC1+(exy_n-exy0)*d2S_dB1dC1+(eyz_n-eyz0)*d2S_dC1dC2
        dsdc2 = dS_dC2 + (eyz_n-eyz0)*d2S_dC22+(exx_n-exx0)*d2S_dA1dC2+(eyy_n-eyy0)*d2S_dB2dC2+(ezz_n-ezz0)*d2S_dC3dC2+(exy_n-exy0)*d2S_dB1dC2+(exz_n-exz0)*d2S_dC1dC2

        # Gibbs free energy
        gibbs = self.energy_guess + e[1,1,1,1,1,1] + pressure*v / abu.eVA3_HaBohr3 + \
         (exx_n-exx0)*dF_dA1 + 0.5*(exx_n-exx0)**2*d2F_dA12 +\
         (eyy_n-eyy0)*dF_dB2 + 0.5*(eyy_n-eyy0)**2*d2F_dB22 +\
         (ezz_n-ezz0)*dF_dC3 + 0.5*(ezz_n-ezz0)**2*d2F_dC32 +\
         (exy_n-exy0)*dF_dB1 + 0.5*(exy_n-exy0)**2*d2F_dB12 +\
         (exz_n-exz0)*dF_dC1 + 0.5*(exz_n-exz0)**2*d2F_dC12 +\
         (eyz_n-eyz0)*dF_dC2 + 0.5*(eyz_n-eyz0)**2*d2F_dC22 +\
         (exx_n-exx0)*(eyy_n-eyy0)*d2F_dA1dB2 +\
         (exx_n-exx0)*(ezz_n-ezz0)*d2F_dA1dC3 +\
         (exx_n-exx0)*(exy_n-exy0)*d2F_dA1dB1 +\
         (exx_n-exx0)*(exz_n-exz0)*d2F_dA1dC1 +\
         (exx_n-exx0)*(eyz_n-eyz0)*d2F_dA1dC2 +\
         (eyy_n-eyy0)*(ezz_n-ezz0)*d2F_dB2dC3 +\
         (eyy_n-eyy0)*(exy_n-exy0)*d2F_dB2dB1 +\
         (eyy_n-eyy0)*(exz_n-exz0)*d2F_dB2dC1 +\
         (eyy_n-eyy0)*(eyz_n-eyz0)*d2F_dB2dC2 +\
         (ezz_n-ezz0)*(exy_n-exy0)*d2F_dC3dB1 +\
         (ezz_n-ezz0)*(exz_n-exz0)*d2F_dC3dC1 +\
         (ezz_n-ezz0)*(eyz_n-eyz0)*d2F_dC3dC2 +\
         (exy_n-exy0)*(exz_n-exz0)*d2F_dB1dC1 +\
         (exy_n-exy0)*(eyz_n-eyz0)*d2F_dB1dC2 +\
         (exz_n-exz0)*(eyz_n-eyz0)*d2F_dC1dC2

        # Compute thermal stresses. Eq (57)
        stress_a1 = -dfda1/v*(exx_n+1) * abu.eVA3_HaBohr3
        stress_b2 = -dfdb2/v*(eyy_n+1) * abu.eVA3_HaBohr3
        stress_c3 = -dfdc3/v*(ezz_n+1) * abu.eVA3_HaBohr3
        stress_b1 = -1.0/v*(dfdb1*(eyy_n+1)+dfda1*exy_n) * abu.eVA3_HaBohr3
        stress_c2 = -1.0/v*(dfdc2*(ezz_n+1)+dfdb2*eyz_n) * abu.eVA3_HaBohr3
        stress_c1 = -1.0/v*(dfdc1*(ezz_n+1)+dfdb1*eyz_n+dfda1*exz_n) * abu.eVA3_HaBohr3

        stress = np.zeros(6)
        stress[0] = stress_a1 -pressure
        stress[1] = stress_b2 -pressure
        stress[2] = stress_c3 -pressure
        stress[3] = stress_c2
        stress[4] = stress_c1
        stress[5] = stress_b1

        if self.verbose:
            print("ax/Ax0, by/By0, cz/Cz0: ", ax/Ax0, by/By0, cz/Cz0)
            print(f"{stress=}")

        dtol = np.zeros(6)
        dtol[0] = abs(stress[0]-self.stress_guess[0,0])
        dtol[1] = abs(stress[1]-self.stress_guess[1,1])
        dtol[2] = abs(stress[2]-self.stress_guess[2,2])
        dtol[3] = abs(stress[3]-self.stress_guess[2,1])
        dtol[4] = abs(stress[4]-self.stress_guess[2,0])
        dtol[5] = abs(stress[5]-self.stress_guess[1,0])
        therm = None

        # Check if the stress has converged (all tolerances below self.dtol_tolerance)
        elastic=None
        if all(dtol[i] < self.dtol_tolerance for i in range(6)):
            if self.bo_elastic_voigt is not None:
                # If elastic constants are requested, compute the second derivatives for C44, C66, and C46
                # Read elastic constants from the output of DFPT obtained by abiopen.py (ELASTIC_RELAXED)
                matrix_elastic = self.bo_elastic_voigt
                # Convert elastic constants to second derivative of BO energy (Eq. 39)
                matrix_elastic = matrix_elastic*v/abu.eVA3_GPa
                scale_xx = (exx0+1)*(exx0+1)
                scale_yy = (eyy0+1)*(eyy0+1)
                scale_zz = (ezz0+1)*(ezz0+1)
                scale_xz = (exx0+1)*(ezz0+1)
                scale_yz = (eyy0+1)*(ezz0+1)
                scale_xy = (exx0+1)*(eyy0+1)
                # Initialize the second derivative matrix M using free energy second derivatives.
                M = np.array([[d2F_dA12 *scale_xx, d2F_dA1dB2*scale_xy ,d2F_dA1dC3*scale_xz ,d2F_dA1dC2*scale_xz ,d2F_dA1dC1*scale_xz ,d2F_dA1dB1*scale_xy],
                              [d2F_dA1dB2*scale_xy, d2F_dB22  *scale_yy ,d2F_dB2dC3*scale_yz ,d2F_dB2dC2*scale_yz ,d2F_dB2dC1*scale_yz ,d2F_dB2dB1*scale_yy],
                              [d2F_dA1dC3*scale_xz, d2F_dB2dC3*scale_yz ,d2F_dC32  *scale_zz ,d2F_dC3dC2*scale_zz ,d2F_dC3dC1*scale_zz ,d2F_dC3dB1*scale_yz],
                              [d2F_dA1dC2*scale_xz, d2F_dB2dC2*scale_yz ,d2F_dC3dC2*scale_zz ,d2F_dC22  *scale_zz ,d2F_dC1dC2*scale_zz ,d2F_dB1dC2*scale_yz],
                              [d2F_dA1dC1*scale_xz, d2F_dB2dC1*scale_yz ,d2F_dC3dC1*scale_zz ,d2F_dC1dC2*scale_zz ,d2F_dC12  *scale_zz ,d2F_dB1dC1*scale_yz],
                              [d2F_dA1dB1*scale_xy, d2F_dB2dB1*scale_yy ,d2F_dC3dB1*scale_yz ,d2F_dB1dC2*scale_yz ,d2F_dB1dC1*scale_yz ,d2F_dB12  *scale_yy]])
                # Contribution of pressure to the second derivative matrix
                P = np.array([[0.0 ,pressure*v,pressure*v  ,0,0,0],
                              [pressure*v ,0.0,pressure*v  ,0,0,0],
                              [pressure*v ,pressure*v,0.0  ,0,0,0],
                              [0,0,0,0,0,0],
                              [0,0,0,0,0,0],
                              [0,0,0,0,0,0]])
                # Compute the final second derivative matrix M, including elastic and pressure contributions
                M = M + matrix_elastic + P / abu.eVA3_HaBohr3
                # Scale the first derivative of entropy
                S1 = dsda1*(exx_n+1)
                S2 = dsdb2*(eyy_n+1)
                S3 = dsdc3*(ezz_n+1)
                S4 = (dsdc2*(ezz_n+1)+dsdb2*eyz_n)
                S5 = (dsdc1*(ezz_n+1)+dsdb1*eyz_n+dsda1*exz_n)
                S6 = (dsdb1*(eyy_n+1)+dsda1*exy_n)

                dSde = np.array([S1,S2,S3,S4,S5,S6])
                # Compute thermal expansion using Eq. (37)
                dstrain_dt = np.linalg.inv(M) @ dSde
                #therm=[dstrain_dt[0]*(exx_n+1) , dstrain_dt[1]*(eyy_n+1),dstrain_dt[2]*(ezz_n+1),0,dstrain_dt[4]*(ezz_n+1)+dstrain_dt[0]*exz_n,0]
                therm = [dstrain_dt[0]*(exx_n+1), dstrain_dt[1]*(eyy_n+1),dstrain_dt[2]*(ezz_n+1),dstrain_dt[3]*(eyz_n),dstrain_dt[4]*(exz_n),dstrain_dt[5]*(exy_n)]
                elastic = M/v*abu.eVA3_GPa

        return dtol, gibbs, stress, therm, elastic

    def stress_ZSISA_slab_1DOF(self, temp: float, pressure: float) -> tuple:
        e, S = self.get_vib_free_energies(temp)

        X0 = self.ave_x[0,0,0,0,0,0]
        X1 = self.ave_x[1,0,0,0,0,0]

        dexx = (X0-X1)/X0
        exx0 = X1/X0-1

        dF_dX = (e[0,0,0,0,0,0]-e[2,0,0,0,0,0])/(2*dexx)
        d2F_dX2 = (e[0,0,0,0,0,0]-2*e[1,0,0,0,0,0]+e[2,0,0,0,0,0])/(dexx)**2

        x = self.ave_x_guess
        v = self.volume_guess

        exx_n = x/X0-1

        dfdx = dF_dX + (exx_n-exx0)*d2F_dX2
        gibbs = self.energy_guess + e[1,0,0,0,0,0] + pressure*v / abu.eVA3_HaBohr3 + \
                  (exx_n-exx0)*dF_dX + 0.5*(exx_n-exx0)**2*d2F_dX2

        stress_xx = -dfdx/v*(exx_n+1)*0.5 * abu.eVA3_HaBohr3
        if self.verbose:
            print("x/X0:", x/X0)
            print(f"{stress_xx=}")

        stress = np.zeros(6)
        stress[0] = stress_xx -pressure
        stress[1] = stress_xx -pressure

        dtol = np.zeros(6)
        dtol[0] = abs(stress[0]-self.stress_guess[0,0])
        dtol[1] = abs(stress[1]-self.stress_guess[1,1])

        return dtol, gibbs, stress

    def stress_ZSISA_slab_2DOF(self, temp: float, pressure: float) -> tuple:
        e, S = self.get_vib_free_energies(temp)

        X0 = self.ave_x[0,1,0,0,0,0]
        Y0 = self.ave_y[1,0,0,0,0,0]

        X1 = self.ave_x[1,1,0,0,0,0]
        Y1 = self.ave_y[1,1,0,0,0,0]

        dexx = (X0-X1)/X0
        deyy = (Y0-Y1)/Y0

        exx0 = X1/X0-1
        eyy0 = Y1/Y0-1

        dF_dX = (e[0,1,0,0,0,0]-e[2,1,0,0,0,0])/(2*dexx)
        dF_dY = (e[1,0,0,0,0,0]-e[1,2,0,0,0,0])/(2*deyy)

        d2F_dX2 = (e[0,1,0,0,0,0]-2*e[1,1,0,0,0,0]+e[2,1,0,0,0,0])/(dexx)**2
        d2F_dY2 = (e[1,0,0,0,0,0]-2*e[1,1,0,0,0,0]+e[1,2,0,0,0,0])/(deyy)**2
        d2F_dXdY = (e[1,1,0,0,0,0] - e[0,1,0,0,0,0] - e[1,0,0,0,0,0] + e[0,0,0,0,0,0]) / (dexx *deyy)

        x = self.ave_x_guess
        y = self.ave_y_guess
        v = self.volume_guess

        exx_n = x/X0-1
        eyy_n = y/Y0-1

        dfdx = dF_dX + (exx_n-exx0)*d2F_dX2+(eyy_n-eyy0)*d2F_dXdY
        dfdy = dF_dY + (eyy_n-eyy0)*d2F_dY2+(exx_n-exx0)*d2F_dXdY

        gibbs = self.energy_guess + e[1,1,0,0,0,0] + pressure*v / abu.eVA3_HaBohr3 + \
                  (exx_n-exx0)*dF_dX + 0.5*(exx_n-exx0)**2*d2F_dX2 +\
                  (eyy_n-eyy0)*dF_dY + 0.5*(eyy_n-eyy0)**2*d2F_dY2 +\
                  (exx_n-exx0)*(ezz_n-eyy0)*d2F_dXdY

        stress_xx = -dfdx/v*(exx_n+1)* abu.eVA3_HaBohr3
        stress_yy = -dfdy/v*(eyy_n+1)* abu.eVA3_HaBohr3
        if self.verbose:
            print("x/X0, y/Y0", x/X0, y/Y0)
            print(f"{stress_xx=}, {stress_yy=}")

        stress = np.zeros(6)
        stress[0] = stress_xx - pressure
        stress[1] = stress_yy - pressure

        dtol = np.zeros(6)
        dtol[0] = abs(stress[0]-self.stress_guess[0,0])
        dtol[1] = abs(stress[1]-self.stress_guess[1,1])

        return dtol, gibbs, stress

    def stress_ZSISA_slab_3DOF(self, temp: float, pressure: float) -> tuple:

        e, S = self.get_vib_free_energies(temp)

        Ax0 = self.ax[0,1,1,0,0,0]
        Bx0 = self.bx[0,1,0,0,0,0]
        By0 = self.by[1,0,1,0,0,0]

        Ax1 = self.ax[1,1,1,0,0,0]
        Bx1 = self.bx[1,1,1,0,0,0]
        By1 = self.by[1,1,1,0,0,0]

        v = self.volume_guess

        dexx = (Ax0-Ax1)/Ax0
        deyy = (By0-By1)/By0
        dexy = (Ax0*(Bx0-Bx1)-(Ax0-Ax1)*Bx0)/(Ax0*By0)

        exx0 = Ax1 / Ax0-1
        eyy0 = By1 / By0-1
        exy0 = (Ax0*Bx1-Ax1*Bx0) / (Ax0*By0)

        dF_dA1 = (e[0,1,1,0,0,0]-e[2,1,1,0,0,0])/(2*dexx)
        dF_dB2 = (e[1,0,1,0,0,0]-e[1,2,1,0,0,0])/(2*deyy)
        dF_dB1 = (e[1,1,0,0,0,0]-e[1,1,2,0,0,0])/(2*dexy)

        d2F_dA12 = (e[0,1,1,0,0,0]-2*e[1,1,1,0,0,0]+e[2,1,1,0,0,0])/(dexx)**2
        d2F_dB22 = (e[1,0,1,0,0,0]-2*e[1,1,1,0,0,0]+e[1,2,1,0,0,0])/(deyy)**2
        d2F_dB12 = (e[1,1,0,0,0,0]-2*e[1,1,1,0,0,0]+e[1,1,2,0,0,0])/(dexy)**2

        d2F_dA1dB2 = (e[1,1,1,0,0,0] - e[0,1,1,0,0,0] - e[1,0,1,0,0,0] + e[0,0,1,0,0,0]) / (dexx *deyy)
        d2F_dA1dB1 = (e[1,1,1,0,0,0] - e[0,1,1,0,0,0] - e[1,1,0,0,0,0] + e[0,1,0,0,0,0]) / (dexx *dexy)

        d2F_dB2dB1 = (e[1,1,1,0,0,0] - e[1,0,1,0,0,0] - e[1,1,0,0,0,0] + e[1,0,0,0,0,0]) / (deyy *dexy)

        a = self.lattice_a_guess
        b = self.lattice_b_guess
        c = self.lattice_c_guess

        cos_ab = math.cos(math.pi*self.angles_guess[2]/180)

        ax = a
        bx = b * cos_ab
        by = b * np.sqrt(1-cos_ab**2)
        cz = c

        exx_n = ax/Ax0-1
        eyy_n = by/By0-1
        exy_n = (Ax0*bx-ax*Bx0)/(Ax0*By0)

        dfda1 = dF_dA1 + (exx_n-exx0)*d2F_dA12+(eyy_n-eyy0)*d2F_dA1dB2+(exy_n-exy0)*d2F_dA1dB1
        dfdb2 = dF_dB2 + (eyy_n-eyy0)*d2F_dB22+(exx_n-exx0)*d2F_dA1dB2+(exy_n-exy0)*d2F_dB2dB1
        dfdb1 = dF_dB1 + (exy_n-exy0)*d2F_dB12+(exx_n-exx0)*d2F_dA1dB1+(eyy_n-eyy0)*d2F_dB2dB1
        gibbs = self.energy_guess + e[1,1,1,0,0,0] + pressure*v / abu.eVA3_HaBohr3 + \
              (exx_n-exx0)*dF_dA1 + 0.5*(exx_n-exx0)**2*d2F_dA12 +\
              (eyy_n-eyy0)*dF_dB2 + 0.5*(eyy_n-eyy0)**2*d2F_dB22 +\
              (exy_n-exy0)*dF_dB1 + 0.5*(exy_n-exy0)**2*d2F_dB12 +\
              (exx_n-exx0)*(eyy_n-eyy0)*d2F_dA1dB2 +\
              (exx_n-exx0)*(exy_n-exy0)*d2F_dA1dB1 +\
              (eyy_n-eyy0)*(exy_n-exy0)*d2F_dB2dB1

        stress_a1 = -dfda1/v*(exx_n+1) * abu.eVA3_HaBohr3
        stress_b2 = -dfdb2/v*(eyy_n+1) * abu.eVA3_HaBohr3
        stress_b1 = -1.0/v*(dfdb1*(eyy_n+1)+dfda1*exy_n) * abu.eVA3_HaBohr3
        if self.verbose:
            print("ax/Ax0, by/By0", ax/Ax0, by/By0)
            print(f"{stress_a1=}, {stress_b2=}, {stress_b1=}")

        stress = np.zeros(6)
        stress[0] = stress_a1 - pressure
        stress[1] = stress_b2 - pressure
        stress[5] = stress_b1

        dtol = np.zeros(6)
        dtol[0] = abs(stress[0] - self.stress_guess[0,0])
        dtol[1] = abs(stress[1] - self.stress_guess[1,1])
        dtol[5] = abs(stress[5] - self.stress_guess[1,0])

        return dtol, gibbs, stress

    def get_tstress(self,
                    temp: float,
                    pressure_gpa: float,
                    structure_guess: Structure,
                    stress_guess,
                    energy_guess,
                    mode: str = "TEC",
                    bo_elastic_voigt=None) -> ThermalData:
                    #elastic_path: str | None = "elastic_constant.txt") -> ThermalData:
        """
        Args
            temp: Temperature in K.
            pressure_gpa: Pressure in GPa
            structure_guess: Structure used as a guess.
            stress_guess: Stress tensor corresponding to the initial guess structure as (3,3) matrix in a.u.
            mode: "TEC" or "ECs"
            bo_elastic_voigt: Elastic tensor in Voigt notation i.e. shape = (6,6) with BO elastic constants in GPa.
        """
        self._set_structure_stress_guess(structure_guess, stress_guess, energy_guess)

        self.bo_elastic_voigt = bo_elastic_voigt
        pressure_au = pressure_gpa / abu.HaBohr3_GPa
        if self.verbose:
            print("Mode: ", mode, "\nPressure: ", pressure_gpa, "(GPa)", "\nTemperature: ", temp, "(K)")

        elastic, therm = None, None

        if self.sym == "v_ZSISA":
            dtol, gibbs, stress = self.stress_v_ZSISA(temp, pressure_au)

        elif self.sym == "cubic" and mode == "TEC":
            dtol, gibbs, stress, therm = self.stress_ZSISA_1DOF(temp, pressure_au)

        elif self.sym in ("trigonal", "hexagonal", "tetragonal") and mode == "TEC":
            dtol, gibbs, stress, therm = self.stress_ZSISA_2DOF(temp, pressure_au)

        elif self.sym in ("cubic", "trigonal", "hexagonal", "tetragonal") and mode == "ECs":
            dtol, gibbs, stress, therm, elastic = self.stress_ZSISA_3DOF(temp, pressure_au, mode)

        elif self.sym == "orthorhombic":
            dtol, gibbs, stress, therm, elastic = self.stress_ZSISA_3DOF(temp, pressure_au, mode)

        elif self.sym == "monoclinic":
            dtol, gibbs, stress, therm, elastic = self.stress_ZSISA_monoclinic(temp, pressure_au, mode)

        elif self.sym == "triclinic":
            dtol, gibbs, stress, therm, elastic = self.stress_ZSISA_triclinic(temp, pressure_au, mode)

        elif self.sym == "ZSISA_slab_1DOF":
            dtol, gibbs, stress = self.stress_ZSISA_slab_1DOF(temp, pressure_au)

        elif self.sym == "ZSISA_slab_2DOF":
            dtol, gibbs, stress = self.stress_ZSISA_slab_2DOF(temp, pressure_au)

        elif self.sym == "ZSISA_slab_3DOF":
            dtol, gibbs, stress = self.stress_ZSISA_slab_3DOF(temp, pressure_au)

        else:
            raise ValueError(f"Unknown {self.sym=}")

        converged = False
        if all(dtol[i] < self.dtol_tolerance for i in range(6)):
            converged = True
            #self.print_data(temp, pressure_gpa, therm, stress, elastic, mode)
            if self.verbose: print("Converged !!!")

        return ThermalData(temperature=temp, pressure_gpa=pressure_gpa,
                           converged=converged, dtol=dtol, stress_au=stress, elastic=elastic, therm=therm, gibbs=gibbs)

    def get_vib_free_energies(self, temp: float) -> tuple:
        """
        Compute vibrational free energy and entropy at the specified temperature `temp` in K.

        Return:
            e: Vibrational free energy (F_vib)
            S: Entropy.
        """
        f = np.zeros((self.dim[0],self.dim[1],self.dim[2],self.dim[3],self.dim[4],self.dim[5]))
        entropy = np.zeros((self.dim[0],self.dim[1],self.dim[2],self.dim[3],self.dim[4],self.dim[5]))

        for i, dim0_list in enumerate(self.phdoses):
            for j, dim1_list in enumerate(dim0_list):
                for k, dim2_list in enumerate(dim1_list):
                    for l, dim3_list in enumerate(dim2_list):
                        for m, dim4_list in enumerate(dim3_list):
                            for n, dos in enumerate(dim4_list):
                               if dos is not None:
                                   f[i,j,k,l,m,n] = dos.get_free_energy(temp, temp, 1).values.item()
                                   entropy[i,j,k,l,m,n] = dos.get_entropy(temp, temp, 1).values.item()
                               else:
                                   f[i,j,k,l,m,n] = None
                                   entropy[i,j,k,l,m,n] = None
        return f, entropy

    #def elastic_constants(self, file_name: str) -> np.ndarray:
    #    print("About to read ELASTIC_RELAXED tensor from:", file_name)
    #    with open(file_name, "rt") as file:
    #        lines = file.readlines()

    #    in_relaxed_section = False
    #    elastic_data = []

    #    for line in lines:
    #        if "[ELASTIC_RELAXED]" in line:
    #            in_relaxed_section = True
    #            continue
    #        if "[" in line and in_relaxed_section:  # Stop if another section starts
    #            break
    #        if in_relaxed_section:
    #            try:
    #                # Extract numbers from each line
    #                numbers = list(map(float, line.split()[1:7]))
    #                #print("numbers:", numbers)
    #                if numbers:  # Add non-empty lines
    #                    elastic_data.append(numbers)
    #            except ValueError as exc:
    #                #raise exc
    #                continue

    #    # Convert the list into a numpy array for easier indexing
    #    return np.array(elastic_data)

    #def print_data(self, temp, pressure, therm, stress, elastic, mode):
    #    M = elastic
    #    filename = f"cell_{temp:04.0f}_{pressure:03.0f}.txt"
    #    with open(filename, "w") as f:
    #        f.write(f"{'#T':<8} {'P':<8} {'lattice_a':<13} {'lattice_b':<13} {'lattice_c':<13} "
    #                f"{'alpha':<10} {'beta':<10} {'gamma':<10} {'volume':<13} {'ave_x':<13} "
    #                f"{'ave_y':<13} {'ave_z':<13}\n")

    #        f.write(f"{temp:<8} {pressure:<8.2f} {self.lattice_a_guess:<13.10f} {self.lattice_b_guess:<13.10f} "
    #                f"{self.lattice_c_guess:<13.10f} {self.angles_guess[0]:<10.5f} {self.angles_guess[1]:<10.5f} "
    #                f"{self.angles_guess[2]:<10.5f} {self.volume_guess:<13.10f} {self.ave_x_guess:<13.10f} "
    #                f"{self.ave_y_guess:<13.10f} {self.ave_z_guess:<13.10f} \n")

    #    if therm is not None:
    #        filename = f"TEC_{temp:04.0f}_{pressure:03.0f}.txt"
    #        with open(filename, "w") as f:
    #            f.write(f"{'#T':<8} {'P':<8} {'alpha_xx':<15} {'alpha_yy':<15} {'alpha_zz':<15} {'alpha_yz':<15} {'alpha_xz':<15} {'alpha_xy':<15}\n")
    #            f.write(f"{temp:<8} {pressure:<8.2f} {therm[0]:<15.8e} {therm[1]:<15.8e} {therm[2]:<15.8e} "
    #                    f"{therm[3]:<15.8e} {therm[4]:<15.8e} {therm[5]:<15.8e}\n")

    #    if elastic is not None:
    #        filename = f"ECs_{temp:04.0f}_{pressure:03.0f}.txt"
    #        with open(filename, "w") as f:
    #            f.write("Elastic [GPa] \n")
    #            if self.sym in ("cubic", "trigonal", "hexagonal", "tetragonal", "orthorhombic"):
    #                if mode == 'ECs':
    #                    if self.sym == "cubic":
    #                        f.write(f"    {'C11':12s} {'C12':12s} {'C44':12s}\n")
    #                        f.write(f"  {M[0,0]:12.6f} {M[0,1]:12.6f} {M[3,3]:12.6f} \n")
    #                    elif self.sym == "hexagonal":
    #                        f.write(f"    {'C11':12s} {'C12':12s} {'C13':12s} {'C33':12s} {'C44':12s}\n")
    #                        f.write(f"  {M[0,0]:12.6f} {M[0,1]:12.6f} {M[0,2]:12.6f} {M[2,2]:12.6f} {M[3,3]:12.6f} \n")
    #                    elif self.sym == "trigonal":
    #                        f.write(f"    {'C11':12s} {'C12':12s} {'C13':12s} {'C33':12s} {'C14':12s} {'C44':12s}\n")
    #                        f.write(f"  {M[0,0]:12.6f} {M[0,1]:12.6f} {M[0,2]:12.6f} {M[2,2]:12.6f} {M[0,3]:12.6f} {M[3,3]:12.6f} \n")
    #                    elif self.sym == "tetragonal":
    #                        f.write(f"    {'C11':12s} {'C12':12s} {'C13':12s} {'C33':12s} {'C44':12s} {'C66':12s}\n")
    #                        f.write(f"  {M[0,0]:12.6f} {M[0,1]:12.6f} {M[0,2]:12.6f} {M[2,2]:12.6f} {M[3,3]:12.6f} {M[5,5]:12.6f} \n")
    #                    if  self.sym == "orthorhombic":
    #                        f.write(f"    {'C11':12s} {'C12':12s} {'C13':12s} {'C22':12s} {'C23':12s} {'C33':12s} {'C44':12s} {'C55':12s} {'C66':12s}\n")
    #                        f.write(f"  {M[0,0]:12.6f} {M[0,1]:12.6f} {M[0,2]:12.6f} {M[1,1]:12.6f} {M[1,2]:12.6f} {M[2,2]:12.6f} {M[3,3]:12.6f} {M[4,4]:12.6f} {M[5,5]:12.6f} \n")
    #                else:
    #                    f.write(f"    {'C11':12s} {'C12':12s} {'C13':12s} {'C22':12s} {'C23':12s} {'C33':12s} \n")
    #                    f.write(f"  {M[0,0]:12.6f} {M[0,1]:12.6f} {M[0,2]:12.6f} {M[1,1]:12.6f} {M[1,2]:12.6f} {M[2,2]:12.6f} \n")

    #            elif self.sym == "monoclinic":
    #                if mode != 'ECs':
    #                    f.write(f" Warning: C44, C46, and C66 do not include the free energy contribution (only BO energy).\n")
    #                f.write( " \t   xx\t\tyy\t\tzz\t\tyz\t\txz\t\txy\n")
    #                f.write(f" xx {M[0,0]:14.8f}  {M[0,1]:14.8f}  {M[0,2]:14.8f}  {M[0,3]:14.8f}  {M[0,4]:14.8f}  {M[0,5]:14.8f}\n")
    #                f.write(f" yy {M[1,0]:14.8f}  {M[1,1]:14.8f}  {M[1,2]:14.8f}  {M[1,3]:14.8f}  {M[1,4]:14.8f}  {M[1,5]:14.8f}\n")
    #                f.write(f" zz {M[2,0]:14.8f}  {M[2,1]:14.8f}  {M[2,2]:14.8f}  {M[2,3]:14.8f}  {M[2,4]:14.8f}  {M[2,5]:14.8f}\n")
    #                f.write(f" yz {M[3,0]:14.8f}  {M[3,1]:14.8f}  {M[3,2]:14.8f}  {M[3,3]:14.8f}  {M[3,4]:14.8f}  {M[3,5]:14.8f}\n")
    #                f.write(f" xz {M[4,0]:14.8f}  {M[4,1]:14.8f}  {M[4,2]:14.8f}  {M[4,3]:14.8f}  {M[4,4]:14.8f}  {M[4,5]:14.8f}\n")
    #                f.write(f" xy {M[5,0]:14.8f}  {M[5,1]:14.8f}  {M[5,2]:14.8f}  {M[5,3]:14.8f}  {M[5,4]:14.8f}  {M[5,5]:14.8f}\n")

    #            elif self.sym == "triclinic":
    #                f.write(f" \t   xx\t\tyy\t\tzz\t\tyz\t\txz\t\txy\n")
    #                f.write(f" xx {M[0,0]:14.8f}  {M[0,1]:14.8f}  {M[0,2]:14.8f}  {M[0,3]:14.8f}  {M[0,4]:14.8f}  {M[0,5]:14.8f}\n")
    #                f.write(f" yy {M[1,0]:14.8f}  {M[1,1]:14.8f}  {M[1,2]:14.8f}  {M[1,3]:14.8f}  {M[1,4]:14.8f}  {M[1,5]:14.8f}\n")
    #                f.write(f" zz {M[2,0]:14.8f}  {M[2,1]:14.8f}  {M[2,2]:14.8f}  {M[2,3]:14.8f}  {M[2,4]:14.8f}  {M[2,5]:14.8f}\n")
    #                f.write(f" yz {M[3,0]:14.8f}  {M[3,1]:14.8f}  {M[3,2]:14.8f}  {M[3,3]:14.8f}  {M[3,4]:14.8f}  {M[3,5]:14.8f}\n")
    #                f.write(f" xz {M[4,0]:14.8f}  {M[4,1]:14.8f}  {M[4,2]:14.8f}  {M[4,3]:14.8f}  {M[4,4]:14.8f}  {M[4,5]:14.8f}\n")
    #                f.write(f" xy {M[5,0]:14.8f}  {M[5,1]:14.8f}  {M[5,2]:14.8f}  {M[5,3]:14.8f}  {M[5,4]:14.8f}  {M[5,5]:14.8f}\n")


    def get_phdos_plotter(self, **kwargs) -> PhononDosPlotter:
        """Build and return a PhononDosPlotter."""
        plotter = PhononDosPlotter()
        #print(f"{self.phdoses.shape=}, {self.dim=}")
        for inds, phdos in np.ndenumerate(self.phdoses):
            if phdos is None: continue
            plotter.add_phdos(str(inds), phdos)

        return plotter


def cmat_inds_names(sym: str, mode: str) -> tuple[list, list]:
    """
    Return list with the numpy indices of the non-null components of
    the elastic tensor as well as list with their names e.g. (0, 0) -> "C_11".

    Args:
        sym: Crystalline system.
        mode:
    """
    if sym in ("cubic", "trigonal", "hexagonal", "tetragonal", "orthorhombic"):
        if mode == 'ECs':
            if sym == "cubic":
                inds_list = [(0,0), (0,1), (3,3)]
            elif sym == "hexagonal":
                inds_list = [(0,0), (0,1), (0,2), (2,2), (3,3)]
            elif sym == "trigonal":
                inds_list = [(0,0), (0,1), (0,2), (2,2), (0,3), (3,3)]
            elif sym == "tetragonal":
                inds_list = [(0,0), (0,1), (0,2), (2,2), (3,3), (5,5)]
            if  sym == "orthorhombic":
                inds_list = [(0,0), (0,1), (0,2), (1,1), (1,2), (2,2), (3,3), (4,4), (5,5)]

        elif mode == 'TEC':
            inds_list = [(0,0), (0,1), (0,2), (1,1), (1,2), (2,2)]

        else:
            raise ValueError(f"Invalid {mode=}")

    elif sym == "monoclinic":
        if mode != 'ECs':
            print("Warning: C44, C46, and C66 do not include the free energy contribution (only BO energy).")
        inds_list = [
            (0,0), (0,1), (0,2), (0,3), (0,4), (0,5),
            (1,0), (1,1), (1,2), (1,3), (1,4), (1,5),
            (2,0), (2,1), (2,2), (2,3), (2,4), (2,5),
            (3,0), (3,1), (3,2), (3,3), (3,4), (3,5),
            (4,0), (4,1), (4,2), (4,3), (4,4), (4,5),
            (5,0), (5,1), (5,2), (5,3), (5,4), (5,5),
        ]

    elif sym == "triclinic":
        inds_list = [
            (0,0), (0,1), (0,2), (0,3), (0,4), (0,5),
            (1,0), (1,1), (1,2), (1,3), (1,4), (1,5),
            (2,0), (2,1), (2,2), (2,3), (2,4), (2,5),
            (3,0), (3,1), (3,2), (3,3), (3,4), (3,5),
            (4,0), (4,1), (4,2), (4,3), (4,4), (4,5),
            (5,0), (5,1), (5,2), (5,3), (5,4), (5,5),
        ]

    # Build names. Note +1
    names = [f"C_{inds[0]+1}{inds[1]+1}" for inds in inds_list]

    return inds_list, names
