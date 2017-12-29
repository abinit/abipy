# coding: utf-8
"""HirshfeldCharges."""
from __future__ import print_function, division, unicode_literals, absolute_import

from abipy.core.mixins import Has_Structure
from abipy.core.fields import Density
from abipy.electrons.denpot import DensityFortranFile
from pymatgen.command_line.bader_caller import BaderAnalysis
from pymatgen.io.abinit.pseudos import Pseudo
from pymatgen.core.units import bohr_to_angstrom
from monty.dev import requires
from monty.os.path import which

import numpy as np
import os
import tempfile
import logging
logger = logging.getLogger(__name__)


__all__ = [
    "HirshfeldCharges",
    "BaderCharges"
]


class Charges(Has_Structure):
    """
    Base charges class, describing the atomic charges in a DFT calculation.
    The charges refer to the total charges for each atom (i.e. Z-n_electrons).
    An eccess of electron has a negative charge

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: Charges
    """

    def __init__(self, electron_charges, structure, reference_charges=None):
        """
        Args:
            electron_charges: Charges coming from the electron for each element of the structure (with negative sign)
            structure: |Structure| object.
            reference_charges: Reference charges associated to each atom of the structure, considered as isolated
                (with negative sign). Should represent either all the electrons or just the valence, depending on the
                charge density used to perform the analysis.
        """
        self.electron_charges = np.array(electron_charges)
        self._structure = structure
        self.reference_charges = np.array(reference_charges)

    @property
    def structure(self):
        return self._structure

    @property
    def net_charges(self):
        """
        List of net charges for each atom of the structure.
        """
        if self.reference_charges is None:
            raise ValueError("reference charges are required to calculate the net transfer")

        return self.electron_charges - self.reference_charges


class HirshfeldCharges(Charges):
    """
    Class representing the charges obtained from the Hirshfeld analysis.
    The charges refer to the total charges for each atom (i.e. Z-n_electrons).
    An eccess of electron has a negative charge

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: HirshfeldCharges
    """
    @classmethod
    def from_cut3d_outfile(cls, filepath, structure):
        """
        Generates a HirshfeldCharges object from the outputfile of cut3d and a structure.
        """

        electron_charges = []
        reference_charges = []
        with open(filepath, 'rt') as f:
            lines = f.readlines()

        start_hirshfeld_i = None
        for i, l in enumerate(lines):
            if "Hirshfeld analysis" in l:
                start_hirshfeld_i = i+3
                break
        else:
            raise RuntimeError('The file does not contain Hirshfeld charges')


        for i in range(start_hirshfeld_i, start_hirshfeld_i+len(structure)):
            l = lines[i]
            electron_charges.append(float(l.split()[2]))
            reference_charges.append(-float(l.split()[1]))

        return cls(electron_charges, structure, reference_charges)


class BaderCharges(Charges):
    """
    Class representing the charges obtained from the Bader analysis.
    TThe charges refer to the total charges for each atom (i.e. Z-n_electrons).
    An eccess of electron has a negative charge

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: BaderCharges
    """

    @classmethod
    @requires(which("bader") or which("bader.exe"),
              "BaderCharges.from_files requires the executable bader to be in the path.")
    def from_files(cls, density_path, pseudopotential_paths, with_core=True, workdir=None, **kwargs):
        """
        Uses the abinit density files and the bader_ executable from Henkelmann et al. to calculate
        the bader charges of the system. If pseudopotentials are given, the atomic charges will be
        extracted as well. See also :cite:`Henkelman2006`.

        The extraction of the core charges may be a time consuming calculation, depending on the
        size of the system. A tuning of the parameters may be required (see Density.ae_core_density_on_mesh).

        Args:
            density_path: Path to the abinit density file. Can be a fortran _DEN or a netCDF DEN.nc file.
                In case of fortran file, requires cut3d (version >= 8.6.1) for the convertion.
            pseudopotential_paths: Dictionary {element: pseudopotential path} for all the elements present
                in the system.
            with_core: Core charges will be extracted from the pseudopotentials with the
                Density.ae_core_density_on_mesh method. Requires pseudopotential_paths.
            workdir: Working directory. If None, a temporary directory is created.
            kwargs: arguments passed to the method ``Density.ae_core_density_on_mesh``

        Returns:
            An instance of :class:`BaderCharges`
        """
        # read the valence density
        # if density is not a netcdf file, convert with cut3d
        if not density_path.endswith('.nc'):
            dff = DensityFortranFile.from_file(density_path)
            density = dff.get_density()
        else:
            density = Density.from_file(density_path)

        structure = density.structure

        atomic_charges = None

        if with_core:
            if not pseudopotential_paths:
                raise ValueError("pseudopotentials should be provided to extract the core densities")

            try:
                from pseudo_dojo.ppcodes.oncvpsp import psp8_get_densities
            except ImportError as exc:
                print("PseudoDojo package required to extract core densities. "
                      "Please install it with `pip install pseudo_dojo`")
                raise exc

            # extract core charge from pseudopotentials on a radial grid in the correct units
            rhoc = {}
            for specie, ppath in pseudopotential_paths.items():
                r = psp8_get_densities(ppath)
                rhoc[specie] = [r.rmesh * bohr_to_angstrom, r.aecore / (4.0 * np.pi) / (bohr_to_angstrom ** 3)]

            workdir = tempfile.mkdtemp() if workdir is None else workdir

            # extrapolate the core density on the density grid
            core_density = Density.ae_core_density_on_mesh(density, structure, rhoc, **kwargs)
            density += core_density
            atomic_charges = [s.specie.Z for s in structure]

        elif pseudopotential_paths:
            pseudos = {k: Pseudo.from_file(p) for k, p in pseudopotential_paths.items()}
            atomic_charges = [pseudos[s.specie.name].Z_val for s in structure]

        chgcar_path = os.path.join(workdir, 'CHGCAR')
        density.to_chgcar(chgcar_path)

        ba = BaderAnalysis(chgcar_path)

        charges = [ba.get_charge(i) for i in range(len(structure))]

        return cls(charges, structure, atomic_charges)