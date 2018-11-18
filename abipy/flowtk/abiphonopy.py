# coding: utf-8
"""Interface between phonopy and abipy workflow model."""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import numpy as np

try:
    from phonopy import Phonopy, file_IO
    from phonopy.interface.vasp import read_vasp_from_strings
    from phonopy.interface.abinit import parse_set_of_forces
except ImportError:
    import warnings
    warnings.warn("phonopy is required by abiphonopy. Install it with conda or pip install phonopy")

from abipy.core.structure import Structure
from abipy.flowtk.works import Work


__all__ = [
    "atoms_from_structure",
    "structure_from_atoms",
    "PhonopyWork",
    "PhonopyGruneisenWork",
]


def atoms_from_structure(structure):
    """
    Convert a pymatgen Structure into a phonopy Atoms object.
    """
    s = structure.to(fmt="poscar", filename=None)
    return read_vasp_from_strings(s, symbols=None)


def structure_from_atoms(atoms):
    """
    Convert a phonopy Atoms object into a pymatgen Structure.
    """
    return Structure(lattice=atoms.cell,
                     species=atoms.symbols,
                     coords=atoms.scaled_positions,
                     validate_proximity=False, to_unit_cell=False,
                     coords_are_cartesian=False, site_properties=None)


class PhonopyWork(Work):
    """
    This work computes the inter-atomic force constants with phonopy.

    .. attribute:: scdims(3)

	numpy arrays with the number of cells in the supercell along the three reduced directions

    .. attribute:: phonon

	:class:`Phonopy` object used to construct the supercells with displaced atoms.

    .. attribute:: phonopy_tasks

	List of :class:`ScfTask`. Each task compute the forces in one perturbed supercell.

    .. attribute:: bec_tasks

    .. attribute:: cpdata2dst

	If not None, the work will copy the output results to the outdir of the flow
	once all_ok is invoked. Note that cpdata2dst must be an absolute path.
    """

    @classmethod
    def from_gs_input(cls, gsinp, scdims, phonopy_kwargs=None, displ_kwargs=None):
        """
        Build the work from an :class:`AbinitInput` object representing a GS calculations.

	Args:
	    gsinp:
		:class:`AbinitInput` object representing a GS calculation in the initial unit cell.
	    scdims:
		Number of unit cell replicas along the three reduced directions.
	    phonopy_kwargs:
		(Optional) dictionary with arguments passed to Phonopy constructor.
	    displ_kwargs:
		(Optional) dictionary with arguments passed to generate_displacements.

	Return:
	    PhonopyWork instance.
        """
        new = cls()
        new.phonopy_tasks, new.bec_tasks = [], []
        new.cpdata2dst = None

        # Initialize phonon. Supercell matrix has (3, 3) shape.
        unitcell = atoms_from_structure(gsinp.structure)
        new.scdims = scdims = np.array(scdims)
        if scdims.shape != (3,):
            raise ValueError("Expecting 3 int in scdims but got %s" % str(scdims))

        supercell_matrix = np.diag(scdims)
        phonopy_kwargs = phonopy_kwargs if phonopy_kwargs is not None else {}
        new.phonon = phonon = Phonopy(unitcell, supercell_matrix, **phonopy_kwargs)

        displ_kwargs = displ_kwargs if displ_kwargs is not None else {}
        phonon.generate_displacements(**displ_kwargs)  # distance=0.01,

        # Obtain supercells containing respective displacements (list of Atoms objects).
        for atoms in phonon.get_supercells_with_displacements():
            sc_struct = structure_from_atoms(atoms)
            sc_gsinp = gsinp.new_with_structure(sc_struct, scdims=new.scdims)
            sc_gsinp.pop_tolerances()
            sc_gsinp.pop_vars(["ionmov", "optcell", "ntime"])
            sc_gsinp.set_vars(toldff=1.e-6)
            sc_gsinp.set_vars_ifnotin(chksymbreak=0, chkprim=0)

            task = new.register_scf_task(sc_gsinp)
            new.phonopy_tasks.append(task)

        return new

    def on_all_ok(self):
        """
        This method is called once the `Work` is completed i.e. when all the tasks
        have reached status S_OK. Here we get the forces from the output files, and
        we call phonopy to compute the inter-atomic force constants.
        """
        phonon = self.phonon

	# Write POSCAR with initial unit cell.
        structure = structure_from_atoms(phonon.get_primitive())
        structure.to(filename=self.outdir.path_in("POSCAR"))

	# Write yaml file with displacements.
        supercell = phonon.get_supercell()
        displacements = phonon.get_displacements()
        #directions = phonon.get_displacement_directions()
        file_IO.write_disp_yaml(displacements, supercell, # directions=directions,
                                filename=self.outdir.path_in('disp.yaml'))

        # Extract forces from the main Abinit output files.
        forces_filenames = [task.output_file.path for task in self.phonopy_tasks]
        num_atoms = supercell.get_number_of_atoms()
        force_sets = parse_set_of_forces(num_atoms, forces_filenames)

	# Write FORCE_SETS file.
        displacements = file_IO.parse_disp_yaml(filename=self.outdir.path_in('disp.yaml'))
        num_atoms = displacements['natom']
        for forces, disp in zip(force_sets, displacements['first_atoms']):
            disp['forces'] = forces
        file_IO.write_FORCE_SETS(displacements, filename=self.outdir.path_in('FORCE_SETS'))

        # Write README and configuration files.
        examples_url = "http://atztogo.github.io/phonopy/examples.html"
        doctags_url = "http://atztogo.github.io/phonopy/setting-tags.html#setting-tags"
        kptbounds = np.array([k.frac_coords for k in structure.hsym_kpoints])
        path_coords = " ".join(str(rc) for rc in kptbounds.flat)
        path_labels = " ".join(k.name for k in structure.hsym_kpoints)
        ngqpt = structure.calc_ngkpt(nksmall=30)

        with open(self.outdir.path_in("band.conf"), "wt") as fh:
            fh.write("#" + doctags_url + "\n")
            fh.write("DIM = %d %d %d\n" % tuple(self.scdims))
            fh.write("BAND = %s\n" % path_coords)
            fh.write("BAND_LABELS = %s\n" % path_labels)
            fh.write("BAND_POINTS = 101\n")
            fh.write("#BAND_CONNECTION = .TRUE.\n")

        with open(self.outdir.path_in("dos.conf"), "wt") as fh:
            fh.write("#" + doctags_url + "\n")
            fh.write("DIM = %d %d %d\n" % tuple(self.scdims))
            fh.write("MP = %d %d %d\n" % tuple(ngqpt))
            fh.write("#GAMMA_CENTER = .TRUE.\n")

        with open(self.outdir.path_in("band-dos.conf"), "wt") as fh:
            fh.write("#" + doctags_url + "\n")
            fh.write("DIM = %d %d %d\n" % tuple(self.scdims))
            fh.write("BAND = %s\n" % path_coords)
            fh.write("BAND_LABELS = %s\n" % path_labels)
            fh.write("BAND_POINTS = 101\n")
            fh.write("#BAND_CONNECTION = .TRUE.\n")
            fh.write("MP = %d %d %d\n" % tuple(ngqpt))
            fh.write("#GAMMA_CENTER = .TRUE.\n")

        with open(self.outdir.path_in("README.md"), "wt") as fh:
            fh.write("To plot bands, use:\n\tphonopy -p band.conf\n\n")
            fh.write("To plot phonon dos, use:\n\tphonopy -p dos.conf\n\n")
            fh.write("To plot bands and dos, use:\n\tphonopy -p band-dos.conf\n\n")
            fh.write("See also:\n")
            fh.write("\t" + examples_url + "\n")
            fh.write("\t" + doctags_url + "\n")

        if self.cpdata2dst:
            self.outdir.copy_r(self.cpdata2dst)

        return super(PhonopyWork, self).on_all_ok()


class PhonopyGruneisenWork(Work):
    """
    This work computes the Grüneisen parameters with phonopy. The workflow is as follows:

    It is necessary to run three phonon calculations.
    One is calculated at the equilibrium volume and the remaining two are calculated
    at the slightly larger volume and smaller volume than the equilibrium volume.
    The unitcells at these volumes have to be fully relaxed under the constraint of each volume.

    .. attribute:: scdims(3)

	numpy arrays with the number of cells in the supercell along the three reduced directions.

    """
    @classmethod
    def from_gs_input(cls, gsinp, voldelta, scdims, phonopy_kwargs=None, displ_kwargs=None):
        """
        Build the work from an :class:`AbinitInput` object representing a GS calculations.

	Args:
	    gsinp:
		:class:`AbinitInput` object representing a GS calculation in the initial unit cell.
	    voldelta:
                Absolute increment for unit cell volume. The three volumes are:
                     [v0 - voldelta, v0, v0 + voldelta] where v0 is taken from gsinp.structure.
	    scdims:
		Number of unit cell replicas along the three reduced directions.
	    phonopy_kwargs:
		(Optional) dictionary with arguments passed to Phonopy constructor.
	    displ_kwargs:
		(Optional) dictionary with arguments passed to generate_displacements.

	Return:
	    `PhonopyGruneisenWork` instance.
        """
        new = cls()

	# Save arguments that will be used to call phonopy for creating
        # the supercells with the displacements once the three volumes have been relaxed.
        new.scdims = np.array(scdims)
        if new.scdims.shape != (3,):
            raise ValueError("Expecting 3 int in scdims but got %s" % str(new.scdims))
        new.phonopy_kwargs = phonopy_kwargs if phonopy_kwargs is not None else {}
        new.displ_kwargs = displ_kwargs if displ_kwargs is not None else {}

        # Build three tasks for structural optimization at constant volume.
        v0 = gsinp.structure.volume
        if voldelta <= 0:
            raise ValueError("voldelta must be > 0 but got %s" % voldelta)
        volumes = [v0 - voldelta, v0, v0  + voldelta]
        if any(v <= 0 for v in volumes):
            raise ValueError("volumes must be > 0 but got %s" % str(volumes))

        for vol in volumes:
            # Build new structure
            new_lattice = gsinp.structure.lattice.scale(vol)
            new_structure = Structure(new_lattice, gsinp.structure.species, gsinp.structure.frac_coords)
            new_input = gsinp.new_with_structure(new_structure)
            # Set variables for structural optimization at constant volume.
            new_input.pop_tolerances()
            new_input.set_vars(optcell=3, ionmov=3, tolvrs=1e-10, toldff=1.e-6)
            new_input.set_vars_ifnotin(ecutsm=0.5, dilatmx=1.05)
            new.register_relax_task(new_input)

        return new

    def on_all_ok(self):
        """
        This method is called once the `Work` is completed i.e. when all the tasks
        have reached status S_OK.
        """
        self.add_phonopy_works_and_build()
        return super(PhonopyGruneisenWork, self).on_all_ok()

    def add_phonopy_works_and_build(self):
        """
        Get relaxed structures from the tasks, build Phonopy works with supercells
        constructed from the new structures, add them to the flow and build new directories.
        """
        for i, task in enumerate(self):
            relaxed_structure = task.get_final_structure()
            gsinp = task.input.new_with_structure(relaxed_structure)

            work = PhonopyWork.from_gs_input(gsinp, self.scdims,
                                             phonopy_kwargs=self.phonopy_kwargs,
                                             displ_kwargs=self.displ_kwargs)

            self.flow.register_work(work)
            # Tell the work to copy the results to e.g. `flow/outdir/w0/minus`
            dst = os.path.join(self.pos_str, {0: "minus", 1: "orig", 2: "plus"}[i])
            work.cpdata2dst = self.flow.outdir.path_in(dst)

        # Allocate new works and update the pickle database.
        self.flow.allocate()
        self.flow.build_and_pickle_dump()
