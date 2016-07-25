# coding: utf-8
"""
Interface betwee phonopy and abipy workflow model.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np

from phonopy import Phonopy, file_IO
from phonopy.structure.atoms import Atoms as PhonopyAtoms
from phonopy.interface.vasp import read_vasp_from_strings
from phonopy.interface.abinit import parse_set_of_forces
from pymatgen.io.abinit.works import Work
from abipy.core.structure import Structure
from abipy.abio.inputs import AbinitInput


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
    This work compute the inter-atomic force constants with phonopy.

    .. attribute:: scdims(3)

	numpy arrays with the number of cells in the supercell along the three reduced directions

    .. attribute:: phonon

	:class:`Phonopy` object used to construct the supercells with displaced atoms.

    .. attribute:: phonopy_tasks

	List of :class:`ScfTask`. Each task compute the forces in one perturbed supercell.

    .. attribute:: bec_tasks
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
		Dictionary with arguments passed to Phonopy constructor.
	    displ_kwargs:
		Dictionary with arguments passed to generate_displacements.

	Return:
	    PhonopyWork instance.
        """
        new = cls()
	new.phonopy_tasks, new.bec_tasks = [], []

        # Initialize phonon. Supercell matrix has (3, 3) shape.
        unitcell = atoms_from_structure(gsinp.structure)
	new.scdims = scdims = np.array(scdims)
	if scdims.shape != (3,):
	    raise ValueError("Expecting 3 int in scdims but got %s" % str(scdims))

	supercell_matrix = np.diag(scdims)
	phonopy_kwargs = phonopy_kwargs if phonopy_kwargs is not None else {}
        new.phonon = phonon = Phonopy(unitcell, supercell_matrix, **phonopy_kwargs)
				#primitive_matrix=settings.get_primitive_matrix(),
				#factor=factor,
				#is_auto_displacements=False,
				#dynamical_matrix_decimals=settings.get_dm_decimals(),
				#force_constants_decimals=settings.get_fc_decimals(),
				#symprec=options.symprec,
				#is_symmetry=settings.get_is_symmetry(),
				#use_lapack_solver=settings.get_lapack_solver(),
				#log_level=log_level

        #supercell, primitive = phonon.get_supercell(), phonon.get_primitive()
	displ_kwargs = displ_kwargs if displ_kwargs is not None else {}
        phonon.generate_displacements(**displ_kwargs)
            #distance=0.01,
            #is_plusminus=settings.get_is_plusminus_displacement(),
            #is_diagonal=settings.get_is_diagonal_displacement(),
            #is_trigonal=settings.get_is_trigonal_displacement()

        # Obtain supercells containing respective displacements by get_supercells_with_displacements,
        # which are given by a list of Atoms objects.
        for atoms in phonon.get_supercells_with_displacements():
            #print("atoms", atoms.cell, atoms.scaled_positions)
            sc_struct = structure_from_atoms(atoms)
            #print("sc_struct", sc_struct)

            #sc_gsinp = gsinp.new_for_supercell(sc_struct, divs)
            sc_gsinp = AbinitInput(sc_struct, gsinp.pseudos)
            for k, v in gsinp.items():
                sc_gsinp[k] = v

            # Rescale nband and k-point sampling
	    iscale = int(np.ceil(len(sc_gsinp.structure) / len(gsinp.structure)))
            if "nband" in sc_gsinp:
                sc_gsinp["nband"] = int(gsinp["nband"] * iscale)
                print("gsinp['nband']", gsinp["nband"], "sc_gsinp['nband']", sc_gsinp["nband"])

	    if "ngkpt" in sc_gsinp:
                sc_gsinp["ngkpt"] = (np.rint(np.array(sc_gsinp["ngkpt"]) / scdims)).astype(int)
		print("new sc_gsinp:", sc_gsinp["ngkpt"])
	    #elif "kptrlatt" in sc_gsinp:
            #    sc_gsinp["kptrlatt"] = (np.rint(np.array(sc_gsinp["kptrlatt"]) / iscale)).astype(int)
	    #else:
	    # 	"""Single k-point"""

            sc_gsinp["chksymbreak"] = 0
            sc_gsinp["chkprim"] = 0

	    # TODO
	    # This is the tricky part because variables whose shape depends on natom
            # must be changed to reflect the new supercell.
	    # To solve this problem, we use the database of abinit variables...

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
        directions = phonon.get_displacement_directions()
        file_IO.write_disp_yaml(displacements, supercell, directions=directions,
                                filename=self.outdir.path_in('disp.yaml'))

	# Extact forces from the main output files.
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

	with open(self.outdir.path_in("README"), "wt") as fh:
	    fh.write("To plot bands, use: `phonopy -p band.conf`\n")
	    fh.write("To plot dos, use: `phonopy -p dos.conf`\n")
	    fh.write("To plot bands and dos, use: `phonopy -p band-dos.conf`\n")
	    fh.write("See also:\n")
	    fh.write("\t" + examples_url + "\n")
	    fh.write("\t" + doctags_url + "\n")

        return dict(returncode=0, message="Delta factor computed")
