"""
Converters between abinit/abipy format and other external tools.
Some portions of the code have been imported from the ConvertDDB.py script
developed by Hu Xe, Eric Bousquet and Aldo Romero.
"""

import os
import itertools
import warnings
import numpy as np
from monty.dev import requires

from pymatgen.io.phonopy import get_phonopy_structure, get_pmg_structure
from pymatgen.io.vasp.inputs import Poscar
import abipy.core.abinit_units as abu
from abipy.dfpt.ddb import DdbFile
from abipy.abio.factories import minimal_scf_input
from abipy.electrons.gsr import GsrFile
from monty.os import makedirs_p
try:
    from phonopy import Phonopy
    from phonopy.file_IO import write_FORCE_CONSTANTS, parse_FORCE_CONSTANTS, parse_BORN, parse_FORCE_SETS
    #from phonopy.interface.phonopy_yaml import PhonopyYaml
    from phonopy.interface.calculator import get_default_physical_units, get_force_constant_conversion_factor
except ImportError:
    Phonopy = None


@requires(Phonopy, "phonopy not installed!")
def abinit_to_phonopy(anaddbnc, supercell_matrix, symmetrize_tensors=False, output_dir_path=None,
                      prefix_outfiles="", symprec=1e-5, set_masses=False):
    """
    Converts the interatomic force constants(IFC), born effective charges(BEC) and dielectric
    tensor obtained from anaddb to the phonopy format. Optionally writes the
    standard phonopy files to a selected directory: FORCE_CONSTANTS, BORN (if BECs are available)
    POSCAR of the unit cell, POSCAR of the supercell.

    The conversion is performed taking the IFC in the Wigner–Seitz supercell with weights
    as produced by anaddb and reorganizes them in a standard supercell multiple of the
    unit cell. Operations are vectorized using numpy. This may lead to the allocation of
    large arrays in case of very large supercells.

    Performs a check to verify if the two codes identify the same symmetries and it gives a
    warning in case of failure. Mismatching symmetries may lead to incorrect conversions.

    Args:
        anaddbnc: an instance of AnaddbNcFile. Should contain the output of the IFC analysis,
            the BEC and the dielectric tensor.
        supercell_matrix: the supercell matrix used for phonopy. Any choice is acceptable, however
            the best agreement between the abinit and phonopy results is obtained if this is set to
            a diagonal matrix with on the diagonal the ngqpt used to generate the anaddb.nc.
        symmetrize_tensors: if True the tensors will be symmetrized in the Phonopy object and
            in the output files. This will apply to IFC, BEC and dielectric tensor.
        output_dir_path: a path to a directory where the phonopy files will be created
        prefix_outfiles: a string that will be added as a prefix to the name of the written files
        symprec: distance tolerance in Cartesian coordinates to find crystal symmetry in phonopy.
            It might be that the value should be tuned so that it leads to the the same symmetries
            as in the abinit calculation.
        set_masses: if True the atomic masses used by abinit will be added to the PhonopyAtoms
            and will be present in the returned Phonopy object. This should improve compatibility
            among abinit and phonopy results if frequencies needs to be calculated.

    Returns:
        An instance of a Phonopy object that contains the IFC, BEC and dieletric tensor data.
    """

    ifc = anaddbnc.ifc
    nac_params = None
    becs = None
    epsinf = None
    if anaddbnc.becs is not None and anaddbnc.epsinf is not None:
        becs = anaddbnc.becs.values
        epsinf = anaddbnc.epsinf

        # according to the phonopy website 14.399652 is not the coefficient for abinit
        # probably it relies on the other conventions in the output.
        nac_params = {"born": becs, "dielectric": epsinf, "factor": 14.399652}

    s = anaddbnc.structure

    phon_at = get_phonopy_structure(s)
    if set_masses:
        phon_at.masses = [anaddbnc.amu[n] for n in phon_at.numbers]

    # use phonopy to get the proper supercell given by the primitive and the matrix
    # and convert it to pymatgen
    phonon = Phonopy(phon_at, supercell_matrix, primitive_matrix=np.eye(3), nac_params=nac_params,
                     symprec=symprec)
    phon_supercell = phonon.get_supercell()
    supercell = get_pmg_structure(phon_supercell)

    abi_hall_num = s.abi_spacegroup.get_spglib_hall_number()
    spglib_hall_num = phonon.symmetry.dataset["hall_number"]
    if abi_hall_num != spglib_hall_num:
        warnings.warn("The hall number obtained based on the DDB symmetries differs "
                      f"from the one calculated with spglib: {abi_hall_num} versus "
                      f"{spglib_hall_num}. The conversion may be incorrect. Try changing symprec.")

    # convert to phonopy units
    at_cart = ifc.atoms_cart_coord * abu.Bohr_Ang
    ifccc = ifc.ifc_cart_coord * abu.Ha_eV / abu.Bohr_Ang ** 2
    weights = ifc.ifc_weights
    latt = supercell.lattice

    ifcph = np.zeros((len(s), len(supercell), 3, 3))

    # loop over the atoms in the primitive cell
    # other operations are vectorized using numpy arrays. Some array may require large allocations
    for i, (site, c_list, w_list) in enumerate(zip(s, at_cart, weights)):

        ind_w = np.where(w_list > 0)
        ifccc_loc = ifccc[i, ind_w[0]]

        w_list = w_list[ind_w]
        c_list = c_list[ind_w]

        # align the coordinates of the first atom in the list (the site under consideration)
        # with the site in the primitive cell.
        c_list = c_list - c_list[0] + site.coords

        # convert to fractional coordinates as needed by the Lattice to get the distances
        f_list = latt.get_fractional_coords(c_list)
        sc_fcoords = supercell.frac_coords

        # construct the list of sites of the supercell that are closer to sites in
        # the primitive cell
        dist_and_img = [latt.get_distance_and_image(f_list[0], fc) for fc in sc_fcoords]
        # the function gives the translation of the image, but it should be applied to the coordinates.
        # Only the positions are needed
        nearest_sc_fcoords = [fc + trasl for (_, trasl), fc in zip(dist_and_img, sc_fcoords)]

        # divide by the corresponding weights. Elements with weights 0 were discarded above
        ifccc_loc = np.transpose(ifccc_loc, (0, 2, 1)) / w_list[:, None, None]

        # create an array with all the possible pairs
        # instantiating this array seems slow but seems still faster than the required loops
        coord_pairs = np.array(list(itertools.product(nearest_sc_fcoords, f_list)))

        # find the pairs that match between the coordinates of the modified supercell and the f_list
        ind_match = np.where(np.abs(coord_pairs[:, 0] - coord_pairs[:, 1]).sum(axis=1) < 1e-6)[0]
        # set the ifc for phonopy in the final array corresponding to the matching indices.
        n_points_f_list = len(f_list)
        ifcph[i, ind_match // n_points_f_list] = ifccc_loc[ind_match % n_points_f_list]

    phonon.set_force_constants(ifcph)
    if symmetrize_tensors:
        phonon.symmetrize_force_constants()

    if output_dir_path:
        makedirs_p(output_dir_path)

        fc_filepath = os.path.join(output_dir_path, prefix_outfiles+"FORCE_CONSTANTS")
        write_FORCE_CONSTANTS(phonon.get_force_constants(), fc_filepath)

        if becs is not None and epsinf is not None:
            born_filepath = os.path.join(output_dir_path, prefix_outfiles+"BORN")
            write_BORN(phon_at, borns=becs, epsilon=epsinf, filename=born_filepath,
                       symmetrize_tensors=symmetrize_tensors)

        poscar_filepath = os.path.join(output_dir_path, prefix_outfiles+"POSCAR")
        poscar = Poscar(s)
        poscar.write_file(poscar_filepath, significant_figures=15)

        supercell_filepath = os.path.join(output_dir_path, prefix_outfiles+"supercell_POSCAR")
        superce_poscar = Poscar(supercell)
        superce_poscar.write_file(supercell_filepath, significant_figures=15)

    return phonon


@requires(Phonopy, "phonopy not installed!")
def phonopy_to_abinit(unit_cell, supercell_matrix, out_ddb_path, ngqpt=None, qpt_list=None,
                      force_constants=None, force_sets=None, born=None,
                      primitive_matrix="auto", symprec=1e-5, tolsym=None, supercell=None,
                      calculator=None, manager=None, workdir=None, pseudos=None, verbose=False):
    """
    Converts the data from phonopy to an abinit DDB file. The data can be provided
    in form of arrays or paths to the phonopy files that should be parsed.
    The minimal input should contains the FORCE_CONSTANTS or FORCE_SETS.
    If BORN is present the Born effective charges (BEC) and dielectric
    tensor will also be added to the DDB.

    The best agreement is obtained with supercell_matrix and ngqpt being
    equivalent (i.e. supercell_matrix a diagonal matrix with ngqpt as diagonal
    elements). Non diagonal supercell_matrix are allowed as well, but the information
    encoded in the DDB will be the result of an interpolation done through phonopy.

    Phonopy is used to convert the IFC to the dynamical matrix. However, in order to
    determine the list of q-points in the irreducible Brillouin zone and to prepare the
    base for the final DDB file, abinit will be called for a very short and inexpensive run.

    Performs a check to verify if the two codes identify the same symmetries and it gives a
    warning in case of failure. Mismatching symmetries may lead to incorrect conversions.

    Args:
        unit_cell: a |Structure| object that identifies the unit cell used for the phonopy
            calculation.
        supercell_matrix: a 3x3 array representing the supercell matrix used to generated the
            forces with phonopy.
        out_ddb_path: a full path to the file where the new DDB will be written
        ngqpt: a list of 3 elements indicating the grid of q points that will be used in the DDB.
        qpt_list: alternatively to ngqpt an explicit list of q-points can be provided here.
            At least one among ngqpt and qpt_list should be defined.
        force_constants: an array with shape (num atoms unit cell, num atoms supercell, 3, 3)
            containing the force constants. Alternatively a string with the path to the
            FORCE_CONSTANTS file. This or force_set should be defined. If both given this
            has precedence.
        force_sets: a dictionary obtained from the force sets generated with phonopy.
            Alternatively a string with the path to the FORCE_SETS file. This or force_constants
            should be defined.
        born: a dictionary with "dielectric" and "born" keywords as obtained from the nac_params
            in phonopy. Alternatively a string with the path to the BORN file. Notice that
            the "factor" attribute is not taken into account, so the values should be in
            default phonopy units.
        primitive_matrix: a 3x3 array with the primitive matrix passed to Phonopy. "auto" will
            use spglib to try to determine it automatically. If the DDB file should contain the
            actual unit cell this should be the identity matrix.
        symprec: distance tolerance in Cartesian coordinates to find crystal symmetry in phonopy.
            It might be that the value should be tuned so that it leads to the the same symmetries
            as in the abinit calculation.
        tolsym: Gives the tolerance to identify symmetries in abinit. See abinit documentation for
            more details.
        supercell: if given it should represent the supercell used to get the force constants,
            without any perturbation. It will be used to match it to the phonopy supercell
            and sort the IFC in the correct order.
        calculator: a string with the name of the calculator. Will be used to set the conversion
            factor for the force constants coming from phonopy.
        manager: |TaskManager| object. If None, the object is initialized from the configuration file
        pseudos: List of filenames or list of |Pseudo| objects or |PseudoTable| object. It will be
            used by abinit to generate the base DDB file. If None the abipy.data.hgh_pseudos.HGH_TABLE
            table will be used.
        verbose: verbosity level. Set it to a value > 0 to get more information
        workdir: path to the directory where the abinit calculation will be executed.

    Returns:
        a DdbFile instance of the file written in out_ddb_path.
    """

    if ngqpt is None and qpt_list is None:
        raise ValueError("at least one among nqgpt and qpt_list should be defined")

    if force_sets is None and force_constants is None:
        raise ValueError("at least one of force_sets and force_constants should be provided")

    phon_at = get_phonopy_structure(unit_cell)

    if isinstance(force_constants, str):
        force_constants = parse_FORCE_CONSTANTS(filename=force_constants)
    elif force_constants is not None:
        force_constants = np.array(force_constants)
        force_sets = None

    if isinstance(force_sets, str):
        force_sets = parse_FORCE_SETS(filename=force_sets)

    # no nac_params here, otherwise they will be used for the interpolation
    phonon = Phonopy(phon_at, supercell_matrix, primitive_matrix=primitive_matrix, nac_params=None,
                     symprec=symprec, calculator=calculator)

    primitive = get_pmg_structure(phonon.primitive)

    if isinstance(born, str):
        born = parse_BORN(phonon.primitive, filename=born)

    if supercell is not None:
        ph_supercell = get_pmg_structure(phonon.supercell)
        if not np.allclose(supercell.lattice.matrix, ph_supercell.lattice.matrix):
            raise RuntimeError("The lattice of the supercells do not match")
        sc_mapping = []
        for i, site_orig in enumerate(supercell):
            for j, site_ph in enumerate(ph_supercell):
                d = supercell.lattice.get_distance_and_image(site_orig.frac_coords, site_ph.frac_coords)[0]
                if d < 1e-5:
                    sc_mapping.append(j)
                    break
            else:
                raise RuntimeError(f"Could not find a match for site {i} with coords "
                                   f"{site_orig.cart_coords} in the supercell.")

        # cross check that the same atom was not matched twice
        n_matches = len(set(sc_mapping))
        if n_matches < len(supercell):
            raise RuntimeError(f"Found matches for {n_matches} different atoms in the supercell: {sc_mapping}")

        force_constants = force_constants[:, sc_mapping]

    if force_constants is not None:
        phonon.set_force_constants(force_constants)
    else:
        phonon.dataset = force_sets
        phonon.produce_force_constants()

    if calculator:
        units = get_default_physical_units(calculator)
        fc_factor = get_force_constant_conversion_factor(units["force_constants_unit"], None)
        phonon.set_force_constants(phonon.force_constants * fc_factor)

    if pseudos is None:
        from abipy.data.hgh_pseudos import HGH_TABLE
        pseudos = HGH_TABLE

    inp = minimal_scf_input(primitive, pseudos)

    # get the qpoints list if not defined
    if qpt_list is None:
        inp["ngkpt"] = ngqpt
        qpt_list = inp.abiget_ibz(verbose=verbose)[0]

    dm_list = get_dm(phonon, qpt_list, primitive)

    if born is not None:
        # for the conversion of the BEC the zion (i.e. the ionic charge of the pseudo)
        # it is an additive factor and should be the same that goes in the header of the DDB,
        # so take it from the pseudos used to generate it.
        zion = inp.valence_electrons_per_atom
        born_data = generate_born_deriv(born, zion, primitive)
    else:
        born_data = None

    inp = minimal_scf_input(primitive, pseudos)
    if tolsym is not None:
        inp["tolsym"] = tolsym
    task = inp.run_in_shell(workdir=workdir, manager=manager, verbose=verbose)

    # use the output of abinit to check that the spacegroup identified by
    # phonopy and abinit are the same.
    with GsrFile(task.opath_from_ext("GSR.nc")) as gsr:
        abi_spg = gsr.structure.abi_spacegroup.spgid
    spglib_spg = phonon.symmetry.dataset["number"]
    if abi_spg != spglib_spg:
        warnings.warn("The space group number obtained based on the DDB symmetries differs "
                      f"from the one calculated with spglib: {abi_spg} versus "
                      f"{spglib_spg}. The convertion may be incorrect. Try changing symprec or tolsym.")

    tmp_ddb_path = task.opath_from_ext("DDB")

    ddb = DdbFile(tmp_ddb_path)
    # remove the blocks generated by the calculation and that are meaningless
    ddb.remove_block(dord=0)
    ddb.remove_block(dord=1)

    add_data_ddb(ddb, dm_list, qpt_list, born_data)

    ddb.write(out_ddb_path)

    new_ddb = DdbFile(out_ddb_path)
    return new_ddb


def generate_born_deriv(born, zion, structure):
    """
    Helper function to generate the portion of the derivatives in the DDB
    that are related to the Born effective charges and dielectric tensor,
    starting from the data available in phonopy format.

    Args:
        born: a dictionary with "dielectric" and "born" keywords as obtained from the nac_params
            in phonopy.
        zion: the ionic charge of each atom in the system. It should be in the same order
            as the one present in the header of the DDB.
        structure: a pymatgen |Structure| of the unit cell.

    Returns:
        a complex numpy array with shape (len(structure)+2, 3, len(structure)+2, 3). Only the
        parts relative to the BECs and dielectric tensors will be filled.
    """
    natoms = len(structure)
    mpert = natoms + 2 # only these perturbations are needed here
    born_data = np.zeros((3, mpert, 3, mpert), dtype=complex)

    eps_e = born["dielectric"]
    bec = np.array(born["born"]).transpose((0, 2, 1))

    # Assume that for the generated DDB acell = [1,1,1] and rprimd = rprim. Should be in Bohr.
    rprim = structure.lattice.matrix * abu.Ang_Bohr
    rprimd = rprim

    volume_bohr = np.dot(rprimd[:, 0], np.cross(rprimd[:, 1], rprimd[:, 2]))
    gprimd = np.linalg.inv(rprimd)
    dij = np.identity(3)
    # BEC
    for ipert1 in range(natoms):  # ipert1 is atom position deriv
        ipert2 = natoms + 1  # E field deriv
        dm1 = np.matmul(rprimd,
                        np.matmul(bec[ipert1, :, :] - dij[:, :] * zion[ipert1], gprimd))
        for idir1 in range(3):
            for idir2 in range(3):
                born_data[idir1, ipert1, idir2, ipert2] = dm1[idir1, idir2] * 2 * np.pi + 0.0j
                born_data[idir2, ipert2, idir1, ipert1] = dm1[idir1, idir2] * 2 * np.pi + 0.0j
    # epsinf
    ipert1 = natoms + 1
    ipert2 = natoms + 1
    dm1 = np.matmul(gprimd.transpose(),
                    np.matmul(dij[:, :] - eps_e[:, :], gprimd))
    for idir1 in range(3):
        for idir2 in range(3):
            born_data[idir1, ipert1, idir2, ipert2] = dm1[idir1, idir2] * np.pi * volume_bohr + 0.0j
    born_data = born_data.transpose((1, 0, 3, 2))
    return born_data


@requires(Phonopy, "phonopy not installed!")
def get_dm(phonon, qpt_list, structure):
    """
    Helper function to generate the dynamical matrix in the abinit conventions
    for a list of q-points from a Phonopy object.

    Args:
        phonon: a Phonopy object with force constants.
        qpt_list: a list of fractional coordinates of q-points for which the
            dynamical matrix should be generated.
        structure: a pymatgen |Structure| of the primitive cell.

    Returns:
        a list of arrays with the dynamical matrices of the selected q-points.
    """
    natom = len(structure)
    rprim = structure.lattice.matrix * abu.Ang_Bohr
    # assume acell is [1., 1., 1.]
    rprimd = rprim
    masses = phonon.masses
    dm_list = []
    mass_tile = np.tile(masses, (3, 1)).T
    mass_matrix = np.sqrt(np.einsum("ij,kl", mass_tile, mass_tile))
    # get the difference in coordinates for each pair of atoms
    # diff_coords[i,j] is equivalent to structure[i].coords - structure[j].coords
    coord_matrix = np.repeat(structure.cart_coords[np.newaxis, :, :], natom, axis=0)
    diff_coords = np.transpose(coord_matrix, (1, 0, 2)) - coord_matrix
    rlatt = structure.lattice.reciprocal_lattice
    for q in qpt_list:
        q_cart = rlatt.get_cartesian_coords(q)
        # the phase exp(-i * q.r)
        phase = np.exp(-1j * np.einsum("ijk,k", diff_coords, q_cart))
        dm = phonon.get_dynamical_matrix_at_q(q).reshape(natom, 3, natom, 3)
        # the following is rprim * dm[ipert1,:,ipert2,:] * rprim.T
        dm = np.einsum("ij,kjlm,nm->kiln", rprimd, dm, rprimd)
        dm = dm * mass_matrix / phase[:, None, :, None]
        dm *= abu.eV_Ha / abu.Ang_Bohr**2

        dm_list.append(dm)

    return dm_list


def add_data_ddb(ddb, dm_list, qpt_list, born_data):
    """
    Helper function to add the blocks for the dynamical matrix and BECs to a DdbFile object.

    Args:
        ddb: a DdbFile object to be modified.
        dm_list: the list of dynamical matrices to be added.
        qpt_list: the list of q-points corresponding to dm_list.
        born_data: the data corresponding to BECs and dielectric tensor. If None
            these part will not be set in the DDB.
    """
    dm_data = {}
    natom = len(ddb.structure)
    for q, dm in zip(qpt_list, dm_list):
        q_data = {}
        for ipert1 in range(natom):
            for idir1 in range(3):
                for ipert2 in range(natom):
                    for idir2 in range(3):
                        q_data[(idir1+1, ipert1+1, idir2+1, ipert2+1)] = dm[ipert1, idir1, ipert2, idir2]

        # for gamma set also the born data if present
        if np.allclose(q, (0, 0, 0)) and born_data is not None:
            ipert2 = natom + 1
            for ipert1 in range(natom):
                for idir1 in range(3):
                    for idir2 in range(3):
                        q_data[(idir1+1, ipert1+1, idir2+1, ipert2+1)] = born_data[ipert1, idir1, ipert2, idir2]
                        q_data[(idir2+1, ipert2+1, idir1+1, ipert1+1)] = born_data[ipert2, idir2, ipert1, idir1]

            for idir1 in range(3):
                for idir2 in range(3):
                    q_data[(idir1+1, ipert2+1, idir2+1, ipert2+1)] = born_data[ipert2, idir1, ipert2, idir2]

        dm_data[tuple(q)] = q_data

    ddb.set_2nd_ord_data(dm_data, replace=True)


@requires(Phonopy, "phonopy not installed!")
def tdep_to_abinit(unit_cell, fc_path, supercell_matrix, supercell, out_ddb_path, ngqpt=None,
                   qpt_list=None, primitive_matrix="auto", lotosplitting_path=None, symprec=1e-5,
                   tolsym=None, manager=None, workdir=None, pseudos=None, verbose=False):
    """
    Converts the files produced by TDEP to an abinit DDB file. If the lotosplitting
    file is provided the BEC and dielectric tensor will also be added to the DDB.

    The conversion is performed by first extracting the force constants in phonopy format
    and then using phonopy_to_abinit_py to generate the DDB file. See the phonopy_to_abinit_py
    docstring for further details about the second step of the conversion.
    Notice that the supercell used by TDEP should be provided in order to properly
    sort the IFC to match the order required by phonopy.

    Args:
        unit_cell: a |Structure| object that identifies the unit cell used for the TDEP
            calculation.
        fc_path: the path to the forceconstants file produced by TDEP.
        supercell_matrix: a 3x3 array representing the supercell matrix used to generated the
            forces with TDEP.
        supercell: the supercell used by TDEP to get the force constants, without any
            perturbation (usually named with extension ssposcar). It will be used to match it to
            the phonopy supercell and sort the IFC in the correct order.
        out_ddb_path: a full path to the file where the new DDB will be written
        ngqpt: a list of 3 elements indicating the grid of q points that will be used in the DDB.
        qpt_list: alternatively to ngqpt an explicit list of q-points can be provided here.
            At least one among ngqpt and qpt_list should be defined.
        primitive_matrix: a 3x3 array with the primitive matrix passed to Phonopy. "auto" will
            use spglib to try to determine it automatically. If the DDB file should contain the
            actual unit cell this should be the identity matrix.
        lotosplitting_path: path to the lotosplitting file produced by TDEP. If None no BEC
            contribution will be set to the DDB.
        symprec: distance tolerance in Cartesian coordinates to find crystal symmetry in phonopy.
            It might be that the value should be tuned so that it leads to the the same symmetries
            as in the abinit calculation.
        tolsym: Gives the tolerance to identify symmetries in abinit. See abinit documentation for
            more details.
        manager: |TaskManager| object. If None, the object is initialized from the configuration file
        pseudos: List of filenames or list of |Pseudo| objects or |PseudoTable| object. It will be
            used by abinit to generate the base DDB file. If None the abipy.data.hgh_pseudos.HGH_TABLE
            table will be used.
        verbose: verbosity level. Set it to a value > 0 to get more information
        workdir: path to the directory where the abinit calculation will be executed.

    Returns:
        a DdbFile instance of the file written in out_ddb_path.
    """

    fc = parse_tdep_fc(fc_path, unit_cell, supercell)

    born = None
    if lotosplitting_path:
        eps, becs = parse_tdep_lotosplitting(lotosplitting_path)
        born = {"dielectric": eps, "born": becs, "factor": 1}

    ddb = phonopy_to_abinit(unit_cell=unit_cell, force_constants=fc, supercell_matrix=supercell_matrix, ngqpt=ngqpt,
                            qpt_list=qpt_list, out_ddb_path=out_ddb_path, born=born, pseudos=pseudos,
                            primitive_matrix=primitive_matrix, supercell=supercell, manager=manager,
                            workdir=workdir, symprec=symprec, verbose=verbose, tolsym=tolsym)

    return ddb


def parse_tdep_fc(fc_path, unit_cell, supercell):
    """
    Parses a forceconstants file produced by TDEP an converts it to an array in the
    phonopy format.

    Args:
        fc_path: path to the forceconstants file
        unit_cell: a |Structure| object with the unit cell used for the calculation
            in TDEP.
        supercell: the supercell used for the calculation in TDEP.

    Returns:
        a comple numpy array with shape (len(unit_cell), len(supercell), 3, 3)
    """
    natoms = len(unit_cell)
    fc = np.zeros((natoms, len(supercell), 3, 3))
    # parse the fc file and put it in the format required for phonopy
    # this will be sorted according to the tdep supercell atoms order, not the phonopy supercell.
    with open(fc_path, "rt") as f:
        lines = f.readlines()
    u_latt = unit_cell.lattice
    sc_latt = supercell.lattice
    iline = 2
    for iat in range(natoms):
        n_neighbours = int(lines[iline].split()[0])
        iline += 1
        for jn in range(n_neighbours):
            iprim = int(lines[iline].split()[0]) - 1
            r = [float(sp) for sp in lines[iline + 1].split()]
            # find matching atom in the supercell
            fcoords = sc_latt.get_fractional_coords(u_latt.get_cartesian_coords(r) + unit_cell[iprim].coords)
            for isc, site_sc in enumerate(supercell):
                d = supercell.lattice.get_distance_and_image(fcoords, site_sc.frac_coords)[0]
                if d < 1e-5:
                    break
            else:
                raise RuntimeError(f"could not find a match for: {lines[iline]}")

            for k in range(3):
                fc[iat, isc, k] = [float(sp) for sp in lines[iline + 2 + k].split()]
            iline += 5
    return fc


def parse_tdep_lotosplitting(filepath):
    """
    Parses the lotosplitting file produced by TDEP and transforms them in
    the phonopy format for Born effective charges and dielectric tensor.

    Args:
        filepath: path to the lotosplitting file.

    Returns:
        a tuple with dielectric tensor and Born effective charges.
    """
    with open(filepath) as f:
        lines = f.readlines()

    values = []
    for l in lines:
        if not l.strip():
            continue
        values.append([float(v) for v in l.split()[:3]])

    if len(values) % 3 != 0:
        raise RuntimeError(f"The data parsed has an unexpected shape: {np.shape(values)}")

    eps = np.array(values[:3])
    natoms = len(values) // 3 - 1
    born = np.reshape(values[3:], (natoms, 3, 3))

    return eps, born


def write_tdep_lotosplitting(eps, born, filepath="infile.lotosplitting", fmt="%14.10f"):
    """
    Writes an lotosplitting file starting from arrays containing dielectric tensor and
    Born effective charges.

    Args:
        eps: a 3x3 array with the dielectric tensor.
        born: an array with the Born effective charges.
        filepath: the path where the lotosplitting file should be written.
        fmt: the format for the float numbers.
    """
    values = np.concatenate((eps[None, :, :], born), axis=0)
    values = values.reshape((3 * len(values), 3))
    np.savetxt(filepath, values, fmt=fmt)


def born_to_lotosplitting(born, lotosplitting_path="infile.lotosplitting"):
    """
    Converted of a file from the BORN file produced from phonopy to the lotosplitting
    file used by TDEP.

    Args:
        born: a dictionary with "dielectric" and "born" keywords as obtained from the nac_params
            in phonopy. Notice that the "factor" attribute is not taken into account, so the
            values should be in default phonopy units.
        lotosplitting_path: the path where the lotosplitting file should be written.
    """

    eps = born["dielectric"]
    becs = born["born"]
    write_tdep_lotosplitting(eps, becs, lotosplitting_path)


@requires(Phonopy, "phonopy not installed!")
def write_BORN(primitive, borns, epsilon, filename="BORN", symmetrize_tensors=False):
    """
    Helper function imported from phonopy.file_IO.
    Contrarily to the original, it does not symmetrize the tensor.
    """
    lines = get_BORN_lines(primitive, borns, epsilon, symmetrize_tensors=symmetrize_tensors)
    with open(filename, 'w') as w:
        w.write('\n'.join(lines))


@requires(Phonopy, "phonopy not installed!")
def get_BORN_lines(unitcell, borns, epsilon,
                   factor=None,
                   primitive_matrix=None,
                   supercell_matrix=None,
                   symprec=1e-5, symmetrize_tensors=False):
    """
    Helper function imported from phonopy.file_IO that exposes the
    option of not symmetrizing the tensor.
    """
    from phonopy.structure.symmetry import elaborate_borns_and_epsilon
    borns, epsilon, atom_indices = elaborate_borns_and_epsilon(
        unitcell, borns, epsilon, symmetrize_tensors=symmetrize_tensors,
        primitive_matrix=primitive_matrix,
        supercell_matrix=supercell_matrix,
        symprec=symprec)

    text = "# epsilon and Z* of atoms "
    text += ' '.join(["%d" % n for n in atom_indices + 1])
    lines = [text, ]
    lines.append(("%13.8f " * 9) % tuple(epsilon.flatten()))
    for z in borns:
        lines.append(("%13.8f " * 9) % tuple(z.flatten()))
    return lines
