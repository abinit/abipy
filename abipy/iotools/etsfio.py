"""ETSF-IO specifications."""
from __future__ import print_function, division

# Version of the ETSF-IO specifications.
__version__ = "3.3"

#: List of valid ETSF-IO keywords.
keywords = []

_dict_dims = {
#"character_string_length",     : "str_len",
"real_or_complex_coefficients"  : "cplex_ug",
"real_or_complex_density"       : "cplex_den",
"real_or_complex_gw_corrections": "cplex_qp",
"real_or_complex_potential"     : "cplex_pot",
"real_or_complex_wavefunctions" : "cplex_ur",
"number_of_symmetry_operations" : "nsym",
"number_of_atoms"               : "natom",
"number_of_atom_species"        : "ntypat",
"max_number_of_states"          : "mband",
"number_of_kpoints"             : "nkpt",
"number_of_spins"               : "nsppol",
"number_of_spinor_components"   : "nspinor",
"number_of_components"          : "nspden",
"max_number_of_coefficients"    : "mpw",
"number_of_grid_points_vector1" : "nfft1",
"number_of_grid_points_vector2" : "nfft2",
"number_of_grid_points_vector3" : "nfft3",
}

#: Dimensions specified by the standard.
keywords.extend(_dict_dims.keys())

_dict_basisdata = {
"basis_set"                          : "basis_set",
"kinetic_energy_cutoff"              : "ecut",
"number_of_coefficients"             : "npwarr",
"reduced_coordinates_of_plane_waves" : "kg",
}

#: Information on the basis set.
keywords.extend(_dict_basisdata.keys())

_dict_main = {
"density"                        : "rhor",
"exchange_potential"             : "vx",
"correlation_potential"          : "vc",
"exchange_correlation_potential" : "vxc",
"coefficients_of_wavefunctions"  : "set_of_ug",
"real_space_wavefunctions"       : "set_of_ur",
}

#: ETSF-IO variables.
keywords.extend(_dict_main.keys())

_dict_geometry = {
"primitive_vectors"             : "rprimd",
"reduced_symmetry_matrices"     : "symrel",
"reduced_symmetry_translations" : "tnons",
"space_group"                   : "space_group_id",
"atom_species"                  : "typat",
"reduced_atom_positions"        : "xred",
"atomic_numbers"                : "znucl_type",
"atom_species_names"            : "atom_species_names",
"chemical_symbols"              : "chemical_symbols",
"pseudopotential_types"         : "pseudopotential_types",
#symmorphic
}

#: Specification of the crystalline structure.
keywords.extend(_dict_geometry.keys())

_dict_kpoints = {
"reduced_coordinates_of_kpoints" : "kpoints",
"kpoint_weights"                 : "wtk",
"kpoint_grid_shift"              : "shiftk",
"kpoint_grid_vectors"            : None,  # FIXME Spec are not clear here!
"monkhorst_pack_folding"         : "ngkpt",
}

#: Specification of BZ sampling.
keywords.extend(_dict_kpoints.keys())


_dict_electrons = {
"number_of_electrons"    : "nelect",
"exchange_functional"    : "exchange_functional",
"correlation_functional" : "correlation_functional",
"fermi_energy"           : "fermie",
"smearing_scheme"        : "smearing_scheme",
"smearing_width"         : "tsmear",
"number_of_states"       : "nband_sk",
"eigenvalues"            : "energies",
"occupations"            : "occfacts",
#k dependent
}

#: Specification of the band structure.
keywords.extend(_dict_electrons.keys())


_dict_optional = {
"valence_charges"  : "valence_charges",
}

#: Optional values.
keywords.extend(_dict_optional.keys())

#: Unofficial extensions.
_dict_geometry_unofficial = {
"cartesian_forces"        : "cart_forces",
"cartesian_stress_tensor" : "cart_stress",
#"reduced_forces"          : "red_forces",
#"reduced_stress_tensor" : "red_stress",
}
keywords.extend(_dict_geometry_unofficial.keys())

# Build a dictionary mapping abipy variables to the ETSF-IO keywords
_dicts = [
  _dict_dims,
  _dict_basisdata,
  _dict_main,
  _dict_geometry,
  _dict_kpoints,
  _dict_electrons,
  _dict_optional,
  #
  _dict_geometry_unofficial,
]

abipy2etsfio = {}

for d in _dicts:
    for k, v in d.iteritems():
        if v in abipy2etsfio:
            raise Exception("value: %s compares more than once\n" % v)
        abipy2etsfio[v] = k


def etsfio2abipy(key):
    """
    Translate an ETSF-IO keyword to an abipy name.

    Raise KeyError if key is not a valid ETSF-IO name.
    """
    for d in _dicts:
        for k, v in d.iteritems():
            if k == key:
                return v
    raise KeyError(str(key))


if __name__ == "__main__":
    import doctest
    doctest.testmod()
