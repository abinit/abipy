#!/usr/bin/env python
r"""
Band structure w/wo magnetization
=================================

Calculation of the band structure of Fe with and without magnetization,
including L-projected (FATBANDS and FATDOS)
See also <~abinit/tutorial/Input/tspin_1.in>
"""
import os
import sys
import abipy.abilab as abilab
import abipy.flowtk as flowtk


CIF_STRING = """
data_BaO
_audit_creation_method           'pos2cif.pl'
_cell_length_a               6.159607710989539
_cell_length_b               6.159607710989539
_cell_length_c               6.159607710989539
_cell_angle_alpha            90
_cell_angle_beta             90
_cell_angle_gamma            90
_symmetry_space_group_H-M        'P1'
_symmetry_Int_Tables_number      '1'
_symmetry_cell_setting           'triclinic'
loop_
_symmetry_equiv_pos_as_xyz
x,y,z

loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_occupancy
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_U_iso_or_equiv
Ba1  Ba  1.0000  0.25000  0.25000  0.25000  0.0000
Ba2  Ba  1.0000  0.75000  0.75000  0.25000  0.0000
Ba3  Ba  1.0000  0.75000  0.25000  0.75000  0.0000
Ba4  Ba  1.0000  0.25000  0.75000  0.75000  0.0000
O1  O  1.0000  0.50000  0.00000  0.00000  0.0000
O2  O  1.0000  0.00000  0.50000  0.00000  0.0000
O3  O  1.0000  0.00000  0.00000  0.50000  0.0000
O4  O  1.0000  0.50000  0.50000  0.00000  0.0000
O5  O  1.0000  0.50000  0.00000  0.50000  0.0000
O6  O  1.0000  0.00000  0.50000  0.50000  0.0000
"""


def build_flow(options):
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(sys.argv[0]).replace(".py", "").replace("run_", "flow_")

    # Create the Flow.
    flow = flowtk.Flow(options.workdir, manager=options.manager)

    structure = abilab.Structure.from_string(CIF_STRING, fmt="cif", primitive=True)

    pseudos_list = [
            ("foo", "bar"),
            ("foo1", "bar1"),
    ]

    for pseudos in pseudos_list:
        scf_input = abilab.AbinitInput(structure=structure, pseudos=pseudos)

        nval = structure.num_valence_electrons(pseudos)
        nband = (nval / 2) + 10

        # Global variables
        scf_input.set_vars(
            ecut=18,
            nband=nband,
            occopt=3,
            tsmear=0.01,
            ngkpt=[4, 4, 4],
            shiftk=[0.5, 0.5, 0.5],
            tolvrs=1e-8,
            #paral_kgb=paral_kgb,
        )

        # Build a BandStructureWork from the scf_input with the given nsppol and add it to the flow
        # L-projection (prtdos 3) is used by default.
        work = flowtk.BandStructureWork.from_scf_input(scf_input, dos_ngkpt=(8, 8, 8), nb_extra=4)
        flow.register_work(work)

    return flow


# This block generates the thumbnails in the AbiPy gallery.
# You can safely REMOVE this part if you are using this script for production runs.
if os.getenv("READTHEDOCS", False):
    __name__ = None
    import tempfile
    options = flowtk.build_flow_main_parser().parse_args(["-w", tempfile.mkdtemp()])
    build_flow(options).graphviz_imshow()


@flowtk.flow_main
def main(options):
    """
    This is our main function that will be invoked by the script.
    flow_main is a decorator implementing the command line interface.
    Command line args are stored in `options`.
    """
    return build_flow(options)


if __name__ == "__main__":
    sys.exit(main())
