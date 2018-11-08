#!/usr/bin/env python
"""
Script to analyze/export/visualize the crystal structure saved in the netcdf files produced by ABINIT.
"""
from __future__ import unicode_literals, division, print_function, absolute_import

import sys
import os
import argparse
import numpy as np

from pprint import pprint
from tabulate import tabulate
from warnings import warn
from monty.string import marquee
from monty.functools import prof_main
from monty.termcolor import cprint
from pymatgen.io.vasp.outputs import Xdatcar
from abipy import abilab
from abipy.core.symmetries import AbinitSpaceGroup
from abipy.core.kpoints import Ktables, Kpoint, IrredZone
from abipy.core.structure import diff_structures
from abipy.iotools.visualizer import Visualizer
from abipy.iotools.xsf import xsf_write_structure
from abipy.abio import factories


#def remove_equivalent_atoms(structure):
#    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
#    spgan = SpacegroupAnalyzer(structure) #, symprec=symprec, angle_tolerance=angle_tolerance)
#    spgdata = spgan.get_symmetry_dataset()
#    equivalent_atoms = spgdata["equivalent_atoms"]
#    mask = np.zeros(len(structure), dtype=np.int)
#    for pos, eqpos in enumerate(equivalent_atoms):
#        mask[eqpos] += 1
#
#    indices = [i for i, m in enumerate(mask) if m == 0]
#    new = structure.copy()
#    new.remove_sites(indices)
#    with open("foo.abi", "wt") as fh:
#        fh.write(new.abi_string)


def save_structure(structure, options):
    """Save structure to file."""
    if not options.savefile: return
    print("Saving structure to %s" % options.savefile)
    if os.path.exists(options.savefile):
        backup = options.savefile + ".bkp"
        print("%s already exists. Saving backup copy to: %s" % (options.savefile, backup))
        os.rename(options.savefile, backup)

    structure.to(filename=options.savefile)


def check_ordered_structure(structure):
    """Print a warning and sys.exit 1 if structure is disordered."""
    if not structure.is_ordered:
        cprint("""
Cannot handle disordered structure with fractional site occupancies.
Use OrderDisorderedStructureTransformation or EnumerateStructureTransformation
to build an appropriate supercell from partial occupancies.""", color="magenta")
        sys.exit(1)


def get_epilog():
    return """\
Usage example:

###################
# Space group tools
###################

  abistruct.py spglib FILE                 => Read structure from FILE and analyze it with spglib.
  abistruct.py abispg FILE                 => Read structure from FILE, and compute ABINIT space group.
  abistruct.py abisanitize FILE            => Read structure from FILE, call abisanitize, compare structures
                                              and save "abisanitized" structure to file.
  abistruct.py conventional FILE           => Read structure from FILE, generate conventional structure
                                              following doi:10.1016/j.commatsci.2010.05.010
  abistruct.py proto FILE                 => Read structure from FILE, find possible crystallographic prototypes:
                                             in the AFLOW library of crystallographic prototypes.
                                             http://doi.org/10.1016/j.commatsci.2017.01.017

##################
# Conversion tools
##################

  abistruct.py convert FILE                => Print the ABINIT variables defining the structure.
  abistruct.py convert FILE -f cif         => Read structure from FILE and output CIF file
                                              (Use convert --help to get list of formats supported)
  abistruct.py convert out_HIST.nc         => Read FINAL structure from the HIST file and
                                              print the corresponding ABINIT variables.
  abistruct.py supercell FILE -s 2 2 1     => Read structure from FILE and build [2, 2, 1] supercell,
                                              print new structure using --format (default abivars).
################
# K-points tools
################

  abistruct.py kpath FILE                  => Read structure from FILE and print ABINIT variables for k-path.
  abistruct.py bz FILE                     => Read structure from FILE, plot BZ with matplotlib.
  abistruct.py ngkpt FILE -n 4             => Compute `ngkpt` and `shiftk` from the number of divisions used to sample
                                              the smallest reciprocal lattice vector.
  abistruct.py abikmesh FILE --ngkpt 2 2 2 => Read structure from FILE, call Abinit to sample the BZ
                                              with a 2, 2, 2 k-mesh, print points in IBZ with weights.
  abistruct.py ktables FILE -m 2 2 2       => Read structure from FILE, call spglib to sample the BZ
                                              with a 2, 2, 2 k-mesh, print points in IBZ with weights.
  abistruct.py lgk FILE -k 0.25 0 0        => Read structure from FILE, find little group of k-point,
                                              print Bilbao character table.
  abistruct.py kstar FILE -k 0.25 0 0      => Read structure from FILE, print star of k-point.
  abistruct.py keq FILE -k 0.5 0 0 0 0.5 0  => Read structure from FILE, test whether k1 and k2 are
                                               symmetry-equivalent k-points.

###############
# Miscelleanous
###############

  abistruct.py neighbors FILE              => Get neighbors for each atom in the unit cell, out to a distance radius.
  abistruct.py interpolate FILE1 FILE2     => Interpolate between two structures. Useful for the construction of NEB inputs.
  abistruct.py xrd FILE                    => X-ray diffraction plot.
  abistruct.py oxistate FILE               => Estimate oxidation states with pymatgen bond valence analysis.

###############
# Visualization
###############

  abistruct.py visualize FILE -a vesta     => Visualize the structure with vesta (default if -a is not given)
                                              Supports also ovito, xcrysden, vtk, mayavi, matplotlib See --help
  abistruct.py ipython FILE                => Read structure from FILE and open it in the Ipython terminal.
  abistruct.py notebook FILE               => Read structure from FILE and generate jupyter notebook.

###########
# Databases
###########

  abistruct.py cod_id 1526507              => Get structure from COD database (http://www.crystallography.net/cod).
  abistruct.py cod_search MgB2             => Search for structures in the COD database.
  abistruct.py mp_id mp-149                => Get structure from materials project database and print
                                              JSON representation. Use e.g. `-f abivars` to change format.
  abistruct.py mp_match FILE               => Read structure from FILE and find matching structures on the
                                              Materials Project site. Use e.g. `-f cif` to change output format.
  abistruct.py mp_search LiF               => Connect to the materials project database. Get structures corresponding
                                              to a chemical system or formula e.g. `Fe2O3` or `Li-Fe-O` or
                                              `Ir-O-*` for wildcard pattern matching.
                                              Print info and Abinit input files. Use e.g. `-f POSCAR`
                                              to change output format. `-f None` to disable structure output.
  abistruct.py mp_pd FILE-or-elements      => Generate phase diagram with entries from the Materials Project.
                                              Accept FILE with structure or list of elements e.g `Li-Fe-O`

`FILE` is any file supported by abipy/pymatgen e.g Netcdf files, Abinit input/output, POSCAR, xsf ...
Use `abistruct.py --help` for help and `abistruct.py COMMAND --help` to get the documentation for `COMMAND`.
Use `-v` to increase verbosity level (can be supplied multiple times e.g -vv).
"""

def get_parser(with_epilog=False):

    # Parent parser for commands that need to know the filepath
    path_selector = argparse.ArgumentParser(add_help=False)
    path_selector.add_argument('filepath', nargs="?",
        help="File with the crystalline structure (Abinit Netcdf files, CIF, Abinit input/output files, POSCAR ...)")

    # Parent parser for commands supporting (jupyter notebooks)
    nb_parser = argparse.ArgumentParser(add_help=False)
    nb_parser.add_argument('-nb', '--notebook', default=False, action="store_true", help='Generate jupyter notebook.')
    nb_parser.add_argument('--foreground', action='store_true', default=False,
        help="Run jupyter notebook in the foreground.")

    parser = argparse.ArgumentParser(epilog=get_epilog() if with_epilog else "",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-V', '--version', action='version', version=abilab.__version__)

    # Parser for commands that need to call spglib.
    spgopt_parser = argparse.ArgumentParser(add_help=False)
    spgopt_parser.add_argument('--symprec', default=1e-3, type=float,
        help="""\
symprec (float): Tolerance for symmetry finding. Defaults to 1e-3,
which is fairly strict and works well for properly refined structures with atoms in the proper symmetry coordinates.
For structures with slight deviations from their proper atomic positions (e.g., structures relaxed with electronic structure
codes), a looser tolerance of 0.1 (the value used in Materials Project) is often needed.""")
    spgopt_parser.add_argument('--angle-tolerance', default=5.0, type=float,
        help="angle_tolerance (float): Angle tolerance for symmetry finding. Default: 5.0")
    spgopt_parser.add_argument("--no-time-reversal", default=False, action="store_true", help="Don't use time-reversal.")
    spgopt_parser.add_argument("--site-symmetry", default=False, action="store_true",
                               help="Show site symmetries i.e. the point group operations that leave the site invariant.")

    # Parent parser for common options.
    copts_parser = argparse.ArgumentParser(add_help=False)
    copts_parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
        help='verbose, can be supplied multiple times to increase verbosity')
    copts_parser.add_argument('--loglevel', default="ERROR", type=str,
        help="Set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    # Parent parser for commands that need to save structure to file.
    savefile_parser = argparse.ArgumentParser(add_help=False)
    savefile_parser.add_argument("-s", "--savefile", default="", type=str,
        help="Save final structure to file. Format is detected from file extensions "
             "e.g. out.abi for Abinit input, out.cif for CIF format.")

    # Helper functions to construct sub-parsers.
    def add_primitive_options(parser):
        """Add --no-primitive and --primitive-standard options to a parser."""
        group = parser.add_mutually_exclusive_group()
        group.add_argument('--no-primitive', default=False, action='store_true', help="Do not enforce primitive cell.")
        group.add_argument('--primitive-standard', default=False, action='store_true',
            help="Enforce primitive standard cell.")

    supported_formats = "(abivars, cif, xsf, poscar, qe, siesta, wannier90, cssr, json, None)"
    def add_format_arg(parser, default, option=True, formats=None):
        """Add --format option to a parser with default value `default`."""
        formats = supported_formats if formats is None else formats
        if option:
            parser.add_argument("-f", "--format", default=default, type=str,
                help="Output format. Default: %s. Accept: %s" % (default, formats))
        else:
            parser.add_argument('format', nargs="?", default=default, type=str,
                help="Output format. Default: %s. Accept: %s" % (default, formats))

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help',
        description="Valid subcommands, use command --help for help")

    # Subparser for spglib command.
    p_spglib = subparsers.add_parser('spglib', parents=[copts_parser, path_selector, spgopt_parser],
        help="Analyze structure with spglib.")

    # Subparser for abispg command.
    p_abispg = subparsers.add_parser('abispg', parents=[copts_parser, path_selector, savefile_parser],
        help="Extract/Compute Abinit space group from file with structure.")
    p_abispg.add_argument("-t", "--tolsym", type=float, default=None, help="""\
Gives the tolerance on the atomic positions (reduced coordinates), primitive vectors, or magnetization,
to be considered equivalent, thanks to symmetry operations. This value is used by ABINIT in the recognition of the set
of symmetries of the system, or the application of the symmetry operations to generate from a reduced set of atoms,
The internal default is 1e-8. Setting tolsym to a value larger than 1e-8 will make Abinit detect the spacegroup within
this tolerance and re-symmetrize the input structure. This option is useful if the structure has been taken from a CIF
file that does not have enough significant digits.""")
    p_abispg.add_argument("-d", "--diff-mode", type=str, default="table", choices=["table", "diff"],
        help="Select diff output format.")

    # Subparser for convert command.
    p_convert = subparsers.add_parser('convert', parents=[copts_parser, path_selector],
        help="Convert structure to the specified format.")
    add_format_arg(p_convert, default="cif")

    # Subparser for supercell command.
    p_supercell = subparsers.add_parser('supercell', parents=[copts_parser, path_selector],
        help="Generate supercell.")
    p_supercell.add_argument("-s", "--scaling_matrix", nargs="+", required=True, type=int,
        help="""\
scaling_matrix: A scaling matrix for transforming the lattice vectors.
Has to be all integers. Several options are possible:

    a. A full 3x3 scaling matrix defining the linear combination
       the old lattice vectors. E.g., -s 2,1,0 0,3,0, 0,0,1 generates a new structure with lattice vectors
       a' = 2a + b, b' = 3b, c' = c where a, b, and c are the lattice vectors of the original structure.
    b. An sequence of three scaling factors. E.g., -s 2, 1, 1 specifies that the supercell should
       have dimensions 2a x b x c.
    c. A number, which simply scales all lattice vectors by the same factor.
    """)
    add_format_arg(p_supercell, default="abivars")

    # Subparser for abisanitize
    p_abisanitize = subparsers.add_parser('abisanitize', parents=[copts_parser, path_selector, spgopt_parser, savefile_parser],
        help="Sanitize structure with abi_sanitize, compare structures and save result to file.")
    add_primitive_options(p_abisanitize)

    # Subparser for irefine
    p_irefine = subparsers.add_parser('irefine', parents=[copts_parser, path_selector, spgopt_parser],
        help="Refine structure with abi_sanitize iteratively, stop if target space group is obtained.")
    p_irefine.add_argument("--target-spgnum", required=True, type=int, help="Target space group number.")
    p_irefine.add_argument("--symprec-step", default=0.05, type=float, help='Increment for symprec.')
    p_irefine.add_argument("--angle-tolerance-step", default=0.0, type=float, help='Increment for angle_tolerance.')
    p_irefine.add_argument("--ntrial", default=50, type=int, help='Number of trials. Default 50.')
    add_primitive_options(p_irefine)

    # Subparser for conventional.
    p_conventional = subparsers.add_parser('conventional', parents=[copts_parser, path_selector, spgopt_parser, savefile_parser],
        help="Gives a structure with a conventional cell according to certain standards. "
             "The standards are defined in doi:10.1016/j.commatsci.2010.05.010")

    # Subparser for proto.
    p_proto = subparsers.add_parser('proto', parents=[copts_parser, path_selector],
        help=("Find prototype in the AFLOW LIBRARY OF CRYSTALLOGRAPHIC PROTOTYPES. "
              "http://doi.org/10.1016/j.commatsci.2017.01.017"))
    p_proto.add_argument("--ltol", default=0.2, type=float, help="fractional length tolerance.")
    p_proto.add_argument("--stol", default=0.3, type=float, help="site tolerance.")
    p_proto.add_argument("--angle-tol", default=5, type=float, help="angle tolerance.")

    # Subparser for wyckoff.
    p_wyckoff = subparsers.add_parser('wyckoff', parents=[copts_parser, spgopt_parser, path_selector],
            help="Print wyckoff positions. WARNING: still under development!")
    p_wyckoff.add_argument("--refine", default=False, action="store_true",
                            help="Use spglib to refine structure before computation")

    # Subparser for tensor_site.
    p_tensor_site = subparsers.add_parser('tensor_site', parents=[copts_parser, spgopt_parser, path_selector],
            help="Print symmetry properties of tensors due to site-symmetries. WARNING: still under development!")
    p_tensor_site.add_argument("--refine", default=False, action="store_true",
                                help="Use spglib to refine structure before computation")

    # Subparser for neighbors.
    p_neighbors = subparsers.add_parser('neighbors', parents=[copts_parser, path_selector],
                                        help="Get neighbors for each atom in the unit cell, out to a distance radius.")
    p_neighbors.add_argument("-r", "--radius", default=2, type=float, help="Radius of the sphere in Angstrom.")

    # Subparser for interpolate.
    p_interpolate = subparsers.add_parser('interpolate', parents=[copts_parser],
        help=("Interpolate between two structures. Useful for the construction of NEB inputs."))
    p_interpolate.add_argument("filepaths", nargs=2, help="Files with initial and final structures.")
    p_interpolate.add_argument("-n", "--nimages", default=10, type=int, help="No. of interpolation images. Defaults to 10.")
    p_interpolate.add_argument("--autosort_tol", default=0.5, type=float, help="""\
A distance tolerance in Angstrom in which to automatically sort end_structure to match to the
closest points in this particular structure. This is usually what you want in a NEB calculation.
0 implies no sorting. Otherwise, a 0.5 value (default) usually works pretty well.""")
    add_format_arg(p_interpolate, default="abivars")

    # Subparser for xrd.
    p_xrd = subparsers.add_parser('xrd', parents=[copts_parser, path_selector], help="X-ray diffraction plot.")
    p_xrd.add_argument("-w", "--wavelength", default="CuKa", type=str, help=(
        "The wavelength can be specified as a string. It must be one of the "
        "supported definitions in the WAVELENGTHS dict declared in pymatgen/analysis/diffraction/xrd.py."
        "Defaults to 'CuKa', i.e, Cu K_alpha radiation."))
    p_xrd.add_argument("-s", "--symprec", default=0, type=float, help=(
        "Symmetry precision for structure refinement. "
        "If set to 0, no refinement is done. Otherwise, refinement is performed using spglib with provided precision."))
    p_xrd.add_argument("-t", "--two-theta-range", default=(0, 90), nargs=2, help=(
        "Tuple for range of two_thetas to calculate in degrees. Defaults to (0, 90)."))
    p_xrd.add_argument("-nap", "--no-annotate-peaks", default=False, action="store_true",
        help="Whether to annotate the peaks with plane information.")

    # Subparser for oxistate.
    p_oxistate = subparsers.add_parser('oxistate', parents=[copts_parser, path_selector],
        help="Estimate oxidation states with pymatgen bond valence analysis.")
    # Subparser for ipython.
    p_ipython = subparsers.add_parser('ipython', parents=[copts_parser, path_selector],
        help="Open IPython shell for advanced operations on structure object.")
    # Subparser for notebook.
    p_notebook = subparsers.add_parser('notebook', parents=[copts_parser, path_selector],
        help="Read structure from file and generate jupyter notebook.")
    p_notebook.add_argument('--foreground', action='store_true', default=False,
        help="Run jupyter notebook in the foreground.")
    # Subparser for kpath.
    p_kpath = subparsers.add_parser('kpath', parents=[copts_parser, path_selector],
        help="Read structure from file, generate k-path for band-structure calculations.")
    add_format_arg(p_kpath, default="abinit", formats=["abinit", "wannier90", "siesta"])
    # Subparser for bz.
    p_bz = subparsers.add_parser('bz', parents=[copts_parser, path_selector],
        help="Read structure from file, plot Brillouin zone with matplotlib.")
    # Subparser for ngkpt.
    p_ngkpt = subparsers.add_parser('ngkpt', parents=[copts_parser, path_selector],
        help="Return the Abinit k-point sampling variables "
             "from the number of divisions used to sample the smallest "
             "lattice vector of the reciprocal lattice.")
    p_ngkpt.add_argument("-n", "--nksmall", required=True, type=int,
        help="Number of divisions used to sample the smallest reciprocal lattice vector.")
    # Subparser for ktables.
    p_ktables = subparsers.add_parser('ktables', parents=[copts_parser, path_selector],
        help=("Read structure from filepath, call spglib to sample the BZ, "
              "print k-points in the IBZ with weights."))
    p_ktables.add_argument("-m", "--mesh", nargs=3, required=True, type=int, help="Mesh divisions e.g. 2 3 4")
    p_ktables.add_argument("-s", "--is_shift", nargs="+", default=None, type=int,
        help=("Three integers (spglib API). The kmesh is shifted along " +
              "the axis in half of adjacent mesh points irrespective of the mesh numbers e.g. 1 1 1 " +
              "Default: Unshifted mesh."))
    p_ktables.add_argument("--no-time-reversal", default=False, action="store_true", help="Don't use time-reversal.")

    # Subparser for abikmesh.
    p_abikmesh = subparsers.add_parser('abikmesh', parents=[copts_parser, path_selector],
        help=("Read structure from file, call Abinit to sample the BZ with ngkpt, shiftk, and kptopt. "
              "Print k-points in the IBZ with weights."))
    p_abikmesh.add_argument("--kppa", default=None, type=int,
            help="Number of k-points per reciprocal atom. Mutually exclusive with ngkpt.")
    p_abikmesh.add_argument("--ngkpt", nargs=3, default=None, type=int,
            help="Mesh divisions e.g. 2 3 4")
    p_abikmesh.add_argument("--shiftk", nargs="+", default=(0.5, 0.5, 0.5), type=float,
       help="Kmesh shifts. Default: 0.5 0.5 0.5")
    p_abikmesh.add_argument("--kptopt", default=1, type=int,
            help="Kptopt input variable. Default: 1")

    # Subparser for kmesh_jhu.
    #p_kmesh_jhu = subparsers.add_parser('kmesh_jhu', parents=[copts_parser, path_selector], help="Foobar ")

    # Subparser for lgk.
    p_lgk = subparsers.add_parser('lgk', parents=[copts_parser, path_selector, spgopt_parser],
        help="Read structure from file, find little group of k-point, print Bilbao character table.")
    p_lgk.add_argument("-k", "--kpoint", nargs=3, required=True, type=float,
        help="K-point in reduced coordinates e.g. 0.25 0 0")

    # Subparser for kstar.
    p_kstar = subparsers.add_parser('kstar', parents=[copts_parser, path_selector, spgopt_parser],
        help="Read structure from file, print star of k-point.")
    p_kstar.add_argument("-k", "--kpoint", nargs=3, required=True, type=float,
        help="K-point in reduced coordinates e.g. 0.25 0 0")

    # Subparser for keq.
    p_keq = subparsers.add_parser('keq', parents=[copts_parser, path_selector, spgopt_parser],
        help="Read structure from file, check whether two k-points are equivalent by symmetry.")
    p_keq.add_argument("-k", "--kpoints", nargs=6, required=True, type=float,
        help="K-points in reduced coordinates e.g. 0.25 0 0 0 0.25 0")

    # Subparser for visualize command.
    p_visualize = subparsers.add_parser('visualize', parents=[copts_parser, path_selector],
        help=("Visualize the structure with the specified application. "
              "Requires external app or optional python modules (mayavi, vtk)."))
    p_visualize.add_argument("-a", "--appname", type=str, default="vesta",
        help=("Application name. Possible options: %s, mpl (matplotlib), mayavi, vtk" % ", ".join(Visualizer.all_visunames())))

    # Options for commands accessing the materials project database.
    mp_rest_parser = argparse.ArgumentParser(add_help=False)
    mp_rest_parser.add_argument("--mapi-key", default=None,
        help="Pymatgen PMG_MAPI_KEY. Use value in .pmgrc.yaml if not specified.")
    mp_rest_parser.add_argument("--endpoint", help="Pymatgen database.", default="https://www.materialsproject.org/rest/v2")
    mp_rest_parser.add_argument("-b", "--browser", default=False, action='store_true',
        help="Open materials-project webpages in browser")

    # Subparser for mp_id command.
    p_mpid = subparsers.add_parser('mp_id', parents=[copts_parser, mp_rest_parser],
        help="Get structure from the pymatgen database. Export to format. Requires internet connection and PMG_MAPI_KEY.")
    p_mpid.add_argument("mpid", type=str, default=None, help="Pymatgen identifier.")
    add_format_arg(p_mpid, default="cif")

    # Subparser for mp_match command.
    p_mpmatch = subparsers.add_parser('mp_match', parents=[path_selector, mp_rest_parser, copts_parser, nb_parser],
        help="Get structure from the pymatgen database. Requires internet connection and PMG_MAPI_KEY.")
    add_format_arg(p_mpmatch, default="abivars")

    # Subparser for mp_search command.
    p_mpsearch = subparsers.add_parser('mp_search', parents=[mp_rest_parser, copts_parser, nb_parser],
        help="Get structure from the pymatgen database. Requires internet connection and PMG_MAPI_KEY")
    p_mpsearch.add_argument("chemsys_formula_id", type=str, default=None,
        help="A chemical system (e.g., Li-Fe-O), or formula (e.g., Fe2O3) or materials_id (e.g., mp-1234).")
    p_mpsearch.add_argument("-s", "--select-spgnum", type=int, default=None,
        help="Select structures with this space group number.")
    add_format_arg(p_mpsearch, default="abivars")

    # Subparser for mp_pd command.
    p_mp_pda = subparsers.add_parser('mp_pd', parents=[mp_rest_parser, copts_parser],
        help=("Generate phase diagram with entries from the Materials Project. "
              "Requires internet connection and PMG_MAPI_KEY"))
    p_mp_pda.add_argument("file_or_elements", type=str, default=None,
        help="FILE with structure or elements e.g., Li-Fe-O).")
    p_mp_pda.add_argument("-u", "--show-unstable", type=int, default=0,
        help="""Whether unstable phases will be plotted as
well as red crosses. If a number > 0 is entered, all phases with
ehull < show_unstable will be shown.""")

    # Subparser for cod_search command.
    p_codsearch = subparsers.add_parser('cod_search', parents=[copts_parser, nb_parser],
        help="Get structure from COD database. Requires internet connection and mysql")
    p_codsearch.add_argument("formula", type=str, default=None, help="formula (e.g., Fe2O3).")
    p_codsearch.add_argument("-s", "--select-spgnum", type=int, default=None,
        help="Select structures with this space group number.")
    p_codsearch.add_argument('--primitive', default=False, action='store_true',
        help="Convert COD cells into primitive cells.")
    add_format_arg(p_codsearch, default="abivars")

    # Subparser for cod_id command.
    p_codid = subparsers.add_parser('cod_id', parents=[copts_parser],
        help="Get structure from COD database. Requires internet connection and mysql")
    p_codid.add_argument("cod_identifier", type=int, default=None, help="COD identifier.")
    p_codid.add_argument('--primitive', default=False, action='store_true', help="Convert COD cell into primitive cell.")
    add_format_arg(p_codid, default="abivars")

    # Subparser for animate command.
    p_animate = subparsers.add_parser('animate', parents=[copts_parser, path_selector],
        help="Read structures from HIST.nc or XDATCAR. Print structures in Xcrysden AXSF format to stdout.")

    return parser

@prof_main
def main():

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(get_epilog())
        if err_msg: sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    parser = get_parser(with_epilog=True)

    # Parse command line.
    try:
        options = parser.parse_args()
    except Exception as exc:
        show_examples_and_exit(error_code=1)

    if not options.command:
        show_examples_and_exit(error_code=1)

    # loglevel is bound to the string value obtained from the command line argument.
    # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
    import logging
    numeric_level = getattr(logging, options.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.loglevel)
    logging.basicConfig(level=numeric_level)

    if options.verbose > 2:
        print(options)

    if options.command == "spglib":
        structure = abilab.Structure.from_file(options.filepath)
        print(structure.spget_summary(symprec=options.symprec, angle_tolerance=options.angle_tolerance,
                                      site_symmetry=options.site_symmetry, verbose=options.verbose))

    elif options.command == "abispg":
        structure = abilab.Structure.from_file(options.filepath)
        check_ordered_structure(structure)
        abi_spg = structure.abi_spacegroup

        if abi_spg is not None and options.tolsym is None:
            print(structure.spget_summary(verbose=options.verbose))
        else:
            # Here we compare Abinit wrt spglib. If abi_spg is None, we create a temporary
            # task to run the code in dry-run mode.
            if abi_spg is None:
                print("FILE does not contain Abinit symmetry operations.")
            cprint("Calling Abinit in --dry-run mode with chkprim = 0 to get space group.")
            if options.tolsym is not None and options.tolsym > 1e-8:
                cprint("Crystal structure will be re-symmetrized by Abinit with tolsym: %s" % options.tolsym, "yellow")

            from abipy.data.hgh_pseudos import HGH_TABLE
            gsinp = factories.gs_input(structure, HGH_TABLE, spin_mode="unpolarized")
            gsinp["chkprim"] = 0
            abistructure = gsinp.abiget_spacegroup(tolsym=options.tolsym)
            print(abistructure.spget_summary(verbose=options.verbose))
            print("")

            diff_structures([structure, abistructure], mode=options.diff_mode,
                            headers=["Input structure", "After Abinit symmetrization"], fmt="abivars")

            # Save file.
            save_structure(abistructure, options)

    elif options.command == "convert":
        fmt = options.format
        if fmt == "cif" and options.filepath.endswith(".cif"): fmt = "abivars"
        print(abilab.Structure.from_file(options.filepath).convert(fmt=fmt))

    elif options.command == "supercell":
        structure = abilab.Structure.from_file(options.filepath)

        options.scaling_matrix = np.array(options.scaling_matrix)
        if len(options.scaling_matrix) == 9:
            options.scaling_matrix.shape = (3, 3)
        if options.verbose:
            print("scaling matrix: ", options.scaling_matrix)

        supcell = structure * options.scaling_matrix
        #supcell = structure.make_supercell(scaling_matrix, to_unit_cell=True)
        print(supcell.convert(fmt=options.format))

    elif options.command == "abisanitize":
        print("\nCalling abi_sanitize to get a new structure in which:")
        print("    * Structure is refined.")
        print("    * Reduced to primitive settings.")
        print("    * Lattice vectors are exchanged if the triple product is negative\n")

        structure = abilab.Structure.from_file(options.filepath)
        sanitized = structure.abi_sanitize(symprec=options.symprec, angle_tolerance=options.angle_tolerance,
                                           primitive=not options.no_primitive, primitive_standard=options.primitive_standard)
        index = [options.filepath, "abisanitized"]
        dfs = abilab.dataframes_from_structures([structure, sanitized], index=index, with_spglib=True)

        abilab.print_dataframe(dfs.lattice, title="Lattice parameters:")
        abilab.print_dataframe(dfs.coords, title="Atomic positions (columns give the site index):")

        if not options.verbose:
            print("\nUse -v for more info")
            #print(sanitized.convert(fmt="cif"))
        else:
            #print("\nDifference between structures:")
            if len(structure) == len(sanitized):
                table = []
                for line1, line2 in zip(str(structure).splitlines(), str(sanitized).splitlines()):
                    table.append([line1, line2])
                print(str(tabulate(table, headers=["Initial structure", "Abisanitized"])))

            else:
                print("\nInitial structure:")
                print(structure)
                print("\nabisanitized structure:")
                print(sanitized)

        # Save file.
        save_structure(sanitized, options)

    elif options.command == "irefine":
        structure = abilab.Structure.from_file(options.filepath)
        sanitized = structure.copy()
        symprec, angle_tolerance = options.symprec, options.angle_tolerance
        print("Calling abi_sanitize with increasing tolerances to reach target space group:", options.target_spgnum)
        print("Using symprec_step: ", options.symprec_step, ", angle_tolerance_step:", options.angle_tolerance_step,
              "ntrial", options.ntrial)
        itrial = 0
        while itrial < options.ntrial:
            print(">>> Trying with symprec: %s, angle_tolerance: %s" % (symprec, angle_tolerance))
            sanitized = sanitized.abi_sanitize(symprec=symprec, angle_tolerance=angle_tolerance,
                primitive=not options.no_primitive, primitive_standard=options.primitive_standard)
            spg_symb, spg_num = sanitized.get_space_group_info(symprec=symprec, angle_tolerance=angle_tolerance)
            print(">>> Space-group number:", spg_symb, ", symbol:", spg_num, "for trial:", itrial)
            if spg_num == options.target_spgnum:
                print(2 * "\n", "# Final structure with space group number:", spg_symb, ", symbol:", spg_num, 2 *"\n")
                print(sanitized.convert(fmt="cif"))
                break

            # Increment counter and tols.
            itrial += 1
            symprec += options.symprec_step
            angle_tolerance += options.angle_tolerance_step
        else:
            print("Cannot find space group number:", options.target_spgnum, "after", options.ntrial, "iterations")
            return 1

        # Save file.
        #save_structure(sanitized, options)

    elif options.command == "conventional":
        print("\nCalling get_conventional_standard_structure to get conventional structure:")
        print("The standards are defined in Setyawan, W., & Curtarolo, S. (2010). ")
        print("High-throughput electronic band structure calculations: Challenges and tools. ")
        print("Computational Materials Science, 49(2), 299-312. doi:10.1016/j.commatsci.2010.05.010\n")

        structure = abilab.Structure.from_file(options.filepath)
        conv = structure.get_conventional_standard_structure(international_monoclinic=True,
                                           symprec=options.symprec, angle_tolerance=options.angle_tolerance)
        index = [options.filepath, "conventional"]
        dfs = abilab.dataframes_from_structures([structure, conv], index=index, with_spglib=True)

        abilab.print_dataframe(dfs.lattice, title="Lattice parameters:")
        if options.verbose:
            abilab.print_dataframe(dfs.coords, title="Atomic positions (columns give the site index):")

        if not options.verbose:
            print("\nUse -v for more info")
        else:
            #print("\nDifference between structures:")
            if len(structure) == len(conv):
                table = []
                for line1, line2 in zip(str(structure).splitlines(), str(conv).splitlines()):
                    table.append([line1, line2])
                print(str(tabulate(table, headers=["Initial structure", "Conventional"])))

            else:
                print("\nInitial structure:\n", structure)
                print("\nConventional structure:\n", conv)

        # Save file.
        save_structure(conv, options)

    elif options.command == "proto":
        structure = abilab.Structure.from_file(options.filepath)
        from pymatgen.analysis.aflow_prototypes import AflowPrototypeMatcher
        m = AflowPrototypeMatcher(initial_ltol=options.ltol, initial_stol=options.stol,
                                  initial_angle_tol=options.angle_tol)
        dlist = m.get_prototypes(structure)
        if not dlist:
            cprint("Cannot find AFLOW prototype for structure.")
            print(structure.to_string(verbose=options.verbose))
            return 1
        else:
            cprint("Found %d matches" % len(dlist), "green")
            for d in dlist:
                if "snl" in d:
                    snl = d.pop("snl")
                    if options.verbose: pprint(snl.as_dict())
                pprint(d)
                url = "http://aflow.org/CrystalDatabase/%s.html" % d["tags"]["aflow"]
                print("AFLOW url: %s\n" % url)
            if not options.verbose:
                print("Use --verbose to increase output level")

    elif options.command == "wyckoff":
        structure = abilab.Structure.from_file(options.filepath)
        if options.refine:
            print("Refining structure with symprec: %s, angle_tolerance: %s" % (options.symprec, options.angle_tolerance))
            structure = structure.refine(symprec=options.symprec, angle_tolerance=options.angle_tolerance)
        print(structure.spget_summary(verbose=options.verbose))
        ss = structure.site_symmetries
        df = ss.get_wyckoff_dataframe(verbose=options.verbose)
        abilab.print_dataframe(df, title="\nWyckoff positions in reduced coordinates.")

    elif options.command == "tensor_site":
        structure = abilab.Structure.from_file(options.filepath)
        if options.refine:
            print("Refining structure with symprec: %s, angle_tolerance: %s" % (options.symprec, options.angle_tolerance))
            structure = structure.refine(symprec=options.symprec, angle_tolerance=options.angle_tolerance)
        print(structure.spget_summary(verbose=options.verbose))
        ss = structure.site_symmetries
        df = ss.get_tensor_rank2_dataframe(verbose=options.verbose)
        abilab.print_dataframe(df, title="\nTensor components in reduced coordinates (rank 2, symmetric)")

    elif options.command == "neighbors":
        abilab.Structure.from_file(options.filepath).print_neighbors(radius=options.radius)

    elif options.command == "interpolate":
        initial_structure = abilab.Structure.from_file(options.filepaths[0])
        end_structure = abilab.Structure.from_file(options.filepaths[1])
        structures = initial_structure.interpolate(end_structure, nimages=options.nimages,
                                                   interpolate_lattices=False, pbc=True,
                                                   autosort_tol=options.autosort_tol)
        structures = list(map(abilab.Structure.as_structure, structures))
        for i, s in enumerate(structures):
            print(marquee("Structure #%d" % i, mark="="))
            print(s.convert(fmt=options.format))
            print(" ")

    elif options.command == "xrd":
        structure = abilab.Structure.from_file(options.filepath)
        two_theta_range = tuple(float(t) for t in options.two_theta_range)
        structure.plot_xrd(wavelength=options.wavelength, two_theta_range=two_theta_range,
                           symprec=options.symprec, annotate_peaks=not options.no_annotate_peaks)

    elif options.command == "oxistate":
        print(abilab.Structure.from_file(options.filepath).get_oxi_state_decorated())

    elif options.command == "ipython":
        structure = abilab.Structure.from_file(options.filepath)
        print("Invoking Ipython, `structure` object will be available in the Ipython terminal")
        import IPython
        IPython.start_ipython(argv=[], user_ns={"structure": structure})

    elif options.command == "notebook":
        structure = abilab.Structure.from_file(options.filepath)
        structure.make_and_open_notebook(nbpath=None, foreground=options.foreground)

    elif options.command == "visualize":
        structure = abilab.Structure.from_file(options.filepath)
        print(structure)
        print("Visualizing structure with:", options.appname)
        structure.visualize(appname=options.appname)

    elif options.command == "kpath":
        structure = abilab.Structure.from_file(options.filepath)
        print(structure.get_kpath_input_string(fmt=options.format, line_density=10))

    elif options.command == "bz":
        abilab.Structure.from_file(options.filepath).plot_bz()

    elif options.command == "ngkpt":
        d = abilab.Structure.from_file(options.filepath).calc_ksampling(options.nksmall)
        print("ngkpt %d %d %d" % (d.ngkpt[0], d.ngkpt[1], d.ngkpt[2]))
        print("nshiftk ", len(d.shiftk), "\nshiftk")
        for s in d.shiftk:
            print("  %s %s %s" % (s[0], s[1], s[2]))

    elif options.command == "ktables":
        structure = abilab.Structure.from_file(options.filepath)
        k = Ktables(structure, options.mesh, options.is_shift, not options.no_time_reversal)
        print(k)
        print("")
        print("NB: These results are obtained by calling spglib with the structure read from file.")
        print("The k-points might differ from the ones expected by Abinit, especially if the space groups differ.")

        if not options.verbose:
            print("\nUse -v to obtain the BZ --> IBZ mapping.")
        else:
            print()
            k.print_bz2ibz()

    elif options.command == "abikmesh":
        structure = abilab.Structure.from_file(options.filepath)
        if options.kppa is None and options.ngkpt is None:
            raise ValueError("Either ngkpt or kppa must be provided")

        if options.kppa is not None:
            print("Calling Abinit to compute the IBZ with kppa:", options.kppa, "and shiftk:", options.shiftk)
            ibz = IrredZone.from_kppa(structure, options.kppa, options.shiftk,
                                      kptopt=options.kptopt, verbose=options.verbose)
        else:
            print("Calling Abinit to compute the IBZ with ngkpt:", options.ngkpt, "and shiftk:", options.shiftk)
            ibz = IrredZone.from_ngkpt(structure, options.ngkpt, options.shiftk,
                                       kptopt=options.kptopt, verbose=options.verbose)

        print(ibz.to_string(verbose=options.verbose))

    #elif options.command == "kmesh_jhu":
    #    structure = abilab.Structure.from_file(options.filepath)
    #    from pymatgen.ext.jhu import get_kpoints
    #    kpoints = get_kpoints(structure, min_distance=0, min_total_kpoints=1,
    #                           kppra=None, gap_distance=7, remove_symmetry=None,
    #                           include_gamma="auto", header="simple", incar=None)
    #    #print(kpoints)

    elif options.command == "lgk":
        structure = abilab.Structure.from_file(options.filepath)
        spgrp = structure.abi_spacegroup
        if spgrp is None:
            cprint("Your file does not contain Abinit symmetry operations.", "yellow")
            cprint("Will call spglib to obtain the space group (assuming time-reversal: %s)" %
                   (not options.no_time_reversal), "yellow")
            spgrp = AbinitSpaceGroup.from_structure(structure, has_timerev=not options.no_time_reversal,
                        symprec=options.symprec, angle_tolerance=options.angle_tolerance)
        print()
        print(marquee("Structure", mark="="))
        print(structure.spget_summary(verbose=options.verbose))
        print("\n")

        print(marquee("Little Group", mark="="))
        ltk = spgrp.find_little_group(kpoint=options.kpoint)
        print(ltk.to_string(verbose=options.verbose))

    elif options.command == "kstar":
        structure = abilab.Structure.from_file(options.filepath)

        # Call spglib to get spacegroup if Abinit spacegroup is not available.
        if structure.abi_spacegroup is None:
            structure.spgset_abi_spacegroup(has_timerev=not options.no_time_reversal)

        kpoint = Kpoint(options.kpoint, structure.reciprocal_lattice)
        kstar = kpoint.compute_star(structure.abi_spacegroup, wrap_tows=True)
        print("Found %s points in the star of %s\n" % (len(kstar), repr(kpoint)))
        for k in kstar:
            print(4 * " ", repr(k))

    elif options.command == "keq":
        structure = abilab.Structure.from_file(options.filepath)

        # Call spglib to get spacegroup if Abinit spacegroup is not available.
        if structure.abi_spacegroup is None:
            structure.spgset_abi_spacegroup(has_timerev=not options.no_time_reversal)

        k1, k2 = options.kpoints[:3], options.kpoints[3:6]
        k1tab = structure.abi_spacegroup.symeq(k1, k2)

        if k1tab.isym != -1:
            print("\nk1:", k1, "and k2:", k2, "are symmetry equivalent k-points\n")
            print("Related by the symmetry operation (reduced coords):\n", k1tab.op)
            print("With umklapp vector Go = TO(k1) - k2 =", k1tab.g0)
        else:
            print(k1, "and", k2, "are NOT symmetry equivalent")

    elif options.command == "mp_id":
        # Get the Structure corresponding to material_id.
        structure = abilab.Structure.from_mpid(options.mpid, final=True,
                                               api_key=options.mapi_key, endpoint=options.endpoint)
        # Convert to format and print it.
        print(structure.convert(fmt=options.format))

    elif options.command == "mp_match":
        mp = abilab.mp_match_structure(options.filepath)
        if not mp.structures:
            cprint("No structure found in database", "yellow")
            return 1

        if options.notebook:
            return mp.make_and_open_notebook(foreground=options.foreground)
        else:
            mp.print_results(fmt=options.format, verbose=options.verbose)

        if options.browser:
            mp.open_browser(limit=None if options.verbose == 2 else 10)

    elif options.command == "mp_search":
        mp = abilab.mp_search(options.chemsys_formula_id)
        if not mp.structures:
            cprint("No structure found in Materials Project database", "yellow")
            return 1
        if options.select_spgnum: mp = mp.filter_by_spgnum(options.select_spgnum)

        if options.notebook:
            return mp.make_and_open_notebook(foreground=options.foreground)
        else:
            mp.print_results(fmt=options.format, verbose=options.verbose)

        if options.browser:
            mp.open_browser(limit=None if options.verbose == 2 else 10)

    elif options.command == "mp_pd":
        if os.path.exists(options.file_or_elements):
            structure = abilab.Structure.from_file(options.file_or_elements)
            elements = structure.symbol_set
        else:
            elements = options.file_or_elements.split("-")

        if options.verbose > 1: print("Building phase-diagram for elements:", elements)
        with abilab.restapi.get_mprester(api_key=options.mapi_key, endpoint=options.endpoint) as rest:
            pdr = rest.get_phasediagram_results(elements)
            pdr.print_dataframes(verbose=options.verbose)
            pdr.plot(show_unstable=options.show_unstable)

    elif options.command == "cod_search":
        cod = abilab.cod_search(options.formula, primitive=options.primitive)
        if not cod.structures:
            cprint("No structure found in COD database", "yellow")
            return 1
        if options.select_spgnum: cod = cod.filter_by_spgnum(options.select_spgnum)

        if options.notebook:
            return cod.make_and_open_notebook(foreground=options.foreground)
        else:
            cod.print_results(fmt=options.format, verbose=options.verbose)

    elif options.command == "cod_id":
        # Get the Structure from COD
        structure = abilab.Structure.from_cod_id(options.cod_identifier, primitive=options.primitive)
        # Convert to format and print it.
        print(structure.convert(fmt=options.format))

    elif options.command == "animate":
        filepath = options.filepath
        if any(filepath.endswith(ext) for ext in ("HIST", "HIST.nc")):
            with abilab.abiopen(filepath) as hist:
                structures = hist.structures

        elif "XDATCAR" in filepath:
            structures = Xdatcar(filepath).structures
            if not structures:
                raise RuntimeError("Your Xdatcar contains only one structure. Due to a bug "
                    "in the pymatgen routine, your structures won't be parsed correctly"
                    "Solution: Add another structure at the end of the file.")
        else:
            raise ValueError("Don't know how to handle file %s" % filepath)

        xsf_write_structure(sys.stdout, structures)

    else:
        raise ValueError("Unsupported command: %s" % options.command)

    return 0


if __name__ == "__main__":
    sys.exit(main())
