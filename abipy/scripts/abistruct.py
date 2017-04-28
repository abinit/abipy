#!/usr/bin/env python
"""
Script to analyze/export/visualize the crystal structure saved in the netcdf files produced by ABINIT.
"""
from __future__ import unicode_literals, division, print_function, absolute_import

import sys
import os
import argparse
import numpy as np

from monty.string import marquee
from monty.functools import prof_main
from monty.termcolor import cprint
from pprint import pprint
from tabulate import tabulate
from pymatgen.io.vasp.outputs import Xdatcar
from abipy import abilab
from abipy.core.symmetries import AbinitSpaceGroup
from abipy.iotools.visualizer import Visualizer
from abipy.iotools.xsf import xsf_write_structure
from abipy.core.kpoints import Ktables


@prof_main
def main():

    def str_examples():
        return """\
Usage example:

    abistruct.py spglib FILE                => Read structure from FILE and analyze it with spglib.
    abistruct.py abispg FILE                => Read structure from FILE, extract ABINIT space-group info.
    abistruct.py convert FILE cif           => Read the structure from FILE and print CIF file.
    abistruct.py convert FILE abivars       => Print the ABINIT variables defining the structure.
    abistruct.py convert out_HIST abivars   => Read the last structure from the HIST file and
                                               print the corresponding ABINIT variables.
    abistruct.py supercell FILE -s 2 2 1    => Read structure from FILE and build [2, 2, 1] supercell,
                                               print new structure using --format (default abivars).
    abistruct.py kpath FILE                 => Read structure from FILE and print ABINIT variables for k-path.
    abistruct.py bz FILE                    => Read structure from FILE, plot BZ with matplotlib.
    abistruct.py kmesh FILE -m 2 2 2        => Read structure from FILE, call spglib to sample the BZ with a 2,2,2 mesh,
                                               print points in IBZ with weights.
    abistruct.py lgk FILE -k 0.25 0 0       => Read structure from FILE, find little group of k-point,
                                               print Bilbao character table.
    abistruct.py abisanitize FILE           => Read structure from FILE, call abisanitize, compare structures and save
                                               "abisanitized" structure to file.
    abistruct.py conventional FILE          => Read structure from FILE, generate conventional structure
                                               following doi:10.1016/j.commatsci.2010.05.010
    abistruct.py visualize FILE vesta       => Visualize the structure with e.g. vesta (xcrysden, ...)
    abistruct.py ipython FILE               => Read structure from FILE and open Ipython terminal.
    abistruct.py notebook FILE              => Read structure from FILE and generate jupyter notebook.
    abistruct.py pmgdata mp-149             => Get structure from pymatgen database and print its JSON representation.
                                               Use e.g. `-f abivars` to change format.

Use `abistruct.py COMMAND --help` for futher information.
"""

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg: sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    # Parent parser for commands that need to know the filepath
    path_selector = argparse.ArgumentParser(add_help=False)
    path_selector.add_argument('filepath', nargs="?", help="File with the crystalline structure (netcdf, cif, input files ...)")

    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-V', '--version', action='version', version=abilab.__version__)

    spgopt_parser = argparse.ArgumentParser(add_help=False)
    spgopt_parser.add_argument('--symprec', default=1e-3, type=float,
        help="""\
symprec (float): Tolerance for symmetry finding. Defaults to 1e-3,
    which is fairly strict and works well for properly refined
    structures with atoms in the proper symmetry coordinates. For
    structures with slight deviations from their proper atomic
    positions (e.g., structures relaxed with electronic structure
    codes), a looser tolerance of 0.1 (the value used in Materials
    Project) is often needed.""")
    spgopt_parser.add_argument('--angle-tolerance', default=5.0, type=float,
        help="angle_tolerance (float): Angle tolerance for symmetry finding. Default: 5.0")
    spgopt_parser.add_argument("--no-time-reversal", default=False, action="store_true", help="Don't use time-reversal.")

    # Parent parser for common options.
    copts_parser = argparse.ArgumentParser(add_help=False)
    copts_parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                              help='verbose, can be supplied multiple times to increase verbosity')
    copts_parser.add_argument('--loglevel', default="ERROR", type=str,
                              help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help',
                                       description="Valid subcommands, use command --help for help")

    # Subparser for spglib command.
    p_spglib = subparsers.add_parser('spglib', parents=[copts_parser, path_selector, spgopt_parser],
                                      help="Analyze structure with spglib.")

    # Subparser for abispg command.
    p_abispg = subparsers.add_parser('abispg', parents=[copts_parser, path_selector],
                                      help="Extract Abinit spacegroup info from fileb.")

    # Subparser for convert command.
    p_convert = subparsers.add_parser('convert', parents=[copts_parser, path_selector],
                                      help="Convert structure to the specified format.")
    p_convert.add_argument('format', nargs="?", default="cif", type=str,
                           help="Format of the output file (cif, cssr, POSCAR, json, mson, abivars).")

    # Subparser for supercell command.
    p_supercell = subparsers.add_parser('supercell', parents=[copts_parser, path_selector],
                                        help="Generate supercell.")
    p_supercell.add_argument("-s", '--scaling_matrix', nargs="+", required=True, type=int,
                             help="""\
scaling_matrix: A scaling matrix for transforming the lattice vectors.
Has to be all integers. Several options are possible:

    a. A full 3x3 scaling matrix defining the linear combination
       the old lattice vectors. E.g., -s 2,1,0 0,3,0, 0,0,1
       generates a new structure with lattice vectors
       a' = 2a + b, b' = 3b, c' = c where a, b, and c are the lattice vectors of the original structure.
    b. An sequence of three scaling factors. E.g., -s 2, 1, 1
       specifies that the supercell should have dimensions 2a x b x c.
    c. A number, which simply scales all lattice vectors by the same factor.
    """)

    p_supercell.add_argument("-f", '--format', default="abivars", type=str,
                             help="Format of the output file (cif, cssr, POSCAR, json, mson, abivars).")

    # Subparser for abisanitize
    p_abisanitize = subparsers.add_parser('abisanitize', parents=[copts_parser, path_selector, spgopt_parser],
                                      help="Sanitize structure with abi_sanitize, compare structures and save result to file.")
    p_abisanitize.add_argument("--savefile", default="", type=str,
                               help='Save final structure to file. Format is detected from file extensions e.g. Si.cif')

    # Subparser for conventional.
    p_conventional = subparsers.add_parser('conventional', parents=[copts_parser, path_selector, spgopt_parser],
                                           help="Gives a structure with a conventional cell according to certain standards. "
                                                "The standards are defined in doi:10.1016/j.commatsci.2010.05.010")
    p_conventional.add_argument("--savefile", default="", type=str,
                                help='Save final structure to file. Format is detected from file extensions e.g. Si.cif')

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

    # Subparser for bz.
    p_bz = subparsers.add_parser('bz', parents=[copts_parser, path_selector],
                                 help="Read structure from file, plot Brillouin zone with matplotlib.")

    # Subparser for kmesh.
    p_kmesh = subparsers.add_parser('kmesh', parents=[copts_parser, path_selector],
                                    help=("Read structure from filepath, call spglib to sample the BZ,"
                                           "print k-points in the IBZ with weights."))
    p_kmesh.add_argument("-m", "--mesh", nargs="+", required=True, type=int, help="Mesh divisions e.g. 2 3 4")
    p_kmesh.add_argument("-s", "--is_shift", nargs="+", default=None, type=int,
                         help=("Three integers (spglib API). The kmesh is shifted along " +
                               "the axis in half of adjacent mesh points irrespective of the mesh numbers e.g. 1 1 1 " +
                               "Default: Unshited mesh."))
    p_kmesh.add_argument("--no-time-reversal", default=False, action="store_true", help="Don't use time-reversal.")

    # Subparser for lgk.
    p_lgk = subparsers.add_parser('lgk', parents=[copts_parser, path_selector, spgopt_parser],
                                  help="Read structure from file, find little group of k-point, print Bilbao character table.")
    p_lgk.add_argument("-k", "--kpoint", nargs="+", required=True, type=float,
                       help="K-point in reduced coordinates e.g. 0.25 0 0")

    # Subparser for visualize command.
    p_visualize = subparsers.add_parser('visualize', parents=[copts_parser, path_selector],
                                        help="Visualize the structure with the specified visualizer")
    p_visualize.add_argument('visualizer', nargs="?", default="vesta", type=str, help=("Visualizer name. "
        "List of visualizer supported: %s" % ", ".join(Visualizer.all_visunames())))

    # Subparser for pmgdata command.
    p_pmgdata = subparsers.add_parser('pmgdata', parents=[copts_parser],
                                      help="Get structure from the pymatgen database. Requires internet connection and MAPI_KEY")
    p_pmgdata.add_argument("pmgid", type=str, default=None, help="Pymatgen identifier")
    p_pmgdata.add_argument("--mapi-key", default=None,
                           help="Pymatgen MAPI_KEY. Use value in .pmgrc.yaml if not specified.")
    p_pmgdata.add_argument("--endpoint", help="Pymatgen database.",
                           default="https://www.materialsproject.org/rest/v2")
    p_pmgdata.add_argument("-f", '--format', default="json", type=str,
                           help="Format of the output file (cif, cssr, POSCAR, json, mson, abivars).")

    # Subparser for animate command.
    p_animate = subparsers.add_parser('animate', parents=[copts_parser, path_selector],
        help="Read structures from HIST or XDATCAR. Print structures in Xrysden AXSF format to stdout")

    # Parse command line.
    try:
        options = parser.parse_args()
    except Exception as exc:
        show_examples_and_exit(error_code=1)

    # loglevel is bound to the string value obtained from the command line argument.
    # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
    import logging
    numeric_level = getattr(logging, options.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.loglevel)
    logging.basicConfig(level=numeric_level)

    if options.command == "spglib":
        structure = abilab.Structure.from_file(options.filepath)
        print(structure.spget_summary(verbose=options.verbose))

    elif options.command == "abispg":
        structure = abilab.Structure.from_file(options.filepath)
        spgrp = structure.abi_spacegroup
        if spgrp is None:
            cprint("Your file does not contain Abinit symmetry operations.", "red")
            return 1
        print(spgrp)

    elif options.command == "convert":
        structure = abilab.Structure.from_file(options.filepath)

        if options.format == "abivars":
            print(structure.abi_string)
        else:
            print(structure.convert(fmt=options.format))

    elif options.command == "supercell":
        structure = abilab.Structure.from_file(options.filepath)

        options.scaling_matrix = np.array(options.scaling_matrix)
        if len(options.scaling_matrix) == 9:
            options.scaling_matrix.shape = (3, 3)
        if options.verbose:
            print("scaling matrix: ", options.scaling_matrix)

        supcell = structure * options.scaling_matrix
        #supcell = structure.make_supercell(scaling_matrix, to_unit_cell=True)

        if options.format == "abivars":
            print(supcell.abi_string)
        else:
            print(supcell.convert(fmt=options.format))

    elif options.command == "abisanitize":
        print("\nCalling abi_sanitize to get a new structure in which:")
        print("    * Structure is refined.")
        print("    * Reduced to primitive settings.")
        print("    * Lattice vectors are exchanged if the triple product is negative\n")

        structure = abilab.Structure.from_file(options.filepath)
        sanitized = structure.abi_sanitize(symprec=options.symprec, angle_tolerance=options.angle_tolerance,
                                           primitive=True, primitive_standard=False)
        index = [options.filepath, "abisanitized"]
        dfs = abilab.frames_from_structures([structure, sanitized], index=index, with_spglib=True)

        abilab.print_frame(dfs.lattice, title="Lattice parameters:")
        abilab.print_frame(dfs.coords, title="Atomic positions (columns give the site index):")

        if not options.verbose:
            print("\nUse -v for more info")
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

        # save file.
        if options.savefile:
            print("Saving abisanitized structure as %s" % options.savefile)
            if os.path.exists(options.savefile):
                raise RuntimeError("%s already exists. Cannot overwrite" % options.savefile)
            sanitized.to(filename=options.savefile)

    elif options.command == "conventional":
        print("\nCalling get_conventional_standard_structure to get conventional structure:")
        print("The standards are defined in Setyawan, W., & Curtarolo, S. (2010). ")
        print("High-throughput electronic band structure calculations: Challenges and tools. ")
        print("Computational Materials Science, 49(2), 299-312. doi:10.1016/j.commatsci.2010.05.010\n")

        structure = abilab.Structure.from_file(options.filepath)
        conv = structure.get_conventional_standard_structure(international_monoclinic=True,
                                           symprec=options.symprec, angle_tolerance=options.angle_tolerance)
        index = [options.filepath, "conventional"]
        dfs = abilab.frames_from_structures([structure, conv], index=index, with_spglib=True)

        abilab.print_frame(dfs.lattice, title="Lattice parameters:")
        #abilab.print_frame(dfs.coords, title="Atomic positions (columns give the site index):")

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
                print("\nInitial structure:")
                print(structure)
                print("\nConventional structure:")
                print(conv)

        # save file.
        if options.savefile:
            print("Saving conventional structure as %s" % options.savefile)
            if os.path.exists(options.savefile):
                raise RuntimeError("%s already exists. Cannot overwrite" % options.savefile)
            conv.to(filename=options.savefile)

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
        print("Visualizing structure with:", options.visualizer)
        structure.visualize(options.visualizer)

    elif options.command == "kpath":
        structure = abilab.Structure.from_file(options.filepath)
        print("# Abinit Structure")
        print(structure.abi_string)
        print("\n# K-path in reduced coordinates:")
        print("# tolwfr 1e-20 iscf -2 getden ??")
        print(" ndivsm 10")
        print(" kptopt", -(len(structure.hsym_kpoints) - 1))
        print(" kptbounds")
        for k in structure.hsym_kpoints:
            print("    %.5f  %.5f  %.5f" % tuple(k.frac_coords), "#", k.name)

    elif options.command == "bz":
        structure = abilab.Structure.from_file(options.filepath)
        structure.show_bz()

    elif options.command == "kmesh":
        structure = abilab.Structure.from_file(options.filepath)
        k = Ktables(structure, options.mesh, options.is_shift, not options.no_time_reversal)
        print(k)
        print("")
        print("NB: These results are obtained by calling spglib with the structure read from file.")
        print("The k-points might differ from the ones expected by Abinit, especially if the space groups differ.")

        if not options.verbose:
            print("\nUse -v to obtain the BZ --> IBZ mapping.")
        else:
            k.print_bz2ibz()

    elif options.command == "lgk":

        structure = abilab.Structure.from_file(options.filepath)
        spgrp = structure.abi_spacegroup

        if spgrp is None:
            cprint("Your file does not contain Abinit symmetry operations.", "yellow")
            cprint("Will call spglib to obtain the spacegroup (assuming time-reversal: %s)" %
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

    elif options.command == "pmgdata":
        # Get the Structure corresponding the a material_id.
        structure = abilab.Structure.from_material_id(options.pmgid, final=True,
                                                      api_key=options.mapi_key, endpoint=options.endpoint)

        if options.format == "abivars":
            print(structure.abi_string)
        else:
            # Convert to json and print it.
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
