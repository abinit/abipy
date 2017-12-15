#!/usr/bin/env python
"""
This script visualizes data with external graphical applications.
It's also possible to generate movies from Abinit output files.

WARNING: This script is still under active development. API will change!
"""
from __future__ import unicode_literals, division, print_function, absolute_import

import sys
import os
import argparse
import numpy as np

from monty.functools import prof_main
from abipy import abilab
from abipy.iotools.visualizer import Visualizer


def abiview_structure(options):
    """
    Visualize the structure with the specified visualizer. Requires external app or optional python modules."
    """
    structure = abilab.Structure.from_file(options.filepath)
    print(structure)
    print("Visualizing structure with:", options.appname)
    structure.visualize(appname=options.appname)


def abiview_hist(options):
    """
    Visualize structural relaxation/molecular-dynamics run
    from data stored in the HIST.nc file. Requires mayavi.
    """
    with abilab.abiopen(options.filepath) as hist:
        print(hist.to_string(verbose=options.verbose))
        if options.appname is not None:
            hist.visualize(appname=options.appname)
        else:
            hist.plot()

        return 0


def abiview_abo(options):
    """
    Plot self-consisten iterations extracted from Abinit output file
    as well as timer data (if present)
    """
    with abilab.abiopen(options.filepath) as abo:
        print(abo.to_string(verbose=options.verbose))
        abo.plot()
        #timer = abo.get_timer()
    return 0


def abiview_log(options):
    """
    Print WARNINGS/COMMENTS/ERRORS found in Abinit log file.
    """
    with abilab.abiopen(options.filepath) as abo:
        print(abo.to_string(verbose=options.verbose))
    return 0


def abiview_ebands(options):
    """
    Plot electronic bands if file contains high-symmetry k-path or DOS if k-sampling.
    Accept any file with ElectronBands e.g. GSR.nc, WFK.nc, ...
    """
    with abilab.abiopen(options.filepath) as abifile:
        if options.xmgrace:
            abifile.ebands.to_xmgrace(sys.stdout)
        elif options.bxsf:
            abifile.ebands.to_bxsf(sys.stdout)
        else:
            print(abifile.to_string(verbose=options.verbose))
            if abifile.ebands.kpoints.is_path:
                abifile.ebands.plot()
            else:
                abifile.ebands.get_edos().plot()
        return 0


def abiview_fatbands(options):
    """
    Plot electronic fatbands bands if file contains high-symmetry k-path
    or PJDOS if k-sampling. Requires FATBANDS.nc file.
    """
    with abilab.abiopen(options.filepath) as abifile:
        print(abifile.to_string(verbose=options.verbose))
        if abifile.ebands.kpoints.is_path:
            abifile.plot_fatbands_lview()
            #abifile.plot_fatbands_typeview()
        else:
            abifile.plot_pjdos_lview()
            #abifile.plot_pjdos_typeview()
    return 0


def abiview_optic(options):
    """
    Plot optical spectra. Requires OPTIC.nc file
    """
    with abilab.abiopen(options.filepath) as abifile:
        print(abifile.to_string(verbose=options.verbose))
        #abifile.plot_linear_epsilon()
        #abifile.plot_linopt(self, select="all", itemp=0, xlims=None, **kwargs):
        #abifile.plot_shg()
        #abifile.plot_leo()
    return 0


def abiview_ddb(options):
    """
    Invoke Anaddb to compute phonon bands and DOS from DDB, plot results.
    """
    with abilab.abiopen(options.filepath) as ddb:
        print(ddb.to_string(verbose=options.verbose))

        # TODO: Autodetect presence of lo_to_splitting data in DDB.
        phbst, phdos = ddb.anaget_phbst_and_phdos_files(
            nqsmall=10, ndivsm=20, asr=2, chneut=1, dipdip=1, dos_method="tetra",
            lo_to_splitting=False, ngqpt=None, qptbounds=None, anaddb_kwargs=None,
            verbose=options.verbose, mpi_procs=1)

        phbands = phbst.phbands

        if options.xmgrace:
            phbands.to_xmgrace(sys.stdout)
            return 0
        #elif options.bxsf:
        #    phbands.to_bxsf(sys.stdout)
        #    return 0
        elif options.phononwebsite:
            return phbands.view_phononwebsite(browser=options.browser)
        else:
            phbands.plot_with_phdos(phdos)

        phbst.close()
        phdos.close()

    return 0


def abiview_phbands(options):
    """Plot phonon bands. Accept any file with PhononBands e.g. PHBST.nc, ..."""
    with abilab.abiopen(options.filepath) as abifile:
        if options.xmgrace:
            abifile.phbands.to_xmgrace(sys.stdout)
        #elif options.bxsf:
        #    abifile.phbands.to_bxsf(sys.stdout)
        #    return 0
        elif options.phononwebsite:
            return abifile.phbands.view_phononwebsite(browser=options.browser)
        else:
            print(abifile.to_string(verbose=options.verbose))
            abifile.phbands.plot()

        return 0

def abiview_phdos(options):
    """Plot phonon DOS. Require PHDOS.nc file."""
    with abilab.abiopen(options.filepath) as abifile:
        print(abifile.to_string(verbose=options.verbose))
        abifile.phdos.plot()
    return 0


def abiview_mdf(options):
    """
    Plot the macroscopic dielectric functions computed by the Bethe-Salpeter code.
    Requires MDF.nc file.
    """
    with abilab.abiopen(options.filepath) as abifile:
        print(abifile.to_string(verbose=options.verbose))
        abifile.plot_mdfs()
    return 0


def abiview_denpot(options):
    """
    Plot DEN/POT files
    """
    with abilab.abiopen(options.filepath) as abifile:
        print(abifile.to_string(verbose=options.verbose))
        abifile.field.visualize(appname="vesta")
    return 0


def abiview_gruns(options):
    """Plot Grunesein parameters. Requires GRUNS.nc file."""
    with abilab.abiopen(options.filepath) as abifile:
        print(abifile.to_string(verbose=options.verbose))
        abifile.plot_phbands_with_gruns()
    return 0


def abiview_eph(options):
    """Plot Eliashberg function. Requires EPH.nc file."""
    with abilab.abiopen(options.filepath) as abifile:
        print(abifile.to_string(verbose=options.verbose))
        abifile.plot_with_a2f()
    return 0


def abiview_sigeph(options):
    """Plot e-ph self-energy. Requires SIGEPH.nc file."""
    with abilab.abiopen(options.filepath) as abifile:
        print(abifile.to_string(verbose=options.verbose))
        abifile.plot()
    return 0


def get_epilog():
    return """\
Usage example:

###########
# Structure
###########

    abiview.py structure FILE                ==> Visualize structure with
    abiview.py hist out_HIST.nc              ==> Plot relaxation/molecular dynamics results with matplotlib.
    abiview.py hist out_HIST.nc  -a ovito    ==> Visualize relaxation/molecular dynamics results with ovito.

############
# Text files
############

    abiview.py abo run.abo          ==> Plot SCF iterations extracted from Abinit output file.
    abiview.py log run.log          ==> Print warnings/comments/errors found in Abinit log file.

###########
# Electrons
###########

    abiview.py ebands out_GSR.nc --xmgrace    ==>   Generate xmgrace file.
    abiview.py fatbands out_FATBANDS.nc       ==>   Plot electron fatbands or PJDOS depending on k-sampling.
    abiview.py denpot out_DEN.nc              ==>   Animate densities on FFT mesh.

#########
# Phonons
#########

    abiview.py ddb out_DDB                ==>  Compute ph-bands and DOS from DDB, plot results.
    abiview.py phbands out_PHBST.nc -web  ==>  Visualize phonon bands and displacements with phononwebsite.
    abiview.py phdos out_PHDOS.nc         ==>  Plot phonon DOS.

#########
# E-PH
#########

  abiview.py eph out_EPH.nc                   => Plot EPH results.
  abiview.py sigeph out_SIGEPH.nc             => Plot Fan-Migdal self-energy.

########
# GW/BSE
########

  abicomp.py sigres out_SIGRES.nc                 => Compare multiple SIGRES files.
  abicomp.py mdf out_MDF.nc --seaborn             => Compare macroscopic dielectric functions.
                                                     Use seaborn settings.


Use `abiview.py --help` for help and `abiview.py COMMAND --help` to get the documentation for `COMMAND`.
Use `-v` to increase verbosity level (can be supplied multiple times e.g -vv).
"""

def get_parser(with_epilog=False):

    # Parent parser for common options.
    copts_parser = argparse.ArgumentParser(add_help=False)
    copts_parser.add_argument('filepath', type=str, help="File to visualize.")

    copts_parser.add_argument('--loglevel', default="ERROR", type=str,
                         help="Set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")
    copts_parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                         help='verbose, can be supplied multiple times to increase verbosity.')
    copts_parser.add_argument('--seaborn', action="store_true", help="Use seaborn settings.")

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=get_epilog() if with_epilog else "",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-V', '--version', action='version', version=abilab.__version__)

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    def add_args(p, *args):
        """Add arguments to subparser `p`."""
        for arg in args:
            if arg == "xmgrace":
                p.add_argument('--xmgrace', default=False, action="store_true",
                    help="Print bands in xmgrace format to stdout and exit.")
            elif arg == "bxsf":
                p.add_argument('--bxsf', default=False, action="store_true",
                    help=("Generate BXSF file suitable for the visualization of isosurfaces with Xcrysden"
                          "(xcrysden --bxsf FILE).\n Requires k-points in IBZ. Print to stdout and exit."))
            elif arg == "phononweb":
                p.add_argument("-web", "--phononwebsite", default=False, action="store_true",
                    help=("Visualize phonon band structure on the phononwebsite. "
                          "http://henriquemiranda.github.io/phononwebsite/"))
            elif arg == "browser":
                p.add_argument("-b", "--browser", default=None,
                    help="Define browser used by python webbrowser. "
                         "See https://docs.python.org/2/library/webbrowser.html#webbrowser.register")
            else:
                raise ValueError("Invalid arg: %s" % arg)

    # Subparser for structure command.
    p_structure = subparsers.add_parser('structure', parents=[copts_parser], help=abiview_structure.__doc__)
    p_structure.add_argument("-a", "--appname", nargs="?", type=str, default="vesta",
        help=("Application name. Default: vesta. "
              "Possible options: %s, mayavi, vtk" % ", ".join(Visualizer.all_visunames())))

    # Subparser for hist command.
    p_hist = subparsers.add_parser('hist', parents=[copts_parser], help=abiview_hist.__doc__)
    p_hist.add_argument("-a", "--appname", nargs="?", default=None, const="ovito",
            help=("Application name. Default: ovito. "
                  "Possible options: `%s`, `mayavi`, `vtk`" % ", ".join(Visualizer.all_visunames())))

    # Subparser for abo command.
    p_abo = subparsers.add_parser('abo', parents=[copts_parser], help=abiview_abo.__doc__)

    # Subparser for log command.
    p_log = subparsers.add_parser('log', parents=[copts_parser], help=abiview_log.__doc__)

    # Subparser for ebands commands.
    p_ebands = subparsers.add_parser('ebands', parents=[copts_parser], help=abiview_ebands.__doc__)
    add_args(p_ebands, "xmgrace", "bxsf")

    # Subparser for fatbands commands.
    p_fatbands = subparsers.add_parser('fatbands', parents=[copts_parser], help=abiview_fatbands.__doc__)

    # Subparser for ddb command.
    p_ddb = subparsers.add_parser('ddb', parents=[copts_parser], help=abiview_ddb.__doc__)
    add_args(p_ddb, "xmgrace", "phononweb", "browser")

    # Subparser for phbands command.
    p_phbands = subparsers.add_parser('phbands', parents=[copts_parser], help=abiview_phbands.__doc__)
    add_args(p_phbands, "xmgrace", "phononweb", "browser")

    # Subparser for phdos command.
    p_phdos = subparsers.add_parser('phdos', parents=[copts_parser], help=abiview_phdos.__doc__)

    # Subparser for gruns command.
    p_gruns = subparsers.add_parser('gruns', parents=[copts_parser], help=abiview_gruns.__doc__)

    # Subparser for mdf command.
    p_mdf = subparsers.add_parser('mdf', parents=[copts_parser], help=abiview_mdf.__doc__)

    # Subparser for denpot command.
    p_denpot = subparsers.add_parser('denpot', parents=[copts_parser], help=abiview_denpot.__doc__)

    return parser


@prof_main
def main():

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(get_epilog())
        if err_msg:
            sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    parser = get_parser(with_epilog=True)

    # Parse the command line.
    try:
        options = parser.parse_args()
    except Exception:
        show_examples_and_exit(error_code=1)

    # loglevel is bound to the string value obtained from the command line argument.
    # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
    import logging
    numeric_level = getattr(logging, options.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.loglevel)
    logging.basicConfig(level=numeric_level)

    if options.seaborn:
        # Use seaborn settings.
        import seaborn as sns
        sns.set(context='article', style='darkgrid', palette='deep',
                font='sans-serif', font_scale=1, color_codes=False, rc=None)

    if options.verbose > 2:
        print(options)

    # Dispatch
    return globals()["abiview_" + options.command](options)

if __name__ == "__main__":
    sys.exit(main())
