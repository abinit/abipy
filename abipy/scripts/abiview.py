#!/usr/bin/env python
"""
This script visualizes results with external graphical applications.
or convert data from Abinit files (usually netcdf) to other formats.
"""
import sys
import os
import argparse

from monty.functools import prof_main
from monty.termcolor import cprint
from abipy import abilab
from abipy.iotools.visualizer import Visualizer
from abipy.tools.plotting import MplExpose, PanelExpose, GenericDataFilePlotter, plotlyfigs_to_browser, push_to_chart_studio


def handle_overwrite(path, options):
    """Exit 1 if file ``path`` exists and not options.force else return path."""
    name_parts = os.path.splitext(path)
    print("Writing %s file:" % name_parts[-1].replace(".", "").upper())
    if os.path.exists(path) and not options.force:
        cprint("Cannot overwrite pre-existent file. Use `-f` options.", "red")
        sys.exit(1)
    return path


def df_to_clipboard(options, df):
    """Copy dataframe to clipboard if options.clipboard."""
    if getattr(options, "clipboard", False):
        cprint("Copying dataframe to the system clipboard.", "green")
        cprint("This can be pasted into Excel, for example", "green")
        df.to_clipboard()


#def abiview_fields(options):
#    """Animate fields with Mayavi. Accept any file with density or potential ..."""
#    from abipy.display.mvtk import MayaviFieldAnimator
#    a = MayaviFieldAnimator(options.filepath)
#    a.volume_animate()
#    return 0


def abiview_structure(options):
    """
    Visualize the structure with the specified visualizer. Requires external app
    or optional python modules (mayavi, vtk)
    """
    structure = abilab.Structure.from_file(options.filepath)
    print(structure.to_string(verbose=options.verbose))
    print("Visualizing structure with:", options.appname)
    structure.visualize(appname=options.appname)
    return 0


def abiview_hist(options):
    """
    Visualize structural relaxation/molecular-dynamics run
    from data stored in the HIST.nc file. Requires mayavi.
    """
    with abilab.abiopen(options.filepath) as hist:
        print(hist.to_string(verbose=options.verbose))
        if options.appname is not None:
            hist.visualize(appname=options.appname, to_unit_cell=options.to_unit_cell)
        elif options.xdatcar:
            xpath = options.filepath + ".XDATCAR"
            handle_overwrite(xpath, options)
            hist.write_xdatcar(filepath=xpath, groupby_type=True, overwrite=True,
                               to_unit_cell=options.to_unit_cell)
        else:
            hist.plot()

        return 0


def abiview_data(options):
    """
    Extract data from a generic text file with results
    in tabular format and plot data with matplotlib.
    Multiple datasets are supported.
    No attempt is made to handle metadata (e.g. column name)
    """
    plotter = GenericDataFilePlotter(options.filepath)
    print(plotter.to_string(verbose=options.verbose))
    plotter.plot(use_index=options.use_index)

    return 0


#def abiview_abo(options):
#    """
#    Plot self-consistent iterations extracted from Abinit output file
#    as well as timer data (if present)
#    """
#    with abilab.abiopen(options.filepath) as abo:
#        print(abo.to_string(verbose=options.verbose))
#        abo.plot()
#        #timer = abo.get_timer()
#    return 0


def abiview_dirviz(options):
    """
    Visualize directory tree with graphviz.
    """
    from pymatgen.io.abinit.utils import Dirviz
    import tempfile
    graph = Dirviz(options.filepath).get_cluster_graph(engine=options.engine)
    directory = tempfile.mkdtemp()
    print("Producing source files in:", directory)
    graph.view(directory=directory)

    return 0


def abiview_ebands(options):
    """
    Plot electronic bands if file contains high-symmetry k-path or DOS if k-sampling.
    Accept any file with ElectronBands e.g. GSR.nc, WFK.nc, ...
    """
    with abilab.abiopen(options.filepath) as abifile:
        if options.xmgrace:
            outpath = options.filepath + ".agr"
            abifile.ebands.to_xmgrace(handle_overwrite(outpath, options))
        elif options.bxsf:
            outpath = options.filepath + ".bxsf"
            abifile.ebands.to_bxsf(handle_overwrite(outpath, options))
        #elif options.plotly:
        #    print(abifile.to_string(verbose=options.verbose))
        else:
            print(abifile.to_string(verbose=options.verbose))
            abifile.expose_ebands(slide_mode=options.slide_mode, slide_timeout=options.slide_timeout,
                                  expose_web=options.expose_web, verbose=options.verbose)

        return 0


def abiview_skw(options):
    """
    Interpolate energies in k-space along a k-path with star-function methods
    Note that the interpolation will likely fail if there are symmetrical k-points in the input set of k-points
    so it's recommended to call this method with energies obtained in the IBZ.
    Accept any file with ElectronBands e.g. GSR.nc, WFK.nc, ...
    """
    ebands = abilab.ElectronBands.as_ebands(options.filepath)
    if not ebands.kpoints.is_ibz:
        cprint("SKW interpolator should be called with energies in the IBZ", "yellow")
    r = ebands.interpolate(lpratio=options.lpratio, line_density=options.line_density, verbose=options.verbose)
    r.ebands_kpath.plot()
    return 0


def abiview_fs(options):
    """
    Extract eigenvalues in the IBZ from file and visualize Fermi surface with --appname
    """
    with abilab.abiopen(options.filepath) as abifile:
        eb3d = abifile.ebands.get_ebands3d()

        if options.appname == "mpl":
            eb3d.plot_isosurfaces()
        elif options.appname == "xsf":
            eb3d.xcrysden_view()
        elif options.appname == "mayavi":
            eb3d.mvplot_isosurfaces()
        else:
            raise ValueError("Unsupported value for --appname: %s" % str(options.appname))

        return 0


def abiview_ddb(options):
    """
    Invoke Anaddb to compute phonon bands and DOS from the DDB, plot the results.
    """
    with abilab.abiopen(options.filepath) as ddb:
        print(ddb.to_string(verbose=options.verbose))

        # Don't need PHDOS if phononwebsite
        nqsmall = 0 if options.phononwebsite else 10
        ndivsm = 20; asr = 2; chneut = 1; dipdip = 1; dos_method = "tetra"; lo_to_splitting = "automatic"
        print("""
Computing phonon bands and DOS from DDB file with:
nqsmall: {nqsmall}, ndivsm: {ndivsm},
asr: {asr}, chneut: {chneut}, dipdip: {dipdip}, lo_to_splitting: {lo_to_splitting}, dos_method: {dos_method}
""".format(**locals()))

        print("Invoking anaddb ...  ", end="")
        phbst_file, phdos_file = ddb.anaget_phbst_and_phdos_files(
            nqsmall=nqsmall, ndivsm=ndivsm, asr=asr, chneut=chneut, dipdip=dipdip, dos_method=dos_method,
            lo_to_splitting=lo_to_splitting, verbose=options.verbose, mpi_procs=1)
        print("Calculation completed.\nResults available in:", os.path.dirname(phbst_file.filepath))

        phbands = phbst_file.phbands
        units = "mev"

        if options.xmgrace:
            outpath = options.filepath + ".agr"
            phbands.to_xmgrace(handle_overwrite(outpath, options))
            return 0

        #elif options.bxsf:
        #    outpath = options.filepath + ".bxsf"
        #    phbands.to_bxsf(handle_overwrite(outpath, options))
        #    return 0

        elif options.phononwebsite:
            return phbands.view_phononwebsite(browser=options.browser, verbose=options.verbose)

        elif options.plotly:
            # Plotly + Panel version.
            phdos = phdos_file.phdos
            with PanelExpose(title=f"Vibrational properties of {phdos_file.structure.formula}") as e:
                e(phbands.qpoints.plotly(show=False))
                e(phbands.plotly_with_phdos(phdos, units=units, show=False))
                e(phdos_file.plotly_pjdos_type(units=units, show=False))
                e(phdos.plotly_harmonic_thermo(tstart=5, tstop=300, num=50, units="eV", formula_units=1,
                                               quantities="all", show=False))
                e(phdos.plotly(units=units, show=False))
                e(phbands.plot_colored_matched(units=units, show=False))
                e(phbands.plot_fatbands(units=units, phdos_file=phdos_file, show=False))
                #push_to_chart_studio(figs) if options.chart_studio else plotlyfigs_to_browser(figs)

        else:
            phdos = phdos_file.phdos

            if not options.expose_web:
                # matplotlib figure and X-server.
                e = MplExpose(slide_mode=options.slide_mode, slide_timeout=options.slide_timeout)
            else:
                # panel version with matplotlib.
                e = PanelExpose(title=f"Vibrational properties of {phdos_file.structure.formula}")

            with e:
                e(phbands.qpoints.plot(show=False))
                e(phbands.plot_with_phdos(phdos, units=units, show=False))
                e(phdos_file.plot_pjdos_type(units=units, show=False))
                e(phbands.plot_colored_matched(units=units, show=False))
                e(phdos.plot_harmonic_thermo(tstart=5, tstop=300, num=50, units="eV", formula_units=1,
                             quantities="all", show=False))
                e(phdos.plot(units=units, show=False))
                e(phbands.plot_fatbands(units=units, phdos_file=phdos_file, show=False))

        phbst_file.close()
        phdos_file.close()

    return 0


def abiview_ddb_vs(options):
    """
    Compute speed of sound by fitting phonon frequencies along selected directions.
    """
    num_points = 20; asr = 2; chneut = 1; dipdip = 1
    print("""
Computing phonon frequencies for linear least-squares with:
num_points: {num_points}, asr: {asr}, chneut: {chneut}, dipdip: {dipdip}
""".format(**locals()))

    print("Invoking anaddb ...  ")
    sv = abilab.SoundVelocity.from_ddb(options.filepath, num_points=num_points,
                                       asr=asr, chneut=chneut, dipdip=dipdip, verbose=options.verbose)
    #print("Calculation completed.\nResults available in:", os.path.dirname(phbst_file.filepath))

    df = sv.get_dataframe()
    abilab.print_dataframe(df, title="Speed of sound for different directions:")
    df_to_clipboard(options, df)
    sv.plotly if options.plotly else sv.plot()

    return 0


def abiview_ddb_ir(options):
    """
    Compute infra-red spectrum from DDB. Plot results.
    """
    asr = 2; chneut = 1; dipdip = 1
    print("""
Computing phonon frequencies for infra-red spectrum with:
asr: {asr}, chneut: {chneut}, dipdip: {dipdip}
""".format(**locals()))

    with abilab.abiopen(options.filepath) as ddb:
        tgen = ddb.anaget_dielectric_tensor_generator(asr=asr, chneut=chneut, dipdip=dipdip, verbose=options.verbose)
        print(tgen.to_string(verbose=options.verbose))

        gamma_ev = 1e-4
        print("Plotting dielectric tensor with constant phonon damping: %s (eV)" % gamma_ev)
        if options.plotly:
            tgen.plotly_all(gamma_ev=gamma_ev)
        else:
            tgen.plot_all(gamma_ev=gamma_ev)

    return 0


def abiview_ddb_asr(options):
    """
    Compute phonon band structure from DDB with/without acoustic sum rule. Plot results.
    """
    print("Computing phonon frequencies with/without ASR")

    with abilab.abiopen(options.filepath) as ddb:
        plotter = ddb.anacompare_asr(asr_list=(0, 2), chneut_list=(1,), dipdip=1, lo_to_splitting="automatic",
                                     nqsmall=10, ndivsm=20, dos_method="tetra", ngqpt=None,
                                     verbose=0, mpi_procs=1)

        title = ddb.structure.formula
        renderer = "browser" if not options.chart_studio else "chart_studio"
        plotter.combiplotly(renderer=renderer, title=title) if options.plotly else plotter.plot(title=title)

    return 0


def abiview_ddb_dipdip(options):
    """
    Compute phonon band structure from DDB with/without dipole-dipole interaction. Plot results.
    """
    print("Computing phonon frequencies with/without dipdip")

    with abilab.abiopen(options.filepath) as ddb:
        plotter = ddb.anacompare_dipdip(chneut_list=(1,), asr=2, lo_to_splitting="automatic",
                                        nqsmall=10, ndivsm=20, dos_method="tetra", ngqpt=None,
                                        verbose=0, mpi_procs=1)

        title = ddb.structure.formula
        renderer = "browser" if not options.chart_studio else "chart_studio"
        plotter.combiplotly(renderer=renderer, title=title) if options.plotly else plotter.plot(title=title)

    return 0


def abiview_ddb_quad(options):
    """
    Compute phonon band structure from DDB with/without quadrupole terms. Plot results.
    """
    print("Computing phonon frequencies with/without dip-quad and quad-quad terms.")

    with abilab.abiopen(options.filepath) as ddb:
        plotter = ddb.anacompare_quad(asr=2, chneut=1, dipdip=1, lo_to_splitting="automatic",
                                      nqsmall=0, ndivsm=20, dos_method="tetra", ngqpt=None,
                                      verbose=0, mpi_procs=1)

        title = ddb.structure.formula
        renderer = "browser" if not options.chart_studio else "chart_studio"
        plotter.combiplotly(renderer=renderer, title=title) if options.plotly else plotter.plot(title=title)

    return 0


def abiview_ddb_isodistort_ph(options):
    """
    Compute ph-freqs for given q-point (default: Gamma), produce CIF files for unperturbed and distorded structure
    that can be used with ISODISTORT (https://stokes.byu.edu/iso/isodistort.php) to analyze the symmetry of phonon modes.
    See README.me file produced in output directory.
    """
    with abilab.abiopen(options.filepath) as ddb:
        print(ddb.to_string(verbose=options.verbose))
        qpoint = options.qpoint
        asr = 2; chneut = 1; dipdip = 1; lo_to_splitting = False
        print("""
Computing phonon frequencies and eigenvectors with:
asr: {asr}, chneut: {chneut}, dipdip: {dipdip}, lo_to_splitting: {lo_to_splitting}
qpoint = {qpoint}
""".format(**locals()))

        print("Invoking anaddb ...  ", end="")
        phbands = ddb.anaget_phmodes_at_qpoint(qpoint=qpoint, asr=asr, chneut=chneut, dipdip=dipdip,
                                               verbose=options.verbose, lo_to_splitting=False)

        phbands.make_isodistort_ph_dir(qpoint, select_modes=None, eta=1, workdir=None)

    return 0


def abiview_ddb_ifc(options):
    """
    Visualize interatomic force constants in real space.
    """
    asr = 2; chneut = 1; dipdip = 1
    print("""
Computing interatomic force constants with
asr: {asr}, chneut: {chneut}, dipdip: {dipdip}
""".format(**locals()))

    with abilab.abiopen(options.filepath) as ddb:
        # Execute anaddb to compute the interatomic force constants.
        ifc = ddb.anaget_ifc(asr=asr, chneut=chneut, dipdip=dipdip)
        #print(ifc)

        if not options.expose_web:
            # matplotlib figure and X-server.
            e = MplExpose(slide_mode=options.slide_mode, slide_timeout=options.slide_timeout)
        else:
            # panel version
            e = PanelExpose(title=f"Interatomic force constants of {phdos_file.structure.formula}")

        with e:
            e(ifc.plot_longitudinal_ifc(title="Longitudinal IFCs", show=False))
            e(ifc.plot_longitudinal_ifc_short_range(title="Longitudinal IFCs short range", show=False))
            e(ifc.plot_longitudinal_ifc_ewald(title="Longitudinal IFCs Ewald", show=False))

    return 0


def abiview_phbands(options):
    """
    Plot phonon bands. Accept any file with PhononBands e.g. PHBST.nc, ...
    """
    with abilab.abiopen(options.filepath) as abifile:
        if options.xmgrace:
            outpath = options.filepath + ".agr"
            abifile.phbands.to_xmgrace(handle_overwrite(outpath, options))
        #elif options.bxsf:
        #    outpath = options.filepath + ".bxsf"
        #    abifile.phbands.to_bxsf(handle_overwrite(outpath, options))
        #    return 0
        elif options.phononwebsite:
            return abifile.phbands.view_phononwebsite(browser=options.browser)
        else:
            print(abifile.to_string(verbose=options.verbose))
            abifile.expose_phbands(slide_mode=options.slide_mode, slide_timeout=options.slide_timeout,
                           verbose=options.verbose, units="mev")

        return 0


def abiview_denpot(options):
    """
    Plot DEN/POT files
    """
    with abilab.abiopen(options.filepath) as abifile:
        print(abifile.to_string(verbose=options.verbose))
        field = abifile.field
        if options.chgcar:
            if not field.is_density_like:
                cprint("Warning: expecting DENSITY file but received %s" % field.__class__.__name__, "yellow")
            chg_path = options.filepath + ".CHGCAR"
            field.to_chgcar(filename=handle_overwrite(chg_path, options))
        #elif options.cube:
        #    outpath = options.filepath + ".cube"
        #    field.export_to_cube(filename=handle_overwrite(outpath, options), spin="total")
        else:
            abifile.field.visualize(appname=options.appname)
    return 0


def abiview_lobster(options):
    """Analyze lobster output files in directory."""
    from abipy.electrons.lobster import LobsterAnalyzer
    lobana = LobsterAnalyzer.from_dir(os.path.dirname(options.filepath), prefix=options.prefix)
    print(lobana.to_string(verbose=options.verbose))

    if options.ipython:
        # Start ipython shell with namespace
        # Use embed because I don't know how to show a header with start_ipython.
        import IPython
        IPython.embed(header="The LobsterAnalyzer is bound to the `lobana` variable.\nTry `print(lobana)`")

    elif options.notebook:
        return lobana.make_and_open_notebook(foreground=options.foreground)

    else:
        lobana.plot()
        #lobana.plot_coxp_with_dos(from_site_index=[0, 1])

    return 0


def get_epilog():
    return """\
Usage example:

###########
# Structure
###########

    abiview.py structure FILE                ==> Visualize structure with Vesta (default)
    abiview.py structure FILE -a xcrysden    ==> Visualize structure with Xcrysden (default)
    abiview.py hist out_HIST.nc              ==> Plot relaxation/molecular dynamics results with matplotlib.
    abiview.py hist out_HIST.nc -a ovito     ==> Visualize relaxation/molecular dynamics results with ovito.
    abiview.py hist out_HIST.nc --xdatcar    ==> Convert HIST.nc into XDATCAR format (caveat: assume fixed unit cell!)

############
# Text files
############

    abiview.py data FILE     ==> Parse text FILE with data in tabular format and plot arrays.

###########
# Electrons
###########

    abiview.py ebands out_WFK.nc              ==>  Plot electrons bands (or DOS) with matplotlib.
    abiview.py ebands out_GSR.nc --xmgrace    ==>  Generate xmgrace file with electron bands.
    abiview.py fs FILE_WITH_KMESH.nc          ==>  Visualize Fermi surface from netcdf file with electron energies
                                                   on a k-mesh. Use -a xsf to change application e.g. Xcrysden.
    abiview.py skw out_GSR.nc                 ==> Interpolate IBZ energies with star-functions and plot
                                                  interpolated bands.

#########
# Phonons
#########

    abiview.py ddb in_DDB                 ==>  Compute ph-bands and DOS from DDB, plot results.
    abiview.py ddb_vs                     ==>  Compute speed of sound from DDB by fitting phonon frequencies.
    abiview.py ddb_ir                     ==>  Compute infra-red spectrum from DDB. Plot results.
    abiview.py ddb_asr                    ==>  Compute ph-bands from DDB with/wo acoustic rule.
                                               Plot results with matplotlib (default) or plotly (--plotly)
    abiview.py ddb_asr --plotly -cs       ==>  Compute ph-bands from DDB with/wo acoustic rule.
                                               Plot results with plotly and push figure to plotly chart studio cloud.
                                               See: https://chart-studio.plotly.com/
    abiview.py ddb_dipdip                 ==>  Compute ph-bands from DDB with/wo dipole-dipole treatment.
                                               Plot results with matplotlib (default) or plotly (--plotly)
    abiview.py ddb_quad                   ==>  Compute ph-bands from DDB with/wo dipole-quadrupole terms.
                                               Plot results with matplotlib (default) or plotly (--plotly)
    abiview.py ddb_ifc                    ==>  Visualize interatomic force constants in real space.
    abiview.py phbands out_PHBST.nc -web  ==>  Visualize ph-bands and displacements with phononwebsite.

###############
# Miscelleanous
###############

  abiview.py dirviz DIRECTORY            ==> Visualize directory tree with graphviz.
  abiview.py lobster DIRECTORY           ==> Visualize Lobster results.

Use `abiview.py --help` for help and `abiview.py COMMAND --help` to get the documentation for `COMMAND`.
Use `-v` to increase verbosity level (can be supplied multiple times e.g -vv).
"""

# TODO
#abiview.py denpot out_DEN.nc              ==>  Visualize density with Vesta.
#abiview.py denpot out_DEN.nc --chgcar     ==>  Convert DEN file into CHGCAR fileformat.


def get_parser(with_epilog=False):

    # Parent parser for common options.
    copts_parser = argparse.ArgumentParser(add_help=False)
    copts_parser.add_argument('filepath', type=str, help="File to visualize.")

    copts_parser.add_argument('--loglevel', default="ERROR", type=str,
        help="Set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")
    copts_parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
        help='verbose, can be supplied multiple times to increase verbosity.')
    copts_parser.add_argument('-sns', "--seaborn", const="paper", default=None, action='store', nargs='?', type=str,
        help='Use seaborn settings. Accept value defining context in ("paper", "notebook", "talk", "poster"). Default: paper')
    copts_parser.add_argument('-mpl', "--mpl-backend", default=None,
        help=("Set matplotlib interactive backend. "
              "Possible values: GTKAgg, GTK3Agg, GTK, GTKCairo, GTK3Cairo, WXAgg, WX, TkAgg, Qt4Agg, Qt5Agg, macosx."
              "See also: https://matplotlib.org/faq/usage_faq.html#what-is-a-backend."))

    # Parent parser for commands supporting MplExpose.
    slide_parser = argparse.ArgumentParser(add_help=False)
    slide_parser.add_argument("-s", "--slide-mode", default=False, action="store_true",
            help="Iterate over figures. Expose all figures at once if not given on the CLI.")
    slide_parser.add_argument("-t", "--slide-timeout", type=int, default=None,
            help="Close figure after slide-timeout seconds (only if slide-mode). Block if not specified.")
    slide_parser.add_argument("-ew", "--expose-web", default=False, action="store_true",
            help='Generate matplotlib plots in $BROWSER instead of X-server. WARNING: Not all the features are supported.')

    # Parent parser for commands that operating on pandas dataframes
    pandas_parser = argparse.ArgumentParser(add_help=False)
    pandas_parser.add_argument("-c", '--clipboard', default=False, action="store_true",
            help="Copy dataframe to the system clipboard. This can be pasted into Excel, for example")

    # Parent parser for commands supporting ipython
    ipy_parser = argparse.ArgumentParser(add_help=False)
    ipy_parser.add_argument('-ipy', '--ipython', default=False, action="store_true", help='Invoke ipython terminal.')

    # Parent parser for commands supporting (jupyter notebooks)
    nb_parser = argparse.ArgumentParser(add_help=False)
    nb_parser.add_argument('-nb', '--notebook', default=False, action="store_true", help='Generate jupyter notebook.')
    nb_parser.add_argument('--foreground', action='store_true', default=False,
        help="Run jupyter notebook in the foreground.")

    # Parent parser for commands supporting plotly plots
    plotly_parser = argparse.ArgumentParser(add_help=False)
    plotly_parser.add_argument("-ply", '--plotly', default=False, action="store_true",
        help='Generate plotly plots in browser instead of matplotlib.')
    plotly_parser.add_argument("-cs", "--chart-studio", default=False, action="store_true",
        help="Push figure to plotly chart studio. " +
             "Requires --plotly and user account at https://chart-studio.plotly.com.")

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=get_epilog() if with_epilog else "",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-V', '--version', action='version', version=abilab.__version__)

    # Create the parsers for the sub-command
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
            elif arg == "force":
                p.add_argument("-f", '--force', default=False, action="store_true",
                    help="Overwrite pre-existent files without prompting for confirmation.")
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
              "Possible options: `%s`, `mpl` (matplotlib) `mayavi`, `vtk`" % ", ".join(Visualizer.all_visunames())))
    p_hist.add_argument("--xdatcar", default=False, action="store_true", help="Convert HIST file into XDATCAR format.")
    p_hist.add_argument("--to-unit-cell", default=False, action="store_true",
            help="Whether to translate sites into the unit cell.")
    add_args(p_hist, "force")

    # Subparser for data command.
    p_data = subparsers.add_parser('data', parents=[copts_parser], help=abiview_data.__doc__)
    p_data.add_argument("-i", "--use-index", default=False, action="store_true",
        help="Use the row index as x-value in the plot. By default the plotter uses the first column as x-values")

    # Subparser for abo command.
    #p_abo = subparsers.add_parser('abo', parents=[copts_parser], help=abiview_abo.__doc__)

    # Subparser for dirviz command.
    p_dirviz = subparsers.add_parser('dirviz', parents=[copts_parser], help=abiview_dirviz.__doc__)
    p_dirviz.add_argument("-e", "--engine", type=str, default="fdp",
        help=("graphviz engine: ['dot', 'neato', 'twopi', 'circo', 'fdp', 'sfdp', 'patchwork', 'osage']. "
            "See http://www.graphviz.org/pdf/dot.1.pdf "
            "Use `conda install python-graphviz` or `pip install graphviz` to install the python package."))

    # Subparser for ebands command.
    p_ebands = subparsers.add_parser('ebands', parents=[copts_parser, slide_parser], help=abiview_ebands.__doc__)
    add_args(p_ebands, "xmgrace", "bxsf", "force")

    # Subparser for ebands command.
    p_skw = subparsers.add_parser('skw', parents=[copts_parser], help=abiview_skw.__doc__)
    p_skw.add_argument("-lp", "--lpratio", type=int, default=5,
        help=("Ratio between the number of star functions and the number of ab-initio k-points. "
              "The default should be OK in many systems, larger values may be required for accurate derivatives."))
    p_skw.add_argument("-ld", "--line-density", type=int, default=20,
                      help="Number of points in the smallest segment of the k-path.")

    # Subparser for fs command.
    p_fs = subparsers.add_parser('fs', parents=[copts_parser], help=abiview_fs.__doc__)
    p_fs.add_argument("-a", "--appname", type=str, default="mpl",
        help="Application name. Possible options: mpl (matplotlib, default), xsf (xcrysden), mayavi.")

    # Subparser for ddb command.
    p_ddb = subparsers.add_parser('ddb', parents=[copts_parser, slide_parser, plotly_parser], help=abiview_ddb.__doc__)
    add_args(p_ddb, "xmgrace", "phononweb", "browser", "force")

    # Subparser for ddb_vs command.
    p_ddb_vs = subparsers.add_parser('ddb_vs', parents=[copts_parser, pandas_parser, slide_parser],
                                     help=abiview_ddb_vs.__doc__)

    # Subparser for ddb_ir command.
    p_ddb_ir = subparsers.add_parser('ddb_ir', parents=[copts_parser, pandas_parser, slide_parser],
                                     help=abiview_ddb_ir.__doc__)

    # Subparser for ddb_asr command.
    p_ddb_asr = subparsers.add_parser('ddb_asr', parents=[copts_parser, pandas_parser, slide_parser, plotly_parser],
                                      help=abiview_ddb_asr.__doc__)

    # Subparser for ddb_dipdip command.
    p_ddb_dipdip = subparsers.add_parser('ddb_dipdip', parents=[copts_parser, pandas_parser, slide_parser, plotly_parser],
                                         help=abiview_ddb_dipdip.__doc__)

    # Subparser for ddb_quad command.
    p_ddb_quad = subparsers.add_parser('ddb_quad', parents=[copts_parser, pandas_parser, slide_parser, plotly_parser],
                                       help=abiview_ddb_quad.__doc__)

    # Subparser for ddb_ph_isodistort command.
    p_ddb_isodistort_ph = subparsers.add_parser('ddb_isodistort_ph', parents=[copts_parser],
                                     help=abiview_ddb_isodistort_ph.__doc__)
    p_ddb_isodistort_ph.add_argument("-q", "--qpoint", nargs=3, type=float,
        help="q-point in reduced coordinates e.g. 0.25 0 0. Default: 0, 0, 0", default=[0, 0, 0])

    # Subparser for ddb_ifc command.
    p_ddb_ifc = subparsers.add_parser('ddb_ifc', parents=[copts_parser, pandas_parser, slide_parser],
                                     help=abiview_ddb_ifc.__doc__)

    # Subparser for phbands command.
    p_phbands = subparsers.add_parser('phbands', parents=[copts_parser, slide_parser], help=abiview_phbands.__doc__)
    add_args(p_phbands, "xmgrace", "phononweb", "browser", "force")

    # Subparser for lobster command.
    p_lobster = subparsers.add_parser('lobster', parents=[copts_parser, ipy_parser, nb_parser],
        help=abiview_lobster.__doc__)
    p_lobster.add_argument("--prefix", type=str, default="", help="Prefix for lobster output files. Default: ''")

    # Subparser for denpot command.
    #p_denpot = subparsers.add_parser('denpot', parents=[copts_parser], help=abiview_denpot.__doc__)
    #p_denpot.add_argument("-a", "--appname", nargs="?", default=None, const="vesta",
    #        help=("Application name. Default: vesta. "
    #              "Possible options: `%s`, `mayavi`, `vtk`" % ", ".join(Visualizer.all_visunames())))
    #p_denpot.add_argument("--chgcar", default=False, action="store_true", "Convert Density to CHGCAR format.")
    #p_denpot.add_argument("--cube", default=False, action="store_true", "Convert Density/Potential to CUBE format.")

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

    if not options.command:
        show_examples_and_exit(error_code=1)

    # loglevel is bound to the string value obtained from the command line argument.
    # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
    import logging
    numeric_level = getattr(logging, options.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.loglevel)
    logging.basicConfig(level=numeric_level)

    if getattr(options, "plotly", False): options.expose = True
    if options.verbose > 2: print(options)

    if options.mpl_backend is not None:
        # Set matplotlib backend
        import matplotlib
        matplotlib.use(options.mpl_backend)

    if options.seaborn:
        # Use seaborn settings.
        import seaborn as sns
        sns.set(context=options.seaborn, style='darkgrid', palette='deep',
                font='sans-serif', font_scale=1, color_codes=False, rc=None)

    # Dispatch
    return globals()["abiview_" + options.command](options)


if __name__ == "__main__":
    sys.exit(main())
