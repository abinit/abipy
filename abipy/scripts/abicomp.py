#!/usr/bin/env python
"""
Script to analyze/compare results stored in multiple netcdf/output files.
By default the script displays the results/plots in the shell.
Use --ipython to start an ipython terminal or -nb to generate an ipython notebook.
"""
from __future__ import unicode_literals, division, print_function, absolute_import

import sys
import os
import argparse
import numpy as np

from collections import OrderedDict
from pprint import pprint
from monty.functools import prof_main
from monty.termcolor import cprint
from abipy import abilab
from abipy.tools.plotting import get_ax_fig_plt, GenericDataFilesPlotter


# Not used but could be useful to analyze densities.
def sort_paths(options):
    """
    Sort input files whose name is in the form `out_TIM2_DEN`
    Files are sorted by TIM index.
    """
    if options.no_sort: return
    names = [os.path.basename(p) for p in options.filepath]
    import re
    # out_TIM2_DEN
    tim = re.compile(r".+_TIM(\d+)_.+")
    l = []
    for p, n in zip(options.filepath, names):
        m = tim.match(n)
        if m:
            l.append((int(m.group(1)), p))
    if not l: return
    if len(l) != len(options.filepath):
        print("Cannot sort input path!")
        return

    options.paths = [t[1] for t in sorted(l, key=lambda t: t[0])]
    print("Input files have been automatically sorted")
    for i, p in enumerate(options.paths):
        print("%d: %s" % (i, p))
    print("Use --no-sort to disable automatic sorting.")


def remove_disordered(structures, paths):
    """Remove disordered structures and print warning message."""
    slist = []
    for s, p in zip(structures, paths):
        if not s.is_ordered:
            cprint("Removing disordered structure: %s found in %s" % (s.formula, p), "magenta")
        else:
            slist.append(s)
    return slist


def df_to_clipboard(options, df):
    """Copy dataframe to clipboard if options.clipboard."""
    if getattr(options, "clipboard", False):
        cprint("Copying dataframe to the system clipboard.", "green")
        cprint("This can be pasted into Excel, for example", "green")
        df.to_clipboard()


def abiview_fields(options):
    """Animate fields with Mayavi. Accept any file with density or potential ..."""
    from abipy.display.mvtk import MayaviFieldAnimator
    a = MayaviFieldAnimator(options.filepath)
    a.volume_animate()
    return 0


def abicomp_structure(options):
    """
    Compare crystalline structures. Use `--group` to compare for similarity.
    """
    if options.group:
        return compare_structures(options)

    paths = options.paths
    index = [os.path.relpath(p) for p in paths]

    if options.notebook:
        import nbformat
        nbv = nbformat.v4
        nb = nbv.new_notebook()

        nb.cells.extend([
            nbv.new_code_cell("""\
from __future__ import print_function, division, unicode_literals, absolute_import
import sys
import os

%matplotlib notebook
from IPython.display import display

from abipy import abilab"""),
            nbv.new_code_cell("dfs = abilab.dataframes_from_structures(%s, index=%s)" % (paths, index)),
            # Analyze dataframes.
            nbv.new_code_cell("dfs.lattice"),
            nbv.new_code_cell("dfs.coords"),
            nbv.new_code_cell("# for structure in dfs.structures: display(structure)"),
        ])

        import io, tempfile # os,
        _, nbpath = tempfile.mkstemp(prefix="abinb_", suffix='.ipynb', dir=os.getcwd(), text=True)

        # Write notebook
        import nbformat
        with io.open(nbpath, 'wt', encoding="utf8") as fh:
            nbformat.write(nb, fh)

        cmd = "jupyter notebook %s" % nbpath
        return os.system(cmd)

    dfs = abilab.dataframes_from_structures(paths, index=index,
            symprec=options.symprec, angle_tolerance=options.angle_tolerance)

    if options.ipython:
        import IPython
        IPython.embed(header="Type `dfs` in the terminal and use <TAB> to list its methods", dfs=dfs)
    else:
        #print("File list:")
        #for i, p in enumerate(paths):
        #    print("%d: %s" % (i, p))
        #print()
        print("Spglib options: symprec=", options.symprec, "angle_tolerance=", options.angle_tolerance)
        abilab.print_dataframe(dfs.lattice, title="Lattice parameters:")
        df_to_clipboard(options, dfs.lattice)

        if options.verbose:
            abilab.print_dataframe(dfs.coords, title="Atomic positions (columns give the site index):")

    return 0


def compare_structures(options):
    """Inspired to a similar function in pmg_structure."""
    paths = options.paths
    if len(paths) < 2:
        print("You need more than one structure to compare!")
        return 1

    try:
        structures = [abilab.Structure.from_file(p) for p in paths]
    except Exception as ex:
        print("Error reading structures from files. Are they in the right format?")
        print(str(ex))
        return 1

    from pymatgen.analysis.structure_matcher import StructureMatcher, ElementComparator
    compareby = "species" if options.anonymous else "element"
    m = StructureMatcher() if compareby == "species" else StructureMatcher(comparator=ElementComparator())
    print("Grouping %s structures by `%s` with `anonymous: %s`" % (len(structures), compareby, options.anonymous))

    for i, grp in enumerate(m.group_structures(structures, anonymous=options.anonymous)):
        print("Group {}: ".format(i))
        for s in grp:
            spg_symbol, international_number = s.get_space_group_info()
            print("\t- {} ({}), vol: {:.2f} A^3, {} ({})".format(
                  paths[structures.index(s)], s.formula, s.volume, spg_symbol, international_number))
        print()

    if options.verbose:
        pprint(m.as_dict())


def abicomp_spg(options):
    """
    Compare the space group found by Abinit with the spglib results
    for a set of crystalline structure(s) read from FILE(s).
    """
    try:
        structures = [abilab.Structure.from_file(p) for p in options.paths]
    except Exception as ex:
        print("Error reading structures from files. Are they in the right format?")
        print(str(ex))
        return 1

    # Remove disordered structures.
    structures = remove_disordered(structures, options.paths)

    rows, index = [], []
    symprec, angle_tolerance = options.symprec, options.angle_tolerance
    for structure in structures:
        index.append(structure.formula)
        # Call Abinit
        row = structure.abiget_spginfo(tolsym=options.tolsym, pre="abi_")
        # Call spglib.
        spglib_symbol, spglib_number = structure.get_space_group_info(symprec=symprec, angle_tolerance=angle_tolerance)
        spglib_lattice_type = structure.spget_lattice_type(symprec=symprec, angle_tolerance=angle_tolerance)
        row.update(spglib_symbol=spglib_symbol, spglib_number=spglib_number, spglib_lattice=spglib_lattice_type)
        rows.append(row)

    import pandas as pd
    df = pd.DataFrame(rows, index=index, columns=list(rows[0].keys()) if rows else None)

    print("Spglib options: symprec=", options.symprec, "angle_tolerance=", options.angle_tolerance)
    print("Abinit options: tolsym=", options.tolsym)
    print("")
    abilab.print_dataframe(df, title="Spacegroup found by Abinit and Spglib:")
    df_to_clipboard(options, df)

    return 0


def abicomp_mp_structure(options):
    """
    Compare the crystalline structure(s) read from FILE with the one(s)
    reported in the materials project database.
    """
    return _compare_with_database(options)


def abicomp_cod_structure(options):
    """
    Compare the crystalline structure(s) read from FILE with the one(s)
    given in the COD database (http://www.crystallography.net/cod).
    """
    return _compare_with_database(options)


def _compare_with_database(options):
    structures = [abilab.Structure.from_file(p) for p in options.paths]
    dbname = {"mp_structure": "materials project", "cod_structure": "COD"}[options.command]
    if options.command == "mp_structure":
        mpres = [abilab.mp_search(struct.composition.formula) for struct in structures]
    elif options.command == "cod_structure":
        mpres = [abilab.cod_search(struct.composition.formula) for struct in structures]
    else:
        raise NotImplementedError(str(options.command))

    # Filter by spglib space group number.
    if getattr(options, "same_spgnum", False):
        spgnums = [struct.get_space_group_info()[1] for struct in structures]
        mpres = [r.filter_by_spgnum(spgnum) for spgnum, r in zip(spgnums, mpres)]

    retcode = 0
    for this_structure, r in zip(structures, mpres):
        if r.structures:
            if options.notebook:
                new = r.add_entry(this_structure, "this")
                retcode += new.make_and_open_notebook(foreground=options.foreground)
            else:
                print()
                dfs = abilab.dataframes_from_structures(r.structures + [this_structure], index=r.ids + ["this"])
                abilab.print_dataframe(dfs.lattice, title="Lattice parameters:", sortby="spglib_num")
                if options.verbose:
                    abilab.print_dataframe(dfs.coords, title="Atomic positions (columns give the site index):")
                else:
                    print("Use --verbose to print atomic positions.")
                print()

        else:
            print("Couldn't find %s database entries with formula `%s`" % (dbname, this_structure.composition.formula))
            retcode += 1

    return retcode


def abicomp_xrd(options):
    """
    Compare X-ray diffraction plots (requires FILES with structure).
    """
    if len(options.paths) < 2:
        print("You need more than one structure to compare!")
        return 1

    structures = [abilab.Structure.from_file(p) for p in options.paths]

    dfs = abilab.dataframes_from_structures(structures, index=[os.path.relpath(p) for p in options.paths])
    abilab.print_dataframe(dfs.lattice, title="Lattice parameters:")
    if options.verbose:
        abilab.print_dataframe(dfs.coords, title="Atomic positions (columns give the site index):")

    from pymatgen.analysis.diffraction.xrd import XRDCalculator
    two_theta_range = tuple(float(t) for t in options.two_theta_range)
    xrd = XRDCalculator(wavelength=options.wavelength, symprec=options.symprec)
    xrd.plot_structures(structures, two_theta_range=two_theta_range, fontsize=6,
                        annotate_peaks=not options.no_annotate_peaks, tight_layout=True)
    return 0

def abicomp_data(options):
    """
    Compare results stored in multiple files with data in tabular format.
    """
    plotter = GenericDataFilesPlotter.from_files(options.paths)
    print(plotter.to_string(verbose=options.verbose))
    plotter.plot(use_index=options.use_index)
    return 0


def abicomp_ebands(options):
    """
    Plot electron bands on a grid.
    """
    paths, e0 = options.paths, options.e0
    plotter = abilab.ElectronBandsPlotter(key_ebands=[(os.path.relpath(p), p) for p in paths])

    if options.ipython:
        import IPython
        IPython.embed(header=str(plotter) + "\n\nType `plotter` in the terminal and use <TAB> to list its methods",
                      plotter=plotter)

    elif options.notebook:
        plotter.make_and_open_notebook(foreground=options.foreground)

    else:
        # Print pandas Dataframe.
        frame = plotter.get_ebands_frame()
        abilab.print_dataframe(frame)
        df_to_clipboard(options, frame)

        # Optionally, print info on gaps and their location
        if not options.verbose:
            print("\nUse --verbose for more information")
        else:
            for ebands in plotter.ebands_list:
                print(ebands)

        # Here I select the plot method to call.
        if options.plot_mode != "None":
            plotfunc = getattr(plotter, options.plot_mode, None)
            if plotfunc is None:
                raise ValueError("Don't know how to handle plot_mode: %s" % options.plot_mode)
            plotfunc(e0=e0)

    return 0


def abicomp_edos(options):
    """
    Plot electron DOSes on a grid.
    """
    paths, e0 = options.paths, options.e0
    plotter = abilab.ElectronDosPlotter(key_edos=[(os.path.relpath(p), p) for p in paths])

    if options.ipython:
        import IPython
        IPython.embed(header=str(plotter) + "\n\nType `plotter` in the terminal and use <TAB> to list its methods",
                      plotter=plotter)

    elif options.notebook:
        plotter.make_and_open_notebook(foreground=options.foreground)

    elif options.expose:
        plotter.expose(slide_mode=options.slide_mode, slide_timeout=options.slide_timeout,
                       verbose=options.verbose)

    else:
        # Optionally, print info on gaps and their location
        if not options.verbose:
            print("\nUse --verbose for more information")
        else:
            for edos in plotter.edos_list:
                print(edos)

        # Here I select the plot method to call.
        if options.plot_mode != "None":
            plotfunc = getattr(plotter, options.plot_mode, None)
            if plotfunc is None:
                raise ValueError("Don't know how to handle plot_mode: %s" % options.plot_mode)
            plotfunc(e0=e0)

    return 0


def abicomp_phbands(options):
    """
    Plot phonon bands on a grid.
    """
    paths = options.paths
    plotter = abilab.PhononBandsPlotter(key_phbands=[(os.path.relpath(p), p) for p in paths])

    if options.ipython:
        import IPython
        IPython.embed(header=str(plotter) + "\n\nType `plotter` in the terminal and use <TAB> to list its methods",
                      plotter=plotter)

    elif options.notebook:
        plotter.make_and_open_notebook(foreground=options.foreground)

    elif options.expose:
        plotter.expose(slide_mode=options.slide_mode, slide_timeout=options.slide_timeout,
                       verbose=options.verbose)
    else:
        # Print pandas Dataframe.
        frame = plotter.get_phbands_frame()
        abilab.print_dataframe(frame)

        # Optionally, print info on gaps and their location
        if not options.verbose:
            print("\nUse --verbose for more information")
        else:
            for phbands in plotter.phbands_list:
                print(phbands)

        # Here I select the plot method to call.
        if options.plot_mode != "None":
            plotfunc = getattr(plotter, options.plot_mode, None)
            if plotfunc is None:
                raise ValueError("Don't know how to handle plot_mode: %s" % options.plot_mode)
            plotfunc()

    return 0


def abicomp_phdos(options):
    """
    Compare multiple PHDOS files.
    """
    paths = options.paths
    plotter = abilab.PhononDosPlotter(key_phdos=[(os.path.relpath(p), p) for p in paths])

    if options.ipython:
        import IPython
        IPython.embed(header=str(plotter) + "\n\nType `plotter` in the terminal and use <TAB> to list its methods",
                      plotter=plotter)

    elif options.expose:
        plotter.expose(slide_mode=options.slide_mode, slide_timeout=options.slide_timeout,
                       verbose=options.verbose)

    elif options.notebook:
        plotter.make_and_open_notebook(foreground=options.foreground)

    else:
        # Optionally, print info on gaps and their location
        if not options.verbose:
            print("\nUse --verbose for more information")
        else:
            for phdos in plotter.phdos_list:
                print(phdos)

        # Here I select the plot method to call.
        if options.plot_mode != "None":
            plotfunc = getattr(plotter, options.plot_mode, None)
            if plotfunc is None:
                raise ValueError("Don't know how to handle plot_mode: %s" % options.plot_mode)
            plotfunc()

    return 0


def abicomp_getattr(options):
    """
    Extract attribute from Abipy object for all files listed on the command line and print them.
    Use `--list` option to list the attributes available (in the first file).
    """
    files = []
    attr_name = options.paths[0]
    values = []
    for p in options.paths[1:]:
        with abilab.abiopen(p) as abifile:
            if options.list:
                print("List of attributes available in %s" % p)
                pprint(dir(abifile))
                return 0

            v = getattr(abifile, attr_name)
            print(v, "   # File: ", p)
            if options.plot:
                try:
                    values.append(float(v))
                except TypeError as exc:
                    print("Cannot plot data. Exception:\n", str(exc))

    if options.plot and len(values) == len(options.paths[1:]):
        # Plot values.

        ax, fig, plt = get_ax_fig_plt()
        xs = np.arange(len(options.paths[1:]))
        ax.plot(xs, values)
        ax.set_ylabel(attr_name)
        ax.set_xticks(xs)
        xlabels = options.paths[1:]
        s = set((os.path.basename(s) for s in xlabels))
        if len(s) == len(xlabels): xlabels = s
        ax.set_xticklabels(xlabels) #, rotation='vertical')
        plt.show()

    return 0


##################
# Robot commands #
##################

def abicomp_gsr(options):
    """
    Compare multiple GSR files.
    """
    return _invoke_robot(options)


def abicomp_hist(options):
    """
    Compare multiple HIST files.
    """
    return _invoke_robot(options)


def abicomp_ddb(options):
    """
    Compare multiple DDB files. Assume DDB files with a list of q-points in the IBZ
    corresponding to homogeneous sampling i.e. files that have been merged with mrgddb.
    """
    return _invoke_robot(options)


def abicomp_anaddb(options):
    """
    Compare multiple anaddb.nc files.
    """
    return _invoke_robot(options)


def abicomp_phbst(options):
    """
    Compare multiple PHBST.nc files.
    """
    return _invoke_robot(options)


def abicomp_sigres(options):
    """
    Compare multiple SIGRES files.
    """
    return _invoke_robot(options)


def abicomp_mdf(options):
    """
    Compare macroscopic dielectric functions stored in multiple MDF files.
    """
    return _invoke_robot(options)


def abicomp_optic(options):
    """
    Compare results stored in OPTIC.nc files.
    """
    return _invoke_robot(options)


def abicomp_a2f(options):
    """
    Compare results stored in A2f.nc files.
    """
    return _invoke_robot(options)


def abicomp_gkq(options):
    """
    Compare multiple GKQ files with EPH matrix elements for a given q-point.
    """
    if options.diff:
        robot = _build_robot(options, trim_paths=True) 

        robot.plot_gkq2_diff()
    else:
        return _invoke_robot(options)


def abicomp_wrmax(options):
    """
    Compare multiple WRmax files with first order potential in real-space.
    """
    return _invoke_robot(options)


def abicomp_v1qavg(options):
    """
    Compare multiple V1QAVG files with the average of the DFPT V1 potentials as function of q-point.
    """
    return _invoke_robot(options)


def abicomp_sigeph(options):
    """
    Compare multiple SIGEPH files storing the e-ph self-energy.
    """
    return _invoke_robot(options)


def abicomp_abiwan(options):
    """
    Compare multiple ABIWAN files.
    """
    return _invoke_robot(options)


def dataframe_from_pseudos(pseudos, index=None):
    """
    Build pandas dataframe with the most important info associated to
    a list of pseudos or a list of objects that can be converted into pseudos.

    Args:
        pseudos: List of objects that can be converted to pseudos.
        index: Index of the dataframe.

    Return:
        pandas Dataframe.
    """
    from abipy.flowtk import PseudoTable
    pseudos = PseudoTable.as_table(pseudos)

    import pandas as pd
    attname = ["Z_val", "l_max", "l_local", "nlcc_radius", "xc", "supports_soc", "type"]
    rows = []
    for p in pseudos:
        row = OrderedDict([(k, getattr(p, k, None)) for k in attname])
        row["ecut_normal"], row["pawecutdg_normal"] = None, None
        if p.has_hints:
            hint = p.hint_for_accuracy(accuracy="normal")
            row["ecut_normal"] = hint.ecut
            if hint.pawecutdg: row["pawecutdg_normal"] = hint.pawecutdg
        rows.append(row)

    return pd.DataFrame(rows, index=index, columns=list(rows[0].keys()) if rows else None)


def abicomp_pseudos(options):
    """"Compare multiple pseudos Print table to terminal."""
    # Make sure entries in index are unique.
    index = [os.path.basename(p) for p in options.paths]
    if len(index) != len(set(index)): index = [os.path.relpath(p) for p in options.paths]
    df = dataframe_from_pseudos(options.paths, index=index)
    abilab.print_dataframe(df, sortby="Z_val")
    return 0


def _build_robot(options, trim_paths=False):
    """Build robot instance from CLI options."""
    robot_cls = abilab.Robot.class_for_ext(options.command.upper())

    # To define an Help action
    # http://stackoverflow.com/questions/20094215/argparse-subparser-monolithic-help-output?rq=1
    paths = options.paths
    if options.verbose > 1:
        print("In _invoke_robot with paths", paths)

    if os.path.isdir(paths[0]):
        # Assume directory.
        robot = robot_cls.from_dir(top=paths[0], walk=not options.no_walk)
    else:
        # Assume file.
        robot = robot_cls.from_files([paths[0]])

    if len(paths) > 1:
        # Handle multiple arguments. Could be either other directories or files.
        for p in paths[1:]:
            if os.path.isdir(p):
                robot.scan_dir(top=p, walk=not options.no_walk)
            elif os.path.isfile(p):
                robot.add_file(os.path.abspath(p), p)
            else:
                cprint("Ignoring %s. Neither file or directory." % str(p), "red")

    if len(robot) == 0:
        raise RuntimeError("Empty robot --> No file associated to this robot has been found")

    if trim_paths: robot.trim_paths()
    return robot


def _invoke_robot(options):
    """
    Analyze multiple files with a robot. Support list of files and/or
    list of directories passed on the CLI..

    By default, the script with call `robot.to_string(options.verbose)` to print info to terminal.
    For finer control, use --ipy to start an ipython console to interact with the robot directly
    or --nb to generate a jupyter notebook.
    """
    robot = _build_robot(options)

    if options.notebook:
        robot.make_and_open_notebook(foreground=options.foreground)

    elif options.print or options.expose:
        robot.trim_paths()
        # Print dataframe if robot provides get_dataframe method.
        if hasattr(robot, "get_dataframe"):
            try:
                df = robot.get_dataframe()
                abilab.print_dataframe(df, title="Output of robot.get_dataframe():")
                df_to_clipboard(options, df)

            except Exception as exc:
                cprint("Exception:\n%s\n\nwhile invoking get_dataframe. Falling back to to_string" % str(exc), "red")
                print(robot.to_string(verbose=options.verbose))

        else:
            cprint("%s does not provide `get_dataframe` method. Using `to_string`" % (
                    robot.__class__.__name__), "yellow")
            print(robot.to_string(verbose=options.verbose))

        if not options.verbose:
            print("\nUse --verbose for more information")

        if options.expose and hasattr(robot, "expose"):
            robot.expose(slide_mode=options.slide_mode, slide_timeout=options.slide_timeout,
                         verbose=options.verbose)

    #elif options.ipython:
    else:
        import IPython
        robot.trim_paths()
        IPython.embed(header=repr(robot) + "\n\nType `robot` in the terminal and use <TAB> to list its methods",
                      robot=robot)

    return 0


def abicomp_gs_scf(options):
    """
    Compare ground-state SCF cycles.
    """
    paths = options.paths
    f0 = abilab.AbinitOutputFile(paths[0])
    figures = f0.compare_gs_scf_cycles(paths[1:])
    if not figures:
        cprint("Cannot find GS-SCF sections in output files.", "yellow")
    return 0


def abicomp_dfpt2_scf(options):
    """
    Compare DFPT SCF cycles.
    """
    paths = options.paths
    f0 = abilab.AbinitOutputFile(paths[0])
    figures = f0.compare_d2de_scf_cycles(paths[1:])
    if not figures:
        cprint("Cannot find DFPT-SCF sections in output files.", "yellow")
    return 0


def abicomp_text(options):
    """
    Compare 2+ text files in the browser
    """
    from abipy.tools.devtools import HtmlDiff
    return HtmlDiff(options.paths).open_browser(diffmode=options.diffmode)


def abicomp_time(options):
    """
    Analyze/plot the timing data of single or multiple runs.
    """
    paths = options.paths
    from abipy.abio.timer import AbinitTimerParser

    if len(options.paths) == 1 and os.path.isdir(paths[0]):
        # Scan directory tree
        top = options.paths[0]
        print("Walking directory tree from top:", top, "Looking for file extension:", options.ext)
        parser, paths_found, okfiles = AbinitTimerParser.walk(top=top, ext=options.ext)

        if not paths_found:
            cprint("Empty file list!", color="magenta")
            return 1

        print("Found %d files\n" % len(paths_found))
        if okfiles != paths_found:
            badfiles = [f for f in paths_found if f not in okfiles]
            cprint("Cannot parse timing data from the following files:", color="magenta")
            for bad in badfiles: print(bad)

    else:
        # Parse list of files.
        parser = AbinitTimerParser()
        okfiles = parser.parse(options.paths)

        if okfiles != options.paths:
            badfiles = [f for f in options.paths if f not in okfiles]
            cprint("Cannot parse timing data from the following files:", color="magenta")
            for bad in badfiles: print(bad)

    if parser is None:
        cprint("Cannot analyze timing data. parser is None", color="magenta")
        return 1

    print(parser.summarize())

    if options.verbose:
        for timer in parser.timers():
            print(timer.get_dataframe())

    if options.ipython:
        cprint("Invoking ipython shell. Use parser to access the object inside ipython", color="blue")
        import IPython
        IPython.start_ipython(argv=[], user_ns={"parser": parser})
    elif options.notebook:
        parser.make_and_open_notebook(foreground=options.foreground)
    else:
        parser.plot_all()

    return 0


def get_epilog():
    return """\
Usage example:

############
# Structures
############

  abicomp.py structure */*/outdata/out_GSR.nc   => Compare structures in multiple files.
                                                   Use `--group` to compare for similarity.
  abicomp.py spg si.cif out_GSR.nc              => Compare spacegroup(s) found by Abinit and spglib.
  abicomp.py hist FILE(s)                       => Compare final structures read from HIST.nc files.
  abicomp.py mp_structure FILE(s)               => Compare structure(s) read from FILE(s) with the one(s)
                                                   given in the materials project database.
  abicomp.py cod_structure FILE(s)              => Compare structure(s) read from FILE(s) with the one(s)
                                                   given in the COD database (http://www.crystallography.net/cod).
  abicomp.py xrd *.cif *.GSR.nc                 => Compare X-ray diffraction plots (requires FILES with structure).

###########
# Electrons
###########

  abicomp.py ebands out1_GSR.nc out2_WFK.nc     => Plot electron bands on a grid (Use `-p` to change plot mode)
  abicomp.py ebands *_GSR.nc -ipy               => Build plotter object and start ipython console.
  abicomp.py ebands *_GSR.nc -nb                => Interact with the plotter in the jupyter notebook.
  abicomp.py ebands `find . -name "*_GSR.nc"` -c = Find all GSR.nc files startign from current working directory
                                                   Copy dataframe to the system clipboard.
                                                   This can be pasted into Excel, for example
  abicomp.py edos *_WFK.nc -nb                  => Compare electron DOS in the jupyter notebook.
  abicomp.py optic DIR -nb                      => Compare optic results in the jupyter notebook.
  abicomp.py abiwan *_ABIWAN.nc --expose        => Compare ABIWAN results, produce matplotlib figures.

#########
# Phonons
#########

  abicomp.py phbands *_PHBST.nc -nb             => Compare phonon bands in the jupyter notebook.
  abicomp.py phbst *_PHBST.nc -ipy              => Compare phonon bands with robot in ipython terminal.
  abicomp.py phdos *_PHDOS.nc -nb               => Compare phonon DOSes in the jupyter notebook.
  abicomp.py ddb outdir1 outdir2 out_DDB -nb    => Analyze all DDB files in directories outdir1, outdir2 and out_DDB file.

###############
# Anaddb netcdf
###############

  abicomp anaddb tutorespfn_telast_2-telast_3/anaddb.nc

#########
# E-PH
#########

  abicomp.py a2f *_A2F.nc -nb                   => Compare A2f results in the jupyter notebook.
  abicomp.py sigeph *_SIGEPH.nc -nb             => Compare Fan-Migdal self-energy in the jupyter notebook.
  abicomp.py gkq out1_GKQ.nc out1_GKQ.nc -d     => Plot difference between matrix elements (supports 2+ files).
  abicomp.py v1qavg out_V1QAVG.nc               => Compare V1QAVG files.

########
# GW/BSE
########

  abicomp.py sigres *_SIGRES.nc                 => Compare multiple SIGRES files.
  abicomp.py mdf *_MDF.nc --seaborn             => Compare macroscopic dielectric functions.
                                                   Use seaborn settings.

###############
# Miscelleanous
###############

  abicomp.py data FILE1 FILE2 ...                 => Read data from files with results in tabular format and
                                                     compare results. Mainly used for text files without any schema.
  abicomp.py getattr energy *_GSR.nc              => Extract the `energy` attribute from a list of GSR files
                                                     and print results. Use `--list` to get list of possible names.
  abicomp.py pseudos PSEUDO_FILES                 => Compare pseudopotential files.

############
# Text files
############

  abicomp.py gs_scf run1.abo run2.abo             => Compare the SCF cycles in two output files.
  abicomp.py dfpt2_scf run1.abo run2.abo          => Compare the DFPT SCF cycles in two output files.
  abicomp.py.py time [OUT_FILES]                  => Parse timing data in files and plot results
  abicomp.py.py time . --ext=abo                  => Scan directory tree from `.`, look for files with extension `abo`
                                                     parse timing data and plot results.
  abicomp.py text run1.abo run2.abo               => Produce diff of 2+ text files in the browser.

TIP:

The python code operates on a list of files/directories passed via the command line interface.
The arguments are interpreted by the shell before invoking the script.
This means that one can use the bash syntax and the unix tools to precompute the list of files/directories.

For example, one can use Unix find to select all files with the a given extension and pass them to abicomp.py.
For command:

    abicomp.py structure `find . -name "*_GSR.nc"`

will compare the structures extracted from all GSR.nc files found within the current working directory (note backticks).

Also, remember that in bash {#..#} generates a sequence of numbers or chars, similarly to range() in Python
For instance:

    {1..5} --> 1 2 3 4 5

and this trick can be used to select files/directories as in:

    abicomp.py structure w1/t{2..4}/outdata/*_GSR.nc

The [!name] syntax can be used to exclude patterns, so

    abicomp.py structure w1/t[!2]*/outdata/*_GSR.nc

excludes all the GSR.nc files in the t2/outdata directory.

See also http://wiki.bash-hackers.org/syntax/pattern

NOTE: The `gsr`, `ddb`, `sigres`, `mdf` commands use robots to analyze files.
In this case, one can provide a list of files and/or list of directories on the command-line interface e.g.:

    $ abicomp.py ddb dir1 out_DDB dir2

Directories will be scanned recursively to find files with the extension associated to the robot, e.g.
`abicompy.py mdf .` will read all *_MDF.nc files inside the current directory including sub-directories (if any).
Use --no-walk to ignore sub-directories when robots are used.

Use `abicomp.py --help` for help and `abicomp.py COMMAND --help` to get the documentation for `COMMAND`.
Use `-v` to increase verbosity level (can be supplied multiple times e.g -vv).
"""


def get_parser(with_epilog=False):

    # Parent parser for common options.
    copts_parser = argparse.ArgumentParser(add_help=False)
    copts_parser.add_argument('paths', nargs="+", help="List of files to compare.")
    copts_parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
        help='Verbose, can be supplied multiple times to increase verbosity.')
    copts_parser.add_argument('-sns', "--seaborn", const="paper", default=None, action='store', nargs='?', type=str,
        help='Use seaborn settings. Accept value defining context in ("paper", "notebook", "talk", "poster"). Default: paper')
    copts_parser.add_argument('--pylustrator', action='store_true', default=False,
        help="Style matplotlib plots with pylustrator. See https://pylustrator.readthedocs.io/en/latest/")
    copts_parser.add_argument('-mpl', "--mpl-backend", default=None,
        help=("Set matplotlib interactive backend. "
              "Possible values: GTKAgg, GTK3Agg, GTK, GTKCairo, GTK3Cairo, WXAgg, WX, TkAgg, Qt4Agg, Qt5Agg, macosx."
              "See also: https://matplotlib.org/faq/usage_faq.html#what-is-a-backend."))
    copts_parser.add_argument('--loglevel', default="ERROR", type=str,
        help="Set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG.")

    # Parent parser for commands calling spglib.
    spgopt_parser = argparse.ArgumentParser(add_help=False)
    spgopt_parser.add_argument('--symprec', default=1e-3, type=float,
        help="""\
symprec (float): Tolerance for symmetry finding. Defaults to 1e-3,
which is fairly strict and works well for properly refined structures with atoms in the proper symmetry coordinates.
For structures with slight deviations from their proper atomic positions (e.g., structures relaxed with electronic structure
codes), a looser tolerance of 0.1 (the value used in Materials Project) is often needed.""")
    spgopt_parser.add_argument('--angle-tolerance', default=5.0, type=float,
        help="angle_tolerance (float): Angle tolerance for symmetry finding. Default: 5.0")
    #spgopt_parser.add_argument("--no-time-reversal", default=False, action="store_true", help="Don't use time-reversal.")

    # Parent parser for commands that operating on pandas dataframes
    pandas_parser = argparse.ArgumentParser(add_help=False)
    pandas_parser.add_argument("-c", '--clipboard', default=False, action="store_true",
            help="Copy dataframe to the system clipboard. This can be pasted into Excel, for example")

    # Parent parser for commands supporting (ipython/jupyter)
    ipy_parser = argparse.ArgumentParser(add_help=False)
    ipy_parser.add_argument('-nb', '--notebook', default=False, action="store_true", help='Generate jupyter notebook.')
    ipy_parser.add_argument('--foreground', action='store_true', default=False,
        help="Run jupyter notebook in the foreground.")
    ipy_parser.add_argument('-ipy', '--ipython', default=False, action="store_true", help='Invoke ipython terminal.')

    # Parent parser for commands supporting (jupyter notebooks)
    nb_parser = argparse.ArgumentParser(add_help=False)
    nb_parser.add_argument('-nb', '--notebook', default=False, action="store_true", help='Generate jupyter notebook.')
    nb_parser.add_argument('--foreground', action='store_true', default=False,
        help="Run jupyter notebook in the foreground.")

    # Parent parser for commands supporting expose
    expose_parser = argparse.ArgumentParser(add_help=False)
    expose_parser.add_argument("-e", '--expose', default=False, action="store_true",
            help='Execute robot.expose to produce a pre-defined list of matplotlib figures.')
    expose_parser.add_argument("-s", "--slide-mode", default=False, action="store_true",
            help="Used if --expose to iterate over figures. Expose all figures at once if not given on the CLI.")
    expose_parser.add_argument("-t", "--slide-timeout", type=int, default=None,
            help="Close figure after slide-timeout seconds (only if slide-mode). Block if not specified.")

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=get_epilog() if with_epilog else "",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-V', '--version', action='version', version=abilab.__version__)

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for structure command.
    p_struct = subparsers.add_parser('structure', parents=[copts_parser, ipy_parser, spgopt_parser, pandas_parser],
            help=abicomp_structure.__doc__)
    p_struct.add_argument("-g", "--group", default=False, action="store_true",
        help="Compare a set of structures for similarity.")
    p_struct.add_argument("-a", "--anonymous", default=False, action="store_true",
        help="Whether to use anonymous mode in StructureMatcher. Default False")

    # Subparser for spg command.
    p_spg = subparsers.add_parser('spg', parents=[copts_parser, spgopt_parser, pandas_parser],
            help=abicomp_spg.__doc__)
    p_spg.add_argument("-t", "--tolsym", type=float, default=None, help="""\
Gives the tolerance on the atomic positions (reduced coordinates), primitive vectors, or magnetization,
to be considered equivalent, thanks to symmetry operations. This is used in the recognition of the set
of symmetries of the system, or the application of the symmetry operations to generate from a reduced set of atoms,
the full set of atoms. Note that a value larger than 0.01 is considered to be unacceptable.""")

    # Subparser for mp_structure command.
    p_mpstruct = subparsers.add_parser('mp_structure', parents=[copts_parser, nb_parser],
        help=abicomp_mp_structure.__doc__)
    p_mpstruct.add_argument("--same-spgnum", default=False, action="store_true",
        help="Select only MP structures with same space group number as input structure.")

    # Subparser for cod_structure command.
    p_codstruct = subparsers.add_parser('cod_structure', parents=[copts_parser, nb_parser],
        help=abicomp_cod_structure.__doc__)
    #p_codstruct.add_argument("--same-spgnum", default=False, action="store_true",
    #    help="Select only COD structures with same space group number as input structure.")

    # Subparser for xrd.
    p_xrd = subparsers.add_parser('xrd', parents=[copts_parser], help="Compare X-ray diffraction plots.")
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

    # Subparser for data command.
    p_data = subparsers.add_parser('data', parents=[copts_parser, expose_parser], help=abicomp_data.__doc__)
    p_data.add_argument("-i", "--use-index", default=False, action="store_true",
        help="Use the row index as x-value in the plot. By default the plotter uses the first column as x-values")

    # Subparser for ebands command.
    p_ebands = subparsers.add_parser('ebands', parents=[copts_parser, ipy_parser, pandas_parser],
            help=abicomp_ebands.__doc__)
    p_ebands.add_argument("-p", "--plot-mode", default="gridplot",
        choices=["gridplot", "combiplot", "boxplot", "combiboxplot", "plot_band_edges", "animate", "None"],
        help="Plot mode e.g. `-p combiplot` to plot bands on the same figure. Default is `gridplot`.")
    p_ebands.add_argument("-e0", default="fermie", choices=["fermie", "None"],
        help="Option used to define the zero of energy in the band structure plot. Default is `fermie`.")

    # Subparser for edos command.
    p_edos = subparsers.add_parser('edos', parents=[copts_parser, ipy_parser, expose_parser],
        help=abicomp_edos.__doc__)
    p_edos.add_argument("-p", "--plot-mode", default="gridplot",
        choices=["gridplot", "combiplot", "None"],
        help="Plot mode e.g. `-p combiplot` to plot DOSes on the same figure. Default is `gridplot`.")
    p_edos.add_argument("-e0", default="fermie", choices=["fermie", "None"],
        help="Option used to define the zero of energy in the DOS plot. Default is `fermie`.")

    # Subparser for phbands command.
    p_phbands = subparsers.add_parser('phbands', parents=[copts_parser, ipy_parser, expose_parser],
        help=abicomp_phbands.__doc__)
    p_phbands.add_argument("-p", "--plot-mode", default="gridplot",
        choices=["gridplot", "combiplot", "boxplot", "combiboxplot", "animate", "None"],
        help="Plot mode e.g. `-p combiplot` to plot bands on the same figure. Default is `gridplot`.")

    # Subparser for phdos command.
    p_phdos = subparsers.add_parser('phdos', parents=[copts_parser, ipy_parser, expose_parser],
        help=abicomp_phdos.__doc__)
    p_phdos.add_argument("-p", "--plot-mode", default="gridplot",
        choices=["gridplot", "combiplot", "None"],
        help="Plot mode e.g. `-p combiplot` to plot DOSes on the same figure. Default is `gridplot`.")

    # Subparser for getattr command.
    p_getattr = subparsers.add_parser('getattr', parents=[copts_parser], help=abicomp_getattr.__doc__)
    p_getattr.add_argument('--plot', default=False, action="store_true", help="Plot data with matplotlib (requires floats).")
    p_getattr.add_argument('--list', default=False, action="store_true", help="Print attributes available in file")

    # Subparser for robot commands
    # Use own version of ipy_parser with different default values.
    robot_ipy_parser = argparse.ArgumentParser(add_help=False)
    robot_ipy_parser.add_argument('-nb', '--notebook', default=False, action="store_true", help='Generate jupyter notebook.')
    robot_ipy_parser.add_argument('--foreground', action='store_true', default=False,
        help="Run jupyter notebook in the foreground.")
    #robot_ipy_parser.add_argument('-ipy', '--ipython', default=True, action="store_true", help='Invoke ipython terminal.')
    robot_ipy_parser.add_argument('-p', '--print', default=False, action="store_true", help='Print robot and return.')

    # Parent parser for *robot* commands
    robot_parser = argparse.ArgumentParser(add_help=False)
    robot_parser.add_argument('--no-walk', default=False, action="store_true", help="Don't enter subdirectories.")

    robot_parents = [copts_parser, robot_ipy_parser, robot_parser, expose_parser, pandas_parser]
    p_gsr = subparsers.add_parser('gsr', parents=robot_parents, help=abicomp_gsr.__doc__)
    p_hist = subparsers.add_parser('hist', parents=robot_parents, help=abicomp_hist.__doc__)
    p_ddb = subparsers.add_parser('ddb', parents=robot_parents, help=abicomp_ddb.__doc__)
    p_anaddb = subparsers.add_parser('anaddb', parents=robot_parents, help=abicomp_anaddb.__doc__)
    p_phbst = subparsers.add_parser('phbst', parents=robot_parents, help=abicomp_phbst.__doc__)
    p_sigres = subparsers.add_parser('sigres', parents=robot_parents, help=abicomp_sigres.__doc__)
    p_mdf = subparsers.add_parser('mdf', parents=robot_parents, help=abicomp_mdf.__doc__)
    p_optic = subparsers.add_parser('optic', parents=robot_parents, help=abicomp_optic.__doc__)
    p_a2f = subparsers.add_parser('a2f', parents=robot_parents, help=abicomp_a2f.__doc__)
    p_sigeph = subparsers.add_parser('sigeph', parents=robot_parents, help=abicomp_sigeph.__doc__)
    p_gkq = subparsers.add_parser('gkq', parents=robot_parents, help=abicomp_gkq.__doc__)
    p_gkq.add_argument('-d', '--diff', default=False, action="store_true", help='Plot difference between eph matrix elements.')
    p_v1qavg = subparsers.add_parser('v1qavg', parents=robot_parents, help=abicomp_v1qavg.__doc__)
    p_wrmax = subparsers.add_parser('wrmax', parents=robot_parents, help=abicomp_wrmax.__doc__)
    p_abiwan = subparsers.add_parser('abiwan', parents=robot_parents, help=abicomp_abiwan.__doc__)

    # Subparser for pseudos command.
    p_pseudos = subparsers.add_parser('pseudos', parents=[copts_parser], help=abicomp_pseudos.__doc__)

    # Subparser for time command.
    p_time = subparsers.add_parser('time', parents=[copts_parser, ipy_parser], help=abicomp_time.__doc__)
    p_time.add_argument("-e", "--ext", default=".abo", help=("File extension for Abinit output files. "
                        "Used when the first argument is a directory. Default is `.abo`"))

    # Subparser for gs_scf command.
    p_gs_scf = subparsers.add_parser('gs_scf', parents=[copts_parser], help=abicomp_gs_scf.__doc__)

    # Subparser for dfpt2_scf command.
    p_dftp2_scf = subparsers.add_parser('dfpt2_scf', parents=[copts_parser], help=abicomp_dfpt2_scf.__doc__)

    # Subparser for text command.
    p_text = subparsers.add_parser('text', parents=[copts_parser], help=abicomp_text.__doc__)
    p_text.add_argument("-d", "--diffmode", default="difflib", help=("Select diff application. "
        "Possible values: difflib (default), pygmentize (requires package)."))

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

    if options.mpl_backend is not None:
        # Set matplotlib backend
        import matplotlib
        matplotlib.use(options.mpl_backend)

    if options.seaborn:
        # Use seaborn settings.
        import seaborn as sns
        sns.set(context=options.seaborn, style='darkgrid', palette='deep',
                font='sans-serif', font_scale=1, color_codes=False, rc=None)

    if options.pylustrator:
        # Start pylustrator to style matplotlib plots 
        import pylustrator
        pylustrator.start()

    if options.verbose > 2:
        print(options)

    # Dispatch
    return globals()["abicomp_" + options.command](options)


if __name__ == "__main__":
    sys.exit(main())
