#!/usr/bin/env python
"""
Script to generate movies from Abinit output files.
"""
from __future__ import unicode_literals, division, print_function, absolute_import

import sys
import os
import argparse
import numpy as np

from monty.functools import prof_main
from abipy import abilab


def sort_paths(options):
    """
    Sort input files whose name is in the form `out_TIM2_DEN`
    Files are sorted by TIM index.
    """
    if options.no_sort: return
    names = [os.path.basename(p) for p in options.paths]
    import re
    # out_TIM2_DEN
    tim = re.compile(r".+_TIM(\d+)_.+")
    l = []
    for p, n in zip(options.paths, names):
        m = tim.match(n)
        if m:
            l.append((int(m.group(1)), p))
    if not l: return
    if len(l) != len(options.paths):
        print("Cannot sort input paths!")
        return

    options.paths = [t[1] for t in sorted(l, key=lambda t: t[0])]
    print("Input files have been automatically sorted")
    for i, p in enumerate(options.paths):
        print("%d: %s" % (i, p))
    print("Use --no-sort to disable automatic sorting.")


def abimovie_hist(options):
    """Visualize structural relaxation/molecular dynamics run from data stored in the HIST.nc file.
    Requires mayavi."""
    for path in options.paths:
        with abilab.abiopen(path) as hist:
            print(hist)
            if options.trajectories:
                hist.mvplot_trajectories()
            else:
                hist.mvanimate()
    return 0


def abimovie_ebands(options):
    """Animate electron bands. Accept any file with ElectronBands e.g. GSR.nc, WFK.nc, ..."""
    #paths, e0 = options.paths, "fermie" # options.e0
    plotter = abilab.ElectronBandsPlotter(key_ebands=[(os.path.relpath(p), p) for p in options.paths])
    plotter.animate()
    return 0

    # DRAW A FIGURE WITH MATPLOTLIB
    ebands_list = [abilab.ElectronBands.from_file(p) for p in options.paths]
    import matplotlib.pyplot as plt
    from moviepy.video.io.bindings import mplfig_to_npimage
    import moviepy.editor as mpy
    # See also http://zulko.github.io/blog/2014/11/29/data-animations-with-python-and-moviepy/
    duration = 4
    #fig_mpl, ax = plt.subplots(1,figsize=(5,3), facecolor='white')
    #xx = np.linspace(-2,2,200) # the x vector
    #zz = lambda d: np.sinc(xx**2)+np.sin(xx+d) # the (changing) z vector
    #ax.set_title("Elevation in y=0")
    #ax.set_ylim(-1.5,2.5)
    #line, = ax.plot(xx, zz(0), lw=3)

    # ANIMATE WITH MOVIEPY (UPDATE THE CURVE FOR EACH t). MAKE A GIF.
    def make_frame_mpl(t):
        t = int(t % len(ebands_list))
        fig_mpl = ebands_list[t].plot(show=False)
        return mplfig_to_npimage(fig_mpl) # RGB image of the figure

    animation = mpy.VideoClip(make_frame_mpl, duration=duration)
    animation.write_gif("sinc_mpl.gif", fps=20)
    return 0


def abimovie_phbands(options):
    """Animate phonon bands. Accept any file with PhononBands e.g. PHBST.nc, ..."""
    plotter = abilab.PhononBandsPlotter(key_phbands=[(os.path.relpath(p), p) for p in options.paths])
    plotter.animate()
    return 0


def abimovie_fields(options):
    """Animate fields with Mayavi. Accept any file with densities, potentials ..."""
    from abipy.display.mvtk import MayaviFieldAnimator
    sort_paths(options)
    a = MayaviFieldAnimator(options.paths)
    a.volume_animate()
    return 0

@prof_main
def main():

    def str_examples():
        return """\
Usage example:

###########
# Structure
###########

    abimovie.py hist HIST_FILE(s)    ==> Visualize structural relaxation/molecular dynamics
                                         run from data stored in the HIST.nc file.

###########
# Electrons
###########

    abimovie.py ebands *_GSR.nc      ==>   Animate electron bands.
    abimovie.py fields *_DEN.nc      ==>   Animate densities on FFT mesh.

#########
# Phonons
#########

    abimovie.py phbands *_PHBST.nc   ==>   Animate phonon bands.

Use `abimovie.py --help` for help and `abimovie.py COMMAND --help` to get the documentation for `COMMAND`.
Use `-v` to increase verbosity level (can be supplied multiple times e.g -vv).
"""

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg:
            sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    # Parent parser for common options.
    copts_parser = argparse.ArgumentParser(add_help=False)
    copts_parser.add_argument('paths', nargs="+", help="List of files to analyze.")

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--loglevel', default="ERROR", type=str,
                         help="Set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")
    parser.add_argument('-V', '--version', action='version', version=abilab.__version__)
    parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                         help='verbose, can be supplied multiple times to increase verbosity.')
    parser.add_argument('--seaborn', action="store_true", help="Use seaborn settings.")

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for hist command.
    p_hist = subparsers.add_parser('hist', parents=[copts_parser], help=abimovie_hist.__doc__)
    p_hist.add_argument("-t", "--trajectories", default=False, action="store_true", help="Plot trajectories.")

    # Subparser for ebands command.
    p_ebands = subparsers.add_parser('ebands', parents=[copts_parser], help=abimovie_ebands.__doc__)

    # Subparser for phbands command.
    p_phbands = subparsers.add_parser('phbands', parents=[copts_parser], help=abimovie_phbands.__doc__)

    # Subparser for fields command.
    p_fields = subparsers.add_parser('fields', parents=[copts_parser], help=abimovie_fields.__doc__)
    p_fields.add_argument('--no-sort', default=False, action="store_true", help="Disable automatic sorting of filepaths.")

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

    if options.verbose > 1:
        print(options)



    # Dispatch
    return globals()["abimovie_" + options.command](options)

    return 0

if __name__ == "__main__":
    sys.exit(main())
