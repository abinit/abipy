#!/usr/bin/env python
"""Script for plotting Bethe-Salpeter results obtained with ABINIT."""
from __future__ import print_function, division

import sys
import matplotlib.pyplot as plt

from matplotlib.gridspec import GridSpec
from abipy.electrons import DielectricFunction

__version__ = "0.1"
__author__ = "Matteo Giantomassi"

#########################################################################################

def plot_dielectric_function(df_list, what):
    """Plot the macroscopic dielectric function.

       :arg df_list: list of strings, file objects or DielectricFunction instances.
       :arg what: ("r","i","b") "r" for plotting only the real part, "i" for the
       imaginary part, "b" for both.
    """

    if what not in ["r", "i", "b"]:
        raise ValueError("Wrong value for what: " + str(what))

    ylabels = {"r": "$\epsilon_1$", "i": "$\epsilon_2$"}
    reim = {"r": "e1q_avg", "i": "e2q_avg"}

    fig = plt.figure()
    #
    # Instanciate the axes for the plots.
    if what == "b":
        kws = ("r", "i")
    else:
        kws = (what)
    nplots = len(kws)

    axes = dict()
    for idx_plot, k in enumerate(kws):
        axes[k] = fig.add_subplot(nplots, 1, idx_plot)
        axes[k].grid(True)
        axes[k].set_xlabel('$\omega$ [eV]')
        axes[k].set_ylabel(ylabels[k])

    for df in df_list:
        if not isinstance(df, DielectricFunction):
            df = DielectricFunction(df)

        for k in kws:
            what = reim[k]
            df.plot_ax(axes[k], what) #, *args, **kwargs):

    for k in kws:
        axes[k].legend(loc="upper right")

    plt.show()

#########################################################################################

def plot_driver(options, args):
    df_files = options.df_files

    if df_files:
        if options.verbose: print
        "Plotting MDF files"
        plot_dielectric_function(df_files, "b")

    else:
        show_examples_and_exit(1)

#########################################################################################

def show_examples_and_exit(error_code=0, err_msg=None):
    """Display the usage of the script."""
    examples = """
    Typical examples:
    \n
    Author %(__author__)s
    \n""" % globals()
    sys.stderr.write(examples)
    if err_msg:  sys.stderr.write(err_msg + "\n")
    sys.exit(error_code)


def main():
    from optparse import OptionParser

    usage = "usage: %prog [options]"
    version = "%prog " + str(__version__)
    parser = OptionParser(usage=usage, version=version)

    parser.add_option("-a", "--action", dest="action", type="string", default="plot",
                      help="action to be performed, defaults to 'plot'")

    parser.add_option("-d", "--dielectric-function", dest="df_files", action="append",
                      help="read dielectric function from FILE", metavar="FILE")

    #parser.add_option("-d", "--dos", dest="dos_file", type="string",
    #                  help="read DOS from FILE", metavar="FILE")

    parser.add_option("-v", "--verbose", help="verbose mode",
                      action="store_true", dest="verbose")

    parser.add_option("-q", "--quiet", help="disable verbose mode",
                      action="store_false", dest="verbose")

    (options, args) = parser.parse_args()

    # One argument is required.
    #if len(args) != 1:
    #  parser.error("incorrect number of arguments.")

    children = dict()
    children["plot"] = plot_driver

    action = options.action

    if children.has_key(action):
        children[action](options, args)
    else:
        sys.stderr.write("Unknown value for action: " + str(action))
        show_examples_and_exit(1)

    return 0

#########################################################################################

if __name__ == "__main__":
    sys.exit(main())
