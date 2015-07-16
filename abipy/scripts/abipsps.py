#!/usr/bin/env python
"""
Script to plot pseudopotential data and/or compare mulitple pseudopotentials.
It invokes Abinit to produced the PSPS.nc files with formactors, model core charges 
and other quanties used by Abinit to apply the non-local part of the KS Hamiltonian.
"""
from __future__ import print_function, division, unicode_literals

import sys
import argparse

from abipy.electrons.psps import PspsFile, compare_pseudos
from abipy import abilab


def main():
    def str_examples():
        return  """\
Usage example:\n

    abipsps.py pseudo            => Visualize data relative to a single pseudo.
    abipsps.py pseudo1 pseudo2   => Compare pseudo1 with pseudo1 (accept arbitrary number of pseudos).

"""
    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg: sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--loglevel', default="ERROR", type=str,
                        help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    parser.add_argument("-e", '--ecut', type=float, default=10, help="Cutoff energy in Hartree (default: 10).")
    parser.add_argument('--seaborn', action="store_true", help="Use seaborn settings")

    parser.add_argument('pseudos', nargs="+", help="Pseudo or list of pseudopotential files")

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

    if options.seaborn:
        import seaborn as sns
        #sns.set(style='ticks', palette='Set2')
        sns.set(style="dark", palette="Set2")
        #And to remove "chartjunk", do:
        #sns.despine()
        #plt.tight_layout()
        #sns.despine(offset=10, trim=True)

    pseudos = options.pseudos
    what = "all"

    if len(pseudos) > 1:
      compare_pseudos(pseudos, options.ecut)

    else:
        pseudo = abilab.Pseudo.from_file(pseudos[0])
        with pseudo.open_pspsfile(ecut=options.ecut) as psps:
            psps.plot(what="all")
            #psps.plot_modelcore_rspace()
            #psps.plot_modelcore_qspace()
            #psps.plot_vlspl()
            #psps.plot_ffspl()

    return 0


if __name__ == "__main__":
  sys.exit(main())
