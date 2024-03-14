#!/usr/bin/env python
"""
This scripts runs all the python scripts located in this directory allowing
the user to change the matplotlib backend.

Usage:
    _runplots.py [backend]
"""
from __future__ import annotations

import sys
import os
import time
import argparse

from subprocess import call, Popen


def str_examples():
    examples = """
      Usage example:\n\n
      _runplots.py               => Run all scripts.
      _runplots.py -m auto -t 5  => Run all scripts, kill the process after 5 seconds.
    """
    return examples


def show_examples_and_exit(err_msg=None, error_code=1):
    """Display the usage of the script."""
    sys.stderr.write(str_examples())
    if err_msg:
        sys.stderr.write("Fatal Error\n" + err_msg + "\n")
    sys.exit(error_code)


def main():
    parser = argparse.ArgumentParser(epilog=str_examples(),formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-b', '--backend', type=str, default="Agg",
                        help="matplotlib backend e.g. Agg for non-graphical display.")
    parser.add_argument('-m', '--mode', type=str, default="automatic",
                        help="execution mode. Either s (sequential) or a (automatic)")
    parser.add_argument('-t', '--time', type=float, default=8,
                        help="wait time seconds before running next demo.")
    #parser.add_argument('-p', '--ply-show', type=bool, default=False,
    #                     help="Show plotly figures in browser.")

    options = parser.parse_args()

    import matplotlib
    if options.backend:
        print("Using matplotlib backend: ", options.backend)
        matplotlib.use(options.backend)
    #change_matplotlib_backend(new_backend=options.backend)

    #from abipy.tools.plotting import set_plotly_default_show
    #print("Setting plotly_default_show to: ", options.ply_show)
    #set_plotly_default_show(options.ply_show)

    # Find scripts.
    root = os.path.join(os.path.dirname(__file__), "plot")
    scripts = []
    for fname in os.listdir(root):
        if fname.endswith(".py") and fname.startswith("plot"):
            scripts.append(os.path.join(root, fname))

    #env = {
    #    "MPLBACKEND":  options.backend,
    #}
    env = None
    print(f"{env=}")

    # Run scripts according to mode.
    if options.mode in ["s", "sequential"]:
        for script in scripts:
            print("About to execute", script)
            args = ["python", script]
            retcode = call(args, shell=False, env=env)
            if retcode != 0:
                print("Failure")
                break

    elif options.mode in ["a", "automatic"]:
        for script in scripts:
            print("About to execute", script)
            args = ["python", script]
            p = Popen(args, shell=False, env=env)
            time.sleep(options.time)
            p.kill()
        retcode = 0

    elif options.mode == "screenshot":
        processes = []
        for script in scripts:
            args = ["python", script]
            p = Popen(args, shell=False, env=env)
            processes.append(p)

        time.sleep(options.time)
        for p in processes:
            p.kill()
        retcode = 0

    else:
        show_examples_and_exit(err_msg="Wrong value for mode: %s" % options.mode)

    print("Final return code:", retcode)

    return retcode


if __name__ == "__main__":
    sys.exit(main())
