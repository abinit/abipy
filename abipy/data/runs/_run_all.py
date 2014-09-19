#!/usr/bin/env python
"""
This script runs all the python scripts located in this directory 
"""
from __future__ import print_function

import sys
import os 
import argparse

from subprocess import call, Popen

def str_examples():
    examples = """
      Usage example:\n\n
      runall.py               => Run all scripts.
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

    parser.add_argument('-m', '--mode', type=str, default="sequential",
                        help="execution mode. Default is sequential.")

    parser.add_argument('-e', '--exclude', type=str, default="",
                        help="Exclude scripts.")

    options = parser.parse_args()

    # Find scripts.
    if options.exclude:
        options.exclude = options.exclude.split()
        print("Will exclude:\n", options.exclude)

    dir = os.path.join(os.path.dirname(__file__))
    scripts = []
    for fname in os.listdir(dir):
        if fname in options.exclude: continue
        if fname.endswith(".py") and fname.startswith("run_"):
            path = os.path.join(dir, fname)
            if path != __file__:
                scripts.append(path)

    # Run scripts according to mode.
    if options.mode in ["s", "sequential"]:
        for script in scripts:
            retcode = call(["python", script])
            if retcode != 0: break

    #elif options.mode in ["a", "automatic"]:
    #    for script in scripts:
    #        p = Popen(["python", script])
    #        time.sleep(options.time)
    #        p.kill()
    #    retcode = 0

    else:
        show_examples_and_exit(err_msg="Wrong value for mode: %s" % options.mode)

    return retcode

if __name__ == "__main__":
    sys.exit(main())
