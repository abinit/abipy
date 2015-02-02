#!/usr/bin/env python
"""
This script regenerate all the reference files located in this directory 
"""
from __future__ import print_function, division, unicode_literals

import sys
import os 
import argparse
import shutil
import tempfile

from subprocess import call, Popen


def main():
    def str_examples():
        examples = """
          Usage example:\n\n
          runall.py               => Regenerate all reference files.
        """
        return examples

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg: 
            sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    parser = argparse.ArgumentParser(epilog=str_examples(),formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-m', '--mode', type=str, default="sequential",
                        help="execution mode. Default is sequential.")

    parser.add_argument('-e', '--exclude', type=str, default="", help="Exclude scripts.")

    parser.add_argument('--keep-dirs', action="store_true", default=False,
                        help="Do not remove flowdirectories.")

    parser.add_argument('-b', '--bail-on-failure', default=False, help="Exit at the first error.")

    #parser.add_argument("scripts", nargs="+",help="List of scripts to be executed")

    options = parser.parse_args()

    # Find scripts.
    if options.exclude:
        options.exclude = options.exclude.split()
        print("Will exclude:\n", options.exclude)

    scripts = []
    for dirpath, dirnames, filenames in os.walk(os.path.dirname(__file__)):
        for fname in filenames:
            if fname in options.exclude: continue
            if fname == "gendata.py":
                path = os.path.join(os.path.abspath(dirpath), fname)
                #if path != __file__:
                scripts.append(path)

    #print(scripts)
    #return 0

    # Run scripts according to mode.
    dirpaths, retcode = [], 0
    if options.mode in ["s", "sequential"]:
        for script in scripts:
            os.chdir(os.path.dirname(script))
            ret = call(["python", script])
            retcode += ret

            if ret != 0: 
                print("retcode %d while running %s" % (ret, script))
                if options.bail_on_failure: break

        # Remove directories.
        #if not options.keep_dirs:
        #    for dirpath in dirpaths:
        #        try:
        #            shutil.rmtree(dirpath, ignore_errors=False)
        #        except OSError:
        #            print("Exception while removing %s" % dirpath)

    else:
        show_examples_and_exit(err_msg="Wrong value for mode: %s" % options.mode)

    print("retcode %d" % retcode)
    return retcode

if __name__ == "__main__":
    sys.exit(main())
