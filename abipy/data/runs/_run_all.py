#!/usr/bin/env python
"""
This script runs all the python scripts located in this directory 
"""
from __future__ import print_function, division, unicode_literals

import sys
import os 
import argparse
import shutil
import tempfile

from subprocess import call, Popen
from abipy.abilab import Flow


def main():
    def str_examples():
        examples = """
          Usage example:\n\n
          runall.py => Run all scripts.
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

    dir = os.path.join(os.path.dirname(__file__))
    scripts = []
    for fname in os.listdir(dir):
        if fname in options.exclude: continue
        if fname.endswith(".py") and fname.startswith("run_"):
            path = os.path.join(dir, fname)
            if path != __file__:
                scripts.append(path)

    # Run scripts according to mode.
    dirpaths, retcode = [], 0
    if options.mode in ["s", "sequential"]:
        for script in scripts:
            # flow will be produced in a temporary workdir.
            workdir = tempfile.mkdtemp(prefix='flow_' + os.path.basename(script))
            ret = call(["python", script, "--workdir", workdir])
            retcode += ret

            if ret != 0: 
                print("retcode %d while running %s" % (ret, script))
                if options.bail_on_failure: break

            dirpaths.append(workdir)

            execute_flow = False
            # Comment this line to execute the flow
            #execute_flow = True
            if execute_flow:
                flow = Flow.pickle_load(workdir)
                print(flow)
                flow.make_scheduler().start()

        # Remove directories.
        if not options.keep_dirs:
            for dirpath in dirpaths:
                try:
                    shutil.rmtree(dirpath, ignore_errors=False)
                except OSError:
                    print("Exception while removing %s" % dirpath)

    else:
        show_examples_and_exit(err_msg="Wrong value for mode: %s" % options.mode)

    print("retcode %d" % retcode)
    return retcode

if __name__ == "__main__":
    sys.exit(main())
