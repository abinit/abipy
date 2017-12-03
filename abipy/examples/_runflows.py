#!/usr/bin/env python
"""
This script runs all the python scripts located in this directory 
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os 
import argparse
import shutil
import tempfile

from subprocess import call, Popen
from abipy import __version__

#root = os.path.abspath(os.path.join(os.path.dirname(__file__)), "flows")


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

    parser.add_argument('-V', '--version', action='version', version="%(prog)s version " + __version__)
    parser.add_argument('--loglevel', default="ERROR", type=str,
                        help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    parser.add_argument('-m', '--mode', type=str, default="sequential", help="execution mode. Default is sequential.")

    parser.add_argument('-e', '--exclude', type=str, default="", help="Exclude scripts. Comma-separated names")

    parser.add_argument('-x', '--execute', default=False, action="store_true", help="Execute flows.")

    parser.add_argument('--keep-dirs', action="store_true", default=False,
                        help="Do not remove flowdirectories.")

    parser.add_argument('-b', '--bail-on-failure', default=False, help="Exit at the first error.")

    #parser.add_argument("scripts", nargs="+",help="List of scripts to be executed")

    options = parser.parse_args()

    # loglevel is bound to the string value obtained from the command line argument. 
    # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
    import logging
    numeric_level = getattr(logging, options.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.loglevel)
    logging.basicConfig(level=numeric_level)

    # Find scripts.
    if options.exclude:
        options.exclude = options.exclude.split(",")
        print("Will exclude:\n", options.exclude)

    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "flows"))
    scripts = []
    for fname in os.listdir(root):
        if fname in options.exclude: continue
        if fname.endswith(".py") and fname.startswith("run_"):
            path = os.path.join(root, fname)
            if path != __file__:
                scripts.append(path)

    # Run scripts according to mode.
    dirpaths, errors, retcode, cnt = [], [], 0, 0
    if options.mode in ["s", "sequential"]:
        for script in scripts:
            # flow will be produced in a temporary workdir.
            workdir = tempfile.mkdtemp(prefix='flow_' + os.path.basename(script))
            ret = call(["python", script, "--workdir", workdir])
            retcode += ret

            if ret != 0: 
                e = "python %s returned retcode !=0" % script
                print(e)
                errors.append(e)
                if options.bail_on_failure: 
                    print("Exiting now since bail_on_failure")
                    break

            dirpaths.append(workdir)

            # Here we execute the flow
            if options.execute:
                cnt += 1
                ret = 0
                try:
                    flow = Flow.pickle_load(workdir)
                    flow.make_scheduler().start()
                    if not flow.all_ok: retcode += 1

                except Exception as exc:
                    ret += 1
                    s = "Exception raised during flow execution: %s\n:%s" % (flow, exc)
                    print(s)
                    errors.append(s)
                    if options.bail_on_failure: 
                        print("Exiting now since bail_on_failure")
                        break
                    retcode += ret

        # Remove directories.
        if not options.keep_dirs:
            for dirpath in dirpaths:
                try:
                    shutil.rmtree(dirpath, ignore_errors=False)
                except OSError:
                    print("Exception while removing %s" % dirpath)

    else:
        show_examples_and_exit(err_msg="Wrong value for mode: %s" % options.mode)

    if errors:
        for i, err in enumerate(errors):
            print(92 * "=")
            print("[%d] %s" % (i, err))
            print(92 * "=")

    print("Number of scripts executed %d/%d" % (cnt, len(scripts)), ", final retcode %d" % retcode)
    return retcode


if __name__ == "__main__":
    sys.exit(main())
