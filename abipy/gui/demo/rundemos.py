#!/usr/bin/env python
from __future__ import print_function

import sys
import os
import time
import argparse

from subprocess import call, Popen

def str_examples():
    examples = (
      "Usage example:\n"
      "run.py               => Run all demo scripts.\n"
      "run.py -m automatic -t 7  => Run all tests, close the demo after 7 seconds.\n"
      )
    return examples


def show_examples_and_exit(err_msg=None, error_code=1):
    """Display the usage of the script."""
    sys.stderr.write(str_examples())
    if err_msg: 
        sys.stderr.write("Fatal Error\n" + err_msg + "\n")
    sys.exit(error_code)


def main():
    parser = argparse.ArgumentParser(epilog=str_examples(),formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-m', '--mode', type=str, default="automatic",
                        help="execution mode. Either s (sequential) or a (automatic)")

    parser.add_argument('-t', '--time', type=float, default=5,
                        help="wait time seconds before running next demo.")

    options = parser.parse_args()

    # Find scripts.
    dirpath = os.path.dirname(__file__)
    scripts = []
    for fname in os.listdir(dirpath):
        if fname.endswith(".py") and fname.startswith("demo_"):
            scripts.append(os.path.join(dirpath, fname))

    python = "pythonw"
    python = "python.app"
    print("Using python:", python)

    # Run scripts depending on mode.
    if options.mode in ["s", "sequential"]:
        for script in scripts:
            retcode = call([python, script])
            if retcode != 0: break

    elif options.mode in ["a", "automatic"]:
        for script in scripts:
            p = Popen([python, script])
            time.sleep(options.time)
            p.kill()
        retcode = 0

    elif options.mode == "screenshot":
        processes = []
        for script in scripts:
            p = Popen([python, script])
            processes.append(p)
            
        time.sleep(options.time)
        for p in processes:
            p.kill()
        retcode = 0

    else:
        show_examples_and_exit(err_msg="Wrong value for mode: %s" % options.mode)

    return retcode


if __name__ == "__main__":
    sys.exit(main())
