#!/usr/bin/env python

import sys
import os
import time

from argparse import ArgumentParser
from subprocess import call, Popen

def str_examples():
    examples = """
      Usage example:\n\n
      run.py               => Run all demo scripts.
      run.py -m auto -t 5  => Run all tests, close the demo after 5 seconds.
    """
    return examples

def show_examples_and_exit(err_msg=None, error_code=1):
    """Display the usage of the script."""
    sys.stderr.write(str_examples())
    if err_msg: 
        sys.stderr.write("Fatal Error\n" + err_msg + "\n")
    sys.exit(error_code)


def main():
    parser = ArgumentParser(epilog=str_examples())

    parser.add_argument('-m', '--mode', type=str, default="sequential",
                        help="Execution mode")

    parser.add_argument('-t', '--time', type=float, default=3,
                        help="Wait time seconds before running next demo.")

    options = parser.parse_args()

    # Find scripts.
    scripts = []
    for fname in os.listdir("."):
        if fname.endswith(".py") and fname.startswith("demo_"):
            scripts.append(fname)

    # Run scripts depending on mode.
    if options.mode == "sequential":
        for script in scripts:
            retcode = call(["python", script])
            if retcode != 0: break

    elif options.mode == "auto":
        for script in scripts:
            p = Popen(["python", script])
            time.sleep(options.time)
            p.kill()
        retcode = 0

    elif options.mode == "screenshot":
        processes = []
        for script in scripts:
            p = Popen(["python", script])
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
