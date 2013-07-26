#!/usr/bin/env python
"""
This scripts runs all the python scripts located in this directory allowing
the user to change the matplotlib backend. 

Usage: 
    run.py [backend]
"""
from __future__ import print_function

import sys
import os 
import time
import shutil

from argparse import ArgumentParser
from subprocess import call, Popen

CONF_FILE = None
BKP_FILE = None

def change_backend(new_backend=""):
    """Change the backend by modifying the matplotlib configuration file."""
    global CONF_FILE, BKP_FILE

    if not new_backend:
        return 

    home = os.environ["HOME"]
    CONF_FILE = conf_file = os.path.join(home, ".matplotlib", "matplotlibrc")

    BKP_FILE = conf_file + ".bkp"

    if os.path.exists(conf_file):
        shutil.copy(conf_file, BKP_FILE)

        with open(conf_file, "r") as f:
            lines = f.readlines()

        for i, line in enumerate(lines):
            if line.strip().startswith("backend"):
                lines[i] = "backend : " + new_backend + "\n"

        with open(conf_file, "w") as f:
            f.writelines(lines)


def revert_backend():
    global CONF_FILE, BKP_FILE
    #print("reverting: BKP_FILE %s --> CONF %s" % (BKP_FILE, CONF_FILE))
    if BKP_FILE is not None:
        shutil.move(BKP_FILE, CONF_FILE)

def str_examples():
    examples = """
      Usage example:\n\n
      run.py               => Run all scripts.
      run.py -m auto -t 5  => Run all scripts, kill the process after 5 seconds.
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

    parser.add_argument('-b', '--backend', type=str, default="",
                        help="matplotlib backend")

    parser.add_argument('-m', '--mode', type=str, default="auto",
                        help="Execution mode")

    parser.add_argument('-t', '--time', type=float, default=3,
                        help="Wait time seconds before running next demo.")

    options = parser.parse_args()

    change_backend(new_backend=options.backend)

    # Find scripts.
    dir = os.path.join(os.path.dirname(__file__), "plot")
    scripts = []
    for fname in os.listdir(dir):
        if fname.endswith(".py") and fname.startswith("plot_"):
            scripts.append(os.path.join(dir, fname))

    # Run scripts according to mode.
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

    revert_backend()
    return retcode

if __name__ == "__main__":
    sys.exit(main())
